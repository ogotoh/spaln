/*****************************************************************************
*
*	fwd2b1.c
*
*	Alignment of two protein or nucleotide sequences.
*	Splicing is NOT considered.
*	Assumes there is no internal gap in either sequence.
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#define	DEBUG	1

#include "aln.h"
#include "vmf.h"
#include "wln.h"
#if __SSE4_1__
#include "fwd2s1_simd.h"
#endif

class Aln2b1 {
const	Seq**	seqs;
const	Seq*&	a;
const	Seq*&	b;
const	PwdB*	pwd;
const	bool	Local;
	Mfile*	mfd = 0;
	Vmf*	vmf = 0;
const	INT	lowestlvl;
public:
	Aln2b1(const Seq** _seqs, const PwdB* _pwd) :
	    seqs(_seqs), a(seqs[0]), b(seqs[1]), pwd(_pwd), 
	    Local(algmode.lcl & 16), lowestlvl(b->wllvl) {} 
	~Aln2b1() {delete vmf; delete mfd;}
	void	initB_ng(RVPD* hh[], const WINDOW& wdw);
	RVPD*	lastB_ng(RVPD* hh[], const WINDOW& wdw);
	VTYPE	forwardB_ng(const WINDOW& wdw, int pp[] = 0);
	void	hinitB_ng(Rvwml* hh[], const WINDOW& wdw);
	Rvwml*	hlastB_ng(Rvwml* hh[], const WINDOW& wdw);
	VTYPE	hirschbergB_ng(int cpos[], 		// coordinates on center row
		const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb, INT& bexgl);
	void	sinitB_ng(VTYPE* hh[], const WINDOW& wdw);
	VTYPE	slastB_ng(VTYPE* hh[], const WINDOW& wdw);
	VTYPE	scorealoneB_ng(const WINDOW& wdw);
	VTYPE	diagonalB_ng();
	VTYPE	trcbkalignB_ng(const WINDOW& wdw);
	void	cinitB_ng(RVDWC* hhc[], const WINDOW& wdw);
	Colonies*	fwdswgB_ng(VTYPE* scr, const WINDOW& wdw);
	VTYPE	backforth(int n, const BOUND& lub);
	VTYPE	lspB_ng(const WINDOW& wdw);
	VTYPE	seededB_ng(INT level, int cmode, const BOUND& lub);
	SKL*	globalB_ng(VTYPE* scr, const WINDOW& wdw);
	VTYPE	creepback(int ovr, VTYPE bscr, const BOUND& lub);
	VTYPE	creepfwrd(int& ovr, VTYPE bscr, const BOUND& lub);
};

void Aln2b1::initB_ng(RVPD* hh[], const WINDOW& wdw)
{
	int	r = b->left - a->left;
	int	rr = b->right - a->left;
	int	origin = vmf->add(a->left, b->left, 0);
	FTYPE	ltgapf = a->left? 1: (a->inex.exgl? 0: alprm.tgapf);

	RVPD*	h = hh[0] + r;
	RVPD*	g = hh[1] + r;
	h->val = 0;
	h->dir = NEWD;
	h->ptr = origin;
	if (wdw.up < rr) rr = wdw.up;
	for (int i = 1; ++r <= rr; ++i) {
	    ++h;
	    h->dir = HORI;
	    h->ptr = origin;
const	    VTYPE	gpn = (i == 1)? pwd->GapPenalty(1): pwd->GapExtPen(i);
	    h->val = h[-1].val + (VTYPE) (gpn * ltgapf);
	}

	ltgapf = b->left? 1: (b->inex.exgl? 0: alprm.tgapf);
	r = b->left - a->left;
	rr = b->left - a->right;
	h = hh[0] + r;
	if (wdw.lw > rr) rr = wdw.lw;
	for (int i = 1; --r >= rr; ++i) {
	    --h;
	    h->dir = VERT;
	    h->ptr = origin;
	    VTYPE	gpn = (i == 1)? pwd->GapPenalty(1): pwd->GapExtPen(i);
	    h->val = h[1].val  + (VTYPE) (gpn * ltgapf);
	    *--g = *h;
	}
}

RVPD* Aln2b1::lastB_ng(RVPD* hh[], const WINDOW& wdw)
{
	RVPD*	h9 = hh[0] + b->right - a->right;
	FTYPE	rtgapf = b->inex.exgr? 0: alprm.tgapf;

	int	dm = 0;
	if (b->right == b->len && rtgapf < 1) {
	    int	rw = wdw.up;
	    int	rf = b->right - a->left;
	    if (rf < rw) rw = rf;
	    RVPD*	h = hh[0] + rw;
	    while (--h >= h9) {
		RVPD*	g = h + 1;
		++dm;
		VTYPE	gpn = isntvert(g)? pwd->GapPenalty(1): pwd->GapExtPen(dm);
		g->val += VTYPE (gpn * rtgapf);
		if (g->val > h->val) {*h = *g; h->dir = VERT;}
		else	dm = 0;
	    }
	}
	int	dn = 0;
	rtgapf = a->inex.exgr? 0: alprm.tgapf;
	if (a->right == a->len && rtgapf < 1) {
	    int	rw = wdw.lw;
	    int	rf = b->left - a->right;
	    if (rf > rw) rw = rf;
	    else	 rf = rw;
	    RVPD*	h = hh[0] + rw;
	    while (++h <= h9) {
		RVPD*	f = h - 1;
		++dn;
		VTYPE	gpn = isnthori(f)? pwd->GapPenalty(1): pwd->GapExtPen(dn);
		f->val += (VTYPE) (gpn * rtgapf);
		if (f->val > h->val) {*h = *f; h->dir = VERT;}
		else	dn = 0;
	    }
	}
	if (dn || dm) {
	    if (dn) dm = 0;
	    h9->ptr = vmf->add(a->right - dm, b->right - dn, h9->ptr);
	}
	h9->ptr = vmf->add(a->right, b->right, h9->ptr);
	return (h9);
}

VTYPE Aln2b1::forwardB_ng(const WINDOW& wdw, int pp[])
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
const	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
	RVPD*	hh[NOL];
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};

	RVPD*	rbuf = new RVPD[pwd->Noll * wdw.width];
	vset(rbuf, black_vpd, pwd->Noll * wdw.width);
	hh[0] = rbuf - wdw.lw + 1;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + wdw.width;
	vmf->add(0, 0, 0);	// Skip 0-th record
	initB_ng(hh, wdw);

	int	m = a->left;
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
const	CHAR*	as = a->at(m);

	for ( ; ++m <= a->right; ++as, ++n1, ++n2) {
	    RVPD	e1 = black_vpd;
	    RVPD	e2 = black_vpd;
	    int		n = max(n1, b->left);
const	    int		n9 = min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int		r = n - m;
	    int		k = 0;
	    RVPD*	h = hh[k] + r;
	    RVPD*	f = hh[++k] + r;
	    RVPD*	f2 = dagp? hh[++k] + r: 0;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; ++n <= n9; ++bs) {
		++h; ++f; if (dagp) ++f2;
//	Diagonal
		RVPD*	from = h;
		RVPD*	mx = h;
		VTYPE	diag = h->val;
		h->val += qprof[*bs];
		h->dir = isdiag(from)? DIAG: NEWD;

//	Vertical
		VTYPE	x = (++from)->val + pwd->BasicGOP;
		if (x >= f[1].val) {
		    f->val = x;
		    f->ptr = from->ptr;
		    f->dir = VERT;
		} else	*f = f[1];
		f->val += pwd->BasicGEP;
		if (f->val > mx->val) mx = f;

//	Vertical2
		if (dagp) {
		  x = from->val + pwd->LongGOP;
		  if (x >= f2[1].val) {
		    f2->val = x;
		    f2->ptr = from->ptr;
		    f2->dir = VERT;
		  } else	*f2 = f2[1];
		  f2->val += pwd->LongGEP;
		  if (f2->val > mx->val) mx = f2;
		}
//	Horizontal
		x = h[-1].val + pwd->BasicGOP;
		if (x >= e1.val) {
		    e1.val = x;
		    e1.ptr = h[-1].ptr;
		    e1.dir = HORI;
		}
		e1.val += pwd->BasicGEP;
		if (e1.val >= mx->val) mx = &e1;

//	Horizontal2
		if (dagp) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= e2.val) {
			e2.val = x;
			e2.ptr = h[-1].ptr;
			e2.dir = HORL;
		    }
		    e2.val += pwd->LongGEP;
		    if (e2.val >= mx->val) mx = &e2;
		}

//	Find optimal path
		if (h != mx) *h = *mx;	// non-diagonal
		else if (Local && h->val > diag) {
		    if (LocalL && diag == 0)
			h->ptr = vmf->add(m - 1, n - 1, 0);
		    else if (LocalR && h->val > maxh.val) {
			maxh.val = h->val;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && h->val <= 0) h->val = h->dir = 0;
		else if (h->dir == NEWD || h->dir == NEWV || h->dir == NEWH)
			h->ptr = vmf->add(m - 1, n - 1, h->ptr);

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d %2d ", m, n, mx->dir);
		    putvar(mx->val); putvar(diag); 
		    putvar(f->val); putvar(e1.val);
		    if (dagp) {
			putvar(f2->val); putvar(e2.val);
		    }
		    putchar('\n');
		}
#endif
	    }	/* end of n-loop */
	}	/* end of m-loop */

	if (LocalR) {
	    *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    RVPD*	mx = lastB_ng(hh, wdw);
	    maxh.val = mx->val;
	    pp[0] = mx->ptr;
	    pp[1] = b->left - a->left + mx - hh[0];
	}
	delete[] rbuf;
	return (maxh.val);
}

VTYPE skl_rngB_ng(const Seq* seqs[], Gsinfo* gsi, const PwdB* pwd)
{
const	Seq*&	a = seqs[0];
const	Seq*&	b = seqs[1];
	VTYPE	scr = 0;
const	SKL*	wsk = trimskl(seqs, gsi->skl);
	int	num = (wsk++)->n;
	FSTAT*	fst = &gsi->fstat;
	Cigar*	cigar = gsi->cigar = 0;
	Vulgar*	vlgar = gsi->vlgar = 0;
	Samfmt*	samfm = gsi->samfm = 0;
	switch (algmode.nsa) {
	    case CIG_FORM: cigar = gsi->cigar = new Cigar(); break;
	    case VLG_FORM: gsi->vlgar = vlgar = new Vulgar(); break;
	    case SAM_FORM: gsi->samfm = samfm = new Samfmt(); break;
	    default:	 break;
	}
	vclear(fst);
	int	m = wsk->m;
	int	n = wsk->n;
const	CHAR*	as = a->at(m);
const	CHAR*	bs = b->at(n);
	if (samfm) {
	    samfm->left = m;
	    if ( m) samfm->push('H', m);	// local alignment
	    if (b->inex.sens == 0) samfm->pos = n;
	}
	FTYPE	tgapf = (m == 0. || n == 0)? alprm.tgapf: 1.;
	Iiinfo*	iif = (SpbFact > 0 && ((a->sigII && a->sigII->pfqnum) || (b->sigII && b->sigII->pfqnum)))?
	    new Iiinfo(seqs, a->left, b->left, true): 0;
	int	span = 0;
	while (--num) {
	    ++wsk;
	    int	mi = wsk->m - m;
	    int	ni = wsk->n - n;
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    span += max(mi, ni);
	    if (d) {
		if (cigar) cigar->push('M', d);
		if (samfm) samfm->push('M', d);
		if (vlgar) vlgar->push('M', d, d);
		n += d;
		m += d;
		if (iif) scr += iif->StoreIIinfo(m, n);
		for ( ; d--; ++as, ++bs) {
		    scr += pwd->sim2(as, bs);
		    if (*as == *bs) ++fst->mch;
		    else	++fst->mmc;
		}
	    }
	    if (i < 0) {
		d = -i;
		n -= i;
		bs += d;
		if (cigar) cigar->push('D', d);
		if (samfm) samfm->push('D', d);
		if (vlgar) vlgar->push('G', 0, d);
	    } else if (i > 0) {
		d = i;
		m += i;
		as += d;
		if (cigar) cigar->push('I', d);
		if (samfm) samfm->push('I', d);
		if (vlgar) vlgar->push('G', d, 0);
	    } else d = 0;
	    if (d) {
		if (wsk->m == a->len || wsk->n == b->len) tgapf = alprm.tgapf;
		fst->gap += tgapf;
		fst->unp += d * tgapf;
		scr +=  (VTYPE) (pwd->GapPenalty(d) * tgapf);
		tgapf = 1.;
		if (iif) {
		    d *= iif->step;
		    if (i > 0)  iif->bgap += d;
		    else	iif->agap += d;
		    scr += iif->StoreIIinfo(m, n);
		}
	    }
	    m = wsk->m;
	    n = wsk->n;
	}
	if (samfm) {
	    samfm->right = m;
	    if (m < a->len) samfm->push('H', a->len - m);
	    samfm->flush();
	    samfm->mapq = 30 + int(100 * (fst->mmc + fst->unp) / a->len);
	    if (b->inex.sens) {
		samfm->flag |= 0x10;
		samfm->pos = n - 1;
		vreverse(samfm->rec, samfm->size());
	    }
	}
	if (cigar) cigar->flush();
	if (vlgar) vlgar->flush();
	fst->val = (FTYPE) scr;
	if (gsi) gsi->SaveGsInfo(iif, span);
	delete iif;
	return (scr);
}

/********************************************************
*
*	unidirectional Hirschberg algorithm
*
********************************************************/

void Aln2b1::hinitB_ng(Rvwml* hh[], const WINDOW& wdw)
{
	int	r = b->left - a->left;
	int	r0 = r;
	int	rr = b->right - a->left;

	Rvwml*	h = hh[0]  + r;
	h->val = 0;
	h->lwr = h->upr = h->ulk = r;
	h->ml = a->left;
	if (a->inex.exgl) {	// semi-global
	    if (wdw.up < rr) rr = wdw.up;
	    while (++r <= rr) {
		++h;
		h->val = 0;
		h->lwr = h->upr = h->ulk = r;
	    }
	}

	r = r0;
	rr = b->left - a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	h = hh[0] + r;
	Rvwml*	f = hh[1] + r;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --f;
	    if (b->inex.exgl) {
		h->val = 0;
		h->lwr = h->upr = h->ulk = r;
	    } else {
		*h = h[1];
		if (i == 1) {
		    h->val += pwd->GapPenalty(1);
		    *f = *h;
		} else {
		    *f = f[1];
		    h->val += pwd->GapExtPen(i);
		    f->val += pwd->BasicGEP;
		}
		h->lwr = f->lwr = r;
		h->ulk = f->ulk = r0;
	    }
	}
}

Rvwml* Aln2b1::hlastB_ng(Rvwml* hh[], const WINDOW& wdw)
{
	Rvwml*	h9 = hh[0] + b->right - a->right;
	Rvwml*	mx = h9;

	if (b->inex.exgr) {
	    int	rw = std::min(wdw.up, b->right - a->left);
	    for (Rvwml* h = hh[0] + rw; h > h9; --h)
		if (h->val > mx->val) mx = h;
	}
	if (a->inex.exgr) {
	    int	rw = std::max(wdw.lw, b->left - a->right);
	    for (Rvwml* h = hh[0] + rw; h < h9; ++h)
		if (h->val > mx->val) mx = h;
	}
	return (mx);
}

VTYPE Aln2b1::hirschbergB_ng(int cpos[], 
	const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb, INT& bexgl)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
const	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
	Rvwml*	hh[NOL];
	Rvwml*	hf[NOD];
	Rvwml	e[2];
	Rvwml*&	e1 = hf[1] = e;
	Rvwml*&	e2 = hf[3] = e + 1;
	int*	center_lnk[NOL];
	int*	center_end[NOL];
	int*	center_lwr[NOL];
	int*	center_upr[NOL];
	Rvwmrmn	maxh = {NEVSEL, 0, 0, 0, 0};

const	size_t	bufsiz = pwd->Noll * wdw.width;
	Rvwml*	wbuf = new Rvwml[bufsiz];
	int	r = b->left - a->right;
const	Rvwml	black_vpwml = {NEVSEL, r, r, 0, end_of_ulk};
	vset(wbuf, black_vpwml, bufsiz);
	Rvwml*	blackvwu = wbuf + bufsiz - 1;	// assume to be const
	hh[0] = wbuf - wdw.lw + 1;
	int*	ibuf = new int[4 * bufsiz];
	center_lnk[0] = ibuf - wdw.lw + 1;
	center_end[0] = center_lnk[0] + bufsiz;
	center_lwr[0] = center_end[0] + bufsiz;
	center_upr[0] = center_lwr[0] + bufsiz;
	for (int k = 1; k < pwd->Noll; ++k) {
	    hh[k] = hh[k-1] + wdw.width;
	    center_lnk[k] = center_lnk[k - 1] + wdw.width;
	    center_end[k] = center_end[k - 1] + wdw.width;
	    center_lwr[k] = center_lwr[k - 1] + wdw.width;
	    center_upr[k] = center_upr[k - 1] + wdw.width;
	}
	vset(ibuf, end_of_ulk, 2 * bufsiz);	// center_lnk

	hinitB_ng(hh, wdw);

	int	m = a->left;
const	int	mm = (a->left + a->right + 1) / 2;
	if (!a->inex.exgl) --m; 	// global
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as, ++n1, ++n2) {
	    int	n = std::max(n1, b->left);
const	    int	n9 = std::min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    r = n - m;
	    Rvwml*&	h = hf[0] = hh[0] + r;
	    Rvwml*&	f = hf[2] = hh[1] + r;
	    Rvwml*&	f2 = hf[4] = dagp? hh[2] + r: blackvwu;
	    vset(e, black_vpwml, 2);
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d ", m, n);
		putvar(h->val); putchar('\n');
	    }
#endif
	    for ( ; ++n <= n9; ++bs) {
		++r;
		for (int k = 0; k < pwd->Noll; ++k) {
		    int	kk = k + k;
		    ++hf[kk];
		    if (m == mm) {
			center_end[k][r] = hf[kk]->ulk;
			hf[kk]->ulk = k? end_of_ulk: r;
		    }
		}

//	Diagonal
		Rvwml*	from = h;
		Rvwml*	mx = h;
		VTYPE	x;
		if (m == a->left) goto HorizonF;
		h->val += qprof[*bs];

//	Vertical
		x = (++from)->val + pwd->BasicGOP;
		if (x >= f[1].val) {
		    *f = *from;
		    f->val = x;
		} else	*f = f[1];
		f->val += pwd->BasicGEP;
		if (f->val >= mx->val) mx = f;

//	Vertical2
		if (dagp) {
		  x = from->val + pwd->LongGOP;
		  if (x >= f2[1].val) {
		    *f2 = *from;
		    f2->val = x;
		  } else 	*f2 = f2[1];
		  f2->val += pwd->LongGEP;
		  if (f2->val >= mx->val) mx = f2;
		}
HorizonF:
//	Horizontal
		x = h[-1].val + pwd->BasicGOP;
		if (x >= e1->val) {
		    *e1 = h[-1];
		    e1->val = x;
		}
		e1->val += pwd->BasicGEP;
		if (e1->val >= mx->val) mx = e1;

//	Horizontal2
		if (dagp) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= e2->val) {
			*e2 = h[-1];
			e2->val = x;
		    }
		    e2->val += pwd->LongGEP;
		    if (e2->val >= mx->val) mx = e2;
		}

//	Find optimal path
#if DEBUG
		VTYPE	y = h->val;
#endif
		int	hd = 0;		// diag: 0, hori: 1, vert: 2
		if (h == mx) {
		    if (LocalR && h->val > maxh.val) {
			maxh.val = h->val;
			maxh.upr = h->upr;
			maxh.lwr = h->lwr;
			maxh.ml = h->ml;
			maxh.ulk = h->ulk;
			maxh.mr = m;
			maxh.nr = n;
		    }
		} else {
		    while (mx != hf[++hd]) ;
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		}
		if (LocalL && h->val <= 0) {
		    h->val = 0;
		    h->ml = m;
		    h->ulk = h->upr = h->lwr = r;
		}

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d ", m, n);
		    putvar(mx->val); putvar(y); 
		    putvar(f->val); putvar(e1->val);
		    if (dagp) {
			putvar(f2->val); putvar(e2->val);
		    }
		    putchar('\n');
		}
#endif

// save variables at the center
		if (m == mm) {
const		    int	vk = (hd == 2 || hd == 4)? hd / 2: 0;
		    for (int k = 0; k < pwd->Noll; ++k) {
			int	kk = k + k;
			center_lwr[k][r] = hf[kk]->lwr;
		        center_upr[k][r] = hf[kk]->upr;
			hf[kk]->upr = r;
			if (!(hd % 2)) hf[kk]->lwr = r;
			if (k) hf[kk]->ulk = r + k * wdw.width;
			else {
			    center_lnk[0][r] = hf[0]->ulk;
			    hf[0]->ulk = r + vk * wdw.width;
			}
		    }
		}	// was center
	    }	// end of n-loop
	}	// end of m-loop

const	int	rl = b->left - a->left;
	if (LocalR) {
	    a->right = maxh.mr;
	    b->right = maxh.nr;
	    wdwb.up = maxh.upr;
	    wdwb.lw = maxh.lwr;
	} else {
const	    int	rr = b->right - a->right;
const	    Rvwml*	mx = hlastB_ng(hh, wdw);
	    maxh.val = mx->val;
	    maxh.ulk = mx->ulk;
	    wdwb.up = mx->upr;
	    wdwb.lw = mx->lwr;
	    r = mx - hh[0];
	    if (a->inex.exgr && rr < r) a->right = b->right - r;
	    if (b->inex.exgr && rr > r) b->right = a->right + r;
	}
	int	c = 0;
	int	d = 0;
	for (r = maxh.ulk; r > wdw.up; r -= wdw.width) ++d;
	if (center_end[d][r] < end_of_ulk) {	// cross center
	    cpos[c++] = mm;
	    for (int rp = center_lnk[d][r]; rp < end_of_ulk && r != rp;
		rp = center_lnk[0][r = rp]) {
		cpos[c++] = r + mm;
	    }
	    cpos[c++] = r + mm;
	    cpos[c] = end_of_ulk;
	    wdwf.up = center_upr[d][r];
	    wdwf.lw = center_lwr[0][r];
	    r = center_end[d][r];
	} else	cpos[0] = end_of_ulk;		// don't cross center
	if (LocalL) {
	    a->left = maxh.ml;
	    b->left = r + maxh.ml;
	} else {
	    if (a->inex.exgl && rl > r) a->left = b->left - r;
	    if (b->inex.exgl && rl < r) b->left = a->left + r;
	}
	bexgl = (d > 0)? 1: 0;
	delete[] wbuf;
	delete[] ibuf;
	return (maxh.val);
}

void Aln2b1::cinitB_ng(RVDWC* hhc[], const WINDOW& wdw)
{
	int	n = b->left;
	int	r = b->left - a->left;
	int	rr = b->right - a->left;
	RVDWC*	h = hhc[0] + r;

	if (wdw.up < rr) rr = wdw.up;
	for ( ; r <= rr; ++h, ++r) {
	    h->val = 0;
	    h->dir = DEAD;
	    h->upr = h->lwr = r;
	    h->mlb = a->left;
	    h->nlb = n++;
	}

	r = b->left - a->left;
	rr = b->left - a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	h = hhc[0] + r;
	while (--r >= rr) {
	    --h;
	    h->val = 0;
	    h->dir = DEAD;
	    h->upr = h->lwr = r;
	    h->mlb = b->left - r;
	    h->nlb = b->left;
	}
}

Colonies* Aln2b1::fwdswgB_ng(VTYPE* scr, const WINDOW& wdw)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	RVDWC*	hhc[NOL];
	RVDWC   hl[NOL];
	RVDWC	e[NOL]; 		// [DIAG, HORI, HORL]
	RVDWC*	f2 = 0;
	Colonies*	cl = new Colonies(0);
	COLONY*	clny = cl->at();

	RVDWC*	cbuf = new RVDWC[pwd->Noll * wdw.width];
	vset(cbuf, black_vdwc, pwd->Noll * wdw.width);
        hhc[0] = cbuf - wdw.lw + 1;
        for (int k = 1; k < pwd->Noll; ++k) hhc[k] = hhc[k-1] + wdw.width;
	cinitB_ng(hhc, wdw);
	int	m  = a->left;
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as, ++n1, ++n2) {
	    int 	n = max(n1, b->left);
	    int 	n9 = min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int 	r = n - m;
	    RVDWC*	hlb = hhc[0] + r;
	    RVDWC*	hrb = hhc[0] + n9 - m;
	    RVDWC*	h = hlb;
	    RVDWC*	f = hhc[1] + r;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    if (dagp) f2 = hhc[2] + r;
	    for (int k = 0; k < pwd->Noll; ++k)
		hl[k] = e[k] = black_vdwc;

#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; ++n <= n9; ++bs) {
		++r; ++h; ++f; if (f2) ++f2;
//	Diagonal
		RVDWC*	mx = h;
		e[0] = *h;	/* diag */
#if DEBUG
		VTYPE	diag = h->val;
#endif
		h->val += qprof[*bs];
		h->dir = isdiag(h)? DIAG: NEWD;

//	Vertical
		VTYPE	x = h[1].val + pwd->BasicGOP;
		if (x >= f[1].val) {
		    *f = h[1];
		    f->val = x;
		    f->dir = VERT;
		} else	*f = f[1];
		f->val += pwd->BasicGEP;
		if (f->val > mx->val) mx = f;

//	Vertical2
		if (dagp) {
		  x = h[1].val + pwd->LongGOP;
		  if (x >= f2[1].val) {
		    *f2 = h[1];
		    f2->val = x;
		    f2->dir = VERT;
		  } else 	*f2 = f2[1];
		  f2->val += pwd->LongGEP;
		  if (f2->val > mx->val) mx = f2;
		}

//	Horizontal
		x = h[-1].val + pwd->BasicGOP;
		if (x >= e[1].val) {
		    e[1] = h[-1];
		    e[1].dir = HORI;
		    e[1].val = x;
		}
		if (e[1].dir) {
		    e[1].val += pwd->BasicGEP;
		    e[1].dir = HORI;
		    if (e[1].val >= mx->val) mx = e + 1;
		}

//	Horizontal2
		if (dagp) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= e[2].val) {
			e[2] = h[-1];
			e[2].val = x;
			e[2].dir = HORL;
		    }
		    if (e[2].dir) {
			e[2].val += pwd->LongGEP;
			e[2].dir = HORL;
			if (e[2].val >= mx->val) mx = e + 2;
		    }
		}

//	Find optimal path
		if (h != mx) {			// non-diagonal
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		} else if (h->val > e->val) {
		    if (e->val == 0) {		// new colony
			h->upr = h->lwr = r;
			h->mlb = m - 1;
			h->nlb = n - 1;
		    }
		    if (h->val > clny->val) {	// max local score
			clny->val = h->val;
			clny->mrb = m;
			clny->nrb = n;
			clny->lwr = h->lwr;
			clny->upr = h->upr;
			clny->mlb = h->mlb;
			clny->nlb = h->nlb;
		    }
		}
		if (h->val < 0) {		// reset to blank
		    *h = e[1] = *f = white_vdwc;
		    if (f2) {e[2] = *f2 = white_vdwc;}
		    h->clny = 0;
		}
		if (algmode.mlt > 1 && h->val >= pwd->Vthr && !h->clny) {
		    if (cl->size() >= (int) OutPrm.MaxOut) {	// Gabage collecton
			for (RVDWC* hw = hlb; hw < hrb; ++hw)	// mark active colony
			    if (hw->clny) hw->clny->mark = 1;
			if (algmode.mlt == 2) cl->removeoverlap();
			if (cl->size() >= (int) OutPrm.MaxOut) cl->removelowscore();
			for (RVDWC* hw = hlb; hw < hrb; ++hw)
			    if (hw->clny) hw->clny = cl->at(hw->clny->clno);
			for (COLONY* cc = cl->at(1); cc <= cl->end(); ++cc) {
			    cc->mark = 0;
			    cc->clno = cc - clny;
			}
		    }
		    h->clny = cl->next();
		}
		COLONY*	cc = h->clny;
		if (cc) {
		    if (h->val > cc->val) {
			cc->val = h->val;
			cc->mrb = m;
			cc->nrb = n;
			cc->lwr = h->lwr;
			cc->upr = h->upr;
			cc->mlb = h->mlb;
			cc->nlb = h->nlb;
		    } else if (h->val <= cc->val - pwd->Vthr) {
			*h = e[1] = *f = white_vdwc;		// X-drop
			if (f2) {e[2] = *f2 = white_vdwc;}
			h->clny = 0;
		    }
		}
#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(diag); 
		putvar(f->val); putvar(e[1].val);
		if (f2) {
		    putvar(f2->val); putvar(e[2].val);
		}
		putchar('\n');
	}
#endif
	    }	// end of n-loop
	}	// end of m-loop

	delete[] cbuf;
	*scr = clny->val;
	if (algmode.mlt == 2) cl->removeoverlap();
	cl->sortcolonies();
	return (cl);
}

/********************************************************
*
*	score only algorithm
*
********************************************************/

void Aln2b1::sinitB_ng(VTYPE* hh[], const WINDOW& wdw)
{
	int	r = b->left - a->left;
	int	rr = b->right - a->left;

	VTYPE*	h = hh[0]  + r;
	*h = 0;
	if (a->inex.exgl) {	// semi-global
	    if (wdw.up < rr) rr = wdw.up;
	    vclear(h + 1, rr - r);
	}

	rr = b->left - a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	if (b->inex.exgl) {
	    vclear(hh[0] + rr, r - rr);
	    return;
	}
	VTYPE*	f = hh[1] + r;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --f;
	    *h = h[1];
	    if (i == 1) {
		*h += pwd->GapPenalty(1);
		*f = *h;
	    } else {
		*f = f[1];
		*h += pwd->GapExtPen(i);
		*f += pwd->BasicGEP;
	    }
	}
}

VTYPE Aln2b1::slastB_ng(VTYPE* hh[], const WINDOW& wdw)
{
	VTYPE*	h9 = hh[0] + b->right - a->right;
	VTYPE	mx = *h9;

	if (b->inex.exgr) {
	    int	rw = std::min(wdw.up, b->right - a->left);
	    for (VTYPE* h = hh[0] + rw; h > h9; --h)
		if (*h > mx) mx = *h;
	}
	if (a->inex.exgr) {
	    int	rw = std::max(wdw.lw, b->left - a->right);
	    for (VTYPE* h = hh[0] + rw; h < h9; ++h)
		if (*h > mx) mx = *h;
	}
	return (mx);
}

VTYPE Aln2b1::scorealoneB_ng(const WINDOW& wdw)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	VTYPE*	hh[NOL];
	VTYPE	e1;
	VTYPE	e2;
	VTYPE*	hf[NOD] = {0, &e1, 0, &e2, 0};
	VTYPE	maxh = NEVSEL;
const	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

const	size_t	bufsiz = pwd->Noll * wdw.width;
	VTYPE*	wbuf = new VTYPE[bufsiz];
	vset(wbuf, NEVSEL, bufsiz);
	VTYPE*	blackv = wbuf + bufsiz - 1;	// assume to be const
	hh[0] = wbuf - wdw.lw + 1;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + wdw.width;
	sinitB_ng(hh, wdw);

	int	m = a->left;
	if (!a->inex.exgl) --m; 	// global
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as, ++n1, ++n2) {
	    int	n = std::max(n1, b->left);
const	    int	n9 = std::min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    VTYPE*&	h = hf[0] = hh[0] + r;
	    VTYPE*&	f = hf[2] = hh[1] + r;
	    VTYPE*&	f2 = hf[4] = dagp? hh[2] + r: blackv;
	    e1 = e2 = NEVSEL;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d ", m, n);
		putvar(*h); putchar('\n');
	    }
#endif
	    for ( ; ++n <= n9; ++bs) {
		++h; ++f; if (dagp) ++f2;

//	Diagonal
		VTYPE*	from = h;
		VTYPE*	mx = h;
		VTYPE	x;
		if (m == a->left) goto HorizonF;
		*h += qprof[*bs];

//	Vertical
		x = *++from + pwd->BasicGOP;
		*f = std::max(x, f[1]) + pwd->BasicGEP;
		if (*f > *mx) mx = f;

//	Vertical2
		if (dagp) {
		  x = *from + pwd->LongGOP;
		  *f2 = std::max(x, f2[1]) + pwd->LongGEP;
		  if (*f2 > *mx) mx = f2;
		}
HorizonF:
//	Horizontal
		x = h[-1] + pwd->BasicGOP;
		e1 = std::max(x, e1) + pwd->BasicGEP;
		if (e1 > *mx) mx = &e1;

//	Horizontal2
		if (dagp) {
		    x = h[-1] + pwd->LongGOP;
		    e2 = std::max(x, e2) + pwd->LongGEP;
		    if (e2 > *mx) mx = &e2;
		}

//	Find optimal path
		VTYPE	y = *h;
		if (h != mx) *h = *mx;
		else if (LocalR && y > maxh) maxh = y;
		if (LocalL && *h < 0) *h = 0;
		int	hd = 0;
		for ( ; mx != hf[hd]; ++hd) ;

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d ", m, n);
		    putvar(*mx); putvar(y); 
		    putvar(*f); putvar(e1);
		    if (dagp) {
			putvar(*f2); putvar(e2);
		    }
		    putchar('\n');
		}
#endif
	    }	// end of n-loop
	}	// end of m-loop
	if (!LocalR) maxh = slastB_ng(hh, wdw);

	delete[] wbuf;
	return (maxh);
}

VTYPE Aln2b1::diagonalB_ng()
{
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(b->left);
	VTYPE	scr = 0;
	VTYPE	maxh = NEVSEL;
	SKL	wskl;
	int	mL = a->left;
	int	mR = a->right;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

	for (int m = a->left; m < a->right; ) {
	    scr += pwd->sim2(as++, bs++);
	    ++m;
	    if (LocalL && scr < 0) {
		scr = 0;
		mL = m;
	    }
	    if (LocalR && scr > maxh) {
		maxh = scr;
		mR = m;
	    }
	}
	int	m = b->left - a->left;
	wskl.m = mL;
	wskl.n = mL + m;
	mfd->write((UPTR) &wskl);
	wskl.m = mR;
	wskl.n = mR + m;
	mfd->write((UPTR) &wskl);
	return (LocalR? maxh: scr);
}

VTYPE Aln2b1::trcbkalignB_ng(const WINDOW& wdw)
{
	int	ptr[2];
	VTYPE	scr;

	vmf = new Vmf();

#if !__SSE4_1__	// scalar version of forward DP 
	scr = forwardB_ng(wdw, ptr);
#else
const	int	nelem = 8;
const	int	m = a->right - a->left;
	if (m < nelem) scr = forwardB_ng(wdw, ptr);
	else {
	    float	cvol = wdw.lw - b->left + a->right;
	    cvol =  (float) (a->right - a->left) * (b->right - b->left)
		 - cvol * cvol / 2;
const	    int	mode = 2 + (cvol < USHRT_MAX? 0: 4);
#if _SIMD_PP_	// vectorized forward DP 
#if __AVX2__
	    SimdAln2s1<short, 16, 
		simdpp::int16<16>, simdpp::mask_int16<16>>
#else
	    SimdAln2s1<short, SIMDPP_FAST_INT16_SIZE, 
		simdpp::int16<8>, simdpp::mask_int16<8>>
#endif	// __AVX2__
#elif __AVX512BW__
	    SimdAln2s1<short, 32, __m512i, __m512i> 
#elif __AVX2__
	    SimdAln2s1<short, 16, __m256i, __m256i>
#else	// __SSE4_1__
	    SimdAln2s1<short, 8, __m128i, __m128i> 
#endif	// _SIMD_PP_
		trbfwd(seqs, pwd, wdw, 0, 0, mode, vmf);
	    scr = trbfwd.forwardS1(ptr);
	}
#endif	// !__SSE4_1__
	if (*ptr) {
	    SKL* lskl = vmf->traceback(*ptr);
	    if (!lskl) {
		scr = NEVSEL;
		goto eop;
	    }
	    SKL* lwsk = lskl;
	    while (lskl->n--) mfd->write((UPTR) ++lwsk);
	    if (!Local) {
		lskl->m = a->left;
		lskl->n = b->left;
		if (lskl->m != lwsk->m || lskl->n != lwsk->n)
		    mfd->write((UPTR) lskl);
	    }
	    delete[] lskl;
	}
eop:
	delete vmf;
	vmf = 0;
	return (scr);
}

VTYPE Aln2b1::lspB_ng(const WINDOW& wdw)  /* recursive */ 
{
const	int	m = a->right - a->left;
const	int	n = b->right - b->left;
	if (!m && !n) return (0);
const	INT	aexgl = b->inex.exgl;	// reserve
const	INT	aexgr = a->inex.exgr;	// reserve
const	INT	bexgl = b->inex.exgl;	// reserve
const	INT	bexgr = b->inex.exgr;	// reserve
	if (!m || !n) {
	    SKL	wskl = {a->left, b->left};
	    mfd->write((UPTR) &wskl);
	    wskl.m = a->right;
	    wskl.n = b->right;
	    mfd->write((UPTR) &wskl);
	    if (m)
		return (aexgl || aexgr)? 
		    pwd->GapExtPen(m): pwd->GapPenalty(m);
	    else
		return (bexgl || bexgr)?
		    pwd->GapExtPen(n): pwd->UnpPenalty(n);
	}
	if (wdw.up == wdw.lw) return(diagonalB_ng());
const	float	k = wdw.lw - b->left + a->right;
const	float	q = b->right - a->left - wdw.up;
	float	cvol =  m * n - (k * k + q * q) / 2;
	if (cvol < MaxVmfSpace || m == 1 || n <= 1)
	    return (trcbkalignB_ng(wdw));

	RANGE	rng[2];			// reserve
	save_range(seqs, rng, 2);
	int	cpos[6];
	VTYPE	scr;
	WINDOW	wdwf = {INT_MAX};
	WINDOW	wdwb = {INT_MAX};
	int	upf = INT_MAX; 
	int	upb = INT_MAX;
	INT	exgl = bexgl;

#if !__SSE4_1__	// scalar version of unidirectional Hirschberg method
	scr = hirschbergB_ng(cpos, wdw, wdwf, wdwb, exgl);
#else
	if (m < 4) scr = hirschbergB_ng(cpos, wdw, wdwf, wdwb, exgl);
	else {
	    int	mode = 1 + ((std::max(abs(wdw.lw), wdw.up) < SHRT_MAX)? 0: 4);
#if _SIMD_PP_	// 2B int vectorized unidirectional Hirschberg method
#if __AVX2__
	    SimdAln2s1<short, 16, 
		simdpp::int16<16>, simdpp::mask_int16<16>>
#else
	    SimdAln2s1<short, SIMDPP_FAST_INT16_SIZE, 
		simdpp::int16<8>, simdpp::mask_int16<8>>
#endif	// __AVX2__
#elif __AVX512BW__
	    SimdAln2s1<short, 32, __m512i, __m512i> 
#elif __AVX2__
	    SimdAln2s1<short, 16, __m256i, __m256i>
#else	// __SSE4_1__
	    SimdAln2s1<short, 8, __m128i, __m128i> 
#endif	// _SIMD_PP_
		sb1(seqs, pwd, wdw, 0, 0, mode);
	    scr = sb1.hirschbergS1(cpos, upf, upb, exgl);
	}
#endif	// !__SSE4_1__
	a->inex.exgl = a->inex.exgr = b->inex.exgl = b->inex.exgr = 0;
	SKL	wskl = {cpos[0], 0};
	int	c = 0;
	if (cpos[c] < end_of_ulk) {	// cross center
	    while (cpos[++c] < end_of_ulk) {
		wskl.n = cpos[c];
		mfd->write((UPTR) &wskl);
	    }
const	    int	aright = a->right;	// reserve
const	    int	bright = b->right;
	    a->right = cpos[0];
	    b->right = cpos[c - 1];
	    if (upf < INT_MAX || wdwf.lw == INT_MAX) {
		stripe(seqs, &wdwf, alprm.sh);
		if (upf < INT_MAX) wdwf.up = upf;
	    }
	    wdwf.width = wdwf.up - wdwf.lw + 3;
	    lspB_ng(wdwf);		// first half
	    a->left = cpos[0];
	    b->left = cpos[1];
	    a->right = aright;		// recover
	    b->right = bright;
	    b->inex.exgl = exgl;
	}
	if (upb < INT_MAX || wdwb.lw == INT_MAX) {		// local
	    stripe(seqs, &wdwb, alprm.sh);
	    if (upb < INT_MAX) wdwb.up = upb;
	}
	wdwb.width = wdwb.up - wdwb.lw + 3;
	lspB_ng(wdwb);			// second half
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexgl;
	a->inex.exgr = aexgr;
	b->inex.exgl = bexgl;
	b->inex.exgr = bexgr;
	return scr;
}

VTYPE Aln2b1::backforth(int ovr, const BOUND& lub)
{
	VTYPE*	bscr = new VTYPE[ovr + 1];
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(b->left);
	VTYPE	scr = bscr[ovr] = 0;
	int	i = ovr;
	int	m = a->left;
	int	n = b->left;
	while (--i >= 0 && --m >= lub.la && --n >= lub.lb)
	    bscr[i] = scr += pwd->sim2(--as, --bs);
	VTYPE	maxscr = scr;
	scr = 0;
	int	ii = ++i;
	m = a->right + i;
	n = b->right + i;
	as = a->at(m);
	bs = b->at(n);
	while (i++ < ovr && m++ < lub.ua && n++ <lub.ub) {
	    scr += pwd->sim2(as++, bs++);
	    if ((bscr[i] += scr) > maxscr) {
		maxscr = bscr[i]; ii = i;
	    }
	}
	SKL	skl = {a->right + ii, n = b->right + ii};
	mfd->write((UPTR) &skl);
	int dr = (b->right - a->right) - (b->left - a->left);
	if (dr >= 0) skl.n -= dr;
	else	skl.m -= (dr = -dr);
	mfd->write((UPTR) &skl);
	delete[] bscr;
	return (maxscr + pwd->GapPenalty(dr));
}

VTYPE Aln2b1::creepback(int ovr, VTYPE bscr, const BOUND& lub)
{
const	CHAR*   as = a->at(a->left);
const	CHAR*   bs = b->at(b->left);
	VTYPE   dscr = 0;
	while (a->left > lub.la && b->left > lub.lb
		&& (ovr < 0 || dscr < bscr)) {
	    dscr += pwd->sim2(--as, --bs);
	    a->left--; b->left--;
	    if (++ovr == 0) bscr += dscr;
	}
	return (dscr);
}
		                                                                
VTYPE Aln2b1::creepfwrd(int& ovr, VTYPE bscr, const BOUND& lub)
{
const	CHAR*   as = a->at(a->right);
const	CHAR*   bs = b->at(b->right);
	VTYPE   dscr = 0;
	while (a->right < lub.ua && b->right < lub.ub
		&& (ovr < 0 || dscr < bscr)) {
	    dscr += pwd->sim2(as++, bs++);
	    a->right++; b->right++;
	    if (++ovr == 0) bscr += dscr;
	}
	return (dscr);
}

/* eimode 1: 5'end, 2: 3' end, 3: internal */
VTYPE Aln2b1::seededB_ng(INT level, int eimode, const BOUND& lub)
{
	INEX	ainex = a->inex;
	INEX	binex = b->inex;
	RANGE	rng[2];
	int	cmode = eimode;
	int	agap = a->right - a->left;
	int	bgap = b->right - b->left;
	int	num = 0;
	VTYPE	scr = 0;
	JUXT*	jxt = 0;
	JUXT*	wjxt = 0;
	Wilip*	wl = 0;
	WLUNIT*	wlu = 0;

	if (level == lowestlvl && b->jxt) {
	    jxt = (JUXT*) b->jxt;
	    num = b->CdsNo;
	} else {
	    wl = new Wilip(seqs, pwd, level);
	    wlu = wl->begin();
	    if (wlu) {
		num = wlu->num;
		jxt = wlu->jxt;
	    } else num = 0;
	}
	save_range(seqs, rng, 2);
	VTYPE	slmt = pwd->Vthr / 2;
	WLPRM*	wlprm = setwlprm(level);
	int	backstep = wlprm->width;
	if (++level < algmode.qck) {
	    wlprm = setwlprm(level);
	    backstep = wlprm->width;
	}
	BOUND   bab = lub;
	if (num) {
	    jxt[num].jx = a->right;
	    jxt[num].jy = b->right;
	    a->inex.exgr = 0;
	    b->inex.exgr = 0;
	    agap = jxt->jx - a->left;
	    bgap = jxt->jy - b->left;
	    for (wjxt = jxt; num--; wjxt++) {
		scr += wjxt->jscr;
		a->right = wjxt->jx;
		b->right = wjxt->jy;
		bab.ua = wjxt->jx + wjxt->jlen;
		bab.ub = wjxt->jy + wjxt->jlen;
		if (wjxt[1].jx < bab.ua) bab.ua = wjxt[1].jx;
		if (wjxt[1].jy < bab.ub) bab.ub = wjxt[1].jy;
		bab.ua -= backstep;
		bab.ub -= backstep;
		int     lx = wjxt->jx + wjxt->jlen / 2;
		int     ly = wjxt->jy + wjxt->jlen / 2;
		if (bab.ua < lx) bab.ua = lx;
		if (bab.ub < ly) bab.ub = ly;
		ly = wjxt->jy + bab.ua - wjxt->jx;
		lx = wjxt->jx + bab.ub - wjxt->jy;
		if (lx < bab.ua) bab.ua = ly;
		if (ly < bab.ub) bab.ub = lx;
		int	ovr = min(agap, bgap);
		if (cmode == 1 && agap == 0)    // 5' end
		    mfd->write((UPTR) wjxt);
		else if (cmode == 3 && ovr <= 0)
		    scr += backforth(-ovr, bab);
		else if (level < algmode.qck && ovr >= (int) wlprm->tpl)
		    scr += seededB_ng(level, cmode, bab);
		else {			      // recursive search
		    scr -= creepback(ovr, slmt, bab);
		    scr -= creepfwrd(ovr, slmt, bab);
		    if (ovr < 0) {          // skip this hsp
			agap = wjxt[1].jx - a->left;
			bgap = wjxt[1].jy - b->left;
			continue;
		    } else {
			WINDOW	wdw;
			stripe(seqs, &wdw, alprm.sh);
			scr += lspB_ng(wdw);
		    }
		}
		cmode = 3;
		a->left = wjxt->jx + wjxt->jlen;
		b->left = wjxt->jy + wjxt->jlen;
		a->inex.exgl = 0;
		b->inex.exgl = 0;
		agap = wjxt[1].jx - a->left;
		bgap = wjxt[1].jy - b->left;
		bab.la = a->right;
		bab.lb = b->right;
	    }
	    a->inex.exgr = ainex.exgr;
	    b->inex.exgr = binex.exgr;
	    a->right = rng[0].right;
	    b->right = rng[1].right;
	    bab.ua = lub.ua;
	    bab.ub = lub.ub;
	    if (eimode == 2 || (level == INT(lowestlvl + 1) && eimode == 1))
		cmode = 2;
	}

	int	ovr = min(agap, bgap);
	if (cmode == 2 && ovr == 0) {
	    SKL	tmp = {a->left, b->left};
	    mfd->write((UPTR) &tmp);
	} else if (cmode == 3 && ovr <= 0)
	    scr += backforth(-ovr, bab);
	else if (level < algmode.qck && ovr > (int) wlprm->tpl)
	    scr += seededB_ng(level, cmode, bab);
	else {
	    scr -= creepback(ovr, slmt, bab);
	    scr -= creepfwrd(ovr, slmt, bab);
	    WINDOW	wdw;
	    stripe(seqs, &wdw, alprm.sh);
	    scr += lspB_ng(wdw);
	}
	rest_range(seqs, rng, 2);
	if (--level == lowestlvl && jxt) {
	    wjxt->jx = a->len;
	    wjxt->jy = b->len;
	}
	a->inex.exgl = ainex.exgl;
	b->inex.exgl = binex.exgl;
	delete wl;
	return (scr);
}

SKL* Aln2b1::globalB_ng(VTYPE* scr, const WINDOW& wdw)
{
	mfd = new Mfile(sizeof(SKL));
	SKL	wsk;
	mfd->write(&wsk);	// make a space record
	if (!a->inex.exgl && algmode.qck) {
	    wsk.m = a->left;
	    wsk.n = b->left;
	    mfd->write((UPTR) &wsk);
	}
	if (!a->inex.exgr && algmode.qck) {
	    wsk.m = a->right;
	    wsk.n = b->right;
	    mfd->write((UPTR) &wsk);
	}
	if (algmode.qck) {
	    BOUND bab = {a->left, b->left, a->right, b->right};
	    *scr = seededB_ng(lowestlvl, 1, bab);
	} else
	    *scr = lspB_ng(wdw);
	wsk.n = (int) mfd->size();
	SKL*	skl = (SKL*) mfd->flush();
	skl->n = wsk.n - 1;
	skl->m = 1;
	if (skl->n == 0) {
	    delete[] skl;
	    return (0);
	}
	return (stdskl(&skl));
}

VTYPE HomScoreB_ng(const Seq* seqs[], const PwdB* pwd)
{
	WINDOW	wdw;
	stripe((const Seq**) seqs, &wdw, alprm.sh);
#if !__SSE4_1__	// scalar version of forward DP 
	Aln2b1 alnv(seqs, pwd);
	return (alnv.scorealoneB_ng(wdw));
#else
const	Seq*&	a = seqs[0];
const	int	m = a->right - a->left;
	if (m < 4) {
	    Aln2b1 alnv(seqs, pwd);
	    return (alnv.scorealoneB_ng(wdw));
	}
#if _SIMD_PP_	// vectorized forward DP 
#if __AVX2__
	SimdAln2s1<short, 16, simdpp::int16<16>, simdpp::mask_int16<16>>
#else
	SimdAln2s1<short, 8, simdpp::int16<8>, simdpp::mask_int16<8>>
#endif	// __AVX2__
#elif __AVX512BW__
	SimdAln2s1<short, 32, __m512i, __m512i> 
#elif __AVX2__
	SimdAln2s1<short, 16, __m256i, __m256i>
#else	// __SSE4_1__
	SimdAln2s1<short, 8, __m128i, __m128i> 
#endif	// _SIMD_PP_
	    fwds(seqs, pwd, wdw, 0, 0, 0);
	return (fwds.scoreonlyS1());
#endif	// !__SSE4_1__
}

Colonies* swg1stB_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr)
{
	Aln2b1	alnv(seqs, pwd);
	WINDOW	wdw;
	stripe((const Seq**) seqs, &wdw, alprm.sh);
	return (alnv.fwdswgB_ng(scr, wdw));
}

SKL* swg2ndB_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr, COLONY* clny)
{
	if (clny->val <= 0) return (0);
const	Seq*	a = seqs[0];
const	Seq*	b = seqs[1];
	a->left = clny->mlb;
	b->left = clny->nlb;
	b->right = clny->nrb;
	a->right = clny->mrb;
	Aln2b1 alnv(seqs, pwd);
	WINDOW	wdw = {clny->upr, clny->lwr, clny->upr - clny->lwr + 3};
	return (alnv.globalB_ng(scr, wdw));
}

SKL* alignB_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr)
{
	Aln2b1 alnv(seqs, pwd);
	WINDOW	wdw;
	stripe((const Seq**) seqs, &wdw, alprm.sh);
	return (alnv.globalB_ng(scr, wdw));
}

Gsinfo* localB_ng(Seq* seqs[], const PwdB* pwd) {
	Wilip*	wl = new Wilip((const Seq**) seqs, pwd, 0);
	WLUNIT*	wlu = wl->begin();

	if (!wlu) {delete wl; return(0);}
	WINDOW	wdw;
	stripe((const Seq**) seqs, &wdw, alprm.sh);
	Seq*	b = seqs[1];
	int	n = min((int) OutPrm.MaxOut, wl->size());
	Gsinfo*	mai = new Gsinfo[n + 1];
	for (Gsinfo* wai = mai; n--; ++wlu, ++wai) {
	    b->jxt = wlu->jxt;
	    b->CdsNo = wlu->num;
	    Aln2b1* alnv = new Aln2b1((const Seq**) seqs, pwd);
	    wai->skl = alnv->globalB_ng(&wai->scr, wdw);
	    delete alnv;
	}
	delete wl;
	b->jxt = 0;
	b->CdsNo = 0;
	return (mai);
}

