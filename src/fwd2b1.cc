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

#define	DEBUG	0

#include "aln.h"
#include "vmf.h"
#include "wln.h"

class Aln2b1 {
const	Seq**	seqs;
const	Seq*&	a;
const	Seq*&	b;
const	PwdB*	pwd;
	Mfile*	mfd;
	Vmf*	vmf;
	INT	lowestlvl;
public:
	Aln2b1(const Seq** _seqs, const PwdB* _pwd) :
	    seqs(_seqs), a(seqs[0]), b(seqs[1]), 
	    pwd(_pwd), mfd(0), vmf(0), lowestlvl(b->wllvl) {} 
	~Aln2b1() {delete vmf; delete mfd;}
	void	initB_ng(RVPD* hh[], const WINDOW& wdw);
	RVPD*	lastB_ng(RVPD* hh[], const WINDOW& wdw);
	VTYPE	forwardB_ng(const WINDOW& wdw, long pp[] = 0);
	void	finitB_ng(RVDWL* hhf[], int mm, const WINDOW& wdw);
	void	binitB_ng(RVDWL* hhb[], int mm, const WINDOW& wdw);
	VTYPE	centerB_ng(int* ml, int* mr, int* nl, int* nr,
		const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb);
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
	long	origin = vmf? vmf->add(a->left, b->left, 0): r;
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
	    h->ptr = vmf? origin: r;
	    VTYPE	gpn = (i == 1)? pwd->GapPenalty(1): pwd->GapExtPen(i);
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
	    h->ptr = vmf? origin: r;
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
	if (vmf && (dn || dm)) {
	    if (dn) dm = 0;
	    h9->ptr = vmf->add(a->right - dm, b->right - dn, h9->ptr);
	}
	if (vmf) h9->ptr = vmf->add(a->right, b->right, h9->ptr);
	return (h9);
}

VTYPE Aln2b1::forwardB_ng(const WINDOW& wdw, long pp[])
{
	RVPD*	hh[NOL];
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

	RVPD*	rbuf = new RVPD[pwd->Noll * wdw.width];
	vset(rbuf, black_vpd, pwd->Noll * wdw.width);
	hh[0] = rbuf - wdw.lw + 1;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + wdw.width;
	if (vmf) vmf->add(0, 0, 0);	// Skip 0-th record
	initB_ng(hh, wdw);
	int	m = a->left;
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
const	CHAR*	as = a->at(m);

	for ( ; m < a->right; ++as, ++m, ++n1, ++n2) {
	    RVPD	f1 = black_vpd;
	    RVPD	f2 = black_vpd;
	    int		n = max(n1, b->left);
	    int		n9 = min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int		r = n - m;
	    int		k = 0;
	    RVPD*	h = hh[k] + r;
	    RVPD*	g = hh[++k] + r;
	    RVPD*	g2 = (pwd->Noll == 3)? hh[++k] + r: 0;
#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m + 1, n + 1, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; n < n9; ++bs, ++n, ++h, ++g) {
		VTYPE	x;
//	Diagonal
		RVPD*	from = h;
		RVPD*	mx = h;
		VTYPE	diag = h->val;
		h->val += pwd->sim2(as, bs);
		h->dir = isdiag(from)? DIAG: NEWD;

//	Vertical
		x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    g->val = x;
		    g->ptr = from->ptr;
		    g->dir = VERT;
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val > mx->val) mx = g;

//	Vertical2
		if (g2) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[1].val) {
		    g2->val = x;
		    g2->ptr = from->ptr;
		    g2->dir = VERT;
		  } else	*g2 = g2[1];
		  g2->val += pwd->LongGEP;
		  if (g2->val > mx->val) mx = g2;
		}
//	Horizontal
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f1.val) {
		    f1.val = x;
		    f1.ptr = h[-1].ptr;
		    f1.dir = HORI;
		}
		f1.val += pwd->BasicGEP;
		if (f1.val >= mx->val) mx = &f1;

//	Horizontal2
		if (g2) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= f2.val) {
			f2.val = x;
			f2.ptr = h[-1].ptr;
			f2.dir = HORL;
		    }
		    f2.val += pwd->LongGEP;
		    if (f2.val >= mx->val) mx = &f2;
		}

//	Find optimal path
		if (h != mx) *h = *mx;	// non-diagonal
		else if (Local && h->val > diag) {
		    if (LocalL && diag == 0)
			h->ptr = vmf? vmf->add(m, n, 0): 0;
		    else if (LocalR && h->val > maxh.val) {
			maxh.val = h->val;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && h->val <= 0) h->val = h->dir = 0;
		else if (vmf) {
		    if (h->dir == NEWD || h->dir == NEWV || h->dir == NEWH)
			h->ptr = vmf->add(m, n, h->ptr);
		} else if (m == a->left)
		    h->ptr = n - m;

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d %2d ", m + 1, n + 1, mx->dir);
		    putvar(mx->val); putvar(diag); 
		    putvar(g->val); putvar(f1.val);
		    if (g2) {
			putvar(g2->val); putvar(f2.val);
		    }
		    putchar('\n');
		}
#endif
		if (g2) ++g2;
	    }	/* end of n-loop */
	}	/* end of m-loop */

	if (LocalR) {
	    if (vmf) *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    RVPD*	mx = lastB_ng(hh, wdw);
	    maxh.val = mx->val;
	    if (vmf) {
		pp[0] = mx->ptr;
		pp[1] = b->left - a->left + mx - hh[0];
	    }
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

void Aln2b1::finitB_ng(RVDWL* hhf[], int mm, const WINDOW& wdw)
{
	int	r = b->left - a->left;
	int	r0 = r;
	int	rr = b->right - a->left;
	FTYPE	ltgapf = a->left? 1: (a->inex.exgl? 0: alprm.tgapf);

	RVDWL*	h = hhf[0] + r;
	RVDWL*	g = hhf[1] + r;
	h->val = 0;
	h->dir = DEAD;
	h->lwr = h->upr = h->lst = r;
	if (wdw.up < rr) rr = wdw.up;
	for (int i = 1; ++r <= rr; ++i) {
	    ++h;
	    if (a->inex.exgl) {
		h->val = 0;
		h->dir = DEAD;
		h->lwr = h->upr = h->lst = r;
	    } else {
		h->dir = HORI;
		h->lwr = h->lst = r0;
		h->upr = r;
		VTYPE	gpn = (i == 1)? pwd->GapPenalty(1): pwd->GapExtPen(i);
		h->val = h[-1].val + (VTYPE) (gpn * ltgapf);
	    }
	}

	ltgapf = b->left? 1: (b->inex.exgl? 0: alprm.tgapf);
	r = r0;
	rr = b->left - mm;
	if (wdw.lw > rr) rr = wdw.lw;
	h = hhf[0] + r;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --g;
	    if (b->inex.exgl) {
		h->val = 0;
		h->dir = DEAD;
		h->lwr = h->upr = h->lst = r;
	    } else {
		*h = h[1];
		if (i == 1) {
		    h->val += (VTYPE) (pwd->GapPenalty(1) * ltgapf);
		    h->dir = VERT;
		    *g = *h;
		} else {
		    *g = g[1];
		    h->val += (VTYPE) (pwd->GapExtPen(i) * ltgapf);
		    g->val += (VTYPE) (pwd->BasicGEP * ltgapf);
		}
		h->lwr = g->lwr = r;
	    }
	}
}

void Aln2b1::binitB_ng(RVDWL* hhb[], int mm, const WINDOW& wdw)
{
	int	r = b->right - a->right;
	int	r9 = r;
	int	rr = b->left - a->right;
	FTYPE	rtgapf = a->right < a->len? 1: (a->inex.exgr? 0: alprm.tgapf);

	RVDWL*	h = hhb[0] + r;
	h->val = 0;
	h->dir = DEAD;
	h->lwr = h->upr = h->lst = r;
	if (wdw.lw > rr) rr = wdw.lw;
	for (int i = 1; --r >= rr; ++i) {
	    --h;
	    if (a->inex.exgr) {
		h->val = 0;
		h->dir = DEAD;
		h->lwr = h->upr = h->lst = r;
	    } else {
		h->dir = HORI;
		h->upr = h->lst = r9;
		h->lwr = r;
		VTYPE	gpn = (i == 1)? pwd->GapPenalty(1): pwd->GapExtPen(i);
		h->val = h[1].val + (VTYPE) (gpn * rtgapf);
	    }
	}

	rtgapf = b->right < b->len? 1: (b->inex.exgr? 0: alprm.tgapf);
	r = r9;
	rr = b->right - mm;
	if (wdw.up < rr) rr = wdw.up;
	h = hhb[0] + r;
	RVDWL*	g = hhb[1] + r;
	for (int i = 1; ++r <= rr; ++i) {
	    ++h; ++g;
	    if (b->inex.exgr) {
		h->val = 0;
		h->dir = DEAD;
		h->lwr = h->upr = h->lst = r;
	    } else {
		*h = h[-1];
		if (i == 1) {
		    h->val += (VTYPE) (pwd->GapPenalty(1) * rtgapf);
		    h->dir = VERT;
		    *g = *h;
		} else {
		    *g = g[-1];
		    h->val += (VTYPE) (pwd->GapExtPen(i) * rtgapf);
		    g->val += (VTYPE) (pwd->BasicGEP * rtgapf);
		}
		h->upr = g->upr = r;
	    }
	}
}

VTYPE Aln2b1::centerB_ng(int* ml, int* mr, int* nl, int* nr, 
	const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb)
{
	int	kk = 0;
	int	rr[2] = {0, 0};
	RVDWL*	hhf[NOL];
	RVDWL*	hhb[NOL];
	RVDWL*	hf[NOL];
	RVDWL*	hb[NOL];
	RVDWL	f[NOL];
	RVDWL*	f1 = f + 1;
	RVDWL*	g2 = 0;
	VTYPE	mxh = NEVSEL;
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

	RVDWL*	wbuf = new RVDWL[2 * pwd->Noll * wdw.width];
	vset(wbuf, black_vdwl, 2 * pwd->Noll * wdw.width);
	hhf[0] = wbuf - wdw.lw + 1;
	hhb[0] = hhf[0] + pwd->Noll * wdw.width;
	for (int k = 1; k < pwd->Noll; ++k) hhf[k] = hhf[k-1] + wdw.width;
	for (int k = 1; k < pwd->Noll; ++k) hhb[k] = hhb[k-1] + wdw.width;
	int	mm = (a->left + a->right + 1) / 2;
	RVDWL*	f2 = (pwd->Noll == 3)? f + 2: 0;

/*******	backward phase	*******/

	binitB_ng(hhb, mm, wdw);
	int	m = a->right - 1;
const	CHAR*	as = a->at(m);
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
	for ( ; m >= mm; --m, --as, --n1, --n2) {
	    int	n = min(n2, b->right) - 1;
	    int	n9 = max(n1, b->left);
const	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    int	d = 0;
	    RVDWL*	h = hhb[d] + r;
	    RVDWL*	g = hhb[++d] + r;
	    *f1 = black_vdwl;
	    if (pwd->Noll == 3) {
		g2 = hhb[++d] + r;
		*f2 = black_vdwl;
	    }
#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m + 1, n + 1, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; n >= n9; --n, --bs, --r, --h, --g) {
/*	Diagonal	*/
		RVDWL*	from = h;
		RVDWL*	mx = h;
		h->val += pwd->sim2(as, bs);
		h->dir = isdiag(from)? DIAG: NEWD;
/*	Vertical	*/
		VTYPE	x = (--from)->val + pwd->BasicGOP;
		if (x >= g[-1].val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT;
		} else *g = g[-1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;

/*	Vertical2	*/
		if (g2) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[-1].val) {
		    *g2 = *from;
		    g2->val = x;
		    g2->dir = VERT;
		  } else *g2 = g2[-1];
		  g2->val += pwd->LongGEP;
		  if (g2->val >= mx->val) mx = g2;
		}
/*	Horizontal	*/
		x = h[1].val + pwd->BasicGOP;
		if (x >= f1->val) {
		    *f1 = h[1];
		    f1->val = x;
		    f1->dir = HORI;
		}
		f1->val += pwd->BasicGEP;
		if (f1->val >= mx->val) mx = f1;

/*	Horizontal2	*/
		if (f2) {
		    x = h[1].val + pwd->LongGOP;
		    if (x >= f2->val) {
			*f2 = h[1];
			f2->val = x;
			f2->dir = HORI;
		    }
		    f2->val += pwd->LongGEP;
		    if (f2->val >= mx->val) mx = f2;
		}
		if (mx->dir == NEWD) mx->lst = r;

/*	Find optimal path	*/
#if DEBUG
		VTYPE	y = h->val;
#endif
		if (h != mx) {
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		}
		if (LocalR && h->val <= 0) h->val = h->dir = 0;
#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m + 1, n + 1, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(f1->val);
		if (g2) {
		    putvar(g2->val); putvar(f2->val);
		}
		putchar('\n');
	}
#endif
		if (g2) --g2;
	    }	/* end of n loop */
	}	/* end of m loop */

/*******	forward phase	*******/

	finitB_ng(hhf, mm, wdw);
	m = a->left;
	n1 = m + wdw.lw;
	n2 = m + wdw.up + 1;
	as = a->at(m);
	int	mm1 = mm - 1;
	for ( ; m < mm; ++m, ++as, ++n1, ++n2) {
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    int	d = 0;
	    RVDWL*	h = hhf[d] + r;
	    RVDWL*	g = hhf[++d] + r;
	    *f1 = black_vdwl;
	    if (g2) {
		g2  = hhf[++d] + r;
		*f2 = black_vdwl;
	    }
	    if (m == mm1)
		for (int k = 0; k < pwd->Noll; k++) hb[k] = hhb[k] + r;

	    for ( ; n < n9; ++n, ++bs, ++h, ++g, ++r) {
		RVDWL*	from = h;
		RVDWL*	mx = h;

/*	Diagonal	*/
		h->val += pwd->sim2(as, bs);
		h->dir = isdiag(from)? DIAG: NEWD;

/*	Vertical	*/
		VTYPE	x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT;
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;

/*	Vertical2	*/
		if (g2) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[1].val) {
		    *g2 = *from;
		    g2->val = x;
		    g2->dir = VERT;
		  } else 	*g2 = g2[1];
		  g2->val += pwd->LongGEP;
		  if (g2->val >= mx->val) mx = g2;
		}
/*	Horizontal	*/
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f1->val) {
		    *f1 = h[-1];
		    f1->val = x;
		    f1->dir = HORI;
		}
		f1->val += pwd->BasicGEP;
		if (f1->val >= mx->val) mx = f1;

/*	Horizontal2	*/
		if (f2) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= f2->val) {
			*f2 = h[-1];
			f2->val = x;
			f2->dir = HORI;
		    }
		    f2->val += pwd->LongGEP;
		    if (f2->val >= mx->val) mx = f2;
		}
		if (mx->dir == NEWD) mx->lst = r;

/*	Find optimal path	*/
		VTYPE	y = h->val;
		if (h != mx)  {
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		}
		if (LocalL && h->val <= 0) h->val = h->dir = 0;

#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m + 1, n + 1, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(f1->val);
		if (g2) {
		    putvar(g2->val); putvar(f2->val);
		}
		putchar('\n');
	}
#endif

/* Find Center */
		if (m == mm1) {
		    x = h->val + (*hb)->val;
		    int	k1 = isvert(h)? 1: 0;
		    int	k2 = isvert(*hb)? 1: 0;
		    y = g->val + hb[1]->val - pwd->BasicGOP;
		    int	l = g->lst - hb[1]->lst - pwd->codonk1;
		    if (l > 0) y += pwd->diffu * l;
		    if (y >= x) {x = y; k1 = k2 = 1;}
		    if (g2) {
			y = g->val + hb[2]->val - pwd->BasicGOP + pwd->diffu * (g->lst - r);
			if (y > x) {x = y; k1 = 1; k2 = 2;}
			y = g2->val + hb[1]->val - pwd->BasicGOP + pwd->diffu * (r - hb[1]->lst);
			if (y > x) {x = y; k1 = 2; k2 = 1;}
			y = g2->val + hb[2]->val - pwd->LongGOP;
			if (y > x) {x = y; k1 = k2 = 2;}
		    }
		    if (gt(x, mxh)) {
			mxh = x;
			rr[0] = rr[1] = r;
			kk = k1 + 4 * k2;
		    } else if (ge(x, mxh)) {
			rr[1] = r;
		    }
		    for (int k = 0; k < pwd->Noll; k++) ++hb[k];
		}
		if (g2) ++g2;
	    }	/* end of n-loop */
	}	/* end of m-loop */

	int	k1 = kk % 4;
	int	k2 = kk / 4;
	int	r = rr[0];
	*nl = r + mm;
	for (int k = 0; k < pwd->Noll; ++k) hf[k] = hhf[k] + r;
	wdwf.up = hf[k1]->upr;
	wdwf.lw = hf[k1]->lwr;
	wdwf.width = wdwf.up - wdwf.lw + 3;
	SKL	wskl;
	if (k1) {
	    wskl.m = *ml = *nl - hf[k1]->lst;
	    wskl.n = *nl;
	    mfd->write((UPTR) &wskl);
	} else {
	    wskl.m = *ml = mm;
	}
	wskl.n = *nr = (r = rr[1]) + mm;
	for (int k = 0; k < pwd->Noll; ++k) hb[k] = hhb[k] + r;
	wdwb.up = hb[k2]->upr;
	wdwb.lw = hb[k2]->lwr;
	wdwb.width = wdwb.up - wdwb.lw + 3;
	if (k2) {
	    mfd->write((UPTR) &wskl);
	    wskl.m = *mr = *nr - hb[k2]->lst;
	    mfd->write((UPTR) &wskl);
	} else {
	    *mr = mm;
	}

#if MONITOR
	stop = time(0);
	printf("(%4d, %4d) (%4d, %4d) %4d %4d %4d %4ld\n",
	a->left,b->left,a->right,b->right,mm,wdw.lw,wdw.up, stop - start);
#endif
	delete[] wbuf;
	return (mxh);
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
	RVDWC*	hhc[NOL];
	RVDWC   hl[NOL];
	RVDWC	f[NOL]; 	/* [DIAG, HORI, HORL] */
	RVDWC*	g2 = 0;
	Colonies*	cl = new Colonies(0);
	COLONY*	clny = cl->at();

	RVDWC*	cbuf = new RVDWC[pwd->Noll * wdw.width];
	vset(cbuf, black_vdwc, pwd->Noll * wdw.width);
        hhc[0] = cbuf - wdw.lw + 1;
        for (int k = 1; k < pwd->Noll; ++k) hhc[k] = hhc[k-1] + wdw.width;
	cinitB_ng(hhc, wdw);
	int	m  = a->left;
const	CHAR*	as = a->at(m);
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
	for ( ; m < a->right; ++m, ++as, ++n1, ++n2) {
	    int 	n = max(n1, b->left);
	    int 	n9 = min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int 	r = n - m;
	    RVDWC*	hlb = hhc[0] + r;
	    RVDWC*	hrb = hhc[0] + n9 - m;
	    RVDWC*	h = hlb;
	    RVDWC*	g = hhc[1] + r;
	    if (pwd->Noll == 3) g2 = hhc[2] + r;
	    for (int k = 0; k < pwd->Noll; ++k)
		hl[k] = f[k] = black_vdwc;

#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m+1, n+1, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; n < n9; ++n, ++r, ++h, ++g, ++bs) {
/*	Diagonal	*/
		RVDWC*	mx = h;
		f[0] = *h;	/* diag */
#if DEBUG
		VTYPE	diag = h->val;
#endif
		h->val += pwd->sim2(as, bs);
		h->dir = isdiag(h)? DIAG: NEWD;

/*	Vertical	*/
		VTYPE	x = h[1].val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    *g = h[1];
		    g->val = x;
		    g->dir = VERT;
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val > mx->val) mx = g;

/*	Vertical2	*/
		if (g2) {
		  x = h[1].val + pwd->LongGOP;
		  if (x >= g2[1].val) {
		    *g2 = h[1];
		    g2->val = x;
		    g2->dir = VERT;
		  } else 	*g2 = g2[1];
		  g2->val += pwd->LongGEP;
		  if (g2->val > mx->val) mx = g2;
		}

/*	Horizontal	*/
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f[1].val) {
		    f[1] = h[-1];
		    f[1].dir = HORI;
		    f[1].val = x;
		}
		if (f[1].dir) {
		    f[1].val += pwd->BasicGEP;
		    f[1].dir = HORI;
		    if (f[1].val >= mx->val) mx = f + 1;
		}

/*	Horizontal2	*/
		if (pwd->Noll == 3) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= f[2].val) {
			f[2] = h[-1];
			f[2].val = x;
			f[2].dir = HORL;
		    }
		    if (f[2].dir) {
			f[2].val += pwd->LongGEP;
			f[2].dir = HORL;
			if (f[2].val >= mx->val) mx = f + 2;
		    }
		}

/*	Find optimal path	*/
		if (h != mx) {		/* non-diagonal */
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		} else if (h->val > f->val) {
		    if (f->val == 0) {	/* new colony */
			h->upr = h->lwr = r;
			h->mlb = m;
			h->nlb = n;
		    }
		    if (h->val > clny->val) {	/* max local score */
			clny->val = h->val;
			clny->mrb = m + 1;
			clny->nrb = n + 1;
			clny->lwr = h->lwr;
			clny->upr = h->upr;
			clny->mlb = h->mlb;
			clny->nlb = h->nlb;
		    }
		}
		if (h->val < 0) {		/* reset to blank */
		    *h = f[1] = *g = white_vdwc;
		    if (g2) {f[2] = *g2 = white_vdwc;}
		    h->clny = 0;
		}
		if (algmode.mlt > 1 && h->val >= pwd->Vthr && !h->clny) {
		    if (cl->size() >= (int) OutPrm.MaxOut) {	/* Gabage collecton */
			for (RVDWC* hw = hlb; hw < hrb; ++hw)	/* mark active colony */
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
			cc->mrb = m + 1;
			cc->nrb = n + 1;
			cc->lwr = h->lwr;
			cc->upr = h->upr;
			cc->mlb = h->mlb;
			cc->nlb = h->nlb;
		    } else if (h->val <= cc->val - pwd->Vthr) {
			*h = f[1] = *g = white_vdwc;		// X-drop
			if (g2) {f[2] = *g2 = white_vdwc;}
			h->clny = 0;
		    }
		}
#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m + 1, n + 1, mx->dir);
		putvar(mx->val); putvar(diag); 
		putvar(g->val); putvar(f[1].val);
		if (g2) {
		    putvar(g2->val); putvar(f[2].val);
		}
		putchar('\n');
	}
#endif
		if (g2) ++g2;
	    }	/* end of n-loop */
	}	/* end of m-loop */

	delete[] cbuf;
	*scr = clny->val;
	if (algmode.mlt == 2) cl->removeoverlap();
	cl->sortcolonies();
	return (cl);
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
	int	Local = algmode.lcl & 16;
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
	long	ptr[2];

	vmf = new Vmf();
	VTYPE	scr = forwardB_ng(wdw, ptr);
	if (*ptr) {
	    SKL* lskl = vmf->traceback(*ptr);
	    if (!lskl) {
		scr = NEVSEL;
		goto eop;
	    }
	    SKL* lwsk = lskl;
	    while (lskl->n--) mfd->write((UPTR) ++lwsk);
	    if (!(algmode.lcl & 16)) {
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
	WINDOW	wdwl;
	WINDOW	wdwr;
	RANGE	rng[2];			/* reserve */
	int	aexg = a->inex.exgr;	/* reserve */
	int	bexg = b->inex.exgr;	/* reserve */

	if (wdw.up == wdw.lw) return(diagonalB_ng());
	int	m = a->right - a->left;
	int	n = b->right - b->left;
	int	k = wdw.lw - b->left + a->right;
	int	q = b->right - a->left - wdw.up;
	long	cvol =  m * n - (k * k + q * q) / 2;
	if (cvol < MaxVmfSpace || m == 1 || n <= 1) {
	    return (trcbkalignB_ng(wdw));
	}
	save_range(seqs, rng, 2);
	VTYPE	scr = centerB_ng(&m, &q, &n, &k, wdw, wdwl, wdwr);
	a->inex.exgr = b->inex.exgr = 0;
	int	r = b->left - a->left;
	if (r < wdwl.lw) b->left = a->left + wdwl.lw;
	if (r > wdwl.up) a->left = b->left - wdwl.up;
	a->right = m;
	b->right = n;
	lspB_ng(wdwl);
	rest_range(seqs, rng, 2);
	a->inex.exgr = aexg;
	b->inex.exgr = bexg;
	aexg = a->inex.exgl;
	bexg = b->inex.exgl;
	a->inex.exgl = b->inex.exgl = 0;
	r = b->right - a->right;
	if (r < wdwr.lw) a->right = b->right - wdwr.lw;
	if (r > wdwr.up) b->right = a->right + wdwr.up;
	a->left = q;
	b->left = k;
	lspB_ng(wdwr);
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexg;
	b->inex.exgl = bexg;
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

VTYPE HomScoreB_ng(const Seq* seqs[], const PwdB* pwd, long rr[])
{
	Aln2b1 alnv(seqs, pwd);
	WINDOW	wdw;
	stripe((const Seq**) seqs, &wdw, alprm.sh);
	return alnv.forwardB_ng(wdw, rr);
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

