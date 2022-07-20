/*****************************************************************************
*
*	Alignment of genomic vs. cDNA/mRNA sequences.
*	5' and 3' splice site signals and intron-length distribution
*	are considered.
*	Assumes there is no internal gap in either sequence.
*
*	Implement unidirectional Hirschberg linear-space algorithm
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
*	Osamu Gotoh, Ph.D.	(2003-2012)
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
#include "boyer_moore.h"
#if __SSE4_1__
#include "fwd2s1_simd.h"
#endif
#define	CigarM	0

static	const	int	expected_max_overlap = 1024;
static	const	int	expected_overlap_ext = 16;
static	const	int	max_dist2ss = 9;
static	const	CHAR	Newd = 8;
static	const	float	gudp = 0.3f;

static int infer_orientation(Seq** sqs, const PwdB* pwd);

class Aln2s1 {
protected:
const	Seq**	seqs;
const	Seq*&	a;
const	Seq*&	b;
const 	PwdB*	pwd;
const	bool	Local;
	Mfile*	mfd = 0;
	Vmf*	vmf = 0;
const	INT	lowestlvl;
const	int	end_margin;
const	VTYPE	slmt;
	bool	is3end = false;
const	int	Nod;
	VTYPE	backward[expected_max_overlap];
public:
	SpJunc*	spjcs;
	Cip_score*	cip;
	Aln2s1(const Seq** _seqs, const PwdB* _pwd);
	~Aln2s1() {delete mfd; delete spjcs; delete cip;}
	void	reset(const Seq** sqs) {
	    delete mfd; mfd = 0;
	    delete spjcs;
	    a = sqs[0];
	    b = sqs[1];
	    spjcs = new SpJunc(b, pwd);
	    is3end = false;
	}
	void	initS_ng(RVP* hh[], CHAR* hdir, const WINDOW& wdw, const RANGE* cutrng = 0);
	RVP*	lastS_ng(RVP* hh[], const WINDOW& wdw, const RANGE* cutrng = 0);
	VTYPE	forwardS_ng(const WINDOW& wdw, int pp[] = 0, 
		const RANGE* cutrng = 0);
	VTYPE	shortcutS_ng(int ovr, const BOUND& bab);
	void	hinitS_ng(Rvwml* hh[], const WINDOW& wdw);
	Rvwml*	hlastS_ng(Rvwml* hh[], const WINDOW& wdw);
	VTYPE	hirschbergS_ng(int cpos[], 		// coordinates on center row
		const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb, INT& bexgl);
	void	sinitS_ng(VTYPE* hh[], const WINDOW& wdw);
	VTYPE	slastS_ng(VTYPE* hh[], const WINDOW& wdw);
	VTYPE	scorealoneS_ng(const WINDOW& wdw);
	void	pfinitS_ng(RVP* hh[], const WINDOW& wdw);
	void	pbinitS_ng(RVP* hh[], const WINDOW& wdw);
	VTYPE	back2ward5endS_ng(int* ptr, const WINDOW& wdw);
	VTYPE	for2ward3endS_ng(int* ptr, const WINDOW& wdw);
	VTYPE	openendS_ng(int cmode);
	bool	indelfreespjS(int agap, VTYPE& iscr, const bool write_skl = true);
	VTYPE	diagonalS_ng(bool r_justt = false);
	VTYPE	trcbkalignS_ng(const WINDOW& wdw, bool spj = true,
		const RANGE* cutrng = 0);
	VTYPE	ordinaryS_ng(const WINDOW& wdw);
	VTYPE	backforth(int n, const BOUND& lub);
	VTYPE	lspS_ng(const WINDOW& wdw);
	int	nearest5ss(const BOUND& bab);
	int	nearest3ss(const BOUND& bab);
	VTYPE	failed(const RANGE* rng) {
	    rest_range(seqs, rng, 2);
	    return (NEVSEL);
	}
	int	first_exon_wmm(VTYPE& maxscr);
	VTYPE	first_exon(const BOUND& bab);
	int	last_exon_wmm(VTYPE& maxscr);
	VTYPE	last_exon(const BOUND& bab);
	VTYPE	micro_exon(const BOUND& bab);
	VTYPE	interpolateS(INT level, const int cmode, 
		const JUXT* wjxt, const BOUND& bab);
	VTYPE	seededS_ng(INT level, int cmode, const BOUND& lub);
	SKL*	globalS_ng(const WINDOW& wdw, VTYPE* scr);
	VTYPE	creepback(int ovr, VTYPE bscr, const BOUND& lub);
	VTYPE	creepfwrd(int& ovr, VTYPE bscr, const BOUND& lub);
	WLUNIT*	bestwlu(WLUNIT* wlu, const int nwlu, const int cmode);
};

Aln2s1::Aln2s1(const Seq* _seqs[], const PwdB* _pwd) :
	seqs(_seqs), a(seqs[0]), b(seqs[1]), pwd(_pwd), 
	Local(algmode.lcl & 16), lowestlvl(b->wllvl),
	end_margin((int) ((pwd->Vthr + pwd->BasicGOP) / pwd->BasicGEP)),
	slmt(pwd->Vthr / 2), Nod(2 * pwd->Noll - 1)
{
	spjcs = new SpJunc(b, pwd);
	cip = a->sigII? new Cip_score(a): 0;
}

void Aln2s1::initS_ng(RVP* hh[], CHAR* hdir, const WINDOW& wdw, const RANGE* cutrng)
{
	int	n = b->left;
	int	r = b->left - a->left;
	int	rr = b->right - a->left;
	int	cutlen = cutrng? cutrng->right - cutrng->left: 0;

	RVP*	h = hh[0] + r;
	CHAR*	hd = hdir? hdir + r: 0;
	h->val = 0;
	if (hd) *hd = 0;
	h->ptr = vmf->add(a->left, n, 0);
	if (a->inex.exgl) {		// semi-global
	    if (wdw.up < rr) rr = wdw.up;
	    while (++r <= rr) {
		if (cutrng && n == cutrng->left) {
		    n += cutlen;
		    r += cutlen;
		    continue;
		}
		++h;
		h->val = 0;
		if (hd) *++hd = 1;
		h->ptr = vmf->add(a->left, n, 0);
	    }
	}

	r = b->left - a->left;
	rr = b->left - a->right;
	h = hh[0] + r;
	if (hdir) hd = hdir + r;
	RVP*	f = hh[1] + r;
	if (wdw.lw > rr) rr = wdw.lw;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --f; if (hd) *--hd = 2;
	    if (b->inex.exgl) {		// local
		h->val = 0;
		h->ptr = vmf->add(a->left + i, b->left, 0);
	    } else {			// (semi) global
		*h = h[1];
		if (i == 1) {
		    h->val += pwd->GapPenalty(1);
		    *f = *h;
		} else {
		    *f = f[1];
		    h->val += pwd->GapExtPen(i);
		    f->val += pwd->BasicGEP;
		}
	    }
	}
}

RVP* Aln2s1::lastS_ng(RVP* hh[], const WINDOW& wdw, const RANGE* cutrng)
{
	int	cutlen = cutrng? cutrng->right - cutrng->left: 0;
	int	rw = wdw.lw;
	int	rf = (cutrng? cutrng->left: b->left) - a->right;
	if (rf > rw) rw = rf;
	RVP*	h = hh[0] + rw - cutlen;
	RVP*	h9 = hh[0] + b->right - a->right - cutlen;
	RVP*	mx = h9;

	if (a->inex.exgr) {
	    for ( ; h <= h9; ++h)
		if (h->val > mx->val) mx = h;
	}
	if (b->inex.exgr) {
	    rw = std::min(wdw.up, b->right - a->left) - cutlen;
	    for (RVP* h = hh[0] + rw; h > h9; --h)
		if (h->val > mx->val) mx = h;
	}
const	int	i = mx - h9;
	int	m9 = a->right;
	int	n9 = b->right;
	if (i > 0) m9 -= i;
	if (i < 0) n9 += i;
	mx->ptr = vmf->add(m9, n9, mx->ptr);
	return (mx);
}

VTYPE Aln2s1::forwardS_ng(const WINDOW& wdw, int* pp, const RANGE* cutrng)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
const	bool	spj = b->inex.intr;
	RVP*	hh[NOL];	// H matrix
	RVP	e1, e2;
	RVP*	hf[NOD] = {0, &e1, 0, &e2, 0};	// [DIAG, HORI, VERT, HORL, VERTL]
	RVPDJ	rcd[NCAND + 1]; // [candidates]
	int	idx[NCAND + 1]; // [candidates]
const	RVPDJ*	maxphl[NOD];
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
const	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	int	cutlen = cutrng? cutrng->right - cutrng->left: 0;
const	VTYPE	longgep = pwd->BasicGEP * cutlen;
const	VTYPE	longgep2 = pwd->LongGEP * cutlen;

const	int	width = wdw.width - cutlen;
const	size_t	bufsiz = pwd->Noll * width;
	RVP*	jbuf = new RVP[bufsiz];
	vset(jbuf, black_vp, bufsiz);
	RVP*	blackvp = jbuf + bufsiz - 1;
	hh[0] = jbuf - wdw.lw + 1;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + width;
	CHAR*	dbuf = new CHAR[width];
	vclear(dbuf, width);
	CHAR*	hdir = dbuf - wdw.lw + 1;
	vmf->add(0, 0, 0);	// Skip 0-th record
	initS_ng(hh, hdir, wdw, cutrng);

	int	m = a->left;
	if (!a->inex.exgl) --m; 	// global
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as, ++n1, ++n2) {
const	    bool	internal = spj && (!a->inex.exgr || m < a->right);
const	    VTYPE	sigB = cip? cip->cip_score(m): 0;
	    int	n = std::max(n1, b->left);
const	    int	n9 = std::min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    int k = 0;
	    RVP*&	h = hf[0] = hh[k] + r;
	    RVP*&	f = hf[2] = hh[++k] + r;
	    RVP*&	f2 = hf[4] = dagp? hh[++k] + r: blackvp;
	    CHAR*	dir = hdir + r;
	    CHAR	psp = 0;
	    e1 = e2 = black_vp;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    vset((RVPDJ*) rcd, black_vpdj, NCAND + 1);
	    for (int l = 0; l <= NCAND; ++l) idx[l] = l;
	    int	ncand = -1;
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d %2d", m, n, *dir);
		putvar(h->val); putchar('\n');
	    }
#endif
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x;
		++dir; ++h; ++f; if (dagp) ++f2;
//	Diagonal
		RVP*	from = h;
		RVP*	mx = h;
const		VTYPE	diag = h->val;
		if (m == a->left) goto HorizonS;
		h->val += qprof[*bs];
		*dir = (*dir % Newd)? Newd: 0;

//	Vertical
		x = (++from)->val + pwd->BasicGOP;
		if (x >= f[1].val) {
		    f->val = x;
		    f->ptr = from->ptr;
		} else	*f = f[1];
		f->val += pwd->BasicGEP;
		if (f->val > mx->val) mx = f;

//	Vertical2
		if (dagp) {
		  x = from->val + pwd->LongGOP;
		  if (x >= f2[1].val) {
		    f2->val = x;
		    f2->ptr = from->ptr;
		  } else	*f2 = f2[1];
		  f2->val += pwd->LongGEP;
		  if (f2->val > mx->val) mx = f2;
		}
HorizonS:
//	Horizontal
		x = h[-1].val + pwd->BasicGOP;
const		CHAR	prev_psp = psp;
		if (x >= e1.val) {
		    e1.val = x;
		    e1.ptr = h[-1].ptr;
		    psp = psp? e1_psp: 0;
		} else
		    psp &= e1_psp;
		e1.val += pwd->BasicGEP;
		if (e1.val >= mx->val) mx = &e1;

//	Horizontal2
		if (dagp) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= e2.val) {
			e2.val = x;
			e2.ptr = h[-1].ptr;
			if (prev_psp) psp |= e2_psp;
		    } else
			psp |= (prev_psp & e2_psp);
		    e2.val += pwd->LongGEP;
		    if (e2.val >= mx->val) mx = &e2;
		}

//	intron 3' boundary, assume no overlapping signals
		if (internal && b->exin->isAccpt(n)) {
		    vclear(maxphl, Nod);
		    for (int l = 0; l <= ncand; ++l) {
const			RVPDJ*	prd = rcd + idx[l];
			if (n - prd->jnc < IntronPrm.llmt) continue;
			x = prd->val + sigB + spjcs->spjscr(prd->jnc, n);
			from = hf[prd->dir];
			if (x >= from->val) {
			    from->val = x;
			    maxphl[prd->dir] = prd;
			}
		    }
		    for (int k = 0; k < Nod; ++k) {
const			RVPDJ*	prd = maxphl[k];
			if (!prd) continue;
			from = hf[k];
			psp |= psp_bit[k];
			from->ptr = vmf->add(m, n, 
			    vmf->add(m, prd->jnc, prd->ptr));
			if (from->val >= mx->val) mx = from;
		    }
		}

//	Find optimal path
#if DEBUG
		VTYPE	y = h->val;
#endif
		int	hd = 0;
		if (h != mx) {		// non-diagonal
		    *h = *mx;
		    while (mx != hf[++hd]) ;
		    *dir = hd;
		} else if (Local && h->val > diag) {
		    if (LocalL && diag == 0)
			h->ptr = vmf->add(m - 1, n - 1, 0);
		    else if (LocalR && h->val > maxh.val) {
			maxh.val = h->val;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && h->val <= 0) {h->val = 0; *dir = 1;}
		else if (*dir == Newd && !(psp & psp_bit[0]))
		    h->ptr = vmf->add(m - 1, n - 1, h->ptr);

//	intron 5' boundary
		if (internal && b->exin->isDonor(n)) {
const		    VTYPE	sigJ = b->exin->sig53(n, 0, IE5);
		    for (int k = hd == 0? 0: 1; k < Nod; ++k) {
			from = hf[k];
			if (psp & psp_bit[k]) continue;	// disallow orphan exon
			if (k != hd) {
			    VTYPE	z = mx->val;
			    if (hd == 0 || (k - hd) % 2) z += pwd->GOP[k / 2];
			    if (from->val <= z) continue;	// prune
			}
			x = from->val + sigJ;
			int	l = ncand < NCAND? ++ncand: NCAND;
			while (--l >= 0) {
			    if (x > rcd[idx[l]].val)
				std::swap(idx[l], idx[l + 1]);
			    else
				break;
			}
			if (++l < NCAND) {
			    RVPDJ*	prd = rcd + idx[l];
			    prd->val = x;
			    prd->jnc = n;
			    prd->dir = k;
			    prd->ptr = from->ptr;
			} else --ncand;
		    }
		}

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d %2d ", m, n, *dir);
		    putvar(mx->val); putvar(y); 
		    putvar(f->val); putvar(e1.val);
		    if (dagp) {
			putvar(f2->val); putvar(e2.val);
		    }
		    if (algmode.lsg) {
			putvar(rcd[idx[0]].val);
			putvar(rcd[idx[1]].val);
		    }
		    putchar('\n');
		}
#endif
		if (cutrng && n == cutrng->left) {	// shortcut
		    e1.val += longgep;
		    if (dagp) e2.val += longgep2;
		    *h = dagp? e2: e1;
		    *f = black_vp;
		    n += cutlen;
		    bs += cutlen;
		}
	    } // end of n-loop
	} // end of m-loop

	if (LocalR) {
	    *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    RVP*	mx = lastS_ng(hh, wdw, cutrng);
	    maxh.val = mx->val;
	    if (pp) *pp = mx->ptr;
	}
	delete[] jbuf;
	delete[] dbuf;
	return (maxh.val);
}

VTYPE skl_rngS_ng(const Seq* seqs[], Gsinfo* gsi, const PwdB* pwd)
{
const	Seq*	a = seqs[0];
const	Seq*	b = seqs[1];
const	Exinon*	exin = b->exin;
const	SKL*	wsk = gsi->skl;
	int	num = (wsk++)->n;
	SpJunc*	spjcs = new SpJunc(b, pwd);
	VTYPE	h = 0;
	VTYPE	ha = 0;
	VTYPE	hb = 0;
	VTYPE	xi = 0;
	VTYPE	sig5 = 0;
	VTYPE	sig3 = 0;
	int	insert = 0;
	int	deletn = 0;
	int	intlen = 0;
	int	preint = 0;
	EISCR	rbuf;
	FSTAT*	fst = &gsi->fstat;
	FSTAT	pst;
	int	psp = 0;		// post splicing position
	Eijnc*	eijnc = gsi->eijnc = new Eijnc(true);
	Cigar*	cigar = gsi->cigar = 0;
	Vulgar*	vlgar = gsi->vlgar = 0;
	Samfmt*	samfm = gsi->samfm = 0;
	switch (algmode.nsa) {
	    case CIG_FORM: cigar = gsi->cigar = new Cigar(); break;
	    case VLG_FORM: vlgar = gsi->vlgar = new Vulgar(); break;
	    case SAM_FORM: samfm = gsi->samfm = new Samfmt(); break;
	    default:	 break;
	}
	vclear(fst);
	vclear(&pst);
	vclear(&rbuf);
	if (wsk[1].n == wsk->n && b->inex.exgl) {++wsk; --num;}
	int	m = wsk->m;
	int	n = wsk->n;
const	CHAR*	as = a->at(m);
const	CHAR*	bs = b->at(n);
	PfqItr	api(a, m);
const	bool	usespb = api.size() && use_spb();
	rbuf.left = n;
	rbuf.rleft = m;
	rbuf.iscr = NEVSEL;
	rbuf.sig3 = 0;
	if (samfm) {
	    samfm->left = m;
	    if (m) samfm->push('H', m);	// local alignment
	    if (b->inex.sens == 0) samfm->pos = n;
	}
	while (--num) {
	    ++wsk;
const	    int	mi = wsk->m - m;
	    if (mi && insert) {
const		bool	j = a->inex.exgl && m == a->left;
		VTYPE	x = j? 0: pwd->GapPenalty(insert);
		VTYPE	xi = NEVSEL;
		if (intlen) {
		    insert -= intlen;
		    xi = rbuf.iscr + pwd->GapPenalty(insert);
		}
		if (xi >= x) {		// intron
		    if (cigar) {
			if (preint) cigar->push('D', preint);
			cigar->push('N', intlen);
		    }
		    if (samfm) {
			if (preint) samfm->push('D', preint);
			samfm->push('N', intlen);
		    }
		    if (vlgar) {
			if (preint) vlgar->push('G', 0, preint);
			vlgar->push('5', 0, 2);
			vlgar->push('I', 0, intlen - 4);
			vlgar->push('3', 0, 2);
		    }
		    hb = ha;
		    if (eijnc && rbuf.right - rbuf.left > 0)
			eijnc->push(&rbuf);
		    rbuf.left = rbuf.right + intlen;
		    rbuf.rleft = m;
		    rbuf.sig3 = sig3;
		    rbuf.iscr = NEVSEL;
		    h += xi;
		    insert -= preint;
		} else {		// gap
		    h += x;
		}
		if (insert) {
		    if (cigar) cigar->push('D', insert);
		    if (samfm) samfm->push('D', insert);
		    if (vlgar) vlgar->push('G', 0, insert);
		    insert = intlen = preint = 0;
		}
	    }
const	    int	ni = wsk->n - n;
	    if (ni && deletn) {
		if (!(b->inex.exgl && n == b->left)) {
		    h += pwd->GapPenalty(deletn);
		    fst->gap += 1;
		}
		as += deletn;
		deletn = 0;
	    }
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    if (d) {
		m += d;
		if (cigar) cigar->push('M', d);
		if (vlgar) vlgar->push('M', d, d);
#if CigarM
		if (samfm) samfm->push('M', d);
#endif
		VTYPE	x = 0;
#if CigarM
		for ( ; d; --d, ++as, ++bs, ++n) {
		    if (eijnc) eijnc->shift(rbuf, *fst, 
			++pst == alprm2.jneibr);
		    x += pwd->sim2(as, bs);
		    if (*as == *bs) ++fst->mch;
		    else ++fst->mmc;
		}
#else
		int	run = 0;
		for ( ; d; --d, ++as, ++bs, ++n) {
		  if (eijnc) eijnc->shift(rbuf, *fst, 
		    psp++ == alprm2.jneibr);
		  x += pwd->sim2(as, bs);
		  if (*as == *bs) {
		    if (run < 0) {
			if (samfm) samfm->push('X', -run);
			run = 0;
		    }
		    ++fst->mch;
		    ++run;
		  } else {
		    if (run > 0) {
			if (samfm) samfm->push('=', run);
			run = 0;
		    }
		    ++fst->mmc;
		    --run;
		  }
		  if (samfm) {
		    if (run > 0) samfm->push('=', run); else
		    if (run < 0) samfm->push('X', -run);
		  }
		}
#endif
		h += x;
		fst->val += x;
	    }
	    if (i > 0) {
		deletn += i;
		if (cigar) cigar->push('I', i);
		if (samfm) samfm->push('I', i);
		if (vlgar) vlgar->push('G', i, 0);
		for (int j = 0; j < i; ++j) {
		    if (eijnc)
			eijnc->shift(rbuf, *fst, psp++ == alprm2.jneibr);
		    ++fst->unp;
		}
	    } else if (i < 0) {
const		int	n3 = n + (i = -i);
		if (algmode.lsg && i > IntronPrm.minl) {
		    sig5 = exin->sig53(n, n3, IE5);
		    sig3 = exin->sig53(n, n3, IE3);
		    xi = sig5 + spjcs->spjscr(n, n3);
		    if (usespb) xi += api.match_score(m);
		} else	xi = NEVSEL;
		if (xi > pwd->GapPenalty(i) && xi > rbuf.iscr) {
		    preint = insert;		// intron
		    intlen = i;
		    rbuf.right = n;
		    rbuf.rright = m;
		    rbuf.iscr = xi;
		    rbuf.escr = h + sig5 - hb;
		    rbuf.sig5 = sig5;
		    ha = h + xi - sig3;
		    if (eijnc) {
			eijnc->store(rbuf, *fst, pst,
			psp < alprm2.jneibr);
			pst = *fst;
			psp = 0;
		    }
		} else if (!a->inex.exgl || m != a->left) {
		    if (!insert) fst->gap += 1;
		    for (int j = 0; j < i; ++j, ++n) {
			if (eijnc) eijnc->shift(rbuf, *fst, 
			    psp++ == alprm2.jneibr);
		    	++fst->unp;
		    }
		}
		bs += i;
		insert += i;
    	    }
	    m = wsk->m;
	    n = wsk->n;
	    if (usespb) while (!api.end() && api < m) ++api;
	}
	if (insert && !(a->inex.exgr && m == a->right)) {
	    h += pwd->GapPenalty(insert);
	    fst->gap += 1;
	    fst->unp += insert;
	    if (cigar) cigar->push('D', insert);
	    if (samfm) samfm->push('D', insert);
	    if (vlgar) vlgar->push('G', 0, insert);
	}
	if (deletn && !(b->inex.exgr && n == b->right)) {
	    h += pwd->GapPenalty(deletn);
	    fst->gap += 1;
	    fst->unp += deletn;
	    if (cigar) cigar->push('I', deletn);
	    if (samfm) samfm->push('I', deletn);
	    if (vlgar) vlgar->push('G', deletn, 0);
	}
	if (eijnc) {
	    rbuf.escr = h - hb;
	    rbuf.iscr = 0;
	    rbuf.sig5 = 0;
	    rbuf.right = n;
	    rbuf.rright = m;
	    eijnc->store(rbuf, *fst, pst, n - rbuf.left <= alprm2.jneibr);
	    eijnc->push(&rbuf);
	    rbuf.left = endrng.left;
	    rbuf.right = endrng.right;
	    eijnc->push(&rbuf);
	    eijnc->flush();
	    gsi->noeij = eijnc->size() - 1;
	}
	if (cigar) cigar->flush();
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
	if (vlgar) vlgar->flush();
	delete spjcs;
	fst->val += pwd->BasicGOP * fst->gap + pwd->BasicGEP * fst->unp;;
	return (h);
}

/********************************************************
*
*	unidirectional Hirschberg algorithm
*
********************************************************/

void Aln2s1::hinitS_ng(Rvwml* hh[], const WINDOW& wdw)
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

Rvwml* Aln2s1::hlastS_ng(Rvwml* hh[], const WINDOW& wdw)
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

VTYPE Aln2s1::hirschbergS_ng(int cpos[], 
		const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb, INT& bexgl)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	Rvwml*	hh[NOL];
	Rvwml*	hf[NOD];
	Rvdwmlj	rcd[NCAND + 1];;		// [candidates]
	int	idx[NCAND + 1];		// [candidates]
const	Rvdwmlj*	maxphl[NOD];
	Rvwml	e[2];
	Rvwml*&	e1 = hf[1] = e;
	Rvwml*&	e2 = hf[3] = e + 1;
	int*	center_lnk[NOD];
	int*	center_end[NOL];
	int*	center_lwr[NOL];
	int*	center_upr[NOL];
	Rvwmrmn	maxh = {NEVSEL, 0, 0, 0, 0};
const	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
#if MONITOR
	long	start = time(0);
	long	stop;
#endif

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
	vset(ibuf, end_of_ulk, 2 * bufsiz);	// center_lnk + center_end

	hinitS_ng(hh, wdw);

	int	m = a->left;
const	int	mm = (a->left + a->right + 1) / 2;
	if (!a->inex.exgl) --m; 	// global
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as, ++n1, ++n2) {
const	    VTYPE	sigB = cip? cip->cip_score(m): 0;
	    int	n = std::max(n1, b->left);
const	    int	n9 = std::min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    CHAR	psp = 0;
	    r = n - m;
	    Rvwml*&	h = hf[0] = hh[0] + r;
	    Rvwml*&	f = hf[2] = hh[1] + r;
	    Rvwml*&	f2 = hf[4] = dagp? hh[2] + r: blackvwu;
	    vset(e, black_vpwml, 2);
	    vset(rcd, black_vdwmlj, NCAND + 1);
	    for (int l = 0; l <= NCAND; ++l) idx[l] = l;
	    int	ncand = -1;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d ", m, n);
		putvar(h->val); putchar('\n');
	    }
#endif
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x;
		++r;
		for (int k = 0; k < pwd->Noll; ++k) {
const		    int	kk = k + k;
		    ++hf[kk];
		    if (m == mm) {
			center_end[k][r] = hf[kk]->ulk;
			hf[kk]->ulk = k? end_of_ulk: r;
		    }
		}

//	Diagonal
		Rvwml*	from = h;
		Rvwml*	mx = h;
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
const		CHAR	prev_psp = psp;
		if (x >= e1->val) {
		    *e1 = h[-1];
		    e1->val = x;
		    psp = psp? e1_psp: 0;
		} else
		    psp &= e_psp;
		e1->val += pwd->BasicGEP;
		if (e1->val >= mx->val) mx = e1;

//	Horizontal2
		if (dagp) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= e2->val) {
			*e2 = h[-1];
			e2->val = x;
			if (prev_psp) psp |= e2_psp;
		    } else
			psp |= (prev_psp & e2_psp);
		    e2->val += pwd->LongGEP;
		    if (e2->val >= mx->val) mx = e2;
		}

//	intron 3' boundary
		if (b->exin->isAccpt(n)) {
		    vclear(maxphl, Nod);
		    for (int l = 0; l <= ncand; ++l) {
const			Rvdwmlj*	prd = rcd + idx[l];
			if (n - prd->jnc < IntronPrm.llmt) continue;
			from = hf[prd->dir];
			x = prd->val + sigB + spjcs->spjscr(prd->jnc, n);
			if (x > from->val) {
			    from->val = x;
			    maxphl[prd->dir] = prd;
			}
		    }
		    for (int k = 0; k < Nod; ++k) {
const			Rvdwmlj*	prd = maxphl[k];
			if (!prd) continue;
			psp |= psp_bit[k];
			from = hf[k];
			from->upr = std::max(prd->upr, r);
			from->lwr = std::min(prd->lwr, r);
			from->ulk = (m == mm)? r: prd->ulk;
			if (from->val > mx->val) mx = from;
		    }
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

//	intron 5' boundary
		if (b->exin->isDonor(n)) {
const		    VTYPE	sigJ = b->exin->sig53(n, 0, IE5);
		    for (int k = mx == h? 0: 1; k < Nod; ++k) {
			from = hf[k];
			if (psp & psp_bit[k]) continue;	// disallow orphan exon
			if (k != hd) {
			    y = mx->val;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (from->val <= y) continue;	// prune
			}
			x = from->val + sigJ;
			int	l = ncand < NCAND? ++ncand: NCAND;
			while (--l >= 0) {
			    if (x > rcd[idx[l]].val)
				std::swap(idx[l], idx[l + 1]);
			    else
				break;
			}
			if (++l < NCAND) {
			    Rvdwmlj*	prd = rcd + idx[l];
			    prd->val = x;
			    prd->jnc = n;
			    prd->dir = k;
			    prd->upr = from->upr;
			    prd->lwr = from->lwr;
			    prd->ulk = (m == mm)? r: from->ulk;
			} else --ncand;
		    }
		}

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d ", m, n);
		    putvar(mx->val); putvar(y); 
		    putvar(f->val); putvar(e1->val);
		    if (dagp) {
			putvar(f2->val); putvar(e2->val);
		    }
		    if (algmode.lsg) {
			putvar(rcd[idx[0]].val);
			putvar(rcd[idx[1]].val);
			putvar(rcd[idx[2]].val);
		    }
		    putchar('\n');
		}
#endif

// save variables at the center
		if (m == mm) {
const		    int	vk = (hd == 2 || hd == 4)? hd / 2: 0;
		    for (int k = 0; k < pwd->Noll; ++k) {
const			int	kk = k + k;
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
const	    Rvwml*	mx = hlastS_ng(hh, wdw);
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
	    r = wdwf.lw = center_lwr[0][r];
	} else	cpos[0] = end_of_ulk;		// don't cross center
	if (LocalL) {
	    a->left = maxh.ml;
	    b->left = r + maxh.ml;
	} else {
	    if (a->inex.exgl && rl > r) a->left = b->left - r;
	    if (b->inex.exgl && rl < r) b->left = a->left + r;
	    r = b->left - a->left;
	    if (r < wdwf.lw) wdwf.lw = r;
	    if (r > wdwf.up) wdwf.up = r;
	}
	bexgl = (d > 0)? 1: 0;
	delete[] wbuf;
	delete[] ibuf;
	return (maxh.val);
}

/********************************************************
*
*	score only algorithm
*
********************************************************/

void Aln2s1::sinitS_ng(VTYPE* hh[], const WINDOW& wdw)
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

VTYPE Aln2s1::slastS_ng(VTYPE* hh[], const WINDOW& wdw)
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

VTYPE Aln2s1::scorealoneS_ng(const WINDOW& wdw)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	VTYPE*	hh[NOL];
	RVDJ	rcd[NCAND + 1];
	int	idx[NCAND + 1];		// [candidates]
const	RVDJ*	maxphl[NOD];
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
	sinitS_ng(hh, wdw);

	int	m = a->left;
	if (!a->inex.exgl) --m; 	// global
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as, ++n1, ++n2) {
const	    VTYPE	sigB = cip? cip->cip_score(m): 0;
	    int	n = std::max(n1, b->left);
const	    int	n9 = std::min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    CHAR	psp = 0;
	    int	r = n - m;
	    VTYPE*&	h = hf[0] = hh[0] + r;
	    VTYPE*&	f = hf[2] = hh[1] + r;
	    VTYPE*&	f2 = hf[4] = dagp? hh[2] + r: blackv;
	    e1 = e2 = NEVSEL;
	    vset(rcd, black_vdj, NCAND + 1);
	    for (int l = 0; l <= NCAND; ++l) idx[l] = l;
	    int	ncand = -1;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d ", m, n);
		putvar(*h); putchar('\n');
	    }
#endif
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x;
		++h; ++f; if (dagp) ++f2;

//	Diagonal
		VTYPE*	from = h;
		VTYPE*	mx = h;
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
const		CHAR	prev_psp = psp;
		if (x > e1) {
		    e1 = x;
		    psp = psp? e1_psp: 0;
		} else
		    psp &= e1_psp;
		e1 += pwd->BasicGEP;
		if (e1 > *mx) mx = &e1;

//	Horizontal2
		if (dagp) {
		    x = h[-1] + pwd->LongGOP;
		    if (x > e2) {
			e2 = x;
			if (prev_psp) psp |= e2_psp;
		    } else
			psp |= (prev_psp & e2_psp);
		    e2 += pwd->LongGEP;
		    if (e2 > *mx) mx = &e2;
		}

//	intron 3' boundary
		if (b->exin->isAccpt(n)) {
		    vclear(maxphl, Nod);
		    for (int l = 0; l <= ncand; ++l) {
const			RVDJ*	prd = rcd + idx[l];
			if (n - prd->jnc < IntronPrm.llmt) continue;
			from = hf[prd->dir];
			x = prd->val + sigB + spjcs->spjscr(prd->jnc, n);
			if (x > *from) {
			    *from = x;
			    maxphl[prd->dir] = prd;
			}
		    }
		    for (int k = 0; k < Nod; ++k) {
const			RVDJ*	prd = maxphl[k];
			if (!prd) continue;
			psp |= psp_bit[k];
			from = hf[k];
			if (*from > *mx) mx = from;
		    }
		}

//	Find optimal path
		VTYPE	y = *h;
		if (h != mx) *h = *mx;
		else if (LocalR && y > maxh) maxh = y;
		if (LocalL && *h < 0) *h = 0;
		int	hd = 0;
		for ( ; mx != hf[hd]; ++hd) ;

//	intron 5' boundary
		if (b->exin->isDonor(n)) {
const		    VTYPE	sigJ = b->exin->sig53(n, 0, IE5);
		    for (int k = hd == 0? 0: 1; k < Nod; ++k) {
			from = hf[k];
			if (psp & psp_bit[k]) continue;
			if (k != hd) {
			    y = *mx;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (*from <= y) continue;	// prune
			}
			x = *from + sigJ;
			int	l = ncand < NCAND? ++ncand: NCAND;
			while (--l >= 0) {
			    if (x >= rcd[idx[l]].val)
				std::swap(idx[l], idx[l + 1]);
			    else
				break;
			}
			if (++l < NCAND) {
			    RVDJ*	prd = rcd + idx[l];
			    prd->val = x;
			    prd->jnc = n;
			    prd->dir = k;
			} else --ncand;
		    }
		}

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d ", m, n);
		    putvar(*mx); putvar(y); 
		    putvar(*f); putvar(e1);
		    if (dagp) {
			putvar(*f2); putvar(e2);
		    }
		    if (algmode.lsg) {
			putvar(rcd[0].val);
			putvar(rcd[1].val);
			putvar(rcd[2].val);
		    }
		    putchar('\n');
		}
#endif
	    }	// end of n-loop
	}	// end of m-loop
	if (!LocalR) maxh = slastS_ng(hh, wdw);

	delete[] wbuf;
	return (maxh);
}

void Aln2s1::pfinitS_ng(RVP* hh[], const WINDOW& wdw)
{
	int	r = b->left - a->left;
	int	rr = b->right - a->left;

	RVP*	h = hh[0] + r;
	RVP*	f = hh[1] + r;
	h->val = 0;
	h->ptr = vmf->add(a->left, b->left, 0);

	r = b->left - a->left;
	rr = b->left - a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --f;
	    *h = h[1];
	    if (i == 1)	h->val += pwd->BasicGOP;
	    h->val += pwd->BasicGEP;
	    *f = *h;
	}
}

void Aln2s1::pbinitS_ng(RVP* hh[], const WINDOW& wdw)
{
	int	r = b->right - a->right;
	int	rr = b->left - a->right;

	RVP*	h = hh[0] + r;
	RVP*	f = hh[1] + r;
	h->val = 0;
	h->ptr = vmf->add(a->right, b->right, 0);

	r = b->right - a->right;
	rr = b->right - a->left;
	if (wdw.up < rr) rr = wdw.up;
	for (int i = 1; ++r <= rr; ++i) {
	    ++h; ++f;
	    *h = h[-1];
	    if (i == 1) h->val += pwd->BasicGOP;
	    h->val += pwd->BasicGEP;
	    *f = *h;
	}
}

// intron-less backward extension

VTYPE Aln2s1::back2ward5endS_ng(int *ptr, const WINDOW& wdw)
{
	RVP*	hh[2];			// H matrix
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	RVP	e1 = black_vp;
	RVP*	hf[3] = {0, &e1, 0};
	RVP*&	h = hf[0];
	RVP*&	f = hf[2];

const	size_t	bufsiz = 2 * wdw.width;
	RVP*	jbuf = new RVP[bufsiz];
	vset(jbuf, black_vp, bufsiz);
	RVP*	blackvp = jbuf + bufsiz - 1;
	hh[0] = jbuf - wdw.lw + 1;
	hh[1] = hh[0] + wdw.width;
	CHAR*	dbuf = new CHAR[wdw.width];
	vset(dbuf, CHAR(1), wdw.width);
	CHAR*	hdir = dbuf - wdw.lw + 1;
	vmf->add(0, 0, 0);		// Skip 0-th record
	pbinitS_ng(hh, wdw);

	int	m = a->right;
	if (!a->inex.exgr) ++m; 	// global
const	CHAR*	as = a->at(m);
	int	n = 0;
	int	n9 = 0;
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
	while (--m >= a->left) {
	    --as; --n1; --n2;
	    n = std::min(n2, b->right);
	    n9 = std::max(n1, b->left);
	    int	r = n - m;
	    int	nr = n + 1;
	    int	peak = 0;
const	    CHAR*	bs = b->at(n);
	    h = hh[0] + r;
	    f = hh[1] + r;
	    CHAR*	dir = hdir + r;
	    e1 = black_vp;
	    RVP*	mxd = ((h->val + pwd->Vthr) < maxh.val)? blackvp: h;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)

#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m+1, n, *dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    while (--n >= n9) {
		VTYPE	x;
		--bs; --r; --dir; --h; --f;

//	Diagonal
		RVP*	from = h;
		RVP*	mx = h;
		if (m == a->right) goto HorizonB;
		h->val += qprof[*bs];
		*dir = *dir % Newd? Newd: 0;

//	Vertical
		x = (--from)->val + pwd->BasicGOP;
		if (x >= f[-1].val) {
		    *f = *from;
		    f->val = x;
		} else *f = f[-1];
		f->val += pwd->BasicGEP;
		if (f->val >= mx->val) mx = f;
HorizonB:
//	Horizontal
		x = h[1].val + pwd->BasicGOP;
		if (x >= e1.val) {
		    e1 = h[1];
		    e1.val = x;
		}
		e1.val += pwd->BasicGEP;
		if (e1.val >= mx->val) mx = &e1;

//	Find optimal path
#if DEBUG
		VTYPE	y = maxh.val;
#endif
		if ((*dir & Newd) && vmf)
		    mx->ptr = vmf->add(m + 1, n + 1, mx->ptr);
		if (mx->val > maxh.val) {
		    maxh.val = mx->val;
		    maxh.p = mx->ptr;
		    maxh.m = m;
		    maxh.n = n;
		}
		if (mx->val + pwd->Vthr < maxh.val) {	// x drop off
		    if (peak) {				// end of block
			n1 = n + 1;
			peak = 0;
		    }
		    nr = n;				// before block
		} else if ((*dir % Newd == 0) && mx->val >= mxd->val) {
		    mxd = mx;
		    if (nr < n2) n2 = nr;
		    peak = 1;
		}
		if (h != mx) *h = *mx;
		for (*dir = 0; mx != hf[*dir]; ++*dir) ;

#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m, n, *dir);
		putvar(mx->val); putvar(y); 
		putvar(f->val); putvar(e1.val);
		putchar('\n');
	}
#endif
	    } // end of n loop
	    if (mxd->val <= NEVSEL) break;	// no peak
	    if (peak) n1 = n;
	} // end of m loop

	*ptr = vmf->add(maxh.m, maxh.n, maxh.p);
	delete[] jbuf;
	delete[] dbuf;
	return (maxh.val);
}

// intron-less forward extension

VTYPE Aln2s1::for2ward3endS_ng(int *ptr, const WINDOW& wdw)
{
	RVP*	hh[2];		// H matrix
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	RVP	e1 = black_vp;
	RVP*	hf[3] = {0, &e1, 0};
	RVP*&	h = hf[0];
	RVP*&	f = hf[2];

const	size_t	bufsiz = 2 * wdw.width;
	RVP*	jbuf = new RVP[bufsiz];
	vset(jbuf, black_vp, bufsiz);
	RVP*	blackvp = jbuf + bufsiz - 1;
	hh[0] = jbuf - wdw.lw + 1;
	hh[1] = hh[0] + wdw.width;
	CHAR*	dbuf = new CHAR[wdw.width];
	vset(dbuf, CHAR(1), wdw.width);
	CHAR*	hdir = dbuf  - wdw.lw + 1;
	vmf->add(0, 0, 0);		// Skip 0-th record
	pfinitS_ng(hh, wdw);

	int	m = a->left;
	if (!a->inex.exgl) --m;		// global
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
	int	n = 0;
	int	n9 = 0;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as, ++n1, ++n2) {
	    n = std::max(n1, b->left);
	    n9 = std::min(n2, b->right);
	    int	nr = n - 1;
	    int	peak = 0;
	    int	r = n - m;
const	    CHAR*	bs = b->at(n);
	    h = hh[0] + r;;
	    f = hh[1] + r;
	    CHAR*	dir = hdir + r;
	    RVP*	mxd = ((h->val + pwd->Vthr) < maxh.val)? blackvp: h;
	    e1 = black_vp;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)

#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m, n, *dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x;
		++r; ++dir, ++h; ++f;

//	Diagonal
		RVP*	from = h;
		RVP*	mx = h;
		if (m == a->left) goto HorizonP;
		h->val += qprof[*bs];
		*dir = (*dir % Newd)? Newd: 0;

//	Vertical
		x = (++from)->val + pwd->BasicGOP;
		if (x >= f[1].val) {
		    *f = *from;
		    f->val = x;
		} else	*f = f[1];
		f->val += pwd->BasicGEP;
		if (f->val >= mx->val) mx = f;
HorizonP:
//	Horizontal
		x = h[-1].val + pwd->BasicGOP;
		if (x >= e1.val) {
		    e1 = h[-1];
		    e1.val = x;
		}
		e1.val += pwd->BasicGEP;
		if (e1.val >= mx->val) mx = &e1;
//	Find optimal path
#if DEBUG
		VTYPE	y = h->val;
#endif
		if (*dir & Newd)
		    mx->ptr = vmf->add(m - 1, n - 1, mx->ptr);
		if (mx->val > maxh.val) {
		    maxh.val = mx->val;
		    maxh.m = m;
		    maxh.n = n;
		    maxh.p = mx->ptr;
		}
		if (mx->val + pwd->Vthr < maxh.val) {	// x drop off
		    if (peak) {				// end of block
			n2 = n - 1;
			peak = 0;
		    }
		    nr = n;
		} else if ((*dir % Newd == 0) && mx->val >= mxd->val) {
		    mxd = mx;
		    if (nr > n1) n1 = nr;
		    peak = 1;
		}
		if (h != mx) *h = *mx;
		for (*dir = 0; mx != hf[*dir]; ++*dir) ;

#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m, n, *dir);
		putvar(mx->val); putvar(y); 
		putvar(f->val); putvar(e1.val);
		putchar('\n');
	}
#endif
	    } // end of n-loop
	    if (peak) n2 = n;
	    if (mxd->val <= NEVSEL) break;
	} // end of m-loop

	*ptr = vmf->add(maxh.m, maxh.n, maxh.p);
	is3end = true;
	delete[] jbuf;
	delete[] dbuf;
	return (maxh.val);
}

VTYPE Aln2s1::diagonalS_ng(bool r_just)
{
const	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	int	dlt = Local? 0: ((b->right - b->left) - (a->right - a->left));
	if (dlt < 0) std::swap(a, b);
	int	mL = a->left;
	int	mR = a->right;
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(b->left);
	if (r_just && !Local) bs += abs(dlt);
	VTYPE	scr = 0;
	VTYPE	maxh = NEVSEL;
	SKL	wskl;

	for (int m = a->left; m++ < a->right; ++as, ++bs) {
	    scr += pwd->sim2(as, bs);
	    if (LocalL && scr < 0) {
		scr = 0;
		mL = m;
	    }
	    if (LocalR && scr > maxh) {
		maxh = scr;
		mR = m;
	    }
	}
	int r = b->left - a->left;
	if (dlt < 0) r -= dlt;
	wskl.m = mL;
	wskl.n = mL + r;
	mfd->write((UPTR) &wskl);
	wskl.m = mR;
	wskl.n = mR + r;
	mfd->write((UPTR) &wskl);
	if (dlt < 0) std::swap(a, b);
	return (LocalR? maxh: scr);
}

VTYPE Aln2s1::trcbkalignS_ng(const WINDOW& wdw, bool spj, const RANGE* mc)
{
	int	ptr = 0;
	VTYPE	scr = 0;
	vmf = new Vmf();

#if !__SSE4_1__	// scalar version of forward DP 
	scr = forwardS_ng(wdw, &ptr, mc);
#else
const	int	nelem = 8;
const	int	m = a->right - a->left;
	if (m < nelem || mc) scr = forwardS_ng(wdw, &ptr, mc);
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
		trbfwd(seqs, pwd, wdw, spjcs, cip, mode, vmf);
	    scr = trbfwd.forwardS1(&ptr);
	}
#endif	// !__SSE4_1__
	if (ptr) {
	    SKL* lskl = vmf->traceback(ptr);
	    if (!lskl) {
		scr = NEVSEL;
		goto eop;
	    }
	    SKL* lwsk = lskl;
	    while (lskl->n--) mfd->write((UPTR) ++lwsk);
	    delete[] lskl;
	}
eop:
	delete vmf; vmf = 0;
	return (scr);
}

VTYPE Aln2s1::lspS_ng(const WINDOW& wdw)  // recursive
{
const	int	m = a->right - a->left;
const	int	n = b->right - b->left;
	if (!m && !n) return(0);
const	INT	aexgl = a->inex.exgl;	// reserve
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
	if (wdw.up == wdw.lw) return(diagonalS_ng());
const	float	k = wdw.lw - b->left + a->right;
const	float	q = b->right - a->left - wdw.up;
	float	cvol =  m * n - (k * k + q * q) / 2;
	if (cvol < MaxVmfSpace || m == 1 || n <= 1)
	    return (trcbkalignS_ng(wdw, b->inex.intr));

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
	scr = hirschbergS_ng(cpos, wdw, wdwf, wdwb, exgl);
#else
	if (m < 4) scr = hirschbergS_ng(cpos, wdw, wdwf, wdwb, exgl);
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
		sb1(seqs, pwd, wdw, spjcs, cip, mode);
	    scr = sb1.hirschbergS1(cpos, upf, upb, exgl);
	}
#endif	// !__SSE4_1__
	SKL	wskl = {cpos[0], 0};
	a->inex.exgr = b->inex.exgr = 0;
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
	    lspS_ng(wdwf);		// first half
	    a->left = cpos[0];
	    b->left = cpos[1];
	    a->right = aright;		// recover
	    b->right = bright;
	    a->inex.exgl = 0;
	    b->inex.exgl = exgl;
	    a->inex.exgr = aexgr;
	    b->inex.exgr = bexgr;
	}
	if (upb < INT_MAX || wdwb.lw == INT_MAX) {	// local
	    stripe(seqs, &wdwb, alprm.sh);
	    if (upb < INT_MAX) wdwb.up = upb;
	}
	wdwb.width = wdwb.up - wdwb.lw + 3;
	lspS_ng(wdwb);			// second half
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexgl;
	b->inex.exgl = bexgl;
	return scr;
}

VTYPE Aln2s1::shortcutS_ng(int ovr, const BOUND& bab)
{
	int	margin = IntronPrm.minl;
	VTYPE	scr = 0;
	ovr = (ovr > 0? 0: ovr) - 3;
	scr -= creepback(ovr, slmt, bab);
	scr -= creepfwrd(ovr, slmt, bab);
const	int	alen = a->right - a->left;
	int	sh = alen / 2;
	if (alprm.sh < 0) {
	    float f = -alprm.sh;
	    if (f > 1.) f /= 100;
	    if (f < 0.5) sh = int(alen * f);
	} else if (alprm.sh < sh)
	    sh = int(alprm.sh);
	int	minsh = alen - margin;
	if (sh < minsh) sh = minsh;
	WINDOW	wdw;
	stripe(seqs, &wdw, sh);

const	int	interval = b->right - b->left - 2 * margin;
	RANGE	rng = {b->left + margin, b->right - margin};

const	int	aexg = a->inex.exgl;
const	int	bexg = b->inex.exgl;
	a->inex.exgl = b->inex.exgl = 0;	// global
	scr += trcbkalignS_ng(wdw, true, (interval > 0)? &rng: 0);
	a->inex.exgr = aexg;
	b->inex.exgr = bexg;
	return (scr);
}

VTYPE Aln2s1::openendS_ng(int cmode)
{
	vmf = new Vmf();
	int	ptr;
	WINDOW	wdw;
	if (cmode == 3) {
	    int	ar = a->len - a->right;
	    if (a->left > ar) {
		cmode = 2;
		a->right = a->len;
	    } else {
		cmode = 1;
		a->left = 0;
		mfd->reset(0);
	    }
	}
	stripe(seqs, &wdw, alprm.sh);
	VTYPE scr = (cmode == 1)? 
	    back2ward5endS_ng(&ptr, wdw):
	    for2ward3endS_ng(&ptr, wdw);

	SKLP    sv = {0, 0, ptr};
	while (sv.p) {
	    vmf->readvmf(&sv, sv.p);
	    mfd->write((UPTR) &sv);
	}
	delete vmf;
	vmf = 0;
	return (scr);
}

VTYPE Aln2s1::backforth(int ovr, const BOUND& lub)
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
	while (i++ < ovr && m++ < lub.ua && n++ < lub.ub) {
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

/**********************************************************
   find spj within indel-free overlapped region
   doubly conted alignment score is cancelled
   assume agap <= 0
   return whether a cannonical or reverse cannonical pair 
   is found or not
**********************************************************/

bool Aln2s1::indelfreespjS(int agap, VTYPE& iscr, const bool write_skl)
{
const	int	ilen = b->right - b->left - agap;	// intron length
	if (ilen < IntronPrm.minl) {
	    iscr = pwd->GapPenalty(ilen);
	    return (true);
	}
	agap = 1 - agap;
	int	n = min3(a->left, b->left, agap + expected_overlap_ext);
	VTYPE*	bw = (n + 1 > expected_max_overlap)? new VTYPE[n + 1]: backward;
	int	i = 0;
const	CHAR*	as = a->at(a->left);	// cancel score of
const	CHAR*	bs = b->at(b->left);	// donor side overlap
const	CHAR*	ds = bs + ilen;		// bl + ilen
	bw[0] = 0;
	for (VTYPE v = 0; i < n && (*--bs == *--ds || i < agap); )
	    bw[++i] = v += pwd->sim2(--as, bs);
	vreverse(bw, i + 1);
	SKL	skl = {a->left - i, b->left - i};
	iscr = NEVSEL;

	for (int retry = 0; retry < 2 && iscr == NEVSEL; ++retry) {
	  int	m = skl.m;
	  int	n = skl.n;
	  as = a->at(m);
	  PfqItr	api(a, m);
const	  bool	usespb = api.size() && use_spb();
	  bs = ds;
	  for (VTYPE v = i = 0; n <= b->left; ++n, ++i, ++m) {
const	    bool	a_in_zone = usespb && (api == m);
const	    int	rc = b->exin->isCanon(n, n + ilen);
	    if (retry || rc) {		// canonical or revserse canonical
		VTYPE	x = b->exin->sig53(n, n + ilen, IE5P3);
		if (a_in_zone && rc > 3) x += api.match_score(m);
		VTYPE	y = x - v - bw[i];
		if (y > iscr) {
		    skl.n = n; 
		    skl.m = m;
		    iscr = y;
		}
	    }
	    v += pwd->sim2(as++, ++bs);	// cancel score of acc side
	    if (a_in_zone) ++api;
	  }
	}
	if (bw != backward) delete[] bw;
	if (iscr <= NEVSEL) return (false);
	if (write_skl) {
	    mfd->write(&skl);
	    EXIN*	bd = b->exin->score(skl.n);
	    bd->phs5 = 0;
	    skl.n += ilen;
	    mfd->write(&skl);
	    bd = b->exin->score(skl.n);
	    bd->phs3 = 0;
	    iscr += pwd->IntPen->Penalty(ilen);
	}
	return (true);
}

VTYPE Aln2s1::creepback(int ovr, VTYPE bscr, const BOUND& lub)
{
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(b->left);
	VTYPE	dscr = 0;
	while (a->left > lub.la && b->left > lub.lb
		&& (ovr < 0 || vabs(dscr) < bscr)) {
	    dscr += pwd->sim2(--as, --bs);
	    --a->left; --b->left;
	    if (++ovr == 0) bscr += dscr;
	}
	return (dscr);
}

VTYPE Aln2s1::creepfwrd(int& ovr, VTYPE bscr, const BOUND& lub)
{
const	CHAR*	as = a->at(a->right);
const	CHAR*	bs = b->at(b->right);
	VTYPE	dscr = 0;
	while (a->right < lub.ua && b->right < lub.ub
		&& (ovr < 0 || vabs(dscr) < bscr)) {
	    dscr += pwd->sim2(as++, bs++);
	    ++a->right; ++b->right;
	    if (++ovr == 0) bscr += dscr;
	}
	return (dscr);
}

//	heuristic search for short exon

int Aln2s1::nearest5ss(const BOUND& bab)
{
const	int&	from = b->left;
	int	retry = 0;
ss5:
	int	nu = from;
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(nu);
const	CHAR*	ts = a->at(std::max(bab.la, a->left - max_dist2ss));
const	EXIN*	bbu = b->exin->score(nu);
	bool	eij = false;
	for ( ; as > ts && nu > bab.lb; --nu, --bbu) {
	    eij = bbu->sig5 > b->exin->gc_sig5 || (retry && bbu->phs5 == 0);
	    if (eij || (!algmode.crs && *--as != *--bs)) break;
	}
	if (nu == from && eij) return (nu);
	int	nd = from;
	as = a->at(a->left);
	bs = b->at(nd);
	ts = a->at(std::min(a->left + max_dist2ss, bab.ua));
	eij = false;
const	EXIN*	bbd = b->exin->score(nd);
	while (as < ts && ++nd < bab.ub) {
	    eij = (++bbd)->sig5 > b->exin->gc_sig5 || (retry && bbd->phs5 == 0);
	    if (eij || (!algmode.crs && *as++ != *bs++)) break;
	}
	if (retry++ == 0 && bbu->sig5 <= 0 && bbd->sig5 <= 0) goto ss5;
	else if (bbu->phs5 && bbd->phs5) return (-1);
	if (bbu->phs5) return (nd);
	if (bbd->phs5) return (nu);
	if ((from - nu) == (nd - from)) return ((bbu->sig5 > bbd->sig5)? nu: nd);
	return ((from - nu < nd - from)? nu: nd);
}

int Aln2s1::nearest3ss(const BOUND& bab)
{
const	int&	from = b->right;
	int	retry = 0;
ss3:
	int	nu = from;
const	CHAR*	as = a->at(a->right);
const	CHAR*	bs = b->at(nu);
const	CHAR*	ts = a->at(std::max(bab.la, a->left - max_dist2ss));
const	EXIN*	bbu = b->exin->score(nu);
	bool	eij = false;
	for ( ; as > ts && nu > bab.lb; --nu, --bbu) {
	    eij = bbu->sig3 > 0 || (retry && bbu->phs3 == 0);
	    if (eij || (!algmode.crs && *--as != *--bs)) break;
	}
	if (nu == from && eij) return (nu);
	int	nd = from;
	as = a->at(a->right);
	bs = b->at(nd);
	ts = a->at(std::min(a->left + max_dist2ss, bab.ua));
const	EXIN*	bbd = b->exin->score(nd);
	while (as < ts && ++nd < bab.ub) {
	    eij = (++bbd)->sig3 > 0 || (retry && bbd->phs3 == 0);
	    if (eij || (!algmode.crs && *as++ != *bs++)) break;
	}
	if (retry++ == 0 && bbu->sig3 <= 0 && bbd->sig3 <= 0) goto ss3;
	else if (bbu->phs3 && bbd->phs3) return (-1);
	if (bbu->phs3) return (nd);
	if (bbd->phs3) return (nu);
	if ((from - nu) == (nd - from)) return ((bbu->sig3 > bbd->sig3)? nu: nd);
	return ((from - nu < nd - from)? nu: nd);
}

// no gap with mismatches

VTYPE Aln2s1::micro_exon(const BOUND& bab)
{
const	int	l = nearest5ss(bab);
	if (l < 0) return (NEVSEL);
const	int	r = nearest3ss(bab);
	if (r < 0) return (NEVSEL);
	int	d5 = l - b->left;
	int	d3 = r - b->right;
	b->left = l;
	a->left += d5;
	b->right = r;
	a->right += d3;
const	int	alen = a->right - a->left;
const	int	blen = r - l;
	VTYPE	maxscr = pwd->IntPen->Penalty(blen);
	int	f = -1;
	if (alen <= 0) {
	    SKL	skl = {a->left, b->left};
	    if (alen < 0) {
		skl.m = a->left + alen;
		mfd->write(&skl);
		skl.m = a->left;
	    }
	    mfd->write(&skl);
	    skl.n = b->right;
	    mfd->write(&skl);
	    return (maxscr + b->exin->sig53(l, r, IE5P3));
	}
	int	n5 = b->left + IntronPrm.minl;
const	int	n9 = b->right - alen - IntronPrm.minl;
const	EXIN*	bb5 = b->exin->score(n5);	// 5' end of internal exon
const	EXIN*	bb3 = b->exin->score(n5 + alen);	// 3' end of internal exon
const	CHAR*	ss = a->at(a->left);
const	CHAR*	ts = a->at(a->right);
	for ( ; n5 < n9; ++n5, ++bb5, ++bb3) {
	    if (bb5->phs3 || bb3->phs5) continue;
const	    int	n3 = n5 + alen;
const	    CHAR*	as = ss;
const	    CHAR*	bs = b->at(n5);
	    VTYPE	scr = 0;
	    while (as < ts)
		scr += pwd->simmtx->mtx[*as++][*bs++];
	    scr = alprm2.w * scr + 
		b->exin->sig53(l, n5, IE5P3) + 
		b->exin->sig53(n3, r, IE5P3) +
		pwd->IntPen->Penalty(n5 - l) +
		pwd->IntPen->Penalty(r - n3);
	    if (scr > maxscr) {
		maxscr = scr;
		f = n5;
	    }
	}
	if (f < 0) {
	    b->left -= d5;
	    a->left -= d5;
	    b->right -= d3;
	    a->right -= d3;
	    return (NEVSEL);
	}
	SKL	skl = {a->left, b->left};
	mfd->write(&skl);
	if (f != b->left) {
	    skl.n = f;
	    mfd->write(&skl);
	    skl.n = f + alen;
	}
	skl.m += alen;
	mfd->write(&skl);
	skl.n = b->right;
	mfd->write(&skl);
	return (maxscr);
}

// no gap with mismatches

int Aln2s1::first_exon_wmm(VTYPE& maxscr)
{
const	int	r = b->right;
const	int	rr = b->right - a->right;
	int	n = rr - IntronPrm.minl;
	int	nd = n + a->right;
const	EXIN*	bbi = b->exin->score(n);
const	EXIN*	bb5 = b->exin->score(nd);
const	CHAR*	ts = a->at(a->right);
const	CHAR*	ss = a->at(a->left);
	VTYPE	pmch = 0;
	for (const CHAR* as = ss; as < ts; ++as)
	    pmch += pwd->simmtx->mtx[*as][*as];
	pmch *= alprm2.w;
	int	f = -1;
	for ( ; n >= b->left; --n, --nd, --bbi, --bb5) {
	    if (!b->exin->isCanon(nd, r)) continue;
const	    VTYPE	ip = pwd->IntPen->Penalty(r - nd);
	    VTYPE	scr = bb5->sig5 + ip;
//	    if (ip + pmch < maxscr) break;	// longer intron is unlikely
const	    CHAR*	bs = b->at(n);
	    VTYPE	mchscr = 0;
	    for (const CHAR* as = ss; as < ts; )
		mchscr += pwd->simmtx->mtx[*as++][*bs++];
	    scr += alprm2.w * mchscr;
	    if (scr > maxscr) {
		f = n;
		if (mchscr == pmch) break;
		maxscr = scr;
	    }
	}
	return (f);
}

// exact matches

VTYPE Aln2s1::first_exon(const BOUND& bab)
{
	RANGE	rng[2];
	save_range(seqs, rng, 2);
const	int	r = nearest3ss(bab);
	if (r < 0) return failed(rng);		// no acceptor site nearby
const	int	d = r - b->right;
	b->right = r;
	a->right += d;
	if (a->left >= a->right || b->left >= b->right) return failed(rng);
	if (a->right == 0 || b->right == 0) {
	    SKL	skl = {a->right, b->right};
	    mfd->write(&skl);
	    return (0);
	}
	BoyerMoore	bm(b, a, -1);
	VTYPE	maxscr = NEVSEL;
	int	maxf = -1;
	while (!bm.finished()) {
const	    int	f = bm.nexthit();
	    if (f >= 0) {
		int	nd = f + a->right;
		if (b->exin->isCanon(nd, r)) {
const		    VTYPE	scr = pwd->IntPen->Penalty(r - nd) +
			b->exin->sig53(nd, r, IE5P3);
		    if (scr > maxscr) {
			maxscr = scr;
			maxf = f;
		    }
		}
	    }
	}

	if (maxf < 0) {
	    if ((maxf = first_exon_wmm(maxscr)) < 0) return failed(rng);
	} else {
const	    CHAR*	ts = a->at(a->right);
	    for (const CHAR* as = a->at(a->left); as < ts; ++as)
		maxscr += pwd->simmtx->mtx[*as][*as];
	}
	b->left = maxf;
	SKL	skl = {a->left, b->left};
	mfd->write(&skl);
	skl.m = a->right;
	skl.n = b->left + a->right;
	mfd->write(&skl);
	skl.n = b->right;
	mfd->write(&skl);
	return (maxscr);
}

// no gap with mismatches

int Aln2s1::last_exon_wmm(VTYPE& maxscr)
{
const	int	l = b->left;
const	int	alen = a->right - a->left;
const	int	rr = b->right - alen;
	int	n = b->left + IntronPrm.minl;
const	EXIN*	bb3 = b->exin->score(n);
const	CHAR*	ts = a->at(a->right);
const	CHAR*	ss = a->at(a->left);
	VTYPE	pmch = 0;			// perfect match score
	for (const CHAR* as = ss; as < ts; ++as)
	    pmch += pwd->simmtx->mtx[*as][*as];
	pmch *= alprm2.w;
	int	f = -1;
	for ( ; n < rr; ++n, ++bb3) {
	    if (!b->exin->isCanon(l, n)) continue;
const	    VTYPE	ip = pwd->IntPen->Penalty(n - l);
	    VTYPE	scr = bb3->sig3 + ip;
//	    if (ip + pmch < maxscr) break;	// longer intron is unlikely
	    const CHAR*	bs = b->at(n);
	    VTYPE	mchscr = 0;
	    for (const CHAR* as = ss; as < ts; )
		mchscr += pwd->simmtx->mtx[*as++][*bs++];
	    scr += alprm2.w * mchscr;
	    if (scr > maxscr) {
		f = n;
		if (mchscr == pmch) break;
		maxscr = scr;
	    }
	}
	return (f);
}

// exact matches

VTYPE Aln2s1::last_exon(const BOUND& bab)
{
	RANGE	rng[2];
	save_range(seqs, rng, 2);
const	int	l = nearest5ss(bab);
	if (l < 0) return failed(rng);		// no splice donor site nearby
const	int	d = l - b->left;
	b->left = l;
	a->left += d;
	if (a->left >= a->right || b->left >= b->right) return failed(rng);
	VTYPE	maxscr = NEVSEL;
	int	maxf = -1;
	BoyerMoore	bm(b, a, 1);
	while (!bm.finished()) {
const	    int	f = bm.nexthit();
	    if (f >= 0) {
		if (b->exin->isCanon(l, f)) {
		    VTYPE	scr = pwd->IntPen->Penalty(f - l) +
			b->exin->sig53(l, f, IE5P3);
		    if (scr > maxscr) {
			maxscr = scr;
			maxf = f;
		    }
		}
	    }
	}
	if (maxf < 0) {
	    if ((maxf = last_exon_wmm(maxscr)) < 0) return failed(rng);
	} else {
const	    CHAR*	ts = a->at(a->right);
	    for (const CHAR* as = a->at(a->left); as < ts; ++as)
		maxscr += pwd->simmtx->mtx[*as][*as];
	}
	SKL	skl = {a->left, b->left};
	mfd->write(&skl);
	skl.n = maxf;
	mfd->write(&skl);
	skl.m = a->right;
	skl.n += (a->right - a->left);
	mfd->write(&skl);
	return (maxscr);
}

VTYPE Aln2s1::interpolateS(INT level, const int cmode, 
	const JUXT* wjxt, const BOUND& bab)
{
	if (is3end) return (0); 
	int	agap = a->right - a->left;
	int	bgap = b->right - b->left;
	int	ovr = std::min(agap, bgap);
const	int	dgap = bgap - agap;
const	bool	cont = agap <= 0;			// connect adjacent HSPs
const	int	wlmt = setwlprm(++level)->width;
const	bool	no_rec = ovr < wlmt;			// no recurrsion
	VTYPE	iscore = NEVSEL;
	VTYPE	scr = 0;
	Mfile*	save_mfd = 0;

	if (dgap == 0) {
	    if (agap == 0) {
		if (wjxt) mfd->write((UPTR) wjxt);
		return (0);
	    }
	    iscore = diagonalS_ng();			// no gap
	} else if (cmode == 1 && no_rec && wjxt) {
	    if (cont && wjxt) {				// 5' end
		mfd->write((UPTR) wjxt);
		iscore = 0;
	    } else if (agap < IntronPrm.elmt) {
		SKL	skl = {wjxt->jx - agap, wjxt->jy - agap};
		if (skl.m >= a->left && skl.n >= b->left)
		    mfd->write(&skl);
		mfd->write((UPTR) wjxt);
		iscore = agap * getsmn(4);
	    } else if (bgap <= 0) {
		SKL	skl = {a->left = a->right, b->left = b->right};
		mfd->write((UPTR) &skl);
		iscore = 0;
	    } else
		iscore = first_exon(bab);
	} else if (cmode == 2 && no_rec) {
	    if  (agap < IntronPrm.elmt) {
		SKL	skl = {a->left, b->left};	// 3' end
		mfd->write((UPTR) &skl);
		if (agap < 0) agap = 0;
		iscore = agap * getsmn(4);
		if (agap) {
		    skl.m += agap;
		    skl.n += agap;
		    if (skl.m < a->right && skl.n < b->right)
			mfd->write(&skl);
		}
	    } else if (bgap <= 0) {
		SKL	skl = {a->left = b->right, b->left = b->right};
		mfd->write((UPTR) &skl);
		iscore = 0;
	    } else
		iscore = last_exon(bab);
	} else if (cmode == 3 && cont && dgap >= IntronPrm.minl && 
	    indelfreespjS(agap, iscore)) {		// indel-free spj
	    scr += iscore;
	    iscore = 0;
	} else if (cmode == 3 && no_rec && dgap >= IntronPrm.minl) {
	    if (algmode.crs == 0) iscore = micro_exon(bab);
	    if (iscore == NEVSEL && agap < IntronPrm.elmt)
		iscore = shortcutS_ng(ovr, bab);	// shortcut
	} else if (ovr <= 0 && dgap < IntronPrm.minl) {
	    iscore = backforth(-ovr, bab);		// ordinary gap
	} else if (abs(dgap) < IntronPrm.minl) {
	    scr -= creepback(ovr, slmt, bab);
	    scr -= creepfwrd(ovr, slmt, bab);
	    WINDOW	wdw;
	    stripe(seqs, &wdw, std::min(alprm.sh, abs(dgap) + 3));
	    iscore = trcbkalignS_ng(wdw, false);	// ordinary alignment
	} else if (level < algmode.qck) {		// recursive search
	    save_mfd = new Mfile(*mfd);
	    iscore = seededS_ng(level, cmode, bab);
	}
	if (iscore == NEVSEL && (no_rec || level == algmode.qck)) {
	    RANGE	rng[2];
	    save_range(seqs, rng, 2);
	    if (cmode & 1) scr -= creepfwrd(ovr, slmt, bab);
	    if (cmode & 2) scr -= creepback(ovr, slmt, bab);
	    agap += rng[0].left - a->left + a->right - rng[0].right;
	    bgap += rng[1].left - b->left + b->right - rng[1].right;
	    bool	try_dp = algmode.crs || agap < a->len * gudp;
	    if (try_dp) {				// try DP
		if (save_mfd) *mfd = *save_mfd;
		else	save_mfd = new Mfile(*mfd);
const		float	dpspace = (float) agap * (float) bgap / MEGA;
		try_dp = algmode.crs || fabs(dpspace) < alprm.maxsp;
		if (!try_dp && cmode < 3) {
		    bgap = IntronPrm.maxl;
		    if (cmode == 1) b->left = std::max(b->right - bgap, 0);
		    else	b->right = std::min(b->left + bgap, b->len);
		    if (b->exin) b->exin->resize();
		    try_dp = true;
		}
		if (try_dp) {
		    WINDOW	wdw;
		    stripe(seqs, &wdw, alprm.sh);
		    iscore = lspS_ng(wdw);		// DP
		}
	    }
	    if (iscore == NEVSEL) {
		if (save_mfd) *mfd = *save_mfd;		// reset
		if (cmode == 1) {
		    if (wjxt) {				// near 5' end
			int	bl = wjxt->jy + end_margin;
			if (bl > b->left) b->left = bl;
		    }
		    iscore = openendS_ng(cmode);	// drop off
		} else if (cmode == 2) {
		    if (wjxt) {				// near 3' end
			int	br = bgap - wjxt->jy - end_margin;
			if (br > b->left && br < b->right) b->right = br;
		    }
		    iscore = openendS_ng(cmode);	// drop off
		} else {
		    iscore = Local?			// give up
			openendS_ng(cmode): shortcutS_ng(ovr, bab);
		}
	    }
	}
	delete save_mfd;
	return (scr + iscore);
}

WLUNIT* Aln2s1::bestwlu(WLUNIT* wlu, const int nwlu, const int cmode)
{
	RANGE	rng[2];
	save_range(seqs, rng, 2);
	VTYPE	wlu_s = NEVSEL;
	int	wlu_n = 0;
	for (int n = 0; n < nwlu; ++n) {
const	    JUXT*	jxt = wlu[n].jxt;
	    a->left = rng[0].left;	// restore
	    b->left = rng[1].left;
	    a->right = jxt->jx;
	    b->right = jxt->jy;
	    int	agap = jxt->jx - a->left;
	    if (agap > 0) continue;
	    VTYPE	iscore = NEVSEL;
	    VTYPE	jscore = 0;
	    if (cmode == 1)
		jscore = wlu[n].scr;
	    else if (indelfreespjS(agap, iscore, false))
		jscore = wlu[n].scr + iscore;
	    else continue;
	    jxt = wlu[n].jxt + wlu[n].num - 1;
	    a->left = jxt->jx + jxt->jlen;
	    b->left = jxt->jy + jxt->jlen;
	    a->right = jxt[1].jx;
	    b->right = jxt[1].jy;
	    agap = jxt[1].jx - a->left;
	    if (agap > 0) continue;
	    if (cmode != 2) {
		if (indelfreespjS(agap, iscore, false))
		    jscore += iscore;
		else	continue;
	    }
	    if (jscore > wlu_s) {
		wlu_s = jscore;
		wlu_n = n;
	    }
	}
	rest_range(seqs, rng, 2);
	wlu = wlu_s > NEVSEL? wlu + wlu_n: 0;
	return (wlu);
}

// eimode 1: 5'end, 2: 3' end, 3: internal
VTYPE Aln2s1::seededS_ng(INT level, int eimode, const BOUND& lub)
{
const	INT	aexgl = a->inex.exgl;
const	INT	aexgr = a->inex.exgr;
const	INT	bexgl = b->inex.exgl;
const	INT	bexgr = b->inex.exgr;
	RANGE	rng[2];
	int	cmode = eimode;
	int	num = 0;
	VTYPE	scr = 0;
	JUXT*	jxt = 0;
	JUXT*	wjxt = 0;
	Wilip*	wl = 0;
	WLUNIT*	wlu = 0;

	save_range(seqs, rng, 2);
const	int	wlmt = setwlprm(level)->width;
	BOUND	bab = lub;
	if (level == lowestlvl && b->jxt) {
	    jxt = b->jxt;
	    num = b->CdsNo;
	} else {
	    wl = new Wilip(seqs, pwd, level);
	    int	nwlu = wl->size();
	    wlu = wl->begin();
const	    int	bgap = b->right - b->left;
	    if (nwlu > 1 && bgap >= IntronPrm.minl)
		wlu = bestwlu(wlu, nwlu, cmode);
	    if (wlu) {
		num = wlu->num;
		jxt = wlu->jxt;
	    } else if (nwlu > 1)
		level = algmode.qck - 1;
	}
	if (num) {
	    jxt[num].jx = a->right;
	    jxt[num].jy = b->right;
	    a->inex.exgr = 0;
	    b->inex.exgr = 0;
	    for (wjxt = jxt; num--; wjxt++) {
		scr += wjxt->jscr;
		a->right = wjxt->jx;
		b->right = wjxt->jy;
		bab.ua = wjxt->jx + wjxt->jlen;
		if (wjxt[1].jx < bab.ua) bab.ua = wjxt[1].jx;
		bab.ua -= wlmt;
		int	lx = wjxt->jx + wjxt->jlen / 2;
		if (bab.ua < lx) bab.ua = lx;
		bab.ub = wjxt->jy + bab.ua - wjxt->jx;
		if (cmode == 2) cmode = 3;
		VTYPE	iscore = 
		    interpolateS(level, cmode, wjxt, bab);
		if (iscore != NEVSEL) {
		    scr += iscore;
		    cmode = 3;
		    a->left = wjxt->jx + wjxt->jlen;
		    b->left = wjxt->jy + wjxt->jlen;
		    a->inex.exgl = 0;
		    b->inex.exgl = 0;
		    bab.la = a->right;
		    bab.lb = b->right;
		}
	    }
	    a->inex.exgr = aexgr;
	    b->inex.exgr = bexgr;
	    a->right = rng[0].right;
	    b->right = rng[1].right;
	    bab.ua = lub.ua;
	    bab.ub = lub.ub;
	    if (eimode == 2 || (level == INT(lowestlvl) && eimode == 1))
		cmode = 2;
	}

const	VTYPE	iscore = interpolateS(level, cmode, wjxt, bab);
	if (iscore > NEVSEL) scr += iscore;
	else	scr = NEVSEL;
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexgl;
	b->inex.exgl = bexgl;
	if (level == lowestlvl && wjxt) {
	    wjxt->jx = a->len;
	    wjxt->jy = b->len;
	}
	delete wl;
	return (scr);
}

SKL* Aln2s1::globalS_ng(const WINDOW& wdw, VTYPE* scr)
{
	mfd = new Mfile(sizeof(SKL));
	SKL	wsk;
	mfd->write((UPTR) &wsk);	// dummy call
	if (algmode.qck) {
	    BOUND bab = {a->left, b->left, a->right, b->right};
	    *scr = seededS_ng(lowestlvl, 1, bab);
	} else
	    *scr = lspS_ng(wdw);
	wsk.n = (int) mfd->size();
	SKL*	skl = (SKL*) mfd->flush();
	skl->n = wsk.n - 1;
	if (skl->n < 2) {
	    delete[] skl;
	    return (0);
	}
	skl->m = AlgnTrb;
	if (a->inex.sens) skl->m |= A_RevCom;
	return (trimskl(seqs, stdskl(&skl)));
}

VTYPE HomScoreS_ng(const Seq* seqs[], const PwdB* pwd)
{
	Aln2s1 alnv(seqs, pwd);
	WINDOW	wdw;
	stripe(seqs, &wdw, alprm.sh);
#if !__SSE4_1__	// scalar version of forward DP 
	return (alnv.scorealoneS_ng(wdw));
#else
const	Seq*&	a = seqs[0];
const	int	m = a->right - a->left;
	if (m < 4) return (alnv.scorealoneS_ng(wdw));
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
	    fwds(seqs, pwd, wdw, alnv.spjcs, alnv.cip, 0);
	return (fwds.scoreonlyS1());
#endif	// !__SSE4_1__
}

static int infer_orientation(Seq** sqs, const PwdB* pwd)
{
	Seq*&	a = sqs[0];
const	VTYPE	scr1 = HomScoreS_ng((const Seq**) sqs, pwd);
	a->comrev();
	antiseq(sqs + 1);
const	VTYPE	scr2 = HomScoreS_ng((const Seq**) sqs, pwd);
	if (scr2 > scr1) return (1);
	a->comrev();
	antiseq(sqs + 1);
	return (0);
}

static void reverse_copy_jxt(Seq** seqs)
{
	Seq*&	b = seqs[0];
	if (!b->jxt) return;
	Seq*	c = b;
	if (b->getanti()) {
	    c = seqs[1];
	    if (c->jxt) delete[] c->jxt;
	    c->CdsNo = b->CdsNo;
	    c->jxt = new JUXT[b->CdsNo + 1];
	    vcopy(c->jxt, b->jxt, b->CdsNo + 1);
	}
	c->revjxt();
}

SKL* alignS_ng(Seq* seqs[], const PwdB* pwd, VTYPE* scr, int ori)
{
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
	if (ori == 2) {
	    a->comrev();
	    antiseq(seqs + 1);
	} else if (ori == 3 && !algmode.qck)
	    ori = infer_orientation(seqs, pwd);
	Aln2s1 alnv((const Seq**) seqs, pwd);
	WINDOW	wdw;
	stripe((const Seq**) seqs, &wdw, alprm.sh);
	SKL*	fskl = alnv.globalS_ng(wdw, scr);
	if (ori != 3) return (fskl);
	VTYPE	rscr = NEVSEL;
	reverse_copy_jxt(seqs + 1);
	a->comrev();
	antiseq(seqs + 1);
	alnv.reset((const Seq**) seqs);
	SKL*	rskl = alnv.globalS_ng(wdw, &rscr);
	if (*scr >= rscr) {
	    delete[] rskl;
	    a->comrev();
	    antiseq(seqs + 1);
	    if (!b->getanti()) b->revjxt();
	} else {
	    delete[] fskl;
	    fskl = rskl;
	    *scr = rscr;
	}
	return (fskl);
}

