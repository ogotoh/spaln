/*****************************************************************************
*
*	Alignment of protein vs. nucleotide sequences.
*	5' and 3' splice site signals and coding potential are considered.
*	Assumes there is no internal gap in the reference protein sequence.
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
*	Osamu Gotoh, Ph.D.      (2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#define	DEBUG	1
#define	TERMGOP	0

#include "aln.h"
#include "vmf.h"
#include "wln.h"
#include "boyer_moore.h"

static	const	int	NQUE = 3;
static	const	int	expected_max_overlap = 1024;
static	const	int	expected_overlap_ext = 16;
static	const	int	max_dist2ss = 5;
static	const	int	n_nearest_ss = 2;

inline	bool	avst_equal(const CHAR& a, const CHAR& b) 
	{return (a == b || (a == SER && b == SER2));}

// alignment parameters independent of each sequence

class Aln2h1 {
protected:
const	Seq**	seqs;
const	Seq*&	a;
const	Seq*&	b;
const	PwdB*	pwd;
	Mfile*	mfd;
	Vmf*	vmf;
	SpJunc*	spjcs;
	VTYPE	pinthr;
	WLPRM*	wlprm;
	INT	lowestlvl;
	int	term;
	VTYPE	slmt;
	bool	island;
	int	ss[n_nearest_ss];
	Cip_score*	cip;
const	int	Nod;
	VTYPE	backward[expected_max_overlap];
public:
	Aln2h1(const Seq** _seqs, const PwdB* _pwd);
	~Aln2h1() {delete mfd; delete spjcs; delete cip;}
	void	initH_ng(RVPD* hf[], const WINDOW& wdw, const RANGE* cutrng = 0);
	RVPD*	lastH_ng(RVPD* hhb[], const WINDOW& wdw, const RANGE* cutrng = 0);
	VTYPE	forwardH_ng(long pp[], const WINDOW& wdw, 
		bool spj = true, const RANGE* cutrng = 0);
	VTYPE	shortcutH_ng(int ovr, BOUND& lub);
	void	hinitH_ng(RVDWU* hhf[], const WINDOW& wdw);
	RVDWU&	hlastH_ng(RVDWU* hhf[], const WINDOW& wdw);
	VTYPE	hirschbergH_ng(int cpos[], 
		const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb);
	void	pfinitH_ng(RVPD* hh[], const WINDOW& wdw);
	void	pbinitH_ng(RVPD* hh[], const WINDOW& wdw);
	VTYPE	back2ward5endH_ng(long* ptr, const WINDOW& wdw);
	VTYPE	for2ward3endH_ng(long* ptr, const WINDOW& wdw);
	VTYPE	openendH_ng(int cmode);
	bool	indelfreespjH(int agap, VTYPE& scr);
	bool	isEIJ(int phs) {return (phs > -2);}
	VTYPE	diagonalH_ng();
	VTYPE	backforth(int n, BOUND& lub);
	VTYPE	trcbkalignH_ng(const WINDOW& wdw, 
		bool spj = true, const RANGE* cutrng = 0);
	VTYPE	lspH_ng(const WINDOW& wdw);
	int	nearest5ss(int from, BOUND& bab);
	int	nearest3ss(int from, BOUND& bab);
	int	first_exon_wmm(int d3, VTYPE& maxscr, bool& pm, int nss);
	VTYPE	first_exon(BOUND& bab);
	int	last_exon_wmm(int d3, VTYPE& maxscr, bool& pm, int nss);
	VTYPE	last_exon(BOUND& bab, SKL* pskl);
	VTYPE	micro_exon(BOUND& bab);
	VTYPE	cds5end(const SKL* pskl);
	VTYPE	cds3end(const SKL* pskl, VTYPE maxscr = 0);
	VTYPE	interpolateH(INT level, int agap, int bgap, int cmode, 
		const JUXT* wjxt, BOUND& bab, int wlmt);
	VTYPE	seededH_ng(INT level, int cmode, BOUND& lub);
	SKL*	globalH_ng(VTYPE* scr, const WINDOW& wdw);
	VTYPE	creepback(int ovr, VTYPE bscr, BOUND& lub);
	VTYPE	creepfwrd(int& ovr, VTYPE bscr, BOUND& lub);
};

static	VTYPE	SumCodePot(const EXIN* bb, int i, const CHAR* cs, const PwdB* pwd);
static	void	addsigEjxt(JUXT* jxt, int num, const Seq* b);

// initialization of alignment process

Aln2h1::Aln2h1(const Seq* _seqs[], const PwdB* _pwd) :
	seqs(_seqs), a(seqs[0]), b(seqs[1]), pwd(_pwd), 
	mfd(0), vmf(0), wlprm(0), lowestlvl(b->wllvl),
	term((int) (pwd->Vthr / pwd->BasicGEP)),
	slmt(pwd->Vthr / 2), island(false), Nod(2 * pwd->Noll - 1)
{
	pinthr = pwd->GapPenalty(1);
	spjcs = new SpJunc(b, pwd);
	cip = new Cip_score(a);
}

void Aln2h1::initH_ng(RVPD* hh[], const WINDOW& wdw, const RANGE* cutrng)
{
	int	n = b->left;
	int	r = b->left - 3 * a->left;
	int	rr = (cutrng? cutrng->left: b->right) - 3 * a->left;
	int	dir = a->inex.exgl? DEAD: DIAG;
	int	jnc[3] = {n};
const	EXIN*	bb = b->exin->score(n);

	RVPD*	h = hh[0] + r;
	h->val = (a->inex.exgl && bb[1].sigS > 0)? bb[1].sigS: 0;
	h->dir = dir;
	h->ptr = vmf? vmf->add(a->left, n, 0): 0;
	if (a->inex.exgl) {		// semi-gloabal
	    if (wdw.up < rr) rr = wdw.up;
	    for (int i = 1; ++r <= rr; ++i) {
		++h; ++bb; ++n;
		if (i < 3) {
 		    h->val = (bb[1].sigS > 0)? bb[1].sigS: 0;
		    h->dir = dir;
		    h->ptr = vmf? vmf->add(a->left, n, 0): 0;
		    jnc[i] = n;
		} else {
		    *h = h[-3];
		    int	k = n - jnc[i % 3];
		    if (k == 3 && !(a->inex.exgl & 1)) h->val += pwd->BasicGOP;
		    if (!(a->inex.exgl & 2)) h->val += pwd->GapExtPen3(k);
		    h->val += bb[-2].sigE;
		    h->dir = HORI;
		    VTYPE	x = h[-1].val + pwd->GapW1;
		    if (x > h->val) {*h = h[-1]; h->val = x; h->dir = HOR1;}
		    x = h[-2].val + pwd->GapW2;
		    if (x > h->val) {*h = h[-2]; h->val = x; h->dir = HOR2;}
		}
		VTYPE	x = (bb[1].sigS > 0)? bb[1].sigS: 0;
	 	if (h->val < x) {
		    h->val = x;
		    h->dir = DEAD;
		    h->ptr = vmf? vmf->add(a->left, n, 0): 0; 
		    jnc[i % 3] = n;
		}
	    }
	}

	r = b->left - 3 * a->left;
	rr = b->left - 3 * a->right;
	h = hh[0] + r;
	RVPD*	g = hh[1] + r;
	if (wdw.lw > rr) rr = wdw.lw;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --g;
	    if (b->inex.exgl == 3) {	// no terminal gap penalty
		h->val = 0;
		h->dir = DEAD;
		h->ptr = vmf? vmf->add(a->left + i / 3, b->left + i % 3, 0): 0;
	    } else if (i <= 3) {
		*h = h[i];
		if (!(b->inex.exgl & 2)) h->val += pwd->BasicGEP;
		if (!(b->inex.exgl & 1)) h->val += pwd->BasicGOP;
		if (i < 3) h->val += pwd->ExtraGOP;
		h->dir = VERT;
		*g = *h;
	    } else {
		*h = h[3];
		*g = g[3];
		if (!(b->inex.exgl & 2)) {
		    h->val += pwd->GapExtPen3(i);
		    g->val += pwd->BasicGEP;
		}
	    }
	}
}

RVPD* Aln2h1::lastH_ng(RVPD* hh[], const WINDOW& wdw, const RANGE* cutrng)
{
	int	cutlen = cutrng? cutrng->right - cutrng->left: 0;
	int	glen[3];
	int	rw = wdw.lw;
	int	m3 = 3 * a->right;
	int	rf = (cutrng? cutrng->left: b->left) - m3;
	if (rf > rw) rw = rf;
	else	rf = rw;
	RVPD*	h = hh[0] + rw - cutlen;
	RVPD*	h9 = hh[0] + b->right - m3 - cutlen;
	RVPD*	f = h9;
	VTYPE	mx = h9->val;
const	EXIN*	bb = b->exin->score(rw + m3);

	for (int p = 0; p < 3; ++p) glen[p] = 0;
	if (a->inex.exgr) {
	    for (int p = 0; h <= h9; ++h, ++bb, ++rf, ++p) {
		if (p == 3) p = 0;
		VTYPE	y = NEVSEL;
		glen[p] += 3;
		if (rf - rw >= 3 && h[-3].dir != DEAD) {
		    VTYPE	x = h[-3].val + bb[-2].sigE;
		    if (!(a->inex.exgr & 2)) x += pwd->GapExtPen3(glen[p]);
		    if (!(a->inex.exgr & 1) && glen[p] == 3) x += pwd->BasicGOP;
		    if ((algmode.lcl & 2) && bb[-2].sigT > 0 && !(h->dir & SPIN))
			y = h[-3].val + bb[-2].sigT;
		    if (x > h->val) {
			h->val = x;
			h->dir = HORI;
			h->ptr = h[-3].ptr;
		    } else if (isnthori(h))
			glen[p] = 0;
		}
		VTYPE	x = h->val;
		if ((algmode.lcl & 16) && bb->sig5 > 0) x += bb->sig5;
		if (x > mx && x >= y) {
		    f = h;
		    mx = x;
		} else if (y > mx) {	// termination codon
		    f = h;
		    *h = h[-3];
		    f->ptr = vmf->add(a->right, rf + m3 - 3, h->ptr);
		    mx = y;
		    h->dir = DEAD;
		}
	    }
	}
	f->val = mx;
	if (b->inex.exgr) {
	    rw = std::min(wdw.up, b->right - 3 * a->left) - cutlen;
	    VTYPE	g[3] = {NEVSEL, NEVSEL, NEVSEL};
	    h = hh[0] + rw - 3;
	    for (int p = 0; h >= h9; --h) {
		VTYPE	x = h[3].val;
		if (!(b->inex.exgr & 1)) x += pwd->BasicGOP;
		if (x > g[p]) g[p] = x;
		if (!(b->inex.exgr & 2)) g[p] += pwd->BasicGEP;
		if (h->val > g[p]) g[p] = NEVSEL;
		else if (g[p] > f->val) {f = h; f->val = g[p];}
		if (++p == 3) p = 0;
	    }
	}
	int p = f - h9;
	rf = a->right;		// m9
	rw = b->right;		// n9
	if (p > 0) {
	    rf -= (p + 2) / 3;
	    if (p %= 3) rw -= (3 - p); 
	} else if (p < 0) rw += p;
	if (vmf) f->ptr = vmf->add(rf, rw, f->ptr);
	return (f);
}

VTYPE Aln2h1::forwardH_ng(long* pp, const WINDOW& wdw, 
	bool spj, const RANGE* cutrng)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	RVPD*	hh[NOL + 1];
	RVPD*	hf[NOD];		// [DIAG, HORI, VERT, HORL, VERTL]
	RVPD	e1[2 * NQUE];
	RVPD*	e2 = e1 + NQUE;
	RVPDJ	hl[3][NCAND + 1];	// [phase][candidates]
	int	nx[3][NCAND + 1];	// [phase][candidates]
	RVPDJ*	maxphl[NOD];
	VSKLP   maxh = {NEVSEL, a->left, b->left, 0L};
	int     Local = algmode.lcl & 16;
	int     LocalL = Local && a->inex.exgl && b->inex.exgl;
	int     LocalR = Local && a->inex.exgr && b->inex.exgr;
	int	cutlen = cutrng? (cutrng->right - cutrng->left): 0;
	VTYPE	longgep = pwd->BasicGEP * cutlen / 3;
	VTYPE	longgep2 = pwd->LongGEP * cutlen / 3;

	int	width = wdw.width - cutlen;
	size_t	bufsiz = pwd->Noll * width;
	RVPD*	buf = new RVPD[bufsiz];
	vset(buf, black_vpd, bufsiz);
	hh[0] = buf - wdw.lw + 3;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + width;
	if (vmf) vmf->add(0, 0, 0);	// Skip 0-th record
	initH_ng(hh, wdw);

	int	m = a->left;
	if (!a->inex.exgl) --m;		// global
	int	n1 = 3 * m + wdw.lw - 1;
	int	n2 = 3 * m + wdw.up;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as) {
	    bool	internal = spj && (!a->inex.exgr || m < a->right);
	    n1 += 3; n2 += 3;
	    int		n0 = std::max(n1, b->left);
	    int		n9 = std::min(n2, b->right);
	    int		n  = n0;
	    int		r = n - 3 * m;
	    vset(e1, black_vpd, 2 * NQUE);
	    if (!b->inex.exgl && m == a->left) {
		e1[2] = e2[2] = hh[0][r];
		e1[2].val = pwd->GapW3;
		e2[2].val = pwd->GapW3L;
	    }
const	    CHAR*	bs = b->at(n - 1);
const	    EXIN*	bb = b->exin->score(n);
	    int		k = 0;
	    RVPD*&	h = hf[0] = hh[k++] + r;
	    RVPD*&	g = hf[2] = hh[k++] + r;
	    RVPD*&	g2 = hf[4] = dagp? hh[k++] + r: 0;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    vset((RVPDJ*) &hl, black_vpdj, 3 * (NCAND + 1));
	    for (int p = 0; p < 3; ++p)
		for (int l = 0; l <= NCAND; ++l)
		    nx[p][l] = l;
	    int	ncand[3] = {-1, -1, -1};
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d %2d", a->SiteNz(m), b->SiteNz(n), h->dir);
		putvar(h->val); putchar('\n');
	    }
#endif
	    for (int q = 0; ++n <= n9; ++bs) {
		VTYPE	x, y;
		++bb; ++h; ++g; if (g2) ++g2;
		RVPD*&	eq1 = hf[1] = e1 + q;
		RVPD*&	eq2 = hf[3] = dagp? e2 + q: 0;
	        RVPD	hq = *h;	// previous state

//	diagonal match
		VTYPE	sigE = bb[-2].sigE;
		RVPD*	from = h;
		RVPD*	mx = h;
		if (m == a->left) goto Horizon;
		if (n > b->left + 2) {
		    h->val += qprof[*bs] + sigE;
		    h->dir = isdiag(from)? DIAG: NEWD;
		} else	*h = black_vpd;

//	vertical gap extention
		y = g[3].val + pwd->BasicGEP;

//	1 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    g->val = x;
		    g->dir = SLA2;
		    g->ptr = from->ptr;
		} else	g->val = y;

//	2 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val) {
		    g->val = x;
		    g->dir = SLA1;
		    g->ptr = from->ptr;
		}

//	normal deletion
		x = (++from)->val + pwd->GapW3;
		if (x >= g->val) {
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		    g->ptr = from->ptr;
		} else if (y >= g->val) {
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		    g->ptr = g[3].ptr;
		}
		if (g->val > mx->val) mx = g;

//	long deletion
		if (dagp) {
		  x = from->val + pwd->GapW3L;
		  y = g2[3].val + pwd->LongGEP;
		  if (x >= y) {
		    g2->val = x;
		    g2->dir = VERL | (g2->dir & SPJC);
		    g2->ptr = from->ptr;
		  } else {
		    *g2 = g2[3];
		    g2->val = y;
		  }
		  if (g2->val > mx->val) mx = g2;
		}
Horizon:
//	nomal insertion
		if (n > n0 + 2) {
		    from = h - 3;
		    x = from->val + pwd->GapW3;
		    y = eq1->val += pwd->BasicGEP;
		    if (x > y) {
			*eq1 = *from;
			eq1->val = x;
		    }
		    eq1->val += sigE;
		    eq1->dir = (eq1->dir & SPIN) + HORI;
//	long insertion
		    if (dagp) {
			x = from->val + pwd->GapW3L;
			y = eq2->val += pwd->LongGEP;
			if (x > y) {
			    *eq2 = *from;
			    eq2->val = x;
			}
			eq2->val += sigE;
			eq2->dir = (eq2->dir & SPIN) + HORL;
			if (eq2->val > mx->val) mx = e2 + q;
		    }
		}
//	2 nt insertion
		if (n > n0 + 1) {
		    from = h - 2;
		    x = from->val + pwd->GapW2;
		    if (x > eq1->val) {
			*eq1 = *from;
			eq1->val = x;
			eq1->dir = (eq1->dir & SPIN) + HOR2;
		    }
		}
//	1 nt insertion
		from = h - 1;
		x = from->val + pwd->GapW1;
		if (x > eq1->val) {
		    *eq1 = *from;
		    eq1->val = x;
		    eq1->dir = (eq1->dir & SPIN) + HOR1;
		}
		if (eq1->val > mx->val) mx = e1 + q;
		if (++q == NQUE) q = 0;

//	intron 3' boundary, assume no overlapping signals
		if (internal && isEIJ(bb->phs3)) {
		    int	phs = (bb->phs3 == 2)? -1: bb->phs3;
Acceptor:
		    int nb = n - phs;
		    VTYPE sigJ = cip->cip_score(3 * m - phs);
		    int*	pnl = nx[phs + 1];
		    vclear(maxphl, Nod);
		    for (int l = 0; l <= ncand[phs + 1]; ++l) {
			RVPDJ*	phl = hl[phs + 1] + pnl[l];
			if (nb - phl->jnc < IntronPrm.minl) continue;
			x = phl->val + sigJ + spjcs->spjscr(phl->jnc, nb);
const			CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
			if (phl->dir == 0 && phs) {
			    if (phs == 1) x += qprof[*cs];
			    else	x += qprof[*++cs];
			}
			from = hf[phl->dir];
			if (x > from->val) {
			    from->val = x;
			    maxphl[phl->dir] = phl;
			}
		    }
		    for (int d = 0; d < Nod; ++d) {
			RVPDJ*	phl = maxphl[d];
			if (!phl) continue;
			from = hf[d];
		 	if (vmf) {
			    from->ptr = vmf->add(m, n, 
			    vmf->add(m, phl->jnc + phs, phl->ptr));
			}
			from->dir |= SPJCI;
			if (from->val > mx->val) mx = from;
		    }
		    if (bb->phs3 - phs == 3) {	// AGAG
			phs = 1;
			goto Acceptor;
		    }
		}

//	Find optimal path
		y = h->val;
		if (h != mx) *h = *mx;	// non-diagonal
		else if (Local && y > hq.val) {
		    if (LocalL && hq.dir == 0 && !(h->dir & SPJC))
			h->ptr = vmf? vmf->add(m - 1, n - 3, 0): 0;
		    else if (LocalR && y > maxh.val) {
			maxh.val = y;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && h->val <= 0) h->val = h->dir = 0;
		else if (vmf && h->dir == NEWD)
		    h->ptr = vmf->add(m - 1, n - 3, h->ptr);

//	intron 5' boundary
		if (internal && isEIJ(bb->phs5)) {
		    int	phs = (bb->phs5 == 2)? -1: bb->phs5;
Donor:
		    int	nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(nb, 0, IE5);
		    int		hd = dir2nod[mx->dir & 15];
		    for (int k = (hd == 0 || phs == 1)? 0: 1; k < Nod; ++k) {
			bool	crossspj = phs == 1 && k == 0;
			from = crossspj? &hq: hf[k];
			// An orphan exon is disallowed
			if (!from->dir || (from->dir & SPIN)) continue;
			if (!crossspj && k != hd && hd >= 0) {
			    y = mx->val;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (from->val <= y) continue;	// prune
			}
			x = from->val + sigJ;
			RVPDJ*	phl = hl[phs + 1];
			int*	pnl = nx[phs + 1];
			int&	nc = ncand[phs + 1];
			int	l = nc < NCAND? ++nc: NCAND;
			while (--l >= 0) {
			    if (x >= phl[pnl[l]].val)
				std::swap(pnl[l], pnl[l + 1]);
			    else
				break;
			}
			if (++l < NCAND) {
			    phl += pnl[l];
			    phl->val = x;
			    phl->jnc = nb;
			    phl->dir = k;
			    phl->ptr = (vmf && crossspj && isntdiag(from))?
				vmf->add(m - 1, n - 3, from->ptr): from->ptr;
			} else --nc;
		    }
		    if (bb->phs5 - phs == 3) {	//	GTGT..
			phs = 1;
			goto Donor;
		    }
		}

#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", a->SiteNz(m), b->SiteNz(n), mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(eq1->val);
		if (dagp) {
		    putvar(g2->val); putvar(eq2->val);
		}
		putvar(hl[2][nx[2][0]].val); 
		putvar(hl[1][nx[1][0]].val);
		putvar(hl[0][nx[0][0]].val);
		printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
		printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		putchar('\n');
	}
#endif

		if (cutrng && n == cutrng->left) {	// shortcut
		    h -= 3; g -= 3;
		    if (dagp) g2 -= 3;
		    for (int p = 0; p < 3; ++p) {
			if (++q == 3) q = 0;
			e1[q].val += longgep;
			if (dagp) e2[q].val += longgep2;
			*++h = dagp? e2[q]: e1[q];
			*++g = black_vpd;
			if (dagp) *++g2 = black_vpd;
		    }
		    n += cutlen;
		    bs += cutlen;
		    bb += cutlen;
		}
	    }	// end of n-loop
	}	// end of m-loop

	RVPD*	mx = lastH_ng(hh, wdw, cutrng);
	if (maxh.val > mx->val) {
	    if (pp) *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    maxh.val = mx->val;
	    if (pp) *pp = mx->ptr;
	}
	delete[] buf;
	return (maxh.val);
}

static VTYPE SumCodePot(const EXIN* bb, int i, const CHAR* cs, const PwdB* pwd)
{
	VTYPE	y = 0;

	for (++bb; i > 0; i -= 3, bb += 3) {
	    VTYPE	hvl;
	    if (cs) {
		hvl = 0;
		cs = 0;
	    } else
		hvl = bb->sigE;
	    y += hvl;
	}
	return (y);
}

VTYPE skl_rngH_ng(const Seq* seqs[], Gsinfo* gsi, const PwdB* pwd)
{
static	const char*	fswarn = "%s %s FrameShift: %d %d\n";
const	Seq*	a = seqs[0];
const	Seq*	b = seqs[1];
	SKL*	wsk = gsi->skl;
	int	num = (wsk++)->n;
	SpJunc*	spjcs = new SpJunc(b, pwd);
	VTYPE	h = 0;		// 
	VTYPE	hi = NEVSEL;		// intron-containing path
	VTYPE	ha = 0;
	VTYPE	hb = 0;
	VTYPE	hvl = 0;
	bool	ivl = false;
	VTYPE	gop = 0;
	VTYPE	sig5 = 0;
	VTYPE	sig3 = 0;
	int	insert = 0;
	int	deletn = 0;
	int	intlen = 0;
	int	preint = 0;
	int	phs = 0;
	int	nb = 0;
	int	maxexon = IntronPrm.rlmt;
	EISCR	rbuf;
	FSTAT*	fst = &gsi->fstat;
	FSTAT	pst;
	Eijnc*	eijnc = gsi->eijnc = new Eijnc();
	Cigar*	cigar = gsi->cigar = 0;
	Vulgar*	vlgar = gsi->vlgar = 0;
	switch (algmode.nsa) {
	    case CIG_FORM: cigar = gsi->cigar = new Cigar(); break;
	    case VLG_FORM: vlgar = gsi->vlgar = new Vulgar(); break;
	    default: break;
	}
	vclear(fst);
	vclear(&pst);
	vclear(&rbuf);
	gsi->noeij = 0;
const	CHAR*	cs = b->at(wsk[num - 1].n - 2);
	int	termcodon = (*cs == TRM || *cs == TRM2);
	cs = 0;
	if (termcodon) wsk[num-1].n -= 3;
	if (wsk[1].n == wsk->n && b->inex.exgl) {++wsk; --num;}
	int	m = wsk->m;
	int	n = wsk->n;
const	CHAR*	as = a->at(m);
const	CHAR*	bs = b->at(n);
const	EXIN*	bb = b->exin->score(n);
const	EXIN*	bb_last = b->exin->score(b->right);
	PfqItr	api(a, m);
	int	api_size = api.size();
	bool	usespb = api_size && use_spb();
	if ((algmode.lcl & (16 + 1)) && bb[1].sigS > h)	h = bb[1].sigS;
	if ((algmode.lcl & (16 + 4)) && bb->sig3 > h)	h = bb->sig3;
	rbuf.left = n;
	rbuf.rleft = m;
	rbuf.iscr = NEVSEL;
	rbuf.sig3 = h;
	if (cigar && m) if (m) cigar->push('H', m);	// local alignment
	while (--num > 0 || hi > NEVSEL) {
	    if (num > 0) ++wsk;
	    bool	term = num == 1;		// tail gap?
	    int	mi = (wsk->m - m) * 3;
	    if (insert && (mi || (h > NEVSEL && hi > NEVSEL) || term)) {
		term =  (a->inex.exgl && m == a->left) ||
			(a->inex.exgr && m == a->right);
		h += term? pwd->UnpPenalty3(insert):	// not intron
			pwd->GapPenalty3(insert);
		if (hi > NEVSEL && insert > intlen)	// pre-intron gap
		    hi += pwd->GapPenalty3(insert - intlen);
		if (hi > NEVSEL && hi >= h) {		// intron
		    if (cigar) {
			if (preint) cigar->push('D', preint);
			cigar->push('N', intlen);
		    }
		    if (vlgar) {
			if (preint) vlgar->push('G', 0, preint);
			if (phs == -1) vlgar->push('S', 0, 1); else
			if (phs ==  1) vlgar->push('S', 0, 2);
			vlgar->push('5', 0, 2);
			vlgar->push('I', 0, intlen - 4);
			vlgar->push('3', 0, 2);
			if (phs == -1) vlgar->push('S', 1, 2); else
			if (phs ==  1) vlgar->push('S', 1, 1);
		    }
		    hb = ha;
		    if (eijnc && rbuf.right - rbuf.left > 1)
			eijnc->push(&rbuf);
		    rbuf.left = rbuf.right + intlen;
		    rbuf.rleft = m;
		    rbuf.sig3 = sig3;
		    h = hi;
		    insert -= (preint + intlen);
		}
		hi = NEVSEL;
		if (insert) {				// post-intron gap
		    if (term && IsTerm(bs[insert - 1])) insert -= 3;
		    if (cigar) cigar->push('D', insert);
		    phs = insert % 3;
		    insert -= phs;
		    if (!((a->inex.exgl && m == a->left) ||
			  (a->inex.exgr && m == a->right))) 
			    fst->gap += gop;
		    fst->unp += insert;
		    if (phs) {	// insertion frame shift
			rbuf.right = n - phs;
			rbuf.rright = m;
			rbuf.iscr = NEVSEL;
		 	prompt(fswarn, a->sqname(), b->sqname(), 
			b->SiteNo(rbuf.right), phs);
			if (eijnc) eijnc->push(&rbuf);
			if (vlgar) vlgar->push('F', 0, phs);
			rbuf.left = n;
			rbuf.rleft = m;
			h += (phs == 1? pwd->GapE1: pwd->GapE2);
			fst->val += (phs == 1? pwd->GapE1: pwd->GapE2);
		    }
		    if (vlgar && insert) vlgar->push('G', 0, insert);
		    gop = insert = intlen = preint = 0;
		}
	    }
	    int	ni = wsk->n - n;
	    if (ni && deletn) {
		if (!(b->inex.exgl && n == b->left)) {
		    h += pwd->GapPenalty3(deletn);
		    fst->gap += 1;
		    fst->unp += deletn;
		}
		as += deletn / 3;
		if ((phs = deletn % 3)) {	// deletion frame shift
		    rbuf.right = n + phs;
		    rbuf.rright = m;
		    rbuf.iscr = NEVSEL;
		    prompt(fswarn, a->sqname(), b->sqname(), b->SiteNo(rbuf.right), -phs);
		    if (eijnc) eijnc->push(&rbuf);
		    if (vlgar) vlgar->push('F', phs, 0);
		    rbuf.left = n;
		    rbuf.rleft = m;
		    h += pwd->ExtraGOP;
		    fst->val += pwd->ExtraGOP;
		    ++as;
		    deletn -= phs;
		    phs = 3 - phs;
		    bs += phs;
		    bb += phs;
		}
		if (vlgar && deletn > 2) vlgar->push('G', deletn / 3, 0);
		deletn = 0;
	    }
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    if (d) {
		if (cigar) cigar->push('M', d);
		if (vlgar) vlgar->push('M', d / 3, d);
		n += d;
		m += d / 3;
		for ( ; d > 2; d -= 3, ++as, bs += 3, bb += 3) {
const		    CHAR*	gs = (cs? cs: bs) + 1;
		    hvl = pwd->sim2(as, gs);
		    fst->val += hvl;
		    h += hvl + (cs? 0: bb[1].sigE);
		    ivl = avst_equal(*as, *gs);
		    if (ivl)	++fst->mch;
		    else	++fst->mmc;
		    cs = 0;
		}
	    }
	    if (i > 0) {
		cs = 0;
		deletn += i;
		if (cigar) cigar->push('I', i);
	    } else if (i < 0) {
const		EXIN*	b3 = bb + (i = -i);
		if (hi <= NEVSEL && i >= IntronPrm.minl && wsk->n < b->right) {	// intron?
const		    CHAR*	cm = 0;
		    VTYPE	sig5m = 0;
		    int	n3 = n + i;
		    int	phs5 = (bb->phs5 == 2)? b3->phs3: bb->phs5;
		    int	phs3 = (b3->phs3 == 2)? bb->phs5: b3->phs3;
		    VTYPE	xm = NEVSEL;
		    VTYPE	xi = NEVSEL;
		 	// potential intron
		    if (phs3 == 2 && phs5 == 2) {	// GTGT....AGAG
			// phs3 = phs5 = -1; 
			nb = n + 1; n3 = nb + i;	// upstream
			sig5m = b->exin->sig53(nb, 0, IE5);
			xm = sig5m + spjcs->spjscr(nb, n3);
			cm = spjcs->spjseq(nb, n3);
			if (usespb) xm += api.match_score(3 * m + 1);
			phs3 = phs5 = 1;
		    }
		    nb = n - phs3; n3 = nb + i;
		    if (isJunct(phs3, phs5)) {
			sig5 = b->exin->sig53(nb, 0, IE5);
			sig3 = b->exin->sig53(nb, n3, IE53);
			xi = sig5 + spjcs->spjscr(nb, n3);
			cs = spjcs->spjseq(nb, n3);
			if (usespb) xi += api.match_score(3 * m - phs3);
			preint = insert;
			if (phs3 == 0) cs = 0;
			if (insert == 0 && phs3 == 1) {
			    VTYPE	hdlt = (pwd->sim2)(as - 1, cs) - hvl;
			    xi += hdlt;
			    fst->val += hdlt;
			    bool match = avst_equal(as[-1], *cs);
			    if (match && !ivl) {
				++fst->mch; --fst->mmc;
			    } else if (!match && ivl) {
				--fst->mch; ++fst->mmc;
			    }
			}
		    }
		    if (xm > xi) {
			xi = xm;
			phs3 = -1;
			nb = n - phs3; n3 = nb + i;
			sig5 = sig5m;
			sig3 = b->exin->sig53(nb, n3, IE53);
			cs = cm;
		    }
		    if (xi > NEVSEL) {
			if (phs3 != -1) cs = 0;
			hi = h + xi;
			intlen = i;
			rbuf.right = nb;		// 5' information
			rbuf.rright = m;
			rbuf.phs = phs3;
			rbuf.iscr = xi;
			rbuf.sig5 = sig5;
			rbuf.escr = h + pwd->GapPenalty3(insert);
			ha = rbuf.escr + xi - sig3;
			rbuf.escr += sig5 - hb;
			rbuf.mch = (int) (fst->mch - pst.mch);
			rbuf.mmc = (int) (fst->mmc - pst.mmc);
			rbuf.gap = (int) (fst->gap - pst.gap);
			rbuf.unp = (int) (fst->unp - pst.unp) / 3;
			pst = *fst;
		    }
		} else	++gop;
		if (i < maxexon) h += SumCodePot(bb, i, 0, pwd);
		else	h = NEVSEL;
		bb = b3;
		bs += i;
		insert += i;
	    }
	    m = wsk->m;
	    n = wsk->n;
	    if (usespb) while (!api.end() && api.lt(m)) ++api;
	    if (num == 0) hi = NEVSEL;
	}
	if (bb + 3 < bb_last) {
	    sig5 = 0;
	    if (algmode.lcl & (16 + 2) && bb[1].sigT > 0)
		sig5 = bb[1].sigT;
	    if (algmode.lcl & (16 + 8) && bb->sig5 > 0
		&& bb->sig5 > bb[1].sigT)
		sig5 = bb->sig5;
	    h += sig5;
	}
	if (eijnc) {
	    rbuf.escr = h - hb;
	    rbuf.iscr = 0;
	    rbuf.sig5 = sig5;
	    rbuf.right = n;
	    rbuf.rright = m;
	    rbuf.mch = (int) (fst->mch - pst.mch);
	    rbuf.mmc = (int) (fst->mmc - pst.mmc);
	    rbuf.gap = (int) (fst->gap - pst.gap);
	    rbuf.unp = (int) (fst->unp - pst.unp) / 3;
	    eijnc->push(&rbuf);
	    rbuf.left = endrng.left;
	    rbuf.right = endrng.right;
	    eijnc->push(&rbuf);
	    eijnc->flush();
	    gsi->noeij = eijnc->size() - 1;
	}
	if (cigar) cigar->flush();
	if (vlgar) {
	    vlgar->push('E', 0, 0);	// dummy
	    vlgar->flush();
	    vlgar->postproc();		// correct match length
	}
	delete spjcs;
	fst->mch /= a->many;
	fst->mmc /= a->many;
	fst->unp /= 3;
	fst->val += pwd->BasicGOP * fst->gap + pwd->BasicGEP * fst->unp;
	if (termcodon) wsk->n += 3;
	return (h);
}

/********************************************************
*
*	unidirectional Hirschberg algorithm
*
********************************************************/

void Aln2h1::hinitH_ng(RVDWU* hhg[], const WINDOW& wdw)
{
	int	n = b->left;
	int	r = b->left - 3 * a->left;
	int	r0 = r;
	int	rr = b->right - 3 * a->left;
	int	dir = a->inex.exgl? DEAD: DIAG;
const	EXIN*	bb = b->exin->score(n);

	RVDWU*	h = hhg[0] + r;
	h->val = (a->inex.exgl && bb[1].sigS > 0)? bb[1].sigS: 0;
	h->dir = dir;
	h->lwr = h->upr = r0;
	if (a->inex.exgl) {
	    if (wdw.up < rr) rr = wdw.up;
	    int	jnc = n;
	    for (int i = 1; ++r <= rr; ++i) {
		++h; ++bb; ++n;
		if (i < 3) {
 		    h->val = (bb[1].sigS > 0)? bb[1].sigS: 0;
		    h->dir = dir;
		    h->lwr = r0;
		} else {
		    *h = h[-3];
		    int	d = n - jnc;
#if TERMGOP
		    if (d == 3) h->val += pwd->BasicGOP;
#else
		    if (d == 3 && !a->inex.exgl) h->val += pwd->BasicGOP;
#endif
		    h->val += pwd->GapExtPen3(d);
	    	    h->val += bb[-2].sigE;
		    h->dir = HORI;
		    VTYPE	x = h[-1].val + pwd->GapW1;
		    if (x > h->val) {*h = h[-1]; h->val = x; h->dir = HOR1;}
		    x = h[-2].val + pwd->GapW2;
		    if (x > h->val) {*h = h[-2]; h->val = x; h->dir = HOR2;}
		}
		h->upr = r;
		VTYPE	x = (bb[1].sigS > 0)? bb[1].sigS: 0;
		if (h->val < x) {
		    h->val = x;
		    h->dir = DEAD;
		    jnc = n;
		    h->lwr = r;
		}
	    }
	}

	r = r0;
	rr = b->left - 3 * a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	h = hhg[0] + r;
	RVDWU*	g = hhg[1] + r;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --g;
	    if (b->inex.exgl) {
		h->val = 0;
		h->dir = DEAD;
		h->upr = h->lwr = r;
	    } else {
		if (i <= 3) {
		    *h = h[i];
		    h->val += pwd->GapPenalty(i);
		    if (i < 3) h->val += pwd->ExtraGOP;
		    h->dir = VERT;
		    *g = *h;
		} else {
		    *h = h[3];
		    *g = g[3];
		    h->val += pwd->GapExtPen3(i);
		    g->val += pwd->BasicGEP;
		}
		h->lwr = g->lwr = r;
	    }
	}
}

RVDWU& Aln2h1::hlastH_ng(RVDWU* hhg[], const WINDOW& wdw)
{
	int	glen[3];
	int	m3 = 3 * a->right;
	int	rw = wdw.lw;
	int	rf = b->left - m3;
	if (rf > rw) rw = rf;
	else	rf = rw;
	RVDWU*	h = hhg[0] + rw;
	RVDWU*	h9 = hhg[0] + b->right - m3;
	RVDWU*	f = h9;
	VTYPE	mx = h9->val;
const	EXIN*	bb = b->exin->score(rw + m3);

	vclear(glen, 3);
	if (a->inex.exgr) {
	    for (int p = 0; h <= h9; ++h, ++bb, ++rf, ++p) {
		if (p == 3) p = 0;
		VTYPE	y = NEVSEL;
		glen[p] += 3;
		if (rf - rw >= 3 && h[-3].dir != DEAD) {
		    VTYPE	x = h[-3].val + bb[-2].sigE 
			+ pwd->GapExtPen3(glen[p]);
#if TERMGOP
		    if (glen[p] == 3) x += pwd->BasicGOP;
#endif
		    if ((algmode.lcl & 2) && bb[-2].sigT > 0 && !(h->dir & SPIN))
			y = h[-3].val + bb[-2].sigT;
		    if (x > h->val) {
			*h = h[-3];
			h->val = x;
			h->dir = HORI;
		    } else if (isnthori(h))
			glen[p] = 0;
		}
		VTYPE	x = h->val;
		if ((algmode.lcl & 16) && bb->sig5 > 0) x += bb->sig5;
		if (x > mx && x >= y) {
		    f = h;
		    mx = x;
		} else if (y > mx) {	// termination codon
		    f = h;
		    *h = h[-3];
		    h->upr += 3;
		    mx = y;
		    h->dir = DEAD;
		}
	    }
	}
	f->val = mx;
	if (b->inex.exgr) {
	    rw = std::min(wdw.up, b->right - 3 * a->left);
	    for (h = hhg[0] + rw; h > h9; --h, --rw) {
		VTYPE	x = h->val + (rw % 3? pwd->ExtraGOP: 0);
		if (f->val < x) {f = h; f->val = x;}
	    }
	}
	return (*f);
}

VTYPE Aln2h1::hirschbergH_ng(int cpos[], 
	const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	RVDWU*	hhg[NOL + 1];
	RVDWU*	hhb[NOL];
	RVDWU*	hb[NOL];		// hhb[k] + r
	RVDWU*	hf[NOD];		// horizontal
	RVDWU	e1[2 * NQUE];
	RVDWU*	e2 = e1 + NQUE;
	RVDWUJ	hl[3][NCAND + 1];	// [phase][candidates]
	int	nx[3][NCAND + 1];	// [phase][candidates]
	RVDWUJ*	maxphl[NOD];
	RVDWU	maxh = black_vdwu;
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

	size_t	bufsiz = pwd->Noll * wdw.width;
	RVDWU*	wbuf = new RVDWU[bufsiz + bufsiz];
	vset(wbuf, black_vdwu, bufsiz);
	vclear(hb, NOL);
	RVDWU*	blackvdwuj = wbuf + bufsiz - 1;	// assume to be const
	hhg[0] = wbuf - wdw.lw + 3;	 
	for (int k = 1; k < pwd->Noll; ++k) hhg[k] = hhg[k-1] + wdw.width;
	for (int k = 0; k < pwd->Noll; ++k) hhb[k] = hhg[k] + bufsiz;
	hinitH_ng(hhg, wdw);

	int	m = a->left;
	int	mm = (a->left + a->right + 1) / 2;
	int	mm3 = 3 * mm;
	if (!a->inex.exgl) --m; // global
	int	n1 = 3 * m + wdw.lw - 1;
	int	n2 = 3 * m + wdw.up;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as) {
	    VTYPE	x, y;
	    bool	internal = (!a->inex.exgr || m < a->right);
	    n1 += 3; n2 += 3;
	    int		n0 = std::max(n1, b->left);
	    int		n9 = std::min(n2, b->right);
	    int		n = n0;
	    int		r = n - 3 * m;
	    vset(e1, black_vdwu, 2 * NQUE);
	    if (!b->inex.exgl && m == a->left) {
		e1[2] = e2[2] = hhg[0][r];
		e1[2].val += pwd->GapW3;
		e2[2].val += pwd->GapW3L;
	    }
const	    CHAR*	bs = b->at(n - 1);
const	    EXIN*	bb = b->exin->score(n);
	    int		k = 0;
	    RVDWU*&	h = hf[0] = hhg[k++] + r;
	    RVDWU*&	g = hf[2] = hhg[k++] + r;
	    RVDWU*&	g2 = hf[4] = dagp? hhg[k++] + r: blackvdwuj;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    if (m == mm) {
		for (int k = 0; k < pwd->Noll; ++k) {
		    hb[k] = hhb[k] + r;
		    hf[2 * k]->ulk = r;
		}
	    }
	    vset((RVDWUJ*) &hl, black_vdwuj, 3 * (NCAND + 1));
	    for (int p = 0; p < 3; ++p)
		for (int l = 0; l <= NCAND; ++l)
		    nx[p][l] = l;
	    int	ncand[3] = {-1, -1, -1};
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d %2d", m, n, h->dir);
		putvar(h->val); putchar('\n');
	    }
#endif
	    int	q = 0;
	    for ( ; ++n <= n9; ++bs) {
		++bb; ++r; ++h; ++g; if (dagp) ++g2;
		RVDWU*&	eq1 = hf[1] = e1 + q;
		RVDWU*&	eq2 = hf[3] = dagp? e2 + q: 0;
		RVDWU	hq = *h;

//	diagonal match
		VTYPE	sigE = bb[-2].sigE;
		RVDWU*	from = h;
		RVDWU*	mx = h;
		if (m == a->left) goto HorizonF;
		if (n > b->left + 2) {
		    h->val += qprof[*bs] + sigE;
		    h->dir = (from->dir & DIAG)? DIAG: NEWD;
		    if (m <= mm) h->ulk = r;
		} else	*h = black_vdwu;

//	vertical gap extention
		y = g[3].val + pwd->BasicGEP;

//	1 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA2;
		} else	g->val = y;

//	2 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA1;
		}

//	nomal deletion
		x = (++from)->val + pwd->GapW3;
		if (x >= g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		} else if (y >= g->val) {
		    *g = g[3];
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		}
		if (g->val >= mx->val) mx = g;

//	long deletion
		if (dagp) {
		  x = from->val + pwd->GapW3L;
		  y = g2[3].val + pwd->LongGEP;
		  if (x >= y) {
		    *g2 = *from;
		    g2->val = x;
		    g2->dir = VERL | (g2->dir & SPJC);
		  } else {
		    *g2 = g2[3];
		    g2->val = y;
		  }
		  if (g2->val >= mx->val) mx = g2;
		}
HorizonF:
//	nomal insertion
		if (n > n0 + 2) {
		    from = h - 3;
		    x = from->val + pwd->GapW3;
		    y = eq1->val += pwd->BasicGEP;
		    if (x > y) {
			*eq1 = *from;
			eq1->val = x;
		    }
		    if (!(eq1->dir & SPF2))	eq1->val += sigE;
		    eq1->dir = (eq1->dir & SPIN) + HORI;
//	long insertion
		    if (dagp) {
			x = from->val + pwd->GapW3L;
			y = eq2->val += pwd->LongGEP;
			if (x > y) {
			    *eq2 = *from;
			    eq2->val = x;
			}
			if (!(eq2->dir & SPF2)) eq2->val += sigE;
			eq2->dir = (eq2->dir & SPIN) + HORL;
			if (eq2->val > mx->val) mx = eq2;
		    }
		}
//	2 nt insertion
		if (n > n0 + 1) {
		    from = h - 2;
		    x = from->val + pwd->GapW2;
		    if (x > eq1->val) {
			*eq1 = *from;
			eq1->val = x;
			eq1->dir = (eq1->dir & SPIN) + HOR2;
		    }
		}
//	1 nt insertion
		from = h - 1;
		x = from->val + pwd->GapW1;
		if (x > eq1->val) {
		    *eq1 = *from;
		    eq1->val = x;
		    eq1->dir = (eq1->dir & SPIN) + HOR1;
		}
		if (eq1->val > mx->val) mx = eq1;
		if (++q == NQUE) q = 0;

//	intron 3' boundary, assume no overlapping signals
		if (internal && isEIJ(bb->phs3)) {
		    int	phs = (bb->phs3 == 2)? -1: bb->phs3;
AccFwd:
		    int		nb = n - phs;
		    VTYPE	sigJ = cip->cip_score(3 * m - phs);
		    int*	pnx = nx[phs + 1];
		    vclear(maxphl, Nod);
		    for (int l = 0; l <= ncand[phs + 1]; ++l) {
			RVDWUJ*	phl = hl[phs + 1] + pnx[l];
			if (nb - phl->jnc < IntronPrm.minl) continue;
			x = phl->val + sigJ + spjcs->spjscr(phl->jnc, nb);
const			CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
			if (phl->dir == 0 && phs) {
			    if (phs == 1) x += qprof[*cs];
			    else	x += qprof[*++cs];
			}
			from = hf[phl->dir];
			if (x > from->val) {
			    from->val = x;
			    maxphl[phl->dir] = phl;
			}
		    }
		    for (int k = 0; k < Nod; ++k) {
			RVDWUJ*	phl = maxphl[k];
			if (!phl) continue;
			from = hf[k];
			from->dir = phl->dir | SPJCI;
			from->upr = std::max(phl->upr, r);
			from->lwr = std::min(phl->lwr, r);
			if (m >= mm) from->ulk = phl->ulk;
			if (ge(from->val, mx->val)) mx = from;
		    }
		    if (bb->phs3 - phs == 3) {	// AGAG
			phs = 1;
			goto AccFwd;
		    }
		}

//	Find optimal path
		y = h->val;
		if (h == mx) {	// diagonal
		    if (LocalR && y > maxh.val)
			maxh = *h;
		} else {	// non-diagonal
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		}
		if (LocalL && h->val <= 0)
		    h->val = h->dir = 0;

//	intron 5' boundary

		if (internal && isEIJ(bb->phs5)) {
		    int		phs = (bb->phs5 == 2)? -1: bb->phs5;
Donor:
		    int	nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(nb, 0, IE5);
		    int		hd = dir2nod[mx->dir & 15];
		    for (int k = (hd == 0 || phs == 1)? 0: 1; k < Nod; ++k) {
			bool	crossspj = phs == 1 && k == 0;
			from = crossspj? &hq: hf[k];
			// An orphan exon is disallowed
			if (!from->dir || (from->dir & SPIN)) continue;
			if (k != hd && !crossspj && hd >= 0) {
			    y = mx->val;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (from->val <= y) continue;	// prune
			}
			x = from->val + sigJ;
			RVDWUJ*	phl = hl[phs + 1];
			int*	pnx = nx[phs + 1];
			int&	nc = ncand[phs + 1];
			int	l = nc < NCAND? ++nc: NCAND;
			while (--l >= 0) {
			    if (x >= phl[pnx[l]].val)
				std::swap(pnx[l], pnx[l + 1]);
			    else
				break;
			}
			if (++l < NCAND) {
			    phl += pnx[l];
			    phl->val = x;
			    phl->jnc = nb;
			    phl->dir = k;
			    phl->upr = from->upr;
			    phl->lwr = from->lwr;
			    if (m > mm) phl->ulk = from->ulk;
			    else if (m == mm) phl->ulk = r;
			} else --nc;
		    }
		    if (bb->phs5 - phs == 3) {	//	GTGT..
			phs = 1;
			goto Donor;
		    }
		}

#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); putvar(g->val); putvar(eq1->val);
		if (dagp) {
		    putvar(g2->val); putvar(eq2->val);
		}
		putvar(hl[2][nx[2][0]].val); 
		putvar(hl[1][nx[1][0]].val);
		putvar(hl[0][nx[0][0]].val);
		printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
		printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		putchar('\n');
	}
#endif

// save variables at the center
		if (m == mm) {
		    for (int k = 0; k < pwd->Noll; ++k) {
			*++hb[k] = *hf[2 * k];
			if (hb[k]->dir & SPJC) hf[2 * k]->ulk = r;
		    }
		    for (int k = 0, j = 0; k < pwd->Noll; ++k, j += wdw.width) {
			hf[2 * k]->upr = hf[2 * k]->lwr = r;
			hf[2 * k]->ulk = r + j;
		    }
		}	// was center
	    }	// end of n-loop
	}	// end of m-loop

	if (!LocalR) maxh = hlastH_ng(hhg, wdw);
	wdwb.up = maxh.upr;
	wdwb.lw = maxh.lwr;
	wdwb.width = wdwb.up - wdwb.lw + 7;
	int	r = maxh.ulk;
	int	k = (r - wdw.lw) / wdw.width;
	if (k)	r -= k * wdw.width;
	RVDWU*	mx = hhb[k] + r;
	b->inex.exgr = k? 1: 0;
	k = k || isvert(mx);
	b->inex.exgl = k? 1: 0;
	if (k) mx->ulk = end_of_ulk;
	int	c = 0;
	if (b->exin->data[r + mm3].phs3 % 3 == 2) {
	    --mm; mm3 -= 3;
	}
	cpos[c++] = mm;
	while (mx->ulk != end_of_ulk && mx->ulk != r) {
	    cpos[c++] = r + mm3;
	    mx = hhb[0] + (r = mx->ulk);
	}
	cpos[c++] = r + mm3;
	cpos[c] = end_of_ulk + k;
	wdwf.up = mx->upr;
	wdwf.lw = std::min(mx->lwr, r);
	wdwf.width = wdwf.up - wdwf.lw + 7;

	delete[] wbuf;
	return (maxh.val);
}

void Aln2h1::pfinitH_ng(RVPD* hh[], const WINDOW& wdw)
{
	int	r = b->left - 3 * a->left;
	int	rr = b->right - 3 * a->left;

	RVPD*	h = hh[0] + r;
	RVPD*	g = hh[1] + r;
	h->val = 0;
	h->dir = DIAG;
	h->ptr = vmf->add(a->left, b->left, 0L);

	r = b->left - 3 * a->left;
	rr = b->left - 3 * a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --g;
	    if (i <= 3) {
		*h = h[i];
		h->val += pwd->GapPenalty(i);
		if (i < 3) h->val += pwd->ExtraGOP;
		h->dir = VERT;
		*g = *h;
	    } else {
		*h = h[3];
		*g = g[3];
		h->val += pwd->GapExtPen3(i);
		g->val += pwd->BasicGEP;
	    }
	}
}

void Aln2h1::pbinitH_ng(RVPD* hh[], const WINDOW& wdw)
{
	int	r = b->right - 3 * a->right;
	int	rr = b->left - 3 * a->right;

	RVPD*	h = hh[0] + r;
	RVPD*	g = hh[1] + r;
	h->val = 0;
	h->dir = DIAG;
	h->ptr = vmf->add(a->right, b->right, 0L);

	r = b->right - 3 * a->right;
	rr = b->right - 3 * a->left;
	if (wdw.up < rr) rr = wdw.up;
	for (int i = 1; ++r <= rr; ++i) {
	    ++h; ++g;
	    if (i <= 3) {
		*h = h[-i];
		h->val += pwd->GapPenalty(i);
		if (i < 3) h->val += pwd->ExtraGOP;
		h->dir = VERT;
		*g = *h;
	    } else {
		*h = h[-3];
		*g = g[-3];
		h->val += pwd->GapExtPen3(i);
		g->val += pwd->BasicGEP;
	    }
	}
}

// intron-less backward extension

VTYPE Aln2h1::back2ward5endH_ng(long *ptr, const WINDOW& wdw)
{
	int	nr[3];
	RVPD*	hh[2];
	RVPD	e1[NQUE];
	VSKLP	maxh = {NEVSEL, 0, 0, 0L};	// start point of traceback
	VTYPE	maxval = NEVSEL;	// for x-drop off
	VTYPE	maxscr = NEVSEL;	// overall

	size_t	bufsiz = 2 * wdw.width;
	RVPD*	jbuf = new RVPD[bufsiz];
	vset(jbuf, black_vpd, bufsiz);
	RVPD*	blackvpd = jbuf + bufsiz - 1;
	hh[0] = jbuf - wdw.lw + 3;
	hh[1] = hh[0] + wdw.width;
	vmf->add(0, 0, 0L);		// Skip 0-th record
	pbinitH_ng(hh, wdw);

	int	m = a->right;
	if (!a->inex.exgr) ++m; 	// global
const	CHAR*	as = a->at(m);
	int	n1 = 3 * m + wdw.lw;
	int	n2 = 3 * m + wdw.up + 1;
	while (--m >= a->left) {
	    --as; n1 -= 3; n2 -= 3;
	    int		n0 = std::min(n2, b->right);
	    int		n9 = std::max(n1, b->left);
	    int		n  = n0;
const	    CHAR*	bs = b->at(n);
const	    EXIN*	bb = b->exin->score(n);
	    int		r = n - 3 * m;
	    vset(e1, black_vpd, NQUE);
	    if (!b->inex.exgr && n == b->right && m == a->right) {
		e1[2] = hh[0][r];
		e1[2].val = pwd->GapW3;
	    }
	    RVPD*	h = hh[0] + r;
	    RVPD*	g = hh[1] + r;
	    RVPD*	mxd = ((h->val + pwd->Vthr) < maxval)? blackvpd: h;
	    int		count3 = 0;
	    for (int p = 0; p < 3; ++p) nr[(n + p) % 3] = n + p;
	    nr[n % 3] = n + 3;

#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m+1, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    int	peak = 0;
	    int	q = 0;		// queue pointer
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    while (--n >= n9) {
		VTYPE	x, y;
		--bs; --r; --h; --g;
		RVPD*	eq1 = e1 + q;
		VTYPE	sigE = (bb--)->sigE;

//	diagonal match
		RVPD*	from = h;
		RVPD*	mx = h;
		if (m == a->right) goto HorizonB;
		if (n < b->right - 2) {
		    h->val += qprof[bs[1]] + sigE;
		    h->dir = isdiag(from)? DIAG: NEWD;
		} else	*h = black_vpd;

//	vertical gap extension
		y = g[-3].val + pwd->BasicGEP;

//	1 nt deletion
		--from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA2;
		} else	g->val = y;

//	2 nt deletion
		--from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val)	{
		    *g = *from;
		    g->val = x;
		    g->dir = SLA1;
		}

//	normal deletion
		x = (--from)->val + pwd->GapW3;
		if (x >= g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		} else if (y >= g->val) {
		    *g = g[-3];
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		}
		if (g->val >= mx->val) mx = g;

HorizonB:
//	normal insertion
		if (n < n0 - 2) {
		    from = h + 3;
		    x = from->val + pwd->GapW3;
		    y = eq1->val += pwd->BasicGEP;
		    if (x > y) {
			*eq1 = *from;
			eq1->val = x;
			if (eq1->dir) eq1->dir = HORI;
		    }
		    eq1->val += sigE;
		}
//	2 nt insertion
		if (n < n0 - 1) {
		    from = h + 2;
		    x = from->val + pwd->GapW2;
		    if (x > eq1->val) {
			*eq1 = *from;
			eq1->val = x;
			eq1->dir = HOR2;
		    }
		}
//	1 nt insertion
		from = h + 1;
		x = from->val + pwd->GapW1;
		if (x > eq1->val) {
		    *eq1 = *from;
		    eq1->val = x;
		    eq1->dir = HOR1;
		}
		if (eq1->val >= mx->val) mx = e1 + q;
		if (++q == NQUE) q = 0;

//	Find optimal path
		if (mx->dir == NEWD)
		    mx->ptr = vmf->add(m + 1, n + 3, mx->ptr);
		y = h->val;
		x = mx->val;
		if (x > maxval) maxval = x;
		if (m == a->left && bb[1].sigS > 0)
		    x += bb[1].sigS;
		if (x > maxh.val) {
		    maxh.val = x;
		    maxh.m = m;
		    maxh.n = n;
		    maxh.p = mx->ptr;
		}
		if (mx->val + pwd->Vthr < maxval) {	// x drop off
		    if  (++count3 == 3 && peak)	{	// end of block
			n1 = n + 3;
			peak = 0;
		    }
		    nr[n % 3] = n;			// before block
		} else {
		   if (isdiag(mx) && mx->val >= mxd->val) {
			mxd = mx;
			if (nr[n % 3] < n2) n2 = nr[n % 3];
			peak = 1;
		    }
		    count3 = 0;
		}
		if (h != mx) *h = *mx;

#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(eq1->val);
		printf(" %6.2f %6.2f", (float) sigE, (float) bb[1].sigS);
		putchar('\n');
	}
#endif
	    } // end of n loop
	    if (!mxd->dir) break;	// no peak
	    if (peak) n1 = n;
	}	// end of m loop
	maxscr = maxh.val;
	if (vmf) *ptr = vmf->add(maxh.m, maxh.n, maxh.p);
	delete[] jbuf;
	return (maxscr);
}

// intron-less forward extension

VTYPE Aln2h1::for2ward3endH_ng(long *ptr, const WINDOW& wdw)
{
	int	nr[3];
	RVPD*	hh[2];
	RVPD	e1[NQUE];
	VSKLP	maxh = {NEVSEL, 0, 0, 0L};	// start point of traceback
	VTYPE	maxval = NEVSEL;	// for x-drop off
	VTYPE	maxscr = NEVSEL;	// overall

	size_t	bufsiz = 2 * wdw.width;
	RVPD*	jbuf = new RVPD[bufsiz];
	vset(jbuf, black_vpd, bufsiz);
	RVPD*	blackvpd = jbuf + bufsiz - 1;
	hh[0] = jbuf - wdw.lw + 3;
	hh[1] = hh[0] + wdw.width;
	vmf->add(0, 0, 0L);		// Skip 0-th record
	pfinitH_ng(hh, wdw);

	int	m = a->left;
	if (!a->inex.exgl) --m; 	// global
	int	n1 = 3 * m + wdw.lw - 1;
	int	n2 = 3 * m + wdw.up;
const	CHAR*	as = a->at(m);
	maxh.val = maxval = NEVSEL;
	for ( ; ++m <= a->right; ++as) {
	    n1 += 3;
	    n2 += 3;
	    int		n0 = std::max(n1, b->left);
	    int		n9 = std::min(n2, b->right);
	    int		n  = n0;
	    int		r = n - 3 * m;
	    int		count3 = 0;
	    vset(e1, black_vpd, NQUE);
	    if (!b->inex.exgl && n == b->left && m == a->left) {
		e1[2] = hh[0][n - 3 * m];
		e1[2].val = pwd->GapW3;
	    }
const	    CHAR*	bs = b->at(n);
const	    EXIN*	bb = b->exin->score(n);
	    RVPD*	h = hh[0] + r;
	    RVPD*	g = hh[1] + r;
	    for (int p = 0; p < 3; ++p) nr[(n - p) % 3] = n - p;
	    RVPD*	mxd = ((h->val + pwd->Vthr) < maxval)? blackvpd: h;
	    nr[n % 3] = n - 3;

#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif

	    int		peak = 0;
	    int		q = 0;		// queue pointer
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x, y;
		++r; ++bb; ++h; ++g;
		RVPD*	eq1 = e1 + q;

//	diagonal match
		VTYPE	sigE = bb[-2].sigE;
		RVPD*	from = h;
		RVPD*	mx = h;
		if (m == a->left) goto HorizonP;
		if (n > b->left + 2) {
		    h->val += qprof[bs[-1]] + sigE;
		    h->dir = isdiag(from)? DIAG: NEWD;
		} else	*h = black_vpd;

//	vertical gap extention
		y = g[3].val + pwd->BasicGEP;

//	1 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA2;
		} else	g->val = y;

//	2 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA1;
		}

//	normal deletion
		x = (++from)->val + pwd->GapW3;
		if (x >= g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		} else	if (y >= g->val) {
		    *g = g[3];
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		}
		if (g->val >= mx->val) mx = g;

HorizonP:
//	normal insertion
		if (n > n0 + 2) {
		    from = h - 3;
		    bool	l = m == a->right && bb[-2].sigT > 0;
		    if (l) {
#if TERMGOP
			x = from->val + (l? bb[-2].sigT: pwd->GapW3);
#else
			x = from->val + (l? bb[-2].sigT: pwd->BasicGEP);
#endif
		    } else
			x = from->val + pwd->GapW3;
		    y = eq1->val += pwd->BasicGEP;
		    if (x > y) {
			*eq1 = *from;
			eq1->val = x;
			if (l) eq1->dir = DEAD;
			else eq1->dir = HORI;
		    }
		    if (!l) eq1->val += sigE;
		}
//	2 nt insertion
		if (n > n0 + 1) {
		    from = h - 2;
		    x = from->val + pwd->GapW2;
		    if (x > eq1->val) {
			*eq1 = *from;
			eq1->val = x;
			eq1->dir = HOR2;
		    }
		}
//	1 nt insertion
		from = h - 1;
		x = from->val + pwd->GapW1;
		if (x > eq1->val) {
		    *eq1 = *from;
		    eq1->val = x;
		    eq1->dir = HOR1;
		}
		if (eq1->val >= mx->val) mx = e1 + q;
		if (++q == NQUE) q = 0;

//	Find optimal path
		if (mx->dir == NEWD)
		    mx->ptr = vmf->add(m - 1, n - 3, mx->ptr);
		y = h->val;
		x = mx->val;
		int	phs = 0;
		if (x > maxh.val) {
		    maxh.val = x;
		    maxh.m = m;
		    maxh.n = n + phs;
		    maxh.p = mx->ptr;
		}
		if (mx->val > maxval) maxval = mx->val;
		if (mx->val + pwd->Vthr < maxval) {	// x drop off
		    if (++count3 == 3 && peak) {	// end of block
			n2 = n - 3;
			peak = 0;
		    }
		    nr[n % 3] = n;
		} else {
		    if (isdiag(mx) && mx->val >= mxd->val) {
			mxd = mx;
			if (nr[n % 3] > n1) n1 = nr[n % 3];
			peak = 1;
		    }
		    count3 = 0;
		}
		if (h != mx) *h = *mx;

#if DEBUG
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y);
		putvar(g->val); putvar(eq1->val);
		printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
		printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		putchar('\n');
	}
#endif

	    } // end of n-loop
	    if (peak) n2 = n; 
	    if (!mxd->dir) break;
	} // end of m-loop
	maxscr = maxh.val;
	if (vmf) *ptr = vmf->add(maxh.m, maxh.n, maxh.p);

	delete[] jbuf;
	return (maxscr);
}

VTYPE Aln2h1::diagonalH_ng()
{
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(b->left + 1);
const	EXIN*	bb = b->exin->score(b->left + 1);
	VTYPE	scr = 0;
	VTYPE	maxh = NEVSEL;
	SKL	wskl;
	int	mL = a->left;
	int	mR = a->right;
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

	for (int m = a->left; m < a->right; ) {
	    scr += pwd->sim2(as, bs) + bb->sigE;
	    ++as; bs += 3; bb += 3; ++m;
	    if (LocalL && scr < 0) {
		scr = 0;
		mL = m;
	    }
	    if (LocalR && scr > maxh) {
		maxh = scr;
		mR = m;
	    }
	}
	wskl.m = mL;
	wskl.n = 3 * (mL - a->left) + b->left;
	mfd->write((UPTR) &wskl);
	wskl.m = mR;
	wskl.n = 3 * (mR - a->left) + b->left;
	mfd->write((UPTR) &wskl);
	return (LocalR? maxh: scr);
}

VTYPE Aln2h1::trcbkalignH_ng(const WINDOW& wdw, bool spj, const RANGE* mc)
{
	long	ptr;

	vmf = new Vmf();
	VTYPE	scr = forwardH_ng(&ptr, wdw, spj, mc);
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

// recursive call 
VTYPE Aln2h1::lspH_ng(const WINDOW& wdw)
{
	int	m = a->right - a->left;
	int	n = b->right - b->left;
	if (!m && !n) return (0);
	INT	aexgl = a->inex.exgl;	// reserve
	INT	aexgr = a->inex.exgr;	// reserve
	INT	bexgl = b->inex.exgl;	// reserve
	INT	bexgr = b->inex.exgr;	// reserve
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
		    pwd->GapExtPen3(n): pwd->UnpPenalty3(n);
	}
	if (wdw.up == wdw.lw) return(diagonalH_ng());
	float	cvol = wdw.lw - b->left + 3 * a->right;
	cvol =  (float) m * n - cvol * cvol / 3;
	if (cvol < MaxVmfSpace || m == 1 || 3 * m + 3 >= n)
	    return (trcbkalignH_ng(wdw));

	WINDOW	wdwl;
	WINDOW	wdwr;
	RANGE	rng[2];			// reserve
	int	cpos[12];

	save_range(seqs, rng, 2);
	VTYPE	scr = hirschbergH_ng(cpos, wdw, wdwl, wdwr);
	SKL	wskl = {cpos[0], 0};
	int	c = 1;
	for ( ; cpos[c] < end_of_ulk; ++c) {
	    wskl.n = cpos[c];
	    mfd->write((UPTR) &wskl);
	}
	a->inex.exgr = 0;
	INT	bexg = b->inex.exgl;	// temporary
	b->inex.exgl = bexgl;
	int	r = b->left - 3 * a->left;
	if (r < wdwl.lw && bexgl) b->left = 3 * a->left + wdwl.lw;
	if (r > wdwl.lw && aexgl) a->left = (b->left - wdwl.lw) / 3;
	a->right = cpos[0];
	b->right = cpos[c - 1];
	lspH_ng(wdwl);
	rest_range(seqs, rng, 2);
	a->inex.exgl = 0;
	a->inex.exgr = aexgr;
	b->inex.exgr = bexgr;
	b->inex.exgl = bexg;
	r = b->right - 3 * a->right;
	if (r > wdwr.up && bexgr) b->right = 3 * a->right + wdwr.up;
	if (r < wdwr.up && aexgr) a->right = (b->right - wdwr.up) / 3;
	a->left = cpos[0];
	b->left = cpos[1];
	r = b->left - 3 * a->left;
	if (r > wdwr.lw && cpos[c] > end_of_ulk) {
	    mfd->write((UPTR) &wskl);
	    a->left += (r - wdwr.lw) / 3;
	}
	lspH_ng(wdwr);
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexgl;
	b->inex.exgl = bexgl;
	return scr;
}

VTYPE Aln2h1::shortcutH_ng(int ovr, BOUND& bab)
{
	int	margin = IntronPrm.minl;
	int	interval = b->right - b->left - 2 * margin;
	interval = (interval > 0)? interval / 3 * 3: 0;
	RANGE	rng = {b->left + margin, b->left + margin + interval};

	VTYPE	scr = 0;
	ovr = (ovr > 0? 0: ovr) - 3;
	scr -= creepback(ovr, 0, bab);
	scr -= creepfwrd(ovr, 0, bab);
	int	alen = a->right - a->left;
	int	sh = alen / 2;
	if (alprm.sh < 0) {
	    float f = -alprm.sh;
	    if (f > 1.) f /= 100;
	    if (f < 0.5) sh = int(alen * f);
	} else if (alprm.sh < sh)
	    sh = int(alprm.sh);
	int	minsh = alen - margin / 3;
	if (sh < minsh) sh = minsh;
	WINDOW	wdw;
	stripe31(seqs, &wdw, sh);

	a->inex.exgl = b->inex.exgl = 0;	// global
	a->inex.exgr = b->inex.exgr = 0;	// global
	scr += trcbkalignH_ng(wdw, true, interval? &rng: 0);
	return (scr);
}

VTYPE Aln2h1::openendH_ng(int cmode)
{
//printf("%s\t%s\n", (*a->sname)[0], (*b->sname)[0]);
	vmf = new Vmf();
	long	ptr;
	WINDOW	wdw;
	stripe31(seqs, &wdw, alprm.sh);
	VTYPE	scr = (cmode == 1)?
	    back2ward5endH_ng(&ptr, wdw):
	    for2ward3endH_ng(&ptr, wdw);

	SKLP    sv = {0, 0, ptr};
	while (sv.p) {
	    vmf->readvmf(&sv, sv.p);
	    mfd->write((UPTR) &sv);
	}
	delete vmf;
	vmf = 0;
	return (scr);
}

VTYPE Aln2h1::backforth(int ovr, BOUND& lub)
{
	int	ocodon = ovr / 3;
	VTYPE*	bscr = new VTYPE[ocodon + 1];
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(b->left + 1);
const	EXIN*	bb = b->exin->score(b->left + 1);
	VTYPE	scr = bscr[ocodon] = 0;
	int	i = ocodon;
	int	m = a->left;
	int	n = b->left;
	while (--i >= 0 && (m -= 3) >= lub.la && (n -= 3) >= lub.lb)
	    bscr[i] = scr += pwd->sim2(--as, bs -= 3) + (bb -= 3)->sigE;
	VTYPE	maxscr = scr;
	scr = 0;
	int	ii = ++i;
	m = a->right + i;
	n = b->right + 3 * i;
	as = a->at(m);
	bs = b->at(n);
	bb = b->exin->score(n + 1);
	for ( ; i++ < ocodon && m++ < lub.ua && n < lub.ub; n += 3, bs += 3) {
	    scr += pwd->sim2(as++, bs);
	    if ((bscr[i] += scr) > maxscr) {
		maxscr = bscr[i]; ii = i;
	    }
	}
	SKL	skl = {a->right + ii, b->right + 3 * ii};
	mfd->write((UPTR) &skl);
	int dr = (b->right - 3 * a->right) - (b->left - 3 * a->left);
	if (dr >= 0) skl.n -= dr;
	else	skl.m -= (dr = -dr) / 3;
	if (skl.n >= 0) mfd->write((UPTR) &skl);
	delete[] bscr;
	return (maxscr + pwd->GapPenalty3(dr));
}

VTYPE Aln2h1::cds5end(const SKL* pskl)
{
const	EXIN*	bb = b->exin->score(pskl->n + 1);
	VTYPE	scr = 0;
	VTYPE	maxscr = 0;
	SKL	skl = *pskl;
	mfd->write(pskl);

	for (int y = pskl->n; y > b->left; y -= 3) {
	    if (bb->sigS > 0) scr += bb->sigS;
	    if (scr > maxscr) {
		maxscr = scr;
		skl.m = y;
	    }
	    if (bb->sigS > 0 || scr + pwd->Vthr < 0) break;
	    bb -= 3;
	    scr += bb->sigE + pwd->BasicGEP;
	}
	if (skl.n != pskl->n) mfd->write(&skl);
	return (maxscr);
}

VTYPE Aln2h1::cds3end(const SKL* pskl, VTYPE maxscr)
{
const	EXIN*	bb = b->exin->score(pskl->n + 1);
	VTYPE	scr = maxscr;
	SKL	skl = *pskl;
	mfd->write(pskl);

	for (int y = pskl->n; y < b->right; y += 3, bb += 3) {
	    if (bb->sigT > 0) scr += bb->sigT;
	    else	scr += bb->sigE + pwd->BasicGEP;
	    if (scr > maxscr) {
		maxscr = scr;
		skl.n = y + 3;
	    }
	    if (bb->sigT > 0 || scr + pwd->Vthr < 0) break;
	}
	if (skl.n != pskl->n) mfd->write(&skl);
	return (maxscr);
}

static void addsigEjxt(JUXT* jxt, int num, const Seq* b)
{
	if (!b->exin) return;
	for ( ; num--; ++jxt) {
const	    EXIN*	bb = b->exin->score(jxt->jy + 1);
	    VTYPE	scr = 0;
	    for (int i = 0; i < jxt->jlen; ++i) {
		scr += bb->sigE;
		bb += 3;
	   }
	   jxt->jscr += scr;
	}
}

VTYPE Aln2h1::creepback(int ovr, VTYPE bscr, BOUND& lub)
{
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(b->left + 1);
const	EXIN*	bb = b->exin->score(b->left + 1);
	VTYPE	dscr = 0;
	while (a->left > lub.la && b->left > lub.lb
		&& (ovr < 0 || dscr < bscr)) {
	    --as; bs -= 3; bb -= 3;
	    dscr += pwd->sim2(as, bs) + bb->sigE;
	    a->left--; b->left -= 3;
	    if ((ovr += 3) == 0) bscr += dscr;
	}
	return (dscr);
}

VTYPE Aln2h1::creepfwrd(int& ovr, VTYPE bscr, BOUND& lub)
{
const	CHAR*	as = a->at(a->right);
const	CHAR*	bs = b->at(b->right + 1);
const	EXIN*	bb = b->exin->score(b->right + 1);
	VTYPE	dscr = 0;
	while (a->right < lub.ua && b->right < lub.ub
		&& (ovr < 0 || dscr < bscr)) {
	    dscr += pwd->sim2(as, bs) + bb->sigE;
	    ++as; bs += 3; bb += 3;
	    a->right++; b->right += 3;
	    if ((ovr += 3) == 0) bscr += dscr;
	}
	return (dscr);
}

// indel-free version of pincsersH within overlapped region
// doubly conted alignment score is cancelled
// assume agap <= 1

bool Aln2h1::indelfreespjH(int agap, VTYPE& iscr)
{
const	int	ilen = b->right - b->left - 3 * agap;
const	int	d5 = b->left + 3 * agap - 2;	// donor 5' end
const	int	d3 = b->left + 2;		// donor 3' end
const	int	a5 = b->right - 2;		// accpt 5' end
	agap = 2 - agap;
	int	n = min3(a->left, b->left / 3, agap + expected_overlap_ext);
	VTYPE*	bw = (n + 1 > expected_max_overlap)? new VTYPE[n + 1]: backward;
const	CHAR*	as = a->at(a->left);		// cancel score of
const	CHAR*	bs = b->at(b->left + 1);	// donor side overlap
const	CHAR*	ds = bs + ilen;
	EXIN*	bd = b->exin->score(b->left + 1);
	int	i = 0;
	bw[0] = 0;
	for (VTYPE v = 0; ++i < n && (*(bs -= 3) == *(ds -= 3) || i < agap); )
	    bw[i] = v += pwd->sim2(--as, bs) + (bd -= 3)->sigE;
	vreverse(bw, i);
	SKL	skl = {0, 0};
	int	phs53 = 0;
	iscr = NEVSEL;

	int	m = a->right - 1;
	bool	all_mch = true;
	for (int nt = 0; nt < 2 && iscr == NEVSEL; ++nt) {
	  int	phs = 1;
	  as = a->at(m = a->right - 1);
	  bs = b->at(a5);
	  bd = b->exin->score(n = d5);
	  EXIN*	ba = b->exin->score(a5);
	  for (VTYPE v = i = 0; n <= d3; ++n, ++bd, ++ba) {
	    if (nt || b->exin->isCanon(n, n + ilen)) {
	  	bool	mch = true;
		VTYPE	x = b->exin->sig53(n, n + ilen, IE5P3)
			 + cip->cip_score(3 * m + phs);
		VTYPE	y = x - v - bw[i];
		if (phs) {
const		    CHAR*	cs = spjcs->spjseq(n, n + ilen);
		    if (phs == 1) ++cs;
		    if ((mch = algmode.crs || avst_equal(*as, *cs)))
			y += pwd->sim2(as, cs);
		    else
			y = NEVSEL;

		}
		if (y > iscr) {
		    skl.n = n - phs;
		    skl.m = m;
		    iscr = y;
		    phs53 = -phs;
		    all_mch = mch;
		}
	    }
	    if (++phs == 0) ++i;
	    else if (phs == 2) {++m; phs = -1;}
	    else v += pwd->sim2(++as, bs += 3) + ba[1].sigE;
	    		// cancel score of acc side
	  }
	}
	if (bw != backward) delete[] bw;
	if (algmode.crs || (all_mch && iscr > NEVSEL)) {
	    mfd->write(&skl);
	    bd = b->exin->score(skl.n);
	    bd->phs5 = phs53;
	    skl.n += ilen;
	    mfd->write(&skl);
	    bd = b->exin->score(skl.n);
	    bd->phs3 = phs53;
	    iscr += pwd->IntPen->Penalty(ilen);
	    return (true);
	} else
	    return (false);
}

//	heuristic search for short exon

int Aln2h1::nearest5ss(int from, BOUND& bab)
{
	int	retry = 0;
const	EXIN*	bbs[n_nearest_ss];
const	EXIN*	bbb = b->exin->score(from);
ss5:
	int	pa = a->left;
	int	ta = std::max(bab.la, pa - max_dist2ss);
	int	nn = from;
const	EXIN*	bb = bbb;
	int	nss = 0;
	for (int p = 0; pa > ta && nn > bab.lb; --nn, --bb) {
	    if (bb->sig5 > 0 || (retry && bb->phs5 == 0)) {
		if (nss < n_nearest_ss) bbs[nss++] = bb;
		else	break;
	    } else if ((++p % 3) == 0) --pa;
	}
	int	dd = nss? bbb - bbs[nss - 1]: max_dist2ss;
	pa = a->left;
	ta = std::min(pa + max_dist2ss, bab.ua);
	nn = from;
	bb = bbb;
	for (int p = 1; pa < ta && ++nn < bab.ub; ) {
	    if ((++bb)->sig5 > 0 || (retry && bb->phs5 == 0)) {
		if (nss < n_nearest_ss) bbs[nss++] = bb;
		else if (p < dd || bb->sig5 > bbs[nss - 1]->sig5)
		    bbs[nss - 1] = bb;
		if (nss == n_nearest_ss) break;
	    } else if ((++p % 3) == 0) ++pa;
	}
	if (nss == 0) {
	    if (!retry++) goto ss5;
	    else	return (0);
	}
	if (nss == 2) {
	    int	d0 = abs(bbs[0] - bbb);
	    int	d1 = abs(bbs[1] - bbb);
	    if (d0 > d1 || (d0 == d1 && bbs[0]->sig5 < bbs[1]->sig5))
		std::swap(bbs[0], bbs[1]);
	    if (bbs[0]->sig5 > bbs[1]->sig5) --nss;
	}
	for (int n = 0; n < nss; ++n)
	    ss[n] = from + bbs[n] - bbb;
	return (nss);
}

int Aln2h1::nearest3ss(int from, BOUND& bab)
{
	int	retry = 0;
const	EXIN*	bbs[n_nearest_ss];
const	EXIN*	bbb = b->exin->score(from);
ss3:
	int	pa = a->right;
	int	nn = from;
	int	ta = std::max(bab.la, pa - max_dist2ss);
const	EXIN*	bb = bbb;
	int	nss = 0;
	for (int p = 0; pa > ta && nn > bab.lb; --nn, --bb) {
	    if (bb->sig3 > 0 || (retry && bb->phs3 == 0)) {
		if (nss < n_nearest_ss) bbs[nss++] = bb;
		else	break;
	    } else if ((++p % 3) == 0) --pa;
	}
	int	dd = nss? bbb - bbs[nss - 1]: max_dist2ss;
	pa = a->right;
	ta = std::min(pa + max_dist2ss, bab.ua);
	nn = from;
	bb = bbb;
	for (int p = 1; pa < ta && ++nn < bab.ub; ) {
	    if ((++bb)->sig3 > 0 || (retry && bb->phs3 == 0)) {
		if (nss < n_nearest_ss) bbs[nss++] = bb;
		else if (p < dd || bb->sig3 > bbs[nss - 1]->sig3)
		    bbs[nss - 1] = bb;
		if (nss == n_nearest_ss) break;
	    } else if ((++p % 3) == 0) ++pa;
	}
	if (nss == 0) {
	    if (!retry++) goto ss3;
	    else	return (0);
	}
	if (nss == 2) {
	    int	d0 = abs(bbs[0] - bbb);
	    int	d1 = abs(bbs[1] - bbb);
	    if (d0 > d1 || (d0 == d1 && bbs[0]->sig3 < bbs[1]->sig3))
		std::swap(bbs[0], bbs[1]);
	    if (bbs[0]->sig3 > bbs[1]->sig3) --nss;
	}
	for (int n = 0; n < nss; ++n)
	    ss[n] = from + bbs[n] - bbb;
	return (nss);
}

// no gap with mismatches

VTYPE Aln2h1::micro_exon(BOUND& bab)
{
	int	nl = nearest5ss(b->left, bab);
	if (nl == 0) return (NEVSEL);
	int	l = ss[0];
	int	nr = nearest3ss(b->right, bab);
	if (nr == 0) return (NEVSEL);
	int	r = ss[0];
	RANGE	arng;
	RANGE	brng;
	a->saverange(&arng);
	b->saverange(&brng);

	int	d5 = b->left - l;
	d5 = (d5 >= 0? d5 + 1: d5 - 1) / 3;
	a->left -= d5;
	b->left -= 3 * d5;
	d5 = b->left - l;	// -1, 0, 1
	int	d3 = b->right - r;
	d3 = (d3 >= 0? d3 + 1: d3 - 1) / 3;
	a->right -= d3;
	b->right -= 3 * d3;
	d3 = b->right - r;	// -1, 0, 1
	int	alen = a->right - a->left;
	if (alen <= 0) {	// no intron in between
	    VTYPE	scr = 0;
	    if (indelfreespjH(alen, scr)) return (scr);
	    a->restrange(&arng);
	    b->restrange(&brng);
	    return (NEVSEL);
	}
const	CHAR*	ts = a->at(a->right);
const	CHAR*	ss = a->at(a->left);
	if (d5 == 1) --ss;
	if (d3 == 1) --ts;
	VTYPE	maxscr = NEVSEL;
	int	f = -1;
	int	n = b->left + IntronPrm.minl;
	int	n9 = b->right - 3 * alen - IntronPrm.minl;
	int	cds = 3 * alen + d5 - d3;
	EXIN*	bb5 = b->exin->score(n);	// 5' end of internal exon
	EXIN*	bb3 = b->exin->score(n + cds);	// 3' end of internal exon
	for ( ; n < n9; ++n, ++bb5, ++bb3) {
	    if (bb5->phs3 || bb3->phs5) continue;
const	    CHAR*	as = ss;
const	    CHAR*	bs = b->at(n + d5 + 1);
	    VTYPE	scr = 0;
	    if (d5) {
const		CHAR*	cs = spjcs->spjseq(l, n);
		if (d5 == -1) {++cs; bs += 3;}
		scr += pwd->simmtx->mtx[*as++][*cs];
	    }
	    if (d3) {
const		CHAR*	cs = spjcs->spjseq(n + cds, r);
		if (d3 == -1) ++cs;
		scr += pwd->simmtx->mtx[*ts][*cs];
	    }
	    for ( ; as < ts; bs += 3)
		scr += pwd->simmtx->mtx[*as++][*bs];
	    scr = alprm2.w * scr + bb5->sig3 + bb3->sig5 +
		pwd->IntPen->Penalty(n - l) +
		pwd->IntPen->Penalty(r - n);
	    if (scr > maxscr) {
		maxscr = scr;
		f = n;
	    }
	}
	if (f < 0) {
	    a->restrange(&arng);
	    b->restrange(&brng);
	    return (NEVSEL);
	}
	SKL	skl = {a->left, b->left};
	mfd->write(&skl);
	skl.n = f + d5;
	mfd->write(&skl);
	skl.n = f + cds + d3;
	skl.m += alen;
	mfd->write(&skl);
	skl.n = b->right;
	mfd->write(&skl);
	return (maxscr);
}

// no gap with mismatches

int Aln2h1::first_exon_wmm(int d3, VTYPE& retscr, bool& pm, int nss)
{
	int	r = b->right - d3;
	int	n = std::max(b->left, b->right - 3 * a->right - IntronPrm.minl);
	int	nd = n + 3 * a->right - d3;
	EXIN*	bbi = b->exin->score(n + 1);
	EXIN*	bb5 = b->exin->score(nd);
const	CHAR*	ts = a->at(a->right);
const	CHAR*	ss = a->at(a->left);
	VTYPE	pmch = 0;
	VTYPE	maxscr = NEVSEL;
	VTYPE	dfact = nss > 1? alprm2.w - 1: 0;
	for (const CHAR* as = ss; as < ts; ++as)
	    pmch += pwd->simmtx->mtx[*as][*as];
	if (d3 == -1) pmch += pwd->simmtx->mtx[*ts][*ts];
	else if (d3 == 1) --ts;
	int	f = -1;
	for ( ; n >= b->left; --n, --nd, --bbi, --bb5) {
	    if (bbi->sigS < 0 || !b->exin->isCanon(nd, r)) continue;
	    const CHAR*	bs = b->at(n + 1);
	    VTYPE	mchscr = 0;
	    for (const CHAR* as = ss; as < ts; bs += 3)
		mchscr += pwd->simmtx->mtx[*as++][*bs];
	    if (d3) {
const		CHAR*	cs = spjcs->spjseq(nd, r);
		if (d3 == -1) ++cs;
		mchscr += pwd->simmtx->mtx[*ts][*cs];
	    }
	    VTYPE	scr = alprm2.w * mchscr + bbi->sigS + 
		bb5->sig5 + pwd->IntPen->Penalty(r - nd);
	    if (scr > maxscr) {
		f = n;
		maxscr = scr;
		retscr = maxscr - dfact * mchscr;
		pm = mchscr == pmch;
		if (pm && ((r - nd) > IntronPrm.mode)) break;
	    }
	    if (((r - nd) % IntronPrm.maxl) == 0 && maxscr > NEVSEL) break;
	}
	return (f);
}

// exact matches

VTYPE Aln2h1::first_exon(BOUND& bab)
{
	int	nss = nearest3ss(b->right, bab);
	if (nss == 0) {			// no acceptor site nearby
	    SKL	skl = {a->right, b->right};
	    return (cds5end(&skl));
	}
	RANGE	arng[2] = {{0, 0}, {0, 0}};
	RANGE	brng[2] = {{0, 0}, {0, 0}};
	a->saverange(arng);
	b->saverange(brng);
	int	maxf = -1;
	VTYPE	maxscr = NEVSEL;
	int	nn = 1;

	for (int n = 0; n < nss; ++n) {
	    if (n) {
	 	a->restrange(arng);
		b->restrange(brng);
	    }
	    int	r = ss[n];
	    int	d3 = b->right - r;
	    d3 = (d3 >= 0? d3 + 1: d3 - 1) / 3;
	    a->right -= d3;
	    b->right -= 3 * d3;
	    if (a->right == 0 || b->right < 3) {
		SKL	skl = {a->right, b->right};
		return (cds5end(&skl));
	    }
	    d3 = b->right - r;	// -1, 0, 1

	    if (algmode.crs || a->right < 2) {
		if (n == 0) {
		    a->saverange(arng + 1);
		    b->saverange(brng + 1);
		}
		VTYPE	scr = NEVSEL;
		bool	pm = false;
		int	f = first_exon_wmm(d3, scr, pm, nss);
	 	if (scr > maxscr) {
		    maxscr = scr;
		    maxf = f;
		    nn = n;
		    if(pm) break;
		}
	    } else {
		int	cds = 3 * a->right - d3;
		if (d3 == 1) --a->right;
		BoyerMoore	bm(b, a, -3);
		int	l = std::max(b->left, b->right - IntronPrm.maxl);
const		CHAR*	as = a->at(a->right);
		while (!bm.finished()) {
		    int	f = bm.nexthit3(l) - 1;
		    if (f >= 0) {
			int	nd = f + cds;
			if (b->exin->isCanon(nd, r)) {
			    if (d3) {
const				CHAR*	cs = spjcs->spjseq(nd, r);
				if (d3 == -1) ++cs;
				if (!avst_equal(*as, *cs)) f = -1;
	    		    }
			    if (f >= 0) {
				EXIN*	bbi = b->exin->score(f + 1);
				VTYPE	scr = pwd->IntPen->Penalty(r - nd) +
				    b->exin->sig53(nd, r, IE5P3) + bbi->sigS;
				if (scr > maxscr) {
				    maxscr = scr;
				    maxf = f;
				}
			    }
			}
	 	    }
	 	    if (bm.scanned(l)) {
			if (maxf < 0) l = std::max(b->left, l - IntronPrm.maxl);
			else	break;
		    }
		}
		if (d3 == 1) ++a->right;
		if (maxf >= 0) {nn = 1; break;}
	      }
	}
	if (maxf < 0) {
	    a->restrange(arng);
	    b->restrange(brng);
	    return (NEVSEL);
	}
	if (nn == 0) {
	    a->restrange(arng + 1);
	    b->restrange(brng + 1);
	}
	b->left = maxf;
	SKL	skl = {a->left, b->left};
	mfd->write(&skl);
	skl.m = a->right;
	skl.n = b->left + 3 * a->right;
	mfd->write(&skl);
	skl.n = b->right;
	mfd->write(&skl);
	return (0);
}

// no gap with mismatches

int Aln2h1::last_exon_wmm(int d5, VTYPE& retscr, bool& pm, int nss)
{
	int	l = b->left - d5;
	int	alen = a->right - a->left;
	int	rr = b->right - 3 * alen - d5 - 1;
	int	n = b->left + IntronPrm.minl;
	EXIN*	bb3 = b->exin->score(n);
	EXIN*	bbt = b->exin->score(n + 3 * alen + d5 + 1);
const	CHAR*	ts = a->at(a->right);
const	CHAR*	ss = a->at(a->left);
	if (d5 == 1) --ss;
	VTYPE	maxscr = NEVSEL;
	VTYPE	dfact = nss > 1? alprm2.w - 1: 0;
	VTYPE	pmch = 0;			// perfect match score
	for (const CHAR* as = ss; as < ts; ++as)
	    pmch += pwd->simmtx->mtx[*as][*as];
const	CHAR*	bss = b->at(n + d5 + 1);
	int	f = INT_MIN / 2;
	for ( ; n < rr; ++n, ++bb3, ++bbt, ++bss) {
	    if (bbt->sigT < 0 || !b->exin->isCanon(l, n)) continue;
	    VTYPE	mchscr = 0;
const	    CHAR*	as = ss;
const	    CHAR*	bs = bss;
	    if (d5) {
const		CHAR*	cs = spjcs->spjseq(l, n);
		if (d5 == -1) {++cs; bs += 3;}
		mchscr += pwd->simmtx->mtx[*as++][*cs];
	    }
	    for ( ; as < ts; bs += 3)
		mchscr += pwd->simmtx->mtx[*as++][*bs];
	    VTYPE	scr = alprm2.w * mchscr + bb3->sig3 + 
		bbt->sigT + pwd->IntPen->Penalty(n - l);
	    if (scr > maxscr) {
		f = n;
		maxscr = scr;
		retscr = maxscr - dfact * mchscr;
		pm = mchscr == pmch;
		if (pm && ((n - l) > IntronPrm.mode)) break;
	    }
	    if (((n - l) % IntronPrm.maxl) == 0 && maxscr > NEVSEL) break;
	}
	return (f + d5);
}

// exact matches

VTYPE Aln2h1::last_exon(BOUND& bab, SKL* pskl)
{
	if (a->right == a->len) ++bab.ua;
	int	nss = nearest5ss(b->left, bab);
	if (nss == 0) return cds3end(pskl);	// no splice donor site nearby

	RANGE	arng[2] = {{0, 0}, {0, 0}};
	RANGE	brng[2] = {{0, 0}, {0, 0}};
	a->saverange(arng);
	b->saverange(brng);
	int	maxf = -1;
	VTYPE	maxscr = NEVSEL;
	int	alen = 0;
	int	alen_0 = 0;
	int	nn = 1;

	for (int n = 0; n < nss; ++n) {
	    if (n) {
		a->restrange(arng);
		b->restrange(brng);
	    }
	    int	l = ss[n];
	    int	d5 = b->left - l;
	    d5 = (d5 >= 0? d5 + 1: d5 - 1) / 3;
	    a->left -= d5;
	    b->left -= 3 * d5;
	    d5 = b->left - l;	// -1, 0, 1
	    alen = a->right - a->left;

	    if (alen < 0) { 		// xxx'tg|gt...ag|a'
const		CHAR*	bs = b->at(brng->left + 1);
		if (acodon[*bs] != 'W') return cds3end(pskl);
		mfd->write(pskl);
		return (0);
	    } else if (alen == 0 && d5 == 0) {
		mfd->write(pskl);		// xxx'|gt
		return (0);
	    } else if (alen == 0 && d5 == -1) {	// XXX'T|gt...ag|AR or |AG
const		CHAR*	bs = b->at(brng->left);
		if (ncodon[*bs] != 'T') return cds3end(pskl);
		mfd->write(pskl);
		return (0);
	    }

	    if (algmode.crs || alen < 2) {
		if (n == 0) {
		    a->saverange(arng + 1);
		    b->saverange(brng + 1);
		    alen_0 = alen;
		}
		VTYPE	scr = NEVSEL;
		bool	pm = false;
		int	f = last_exon_wmm(d5, scr, pm, nss);
		if (scr > maxscr) {
		    maxscr = scr;
		    maxf = f;
		    nn = n;
		    if (pm) break;
		}
	    } else {
const		CHAR*	as = a->at(a->left);
		if (d5 < 0) {
		    ++a->left; --alen; d5 += 3;
		}
		BoyerMoore	bm(b, a, 3);
		int	r = std::min(b->right, l + IntronPrm.maxl);
		if (d5 == 1) --as;
		while (!bm.finished()) {
		    int	f = bm.nexthit3(-1, r) - 1;
		    if (f >= 0) {
			int	na = f - d5;	// acceptor
			if (b->exin->isCanon(l, na)) {
			    if (d5) {
const				CHAR*	cs = spjcs->spjseq(l, na);
				if (d5 != 1) ++cs;
				if (!avst_equal(*as, *cs)) f = -1;
			    }
			    if (f >= 0) {
				EXIN*	bb = b->exin->score(f + 3 * alen + 1);
				if (bb->sigT > 0) {
				    maxscr = pwd->IntPen->Penalty(na - l) +
					b->exin->sig53(l, na, IE5P3);
				    maxf = f;
				    break;
				}
			    }
			}
	 	    }
		    if (bm.scanned(r)) {
			if (maxf < 0) r = std::min(b->right, r + IntronPrm.maxl);
			else	break;
		    }
		}	// end of BM scan
		if (d5 == 2) {--a->left; ++alen; d5 -= 3; maxf -= 3;}
		if (maxf >= 0) {nn = 1; break;}
	    }	// BM
	}

	if (maxf < 0) {
	    a->restrange(arng);
	    b->restrange(brng);
	    return (NEVSEL);
	}
	if (nn == 0) {
	    a->restrange(arng + 1);
	    b->restrange(brng + 1);
	    alen = alen_0;
	}
	SKL	skl = {a->left, b->left};
	mfd->write(&skl);
	skl.n = maxf;
	mfd->write(&skl);
	skl.m = a->right;
	skl.n = maxf + 3 * alen;
	return cds3end(&skl, maxscr);
}

VTYPE Aln2h1::interpolateH(INT level, int agap, int bgap, int cmode, 
	const JUXT* wjxt, BOUND& bab, int wlmt)
{
	int	ovr = std::min(3 * agap, bgap);
	int	dgap = bgap - 3 * agap;
	bool	no_rec = ovr <= (cmode == 3? wlmt: 3 * wlmt);// no recurrsion
	VTYPE	iscore = NEVSEL;
	VTYPE	scr = 0;
	Mfile*	save_mfd = 0;

	if (dgap == 0) {
	    if (agap == 0) mfd->write((UPTR) wjxt);
	    else	scr += diagonalH_ng();
	    iscore = 0;
	} if (cmode == 1 && no_rec) {			// 5' end
	    iscore = (wjxt && (agap <= 0 || algmode.crs > 1))?
		cds5end((SKL*) wjxt): first_exon(bab);
	} else if (cmode == 2 && no_rec) {
	    SKL	skl = {a->left, b->left};
	    if (agap <= 0) {
		mfd->write(&skl);
		EXIN*	bb = b->exin->score(b->left + 1);
		if (bb->sigT > 0) {
		    skl.n += 3;
		    mfd->write(&skl);
		}
		iscore = 0;
	    } else
		iscore = algmode.crs > 1? 
		    cds3end(&skl): last_exon(bab, &skl);
	} else if (cmode == 3 && agap <= 1 && dgap >= IntronPrm.minl &&
	    indelfreespjH(agap, iscore)) {		// indel-free spj
	} else if (cmode == 3 && no_rec && dgap >= IntronPrm.minl) {
	    iscore = micro_exon(bab);
	    if (iscore == NEVSEL && agap < IntronPrm.elmt)
		iscore = shortcutH_ng(ovr, bab);	// shortcut
	} else if (agap <= 0 && dgap < IntronPrm.minl) {
	    iscore = backforth(-ovr, bab);		// ordinary gap
	} else if (dgap < IntronPrm.minl) {
	    scr -= creepback(ovr, slmt, bab);
	    scr -= creepfwrd(ovr, slmt, bab);
	    WINDOW	wdw;
	    bool	abnormal = b->left > b->right;
	    if (abnormal) std::swap(b->left, b->right);
	    stripe31(seqs, &wdw, std::min(alprm.sh, abs(dgap) + 3));
	    iscore = trcbkalignH_ng(wdw, false);	// ordinary alignment
	    if (abnormal) std::swap(b->left, b->right);
	} else if (level + 1 < algmode.qck) {
	    save_mfd = new Mfile(*mfd);			// recursive search
	    island = false;
	    iscore = seededH_ng(level + 1, cmode, bab);
	    if (iscore < 0 && island) iscore = NEVSEL;
	}
	if (iscore > NEVSEL) scr += iscore;
	else {
	    if (save_mfd) std::swap(mfd, save_mfd);
	    scr -= creepback(ovr, slmt, bab);
	    scr -= creepfwrd(ovr, slmt, bab);
	    if (ovr < 0 && wjxt) return (NEVSEL);	// excessive overlap skip this hsp
	    float	dpspace = (float) agap * (float) bgap / MEGA;
	    VTYPE	kscore = NEVSEL;
	    if (algmode.crs < 2 && dpspace >= alprm.maxsp && cmode < 3) {
		bgap = int(alprm.maxsp * MEGA / agap);
		if (cmode == 1) b->left = b->right - bgap;
		else	b->right = b->left + bgap;
		dpspace = (float) agap * (float) bgap / MEGA;
	    }
	    if (dpspace < alprm.maxsp) {
		WINDOW	wdw;
		stripe31(seqs, &wdw, alprm.sh);
		kscore = lspH_ng(wdw);			// DP
	    } else if (cmode == 1) {
		if (wjxt) {				// near 5' end
		    int	bl = wjxt->jy + term;
		    if (bl > b->left) b->left = bl;
		}
		kscore = openendH_ng(cmode);		// drop off
	    } else if (cmode == 2) {
		if (wjxt) {				// near 3' end
		    int	br = bgap - wjxt->jy - term;
		    if (br > b->left && br < b->right) b->right = br;
		}
		kscore = openendH_ng(cmode);		// drop off
	    } else { 		// cmode == 3, space-saving DP
		WINDOW	wdw;
		stripe31(seqs, &wdw, alprm.sh);
		kscore = lspH_ng(wdw);
	    }
	    if (iscore > kscore) {
		if (save_mfd) std::swap(mfd, save_mfd);
		scr += iscore;
	    } else if (kscore > NEVSEL)
		scr += kscore;
	    else	scr = NEVSEL;
	}
	delete save_mfd;
	return (scr);
}

VTYPE Aln2h1::seededH_ng(INT level, int eimode, BOUND& lub)
{
	INT	aexgl = a->inex.exgl;
	INT	aexgr = a->inex.exgr;
	INT	bexgl = b->inex.exgl;
	INT	bexgr = b->inex.exgr;
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
	    jxt = b->jxt;
	    num = b->CdsNo;
	    addsigEjxt(jxt, num, b);
	} else {
	    wl = new Wilip(seqs, pwd, level);
	    wlu = wl->begin();
	    int	nwlu = wl->size();
	    if (nwlu >= 1 && wlu[0].scr > wlu[1].scr) {
		num = wlu->num;
		jxt = wlu->jxt;
		island = true;
const		CHAR*	as = a->at(jxt->jx - 1);
const		CHAR*	bs = b->at(jxt->jy - 2);

		while (avst_equal(*as, *bs)) {
		    --as; bs -= 3; --jxt->jx; jxt->jy -= 3; ++jxt->jlen;
		}
		as += jxt->jlen + 1;
		bs += 3 * (jxt->jlen + 1);
		wjxt = jxt + num - 1;
		while (avst_equal(*as, *bs)) {
		    ++as; bs += 3; ++wjxt->jlen;
		}
	    } else if (nwlu > 1)
		level = algmode.qck - 1;
	}
	save_range(seqs, rng, 2);
	wlprm = setwlprm(level);
	int	backstep = wlprm->width;
	BOUND	bab = lub;
	if (num) {
	    jxt[num].jx = a->right;
	    jxt[num].jy = b->right;
	    a->inex.exgr = 0;
	    b->inex.exgr = 0;
	    agap = jxt->jx - a->left;
	    bgap = jxt->jy - b->left;
	    for (wjxt = jxt; num--; ++wjxt) {
		scr += wjxt->jscr;
		a->right = wjxt->jx;
		b->right = wjxt->jy;
		bab.ua = wjxt->jx + wjxt->jlen;
		if (wjxt[1].jx < bab.ua) bab.ua = wjxt[1].jx;
		bab.ua -= backstep;
		int	lx  = wjxt->jx + wjxt->jlen / 2;
		if (bab.ua < lx) bab.ua = lx;
		bab.ub = wjxt->jy + 3 * (bab.ua - wjxt->jx);
		if (cmode == 2) cmode = 3;
		VTYPE	iscore = interpolateH(level, agap, bgap, cmode, wjxt, bab, backstep);
		if (iscore != NEVSEL) {
		    cmode = 3;
		    scr += iscore;
		    a->left = wjxt->jx + wjxt->jlen;
		    b->left = wjxt->jy + 3 * wjxt->jlen;
		    a->inex.exgl = 0;
		    b->inex.exgl = 0;
		    bab.la = a->right;
		    bab.lb = b->right;
		}
		agap = wjxt[1].jx - a->left;
		bgap = wjxt[1].jy - b->left;
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

	VTYPE	iscore = (bgap < 0)? NEVSEL:
	    interpolateH(level, agap, bgap, cmode, wjxt, bab, backstep);
	if (iscore > NEVSEL) scr += iscore;
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexgl;
	b->inex.exgl = bexgl;
	if (level == lowestlvl && jxt) {
	    wjxt->jx = a->len;
	    wjxt->jy = b->len;
	}
	delete wl;
	return (scr);
}

SKL* Aln2h1::globalH_ng(VTYPE* scr, const WINDOW& wdw)
{
	mfd = new Mfile(sizeof(SKL));
	SKL	wsk;
	mfd->write((UPTR) &wsk);	// dummy call
	if (algmode.qck) {
	    BOUND bab = {a->left, b->left, a->right, b->right};
	    *scr = seededH_ng(lowestlvl, 1, bab);
	} else
	    *scr = lspH_ng(wdw);
	wsk.n = (int) mfd->size();
	SKL*	skl = (SKL*) mfd->flush();
	skl->n = wsk.n - 1;
	skl->m = 1;
	if (skl->n == 0) {
	    delete[] skl;
	    return (0);
	}
	return stdskl3(&skl);
}

VTYPE HomScoreH_ng(const Seq* seqs[], const PwdB* pwd)
{
	Aln2h1 alnv(seqs, pwd);
	WINDOW	wdw;
	stripe31(seqs, &wdw, alprm.sh);
	return alnv.forwardH_ng(0, wdw);
}

SKL* alignH_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr)
{
	Aln2h1 alnv(seqs, pwd);
	WINDOW  wdw;
	stripe31(seqs, &wdw, alprm.sh); // Recacl. window boundaries
	return (alnv.globalH_ng(scr, wdw));
}
