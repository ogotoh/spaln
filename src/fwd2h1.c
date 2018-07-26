/*****************************************************************************
*
*	Alignment of protein vs. nucleotide sequences.
*	5' and 3' splice site signals and coding potential are considered.
*	Assumes there is no internal gap in the reference protein sequence.
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

#define	DEBUG	0
#define	TERMGOP	0

#include "aln.h"
#include "vmf.h"
#include "wln.h"

#define isEIJ(phs) (algmode.lsg && (phs) > -2)

#define	INTR	2
#define	NCANDS	(INTR * NOL)
#define	NQUE	3
#define	MAXR	4

struct RECD {
	VTYPE	val;
	int	dir;
	long	ptr;
	int	jnc;
};

struct LRECD {
	VTYPE	val;
	int	dir;
	long	ptr;
	int	jnc;
	int	lst;
};

struct WRECD {
	VTYPE	val;
	int	dir;
	int	upr;
	int	lwr;
	int	lst;
	int	jnc;
	VTYPE	sig;
};

struct CRECD {
	VTYPE	val;
	int	dir;
	int	upr;
	int	lwr;
	int	mlb;
	int	nlb;
	int	jnc;
	COLONY* clny;
};

typedef LRECD	EIJUNCP[3][NCANDS+1];
typedef int	EIJODRP[3][NCANDS+1];

struct BOUND {int la, lb, ua, ub;};

// alignment parameters independent of each sequence
// the base class is described in fwd2b1.C

class Aln2h1 {
protected:
	Seq**	seqs;
	Seq*&	a;
	Seq*&	b;
	WINDOW*	wdw;
	PwdB*	pwd;
	Mfile*	mfd;
	Vmf*	vmf;
	SpJunc*	spjcs;
	VTYPE	pinthr;
	VTYPE	spjthr;
	INT	lowestlvl;
public:
	Aln2h1(Seq** _seqs, PwdB* _pwd, WINDOW* _wdw);
	~Aln2h1() {delete mfd; delete spjcs;}
	void	initH_ng(RECD* hf[]);
	RECD*	lastH_ng(RECD* hhb[]);
	VTYPE	forwardH_ng(long pp[]);
	void	finitH_ng(WRECD* hhf[], int mm);
	void	binitH_ng(WRECD* hhf[], int mm);
	VTYPE	centerH_ng(int* ml, int* mr, int* nl, int* nr, 
		WINDOW* wdwf, WINDOW* wdwb);
	void	cinitH_ng(CRECD* hhc[]);
	void	pfinitH_ng(LRECD* hhg[]);
	void	pbinitH_ng(LRECD* hhg[]);
	VTYPE	pincersH_ng(long* ptr, int cmode);
	bool	pincerTrcbkH_ng(int cmode, VTYPE& iscore);
	bool	indelfreespjH(int agap, VTYPE& scr);
	Colonies*	fwdswgH_ng(VTYPE* scr);
	VTYPE	diagonalH_ng();
	VTYPE	backforth(int n, BOUND& lub);
	VTYPE	trcbkalignH_ng();
	VTYPE	lspH_ng();
	VTYPE	cds5end(int x, int y);
	VTYPE	cds3end(int x, int y);
	VTYPE	seededH_ng(INT level, int cmode, BOUND& lub);
	SKL*	globalH_ng(VTYPE* scr);
	VTYPE	creepback(int ovr, VTYPE bscr, BOUND& lub);
	VTYPE	creepfwrd(int& ovr, VTYPE bscr, BOUND& lub);
};

static	RECD	oom = {NEVSEL, 0, 0L, 0};
static	LRECD	ooml = {NEVSEL, 0, 0L, 0, 0};
static	WRECD	oomW = {NEVSEL, 0, INT_MIN, INT_MAX, 0, 0};
static	CRECD	oomC = {NEVSEL, 0, INT_MIN, INT_MAX, 0, 0, 0, 0};

static	VTYPE	SumCodePot(EXIN* bb, int i, CHAR* cs, PwdB* pwd);
static	void	addsigEjxt(JUXT* jxt, int num, Seq* b);

// initialization of alignment process

Aln2h1::Aln2h1(Seq* _seqs[], PwdB* _pwd, WINDOW* _wdw) :
	seqs(_seqs), a(seqs[0]), b(seqs[1]), wdw(_wdw), pwd(_pwd), 
	mfd(0), vmf(0), lowestlvl(b->wllvl)
{
	pinthr = pwd->GapPenalty(1);
	spjthr = pinthr + pwd->IntPen->Penalty();
	spjcs = new SpJunc(b);
}

void Aln2h1::initH_ng(RECD* hh[])
{
	int	n = b->left;
	int	r = b->left - 3 * a->left;
	int	rr = b->right - 3 * a->left;
	int	dir = a->inex.exgl? DEAD: DIAG;
	RECD*	hf[NOL + 1];
	EXIN*	bb = b->exin->score(n);

	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hh[k] + r;
	if (wdw->up < rr) rr = wdw->up;
	for (int i = 0; r <= wdw->up + 3; ++r, ++bb, ++n, ++i) {
	    RECD* h = hf[0]++;
 	    for (int k = 1; k < pwd->Nrow; ++k) *hf[k]++ = oom;
	    if (i == 0 || (a->inex.exgl && i < 3)) {
 		h->val = 0;
		h->dir = dir;
		h->ptr = vmf? vmf->add(a->left, n, 0): 0;
		h->jnc = n;
		if (!a->inex.exgl) continue;
	    } else if (!a->inex.exgl || r > rr) {       // global
		*h = oom;
		continue;
	    } else {
		*h = h[-3];
		int	k = n - h->jnc;
#if TERMGOP
		if (k == 3) h->val += pwd->BasicGOP;
#else
		if (k == 3 && !a->inex.exgl) h->val += pwd->BasicGOP;
#endif
		h->val += pwd->GapExtPen3(k) + bb[-2].sigE;
		h->dir = HORI;
		VTYPE	x = h[-1].val + pwd->GapW1;
		if (x > h->val) {*h = h[-1]; h->val = x; h->dir = HOR1;}
		x = h[-2].val + pwd->GapW2;
		if (x > h->val) {*h = h[-2]; h->val = x; h->dir = HOR2;}
	    }
	    VTYPE	x = 0;
	    if ((algmode.lcl & 1) && bb[1].sigS > x) x = bb[1].sigS;
	    if ((algmode.lcl & 4) && bb->sig3 > x) x = bb->sig3;
	    if (h->val < x) {
		h->val = x;
		h->dir = dir;
		h->ptr = vmf? vmf->add(a->left, n, 0): 0; 
		h->jnc = n;
	    }
	}

	r = b->left - 3 * a->left;
	rr = b->left - 3 * a->right;
	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hh[k] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	for (int i = 1; --r >= wdw->lw - 3; ++i) {
	    RECD*	h = --hf[0];
	    RECD*	g = --hf[1];
	    for (int k = 2; k < pwd->Nrow; ++k) *--hf[k] = oom;
	    if (r < rr) *h = *g = oom;
	    else if (b->inex.exgl) {
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->left + i % 3;
		h->ptr = vmf? vmf->add(a->left + i / 3, h->jnc, 0): 0;
		*g = oom;
	    } else if (i <= 3) {
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

RECD* Aln2h1::lastH_ng(RECD* hh[])
{
	int	glen[3];
	int	rw = wdw->lw;
	int	m3 = 3 * a->right;
	int	rf = b->left - m3;
	if (rf > rw) rw = rf;
	else	rf = rw;
	RECD*	h = hh[0] + rw;
	RECD*	h9 = hh[0] + b->right - m3;
	RECD*	f = h9;
	VTYPE	mx = h9->val;
	EXIN*	bb = b->exin->score(rw + m3);

	for (int p = 0; p < 3; ++p) glen[p] = 0;
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
			h->val = x;
			h->dir = HORI;
			h->ptr = h[-3].ptr;
			h->jnc = h[-3].jnc;
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
	    rw = wdw->up;
	    rf = b->right - 3 * a->left;
	    if (rf < rw) rw = rf;
	    h = hh[0] + rw;
	    for ( ; h > h9; --h, --rw) {
		VTYPE	x = h->val + (rw % 3? pwd->ExtraGOP: 0);
		if (f->val < x) {f = h; f->val = x;}
	    }
	}
	int p = f - h9;
	rf = a->right;		// m9
	rw = b->right;		// n9
	if (p > 0) {
	    rf -= (p + 2) / 3;
	    if (p %= 3) rw -= (3 - p); 
	} else if (p < 0) rw += p;
	if (vmf && a->inex.exgr) f->ptr = vmf->add(rf, rw, f->ptr);
	return (f);
}

VTYPE Aln2h1::forwardH_ng(long* pp)
{
	RECD*	hh[NOL];
	RECD*	hf[NOL];	// [DIAG, HORI, HORL]
	RECD	e1[NQUE];
	RECD	e2[NQUE];
#if INTR > 1
	RECD	hl[NOL][3][INTR+1]; // [DIA, HORI, HORL][phase][candidates]
	int	nl[NOL][3][INTR+1]; // [DIA, HORI, HORL][phase][candidates]
#else
	RECD    hl[NOL][3];
#endif
	VSKLP   maxh = {NEVSEL, a->left, b->left, 0L};
	int     Local = algmode.lcl & 16;
	int     LocalL = Local && a->inex.exgl && b->inex.exgl;
	int     LocalR = Local && a->inex.exgr && b->inex.exgr;

	RECD*	buf = new RECD[pwd->Nrow * wdw->width];
	hh[0] = buf - wdw->lw + 3;
	for (int k = 1; k < pwd->Nrow; ++k) hh[k] = hh[k-1] + wdw->width;
	if (vmf) vmf->add(0, 0, 0);	// Skip 0-th record
	initH_ng(hh);

	int	m = a->left;
	PfqItr	api(a, m);		// iterator
	int	api_size = api.size();
	if (!a->inex.exgl) --m;		// global
	int	n1 = 3 * m + wdw->lw - 1;
	int	n2 = 3 * m + wdw->up;
	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as) {
	    bool	internal = !a->inex.exgr || m < a->right;
	    n1 += 3;
	    n2 += 3;
	    int		n0 = max(n1, b->left);
	    int		n9 = min(n2, b->right);
	    int		n  = n0;
	    int		r = n - 3 * m;
	    for (int q = 0; q < NQUE; ++q) e1[q] = e2[q] = oom;
	    if (!b->inex.exgl && m == a->left && n == b->left) {
		e1[2] = e2[2] = hh[0][r];
		e1[2].val = pwd->GapW3;
		e2[2].val = pwd->GapW3L;
	    }
	    CHAR*	bs = b->at(n);
	    EXIN*	bb = b->exin->score(n);
	    RECD*	h = hh[0] + r;
	    RECD*	g = hh[1] + r;
	    RECD*	g2 = (pwd->Noll == 3)? hh[2] + r: 0;
	    RECD*	sj = algmode.lsg? hh[pwd->Noll] + r: &oom;
	    RECD	hq[NOL] = {oom};	// previous state
	    for (int k = 0; k < pwd->Noll; ++k) {
	      for (int p = 0; p < 3; ++p) {
#if INTR > 1
		RECD*	phl = hl[k][p];
		for (int l = 0; l <= INTR; ++l, ++phl) {
		  nl[k][p][l] = l;
		  *phl = oom;
		}
#else
		hl[k][p] = oom;
#endif
	      }
	    }
#if DEBUG
	    if (algmode.nsa & 8) {
		printf("%2d %2d %2d", m, n, h->dir);
		putvar(h->val); putchar('\n');
	    }
#endif
	    bool	a_in_zone = api_size && api.eq(m);
	    int	q = 0;		// queue pointer
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x, y;
		++bb; ++g;
		if (algmode.lsg) ++sj;
		hf[0] = ++h;
		RECD*	eq1 = e1 + q;
		RECD*	eq2 = 0;
		hf[1] = eq1;
		if (g2) {
		    g2++;
		    eq2 = e2 + q;
		    hf[2] = eq2;
		}

//	diagonal match
		VTYPE	sigE = bb[-2].sigE;
		RECD*	from = h;
		RECD*	mx = h;
		if (m == a->left) goto Horizon;
		if (n > b->left + 2) {
		    *hq = *h;
		    if (sj->dir) {
			*h = *sj;
			sj->dir = 0;
		    } else     h->val += pwd->sim2(as, bs - 1) + sigE;
		    h->dir = isdiag(from)? DIAG: NEWD;
		}   else	*h = oom;

//	vertical gap extention
		y = g[3].val + pwd->BasicGEP;

//	1 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    g->val = x;
		    g->dir = SLA2;
		    g->jnc = from->jnc;
		    g->ptr = from->ptr;
		} else	g->val = y;

//	2 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val) {
		    g->val = x;
		    g->dir = SLA1;
		    g->jnc = from->jnc;
		    g->ptr = from->ptr;
		}

//	normal deletion
		x = (++from)->val + pwd->GapW3;
		if (x >= g->val) {
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		    g->jnc = from->jnc;
		    g->ptr = (vmf && from->dir & SPIN)?
			vmf->add(m - 1, n, from->ptr): from->ptr;
		} else if (y >= g->val) {
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		    g->jnc = g[3].jnc;
		    g->ptr = g[3].ptr;
		}
		if (g->val > mx->val) mx = g;

//	long deletion
		if (g2) {
		  x = from->val + pwd->GapW3L;
		  y = g2[3].val + pwd->LongGEP;
		  if (x >= y) {
		    g2->val = x;
		    g2->dir = VERL | (g2->dir & SPJC);
		    g2->jnc = from->jnc;
		    g2->ptr = (vmf && from->dir & SPIN)? 
			vmf->add(m - 1, n, from->ptr): from->ptr;
		  } else {
		    *g2 = g2[3];
		    g2->val = y;
		  }
		  if (g2->val > mx->val) mx = g2;
		}
Horizon:
//	nomal insertion
		hq[1] = *eq1;
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
		    if (eq2) {
			hq[2] = *eq2;
			x = from->val + pwd->GapW3L;
			y = eq2->val += pwd->LongGEP;
			if (x > y) {
			    *eq2 = *from;
			    eq2->val = x;
			}
			if (!(eq2->dir & SPF2)) eq2->val += sigE;
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
		if (isEIJ(bb->phs3) && internal) {
		    int phs = (bb->phs3 == 2)? 1: bb->phs3;
Acceptor:
		    if (phs == -1)
			sj->val =  mx->val + pwd->sim2(as + 1, bs + 2);
		    int nb = n - phs;
		    VTYPE sigJ = a_in_zone? api.match_score(3 * m - phs): 0;
		    RECD*	maxdhl = 0;
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = hf[d];
#if INTR > 1
			RECD*	maxphl = 0;
			int*	pnl = nl[d][phs + 1];
			for (int l = 0; l < INTR; ++l) {
			    RECD*	phl = hl[d][phs + 1] + pnl[l];
			    if (!phl->dir) break;
			    y = phl->val + sigJ +
				pwd->IntPen->Penalty(nb - phl->jnc) +
				b->exin->sig53(phl->jnc, nb, IE53);
			    if (phs == 1) {
				CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
				y += pwd->pmt->prematT(cs);
				if (isdiag(phl)) y += pwd->sim2(as, cs);
			    } else if (phs == -1) {
				CHAR*	cs = spjcs->spjseq(phl->jnc, nb) + 1;
				y += pwd->pmt->prematT(cs);
				x = y + pwd->sim2(as + 1, cs);
				if (x > sj->val) {
				    maxdhl = phl;
				    sj->val = x;
				}
			    }
			    if (y > from->val) {
				from->val = y;
				maxphl = phl;
		    	    }
			}
			if (maxphl) {
		 	    if (vmf)
				from->ptr = vmf->add(m, n, 
				vmf->add(m, maxphl->jnc + phs, maxphl->ptr));
		 	    from->jnc = nb;
			    from->dir = maxphl->dir | SPJCI;
			    if (phs == -1) from->dir |= SPF2;
			    if (from->val > mx->val) mx = from;
			}
#else
			RECD*	phl = hl[d] + phs + 1;
			if (!phl->dir) break;
			y = phl->val + sigJ + 
			    pwd->IntPen->Penalty(nb - phl->jnc) + 
			    b->exin->sig53(phl->jnc, nb, IE53);
			if (phs == 1) {
			    CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
			    y += pwd->pmt->prematT(cs);
			    if (isdiag(phl)) y += pwd->sim2(as, cs);
			} else if (phs == -1) {
			    CHAR*	cs = spjcs->spjseq(phl->jnc, nb) + 1;
			    y += pwt->pmt->prematT(cs);
			    VTYPE	x = y + pwd->sim2(as + 1, cs);
			    if (x > sj->val) {
				maxdhl = phl;
				sj->val = x;
			    }
			}
			if (y > from->val) {
			    from->val = y;
		 	    if (vmf)
				from->ptr = vmf->add(m, n, 
				vmf->add(m, phl->jnc + phs, phl->ptr));
		 	    from->jnc = nb;
			    from->dir = phl->dir | SPJCI;
			    if (phs == -1) from->dir |= SPF2;
			    if (from->val > mx->val) mx = from;
			}
#endif
		    }
		    if (maxdhl) {
			sj->jnc = nb;
			sj->dir = maxdhl->dir | SPALL;
		 	if (vmf)
			    sj->ptr = vmf->add(m, n, 
				vmf->add(m, maxdhl->jnc + phs, maxdhl->ptr));
		    }
		    if (bb->phs3 - phs == 1) {	// AGAG
			phs = -1;
			goto Acceptor;
		    }
		}

//	Find optimal path
		y = h->val;
		if (h != mx) *h = *mx;	// non-diagonal
		else if (Local && y > hq->val) {
		    if (LocalL && hq->dir == 0 && !(h->dir & SPJC))
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
		if (isEIJ(bb->phs5) && internal) {
		    int	phs = (bb->phs5 == 2)? 1: bb->phs5;
Donor:
		    int	nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(nb, 0, IE5);
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = (phs == 1)? hq + d: hf[d];
			// An orphan exon is disallowed
			// vertical path is not transfered
			// not rigorous when phase != 0
			if (!from->dir || isvert(from) || (from->dir & SPIN))
			    continue;
			if (phs == 1 && nb - from->jnc == 2 && bb[-2].sigS > 0)
			    from->val -= bb[-2].sigS;
			x = from->val + sigJ;
#if INTR > 1
			RECD*	phl = hl[d][phs + 1];
			int*	pnl = nl[d][phs + 1];
			int	l = INTR;
			while (--l >= 0) {
			    if (x > phl[pnl[l]].val)
				gswap(pnl[l], pnl[l + 1]);
			    else
				break;
			}
			if (++l < INTR) {
			    phl[pnl[INTR]].dir = 0;
			    phl += pnl[l];
			    phl->val = x;
			    phl->jnc = nb;
			    phl->dir = from->dir;
			    phl->ptr = from->ptr;
			}
#else
			phl = hl[d] + phs + 1;
			if (x > phl->val) {
			    phl->val = x;
			    phl->jnc = nb;
			    phl->dir = from->dir;
			    phl->ptr = from->ptr;
			}
#endif
		    }
		    if (bb->phs5 - phs == 1) {	//	GTGT..
			phs = -1;
			goto Donor;
		    }
		}

#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(eq1->val);
		if (g2) {
		    putvar(g2->val); putvar(eq2->val);
		}
		if (algmode.lsg) {
#if INTR > 1
		    putvar(hl[0][2][nl[0][2][0]].val); 
		    putvar(hl[0][1][nl[0][1][0]].val);
		    putvar(hl[0][0][nl[0][0][0]].val);
#else
		    putvar(hl[0][2].val); 
		    putvar(hl[0][1].val);
		    putvar(hl[0][0].val);
#endif
		    printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
		    printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		}
		putchar('\n');
	}
#endif

	    }	// end of n-loop
	    if (a_in_zone) ++api;	// has exon-exon junctions
	}	// end of m-loop

	RECD*	mx = lastH_ng(hh);
	if (maxh.val > mx->val) {
	    if (pp) *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    maxh.val = mx->val;
	    if (pp) *pp = mx->ptr;
	}
	delete[] (hh[0] + wdw->lw - 3);
	return (maxh.val);
}

static VTYPE SumCodePot(EXIN* bb, int i, CHAR* cs, PwdB* pwd)
{
	VTYPE	y = 0;

	for (++bb; i > 0; i -= 3, bb += 3) {
	    VTYPE	hvl;
	    if (cs) {
		hvl = pwd->pmt->prematT(cs);
		cs = 0;
	    } else
		hvl = bb->sigE;
	    y += hvl;
	}
	return (y);
}

VTYPE skl_rngH_ng(Seq* seqs[], Gsinfo* gsi, PwdB* pwd)
{
static	const char*	fswarn = "%s %s FrameShift: %d %d\n";
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	SKL*	wsk = gsi->skl;
	int	num = (wsk++)->n;
	SpJunc*	spjcs = new SpJunc(b);
	VTYPE	h = 0;		// 
	VTYPE	hi = NEVSEL;		// intron-containing path
	VTYPE	ha = 0;
	VTYPE	hb = 0;
	VTYPE	hvl = 0;
	bool	ivl = false;
	VTYPE	sig5 = 0;
	VTYPE	sig3 = 0;
	int	insert = 0;
	int	deletn = 0;
	int	intlen = 0;
	int	preint = 0;
	int	phs = 0;
	int	nb = 0;
	EISCR	rbuf;
	FSTAT*	fst = &gsi->fstat;
	FSTAT	pst;
	Eijnc*	eijnc = gsi->eijnc = new Eijnc();
	Cigar*	cigar = gsi->cigar = 0;
	Vulgar*	vlgar = gsi->vlgar = 0;
	switch (algmode.nsa) {
	    case EXN_FORM: 
	    case CIG_FORM: cigar = gsi->cigar = new Cigar(); break;
	    case VLG_FORM: vlgar = gsi->vlgar = new Vulgar(); break;
	    default: break;
	}
	vclear(fst);
	vclear(&pst);
	vclear(&rbuf);
	gsi->noeij = 0;
	CHAR*	cs = b->at(wsk[num-1].n - 2);
	int	termcodon = (*cs == TRM || *cs == TRM2);
	cs = 0;
	if (termcodon) wsk[num-1].n -= 3;
	if (wsk[1].n == wsk->n && b->inex.exgl) {++wsk; --num;}
	int	m = wsk->m;
	int	n = wsk->n;
	CHAR*	as = a->at(m);
	CHAR*	bs = b->at(n);
	EXIN*	bb = b->exin->score(n);
	EXIN*	bb_last = b->exin->score(b->right);
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
	while (--num) {
	    ++wsk;
	    bool	term = num == 1;		// tail gap?
	    int	mi = (wsk->m - m) * 3;
	    if (insert && (mi || term)) {
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
		if (insert) {				// post-intron gap
		    if (term && IsTerm(bs[insert - 1])) insert -= 3;
		    if (cigar) cigar->push('D', insert);
		    phs = insert % 3;
		    insert -= phs;
		    fst->gap += 1;
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
		    }
		    if (vlgar && insert) vlgar->push('G', 0, insert);
		}
		insert = intlen = preint = 0; hi = NEVSEL;
	    }
	    if (deletn) {
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
		    ++as;
		    deletn -= phs;
		    phs = 3 - phs;
		    bs += phs;
		    bb += phs;
		}
		if (vlgar && deletn > 2) vlgar->push('G', deletn / 3, 0);
		deletn = 0;
	    }
	    int	ni = wsk->n - n;
	    int	i = mi - ni;
	    int	d = (i >= 0)? ni: mi;
	    if (d) {
		if (cigar) cigar->push('M', d);
		if (vlgar) vlgar->push('M', d / 3, d);
		n += d;
		m += d / 3;
		for ( ; d > 2; d -= 3, ++as, bs += 3, bb += 3) {
		   CHAR*	gs = (cs? cs: bs) + 1;
		   hvl = cs?
			pwd->sim2(as, gs) + pwd->pmt->prematT(gs):
			pwd->sim2(as, gs) + bb[1].sigE;
		    h += hvl;
		    ivl = (*as == *gs || (*as == SER && *gs == SER2));
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
		EXIN*	b3 = bb + (i = -i);
		if (hi > NEVSEL) {			// after intron
		    hi += SumCodePot(bb, i, cs? cs + 1: 0, pwd);
		} else if (i >= IntronPrm.llmt && wsk->n < b->right) {	// intron?
		    int	n3 = n + i;
		    int	phs5 = (bb->phs5 == 2)? b3->phs3: bb->phs5;
		    int	phs3  = (b3->phs3 == 2)? bb->phs5: b3->phs3;
		    VTYPE	xm = NEVSEL;
		    VTYPE	xi = NEVSEL;
		 	// potential intron
		    if (phs3 == 2 && phs5 == 2) {	// GTGT....AGAG
			// phs3 = phs5 = -1; 
			nb = n + 1; n3 = nb + i;
			if (isJunct(phs3, phs5)) {
			    xm = b->exin->sig53(nb, n3, IE5P3);
			    if (usespb) xm += api.match_score(3 * m + 1);
			}
			phs3 = phs5 = 1;
		    }
		    nb = n - phs3; n3 = nb + i;
		    if (isJunct(phs3, phs5)) {
			sig5 = b->exin->sig53(nb, 0, IE5);
			sig3 = b->exin->sig53(nb, n3, IE53);
			xi = b->exin->sig53(nb, n3, IE5P3);
			if (usespb) xi += api.match_score(3 * m - phs3);
			preint = insert;
			if (phs3 != 0) cs = spjcs->spjseq(nb, n3);
			if (insert == 0 && phs3 == 1) {
			    xi += (pwd->sim2)(as - 1, cs) - hvl 
				+ pwd->pmt->prematT(cs);
			    bool match = (as[-1] == *cs || 
				(as[-1] == SER && *cs == SER2));
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
			sig5 = b->exin->sig53(nb, 0, IE5);
			sig3 = b->exin->sig53(nb, n3, IE53);
			cs = spjcs->spjseq(nb, n3);
		    }
		    if (xi > NEVSEL) {
			xi += pwd->IntPen->Penalty(i);
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
		}
		h += SumCodePot(bb, i, 0, pwd);
		bb = b3;
		bs += i;
		insert += i;
	    }
	    m = wsk->m;
	    n = wsk->n;
	    if (usespb) while (api.lt(m)) ++api;
	}
	if (bb + 3 < bb_last) {
	    if (algmode.lcl & (16 + 2) && bb[3].sigT > 0)
		sig5 = bb[3].sigT;
	    if (algmode.lcl & (16 + 8) && bb[2].sig5 > 0 && bb[2].sig5 > bb[3].sigT)
		sig5 = bb[2].sig5;
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
	fst->val = (FTYPE) h;
	if (termcodon) wsk->n += 3;
	return (h);
}

void Aln2h1::finitH_ng(WRECD* hhf[], int mm)
{
	int	n = b->left;
	int	r = b->left - 3 * a->left;
	int	r0 = r;
	int	rr = b->right - 3 * a->left;
	int	dir = a->inex.exgl? DEAD: DIAG;
	WRECD*	hf[NOL + 1];
	EXIN*	bb = b->exin->score(n);

	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hhf[k] + r;
	if (wdw->up < rr) rr = wdw->up;
	for (int i = 0; r <= wdw->up + 3; ++r, ++bb, ++n, ++i) {
	    WRECD*	h = hf[0]++;
 	    for (int d = 1; d < pwd->Nrow; ++d) *hf[d]++ = oomW;
	    if (i == 0 || (a->inex.exgl && i < 3)) {
 		h->val = 0;
		h->dir = dir;
		h->jnc = n;
		h->lwr = r0;
	    } else if (!a->inex.exgl || r > rr) {       // global
		*h = oomW;
		continue;
	    } else {
		*h = h[-3];
		int	d = n - h->jnc;
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
	    h->upr = h->lst = r;
	    if (!a->inex.exgl) continue;
	    VTYPE	x = 0;
	    if ((algmode.lcl & 1) && bb[1].sigS > x) x = bb[1].sigS;
	    if ((algmode.lcl & 4) && bb->sig3 > x) x = bb->sig3;
	    if (h->val < x) {
		h->val = x;
		h->dir = dir;
		h->jnc = n;
		h->lwr = r;
	    }
	}

	r = r0;
	rr = b->left - 3 * mm;
	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hhf[k] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	for (int i = 1; --r >= wdw->lw - 3; ++i) {
	    WRECD*	h = --hf[0];
	    WRECD*	g = --hf[1];
	    for (int d = 2; d < pwd->Nrow; ++d) *--hf[d] = oomW;
	    if (r < rr) *h = *g = oomW;
	    else if (b->inex.exgl) {
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->left + i % 3;
		h->upr = h->lwr = h->lst = r;
		*g = oomW;
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

void Aln2h1::binitH_ng(WRECD* hhb[], int mm)
{
	int	n = b->right;
	int	r = b->right - 3 * a->right;
	int	r9 = r;
	int	rr = b->left - 3 * a->right;
	int	dir = a->inex.exgr? DEAD: DIAG;
	WRECD*	hb[NOL + 1];
	EXIN*	bb = b->exin->score(b->right);

	for (int k = 0; k < pwd->Nrow; ++k) hb[k] = hhb[k] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	for (int i = 0; r >= wdw->lw - 3; --r, --n, --bb, ++i) {
	    WRECD*	h = hb[0]--;
 	    for (int d = 1; d < pwd->Nrow; ++d) *hb[d]-- = oomW;
	    if (i == 0 || (a->inex.exgr && i < 3)) {
 		h->val = 0;
		h->dir = dir;
		h->jnc = n;
		h->upr = r;
	    } else if (!a->inex.exgr || r < rr) {       // global
		*h = oomW;
		continue;
 	    } else {
		*h = h[3];
		int	k = h->jnc - n;
#if TERMGOP
		if (k == 3) h->val += pwd->BasicGOP;
#else
		if (k == 3 && !a->inex.exgr) h->val += pwd->BasicGOP;
#endif
		h->val += pwd->GapExtPen3(k);
	    	h->val += bb[1].sigE;
		h->dir = HORI;
		VTYPE	x = h[1].val + pwd->GapW1;
		if (x > h->val) {*h = h[1]; h->val = x; h->dir = HOR1;}
		x = h[2].val + pwd->GapW2;
		if (x > h->val) {*h = h[2]; h->val = x; h->dir = HOR2;}
	    }
	    h->lwr = h->lst = r;
	    if (!a->inex.exgr) continue;
	    VTYPE	x = 0;
	    if ((algmode.lcl & 2) && bb[1].sigT > x)	x = bb[1].sigT;
	    if ((algmode.lcl & 8) && bb->sig5 > x) 	x = bb->sig5;
	    if (h->val < x) {
		h->val = x;
		h->dir = dir;
		h->jnc = n;
		h->upr = r;
		if (x == bb[1].sigT) h->upr += 3;
	    }
	}

	r = r9;
	rr = b->right - 3 * mm;
	for (int d = 0; d < pwd->Nrow; ++d) hb[d] = hhb[d] + r;
	if (wdw->up < rr) rr = wdw->up;
	for (int i = 1; ++r <= wdw->up + 3; ++i) {
	    WRECD*	h = ++hb[0];
	    WRECD*	g = ++hb[1];
	    for (int d = 2; d < pwd->Nrow; ++d) *++hb[d] = oomW;
	    if (r > rr) *h = *g = oomW;
	    else if (b->inex.exgr) {
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->right - i % 3;
		h->lwr = h->upr = h->lst = r;
		*g = oomW;
	    } else {
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
		h->upr = g->upr = r;
	    }
	}
}

VTYPE Aln2h1::centerH_ng(int* ml, int* mr, int* nl, int* nr, 
	WINDOW* wdwf, WINDOW* wdwb)
{
	int	kk = 0;
	int	jtigh = 0;
	int	kb = 0;
	int	jb = 0;
	int	rr[MAXR+1];
	int	pp[MAXR+1];
	int	phs3 = 0;
#if INTR > 1
	WRECD	hl[NOL][3][INTR+1];
	int	nx[NOL][3][INTR+1]; // [DIA, HORI, HORL][phase][candidates]
#else
	WRECD   hl[NOL][3];
#endif
	WRECD*	hhf[NOL+1];
	WRECD*	hhb[2 * NOL];
	WRECD*	hf[NOL+1];
	WRECD*	hb[NOL+1];
	WRECD*	f1 = 0;
	WRECD*	f2 = 0;
	WRECD	e1[NQUE];
	WRECD	e2[NQUE];
	VTYPE	mxh = NEVSEL;
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

	hhf[0] = new WRECD[pwd->Nrow * wdw->width] - wdw->lw + 3; // +1: sp junc
	hhb[0] = new WRECD[pwd->Nrwb * wdw->width] - wdw->lw + 3; // +1: sp junc
	for (int k = 1; k < pwd->Nrow; ++k) hhf[k] = hhf[k-1] + wdw->width;
	for (int k = 1; k < pwd->Nrwb; ++k) hhb[k] = hhb[k-1] + wdw->width;
	int	mm = (a->left + a->right + 1) / 2;

/*******	backward phase	*******/

	binitH_ng(hhb, mm);
	int	m = a->right;
	PfqItr	api(a, m - 1);
	int	api_size = api.size();
	if (!a->inex.exgr) ++m; // global
	CHAR*	as = a->at(m);
	int	n1 = 3 * m + wdw->lw;
	int	n2 = 3 * m + wdw->up + 1;
	for ( ; --m >= mm; ) {
	    --as;
	    n1 -= 3; n2 -= 3;
	    int		n0 = min(n2, b->right);
	    int		n9 = max(n1, b->left);
	    int		n  = n0;
	    CHAR*	bs = b->at(n);
	    EXIN*	bb = b->exin->score(n);
	    int		r = n - 3 * m;
	    for (int q = 0; q < NQUE; ++q) e1[q] = e2[q] = oomW;
	    if (!b->inex.exgr && n == b->right && m == a->right) {
		e1[2] = e2[2] = hhb[0][r];
		e1[2].val = pwd->GapW3;
		e2[2].val = pwd->GapW3L;
	    }
	    int		d = 0;
	    WRECD*	h = hhb[d] + r;
	    WRECD*	g = hhb[++d] + r;
	    WRECD*	g2 = (pwd->Noll == 3)? hhb[++d] + r: 0;
	    WRECD*	sj = algmode.lsg? hhb[++d] + r: &oomW;
	    WRECD	hq[NOL] = {oomW};	// previous state
	    if (m == mm) {
		*(f1 = hhb[++d] + r) = oomW;
		if (pwd->Noll == 3) *(f2 = hhb[++d] + r) = oomW;
	    }
	    for (int k = 0; k < pwd->Noll; ++k) {
	      for (int p = 0; p < 3; ++p) {
#if INTR > 1
		WRECD*	phl = hl[k][p];
		for (int l = 0; l <= INTR; ++l, ++phl) {
		  nx[k][p][l] = l;
		  *phl = oomW;
		}
#else
		hl[k][p] = oomW;
#endif
	      }
	    }
#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    bool	a_in_zone = api_size && api.eq(m);
	    int	q = 0;
	    for ( ; --n >= n9; ) {
		VTYPE	x, y;
		--bs; --g; --r;
		WRECD*	eq1 = e1 + q;
		WRECD*	eq2 = 0;
		hb[0] = --h;
		hb[1] = eq1;
		WRECD*	from = h;
		WRECD*	mx = h;
		if (algmode.lsg) --sj;
		if (g2) {
		    --g2;
		    hb[2] = eq2 = e2 + q;
		}
		VTYPE	sigE = (bb--)->sigE;

//	diagonal match
		if (m == a->right) goto HorizonB;
		if (n < b->right - 2) {
		    *hq = *h;
		    if (sj->dir) {
			*h = *sj;
			sj->dir = 0;
		    } else	h->val += pwd->sim2(as, bs + 1) + sigE;
		    h->dir = isdiag(from)? DIAG: NEWD;
		} else	*h = oomW;

//	vertical gap extention
		y = g[-3].val + pwd->BasicGEP;

//	1 nt deletion
		--from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA2;
		    g->lst = r - 1;
		} else	g->val = y;

//	2 nt deletion
		--from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val)	{
		    *g = *from;
		    g->val = x;
		    g->dir = SLA1;
		    g->lst = r - 2;
		}
//	normal deletion
		x = (--from)->val + pwd->GapW3;
		if (x >= g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		    g->lst = r - 3;
		} else if (y >= g->val) {
		    *g = g[-3];
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		}
		if (g->val >= mx->val) mx = g;

//	long deletion
		if (g2) {
		  x = from->val + pwd->GapW3L;
		  y = g2[-3].val + pwd->LongGEP;
		  if (x >= y) {
		    *g2 = *from;
		    g2->val = x;
		    g2->dir = VERL | (g2->dir & SPJC);
		    g2->lst = r - 3;
		  } else {
		    *g2 = g2[-3];
		    g2->val = y;
		  }
		  if (g2->val >= mx->val) mx = g2;
		}
HorizonB:
//	nomal insertion
		hq[1] = *eq1;
		if (n < n0 - 2) {
		    from = h + 3;
		    x = from->val + pwd->GapW3;
		    y = eq1->val += pwd->BasicGEP;
		    if (x > y) {
			*eq1 = *from;
			eq1->val = x;
			eq1->lst = r + 3;
		    }
		    if (!(eq1->dir & SPF2))	eq1->val += sigE;
		    eq1->dir = (eq1->dir & SPIN) + HORI;
//	long insertion
		    if (eq2) {
			hq[2] = *eq2;
			x = from->val + pwd->GapW3L;
			y = eq2->val += pwd->LongGEP;
			if (x > y) {
			    *eq2 = *from;
			    eq2->val = x;
			    eq2->lst = r + 3;
			    eq2->dir = (from->dir & SPALL) + HORL;
			}
			if (!(eq2->dir & SPF2)) eq2->val += sigE;
			eq2->dir = (eq2->dir & SPIN) + HORL;
			if (eq2->val > mx->val) mx = e2 + q;
		    }
		}
//	2 nt insertion
		if (n < n0 - 1) {
		    from = h + 2;
		    x = from->val + pwd->GapW2;
		    if (x > eq1->val) {
			*eq1 = *from;
			eq1->val = x;
			eq1->lst = r + 2;
			eq1->dir = (eq1->dir & SPIN) + HOR2;
		    }
		}
//	1 nt insertion
		from = h + 1;
		x = from->val + pwd->GapW1;
		if (x > eq1->val) {
		    *eq1 = *from;
		    eq1->val = x;
		    eq1->lst = r + 1;
		    eq1->dir = (from->dir & SPIN) + HOR1;
		}
		if (eq1->val >= mx->val) mx = e1 + q;
		if (++q == NQUE) q = 0;

//	intron 5' boundary reverse phase
		if (isEIJ(bb->phs5)) {
		    int	phs = (bb->phs5 == 2)? 1: bb->phs5;
		    if (phs == 1)
			sj->val = mx->val + pwd->sim2(as - 1, bs - 2);
DonRev:
		    int	nb = n - phs;
		    VTYPE	sigJ = a_in_zone? api.match_score(3 * m - phs): 0;
		    WRECD*	maxdhl = 0;
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = hb[d];
#if INTR > 1
			WRECD*	maxphl = 0;
			int*	pnx = nx[d][phs + 1];
			for (int l = 0; l < INTR; ++l) {
			    WRECD*	phl = hl[d][phs + 1] + pnx[l];
			    if (!phl->dir) break;
			    y = phl->val + sigJ + 
				b->exin->sig53(nb, phl->jnc, IE53) +
				(b->exin->sig53(nb, 0, IE5) - phl->sig) +
				pwd->IntPen->Penalty(phl->jnc - nb);
			    if (phs == -1) {
				CHAR*	cs = spjcs->spjseq(nb, phl->jnc) + 1;
				y += pwd->pmt->prematT(cs);
				if (isdiag(phl)) y += pwd->sim2(as, cs);
			    } else if (phs == 1) {
				CHAR*	cs = spjcs->spjseq(nb, phl->jnc);
				y += pwd->pmt->prematT(cs);
				x = y + pwd->sim2(as - 1, cs);
				if (x > sj->val) {
				    maxdhl = phl;
				    sj->val = x;
				}
			    }
			    if (y > from->val) {
				from->val = y;
				maxphl = phl;
		    	    }
			}
			if (maxphl) {
	 		    from->jnc = maxphl->jnc;
			    from->dir = maxphl->dir | SPJCI;
			    if (phs == 1) from->dir |= SPF2;
			    from->upr = max(maxphl->upr, r);
			    from->lwr = min(maxphl->lwr, r);
			    from->lst = maxphl->lst + nb - maxphl->jnc;
			    if (from->val > mx->val) mx = from;
			}
#else
			WRECD*	phl = hl[d] + phs + 1;
			if (!phl->dir) break;
			y = phl->val + sigJ +
			    b->exin->sig53(nb, phl->jnc, IE53) + 
			    (b->exin->sig53(nb, 0, IE5) - phl->sig) +
			    pwd->IntPen->Penalty(phl->jnc - nb);
			if (phs == -1) {
			    CHAR*	cs = spjcs->spjseq(nb, phl->jnc) + 1;
			    y += pwd->pmt->prematT(cs);
			    if (isdiag(phl)) y += pwd->sim2(as, cs);
			} else if (phs == 1) {
			    CHAR*	cs = spjcs->spjseq(nb, phl->jnc);
			    y += pwd->pmt->prematT(cs);
			    x = y + pwd->sim2(as - 1, cs);
			    if (x > sj->val) {
				maxdhl = phl;
				sj->val = x;
			    }
			}
			if (y > from->val) {
			    from->val = y;
	 		    from->jnc = phl->jnc;
			    from->dir = phl->dir | SPJCI;
			    if (phs == 1) from->dir |= SPF2;
			    from->upr = max(phl->upr, r);
			    from->lwr = min(phl->lwr, r);
			    from->lst = phl->lst + nb - phl->jnc;
			    if (from->val > mx->val) mx = from;
			}
#endif
		    }
		    if (maxdhl) {
			sj->jnc = maxdhl->jnc;
			sj->dir = maxdhl->dir | SPALL;
			sj->upr = max(maxdhl->upr, r);
			sj->lwr = min(maxdhl->lwr, r);
			sj->lst = maxdhl->lst + nb - maxdhl->jnc;
		    }
		    if (bb->phs5 - phs == 1) {
			phs = -1;
			goto DonRev;
		    }
		}

//	Find optimal path
		y = h->val;
		if (h != mx) {	// non-diagonal
		    if (h->val == mx->val) {
			if (mx->upr < h->upr) mx->upr = h->upr;
			if (mx->lwr > h->lwr) mx->lwr = h->lwr;
		    }
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		}
		if (LocalR && h->val <= 0) h->val = h->dir = 0;

//	intron 3' boundary reverse phase

		if (isEIJ(bb->phs3)) {
		    int	phs = (bb->phs3 == 2)? 1: bb->phs3;
AccRev:
		    int	nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(0, nb, IE3);
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = (phs == -1)? hq + d: hb[d];
			// An orphan exon is disallowed
			if (!from->dir || from->dir & SPIN) continue;
			x = from->val + sigJ;
#if INTR > 1
			WRECD*	phl = hl[d][phs + 1];
			int*	pnx = nx[d][phs + 1];
			int	l = INTR;
			while (--l >= 0) {
			    if (x > phl[pnx[l]].val)
				gswap(pnx[l], pnx[l + 1]);
			    else
				break;
			}
			if (++l < INTR) {
			    phl[pnx[INTR]].dir = 0;
			    phl += pnx[l];
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = nb;
			    phl->sig = sigJ;
			}
#else
			phl = hl[d] + phs + 1;
			if (x > phl->val) {
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = nb;
			    phl->sig = sigJ;
			}
#endif
		    }
		    if (bb->phs3 - phs == 1) {
			phs = -1;
			goto AccRev;
		    }
		}

#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y);
		putvar(g->val); putvar(eq1->val);
		if (g2) {
		    putvar(g2->val); putvar(eq2->val);
		}
		if (algmode.lsg) {
#if INTR > 1
		    putvar(hl[0][2][nx[0][2][0]].val); 
		    putvar(hl[0][1][nx[0][1][0]].val);
		    putvar(hl[0][0][nx[0][0][0]].val);
#else
		    putvar(hl[0][2].val); 
		    putvar(hl[0][1].val);
		    putvar(hl[0][0].val);
#endif
		    printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
		    printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		}
		putchar('\n');
	}
#endif
		if (m == mm) {
		    if (isvert(mx) && g->val < mx->val) *g = *mx;
		    (*--f1) = *eq1;
		    if (f2) (*--f2) = *eq2;
		}
	    }	// end of n loop
	    if (a_in_zone) --api;		// has exon-exon junctions
	}	// end of m loop

/*******	forward phase	*******/

	finitH_ng(hhf, mm);
	m = a->left;
	if (!a->inex.exgl) --m; // global
	n1 = 3 * m + wdw->lw - 1;
	n2 = 3 * m + wdw->up;
	as = a->at(m);
	if (api_size) api.reset(m);
	for ( ; ++m <= mm; ++as) {
	    n1 += 3;
	    n2 += 3;
	    int		n0 = max(n1, b->left);
	    int		n9 = min(n2, b->right);
	    int		n = n0;
	    CHAR*	bs = b->at(n);
	    EXIN*	bb = b->exin->score(n);
	    int		r = n - 3 * m;
	    int		l = 0;
	    WRECD*	h = hhf[l] + r;
	    WRECD*	g = hhf[++l] + r;
	    WRECD*	g2 = (pwd->Noll == 3)? hhf[++l] + r: 0;
	    WRECD*	sj = algmode.lsg? hhf[++l] + r: 0;
	    for (int q = 0; q < NQUE; ++q) e1[q] = e2[q] = oomW;
	    if (m == mm) {
		f1 = hhb[++l] + r;
		if (pwd->Noll == 3) f2 = hhb[++l] + r;
		for (int d = 0; d < pwd->Nrow; ++d) hb[d] = hhb[d] + r;
	    }
	    WRECD	hq[NOL] = {oomW};
	    for (int d = 0; d < pwd->Noll; ++d) {
	      for (int p = 0; p < 3; ++p) {
#if INTR > 1
		WRECD*	phl = hl[d][p];
		for (int l = 0; l <= INTR; ++l, ++phl) {
		  nx[d][p][l] = l;
		  *phl = oomW;
		}
#else
		hl[d][p] = oomW;
#endif
	      }
	    }
#if DEBUG
	    if (algmode.nsa & 8) {
		printf("%2d %2d %2d", m, n, h->dir);
		putvar(h->val); putchar('\n');
	    }
#endif
	    bool	a_in_zone = api.eq(m);
	    int	q = 0;
	    for ( ; n <= n9; ++n, ++h, ++g, ++r, ++bb, ++bs) {
		VTYPE	x, y;
		WRECD*	eq1 = e1 + q;
		WRECD*	eq2 = 0;
		hf[0] = h;
		hf[1] = eq1;
		if (g2) hf[2] = eq2 = e2 + q;
		VTYPE	sigE = bb[-2].sigE;
		WRECD*	from = h;
		WRECD*	mx = h;
		if (n == n0) goto FindCenter;
		if (algmode.lsg) ++sj;
		if (g2) ++g2;

//	diagonal match
		if (m == a->left) goto HorizonF;
		if (n > b->left + 2) {
		    *hq = *h;
		    if (sj->dir) {
			*h = *sj;
			sj->dir = 0;
		    }
		    else h->val += pwd->sim2(as, bs - 2) + sigE;
		    h->dir = (from->dir & DIAG)? DIAG: NEWD;
		} else	*h = oomW;

//	vertical gap extention
		y = g[3].val + pwd->BasicGEP;

//	1 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA2;
		    g->lst = r + 1;
		} else	g->val = y;

//	2 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA1;
		    g->lst = r + 2;
		}

//	nomal deletion
		x = (++from)->val + pwd->GapW3;
		if (x >= g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		    g->lst = r + 3;
		} else if (y >= g->val) {
		    *g = g[3];
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		}
		if (g->val >= mx->val) mx = g;

//	long deletion
		if (g2) {
		  x = from->val + pwd->GapW3L;
		  y = g2[3].val + pwd->LongGEP;
		  if (x >= y) {
		    *g2 = *from;
		    g2->val = x;
		    g2->dir = VERL | (g2->dir & SPJC);
		    g2->lst = r + 3;
		  } else {
		    *g2 = g2[3];
		    g2->val = y;
		  }
		  if (g2->val >= mx->val) mx = g2;
		}
HorizonF:
//	nomal insertion
		hq[1] = *eq1;
		if (n > n0 + 2) {
		    from = h - 3;
		    x = from->val + pwd->GapW3;
		    y = eq1->val += pwd->BasicGEP;
		    if (x > y) {
			*eq1 = *from;
			eq1->val = x;
			eq1->lst = r - 3;
		    }
		    if (!(eq1->dir & SPF2))	eq1->val += sigE;
		    eq1->dir = (eq1->dir & SPIN) + HORI;
//	long insertion
		    if (eq2) {
			hq[2] = *eq2;
			x = from->val + pwd->GapW3L;
			y = eq2->val += pwd->LongGEP;
			if (x > y) {
			    *eq2 = *from;
			    eq2->val = x;
			    eq2->lst = r - 3;
			}
			if (!(eq2->dir & SPF2)) eq2->val += sigE;
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
			eq1->lst = r - 2;
			eq1->dir = (eq1->dir & SPIN) + HOR2;
		    }
		}
//	1 nt insertion
		from = h - 1;
		x = from->val + pwd->GapW1;
		if (x > eq1->val) {
		    *eq1 = *from;
		    eq1->val = x;
		    eq1->lst = r - 1;
		    eq1->dir = (eq1->dir & SPIN) + HOR1;
		}
		if (eq1->val > mx->val) mx = e1 + q;
		if (++q == NQUE) q = 0;

//	intron 3' boundary, assume no overlapping signals
		if (isEIJ(bb->phs3)) {
		    int	phs = (bb->phs3 == 2)? 1: bb->phs3;
AccFwd:
		    if (phs == -1)
			sj->val = mx->val + pwd->sim2(as + 1, bs + 1); 
		    int		nb = n - phs;
		    VTYPE	sigJ = a_in_zone? api.match_score(3 * m - phs): 0;
		    WRECD*	maxdhl = 0;
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = hf[d];
#if INTR > 1
			WRECD*	maxphl = 0;
			int*	pnx = nx[d][phs + 1];
			for (l = 0; l < INTR; ++l) {
			    WRECD*	phl = hl[d][phs + 1] + pnx[l];
			    if (!phl->dir) break;
			    y = sigJ + phl->val + 
				b->exin->sig53(phl->jnc, nb, IE53) +
				pwd->IntPen->Penalty(nb - phl->jnc);
			    if (phs == 1) {
				CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
				y += pwd->pmt->prematT(cs);
				if (phl->dir & DIAG) y += pwd->sim2(as, cs);
			    } else if (phs == -1) {
				CHAR*	cs = spjcs->spjseq(phl->jnc, nb) + 1;
				y += pwd->pmt->prematT(cs);
				x = y + pwd->sim2(as + 1, cs);
				if (x > sj->val) {
				    maxdhl = phl;
				    sj->val = x;
				}
			    }
			    if (y > from->val) {
				from->val = y;
				maxphl = phl;
		    	    }
			}
			if (maxphl) {
			    phs3 = phs;
		 	    from->jnc = maxphl->jnc;
			    from->dir = maxphl->dir | SPJCI;
			    if (phs == -1) from->dir |= SPF2;
			    from->upr = max(maxphl->upr, r);
			    from->lwr = min(maxphl->lwr, r);
			    from->lst = maxphl->lst + nb - maxphl->jnc;
			    if (ge(from->val, mx->val)) mx = from;
			}
#else
			WRECD*	phl = hl[d] + phs + 1;
			if (!phl->dir) break;
			y = sigJ + phl->val + 
			    b->exin->sig53(phl->jnc, nb, IE53) +
			    pwd->IntPen->Penalty(nb - phl->jnc);
			if (phs == 1) {
			    CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
			    y += pwd->pmt->prematT(cs);
			    if (isdiag(phl)) y += pwd->sim2(as, cs);
			} else if (phs == -1) {
			    CHAR*	cs = spjcs->spjseq(phl->jnc, nb) + 1;
			    y += pwd->pmt->prematT(cs);
			    x = y + pwd->sim2(as + 1, cs);
			    if (x > sj->val) {
				maxdhl = phl;
				sj->val = x;
			    }
			}
			if (y > from->val) {
			    phs3 = phs;
			    from->val = y;
		 	    from->jnc = phl->jnc;
			    from->dir = phl->dir | SPJCI;
			    if (phs == -1) from->dir |= SPF2;
			    from->upr = max(phl->upr, r);
			    from->lwr = min(phl->lwr, r);
			    from->lst = phl->lst + nb - phl->jnc;
			    if (ge(from->val, mx->val)) mx = from;
			}
#endif
		    }
		    if (maxdhl) {
			sj->jnc = maxdhl->jnc;
			sj->dir = maxdhl->dir | SPALL;
			sj->upr = max(maxdhl->upr, r);
			sj->lwr = min(maxdhl->lwr, r);
			sj->lst = maxdhl->lst + nb - maxdhl->jnc;
		    }
		    if (bb->phs3 - phs == 1) {	// AGAG
			phs = -1;
			goto AccFwd;
		    }
		}

//	Find optimal path
		y = h->val;
		if (h != mx) {	// non-diagonal
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		}
		if (LocalL && h->val <= 0) h->val = h->dir = 0;

//	intron 5' boundary

		if (isEIJ(bb->phs5)) {
		    int		phs = (bb->phs5 == 2)? 1: bb->phs5;
Donor:
		    int	nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(nb, 0, IE5);
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = (phs == 1)? hq + d: hf[d];
			// An orphan exon is disallowed
			if (!from->dir || (from->dir & SPIN) || isvert(from))
			    continue;
			x = from->val + sigJ;
			if (phs == 1 && nb - from->jnc == 2 && bb[-2].sigS > 0)
			    x -= bb[-2].sigS;
#if INTR > 1
			WRECD*	phl = hl[d][phs + 1];
			int*	pnx = nx[d][phs + 1];
			int	l = INTR;
			while (--l >= 0) {
			    if (x > phl[pnx[l]].val)
				gswap(pnx[l], pnx[l + 1]);
			    else
				break;
			}
			if (++l < INTR) {
			    phl[pnx[INTR]].dir = 0;
			    phl += pnx[l];
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = nb;
			}
#else
			phl = hl[d] + phs + 1;
			if (x > phl->val) {
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = nb;
			}
#endif
		    }
		    if (bb->phs5 == 2 && phs == 1) {	//	GTGT..
			phs = -1;
			goto Donor;
		    }
		}

#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); putvar(g->val); putvar(eq1->val);
		if (g2) {
		    putvar(g2->val); putvar(eq2->val);
		}
		if (algmode.lsg) {
#if INTR > 1
		    putvar(hl[0][2][nx[0][2][0]].val); 
		    putvar(hl[0][1][nx[0][1][0]].val);
		    putvar(hl[0][0][nx[0][0][0]].val);
#else
		    putvar(hl[0][2].val); 
		    putvar(hl[0][1].val);
		    putvar(hl[0][0].val);
#endif
		    printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
		    printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		}
		putchar('\n');
	}
#endif

// Find Center
FindCenter:
		if (m == mm) {
		    int	acc = (*hf)->dir & SPJC;
		    int	dnr = (*hb)->dir & SPJC;
		    int	j1, k1, l;
		    WRECD*	sjb = hb[pwd->Noll];
		    if (((*hf)->dir & SPIN) && ((*hb)->dir & SPIN))
			goto NextCycle;
		    x = (*hf)->val + (*hb)->val;
		    if (((*hf)->dir & SPF2) && isntvert(*hb)) x -= bb[1].sigE;
		    if (((*hb)->dir & SPF2) && isntvert(*hf)) x -= bb[-2].sigE;
		    if (ishori(*hf) && ishori(*hb)) x -= pwd->BasicGOP;
		    if (sj->dir && isdiag(*hb)) { 	// phs == -1
			y = sj->val + (*hb)->val - pwd->sim2(as + 1, bs + 1) - bb[1].sigE;
			if (y > x || (*hf)->dir & SPF2) {
			    x = y;
			    acc = SPJC;
			}
		    }
		    if (sjb->dir && isdiag(*hf)) {	// phs === 1
			y = sjb->val + (*hf)->val - pwd->sim2(as, bs - 2) - bb[-2].sigE;
			if (y > x || (*hb)->dir & SPF2) {
			    x = y;
			    dnr = SPJC;
			}
		    }
		    j1 = eq1->dir & SPJC;	// 3' boundary
		    k1 = f1->dir & SPJC;	// 5' boundary
		    VTYPE	z;
		    if ((j1 && !k1) || (!j1 && k1)) {
			z = eq1->val;
			y = f1->val;
			if (eq1->dir & SPF2) z -= bb[1].sigE;
			if (f1->dir & SPF2) y -= bb[-2].sigE;
			y += z - pwd->BasicGOP;
			l = f1->lst - eq1->lst - pwd->codonk1;	// gap length
			if (l > 0) y += pwd->diffu * l;		// > codonK1 ?
			if (y > x) {
			    x = y;
			    dnr = k1;
			    acc = j1;
			}
		    }
		    if (f2) {
			int	j2 = eq2->dir & SPJC;
			int	k2 = f2->dir & SPJC;
			if ((j1 && !k2) || (!j1 && k2)) {
			    z = eq1->val;
			    y = f2->val;
			    if (eq1->dir & SPF2) z -= bb[1].sigE;
			    if (f2->dir & SPF2)  y -= bb[-2].sigE;
			    l = r - eq1->lst;
			    y += z - pwd->BasicGOP + pwd->diffu * l;
			    if (y > x) {
				x = y;
				dnr = k2;
				acc = j1;
			    }
			}
			if ((j2 && !k1) || (!j2 && k1)) {
			    z = eq2->val;
			    y = f2->val;
			    if (eq2->dir & SPF2) z -= bb[1].sigE;
			    if (f1->dir & SPF2)   y -= bb[-2].sigE;
			    l = f1->lst - r;
			    y += z - pwd->BasicGOP + pwd->diffu * l;
			    if (y > x) {
				x = y;
				dnr = k1;
				acc = j2;
			    }
			}
			if ((j2 && !k2) || (!j2 && k2)) {
			    z = eq2->val;
			    y = f2->val;
			    if (eq2->dir & SPF2) z -= bb[1].sigE;
			    if (f2->dir & SPF2)  y -= bb[-2].sigE;
			    y += z - pwd->LongGOP;
			    if (y > x) {
				x = y;
				dnr = k2;
				acc = j2;
			    }
			}
		    }
		    k1 = isvert(*hf)? 1: 0;
		    j1 = isvert(*hb)? 1: 0;
		    y = g->val + hb[1]->val - pwd->BasicGOP;
		    l = g->lst - hb[1]->lst - pwd->codonk1;
		    if (l > 0) y += pwd->diffu * l;
		    if (y >= x) {x = y; k1 = j1 = 1;}
		    if (g2) {
			y = g->val + hb[2]->val - pwd->BasicGOP + pwd->diffu * (g->lst - r);
			if (y > x) {x = y; k1 = 1; j1 = 2;}
			y = g2->val + hb[1]->val - pwd->BasicGOP + pwd->diffu * (r - hb[1]->lst);
			if (y > x) {x = y; k1 = 2; j1 = 1;}
			y = g2->val + hb[2]->val - pwd->LongGOP;
			if (y > x) {x = y; k1 = j1 = 2;}
		    }
		    if (gt(x, mxh)) {
			mxh = x;
			kk = k1 + 4 * j1;
			rr[jtigh = 0] = r;
			jb = ((*hb)->dir & SPJC)? 1: 0;
			kb = ishori(*hb);
			pp[0] = (jb && !vmf);
		    } else if (ge(x, mxh)) {
			if (j1 && dnr) rr[jtigh = 0] = r;
			else if (j1 && acc) rr[jtigh = 1] = r;
			else if 
			(((kb || jb == 2) && ishori(*hf)) ||		// --*
			(jb == 1 && (acc || ((*hf)->dir & SPJC))) ||	// ..|
			(kb && dnr)) {					// *-|
			    if (++jtigh < MAXR) rr[jtigh] = r;
			    if (acc || dnr) {++jb; pp[jtigh] = 1;}
			    else	jb = pp[jtigh] = 0;
			    kb = 0;
			}
		    }
NextCycle:
		    for (int d = 0; d < pwd->Nrow; ++d) ++hb[d];
		    ++f1;
		    if (f2) ++f2;
		}	// was center
	    }	// end of n-loop
	    if (a_in_zone) ++api;		// has exon-exon junctions
	}	// end of m-loop

	int	k1 = kk % 4;
	int	k2 = kk / 4;
	int	r = rr[0];
	*nl = r + 3 * mm;
	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hhf[k] + r;
	wdwf->up = max(hf[k1]->upr, hf[k1]->lst);
	wdwf->lw = min(hf[k1]->lwr, hf[k1]->lst);
	wdwf->width = wdwf->up - wdwf->lw + 7;
	SKL	wskl;
	if (k1) {
	    wskl.m = *ml = (*nl - hf[k1]->lst) / 3;
	    wskl.n = *nl;
	    mfd->write((UPTR) &wskl);
	} else {
	    wskl.m = *ml = mm;
	}
	*nr = (r = rr[jtigh]) + 3 * mm;
	for (int k = 0; k < pwd->Nrow; ++k) hb[k] = hhb[k] + r;
	wdwb->up = max(hb[k2]->upr, hb[k2]->lst);
	wdwb->lw = min(hb[k2]->lwr, hb[k2]->lst);
	wdwb->width = wdwb->up - wdwb->lw + 7;
	if (k2) {
	    wskl.n = *nr;
	    mfd->write((UPTR) &wskl);
	    wskl.m = *mr = (*nr - hb[k2]->lst) / 3;
	    mfd->write((UPTR) &wskl);
	} else {
	    *mr = mm;
	    for ( ; jtigh >= 0; --jtigh) {
		wskl.m = mm;
		wskl.n = rr[jtigh] + 3 * mm;
		mfd->write((UPTR) &wskl);
	    }
	}

	delete[] (hhf[0] + wdw->lw - 3);
	delete[] (hhb[0] + wdw->lw - 3);
	return (mxh);
}

static const int pNoll = 2;
static const int pNrow = 3;
static const int Nedge = 3;

void Aln2h1::pfinitH_ng(LRECD* hhg[])
{
	int	r = b->left - 3 * a->left;
	int	rr = b->right - 3 * a->left;
	LRECD*	hf[pNrow];

	for (int k = 0; k < pNrow; ++k) hf[k] = hhg[k] + r;
	if (wdw->up < rr) rr = wdw->up;
	LRECD*	h = hf[0]++;
	h->val = 0;
	h->dir = DIAG;
	h->jnc = b->left;
	h->lst = r;
	h->ptr = vmf->add(a->left, b->left, 0L);
 	for (int k = 1; k < pNrow; ++k) *hf[k]++ = ooml;
	while (++r <= wdw->up + 3) {
 	    for (int k = 0; k < pNrow; ++k) *hf[k]++ = ooml;
	}

	r = b->left - 3 * a->left;
	rr = b->left - 3 * a->right;
	for (int k = 0; k < pNrow; ++k) hf[k] = hhg[k] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	for (int i = 1; --r >= wdw->lw - 3; ++i) {
	    h = --hf[0];
	    if (r < rr) *h = ooml;
	    else {			// (semi) global
		if (i <= 3) {
		    *h = h[i];
		    h->val += pwd->GapPenalty(i);
		    if (i < 3) h->val += pwd->ExtraGOP;
		    h->dir = VERT;
		} else {
		    *h = h[3];
		    h->val += pwd->BasicGEP;
		}
	    }
	}
}

void Aln2h1::pbinitH_ng(LRECD* hhg[])
{
	int	r = b->right - 3 * a->right;
	int	rr = b->left - 3 * a->right;
	LRECD*	hb[pNrow];

	for (int k = 0; k < pNrow; ++k) hb[k] = hhg[k] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	LRECD*	h = hb[0]--;
	h->val = 0;
	h->dir = DIAG;
	h->jnc = b->right;
	h->lst = r;
	h->ptr = vmf->add(a->right, b->right, 0L);
	for (int d = 1; d < pNrow; ++d) *hb[d]-- = ooml;
	while (--r >= wdw->lw - 3) {
	    for (int d = 0; d < pNrow; ++d) *hb[d]-- = ooml;
	}

	r = b->right - 3 * a->right;
	rr = b->right - 3 * a->left;
	if (wdw->up < rr) rr = wdw->up;
	for (int k = 0; k < pNrow; ++k) hb[k] = hhg[k] + r;
	for (int i = 1; ++r <= wdw->up + 3; ++i) {
	    h = ++hb[0];
	    if (r > rr) *h = ooml;
	    else {
		if (i <= 3) {
		    *h = h[-i];
		    h->val += pwd->GapPenalty(i);
		    if (i < 3) h->val += pwd->ExtraGOP;
		    h->dir = VERT;
		} else {
		    *h = h[-3];
		    h->val += pwd->BasicGEP;
		}
	    }
	}
}

VTYPE Aln2h1::pincersH_ng(long *ptr, int cmode)
{
	int	nr[3];
	LRECD*	hhg[pNrow];
	LRECD*	hf[Nedge];
	LRECD	e1[NQUE];
	LRECD	maxl = ooml;		// left  side
	LRECD	maxr = ooml;		// right side
	VSKLP	maxh = {NEVSEL, 0, 0, 0L};	// start point of traceback
	VTYPE	maxval = NEVSEL;	// for x-drop off
	VTYPE	maxscr = NEVSEL;	// overall
	int	dd = a->right - a->left + 1;
	EIJUNCP* hhl = new EIJUNCP[dd] - a->left;
	EIJODRP* nnl = new EIJODRP[dd] - a->left;
	int	n1, n2;
	CHAR*	as;

	hhg[0] = new LRECD[pNrow * wdw->width] - wdw->lw + 3;
	for (int d = 1; d < pNrow; ++d) hhg[d] = hhg[d-1] + wdw->width;
	vmf->add(0, 0, 0L);	// Skip 0-th record
	int	m = a->left;
	PfqItr	api(a, m);
	int	api_size = api.size();
	if (cmode == 2) goto forward;	// forward only

/*******	backward phase	*******/

	pbinitH_ng(hhg);
	m = a->right;
	if (!a->inex.exgr) ++m; // global
	else {
	    for (int p = 0; p < 3; ++p) {
		for (int l = 0; l <= NCANDS; ++l) {
		    hhl[m][p][l] = ooml;
		    nnl[m][p][l] = l;
		}
	    }
	}
	as = a->at(m);
	n1 = 3 * m + wdw->lw;
	n2 = 3 * m + wdw->up + 1;
	while (--m >= a->left) {
	    for (int p = 0; p < 3; ++p) {
		for (int l = 0; l <= NCANDS; ++l) {
		    hhl[m][p][l] = ooml;
		    nnl[m][p][l] = l;
		}
	    }
	    --as; n1 -= 3; n2 -= 3;
	    int		n0 = min(n2, b->right);
	    int		n9 = max(n1, b->left);
	    int		n  = n0;
	    CHAR*	bs = b->at(n);
	    EXIN*	bb = b->exin->score(n);
	    int		r = n - 3 * m;
	    for (int q = 0; q < NQUE; ++q) e1[q] = ooml;
	    if (!b->inex.exgr && n == b->right && m == a->right) {
		e1[2] = hhg[0][r];
		e1[2].val = pwd->GapW3;
	    }
	    LRECD*	h = hhg[0] + r;
	    LRECD*	g = hhg[1] + r;
	    LRECD*	mxd = ((h->val + pwd->Vthr) < maxval)? &ooml: h;
	    LRECD	hq[NOL] = {ooml};
	    int		count3 = 0;
	    for (int p = 0; p < 3; ++p) nr[(n + p) % 3] = n + p;
	    nr[n % 3] = n + 3;

#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m+1, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    int	peak = 0;
	    int	q = 0;		// queue pointer
	    while (--n >= n9) {
		VTYPE	x, y;
		--bs; --r;
		hf[0] = --h;
		hf[1] = --g;
		LRECD*	eq1 = e1 + q;
		hf[2] = eq1;
		VTYPE	sigE = (bb--)->sigE;

//	diagonal match
		LRECD*	from = h;
		LRECD*	mx = h;
		if (m == a->right) goto HorizonB;
		if (n < b->right - 2) {
		    *hq = *h;
		    h->val += pwd->sim2(as, bs + 1) + sigE;
		    h->dir = isdiag(from)? DIAG: NEWD;
		} else	*h = ooml;

//	vertical gap extension
		y = g[-3].val + pwd->BasicGEP;

//	1 nt deletion
		--from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA2;
		    g->lst = r - 1;
		} else	g->val = y;

//	2 nt deletion
		--from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val)	{
		    *g = *from;
		    g->val = x;
		    g->dir = SLA1;
		    g->lst = r - 2;
		}

//	normal deletion
		x = (--from)->val + pwd->GapW3;
		if (x >= g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		    g->lst = r - 3;
		} else if (y >= g->val) {
		    *g = g[-3];
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		}
		if (g->val >= mx->val) mx = g;

HorizonB:
//	normal insertion
		hq[1] = *eq1;
		if (n < n0 - 2) {
		    from = h + 3;
		    x = from->val + pwd->GapW3;
		    y = eq1->val += pwd->BasicGEP;
		    if (x > y) {
			*eq1 = *from;
			eq1->val = x;
			eq1->lst = r + 3;
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
		if (cmode == 1 && m == a->left && bb[1].sigS > 0)
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

//	intron 3' boundary reverse phase

		if (cmode == 3 && isEIJ(bb->phs3)) {
		    int	phs = (bb->phs3 == 2)? 1: bb->phs3;
AccRev:
		    int		nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(0, nb, IE3);
		    LRECD*	phl = hhl[m][phs + 1];
		    int*	pnx = nnl[m][phs + 1];
		    int 	p = 0;
		    for ( ; p < NCANDS && phl[pnx[p]].dir; ++p) ;
		    for (int d = 0; d < Nedge; ++d) {
			from = (phs == -1)? hq + d: hf[d];
			int	l = (d + 1) / 2;
			if (!from->dir || (d && 
			(mx == from || from->val <= mx->val + pwd->GOP[l])))
			    continue;
			x = from->val + sigJ;
			if (p == NCANDS) --p;
			for (l = p++; --l >= 0; ) {
			    if (x > phl[pnx[l]].val)
				gswap(pnx[l], pnx[l + 1]);
			    else
				break;
			}
			if (++l < NCANDS) {
			    phl[pnx[NCANDS]].dir = 0;
			    LRECD*	phr = phl + pnx[l];
			    *phr = *from;
			    phr->val = x;
			    phr->jnc = nb;
			}
		    }
		    if (bb->phs3 - phs == 1) {
			phs = -1;
			goto AccRev;
		    }
		}
#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(eq1->val);
		if (algmode.lsg) {
		    putvar(hhl[m][2][nnl[m][2][0]].val);
		    putvar(hhl[m][1][nnl[m][1][0]].val);
		    putvar(hhl[m][0][nnl[m][0][0]].val);
		    printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
		    printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		}
		putchar('\n');
	}
#endif
	    }	// end of n loop
	    if (!mxd->dir) break;	// no peak
	    if (peak) n1 = n;
	    LRECD*	phr = hhl[m][0] + nnl[m][0][NCANDS];
	    *phr = *mxd;	// maxscr in row
	}	// end of m loop
	if (cmode == 1) {
	    maxscr = maxh.val;
	    ptr[0] = 0L;
	    ptr[1] = vmf->add(maxh.m, maxh.n, maxh.p);
	    goto freeprec;
	} else	++m;

/*******	forward phase	*******/

forward:
	pfinitH_ng(hhg);
	if (!a->inex.exgl) --m; // global
	n1 = 3 * m + wdw->lw - 1;
	n2 = 3 * m + wdw->up;
	as = a->at(m);
	if (api_size) api.reset(m);
	maxh.val = maxval = NEVSEL;
	for ( ; ++m <= a->right; ++as) {
	    n1 += 3;
	    n2 += 3;
	    int		n0 = max(n1, b->left);
	    int		n9 = min(n2, b->right);
	    int		n  = n0;
	    int		r = n - 3 * m;
	    int		count3 = 0;
	    for (int q = 0; q < NQUE; ++q) e1[q] = ooml;
	    if (!b->inex.exgl && n == b->left && m == a->left) {
		e1[2] = hhg[0][n - 3 * m];
		e1[2].val = pwd->GapW3;
	    }
	    CHAR*	bs = b->at(n);
	    EXIN*	bb = b->exin->score(n);
	    LRECD*	h = hhg[0] + r;
	    LRECD*	g = hhg[1] + r;
	    LRECD	hq[NOL] = {ooml};
	    for (int p = 0; p < 3; ++p) nr[(n - p) % 3] = n - p;
	    LRECD*	mxd = ((h->val + pwd->Vthr) < maxval)? &ooml: h;
	    nr[n % 3] = n - 3;

#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif

	    int		peak = 0;
	    int		q = 0;		// queue pointer
	    bool	a_in_zone = api_size && api.eq(m);
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x, y;
		++r; ++bb;
		LRECD*	eq1 = e1 + q;
		hf[0] = ++h;
		hf[1] = ++g;
		hf[2] = eq1;

//	diagonal match
		VTYPE	sigE = bb[-2].sigE;
		LRECD*	from = h;
		LRECD*	mx = h;
		if (m == a->left) goto HorizonP;
		if (n > b->left + 2) {
		    *hq = *h;
		    h->val += pwd->sim2(as, bs - 1) + sigE;
		    h->dir = isdiag(from)? DIAG: NEWD;
		} else	*h = ooml;

//	vertical gap extention
		y = g[3].val + pwd->BasicGEP;

//	1 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
		if (x > y) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA2;
		    g->lst = r + 1;
		} else	g->val = y;

//	2 nt deletion
		++from;
		x = from->val + (isvert(from)? pwd->GapE2: pwd->GapW2);
		if (x > g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = SLA1;
		    g->lst = r + 2;
		}

//	normal deletion
		x = (++from)->val + pwd->GapW3;
		if (x >= g->val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		    g->lst = r + 3;
		} else	if (y >= g->val) {
		    *g = g[3];
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		}
		if (g->val >= mx->val) mx = g;

HorizonP:
//	normal insertion
		hq[1] = *eq1;
		if (n > n0 + 2) {
		    from = h - 3;
		    bool	l = cmode == 2 && m == a->right;
		    if (l) {
			l = l && bb[-2].sigT > 0;
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
			eq1->lst = l? r: r - 3;
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
			eq1->lst = r - 2;
		    }
		}
//	1 nt insertion
		from = h - 1;
		x = from->val + pwd->GapW1;
		if (x > eq1->val) {
		    *eq1 = *from;
		    eq1->val = x;
		    eq1->dir = HOR1;
		    eq1->lst = r - 1;
		}
		if (eq1->val >= mx->val) mx = e1 + q;
		if (++q == NQUE) q = 0;

//	Find optimal path
		if (mx->dir == NEWD) {
		    mx->lst = r;
		    mx->ptr = vmf->add(m - 1, n - 3, mx->ptr);
		}
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

//	intron 5' boundary
		if (cmode == 3 && isEIJ(bb->phs5)) {
		    int	phs = (bb->phs5 == 2)? 1: bb->phs5;
Donor:
		    int		nb = n - phs;
		    VTYPE sigJ = a_in_zone? api.match_score(3 * m - phs): 0;
		    for (int d = 0; d < Nedge; ++d) {
			int	l = (d + 1) / 2;
			from = (phs == 1)? hq + d: hf[d];
			if (!from->dir || (d &&
			    (mx == from || from->val <= mx->val + pwd->GOP[l])))
			    continue;
			for (l = 0; l < NCANDS; ++l) {
			    LRECD*	phr = hhl[m][phs + 1] + nnl[m][phs + 1][l];
			    if (!phr->dir) break;
			    x = from->val + phr->val + sigJ +
				pwd->IntPen->Penalty(phr->jnc - nb) +
				b->exin->sig53(nb, phr->jnc, IE35);
			    if (phs == 1) {
				CHAR*	cs = spjcs->spjseq(nb, phr->jnc);
				x += pwd->pmt->prematT(cs) +
					pwd->sim2(as, cs);
				if (nb - from->jnc == 2 && bb[-2].sigS > 0)
				    x -= bb[-2].sigS;
			    } else if (phs == -1) {
				CHAR*	cs = spjcs->spjseq(nb, phr->jnc) + 1;
				x += pwd->pmt->prematT(cs) +
					pwd->sim2(as + 1, cs);
			    }
			    if  ((isvert(from) && isvert(phr)) ||
				(ishori(from) && ishori(phr)))
				    x -= pwd->BasicGOP;
			    if (x > maxscr) {
				maxscr = x;
				maxl = *from;
				maxr = *phr;
				maxl.ptr = vmf->add(m, n, from->ptr);
				maxl.ptr = vmf->add(m, phr->jnc + phs, maxl.ptr);
			    }
			}
		    }
		    if (bb->phs5 == 2 && phs == 1) {	//	GTGT..
			phs = -1;
			goto Donor;
		    }
		}
#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y);
		putvar(g->val); putvar(eq1->val);
		if (algmode.lsg) {
		    printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
		    printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		}
		putchar('\n');
	}
#endif

	    }	// end of n-loop
	    if (peak) n2 = n; 
	    if (cmode == 3) {
		LRECD*	phr = hhl[m][0] + nnl[m][0][NCANDS];
		int	l = phr->lst - mxd->lst;
		if (l >= 0 && l < pwd->MaxGapL && l % 3 == 0) {
		    VTYPE	x = mxd->val + phr->val +
			+ SumCodePot(bb - r + mxd->lst, l, 0, pwd) + pwd->GapPenalty3(l);
		    if (x > maxscr) {
			maxscr = x;
			maxl = *mxd;
			maxr = *phr;
			maxl.ptr = vmf->add(m, 3 * m + mxd->lst, mxd->ptr);
			maxl.ptr = vmf->add(m, 3 * m + phr->lst, maxl.ptr);
		    }
		}
	    }
	    if (!mxd->dir) break;
	    if (a_in_zone) ++api;
	}	// end of m-loop
	if (cmode == 2) {
	    maxscr = maxh.val;
	    ptr[0] = vmf->add(maxh.m, maxh.n, maxh.p);
	    ptr[1] = 0L;
	} else {
	    ptr[0] = maxl.ptr;
	    ptr[1] = maxr.ptr;
	}

freeprec:
	delete[] (hhg[0] + wdw->lw - 3);
	delete[] (hhl + a->left);
	delete[] (nnl + a->left);
	return (maxscr);
}

void Aln2h1::cinitH_ng(CRECD* hhc[])
{
	int	n = b->left;
	int	r = b->left - 3 * a->left;
	int	r0 = r;
	int	rr = b->right - 3 * a->left;
	CRECD*	hf[NOL + 1];
	EXIN*	bb = b->exin->score(n);

	for (int d = 0; d < pwd->Nrow; ++d) hf[d] = hhc[d] + r;
	if (wdw->up < rr) rr = wdw->up;
	for (int i = 0; r <= wdw->up + 3; r++, bb++, n++, i++) {
	    for (int d = 1; d < pwd->Nrow; ++d) *hf[d]++ = oomC;
	    if (r > rr) {
		*hf[0]++ = oomC;
		continue;
	    }
	    CRECD*	h = hf[0]++;
	    if (i >= 3) {
		*h = h[-3];
		h->dir = HORI;
		h->val += bb[-2].sigE + pwd->GapExtPen3(n - h->nlb);
	    }
	    VTYPE	x = 0;
	    if ((algmode.lcl & 1) && bb[1].sigS > x) x = bb[1].sigS;
	    if ((algmode.lcl & 4) && bb->sig3 > x) x = bb->sig3;
	    if (h->val <= x) {
		h->val = x;
		h->lwr = r;
		h->mlb = a->left;
		h->nlb = n;
		h->jnc = n;
		h->dir = DEAD;
	    }
	    h->upr = r;
	}

	r = r0;
	rr = b->left - 3 * a->right;
	for (int d = 0; d < pwd->Nrow; ++d) hf[d] = hhc[d] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	while (--r >= wdw->lw - 3) {
	    for (int d = 0; d < pwd->Nrow; ++d) *--hf[d] = oomC;
	    if (r < rr) continue;
	    CRECD*	h = hf[0];
	    h->val = 0;
	    h->dir = DEAD;
	    h->upr = h->lwr = r;
	    h->mlb = (b->left - r) / 3;
	    h->nlb = b->left;
	    h->jnc = b->left;
	}
}

Colonies* Aln2h1::fwdswgH_ng(VTYPE* scr)
{
	CRECD*	hhc[NOL+1];
	CRECD*	hf[NOL];	// [DIAG, HORI, HORL]
	CRECD	e1[NQUE];
	CRECD	e2[NQUE];
#if INTR > 1
	CRECD	hl[NOL][3][INTR+1]; // [DIA, HORI, HORL][phase][candidates]
	int	nl[NOL][3][INTR+1]; // [DIA, HORI, HORL][phase][candidates]
#else
	CRECD   hl[NOL][3];
#endif
	Colonies*	cl = new Colonies(0);
	COLONY*	clny = cl->at();

	hhc[0] = new CRECD[pwd->Nrow * wdw->width] - wdw->lw + 3; // +1: sp junc
	for (int k = 1; k < pwd->Nrow; ++k) hhc[k] = hhc[k-1] + wdw->width;
	cinitH_ng(hhc);
	int	m = a->left;
	CHAR*	as = a->at(m);
	PfqItr	api(a, m);
	int	api_size = api.size();
	int	n1 = 3 * m + wdw->lw - 1;
	int	n2 = 3 * m + wdw->up;
	for ( ; ++m <= a->right; ++as) {
	    n1 += 3;
	    n2 += 3;
	    int 	n0 = max(n1, b->left);
	    int 	n9 = min(n2, b->right);
	    int		n = n0;
	    int 	r = n - 3 * m;
	    CHAR*	bs = b->at(n);
	    EXIN*	bb = b->exin->score(n);
	    CRECD*	h = hhc[0] + r;
	    CRECD*	g = hhc[1] + r;
	    CRECD*	g2 = (pwd->Noll == 3)? hhc[2] + r: 0;
	    CRECD*	sj = algmode.lsg? hhc[pwd->Noll] + r: 0;
	    CRECD	hq[NOL] = {oomC};
	    for (int q = 0; q < NQUE; ++q) e1[q] = e2[q] = oomC;
	    for (int d = 0; d < pwd->Noll; ++d) {
	      for (int p = 0; p < 3; ++p) {
#if INTR > 1
		CRECD*	phl = hl[d][p];
		for (int l = 0; l <= INTR; ++l, ++phl) {
		  nl[d][p][l] = l;
		  *phl = oomC;
		}
#else
		hl[d][p] = oomC;
#endif
	      }
	    }
#if DEBUG
	    if (algmode.nsa & 8) {
		printf("%2d %2d %2d", m, n, h->dir);
		putvar(h->val); putchar('\n');
	    }
#endif
	    int	q = 0;		// queue pointer
	    bool	a_in_zone = api_size && api.eq(m);
	    for ( ; ++n <= n9; ++bs) {
		++g; ++r; ++bb;
		if (algmode.lsg) ++sj;
		CRECD*	eq1 = e1 + q;
		CRECD*	eq2 = 0;
		hf[0] = ++h;
		hf[1] = eq1;
		if (g2) {
		    ++g2;
		    hf[2] = eq2 = e2 + q;
		}

//	diagonal match
		VTYPE	sigE = bb[-2].sigE;
		CRECD*	from = h;
		CRECD*	mx = h;
		if (n > b->left + 2) {
		    *hq = *h;
		    if (sj->dir) {
			*h = *sj;
			sj->dir = 0;
		    } else	h->val += pwd->sim2(as, bs - 1) + sigE;
		    h->dir = isdiag(from)? DIAG: NEWD;
		} else	*h = oomC;

//	vertical gap extention
		VTYPE	y = g[3].val + pwd->BasicGEP;

//	1 nt deletion
		++from;
		VTYPE	x = from->val + (isvert(from)? pwd->GapE1: pwd->GapW1);
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
		} else if (y >= g->val) {
		    *g = g[3];
		    g->val = y;
		    g->dir = VERT | (g->dir & SPJC);
		}
		if (g->val > mx->val) mx = g;
//	long deletion
		if (g2) {
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
		  if (g2->val > mx->val) mx = g2;
		}
//	Horizongal
//	normal insertion
		hq[1] = *eq1;
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
		    if (eq2) {
			hq[2] = *eq2;
			x = from->val + pwd->GapW3L;
			y = eq2->val += pwd->LongGEP;
			if (x > y) {
			    *eq2 = *from;
			    eq2->val = x;
			    eq2->dir = (from->dir & SPALL) + HORL;
			}
			if (!(eq2->dir & SPF2)) eq2->val += sigE;
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
		if (isEIJ(bb->phs3)) {
		    int	phs = (bb->phs3 == 2)? bb->phs3: 1;
Acceptor:
		    if (phs == -1)
			sj->val = mx->val + pwd->sim2(as + 1, bs + 2);
		    int	nb = n - phs;
		    VTYPE	sigJ = a_in_zone? api.match_score(3 * m - phs): 0;
		    CRECD*	maxdhl = 0;
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = hf[d];
#if INTR > 1
			CRECD*	maxphl = 0;
			int*	pnl = nl[d][phs + 1];
			for (int l = 0; l < INTR; ++l) {
			    CRECD*	phl = hl[d][phs + 1] + pnl[l];
			    if (!phl->dir) break;
			    y = sigJ + phl->val 
				+ pwd->IntPen->Penalty(nb - phl->jnc)
				+ b->exin->sig53(phl->jnc, nb, IE53);
			    if (phs == 1) {
				CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
				y += pwd->pmt->prematT(cs);
				if (isdiag(phl)) y += pwd->sim2(as, cs);
			    } else if (phs == -1) {
				CHAR*	cs = spjcs->spjseq(phl->jnc, nb) + 1;
				y += pwd->pmt->prematT(cs);
				x = y + pwd->sim2(as + 1, cs);
				if (x > sj->val) {
				    maxdhl = phl;
				    sj->val = x;
				}
			    }
			    if (y > from->val) {
				from->val = y;
				maxphl = phl;
		    	    }
			}
			if (maxphl) {
		 	    from->jnc = nb;
			    from->dir = maxphl->dir | SPJCI;
			    if (phs == -1) from->dir |= SPF2;
			    from->upr = max(maxphl->upr, r);
			    from->lwr = min(maxphl->lwr, r);
			    from->mlb = maxphl->mlb;
			    from->nlb = maxphl->nlb;
			    if (from->val > mx->val) mx = from;
			}
#else
			CRECD*	phl = hl[d] + phs + 1;
			if (!phl->dir) break;
			y = sigJ + phl->val + 
			    pwd->IntPen->Penalty(nb - phl->jnc) +
			    b->exin->sig53(phl->jnc, nb, IE53);
			if (phs == 1) {
			    CHAR*	cs = spjcs->spjseq(phl->jnc, nb);
			    y += pwd->pmt->prematT(cs);
			    if (isdiag(phl)) y += pwd->sim2(as, cs);
			} else if (phs == -1) {
			    CHAR*	cs = spjcs->spjseq(phl->jnc, nb) + 1;
			    y += pwd->pmt->prematT(cs);
			    x = y + pwd->sim2(as + 1, cs);
			    if (x > sj->val) {
				maxdhl = phl;
				sj->val = x;
			    }
			}
			if (y > from->val) {
			    from->val = y;
		 	    from->jnc = nb;
			    from->dir = phl->dir | SPJCI;
			    if (phs == -1) from->dir |= SPF2;
			    from->upr = max(phl->upr, r);
			    from->lwr = min(phl->lwr, r);
			    from->mlb = phl->mlb;
			    from->nlb = phl->nlb;
			    if (from->val > mx->val) mx = from;
			}
#endif
		    }
		    if (maxdhl) {
			sj->jnc = nb;
			sj->dir = maxdhl->dir | SPALL;
			sj->upr = max(maxdhl->upr, r);
			sj->lwr = min(maxdhl->lwr, r);
			sj->mlb = maxdhl->mlb;
			sj->nlb = maxdhl->nlb;
		    }
		    if (bb->phs3 - phs == 1) {	// AGAG
			phs = -1;
			goto Acceptor;
		    }
		}

//	Find optimal path
		y = h->val;
		if (h != mx) {		// non-diagonal
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		} else if (h->val > hq->val) {
		    if (hq->val == 0 && !(h->dir & SPJC)) {
			h->upr = h->lwr = r;	// new colony
			h->mlb = m - 1;
			h->nlb = n - 3;
		    }
		    x = h->val;
		    bool	k = (algmode.lcl & 2) && n < b->right && bb[2].sigT > 0;
		    if (k) x = h->val + bb[2].sigT; // termination codon
		    else if ((algmode.lcl & 8) && bb->sig5 > 0) x += bb->sig3;
		    if (x > clny->val) {
			clny->val = x;
			clny->mrb = m;
			clny->nrb = n + (k? 3: 0);
			clny->lwr = h->lwr;
			clny->upr = h->upr + (k? 3: 0);
			clny->mlb = h->mlb;
			clny->nlb = h->nlb;
		    }
		}
		if (h->val < 0) {
		    *h = *eq1 = *g = oomC;
		    h->val = 0;
		    if (g2) *eq2 = *g2 = oomC;
		    h->clny = 0;
		}
		if (h->val >= pwd->Vthr && !h->clny)
		    h->clny = cl->next();
		if (h->clny) {
		    COLONY*	cc = h->clny;
		    if (h->val > cc->val) {
			cc->val = h->val;
			cc->mrb = m;
			cc->nrb = n;
			cc->lwr = h->lwr;
			cc->upr = h->upr;
			cc->mlb = h->mlb;
			cc->nlb = h->nlb;
		    } else if (algmode.mlt > 1 && h->val <= cc->val - pwd->Vthr) {
			*h = *eq1 = *g = oomC;	 // X-drop
			h->val = 0;
			if (g2) *eq2 = *g2 = oomC;
			h->clny = 0;
		    }
		}

//	intron 5' boundary
		if (isEIJ(bb->phs5)) {
		    int	phs = (bb->phs5 == 2)? 1: bb->phs5;
Donor:
		    int 	nb = n - phs;
		    VTYPE	sigJ = b->exin->sig53(nb, 0, IE5);
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = (phs == 1)? hq + d: hf[d];
			// An orphan exon is disallowed
			if (!from->dir || (from->dir & SPIN) || isvert(from))
			    continue;
			x = from->val + sigJ;
			if (phs == 1 && nb - from->jnc == 2 && bb[-2].sigS > 0)
			    x -= bb[-2].sigS;
#if INTR > 1
			CRECD*	phl = hl[d][phs + 1];
			int*	pnl = nl[d][phs + 1];
			int	l = INTR;
			while (--l >= 0) {
			    if (x > phl[pnl[l]].val)
				gswap(pnl[l], pnl[l + 1]);
			    else
				break;
			}
			if (++l < INTR) {
			    phl[pnl[INTR]].dir = 0;
			    phl += pnl[l];
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = nb;
			}
#else
			CRECD*	phl = hl[d] + phs + 1;
			if (x > phl->val) {
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = nb;
			}
#endif
		    }
		    if (bb->phs5 == 2 && phs == 1) {	//	GTGT..
			phs = -1;
			goto Donor;
		    }
		}

#if DEBUG
		if (algmode.nsa & 8) {
		    printf("%2d %2d %2d ", m, n, mx->dir);
		    putvar(mx->val); putvar(y); 
		    putvar(g->val); putvar(eq1->val);
		    if (g2) {
			putvar(g2->val); putvar(eq2->val);
		    }
		    if (algmode.lsg) {
#if INTR > 1
			putvar(hl[0][2][nl[0][2][0]].val); 
			putvar(hl[0][1][nl[0][1][0]].val);
		 	putvar(hl[0][0][nl[0][0][0]].val);
#else
			putvar(hl[0][2].val); 
			putvar(hl[0][1].val);
			putvar(hl[0][0].val);
#endif
			printf(" %2d %2d %6.2f", bb->phs5, bb->phs3, (double) sigE);
			printf(" %6.1f %6.1f", (double) bb->sig5, (double) bb->sig3);
		    }
		    putchar('\n');
		}
#endif
	    }		// end of n loop
	    if (a_in_zone) ++api;	// has exon-exon junctions
	}		// end of m loop

	delete[] (hhc[0] + wdw->lw - 3);
	*scr = clny->val;
	cl->sortcolonies();
	return (cl);
}

VTYPE Aln2h1::diagonalH_ng()
{
	CHAR*	as = a->at(a->left);
	CHAR*	bs = b->at(b->left + 1);
	EXIN*	bb = b->exin->score(b->left + 1);
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

VTYPE Aln2h1::trcbkalignH_ng()
{
	long	ptr;

	vmf = new Vmf();
	VTYPE	scr = forwardH_ng(&ptr);
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
VTYPE Aln2h1::lspH_ng()
{
	WINDOW*	wrsv = wdw;
	WINDOW	wdwl;
	WINDOW	wdwr;
	RANGE	rng[2];			// reserve
	int	aexg = a->inex.exgr;	// reserve
	int	bexg = b->inex.exgr;	// reserve

	if (wdw->up == wdw->lw) return(diagonalH_ng());
	int	m = a->right - a->left;
	int	n = b->right - b->left;
	int	k = wdw->lw - b->left + 3 * a->right;
	int	q = b->right - 3 * a->left - wdw->up;
	long	cvol =  m * n - (k * k + q * q) / 6;
	if (cvol < MaxVmfSpace || m == 1 || n < 3) {
	    return (trcbkalignH_ng());
	}
	save_range(seqs, rng, 2);
	VTYPE	scr = centerH_ng(&m, &q, &n, &k, &wdwl, &wdwr);
	a->inex.exgr = b->inex.exgr = 0;
	int	r = b->left - 3 * a->left;
	if (r < wdwl.lw) b->left = 3 * a->left + wdwl.lw;
	if (r > wdwl.up) a->left = (b->left - wdwl.up) / 3;
	a->right = m;
	b->right = n;
	wdw = &wdwl;
	lspH_ng();
	rest_range(seqs, rng, 2);
	a->inex.exgr = aexg;
	b->inex.exgr = bexg;
	aexg = a->inex.exgl;
	bexg = b->inex.exgl;
	a->inex.exgl = b->inex.exgl = 0;
	r = b->right - 3 * a->right;
	if (r > wdwr.up) b->right = 3 * a->right + wdwr.up;
	if (r < wdwr.lw) a->right = (b->right - wdwr.lw) / 3;
	a->left = q;
	b->left = k;
	wdw = &wdwr;
	lspH_ng();
	wdw = wrsv;
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexg;
	b->inex.exgl = bexg;
	return scr;
}

bool Aln2h1::pincerTrcbkH_ng(int cmode, VTYPE& scr)
{
	vmf = new Vmf();
	long	ptr[2];
	stripe31(seqs, wdw, alprm.sh);
	scr = pincersH_ng(ptr, cmode);
	SKLP    sv;

	if (scr <= pinthr) {
	    delete vmf;
	    return (false);
	}
	sv.p = ptr[0];
	while (sv.p) {
	    vmf->readvmf(&sv, sv.p);
	    mfd->write((UPTR) &sv);
	}
	sv.p = ptr[1];	// restore
	while (sv.p) {
	    vmf->readvmf(&sv, sv.p);
	    mfd->write((UPTR) &sv);
	}
	delete vmf;
	vmf = 0;
	return (true);
}

VTYPE Aln2h1::backforth(int ovr, BOUND& lub)
{
	int	ocodon = ovr / 3;
	VTYPE*	bscr = new VTYPE[ocodon + 1];
	CHAR*	as = a->at(a->left);
	CHAR*	bs = b->at(b->left + 1);
	EXIN*	bb = b->exin->score(b->left + 1);
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

VTYPE Aln2h1::cds5end(int x, int y)
{
	EXIN*	bb = b->exin->score(y + 1);
	VTYPE	scr = 0;
	VTYPE	maxv = 0;
	int	maxy = y;

	JUXT	jxt = {x, y};
	mfd->write(&jxt);
	for ( ; y > b->left; y -= 3) {
	    if (bb->sigS > 0) scr += bb->sigS;
	    if (scr > maxv) {
		maxv = scr;
		maxy = y;
	    }
	    if (bb->sigS > 0 || scr + pwd->Vthr < 0) break;
	    bb -= 3;
	    scr += bb->sigE + pwd->BasicGEP;
	}
	if (maxy != jxt.jy) {
	    jxt.jy = maxy;
	    mfd->write(&jxt);
	}
	return (maxv);
}

VTYPE Aln2h1::cds3end(int x, int y)
{
	EXIN*	bb = b->exin->score(y + 1);
	VTYPE	scr = 0;
	VTYPE	maxv = 0;
	int	maxy = y;

	JUXT	jxt = {x, y};
	mfd->write(&jxt);
	for ( ; y < b->right; y += 3, bb += 3) {
	    if (bb->sigT > 0) scr += bb->sigT;
	    else	scr += bb->sigE + pwd->BasicGEP;
	    if (scr > maxv) {
		maxv = scr;
		maxy = y + 3;
	    }
	    if (bb->sigT > 0 || scr + pwd->Vthr < 0) break;
	}
	if (maxy != jxt.jy) {
	    jxt.jy = maxy;
	    mfd->write(&jxt);
	}
	return (maxv);
}

static void addsigEjxt(JUXT* jxt, int num, Seq* b)
{
	if (!b->exin) return;
	for ( ; num--; ++jxt) {
	    EXIN*	bb = b->exin->score(jxt->jy + 1);
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
	CHAR*	as = a->at(a->left);
	CHAR*	bs = b->at(b->left + 1);
	EXIN*	bb = b->exin->score(b->left + 1);
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
	CHAR*	as = a->at(a->right);
	CHAR*	bs = b->at(b->right + 1);
	EXIN*	bb = b->exin->score(b->right + 1);
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
	int	play = -b->exin->lplay(b->left);
	if (play > agap) agap = play;
	int	d5 = b->left + 3 * agap - 2;	// donor 5' end
	int	d3 = b->left + 2;		// donor 3' end
	int	a5 = b->right - 2;		// accpt 5' end
	int	ilen = a5 - d5;			// intron length
	int	i = 2 - agap;
	VTYPE*	backward = new VTYPE[i];
	backward[--i] = 0;
	CHAR*	as = a->at(a->left - 1);	// cancel score of
	CHAR*	bs = b->at(b->left - 2);	// donor side overlap
	EXIN*	bd = b->exin->score(b->left - 2);
	for (VTYPE v = 0; --i >= 0; --as, bs -= 3, bd -= 3)
	    backward[i] = v += pwd->sim2(as, bs) + bd->sigE;
	int	phs = 1;
	int	m = a->right - 1;
	int	mm = m;
	int	n = d5;
	int	nn = 0;
	as = a->at(m);
	PfqItr	api(a, m);
	int	api_size = api.size();
	bs = b->at(a5);
	bd = b->exin->score(d5);
	EXIN*	ba = b->exin->score(a5);
	VTYPE	maxspj = iscr = NEVSEL;
	VTYPE	ip = pwd->IntPen->Penalty(ilen);
	bool	usespb = api_size && use_spb();
	for (VTYPE v = i = 0; n <= d3; ++n, ++bd, ++ba) {
	    bool	a_in_zone = usespb && api == (3 * m + phs);
	    if (b->exin->isCanon(n, n + ilen)) {
		VTYPE	x = ip + b->exin->sig53(n, n + ilen, IE5P3);
		if (a_in_zone) x += api.match_score(3 * m + phs);
		VTYPE	y = x - v - backward[i];
		if (phs) {
		    CHAR*	cs = spjcs->spjseq(n, n + ilen);
		    if (phs == 1) cs++;
		    y += pwd->pmt->prematT(cs) + pwd->sim2(as, cs);
		}
		if (y > iscr) {nn = n - phs; mm = m; iscr = y; maxspj = x;}
	    }
	    if (++phs == 0) ++i;
	    else if (phs == 2) {++m; phs = -1;}
	    else v += pwd->sim2(++as, bs += 3) + ba[1].sigE;
	    		// cancel score of acc side
	    if (a_in_zone) ++api;
	}
	delete[]	backward;
	if (maxspj > spjthr) {
	    SKL	tmp = {mm, nn};
	    mfd->write(&tmp);
	    tmp.n += ilen;
	    mfd->write(&tmp);
	    return (true);
	} else
	    return (false);
}

VTYPE Aln2h1::seededH_ng(INT level, int eimode, BOUND& lub)
{
	INEX	ainex = a->inex;
	INEX	binex = b->inex;
	RANGE	rng[2];
	int	cmode = eimode;
	bool	qck = algmode.blk && algmode.crs && level == algmode.qck && algmode.lcl < 16;
	int	elmt = algmode.crs? (IntronPrm.elmt + 2) / 3: 1;
	int	agap = a->right - a->left;
	int	bgap = b->right - b->left;
	VTYPE	scr = 0;
	JUXT*	wjxt = 0;
	JUXT*	jxt = 0;
	int	num;
	Wilip*	wl = 0;
	WLUNIT*	wlu = 0;

	if (level == lowestlvl && b->jxt) {
	    jxt = b->jxt;
	    num = b->CdsNo;
	    addsigEjxt(jxt, num, b);
	} else {
	    wl = new Wilip(seqs, pwd, level);
	    wlu = wl->begin();
	    if (wlu) {
		num = wlu->num;
		jxt = wlu->jxt;
	    } else	num = 0;
	}
	save_range(seqs, rng, 2);
	WLPRM*	wlprm = setwlprm(level);
	VTYPE	slmt = pwd->Vthr / 2;
	int	term = (int) (pwd->Vthr / pwd->BasicGEP);
	int	backstep = wlprm->width;
	if (++level < algmode.qck) {
	    wlprm = setwlprm(level);
	    backstep = wlprm->width;
	}
	int	wlmt = (int) (3 * wlprm->tpl);
	BOUND	bab = lub;
	if (num) {
	    jxt[num].jx = a->right;
	    jxt[num].jy = b->right;
	    a->inex.exgr = 0;
	    b->inex.exgr = 0;
	    agap = jxt->jx - a->left;
	    bgap = jxt->jy - b->left;
	    for (wjxt = jxt; num--; ++wjxt) {
		a->right = wjxt->jx;
		b->right = wjxt->jy;
		bab.ua = wjxt->jx + wjxt->jlen;
		if (wjxt[1].jx < bab.ua) bab.ua = wjxt[1].jx;
		bab.ua -= backstep;
		int	lx  = wjxt->jx + wjxt->jlen / 2;
		if (bab.ua < lx) bab.ua = lx;
		bab.ub = wjxt->jy + 3 * (bab.ua - wjxt->jx);
		int	ovr = min(3 * agap, bgap);
		VTYPE	iscore;
		scr += wjxt->jscr;
		if (cmode == 1 && agap == 0) scr += cds5end(wjxt->jx, wjxt->jy);
		else if (cmode == 3 && agap <= 0 && bgap >= IntronPrm.llmt &&
			indelfreespjH(agap, iscore))
		    scr += iscore;
		else if (level < algmode.qck && ovr >= wlmt)
		    scr += seededH_ng(level, cmode, bab);	// recursive search
		else if (ovr <= 0 && bgap < IntronPrm.llmt)
		    scr += backforth(-ovr, bab);
		else {
		    scr -= creepback(ovr, slmt, bab);
		    scr -= creepfwrd(ovr, slmt, bab);
		    if (ovr < 0) {	// excessive overlap skip this hsp
			agap = wjxt[1].jx - a->left;
			bgap = wjxt[1].jy - b->left;
			continue;
		    }
		    float	dpspace = (float) agap * (float) bgap / MEGA;
		    if (cmode == 1 && (qck || dpspace >= alprm.maxsp)) {
			int	bl = wjxt->jy + term;	// 5' end
			if (bl > b->left) b->left = bl;
			pincerTrcbkH_ng(cmode, iscore);
			scr += iscore;
		    } else if (cmode == 3 && agap < elmt && bgap >= IntronPrm.llmt
			&& pincerTrcbkH_ng(cmode, iscore)) {	// Sandwitch
			scr += iscore;
		    } else if (dpspace < alprm.maxsp) {
			stripe31(seqs, wdw, alprm.sh);
			if (wdw->up == wdw->lw) diagonalH_ng();
			else	trcbkalignH_ng();
		    } else { 		// give up alignment
//			scr += pincerTrcbkH_ng(cmode, iscore);
			mfd->write((UPTR) wjxt);
			mfd->write((UPTR) (wjxt + 1));
		    }
		}
		cmode = 3;
		a->left = wjxt->jx + wjxt->jlen;
		b->left = wjxt->jy + 3 * wjxt->jlen;
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

	int	ovr = min(3 * agap, bgap);
	VTYPE	iscore;
	if (cmode == 2 && agap == 0)		// 3' end
	    scr += cds3end(wjxt->jx, wjxt->jy - bgap);
	else if (cmode == 3 && agap <= 0 && bgap >= IntronPrm.llmt &&
	    indelfreespjH(agap, iscore))
	    scr += iscore;
	else if (level < algmode.qck && ovr >= wlmt)
	    scr += seededH_ng(level, cmode, bab);	// recursive search
	else if (ovr <= 0 && bgap < IntronPrm.llmt)
	    scr += backforth(-ovr, bab);
	else {
	    scr -= creepback(ovr, slmt, bab);
	    scr -= creepfwrd(ovr, slmt, bab);
	    float	dpspace = (float) agap * (float) bgap / MEGA;
	    if (cmode == 1 && (qck || dpspace >= alprm.maxsp)) {
		pincerTrcbkH_ng(cmode, iscore);
		scr += iscore;
	    } else if (cmode == 2 && (qck || dpspace >= alprm.maxsp)) {
		if (wjxt) {
		    int	br = bgap - wjxt->jy - term;	// near 3' end
		    if (br > b->left && br < b->right) b->right = br;
		}
		pincerTrcbkH_ng(cmode, iscore);
		scr += iscore;
	    } else if (cmode == 3 && agap < elmt && bgap >= IntronPrm.llmt
		&& pincerTrcbkH_ng(cmode, iscore)) {	// Sandwitch
		scr += iscore;
	    } else if (dpspace < alprm.maxsp) {
		stripe31(seqs, wdw, alprm.sh);
		if (wdw->up == wdw->lw) diagonalH_ng();
		else	trcbkalignH_ng();
	    } else { 	// cmode == 3, give up alignment
		SKL	tmp = {a->left, b->left};
		mfd->write((UPTR) &tmp);
		tmp.m = a->right; tmp.n = b->right;
		mfd->write((UPTR) &tmp);
	    }
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

SKL* Aln2h1::globalH_ng(VTYPE* scr)
{
	mfd = new Mfile(sizeof(SKL));
	SKL	wsk;
	mfd->write((UPTR) &wsk);	// dummy call
	if (algmode.qck) {
	    BOUND bab = {a->left, b->left, a->right, b->right};
	    *scr = seededH_ng(lowestlvl, 1, bab);
	} else
	    *scr = lspH_ng();
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

VTYPE HomScoreH_ng(Seq* seqs[], PwdB* pwd)
{
	WINDOW	wdw;

	stripe31(seqs, &wdw, alprm.sh);
	Aln2h1 alnv(seqs, pwd, &wdw);
	return alnv.forwardH_ng(0);
}

Colonies* swg1stH_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr)
{
	WINDOW  wdw;

	stripe31(seqs, &wdw, alprm.sh); // Recacl. window boundaries
	Aln2h1	alnv(seqs, pwd, &wdw);
	return alnv.fwdswgH_ng(scr);
}

SKL* swg2ndH_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr, const COLONY* clny)
{
	if (clny->val <= 0) return (0);
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	a->left = clny->mlb;
	b->left = clny->nlb;
	b->right = clny->nrb;
	a->right = clny->mrb;
	WINDOW	wdw = {clny->upr, clny->lwr, clny->upr - clny->lwr + 7};
	Aln2h1 alnv(seqs, pwd, &wdw);
	return (alnv.globalH_ng(scr));
}

SKL* alignH_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr)
{
	WINDOW  wdw;

	stripe31(seqs, &wdw, alprm.sh); // Recacl. window boundaries
	Aln2h1 alnv(seqs, pwd, &wdw);
	return (alnv.globalH_ng(scr));
}
