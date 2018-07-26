/*****************************************************************************
*
*	Alignment of genomic vs. cDNA/mRNA sequences.
*	5' and 3' splice site signals are considered.
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
*	Osamu Gotoh, Ph.D.	(2003-2012)
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

#define	INTR	2
#define	NCANDS	(INTR * NOL)
#define MAXR	4
#define	CigarM	0

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

typedef LRECD	EIJUNC[NCANDS+1];
typedef int	EIJODR[NCANDS+1];

struct BOUND {int la, lb, ua, ub;};

class Aln2s1 {
protected:
	Seq**	seqs;
	Seq*&	a;
	Seq*&	b;
	WINDOW*	wdw;
	PwdB*	pwd;
	Mfile*	mfd;
	Vmf*	vmf;
	Exinon*	exin;
	INT	lowestlvl;
public:
	Aln2s1(Seq** _seqs, PwdB* _pwd, WINDOW* _wdw);
	~Aln2s1() {delete mfd;}
	void	initS_ng(RECD* hf[]);
	RECD*	lastS_ng(RECD* hhb[]);
	VTYPE	forwardS_ng(long pp[]);
	void	finitS_ng(WRECD* hhf[], int mm);
	void	binitS_ng(WRECD* hhf[], int mm);
	VTYPE	centerS_ng(int* ml, int* mr, int* nl, int* nr, 
		WINDOW* wdwf, WINDOW* wdwb);
	void	cinitS_ng(CRECD* hhc[]);
	void	pfinitS_ng(LRECD* hhg[]);
	void	pbinitS_ng(LRECD* hhg[]);
	VTYPE	pincerTrcbkS_ng(int cmode);
	VTYPE	pincersS_ng(long* ptr, int cmode);
	Colonies*	fwdswgS_ng(VTYPE* scr);
	VTYPE	diagonalS_ng();
	VTYPE	trcbkalignS_ng();
	VTYPE	backforth(int n, BOUND& lub);
	VTYPE	lspS_ng();
	VTYPE	seededS_ng(INT level, int cmode, BOUND& lub);
	SKL*	globalS_ng(VTYPE* scr);
	VTYPE	creepback(int ovr, VTYPE bscr, BOUND& lub);
	VTYPE	creepfwrd(int& ovr, VTYPE bscr, BOUND& lub);
	bool	indelfreespjS(int agap, VTYPE& iscr);
};

extern	Colonies* swg1stS_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr);
extern	VTYPE	HomScoreS_ng(Seq* seqs[], PwdB* pwd, long ptr[]);
extern	SKL*	swg2ndS_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr, COLONY* clny);
extern	SKL*	alignS_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr);

static	RECD	oom = {NEVSEL, 0, 0, 0};
static	LRECD	ooml = {NEVSEL, 0, 0, 0, 0};
static	WRECD	oomW = {NEVSEL, 0, INT_MIN, INT_MAX, 0, 0};
static	CRECD	oomC = {NEVSEL, 0, INT_MIN, INT_MAX, 0, 0, 0, 0};

Aln2s1::Aln2s1(Seq* _seqs[], PwdB* _pwd, WINDOW* _wdw) :
	seqs(_seqs), a(seqs[0]), b(seqs[1]), wdw(_wdw), pwd(_pwd),
	mfd(0), vmf(0), exin(b->exin), lowestlvl(b->wllvl)
{ }

void Aln2s1::initS_ng(RECD* hh[])
{
	int	n = b->left;
	int	r = b->left - a->left;
	int	rr = b->right - a->left;
	RECD*	hf[NOL];

	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hh[k] + r;
	if (wdw->up < rr) rr = wdw->up;
	for (int i = 0; r <= wdw->up + 1; ++r, ++n, ++i) {
	    RECD*	h = hf[0]++;
 	    for (int k = 1; k < pwd->Nrow; ++k) *hf[k]++ = oom;
	    if (i == 0) {
		h->val = 0;
	    } else if (!a->inex.exgl || r > rr) {
		*h = oom;
		continue;
	    } else {
		*h = h[-1];
		int	k = n - h->jnc;
		if (k == 1) h->val += pwd->BasicGOP;
		h->val += pwd->GapExtPen(k);
		h->dir = HORI;
	    }
	    VTYPE	x = 0;
	    if (h->val <= x) {
		h->val = x;
		h->dir = i? DEAD: DIAG;
		h->jnc = n;
		h->ptr = vmf? vmf->add(a->left, n, 0): 0;
	    }
	}

	r = b->left - a->left;
	rr = b->left - a->right;
	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hh[k] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	for (int i = 1; --r >= wdw->lw - 1; ++i) {
	    RECD*	h = --hf[0];
	    RECD*	g = --hf[1];
	    for (int k = 2; k < pwd->Nrow; ++k) *--hf[k] = oom;
	    if (r < rr) *h = *g = oom;
	    else if (b->inex.exgl) {	/* local */
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->left;
		h->ptr = vmf? vmf->add(a->left + i, b->left, 0): 0;
		*g = oom;
	    } else {			/* (semi) global */
		*h = h[1];
		if (i == 1) {
		    h->val += pwd->GapPenalty(1);
		    h->dir = VERT;
		    *g = *h;
		} else {
		    *g = g[1];
		    h->val += pwd->GapExtPen(i);
		    g->val += pwd->BasicGEP;
		}
	    }
	}
}

RECD* Aln2s1::lastS_ng(RECD* hh[])
{
	RECD*	h9 = hh[0] + b->right - a->right;
	RECD*	mx = h9;

	if (b->inex.exgr) {
	    int	rw = wdw->up;
	    int	rf = b->right - a->left;
	    if (rf < rw) rw = rf;
	    for (RECD* h = hh[0] + rw; h > h9; --h)
		if (h->val > mx->val) mx = h;
	}
	int	rw = wdw->lw;
	int	rf = b->left - a->right;
	if (rf > rw) rw = rf;
	else	 rf = rw;
	RECD*	h = hh[0] + rw;
	if (a->inex.exgr)
	    for (int i = 0; h <= h9; ++h, ++i)
		if (h->val > mx->val) mx = h;
	int	i = mx - h9;
	rf = a->right;	/* m9 */
	rw = b->right;	/* n9 */
	if (i > 0) rf -= i;
	if (i < 0) rw += i;
	if (vmf || a->inex.exgr) mx->ptr = vmf->add(rf, rw, mx->ptr);
	return (mx);
}

VTYPE Aln2s1::forwardS_ng(long* pp)
{
	RECD*	hh[NOL];	/* H matrix	 	*/
	RECD*	hf[NOL];
	RECD	f1, f2;
	RECD*	g2 = 0;
#if INTR > 1
	RECD	hl[NOL][INTR+1]; /* [DIA, HORI, HORL][candidates] */
	int	nl[NOL][INTR+1]; /* [DIA, HORI, HORL][candidates] */
#else
	RECD    hl[NOL];
#endif
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

	hh[0] = new RECD[pwd->Nrow * wdw->width] - wdw->lw + 1;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + wdw->width;
	if (vmf) vmf->add(0, 0, 0);	/* Skip 0-th record */
	initS_ng(hh);

	int	m = a->left;
	PfqItr	api(a, m);
	int	api_size = api.size();
	if (!a->inex.exgl) --m; /* global */
	int	n1 = m + wdw->lw - 1;
	int	n2 = m + wdw->up;
	CHAR*	as = a->at(m);
	hf[1] = &f1;
	hf[2] = &f2;
	for ( ; ++m <= a->right; ++as) {
	    bool	internal = !a->inex.exgr || m < a->right;
	    n1++; n2++;
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    int k = 0;
	    RECD*	h = hh[k] + r;
	    RECD*	g = hh[++k] + r;
	    f1 = oom;
	    if (pwd->Noll == 3) {
		g2 = hh[++k] + r;
		f2 = oom;
	    }
	    for (int k = 0; k < pwd->Noll; ++k) {
#if INTR > 1
		RECD*	phl = hl[k];
		for (int l = 0; l <= INTR; ++l, ++phl) {
		    nl[k][l] = l;
		    *phl = oom;
	  	}
#else
		hl[k] = oom;
#endif
	    }
#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    bool	a_in_zone = api_size && (api == m);
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x;
		++g;
		hf[0] = ++h;
		if (g2) ++g2;
/*	Diagonal	*/
		RECD*	from = h;
		RECD*	mx = h;
		VTYPE	diag = h->val;
		if (m == a->left) goto HorizonS;
		h->val += pwd->sim2(as, bs);
		h->dir = (from->dir & DIAG)? DIAG: NEWD;

/*	Vertical	*/
		x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		    g->ptr = (vmf && from->dir & SPIN)?
			vmf->add(m - 1, n, from->ptr): from->ptr;
		    g->jnc = from->jnc;
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val > mx->val) mx = g;

/*	Vertical2	*/
		if (g2) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[1].val) {
		    g2->val = x;
		    g2->dir = VERL | (g->dir & SPJC);
		    g2->ptr = (vmf && from->dir & SPIN)? 
			vmf->add(m - 1, n, from->ptr): from->ptr;
		    g2->jnc = from->jnc;
		  } else	*g2 = g2[1];
		  g2->val += pwd->LongGEP;
		  if (g2->val > mx->val) mx = g2;
		}
HorizonS:
/*	Horizontal	*/
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f1.val) {
		    f1 = h[-1];
		    f1.val = x;
		}
		if (f1.dir) {
		    f1.val += pwd->BasicGEP;
		    f1.dir = (f1.dir & SPIN) + HORI;
		    if (f1.val >= mx->val) mx = &f1;
		}

/*	Horizontal2	*/
		if (g2) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= f2.val) {
			f2 = h[-1];
			f2.val = x;
		    }
		    if (f2.dir) {
			f2.val += pwd->LongGEP;
			f2.dir = (f2.dir & SPIN) + HORL;
			if (f2.val >= mx->val) mx = &f2;
		    }
		}

/*	intron 3' boundary, assume no overlapping signals	*/
		if (exin->isAccpt(n) && internal) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    for (int k = 0; k < pwd->Noll; ++k) {
			from = hf[k];
#if INTR > 1
			RECD*	maxphl = 0;
			int*	pnl = nl[k];
			for (int l = 0; l < INTR; ++l) {
			    RECD*	phl = hl[k] + pnl[l];
			    if (!phl->dir) break;
			    x = phl->val + sigJ +
				pwd->IntPen->Penalty(n - phl->jnc) +
				exin->sig53(phl->jnc, n, IE53);
			    if (x > from->val) {
				from->val = x;
				maxphl = phl;
			    }
			}
			if (maxphl) {
			    if (vmf) from->ptr = vmf->add(m, n, 
				vmf->add(m, maxphl->jnc, maxphl->ptr));
			    from->jnc = n;
			    from->dir = maxphl->dir | SPJCI;
			    if (from->val > mx->val) mx = from;
			}
#else
			RECD*	phl = hl + k;
			if (!phl->dir) break;
			x = phl->val + pwd->IntPen->Penalty(n - phl->jnc)
			    + exin->sig53(phl->jnc, n, IE53);
			if (x > from->val) {
			    from->val = x;
			    if (vmf) from->ptr = vmf->add(m, n, 
				vmf->add(m, phl->jnc, phl->ptr));
			    from->jnc = n;
			    from->dir = phl->dir | SPJCI;
			    if (from->val > mx->val) mx = from;
			}
#endif
		    }
		}

/*	Find optimal path	*/
		if (h != mx) *h = *mx;	/* non-diagonal */
		else if (Local && h->val > diag) {
		    if (LocalL && diag == 0 && !(h->dir & SPJC))
			h->ptr = vmf? vmf->add(m - 1, n - 1, 0): 0;
		    else if (LocalR && h->val > maxh.val) {
			maxh.val = h->val;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && h->val <= 0) h->val = h->dir = 0;
		else if (vmf && h->dir == NEWD)
		    h->ptr = vmf->add(m - 1, n - 1, h->ptr);


/*	intron 5' boundary	*/
		if (exin->isDonor(n) && internal) {
		    VTYPE	sigJ = exin->sig53(n, 0, IE5);
		    for (int k = 0; k < pwd->Noll; ++k) {
				/* An orphan exon is disallowed */
				/* Vertical path is not transferred */
			from = hf[k];
			if (!from->dir || (from->dir & SPIN) || 
				from->dir == VERT || from->dir == VERL)
			    continue;
			x = from->val + sigJ;
#if INTR > 1
			RECD*	phl = hl[k];
			int*	pnl = nl[k];
			int	l = INTR;
			for ( ; --l >= 0; ) {
			    if (x > phl[pnl[l]].val)
				gswap(pnl[l], pnl[l + 1]);
			    else
				break;
			}
			if (++l < INTR) {
			    phl[pnl[INTR]].dir = 0;
			    phl += pnl[l];
			    phl->val = x;
			    phl->jnc = n;
			    phl->dir = from->dir;
			    phl->ptr = from->ptr;
			}
#else
			RECD*	phl = hl + k;
			if (x > phl->val) {
			    phl->val = x;
			    phl->jnc = n;
			    phl->dir = from->dir;
			    phl->ptr = from->ptr;
			}
#endif
		    }
		}

#if DEBUG
		if (algmode.nsa & 8) {
		    printf("%2d %2d %2d ", m, n, mx->dir);
		    putvar(mx->val); putvar(y); 
		    putvar(g->val); putvar(f1.val);
		    if (g2) {
			putvar(g2->val); putvar(f2.val);
		    }
		    if (algmode.lsg) {
#if INTR > 1
			putvar(hl[0][nl[0][0]].val);
			putvar(hl[1][nl[1][0]].val);
			if (g2) putvar(hl[2][nl[2][0]].val);
#else
			putvar(hl[0].val);
			putvar(hl[1].val);
			if (g2) putvar(hl[2].val);
#endif
		    }
		    putchar('\n');
		}
#endif
	    }
	    if (a_in_zone) ++api;		// has exon-exon junctions
	}

	if (LocalR) {
	    if (vmf) *pp = vmf->add(maxh.m, maxh.n, maxh.p);
	} else {
	    RECD*	mx = lastS_ng(hh);
	    maxh.val = mx->val;
	    if (pp) *pp = mx->ptr;
	}
	delete[] (hh[0] + wdw->lw - 1);
	return (maxh.val);
}

VTYPE skl_rngS_ng(Seq* seqs[], Gsinfo* gsi, PwdB* pwd)
{
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];
	Exinon*	exin = b->exin;
	SKL*	wsk = gsi->skl;
	int	num = (wsk++)->n;
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
	Eijnc*	eijnc = gsi->eijnc = new Eijnc(true);
	Cigar*	cigar = gsi->cigar = 0;
	Vulgar*	vlgar = gsi->vlgar = 0;
	Samfmt*	samfm = gsi->samfm = 0;
	switch (algmode.nsa) {
	    case EXN_FORM: 
	    case CIG_FORM: cigar = gsi->cigar = new Cigar(); break;
	    case VLG_FORM: gsi->vlgar = vlgar = new Vulgar(); break;
	    case SAM_FORM: gsi->samfm = samfm = new Samfmt(); break;
	    default:	 break;
	}
	vclear(fst);
	vclear(&pst);
	vclear(&rbuf);
	if (wsk[1].n == wsk->n && b->inex.exgl) {++wsk; --num;}
	int	m = wsk->m;
	int	n = wsk->n;
	CHAR*	as = a->at(m);
	CHAR*	bs = b->at(n);
	PfqItr	api(a, m);
	bool	usespb = api.size() && use_spb();
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
	    int	mi = wsk->m - m;
	    if (mi && insert) {
		bool	j = a->inex.exgl && m == a->left;
		VTYPE	x = j? 0: pwd->GapPenalty(insert);
		VTYPE	xi = NEVSEL;
		if (intlen) {
		    insert -= intlen;
		    xi = rbuf.iscr + pwd->GapPenalty(insert);
		}
		if (xi >= x) {  /* intron */
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
		    h += xi;
		    insert -= preint;
		} else {	/* gap */
		    h += x;
		}
		if (insert) {
		    fst->gap += 1;
		    fst->unp += insert;
		    if (cigar) cigar->push('D', insert);
		    if (samfm) samfm->push('D', insert);
		    if (vlgar) vlgar->push('G', 0, insert);
		}
		insert = intlen = preint = 0;
		rbuf.iscr = NEVSEL;
	    }
	    int	ni = wsk->n - n;
	    if (ni && deletn && mi) {
		if (!(b->inex.exgl && n == b->left)) {
		    h += pwd->GapPenalty(deletn);
		    fst->gap += 1;
		    fst->unp += deletn;
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
		x += pwd->sim2(as, bs);
		if (*as == *bs) ++fst->mch;
		else ++fst->mmc;
		if (eijnc) eijnc->shift(rbuf, *fst, pst,
		    n - rbuf.left == alprm2.jneibr);
		}
#else
		int	run = 0;
		for ( ; d; --d, ++as, ++bs, ++n) {
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
		if (eijnc) eijnc->shift(rbuf, *fst, pst,
		    n - rbuf.left == alprm2.jneibr);
		}
		if (samfm) {
		if (run > 0) samfm->push('=', run); else
		if (run < 0) samfm->push('X', -run);
		}
#endif
		h += x;
	    }
	    if (i > 0) {
		deletn += i;
		if (cigar) cigar->push('I', i);
		if (samfm) samfm->push('I', i);
		if (vlgar) vlgar->push('G', i, 0);
	    } else if (i < 0) {
		int	n3 = n + (i = -i);
		if (algmode.lsg && i > IntronPrm.llmt) {
		    sig5 = exin->sig53(n, n3, IE5);
		    sig3 = exin->sig53(n, n3, IE3);
		    xi = exin->sig53(n, n3, IE5P3) + pwd->IntPen->Penalty(i);
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
			 n - rbuf.left < alprm2.jneibr);
			pst = *fst;
		    }
		} else if (!a->inex.exgl || m != a->left) {
		    if (!insert) fst->gap += 1;
		    for (int j = 0; j < i; ++j, ++n) {
		    	++fst->unp;
			if (eijnc) eijnc->shift(rbuf, *fst, pst,
			    n - rbuf.left == alprm2.jneibr);
		    }
		}
		bs += i;
		insert += i;
    	    }
	    m = wsk->m;
	    n = wsk->n;
	    if (usespb) while (api < m) ++api;
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
	    eijnc->unshift();
	    eijnc->store(rbuf, *fst, pst, n - rbuf.left <= alprm2.jneibr);
	    eijnc->push(&rbuf);
	    rbuf.left = endrng.left;
	    rbuf.right = endrng.right;
	    eijnc->push(&rbuf);
	    eijnc->flush();
	    gsi->noeij = eijnc->size() - 1;
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
	fst->val = (FTYPE) h;
	return (h);
}

void Aln2s1::finitS_ng(WRECD* hhf[], int mm)
{
	int	n = b->left;
	int	r = b->left - a->left;
	int	r0 = r;
	int	rr = b->right - a->left;
	WRECD*	hf[NOL];

	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hhf[k] + r;
	if (wdw->up < rr) rr = wdw->up;
	for (int i = 0; r <= wdw->up + 1; ++r, ++n, ++i) {
	    WRECD*	h = hf[0]++;
 	    for (int k = 1; k < pwd->Nrow; ++k) *hf[k]++ = oomW;
	    if (i == 0) {
		h->val = 0;
	    } else if (!a->inex.exgl || r > rr) {
		*h = oomW;
		continue;
	    } else {
		*h = h[-1];
		int	k = n - h->jnc;
		if (k == 1) h->val += pwd->BasicGOP;
		h->val += pwd->GapExtPen(k);
		h->dir = HORI;
	    }
	    h->upr = h->lst = r;
	    if (h->val <= 0) {
		h->val = 0;
		h->dir = i? DEAD: DIAG;
		h->jnc = n;
		h->lwr = r;
	    }
	}

	r = r0;
	rr = b->left - mm;
	for (int k = 0; k < pwd->Nrow; ++k) hf[k] = hhf[k] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	for (int i = 1; --r >= wdw->lw - 1; ++i) {
	    WRECD*	h = --hf[0];
	    WRECD*	g = --hf[1];
	    for (int k = 2; k < pwd->Nrow; ++k) *--hf[k] = oomW;
	    if (r < rr) *h = *g = oomW;
	    else if (b->inex.exgl) {
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->left;
		h->lwr = h->upr = h->lst = r;
		*g = oomW;
	    } else {
		*h = h[1];
		if (i == 1) {
		    h->val += pwd->GapPenalty(1);
		    h->dir = VERT;
		    *g = *h;
		} else {
		    *g = g[1];
		    h->val += pwd->GapExtPen(i);
		    g->val += pwd->BasicGEP;
		}
		h->lwr = g->lwr = r;
	    }
	}
}

void Aln2s1::binitS_ng(WRECD* hhb[], int mm)
{
	int	n = b->right;
	int	r = b->right - a->right;
	int	r9 = r;
	int	rr = b->left - a->right;
	WRECD*	hb[2 * NOL];

	for (int k = 0; k < pwd->Nrow; ++k) hb[k] = hhb[k] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	for (int i = 0; r >= wdw->lw - 1; --r, --n, ++i) {
	    WRECD*	h = hb[0]--;
	    for (int k = 1; k < pwd->Nrow; ++k) *hb[k]-- = oomW;
	    if (i == 0) {
		h->val = 0;
	    } else if (!a->inex.exgr || r < rr) {	/* global */
		*h = oomW;
		continue;
	    } else {
		*h = h[1];
		int	k = h->jnc - n;
		if (k == 1) h->val += pwd->BasicGOP;
		h->val += pwd->GapExtPen(k);
		h->dir = HORI;
	    }
	    h->lwr = h->lst = r;
	    if (h->val <= 0) {
		h->val = 0;
		h->dir = i? DEAD: DIAG;
		h->jnc = n;
		h->upr = r;
	    }
	}

	r = r9;
	rr = b->right - mm;
	if (wdw->up < rr) rr = wdw->up;
	for (int k = 0; k < pwd->Nrow; ++k) hb[k] = hhb[k] + r;
	for (int i = 1; ++r <= wdw->up + 1; ++i) {
	    WRECD*	h = ++hb[0];
	    WRECD*	g = ++hb[1];
	    for (int k = 2; k < pwd->Nrow; ++k) *++hb[k] = oomW;
	    if (r > rr) *h = *g = oomW;
	    else if (b->inex.exgr) {
		h->val = 0;
		h->dir = DEAD;
		h->jnc = b->right;
		h->lwr = h->upr = h->lst = r;
		*g = oomW;
	    } else {
		*h = h[-1];
		if (i == 1) {
		    h->val = pwd->GapPenalty(1);
		    h->dir = VERT;
		    *g = *h;
		} else {
		    *g = g[-1];
		    h->val += pwd->GapExtPen(i);
		    g->val += pwd->BasicGEP;
		}
		h->upr = g->upr = r;
	    }
	}
}

VTYPE Aln2s1::centerS_ng(int* ml, int* mr, int* nl, int* nr, 
		WINDOW* wdwf, WINDOW* wdwb)
{
	int	kk = 0;
	int	jtigh = 0;
	int	kb = 0;
	int	jb = 0;
	int	rr[MAXR+1];
#if INTR > 1
	WRECD	hl[NOL][INTR+1];
	int	nx[NOL][INTR+1]; /* [DIA, HORI, HORL][candidates] */
#else
	WRECD   hl[NOL];
#endif
	WRECD*	hhf[NOL];
	WRECD*	hhb[2 * NOL];
	WRECD*	hf[NOL];
	WRECD*	hb[2 * NOL];
	WRECD	f[NOL];
	WRECD*	f1 = f + 1;
	WRECD*	f2 = 0;
	WRECD*	b1 = f1;
	WRECD*	b2 = f2;
	WRECD*	g2 = 0;
	VTYPE	mxh = NEVSEL;
	SKL	wskl;
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
#if MONITOR
	long	start = time(0);
	long	stop;
#endif

	hhf[0] = new WRECD[pwd->Nrow * wdw->width] - wdw->lw + 1;
	hhb[0] = new WRECD[pwd->Nrwb * wdw->width] - wdw->lw + 1;
	for (int k = 1; k < pwd->Nrow; ++k) hhf[k] = hhf[k-1] + wdw->width;
	for (int k = 1; k < pwd->Nrwb; ++k) hhb[k] = hhb[k-1] + wdw->width;
	int	mm = (a->left + a->right + 1) / 2;
	if (pwd->Noll == 3) b2 = f2 = f + 2;
	for (int k = 1; k < pwd->Noll; ++k) hf[k] = hb[k] = f + k;

/*******	backward phase	*******/

	binitS_ng(hhb, mm);
	int	m = a->right;
	PfqItr	api(a, m - 1);
	int	api_size = api.size();
	if (!a->inex.exgr) ++m; /* global */
	CHAR*	as = a->at(m);
	int	n1 = m + wdw->lw;
	int	n2 = m + wdw->up + 1;
	for ( ; --m >= mm; ) {
	    --as; --n1; --n2;
	    int	n = min(n2, b->right);
	    int	n9 = max(n1, b->left);
	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    int	d = 0;
	    WRECD*	h = hhb[d] + r;
	    WRECD*	g = hhb[++d] + r;
	    *b1 = oomW;
	    if (pwd->Noll == 3) {
		g2 = hhb[++d] + r;
		*b2 = oomW;
	    }
	    if (m == mm) 
		while (d++ < 2 * pwd->Noll - 1) hb[d] = hhb[d] + r;
	    for (int k = 0; k < pwd->Noll; ++k) {
#if INTR > 1
		WRECD*	phl = hl[k];
		for (int l = 0; l <= INTR; ++l, ++phl) {
		    nx[k][l] = l;
		    *phl = oomW;
		}
#else
		hl[k] = oomW;
#endif
	    }

#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m+1, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    bool	a_in_zone = api_size && (api == m);
	    for ( ; --n >= n9; ) {
		VTYPE	x, y;
		--bs; --g; --r;
		hb[0] = --h;
		if (g2) --g2;

/*	Diagonal	*/
		WRECD*	from = h;
		WRECD*	mx = h;
		if (m == a->right) goto HorizonB;
		h->val += pwd->sim2(as, bs);
		h->dir = (from->dir & DIAG)? DIAG: NEWD;

/*	Vertical	*/
		x = (--from)->val + pwd->BasicGOP;
		if (x >= g[-1].val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		} else *g = g[-1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;

/*	Vertical2	*/
		if (g2) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[-1].val) {
		    *g2 = *from;
		    g2->val = x;
		    g2->dir = VERL | (g2->dir & SPJC);
		  } else *g2 = g2[-1];
		  g2->val += pwd->LongGEP;
		  if (g2->val >= mx->val) mx = g2;
		}
HorizonB:
/*	Horizontal	*/
		x = h[1].val + pwd->BasicGOP;
		if (x >= b1->val) {
		    *b1 = h[1];
		    b1->val = x;
		}
		if (b1->dir ) {
		    b1->val += pwd->BasicGEP;
		    b1->dir = (b1->dir & SPIN) + HORI;
		    if (b1->val >= mx->val) mx = b1;
		}

/*	Horizontal2	*/
		if (b2) {
		    x = h[1].val + pwd->LongGOP;
		    if (x >= b2->val) {
			*b2 = h[1];
			b2->val = x;
		    }
		    if (b2->dir) {
			b2->val += pwd->LongGEP;
			b2->dir = (b2->dir & SPIN) + HORI;
			if (b2->val >= mx->val) mx = b2;
		    }
		}
		if (mx->dir == NEWD) mx->lst = r;

/*	intron 5' boundary reverse phase */
		if (exin->isDonor(n)) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = hb[d];
#if INTR > 1
			WRECD*	maxphl = 0;
			int*	pnx = nx[d];
			for (int l = 0; l < INTR; ++l) {
			    WRECD*	phl = hl[d] + pnx[l];
			    if (!phl->dir) break;
/*			    x = exin->sig53(n, phl->jnc, IE35) + sigJ +	*/
			    x = exin->sig53(n, phl->jnc, IE53) + sigJ +
				(exin->sig53(n, 0, IE5) - phl->sig) +
				phl->val + pwd->IntPen->Penalty(phl->jnc - n);
			    if (x > from->val) {
				from->val = x;
				maxphl = phl;
		    	    }
			}
			if (maxphl) {
	 		    from->jnc = maxphl->jnc;
			    from->dir = maxphl->dir | SPJCI;
			    from->upr = max(maxphl->upr, r);
			    from->lwr = min(maxphl->lwr, r);
			    from->lst = maxphl->lst + n - maxphl->jnc;
			    if (from->val > mx->val) mx = from;
			}
#else
			WRECD*	phl = hl + d;
			if (!phl->dir) break;
/*			x = exin->sig53(n, phl->jnc, IE35) +	*/
			x = exin->sig53(n, phl->jnc, IE53) + sigJ +
			    (exin->sig53(n, 0, IE5) - phl->sig) +
			    phl->val + pwd->IntPen->Penalty(phl->jnc - n);
			if (x > from->val) {
			    from->val = x;
	 		    from->jnc = phl->jnc;
			    from->dir = phl->dir | SPJCI;
			    from->upr = max(phl->upr, r);
			    from->lwr = min(phl->lwr, r);
			    from->lst = phl->lst + n - phl->jnc;
			    if (from->val > mx->val) mx = from;
			}
#endif
		    }
		}

/*	Find optimal path	*/
		y = h->val;
		if (h != mx) {
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		}
		if (LocalR && h->val <= 0) h->val = h->dir = 0;

/*	intron 3' boundary reverse phase	*/

		if (exin->isAccpt(n)) {
		    VTYPE	sigJ = exin->sig53(0, n, IE3);
		    for (int d = 0; d < pwd->Noll; ++d) {
				/* An orphan exon is disallowed */
			from = hb[d];
			if (!from->dir || from->dir & SPIN) continue;
			x = from->val + sigJ;
#if INTR > 1
			WRECD*	phl = hl[d];
			int*	pnx = nx[d];
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
			    phl->jnc = n;
			    phl->sig = sigJ;
			}
#else
			WRECD*	phl = hl + d;
			if (x > phl->val) {
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = n;
			    phl->sig = sigJ;
			}
#endif
		    }
		}

#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(b1->val);
		if (g2) {
		    putvar(g2->val); putvar(b2->val);
		}
		if (algmode.lsg) {
#if INTR > 1
		    putvar(hl[0][nx[0][0]].val);
		    putvar(hl[1][nx[1][0]].val);
		    if (g2) putvar(hl[2][nx[2][0]].val);
#else
		    putvar(hl[0].val);
		    putvar(hl[1].val);
		    if (g2) putvar(hl[2].val);
#endif
		}
		putchar('\n');
	}
#endif

/*	Prepare for finding the center */

		if (m == mm) {
		    if (isvert(mx) && g->val < mx->val)
			*g = *mx;
		    *--hb[d = pwd->Noll] = *b1;
		    if (b2) *--hb[++d] = *b2;
		}

	    }	/* end of n-loop */
	    if (a_in_zone) --api;		// has exon-exon junctions
	}	/* end of m-loop */

/*******	forward phase	*******/

	finitS_ng(hhf, mm);
	m = a->left;
	if (api_size) api.reset(m);
	if (!a->inex.exgl) --m; /* global */
	n1 = m + wdw->lw - 1;
	n2 = m + wdw->up;
	as = a->at(m);
	for ( ; ++m <= mm; ++as) {
	    ++n1; ++n2;
	    int	n = max(n1, b->left);
	    int	n0 =n;
	    int	n9 = min(n2, b->right);
	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    WRECD*	h = hhf[0] + r;
	    WRECD*	g = hhf[1] + r;
	    *f1 = oomW;
	    if (g2) {
		g2  = hhf[2] + r;
		*f2 = oomW;
	    }
	    if (m == mm) {
		int	k = 0;
		for ( ; k < pwd->Nrow; ++k) hb[k] = hhb[k] + r;
		b1 = hhb[k++] + r;
		if (b2) b2 = hhb[k] + r;
	    }
	    for (int d = 0; d < pwd->Noll; ++d) {
#if INTR > 1
		WRECD*	phl = hl[d];
		for (int l = 0; l <= INTR; ++l, ++phl) {
		    nx[d][l] = l;
		    *phl = oomW;
		}
#else
		hl[d] = oomW;
#endif
	    }
	    bool	a_in_zone = api == m;
	    for ( ; n <= n9; ++n, ++h, ++g, ++r, ++bs) {
		VTYPE	x, y;
		hf[0] = h;
		WRECD*	from = h;
		WRECD*	mx = h;
		if (n == n0) goto FindCenter;
		if (g2) ++g2;

/*	Diagonal	*/
		if (m == a->left) goto HorizonF;
		h->val += pwd->sim2(as, bs);
		h->dir = (from->dir & DIAG)? DIAG: NEWD;

/*	Vertical	*/
		x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;

/*	Vertical2	*/
		if (g2) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[1].val) {
		    *g2 = *from;
		    g2->val = x;
		    g2->dir = VERL | (g2->dir & SPJC);
		  } else 	*g2 = g2[1];
		  g2->val += pwd->LongGEP;
		  if (g2->val >= mx->val) mx = g2;
		}
HorizonF:
/*	Horizontal	*/
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f1->val) {
		    *f1 = h[-1];
		    f1->val = x;
		}
		if (f1->dir) {
		    f1->val += pwd->BasicGEP;
		    f1->dir =(f1->dir & SPIN) + HORI;
		    if (f1->val >= mx->val) mx = f1;
		}

/*	Horizontal2	*/
		if (f2) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= f2->val) {
			*f2 = h[-1];
			f2->val = x;
			f2->dir = HORI;
		    }
		    if (f2->dir) {
			f2->val += pwd->LongGEP;
			f2->dir =(f2->dir & SPIN) + HORL;
			if (f2->val >= mx->val) mx = f2;
		    }
		}
		if (mx->dir == NEWD) mx->lst = r;

/*	intron 3' boundary	*/
		if (exin->isAccpt(n)) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = hf[d];
#if INTR > 1
			WRECD*	maxphl = 0;
			int*	pnx = nx[d];
			for (int l = 0; l < INTR; ++l) {
			    WRECD*	phl = hl[d] + pnx[l];
			    if (!phl->dir) break;
			    x = sigJ + phl->val + exin->sig53(phl->jnc, n, IE53) +
				pwd->IntPen->Penalty(n - phl->jnc);
			    if (x > from->val) {
				from->val = x;
				maxphl = phl;
		    	    }
			}
			if (maxphl) {
		 	    from->jnc = maxphl->jnc;
			    from->dir = maxphl->dir | SPJCI;
			    from->upr = max(maxphl->upr, r);
			    from->lwr = min(maxphl->lwr, r);
			    from->lst = maxphl->lst + n - maxphl->jnc;
			    if (from->val > mx->val) mx = from;
			}
#else
			WRECD*	phl = hl + d;
			if (!phl->dir) break;
			x = phl->val + exin->sig53(phl->jnc, n, IE53) +
			    pwd->IntPen->Penalty(n - phl->jnc);
			if (x > from->val) {
			    from->val = x;
		 	    from->jnc = phl->jnc;
			    from->dir = phl->dir | SPJCI;
			    from->upr = max(phl->upr, r);
			    from->lwr = min(phl->lwr, r);
			    from->lst = phl->lst + n - phl->jnc;
			    if (from->val > mx->val) mx = from;
			}
#endif
		    }
		}

/*	Find optimal path	*/
		y = h->val;
		if (h != mx)  {
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		}
		if (LocalL && h->val <= 0) h->val = h->dir = 0;

/*	intron 5' boundary	*/
		if (exin->isDonor(n)) {
		    VTYPE	sigJ = exin->sig53(n, 0, IE5);
		    for (int d = 0; d < pwd->Noll; ++d) {
				/* An orphan exon is disallowed */
			from = hf[d];
			if (!from->dir || (from->dir & SPIN) || isvert(from))
			    continue;
			x = from->val + sigJ;
#if INTR > 1
			WRECD*	phl = hl[d];
			int*	pnx = nx[d];
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
			    phl->jnc = n;
			}
#else
			WRECD*	phl = hl + d;
			if (x > phl->val) {
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = n;
			}
#endif
		    }
		}

#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(f1->val);
		if (g2) {
		    putvar(g2->val); putvar(f2->val);
		}
		if (algmode.lsg) {
#if INTR > 1
		    putvar(hl[0][nx[0][0]].val);
		    putvar(hl[1][nx[1][0]].val);
		    if (g2) putvar(hl[2][nx[2][0]].val);
#else
		    putvar(hl[0].val);
		    putvar(hl[1].val);
		    if (g2) putvar(hl[2].val);
#endif
		}
		putchar('\n');
	}
#endif

/* Find Center */
FindCenter:
		if (m == mm) {
		    int	acc = (*hf)->dir & SPJC;
		    int	dnr = (*hb)->dir & SPJC;
		    int	j1, k1, l;
		    if (((*hf)->dir & SPIN) && ((*hb)->dir & SPIN))
			goto NextCycle;
		    x = (*hf)->val + (*hb)->val;
		    j1 = f1->dir & SPJC;	/* 3' boundary */
		    k1 = b1->dir & SPJC;	/* 5' boundary */
		    if ((j1 && !k1) || (!j1 && k1)) {
			y = f1->val + b1->val - pwd->BasicGOP;
			l = b1->lst - f1->lst - pwd->codonk1;/* gap length */
			if (l > 0) y += pwd->diffu * l;		/* > codonK1 ? */
			if (y > x) {
			    x = y;
			    dnr = k1;
			    acc = j1;
			}
		    }
		    if (f2) {
			int	j2 = f2->dir & SPJC;
			int	k2 = b2->dir & SPJC;
			if ((j1 && !k2) || (!j1 && k2)) {
			    l = r - f1->lst;
			    y = f1->val + b2->val - pwd->BasicGOP + pwd->diffu * l;
			    if (y > x) {
				x = y;
				dnr = k2;
				acc = j1;
			    }
			}
			if ((j2 && !k1) || (!j2 && k1)) {
			    int l = b1->lst - r;
			    y = f2->val + b1->val - pwd->BasicGOP + pwd->diffu * l;
			    if (y > x) {
				x = y;
				dnr = k1;
				acc = j2;
			    }
			}
			if ((j2 && !k2) || (!j2 && k2)) {
			    y = f2->val + b2->val - pwd->LongGOP;
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
			rr[jtigh = 0] = r;
			kk = k1 + 4 * j1;
			jb = ((*hb)->dir & SPJC)? 1: 0;
			kb = (*hb)->dir & HORI;
		    } else if (ge(x, mxh)) {
			if (j1 && dnr) {
			    rr[jtigh = 0] = r;
			} else if (j1 && acc) {
			    rr[jtigh = 1] = r;
			} else if
			(((kb || jb == 2) && ((*hf)->dir & HORI)) ||	/* --* */
			(jb == 1 && (acc || ((*hf)->dir & SPJC))) ||	/* ..| */
			(kb && dnr)) {					/* *-| */
			    if (++jtigh < MAXR) rr[jtigh] = r;
			    if (acc || dnr) ++jb;
			    else	jb = 0;
			    kb = 0;
			}
		    }
NextCycle:
		    for (int d = 0; d < pwd->Nrow; d++) ++hb[d];
		    b1++;
		    if (b2) b2++;
		}
	    }	/* end of n-loop */
	    if (a_in_zone) ++api;		// has exon-exon junctions
	}	/* end of m-loop */

	int	k1 = kk % 4;
	int	k2 = kk / 4;
	int	r = rr[0];
	*nl = r + mm;
	for (int d = 0; d < pwd->Nrow; ++d) hf[d] = hhf[d] + r;
	wdwf->up = max(hf[k1]->upr, hf[k1]->lst);
	wdwf->lw = min(hf[k1]->lwr, hf[k1]->lst);;
	wdwf->width = wdwf->up - wdwf->lw + 3;
	if (k1) {
	    wskl.m = *ml = *nl - hf[k1]->lst;
	    wskl.n = *nl;
	    mfd->write((UPTR) &wskl);
	} else {
	    wskl.m = *ml = mm;
	}
	*nr = (r = rr[jtigh]) + mm;
	for (int d = 0; d < pwd->Nrow; ++d) hb[d] = hhb[d] + r;
	wdwb->up = max(hb[k2]->upr, hb[k2]->lst);
	wdwb->lw = min(hb[k2]->lwr, hb[k2]->lst);
	wdwb->width = wdwb->up - wdwb->lw + 3;
	if (k2) {
	    wskl.n = *nr;
	    mfd->write((UPTR) &wskl);
	    wskl.m = *mr = *nr - hb[k2]->lst;
	    mfd->write((UPTR) &wskl);
	} else {
	    *mr = mm;
	    for ( ; jtigh >= 0; --jtigh) {
		wskl.m = mm;
		wskl.n = rr[jtigh] + mm;
		mfd->write((UPTR) &wskl);
	    }
	}
	delete[] (hhf[0] + wdw->lw - 1);
	delete[] (hhb[0] + wdw->lw - 1);
	return (mxh);
}

static const int pNoll = 2;
static const int pNrow = 2;
static const int Nedge = 3;

void Aln2s1::pfinitS_ng(LRECD* hhg[])
{
	int	r = b->left - a->left;
	int	rr = b->right - a->left;
	LRECD*	hf[pNrow];

	for (int d = 0; d < pNrow; ++d) hf[d] = hhg[d] + r;
	if (wdw->up < rr) rr = wdw->up;
	LRECD*	h = hf[0]++;
	h->val = 0;
	h->dir = DIAG;
	h->jnc = b->left;
	h->lst = r;
	h->ptr = vmf? vmf->add(a->left, b->left, 0): 0;
 	for (int d = 1; d < pNrow; ++d) *hf[d]++ = ooml;
	while (++r <= wdw->up + 1) {
 	    for (int d = 0; d < pNrow; ++d) *hf[d]++ = ooml;
	}

	r = b->left - a->left;
	rr = b->left - a->right;
	for (int d = 0; d < pNrow; ++d) hf[d] = hhg[d] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	for (int i = 1; --r >= wdw->lw - 1; ++i) {
	    h = --hf[0];
	    if (r < rr) *h = ooml;
	    else {			/* (semi) global */
		*h = h[1];
		if (i == 1) {
		    h->val += pwd->BasicGOP;
		    h->dir = VERT;
		}
		h->val += pwd->BasicGEP;
	    }
	    *--hf[1] = *h;
	}
}

void Aln2s1::pbinitS_ng(LRECD* hhg[])
{
	int	r = b->right - a->right;
	int	rr = b->left - a->right;
	LRECD*	hb[pNrow];

	for (int d = 0; d < pNrow; ++d) hb[d] = hhg[d] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	LRECD*	h = hb[0]--;
	h->val = 0;
	h->dir = DIAG;
	h->jnc = b->right;
	h->lst = r;
	h->ptr = vmf? vmf->add(a->right, b->right, 0): 0;
	for (int d = 1; d < pNrow; ++d) *hb[d]-- = ooml;
	while (--r >= wdw->lw - 1) {
	    for (int d = 0; d < pNrow; ++d) *hb[d]-- = ooml;
	}

	r = b->right - a->right;
	rr = b->right - a->left;
	if (wdw->up < rr) rr = wdw->up;
	for (int d = 0; d < pNrow; ++d) hb[d] = hhg[d] + r;
	for (int i = 1; ++r <= wdw->up + 1; ++i) {
	    h = ++hb[0];
	    if (r > rr) *h = ooml;
	    else {
		*h = h[-1];
		if (i == 1) {
		    h->val += pwd->BasicGOP;
		    h->dir = VERT;
		}
		h->val += pwd->BasicGEP;
	    }
	    *++hb[1] = *h;
	}
}

VTYPE Aln2s1::pincersS_ng(long *ptr, int cmode)
{
	int	enda = a->inex.exgl && algmode.lcl < 16;
	int	endb = b->inex.exgl && algmode.lcl < 16;
	int	nx = 0;
	int	n = 0, n1, n2, n9 = 0;
	LRECD*	hhg[pNrow];
	LRECD*	hf[Nedge];
	LRECD	f[pNoll];
	LRECD*	h = 0;
	LRECD*	f1 = f + 1;
	LRECD*	b1 = f1;
	LRECD	maxl = ooml;	/* left  side */
	LRECD	maxr = ooml;	/* right side */
	int	dd = a->right - a->left + 1;
	EIJUNC*	hhl = new EIJUNC[dd] - a->left;	/* intron junc	*/
	EIJODR*	nnl = new EIJODR[dd] - a->left;	/* intron junc	*/
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	VSKLP	lsth = maxh;
	VTYPE	maxscr = NEVSEL;	/* overall */

	hhg[0] = new LRECD[pNrow * wdw->width] - wdw->lw + 1;
	for (int d = 1; d < pNrow; ++d) hhg[d] = hhg[d-1] + wdw->width;
	for (int d = 1; d < pNoll; ++d) hf[2 * d] = f + d;
	vmf->add(0, 0, 0);	/* Skip 0-th record */
	int	m = a->left;
	CHAR*	as = a->at(m);
	PfqItr	api(a, m);
	int	api_size = api.size();
	if (cmode == 2) goto forward;	// forward only

/*******	backward phase	*******/

	pbinitS_ng(hhg);
	m = a->right;
	if (!a->inex.exgr) ++m; /* global */
	else {
	    for (int l = 0; l <= NCANDS; ++l) {
		hhl[m][l] = ooml;
		nnl[m][l] = l;
	    }
	}
	as = a->at(m);
	n1 = m + wdw->lw;
	n2 = m + wdw->up + 1;
	while (--m >= a->left) {
	    for (int l = 0; l <= NCANDS; ++l) {
		hhl[m][l] = ooml;
		nnl[m][l] = l;
	    }
	    --as; --n1; --n2;
	    n  = min(n2, b->right);
	    n9 = max(n1, b->left);
	    int	r = n - m;
	    int	nr = n + 1;
	    int	peak = 0;
	    CHAR*	bs = b->at(n);
	    h = hhg[0] + r;
	    LRECD*	g = hhg[1] + r;
	    LRECD*	mxd = ((h->val + pwd->Vthr) < maxh.val)? &ooml: h;
	    *b1 = ooml;

#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m+1, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    while (--n >= n9) {
		VTYPE	x;
		--bs; --r;
		hf[0] = --h;
		hf[1] = --g;

/*	Diagonal	*/
		LRECD*	from = h;
		LRECD*	mx = h;
		if (m == a->right) goto HorizonB;
		h->val += pwd->sim2(as, bs);
		h->dir = (from->dir & DIAG)? DIAG: NEWD;

/*	Vertical	*/
		x = (--from)->val + pwd->BasicGOP;
		if (x >= g[-1].val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT;
		} else *g = g[-1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;
HorizonB:
/*	Horizontal	*/
		x = h[1].val + pwd->BasicGOP;
		if (x >= b1->val) {
		    *b1 = h[1];
		    b1->val = x;
		}
		if (b1->dir) {
		    b1->val += pwd->BasicGEP;
		    b1->dir = HORI;
		    if (b1->val >= mx->val) mx = b1;
		}

		if (mx->dir == NEWD) {
		    mx->lst = r;
		    if (vmf) mx->ptr = vmf->add(m + 1, n + 1, mx->ptr);
		}

/*	Find optimal path	*/
		if (mx->val > maxh.val) {
		    maxh.val = mx->val;
		    maxh.p = mx->ptr;
		    maxh.m = m;
		    maxh.n = n;
		}
		if (mx->val + pwd->Vthr < maxh.val) {	/* x drop off */
		    if (peak) {		/* end of block */
			n1 = n + 1;
			peak = 0;
		    }
		    nr = n;		/* before block	*/
		} else if ((mx->dir & DIAG) && mx->val >= mxd->val) {
		    mxd = mx;
		    if (nr < n2) n2 = nr;
		    nx = n;
		    peak = 1;
		}
		if (h != mx) *h = *mx;

/*	intron 3' boundary reverse phase	*/

		if (cmode == 3 && exin->isAccpt(n)) {
		    LRECD*	phl = hhl[m];
		    int*	pnx = nnl[m];
		    int		p = 0;
		    for ( ; p < NCANDS && phl[pnx[p]].dir; ++p) ;
		    VTYPE	sigJ = exin->sig53(0, n, IE3);
		    for (int d = 0; d < Nedge; ++d) {
			int	l = (d + 1) / 2;
			from = hf[d];
			if (!from->dir || (d && (mx == from || 
			from->val <= mx->val + pwd->GOP[l])))
			    continue;
			x = from->val + sigJ;
			if (p == NCANDS) --p;
			for (l = p++; l >= 0; --l) {
			    if (x <= phl[pnx[l]].val)
				break;
			    else if (l < NCANDS - 1)
				gswap(pnx[l], pnx[l + 1]);
			}
			if (++l < NCANDS) {
			    LRECD*	phr = phl + pnx[l];
			    *phr = *from;
			    phr->val = x;
			    phr->jnc = n;
			}
		    }
		}

#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(b1->val);
		if (algmode.lsg) {
		    putvar(hhl[m][nnl[m][0]].val);
		    putvar(hhl[m][nnl[m][1]].val);
		    putvar(hhl[m][nnl[m][2]].val);
		}
		putchar('\n');
	}
#endif
		if (m == a->left && enda && h->val > lsth.val)
		    {lsth.val = h->val; lsth.p = h->ptr; lsth.m = m; lsth.n = n;}
	    }	// end of n loop
	    if (!mxd->dir) break;	// no peak
	    if (peak) n1 = n;
	    LRECD*	phr = hhl[m] + nnl[m][NCANDS];
	    *phr = *mxd;		// maxscr in row
	    phr->jnc = nx;
	    if (++n == b->left && endb && h->val > lsth.val)
		{lsth.val = h->val; lsth.p = h->ptr; lsth.m = m; lsth.n = n;}
	}	// end of m loop
	if (cmode == 1) {
	    if (!a->inex.exgl && !b->inex.exgl && h)		// global
		{lsth.val = h->val; lsth.p = h->ptr; lsth.m = ++m; lsth.n = max(n, n9);}
	    else if (lsth.val == NEVSEL) lsth = maxh;
	    maxscr = lsth.val;
	    ptr[0] = 0;
	    if (vmf) ptr[1] = vmf->add(lsth.m, lsth.n, lsth.p);
	    goto freeprec;
	} else	++m;

/*******	forward phase	*******/

forward:
	pfinitS_ng(hhg);
	if (api_size) api.reset(m);
	if (!a->inex.exgl) --m; /* global */
	n1 = m + wdw->lw - 1;
	n2 = m + wdw->up;
	as = a->at(m);
	maxh.val = lsth.val = NEVSEL;
	enda = a->inex.exgr && algmode.lcl < 16;
	endb = b->inex.exgr && algmode.lcl < 16;
	for ( ; ++m <= a->right; ++as) {
	    ++n1; ++n2;
	    n = max(n1, b->left);
	    n9 = min(n2, b->right);
	    int	nr = n - 1;
	    int	peak = 0;
	    int	r = n - m;
	    CHAR*	bs = b->at(n);
	    LRECD*	phr = hhl[m] + nnl[m][NCANDS];
	    h = hhg[0] + r;
	    LRECD*	g = hhg[1] + r;
	    LRECD*	mxd = ((h->val + pwd->Vthr) < maxh.val)? &ooml: h;
	    *f1 = ooml;

#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    bool	a_in_zone = api_size && (api == m);
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x;
		++r;
		hf[0] = ++h;
		hf[1] = ++g;

/*	Diagonal	*/
		LRECD*	from = h;
		LRECD*	mx = h;
		if (m == a->left) goto HorizonP;
		h->val += pwd->sim2(as, bs);
		h->dir = (from->dir & DIAG)? DIAG: NEWD;

/*	Vertical	*/
		x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT;
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;
HorizonP:
/*	Horizontal	*/
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f1->val) {
		    *f1 = h[-1];
		    f1->val = x;
		}
		if (f1->dir) {
		    f1->val += pwd->BasicGEP;
		    f1->dir = HORI;
		    if (f1->val >= mx->val) mx = f1;
		}

		if (mx->dir == NEWD) {
		    mx->lst = r;
		    mx->ptr = vmf->add(m - 1, n - 1, mx->ptr);
		}

/*	Find optimal path	*/
		if (mx->val > maxh.val) {
		    maxh.val = mx->val;
		    maxh.m = m;
		    maxh.n = n;
		    maxh.p = mx->ptr;
		}
		if (mx->val + pwd->Vthr < maxh.val) {	/* x drop off */
		    if (peak) {		/* end of block */
			n2 = n - 1;
			peak = 0;
		    }
		    nr = n;
		} else if ((mx->dir & DIAG) && mx->val >= mxd->val) {
		    mxd = mx;
		    if (nr > n1) n1 = nr;
		    peak = 1;
		}
		if (h != mx) *h = *mx;

/*	intron 5' boundary	*/
		if (cmode == 3 && exin->isDonor(n)) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    for (int k = 0; k < Nedge; ++k) {
			int	l = (k + 1) / 2;
			LRECD*	phl = hf[k];
			if (!phl->dir || (k && (mx == phl ||
				phl->val <= mx->val + pwd->GOP[l])))
			    continue;
			for (l = 0; l < NCANDS; ++l) {
			    phr = hhl[m] + nnl[m][l];
			    if (!phr->dir) break;
			    x = sigJ + phl->val + 
				phr->val + pwd->IntPen->Penalty(phr->jnc - n) +
				exin->sig53(n, phr->jnc, IE35);
			    if ((phl->dir == VERT && phr->dir == VERT) ||
				((phl->dir & HORI) && (phr->dir & HORI))) {
				    x -= pwd->BasicGOP;
			    }
			    if (x > maxscr) {
				maxscr = x;
				maxl = *phl;
				maxr = *phr;
				maxl.ptr = vmf->add(m, n, phl->ptr);
				maxl.ptr = vmf->add(m, phr->jnc, maxl.ptr);
			    }
			}
		    }
		}
#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(f1->val);
		putchar('\n');
	}
#endif
		if (m == a->right && enda && h->val > lsth.val)
		    {lsth.val = h->val; lsth.p = h->ptr; lsth.m = m; lsth.n = n;}
	    }	/* end of n-loop */
	    if (peak) n2 = n;
	    if (cmode == 3) {
		LRECD*	phr = hhl[m] + nnl[m][NCANDS];
		int	l = phr->lst - mxd->lst;
		if (l < 0) continue;
		VTYPE	x = mxd->val + phr->val + pwd->GapPenalty(l);
		if (x > maxscr) {
		    maxscr = x;
		    maxl = *mxd;
		    maxr = *phr;
		    maxl.ptr = vmf->add(m, m + mxd->lst, mxd->ptr);
		    maxl.ptr = vmf->add(m, m + phr->lst, maxl.ptr);
		}
	    }
	    if (!mxd->dir) break;
	    if (a_in_zone) ++api;		// has exon-exon junctions
	    if (--n == b->right && endb && h->val > lsth.val)
		{lsth.val = h->val; lsth.p = h->ptr; lsth.m = m; lsth.n = n;}
	}	/* end of m-loop */
	if (cmode == 2) {
	    if (!a->inex.exgr && !b->inex.exgr && h)	// global
		{lsth.val = h->val; lsth.p = h->ptr; lsth.m = --m; lsth.n = min(n, n9);}
	    else if (lsth.val == NEVSEL) lsth = maxh;	// local
	    maxscr = lsth.val;
	    ptr[0] = vmf->add(lsth.m, lsth.n, lsth.p);
	    ptr[1] = 0;
	} else {
	    ptr[0] = maxl.ptr;
	    ptr[1] = maxr.ptr;
	}

freeprec:
	delete[] (hhg[0] + wdw->lw - 1);
	delete[] (hhl + a->left);
	delete[] (nnl + a->left);
	return (maxscr);
}

void Aln2s1::cinitS_ng(CRECD* hhc[])
{
	int	n = b->left;
	int	r = b->left - a->left;
	int	rr = b->right - a->left;
	CRECD*	hf[NOL];

	for (int d = 0; d < pwd->Nrow; ++d) hf[d] = hhc[d] + r;
	if (wdw->up < rr) rr = wdw->up;
	CRECD*	h = hhc[0] + r;
	for ( ; r <= wdw->up + 1; ++h, ++r) {
	    for (int d = 1; d < pwd->Nrow; ++d) *hf[d]++ = oomC;
	    if (r > rr) {
		*hf[0]++ = oomC;
		continue;
	    }
	    h = hf[0]++;
	    h->val = 0;
	    h->dir = DEAD;
	    h->upr = h->lwr = r;
	    h->jnc = n;
	    h->mlb = a->left;
	    h->nlb = n++;
	}

	r = b->left - a->left;
	rr = b->left - a->right;
	for (int d = 0; d < pwd->Nrow; ++d) hf[d] = hhc[d] + r;
	if (wdw->lw > rr) rr = wdw->lw;
	while (--r >= wdw->lw - 1) {
	    for (int d = 0; d < pwd->Nrow; ++d) *--hf[d] = oomC;
	    if (r < rr) continue;
	    h = hf[0];
	    h->val = 0;
	    h->dir = DEAD;
	    h->upr = h->lwr = r;
	    h->mlb = b->left - r;
	    h->nlb = b->left;
	    h->jnc = b->left;
	}
}

Colonies* Aln2s1::fwdswgS_ng(VTYPE* scr)
{
	CRECD*	hhc[NOL];	/* local H matrix	*/
	CRECD*	hf[NOL]; 	/* [DIAG, HORI, HORL] */
	CRECD	f[NOL]; 	/* [DIAG, HORI, HORL] */
#if INTR > 1
	CRECD	hl[NOL][INTR+1];
	int	nl[NOL][INTR+1]; /* [DIA, HORI, HORL][candidates] */
#else
	CRECD   hl[NOL];
#endif
	CRECD*	g2 = 0;
	Colonies*	cl = new Colonies(0);
	COLONY*	clny = cl->at();

	hhc[0] = new CRECD[pwd->Noll * wdw->width] - wdw->lw + 1; /* +1: sp junc */
	for (int d = 1; d < pwd->Noll; ++d) hhc[d] = hhc[d-1] + wdw->width;
	cinitS_ng(hhc);
	for (int d = 1; d < pwd->Noll; ++d) hf[d] = f + d;
	int	m = a->left;
	CHAR*	as = a->at(m);
	PfqItr	api(a, m);
	int	api_size = api.size();
	int	n1 = m + wdw->lw - 1;
	int	n2 = m + wdw->up;
	for ( ; ++m <= a->right; ++as) {
	    ++n1; ++n2;
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
	    int	r = n - m;
	    CHAR*	bs = b->at(n);
	    CRECD*	h = hhc[0] + r;
	    CRECD*	 g = hhc[1] + r;
	    if (pwd->Noll == 3) g2 = hhc[2] + r;
	    for (int d = 0; d < pwd->Noll; ++d) {
#if INTR > 1
		CRECD*	phl = hl[d];
		for (int l = 0; l <= INTR; ++l, ++phl) {
		    nl[d][l] = l;
		    *phl = oomC;
		}
#else
		hl[d] = oomC;
#endif
		f[d] = oomC;
	    }

#if DEBUG
	if (algmode.nsa & 8) {
	    printf("%2d %2d %2d", m+1, n, h->dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    bool	a_in_zone = api_size && (api == m);
	    for ( ; ++n <= n9; ++bs) {
		++g; ++r;
		hf[0] = ++h;
		if (g2) ++g2;
/*	Diagonal	*/
		CRECD*	from = h;
		CRECD*	mx = h;
		f[0] = *h;	/* diag */
		h->val += pwd->sim2(as, bs);
		h->dir = (from->dir & DIAG)? DIAG: NEWD;

/*	Vertical	*/
		VTYPE	x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    *g = *from;
		    g->val = x;
		    g->dir = VERT | (g->dir & SPJC);
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val > mx->val) mx = g;

/*	Vertical2	*/
		if (g2) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[1].val) {
		    *g2 = *from;
		    g2->val = x;
		    g2->dir = VERL | (g2->dir & SPJC);
		  } else 	*g2 = g2[1];
		  g2->val += pwd->LongGEP;
		  if (g2->val > mx->val) mx = g2;
		}

/*	Horizontal	*/
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f[1].val) {
		    f[1] = h[-1];
		    f[1].val = x;
		}
		if (f[1].dir) {
		    f[1].val += pwd->BasicGEP;
		    f[1].dir = (f[1].dir & SPIN) + HORI;
		    if (f[1].val >= mx->val) mx = f + 1;
		}

/*	Horizontal2	*/
		if (pwd->Noll == 3) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= f[2].val) {
			f[2] = h[-1];
			f[2].val = x;
		    }
		    if (f[2].dir) {
			f[2].val += pwd->LongGEP;
			f[2].dir = (f[2].dir & SPIN) + HORL;
			if (f[2].val >= mx->val) mx = f + 2;
		    }
		}

/*	intron 3' boundary, assume no overlapping signals	*/
		if (exin->isAccpt(n)) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    for (int d = 0; d < pwd->Noll; ++d) {
			from = hf[d];
#if INTR > 1
			CRECD*	maxphl = 0;
			int*	pnl = nl[d];
			VTYPE	y = NEVSEL;
			for (int l = 0; l < INTR; ++l) {
			    CRECD*	phl = hl[d] + pnl[l];
			    if (!phl->dir) break;
			    x = sigJ + phl->val + pwd->IntPen->Penalty(n - phl->jnc)
				+ exin->sig53(phl->jnc, n, IE53);
			    if (x > from->val) {
				y = x;
				maxphl = phl;
		    	    }
			}
			if (maxphl) {
			    *from = *maxphl;
			    from->val = y;
		 	    from->jnc = n;
			    from->dir |= SPJCI;
			    if (r > from->upr) from->upr = r;
			    if (r < from->lwr) from->lwr = r;
			    if (from->val > mx->val) mx = from;
			}
#else
			CRECD*	phl = hl + d;
			if (!phl->dir) break;
			x = phl->val + pwd->IntPen->Penalty(n - phl->jnc)
			    + exin->sig53(phl->jnc, n, IE53);
			if (x > from->val) {
			    *from = *phl;
			    from->val = x;
		 	    from->jnc = n;
			    from->dir |= SPJCI;
			    if (r > from->upr) from->upr = r;
			    if (r < from->lwr) from->lwr = r;
			    if (from->val > mx->val) mx = from;
			}
#endif

		    }
		}

/*	intron 5' boundary	*/

		if (exin->isDonor(n)) {
		    VTYPE	sigJ = exin->sig53(n, 0, IE5);
		    for (int d = 0; d < pwd->Noll; ++d) {
				/* An orphan exon is disallowed */
			from = hf[d];
			if (!from->dir || (from->dir & SPIN) || 
				from->dir == VERT || from->dir == VERL)
			    continue;
			x = from->val + sigJ;
#if INTR > 1
			CRECD*	phl = hl[d];
			int*	pnl = nl[d];
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
			    phl->jnc = n;
			}
#else
			CRECD*	phl = hl + d;
			if (x > phl->val) {
			    *phl = *from;
			    phl->val = x;
			    phl->jnc = n;
			}
#endif
		    }
		}

/*	Find optimal path	*/
		if (h != mx) {		/* non-diagonal */
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		} else if (h->val > f->val) {
		    if (f->val == 0 && !(h->dir & SPJC)) {	/* new colony */
			h->upr = h->lwr = r;
			h->mlb = m - 1;
			h->nlb = n - 1;
		    }
		    if (h->val > clny->val) {	/* max local score */
			clny->val = h->val;
			clny->mrb = m;
			clny->nrb = n;
			clny->lwr = h->lwr;
			clny->upr = h->upr;
			clny->mlb = h->mlb;
			clny->nlb = h->nlb;
		    }
		}
		if (h->val < 0) {		/* reset to blank */
		    h->val = f[1].val = g->val = 0;
		    if (g2) f[2].val = g2->val = 0;
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
			h->val = f[1].val = g->val = 0;	/* X-drop */
			if (g2) f[2].val = g2->val = 0;
			h->clny = 0;
		    }
		}
#if DEBUG
	if (algmode.nsa & 8) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(f[1].val);
		if (g2) {
		    putvar(g2->val); putvar(f[2].val);
		}
		if (algmode.lsg) {
#if INTR > 1
		    putvar(hl[0][nl[0][0]].val); 
		    putvar(hl[1][nl[1][0]].val);
		    if (g2) putvar(hl[2][nl[2][0]].val);
#else
		    putvar(hl[0].val); 
		    putvar(hl[1].val);
		    if (g2) putvar(hl[2].val);
#endif
		}
		putchar('\n');
	}
#endif

	    }
	    if (a_in_zone) ++api;		// has exon-exon junctions
	}

	delete[] (hhc[0] + wdw->lw - 1);
	*scr = clny->val;
	cl->sortcolonies();
	return (cl);
}

VTYPE Aln2s1::diagonalS_ng()
{
	CHAR*	as = a->at(a->left);
	CHAR*	bs = b->at(b->left);
	VTYPE	scr = 0;
	VTYPE	maxh = NEVSEL;
	SKL	wskl;
	int	mL = a->left;
	int	mR = a->right;
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

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
	int m = b->left - a->left;
	if (vmf) {
	    wskl.m = mL;
	    wskl.n = mL + m;
	    mfd->write((UPTR) &wskl);
	    wskl.m = mR;
	    wskl.n = mR + m;
	    mfd->write((UPTR) &wskl);
	} else if (a->inex.exgl && a->left == 0) {
	    wskl.m = a->left;
	    wskl.n = b->left;
	    mfd->write((UPTR) &wskl);
	} else if (a->inex.exgr && a->right == a->len) {
	    wskl.m = a->right;
	    wskl.n = b->right;
	    mfd->write((UPTR) &wskl);
	}
	return (LocalR? maxh: scr);
}

VTYPE Aln2s1::trcbkalignS_ng()
{
	long	ptr = 0;

	vmf = new Vmf();
	VTYPE	scr = forwardS_ng(&ptr);
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
	delete vmf;
	vmf = 0;
	return (scr);
}

VTYPE Aln2s1::lspS_ng()  /* recursive */ 
{
	WINDOW*	wrsv = wdw;
	WINDOW	wdwl;
	WINDOW	wdwr;
	RANGE	rng[2];			/* reserve */
	int	aexg = a->inex.exgr;	/* reserve */
	int	bexg = b->inex.exgr;	/* reserve */

	if (wdw->up == wdw->lw) return(diagonalS_ng());
	int	m = a->right - a->left;
	int	n = b->right - b->left;
	int	k = wdw->lw - b->left + a->right;
	int	q = b->right - a->left - wdw->up;
	long	cvol =  m * n - (k * k + q * q) / 2;
	if (cvol < MaxVmfSpace || m == 1 || n <= 1) {
	    return (trcbkalignS_ng());
	}
	save_range(seqs, rng, 2);
	VTYPE	scr = centerS_ng(&m, &q, &n, &k, &wdwl, &wdwr);
	a->inex.exgr = b->inex.exgr = 0;
	int	r = b->left - a->left;
	if (r < wdwl.lw) b->left = a->left + wdwl.lw;
	if (r > wdwl.up) a->left = b->left - wdwl.up;
	a->right = m;
	b->right = n;
	wdw = &wdwl;
	lspS_ng();
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
	wdw = &wdwr;
	lspS_ng();
	wdw = wrsv;
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexg;
	b->inex.exgl = bexg;
	return scr;
}

VTYPE Aln2s1::pincerTrcbkS_ng(int cmode)
{
	vmf = new Vmf();
	long	ptr[2];
	VTYPE	scr = pincersS_ng(ptr, cmode);;
	SKLP    sv;

	sv.p = ptr[0];
	while (sv.p) {
	    vmf->readvmf(&sv, sv.p);
	    mfd->write((UPTR) &sv);
	}
	sv.p = ptr[1];	/* restore */
	while (sv.p) {
	    vmf->readvmf(&sv, sv.p);
	    mfd->write((UPTR) &sv);
	}
	delete vmf;
	vmf = 0;
	return (scr);
}

VTYPE Aln2s1::backforth(int ovr, BOUND& lub)
{
	VTYPE*	bscr = new VTYPE[ovr + 1];
	CHAR*	as = a->at(a->left);
	CHAR*	bs = b->at(b->left);
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

// indel-free version of pincsersS within overlapped region
// doubly conted alignment score is cancelled
// assume agap <= 0
// return whether a cannonical pair is found or not

bool Aln2s1::indelfreespjS(int agap, VTYPE& iscr)
{
	int	play = -b->exin->lplay(b->left);
	if (play > agap) agap = play;
	int	d5 = b->left + agap;	// donor 5' end
	int	d3 = b->left;		// donor 3' end
	int	ilen = b->right - d5;	// intron length
	int	i = 1 - agap;
	VTYPE*	backward = new VTYPE[i];
	CHAR*	as = a->at(a->left);	// cancel score of
	CHAR*	bs = b->at(b->left);	// donor side overlap
	for (VTYPE v = 0; --i >= 0; v += pwd->sim2(--as, --bs))
	    backward[i] = v;
	int	m = a->right;
	int	mm = m;
	int	n = d5;
	int	nn = n;
	as = a->at(a->right);
	PfqItr	api(a, m);
	bs = b->at(b->right);
	VTYPE	maxspj = iscr = NEVSEL;
	VTYPE	ip = pwd->IntPen->Penalty(ilen);
	bool	usespb = api.size() && use_spb();
	for (VTYPE v = i = 0; n <= d3; ++n, ++i, ++m) {
	    bool	a_in_zone = usespb && (api == m);
	    if (exin->isCanon(n, n + ilen)) {
		VTYPE	x = ip + exin->sig53(n, n + ilen, IE5P3);
		if (a_in_zone) x += api.match_score(m);
		VTYPE	y = x - v - backward[i];
		if (y > iscr) {nn = n; mm = m; iscr = y; maxspj = x;}
	    }
	    v += pwd->sim2(as++, bs++);	// cancel score of acc side
	    if (a_in_zone) ++api;
	}
	delete[]	backward;
	if (maxspj > NEVSEL) {
	    SKL	tmp = {mm, nn};
	    mfd->write(&tmp);
	    tmp.n += ilen;
	    mfd->write(&tmp);
	    return (true);
	} else
	    return (false);
}

VTYPE Aln2s1::creepback(int ovr, VTYPE bscr, BOUND& lub)
{
	CHAR*	as = a->at(a->left);
	CHAR*	bs = b->at(b->left);
	VTYPE	dscr = 0;
	while (a->left > lub.la && b->left > lub.lb
		&& (ovr < 0 || dscr < bscr)) {
	    dscr += pwd->sim2(--as, --bs);
	    --a->left; --b->left;
	    if (++ovr == 0) bscr += dscr;
	}
	return (dscr);
}

VTYPE Aln2s1::creepfwrd(int& ovr, VTYPE bscr, BOUND& lub)
{
	CHAR*	as = a->at(a->right);
	CHAR*	bs = b->at(b->right);
	VTYPE	dscr = 0;
	while (a->right < lub.ua && b->right < lub.ub
		&& (ovr < 0 || dscr < bscr)) {
	    dscr += pwd->sim2(as++, bs++);
	    ++a->right; ++b->right;
	    if (++ovr == 0) bscr += dscr;
	}
	return (dscr);
}

/* eimode 1: 5'end, 2: 3' end, 3: internal */
VTYPE Aln2s1::seededS_ng(INT level, int eimode, BOUND& lub)
{
	INEX	ainex = a->inex;
	INEX	binex = b->inex;
	RANGE	rng[2];
	int	cmode = eimode;
	int	qck = algmode.blk && algmode.crs && level == algmode.qck && algmode.lcl < 16;
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
	int	term = (int) ((pwd->Vthr + pwd->BasicGOP) / pwd->BasicGEP);
	int	backstep = wlprm->width;
	if (++level < algmode.qck) {
	    wlprm = setwlprm(level);
	    backstep = wlprm->width;
	}
	BOUND	bab = lub;
	if (num) {
	    jxt[num].jx = a->right;
	    jxt[num].jy = b->right;
	    a->inex.exgr = 0;
	    b->inex.exgr = 0;
	    agap = jxt->jx - a->left;
	    bgap = jxt->jy - b->left;
	    for (wjxt = jxt; num--; wjxt++) {
		a->right = wjxt->jx;
		b->right = wjxt->jy;
		bab.ua = wjxt->jx + wjxt->jlen;
		if (wjxt[1].jx < bab.ua) bab.ua = wjxt[1].jx;
		bab.ua -= backstep;
		int	lx = wjxt->jx + wjxt->jlen / 2;
		if (bab.ua < lx) bab.ua = lx;
		bab.ub = wjxt->jy + bab.ua - wjxt->jx;
		int	ovr = min(agap, bgap);
		float	dpspace = (float) agap * (float) bgap / MEGA;
		VTYPE	iscore;
		scr += wjxt->jscr;
		if (cmode == 1 && agap == 0)	// 5' end
		    mfd->write((UPTR) wjxt);
		else if (cmode == 3 && agap <= 0 &&
		    bgap >= IntronPrm.llmt && indelfreespjS(agap, iscore))
			scr += iscore;
		else if (level < algmode.qck && ovr >= (int) wlprm->tpl)
		    scr += seededS_ng(level, cmode, bab);	// recursive search
		else if (ovr <= 0 && bgap < IntronPrm.llmt)
		    scr += backforth(-ovr, bab);
		else {
		    scr -= creepback(ovr, slmt, bab);
		    scr -= creepfwrd(ovr, slmt, bab);
		    if (ovr < 0) {		// skip this hsp
			agap = wjxt[1].jx - a->left;
			bgap = wjxt[1].jy - b->left;
			continue;
		    }
		    if (cmode == 1 && (qck || ovr < IntronPrm.tlmt
			|| dpspace >= alprm.maxsp)) {
			int	bl = wjxt->jy + term;
			if (bl > b->left) b->left = bl;
			stripe(seqs, wdw, alprm.sh);
			scr += pincerTrcbkS_ng(cmode);
		    } else if (cmode == 3 && agap < IntronPrm.elmt
			&& bgap >= IntronPrm.llmt) {
			stripe(seqs, wdw, alprm.sh); // intron ?
			scr += pincerTrcbkS_ng(cmode);
		    } else if (dpspace < alprm.maxsp) {
			stripe(seqs, wdw, alprm.sh);
			if (wdw->up == wdw->lw) diagonalS_ng();
			else	trcbkalignS_ng();
		    } else { 			// give up alignment
//			scr += pincerTrcbkS_ng(seqs, alnv, cmode);
			mfd->write((UPTR) wjxt);
			mfd->write((UPTR) (wjxt + 1));
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
	VTYPE	iscore;
	if (cmode == 2 && ovr == 0) {
	    SKL	tmp = {a->left, b->left};
	    mfd->write((UPTR) &tmp);
	} else if (cmode == 3 && agap <= 0 && bgap >= IntronPrm.llmt && 
	    indelfreespjS(agap, iscore)) {
	    scr += iscore;
	} else if (level < algmode.qck && ovr > (int) wlprm->tpl) {
	    scr += seededS_ng(level, cmode, bab);
	} else if (ovr <= 0 && bgap < IntronPrm.llmt) {
	   scr += backforth(-ovr, bab);
	} else {
	    scr -= creepback(ovr, slmt, bab);
	    scr -= creepfwrd(ovr, slmt, bab);
	    float	dpspace = (float) agap * (float) bgap / MEGA;
	    if (cmode == 1 && (qck || ovr < IntronPrm.tlmt
		|| dpspace >= alprm.maxsp)) {
		stripe(seqs, wdw, alprm.sh);
		scr += pincerTrcbkS_ng(cmode);
	    } else if (cmode == 2 && (qck || ovr < IntronPrm.tlmt
		|| dpspace >= alprm.maxsp)) {
		if (wjxt) {
		    int	br = wjxt->jy - term;
		    if (br > b->left && br < b->right) b->right = br;
		}
		stripe(seqs, wdw, alprm.sh);
		scr += pincerTrcbkS_ng(cmode);
	    } else if (cmode == 3 && agap < IntronPrm.elmt && 
		bgap >= IntronPrm.llmt) {
		stripe(seqs, wdw, alprm.sh);
		scr += pincerTrcbkS_ng(cmode);
	    } else if (dpspace < alprm.maxsp) {
		stripe(seqs, wdw, alprm.sh);
		if (wdw->up == wdw->lw) diagonalS_ng();
		else	trcbkalignS_ng();
	    } else {	// cmode == 3, give up alignment
		SKL	tmp = {a->left, b->left};
		mfd->write((UPTR) &tmp);
		tmp.m = a->right; tmp.n = b->right;
		mfd->write((UPTR) &tmp);
	    }
	}
	rest_range(seqs, rng, 2);
	a->inex.exgl = ainex.exgl;
	b->inex.exgl = binex.exgl;
	if (--level == lowestlvl && wjxt) {
	    wjxt->jx = a->len;
	    wjxt->jy = b->len;
	}
	delete wl;
	return (scr);
}
SKL* Aln2s1::globalS_ng(VTYPE* scr)
{
	mfd = new Mfile(sizeof(SKL));
	SKL	wsk;
	mfd->write((UPTR) &wsk);	// dummy call
	if (algmode.qck) {
	    BOUND bab = {a->left, b->left, a->right, b->right};
	    *scr = seededS_ng(lowestlvl, 1, bab);
	} else
	    *scr = lspS_ng();
	wsk.n = (int) mfd->size();
	SKL*	skl = (SKL*) mfd->flush();
	skl->n = wsk.n - 1;
	skl->m = 1;
	if (skl->n == 0) {
	    delete[] skl;
	    return (0);
	}
	return (trimskl(seqs, stdskl(&skl)));
}


VTYPE HomScoreS_ng(Seq* seqs[], PwdB* pwd)
{
	WINDOW	wdw;

	stripe(seqs, &wdw, alprm.sh);
	Aln2s1 alnv(seqs, pwd, &wdw);
	return alnv.forwardS_ng(0);
}

Colonies* swg1stS_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr)
{
	WINDOW	wdw;

	stripe(seqs, &wdw, alprm.sh);
	Aln2s1	alnv(seqs, pwd, &wdw);
	return (alnv.fwdswgS_ng(scr));
}

SKL* swg2ndS_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr, COLONY* clny)
{
	if (clny->val <= 0) return (0);
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];

	a->left = clny->mlb;
	b->left = clny->nlb;
	b->right = clny->nrb;
	a->right = clny->mrb;
	WINDOW	wdw = {clny->upr, clny->lwr, clny->upr - clny->lwr + 3};
	Aln2s1 alnv(seqs, pwd, &wdw);
	return (alnv.globalS_ng(scr));
}

SKL* alignS_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr)
{
	WINDOW	wdw;

	stripe(seqs, &wdw, alprm.sh);
	Aln2s1 alnv(seqs, pwd, &wdw);
	return alnv.globalS_ng(scr);
}
