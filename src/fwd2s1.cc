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

#define	DEBUG	1

#include "aln.h"
#include "vmf.h"
#include "wln.h"
#include "boyer_moore.h"

#define	CigarM	0

static	const	int	NCAND = 4;
static	const	int	Newd = 8;
static	const	int	expected_max_overlap = 1024;
static	const	int	expected_overlap_ext = 16;
static	const	int	max_dist2ss = 9;

static int infer_orientation(Seq** sqs, const PwdB* pwd);

class Aln2s1 {
protected:
const	Seq**	seqs;
const	Seq*	a;
const	Seq*	b;
const 	PwdB*	pwd;
	Mfile*	mfd;
	Vmf*	vmf;
const	INT	lowestlvl;
const	int	term;
const	VTYPE	slmt;
	int	wlmt;
	bool	island;
const	int	Nod;
	VTYPE	backward[expected_max_overlap];
public:
	Aln2s1(const Seq** _seqs, const PwdB* _pwd);
	~Aln2s1() {delete mfd;}
	void	reset(const Seq** sqs) {
	    delete mfd; mfd = 0;
	    a = sqs[0];
	    b = sqs[1];
	}
	void	initS_ng(RVP* hh[], const WINDOW& wdw, CHAR* hdir, const RANGE* cutrng = 0);
	RVP*	lastS_ng(RVP* hh[], const WINDOW& wdw, const RANGE* cutrng = 0);
	VTYPE	forwardS_ng(const WINDOW& wdw, long pp[] = 0, 
		bool spj = true, bool align = true, const RANGE* cutrng = 0);
	VTYPE	shortcutS_ng(int ovr, const BOUND& bab);
	void	hinitS_ng(RVWU* hhg[], const WINDOW& wdw);
	RVWU*	hlastS_ng(RVWU* hh[], const WINDOW& wdw);
	VTYPE	hirschbergS_ng(int cpos[], 		// coordinates on center row
		const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb);
	void	sinitS_ng(VTYPE* hh[], const WINDOW& wdw);
	VTYPE	slastS_ng(VTYPE* hh[], const WINDOW& wdw);
	VTYPE	scorealoneS_ng(const WINDOW& wdw);
	void	pfinitS_ng(RVP* hh[], const WINDOW& wdw);
	void	pbinitS_ng(RVP* hh[], const WINDOW& wdw);
	VTYPE	back2ward5endS_ng(long* ptr, const WINDOW& wdw);
	VTYPE	for2ward3endS_ng(long* ptr, const WINDOW& wdw);
	VTYPE	openendS_ng(int cmode);
	void	cinitS_ng(RVDWJC* hhc[], const WINDOW& wdw);
	Colonies*	fwdswgS_ng(const WINDOW& wdw, VTYPE* scr);
	bool	indelfreespjS(int agap, VTYPE& iscr);
	VTYPE	diagonalS_ng();
	VTYPE	trcbkalignS_ng(const WINDOW& wdw, bool spj = true,
		const RANGE* cutrng = 0);
	VTYPE	ordinaryS_ng(const WINDOW& wdw);
	VTYPE	backforth(int n, const BOUND& lub);
	VTYPE	lspS_ng(const WINDOW& wdw);
	int	nearest5ss(int from, const BOUND& bab);
	int	nearest3ss(int from, const BOUND& bab);
	int	first_exon_wmm(VTYPE& maxscr);
	VTYPE	first_exon(const BOUND& bab);
	int	last_exon_wmm(VTYPE& maxscr);
	VTYPE	last_exon(const BOUND& bab);
	VTYPE	micro_exon(const BOUND& bab);
	VTYPE	interpolateS(INT level, int agap, int bgap, int cmode, 
		JUXT* wjxt, const BOUND& bab);
	VTYPE	seededS_ng(INT level, int cmode, const BOUND& lub);
	SKL*	globalS_ng(const WINDOW& wdw, VTYPE* scr);
	VTYPE	creepback(int ovr, VTYPE bscr, const BOUND& lub);
	VTYPE	creepfwrd(int& ovr, VTYPE bscr, const BOUND& lub);
};

extern	Colonies* swg1stS_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr);
extern	SKL*	swg2ndS_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr, COLONY* clny);

Aln2s1::Aln2s1(const Seq* _seqs[], const PwdB* _pwd) :
	seqs(_seqs), a(seqs[0]), b(seqs[1]), pwd(_pwd),
	mfd(0), vmf(0), lowestlvl(b->wllvl),
	term((int) ((pwd->Vthr + pwd->BasicGOP) / pwd->BasicGEP)),
	slmt(pwd->Vthr / 2), island(false), Nod(2 * pwd->Noll - 1)
{}

void Aln2s1::initS_ng(RVP* hh[], const WINDOW& wdw, CHAR* hdir, const RANGE* cutrng)
{
	int	n = b->left;
	int	r = b->left - a->left;
	int	rr = b->right - a->left;
	int	cutlen = cutrng? cutrng->right - cutrng->left: 0;

	RVP*	h = hh[0] + r;
	CHAR*	hd = hdir? hdir + r: 0;
	h->val = 0;
	if (hd) *hd = 0;
	h->ptr = vmf? vmf->add(a->left, n, 0): 0;
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
		h->ptr = vmf? vmf->add(a->left, n, 0): 0;
	    }
	}

	r = b->left - a->left;
	rr = b->left - a->right;
	h = hh[0] + r;
	if (hdir) hd = hdir + r;
	RVP*	g = hh[1] + r;
	if (wdw.lw > rr) rr = wdw.lw;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --g; if (hd) *--hd = 2;
	    if (b->inex.exgl) {		/* local */
		h->val = 0;
		h->ptr = vmf? vmf->add(a->left + i, b->left, 0): 0;
	    } else {			/* (semi) global */
		*h = h[1];
		if (i == 1) {
		    h->val += pwd->GapPenalty(1);
		    *g = *h;
		} else {
		    *g = g[1];
		    h->val += pwd->GapExtPen(i);
		    g->val += pwd->BasicGEP;
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
	    rw = min(wdw.up, b->right - a->left) - cutlen;
	    for (RVP* h = hh[0] + rw; h > h9; --h)
		if (h->val > mx->val) mx = h;
	}
	if (vmf) {
	    int	i = mx - h9;
	    int	m9 = a->right;
	    int	n9 = b->right;
	    if (i > 0) m9 -= i;
	    if (i < 0) n9 += i;
	    mx->ptr = vmf->add(m9, n9, mx->ptr);
	}
	return (mx);
}

VTYPE Aln2s1::forwardS_ng(const WINDOW& wdw, long* pp, 
	bool spj, bool align, const RANGE* cutrng)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	RVP*	hh[NOL];	// H matrix
	RVP	f1, f2;
	RVP*	hf[NOD] = {0, &f1, 0, &f2, 0};	// [DIAG, HORI, VERT, HORL, VERTL]
	RVPDJ	hl[NCAND + 1]; // [candidates]
	int	nx[NCAND + 1]; // [candidates]
	RVPDJ*	maxphl[NOD];
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
	int	cutlen = cutrng? cutrng->right - cutrng->left: 0;
	VTYPE	longgep = pwd->BasicGEP * cutlen;
	VTYPE	longgep2 = pwd->LongGEP * cutlen;

	int	width = wdw.width;
	if (cutrng) width -= cutlen;
	size_t	bufsiz = pwd->Noll * width;
	RVP*	jbuf = new RVP[bufsiz];
	vset(jbuf, black_vp, bufsiz);
	RVP*	blackvp = jbuf + bufsiz - 1;
	hh[0] = jbuf - wdw.lw + 1;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + width;
	CHAR*	dbuf = new CHAR[width];
	vclear(dbuf, width);
	CHAR*	hdir = dbuf - wdw.lw + 1;
	if (vmf) vmf->add(0, 0, 0);	/* Skip 0-th record */
	initS_ng(hh, wdw, hdir, cutrng);

	int	m = a->left;
	PfqItr	api(a, m);
	int	api_size = api.size();
	if (!a->inex.exgl) --m; /* global */
	int	n1 = m + wdw.lw - 1;
	int	n2 = m + wdw.up;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as) {
	    bool	internal = spj && (!a->inex.exgr || m < a->right);
	    ++n1; ++n2;
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    int k = 0;
	    RVP*&	h = hf[0] = hh[k] + r;
	    RVP*&	g = hf[2] = hh[++k] + r;
	    RVP*&	g2 = hf[4] = dagp? hh[++k] + r: blackvp;
	    CHAR*	dir = hdir + r;
	    f1 = f2 = black_vp;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    vset((RVPDJ*) hl, black_vpdj, NCAND + 1);
	    for (int l = 0; l <= NCAND; ++l) nx[l] = l;
	    int	ncand = -1;
#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m, n, *dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    bool	a_in_zone = api_size && (api == m);
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x;
		++dir; ++h; ++g; if (dagp) ++g2;
/*	Diagonal	*/
		RVP*	from = h;
		RVP*	mx = h;
		VTYPE	diag = h->val;
		if (m == a->left) goto HorizonS;
		h->val += qprof[*bs];
		*dir = (*dir % Newd)? Newd: 0;

/*	Vertical	*/
		x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    g->val = x;
		    g->ptr = from->ptr;
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val > mx->val) mx = g;

/*	Vertical2	*/
		if (dagp) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[1].val) {
		    g2->val = x;
		    g2->ptr = from->ptr;
		  } else	*g2 = g2[1];
		  g2->val += pwd->LongGEP;
		  if (g2->val > mx->val) mx = g2;
		}
HorizonS:
/*	Horizontal	*/
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f1.val) {
		    f1.val = x;
		    f1.ptr = h[-1].ptr;
		}
		f1.val += pwd->BasicGEP;
		if (f1.val >= mx->val) mx = &f1;

/*	Horizontal2	*/
		if (dagp) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= f2.val) {
			f2.val = x;
			f2.ptr = h[-1].ptr;
		    }
		    f2.val += pwd->LongGEP;
		    if (f2.val >= mx->val) mx = &f2;
		}

/*	intron 3' boundary, assume no overlapping signals	*/
		if (internal && b->exin->isAccpt(n)) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    vclear(maxphl, Nod);
		    for (int l = 0; l <= ncand; ++l) {
			RVPDJ*	phl = hl + nx[l];
			x = phl->val + sigJ +
			    pwd->IntPen->Penalty(n - phl->jnc) +
			    b->exin->sig53(phl->jnc, n, IE53);
			from = hf[phl->dir];
			if (x >= from->val) {
			    from->val = x;
			    maxphl[phl->dir] = phl;
			}
		    }
		    for (int d = 0; d < Nod; ++d) {
			RVPDJ*	phl = maxphl[d];
			if (!phl) continue;
			from = hf[d];
			if (vmf) from->ptr = vmf->add(m, n, 
			    vmf->add(m, phl->jnc, phl->ptr));
			if (from->val >= mx->val) {
			    mx = from;
			    *dir |= SPJCI;
			}
		    }
		}

/*	Find optimal path	*/
#if DEBUG
		VTYPE	y = h->val;
#endif
		int	hd = 0;
		if (h != mx) {		// non-diagonal
		    *h = *mx;
		    while (mx != hf[++hd]) ;
		    if (hd % 2 && !(*dir & SPJC)) *dir = dir[-1] & SPJC;
		    else	*dir &= ~15;
		    *dir += hd;
		} else if (Local && h->val > diag) {
		    if (LocalL && diag == 0 && !(*dir & SPJC))
			h->ptr = vmf? vmf->add(m - 1, n - 1, 0): 0;
		    else if (LocalR && h->val > maxh.val) {
			maxh.val = h->val;
			maxh.p = h->ptr;
			maxh.m = m;
			maxh.n = n;
		    }
		}
		if (LocalL && h->val <= 0) {h->val = 0; *dir = 1;}
		else if (vmf && align && *dir == Newd)
		    h->ptr = vmf->add(m - 1, n - 1, h->ptr);

/*	intron 5' boundary	*/
		if (internal && b->exin->isDonor(n)) {
		    VTYPE	sigJ = b->exin->sig53(n, 0, IE5);
		    for (int k = hd == 0? 0: 1; k < Nod; ++k) {
			from = hf[k];
			if (k == hd) {			// disallow orphan exon
			    int	vert = k == 2 || k == 4;
			    if (dir[vert] & SPIN) continue;
			} else if (hd >= 0) {
			    y = mx->val;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (from->val <= y) continue;	// prune
			}
			x = from->val + sigJ;
			int	l = ncand < NCAND? ++ncand: NCAND;
			while (--l >= 0) {
			    if (x > hl[nx[l]].val +
				pwd->IntPen->PenaltyDrop(n - hl[nx[l]].jnc))
				swap(nx[l], nx[l + 1]);
			    else
				break;
			}
			if (++l < NCAND) {
			    RVPDJ*	phl = hl + nx[l];
			    phl->val = x;
			    phl->jnc = n;
			    phl->dir = k;
			    phl->ptr = from->ptr;
			} else --ncand;
		    }
		}

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d %2d ", m, n, *dir);
		    putvar(mx->val); putvar(y); 
		    putvar(g->val); putvar(f1.val);
		    if (dagp) {
			putvar(g2->val); putvar(f2.val);
		    }
		    if (algmode.lsg) {
			putvar(hl[nx[0]].val);
			putvar(hl[nx[1]].val);
		    }
		    putchar('\n');
		}
#endif
		if (cutrng && n == cutrng->left) {	// shortcut
		    f1.val += longgep;
		    if (dagp) f2.val += longgep2;
		    *h = dagp? f2: f1;
		    *g = black_vp;
		    n += cutlen;
		    bs += cutlen;
		}
	    } // end of n-loop
	    if (a_in_zone) ++api;		// has exon-exon junctions
	} // end of m-loop

	if (LocalR) {
	    if (vmf) *pp = vmf->add(maxh.m, maxh.n, maxh.p);
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
		    rbuf.iscr = NEVSEL;
		    h += xi;
		    insert -= preint;
		} else {	/* gap */
		    h += x;
		}
		if (insert) {
		    if (cigar) cigar->push('D', insert);
		    if (samfm) samfm->push('D', insert);
		    if (vlgar) vlgar->push('G', 0, insert);
		    insert = intlen = preint = 0;
		}
	    }
	    int	ni = wsk->n - n;
	    if (ni && deletn) {
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
	    } else if (i < 0) {
		int	n3 = n + (i = -i);
		if (algmode.lsg && i > IntronPrm.minl) {
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
	    eijnc->unshift();
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
	fst->val += pwd->BasicGOP * fst->gap + pwd->BasicGEP * fst->unp;;
	return (h);
}

/********************************************************
*
*	unidirectional Hirschberg algorithm
*
********************************************************/

void Aln2s1::hinitS_ng(RVWU* hhg[], const WINDOW& wdw)
{
	int	r = b->left - a->left;
	int	r0 = r;
	int	rr = b->right - a->left;

	RVWU*	h = hhg[0]  + r;
	h->val = 0;
	h->lwr = h->upr = r;
	if (a->inex.exgl) {	// semi-global
	    if (wdw.up < rr) rr = wdw.up;
	    while (++r <= rr) {
		++h;
		h->val = 0;
		h->lwr = h->upr = r;
	    }
	}

	r = r0;
	rr = b->left - a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	h = hhg[0] + r;
	RVWU*	g = hhg[1] + r;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --g;
	    if (b->inex.exgl) {
		h->val = 0;
		h->lwr = h->upr = r;
	    } else {
		*h = h[1];
		if (i == 1) {
		    h->val += pwd->GapPenalty(1);
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

RVWU* Aln2s1::hlastS_ng(RVWU* hh[], const WINDOW& wdw)
{
	RVWU*	h9 = hh[0] + b->right - a->right;
	RVWU*	mx = h9;

	if (b->inex.exgr) {
	    int	rw = min(wdw.up, b->right - a->left);
	    for (RVWU* h = hh[0] + rw; h > h9; --h)
		if (h->val > mx->val) mx = h;
	}
	if (a->inex.exgr) {
	    int	rw = max(wdw.lw, b->left - a->right);
	    for (RVWU* h = hh[0] + rw; h < h9; ++h)
		if (h->val > mx->val) mx = h;
	}
	return (mx);
}


VTYPE Aln2s1::hirschbergS_ng(int cpos[], 
		const WINDOW& wdw, WINDOW& wdwf, WINDOW& wdwb)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	RVWU*	hhg[NOL];
	RVWU*	hhb[NOL];
	RVWU*	hb[NOL];		// hhb[k] + r
	RVWU*	hf[NOD];
	RVDWUJ	hl[NCAND + 1];;		// [candidates]
	int	nx[NCAND + 1];		// [candidates]
	RVDWUJ*	maxphl[NOD];
	RVWU	e[2];
	RVWU*&	f1 = hf[1] = e;
	RVWU*&	f2 = hf[3] = e + 1;
	RVWU	maxh = black_vwu;
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;
#if MONITOR
	long	start = time(0);
	long	stop;
#endif

	size_t	bufsiz = pwd->Noll * wdw.width;
	RVWU*	wbuf = new RVWU[bufsiz + bufsiz];
	vset(wbuf, black_vwu, bufsiz + wdw.width);
	vclear(hb, NOL);
	RVWU*	blackvwu = wbuf + bufsiz - 1;	// assume to be const
	hhg[0] = wbuf - wdw.lw + 1;
	CHAR*	from_vert = new CHAR[wdw.width] - wdw.lw + 1;
	vclear(from_vert + wdw.lw - 1, wdw.width);
	for (int k = 1; k < pwd->Noll; ++k) hhg[k] = hhg[k-1] + wdw.width;
	for (int k = 0; k < pwd->Noll; ++k) hhb[k] = hhg[k] + bufsiz;
	hinitS_ng(hhg, wdw);

	int	m = a->left;
	int	mm = (a->left + a->right + 1) / 2;
	PfqItr	api(a, m);		// iterator
	int	api_size = api.size();
	if (!a->inex.exgl) --m; 	// global
	int	n1 = m + wdw.lw - 1;
	int	n2 = m + wdw.up;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as) {
	    ++n1; ++n2;
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    int	lst = r;
	    int	hdir = 0;
	    RVWU*&	h = hf[0] = hhg[0] + r;
	    RVWU*&	g = hf[2] = hhg[1] + r;
	    RVWU*&	g2 = hf[4] = dagp? hhg[2] + r: blackvwu;
	    vset(e, black_vwu, 2);
	    if (m == mm)
		for (int k = 0; k < pwd->Noll; ++k)
		    hb[k] = hhb[k] + r;
	    vset(hl, black_vdwuj, NCAND + 1);
	    for (int l = 0; l <= NCAND; ++l) nx[l] = l;
	    int	ncand = -1;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    bool	a_in_zone = api_size && api == m;
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d ", m, n);
		putvar(h->val); putchar('\n');
	    }
#endif
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x, y;
		++r; ++h; ++g; if (dagp) ++g2;

/*	Diagonal	*/
		RVWU*	from = h;
		RVWU*	mx = h;
		if (m == a->left) goto HorizonF;
		h->val += qprof[*bs];

/*	Vertical	*/
		x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    *g = *from;
		    g->val = x;
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;

/*	Vertical2	*/
		if (dagp) {
		  x = from->val + pwd->LongGOP;
		  if (x >= g2[1].val) {
		    *g2 = *from;
		    g2->val = x;
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
		f1->val += pwd->BasicGEP;
		if (f1->val >= mx->val) mx = f1;

/*	Horizontal2	*/
		if (dagp) {
		    x = h[-1].val + pwd->LongGOP;
		    if (x >= f2->val) {
			*f2 = h[-1];
			f2->val = x;
		    }
		    f2->val += pwd->LongGEP;
		    if (f2->val >= mx->val) mx = f2;
		}

/*	intron 3' boundary	*/
		if (b->exin->isAccpt(n)) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    vclear(maxphl, Nod);
		    for (int l = 0; l <= ncand; ++l) {
			RVDWUJ*	phl = hl + nx[l];
			from = hf[phl->dir];
			x = phl->val + sigJ + b->exin->sig53(phl->jnc, n, IE53) +
			    pwd->IntPen->Penalty(n - phl->jnc);
			if (x > from->val) {
			    from->val = x;
			    maxphl[phl->dir] = phl;
			}
		    }
		    for (int k = 0; k < Nod; ++k) {
			RVDWUJ*	phl = maxphl[k];
			if (!phl) continue;
			from = hf[k];
			from->upr = max(phl->upr, r);
			from->lwr = min(phl->lwr, r);
			if (m > mm) from->ulk = phl->ulk;
			else if (m == mm) 
			    from->ulk = phl->jnc - m;
			if (from->val > mx->val) {
			    mx = from;
			    hdir = SPJCI;
			}
		    }
		}

/*	Find optimal path	*/
		y = h->val;
		if (h != mx) {
		    *h = *mx;
		    if (h->upr < r) h->upr = r;
		    if (h->lwr > r) h->lwr = r;
		} else if (LocalR && y > maxh.val)
		    maxh = *h;
		if (LocalL && h->val <= 0) h->val = 0;
		int		hd = 0;
		for ( ; mx != hf[hd]; ++hd) ;

/*	intron 5' boundary	*/
		if (b->exin->isDonor(n)) {
		    VTYPE	sigJ = b->exin->sig53(n, 0, IE5);
		    for (int k = mx == h? 0: 1; k < Nod; ++k) {
			from = hf[k];
			if (k % 2 && (hdir & SPIN)) continue;
			if (k != hd) {
			    y = mx->val;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (from->val <= y) continue;	// prune
			}
			x = from->val + sigJ;
			int	l = ncand < NCAND? ++ncand: NCAND;
			while (--l >= 0) {
			    if (x >= hl[nx[l]].val)
				swap(nx[l], nx[l + 1]);
			    else
				break;
			}
			if (++l < NCAND) {
			    RVDWUJ*	phl = hl + nx[l];
			    phl->val = x;
			    phl->jnc = n;
			    phl->dir = k;
			    phl->upr = from->upr;
			    phl->lwr = from->lwr;
			    if (m > mm) phl->ulk = from->ulk;
			    else if (m == mm)	// hori?
				phl->ulk = (k % 2)? lst: end_of_ulk;
			} else --ncand;
		    }
		}

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d ", m, n);
		    putvar(mx->val); putvar(y); 
		    putvar(g->val); putvar(f1->val);
		    if (dagp) {
			putvar(g2->val); putvar(f2->val);
		    }
		    if (algmode.lsg) {
			putvar(hl[nx[0]].val);
			putvar(hl[nx[1]].val);
			putvar(hl[nx[2]].val);
		    }
		    putchar('\n');
		}
#endif

// save variables at the center
		if (m == mm) {
		    for (int k = 0; k < pwd->Noll; ++k)
			*++hb[k] = *hf[2 * k];
		    if (mx == h || (hdir & SPJC)) lst = r;
		    else if (lst < r) hb[0]->ulk = lst;
		    else from_vert[r] = 1;
		    for (int k = 0, j = 0; k < pwd->Noll; ++k, j += wdw.width) {
			hf[2 * k]->upr = hf[2 * k]->lwr = r;
			hf[2 * k]->ulk = r + j;
		    }
		}	// was center
		hdir &= (hd % 2? SPIN: 0);
	    }	// end of n-loop
	    if (a_in_zone) ++api;		// has exon-exon junctions
	}	// end of m-loop

	if (!LocalR) maxh = *hlastS_ng(hhg, wdw);
	wdwb.up = maxh.upr;
	wdwb.lw = maxh.lwr;
	wdwb.width = wdwb.up - wdwb.lw + 3;
	int	r = maxh.ulk;
	int	k = (r - wdw.lw) / wdw.width;
	if (k)	r -= k * wdw.width;
	RVWU*	mx = hhb[k] + r;
	b->inex.exgr = k? 1: 0;
	b->inex.exgl = (k == from_vert[r])? 1: 0;
	int	c = 0;
	cpos[c++] = mm;
	while (mx->ulk != end_of_ulk) {
	    cpos[c++] = r + mm;
	    mx = hhb[0] + (r = mx->ulk);
	}
	cpos[c++] = r + mm;
	cpos[c] = end_of_ulk;
	wdwf.up = mx->upr;
	wdwf.lw = mx->lwr;
	wdwf.width = wdwf.up - wdwf.lw + 3;

	delete[] wbuf;
	delete[] (from_vert + wdw.lw - 1);
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
	VTYPE*	g = hh[1] + r;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --g;
	    *h = h[1];
	    if (i == 1) {
		*h += pwd->GapPenalty(1);
		*g = *h;
	    } else {
		*g = g[1];
		*h += pwd->GapExtPen(i);
		*g += pwd->BasicGEP;
	    }
	}
}

VTYPE Aln2s1::slastS_ng(VTYPE* hh[], const WINDOW& wdw)
{
	VTYPE*	h9 = hh[0] + b->right - a->right;
	VTYPE	mx = *h9;

	if (b->inex.exgr) {
	    int	rw = min(wdw.up, b->right - a->left);
	    for (VTYPE* h = hh[0] + rw; h > h9; --h)
		if (*h > mx) mx = *h;
	}
	if (a->inex.exgr) {
	    int	rw = max(wdw.lw, b->left - a->right);
	    for (VTYPE* h = hh[0] + rw; h < h9; ++h)
		if (*h > mx) mx = *h;
	}
	return (mx);
}

VTYPE Aln2s1::scorealoneS_ng(const WINDOW& wdw)
{
const	bool	dagp = pwd->Noll == 3;	// double affine gap penalty
	VTYPE*	hh[NOL];
	RVDJ	hl[NCAND + 1];
	int	nx[NCAND + 1];		// [candidates]
	RVDJ*	maxphl[NOD];
	VTYPE	f1;
	VTYPE	f2;
	VTYPE*	hf[NOD] = {0, &f1, 0, &f2, 0};
	VTYPE	maxh = NEVSEL;
	int	Local = algmode.lcl & 16;
	int	LocalL = Local && a->inex.exgl && b->inex.exgl;
	int	LocalR = Local && a->inex.exgr && b->inex.exgr;

	size_t	bufsiz = pwd->Noll * wdw.width;
	VTYPE*	wbuf = new VTYPE[bufsiz];
	vset(wbuf, NEVSEL, bufsiz);
	VTYPE*	blackv = wbuf + bufsiz - 1;	// assume to be const
	hh[0] = wbuf - wdw.lw + 1;
	for (int k = 1; k < pwd->Noll; ++k) hh[k] = hh[k-1] + wdw.width;
	sinitS_ng(hh, wdw);

	int	m = a->left;
	if (!a->inex.exgl) --m; 	// global
	int	n1 = m + wdw.lw - 1;
	int	n2 = m + wdw.up;
const	CHAR*	as = a->at(m);
	for ( ; ++m <= a->right; ++as) {
	    ++n1; ++n2;
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
const	    CHAR*	bs = b->at(n);
	    int	r = n - m;
	    VTYPE*&	h = hf[0] = hh[0] + r;
	    VTYPE*&	g = hf[2] = hh[1] + r;
	    VTYPE*&	g2 = hf[4] = dagp? hh[2] + r: blackv;
	    f1 = f2 = NEVSEL;
	    vset(hl, black_vdj, NCAND + 1);
	    for (int l = 0; l <= NCAND; ++l) nx[l] = l;
	    int	ncand = -1;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
#if DEBUG
	    if (OutPrm.debug) {
		printf("%2d %2d ", m, n);
		putvar(*h); putchar('\n');
	    }
#endif
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x, y;
		++h; ++g; if (dagp) ++g2;

/*	Diagonal	*/
		VTYPE*	from = h;
		VTYPE*	mx = h;
		if (m == a->left) goto HorizonF;
		*h += qprof[*bs];

/*	Vertical	*/
		x = *++from + pwd->BasicGOP;
		*g = max(x, g[1]) + pwd->BasicGEP;
		if (*g > *mx) mx = g;

/*	Vertical2	*/
		if (dagp) {
		  x = *from + pwd->LongGOP;
		  *g2 = max(x, g2[1]) + pwd->LongGEP;
		  if (*g2 > *mx) mx = g2;
		}
HorizonF:
/*	Horizontal	*/
		x = h[-1] + pwd->BasicGOP;
		f1 = max(x, f1) + pwd->BasicGEP;
		if (f1 > *mx) mx = &f1;

/*	Horizontal2	*/
		if (dagp) {
		    x = h[-1] + pwd->LongGOP;
		    f2 = max(x, f2) + pwd->LongGEP;
		    if (f2 > *mx) mx = &f2;
		}

/*	intron 3' boundary	*/
		if (b->exin->isAccpt(n)) {
		    vclear(maxphl, Nod);
		    for (int l = 0; l <= ncand; ++l) {
			RVDJ*	phl = hl + nx[l];
			from = hf[phl->dir];
			x = phl->val + b->exin->sig53(phl->jnc, n, IE53) +
			    pwd->IntPen->Penalty(n - phl->jnc);
			if (x > *from) {
			    *from = x;
			    maxphl[phl->dir] = phl;
			}
		    }
		    for (int k = 0; k < Nod; ++k) {
			RVDJ*	phl = maxphl[k];
			if (!phl) continue;
			from = hf[k];
			if (*from > *mx) mx = from;
		    }
		}

/*	Find optimal path	*/
		y = *h;
		if (h != mx) *h = *mx;
		else if (LocalR && y > maxh) maxh = y;
		if (LocalL && *h < 0) *h = 0;
		int	hd = 0;
		for ( ; mx != hf[hd]; ++hd) ;

/*	intron 5' boundary	*/
		if (b->exin->isDonor(n)) {
		    VTYPE	sigJ = b->exin->sig53(n, 0, IE5);
		    for (int k = hd == 0? 0: 1; k < Nod; ++k) {
			from = hf[k];
			if (k != hd) {
			    y = *mx;
			    if (hd == 0 || (k - hd) % 2) y += pwd->GOP[k / 2];
			    if (*from <= y) continue;	// prune
			}
			x = *from + sigJ;
			int	l = ncand < NCAND? ++ncand: NCAND;
			while (--l >= 0) {
			    if (x >= hl[nx[l]].val)
				swap(nx[l], nx[l + 1]);
			    else
				break;
			}
			if (++l < NCAND) {
			    RVDJ*	phl = hl + nx[l];
			    phl->val = x;
			    phl->jnc = n;
			    phl->dir = k;
			} else --ncand;
		    }
		}

#if DEBUG
		if (OutPrm.debug) {
		    printf("%2d %2d ", m, n);
		    putvar(*mx); putvar(y); 
		    putvar(*g); putvar(f1);
		    if (dagp) {
			putvar(*g2); putvar(f2);
		    }
		    if (algmode.lsg) {
			putvar(hl[0].val);
			putvar(hl[1].val);
			putvar(hl[2].val);
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
	RVP*	g = hh[1] + r;
	h->val = 0;
	h->ptr = vmf? vmf->add(a->left, b->left, 0): 0;

	r = b->left - a->left;
	rr = b->left - a->right;
	if (wdw.lw > rr) rr = wdw.lw;
	for (int i = 1; --r >= rr; ++i) {
	    --h; --g;
	    *h = h[1];
	    if (i == 1)	h->val += pwd->BasicGOP;
	    h->val += pwd->BasicGEP;
	    *g = *h;
	}
}

void Aln2s1::pbinitS_ng(RVP* hh[], const WINDOW& wdw)
{
	int	r = b->right - a->right;
	int	rr = b->left - a->right;

	RVP*	h = hh[0] + r;
	RVP*	g = hh[1] + r;
	h->val = 0;
	h->ptr = vmf? vmf->add(a->right, b->right, 0): 0;

	r = b->right - a->right;
	rr = b->right - a->left;
	if (wdw.up < rr) rr = wdw.up;
	for (int i = 1; ++r <= rr; ++i) {
	    ++h; ++g;
	    *h = h[-1];
	    if (i == 1) h->val += pwd->BasicGOP;
	    h->val += pwd->BasicGEP;
	    *g = *h;
	}
}

// intron-less backward extension

VTYPE Aln2s1::back2ward5endS_ng(long *ptr, const WINDOW& wdw)
{
	RVP*	hh[2];			// H matrix
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	VSKLP	lsth = maxh;
	RVP	f1 = black_vp;
	RVP*	hf[3] = {0, &f1, 0};
	RVP*&	h = hf[0];
	RVP*&	g = hf[2];

	size_t	bufsiz = 2 * wdw.width;
	RVP*	jbuf = new RVP[bufsiz];
	vset(jbuf, black_vp, bufsiz);
	RVP*	blackvp = jbuf + bufsiz - 1;
	hh[0] = jbuf - wdw.lw + 1;
	hh[1] = hh[0] + wdw.width;
	CHAR*	dbuf = new CHAR[wdw.width];
	vset(dbuf, CHAR(1), wdw.width);
	CHAR*	hdir = dbuf  - wdw.lw + 1;
	vmf->add(0, 0, 0);		// Skip 0-th record
	pbinitS_ng(hh, wdw);

	int	m = a->right;
	if (!a->inex.exgr) ++m; 	// global
const	CHAR*	as = a->at(m);
	int	n = 0;
	int	n9 = 0;
	int	n1 = m + wdw.lw;
	int	n2 = m + wdw.up + 1;
	bool	enda = a->inex.exgl && algmode.lcl < 16;
	bool	endb = b->inex.exgl && algmode.lcl < 16;
	while (--m >= a->left) {
	    --as; --n1; --n2;
	    n = min(n2, b->right);
	    n9 = max(n1, b->left);
	    int	r = n - m;
	    int	nr = n + 1;
	    int	peak = 0;
const	    CHAR*	bs = b->at(n);
	    h = hh[0] + r;
	    g = hh[1] + r;
	    CHAR*	dir = hdir + r;
	    f1 = black_vp;
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
		--bs; --r; --h; --g;

/*	Diagonal	*/
		RVP*	from = h;
		RVP*	mx = h;
		if (m == a->right) goto HorizonB;
		h->val += qprof[*bs];
		*dir = *dir % Newd? Newd: 0;

/*	Vertical	*/
		x = (--from)->val + pwd->BasicGOP;
		if (x >= g[-1].val) {
		    *g = *from;
		    g->val = x;
		} else *g = g[-1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;
HorizonB:
/*	Horizontal	*/
		x = h[1].val + pwd->BasicGOP;
		if (x >= f1.val) {
		    f1 = h[1];
		    f1.val = x;
		}
		f1.val += pwd->BasicGEP;
		if (f1.val >= mx->val) mx = &f1;

/*	Find optimal path	*/
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
		putvar(g->val); putvar(f1.val);
		putchar('\n');
	}
#endif
		if (m == a->left && enda && h->val > lsth.val) {
		    lsth.val = h->val;
		    lsth.p = h->ptr;
		    lsth.m = m;
		    lsth.n = n;
		}
	    } // end of n loop
	    if (mxd->val <= NEVSEL) break;	// no peak
	    if (peak) n1 = n;
	    if (++n == b->left && endb && h->val > lsth.val) {
		lsth.val = h->val;
		lsth.p = h->ptr;
		lsth.m = m;
		lsth.n = n;
	    }
	} // end of m loop

	if (!a->inex.exgl && !b->inex.exgl && h) {		// global
	    lsth.val = h->val;
	    lsth.p = h->ptr;
	    lsth.m = ++m;
	    lsth.n = max(n, n9);
	} else if (lsth.val == NEVSEL) lsth = maxh;
	if (vmf) *ptr = vmf->add(lsth.m, lsth.n, lsth.p);
	delete[] jbuf;
	delete[] dbuf;
	return (lsth.val);
}

// intron-less forward extension

VTYPE Aln2s1::for2ward3endS_ng(long *ptr, const WINDOW& wdw)
{
	RVP*	hh[2];		// H matrix
	VSKLP	maxh = {NEVSEL, a->left, b->left, 0};
	VSKLP	lsth = maxh;
	RVP	f1 = black_vp;
	RVP*	hf[3] = {0, &f1, 0};
	RVP*&	h = hf[0];
	RVP*&	g = hf[2];

	size_t	bufsiz = 2 * wdw.width;
	RVP*	jbuf = new RVP[bufsiz];
	vset(jbuf, black_vp, bufsiz);
	RVP*	blackvp = jbuf + bufsiz - 1;
	hh[0] = jbuf - wdw.lw + 1;
	hh[1] = hh[0] + wdw.width;
	CHAR*	dbuf = new CHAR[wdw.width];
	vset(dbuf, CHAR(1), wdw.width);
	CHAR*	hdir = dbuf  - wdw.lw + 1;
	vmf->add(0, 0, 0);		// Skip 0-th record
	pbinitS_ng(hh, wdw);

	int	m = a->left;
	if (!a->inex.exgl) --m;		// global
	int	n1 = m + wdw.lw - 1;
	int	n2 = m + wdw.up;
	int	n = 0;
	int	n9 = 0;
const	CHAR*	as = a->at(m);
	bool	enda = a->inex.exgr && algmode.lcl < 16;
	bool	endb = b->inex.exgr && algmode.lcl < 16;
	for ( ; ++m <= a->right; ++as) {
	    ++n1; ++n2;
	    n = max(n1, b->left);
	    n9 = min(n2, b->right);
	    int	nr = n - 1;
	    int	peak = 0;
	    int	r = n - m;
const	    CHAR*	bs = b->at(n);
	    h = hh[0] + r;;
	    g = hh[1] + r;
	    CHAR*	dir = hdir + r;
	    RVP*	mxd = ((h->val + pwd->Vthr) < maxh.val)? blackvp: h;
	    f1 = black_vp;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)

#if DEBUG
	if (OutPrm.debug) {
	    printf("%2d %2d %2d", m, n, *dir);
	    putvar(h->val); putchar('\n');
	}
#endif
	    for ( ; ++n <= n9; ++bs) {
		VTYPE	x;
		++r; ++h; ++g;

/*	Diagonal	*/
		RVP*	from = h;
		RVP*	mx = h;
		if (m == a->left) goto HorizonP;
		h->val += qprof[*bs];
		*dir = (*dir % Newd)? Newd: 0;

/*	Vertical	*/
		x = (++from)->val + pwd->BasicGOP;
		if (x >= g[1].val) {
		    *g = *from;
		    g->val = x;
		} else	*g = g[1];
		g->val += pwd->BasicGEP;
		if (g->val >= mx->val) mx = g;
HorizonP:
/*	Horizontal	*/
		x = h[-1].val + pwd->BasicGOP;
		if (x >= f1.val) {
		    f1 = h[-1];
		    f1.val = x;
		}
		f1.val += pwd->BasicGEP;
		if (f1.val >= mx->val) mx = &f1;
/*	Find optimal path	*/
#if DEBUG
		VTYPE	y = h->val;
#endif
		if ((*dir & Newd) && vmf)
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
		putvar(g->val); putvar(f1.val);
		putchar('\n');
	}
#endif
		if (m == a->right && enda && h->val > lsth.val) {
		    lsth.val = h->val;
		    lsth.p = h->ptr;
		    lsth.m = m;
		    lsth.n = n;
		}
	    } // end of n-loop
	    if (peak) n2 = n;
	    if (mxd->val <= NEVSEL) break;
	    if (--n == b->right && endb && h->val > lsth.val) {
		lsth.val = h->val;
		lsth.p = h->ptr;
		lsth.m = m;
		lsth.n = n;
	    }
	} // end of m-loop

	if (!a->inex.exgr && !b->inex.exgr && h) {	// global
	    lsth.val = h->val;
	    lsth.p = h->ptr;
	    lsth.m = --m;
	    lsth.n = min(n, n9);
	} else if (lsth.val == NEVSEL) lsth = maxh;	// local
	*ptr = vmf->add(lsth.m, lsth.n, lsth.p);
	delete[] jbuf;
	delete[] dbuf;
	return (lsth.val);
}

void Aln2s1::cinitS_ng(RVDWJC* hhc[], const WINDOW& wdw)
{
	int	n = b->left;
	int	r = b->left - a->left;
	int	rr = b->right - a->left;
	RVDWJC*	h = hhc[0] + r;
	RVDWJC*	g = hhc[1] + r;

	if (wdw.up < rr) rr = wdw.up;
	for ( ; r <= rr; ++h, ++r) {
	    h->val = 0;
	    h->dir = DEAD;
	    h->upr = h->lwr = r;
	    h->jnc = n;
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
	    h->jnc = b->left;
	    *--g = *h;
	}
}

Colonies* Aln2s1::fwdswgS_ng(const WINDOW& wdw, VTYPE* scr)
{
	RVDWJC*	hhc[NOL];	/* local H matrix	*/
	RVDWJC*	hf[NOL]; 	/* [DIAG, HORI, HORL] */
	RVDWJC	f[NOL]; 	/* [DIAG, HORI, HORL] */
	RVDWJC	hl[NOL][INTR+1];
	int	nx[NOL][INTR+1]; /* [DIA, HORI, HORL][candidates] */
	RVDWJC*	g2 = 0;
	Colonies*	cl = new Colonies(0);
	COLONY*	clny = cl->at();

	RVDWJC*	cbuf = new RVDWJC[pwd->Noll * wdw.width];
	vset(cbuf, black_vdwjc, pwd->Noll * wdw.width);
	hhc[0] = cbuf - wdw.lw + 1; /* +1: sp junc */
	for (int k = 1; k < pwd->Noll; ++k) hhc[k] = hhc[k-1] + wdw.width;
	cinitS_ng(hhc, wdw);
	for (int k = 1; k < pwd->Noll; ++k) hf[k] = f + k;
	int	m = a->left;
const	CHAR*	as = a->at(m);
	PfqItr	api(a, m);
	int	api_size = api.size();
	int	n1 = m + wdw.lw - 1;
	int	n2 = m + wdw.up;
	for ( ; ++m <= a->right; ++as) {
	    ++n1; ++n2;
	    int	n = max(n1, b->left);
	    int	n9 = min(n2, b->right);
	    int	r = n - m;
const	    CHAR*	bs = b->at(n);
	    RVDWJC*	h = hhc[0] + r;
	    RVDWJC*	 g = hhc[1] + r;
	    if (pwd->Noll == 3) g2 = hhc[2] + r;
const	    VTYPE*	qprof = pwd->simmtx->mtx[*as];		// sim2(as, .)
	    vset((RVDWJC*) hl, black_vdwjc, NOL * (INTR + 1));
	    for (int k = 0; k < pwd->Noll; ++k) {
		f[k] = black_vdwjc;
		for (int l = 0; l <= INTR; ++l)
		    nx[k][l] = l;
	    }

#if DEBUG
	if (OutPrm.debug) {
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
		RVDWJC*	from = h;
		RVDWJC*	mx = h;
		f[0] = *h;	/* diag */
		h->val += qprof[*bs];
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
		if (b->exin->isAccpt(n)) {
		    VTYPE	sigJ = a_in_zone? api.match_score(m): 0;
		    for (int k = 0; k < pwd->Noll; ++k) {
			from = hf[k];
			RVDWJC*	maxphl = 0;
			int*	pnl = nx[k];
			VTYPE	y = NEVSEL;
			for (int l = 0; l < INTR; ++l) {
			    RVDWJC*	phl = hl[k] + pnl[l];
			    if (!phl->dir) break;
			    x = sigJ + phl->val + pwd->IntPen->Penalty(n - phl->jnc)
				+ b->exin->sig53(phl->jnc, n, IE53);
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

		    }
		}

/*	intron 5' boundary	*/

		if (b->exin->isDonor(n)) {
		    VTYPE	sigJ = b->exin->sig53(n, 0, IE5);
		    for (int k = 0; k < pwd->Noll; ++k) {
				/* An orphan exon is disallowed */
			from = hf[k];
			if (!from->dir || (from->dir & SPIN) || 
				from->dir == VERT || from->dir == VERL)
			    continue;
			x = from->val + sigJ;
			RVDWJC*	phl = hl[k];
			int*	pnl = nx[k];
			int	l = INTR;
			while (--l >= 0) {
			    if (x > phl[pnl[l]].val)
				swap(pnl[l], pnl[l + 1]);
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
		    }
		}

/*	Find optimal path	*/
#if DEBUG
		VTYPE	y = h->val;
#endif
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
	if (OutPrm.debug) {
		printf("%2d %2d %2d ", m, n, mx->dir);
		putvar(mx->val); putvar(y); 
		putvar(g->val); putvar(f[1].val);
		if (g2) {
		    putvar(g2->val); putvar(f[2].val);
		}
		if (algmode.lsg) {
		    putvar(hl[0][nx[0][0]].val); 
		    putvar(hl[1][nx[1][0]].val);
		    if (g2) putvar(hl[2][nx[2][0]].val);
		}
		putchar('\n');
	}
#endif

	    }
	    if (a_in_zone) ++api;		// has exon-exon junctions
	}

	delete[] cbuf;
	*scr = clny->val;
	cl->sortcolonies();
	return (cl);
}

VTYPE Aln2s1::diagonalS_ng()
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
	wskl.m = mL;
	wskl.n = mL + m;
	mfd->write((UPTR) &wskl);
	wskl.m = mR;
	wskl.n = mR + m;
	mfd->write((UPTR) &wskl);
	return (LocalR? maxh: scr);
}

VTYPE Aln2s1::trcbkalignS_ng(const WINDOW& wdw, bool spj, const RANGE* cutrng)
{
	long	ptr = 0;

	vmf = new Vmf();
	VTYPE	scr = forwardS_ng(wdw, &ptr, spj, true, cutrng);
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

VTYPE Aln2s1::lspS_ng(const WINDOW& wdw)  /* recursive */ 
{
	if (wdw.up == wdw.lw) return(diagonalS_ng());
	int	m = a->right - a->left;
	int	n = b->right - b->left;
	float	cvol = wdw.lw - b->left + a->right;
	cvol = (float) m * n - cvol * cvol / 2;
	if (cvol < MaxVmfSpace || m == 1 || n <= 1)
	    return (trcbkalignS_ng(wdw, b->inex.intr));

	WINDOW	wdwl;
	WINDOW	wdwr;
	RANGE	rng[2];			/* reserve */
	INT	aexg = a->inex.exgr;	/* reserve */
	INT	bexg = b->inex.exgr;	/* reserve */
	int	cpos[6];

	save_range(seqs, rng, 2);
	VTYPE	scr = hirschbergS_ng(cpos, wdw, wdwl, wdwr);
	SKL	wskl = {cpos[0], 0};
	int	c = 1;
	for ( ; cpos[c] != end_of_ulk; ++c) {
	    wskl.n = cpos[c];
	    mfd->write((UPTR) &wskl);
	}
	INT	bexgl = b->inex.exgl;
	a->inex.exgr = b->inex.exgr = 0;
	int	r = b->left - a->left;
	if (r < wdwl.lw) b->left = a->left + wdwl.lw;
	if (r > wdwl.lw) a->left = b->left - wdwl.lw;
	a->right = cpos[0];
	b->right = cpos[--c];
	lspS_ng(wdwl);
	rest_range(seqs, rng, 2);
	a->inex.exgr = aexg;
	b->inex.exgr = bexg;
	aexg = a->inex.exgl;
	bexg = b->inex.exgl;
	a->inex.exgl = 0;
	b->inex.exgl = bexgl;
	r = b->right - a->right;
	if (r < wdwr.up) a->right = b->right - wdwr.up;
	if (r > wdwr.up) b->right = a->right + wdwr.up;
	a->left = cpos[0];
	b->left = cpos[1];
	lspS_ng(wdwr);
	rest_range(seqs, rng, 2);
	a->inex.exgl = aexg;
	b->inex.exgl = bexg;
	return scr;
}

VTYPE Aln2s1::shortcutS_ng(int ovr, const BOUND& bab)
{
	int	margin = IntronPrm.minl;
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
	int	minsh = alen - margin;
	if (sh < minsh) sh = minsh;
	WINDOW	wdw;
	stripe(seqs, &wdw, sh);

	int	interval = b->right - b->left - 2 * margin;
	RANGE	rng = {b->left + margin, b->right - margin};

	int	aexg = a->inex.exgl;
	int	bexg = b->inex.exgl;
	a->inex.exgl = b->inex.exgl = 0;	// global
	scr += trcbkalignS_ng(wdw, true, (interval > 0)? &rng: 0);
	a->inex.exgr = aexg;
	b->inex.exgr = bexg;
	return (scr);
}

VTYPE Aln2s1::openendS_ng(int cmode)
{
	vmf = new Vmf();
	long	ptr;
	WINDOW	wdw;
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

/**********************************************************
   find spj within indel-free overlapped region
   doubly conted alignment score is cancelled
   assume agap <= 0
   return whether a cannonical or reverse cannonical pair 
   is found or not
**********************************************************/

bool Aln2s1::indelfreespjS(int agap, VTYPE& iscr)
{
	int	ilen = b->right - b->left - agap;	// intron length
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
	  bool	usespb = api.size() && use_spb();
	  bs = ds;
	  for (VTYPE v = i = 0; n <= b->left; ++n, ++i, ++m) {
	    bool	a_in_zone = usespb && (api == m);
	    int	rc = b->exin->isCanon(n, n + ilen);
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
	if (iscr > NEVSEL) {
	    mfd->write(&skl);
	    EXIN*	bd = b->exin->score(skl.n);
	    bd->phs5 = 0;
	    skl.n += ilen;
	    mfd->write(&skl);
	    bd = b->exin->score(skl.n);
	    bd->phs3 = 0;
	    iscr += pwd->IntPen->Penalty(ilen);
	    return (true);
	} else
	    return (false);
}

VTYPE Aln2s1::creepback(int ovr, VTYPE bscr, const BOUND& lub)
{
const	CHAR*	as = a->at(a->left);
const	CHAR*	bs = b->at(b->left);
	VTYPE	dscr = 0;
	while (a->left > lub.la && b->left > lub.lb
		&& (ovr < 0 || dscr < bscr)) {
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
		&& (ovr < 0 || dscr < bscr)) {
	    dscr += pwd->sim2(as++, bs++);
	    ++a->right; ++b->right;
	    if (++ovr == 0) bscr += dscr;
	}
	return (dscr);
}

//	heuristic search for short exon

int Aln2s1::nearest5ss(int from, const BOUND& bab)
{
	int	retry = 0;
ss5:
	int	pa = a->left;
	int	nu = from;
	int	ta = max(bab.la, pa - max_dist2ss);
const	EXIN*	bbu = b->exin->score(nu);
	for ( ; pa > ta && nu > bab.lb; --nu, --bbu, --pa) 
	    if (bbu->sig5 > 0 || (retry && bbu->phs5 == 0)) break;
	if (nu == from) return (nu);
	pa = a->left;
	int	nd = from;
	ta = min(pa + max_dist2ss, bab.ua);
const	EXIN*	bbd = b->exin->score(nd);
	while (++pa < ta && ++nd < bab.ub)
	    if ((++bbd)->sig5 > 0 || (retry && bbd->phs5 == 0)) break;
	if (retry++ == 0 && bbu->sig5 <= 0 && bbd->sig5 <= 0) goto ss5;
	else if (bbu->phs5 || bbd->phs5) return (-1);
	if ((from - nu) == (nd - from)) return ((bbu->sig5 > bbd->sig5)? nu: nd);
	return ((from - nu < nd - from)? nu: nd);
}

int Aln2s1::nearest3ss(int from, const BOUND& bab)
{
	int	retry = 0;
ss3:
	int	pa = a->right;
	int	nu = from;
	int	ta = max(bab.la, pa - max_dist2ss);
const	EXIN*	bbu = b->exin->score(nu);
	for ( ; pa > ta && nu > bab.lb; --nu, --bbu, --pa)
	    if (bbu->sig3 > 0 || (retry && bbu->phs3 == 0)) break;
	if (nu == from) return (nu);
	pa = a->right;
	int	nd = from;
	ta = min(pa + max_dist2ss, bab.ua);
const	EXIN*	bbd = b->exin->score(nd);
	while (++pa < ta && ++nd < bab.ub)
	    if ((++bbd)->sig3 > 0 || (retry && bbd->phs3 == 0)) break;
	if (retry++ == 0 && bbu->sig3 <= 0 && bbd->sig3 <= 0) goto ss3;
	else if (bbu->phs3 || bbd->phs3) return (-1);
	if ((from - nu) == (nd - from)) return ((bbu->sig3 > bbd->sig3)? nu: nd);
	return ((from - nu < nd - from)? nu: nd);
}

// no gap with mismatches

VTYPE Aln2s1::micro_exon(const BOUND& bab)
{
	int	l = nearest5ss(b->left, bab);
	if (l < 0) return (NEVSEL);
	int	r = nearest3ss(b->right, bab);
	if (r < 0) return (NEVSEL);
	int	d5 = l - b->left;
	int	d3 = r - b->right;
	int	alen = r - l;
	b->left = l;
	a->left += d5;
	b->right = r;
	a->right += d3;
	VTYPE	maxscr = pwd->IntPen->Penalty(r - l);
	int	f = b->left;
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
	int	n = b->left + IntronPrm.minl;
	int	n9 = b->right - alen - IntronPrm.minl;
	EXIN*	bb5 = b->exin->score(n);		// 5' end of internal exon
	EXIN*	bb3 = b->exin->score(n + alen);	// 3' end of internal exon
const	CHAR*	ts = a->at(a->right);
const	CHAR*	ss = a->at(a->left);
	for ( ; n < n9; ++n, ++bb5, ++bb3) {
	    if (bb5->phs3 || bb3->phs5) continue;
	    VTYPE	scr = bb5->sig3 + bb3->sig5 +
		pwd->IntPen->Penalty(n - l) +
		pwd->IntPen->Penalty(r - n);
const	    CHAR*	as = ss;
const	    CHAR*	bs = b->at(n);
	    while (as < ts)
		scr += pwd->simmtx->mtx[*as++][*bs++];
	    if (scr > maxscr) {
		maxscr = scr;
		f = n;
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
	int	r = b->right;
	int	rr = b->right - a->right;
	int	n = rr - IntronPrm.minl;
	int	nd = n + a->right;
	EXIN*	bbi = b->exin->score(n);
	EXIN*	bb5 = b->exin->score(nd);
const	CHAR*	ts = a->at(a->right);
const	CHAR*	ss = a->at(a->left);
	VTYPE	pmch = 0;
	for (const CHAR* as = ss; as < ts; ++as)
	    pmch += pwd->simmtx->mtx[*as][*as];
	pmch *= alprm2.w;
	int	f = -1;
	for ( ; n >= b->left; --n, --nd, --bbi, --bb5) {
	    if (!b->exin->isCanon(nd, r)) continue;
	    VTYPE	ip = pwd->IntPen->Penalty(r - nd);
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
	int	r = nearest3ss(b->right, bab);
	if (r < 0) {		// no acceptor site nearby
	    b->left = max(b->right - a->right - wlmt, 0);
	    return (NEVSEL);
	}
	int	d = r - b->right;
	b->right = r;
	a->right += d;
	if (a->right == 0 || b->right == 0) {
	    SKL	skl = {a->right, b->right};
	    mfd->write(&skl);
	    return (0);
	}
	BoyerMoore	bm(b, a, -1);
	VTYPE	maxscr = NEVSEL;
	int	maxf = -1;
	while (!bm.finished()) {
	    int	f = bm.nexthit();
	    if (f >= 0) {
		int	nd = f + a->right;
		if (b->exin->isCanon(nd, r)) {
		    VTYPE	scr = pwd->IntPen->Penalty(r - nd) +
			b->exin->sig53(nd, r, IE5P3);
		    if (scr > maxscr) {
			maxscr = scr;
			maxf = f;
		    }
		}
	    }
	}

	if (maxf < 0) {
	    if ((maxf = first_exon_wmm(maxscr)) < 0) return (NEVSEL);
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
	int	l = b->left;
	int	alen = a->right - a->left;
	int	rr = b->right - alen;
	int	n = b->left + IntronPrm.minl;
	EXIN*	bb3 = b->exin->score(n);
const	CHAR*	ts = a->at(a->right);
const	CHAR*	ss = a->at(a->left);
	VTYPE	pmch = 0;			// perfect match score
	for (const CHAR* as = ss; as < ts; ++as)
	    pmch += pwd->simmtx->mtx[*as][*as];
	pmch *= alprm2.w;
	int	f = -1;
	for ( ; n < rr; ++n, ++bb3) {
	    if (!b->exin->isCanon(l, n)) continue;
	    VTYPE	ip = pwd->IntPen->Penalty(n - l);
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
	int	l = nearest5ss(b->left, bab);
	if (l < 0) {		// no splice donor site nearby
	    b->right = min(b->left + a->right + wlmt, b->len);
	    return (NEVSEL);
	}
	int	d = l - b->left;
	b->left = l;
	a->left += d;
	VTYPE	maxscr = NEVSEL;
	int	maxf = -1;
	BoyerMoore	bm(b, a, 1);
	while (!bm.finished()) {
	    int	f = bm.nexthit();
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
	    if ((maxf = last_exon_wmm(maxscr)) < 0) return (NEVSEL);
	} else {
const	    CHAR*	ts = a->at(a->right);
	    for (const CHAR* as = a->at(a->left); as < ts; ++as)
		maxscr += pwd->simmtx->mtx[*as][*as];
	}
	SKL	skl = {a->left, b->left};
	mfd->write(&skl);
	skl.n = maxf;
	mfd->write(&skl);
	return (maxscr);
}

VTYPE Aln2s1::interpolateS(INT level, int agap, int bgap, int cmode, 
		JUXT* wjxt, const BOUND& bab)
{
	int	ovr = min(agap, bgap);
	int	dgap = bgap - agap;
	bool	cont = agap <= 0;			// connect adjacent HSPs
	bool	no_rec = ovr < wlmt;			// no recurrsion
	VTYPE	iscore = NEVSEL;
	VTYPE	scr = 0;
	Mfile*	save_mfd = 0;
	if (dgap == 0) {
	    if (agap == 0) {
		mfd->write((UPTR) wjxt);
		return (0);
	    }
	    scr += diagonalS_ng();			// no gap
	    iscore = 0;
	} else if (cmode == 1 && no_rec && wjxt) {
	    if (cont) {		// 5' end
		mfd->write((UPTR) wjxt);
		iscore = 0;
	    } else if (agap < IntronPrm.elmt) {
		SKL	skl = {wjxt->jx - agap, wjxt->jy - agap};
		mfd->write(&skl);
		mfd->write((UPTR) wjxt);
		iscore = agap * getsmn(4);
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
		    mfd->write(&skl);
		}
	    } else
		iscore = last_exon(bab);
	} else if (cmode == 3 && cont && dgap >= IntronPrm.minl && 
	    indelfreespjS(agap, iscore)) {		// indel-free spj
	    scr += iscore;
	    iscore = 0;
	} else if (cmode == 3 && no_rec && dgap >= IntronPrm.minl) {
	    if (!cont) iscore = micro_exon(bab);
	    if (iscore == NEVSEL && agap < IntronPrm.elmt)
		iscore = shortcutS_ng(ovr, bab);	// shortcut
	} else if (cont && dgap < IntronPrm.minl) {
	    scr += backforth(-ovr, bab);		// ordinary gap
	    iscore = 0;
	} else if (dgap < IntronPrm.minl) {
	    scr -= creepback(ovr, slmt, bab);
	    scr -= creepfwrd(ovr, slmt, bab);
	    WINDOW	wdw;
	    stripe(seqs, &wdw, min(alprm.sh, abs(dgap) + 3));
	    scr += trcbkalignS_ng(wdw, false);		// ordinary alignment
	    iscore = 0;
	} else if (level + 1 < algmode.qck) {
	    save_mfd = new Mfile(*mfd);			// recursive search
	    island = false;
	    iscore = seededS_ng(level + 1, cmode, bab);
	    if (iscore >= 0 || !island) {		// reached bottom
		scr += iscore;
		iscore = 0;
	    }
	}
	if (iscore < 0) {
	    if (save_mfd) swap(mfd, save_mfd);
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
		stripe(seqs, &wdw, alprm.sh);
		kscore = trcbkalignS_ng(wdw, b->inex.intr);	// DP
	    } else if (cmode == 1) {
		if (wjxt) {
		    int	bl = wjxt->jy + term;		// near 5' end
		    if (bl > b->left) b->left = bl;
		}
		kscore = openendS_ng(cmode);		// drop off
	    } else if (cmode == 2) {
		if (wjxt) {
		    int	br = wjxt->jy - term;		// near 3' end
		    if (br > b->left && br < b->right) b->right = br;
		}
		kscore = openendS_ng(cmode);		// drop off
	    } else {		// cmode == 3, space-saving DP
		WINDOW	wdw;
		stripe(seqs, &wdw, alprm.sh);
		kscore = lspS_ng(wdw);
	    }
	    if (iscore > kscore) {
		if (save_mfd) swap(mfd, save_mfd);
		scr += iscore;
	    } else if (kscore > NEVSEL)
		scr += kscore;
	    else	scr = NEVSEL;
	}
	delete save_mfd;
	return (scr);
}

/* eimode 1: 5'end, 2: 3' end, 3: internal */
VTYPE Aln2s1::seededS_ng(INT level, int eimode, const BOUND& lub)
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
	} else {
	    wl = new Wilip(seqs, pwd, level);
	    wlu = wl->begin();
	    if (wlu) {
		num = wlu->num;
		jxt = wlu->jxt;
		island = true;
	    } else	num = 0;
	}
	save_range(seqs, rng, 2);
	WLPRM*	wlprm = setwlprm(level);
	wlmt = wlprm->width;
	BOUND	bab = lub;
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
		if (wjxt[1].jx < bab.ua) bab.ua = wjxt[1].jx;
		bab.ua -= wlmt;
		int	lx = wjxt->jx + wjxt->jlen / 2;
		if (bab.ua < lx) bab.ua = lx;
		bab.ub = wjxt->jy + bab.ua - wjxt->jx;
		if (cmode == 2) cmode = 3;
		VTYPE	iscore = interpolateS(level, agap, bgap, cmode, wjxt, bab);
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

	VTYPE	iscore = interpolateS(level, agap, bgap, cmode, wjxt, bab);
	if (iscore != NEVSEL) scr += iscore;
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
	if (skl->n == 0) {
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
	return alnv.scorealoneS_ng(wdw);
}

static int infer_orientation(Seq** sqs, const PwdB* pwd)
{
	Seq*&	a = sqs[0];
	VTYPE	scr1 = HomScoreS_ng((const Seq**) sqs, pwd);
	a->comrev();
	antiseq(sqs + 1);
	VTYPE	scr2 = HomScoreS_ng((const Seq**) sqs, pwd);
	if (scr2 > scr1) return (1);
	a->comrev();
	antiseq(sqs + 1);
	return (0);
}

Colonies* swg1stS_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr)
{
	Aln2s1	alnv(seqs, pwd);
	WINDOW	wdw;
	stripe(seqs, &wdw, alprm.sh);
	return (alnv.fwdswgS_ng(wdw, scr));
}

SKL* swg2ndS_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr, COLONY* clny)
{
	if (clny->val <= 0) return (0);
const	Seq*	a = seqs[0];
const	Seq*	b = seqs[1];

	a->left = clny->mlb;
	b->left = clny->nlb;
	b->right = clny->nrb;
	a->right = clny->mrb;
	WINDOW	wdw = {clny->upr, clny->lwr, clny->upr - clny->lwr + 3};
	Aln2s1 alnv(seqs, pwd);
	return (alnv.globalS_ng(wdw, scr));
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

