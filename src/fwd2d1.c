/*****************************************************************************
*
*	fwd2d1.c
*
*	Calculate Pairwise alignment score of two protein/nucleotide sequences 
*	with affine gap-penalty function
*	> scan DP matrix in aniti-diagnal order
*	> may be useful for fine-grain parallell processing
*	> score variables are of integer type
*	> must run prePwd(molc) beforehand
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

#include "seq.h"
#include "aln.h"

#define	D_DEBUG		0

struct ValDia {VTYPE val; int r;};

class Fwd2d {
	Seq**	seqs;
	Seq*&	a;
	Seq*&	b;
	WINDOW	wdw;
	VTYPE*	hh;
	VTYPE*	ff;
	VTYPE*	gg;
	VTYPE	lastD();
	CHAR*	as;
	CHAR*	bs;
	VTYPE**	mtx;
	VTYPE	uu;
	VTYPE	vv;
public:
	Fwd2d(Seq** sqs, Simmtx* sm);
	VTYPE	swgforwardD();
	~Fwd2d();
	VTYPE	forwardD();
};

Fwd2d::Fwd2d(Seq** sqs, Simmtx* sm)
	: seqs(sqs), a(seqs[0]), b(seqs[1]), mtx(sm->mtx)
{
	as = a->at(a->left);
	bs = b->at(b->left);
	uu = (VTYPE) (alprm.u * alprm.scale);
	vv = (VTYPE) (alprm.v * alprm.scale);
	stripe(seqs, &wdw, alprm.sh);
	hh = new VTYPE[3 * wdw.width] - wdw.lw + 1;
	ff = hh + wdw.width;
	gg = ff + wdw.width;
	int	r0 = b->left - a->left;
	if (a->inex.exgl) vclear(hh + r0, wdw.up - r0);
	else { 
	    hh[r0] = 0;
	    FTYPE	ltgapf = a->left? 1: alprm.tgapf;
	    VTYPE	gp = (VTYPE) (-vv * ltgapf);
	    VTYPE	ge = (VTYPE) (-uu * ltgapf);
	    for (int r = r0 + 1; r <= wdw.up; ++r)
		hh[r] = gp += ge;
	}
	hh[wdw.up + 1] = NEG_INT;
	if (b->inex.exgl) vclear(hh + wdw.lw, r0 - wdw.lw - 1);
	else {
	    FTYPE	ltgapf = b->left? 1: alprm.tgapf;
	    VTYPE	gp = (VTYPE) (-vv * ltgapf);
	    VTYPE	ge = (VTYPE) (-uu * ltgapf);
	    for (int r = r0 - 1; r >= wdw.lw; --r)
		hh[r] = gp += ge;
	}
	hh[wdw.lw - 1] = NEG_INT;
	vset(ff + wdw.lw - 1, NEVSEL, wdw.width);
	vset(gg + wdw.lw - 1, NEVSEL, wdw.width);
}

Fwd2d::~Fwd2d()
{
	delete[] (hh + wdw.lw - 1);
}

VTYPE Fwd2d::lastD()
{
	VTYPE*	h9 = hh + b->right - a->right;

	FTYPE	rtgapf = b->inex.exgr? 0: alprm.tgapf;
	if (b->right == b->len && rtgapf < 1) {
	    int	dm = 0;
	    int	rw = wdw.up + 1;
	    int	rf = b->right - a->left;
	    if (rf < rw) rw = rf;
	    VTYPE*	h = hh + rw;
	    while (--h >= h9) {
		++dm;
		VTYPE	gpn = dm == 1? vv + uu: uu;
		VTYPE*	g = h + 1;
		*g += (VTYPE) (gpn * rtgapf);
		if (*h < *g) *h = *g;
		else	dm = 0;
	    }
	}
	rtgapf = a->inex.exgr? 0: alprm.tgapf;
	if (a->right == a->len && rtgapf < 1) {
	    int	dn = 0;
	    int rw = wdw.lw;
	    int rf = b->left - a->right + 1;
	    if (rf > rw) rw = rf;
	    VTYPE*	h = hh + rw;
	    while (++h <= h9) {
		++dn;
		VTYPE   gpn = dn == 1? vv + uu: uu;
		VTYPE*  f = h - 1;
		*f += (VTYPE) (gpn * rtgapf);
		if (*h < *f) *h = *f;
		else	dn = 0;
	    }
	}
	return (*h9);
}

VTYPE Fwd2d::forwardD()
{
	for (int d = a->left + b->left; d < a->right + b->right - 1; ++d) {
	    int	n = max(b->left, max(d - a->right + 1, (d + wdw.lw + 1) / 2));
	    int	m = d - n;
	    int m9 = max(a->left, max(d - b->right + 1, (d - wdw.up + 1) / 2));
	    int	n9 = d - m9 + 1;
	    int	r0 = n - m;
	    int	r9 = n9 - m9;

	    for (int r = r0; r < r9; r += 2, --m, ++n) {
		ff[r] = max(hh[r-1] - vv, ff[r-1]) - uu;
	        gg[r] = max(hh[r+1] - vv, gg[r+1]) - uu;
		hh[r] += (VTYPE) mtx[as[m]][bs[n]];
		hh[r] = max(max(hh[r], ff[r]), gg[r]);
#if D_DEBUG
		if (algmode.nsa & 8)
		    printf("%5d %5d %5d %5d: %7.1f %7.1f %7.1f\n",
			m, n, d, r, (float) hh[r], 
			(float) ff[r], (float) gg[r]);
#endif	// D_DEBUG
	    }
	}
	return lastD();
}

VTYPE Fwd2d::swgforwardD()
{
	VTYPE	maxh = NEVSEL;
	for (int d = a->left + b->left; d < a->right + b->right - 1; ++d) {
	    int	n = max(b->left, max(d - a->right + 1, (d + wdw.lw + 1) / 2));
	    int	m = d - n;
	    int m9 = max(a->left, max(d - b->right + 1, (d - wdw.up + 1) / 2));
	    int	n9 = d - m9 + 1;
	    int	r0 = n - m;
	    int	r9 = n9 - m9;
	    VTYPE	vzero = 0;

	    for (int r = r0; r < r9; r += 2, --m, ++n) {
		ff[r] = max(hh[r-1] - vv, ff[r-1]) - uu;
	        gg[r] = max(hh[r+1] - vv, gg[r+1]) - uu;
		hh[r] += (VTYPE) mtx[as[m]][bs[n]];
		hh[r] = max(max(max(hh[r], ff[r]), gg[r]), vzero);
		maxh = max(maxh, hh[r]);
#if D_DEBUG
		if (algmode.nsa & 8)
		    printf("%5d %5d %5d %5d: %7.1f %7.1f %7.1f\n",
			m, n, d, r, (float) hh[r], 
			(float) ff[r], (float) gg[r]);
#endif	// D_DEBUG
	    }
	}
	return maxh;
}

class Fwd2d_vd {
	Seq**	seqs;
	Seq*&	a;
	Seq*&	b;
	WINDOW	wdw;
	ValDia*	hh;
	ValDia*	ff;
	ValDia*	gg;
	VTYPE	lastD(int* ends);
	CHAR*	as;
	CHAR*	bs;
	VTYPE**	mtx;
	VTYPE	uu;
	VTYPE	vv;
public:
	Fwd2d_vd(Seq** sqs, Simmtx* sm);
	~Fwd2d_vd();
	VTYPE	forwardD(int* ends);
};

Fwd2d_vd::Fwd2d_vd(Seq** sqs, Simmtx* sm)
	: seqs(sqs), a(seqs[0]), b(seqs[1]), mtx(sm->mtx)
{
	as = a->at(a->left);
	bs = b->at(b->left);
	uu = (VTYPE) (alprm.u * alprm.scale);
	vv = (VTYPE) (alprm.v * alprm.scale);
	stripe(seqs, &wdw, alprm.sh);
	hh = new ValDia[wdw.width] - wdw.lw + 1;
	ff = new ValDia[wdw.width] - wdw.lw + 1;
	gg = new ValDia[wdw.width] - wdw.lw + 1;
	int	r0 = b->left - a->left;
	ValDia	black = {NEG_INT, 0};
	vclear(hh + r0);
	FTYPE	ltgapf = a->left? 1: alprm.tgapf;
	VTYPE	gp = (VTYPE) (-vv * ltgapf);
	VTYPE	ge = (VTYPE) (-uu * ltgapf);
	for (int r = r0 + 1; r <= wdw.up; ++r) {
	    hh[r].val = gp += ge;
	    hh[r].r = r;
	}
	ltgapf = b->left? 1: (b->inex.exgl? 0: alprm.tgapf);
	gp = (VTYPE) (-vv * ltgapf);
	ge = (VTYPE) (-uu * ltgapf);
	for (int r = r0 - 1; r >= wdw.lw; --r) {
	    hh[r].val = gp += ge;
	    hh[r].r = r;
	}
	hh[wdw.lw - 1] = hh[wdw.up + 1] = black;
	vset(ff + wdw.lw - 1, black, wdw.width);
	vset(gg + wdw.lw - 1, black, wdw.width);
}

Fwd2d_vd::~Fwd2d_vd()
{
	delete[] (hh + wdw.lw - 1);
	delete[] (ff + wdw.lw - 1);
	delete[] (gg + wdw.lw - 1);
}

VTYPE Fwd2d_vd::lastD(int* ends)
{
	ValDia*	h9 = hh + b->right - a->right;

	FTYPE   rtgapf = b->inex.exgr? 0: alprm.tgapf;
	int	dm = 0;
	if (b->right == b->len && rtgapf < 1) {
	    int	rw = wdw.up + 1;
	    int	rf = b->right - a->left;
	    if (rf < rw) rw = rf;
	    ValDia*	h = hh + rw;
	    while (--h >= h9) {
		++dm;
		VTYPE	gpn = dm == 1? vv + uu: uu;
		ValDia*	g = h + 1;
		g->val += (VTYPE) (gpn * rtgapf);
		if (h->val < g->val) *h = *g;
		else	dm = 0;
	    }
	}
	rtgapf = a->inex.exgr? 0: alprm.tgapf;
	int	dn = 0;
	if (a->right == a->len && rtgapf < 1) {
	    int rw = wdw.lw;
	    int rf = b->left - a->right + 1;
	    if (rf > rw) rw = rf;
	    ValDia*	h = hh + rw;
	    while (++h <= h9) {
		++dn;
		VTYPE   gpn = dn == 1? vv + uu: uu;
		ValDia*  f = h - 1;
		f->val += (VTYPE) (gpn * rtgapf);
		if (h->val < f->val) *h = *f;
		else	dn = 0;
	    }
	}
	int	r0 = b->left - a->left;
	ends[0] = h9->r - r0;	// +? n0 - b->left: a->left - m0
	ends[1] = dn? dn: -dm;	// +? a->right - m9: n9 - b->right
	return (h9->val);
}

VTYPE Fwd2d_vd::forwardD(int* ends)
{
	for (int d = a->left + b->left; d < a->right + b->right - 1; ++d) {
	    int	n = max(b->left, max(d - a->right + 1, (d + wdw.lw + 1) / 2));
	    int	m = d - n;
	    int m9 = max(a->left, max(d - b->right + 1, (d - wdw.up + 1) / 2));
	    int	n9 = d - m9 + 1;
	    int	r0 = n - m;
	    int	r9 = n9 - m9;

	    for (int r = r0; r < r9; r += 2, --m, ++n) {
		ValDia	ng = {hh[r-1].val - vv - uu, hh[r-1].r};
		ValDia	eg = {ff[r-1].val - uu, ff[r-1].r};
		ff[r] = ng.val > eg.val? ng: eg;
		ng = hh[r+1]; ng.val -= (uu + vv);
		eg = gg[r+1]; eg.val -= uu;
	        gg[r] = ng.val > eg.val? ng: eg;
		hh[r].val += (VTYPE) mtx[as[m]][bs[n]];
		ng = ff[r].val > gg[r].val? ff[r]: gg[r];
		hh[r] = hh[r].val > ng.val? hh[r]: ng;
#if D_DEBUG
		if (algmode.nsa & 8)
		    printf("%5d %5d %5d %5d %5d: %7.1f %7.1f %7.1f\n",
			m, n, d, r, hh[r].r,(float) hh[r].val, 
			(float) ff[r].val, (float) gg[r].val);
#endif	// D_DEBUG
	    }
	}
	return lastD(ends);
}

VTYPE alnScoreD(Seq* seqs[], Simmtx* sm, int* ends)
{
	if (!sm) sm = getSimmtx(0);
	if (algmode.lcl & 16) {		// SWG local alignment
	    Fwd2d	fwd2d(seqs, sm);
	    return (fwd2d.swgforwardD());
	} else if (ends) {		// semi-global with r and l ends
	    Fwd2d_vd	fwd2d(seqs, sm);
	    return (fwd2d.forwardD(ends));
	} else {			// (semi-) global
	    Fwd2d	fwd2d(seqs, sm);
	    return (fwd2d.forwardD());
	}
}

