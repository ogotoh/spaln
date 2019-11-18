/*****************************************************************************
*
*	abstract structures of alingmnet
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

#include "aln.h"

#define CDEBUG	1
#define	UNITE_INDEL_FS	0

typedef	int	DIM2[2];

static char memerror[] = "Memory error at %s, %d > %d";

static	int	scmpf(const SKL* a, const SKL* b);

/*	Each GAPS* variable has a header containing the length of the
	alignment and the number of elements in that array.  These values
	are returned by the macros defined in "gaps.h".
*/

#if CDEBUG

void putskl(const SKL* skl)
{
	int	i = 0;
	int	w = 0;
	int	prv[2];
	int	nn = (skl++)->n;

	prv[0] = skl->m;
	prv[1] = skl->n;
	while (i++ < nn) {
	    printf("%3d %3d: ", skl->m, skl->n);
	    if (!(i % 8)) putchar('\n');
	    prv[0] -= skl->m;
	    prv[1] -= skl->n;
	    if (prv[0] && prv[1] && prv[0] != prv[1]) w = 1;
	    prv[0] = skl->m;
	    prv[1] = skl->n;
	    skl++;
	}
	if (w) putchar('*');
	putchar('\n');
}

#endif


/*	A folded form of GAPS record differs from an unfolded form
	in that the 'gps' field indicates the sequence position
	after removal of all preceding gaps.
*/

void toimage(GAPS* gaps[], int numseq)
{
	for (int j = 0; j < numseq; j++) {
	    int 	gap = 0;
	    GAPS*	gp = gaps[j];
	    for ( ; gaps_intr(gp); gp++) {
		gp->gps += gap;
		gap += gp->gln;
	    }
	    (gp--)->gps += gap;
	    if (gp[1].gps == gp->gps + gp->gln)
		gp->gln = gaps_end;
	}
}

void unfoldgap(GAPS* gp, int step, bool hl)
{
	GAPS*	gg = gp;
	if (hl) ++gp;
	int	bas = gp->gps;
	int 	gap = gp->gln + bas;

	while (gaps_intr(++gp)) {
	    gp->gps = step * (gp->gps - bas) + gap;
	    gap += gp->gln;
	}
	gp->gps = step * (gp->gps - bas) + gap;
	if (hl) {
	    gg->gps = gp->gps;
	    return;
	}
/*
	--gp;
	if (gp[1].gps == gp->gps + gp->gln)
	    gp->gln = gaps_end;
*/
}

void swapskl(SKL* skl)
{
	if (!skl) return;
	int	num = (skl++)->n;
	for ( ; num--; ++skl) swap(skl->m, skl->n);
}

static int scmpf(const SKL* a,const SKL* b)
{
	int	d = a->m - b->m;
	if (d) return d;
	return (a->n - b->n);
}

bool badskl(const SKL* skl, const Seq** sqs)
{
	int	num = skl->n;
const	SKL*	prv = ++skl;
const	SKL*	trm = prv + num;
	int	m0 = sqs? sqs[0]->left: 0;
	int	n0 = sqs? sqs[1]->left: 0;
	if (skl->m != m0 || skl->n != n0) return (true);
	for (++skl; skl < trm; prv = skl++) {
	    int	dm = skl->m - prv->m;
	    int	dn = skl->n - prv->n;
	    if (dm != dn && dm && dn) return (true);
	}
	return (false);
}

SKL* stdskl(SKL** pskl)
{
	int	num = (*pskl)->n;
	SKL*	org = *pskl + 1;
	int	pr = 2;

	if (num < 2) return (*pskl);
	SKL*	std = new SKL[2 * num + 1];
	SKL*	wrk = std + 1;
	qsort((UPTR) org, (INT) num, sizeof(SKL), (CMPF) scmpf);
	SKL*	prv = org;
	for (int i = 1; i < num; i++) {
	    org++;
	    int	dm = org->m - prv->m;
	    int	dn = org->n - prv->n;
	    if (!dm && !dn) continue;	/* no increment */
	    if (dm < 0 || dn < 0) continue;	/* skip inconsistent */
	    int	dd = min(dm, dn);
	    int	df = dn - dm;
	    if (df) df = (df > 0)? 1: -1;
	    if (dd && df) {		/* interpolate */
	        if (pr) *wrk++ = *prv;
		wrk->m = prv->m + dd;
		wrk->n = prv->n + dd;
		wrk++;
	    } else if (df != pr || !dm)
		*wrk++ = *prv;
	    pr = df;
	    prv = org;
	}
	*wrk++ = *prv;
	wrk->m = wrk->n = EOS;
	std->n = wrk - std - 1;
	std->m = (*pskl)->m;
	delete[] *pskl;
	return (*pskl = std);
}

SKL* stdskl3(SKL** pskl)
{
	int	num = (*pskl)->n;
	SKL*	org = *pskl + 1;
	int	pr = -2;

	if (num < 2) return (*pskl);
	qsort((UPTR) org, (INT) num, sizeof(SKL), (CMPF) scmpf);
	Mfile	mfd(sizeof(SKL));
	mfd.write(*pskl);			//	make 1st line
	SKL*	prv = org;
	SKL	sklbuf = *prv;
	for (int i = 1; i < num; ++i) {
	    ++org;
	    int	dm = (org->m - prv->m) * 3;
	    int	dn = org->n - prv->n;
	    if (!dm && !dn) continue;		// no increment
	    if (dn < 0) continue;		// skip inconsistency
	    int	dd = min(dm, dn);
	    int	df = dn - dm;
	    int	dr = df? ((df > 0)? 1: -1): 0;	// direction
#if UNITE_INDEL_FS
	    if (dd && (df > 0 || !(df % 3))) {	// unite deletion + FS
#else
	    if (dd && df) {			// interpolate
#endif
		if (pr) mfd.write(prv);
		sklbuf.n = prv->n + dd;
		if (df < 0 && df % 3) dd += 2;	// deletion frame shift
		sklbuf.m = prv->m + dd / 3;
		mfd.write(&sklbuf);
#if !UNITE_INDEL_FS
		if (df > 0 && df % 3) {
		    sklbuf.n += df % 3;
		    mfd.write(&sklbuf);
		}
#endif
	    } else if (dr != pr || !dm)	{	// change direction
		mfd.write(prv);			// or retain spj
	    }
	    pr = dr;
	    prv = org;
	}
	mfd.write(prv);
	num = mfd.size();
	sklbuf.m = sklbuf.n = EOS;
	mfd.write(&sklbuf);
	SKL*	std = (SKL*) mfd.flush();
	std->n = num - 1;
	std->m = (*pskl)->m;
	delete[] *pskl;
	return (*pskl = std);
}

int sklpartner(const SKL* skl, int m, int given)
{
	int	partner = 1 - given;
const 	DIM2*	coo = (const DIM2*) skl;
const 	DIM2*	boo = coo + skl->n - 2;

	if (given && partner) return (ERROR);
	int	step[2] = {skl->m & STEP3m? 3: 1, skl->m & STEP3n? 3: 1};
	while (++coo < boo) {
	    int	offset = m - coo[0][given];
	    if (offset >= 0 && m < coo[1][given]) {
		offset *= step[partner];
		offset /= step[given];
		return (coo[0][partner] + offset);
	    }
	}
	return (coo[0][partner]);
}

//	delete terminal gaps
//	assume already standarized

SKL* trimskl(const Seq* seqs[], SKL* skl)
{
	int&	pn = skl->n;
	SKL*	wsk = skl + 1;
	SKL*	boo = skl + pn;
	int	i = wsk[1].m - wsk->m;
	int	j = wsk[1].n - wsk->n;

	if ((seqs[0]->inex.exgl && !i) || (seqs[1]->inex.exgl && !j)) {
	    for ( ; wsk <= boo; ++wsk) *wsk = wsk[1];
	    --pn; --boo;
	}
	i = boo->m - boo[-1].m;
	j = boo->n - boo[-1].n;
	if ((seqs[0]->inex.exgr && !i) || (seqs[1]->inex.exgr && !j)) {
	    --pn;
	    *boo = boo[1];	// sentinel
	}
	return (skl);
}

SKL* gap2skl(const GAPS* gga, const GAPS* ggb)
{
const 	GAPS*	ga = gga;
const 	GAPS*	gb = ggb;
	int	maxskl = 2 * (gaps_size(ga) + gaps_size(gb));
const 	GAPS*	gaps[2] = {++ga, ++gb};
	int	mn[2] = {ga->gps - gb->gps, gb->gps - ga->gps};
	int	ndel[2] = {0, 0};
	int	node[2] = {gaps[0]->gps, gaps[1]->gps};
	bool	parity[2] = {false, false};
	SKL	pinc = {1, -1};
	SKL*	skl = new SKL[maxskl + 1];
	SKL*	wk = skl + 1;

	wk->m = node[0];
	wk->n = node[1];
	wk++;

	do {
	    int	i = (node[0] > node[1] || !gaps_intr(gaps[0]))? 1: 0;
	    int	j = 1 - i;
	    if (parity[i]) ndel[i] += gaps[i]->gln;
	    wk->m = node[i] - ndel[i];
	    wk->n = min(node[i] - mn[i], gaps[j]->gps) - ndel[j];
	    if (i) swap(wk->m, wk->n);
	    SKL	cinc = {wk->m - (wk-1)->m, wk->n - (wk-1)->n};
	    if (cinc.m < 0 || cinc.n < 0) {
		delete[] skl;
		return (0);
	    }
	    if (pinc.m * cinc.n != pinc.n * cinc.m) {
		wk++;
		pinc = cinc;
	    } else
		*(wk - 1) = *wk;
	    if (parity[i])
		node[i] = (++gaps[i])->gps;
	    else
		node[i] += gaps[i]->gln;
	    parity[i] = !parity[i];
	} while (gaps_intr(gaps[0]) || gaps_intr(gaps[1]));
	wk->m = node[0] - ndel[0];
	wk->n = node[1] - ndel[1];
	SKL	cinc = {wk->m - (wk-1)->m, wk->n - (wk-1)->n};
	if (cinc.m < 0 || cinc.n < 0) {
	    delete[] skl;
	    return (0);
	} else if (cinc.m || cinc.n) ++wk;
	wk->m = wk->n = EOS;
	int	i = wk - skl;
	skl->n = i - 1;
	skl->m = 1;
	if (++i > maxskl) fatal(memerror, "gap2skl()", i, maxskl);
	if (badskl(skl)) {
	    delete[] skl; skl = 0;
	}
	return (skl);
}

void skl2gaps(GAPS* gaps[], const SKL* skl, bool hl)
{
	int	num = (skl++)->n;
	GAPS*	wga = gaps[0] = new GAPS[num + 2];
	GAPS*	wgb = gaps[1] = new GAPS[num + 2];
	GAPS	gpa = {skl->m, 0};
	GAPS	gpb = {skl->n, 0};

	*wga++ = gpa;
	*wgb++ = gpb;
	if (hl) {
	    *wga++ = gpa;
	    *wgb++ = gpb;
	}

	while (--num) {
	    skl++;
	    int	i = skl->m - skl->n - gpa.gps + gpb.gps;
	    gpa.gps = skl->m;
	    gpb.gps = skl->n;
	    if (i < 0) {
		gpa.gln = -i;
		*wga++ = gpa;
	    } else if (i > 0 ) {
		gpb.gln = i;
		*wgb++ = gpb;
	    }
	} 
	gpa.gln = gpb.gln = gaps_end;
	*wga++ = gpa;
	*wgb++ = gpb;
	if (hl) {
	    gaps[0]->gln = wga - gaps[0];
	    gaps[1]->gln = wgb - gaps[1];
	    gaps[0]->gps = skl->m - gaps[0][1].gps;
	    gaps[1]->gps = skl->n - gaps[1][1].gps;
	}
}

void skl2gaps3(GAPS* gaps[], const SKL* skl, int pro)
{
	int	num = (skl++)->n;
	GAPS*	wga = gaps[0] = new GAPS[num + 2];
	GAPS*	wgb = gaps[1] = new GAPS[num + 2];
	GAPS	gpa = {skl->m, 0};
	GAPS	gpb = {skl->n, 0};

	*wga++ = gpa;
	*wgb++ = gpb;
	while (--num) {
	    skl++;
	    int	m = skl->m;
	    int	n = skl->n;
	    int	i;
	    if (pro == 1) {
		i = 3 * (m - gpa.gps) - (n - gpb.gps);
	    } else {
		i = m - gpa.gps - 3 * (n - gpb.gps);
	    }
	    gpa.gps = m;
	    gpb.gps = n;
	    if (i < 0) {
		gpa.gln = -i;
		*wga++ = gpa;
	    } else if (i > 0 ) {
		gpb.gln = i;
		*wgb++ = gpb;
	    }
	}
	gpa.gln = gpb.gln = gaps_end;
	*wga++ = gpa;
	*wgb++ = gpb;
}

bool sameskl(const SKL* a, const SKL* b)
{
	if (!a || !b || a->n != b->n) return (false);
	const	SKL*	tsk = a + a->n + 1;
	while (++a < tsk) {
	    ++b;
	    if (a->m != b->m || a->n != b->n) return (false);
	}
	return (true);
}

