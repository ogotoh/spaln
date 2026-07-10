/*****************************************************************************
*
*	abstract structures of alingmnet
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-2023)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
*
*****************************************************************************/

#include "aln.h"

#define CDEBUG	1
#define	UNITE_INDEL_FS	0

typedef	int	DIM2[2];

static	int	scmpf(const SKL* a, const SKL* b);

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

/*	A folded form of Gaps record differs from an unfolded form
	in that the 'gps' field indicates the sequence position
	after removal of all preceding gaps.
*/

void toimage(Gaps* gaps[], int numseq)
{
	for (int j = 0; j < numseq; j++) {
	    int 	gap = 0;
	    GAPS*	gp = gaps[j]->begin();
	    GAPS*	gt = gaps[j]->end();
	    for ( ; gp < gt; ++gp) {
		gp->gps += gap;
		gap += gp->gln;
	    }
	}
}

void unfoldgap(Gaps* gg, int step)
{
	GAPS*	gp = gg->begin();
	GAPS*	gt = gg->end();
	int	base = gp->gps;
	int 	glen = gp->gln + base;

	while (++gp < gt) {
	    gp->gps = step * (gp->gps - base) + glen;
	    glen += gp->gln;
	}
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
const	int	num = skl->n;
const	SKL*	prv = ++skl;
const	SKL*	trm = prv + num;
	if (sqs) {
const	    Seq*&	a = sqs[0];
const	    Seq*&	b = sqs[1];
const	    int	m0 = sqs? sqs[0]->left: 0;
const	    int	n0 = sqs? sqs[1]->left: 0;
	    if ((!b->inex.exgl && skl->m != m0) || 
		(!a->inex.exgl && skl->n != n0)) return (true);
	}
	for (++skl; skl < trm; prv = skl++) {
const	    int	dm = skl->m - prv->m;
const	    int	dn = skl->n - prv->n;
	    if (dm != dn && dm && dn) return (true);
	}
	return (false);
}

SKL* stdskl(SKL** pskl)
{
const	int	num = (*pskl)->n;
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
const	int	partner = 1 - given;
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

SKL* gap2skl(const Gaps* gga, const Gaps* ggb, const Seq** sqs)
{
const 	GAPS*	gaps[2] = {gga->begin(), ggb->begin()};
const 	GAPS*&	ga = gaps[0];
const 	GAPS*&	gb = gaps[1];
const 	GAPS*	gat = gga->last();
const 	GAPS*	gbt = ggb->last();
	int	maxskl = 2 * (gga->size() + ggb->size());
	int	mn[2] = {ga->gps - gb->gps, gb->gps - ga->gps};
	int	ndel[2] = {0, 0};
	int	node[2] = {ga->gps, gb->gps};
	bool	parity[2] = {false, false};
	SKL	pinc = {1, -1};
	SKL*	skl = new SKL[maxskl + 2];
	SKL*	wsk = skl + 1;

	wsk->m = node[0];
	wsk->n = node[1];
	++wsk;

	do {
	    int	i = (node[0] > node[1] || ga >= gat)? 1: 0;
	    int	j = 1 - i;
	    if (parity[i]) ndel[i] += gaps[i]->gln;
	    wsk->m = node[i] - ndel[i];
	    wsk->n = std::min(node[i] - mn[i], gaps[j]->gps) - ndel[j];
	    if (i) swap(wsk->m, wsk->n);
	    SKL	cinc = {wsk->m - (wsk-1)->m, wsk->n - (wsk-1)->n};
	    if (cinc.m < 0 || cinc.n < 0) {
		delete[] skl;
		return (0);
	    }
	    if (pinc.m * cinc.n != pinc.n * cinc.m) {
		++wsk;
		pinc = cinc;
	    } else
		*(wsk - 1) = *wsk;
	    if (parity[i])
		node[i] = (++gaps[i])->gps;
	    else
		node[i] += gaps[i]->gln;
	    parity[i] = !parity[i];
	} while (ga < gat || gb < gbt);
	wsk->m = node[0] - ndel[0];
	wsk->n = node[1] - ndel[1];
	SKL	cinc = {wsk->m - (wsk-1)->m, wsk->n - (wsk-1)->n};
	if (cinc.m < 0 || cinc.n < 0) {
	    delete[] skl;
	    return (0);
	} else if (cinc.m || cinc.n) 
	    ++wsk;
	wsk->m = wsk->n = EOS;
const	int	i = wsk - skl;
	skl->n = i - 1;
	skl->m = 1;
	if (badskl(skl, sqs)) {
	    delete[] skl; skl = 0;
	}
	return (skl);
}

void skl2gaps(Gaps* gaps[], const SKL* skl)
{
	int	num = (skl++)->n;
	delete gaps[0];
	delete gaps[1];
	gaps[0] = new Gaps(num);
	gaps[1] = new Gaps(num);
	GAPS*	wga = gaps[0]->begin();
	GAPS*	wgb = gaps[1]->begin();
	GAPS	gpa = {skl->m, 0};
	GAPS	gpb = {skl->n, 0};
	int	span = 0;

	*wga++ = gpa;
	*wgb++ = gpb;
	while (--num) {
	    ++skl;
const	    int	d[2] = {skl->m - gpa.gps, skl->n - gpb.gps};
const	    int	i = d[0] - d[1];
	    gpa.gps = skl->m;
	    gpb.gps = skl->n;
	    if (i < 0) {
		gpa.gln = -i;
		span -= i;
		*wga++ = gpa;
	    } else if (i > 0 ) {
		gpb.gln = i;
		span += i;
		*wgb++ = gpb;
	    } else {
		gpa.gln = gpb.gln = 0;
		span += d[0];
	    }
	}
	gpa.gln = gpb.gln = 0;
	*wga++ = gpa;
	*wgb++ = gpb;
	gaps[0]->elms = wga - gaps[0]->begin();
	gaps[1]->elms = wgb - gaps[1]->begin();
	gaps[0]->alen = gaps[1]->alen = span;
}

void skl2gaps3(Gaps* gaps[], const SKL* skl, const int& pro)
{
	int	num = (skl++)->n;
	delete gaps[0];
	delete gaps[1];
	gaps[0] = new Gaps(num);
	gaps[1] = new Gaps(num);
	GAPS*	wga = gaps[0]->begin();
	GAPS*	wgb = gaps[1]->begin();
	GAPS	gpa = {skl->m, 0};
	GAPS	gpb = {skl->n, 0};
	int	span = 0;

	*wga++ = gpa;
	*wgb++ = gpb;
	while (--num) {
	    ++skl;
const	    int	d[2] = {skl->m - gpa.gps, skl->n - gpb.gps};
const	    int	i = (pro == 1)? 3 * d[0] - d[1]: d[0] - 3 * d[1];
	    gpa.gps = skl->m;
	    gpb.gps = skl->n;
	    if (i < 0) {
		gpa.gln = -i;
		span -= i;
		*wga++ = gpa;
	    } else if (i > 0 ) {
		gpb.gln = i;
		span += i;
		*wgb++ = gpb;
	    } else {
		span += (pro == 1)? d[1]: d[0];
	    }
	}
	gpa.gln = gpb.gln = 0;
	*wga++ = gpa;
	*wgb++ = gpb;
	gaps[0]->elms = wga - gaps[0]->begin();
	gaps[1]->elms = wgb - gaps[1]->begin();
	gaps[0]->alen = gaps[1]->alen = span;
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

