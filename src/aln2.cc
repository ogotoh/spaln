/*****************************************************************************
*
*	Alignment of two protein/nucleotide sequences 
*	Major subroutines for obtaining a global alignment
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
#include "utilseq.h"
#include <math.h>

static	int	estimlen(int na, int nb, const SKL* skl);
static	void	synthi(CHAR* cs, const CHAR* sss[], int an, int bn, int mi, int ni);

SKL* nogap_skl(const Seq* a, const Seq* b)
{
	SKL*	skl = new SKL[4];
	if (!b) b = a;
	skl->m = 0; skl->n = 2;
	skl[1].m = a->left; skl[1].n = b->left;
	skl[2].m = a->right; skl[2].n = b->right;
	skl[3].m = skl[3].n = EOS;
	return (skl);
}

VTYPE PwdB::GapPenalty3(int i, VTYPE bgop) const {
	if (i == 0) return 0;
	int	d = i / 3;
	VTYPE	x = 0;
	switch (i % 3) {
	    case 1: x = GapE1; break;
	    case 2: x = GapE2; break;
	    default: break;
	}
	return x + ((i > codonk1)? LongGOP * bgop / BasicGOP + d * LongGEP:
	    bgop + d * BasicGEP);
}

VTYPE selfAlnScr(const Seq* sd, const Simmtx* sm)
{
	if (sd->many > 1) prompt("Incorrect self alignment score value !\n");
	VTYPE	scr = 0;
const 	CHAR*	ss = sd->at(sd->left);
const 	CHAR*	tt = sd->at(sd->right);

	for ( ; ss < tt; ++ss)
	    scr += sm->mtx[*ss][*ss];
	return (scr / sd->many);
}

int prePwd(int molc, bool use_mdm)
{
	int	dvsp;
	switch (molc) {
	  case PROTEIN:	dvsp = 3; break;
	  case TRON:	dvsp = 4; break;
	  default:	dvsp = 1; break;
	}
	if (OutPrm.SkipLongGap == 3) OutPrm.SkipLongGap = 0;
	setSimmtxes((ComPmt) dvsp, use_mdm);		// default
	return (dvsp);
}

int prePwd(const Seq* sd, bool use_mdm)
{
	return prePwd(sd->inex.molc, use_mdm);
}

int prePwd(const Seq** seqs, bool use_mdm)
{
	int	dvsp = seqs[0]->isprotein() + 2 * seqs[1]->isprotein();
	if (dvsp == 2) dvsp = 1;		// swap later
	if (!dvsp && seqs[0]->istron() && seqs[1]->istron()) dvsp = 4;
	if (algmode.crs == 3) algmode.crs = 1;	// cross-species
	if (OutPrm.SkipLongGap == 3) OutPrm.SkipLongGap = algmode.lsg? 1: 0;
	setSimmtxes((ComPmt) dvsp, use_mdm);		// default
	if (alprm2.sss < 0.) alprm2.sss = defSss[algmode.crs];
	if (alprm2.z < 0) alprm2.z = dvsp? def_alprm2z: 0;
	return (dvsp);
}

PwdB::PwdB(const Seq** seqs, const ALPRM* alp) :
	simmtx(getSimmtx((alp = alp? alp: &alprm)->mtx_no)), 
	DvsP(seqs[0]->isprotein() + 2 * seqs[1]->isprotein()),
	Noll(max(2, min(NOL, alp->ls))),
	Nrow(Noll), Nrwb(2 * Noll - 1),
	Vab(axbscale(seqs)),
	Vthr((VTYPE) (alp->thr * Vab)),
	BasicGOP((VTYPE) (-alp->v * Vab)),
	BasicGEP((VTYPE) (-alp->u * Vab)),
	LongGEP((VTYPE) (-alp->u1 * Vab)),
	diffu(LongGEP - BasicGEP),
	LongGOP(BasicGOP - diffu * alp->k1),
	MaxGapL(0), ExtraGOP(0), GapE1(0), GapE2(0),
	GapW1(0), GapW2(0), GapW3(0), GapW3L(0),
	pmt(0), codepot(0), exinpot(0), eijpat(0), IntPen(0)
{
	GOP[0] = 0;
	GOP[1] = BasicGOP;
	GOP[2] = LongGOP;
	int	step = (DvsP == 3)? 1: 3;
	codonk1 = alp->ls == 3? step * alp->k1: LARGEN;
	if (!seqs[0]->inex.intr && !seqs[1]->inex.intr) return;	// without splice
	if (DvsP || alprm2.sss > 0.) eijpat = new EijPat(DvsP);	// boundary signal
	if (DvsP || algmode.mns == 0) codepot = new CodePot();	// coding potential
	if (!DvsP) {	// C vs G
	    int	zZ = (alprm2.z > 0) + 2 * (alprm2.Z > 0);
	    if (zZ) exinpot = new ExinPot(zZ);
	} else {	// A vs G
	    ++Nrow; ++Nrwb;
	    MaxGapL = 100;
	    ExtraGOP = (VTYPE) (-alprm2.x * Vab);
	    GapE1 = BasicGEP + ExtraGOP;
	    GapE2 = GapE1 + BasicGEP;
	    GapW1 = GapE1 + BasicGOP;
	    GapW2 = GapE2 + BasicGOP;
	    GapW3 = BasicGOP + BasicGEP;
	    GapW3L = LongGOP + LongGEP;
	    pmt = new Premat(seqs);
	    if (alprm2.Z > 0) exinpot = new ExinPot(2);
	}
	IntPen = new IntronPenalty(Vab, DvsP, eijpat, exinpot);
}

PwdB::~PwdB()
{
	delete pmt;
	delete IntPen;
	delete codepot;
	delete exinpot;
	delete eijpat;
}

void putvar(VTYPE x)
{
	float	fx = (float) x;
	if (fx > NEVSEL/2) printf(" %6.1f", fx);
	else	fputs("  *****", stdout);
}

void stripe(const Seq* seqs[], WINDOW* wdw, int sh)
{
const 	Seq*	a = seqs[0];
const 	Seq*	b = seqs[1];

	if (sh < 0) {
	    int	shorter = min(a->right - a->left, b->right - b->left);
	    sh = -sh * shorter / 100;
	}
	wdw->up = b->right - a->right;
	wdw->lw = b->left - a->left;
	if (wdw->up < wdw->lw) swap(wdw->up, wdw->lw);
	wdw->up += sh;
	wdw->lw -= sh;
	int	p;
	if ((p = b->right - a->left) < wdw->up) wdw->up = p;
	if ((p = b->left - a->right) > wdw->lw) wdw->lw = p;
	wdw->width = wdw->up - wdw->lw + 3;
}

void stripe31(const Seq* seqs[], WINDOW* wdw, int shld)
{
const 	Seq*	a = seqs[0];
const 	Seq*	b = seqs[1];

	if (shld < 0) {
	    int	shorter = min(a->right - a->left, b->right - b->left);
	    shld = -shld * shorter / 100;
	}
	shld *= 3;
	wdw->up = b->right - 3 * a->right;
	wdw->lw = b->left - 3 * a->left;
	if (wdw->up < wdw->lw) swap(wdw->up, wdw->lw);
	wdw->up += shld;
	wdw->lw -= shld;
	int	p;
	if ((p = b->right - 3 * a->left) < wdw->up) wdw->up = p;
	if ((p = b->left - 3 * a->right) > wdw->lw) wdw->lw = p;
	wdw->width = wdw->up - wdw->lw + 7;
}

static int estimlen(int na, int nb, const SKL* skl)
{
	int	num = (skl++)->n;
	int	m = skl->m;
	int	n = skl->n;

	while (--num) {
		++skl;
		int	mi = skl->m - m;
		int	ni = skl->n - n;
		int	i = mi - ni;
		if (i > 0) nb += i;
		if (i < 0) na -= i; 
		m = skl->m;
		n = skl->n;
	}
	return (max(na, nb));
}

static void synthi(CHAR* cs, const CHAR* sss[], int an, int bn, int mi, int ni)
{
	if (mi == ni)
	    for (int i = 0; i < mi; i++) {
		for (int j = 0; j < an; j++)
		    *cs++ = *sss[0]++;
		for (int j = 0; j < bn; j++)
		    *cs++ = *sss[1]++;
	    }
	else if (mi)
	    for (int i = 0; i < mi; i++) {
		for (int j = 0; j < an; j++)
		    *cs++ = *sss[0]++;
		for (int j = 0; j < bn; j++)
		    *cs++ = gap_code;
	    }
	else
	    for (int i = 0; i < ni; i++) {
		for (int j = 0; j < an; j++)
		    *cs++ = gap_code;
		for (int j = 0; j < bn; j++)
		    *cs++ = *sss[1]++;
	    }
}

Seq* synthseq(Seq* c, const Seq* a, const Seq* b, const SKL* skl)
{
	int	mi = a->right - a->left + 2;
	int	ni = b->right - b->left + 2;
	int	estlen = estimlen(mi, ni, skl);
	int	num = (skl++)->n;
	int	m = skl->m;
	int	n = 0;

	if (c)	c->refresh(a->many + b->many, estlen);
	else	c = new Seq(a->many + b->many, estlen);
	a->copyattr(c);
	for ( ; n < a->many; n++) {
	    c->nbr[n] = a->nbr[n];
	    c->sname->push((*a->sname)[n]);
	}
	for (int i = 0; i < b->many; i++, n++) {
	    c->nbr[n] = b->nbr[i];
	    c->sname->push((*b->sname)[i]);
	}
	CHAR*	cs = c->at(0);
const	CHAR*	sss[2] = {a->at(a->left), b->at(b->left)};
	c->inex.dels = a->inex.dels | b->inex.dels;
	c->inex.ambs = a->inex.ambs | b->inex.ambs;

	n = skl->n;
	while (--num) {
	    skl++;
	    mi = skl->m - m;
	    ni = skl->n - n;
	    int	i = mi - ni;
	    if (i) c->inex.dels = 1;
	    if (!i || !mi || !ni) {
		synthi(cs, sss, a->many, b->many, mi, ni);
	    } else if (i > 0) {
		synthi(cs, sss, a->many, b->many, ni, ni);
		synthi(cs, sss, a->many, b->many, i, 0);
	    } else {
		synthi(cs, sss, a->many, b->many, mi, mi);
		synthi(cs, sss, a->many, b->many, 0, -i);
	    }
	    m = skl->m;
	    n = skl->n;
	}
	return (c->postseq(cs));
}


FTYPE alnscore2dist(Seq* sqs[], const PwdB* pwd, int* end, FTYPE denome)
{
	Seq*&   a = sqs[0];
	Seq*&   b = sqs[1];
	FTYPE   scr;
	int	dlen;
	if (algmode.lcl) {
	    int*	ends = end? end: new int[2];
	    a->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	    b->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	    scr = alnScoreD((const Seq**) sqs, pwd->simmtx, ends);
	    int	al = a->left;
	    int	bl = b->left;
	    int ar = a->right;
	    int br = b->right;
	    if (ends[0] > 0) {
		bl += ends[0];
		swap(bl, b->left);
	    } else if (ends[0] < 0) {
		al -= ends[0];
		swap(al, a->left);
	    }
	    if (ends[1] > 0) {
		br -= ends[1];
		swap(br, b->right);
	    } else if (ends[1] < 0) {
		ar += ends[1];
		swap(ar, a->right);
	    }
	    denome = sqrt(selfAlnScr(a, pwd->simmtx) * selfAlnScr(b, pwd->simmtx));
	    if (ends[0] > 0) swap(bl, b->left);
	    else if (ends[0] < 0) swap(al, a->left);
	    if (ends[1] > 0) swap(br, b->right);
	    else if (ends[1] < 0) swap(ar, a->right);
	    dlen = ar - al - br + bl;
	    if (end != ends) delete[] ends;
	} else {
	    int	al = a->right - a->left;
	    int	bl = b->right - b->left;
	    scr = alnScoreD((const Seq**) sqs, pwd->simmtx, end);
	    if (!denome) denome = selfAlnScr(al < bl? a: b, pwd->simmtx);
	    dlen = al - bl;
	}
	scr += alprm.u * abs(dlen) / 2;
	return (1. - scr / denome);
}

static int bymrb(const COLONY* a, const COLONY* b)
{
	return (a->mrb - b->mrb);
}

static int byclno(const COLONY* a, const COLONY* b)
{
	return (a->clno - b->clno);
}

static int byscr(const COLONY* a, const COLONY* b)
{
	if (a->val == b->val) return (0);
	return (b->val > a->val)? 1: -1;
}

void Colonies::detectoverlap(COLONY* cc)
{
	COLONY* cw = cc;
	int	mlb = cc->mlb + OutPrm.AllowdOverlap;

	while ((--cw)->mrb > mlb) {
	    if (cw->mark > 0) continue;		/* active */
	    if (cc->mrb - cw->mlb > OutPrm.AllowdOverlap &&
		cc->nrb - cw->nlb > OutPrm.AllowdOverlap &&
		cw->nrb - cc->nlb > OutPrm.AllowdOverlap) {
		if (cc->val < cw->val)  cc->mark = -1;
		else    cw->mark = -1;		/* delete */
	    }
	}
}

void Colonies::sortcolonies()
{
	*clny = clny[no_clny];
	qsort((UPTR) clny, (INT) no_clny, sizeof(COLONY), (CMPF) byscr);
}

void Colonies::removeoverlap()
{
	COLONY*	cc;
	COLONY* cix = clny + no_clny;

	clny->mlb = clny->mrb = 0;
	qsort((UPTR) (clny + 1), (INT) no_clny, sizeof(COLONY), (CMPF) bymrb);
	for (cc = cix; cc > clny + 1; --cc)
	    if (cc->mark == 0) detectoverlap(cc);	/* non-active */
	qsort((UPTR) (clny + 1), (INT) no_clny, sizeof(COLONY), (CMPF) byclno);
	for (COLONY* cw = cc = clny + 1; cw <= cix; ++cw) {
	    int	nc = 0;
	    if (cw->mark >= 0) {			/* retain */
		nc = cc - clny;
		*cc++ = *cw;
	    }
	    cw->clno = nc;
	}
	vclear(cc, cix - cc + 1);
	no_clny = --cc - clny;
}

void Colonies::removelowscore()
{
	COLONY*	cc;
	COLONY* cix = clny + no_clny;

	clny->mlb = clny->mrb = 0;
	qsort((UPTR) (clny + 1), (INT) no_clny, sizeof(COLONY), (CMPF) byscr);
	for (cc = cix; cc > clny + OutPrm.NoOut; --cc)
	    if (cc->mark == 0) cc->mark = -1;	/* non-active */
	qsort((UPTR) (clny + 1), (INT) no_clny, sizeof(COLONY), (CMPF) byclno);
	for (COLONY* cw = cc = clny + 1; cw <= cix; ++cw) {
	    int	nc = 0;
	    if (cw->mark >= 0) {			/* retain */
		nc = cc - clny;
		*cc++ = *cw;
	    }
	    cw->clno = nc;
	}
	vclear(cc, cix - cc + 1);
	no_clny = cc - clny;
}

Colonies::Colonies(int n) : no_clny(0)
{
	if (n == 0) n = OutPrm.NoOut;
	OutPrm.MaxOut = n + MAX_COLONY;
	n = OutPrm.MaxOut + 1;
	clny = new COLONY[n];
	vclear(clny, n);
}
