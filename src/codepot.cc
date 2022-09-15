/*****************************************************************************
*
*	Initiation, termination, and splicing signals
*	and Coding potentials
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

#define IDEBUG	0

#include "aln.h"
#include "utilseq.h"
#include <math.h>

#define	FoldSpjunc	4
#define	FoldStore	2
#define	INTPEN	(INT_MIN / 8 * 6)
#define	BpQueSize	3

enum {Ad, Cy, Gu, Ty};
enum {AA, AC, AG, AT, CA, CC, CG, CT, GA, GC, GG, GT, TA, TC, TG, TT};

INTRONPEN IntronPrm = {FQUERY, FQUERY, -2.767, 20, 224, 
//		     	ip, 	fact,	mean,llmt,  mu, 
	825, 2, 5, 20, 0, 0, 0, 0, 0,
//	rlmt, minexon, tlmt, maxl, array, table, sip
	0.2767, -22.80, 83.35, 5.488, 21.870, 223.95, 0.7882,
//	a1      m1      t1     k1     m2      t2      k2
	0, 0, 0, 0};
//	a2 m3 t3 k3

struct ALSP {
	VTYPE	scr;
	RANGE*	rng;
	INT	eij;
};

struct DefPrm2 {float y, Y;};

// exon-intron junction signals

struct Sig53 {
	STYPE**	sig53tab;
	Sig53(const FTYPE fS = 1, const char* fname = INT53PAT);
	~Sig53();
};

static	double	ProbDist(int i, double mu, double th, double kk);

// [Dvsp][crs]
static	DefPrm2	defprm2[2] = {{4., 4.}, {8., 8.}};
static  float	avrsig53[2] = {2.446, 4.807};

static	Sig53*	stdSig53 = 0;
static	FTYPE	stdfS = 0.;

inline	VTYPE	GapPenalty(int k) {return VTYPE(alprm.u * k + alprm.v);}

VTYPE SpJunc::spjscr(int n5, int n3) {
	return (pwd->IntPen->Penalty(n3 - n5) + 
		b->exin->sig53(n5, n3, IE53));
}

const CHAR* SpJunc::spjseq(int n5, int n3)
{
	if (n5 == b->left || n3 == b->right) return (spj_tron_tab[256]);
static	const	CHAR	pyrim[2] = {PHE, PHE};
const	CHAR*	b5 = b->at(n5 - 2);
const	CHAR*	b3 = n3? b->at(n3): pyrim;
	int	amb = 0;
	int	c = ncredctab[aa2nuc[b5[0]]];
	if (c >= 4) {amb = 1; c = 0;}
	INT	w = c;
	if ((c = ncredctab[aa2nuc[b5[1]]]) < 4) {
	    w = 4 * w + c;
	    if ((c = ncredctab[aa2nuc[b3[0]]]) < 4) {
		w = 4 * w + c;
		if ((c = ncredctab[aa2nuc[b3[1]]]) < 4)
		    w = 4 * w + c;
		else if (amb) w = 256;
		else amb = 2;
	    } else w = 256;
	} else w = 256;
	if (amb == 0 || w == 256) return (spj_tron_tab[w]);
	if (amb == 1) return (spj_amb_tron_tab[w]);
	return (spj_tron_amb_tab[w]);
}

STYPE	Premat::prematT(const CHAR* ps) const
{
	if (bn == 1)
	    return (STYPE) ((*ps == TRM || *ps == TRM2)? fO: 0);
	int	tp = 0;
	for (int n = bn; n--; ++ps) {
	    if (*ps == TRM || *ps == TRM2) ++tp;
	}
	return (STYPE) (fO * tp);
}

Premat::Premat(const Seq* seqs[])
{
	int	swp = !seqs[1]->inex.intr;
	if (swp) swap(seqs[0], seqs[1]);
const	Seq*&	a = seqs[0];
const	Seq*&	b = seqs[1];

	fO = -alprm2.o * a->many * b->many;
	bn = b->many;
	if (swp) swap(seqs[0], seqs[1]);
}

// Frechet Distribution

static double ProbDist(int i, double mu, double th, double kk)
{
	double	z, zz;

	if (i <= mu) return (0.);
	z = th / (i - mu);
	zz = pow(z, kk);
	return (kk / th * z * zz * exp(-zz));
}

IntronPenalty::IntronPenalty(VTYPE f, int dvsp, EijPat* eijpat, ExinPot* exinpot)
	: optip(SHRT_MIN), array(0), table(0)
{
	if (IntronPrm.fact == FQUERY) IntronPrm.fact = defprm2[dvsp > 0].Y;
	if (alprm2.y == FQUERY) alprm2.y = defprm2[dvsp > 0].y;
	if (IntronPrm.ip == FQUERY) IntronPrm.ip = dvsp? 15: 12;
	IntronPrm.sip = (STYPE) (f * IntronPrm.ip);

	if (IntronPrm.maxl <= 0)
	    IntronPrm.maxl = max_intron_len(0.99, 0);
	FTYPE	expsig = 0;
	FTYPE	fy = f * alprm2.y;
	FTYPE	fY = f * IntronPrm.fact;
	if (fy > 0.) {
	    expsig = fy * (1. - alprm2.sss) * avrsig53[0];
	    fy *= alprm2.sss;
	    if (eijpat) {
		FTYPE	fB = f * alprm2.sss * bpprm.factor;
		expsig += fy * (eijpat->pattern5->mmm.mean + eijpat->pattern3->mmm.mean);
		if (eijpat->patternB) expsig += fB * eijpat->patternB->mmm.mean;
	    } else	expsig += fy * avrsig53[1];
	}
	if (exinpot) expsig += exinpot->avrpot(f * alprm2.Z);
	AvrSig = (STYPE) expsig;
	STYPE	IntPen = (STYPE) (expsig + fY * IntronPrm.mean + f * IntronPrm.ip);
	GapWI = (STYPE) (fY * IntronPrm.mean) - IntPen;
	int	i = IntronPrm.rlmt - IntronPrm.llmt + 1;
	array = new STYPE[i];
	table = array - IntronPrm.llmt;
	STYPE*	pentab = array;
	float	a2 = IntronPrm.a2? IntronPrm.a2: 1 - IntronPrm.a1;
	float	a3 = IntronPrm.a2? 1 - IntronPrm.a1 - IntronPrm.a2: 0;

	VTYPE	gep = f * alprm.u;
	VTYPE	gappen = -(f * alprm.v + IntronPrm.llmt * gep);
	IntronPrm.minl = 0;
	for (i = IntronPrm.llmt; i <= IntronPrm.rlmt; ++i) {
	    double z = ProbDist(i, IntronPrm.m1, IntronPrm.t1, IntronPrm.k1);
	    if (a2 > 0.) {
		z *= IntronPrm.a1;
		z += a2 * ProbDist(i, IntronPrm.m2, IntronPrm.t2, IntronPrm.k2);
		if (a3) 
		    z += a3 * ProbDist(i, IntronPrm.m3, IntronPrm.t3, IntronPrm.k3);
	    }
	    STYPE	gp = (STYPE) (fY * log10(z)) - IntPen;
	    *pentab++ = gp;
	    if (gp > optip) {
		optip = gp;
		IntronPrm.mode = i;
	    }
	    if (!IntronPrm.minl) {
		if (gp > gappen) IntronPrm.minl = i;
		else	gappen -= gep;
	    }
	}
	if (!IntronPrm.minl) IntronPrm.minl = IntronPrm.llmt;
	int	k = 1;
	double	z = ProbDist(IntronPrm.rlmt, IntronPrm.m1, IntronPrm.t1, IntronPrm.k1);
	double	zz = a2? ProbDist(IntronPrm.rlmt, IntronPrm.m2, IntronPrm.t2, IntronPrm.k2): 0;
	if (zz > z) {k = 2; z = zz;}
	zz = a3? ProbDist(IntronPrm.rlmt, IntronPrm.m3, IntronPrm.t3, IntronPrm.k3): 0;
	if (zz > z) {k = 3; z = zz;}
	if (k == 1) {
	    IntronPrm.mu = (int) IntronPrm.m1;
	    z = IntronPrm.k1;
	} else if (k == 2) {
	    IntronPrm.mu = (int) IntronPrm.m2;
	    z = IntronPrm.k2;
	} else {
	    IntronPrm.mu = (int) IntronPrm.m3;
	    z = IntronPrm.k3;
	}
	IntEp = (FTYPE) (-(z + 1) * fY / log(10.));
	IntFx = *--pentab - IntEp * log((double) (IntronPrm.rlmt - IntronPrm.mu));
#if IDEBUG
	for (i = IntronPrm.llmt; i < IntronPrm.rlmt + 10; ++i)
	    printf("%d\t%7.2f\n", i, (double) Penalty(i));
#endif
}

EijPat::EijPat(int dvsp)
{
const	char*	fn;

	if (alprm2.y == FQUERY) alprm2.y = defprm2[dvsp > 0].y;
	if (alprm2.y > 0.) {
	    pattern5 = new PatMat(fn = SPLICE5PAT);
	    if (!pattern5->mtx) fatal("Can't open %s file!", fn);
	    pattern3 = new PatMat(fn = SPLICE3PAT);
	    if (!pattern3->mtx) fatal("Can't open %s file!", fn);
	    tonic3 = (FTYPE) pattern3->tonic;
	    tonic5 = (FTYPE) pattern5->tonic;
	    pattern3->tonic = pattern5->tonic = 0;
	} else {
	    pattern5 = pattern3 = 0;
	    tonic3 = tonic5 = 0;
	}
	if (alprm2.bti > 0.) {
	    patternI = new PatMat(fn = INITIATPAT);
	    if (!patternI->mtx) fatal("Can't open %s file!", fn);
	    patternT = new PatMat(fn = TERMINPAT);
	    if (!patternT->mtx) fatal("Can't open %s file!", fn);
	} else {
	    patternI = patternT = 0;
	}
	if (bpprm.factor > 0.) {
	    patternB = new PatMat(fn = BRANCHPAT);
	    if (!patternB->mtx) fatal("Can't open %s file!", fn);
	    tonicB = (FTYPE) patternB->tonic;
	} else {
	    patternB = 0;
	    tonicB = 0;
	}
}

EijPat::~EijPat()
{
	delete pattern5;
	delete pattern3;
	delete patternI;
	delete patternT;
	delete patternB;
}

Sig53::Sig53(const FTYPE fS, const char* fname)
{
	sig53tab = new STYPE*[4];
	sig53tab[0] = new STYPE[544];
	sig53tab[1] = sig53tab[0] + 16;
	sig53tab[2] = sig53tab[1] + 16;
	sig53tab[3] = sig53tab[2] + 256;
	PatMat*	ptmt = 0;
const	float	fs = fS * (1. - alprm2.sss);

	FILE*	fd = ftable.fopen(fname, "r");
	if (!fd) goto abort;

	ptmt = new PatMat(fd);	/* INT5PAT */
	if (!ptmt->mtx) goto abort;
	for (int i = 0; i < 16; ++i)
	    sig53tab[0][i] = (STYPE) (fs * ptmt->mtx[i]);
	delete ptmt;

	ptmt = new PatMat(fd);	/* INT3PAT */
	if (!ptmt->mtx) goto abort;
	for (int i = 0; i < 16; ++i)
	    sig53tab[1][i] = (STYPE) (fs * ptmt->mtx[i]);
	delete ptmt;

	ptmt = new PatMat(fd);	/* INT53PAT */
	if (!ptmt->mtx) goto abort;
	for (int i = 0; i < 256; ++i)
	    sig53tab[2][i] = (STYPE) (fs * ptmt->mtx[i]);
	delete ptmt;

	ptmt = new PatMat(fd);	/* INT35PAT */
	if (!ptmt->mtx) goto abort;
	for (int i = 0; i < 256; ++i)
	    sig53tab[3][i] = (STYPE) (fs * ptmt->mtx[i]);
	delete ptmt;
	fclose(fd);
	return;
abort:
	if (fd) fclose(fd);
	delete[] *sig53tab;
	delete[] sig53tab;
	delete ptmt;
}

STYPE** new53tab(const FTYPE fS, const Sig53* srcSig53)
{
	STYPE**	newtab = new STYPE*[4];
	newtab[0] = new STYPE[544];
	newtab[1] = newtab[0] + 16;
	newtab[2] = newtab[1] + 16;
	newtab[3] = newtab[2] + 256;

	STYPE*	s = srcSig53->sig53tab[0];
	STYPE*	d = newtab[0];
	STYPE*	t = d + 544;
	while (d < t) *d++ = (STYPE) (fS * *s++);
	return (newtab);
}

Sig53::~Sig53()
{
	if (sig53tab) {
	    delete[] *sig53tab;
	    delete[] sig53tab;
	}
}

void EraStdSig53() {delete stdSig53;}

void makeStdSig53()
{
	if (!stdSig53) {
	    stdfS = (FTYPE) (alprm2.y * alprm.scale);
	    stdSig53 = new Sig53(stdfS);
	}
}

Exinon::Exinon(Seq* sd_, const PwdB* pwd_, bool bo) :
	sd(sd_), pwd(pwd_), both_ori(bo),
	size(sd->right - sd->left + 2), bias(sd->left - 1), 
	intnpot(pwd? pwd->intnpot: 0), 
	exonpot(pwd? pwd->exonpot: 0),
	fact((FTYPE) (pwd->Vab / sd->many)),
	fS(alprm2.y * fact)
{
	if (!stdSig53) fatal("call makeStdSig53 beforehand !\n");
const	INT53	zero53 = {0, 0, 0, 0};
	int53 = new INT53[size + 1];
	int53[0] = int53[size] = zero53;
	int53 -= bias;
	data = new EXIN[size + 1] - bias;

	if (fS == stdfS) {
	    sig53tab = stdSig53->sig53tab;
	    statictab = true;
	} else {
	    sig53tab = new53tab(fS / stdfS, stdSig53);
	    statictab = false;
	}
	intron53N();
	intron53();
}

Exinon::~Exinon()
{
	if (data) 	delete[] (data + bias);
	if (int53)	delete[] (int53 + bias);
	if (!statictab && sig53tab) {
	    delete[] *sig53tab;
	    delete[] sig53tab;
	}
}

STYPE Exinon::sig53(int m, int n, INTENDS c) const
{
	STYPE	sig = 0;
	switch (c) {		// acanonical signal
	    case IE5:	sig = data[m].sig5; break;
	    case IE3:	sig = data[n].sig3; break;
	    case IE53:	sig = data[n].sig3 - sig53tab[1][int53[n].dinc3]
		+ sig53tab[2][16 * int53[m].dinc5 + int53[n].dinc3];
		break;
	    case IE35:	sig = data[m].sig5 -  sig53tab[0][int53[m].dinc5]
		+ sig53tab[3][16 * int53[m].dinc5 + int53[n].dinc3];
		break;
	    case IE5P3: sig = data[m].sig5 + data[n].sig3 - sig53tab[1][int53[n].dinc3]
		+ sig53tab[2][16 * int53[m].dinc5 + int53[n].dinc3];
		break;
	}
	if (alprm2.Z > 0) {
	    if (c == IE53 || c == IE5P3) sig += intnpot->intpot(data + m, data + n);
	    else if (c == IE35) sig += intnpot->intpot(data + n, data + m);
	}
	return (sig);
}

void Exinon::intron53N()
{
const	static	INT	jlevelac[4] = {0, 2, 3, 1};
const	static	INT	jlevelgt[4] = {0, 0, 3, 1};
	int	nc = 1;		/* reduced numeric code for 'C' */
	CHAR*	redctab = sd->istron()? tnredctab: ncredctab;

	INT53*	wk5 = int53 + bias;
	INT53*	wk3 = wk5 + 2;
	CHAR*	ss = sd->at(sd->left);
	CHAR*	ts = sd->at(sd->right);
	for ( ; ss < ts; ++ss, ++wk5, ++wk3) {
	    int	c = redctab[*ss];
	    if (c >= 4) c = 1;	// 'C'
	    nc = ((nc << 2) + c) & 0xf;
	    wk5->dinc5 = wk3->dinc3 = nc;
	    wk5->cano5 = wk3->cano3 = algmode.any == 3? 1: 0;
	    switch (nc) {
		case AA: wk3->cano3 = jlevelac[algmode.any]; break;
		case AC: wk3->cano3 = 2;
			 if (both_ori) wk5->cano5 = 1;
			 break;
		case AG: wk3->cano3 = 3; break;
		case AT: wk5->cano5 = 2; wk3->cano3 = jlevelac[algmode.any]; break;
		case CG: wk3->cano3 = jlevelgt[algmode.any]; break;
		case CT: wk5->cano5 = jlevelgt[algmode.any];
			 if (both_ori) wk3->cano3 = 1;
			 break;
		case GA: wk5->cano5 = jlevelgt[algmode.any]; break;
		case GC: wk5->cano5 = 3; break;
		case GG: wk5->cano5 = jlevelgt[algmode.any];
			 wk3->cano3 = jlevelgt[algmode.any]; break;
		case GT: wk5->cano5 = 3;
			 if (both_ori) wk3->cano3 = 1;
			 break;
		case TG: wk3->cano3 = jlevelgt[algmode.any]; break;
		case TT: wk5->cano5 = jlevelgt[algmode.any]; break;
		default: break;
	    }
	}
}

void Exinon::intron53()
{
	FTYPE	fE = alprm2.z * fact;
	FTYPE	fI = alprm2.Z * fact;
	FTYPE	fT = alprm2.bti * fact;
	FTYPE	fB = bpprm.factor * fact;
	FTYPE	fO = -alprm2.o * fact;

	--sd->left; ++sd->right;
	STYPE	th3 = (STYPE) (fS * pwd->eijpat->tonic3);
	STYPE	th5 = (STYPE) (fS * pwd->eijpat->tonic5);
	STYPE	thB = (STYPE) (pwd->eijpat->tonicB);
	float*	pref5 = pwd->eijpat->pattern5? pwd->eijpat->pattern5->calcPatMat(sd): 0;
	float*	pref3 = pwd->eijpat->pattern3? pwd->eijpat->pattern3->calcPatMat(sd): 0;
	float*	prefS = pwd->eijpat->patternI? pwd->eijpat->patternI->calcPatMat(sd): 0;
	float*	prefT = pwd->eijpat->patternT? pwd->eijpat->patternT->calcPatMat(sd): 0;
	float*	prefB = pwd->eijpat->patternB? pwd->eijpat->patternB->calcPatMat(sd): 0;
	float*	prefE = pwd->codepot? pwd->codepot->calcScr(sd): 0;
	float*	prefI = pwd->intnpot? pwd->intnpot->calcScr(sd): 0;
	if (!prefE && pwd->exonpot) prefE = pwd->exonpot->calcScr(sd);
	float*	prf5 = pref5;	// 5'ss signal
	float*	prf3 = pref3;	// 3'ss signal
	float*	prfS = prefS;	// start codon
	float*	prfT = prefT;	// termination codon
	float*	prfB = prefB;	// branch site signal
	float*	prfE = prefE;	// coding potential
	float*	prfI = prefI;	// intron potential
	++sd->left; --sd->right;
	if (pref5)	++prf5;
	if (pref3)	++prf3;
	if (prefS)	++prfS;
	if (prefT)	++prfT;
	if (prefB)	++prfB;
	if (prefE)	++prfE;
	if (prefI)	++prfI;

	vset(data + bias, ZeroExin, size + 1);
	at_sig5 = sig53tab[0][3];
	gc_sig5 = sig53tab[0][9];
	CHAR*	ss = sd->at(sd->left);
	STYPE	sigB = 0;
	CHAR*	posB = 0;
	EXIN*	last = end();
	INT53*	wk53 = int53 + bias + 1;
const	float	fs = fS * alprm2.sss;
	for (EXIN* wkb = begin(); ++wkb < last; ++ss, ++wk53) {
	    wkb->sigS = prefS? (STYPE) (fT * *prfS++): 0;
	    wkb->sigT = prefT? (STYPE) (fT * *prfT++): 0;
	    wkb->sigI = prefI? (STYPE) (fI * *prfI++): 0;
	    if (prfE) {
		FTYPE	sigE = fE * *prfE++;
		if (*ss == TRM || *ss == TRM2) sigE += fO;
		else if (wkb + 3 < last && (ss[3] == TRM || ss[3] == TRM2)) sigE = 0;
		wkb->sigE = (STYPE) sigE;
	    } else	wkb->sigE = 0;
	    STYPE	sig5 = pref5? (STYPE) (fs * *prf5++): 0;
	    STYPE	sig3 = pref3? (STYPE) (fs * *prf3++): 0;
	    sig5 += sig53tab[0][wk53->dinc5];
	    sig3 += sig53tab[1][wk53->dinc3];
	    if (prefB) {
		sig3 += sigB;
		FTYPE	sigb = *prfB++;
		if (sigb > thB) {	// stronger than given threshold
		    sigB = (STYPE) (fB * sigb);
		    posB = ss;
		}	// reset bp signal if bp-3'ss distance exeeds the threshold
		if (posB && ss - posB > bpprm.maxb3d) {
		    sigB = 0;
		    posB = 0;
		}
	    }
	    wkb->sig5 = sig5;
	    if (wkb->phs5 == -2 && ((algmode.any == 2 && sig5 > th5) || wk53->cano5)) {
		wkb->phs5 = 0;
		if (wk53->cano5 > 1) {
		    wkb[1].phs5 = 1;
		    if (wkb[-1].phs5 == 1) wkb[-1].phs5 = 2;	// GTGT
		    else wkb[-1].phs5 = -1;
		}
	    }
	    wkb->sig3 = sig3;
	    if (wkb->phs3 == -2 && ((algmode.any == 2 && sig3 > th3) || wk53->cano3)) {
		wkb->phs3 = 0;
		if (wk53->cano3 > 1) {
		    wkb[1].phs3 = 1;
		    if (wkb[-1].phs3 == 1) wkb[-1].phs3 = 2;	// AGAG
		    else wkb[-1].phs3 = -1;
		}
	    }
	}
	delete[] pref5; delete[] pref3; delete[] prefS; delete[] prefT;
	delete[] prefB; delete[] prefE; delete[] prefI;
}

void Exinon::resize() {
	if (bias < sd->left && sd->right < (bias + size)) return;
	RANGE	rng;
	sd->saverange(&rng);
	if (data) 	delete[] (data + bias);
	if (int53)	delete[] (int53 + bias);
	if (bias + 1 < sd->left) sd->left = bias + 1;
	if (sd->left + size - 2 > sd->right) sd->right = sd->left + size - 2;
	bias = sd->left - 1;
	size = sd->right - sd->left + 2;
	int53 = new INT53[size + 1];
const	INT53	zero53 = {0, 0, 0, 0};
	int53[0] = int53[size] = zero53;
	int53 -= bias;
	data = new EXIN[size + 1] - bias;
	intron53N();
	intron53();
	sd->restrange(&rng);
}

STYPE IntronPenalty::Penalty(int n) const 
{
	if (n < 0 || !array) return (GapWI);
	if (n < IntronPrm.llmt) return (SHRT_MIN);
	if (n >= IntronPrm.rlmt)
	    return (STYPE) (IntFx + IntEp * log((double)(n - IntronPrm.mu)));
	return (table[n]);
}

STYPE IntronPenalty::Penalty(int n, bool addsig) const 
{
	if (n < IntronPrm.llmt) return (SHRT_MIN);
	if (!array) return (GapWI + AvrSig);
	if (n >= IntronPrm.rlmt)
	    return (STYPE) (IntFx + IntEp * log((double)(n - IntronPrm.mu)) + AvrSig);
	return (table[n] + AvrSig);
}

INT max_intron_len(float p, const char* fn)
{
	INT	observed = 0;
	double	mu = IntronPrm.m2;
	double	th = IntronPrm.t2;
	double	ki = IntronPrm.k2;
	if (IntronPrm.a2 > 0.) {
	    mu = IntronPrm.m3;
	    th = IntronPrm.t3;
	    ki = IntronPrm.k3;
	} else if (IntronPrm.a1 == 0) {
	    mu = IntronPrm.m1;
	    th = IntronPrm.t1;
	    ki = IntronPrm.k1;
	}
	if (fn) {
	    FILE*	fd = ftable.fopen(ipstat, "r");
	    if (!fd) goto ild_model;
	    const	char*	sl = strrchr(fn, '/');
	    if (sl) fn = sl + 1;
	    char	str[LINE_MAX];
	    while (fgets(str, LINE_MAX, fd)) {
		if (strncmp(str, fn, 8)) continue;
		Strlist	stl(str, stddelim);
		if (atoi(stl[2]) < 1000) break;
		observed = atoi(stl[5]);
		int	n = stl.size() - 6;
		mu = atof(stl[n]);
		th = atof(stl[n + 1]);
		ki = atof(stl[n + 2]);
		break;
	    }
	    fclose(fd);
	}
ild_model:
	INT	q99 = (INT) frechet_quantile(p, mu, th, ki);
	return (max(observed, q99));
}
