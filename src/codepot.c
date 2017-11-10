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
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#define EDEBUG	0
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

//		     ip, fact, mean, llmt, mu, rlmt, minexon, tlmt, maxl, array, table
INTRONPEN IntronPrm = {DQUERY, DQUERY, -2.767, 20, 224, 825, 2, 5, 0, 0, 0,
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
	VTYPE**	sig53tab;
	Sig53(const FTYPE fS = 1, const char* fname = INT53PAT);
	~Sig53();
};

static	double	ProbDist(int i, double mu, double th, double kk);

// [Dvsp][crs]
static	DefPrm2	defprm2[2][2] = {{{4., 4.}, {4., 4.}}, {{5., 5.}, {8., 8.}}};
static  float	avrsig53[2] = {2.446, 4.807};
static	const	float	imaxfact = 5.;

static	Sig53*	stdSig53 = 0;
static	FTYPE	stdfS = 0.;

inline	VTYPE	GapPenalty(int k) {return VTYPE(alprm.u * k + alprm.v);}

SpJunc::SpJunc(Seq* sd)
{
static	SPJ	nullspj = {0, 0, 0, 0, 0};

	b = sd;
	nent = (sd->right - sd->left) / FoldSpjunc;
	nent = (nent / 6) * 6 + 5;
	max_nent = FoldStore * nent;
	hashent = new SPJ*[nent];
	spjunc = new SPJ[max_nent];
	juncseq = new CHAR[2 * max_nent];
	spliced = new CHAR[6];
		 /* 2 trons  + 4 bases */
	pyrim = new CHAR[2];
	for (int i = 0; i < nent; ++i) hashent[i] = 0;
	for (int i = 0; i < max_nent; ++i) spjunc[i] = nullspj;
	for (int i = 0; i < 2; ) pyrim[i++] = PHE;	/* TT <- YY */
	spp = spjunc;
	jsp = juncseq;
}

SpJunc::~SpJunc()
{
	delete[] hashent;
	delete[] spjunc;
	delete[] juncseq;
	delete[] spliced;
	delete[] pyrim;
}

CHAR* SpJunc::spjseq(int n5, int n3)
{
	int     cc = n5 + n3;
	SPJ*    find;
	SPJ**   htop = hashent + (cc % nent);
	CHAR*   b5 = b->at(n5 - 2);
	CHAR*   b3 = b->at(n3);

	for (find = *htop; find; find = find->dlnk)
	    if (find->n5 == n5 && find->n3 == n3)
		return(find->junc);

	if (spp - spjunc == max_nent) {
	    spp = spjunc;
	    jsp = juncseq;
	}

	find = spp++;
	if (find->ulnk) *find->ulnk = 0;        /* cut the previous link */
	find->n5 = n5;
	find->n3 = n3;
	spliceTron(spliced, b5, n3? b3: pyrim, 1);
	memcpy(jsp, spliced, 2);
	find->junc = jsp;
	jsp += 2;
	find->dlnk = *htop;     /* old 1st link */
	find->ulnk = htop;      /* put the new record */
	*htop = find;           /* at the top of list */
	return (find->junc);
}

inline	bool	Exinon::within(Seq* sd)
{
	return (bias < sd->left && sd->right < (bias + size));
}

VTYPE	Premat::prematT(const CHAR* ps)
{
	if (bn == 1)
	    return (VTYPE) ((*ps == TRM || *ps == TRM2)? fO: 0);
	int	tp = 0;
	for (int n = bn; n--; ++ps) {
	    if (*ps == TRM || *ps == TRM2) ++tp;
	}
	return (VTYPE) (fO * tp);
}

Premat::Premat(Seq* seqs[])
{
	int	swp = !seqs[1]->inex.intr;
	if (swp) swapseq(seqs, seqs + 1);
	Seq*	a = seqs[0];
	Seq*	b = seqs[1];

	fO = -alprm2.o * a->many * b->many;
	bn = b->many;
	if (swp) swapseq(seqs, seqs + 1);
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

IntronPenalty::IntronPenalty(VTYPE f, int hh, EijPat* eijpat, ExinPot* exinpot)
	: array(0), table(0), optlen(0)
{
	if (IntronPrm.fact == DQUERY) IntronPrm.fact = defprm2[hh > 0][algmode.crs].Y;
	if (alprm2.y == DQUERY) alprm2.y = defprm2[hh > 0][algmode.crs].y;

	if (IntronPrm.maxl <= 0)
	    IntronPrm.maxl = int (imaxfact * (IntronPrm.t2 > 0.? IntronPrm.t2: IntronPrm.t1));
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
	if (exinpot) expsig +=  exinpot->avrpot(f * alprm2.Z);
	AvrSig = (VTYPE) expsig;
	VTYPE	IntPen = (VTYPE) (expsig + fY * IntronPrm.mean +
		f * ((IntronPrm.ip == DQUERY)? GapPenalty(Ip_equ_k): IntronPrm.ip));
	GapWI = (VTYPE) (fY * IntronPrm.mean) - IntPen;
	int	i = IntronPrm.rlmt - IntronPrm.llmt + 1;
	array = new VTYPE[i];
	table = array - IntronPrm.llmt;
	VTYPE*	pentab = array;

	float	a2 = IntronPrm.a2? IntronPrm.a2: 1 - IntronPrm.a1;
	float	a3 = IntronPrm.a2? 1 - IntronPrm.a1 - IntronPrm.a2: 0;
	VTYPE	optgp = NEVSEL;
	for (i = IntronPrm.llmt; i <= IntronPrm.rlmt; ++i) {
	    double z = ProbDist(i, IntronPrm.m1, IntronPrm.t1, IntronPrm.k1);
	    if (a2 > 0.) {
		z *= IntronPrm.a1;
		z += a2 * ProbDist(i, IntronPrm.m2, IntronPrm.t2, IntronPrm.k2);
		if (a3) 
		    z += a3 * ProbDist(i, IntronPrm.m3, IntronPrm.t3, IntronPrm.k3);
	    }
	    VTYPE	gp = (VTYPE) (fY * log10(z)) - IntPen;
	    *pentab++ = gp;
	    if (gp > optgp) {
		optgp = gp;
		optlen = i;
	    }
	}
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

EijPat::EijPat(int hh)
{
const	char*	fn;

	if (alprm2.y == DQUERY) alprm2.y = defprm2[hh > 0][algmode.crs].y;
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

Exinon::Exinon(Seq* sd, FTYPE ff, PwdB* pwd)
{
static	const	INT53	zero53 = {0, 0, 0, 0};
	size = sd->right - sd->left + 2;
	bias = sd->left - 1;
	fact = ff;
	sig53tab = 0;		// +1 for dummy
	int53 = new INT53[size + 1];
	int53[0] = int53[size] = zero53;
	int53 -= bias;
	data = (pwd->DvsP || alprm2.sss > 0.)? new EXIN[size + 1] - bias: 0;
	exinpot = pwd? pwd->exinpot: 0;
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

VTYPE Exinon::sig53(int m, int n, INTENDS c)
{
	VTYPE	sig = 0;
	if (alprm2.sss < 1.) {
	  switch (c) {		// canonical signal
	    case IE5:	sig = sig53tab[0][int53[m].dinc5]; break;
	    case IE3:	sig = sig53tab[1][int53[n].dinc3]; break;
	    case IE53:	sig = sig53tab[2][16 * int53[m].dinc5 + int53[n].dinc3];  break;
	    case IE35:	sig = sig53tab[3][16 * int53[m].dinc5 + int53[n].dinc3]; break;
	    case IE5P3:	sig = sig53tab[0][int53[m].dinc5] +
			sig53tab[2][16 * int53[m].dinc5 + int53[n].dinc3]; break;
	  }
	}
	if (!data) return sig;
	VTYPE	sss = 0;	// species specific signal
	switch (c) {
	    case IE5:	sss = data[m].sig5; break;
	    case IE35:	sss = data[m].sig5; break;
	    case IE3:	sss = data[n].sig3; break;
	    case IE53:	sss = data[n].sig3; break;
	    case IE5P3:	sss = data[m].sig5 + data[n].sig3; break;
	}
	sig = (VTYPE) ((1. - alprm2.sss) * sig + alprm2.sss * sss);
	if (alprm2.Z > 0) {
	    if (c == IE53 || c == IE5P3) sig += exinpot->intpot(data + m, data + n);
	    else if (c == IE35) sig += exinpot->intpot(data + n, data + m);
	}
	return (sig);
}

Sig53::Sig53(const FTYPE fS, const char* fname)
{
	sig53tab = new VTYPE*[4];
	sig53tab[0] = new VTYPE[544];
	sig53tab[1] = sig53tab[0] + 16;
	sig53tab[2] = sig53tab[1] + 16;
	sig53tab[3] = sig53tab[2] + 256;
	PatMat*	ptmt = 0;

	FILE*	fd = ftable.fopen(fname, "r");
	if (!fd) goto abort;

	ptmt = new PatMat(fd);	/* INT5PAT */
	if (!ptmt->mtx) goto abort;
	for (int i = 0; i < 16; ++i)
	    sig53tab[0][i] = (VTYPE) (fS * ptmt->mtx[i]);
	delete ptmt;

	ptmt = new PatMat(fd);	/* INT3PAT */
	if (!ptmt->mtx) goto abort;
	for (int i = 0; i < 16; ++i)
	    sig53tab[1][i] = (VTYPE) (fS * ptmt->mtx[i]);
	delete ptmt;

	ptmt = new PatMat(fd);	/* INT53PAT */
	if (!ptmt->mtx) goto abort;
	for (int i = 0; i < 256; ++i)
	    sig53tab[2][i] = (VTYPE) (fS * ptmt->mtx[i]);
	delete ptmt;

	ptmt = new PatMat(fd);	/* INT35PAT */
	if (!ptmt->mtx) goto abort;
	for (int i = 0; i < 256; ++i)
	    sig53tab[3][i] = (VTYPE) (fS * ptmt->mtx[i]);
	delete ptmt;
	fclose(fd);
	return;
abort:
	if (fd) fclose(fd);
	delete[] *sig53tab;
	delete[] sig53tab;
	delete ptmt;
}

VTYPE** new53tab(const FTYPE fS, const Sig53* srcSig53)
{
	VTYPE**	newtab = new VTYPE*[4];
	newtab[0] = new VTYPE[544];
	newtab[1] = newtab[0] + 16;
	newtab[2] = newtab[1] + 16;
	newtab[3] = newtab[2] + 256;

	VTYPE*	s = srcSig53->sig53tab[0];
	VTYPE*	d = newtab[0];
	VTYPE*	t = d + 544;
	while (d < t) *d++ = (VTYPE) (fS * *s++);
	return (newtab);
}

Sig53::~Sig53()
{
	if (sig53tab) {
	    delete[] *sig53tab;
	    delete[] sig53tab;
	}
}

// Intron53N

void EraStdSig53() {delete stdSig53;}

void makeStdSig53()
{
	if (!stdSig53) {
	    stdfS = (FTYPE) (alprm2.y * alprm.scale);
	    stdSig53 = new Sig53(stdfS);
	}
}

void Intron53N(Seq* sd, FTYPE ff, PwdB* pwd)
{
const	static	INT	jlevelac[4] = {0, 2, 3, 1};
const	static	INT	jlevelgt[4] = {0, 0, 3, 1};
	int	nc = 1;		/* reduced numeric code for 'C' */
	CHAR*	redctab = sd->istron()? tnredctab: ncelements;
	FTYPE	fS = alprm2.y * ff;
	Exinon*	exin = sd->exin;

	if (!stdSig53) fatal("call makeStdSig53 beforehand !\n");
	delete exin;
	exin = sd->exin = new Exinon(sd, ff, pwd);
	if (fS == stdfS) {
	    exin->sig53tab = stdSig53->sig53tab;
	    exin->statictab = true;
	} else {
	    exin->sig53tab = new53tab(fS / stdfS, stdSig53);
	    exin->statictab = false;
	}

	INT53*	wk5 = exin->int53 + exin->bias;
	INT53*	wk3 = wk5 + 2;
	CHAR*	ss = sd->at(sd->left);
	CHAR*	ts = sd->at(sd->right);
	for ( ; ss < ts; ++ss, wk5++, wk3++) {
	    nc = ((nc << 2) + redctab[*ss]) & 0xf;
	    wk5->dinc5 = wk3->dinc3 = nc;
	    wk5->cano5 = wk3->cano3 = algmode.any == 3? 1: 0;
	    switch (nc) {
		case AA: wk3->cano3 = jlevelac[algmode.any]; break;
		case AC: wk3->cano3 = 2; break;
		case AG: wk3->cano3 = 3; break;
		case AT: wk5->cano5 = 2; wk3->cano3 = jlevelac[algmode.any]; break;
		case CG: wk3->cano3 = jlevelgt[algmode.any]; break;
		case CT: wk5->cano5 = jlevelgt[algmode.any]; break;
		case GA: wk5->cano5 = jlevelgt[algmode.any]; break;
		case GC: wk5->cano5 = 3; break;
		case GG: wk5->cano5 = jlevelgt[algmode.any];
			 wk3->cano3 = jlevelgt[algmode.any]; break;
		case GT: wk5->cano5 = 3; break;
		case TG: wk3->cano3 = jlevelgt[algmode.any]; break;
		case TT: wk5->cano5 = jlevelgt[algmode.any]; break;
		default: break;
	    }
	}
}

void Intron53(Seq* sd, PwdB* pwd)
{
	FTYPE	ff = (FTYPE) (pwd->Vab / sd->many);
	if (sd->exin && sd->exin->int53 && ff == sd->exin->fact && sd->exin->within(sd))
	    return;
	Intron53N(sd, ff, pwd);
	if (pwd->DvsP == 0 && alprm2.sss <= 0.) return;

	FTYPE	fE = alprm2.z * ff;
	FTYPE	fI = alprm2.Z * ff;
	FTYPE	fS = alprm2.y * ff;
	FTYPE	fT = alprm2.bti * ff;
	FTYPE	fB = bpprm.factor * ff;
	FTYPE	fO = -alprm2.o * ff;
#if EDEBUG
	int	i;
	char*	ps;
	int	trn = !isdrna(sd);
#endif

	--sd->left; ++sd->right;
	VTYPE	th3 = (VTYPE) (-fS * pwd->eijpat->tonic3);
	VTYPE	th5 = (VTYPE) (-fS * pwd->eijpat->tonic5);
	float*	pref5 = pwd->eijpat->pattern5? pwd->eijpat->pattern5->calcPatMat(sd): 0;
	float*	pref3 = pwd->eijpat->pattern3? pwd->eijpat->pattern3->calcPatMat(sd): 0;
	float*	prefS = pwd->eijpat->patternI? pwd->eijpat->patternI->calcPatMat(sd): 0;
	float*	prefT = pwd->eijpat->patternT? pwd->eijpat->patternT->calcPatMat(sd): 0;
	float*	prefB = pwd->eijpat->patternB? pwd->eijpat->patternB->calcPatMat(sd): 0;
	float*	prefE = pwd->codepot? pwd->codepot->calcPrefCodePot(sd, 0): 0;
	float*	prefI = pwd->exinpot? pwd->exinpot->calcExinPot(sd, false): 0;
	if (!prefE && pwd->exinpot) prefE = pwd->exinpot->calcExinPot(sd, true);
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

	EXIN*	wkb = sd->exin->data + sd->exin->bias;
	EXIN*	last = wkb + sd->exin->size - 1;
	memset(wkb, '\0', sizeof(EXIN));
	memset(last, '\0', sizeof(EXIN)); 
	for ( ; wkb <= last; ++wkb) wkb->phs5 = wkb->phs3 = -2;
	CHAR*	ss = sd->at(sd->left);
	VTYPE	sigB = 0;
	CHAR*	posB = 0;
	INT53*	wk53 = sd->exin->int53 + sd->exin->bias + 1;
	for (wkb = sd->exin->data + sd->exin->bias; ++wkb < last; ++ss, ++wk53) {
	    wkb->sigS = (VTYPE) (fT * *prfS++);
	    wkb->sigT = (VTYPE) (fT * *prfT++);
	    if (prfE) {
		FTYPE	sigE = fE * *prfE++;
		if (*ss == TRM || *ss == TRM2) sigE += fO;
		else if (wkb + 3 < last && (ss[3] == TRM || ss[3] == TRM2)) sigE = 0;
		wkb->sigE = (VTYPE) sigE;
	    } else	wkb->sigE = 0;
	    wkb->sigI = prfI? (VTYPE) (fI * *prfI++): 0;
	    VTYPE	sig5 = pref5? (VTYPE) (fS * *prf5++): 0;
	    VTYPE	sig3 = pref3? (VTYPE) (fS * *prf3++): 0;
	    if (prefB) {
		sig3 += sigB;
		FTYPE	sigb = *prfB++;
		if (sigb + pwd->eijpat->tonicB > 0) {	// stronger than given threshold
		    sigB = (VTYPE) (fB * sigb);
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
		    if (wkb[-1].phs5 == -2) wkb[-1].phs5 = -1;
		    else wkb[-1].phs5 = 2;
		}
	    }
	    wkb->sig3 = sig3;
	    if (wkb->phs3 == -2 && ((algmode.any == 2 && sig3 > th3) || wk53->cano3)) {
		wkb->phs3 = 0;
		if (wk53->cano3 > 1) {
		    wkb[1].phs3 = 1;
		    if (wkb[-1].phs3 == 1) wkb[-1].phs3 = 2;
		    else wkb[-1].phs3 = -1;
/* modified to accomodate the case of ..AGAG.. */
		}
	    }
	}
	delete[] pref5; delete[] pref3; delete[] prefS; delete[] prefT;
	delete[] prefB; delete[] prefE; delete[] prefI;
}

VTYPE IntronPenalty::Penalty(int n)
{
	if (n < 0 || !array) return (GapWI);
	if (n < IntronPrm.llmt) return (NEVSEL);
	if (n >= IntronPrm.rlmt)
	    return (VTYPE) (IntFx + IntEp * log((double)(n - IntronPrm.mu)));
	return (table[n]);
}

VTYPE IntronPenalty::Penalty(int n, bool addsig)
{
	if (n < IntronPrm.llmt) return (NEVSEL);
	if (!array) return (GapWI + AvrSig);
	if (n >= IntronPrm.rlmt)
	    return (VTYPE) (IntFx + IntEp * log((double)(n - IntronPrm.mu)) + AvrSig);
	return (table[n] + AvrSig);
}

