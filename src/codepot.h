/*****************************************************************************
*
*	Header to codepot.c
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

#ifndef  _CODEPOT_H_
#define  _CODEPOT_H_

class	ExinPot;

struct EXIN {
	STYPE   sig5;
	STYPE   sig3;
	STYPE   sigS;
	STYPE   sigT;
	STYPE   sigE;
	STYPE   sigI;
	char   phs5;
	char   phs3;
};

static	const	EXIN	ZeroExin = {0, 0, 0, 0, 0, 0, -2, -2};

struct INT53 {
	INT	dinc5:	4;
	INT	dinc3:	4;
	INT	cano5:	4;
	INT	cano3:	4;
};

class Premat {
	FTYPE   fO;
	int     bn;
public:
	Premat(const Seq* seqs[]);
	STYPE	prematT(const CHAR* ps) const;
	STYPE	prematT1(const CHAR* ps) const {
	    return (STYPE) ((*ps == TRM || *ps == TRM2)? fO: 0);
	}
};

struct PwdB;
struct EijPat;

enum INTENDS {IE5, IE3, IE53, IE35, IE5P3};

struct SPJ {
	int	n5;
	int	n3;
	INT	junc;	// CHAR* change from ptr
	INT	dlnk;	// SPJ*  to index
	INT	ulnk;	// SPJ** to save space
};

class Exinon {
	int	size;
	int	bias;
	INT53*	int53;
	STYPE**	sig53tab;
	bool	statictab;
	ExinPot* exinpot;
public:
	double	fact;
	EXIN*	data;
	Exinon(const Seq* sd, FTYPE f, const PwdB* pwd);
	~Exinon();
	STYPE	sig53(int m, int n, INTENDS c) const;
	STYPE	sigST(int n, bool init) const {
	    return (STYPE) (data? (init? data[n].sigS: data[n].sigT) / fact : 0);
	}
	EXIN*	score(int n) const {return (data + n);}
	bool	isDonor(int n) const {return (int53[n].cano5);}
	bool	isAccpt(int n) const {return (int53[n].cano3);}
	int	isCanon(int d, int a) const {return
	    (((int53[d].cano5 == 3) && int53[a].cano3 == 3) ||
	    (int53[d].cano5 == 2 && int53[a].cano3 == 2) ||
	    (int53[d].cano5 == 1 && int53[a].cano3) ||
	    (int53[d].cano5 && int53[a].cano3 == 1))?
		int53[d].cano5 + int53[a].cano3: 0;}
	bool	within(const Seq* sd) const;
	int	lplay(int n) const {return (n - bias);}
	int	rplay(int n) const {return (bias + size - n);}
	void	clear() {vset(data + bias, ZeroExin, size + 1);}
	EXIN*	begin() const {return (data + bias);}
	EXIN*	end()	const {return (data + bias + size - 1);}
friend	void Intron53(Seq* sd, const PwdB* pwd, bool both_ori);
friend	void Intron53N(Seq* sd, FTYPE ff, const PwdB* pwd, bool both_ori);
};

inline bool isJunct(int phs5, int phs3) {
	return (phs5 == phs3 && phs5 > -2);
}

class SpJunc {
const	Seq*	b;
	SPJ**	hashent;
	SPJ*	spjunc;
	SPJ*	spp;
	CHAR*	juncseq;
	INT	jsp;
	int	nent;
	int	max_nent;
	CHAR	spliced[6];	// 2 trons  + 4 bases
	CHAR	pyrim[2];
public:
	SpJunc(const Seq* sd);
	~SpJunc();
	CHAR*	spjseq(int n5, int n3);
};

// intron penalty

struct INTRONPEN {
	float	ip, fact, mean;
	int	llmt, mu, rlmt, elmt, tlmt, minl, maxl, mode;
	STYPE*	array, *table;
	float	a1, m1, t1, k1, m2, t2, k2, a2, m3, t3, k3;
};

extern	INTRONPEN IntronPrm;

class IntronPenalty {
	STYPE	GapWI;
	STYPE	AvrSig;
	STYPE	optip;
	FTYPE	IntEp;
	FTYPE	IntFx;
	STYPE*	array;
	STYPE*	table;
public:
	IntronPenalty(VTYPE f, int hh, EijPat* eijpat, ExinPot* exinpot);
	~IntronPenalty() {delete[] array;}
	STYPE	Penalty(int n = -1) const;
	STYPE	Penalty(int n, bool addsig53) const;
	STYPE	PenaltyDrop(int n) const {	// <= 0
	    return (n > IntronPrm.mode? Penalty(n) - optip: 0);
	}
};

/**********
  ip:  intron penalty const. term; fact: factor to length-dep. term; 
  llmt: lower limt; rlmt upper limit to table;
  table: length-dependent term
**********/

extern	void	makeStdSig53();
extern	void	EraStdSig53();
extern	INT	max_intron_len(float p, const char* fn = 0);

/*****************************************************
	Frechet Distribution
*****************************************************/

inline	double frechet_quantile(double p, double mu, double th, double ki) {
	    return (mu + (p > 0.? th / pow(-log(p), 1. / ki): 0));
};

static	const	char	INITIATPAT[]	= "TransInit";
static	const	char	TERMINPAT[]	= "TransTerm";
static	const	char	SPLICE3PAT[]	= "Splice3";
static	const	char	SPLICE5PAT[]	= "Splice5";
static	const	char	BRANCHPAT[]	= "Branch";
static	const	char	INT5PAT[]	= "Intron5";
static	const	char	INT3PAT[]	= "Intron3";
static	const	char	INT53PAT[]	= "Intron53";
static	const	char	INT35PAT[]	= "Intron35";
static	const	char	GNM2TAB[]	= "gnm2tab";
static	const	int	BoundRng = 20;
static	const	int	Ip_equ_k = 3;	// gap length equivalen to Intron Penalty
static	const	char	ipstat[] = "IldModel.txt";

#endif
