/*****************************************************************************
*
*	Header to codepot.c
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

#ifndef  _CODEPOT_H_
#define  _CODEPOT_H_

class	ExinPot;

struct SGPT2 {
	STYPE   sig5;
	STYPE   sig3;
	char   phs5;
	char   phs3;
};

struct SGPT6 {
	STYPE   sig5;
	STYPE   sig3;
	STYPE   sigS;
	STYPE   sigT;
	STYPE   sigE;
	STYPE   sigI;
	char   phs5;
	char   phs3;
};

static	const	SGPT2	ZeroSGPT2 = {0, 0, -2, -2};
static	const	SGPT6	ZeroSGPT6 = {0, 0, 0, 0, 0, 0, -2, -2};
static	const	float	rlmt_quant = 0.8;

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

class Exinon {
	Seq*	sd;
const	PwdB*	pwd;
const	bool	both_ori;
	int	size;
	int	bias;
	INT53*	int53 = 0;
	STYPE**	sig53tab = 0;
	bool	statictab;
	ExinPot* intnpot = 0;
public:
	FTYPE	fact;
	FTYPE	fS;
	STYPE	at_sig5 = 0;
	STYPE	gc_sig5 = 0;
	SGPT2*	data_n = 0;
	SGPT6*	data_p = 0;
	Exinon(Seq* sd_, const PwdB* pwd_, const bool bo);
	~Exinon() {
	    if (data_n) 	delete[] (data_n + bias);
	    if (data_p) 	delete[] (data_p + bias);
	    if (int53)	delete[] (int53 + bias);
	    if (!statictab && sig53tab) {
		delete[] *sig53tab;
		delete[] sig53tab;
	    }
	}
	STYPE	sig53(int m, int n, INTENDS c) const;
	STYPE	sigST(int n, bool init) const {
	    return (STYPE) (data_p? 
		(init? data_p[n].sigS: data_p[n].sigT) / fact : 0);
	}
	SGPT2*	score_n(const int& n) const {return (data_n + n);}
	SGPT6*	score_p(const int& n) const {return (data_p + n);}
	bool	isDonor(const int& n) const {return (int53[n].cano5);}
	bool	isAccpt(const int& n) const {return (int53[n].cano3);}
	int	isCanon(const int& d, const int& a) const {return
	    (((int53[d].cano5 == 3) && int53[a].cano3 == 3) ||
	    (int53[d].cano5 == 2 && int53[a].cano3 == 2) ||
	    (int53[d].cano5 == 1 && int53[a].cano3) ||
	    (int53[d].cano5 && int53[a].cano3 == 1))?
		int53[d].cano5 + int53[a].cano3: 0;}
	int	lplay(const int& n) const {return (n - bias);}
	int	rplay(const int& n) const {return (bias + size - n);}
	SGPT2*	begin_n() const {return (data_n + bias);}
	SGPT2*	end_n()	const {return (data_n + bias + size - 1);}
	bool	good(const SGPT2* bb) const 
		{return (data_n && begin_n() <= bb && bb < end_n());}
	SGPT6*	begin_p() const {return (data_p + bias);}
	SGPT6*	end_p()	const {return (data_p + bias + size - 1);}
	bool	good(const SGPT6* bb) const 
		{return (data_p && begin_p() <= bb && bb < end_p());}
	void	intron53_c();
	void	intron53_n();
	void	intron53_p(const bool dvsp);
	void	resize();
};

static const CHAR spj_tron_tab[257][2] = {
 {14, 14}, {14,  5}, {14, 14}, {14,  5}, { 5, 19}, { 5, 19}, { 5, 19}, { 5, 19}, 
 {14,  4}, {14, 23}, {14,  4}, {14, 23}, { 5, 12}, { 5, 12}, { 5, 15}, { 5, 12}, 
 {19,  8}, {19, 11}, {19,  8}, {19, 11}, {19, 17}, {19, 17}, {19, 17}, {19, 17}, 
 {19,  4}, {19,  4}, {19,  4}, {19,  4}, {19, 13}, {19, 13}, {19, 13}, {19, 13}, 
 { 4,  9}, { 4,  6}, { 4,  9}, { 4,  6}, {23,  3}, {23,  3}, {23,  3}, {23,  3}, 
 { 4, 10}, { 4, 10}, { 4, 10}, { 4, 10}, {23, 22}, {23, 22}, {23, 22}, {23, 22}, 
 {12, 25}, {12, 21}, {12, 25}, {12, 21}, {12, 18}, {12, 18}, {12, 18}, {12, 18}, 
 {15, 24}, {15,  7}, {15, 20}, {15,  7}, {12, 13}, {12, 16}, {12, 13}, {12, 16}, 
 { 8, 14}, { 8,  5}, { 8, 14}, { 8,  5}, {11, 19}, {11, 19}, {11, 19}, {11, 19}, 
 { 8,  4}, { 8, 23}, { 8,  4}, { 8, 23}, {11, 12}, {11, 12}, {11, 15}, {11, 12}, 
 {17,  8}, {17, 11}, {17,  8}, {17, 11}, {17, 17}, {17, 17}, {17, 17}, {17, 17}, 
 {17,  4}, {17,  4}, {17,  4}, {17,  4}, {17, 13}, {17, 13}, {17, 13}, {17, 13}, 
 { 4,  9}, { 4,  6}, { 4,  9}, { 4,  6}, { 4,  3}, { 4,  3}, { 4,  3}, { 4,  3}, 
 { 4, 10}, { 4, 10}, { 4, 10}, { 4, 10}, { 4, 22}, { 4, 22}, { 4, 22}, { 4, 22}, 
 {13, 25}, {13, 21}, {13, 25}, {13, 21}, {13, 18}, {13, 18}, {13, 18}, {13, 18}, 
 {13, 24}, {13,  7}, {13, 20}, {13,  7}, {13, 13}, {13, 16}, {13, 13}, {13, 16}, 
 { 9, 14}, { 9,  5}, { 9, 14}, { 9,  5}, { 6, 19}, { 6, 19}, { 6, 19}, { 6, 19}, 
 { 9,  4}, { 9, 23}, { 9,  4}, { 9, 23}, { 6, 12}, { 6, 12}, { 6, 15}, { 6, 12}, 
 { 3,  8}, { 3, 11}, { 3,  8}, { 3, 11}, { 3, 17}, { 3, 17}, { 3, 17}, { 3, 17}, 
 { 3,  4}, { 3,  4}, { 3,  4}, { 3,  4}, { 3, 13}, { 3, 13}, { 3, 13}, { 3, 13}, 
 {10,  9}, {10,  6}, {10,  9}, {10,  6}, {10,  3}, {10,  3}, {10,  3}, {10,  3}, 
 {10, 10}, {10, 10}, {10, 10}, {10, 10}, {10, 22}, {10, 22}, {10, 22}, {10, 22}, 
 {22, 25}, {22, 21}, {22, 25}, {22, 21}, {22, 18}, {22, 18}, {22, 18}, {22, 18}, 
 {22, 24}, {22,  7}, {22, 20}, {22,  7}, {22, 13}, {22, 16}, {22, 13}, {22, 16}, 
 {25, 14}, {25,  5}, {25, 14}, {25,  5}, {21, 19}, {21, 19}, {21, 19}, {21, 19}, 
 {25,  4}, {25, 23}, {25,  4}, {25, 23}, {21, 12}, {21, 12}, {21, 15}, {21, 12}, 
 {18,  8}, {18, 11}, {18,  8}, {18, 11}, {18, 17}, {18, 17}, {18, 17}, {18, 17}, 
 {18,  4}, {18,  4}, {18,  4}, {18,  4}, {18, 13}, {18, 13}, {18, 13}, {18, 13}, 
 {24,  9}, {24,  6}, {24,  9}, {24,  6}, { 7,  3}, { 7,  3}, { 7,  3}, { 7,  3}, 
 {20, 10}, {20, 10}, {20, 10}, {20, 10}, { 7, 22}, { 7, 22}, { 7, 22}, { 7, 22}, 
 {13, 25}, {13, 21}, {13, 25}, {13, 21}, {16, 18}, {16, 18}, {16, 18}, {16, 18}, 
 {13, 24}, {13,  7}, {13, 20}, {13,  7}, {16, 13}, {16, 16}, {16, 13}, {16, 16}, 
 { 2,  2}
};

static const CHAR spj_amb_tron_tab[64][2] = {
 { 2, 14}, { 2,  5}, { 2, 14}, { 2,  5}, { 2, 19}, { 2, 19}, { 2, 19}, { 2, 19}, 
 { 2,  4}, { 2, 23}, { 2,  4}, { 2, 23}, { 2, 12}, { 2, 12}, { 2, 15}, { 2, 12}, 
 { 2,  8}, { 2, 11}, { 2,  8}, { 2, 11}, { 2, 17}, { 2, 17}, { 2, 17}, { 2, 17}, 
 { 2,  4}, { 2,  4}, { 2,  4}, { 2,  4}, { 2, 13}, { 2, 13}, { 2, 13}, { 2, 13}, 
 { 2,  9}, { 2,  6}, { 2,  9}, { 2,  6}, { 2,  3}, { 2,  3}, { 2,  3}, { 2,  3}, 
 { 2, 10}, { 2, 10}, { 2, 10}, { 2, 10}, { 2, 22}, { 2, 22}, { 2, 22}, { 2, 22}, 
 { 2, 25}, { 2, 21}, { 2, 25}, { 2, 21}, { 2, 18}, { 2, 18}, { 2, 18}, { 2, 18}, 
 { 2, 24}, { 2,  7}, { 2, 20}, { 2,  7}, { 2, 13}, { 2, 16}, { 2, 13}, { 2, 16}
};

static const CHAR spj_tron_amb_tab[64][2] = {
 {14, 2}, { 5, 2}, {14, 2}, { 5, 2}, {19, 2}, {19, 2}, {19, 2}, {19, 2}, 
 { 4, 2}, {23, 2}, { 4, 2}, {23, 2}, {12, 2}, {12, 2}, {15, 2}, {12, 2}, 
 { 8, 2}, {11, 2}, { 8, 2}, {11, 2}, {17, 2}, {17, 2}, {17, 2}, {17, 2}, 
 { 4, 2}, { 4, 2}, { 4, 2}, { 4, 2}, {13, 2}, {13, 2}, {13, 2}, {13, 2}, 
 { 9, 2}, { 6, 2}, { 9, 2}, { 6, 2}, { 3, 2}, { 3, 2}, { 3, 2}, { 3, 2}, 
 {10, 2}, {10, 2}, {10, 2}, {10, 2}, {22, 2}, {22, 2}, {22, 2}, {22, 2}, 
 {25, 2}, {21, 2}, {25, 2}, {21, 2}, {18, 2}, {18, 2}, {18, 2}, {18, 2}, 
 {24, 2}, { 7, 2}, {20, 2}, { 7, 2}, {13, 2}, {16, 2}, {13, 2}, {16, 2}
};

inline bool isJunct(int phs5, int phs3) {
	return (phs5 == phs3 && phs5 > -2);
}

// Frechet Distribution

class SpJunc {
const	Seq*	b;
const	PwdB*	pwd;
public:
	SpJunc(const Seq* sd, const PwdB* pwd_) :
	    b(sd), pwd(pwd_) {}
	~SpJunc() {}
	VTYPE	spjscr(int n5, int n3);
const	CHAR*	spjseq(int n5, int n3);
};

// intron penalty

struct INTRONPEN {
	float	ip, fact, mean;
	int	llmt, mu, rlmt, elmt, tlmt, minl, maxl, mode, nquant;
	STYPE	sip;
	float	a1, m1, t1, k1, m2, t2, k2, a2, m3, t3, k3;
};

extern	INTRONPEN IntronPrm;

struct LenPen {
	SHORT	len;	// intron length
	short	pen;	// penalty
};

class IntronPenalty {
	STYPE	GapWI;
	STYPE	AvrSig;
	STYPE	optip = SHRT_MIN;
	FTYPE	IntEp = 0;
	FTYPE	IntFx = 0;
	STYPE*	array = 0;
	STYPE*	table = 0;
public:
	LenPen*	qm = 0;
	IntronPenalty(VTYPE f, int hh, EijPat* eijpat, ExinPot* exinpot);
	~IntronPenalty() {delete[] array; delete[] qm;}
	double ProbDist(int i, double mu, double th, double kk) {
	    if (i <= mu) return (0.);
	    double	z = th / (i - mu);
	    double	zz = pow(z, kk);
	    return (kk / th * z * zz * exp(-zz));
	}
	STYPE	Penalty() const {return (GapWI);}
	STYPE	Penalty(const int& n) const {
	    if (n < IntronPrm.llmt) return (SHRT_MIN); else
	    if (n < IntronPrm.rlmt && table)
		return (table[n]); else
	    return (STYPE) (IntFx + IntEp * log((double)(n - IntronPrm.mu)));
	}
	STYPE	PenaltyPlus(const int& n) const {
	    if (n < IntronPrm.llmt) return (SHRT_MIN); else
	    if (n < IntronPrm.rlmt && table)
		return (table[n] + AvrSig); else
	    return (STYPE) (IntFx + IntEp * log((double)(n - IntronPrm.mu)) + AvrSig);
	}
	STYPE	PenaltyDrop(const int& n) const {	// <= 0
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
static	const	int	BoundRng = 20;
static	const	int	Ip_equ_k = 3;	// gap length equivalen to Intron Penalty
static	const	char	ipstat[] = "IldModel.txt";

#endif
