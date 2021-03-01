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
