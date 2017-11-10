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
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#ifndef  _CODEPOT_H_
#define  _CODEPOT_H_

struct ExinPot;

struct EXIN {
	VTYPE   sig5;
	VTYPE   sig3;
	VTYPE   sigS;
	VTYPE   sigT;
	VTYPE   sigE;
	VTYPE   sigI;
	short   phs5;
	short   phs3;
};

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
	Premat(Seq* seqs[]);
	VTYPE	prematT(const CHAR* ps);
};

struct PwdB;
struct EijPat;

enum INTENDS {IE5, IE3, IE53, IE35, IE5P3};

struct SPJ {
	int	n5;
	int	n3;
	CHAR*	junc;
struct	SPJ*	dlnk;
struct	SPJ**	ulnk;
};

class Exinon {
	int	size;
	int	bias;
	INT53*	int53;
	VTYPE**	sig53tab;
	bool	statictab;
	ExinPot* exinpot;
public:
	double	fact;
	EXIN*	data;
	Exinon(Seq* sd, FTYPE f, PwdB* pwd);
	~Exinon();
	VTYPE	sig53(int m, int n, INTENDS c);
	VTYPE	sigST(int n, bool init)
	{return (VTYPE) (data? (init? data[n].sigS: data[n].sigT) / fact : 0);}
	EXIN*	score(int n) {return (data + n);}
	bool	isDonor(int n) {return (int53[n].cano5);}
	bool	isAccpt(int n) {return (int53[n].cano3);}
	bool	isCanon(int d, int a) {return
	    ((int53[d].cano5 == 3) && int53[a].cano3 == 3) ||
	    (int53[d].cano5 == 2 && int53[a].cano3 == 2) ||
	    (int53[d].cano5 == 1 && int53[a].cano3) ||
	    (int53[d].cano5 && int53[a].cano3 == 1);}
	bool	within(Seq* sd);
friend	void Intron53(Seq* sd, PwdB* pwd);
friend	void Intron53N(Seq* sd, FTYPE ff, PwdB* pwd);
friend	void BoundarySeq(FILE* fd, Seq* sd, RANGE* rng);
friend	void BoundaryScr(FILE* fd, Seq* sd, Gsinfo* gsinf);
};

inline bool isJunct(int ph5, int ph3) {return (ph5 == ph3 && ph5 > -2);}

class SpJunc {
	Seq*	b;
	SPJ**	hashent;
	SPJ*	spjunc;
	SPJ*	spp;
	CHAR*	juncseq;
	CHAR*	jsp;
	int	nent;
	int	max_nent;
	CHAR*	spliced;
	CHAR*	pyrim;
public:
	SpJunc(Seq* sd);
	~SpJunc();
	CHAR*	spjseq(int n5, int n3);
};

// intron penalty

class IntronPenalty {
	VTYPE	GapWI;
	VTYPE	AvrSig;
	FTYPE	IntEp;
	FTYPE	IntFx;
	VTYPE*	array;
	VTYPE*	table;
	int	optlen;
public:
	IntronPenalty(VTYPE f, int hh, EijPat* eijpat, ExinPot* exinpot);
	~IntronPenalty() {delete[] array;}
	VTYPE	Penalty(int n = -1);
	VTYPE	Penalty(int n, bool addsig53);
	int	mode() {return (optlen);}
	VTYPE	maxpenalty() {return Penalty(optlen);}
};

struct INTRONPEN {
	float	ip, fact, mean;
	int	llmt, mu, rlmt, elmt, tlmt, maxl;
	VTYPE*	array, *table;
	float	a1, m1, t1, k1, m2, t2, k2, a2, m3, t3, k3;
};

/**********
  ip:  intron penalty const. term; fact: factor to length-dep. term; 
  llmt: lower limt; rlmt upper limit to table;
  table: length-dependent term
**********/

extern	void	makeStdSig53();
extern	void	EraStdSig53();

extern	INTRONPEN IntronPrm;

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

#endif
