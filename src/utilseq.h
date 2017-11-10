/*****************************************************************************
*
*	Micelleious utilities subroutines
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

#ifndef _UTILSEQ_H_
#define _UTILSEQ_H_ 1

#define	CODEPOT		"CodePotTab"
#define	INTRONPOT	"IntronPotTab"
#define	EXONPOT		"ExonPotTab"

extern	CHAR	gencode[];

static	const	int	CP_NTERM = 4;
static	const	int	STOP = INT_MAX;
static	const	float	maxtonic = 5.;

struct	MATRIX	{int rows, cols, offset; int maxidx, minidx; float tonic; double* mtx;};
struct	ExpectMmm {float min, mean, max;};

class PatMat {
	void	readPatMat(FILE* fd);
protected:
	int	rows, cols, offset;
	int	maxidx, minidx, nsupport;
public:
	float	tonic;
	ExpectMmm	mmm;	// min, mean, max
	double* mtx;
	PatMat(FILE* fd);
	PatMat(const char* fname = 0);
	PatMat(const int r, const int c, const int o = 0, float* m = 0);
	~PatMat() {delete[] mtx;}
	float*	calcPatMat(Seq* sd);
};

class CodePot {
	FTYPE*	CodePotTab[CP_NTERM];
	FTYPE*	CodePotBuf;
	int	CodePotType;
	int	CodePotClmn;
	float*	calc5MMCodePot(Seq* sd, int phase);
	float*	calcDitCodePot(Seq* sd, int phase);
public:
	CodePot();
	~CodePot() {delete[] CodePotBuf;}
	float*	calcPrefCodePot(Seq* sd, int phase);
};

class ExinPot {
	int	morder;
	int	size;
	int	lm;
	int	rm;
	float	avpot;
	float	avlen;
	FTYPE*	ExonPot;
	FTYPE*	IntronPot;
	bool	readExinPot(const char* fname);
public:
	ExinPot(int zZ, const char* fname = 0);
	~ExinPot() {delete[] ExonPot; delete[] IntronPot;}
	float*	calcExinPot(Seq* sd, bool exon);
	VTYPE	intpot(EXIN* b5, EXIN* b3);
	VTYPE	avrpot(float f) {return (VTYPE) (f * avpot);}
};

struct EijPat {
	PatMat* pattern5;
	PatMat* pattern3;
	PatMat* patternI;
	PatMat* patternT;
	PatMat*	patternB;
	FTYPE   tonic3;
	FTYPE   tonic5;
	FTYPE   tonicB;
	EijPat(int hh);
	~EijPat();
};

struct SpbInfo {
	CodePot*	codepot;
	EijPat*	eijpat;
	SpbInfo();
	~SpbInfo() {delete codepot; delete eijpat;}
};

extern	void	setorf(int len, int ic);
extern	int	setcodon(int c);
extern	int	codon_id(CHAR* s, int byte);
extern	void	de_codon_4(CHAR* ncs, int n);
extern	int	toaa(CHAR* ns);
extern	int	toaa3(CHAR* ns, int inc);
extern	int	nuc2tron3(CHAR* ns, int inc);
extern	int	initcodon(int code);
extern	int	initcodon(const char* genspc);
extern	void	mkinvtab();
extern	int	getCodonUsage(char* fname);
extern	int	setCodonUsage(int gc);

#endif
