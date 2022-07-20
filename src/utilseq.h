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
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef _UTILSEQ_H_
#define _UTILSEQ_H_ 1

extern	CHAR	gencode[];
class	EiJuncSeq;

enum class Iefp {IF, IP, IB, EF, EP, EB, CF, CP, CB, GF, NG};

static	const	int	itn_lm = 6;	// intron 5' immune margin
static	const	int	itn_rm = 16;	// intron 3' immune margin
static	const	char	CODEPOT[] = "CodePotTab";
static	const	char	INTRONPOT[] = "IntronPotTab";
static	const	char	EXONPOT[] = "ExonPotTab";
static	const	int	CP_NTERM = 4;
static	const	int	STOP = INT_MAX;
static	const	float	maxtonic = 5.;
static	const	CHAR	max_add_size = 24;	// capacity of additional bytes
static	const	char	text_ext[] = ".txt";	// text file
static	const	char	data_ext[] = ".dat";	// binary data file
static	const	char	patmat_ext[] = ".psm";	// extension to binary pssm file
static	const	char	iefp_ext[11][8] = 
	{".ifq", ".ipt", ".ifp", ".efq", ".ept", ".efp", 
	 ".cfq", ".cdp", ".cfp", ".gfq", ""};
static	const	char	iefp_tid[11][16] = 
	{"IntronFrqTab", "IntronPotTab", "IntronFpPTab", 
	 "ExonFrqTab", "ExonPotTab", "ExonFpPTab", 
	 "CodeFrqTab", "CodePotTab", "CodeFpPTab", "GenomeFrqTab", ""};
static	const	char	incompatible[] = "%s is incompatible !\n";

struct	PAT_MATRIX {
	int	rows, cols, offset;
	float*	mtx;
};

struct	ExpectMmm {float min, mean, max;};

struct PatMatHead {
	CHAR	vtype, vsize, nelm, add;	// 0/1:int/float, sizeof, 4/20/23
	int	rows, cols;
};

class PatMat {
public:
	int	rows = 0, cols = 0, offset = 0;
	float	tonic = 0, min_elem = 0; 
	int	transvers = 0, skip = 0;	// used in reading
	int	nsupport = 0, nalpha = 0, morder = 0;
	ExpectMmm	mmm;		// min, mean, max
	float* mtx = 0;

	PatMat(const PatMat& src);	// copy constructor
	PatMat(FILE* fd, bool binary = false);
	PatMat(const char* fname = 0);
	PatMat(const int r, const int c, const int o = 0, float* m = 0);
	~PatMat() {delete[] mtx;}
	void	readPatMat(FILE* fd);
	void	readBinPatMat(FILE* fd);
	float	pwm_score(const Seq* sd, const CHAR* ps, const CHAR* redctab = 0) const;
	float*	calcPatMat(const Seq* sd) const;
	int	order() const {return (morder);}
	int	columns() const {return cols;}
	void	clearmtx() {vclear(mtx, rows * cols);}
	CHAR*	setredctab(const Seq* sd) const;
	void	increment(const Seq* sd, int pos, const CHAR* redctab = 0);
	PatMat&	operator=(const PatMat& src);
};

class ExinPot {
protected:
	int	exin = static_cast<int>(Iefp::IP);
	int	morder = 0;
	int	nphase = 1;	// 1 or 3
	int	ndata = 0;
	int	nsupport = 0;
	int	lm = 0;
	int	rm = 0;
	float	total = 0;	// total # of foreground kmers
	float	avpot = 0;	// average of nsupport self scores
	float	avlen = 0;	// average length of nsupport introns
	float	ess = 0;	// expected mean of self scores
	float*	data = 0;
	float*	read_wdfq(const char* wdfq, const bool& conditional);
	void	count_kmers_1(const Seq* sd, float* fq);
	void	count_kmers_3(const Seq* sd, float* fq);
	void	reform_1(float* bkg = 0);
	void	reform_3();
	float*	calcScr_1(const Seq* sd, float* scr = 0) const;
	float*	calcScr_3(const Seq* sd, float* scr = 0) const;
public:
	ExinPot(const int& ein, const int& mo, const int& nf = 1) :
	    exin(ein), morder(mo), nphase(nf), ndata(ipower(4, mo + 1)), 
	    lm(ein / 3? 0: itn_lm), rm(ein / 3? 0: itn_rm) {
	    if (static_cast<Iefp>(exin) != Iefp::NG)
		readFile(iefp_tid[exin]);
	    else
		fatal("Invalid exin code (%d) !\n", exin);
	}
	ExinPot(const char*& fname) {
	    readFile(fname);
	}
	ExinPot(const int& exn) : exin(exn) {
	    if (static_cast<Iefp>(exin) != Iefp::NG)
		readFile(iefp_tid[exin]);
	    else
		fatal("Invalid exin code (%d) !\n", exin);
	}
	~ExinPot() {delete[] data;}
	bool	ispot() const {
const	    Iefp	iefp = static_cast<Iefp>(exin);
	    return (iefp == Iefp::IP || iefp == Iefp::IB ||
		iefp == Iefp::EP || iefp == Iefp::EB ||
		iefp == Iefp::CP || iefp == Iefp::CB);
	}
	bool	isfpp() const {
const	    Iefp	iefp = static_cast<Iefp>(exin);
	    return (iefp == Iefp::IB || iefp == Iefp::EB || 
		    iefp == Iefp::CB);
	}
	int	size() const {return (ndata);}
	int	dsize() const {
	    return (nphase * ndata + (isfpp()? nphase * ndata: 0));
	}
	bool	readFile(const char* fname);
	float*	getKmers(const char* wdfq, const bool foregrd = true);
	float*	getKmers(EiJuncSeq* eijseq);
	float*	getKmers(int argc, const char** argv);
	void	reform(float* background = 0);
	bool	makeExinPot(const float* gfq);
	bool	readBinary(const char* fname, FILE* fd = 0);
	bool	writeBinary(const char* oname);
	float*	begin() const {return (data);}
	float*	end() const {return (data + nphase * ndata);}
	float*	fbegin() const {
	    return (data + nphase * (isfpp()? ndata: 0));
	}
	float*	fend() const {
	    return (data + dsize());
	}
	float*	calcScr(const Seq* sd, float* scr = 0) const {
	    return ((nphase == 1)? calcScr_1(sd, scr): calcScr_3(sd, scr));
	}
	VTYPE	intpot(const EXIN* b5,const  EXIN* b3) const;
	VTYPE	avrpot(float f = 1.) const {return (VTYPE) (f * avpot);}
	VTYPE	self_score() const {return (ess);}
};

struct EijPat {
	PatMat* pattern5;
	PatMat* pattern3;
	PatMat* patternI;
	PatMat* patternT;
	PatMat*	patternB;
	float   tonic3;
	float   tonic5;
	float   tonicB;
	EijPat(int hh);
	~EijPat();
};

extern	void	setorf(int len, int ic = SILENT);
extern	int	setcodon(int c);
extern	int	codon_id(const CHAR* s, int byte);
extern	void	de_codon_4(CHAR* ncs, int n);
extern	int	toaa(const CHAR* ns);
extern	int	toaa3(const CHAR* ns, int inc);
extern	int	nuc2tron3(const CHAR* ns, int inc);
extern	int	initcodon(int code);
extern	int	initcodon(const char* genspc);
extern	void	mkinvtab();
extern	int	getCodonUsage(const char* fname);
extern	int	setCodonUsage(int gc);
extern	int	fname2exin(const char* fname, int& file_type);

#endif
