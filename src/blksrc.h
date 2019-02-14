/*****************************************************************************
*
*	Quickly search for similar blocks
*	Using oligomer compositions
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

#ifndef	_BLKSRC_H_
#define	_BLKSRC_H_

#include <math.h>
#include "bitpat.h"

#define	TESTRAN	0

static	const	int	BsVersion = 25;
static	const	float	RbsFactLog = 0.4;
static	const	float	RbsFactSqr = 0.303 * 2;
static	const	float	RlnFact = 0.00520 * 2;
static	const	INT	NRTAB = 128;
static	const	int	ENTSIZE = 23;	// for back compatibility
static	const	char	ipstat[] = "IldModel.txt";
static  const   CHAR    ntconv[26] =
	{0,4,1,4,7,7,2,4,7,7,4,7,4,4,7,4,7,4,4,3,3,4,4,7,4,7};
//	 a b c d e f g h i j k l m n o p q r s t u v w x y z

union	BYTE4	{INT i; CHAR c[4];};
union	BYTE2	{SHORT s; CHAR c[2];};
typedef	INT	BLKTYPE;

struct Block2Chr {
	double	BClw;	// lowest off-diagonal
	double	BCup;	// highest off-diagonal
	double	BCce;	// Block->ChrId coefficient
};

struct BPAIR {
	int	bscr, chr;
	VTYPE	jscr;
	SHORT	d, e;
	BLKTYPE	lb, rb;
};

struct BlkScr {
	BLKTYPE key;
	int	bscr;
	bool	operator<(const BlkScr& b) const {return (bscr < b.bscr);}
	bool	operator==(const BlkScr& b) const {return (key == b.key);}
};

class Bhit4 {
	BlkScr*	blkscr;
public:
	PrQueue<BlkScr>*	prqueue_a[4];
	PrQueue<BlkScr>*	prqueue_b[4];
	INT	sigm;
	BPAIR*	bpair;
	CHAR**	as[4];
	BLKTYPE*	sigb[4];	// top blocks with high bscr values
	BLKTYPE*	sigw;		// working 
	SHORT	sign[4];
	SHORT	maxs[4];
	SHORT	nhit[4];
	SHORT	mmct[4];
	int*	ascr[4];	// count all word hits
	int*	bscr[4];	// direct / reverse x left / right
	CHAR*	bend;		// indicator containing gene end
	INT	testword[4];
	void	update_a(int d, BLKTYPE blk);
	int	update_b(int d, BLKTYPE blk);
	SHORT	extract_to_work(int d);
	Bhit4(int nseg);
	~Bhit4();
};

class Bhit2 {
	BlkScr*	blkscr;
public:
	PrQueue<BlkScr>*	prqueue_b;
	INT	sigm;
	CHAR**	as[2];
	int*	rscr;
	int*	hpos;
	Bhit2(int nseg);
	~Bhit2();
	INT	hsort();
	int	update(BLKTYPE blk);
};

struct SeqList {
	Seq**	seqs;
	Seq**	lstgr;
	Seq**	curgr;
};

class Randbs {
	float	RbsCoef;
	float	RlnCoef;
	float	RbsCons;
	double	(*trans_form)(double x);
	int	rscrtab[NRTAB];
public:
	int	Phase1T;
	Randbs(double avr, bool gdb);
	~Randbs() {};
	int	randbs(INT mmc);
	int	lencor(int len) {return (int) (RlnCoef * sqrt(double(len)));}
	int	base() {return rscrtab[0];}
};

struct GeneRng {
	int	sid;
	int	lend;
	int	rend;
	int	lgap;
	int	rgap;
	int	len;
	int	num;
	VTYPE	scr;
	JUXT*	jxt;
};

struct CHROMO {
	long	fpos;
	INT	spos;
	INT	segn;
};

struct CHROMO21 {
	long	fpos;
	INT	spos;
	INT	segn;
	char	entry[ENTSIZE + 1];
};

struct ContBlk22 {
	INT	ConvTS;	// Convertor table size
	INT	WordNo;	// total number of proper blocks
	INT	WordSz;	// block data size in 2 * bytes
	INT	ChrNo;	// # of chromosomes/contigs
	INT	glen;	// genome length
	SHORT	AvrScr;	// average word score
	SHORT	MaxBlk;	// max # of hits in a block
	SHORT	BytBlk;	// bytes to present a block
	SHORT	VerNo;	// version no
	SHORT*	Nblk;	// # of blocks containing each kmer
	INT**	blkp;	// list of blocks containing each kmer
	INT*	blkb;	// buffer of blkp
	short*	wscr;	// score of the kmer
	CHROMO*	ChrID;	// Info of chromosomes/contigs
};

struct ContBlk {
	INT	ConvTS;	// Convertor table size
	size_t	WordNo;	// total number of proper blocks
	size_t	WordSz;	// block data size in 2 * bytes
	size_t	ChrNo;	// # of chromosomes/contigs
	size_t	glen;	// genome length
	SHORT	AvrScr;	// average word score
	SHORT	MaxBlk;	// max # of hits in a block
	SHORT	BytBlk;	// bytes to present a block
	SHORT	VerNo;	// version no
	SHORT*	Nblk;	// # of blocks containing each kmer
	INT**	blkp;	// list of blocks containing each kmer
	INT*	blkb;	// buffer of blkp
	short*	wscr;	// score of the kmer
	CHROMO*	ChrID;	// Info of chromosomes/contigs
	void	clean() {Nblk = 0; blkp = 0; blkb = 0; wscr = 0; ChrID = 0;}
	ContBlk();
	~ContBlk();
};

class Chash : public Dhash<INT, int> {
	ContBlk*	pwc;
public:
	Chash(int n, ContBlk* pwc_) : Dhash<INT, int>(n), pwc(pwc_) {}
	~Chash() {}
	void	countBlk();
	void	registBlk(INT m);
};

class WordTab : public ReducWord {
protected:
	INT	Nbitpat;
	INT	BitPat;
	INT	BitPat2;
	INT	Nshift;
	INT	Nss;
	size_t	total;
	Chash*	ch;
	Bitpat_wq**	bpp;
	SHORT*	ss;
	SHORT**	ss6;
	INT*	ww;
	SHORT	cc[6];
	INT	xx[3];
	INT**	idx;
	INT*	count;
	INT*	tcount;
	void	reset();
	void	MakeAaWc(const char* argv[]);
	void	MakeNucWc(const char* argv[]);
	void	MakeTronWc(const char* argv[]);
public:
	WordTab(Seq* sd, INT kmer, INT nsft = 1, INT elms = 0, const char* ap = 0, 
	    INT bp = 0, INT bp2 = 0, INT nbt = 1);
	~WordTab();
	void	maktab(const char* argv[]);
	bool	c2w(int c, int mode = 0, INT i = UINT_MAX);
	bool	c2w6(int c, int& p, int mode = 0, INT i = UINT_MAX);
	void	countwd(const char* argv[]);
	int	findword(INT* locs, INT word, INT k, INT from, INT to);
};

class MakeBlk : WordTab {
friend	class	SrchBlk;
friend	class	AdjacentMat;
	ContBlk*	pwc;
	Block2Chr*	pbc;
	double	deltaa;
	double*	acomp;
	void	prepacomp();
	void	blkscrtab(size_t segn);
	void	blkscrtab(size_t segn, INT blksz);
	INT	chrblk(int m) {
	    return (pwc->ChrID? pwc->ChrID[m].segn: m + 1);
	}
	DbsDt*	wdbf;
	bool	isaa;
	bool	istron;
	SEQ_CODE*	defcode;
public:
	MakeBlk(Seq* sd, DbsDt* dd = 0);
	~MakeBlk() {delete pwc; delete pbc; delete[] acomp;}
	void	findChrBbound();
	void	idxblk(int argc, const char** argv, int molc);
	void	idxblk(Seq* sd, SeqServer* svr);
	void	idxblk(Seq* sd);
	void	WriteBlkInfo();
	int	dbsmolc() {return (molc);}
	int	no_entry() {return pwc->ChrNo;}
	void	delete_dbf() {delete wdbf; wdbf = 0;}
template <typename file_t>
	void	first_phase(file_t fd, int& entlen, CHROMO& chrbuf);
template <typename file_t>
	int	second_phase(file_t fd, CHROMO& chrbuf, int m,
		DbsRec*& pr, CHAR*& ps, char*& pe,
		DbsRec& prvrec, CHROMO*& pchrid);
template <typename file_t>
	void	writeBlkInfo(file_t fd, const char* fn);
};

#if TESTRAN
struct Testran {
	INT	trnbr;
	INT	trcnt[TESTRAN];
	INT	trmax[TESTRAN];
	INT	trmin[TESTRAN];
	INT	trwrd[TESTRAN];
	double	travr[TESTRAN];
	double	trvar[TESTRAN];
	Testran();
	void	out(int avrscr);
	void	add(const Testran* tr);
}
#endif

class SrchBlk {
	int	kk;
	int	DRNA;
	bool	gnmdb;
	VTYPE	vthr;
	int	DeltaPhase2;
	SHORT	ptpl;
	int	MaxBlock;
	INT	ExtBlock;
	INT	maxmmc;
	INT*	ConvTab;
	ContBlk*	pwc;
	Block2Chr*	pbc;
	INT*	SegLen;
	INT	MinGeneLen;
	short	poslmt;
	Bhit2*	bh2;	// thread specific
	Bhit4*	bh4;
	SrchBlk*	master;
	void	initialize(Seq* seqs[], const char* fn = 0);
	void	init2(Seq* a);
	void	init4(Seq* a);
	VTYPE	FindHsp(BPAIR* wrkbp, VTYPE maxjscr, SeqList* sl);
	int	TestOutput(Seq* seqs[], int force);
	INT	bestref(Seq* seqs[], KVpair<INT, int>* sh, int n);
	Seq*	setgnmrng(BPAIR* wrkbp, SeqList* sl);
	void	setaaseq(Seq* a, int chn);
	INT	chrblk(int m) {
	    return (pwc->ChrID? pwc->ChrID[m].segn: m + 1);
	}
	int	chrsize(int m) {
	    return (pwc->ChrID? pwc->ChrID[m + 1].spos - pwc->ChrID[m].spos:
		 dbf->recidx[m].seqlen);
	}
#if TESTRAN
	Testran*	tstrn;
#endif
public:
	PwdB*	pwd;
	int	bbt;
	INT	nseg;
	DbsDt*	dbf;
	Randbs*	rdbt;
	Bitpat** bpp;
	void	reset(DbsDt* df) {
		dbf = df;
#if TESTRAN
		tstrn = new Testran();
#endif
	}
	SrchBlk(Seq* seqs[], const char* fn, bool gdb = false);
	SrchBlk(Seq* seqs[], MakeBlk* mb, bool gdb = false);
	SrchBlk(SrchBlk* sbk, DbsDt* df);
	~SrchBlk();
	void	setSegLen();
	int	MinQuery();
	int	MaxGene();
template <typename file_t>
	int	read_blk_dt(file_t fd);
template <typename file_t>
	void	ReadBlkInfo(file_t fd, const char* fname);
	int	findChrNo(INT n);
	int	findblock(Seq* seqs[]);
	int	findh(Seq* seqs[]);
	int	finds(Seq* seqs[]);
	bool	grngoverlap(GeneRng* a, GeneRng* b);
	bool	incompatible(Seq* a) {return (DRNA ^ a->isdrna());}
};

class Qwords {
	int	kk;
	int	DRNA;
	INT*	ww;
	int*	xx;
	CHAR**	endss;
	BLKTYPE* front;
	INT*	ConvTab;
	ContBlk*	wc;
	Bitpat** bpp;
public:
	void	init_mrglist();
	BLKTYPE	next_mrglist();
	int	querywords(CHAR* ss);
	int	querywords(CHAR* ss, int d, bool rvs);
	void	reset(Seq* a);
	Qwords(int k, int nc, INT* ct, ContBlk* pwc, Bitpat** bp, Seq *a = 0);
	~Qwords() {delete[] ww; delete[] xx; delete[] endss; delete[] front;}
};

extern	MakeBlk*	makeblock(int argc, const char** argv, int molc = UNKNOWN);
extern	MakeBlk*	makeblock(SeqServer* svr);
extern	MakeBlk*	makeblock(Seq* sd);
extern	int	setQ4prm(const char* s, const char* ss = 0);
extern	void	ReportBlkInfo(const char* fn);
extern	int	genemergin(int agap, int mingap, Seq* sd, bool rend);
extern	bool	extend_gene_rng(Seq* sqs[], PwdB* pwd, DbsDt* dbf);
extern	void	set_max_extend_gene_rng(int n, bool forced = false);
extern	int	get_max_extend_gene_rng();

#endif

