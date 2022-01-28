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
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include "bitpat.h"
#include "wln.h"

#define	TESTRAN	0

class Block;
class MakeBlk;

typedef	std::vector<INT>::iterator IntItr;
struct	SeqLenStat {int num; INT maxv; LONG total;};

extern	CHAR	gencode[];
extern	int	thread_num;
extern	bool	ignoreamb;

static	const	int	BsVersion = 26;
static	const	int	ENTSIZE = 23;	// for back compatibility

/***********************************************************
	multi-threads
************************************************************/

#if M_THREAD

#include <pthread.h>
#include <sched.h>

struct BlkQueue {
int	qid;
	int	rp, wp, remain;
	int	qsize;
	pthread_cond_t	not_full;
	pthread_cond_t	not_empty;
	pthread_cond_t	not_busy;
	Block**	blkque;
	pthread_mutex_t	mutex;
	BlkQueue(int qs);
	~BlkQueue() {delete[] blkque;}
	void	enqueue(Block* blk);
	Block*	dequeue();
	bool	is_empty() const {return (remain == 0);}
	void	reset()	{rp = wp = remain = 0;}
};

struct Targ {
	int		tid;
	INT*		tc;
	BlkQueue*	prev_q;
	BlkQueue*	next_q;
};

#endif

/*****************************************************
	MakeDbs
*****************************************************/

class MakeDbs {
	DbsRec	rec;
	size_t	recnbr;
	bool	cridxf;
const	bool	isaa;
	int	b;
	INT	maxlen;
	char	entrystr[MAXL];
	char	entryprv[MAXL];
	FILE*	fgrp;
	FILE*	fseq;
	FILE*	fidx;
	FILE*	fent;
#if USE_ZLIB
	gzFile	gzseq;
#endif
	char*	dbname;
	void	putsq(int c) {
		if (fseq)	fputc(c, fseq);
#if USE_ZLIB
		else	 fputc(c, gzseq);
#endif
	}
public:
	MakeDbs(const char* dbname, int molc);
	~MakeDbs() {
		fclose(fgrp);
		if (fseq) fclose(fseq);
#if USE_ZLIB
		if (gzseq) fclose(gzseq);
#endif
		fclose(fidx);
		fclose(fent);
		delete[] dbname;
	}
	INT	max_dbseq_len() {return maxlen;}
	void	mkidx();
template <typename file_t>
	int	write_recrd(file_t fd, int c = 0);
	void	putseq(int c) {
		if (isaa) putsq(c);
		else if (rec.seqlen & 1) putsq(c + b);
		else {
			b = c << 4;
			if (c == SEQ_DELIM) putsq(c + b);	// double delimiter
		}
		if (c != SEQ_DELIM) ++rec.seqlen;
	}
	void	wrtgrp(const char* ps) {
		long	 fpos = 0L;
		if (fseq)	fpos = ftell(fseq);
#if USE_ZLIB
		else	 fpos = ftell(gzseq);
#endif
		fprintf(fgrp, "%8ld %u %s\n", fpos, (INT) recnbr, ps);
	}
	void stamp21() {
		DbsRec	rec21 = {magicver21, false, 0};
		fwrite(&rec21, sizeof(DbsRec), 1, fidx);	// header record
	}

};

/*****************************************************
	Pre Scan
*****************************************************/

class PreScan : ReducWord {
template <typename file_t>
	void	scan_genome(file_t fd, SeqLenStat& as);
public:
	PreScan(const Seq* sd, INT elem) : ReducWord(sd, elem) {}
	void	lenStat(int argc, const char** av, SeqLenStat& as);
};

/*****************************************************
	Block
*****************************************************/

struct BlkWcPrm {
	INT	Nalpha;
	INT	Ktuple;
	INT	Bitpat2;
	INT	TabSize;
	INT	BitPat;
	INT	Nshift;
	INT	blklen;
	INT	MaxGene;
	short	Nbitpat;
	SHORT	afact;	  // abundant / even
};

struct Block2Chr {
	double	BClw;	// lowest off-diagonal
	double	BCup;	// highest off-diagonal
	double	BCce;	// Block->ChrId coefficient
};

struct CHROMO {
	INT	spos;
	INT	segn;
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
public:
	Chash(int n) 
	  : Dhash<INT, int>(n)  {}
	~Chash() {}
	void	countBlk(ContBlk& cntblk);
	void	registBlk(INT m, ContBlk& cntblk);
};

class Block : public WordTab {
protected:
const	Seq*	sd;
const	bool	isaa;
const	bool	istron;
const	bool	mkblk;
	Chash*	ch;
public:
	void	c2w(INT uc, INT* tcc = 0);
	void	c2w6(INT uc, INT* tcc = 0);
	void	c2w6_pp(int n);
	Block(const Seq* sq, INT blklen, bool mk_blk = true);
	~Block() {delete ch;}
#if M_THREAD
protected:
	int	cps;	// chromosomal position
	INT	prelude;
	INT	margin;
	char*	seq;
	char*	esq;
	char*	tsq;
public:
	int	bid;	// block ID
	void	setbid(int b) {bid = b;}
	int	getbid() const {return(bid);}
	int	cpos() const {return (cps);}
	int	cpos(int p) {return (cps = p);}
	char*	begin(char* s = 0) {if (s) seq = s; return seq;}
	char*	right() const {return (tsq == esq? esq - margin: tsq);}
	char*	end() const {return (tsq);}
	void	end(char* t) {tsq = t;}
	void	restore() {tsq = esq;}
	Block(Block& src, char* s);
friend	class	MakeBlk;
friend	void*	worker(void* arg);
#endif	// M_THREAD
};

/*****************************************************
	MakeBlk
*****************************************************/

class MakeBlk : public Block {
const	Seq*	sd;
	MakeDbs*	mkdbs;
	DbsDt*	wdbf;
	SEQ_CODE*	defcode;
	INT*	tcount;
	ContBlk	cntblk;
	Block2Chr	b2c;
	CHROMO	chrbuf;
	CHROMO*	pchrid;
	CHAR*   s2r;    // seq to reduced word
	int	bias;
const	int	s_size;
	double	deltaa;
	double* acomp;
	void	prepacomp();
	void	blkscrtab(size_t segn);
	void	blkscrtab(size_t segn, INT blksz);
	INT	chrblk(int m) {
	    return (cntblk.ChrID? cntblk.ChrID[m].segn: m + 1);
	}
	int	encode(int c) {
	    return (defcode->encode[toupper(c) - 'A'] - bias);
	}
friend	class	Chash;
friend	class	SrchBlk;
friend	class	AdjacentMat;
#if M_THREAD
	int	did;
	int	c_qsize;	// circular queue size
	BlkQueue**	bq;
	Block**	blks;
	char*	seqbuf;
	void	harvest(Block* blk, bool first);
friend	void*	worker(void* arg);
#endif	// M_THREAD

public:
	MakeBlk(const Seq* sd, DbsDt* dd = 0, MakeDbs* mkdbs = 0, bool mk_blk = true);
	~MakeBlk();
	void	findChrBbound(Block2Chr* pb2c);
	void	idxblk(int argc, const char** argv);
	void	idxblk(Seq* sd, SeqServer* svr);
	void	idxblk(Seq* sd);
	void	WriteBlkInfo();
	int	dbsmolc() const {return (molc);}
	int	no_entry() const {return cntblk.ChrNo;}
	DbsDt*	get_dbf() const {return wdbf;}
	void	delete_dbf() {delete wdbf; wdbf = 0;}
	void	store_blk(bool first, Block* blk = 0);
template <typename file_t>
	void	scan_genome(file_t fd, INT* tc = 0);
template <typename file_t>
	void	writeBlkInfo(file_t fd, const char* fn);

#if M_THREAD
template <typename file_t>
	void	m_scan_genome(file_t fd, bool first);
	void	m_idxblk(int argc, const char** argv);
#endif	// M_THREAD
};

/*****************************************************
	For back compatibility
*****************************************************/

struct CHROMO21 {	
	long	fpos;
	INT	spos;
	INT	segn;
	char	entry[ENTSIZE + 1];
};

struct CHROMO25 {
	long	fpos;
	INT	spos;
	INT	segn;
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

/*****************************************************
	SrchBlk
*****************************************************/

static	const	float	RbsFactLog = 0.4;
static	const	float	RbsFactSqr = 0.303 * 2;
static	const	float	RlnFact = 0.00520 * 2;
static	const	INT	NRTAB = 128;
union	BYTE4	{INT i; CHAR c[4];};
union	BYTE2	{SHORT s; CHAR c[2];};
typedef	INT	BLKTYPE;

struct BPAIR {
	int	bscr, chr;
	VTYPE	jscr;
	BLKTYPE	lb, rb;		// left-most & right-most blocks
	BLKTYPE	ub, db;		// nearest upstream & downstream bolck pairs
	BLKTYPE zl, zr;		// chromosomal boundaries
	SHORT	rvs;		// reverse strnad
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
	PrQueue_wh<BlkScr>*	prqueue_a[4];
	PrQueue_wh<BlkScr>*	prqueue_b[4];
	INT	sigm;
	BPAIR*	bpair;
const	CHAR**	as[4];
	BLKTYPE*	sigw[2];	// working 
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
	Bhit4(int nseg);
	~Bhit4();
};

class Bhit2 {
	BlkScr*	blkscr;
public:
	PrQueue<BlkScr>*	prqueue_b;
	INT	sigm;
const	CHAR**	as[2];
	int*	rscr;
	int*	hpos;
	Bhit2(int nseg);
	~Bhit2();
	INT	hsort();
	int	update(BLKTYPE blk);
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
	int	lencor(int len) const {return (int) (RlnCoef * sqrt(double(len)));}
	int	base() const {return rscrtab[0];}
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
	Seq**	seqs;
	Seq*	query;
	Seq**	gener;
	Seq**	lstgr;
	Seq**	curgr;
	int	kk;
	int	DRNA;
	bool	gnmdb;
	VTYPE	vthr;	// each gene
	VTYPE	critjscr;
	int	DeltaPhase2;
	SHORT	ptpl;
	INT	maxmmc;
	CHAR*	ConvTab = 0;
	ContBlk*	pbwc = 0;
	Block2Chr*	pb2c = 0;
	INT*	SegLen = 0;
	INT	MinGeneLen;
	int	min_agap;	// codons
	Bhit2*	bh2 = 0;	// thread specific
	Bhit4*	bh4 = 0;
	SrchBlk*	master = 0;
	void	initialize(Seq** sqs, const char* fn = "");
	void	init2(const Seq* sd);
	void	init4(const Seq* sd);
	int	FindHsp(BPAIR* wrkbp);
	int	TestOutput(int force);
	INT	bestref(Seq* sqs[], KVpair<INT, int>* sh, int n);
	Seq*	setgnmrng(const BPAIR* wrkbp);
	void	setaaseq(Seq* sd, int chn);
	SHORT	extract_to_work(int d);
	INT	chrblk(int m) {
	    return (pbwc->ChrID? pbwc->ChrID[m].segn: m + 1);
	}
	int	chrsize(int m) {
	    return (pbwc->ChrID? pbwc->ChrID[m + 1].spos - pbwc->ChrID[m].spos:
		dbf->recidx[m].seqlen);
	}
#if TESTRAN
	Testran*	tstrn = 0;
#endif
public:
	PwdB*	pwd = 0;
	int	bbt;
	INT	nseg;
	int	NoWorkSeq;
	DbsDt*	dbf;
	Randbs*	rdbt = 0;
	Bitpat** bpp = 0;
	void	reset(DbsDt* df) {
		dbf = df;
#if TESTRAN
		tstrn = new Testran();
#endif
	}
	void	setseqs(Seq** sqs) {
	    seqs = sqs;
	    query = seqs[0];
	    gener = seqs + 1 + NoWorkSeq;
	    lstgr = gener + OutPrm.MaxOut2;
	}
	Seq**	get_gener() const {return gener;}
	SrchBlk(Seq* sqs[], const char* fn, bool gdb = false);
	SrchBlk(Seq* sqs[], MakeBlk* mb, bool gdb = false);
	SrchBlk(SrchBlk* sbk, DbsDt* df);
	~SrchBlk();
	void	setSegLen();
	int	MinQuery() const;
	int	MaxGene() const;
template <typename file_t>
	int	read_blk_dt(file_t fd);
template <typename file_t>
	void	ReadBlkInfo(file_t fd, const char* fname);
	int	findChrNo(const INT& n);
	int	findblock(Seq** sqs);
	int	findh(Seq** sqs);
	int	finds(Seq** sqs);
	bool	grngoverlap(const GeneRng* a, const GeneRng* b);
	bool	incompatible(Seq* a) {return (DRNA ^ a->isdrna());}
};

class Qwords {
	int	kk;
	int	DRNA;
	INT*	ww;
	int*	xx;
	CHAR**	endss;
	BLKTYPE* front;
	CHAR*	ConvTab;
	ContBlk*	wc;
	Bitpat** bpp;
public:
	void	init_mrglist();
	BLKTYPE	next_mrglist();
	int	querywords(const CHAR* ss);
	int	querywords(const CHAR* ss, int d, bool rvs);
	void	reset(Seq* a);
	Qwords(int k, int nc, CHAR* ct, ContBlk* pwc, Bitpat** bp, Seq *a = 0);
	~Qwords() {delete[] ww; delete[] xx; delete[] endss; delete[] front;}
};

extern	MakeBlk*	makeblock(int argc, const char** argv, int molc, bool mkblk = true);
extern	MakeBlk*	makeblock(SeqServer* svr);
extern	MakeBlk*	makeblock(Seq* sd);
extern	int	setQ4prm(const char* ps, const char* ss = 0);
extern	void	ReportBlkInfo(const char* fn);
extern	int	genemergin(int agap, int mingap, Seq* sd, bool rend);
extern	bool	extend_gene_rng(Seq* sqs[], const PwdB* pwd, DbsDt* dbf);
extern	void	set_max_extend_gene_rng(int n, bool forced = false);
extern	int	get_max_extend_gene_rng();

#endif

