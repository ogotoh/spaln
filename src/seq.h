/*****************************************************************************
*
*	Header for sequence analysis
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

#ifndef  _BSEQ_H_
#define  _BSEQ_H_

struct	ALPRM {float u, v, u0, u1, v0, tgapf, thr, scale, maxsp, gamma; int k1, ls, sh, mtx_no;};
struct	ALPRM2 {float x, y, z, o, bti, spb, Z, sss; int jneibr, termk1;};
struct	ALPRM3 {float scnd, hydr, hpmt; int hpwing, no_angle;};

extern	ALPRM2	alprm2;
extern	ALPRM3	alprm3;

#include "cmn.h"
#include "dbs.h"
#include "gsinfo.h"
#include "gaps.h"
#include "codepot.h"
#include "utilseq.h"

#define	USE_ETHER 	FVAL

#ifndef _COMM
#define	_COMM	';'
#endif
static	const	char	_NHEAD = '>';		/* Normal Strand	*/
static	const	char	_CHEAD = '<';		/* Compl. Strand	*/
static	const	char	_APPN = '+';		/* Append on 		*/
static	const	char	_UNP = '-';		/* Gap			*/
static	const	char	_TRM = '*';		/* Termination Code	*/
static	const	char	_EOS = '/';		/* End of Seq.		*/
static	const	char	_LCOM = '#';		/* Comment of list	*/
static	const	char	_SAME = '~';		/* Same as 1st		*/
static	const	char	_IBID = '^';		/* Same as above	*/
static	const	char	_PEPT = '@';		/* Pept. Coding 	*/
static	const	char	E_DLM = '.';		/* Exon Number		*/
static	const	char	I_DLM = '/';		/* Intron Number	*/
static	const	char	_PROF = '%';		/* Profile expres.	*/
static	const	char	_DELG = '#';		/* Eliminate del sites	*/
static	const	char	_LABL = '|';		/* Label to seq.	*/
static	const	char	_WGHT = '%';		/* Weight	       */
static	const	char	GBKID = '$';		/* GenBank Seq.		*/

static	const	int	DEFSEQLEN = 1024;
static	const	int	MAXCR = 128;
static	const	int	BUFLEN = MAXL;
static	const	int	CPY_SEQ = 1;	/* Copy seq		*/
static	const	int	CPY_NBR	= 2;	/* Copy nbr		*/
static	const	int	CPY_LBL = 4;	/* Copy lbl		*/
static	const	int	CPY_SBI = 8;	/* copy sp. bound. info	*/
static	const	int	CPY_ALL	= 15;	/* Copy all		*/
static	const	int	RLT_SBI	= 16;	/* renumber lst		*/

static	const	int	nil_code = 0;
static	const	int	gap_code = 1;

static	const	int	USE_EXG = 1;
static	const	int	NSIMD = NTS+USE_EXG;
static	const	int	ASIMD = AAS+USE_EXG;
static	const	int	CSIMD = 4+2+USE_EXG;
static	const	int	TSIMD = 23+2+USE_EXG;
static	const	int	RSIMD = 20+2+USE_EXG;
static	const	int	NAMB = A + 4;
static	const	int	NUCDIM = 4+2+USE_EXG;

static	const	CHAR	t2atab[] = 
{NIL,UNP,AMB,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,
PHE,PRO,SER,THR,TRP,TYR,VAL,SER,AMB,AMB};

enum {DEL_GAP, DEL_AMB, DEL_GRC, DEL_ARC, RET_GAP};

enum {	LABEL_None,		/* No r-field label	*/
	LABEL_Once,		/* Unsed for phylip	*/
	LABEL_Numb,		/* Number r-field label	*/
	LABEL_Name		/* Name r-field label	*/
};

enum Row_Mode {
	Row_None, Row_Last, Row_Every, Row_First, Row_Ditto, Form_Native = 4,
	Form_Bare, Form_Phylp, Form_GCG, Form_CLW, Form_Serial, Form_GDE
};

struct	SEQ_CODE {
	int	max_code;	// # of alphabets
	int	amb_code;	// ambigity code
	int	base_code;	// elementary code
	int	ceil_code;	// upper bound of
	int	gap_prof;	// deletion profile
	CHAR*	encode;
	char*	decode;
	CHAR*	redctab;
};

struct INEX {
	INT	molc:	3;
	INT	sens:	2;
	INT	form:	1;
	INT	vect:	2;
	INT	cmpc:	1;
	INT	dels:	1;
	INT	ambs:	2;
	INT	gfrq:	1;
	INT	prof:	1;
	INT	sngl:	1;
	INT	trcv:	1;
	INT	algn:	2;
	INT	exgl:	1;
	INT	exgr:	1;
	INT	intr:	1;
	INT	dela:	1;
	INT	polA:	2;
	INT	est:	1;
	INT	ori:	2;
	INT	vtl:	1;
	INT	nils:	2;
	INT	sshp:	1;
	INT	cmps:	1;
	INT	vtwt:	1;
};

struct JUXT {
	int	jx;
	int	jy;
	int	jlen;
	int	nid;
	VTYPE	jscr;
};

inline	int	isDRNA(int a) {
	return (a == DNA || a == RNA || a == GENOME);
}

extern	ALPRM	alprm;

static	const	char	SensChar[] = {NM, CM, RV, CR, '\0'};
static	const	INEX	def_inex = {0};

/******** cDNA sequence *********
	TRANSCRIBEDaaaaaaaaaaa
	<--tlen---><---tail-->
	<---------len-------->
	tttttttttttTRANSCRIBED
	<--lead---><---tlen-->
	<---------len-------->
*********************************/

// basic single sequence class

struct OUTPRM {
	int	lpw;
	INT	BlkSz;
	int	NoOut;
	INT	MaxOut;
	int	EijMergin;
	int	AllowdOverlap;
const	char*	out_file;
	INT	RemoveTmp:	1;
	INT	trimend:	1;
	INT	SkipLongGap:	2;	// 0: don't skip, 1: skip, 3: run time
	INT	deflbl:		2;	// 0: sname, 1: fname, 2: tally, 3: new
	INT	fastanno:	1;	// 0: add annotation in fasta output
	INT	descrp:		1;	// output description
	INT	sortodr:	2;	// 0: input, 1: bytree, 2: shorter, 3: longer
	INT	spjinf:		1;
	INT	olrsum:		1;
	INT	ColorEij:	2;	// color mark intron position
	INT	supself:	2;	// suppress result of self comparison
	INT	noseqline:	1;	// additional prefix line in m-fasta
	INT	asterisk:	1;	// add asterisk as the eos mark
	INT	trimendgap:	1;	// suppress tail gap characters
	INT	taxoncode:	3;	// add taxinomic code field X in gnm2tab
};

extern	OUTPRM	OutPrm;

#define strmatch(str, pat) (!strncmp(str, pat, sizeof(pat) - 1))

inline	bool	IsGap(CHAR x) {return (x <= gap_code);}
inline	bool	IsTrueGap(CHAR x) {return (x == gap_code);}
inline	bool	IsntGap(CHAR x) {return (x > gap_code);}
inline	bool	IsAA(CHAR x) {return (x >= ALA && x <= VAL);}
inline	bool	IsTerm(CHAR x)	{return (x == TRM || x == TRM2);}
inline	int	getlpw()	{return (OutPrm.lpw);}

class Seq {
protected:
	int	base_;		// global numering base
	int	area_;		// memory size ocupied
	CHAR*	seq_;		// sequence body
	CHAR*	end_;		// end of seq_ memory
	Seq**	anti_;		// rlocation of everse-complement of this
	StrHash<int>*	mnhash;
	void	fillpad();
	RANGE*	setrange(const char* pa, int* ncr = 0);
	int	getcds(FILE* fd, char* str, int cdscolumn);
	char*	readanno(FILE* fd, char* str, SeqDb* db, Mfile& gapmfd);
	void	estimate_len(FILE* fd, int nos);
	void	header_nat_aln(int n, FTYPE sumwt);
	CHAR*	seq_realloc();
	int	calcResNum(int i);
public:
	int	sid;		// seq id.
	int	did;		// databse id number
	int	vrtl;		// hare seuence memory with other instance
	int	many;		// number of sequences = 1
	int	byte;		// = many
	int	len;		// seq length
	int	tlen;		// transcribed length without poly A or T
	int	left;		// left boundary to be operated
	int	right;		// right boundary to be operated
	int	CdsNo;		// number of exons or hsps
	int	wllvl;		// reserved
	VTYPE	jscr;		// sum of weight2s
	int*	nbr;		// position of the first residue
	SEQ_CODE*	code;	// code table
	Strlist*	spath;	// path name read from
	Strlist*	sname;	// sequence name
	Strlist*	descr;	// descrption
	FTYPE*	cmps;		// composition
	JUXT*	jxt;		// HSPs
	RANGE*	exons;		// CDS/exons in a gene seq
	Exinon*	exin;		// statistical scores related to gene organization
	SigII*	sigII;		// eson-intron junction signal
	INEX	inex;		// internally used flags
#if USE_WEIGHT
	FTYPE*	weight;		// relative weight to each member
	FTYPE*	pairwt;		// weight to each pair of members
	void	fpweight();	// print weight
	FTYPE*	saveseqwt();	// save weight vector
	void	restseqwt(FTYPE* tmpwt);// retore weight vector
	void	copyweight(Seq* dest);
#endif
	int	isAmb(CHAR r) const;
	bool	isGap(CHAR r)	{return r == gap_code || r == nil_code;}
	bool	isGap(CHAR* s) {
	    int n = 0; int i = 0;
	    for ( ; i < many; ++i) if (!IsGap(*s++)) break;
	    return (i == n);
	}
	bool	isdrna()	{return inex.molc == DNA ||
		inex.molc == RNA || inex.molc == GENOME;}
	bool	isprotein()	{return inex.molc == PROTEIN;}
	bool	istron()	{return inex.molc == TRON;}
	bool	empty()		{return left == right;}
	INT	r2s(INT r)	{return (isprotein()? r + ALA: ((1 << r) + _));}
	void	setmolc(int molc);
const	char	Strand()	{return inex.sens? '-': '+';}
const	char*	path2fn(const char* pname);
const	char*	sqname(bool fpri = false) {
		if (!fpri) fpri = many > 1;
		if (sname && *(*sname)[0] && !(fpri && spath)) return (*sname)[0];
		if (spath) return path2fn((*spath)[0]);
		return ("");
	}
	Seq*	attrseq(const char* pa);
	Seq*	splice(Seq* dest, RANGE* rng, int edit);
	int	siteno(int n)	{return base_ + n;}
	int	SiteNm(int n)	{return base_ + ((inex.sens & REVERS)? len - 1 - n: n);}
	int	SiteNz(int n)	{return base_ + ((inex.sens & REVERS)? len - n: n);}
	int	SiteLe()	{return base_ + ((inex.sens & REVERS)? len - right: left);}
	int	SiteRe()	{return base_ + ((inex.sens & REVERS)? len - left: right);}
	int	SiteNoBO(int n)	{return (inex.sens & REVERS)? len - n: n + 1;}
	int	SiteNo(int n)	{return base_ + SiteNoBO(n);}
	int	SiteOn(int n)	{return (inex.sens & REVERS)? len - n + base_: n - 1 - base_;}
	int	senschar()	{return SensChar[inex.sens];}
	int	index(CHAR* ss)	{return (ss - seq_) / many;}
	void	seqalloc(int num, int length, bool keep = false);
	void	fullrange()	{left = 0; right = len;}
	void	saverange(RANGE* rng)	{rng->left = left; rng->right = right;}
	void	restrange(RANGE* rng)	{left = rng->left; right = rng->right;}
	void	refresh(const int& num = 0, const int& length = 0);
	void	exg_seq(int gl, int gr);
	CHAR*	at(const int n)	{return seq_ + many * n;}
	void	rev_attr();
	void	reverse();
	void	comple();
	void	comrev() {reverse(); comple();}
	void	comrev(Seq** sqs);
	void	setanti(Seq** cmpl)	{anti_ = cmpl;}
	void	copyattr(Seq* dest) const;
	Seq*	cutseq(Seq* dest, int snl = CPY_ALL);
	Seq*	copyseq(Seq* dest, int snl = CPY_ALL);
	Seq*	duplseq(Seq* dest);
	Seq*	extseq(Seq* dest, int* which, int snl = CPY_ALL, FTYPE nfact = 1.);
	Seq*	aliaseq(Seq* dest, bool thisisalias = false);
	Seq*	catseq(Seq* tail, int cushion = 0);
	Seq&	operator=(Seq& src) {return (*src.copyseq(this));}
	Seq*	rndseq();
	Seq*	deamb(int bzx = 3);
	void	fillnbr(Seq* dest);
	void	copynbr(Seq* dest);
	void	copylbl(Seq* dest);
	Seq*	postseq(const CHAR* last);
	void	nuc2tron();
	void	tron2nuc(bool rev);
	CHAR	tron2aa(const CHAR res) {return inex.molc == TRON? t2atab[res]: res;}
	void	elim_amb();
	GAPS*	elim_column(int which, float frac = 0.);
	void	push(CHAR c, CHAR*& ps, int step = 1) {
		    if (ps >= end_) ps = seq_realloc();
		    *ps = c; ps += step;
		}
	CHAR*	seq_readin(FILE* fd, char* str, int mem, RANGE* pcr, Mfile* pfqmfd = 0);
	CHAR*	get_seq_aln(FILE* fd, char* str, RANGE* pcr);
	CHAR*	get_mfasta(FILE* fd, long fpos, char* str, RANGE* pcr, SeqDb* dbf);
	int	infermolc(FILE* fd, char* str, bool msf = false);
	CHAR*	ToInferred(CHAR* src, CHAR* lastseq, int step);
	Seq*	fgetseq(FILE* fd, const char* attr = 0, const char* attr2 = 0);
	int	fget(FILE* fd, const char* fn = 0) {return (fgetseq(fd)? 1: EOF);}	// alias 
	Seq*	getseq(const char* str, DbsDt* dbf = 0);
	char*	sqline(int i, char* ps);
	CHAR*	get_nat_aln(FILE* fd, RANGE* pcr);
	CHAR*	get_msf_aln(FILE* fd, RANGE* pcr);
	CHAR*	read_dbres(FILE* fd, RANGE* rng);
	CHAR*	read_dbres(CHAR* dbs, RANGE* rng);
	Seq*	read_dbseq(DbsDt* dbf, DbsRec* rec, RANGE* rng);
	Seq*	read_dbseq(DbsDt* dbf, long pos);
	Seq*	getdbseq(DbsDt* dbf, const char* code, int c = -1, bool readin = true);
	Seq*	apndseq(char* aname);
	void	fphseq(int n, FILE* fd = 0);
	FTYPE*	composition();
	void	printseq(FILE* fdi, int);
	void	fpmem_len(FILE* fd);
	FILE*	openseq(const char* str);// read from named file
	int	calcnbr(int gp, int i);
	void	listseq(int j, bool sub = false);
	void	num2pos(int which, int* array);
	void	pos2num(int which, int* array);
	bool	intronless(Seq* a);
	void 	putSigII();
	bool	findGate(RANGE* gate);
	void	typeseq(FILE* fd = 0, bool in_line = false);
	ORF*	getorf();
	void	passcom(FILE* fo);
	Seq*	translate(Seq* aas, ORF& orf);
	void	ftranslate(FILE* fd, int at_term, SeqDb* form, int orfn);
	void	fmtranslate(FILE* fd, int at_term, int orfn);
	int	transout(FILE* fd, int at_term, int orfn);
	void	setstrand(const int id_or_com, const char* text);
	void	test_gap_amb();
	bool	test_aligned();
	void	pregap(int* gl);
	bool	isgap(CHAR* ps);
	bool	nogap(CHAR* ps);
	VTYPE	countunps();
	int	pfqPos(int n) {return isprotein()? 3 * n: n;}
	int	sname2memno(const char* memid);
	Seq*	randseq(double* pcmp);
	Seq*	substseq(int n, int which);
	void	initialize();
	Seq(const int& num = 1, const int& length = 0);
	Seq(const char* fname);
	Seq(Seq& sd, int* which = 0, int snl = CPY_ALL);
	Seq(Seq* sd, int* which = 0, int snl = CPY_ALL);
				// copy or extract
	~Seq();
friend	void 	antiseq(Seq** seqs);
};

// supplier of one or two sets of sequences in series

class SeqServer {
	int	argc;
	int	argc0;
const	char**	argv;
const	char**	argv0;
	FILE*	fd[2];
	FILE*	fc;
	int	nfrom;
	int	nto;
	char*	cfrom;
	char*	cto;
	int	counter;
	bool	sw;
	int	molc[2];
	char*	attr[2];
	int	atsz[2];
public:
	InputMode	input_form;
	int	input_ns;
	DbsDt*	target_dbf;
	DbsDt*	query_dbf;
	SeqServer(int ac, const char** av, InputMode infm, 
	    const char* catalog = 0, int mq = UNKNOWN, int mt = UNKNOWN);
	~SeqServer() {
	    if (fd[0]) fclose(fd[0]);
	    if (fd[1]) fclose(fd[1]);
	    if (fc) fclose(fc);
	    delete[] attr[0];
	    delete[] attr[1];
	}
	InSt nextseq(Seq* sd, int which = 0);
	void	reset();
	int	getmolc(int i = 0) {return (molc[i]);}
	size_t	total_seq_len(Seq* sd);
};

// phrases in headline comments that specify the orientation of the seq

struct STRPHRASE {int sens; char* phrase;};

class StrPhrases {
	STRPHRASE*	strphrase[2];
	char*	strprhase_buf;
public:
	StrPhrases(const char* aname);
	~StrPhrases() {delete[] strprhase_buf; delete[] strphrase[0];}
friend	void	Seq::setstrand(const int idorcom, const char* text);
};

struct ExonRecord {
	int	Elen;
	int	Nmmc;
	int	Nunp;
	int	Rleft;
	int	Rright;
	int	Gleft;
	int	Gright;
	int	Ilen;
	int	Bmmc;
	int	Bunp;
	int	miss;
	int	phase;
	float	Pmatch;
	float	Escore;
	float	Iscore;
	float	Sig3;
	float	Sig5;
	char	Iends[4];
};

struct VulgRecord {short type; short len;};

struct GeneRecord {
	int	Cid;
	int	Gstart;
	int	Gend;
	INT	Nrecord;
	INT	nexn;
	int	Rid;
	int	Rlen;
	int	Rstart;
	int	Rend;
	int	mmc;
	int	unp;
	int	bmmc;
	int	bunp;
	int	ng;
	float	Gscore;
	float	Pmatch;
	float	Pcover;
	short	Csense;
	short	Rsense;
};

class SeqItr {
public:
	int	many;
	int	pos;
	CHAR*	res;
	EXIN*	bb;
	SeqItr& operator++()
	{
	    ++pos;
	    res += many;
	    if (bb) ++bb;
	    return (*this);
	}
	SeqItr& operator--()
	{
	    --pos;
	    res -= many;
	    if (bb) --bb;
	    return (*this);
	}
	SeqItr& operator+=(int n)
	{
	    pos += n;
	    res += n * many;
	    if (bb) bb += n;
	    return (*this);
	}
	SeqItr& operator-=(int n)
	{
	    pos -= n;
	    res -= n * many;
	    if (bb) bb -= n;
	    return (*this);
	}
	SeqItr operator+(int n)
	{
	    SeqItr temp(*this);
	    temp.pos += n;
	    temp.res += n * many;
	    if (bb) temp.bb += n;
	    return (*this);
	}
	SeqItr& operator-(int n)
	{
	    SeqItr temp(*this);
	    temp.pos -= n;
	    temp.res -= n * many;
	    if (bb) temp.bb -= n;
	    return (*this);
	}
	SeqItr&	sqset(CHAR* ss, int n = 0)
	{
	    res = ss + many * n;
	    pos = n;
	    return (*this);
	}
	bool	operator==(SeqItr& b) {return res == b.res;}
	bool	operator<(SeqItr& b) {return res < b.res;}
	bool	operator>(SeqItr& b) {return res > b.res;}
	bool	operator<=(SeqItr& b) {return res <= b.res;}
	bool	operator>=(SeqItr& b) {return res >= b.res;}
	bool	dullend() {
	    return (res[0] == nil_code || res[-many] == nil_code);
	}
	SeqItr&	reset(int n = 0, Seq* sd = 0) {
	    if (sd) {
		many = sd->many; res = sd->at(n);
		bb = sd->exin? sd->exin->score(n): 0;
	    } else if (res) {
		int shft = n - pos;
		res += shft * many;
		if (bb) bb += shft;
	    } else {
		many = 0; res = 0; bb = 0;
	    }
	    pos = n;
	    return (*this);
	}
	SeqItr(Seq* sd = 0, int n = 0) {res = 0; reset(n, sd);}
	SeqItr(SeqItr& src) {*this = src;}
};

class PrintMember {
	Strlist*	sname;
	Gnm2tab*	g2t;
	char    fmt[MAXL];
public:
	PrintMember(Strlist* sn, bool pad_space = false, const char* tl = 0);
	void    put_member(FILE* fd, int i);
	char*	operator[](int i) {
	    if (!g2t) return (0);
	    char*	taxon;
	    g2t->taxon_code((*sname)[i], &taxon);
	    return (taxon);
	}
};

class PrintAln {
	GAPS**  gaps;
	Seq**   seqs;
	int     seqnum;
	CHAR**  seq;
	GAPS**  gp;
	CHAR**  wkr;
	int*    nbr;
	int*    left;
	CHAR**  image;
	CHAR*   wbuf;
	int     rows;
	int     globpt;
	Row_Mode	cpm;
	int     markeij;
	EscCharCtl*     ecc;
	HtmlCharCtl*    hcc;
	PFQ**   pfqs;
	PFQ**   tfqs;
	int**   lsts;
	int*    agap;
	int     htl;
	int     pro;
	int     gene;
	void    markiis(int k, int j, int clm);
	void    prnt_aln(const int* lft, const int* rght);
	void    seqline(const CHAR* qry, const CHAR* brc,const  CHAR* cur,
		int rpl, const char decode[]);
	void    putclmark();
	void    pair_mrk(const CHAR* p1, const CHAR* p2, int clms);
	void    calc_mrk(int rows, int clms,
		const SEQ_CODE* defcode);
public:
	void    printaln();
	PrintAln(GAPS* _gaps[], Seq* _seqs[], int _seqnum);
	~PrintAln();
};

extern	int 	noseq;
extern	CHAR	ncredctab[];
extern	CHAR	ncelements[];
extern	CHAR	tnredctab[];
extern	CHAR	nccmpctab[];
extern	CHAR	aacmpctab[];
extern	CHAR	aaredctab[];
extern	CHAR	trelements[];
extern	CHAR	nccode[];
extern	CHAR	aacode[];
extern	CHAR	trccode[];
extern	CHAR	rdcode[];
extern	char	sensseq[];
extern	char	nucl[];
extern	char	nucr[];
extern	char	amino[];
extern	char	acodon[];
extern	char	ncodon[];
extern	char	amino3[][4];
extern	char*	Nucl;
extern	char*	Amin;
extern	char*	Tron;
extern	char	rdnucl[];
extern	CHAR	red2nuc[];
extern	CHAR	complcod[];
extern	CHAR	aa2nuc[];

/*	Headers to seq.c	*/

extern	int	setdefmolc(int molc = QUERY);
extern	void	setthr(double thr);
extern	void	setminus(int yon);
extern	int	setoutmode(int nsa);
extern	void	setalgmode(int nsa, int reg);
extern	void	setstrip(int sh);
extern	void	initseq(Seq** seqs, int n);
extern	void	clearseq(Seq** seqs, int n);
extern	void	cleanseq(Seq** seqs, int n);
extern	void	setdfn(const char* newdfn);
extern	void	antiseq(Seq** seqs);
extern	SEQ_CODE* setSeqCode(Seq* sd, int molc);
extern	CHAR*	spliceTron(CHAR* spliced, CHAR* b5, CHAR* b3, int n);
extern	int	Nprim_code(int c);
extern	int	en_code(int c, SEQ_CODE* code);
extern	CHAR*	tosqcode(CHAR* ns, SEQ_CODE* code);
extern	Seq*	inputseq(Seq* seqs[], char* str);
extern	int	samerange(Seq* a, Seq* b);
extern	void	swapseq(Seq** x, Seq** y);
extern	void	save_range(Seq* seqs[], RANGE* temp, int n);
extern	void	rest_range(Seq* seqs[], RANGE* temp, int n);
extern	void	eraStrPhrases();
extern	int	infermolc(const char* fname);
extern	InputMode	get_def_input_mode();

/*	Headers to sqpr.c	*/

extern	int	setlpw(int lpwd);
extern	INT	setdeflbl(int msf);
extern	FILE*	setup_output(int omode = 0, const char* def_fn = 0);
extern	void	close_output();
extern	void	GBcdsForm(RANGE* rng, Seq* sd);
extern	void	closeGeneRecord();
extern	void	setprmode(int pmd, int lorn, int trc);
extern	void	fphseqs(Seq* seqs[], int n);
extern	void	fprint_seq_mem(Seq* seqs[], int n);
extern	void	pralnseq(GAPS* gaps[], Seq* seqs[], int seqnum);

inline	VTYPE axbscale(Seq* seqs[])
{
	return (VTYPE) (alprm.scale * seqs[0]->many * seqs[1]->many);
}

#endif	// _BSEQ_H_
