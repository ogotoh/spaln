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
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef  _BSEQ_H_
#define  _BSEQ_H_

#include "cmn.h"

struct	ALPRM {float u, v, u0, u1, v0, tgapf, thr, scale, maxsp, gamma; 
	    int k1, ls, sh, ubh, mtx_no;};
struct	ALPRM2 {float x, y, z, o, w, bti, spb, Z, sss; int jneibr, termk1;};
struct	ALPRM3 {float scnd, hydr, hpmt; int hpwing, no_angle;};
struct	DefSetup {int defmolc, delamb; InputMode def_input_mode;};

enum {SINGLE_SQ, NATIVE_MF, SEQUENTIAL_MF, SEQUENTIAL_FA, GCG_MSF};

extern	int 	noseq;
extern	ALPRM2	alprm2;
extern	ALPRM3	alprm3;
extern	ALPRM	alprm;
extern	RANGE	fullrng;

#define	USE_ETHER 	FVAL

static	const	char	_NHEAD = '>';		/* Normal Strand	*/
static	const	char	_CHEAD = '<';		/* Compl. Strand	*/
static	const	char	_APPN = '+';		/* Append on 		*/
static	const	char	_UNP = '-';		/* Gap			*/
static	const	char	_TRM = '*';		/* Termination Code	*/
static	const	char	_EOS = '/';		/* End of Seq.		*/
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

#include "dbs.h"
#include "gsinfo.h"
#include "gaps.h"
#include "codepot.h"
#include "utilseq.h"

static	const	int	DEFSEQLEN = 1024;
static	const	int	MAXCR = 512;
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

/* attribute characters	*/

inline	bool isnbr(const char* ps)
{
	return (isdigit(*ps) || ((*ps == '-' || *ps == '+') && isdigit(ps[1])));
}
inline	bool isseqchar(int c)
{
	return (isalpha(c) || c == _UNP || c == _TRM);
}
inline	int aton(int c)
{
	return (c - 'A' + A);
}

static	const	char	StrandPhrase[] = "StrandPhrase";
static	const	int	def_polya_thr = 12;
static	const	int	tribial_area = 64;
static	const	int	DelFrac	= 90;

	/* GenBank aa coding region */
static	const	int	typecolumn = 5;
static	const	int	datacolumn = 21;
static	const	int	MinPctNucChar = 75;
static	const	int	MinPctTronChar = 5;
static	const	INT	MaxTestChar = 1000;
static  const   char   NoSeqSpace[] = "No space for sequence!\n";

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
mutable	INT	exgl:	2;
mutable	INT	exgr:	2;
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

static	const	char	SensChar[] = {NM, CM, RV, CR, '\0'};

extern	int	cmpPfqPos(const PFQ* a, const PFQ* b);
extern	int	en_code(int c, const SEQ_CODE* code);

struct OUTPRM {
	int	lpw;
	INT	BlkSz;
	int	NoOut;
	INT	MaxOut;
	INT	MaxOut2;		// intermediate state
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
	INT	printweight:	1;	// output seq weights in MSA
	INT	gzipped:	1;	// gzipped output for spaln
	INT	debug:		1;	// print debug lines
};

extern	OUTPRM	OutPrm;
extern	DefSetup def_setup;
extern	SEQ_CODE* setSeqCode(Seq* sd, int molc);

inline	int	isDRNA(int a)
	{ return (a == DNA || a == RNA || a == GENOME); }
inline	bool	IsGap(const CHAR x) {return (x <= gap_code);}
inline	bool	IsTrueGap(const CHAR x) {return (x == gap_code);}
inline	bool	IsntGap(const CHAR x) {return (x > gap_code);}
inline	bool	IsAA(const CHAR x) {return (x >= ALA && x <= VAL);}
inline	bool	IsTerm(const CHAR x) {return (x == TRM || x == TRM2);}
inline	int	getlpw() {return (OutPrm.lpw);}
inline	bool	strmatch(const char* str, const char* pat) {
	while (*str && *pat) {if (*str++ != *pat++) return (false);}
	return (true);
}

/******** cDNA sequence *********
	TRANSCRIBEDaaaaaaaaaaa
	<--tlen---><---tail-->
	<---------len-------->
	tttttttttttTRANSCRIBED
	<--lead---><---tlen-->
	<---------len-------->
*********************************/

class PolyA {
	int	thr;
public:
	PolyA(int s = 0) : thr(s) {}	// remove polyA
	void	setthr(const char* s) {thr = *s? atoi(s): 0;}
	void	setthr(int v) {thr = v;}
	int	getthr() {return (thr);}
	SHORT	rmpolyA(Seq* sd, int q_mns = 1) const;
};

static	PolyA	polyA;

// basic single sequence class

class Seq {
protected:
	int	base_;		// global numering base
	int	area_;		// memory size ocupied
	CHAR*	seq_;		// sequence body
	CHAR*	end_;		// end of seq_ memory
	Seq**	anti_;		// rlocation of everse-complement of this
	StrHash<int>*	mnhash;
	void	fillpad();
	RANGE*	setrange(const char* pa, int* ncr = 0) const;
template <typename file_t>
	char*	readanno(file_t fd, char* str, SeqDb* db, Mfile& gapmfd);
	void	estimate_len(FILE* fd, int nos);
#if USE_ZLIB
	void	estimate_len(gzFile fd, int nos);
#endif
	void	header_nat_aln(int n, FTYPE sumwt);
	CHAR*	seq_realloc();
	int	calcResNum(int i);
	char*	onecds(RANGE& wexon, char* ps, int& par);
public:
	int	sid;		// seq id.
	int	did;		// databse id number
	int	vrtl;		// hare seuence memory with other instance
	int	many;		// number of sequences = 1
	int	byte;		// = many
	int	len;		// seq length
	int	tlen;		// transcribed length without poly A or T
mutable	int	left;		// left boundary to be operated
mutable	int	right;		// right boundary to be operated
	int	CdsNo;		// number of exons or hsps
	int	wllvl;		// reserved
	VTYPE	jscr;		// sum of weight2s
	int*	lens;		// seq length of each member
	int*	nbr;		// position of the first residue
	SEQ_CODE*	code;	// code table
	char*	spath;		// path name read from
	char*	msaname;	// msa name
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
	void	fpweight(FILE* fd = 0) const;	// print weight
	FTYPE*	saveseqwt() const;	// save weight vector
	void	restseqwt(FTYPE* tmpwt);// retore weight vector
	void	copyweight(Seq* dest) const;
#endif
	int	isAmb(CHAR r) const;
	bool	isGap(CHAR r) const {return r == gap_code || r == nil_code;}
	bool	isGap(CHAR* s) const {
	    int n = 0; int i = 0;
	    for ( ; i < many; ++i) if (!IsGap(*s++)) break;
	    return (i == n);
	}
	bool	isdrna() const {return inex.molc == DNA ||
		inex.molc == RNA || inex.molc == GENOME;}
	bool	isprotein() const {return inex.molc == PROTEIN;}
	bool	istron() const {return inex.molc == TRON;}
	bool	empty()	const {return left == right;}
	INT	r2s(INT r) const {return (isprotein()? r + ALA: ((1 << r) + _));}
	void	setmolc(int molc);
const	char	Strand() const {return inex.sens? '-': '+';}
const	char*	path2fn(const char* pname) const;
const	char*	sqname(bool fpri = false) const {
		if (!fpri) fpri = many > 1;
		if (sname && *(*sname)[0] && !(fpri && msaname)) return (*sname)[0];
		if (msaname && *msaname) return (msaname);
		if (spath) return path2fn(spath);
		return ("");
	}
	Seq*	attrseq(const char* pa);
	Seq*	splice(Seq* dest, RANGE* rng, const int edit = 0);
	int	siteno(int n) const {return base_ + n;}
	int	SiteNm(int n, int m = 0) const {return base_ + ((inex.sens & REVERS)? lens[m] - 1 - n: n);}
	int	SiteNz(int n, int m = 0) const {return base_ + ((inex.sens & REVERS)? lens[m] - n: n);}
	int	SiteLe(int m = 0) const {return base_ + ((inex.sens & REVERS)? lens[m] - right: left);}
	int	SiteRe(int m = 0) const {return base_ + ((inex.sens & REVERS)? lens[m] - left: right);}
	int	SiteNoBO(int n, int m = 0) const {return (inex.sens & REVERS)? lens[m] - n: n + 1;}
	int	SiteNo(int n, int m = 0) const {return base_ + SiteNoBO(n, m);}
	int	SiteOn(int n, int m = 0) const {return (inex.sens & REVERS)? lens[m] - n + base_: n - 1 - base_;}
	int	senschar() const {return SensChar[inex.sens];}
	int	index(CHAR* ss) const {return (ss - seq_) / many;}
	void	seqalloc(int num, int length, bool keep = false);
	void	fullrange() {left = 0; right = len;}
	void	saverange(RANGE* rng) const {rng->left = left; rng->right = right;}
	void	restrange(RANGE* rng) const {left = rng->left; right = rng->right;}
	void	refresh(const int& num = 0, const int& length = 0);
	void	exg_seq(int gl, int gr);
	CHAR*	at(const int n) const {return seq_ + many * n;}
	void	rev_attr();
	void	reverse();
	void	comple();
	void	comrev() {reverse(); if (!isprotein()) comple();}
	void	comrev(Seq** sqs);
	JUXT*	revjxt();
	void	setanti(Seq** cmpl) {anti_ = cmpl;}
	Seq**	getanti() {return (anti_);}
	void	copyattr(Seq* dest) const;
	Seq*	cutseq(Seq* dest, int snl = CPY_ALL) const;
	Seq*	copyseq(Seq* dest, int snl = CPY_ALL) const;
	Seq*	duplseq(Seq* dest) const;
	Seq*	extseq(Seq* dest, int* which, int snl = CPY_ALL, FTYPE nfact = 1.);
	Seq*	aliaseq(Seq* dest, bool thisisalias = false);
	Seq*	catseq(Seq* tail, int cushion = 0);
	Seq&	operator=(Seq& src) {return (*src.copyseq(this));}
	Seq*	rndseq();
	Seq*	deamb(int bzx = 3);
	void	fillnbr(Seq* dest) const;
	void	copynbr(Seq* dest) const;
	void	copylbl(Seq* dest) const;
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
template <typename file_t>
	CHAR*	seq_readin(file_t fd, char* str, int mem, RANGE* pcr, Mfile* pfqmfd = 0);
template <typename file_t>
	CHAR*	get_seq_aln(file_t fd, char* str, RANGE* pcr);
template <typename file_t>
	CHAR*	get_mfasta(file_t fd, long fpos, char* str, RANGE* pcr, const SeqDb* dbf);
template <typename file_t>
	int	infermolc(file_t fd, char* str, bool msf = false);
	CHAR*	ToInferred(CHAR* src, CHAR* lastseq, int step);
template <typename file_t>
	Seq*	fgetseq(file_t fd, const char* attr = 0, const char* attr2 = 0);
template <typename file_t>
	int	fget(file_t fd, const char* fn = 0) {return (fgetseq(fd)? 1: EOF);}	// alias 
	Seq*	getseq(const char* str, DbsDt* dbf = 0);
template <typename file_t>
	int	getcds(file_t fd, char* str, int cdscolumn);
	char*	sqline(int i, char* ps);
template <typename file_t>
	CHAR*	get_nat_aln(file_t fd, char* str, RANGE* pcr);
template <typename file_t>
	CHAR*	get_msf_aln(file_t fd, char* str, RANGE* pcr);
template <typename file_t>
	CHAR*	read_dbres(file_t fd, const RANGE* rng);
	CHAR*	read_dbres(CHAR* dbs, const RANGE* rng);
	Seq*	read_dbseq(DbsDt* dbf, DbsRec* rec, const RANGE* rng);
	Seq*	read_dbseq(DbsDt* dbf, long pos);
	Seq*	getdbseq(DbsDt* dbf, const char* code, int c = -1, bool readin = true);
	Seq*	apndseq(char* aname);
	void	fphseq(FILE* fd = 0, int n = 3) const;
	FTYPE*	composition();
	FTYPE*	composition(FTYPE* cmps) const;
	void	printseq(FILE* fdi, int);
	void	fpmem_len(FILE* fd);
#if USE_ZLIB
	FILE*	openseq(const char* str, gzFile* gzfd = 0);// read from named file
#else
	FILE*	openseq(const char* str);
#endif
	int	calcnbr(int gp, int i);
	void	listseq(FILE* fd = 0, int j = 0, bool sub = false);
	void	num2pos(int which, int* array) const;
	void	pos2num(int which, int* array) const;
	void 	putSigII(FILE* fd) const;
	bool	findGate(RANGE* gate);
	void	typeseq(FILE* fd = 0, bool in_line = false);
	ORF*	getorf() const;
	void	passcom(FILE* fo) const;
	Seq*	translate(Seq* aas, ORF& orf) const;
	void	ftranslate(FILE* fd, int at_term, SeqDb* form, int orfn) const;
	void	fmtranslate(FILE* fd, int at_term, int orfn) const;
	int	transout(FILE* fd, int at_term, int orfn) const;
	void	setstrand(const int id_or_com, const char* text);
	void	test_gap_amb();
	bool	test_aligned();
	void	pregap(int* gl) const;
	bool	isgap(const CHAR* ps) const;
	bool	nogap(const CHAR* ps) const;
	int	countgap(int mem = 0, int from = 0, int to = INT_MAX) const;
	VTYPE	countunps() const;
	int	pfqPos(int n) const {return isprotein()? 3 * n: n;}
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

template <typename file_t>
SeqDb* seq_NandL(int& num, int& len, int& mode, char* str, file_t fd, int molc)
{
// Force to single sequence

	if (mode == IM_SNGL) {
	    num = 1; mode = SINGLE_SQ;
	    return (0);
	}
// Givien Number & Length
	if ((num = atoi(str))) {
	    char*       ps = cdr(str);
	    if (ps && isdigit(*ps)) len = atoi(ps);	     // Phylip like
	    else if (ps && *ps) fatal("Unsupported format:\n%s\n", str);
	    if (num > 1) mode = SEQUENTIAL_MF;
	    else	mode = SINGLE_SQ;
	    return (whichdb(str, fd));
	}
// spaln output?
	int     m;
	char    c;
	int     n = sscanf(str, "%*s %*s %c [%d", &c, &m);
	if (n == 2 && m == 1 && (c == '+' || c == '-')) {
	    num = 1; mode = SINGLE_SQ;
	    Strlist     stl(str, stddelim);
	    int k = stl.size();
	    len = (stl[k - 2][0] == 'N' || stl[k - 2][0] == 'Q')?
		atoi(stl[k - 4]): 0;
	    return (0);
	}
// native mfa ?
	for (const char* ps = str; (ps = strchr(ps, '[')); ) {
	    num += atoi(++ps);
	    const char* qs = strchr(ps, ':');
	    int tl = qs? len += atoi(ps = ++qs): 0;
	    if (tl > len) len = tl;
	}
	if (num > 1) {mode = NATIVE_MF; return (0);}
	SeqDb*  dbf = whichdb(str, fd);
// MSF format ?
	if (dbf->FormID == MSF) {
	    mode = GCG_MSF;
	    return (dbf);
	}
// sequential mfa unkonw # of seqs
	if (mode == IM_MULT && dbf && dbf->FormID <= FASTA) {
	    num = 0;
	    mode = SEQUENTIAL_FA;
	    return (dbf);
	}
// single seq
	num = 1; mode = SINGLE_SQ;
	return (0);
}

template <typename file_t>
Seq* Seq::fgetseq(file_t fd, const char* attr, const char* attr2)
{
	int     dm = def_setup.defmolc;
	int     dela = def_setup.delamb;
	int     mode = def_setup.def_input_mode;
const   char*   attrs[3] = {attr, attr2, 0};
	for (const char** ars = attrs; *ars; ++ars) {
	    for (const char* as = *ars; as && *as; ++as) {
	      switch (toupper(*as)) {
		case 'N': dela = 1;
		case 'D': dm = DNA; break;
		case 'R': dm = RNA; break;
		case 'X': dela = 1;
		case 'A':
		case 'P': dm = PROTEIN; break;
		case 'T': dm = TRON; break;
		case 'G': dm = GENOME;
		case 'S': mode = IM_SNGL; break;
		case 'M': mode = IM_MULT; break;
		default:  break;
	      }
	    }
	}

	int     nos = len = 0;
	char    str[MAXL];
	long    fpos = 0L;
// skip comment lines
	for (;;) {
	    fpos = ftell(fd);
	    if (!fgets(str, MAXL, fd)) return (0);      // empty
	    if (*str != _LCOMM && !isBlankLine(str)) break;
	    if ((strlen(str) + 1) == MAXL) flush_line(fd);
	}

// infer input sequence format
	SeqDb*  dbf = seq_NandL(nos, len, mode, str, fd, dm);
	if (nos) {
	    if (!len) estimate_len(fd, nos);
	    refresh(nos, len);
	}
	setSeqCode(this, dm);
	inex.dela = dela;
	CHAR*   lastseq;
	RANGE   rng = {1, 0};
	int     ncr = 0;
	RANGE*  slices = setrange(attr, &ncr);
	RANGE*  pcr = slices? slices: &fullrng;
	switch (mode) {   // Multiple seqs ?
	    case SEQUENTIAL_FA:  // Sequential: unknown # of seqs
		lastseq = get_mfasta(fd, fpos, str, pcr, dbf); break;
	    case SEQUENTIAL_MF:  // Sequential: given # of seqs
		lastseq = get_seq_aln(fd, str, pcr); break;
	    case NATIVE_MF:       // Native
		lastseq = get_nat_aln(fd, str, pcr); break;
	    case GCG_MSF:	     // GCG MSF format
		lastseq = get_msf_aln(fd, str, pcr); break;
	    default:	    // Single seq
		lastseq = seq_readin(fd, str, 0, pcr, 0);
		while (lastseq > seq_ && IsGap(lastseq[-1]))
		    --lastseq;    // trim the end gaps
		break;
	}
	if (lastseq <= seq_) return (0);
	base_ = left;
	postseq(lastseq);
	if (sigII) {
	    sigII->step = isprotein()? 3: 1;
	    sigII->pfqrepos(slices);
	}
	delete[] slices;
	if (rng.left < rng.right) {
	    dm = --rng.left - base_;
	    if (dm > 0) left = dm;
	    dm = rng.right - base_;
	    if (dm < right) right = dm;
	}
	if (spath) {
	    if (!sname)  sname = new Strlist(path2fn(spath));
	    else if (sname->empty()) sname->assign(path2fn(spath));
	}
	return attrseq(attr);
}

template <typename file_t>
int Seq::getcds(file_t fd, char* str, int cdscolumn)
{
	Mfile	exonrng(sizeof(RANGE));
	RANGE	ebuf = fullrng;
	bool	rv = false;
	int	par = 0;
static	const	char* readthru = "CDS may be read through !\n";

	char*	ps = str + cdscolumn;
	if (*ps == 'j' && strmatch(ps, "join(")) {
	    ps += 5; ++par;
	} else if (*ps == 'c' && strmatch(ps, "complement(join(")) {
	    rv = true;
	    ps += 16;
	    par = 2;
	} else if (*ps == '-') {
	    rv = true; ++ps;
	} 
	char*	qs = strchr(ps, '(');
	if (qs) {
	    ++par; ps = qs + 1;
	} else if (!par) {
	    qs = onecds(ebuf, ps, par);
	    exonrng.write(&ebuf);
	}
	while (par > 0) {
	    int	n = 0;
	    while ((qs = onecds(ebuf, ps, par))) {
		if (ebuf.right > ebuf.left) exonrng.write(&ebuf);
		if (par < 0) prompt("Bad CDS inf.!\n");
		if (par <= 0) goto eoj;
		ps = qs;
		++n;
	    }
	    if (!n) {
		prompt(readthru);
		break;
	    }
	    ps = strcpy(str, ps);
	    int	len = strlen(ps);
	    if (!fgets(str + len, MAXL - len, fd)) {
		prompt(readthru);
		break;
	    }
	}
eoj:
	if (exons) delete[] exons;
	CdsNo = exonrng.size();
	exonrng.write(&zerorng);
	exons = CdsNo? (RANGE*) exonrng.flush(): 0;
	return (rv? -CdsNo: CdsNo);
}

template <typename file_t>
char* Seq::readanno(file_t fd, char* str, SeqDb* db, Mfile& gapmfd)
{
	char*	ps = str;
	int	feature = 0;
	int	nseg = 0;

	for ( ; ps && ((db->FormID >= FASTA &&
	    (*str == _COMM || *str == _NHEAD || *str == _CHEAD)) ||
	    (db->SeqLabel && wordcmp(str, db->SeqLabel)));
	    ps = fgets(str, MAXL, fd)) {
	    if (db->FormID == GenBank && 
		*str == 'F' && strmatch(str, "FEATURES")) {
		    feature = 1;
	    } else if ((db->FormID == GenBank && *str == ' ' && feature) ||
		     (db->FormID == EMBL && strmatch(str, "FT"))) {
		ps = str + typecolumn;
		if (*ps == 'C' && strmatch(ps, "CDS")) {
		    nseg = getcds(fd, str, datacolumn);
		}
	    } else if (db->FormID >= FASTA) {
		if (strmatch(str, ";C")) {
		    nseg = getcds(fd, str, 3);
		} else if (strmatch(str, ";M")) {
		    ps = cdr(str);
		    bool	delfs = strmatch(ps, "Deleted");
		    GAPS	gaps;
		    gaps.gln = atoi(ps = cdr(ps));
		    if (delfs) gaps.gln = -gaps.gln;
		    ps = cdr(ps); ps = cdr(ps);
		    gaps.gps = atoi(cdr(ps));
		    gapmfd.write(&gaps);
		}
	    }
	}
	if (exons && ((exons->left > exons[CdsNo - 1].left) ^ (nseg < 0)))
	    vreverse(exons, CdsNo);
	if (ps && db->FormID < FASTA) ps = fgets(str, MAXL, fd);
	return (ps);
}

template <typename file_t>
CHAR* Seq::get_nat_aln(file_t fd, char* str, RANGE* qcr)
{
	char*	ps = str + 1;
	msaname = strrealloc(msaname, car(ps));

	CHAR*	ss = seq_;
	long	fpos = ftell(fd);
	if (inex.molc == UNKNOWN) {
//	Infer molecular type
	    if (!fgets(str, MAXL, fd)) return (ss);
	    infermolc(fd, str, true);
	    fseek(fd, fpos, SEEK_SET);
	}

	CHAR**	wrk = new CHAR*[many];
	if (sigII) {delete sigII; sigII = 0;}
	int	column = 0;
	FTYPE	sumwt = 0.;
	inex.algn = 1;
	int	n = 0;
	int	blkno = 0;
	for (int i = 0; i < many; ) wrk[i++] = ss++;
	int	i = 0;
	for ( ; (ps = fgets(str, MAXL, fd)); fpos = ftell(fd)) {
	    if (*ps == _NHEAD || *ps == _CHEAD) {
		fseek(fd, fpos, SEEK_SET);	// undo 1 line
		break;
	    }
	    if (strmatch(ps, ";B") && (alprm2.spb > 0)) {
#if USE_WEIGHT
		sigII = new SigII(fd, ps, weight);
#else
		sigII = new SigII(fd, ps);
#endif
		if (!sigII->pfqnum) {
		    delete sigII; sigII = 0;
		}
	    }
	    if (*ps == _COMM) {
		while (withinline(ps, MAXL, fd)) ;
		continue;
	    }
	    chop(ps);
	    while (*ps && isspace(*ps)) ps++;
	    if (isnbr(ps) || *ps == _LABL) {
		ps = sqline(i, ps);
#if USE_WEIGHT
	    } else if (*ps++ == _WGHT) {
		if (!weight) {
		    n = 0;
		    weight = new FTYPE[many];
		}
		while (*ps && n < many) {
		    sumwt += weight[n++] = atof(ps);
		    ps = cdr(ps);
		}
		continue;
#endif
	    } else	continue;
	    RANGE*	pcr = qcr;
	    int 	clm = column;
	    for ( ; *ps; ++ps) {
		int	res_code = -1;
		switch (*ps) {
		  case '\0':
		  case _COMM:
		  case _NHEAD:
		  case _CHEAD:
		  case _LABL:	goto eob;
		  case ESC:
			while (!isalpha(*++ps)) ;
			if (*ps) ++ps;
			break;
		  case _SAME:	res_code = *wrk[0]; break;
		  case _IBID:	res_code = *wrk[i-1]; break;
		  default:
			res_code = en_code(*ps, code);
			if (res_code < 0) continue;
			if (IsGap(res_code)) inex.algn |= 2;
			break;
		}
		if (res_code >= 0) {
		    if (clm == pcr->right) ++pcr;
		    else if (pcr->left <= clm) {
			*wrk[i] = res_code;
			wrk[i] += many;
		    }
		    ++clm;
		}
	    }
eob:	    ; /* end of block --- this line is a dummy */
	    if (++i == many) {
		i = 0;
		qcr = pcr;
		column = clm;
		if (blkno++ == 0) header_nat_aln(n, sumwt);
	    }
	}
/* end of input */
	ss = wrk[0];
	delete[] wrk;
	return (ss);
}

template <typename file_t>
CHAR* Seq::get_msf_aln(file_t fd, char* str, RANGE* pcr)
{
	int	num  = 0;
	if (!sname) sname = new Strlist;
	while (fgets(str, MAXL, fd)) {
	    if (!wordcmp(str, "//")) break;
	    char*	qs = str;
	    char*	ps = car(qs);
	    if (!strcmp(ps, "MSF:")) {
		len = atoi(++qs);
		ps = cdr(qs); ps = cdr(ps);
		setSeqCode(this, *ps == 'P'? PROTEIN: DNA);
	    } else if (!strcmp(ps, "Name:")) {
		++qs;
		ps = car(qs);
		sname->push(ps);
		++num;
	    }
	}

	seqalloc(num, len);
	delete[]	lens;
	lens = new int[many];
	vclear(lens, many);
	delete[]	nbr;
	nbr = new int[many];
	vclear(nbr, many);
	CHAR*	ss = at(0);
	int	m = 0;
	inex.algn = 1;
	while (fgets(str, MAXL, fd)) {
	    char*	ps = cdr(str);
	    if (!*ps) continue;	// blank line
	    CHAR*	rr = ss + m;
	    for ( ; *ps; ++ps) {
		int	c = -1;
		if (isalpha(*ps)) c = en_code(*ps, code);
		else if (*ps == '.') {c = gap_code; inex.algn = 3;}
		if (c >= 0) {
		    *rr = c;
		    rr += many;
		}
	    }
	    if (++m == many) {
		ss = rr - many + 1;
		m = 0;
	    }
	}
	return (ss);
}

template <typename file_t>
CHAR* Seq::seq_readin(file_t fd, char* str, int mem, RANGE* pcr, Mfile* pfqmfd)
{
	SeqDb*	db = whichdb(str, fd);
	if (!db)	return (0);		// bad format
	GAPS gend = {INT_MAX, EOS};

/*************************
* read annotation
*************************/

	char*	ps = 0;
	Mfile	gapmfd(sizeof(GAPS));
	if (db->FormID < FASTA)	{		// Public 1-seq Format
	    char*	ss = cdr(str);
	    char*	ct = 0;
	    char	cc;
	    sname->push(car(ss, ct, cc));	// seq label
	    if (ct) *ct = cc;
	    ps = readanno(fd, str, db, gapmfd);
	} else if (db->FormID == FASTA)	{	// FASTA format
	    char*	ss = str;		// process header line
	    char*	sp = str;
	    while (*sp && !isspace(*sp)) ++sp;
	    if (*sp == '\n') *sp = 0;
	    else	*sp = '\0';
	    if (algmode.mns == 0) setstrand(0, ss);
	    while (char* sb = strchr(++ss, '|')) {
		if (sb[1] == '\0' || isspace(sb[1])) {
		    *sb = ' ';
		    break;
		}
		ss = sb;
	    }
	    sname->push(ss);
	    if (sp) {				// process comments in header
		*sp++ = ' ';
		do {
		    if (algmode.mns == 0 && inex.ori == 0)
			setstrand(1, sp); 
		} while ((sp = withinline(str, MAXL, fd)));
	    }
	    if (!fgets(str, MAXL, fd)) return (0);	// empty seq
	    ps = readanno(fd, str, db, gapmfd);
	} else
	    ps = readanno(fd, str, db, gapmfd);
	if (!ps)  return (0);

	int	ngaps = gapmfd.size();
	GAPS*	gg = 0;
	if (ngaps) {
	    gapmfd.write(&gend);
	    gg = (GAPS*) gapmfd.flush();
	}

/*************************
* read body of seq
*************************/

	long	fpos = ftell(fd);
	if (inex.molc == UNKNOWN) {
	    char	buf[MAXL];
	    strcpy(buf, str);
	    infermolc(fd, buf);
	    fseek(fd, fpos, SEEK_SET);
	}

	int	step = isprotein()? 3: 1;
	CHAR*	ns = at(0) + mem;
	int 	flag = -1;
	int	total = 0;
	int	unps = 0;
	int	nres = 0;
	int	eij = 0;	// cumulated exon length
	int	cds = 0;	// virtual CDS length != eij if frameshifts
	bool	rev = exons && (exons->left > exons[CdsNo - 1].left);
	RANGE*	bcr = 0;

	if (sscanf(ps, "%d", nbr + mem) > 0)	nbr[mem]--;
	else					nbr[mem] = 0;
	nbr[mem] -= pcr->left;
	GAPS*	gwk = gg;
	bool	readspb = alprm2.spb > 0 && inex.molc != GENOME;
	RANGE*&	cr = exons;
	if ((readspb = readspb && exons)) {
	    if (rev && ngaps) vreverse(gg, ngaps);
	    bcr = (inex.molc != UNKNOWN)? cr: 0;
	    eij = cds = cr->right - cr->left;
	    if (gg) {
		if (rev) {
		    while (gaps_intr(gwk) && gwk->gps > cr->left)
			cds += (gwk++)->gln;
		} else {
		    while (gaps_intr(gwk) && gwk->gps < cr->right)
			cds += (gwk++)->gln;
		}
	    }
	}

	for ( ; ps; fpos = ftell(fd), ps = fgets(str, MAXL, fd)) {
	    if (db->is_DbEntry(ps)) {
		fseek(fd, fpos, SEEK_SET);	// undo 1 line
		break;
	    }
	    if (db->EndLabel && !wordcmp(str, db->EndLabel)) flag = -2;
	    if (flag == -2) continue;
	    while (int c = *ps++) {
		if (c == _COMM || c == _LCOMM) {
		    while (withinline(str, MAXL, fd)) ;
		    break;
		}
		if (c == _EOS && *ps == _EOS) {	/* // */
		    flag = -2; break;
		}
		if (c == ESC) {			// skip esc cntl
		    while (!isalpha(*++ps)) ;
		    if (*ps && (c = *++ps)) ++ps;
		    else	break;
		}
		if ((c = en_code(c, code)) >= 0) {
		    if (total == pcr->left) flag = 1;
		    else if (total == pcr->right) {
			flag = 0;
			pcr++;
		    }
		    ++total;
		    if (IsntGap(c)) ++lens[mem];
		    if (flag == 1) push(c, ns, many);
		    else if (flag == -1 && IsntGap(c)) ++nbr[mem];
		    if (bcr && bcr < cr + CdsNo) {
			if (IsGap(c))	unps += step;
			else		nres += step;
			if ((step == 3 && cds < nres + 2) || cds == nres) {
			    bcr->left = eij + unps;
			    ++bcr;
			    if (gg) {
				if (rev) {
				    while (gaps_intr(gwk) && gwk->gps > bcr->left)
					cds += (gwk++)->gln;
				} else {
				    while (gaps_intr(gwk) && gwk->gps < bcr->right)
					cds += (gwk++)->gln;
				}
			    }
			    eij += bcr->right - bcr->left;
			    cds += bcr->right - bcr->left;
			}
		    } else if (IsGap(c))	unps++;
		}
	    }
	}
	if (bcr && bcr < cr + CdsNo) bcr->left = eij + unps;
	inex.algn = unps? 3: 0;
	if (!readspb) {delete[] gg; return (ns);}

/***************************************
* store exon-intron junction information
***************************************/

	if (!sigII) sigII = new SigII;
	unps = 0;
	gwk = gg;
	if (pfqmfd) {
	    bcr = cr + CdsNo;
	    for (pcr = cr; pcr < bcr; ++pcr) {
		if (gg) {
		    while (gaps_intr(gwk) && (rev ^ (gwk->gps < pcr[rev].right)))
			unps += (gwk++)->gln;
		}
		PFQ	pfqbf = {pcr->left + unps, mem, pcr->right};
		pfqmfd->write(&pfqbf);
	    }
	} else {
	    PFQ*	pfq = sigII->pfq = new PFQ[CdsNo];
	    sigII->pfqnum = sigII->lstnum = CdsNo - 1;
	    bcr = cr + CdsNo;
	    for (pcr = cr; pcr < bcr; ++pcr, ++pfq) {
		if (gg) {
		    while (gaps_intr(gwk) && (rev ^ (gwk->gps < pcr[rev].right)))
			unps += (gwk++)->gln;
		}
		pfq->pos = pcr->left + unps;
		pfq->num = 1;
		pfq->gps = pcr->right;
#if USE_WEIGHT
		pfq->dns = weight? weight[mem]: 1;
#endif
	    }
	    (--pfq)->num = 0;
#if USE_WEIGHT
	    pfq->dns = 0;
#endif
	}
	delete[] gg;
	delete[] exons; exons = 0; CdsNo = 0;
	return (ns);
}

/*	Get multiple sequences in sequential format	*/

template <typename file_t>
CHAR* Seq::get_seq_aln(file_t fd, char* str, RANGE* pcr)
{
	size_t*	slen = new size_t[many];
	Mfile*	pfqmfd = new Mfile(sizeof(PFQ));

	size_t	span = 0;
	for (int i = 0; i < many; ++i) {
	    if (*str != _NHEAD && *str != _CHEAD) *str = '\0';
	    CHAR* ss = seq_readin(fd, str, i, pcr, pfqmfd);
	    if (!ss)
		fatal("Insufficient # of sequences %d < %d", i, many);
	    slen[i] = ss - i - seq_;
	    if (slen[i] > span) span = slen[i];
	}
	CHAR*	last = seq_ + span;
	for (int i = 0; i < many; ++i) {
	    CHAR*	 ss = seq_ + slen[i];
	    if (ss >= last) continue;
	    for ( ; ss < last; ss += many) *ss = UNP;
	}
	delete[] slen;
	int	n = pfqmfd->size();
	pfqmfd->write(&pfqend);
	PFQ*	pfq = (PFQ*) pfqmfd->flush();
	delete pfqmfd;
	if (sigII) {
	    int*	lst = new int[n];
	    sigII->lst = lst;
	    qsort((UPTR) pfq, (INT) n, sizeof(PFQ), (CMPF) cmpPfqPos);
	    PFQ*	wfq = pfq;
	    sigII->pfq = pfq;
	    sigII->lstnum = n;
	    for (int i = 0; i < n; ++i, ++wfq) {
		*lst++ = wfq->num;
#if USE_WEIGHT
		wfq->dns = weight? weight[wfq->num]: 1;
#endif
		wfq->num = 0;
	    }
	    wfq = pfq = sigII->pfq;
	    for (int i = 0; i < n; ++i, ++wfq) {
		if (wfq->pos != pfq->pos) {
		    *++pfq = *wfq;
		} else if ((wfq->gps < 0) ^ (wfq->gps > pfq->gps)) {
		   pfq->gps = wfq->gps;
		    ++pfq->num;
#if USE_WEIGHT
		    pfq->dns += wfq->dns;
#endif
		}
	    }
	    sigII->pfqnum = pfq - sigII->pfq;
	    *++pfq = *wfq;
	} else if (pfq) delete[] pfq;
	return (last);
}

template <typename file_t>
int Seq::infermolc(file_t fd, char* str, bool msf)
{
	INT	cmp[28];
	vclear(cmp, 28);
	INT	total = 0;
	bool	pending = false;

	do {
	    char*	ps = str;
	    if (*ps == _NHEAD || *ps == _CHEAD || *ps == _WGHT || 
		*ps == _COMM || *ps == _LCOMM || pending) {
		    pending = (strlen(str) + 1) == (INT) MAXL;
		    continue;
	    }
	    if (*ps == _EOS) break;
	    if (msf) {
		while (*ps && isspace(*ps)) ++ps;
		if (!isdigit(*ps)) continue;
	    }
	    for ( ; *ps; ++ps) {
		if (isalpha(*ps)) {
		    ++cmp[aton(toupper(*ps))];
		    if (++total >= MaxTestChar) goto eol;
		} else if (msf && *ps == _LABL) break;
		else if (*ps == _COMM || *ps == _LCOMM) break;
	    }
	} while (fgets(str, MAXL, fd));
eol:
	FTYPE	ncode = cmp[aton('A')] + cmp[aton('C')] + cmp[aton('G')]
		+ cmp[aton('T')] + cmp[aton('U')];
	int	molc = PROTEIN;
	INT	tmark = cmp[aton('J')] + cmp[aton('U')];
	if (total == cmp[aton('N')]) molc = DNA;
	else if (100. * ncode/(total - cmp[aton('N')] - cmp[aton('X')]) > MinPctNucChar)
	    molc = (cmp[aton('T')] >= cmp[aton('U')])? DNA: RNA;
	else if (tmark && 100. * (tmark + cmp[aton('O')]) / total > MinPctTronChar) {
	    molc = TRON;
	    prompt("Warning: %s is regarded as TRON sequence!\n", sqname());
	}
	setSeqCode(this, molc);
	return molc;
}

template <typename file_t>
CHAR* Seq::get_mfasta(file_t fd, long fpos, char* str, RANGE* pcr, const SeqDb* dbf)
{
	if (!dbf) fatal("Unknown sequence format !\n");
	int	num = 0;
	do {
	    if (dbf->is_DbEntry(str)) ++num;
	} while (fgets(str, MAXL, fd));
	fseek(fd, fpos, SEEK_SET);
	if (!num) {
	    prompt("Empty or not MSA !\n");
	    exit(0);
	}
	refresh(num);
	*str = '\0';
	return (many == 1?
	    seq_readin(fd, str, 0, pcr, 0):
	    get_seq_aln(fd, str, pcr));
}

// supplier of one or two sets of sequences in series

class SeqServer {
	int	argc;
	int	argc0;
const	char**	argv;
const	char**	argv0;
	FILE*	fd[2];
#if USE_ZLIB
	gzFile	gzfd[2];
#endif
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
	    delete[] cfrom;
	    delete[] cto;
	}
	InSt nextseq(Seq* sd, int which = 0);
	void	reset();
	int	getmolc(int i = 0) {return (molc[i]);}
	size_t	total_seq_len(Seq* sd, int* many = 0);
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
const	Strlist*	sname;
	Gnm2tab*	g2t;
	char    fmt[MAXL];
public:
	PrintMember(const Strlist* sn, bool pad_space = false, const char* tl = 0);
	void    put_member(FILE* fd, int i) const;
	char*	operator[](int i) const {
	    if (!g2t) return (0);
	    char*	taxon;
	    g2t->taxon_code((*sname)[i], &taxon);
	    return (taxon);
	}
};

class PrintAln {
const	GAPS**  gaps;
	Seq**   seqs;
	int     seqnum;
	CHAR**  seq;
const	GAPS**  gp;
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
	FILE*	fd;
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
	PrintAln(const GAPS** _gaps, Seq* _seqs[], int _seqnum, FILE* _fd = 0);
	~PrintAln();
};

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
extern	CHAR*	spliceTron(CHAR* spliced, const CHAR* b5, const CHAR* b3, int n);
extern	int	Nprim_code(int c);
extern	CHAR*	tosqcode(CHAR* ns, const SEQ_CODE* code);
extern	Seq*	inputseq(Seq* seqs[], char* str);
extern	int	samerange(const Seq* a, const Seq* b);
extern	void	swapseq(Seq** x, Seq** y);
extern	void	save_range(const Seq* seqs[], RANGE* temp, int n);
extern	void	rest_range(const Seq* seqs[], RANGE* temp, int n);
extern	void	eraStrPhrases();
extern	int	infermolc(const char* fname);
extern	InputMode	get_def_input_mode();

/*	Headers to sqpr.c	*/

extern	int	setlpw(int lpwd);
extern	INT	setdeflbl(int msf);
extern	FILE*	setup_output(int omode = 0, const char* def_fn = 0, bool setup_out_fd = true);
extern	void	close_output();
extern	void	closeGeneRecord();
extern	void	setprmode(int pmd, int lorn, int trc);
extern	void	fphseqs(const Seq* seqs[], int n = 3, FILE* fd = 0);
extern	void	GBcdsForm(const RANGE* rng, const Seq* sd, FILE* fd = 0);
extern	void	fprint_seq_mem(const Seq* seqs[], int n, FILE* fd = 0);
extern	void	pralnseq(const GAPS** gaps, Seq* seqs[], int seqnum, FILE* fd = 0);

inline	VTYPE axbscale(const Seq** seqs)
{
	return (VTYPE) (alprm.scale * seqs[0]->many * seqs[1]->many);
}

template <typename file_t>
inline	char*	withinline(char* str, INT maxl, file_t fd)
{
	return ((strlen(str) + 1) == maxl? fgets(str, maxl, fd): 0);
}

#endif	// _BSEQ_H_
