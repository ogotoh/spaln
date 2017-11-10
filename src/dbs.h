/*****************************************************************************
*
*	Type definitions for various database formats
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-2013)
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

#ifndef  _DBS_
#define  _DBS_

#ifndef _COMM
#define _COMM	';'
#endif 

#define	SRC_EXT	".dt"
#define	SEQ_EXT	".seq"
#define	IDX_EXT	".idx"
#define	HED_EXT	".hed"
#define	GRP_EXT	".grp"
#define	ENT_EXT	".ent"
#define	ODR_EXT	".odr"
#define	BKA_EXT	".bka"
#define	BKN_EXT	".bkn"
#define	BKP_EXT	".bkp"

#ifndef	DBS_DIR
#define	DBS_DIR		"/home/gotoh/pub/spaln2.3.1/src/../seqdb"
#endif

#ifndef	DBS_SDIR
#define	DBS_SDIR	"/home/gotoh/pub/spaln2.3.1/src/../seqdb"
#endif

#ifndef ALN_DBS
#define	ALN_DBS		"ALN_DBS"
#endif

#ifndef ALN_SDBS
#define	ALN_SDBS	"ALN_SDBS"
#endif

static	const	char	DBSID = '$';
static	const	int	MAX_DBS = 2;
static	const	int	ENTLEN = 14;
static	const	CHAR	SEQ_DELIM = 0x00;
static	const	int	MAXCODE	= 21;
static	const	int	NTESTC = 1000;
static	const	long	magicver21 = 1117114721;

struct SeqDb {
	int	FormID;		/* ID of Database Format    */
	int	defmolc;	/* Default moleclur type    */
const	char*	DbName;
const	char*	EntLabel;	/* Entry	*/
const	char*	DefLabel;	/* Definition	*/
const	char*	AccLabel;	/* Accession	*/
const	char*	KeyLabel;	/* Key Word	*/
const	char*	SouLabel;	/* Source	*/
const	char*	RefLabel;	/* Reference	*/
const	char*	AutLabel;	/* Reference Authors	*/
const	char*	TitLabel;	/* Reference Title	*/
const	char*	JouLabel;	/* Reference Journal	*/
const	char*	ComLabel;	/* Comments	*/
const	char*	FeaLabel;	/* Feature Table	*/
const	char*	SeqLabel;	/* Begining of Sequence	*/
const	char*	EndLabel;	/* End of an Entry	*/
const	char*	SeqHead;	/* Head line on Seqence	*/
const	char*	SeqForm;	/* Format for Left Margin*/
	int	SeqBlkNo;	/* No. of Blocks	*/
	int	SeqBlkSz;	/* Size of a Block	*/
	int	SeqBlkSp;	/* Size of an Inter-Block space	*/
	INT	ContSpc;	/* Continue Space	*/
	int	is_DbEntry(char* str);
	int	is_DbEnd(char* str);
	int	is_DbOrigin(char* str);
	long	dbnextentry(FILE* fd, char* ps);
};

struct DbsRec20 {
	char	entry[ENTLEN];
	CHAR	master;
	CHAR	subset;
	long	recnbr;
	long	srcptr;
	long	seqptr;
	long	seqlen;
};

struct DbsRec {
	long	seqptr;
	size_t	seqlen;
	size_t	entptr;
};

struct DbsGrp {
	long	seqptr;
	long	recnbr;
};

class DbsDt {
friend	class	MakeBlk;
friend	class	SrchBlk;
	char*	pseq;
	CHAR*	dbsseq;
	Strlist*	grplbl;
	bool	comment;
	DbsRec*	recidx;
	int*	gsiidx;
	DbsRec20*	recidx20;
	INT*	recodr;
	char*	entry;
	int*	gsipool;
	DbsRec*	bisearch(const char* key);
	void	readseq();
	void	readgrp(FILE* fgrp);
	void	readidx20(FILE* fd);
	DbsRec*	readidx(FILE* fidx);
	void	readentry(FILE* fent);
	void	readodr(FILE* fodr);
public:
const	char*	dbsid;
	SeqDb*  curdb;
	FILE*	fseq;
	DbsGrp* dbsgrp;
	int     numgrp;
	INT	numidx;
	void	clean();
	DbsDt(int c = 0, int molc = UNKNOWN);
	DbsDt(const char* form);
	~DbsDt();
	int	grpno(DbsGrp* grp) {return (grp - dbsgrp);}
	int	recno(const DbsRec* rec) {return (rec - recidx);}
	DbsGrp*	finddbsgrp(const char* name);
	DbsRec*	findcode(const char* code);
	DbsRec*	dbsrec(INT pos) {return (pos < numidx? recidx + pos: 0);}
	DbsRec*	dbsrec(const DbsRec* rec = 0) {return (recidx + (rec? recodr[recno(rec)]: 0));}
	int	guessmolc();
	CHAR*	dbseq(DbsRec* rec = 0) {return (dbsseq + (rec? rec->seqptr: 0));}
	char*	entname(const DbsRec* rec = 0) {return (entry + (rec? rec->entptr: 0));}
	char*	entname(int recno) {return (entry + recidx[recno].entptr);}
	int*	gsient(int recno) {return (gsiidx? gsipool + gsiidx[recno]: 0);}
	int	gsisize(int recno) {return (gsiidx? gsiidx[recno + 1] - gsiidx[recno]:0);}
	char*	fsrcname(int sub) {return ((*grplbl)[sub]);}
	FILE*	dbsfopen() {return (dbsseq? 0: fopen(pseq, "r"));}
	void	prepare(size_t entry_space, size_t num, size_t seq_space, size_t gsi_space = 0);
	void	prepare(Strlist& sname, int num, size_t space, size_t gsi_space = 0);
	int	seqloc(CHAR* ps) {return (ps - dbsseq);}
	int	entloc(char* pe) {return (pe - entry);}
	int	gsiloc(int* pg) {return (pg - gsipool);}
	void	makodr();
};

enum DBs {GenBank, EMBL, Swiss, TFDS, NBRF, ProDB, 
		FASTA, MSF, PIR, NEXUS, Bare, EndOfDB = Bare};

extern	DbsDt*	dbs_dt[];
extern	SeqDb	SeqDBs[];

extern	SeqDb*	whichdb(char* str, FILE* fin = 0);
extern	SeqDb*	setform(int c);
extern	void	setdefdbf(DbsDt* dbf);
extern	void	EraDbsDt();
extern	char*	path2dbf(char* str, const char* fn, const char* ext = 0);
extern	const	char*	finddbfpath(const char* fn, const char* ext = 0);

#endif
