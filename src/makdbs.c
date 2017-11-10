/*****************************************************************************
*
*	Make unformat and index files from flat database files
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-4-7 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include <new>
#include "seq.h"

static	const	size_t	MAXL2K = 2048;
static	const	INT	ddelim = SEQ_DELIM + (SEQ_DELIM << 4);
static	bool	comment = false;
static	INT	max_memory = 1024 * 1024 * 1024;

class Makdbs {
	int	bias;
	int	ceil;
	SeqDb*	db;
	SEQ_CODE*	defcode;
	long	memsize;
	size_t*	fsize;
	int	memmode;
	FILE*	fsrc;
	char*	msrc;
	char*	mptr;
	FILE*	fgrp;
	FILE*	fseq;
	FILE*	fidx;
	FILE*	fent;
	size_t	recnbr;
	char	str[MAXL];
	char	prv[MAXL];
	char*	mss;
	char*	tss;
	bool	halfway;
	int	encode(int c);
	void	examem(int arc, const char** argv);
	int	convert();
	int	convert2();
	void	initialize(const char* av);
	char*	get_str();
	void	skip_till_nl();
	char*	getDbEntry(DbsRec* rec, int idf);
public:
	Makdbs(int ac, const char** av);
	~Makdbs();
	void	resource(const char* arg, int i);
	void	mkdbs();
	void	mkidx();
	void	stamp21() {
	    DbsRec	rec21 = {magicver21, comment, 0};
	    fwrite(&rec21, sizeof(DbsRec), 1, fidx);	// header record
	}
	void	wrtgrp(const char* ps) {
	    fprintf(fgrp, "%8ld %u %s\n", ftell(fseq), (INT) recnbr, ps);
	}
};

static	void	usage(const char* fmt = 0, const char* arg = 0);
static	int	cmpkey(INT* a, INT* b);

static	int	defmolc = 0;
static	int	idfy = 0;
static	int	monit = 0;
static	int	ignoreamb = 0;
static	bool	cridxf = false;
static	int	dbsch = 0;
static	char	dbname[MAXL] = "";
static	const	char*	srcpath = 0;
static	const	char*	dstpath = "";

static void usage(const char* fmt, const char* arg)
{
	if (fmt && arg) fprintf(stderr, fmt, arg);
	fputs("makdbs Version 2.0.0: format a set of sequences to be quickly addressed\n", stderr);
	fputs("Usage:\n", stderr);
	fputs("\tmakdbs [-d[Name]] [-sSrcdir] [-pDstdir] source_files\n",
		stderr);
	fputs("Option:\n", stderr);
	fputs("\t-a\t: Ignore abmiguity code\n", stderr);
	fputs("\t-bC\t: C = [e]mbl|[g]enbank|[n]brf|[p]rodb|[s]wiss\n", stderr);
	fputs("\t-c\t: include comments/descriptions\n", stderr);
	fputs("\t-dName\t: Body of file names to be created\n", stderr);
	fputs("\t\tIf omitted, default name depending on DB type\n", stderr);
	fputs("\t\tIf Name is omitted, taken from 1st source file name\n", stderr);
	fputs("\t-sSrcdir\t: Directory where source files are located\n", stderr);
	fprintf(stderr, "\t\tIf omitted, srcdir <= {., %s, %s}\n",
		getenv(ALN_DBS), DBS_SDIR);
	fputs("\t-pdstdir\t: Directory where destination files are located\n", stderr);
	fprintf(stderr, "\t\tIf omitted, dstdir = current dir\n");
	fputs("\t\t*.ent, *.grp, *idx, (*.odr), and *.seq files ", stderr);
	fputs("are made in this directory\n", stderr);
	fputs("\t-K[D|P]\t: DNA or Protein sequence\n", stderr);
	exit(1);
}

int Makdbs::encode(int c)
{
	if (isalpha(c)) {
	    return (defcode->encode[toupper(c) - 'A'] - bias);
	} else switch (c) {
	    case _UNP:
	    case _TRM:	return (gap_code - bias);
	    default:	return (IGNORE);
	}
}

char* Makdbs::get_str()
{
	if (memmode != 2) {
	    char* ps = fgets(str, MAXL, fsrc);
	    halfway = ps && (strlen(ps) + 1) == MAXL;
	    return (ps);
	}
	if (!*mptr) return (0);
	char*	qs = str;
	char*	ts = qs + MAXL - 1;
	while (*mptr && qs < ts)
	    if ((*qs++ = *mptr++) == '\n') break;
	*qs = '\0';
	halfway = ts == qs;
	return (str);
}

void Makdbs::skip_till_nl()
{
	if (memmode != 2) {
	    int	c;
	    while ((c = fgetc(fsrc)) != EOF && c != '\n') ;
	} else {
	    while (*mptr && *mptr++ != '\n') ;
	}
}

int Makdbs::convert()
{
	int	n = 0;

	while (get_str()) {
	    if (db->is_DbEnd(str) || db->is_DbEntry(str)) break;
	    for (char* pc = str; *pc; pc++) {
		if (*pc == _COMM) break;
		int	c = encode(*pc);
		if (0 <= c && c < ceil &&
		!(ignoreamb && c == defcode->amb_code)) {
		    putc(c, fseq);
		    ++n;
		}
	    }
	}
	putc(SEQ_DELIM, fseq);
	return (n);
}

int Makdbs::convert2()
{
	int	n = 0;
	int	b = 0;

	while (get_str()) {
	    if (db->is_DbEnd(str) || db->is_DbEntry(str)) break;
	    for (char* pc = str; *pc; pc++) {
		if (*pc == _COMM) break;
		int	c = encode(*pc);
		if (0 <= c && c < ceil &&
		!(ignoreamb && c == defcode->amb_code)) {
		    if (n++ & 1)	putc(b + c, fseq);
		    else	b = c << 4;
		}
	    }
	}
	if (n & 1)	putc(b + SEQ_DELIM, fseq);
	else	putc(ddelim, fseq);
	return (n);
}

char* Makdbs::getDbEntry(DbsRec* rec, int idf)
{
	char*	ps = str;
	char*	sp;
	if (db->FormID == FASTA) {
	    char*	sb = 0;
	    sp = ps++;
	    while (*++sp && !isspace(*sp)) {
		if (*sp == '|') {
		    if (sp[1] == '\0' || isspace(sp[1])) {  // ...| ..
			*sp = ' ';
		    } else {
			sb = sp + 1;
		   }
		}
	    }
	    if (sb) ps = sb;
	} else if (db->EntLabel) {
	    sp = ps = cdr(str);
	    while (*sp && !isspace(*sp)) ++sp;
	}
	if (!comment) *sp = '\0';
	int	sl = strlen(ps);
	if (idf && sl > idf) sl = idf;
	if (ps[sl - 1] == '\n') ps[sl - 1] = '\0';	// remove CR
	fputs(ps, fent);
	fputc('\0', fent);
	return (ps);
}

void Makdbs::mkdbs()
{
	DbsRec	rec;
	while (get_str()) {
	    if (db->is_DbEntry(str)) {
newentry:
		rec.seqptr = ftell(fseq);
		rec.entptr = ftell(fent);
		char*	ps = getDbEntry(&rec, idfy);
		if (wordcmp(ps, prv) < 0) cridxf = true;
		car(prv, ps);
		if (halfway) skip_till_nl();
	    }
	    halfway = false;
	    if (db->is_DbOrigin(str)) {
		rec.seqlen = isDRNA(defmolc)? 
		    convert2(): convert();
		fwrite(&rec, sizeof(DbsRec), 1, fidx);
		++recnbr;
		if (db->is_DbEntry(str)) goto newentry;
	    }
	}
	if (fsrc) {fclose(fsrc); fsrc = 0;}
}

static	DbsRec*	rbuf;
static	char*	cbuf;

static int cmpkey(INT* a, INT* b)
{
	return(strcmp(cbuf + rbuf[*a].entptr, cbuf + rbuf[*b].entptr));
}

void Makdbs::mkidx()
{
	size_t	flen = ftell(fent);
	cbuf = new char[flen];
	rewind(fent);
	if (fread(cbuf, sizeof(char), flen, fent) != flen)
	    fatal("Corrupted entry file !\n");
	rbuf = new DbsRec[recnbr];
	rewind(fidx);
	if (fread(rbuf, sizeof(DbsRec), recnbr, fidx) != recnbr)
	    fatal("Corrupted index file !\n");
	INT*	order = new INT[recnbr];
	for (INT i = 0; i < recnbr; ++i) order[i] = i;
	qsort((UPTR) order, recnbr, sizeof(INT), (CMPF) cmpkey);
	FILE*	fodr = fopenpbe(dstpath, dbname, ODR_EXT, "w", 2);
	fwrite(order, sizeof(INT), recnbr, fodr);
	fclose(fodr);
	delete[] rbuf;
	delete[] cbuf;
	delete[] order;
}

void Makdbs::resource(const char* av, int i)
{
	if (memmode != 2) {
	    fsrc = fopenpbe(srcpath, av, 0, "r", -1);
	    if (!fsrc) usage(not_found, av);
	    SeqDb*	cdb = dbsch? setform(dbsch): whichdb(str, fsrc);
	    if (db != cdb)
		fatal("Inconsistent DB files %d-%d!\n",
		    db->FormID, cdb->FormID);
	    rewind(fsrc);
	}
	if (memmode == 1) {
	    delete[] msrc;
	    try {
		mptr = msrc = new char[fsize[i] + 1];
	    } catch (std::bad_alloc ba) {
		fatal("no memory !\n");
	    }
	    if (fread(msrc, sizeof(char), fsize[i], fsrc) != fsize[i])
		fatal("Fail to read %s!\n", av);
	    msrc[fsize[i]] = '\0';
	    mptr = msrc;
	    fclose(fsrc);
	    fsrc = 0;
	}
}

void Makdbs::initialize(const char* av)
{
	db = dbsch? setform(dbsch): whichdb(str, fsrc);
	if (!db) fatal("Bad DB sorce: %s/%s!\n", srcpath, av);
	if (!*dbname && db->DbName) strcpy(dbname, db->DbName);
	if (!*dbname) partfnam(dbname, av, "b");
	if (defmolc == UNKNOWN) {
	    Seq	sd;
	    char	str[MAXL] = "";
	    defmolc = sd.infermolc(fsrc, str);
	}
	if (!defcode) defcode = setSeqCode(0, defmolc);
	fseq = fopenpbe(dstpath, dbname, SEQ_EXT, "w", 2);
	fidx = fopenpbe(dstpath, dbname, IDX_EXT, "w+", 2);
	fgrp = fopenpbe(dstpath, dbname, GRP_EXT, "w", 2);
	fent = fopenpbe(dstpath, dbname, ENT_EXT, "w+", 2);
	putc(SEQ_DELIM, fseq);	// for compatibility with blast
}

void Makdbs::examem(int argc, const char** argv)
{
	memsize = 0;
	fsize = new size_t[argc];
	if (!srcpath) srcpath = "";
	for (int i = 0; i < argc; ++i) {
	    fsrc = fopenpbe(srcpath, argv[i], 0, "r", -1);
	    if (!fsrc)
		fsrc = fopenpbe(dstpath, argv[i], 0, "r", 1);
	    if (!fsrc) usage(not_found, argv[0]);
	    if (i == 0) initialize(*argv);
	    fseek(fsrc, 0L, SEEK_END);
	    fsize[i] = ftell(fsrc);
	    fclose(fsrc); fsrc = 0;
	    if (fsize[i] + 1 < max_memory) {
		try {
		    mptr = msrc = new char[fsize[i] + 1];
		    memsize += fsize[i];
		} catch (std::bad_alloc ba) {
		    memsize = memmode = 0;
		    delete[] fsize;
		    fsize = 0;
		    return;
		}
	    } else {
		    memsize = memmode = 0;
		    delete[] fsize;
		    fsize = 0;
		    return;
	    }
	    if (argc > 1) delete[] msrc;
	}
	if (argc > 1) {
	    try {
		mptr = msrc = new char[memsize + 1];
	    } catch (std::bad_alloc ba) {
		memsize = 0;
		memmode = 1;
		delete[] fsize;
		fsize = 0;
		return;
	    }
	}
	memsize = 0;
	for (int i = 0; i < argc; ++i) {
	    FILE*	fsrc = fopenpbe(srcpath, argv[i], 0, "r", -1);
	    if (!fsrc)
		fsrc = fopenpbe(dstpath, argv[i], 0, "r", 1);
	    if (fread(msrc + memsize, sizeof(char), fsize[i], fsrc) != fsize[i])
		usage("Fail to read %s!\n", argv[i]);
	    memsize += fsize[i];
	}
	msrc[memsize] = '\0';
	mptr = msrc;
	memmode = 2;
	return;
}

Makdbs::Makdbs(int argc, const char** argv)
{
	db = 0;
	recnbr = 0;
	defcode = 0;
	msrc = mptr = 0;
	fsrc = fgrp = fseq = fidx = fent = 0;
	mss = tss = 0;
	halfway = false;
	*str = *prv = '\0';
	examem(argc, argv);
	if (defmolc == PROTEIN)
	    bias = ALA - 1;
	else
	    bias = A - 1;
	ceil = defcode->max_code - bias;
}

Makdbs::~Makdbs()
{
	delete[] msrc;
	delete[] fsize;
	fclose(fgrp);
	fclose(fseq);
	fclose(fidx);
	fclose(fent);
}

int main(int argc, const char** argv)
{
	const	char* val = 0;
	while (--argc > 0 && **++argv == OPTCHAR) {
	  switch (argv[0][1]) {
	    case 'a': ignoreamb = 1; break;
	    case 'b': dbsch = atoi(getarg(argc, argv, true)); break;
	    case 'c': comment = true; break;
	    case 'd': strcpy(dbname, getarg(argc, argv)); break;
	    case 'i': idfy = atoi(getarg(argc, argv, true)); break;
	    case 'p':
		dstpath = getarg(argc, argv);
		if (!dstpath) dstpath = "";
		break;
	    case 's':
		srcpath = getarg(argc, argv);
		if (!srcpath) srcpath = "";
		break;
	    case 'O':
		val = getarg(argc, argv);
		monit = val? atoi(val): 1;
		break;
	    case 'K':
		val = getarg(argc, argv);
		if (val) defmolc = setdefmolc(*val);
		break;
	    default:	break;
	  }
	}
	if (argc < 1) usage();

	Makdbs	md(argc, argv);
	for (int i = 0; i < argc; ++i) {
	    md.resource(argv[i], i);
	    md.wrtgrp(argv[i]);
	    md.mkdbs();
	}
	md.stamp21();
	md.wrtgrp("E_O_F");
	if (cridxf) md.mkidx();
	return (0);
}
