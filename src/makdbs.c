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
static	void	usage(const char* fmt = 0, const char* arg = 0);
static	int	cmpkey(INT* a, INT* b);

static	int	defmolc = 0;
static	int	idfy = 0;
static	int	monit = 0;
static	int	ignoreamb = 0;
static	int	dbsch = 0;
static	char	dbname[MAXL] = "";
static	const	char*	srcpath = 0;
static	const	char*	dstpath = "";
static	bool	gzout = false;

class Makdbs {
	int	molc;
	int	bias;
	int	ceil;
	bool	cridxf;
	SeqDb*	db;
	SEQ_CODE*	defcode;
	size_t	recnbr;
	FILE*	fgrp;
	FILE*	fseq;
	FILE*	fidx;
	FILE*	fent;
	char	str[MAXL];
	char	prv[MAXL];
	bool	halfway;
#if USE_ZLIB
	gzFile	gzseq;
#endif
	int	encode(int c);
	char*	getDbEntry(DbsRec* rec, int idf);
template <typename file_t, typename ofile_t>
	int	convert(file_t fsrc, ofile_t);
template <typename file_t, typename ofile_t>
	int	convert2(file_t fsrc, ofile_t);
	void	initialize(const char* av);
template <typename file_t>
	char*	get_str(file_t fsrc);
template <typename file_t>
	void	skip_till_nl(file_t fsrc);
template <typename file_t, typename ofile_t>
	void	makdbs(file_t fsrc, ofile_t);

public:
	Makdbs(int ac, const char** av, int mlc);
	~Makdbs();
	void	makdbs(const char* fn);
	void	mkidx();
	void	stamp21() {
	    DbsRec	rec21 = {magicver21, comment, 0};
	    fwrite(&rec21, sizeof(DbsRec), 1, fidx);	// header record
	}
	void	wrtgrp(const char* ps) {
	    long	fpos = 0L;
	    if (fseq)	fpos = ftell(fseq);
#if USE_ZLIB
	    else	fpos = ftell(gzseq);
#endif
	    fprintf(fgrp, "%8ld %u %s\n", fpos, (INT) recnbr, ps);
	}
};

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

template <typename file_t>
char* Makdbs::get_str(file_t fsrc)
{
	char* ps = fgets(str, MAXL, fsrc);
	halfway = ps && (strlen(ps) + 1) == MAXL;
	return (ps);
}

template <typename file_t>
void Makdbs::skip_till_nl(file_t fsrc)
{
	int	c;
	while ((c = fgetc(fsrc)) != EOF && c != '\n') ;
}

template <typename file_t, typename ofile_t>
int Makdbs::convert(file_t fsrc, ofile_t f_seq)
{
	int	n = 0;

	while (get_str(fsrc)) {
	    if (db->is_DbEnd(str) || db->is_DbEntry(str)) break;
	    for (char* pc = str; *pc; pc++) {
		if (*pc == _COMM) break;
		int	c = encode(*pc);
		if (0 <= c && c < ceil &&
		!(ignoreamb && c == defcode->amb_code)) {
		    fputc(c, f_seq);
		    ++n;
		}
	    }
	}
	fputc(SEQ_DELIM, f_seq);
	return (n);
}

template <typename file_t, typename ofile_t>
int Makdbs::convert2(file_t fsrc, ofile_t f_seq)
{
	int	n = 0;
	int	b = 0;

	while (get_str(fsrc)) {
	    if (db->is_DbEnd(str) || db->is_DbEntry(str)) break;
	    for (char* pc = str; *pc; pc++) {
		if (*pc == _COMM) break;
		int	c = encode(*pc);
		if (0 <= c && c < ceil &&
		!(ignoreamb && c == defcode->amb_code)) {
		    if (n++ & 1)	fputc(b + c, f_seq);
		    else	b = c << 4;
		}
	    }
	}
	if (n & 1)	fputc(b + SEQ_DELIM, f_seq);
	else	fputc(ddelim, f_seq);
	return (n);
}

template <typename file_t, typename ofile_t>
void Makdbs::makdbs(file_t fsrc, ofile_t f_seq)
{
	DbsRec	rec;
	while (get_str(fsrc)) {
	    if (db->is_DbEntry(str)) {
newentry:
		rec.seqptr = ftell(f_seq);
		rec.entptr = ftell(fent);
		char*	ps = getDbEntry(&rec, idfy);
		if (wordcmp(ps, prv) < 0) cridxf = true;
		car(prv, ps);
		if (halfway) skip_till_nl(fsrc);
	    }
	    halfway = false;
	    if (db->is_DbOrigin(str)) {
		rec.seqlen = isDRNA(molc)? 
		    convert2(fsrc, f_seq): convert(fsrc, f_seq);
		fwrite(&rec, sizeof(DbsRec), 1, fidx);
		++recnbr;
		if (db->is_DbEntry(str)) goto newentry;
	    }
	}
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

static	DbsRec*	rbuf;
static	char*	cbuf;

static int cmpkey(INT* a, INT* b)
{
	return(strcmp(cbuf + rbuf[*a].entptr, cbuf + rbuf[*b].entptr));
}

void Makdbs::mkidx()
{
	if (!cridxf) return;		// has been sorted
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

void Makdbs::makdbs(const char* av)
{
	if (is_gz(av)) {
#if USE_ZLIB
	    gzFile	gzfd = gzopenpbe(srcpath, av, 0, "r", -1);
	    if (!gzfd) usage(not_found, av);
	    if (fseq)	makdbs(gzfd, fseq);
	    else	makdbs(gzfd, gzseq); 
	    fclose(gzfd);
#else
	    fatal(gz_unsupport, av);
#endif
	} else {
	    FILE* fsrc = fopenpbe(srcpath, av, 0, "r", -1);
	    if (!fsrc) usage(not_found, av);
	    makdbs(fsrc, fseq); 
	    fclose(fsrc);
	}
}

Makdbs::Makdbs(int argc, const char** argv, int mlc) 
	: molc(mlc), bias(0), cridxf(false), db(0), recnbr(0),
	  fgrp(0), fseq(0), fidx(0), fent(0), halfway(false)
#if USE_ZLIB
	 , gzseq(0)
#endif
{
	*str = *prv = '\0';
	if (dbsch) db = setform(dbsch);
	else {
	    if (is_gz(*argv)) {
#if USE_ZLIB
		gzFile	gzfd = gzopenpbe(srcpath, *argv, 0, "r", -1);
		if (!gzfd) usage(not_found, *argv);
		db = whichdb(dbs_header(str, gzfd));
		fclose(gzfd);
#else
		fatal(gz_unsupport, *argv);
#endif
	    } else {
		FILE* fsrc = fopenpbe(srcpath, *argv, 0, "r", -1);
		if (!fsrc) usage(not_found, *argv);
		db = whichdb(dbs_header(str, fsrc));
	        fclose(fsrc);
	    }
	}
	if (!db) fatal("Bad DB sorce: %s/%s!\n", srcpath, *argv);
	if (!molc) molc = infermolc(*argv);
	bias = A - 1;
	if (molc == PROTEIN) bias = ALA - 1;
	defcode = setSeqCode(0, molc);
	*str = '\0';
	ceil = defcode->max_code - bias;

	if (!*dbname && db->DbName) strcpy(dbname, db->DbName);
	if (!*dbname) partfnam(dbname, *argv, "b");
#if USE_ZLIB
	if (gzout) {
	    gzseq = gzopenpbe(dstpath, dbname, SGZ_EXT, "w", 2);
	    fputc(SEQ_DELIM, gzseq);	// for compatibility with blast
	} else
#endif
	{
	    fseq = fopenpbe(dstpath, dbname, SEQ_EXT, "w", 2);
	    fputc(SEQ_DELIM, fseq);	// for compatibility with blast
	}
	fidx = fopenpbe(dstpath, dbname, IDX_EXT, "w+", 2);
	fgrp = fopenpbe(dstpath, dbname, GRP_EXT, "w", 2);
	fent = fopenpbe(dstpath, dbname, ENT_EXT, "w+", 2);
}

Makdbs::~Makdbs()
{
	fclose(fgrp);
	if (fseq) fclose(fseq);
#if USE_ZLIB
	if (gzseq) fclose(gzseq);
#endif
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
	    case 'g': gzout = 1; break;
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

	Makdbs	md(argc, argv, defmolc);
	for (int i = 0; i < argc; ++i) {
	    md.wrtgrp(argv[i]);
	    md.makdbs(argv[i]);
	}
	md.stamp21();
	md.wrtgrp("E_O_F");
	md.mkidx();
	return (0);
}
