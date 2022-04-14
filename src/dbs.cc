/******************************************************************************
*
*	Type definitions for various database formats
*
**	Osamu Gotoh, ph.D.	(-2001)
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

#include "seq.h"

static	int	cmpodr_r20(const int* a, const int* b);

DbsDt*	dbs_dt[MAX_DBS];

static	DbsDt*	defdbf;

SeqDb SeqDBs[] = {
	{
		GenBank,
		DNA,
		"GenBank",
		"LOCUS",
		"DEFINITION",
		"ACCESSION",
		"KEYWORDS",
		"SOURCE",
		"REFERENCE",
		"  AUTHORS",
		"  TITLE",
		"  JOURNAL",
		"COMMENT",
		"FEATURES",
		"ORIGIN",
		"//",
		0,
		"%9d ",
		6, 10, 1, 12
	},
	{
		EMBL,
		DNA, 
		"EMBL",
		"ID",
		"DE",
		"AC",
		"KW",
		"OS",
		"RN",
		"RA",
		"RT",
		"RL",
		"CC",
		"FT",
		"SQ",
		"//",
		0,
		"     ",
		6, 10, 1, 5
	},
	{
		Swiss,
		PROTEIN,
		"Swiss",
		"ID",
		"DE",
		"AC",
		"KW",
		"OS",
		"RN",
		"RA",
		"RT",
		"RL",
		"CC",
		"FT",
		"SQ",
		"//",
		0,
		"     ",
		6, 10, 1, 5
	},
	{
		TFDS,
		DNA, 
		"TFDS",
		"ID",
		"DE",
		"AC",
		"KW",
		"OS",
		"RN",
		"RA",
		"RT",
		"RL",
		"CC",
		"FT",
		"SQ",
		"//",
		0,
		"     ",
		6, 10, 1, 5
	},
	{
		NBRF,
		PROTEIN,
		"NBRF",
		"ENTRY",
		"TITLE",
		"ACCESSION",
		"KEYWORDS",
		"SOURCE",
		"REFERENCE",
		"   #Authors",
		"   #Title",
		"   #Journal",
		"COMMENT",
		"FEATURE",
		"SEQUENCE",
		"///",
		"		5	10	15	20	25	30\n",
		"%7d ",
		30, 1, 1, 18
	},
	{
		ProDB,
		PROTEIN,
		"ProDB",
		"CODE",
		"NAME",
		0,
		"KEYWORD",
		"SOURCE",
		0,
		"AUTHOR",
		"TITLE",
		"JOURNAL",
		"COMMENT",
		0,
		"SEQUENCE",
		"END",
		0,
		"	  ",
		6, 10, 1, 5
	},
	{
		FASTA,
		UNKNOWN,
		0,
		">",
		";d",
		";n",
		";k",
		";s",
		";r",
		";a",
		";t",
		";j",
		";c",
		";f",
		0,
		0,
		0,
		"",
		1, 60, 0, 0
	},
	{
		MSF,
		PROTEIN,
		0,
		"PileUp",
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		"//",
		0,
		0,
		"",
		5, 10, 2, 0
	},
	{
		PIR,
		PROTEIN,
		0,
		">P1;",
		"",
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		"",
		1, 60, 0, 0
	},
	{
		NEXUS,
		PROTEIN,
		0,
		"",
		"",
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		"",
		1, 0, 0, 0
	},
	{
		Bare,
		UNKNOWN,
		0,
		0,
		";d",
		";n",
		";k",
		";s",
		";r",
		";a",
		";t",
		";j",
		";c",
		";f",
		0,
		0,
		0,
		"", 
		1, 60, 0, 0
	}
};

int SeqDb::is_DbEntry(const char* str) const
{
	if (FormID == FASTA) 
	    return (*str == *EntLabel);
	else if (FormID == PIR)
	    return (!strncmp(str, EntLabel, 4));
	else if (EntLabel)
	    return (!wordcmp(str, EntLabel));
	else	return (0);
}

int SeqDb::is_DbEnd(const char* str) const
{
	return (*str == '/');
}

int SeqDb::is_DbOrigin(const char* str) const
{
	if (SeqLabel)
	    return (!wordcmp(str, SeqLabel));
	if (FormID == FASTA) 
	    return (*str == *EntLabel);
	if (FormID == PIR)
	    return (!strncmp(str, EntLabel, 4));
	return (0);
}

bool space_digit(const char* ps)
{
	while (int c = *ps++)
	    if (!(isspace(c) || isdigit(c))) return (false);
	return (true);
}

SeqDb* whichdb(const char* ps)
{
	SeqDb*	db = SeqDBs + FASTA;
	if (*ps == *db->EntLabel) return (db);	// FASTA
	db = SeqDBs + PIR;
	if (!strncmp(ps, db->EntLabel, 4)) return (db);	// PIR
	db = SeqDBs + MSF;
	if (!strncmp(ps, db->EntLabel, 6)) return (db);	// MSF
	for (db = SeqDBs; db->FormID < FASTA; db++)
	    if (db->is_DbEntry(ps)) break;
	if (db->FormID == FASTA) db = SeqDBs + Bare;	// non specific
	if (db->FormID == EMBL && !strcmp(ps + strlen(ps) - 4, "AA.\n"))
	    db = SeqDBs + Swiss;
	return (db);
}

long SeqDb::dbnextentry(FILE* fd, char* ps) const
{
	long	pos = ftell(fd);

	while ((ps = fgets(ps, MAXL, fd)) && !is_DbEntry(ps))
	    pos = ftell(fd);
	return (ps? pos: EOF);
}

SeqDb* setform(int c)
{
static	char frmmsg[] = 
	    "Form [B]are/[F]asta/[E]mbl/[S]wiss/[N]brf/[G]bk/[P]db :";
static	SeqDb* curform = SeqDBs + FASTA;

	if (c < 0) c = progetc(frmmsg);
	switch (tolower(c)) {
	    case 'b':	c = Bare; break;
	    case 'e':	c = EMBL; break;
	    case 'd':	c = NBRF; break;
	    case 'f':	c = FASTA; break;
	    case 'g':	c = GenBank; break;
	    case 'm':	c = MSF; break;
	    case 'n':	c = NEXUS; break;
	    case 'p':	c = PIR; OutPrm.asterisk = 1; break;
	    case 'q':	c = ProDB; break;
	    case 's':	c = Swiss; break;
	    case 't':	c = TFDS; break;
	    default:	return (curform);
	}
	return (curform = SeqDBs + c);
}

void setdefdbf(DbsDt* dbf)
{
	defdbf = dbf;
}

void DbsDt::readseq(const char* fn)
{
	fseek(fseq, 0L, SEEK_END);
	long	fp = ftell(fseq);
	dbsseq	= new CHAR[fp];
	if (!dbsseq) fatal("No memory to accomodate seqdb!\n");
	rewind(fseq);
	if (fread(dbsseq, sizeof(CHAR), fp, fseq) != (size_t) fp)
	    fatal("%s: Bad seq file !\n", fn);
	fclose(fseq);
	fseq = 0;
}

size_t DbsDt::readgrp(FILE* fgrp)
{
static	const	char	frmt[] = "%ld %ld %s";
	char	str[MAXL];
	Mfile	mfd(sizeof(DbsGrp));
	DbsGrp	grp;
	grplbl = new Strlist();

	size_t	ress = 0;
	while (fscanf(fgrp, frmt, &grp.seqptr, &grp.recnbr, str) == 3) {
	    grplbl->push(str);
	    mfd.write(&grp);
	    ress += grp.seqptr;
	}
	numgrp = (INT) mfd.size() - 1;
	dbsgrp = (DbsGrp*) mfd.flush();
	return (ress);
}

static int cmpodr_r20(const int* a, const int* b)
{
	return (wordcmp(defdbf->entname(*a), defdbf->entname(*b), ENTLEN));
}

void DbsDt::readidx20(FILE* fd, const char* fn)
{
	fseek(fd, 0L, SEEK_END);
	long	fp = ftell(fd);
	if (fp % sizeof(DbsRec20)) fatal("%s: Index file may be corrupted!\n", fn);
	numidx = (INT) (fp / sizeof(DbsRec20));
	DbsRec20*	reched20 = new DbsRec20[numidx];
	rewind(fd);
	if (fread(reched20, sizeof(DbsRec20), numidx, fd) != numidx)
	    fatal("%s: Bad index file!", fn);
	DbsRec20*	rec20 = reched20;
	DbsRec* rec = recidx = new DbsRec[numidx];
	entry = new char[numidx * (ENTLEN + 1)];
	recodr = 0;
	size_t	entidx = 0;
	bool	sorted = true;
	for (INT i = 0; i < numidx; ++i, ++rec, ++rec20) {
	    rec->seqptr = rec20->seqptr;
	    if (i && strncmp(rec20[-1].entry, rec20->entry, ENTLEN) > 0)
		sorted = false;
	    rec->seqlen = (INT) rec20->seqlen;
	    rec->entptr = entidx;
	    entry[entidx + ENTLEN] = '\0';
	    strncpy(entry + entidx, rec20->entry, ENTLEN);
	    entidx += strlen(entry + entidx) + 1;
	}
	delete[] reched20;
	if (sorted) return;
	recodr = new INT[numidx];
	for (INT i = 0; i < numidx; ++i) recodr[i] = i;
	qsort((UPTR) recodr, numidx, sizeof(INT), (CMPF) cmpodr_r20);
}

void DbsDt::readentry(FILE* fent, const char* fn)
{
	fseek(fent, 0L, SEEK_END);
	size_t	fp = ftell(fent);
	entry = new char[fp];
	rewind(fent);
	if (fread(entry, sizeof(char), fp, fent) != fp)
	    fatal("%s: Bad entry file!", fn);
}

void DbsDt::readodr(FILE* fodr, const char* fn)
{
	fseek(fodr, 0L, SEEK_END);
	long	fp = ftell(fodr);
	recodr = new INT[fp / sizeof(INT)];
	rewind(fodr);
	if (fread(recodr, sizeof(INT), numidx, fodr) != numidx)
	    fatal("%s: Bad order file!", fn);
}

DbsRec* DbsDt::readidx(FILE* fidx, const char* fn)
{
	fseek(fidx, 0L, SEEK_END);
	long	fp = ftell(fidx);
	if (fp % sizeof(DbsRec))
	    fatal("%s: Index file may be corrupted!\n", fn);
	numidx = (INT) (fp / sizeof(DbsRec));
	recidx = new DbsRec[numidx];
	rewind(fidx);
	if (fread(recidx, sizeof(DbsRec), numidx, fidx) != numidx)
	    fatal("%s: Index file may be corrupted!\n", fn);
	if (recidx[--numidx].seqptr != magicver21) {
	    delete[] recidx; recidx = 0;
	    return (0);
	}
	return (recidx);
}

DbsGrp* DbsDt::finddbsgrp(const char* name) const
{
	DbsGrp*	dgrp = dbsgrp;

	for (int i = 0; i < numgrp; ++i, ++dgrp)
	    if (!strcmp(name, fsrcname(i))) return (dgrp);
	return (0);
}

int DbsDt::guessmolc() const
{
	int	n = 0;
	int	i = 0;
const	CHAR*	ps = dbsseq;

	for ( ; i < NTESTC; ++i) {
	    int	c = dbsseq? *ps++: getc(fseq);
	    if (c < 0) break;
	    if (c < MAXCODE) ++n;
	}
	if (fseq) rewind(fseq);
	return (n < i * 3 / 4)? DNA: PROTEIN;
}

void DbsDt::clean()
{
	curdb = 0; pseq = 0; fseq = 0; numgrp = 0;
	dbsseq = 0; dbsgrp = 0; grplbl = 0;
	numidx = 0; comment = false;
	recidx = 0; recodr = 0; entry = 0;
	gsiidx = 0; gsipool = 0;
}

void DbsDt::prepare(size_t entry_space, size_t num, size_t seq_space, size_t gsi_space)
{
	dbsseq = new CHAR[seq_space];
	recidx = new DbsRec[numidx = num];
	entry = new char[entry_space];
	if (gsi_space) {
	    gsipool = new int[gsi_space + 1];
	    gsiidx = new int[num + 1];
	}
}

void DbsDt::prepare(Strlist& sname, int num, size_t space, size_t gsi_space)
{
	dbsseq = new CHAR[space];
	recidx = new DbsRec[numidx = num];
	entry = sname.squeeze();
	if (gsi_space) {
	    gsipool = new int[gsi_space + 1];
	    gsiidx = new int[num + 1];
	}
}

DbsDt::DbsDt(int c, int molc)
{
	clean();
	switch (tolower(c)) {
	    case 'e':	curdb = SeqDBs + EMBL; break;
	    case 'd':	curdb = SeqDBs + NBRF; break;
	    case 'g':	curdb = SeqDBs + GenBank; break;
	    case 'm':	curdb = SeqDBs + PIR; break;
	    case 'n':	curdb = SeqDBs + NEXUS; break;
	    case 'p':	curdb = SeqDBs + ProDB; break;
	    case 's':	curdb = SeqDBs + Swiss; break;
	    case 't':	curdb = SeqDBs + TFDS; break;
	    case 'f':	curdb = SeqDBs + FASTA; break;
	    default:	break;
	}
	if (curdb) curdb->defmolc = molc;
}

static	const	char*	dbstab[3] = {".", getenv(ALN_DBS), DBS_DIR};

char* path2dbf(char* str, const char* fn, const char* ext)
{
#if USE_ZLIB
	char	extgz[10];
	strcpy(extgz, ext);
	strcat(extgz, gz_ext);
#endif
	for (int i = 0; i < 3; ++i) {
	    const char*	path = dbstab[i];
	    if (!path) continue;
	    if (FILE* fd = fopenpbe(path, fn, ext, "r", -1, str)) {
		fclose(fd);
		return (str);
	    }
#if USE_ZLIB
	    if (FILE* fd = fopenpbe(path, fn, extgz, "r", -1, str)) {
		fclose(fd);
		return (str);
	    }
#endif
	}
	return (0);
}

DbsDt::DbsDt(const char* form) : dbsid(form)
{
static 	const	char	openmsg[] = 
	    "Form [F]asta/[E]mbl/[S]wiss/[N]brf/[G]bk/[P]db/[C]lose : ";
	char	str[LINE_MAX];
	char	buf[LINE_MAX];

	clean();
	if (form == 0) form = progets(buf, openmsg);
	if (strlen(form) <= 1) {
	  switch (tolower(*form)) {
	    case 'e':	curdb = SeqDBs + EMBL; break;
	    case 'd':	curdb = SeqDBs + NBRF; break;
	    case 'g':	curdb = SeqDBs + GenBank; break;
	    case 'm':	curdb = SeqDBs + PIR; break;
	    case 'n':	curdb = SeqDBs + NEXUS; break;
	    case 'p':	curdb = SeqDBs + ProDB; break;
	    case 's':	curdb = SeqDBs + Swiss; break;
	    case 't':	curdb = SeqDBs + TFDS; break;
	    case 'f':
	    default:	curdb = SeqDBs + FASTA; break;
	  }
	} else		curdb = SeqDBs + FASTA;

	if (curdb->DbName) form = curdb->DbName;
	for (int i = 0; i < 3; ++i) {
	    const	char*	path = dbstab[i];
	    if (!path) continue;
	    FILE* fd = fopenpbe(path, form, IDX_EXT, "r", -1, str);
	    if (fd) {
// read "index" file
		DbsRec*	ridx = readidx(fd, str);
		fclose(fd);
		if (ridx) {
// read "entry" and "order" files
		    fd = fopenpbe(path, form, ENT_EXT, "r", 1, str);
		    readentry(fd, str);
		    fclose(fd);
		    fd = fopenpbe(path, form, ODR_EXT, "r", -1, str);
		    if (fd) {
			readodr(fd, str);
			fclose(fd);
		    }
		} else {
// read "header" file
		    setdefdbf(this);
		    fd = fopenpbe(path, form, HED_EXT, "r", -1, str);
		    if (!fd) fd = fopenpbe(path, form, IDX_EXT, "r", 1, str);
		    readidx20(fd, str);
		    fclose(fd);
		}
// read "group" file
		fd = fopenpbe(path, form, GRP_EXT, "r", 1, str);
#if USE_ZLIB
		size_t	rss = readgrp(fd);
#endif
		fclose(fd);
// assign "seq" file
 		if ((fseq = fopenpbe(path, form, SEQ_EXT, "r", -1, str))) {
		    pseq = strrealloc(0, str);
		    if (algmode.dim) readseq(str);
		} else {
#if USE_ZLIB
		    gzFile	gzfd = gzopenpbe(path, form, SGZ_EXT, "r", 1, str);
		    if (!gzfd) continue;
		    pseq = strrealloc(0, str);
		    dbsseq = new CHAR[rss];
		    if (!dbsseq) fatal(no_space);
		    if (fread(dbsseq, sizeof(CHAR), rss, gzfd) <= 0)
			fatal("Fail to read .seq.gz file !\n");
		    fclose(gzfd);
#else
		    continue;
#endif
		}
		if (curdb->defmolc == UNKNOWN)
		    curdb->defmolc = guessmolc();
		return;
	    }
	}
	fatal("%s was not found !\n", str);
}

DbsDt::~DbsDt()
{
	if (fseq) fclose(fseq);
	delete[] pseq; delete[] dbsseq;
	delete[] dbsgrp; delete grplbl;
	delete[] recidx;
	delete[] recodr; delete[] entry;
	delete[] gsiidx; delete[] gsipool;
}

template <typename file_t>
CHAR* Seq::read_dbres(file_t fd, const RANGE* rng)
{
	int	c;
	CHAR*	ps = at(0);
	long	fbase = ftell(fd);

	if (inex.molc == PROTEIN) {
	    int	bias = ALA - 1;
	    for ( ; neorng(rng); ++rng) {
		fseek(fd, rng->left + fbase, 0); // move to start site
		int	n = rng->right - rng->left;
		while (n-- && (c = getc(fd)) != EOF) {
		    *ps = c? c + bias: AMB;
		    if (isAmb(*ps++)) inex.ambs = 1;
		}
	    }
	} else {
	    int	bias = A - 1;
	    for ( ; neorng(rng); ++rng) {
		fseek(fd, rng->left / 2 + fbase, 0);	// move to start site
		int	n = rng->right - rng->left;
		if (n && rng->left % 2) {	// odd start
		    c = getc(fd);
		    int	h = c & 15;
		    *ps = h? h + bias: N;
		    if (isAmb(*ps++)) inex.ambs = 1;
		    --n;
		}
		while (n-- && (c = getc(fd)) != EOF) {
		    int	h = (c >> 4) & 15;
		    *ps = h? h + bias: N;
		    if (isAmb(*ps++)) inex.ambs = 1;
		    if (!n--) break;
		    h = c & 15;
		    *ps = h? h + bias: N;
		    if (isAmb(*ps++)) inex.ambs = 1;
		}
	    }
	}
	return (ps);
}

CHAR* Seq::read_dbres(CHAR* dbs, const RANGE* rng)
{
	CHAR*	ps = at(0);

	if (inex.molc == PROTEIN) {
	    int	bias = ALA - 1;
	    for ( ; neorng(rng); ++rng) {
		dbs += rng->left; 	// set the start site
		int	n = rng->right - rng->left;
		for ( ; n--; ++dbs) {
		    *ps = *dbs? *dbs + bias: AMB;
		    if (isAmb(*ps++)) inex.ambs = 1;
		}
	    }
	} else {
	    int	bias = A - 1;
	    for ( ; neorng(rng); ++rng) {
		dbs += rng->left / 2;	// set the start site
		int	n = rng->right - rng->left;
		if (n && rng->left % 2) {	// odd start
		    int	h = *dbs & 15;
		    *ps = h? h + bias: N;
		    if (isAmb(*ps++)) inex.ambs = 1;
		    --n;
		}
		for ( ; n--; ++dbs) {
		    int	h = (*dbs >> 4) & 15;
		    *ps = h? h + bias: N;
		    if (isAmb(*ps++)) inex.ambs = 1;
		    if (!n--) break;
		    h = *dbs & 15;
		    *ps = h? h + bias: N;
		    if (isAmb(*ps++)) inex.ambs = 1;
		}
	    }
	}
	return (ps);
}

Seq* Seq::read_dbseq(DbsDt* dbf, DbsRec* rec, const RANGE* rng)
{
	if (dbf->fseq && fseek(dbf->fseq, rec->seqptr, 0) == ERROR) return(0);
	int	n = 0;
	int	area = rng->left;
const 	RANGE*	r = rng;

	for ( ; neorng(r); ++r) n += r->right - r->left;
	left = 0;
	right = n;
	refresh(1, n);
	base_ = area;
	area = area_ - many;
	if (n > area) {
	    prompt("%s is truncated %d --> %d\n", sqname(), n, area);
	    n = area;
	}
	setSeqCode(this, dbf->curdb->defmolc);
	CHAR*	ps = dbf->fseq? read_dbres(dbf->fseq, rng):
		read_dbres(dbf->dbseq(rec), rng);
	postseq(ps);
	sname->assign(dbf->entname(rec));
	did = dbf->recno(rec);
	if (int pfqnum = dbf->gsisize(did)) {
	    delete sigII;
	    sigII = new SigII(dbf->gsient(did), pfqnum, isprotein()? 3: 1);
	    sigII->resetend(len);
	}
	return this;
}

Seq* Seq::read_dbseq(DbsDt* dbf, long pos)
{
	DbsRec*	rec = dbf->dbsrec(pos);

	if (!rec) return (0);
	RANGE	rng[2] = {{0, int(rec->seqlen)}, endrng};
	return (read_dbseq(dbf, rec, rng));
}

DbsRec* DbsDt::bisearch(const char* key) const
{
	int	left = 0;
	int	right = numidx;
	while (right - left >= 0) {
	    int	mid = (right + left) / 2;
	    DbsRec*	found = recidx + mid;
	    int comp = wordcmp(key, entname(recodr? dbsrec(found): found));
	    if (comp < 0)	right = mid - 1;
	    else if (comp > 0)	left = mid + 1;
	    else	return (found);
	}
	return (0);
}

DbsRec* DbsDt::findcode(const char* code) const
{
	DbsRec*	rec = bisearch(code);
	return recodr? dbsrec(rec): rec;
}

Seq* Seq::getdbseq(DbsDt* dbf, const char* code, int c, bool readin)
{
	int	negative = 0;
	RANGE	retained[MAXCR];
	RANGE*	rng = retained;
	RANGE*	maxrng = rng + MAXCR - 1;

	if (!dbf) dbf = defdbf;
	if (!dbf) return (0);
	refresh();
	DbsRec*	record = (0 <= c && c < int(dbf->numidx))? dbf->dbsrec(c): 0;
	if (!record) {
	    if (!code) return (0);
	    if (*code == DBSID) ++code;
	    char	token[MAXL];
	    car(token, code);
	    if (sname)	sname->assign(token);
	    else	sname = new Strlist(token);
	    record = dbf->findcode(token);
	    if (!record) return(0);
	}
	int 	ir = 0;
	char*	term = cdr(code);
	for ( ; term && *term && rng < maxrng; term = cdr(term)) {
	    if (isdigit(*term)) {
		INT	n = atoi(term);
		if (n > record->seqlen) n = record->seqlen;
		if (ir++ % 2) {
		    if (n < INT(rng->left)) ++negative;
		    (rng++)->right = n;
		} else	rng->left = n - 1;
		if (n == record->seqlen) break;
	    }
	}
	if (negative) {
	    fprintf(stderr, "Bad ranges: %s\n", code);
	    return (0);
	}
	if (ir == 0) {	// range not specified
	    rng->left = 0;
	    ++ir;
	} 
	if (ir % 2) (rng++)->right = record->seqlen;
	*rng = endrng;
	if (readin) {
	    return (read_dbseq(dbf, record, retained)?
		attrseq(cdr(code)): 0);
	}
	len = 0;
	for (rng = retained; neorng(rng); ++rng)
	    len += rng->right - rng->left;
	left = retained->left;
	right = (--rng)->right;
	return (this);
}

void EraDbsDt()
{
	for (int n = 0; n < MAX_DBS; ++n) delete dbs_dt[n];
}

