/*****************************************************************************
*
*	biological sequences
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

#include <new>
#include "seq.h"

DefSetup	def_setup = {UNKNOWN, 0, IM_NONE};
RANGE	fullrng = {0, INT_MAX};
int	noseq = 0;

/*		     -   - A C M G R S  V T W Y  H K  D  B  N 	*/
CHAR ncredctab[]  = {15,15,0,1,4,2,5,6,10,3,7,8,10,9,12,13,14};
CHAR nucl2rnucl[] = { 5, 5,0,1,4,2,4,4, 4,3,4,4, 4,4, 4, 4, 4};
CHAR ncelements[] = { 0, 0,0,1,2,2,0,2, 0,3,3,3, 1,1, 2, 3, 0};
CHAR nccmpctab[]  = { 0, 1,2,3,6,4,6,6, 6,5,6,6, 6,6, 6, 6, 6};
/*	 _ _ X A R N D C Q E  G  H  I  L  K  M  F  P  S  T  W  Y  V B Z	*/
CHAR aacmpctab[]  = 
	{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,5,8};
CHAR aaredctab[] = 
     {21,21,20,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,2,5};
/*	{_,_,X,C,G,A,A,G,A,A,G,A,T,T,A,T,T,C,C,C,G,A,T,G,G,A};*/
CHAR tnredctab[] = 
	{4,4,4,1,2,0,0,2,0,0,2,0,3,3,0,3,3,1,1,1,2,0,3,2,2,0};
CHAR trelements[] = 
	{4,4,0,1,2,0,0,2,0,0,2,0,3,3,0,3,3,1,1,1,2,0,3,2,2,0};
CHAR nccode[26] = 
	{A,B,C,D,Z,Z,G,H,Z,Z,K,Z,M,N,Z,Z,Z,R,S,T,U,V,W,N,Y,Z};
CHAR aacode[26] = 
	{ALA,ASX,CYS,ASP,GLU,PHE,GLY,HIS,ILE,ZZZ,LYS,LEU,MET,
	 ASN,AMB,PRO,GLN,ARG,SER,THR,SEC,VAL,TRP,AMB,TYR,GLX};
CHAR trccode[26] =
	{ALA,ASX,CYS,ASP,GLU,PHE,GLY,HIS,ILE,SER2,LYS,LEU,MET,
	 ASN,TRM,PRO,GLN,ARG,SER,THR,TRM2,VAL,TRP,AMB,TYR,GLX};
CHAR rdcode[26] =
	{A,N,C,N,Z,Z,G,N,Z,Z,N,Z,N,N,Z,Z,Z,N,N,T,U,N,N,N,N,Z};

char nucl[] = "--ACMGRSVTWYHKDBN";
char nucr[] = "--ACMGRSVUWYHKDBN";
char amino[] =  "--XARNDCQEGHILKMFPSTWYVBUO";
char acodon[] = "--XARNDCQEGHILKMFPSTWYVJUO";
char ncodon[] = "--NCGAAGAAGATTATTCCCGATGGA";
char alphab[] = "--ABCDEFGHIJKLMNOPQRSTUVWXYZ";
char amino3[][4] = {
	"---", "---", "???", "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", 
	"Glu", "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe", "Pro", 
	"Ser", "Thr", "Trp", "Tyr", "Val", "Asx", "Sec", "***"};
char rdnucl[] = "--ACGTN";
CHAR red2nuc[] = {A,C,G,T};
char *Nucl = rdnucl + 2;
char *Amin = amino + 3;
char *Tron = acodon + 3;

static	char	seqdfn[MAXL] = "";

CHAR complcod[] = {___,_,T,G,K,C,Y,S,B,A,W,R,D,M,H,V,N};
/*		  -  - X A R N D C Q E G H I L K M F P S T W Y V S * */
CHAR aa2nuc[] = {___,_,N,C,G,A,A,G,A,A,G,A,T,T,A,T,T,C,C,C,G,A,T,G,G,A};

/* max_code, amb_code, base_code, ceil_code, gap_prof, encode, decode, redctab */
static	SEQ_CODE str_code = {28, 1, 2, 28, 28, 0, alphab, 0};
static	SEQ_CODE nts_code = {NSIMD, N, A, N, NSIMD, nccode, nucl, ncredctab};
static	SEQ_CODE trc_code = {TSIMD, AMB, ALA, TSIMD, TSIMD, trccode, acodon, 0};
static	SEQ_CODE aas_code = {ASIMD, AMB, ALA, ASIMD, ASIMD, aacode, amino, aaredctab};
static	const	INEX	def_inex = {0};

void setdfn(const char* newdfn) {topath(seqdfn, newdfn);}

int setdefmolc(const int& molc)
{
	switch (molc) {
	    case 'X': case 'x': def_setup.delamb = 1;
	    case PROTEIN: case 'A': case 'a': case 'P': case 'p':
		def_setup.defmolc = PROTEIN;
		break;
	    case 'N': case 'n':	def_setup.delamb = 1;
	    case DNA: case 'D': case 'd':
		def_setup.defmolc = DNA;
		break;
	    case RNA: case 'R': case 'r':
		def_setup.defmolc = RNA;
		break;
	    case TRON: case 'T': case 't':
		def_setup.defmolc = TRON;
		break;
	    case GENOME: case 'G': case 'g':
		def_setup.defmolc = GENOME;
		break;
	    case 'M': case 'm':
		def_setup.def_input_mode = IM_MULT;
		break;
	    case 'S': case 's':
		def_setup.def_input_mode = IM_SNGL;
		break;
	    case UNKNOWN:
		def_setup.defmolc = UNKNOWN;
		break;
	    default: break;
	}
	return (def_setup.defmolc);
}

InputMode get_def_input_mode() {return (def_setup.def_input_mode);}

static bool isattrib(const char* s)
{
	if (strchr(SensChar, *s) || *s == ',' || 
	    *s == _PROF || *s == _TRM) return (true);
	if (*s == '-' || *s == '+') ++s;
	while (isdigit(*s)) ++s;
	return (!*s || isspace(*s));
}

SEQ_CODE* setSeqCode(Seq* sd, const int& molc)
{
	SEQ_CODE* sqcode;

	switch (molc) {
	    case DNA: case RNA:	sqcode = &nts_code; break;
	    case PROTEIN:	sqcode = &aas_code; break;
	    case TRON:		sqcode = &trc_code; break;
	    case GENOME:	sqcode = &nts_code; break;
	    case UNKNOWN:
	    default:		sqcode = &str_code; break;
	}
	if (sd) {
	    sd->code = sqcode;
	    sd->inex.molc = molc;
	}
	return (sqcode);
}

/* examine and get the first query sequence */

SeqServer::SeqServer(int ac, const char** av, InputMode infm, 
	const char* catalog, int mq, int mt)
	: argc(ac), argc0(ac), argv(av), argv0(av), 
	  cfrom(0), cto(0), input_form(infm)
{
	molc[0] = mq;
	molc[1] = mt == TRON? DNA: mt;
	input_ns = (infm == IM_SNGL)? 1: 2;
	if (algmode.blk) {
	    target_dbf = dbs_dt[0];
	    query_dbf = dbs_dt[1];
	} else {
	    if (dbs_dt[1]) {
		target_dbf = dbs_dt[0];
		query_dbf = dbs_dt[1];
	    } else
		target_dbf = query_dbf = dbs_dt[0];
	}
	vclear(fd, 2);
#if USE_ZLIB
	vclear(gzfd, 2);
#endif
	fc = catalog? fopen(catalog, "r"): 0;
	nfrom = counter = 0;
	nto = INT_MAX;
	attr[0] = attr[1] = 0;
	atsz[0] = atsz[1] = 0;
	sw = false;
}

void SeqServer::reset()
{
	argc = argc0;
	argv = argv0;
	if (fd[0]) {fclose(fd[0]); fd[0] = 0;}
	if (fd[1]) {fclose(fd[1]); fd[1] = 0;}
#if USE_ZLIB
	if (gzfd[0]) {fclose(gzfd[0]); gzfd[0] = 0;}
	if (gzfd[1]) {fclose(gzfd[1]); gzfd[1] = 0;}
#endif
	if (fc) rewind(fc);
	nfrom = counter = 0;
	nto = INT_MAX;
	delete[] cfrom; cfrom = 0;
	delete[] cto; cto = 0;
	sw = false;
}

InSt SeqServer::nextseq(Seq* sd, int which)
{
	char	str[MAXL];

// read from seq listed in catalog
	if (fc) {
	    if (fgets(str, MAXL, fc)) {
		if (!sd->getseq(str, target_dbf)) return IS_ERR;
		return IS_OK;
	    }
	    return IS_END;
	}

	while (true) {
#if USE_ZLIB
	  if (!fd[which] && !gzfd[which]) 
#else
	  if (!fd[which])
#endif
	  {
	    const	char*	fn = *argv++;
	    if (argc-- <= 0) return IS_END;

// read from database
	    if (*fn == DBSID) {
		DbsDt*	dbf = which == 0? query_dbf: target_dbf;
		if (!dbf || !sd->getdbseq(dbf, fn)) {
		    prompt("%s not in db!\n", fn + 1);
		    return IS_ERR;
		}
		return IS_OK;
	    }

// prepare to read file
	    strcpy(str, fn);
	    char*	ps = strchr(str, '(');
	    if (ps) {	/* exam. seqs in range */
		*ps = '\0';
		while (*++ps && isspace(*ps)) ;
		if (*ps == '$') ++ps;
		char*	qs = ps;
		int	c = 0;
		while ((c = *ps) && !isspace(c) && !strchr(",-)", c)) ++ps;
		if (c) *ps++ = '\0';
		if (isdigit(*qs))	nfrom = atoi(qs) - 1;
		else if (isalpha(*qs))	cfrom = strrealloc(0, qs);

		while (*ps && isspace(*ps)) ++ps;

		if (*ps == '$') ++ps;
		qs = ps;
		while (*ps && !isspace(*ps) && *ps != ')') ++ps;
		*ps = '\0';
		if (isdigit(*qs))	nto = atoi(qs);
		else if (isalpha(*qs))	cto = strrealloc(0, qs);
		else if (c != ',' && c != '-') {
		    if (cfrom) cto = strrealloc(0, "");
		    else nto = nfrom + 1;
		}
	    } else if ((ps = strchr(str, '$'))) {
		char* qs = ++ps;
		while (*ps && !isspace(*ps)) ++ps;
		*ps = '\0';
		if (isdigit(*qs)) {
		    nfrom = atoi(qs);
		    nto = nfrom + 1;
		} else if (isalpha(*qs)) {
		    cfrom = strrealloc(0, qs);
		    cto = strrealloc(0, "");
		}
	    }
#if USE_ZLIB
	    fd[which] = sd->openseq(str, gzfd + which);
	    if (!fd[which] && !gzfd[which]) fatal(not_found, str);
#else
	    fd[which] = sd->openseq(str);
	    if (!fd[which]) fatal(not_found, str);
#endif
	    int	attrsize = strlen(cdr(fn)) + 8;
	    if (attrsize > atsz[which]) {
		delete[] attr[which];
	  	attr[which] = new char[attrsize];
		atsz[which] = attrsize;
	    }
	    strcpy(attr[which], cdr(fn));
	    strcat(attr[which], " S");
	  }

// read from mulitple seq. files
 	  if (fd[which]) {
	   while (sd->fgetseq(fd[which], attr[which])) {
	    if (molc[which] == UNKNOWN) {
		molc[which] = sd->inex.molc;
		strcat(attr[which], molc[which] == PROTEIN? "P": "D");
	    }
	    if (sw && cto && (!*cto || (sd->sname && !wordcmp(cto, (*sd->sname)[0])))) {
		sw = false; delete[] cto; cto = 0;
	    }
	    if (cfrom && sd->sname && !wordcmp(cfrom, (*sd->sname)[0])) {
		sw = true; delete[] cfrom; cfrom = 0;
	    }
	    if (!cfrom && nfrom == counter) sw = true;
	    ++counter;
	    if (sw && counter <= nto) return IS_OK;
	    if (sw && counter > nto) return IS_END;
	    else	continue;
	   }
	   fclose(fd[which]);
	   fd[which] = 0;
	  }
#if USE_ZLIB
	  else if (gzfd[which]) {
	   while (sd->fgetseq(gzfd[which], attr[which])) {
	    if (molc[which] == UNKNOWN) {
		molc[which] = sd->inex.molc;
		strcat(attr[which], molc[which] == PROTEIN? "P": "D");
	    }
	    if (sw && cto && (!*cto || (sd->sname && !wordcmp(cto, (*sd->sname)[0])))) {
		sw = false; delete[] cto; cto = 0;
	    }
	    if (cfrom && sd->sname && !wordcmp(cfrom, (*sd->sname)[0])) {
		sw = true; delete[] cfrom; cfrom = 0;
	    }
	    if (!cfrom && nfrom == counter) sw = true;
	    ++counter;
	    if (sw && counter <= nto) return IS_OK;
	    if (sw && counter > nto) return IS_END;
	    else	continue;
	   }
	   fclose(gzfd[which]);
	   gzfd[which] = 0;
	  }
#endif
	  delete[] cfrom; cfrom = 0;
	  delete[] cto; cto = 0;
	  nfrom = 0; nto = INT_MAX;
	}
	return IS_END;
}

size_t SeqServer::total_seq_len(Seq* sd, int* many)
{
	reset();
	size_t	tsz = 0;
	InSt	ist;
	if (many) *many = 0;
	while ((ist = nextseq(sd)) != IS_END)
	    if (ist == IS_OK) {
		tsz += sd->len;
		if (many) ++(*many);
	    }
	reset();
	return (tsz);
}

void Seq::refresh(const int& num, const int& length)
{
	if (!num && !seq_) return;
	if (!is_eldest()) {
	    sid = vrtl;
	    seq_ = end_ = 0; cmps = 0; lens = 0; nbr = 0;
	    sname = 0; spath = 0; exons = 0; msaname = 0; descr = 0;
#if USE_WEIGHT
	    if (inex.vtwt) delete[] weight;
	    weight = pairwt = 0;
#endif
	} else {
	    delete sigII;
#if USE_WEIGHT
	    delete[] weight; weight = 0;
	    delete[] pairwt; pairwt = 0;
#endif
	    if (num != many) {
		delete[] lens; lens = 0;
		delete[] nbr; nbr = 0;
	    }
	    if (!num) {
		if (seq_) {delete[] (seq_ - many); seq_ = 0;}
		delete[] cmps; cmps = 0;
		delete[] exons; exons = 0;
		delete[] spath; spath = 0;
		delete[] msaname; msaname = 0;
		delete sname; sname = 0;
		delete descr; descr = 0;
	    }
	}
	delete[] jxt; jxt = 0;
	delete exin; exin = 0; sigII = 0; 
	left = right = base_ = 0;
	jscr = 0;  vrtl = 0; anti_ = 0;
	inex = def_inex;
	if (num) {
	    seqalloc(num, length);
	    if (!lens) lens = new int[num];
	    vclear(lens, num);
	    if (!nbr) nbr = new int[num];
	    vclear(nbr, num);
	    if (sname)	sname->reset(num);
	    else	sname = new Strlist(num);
	} else {
	    many = len = 0;
	}
}

// alias of this

Seq* Seq::aliaseq(Seq* dest, const bool& this_is_alias)
{
	if (!dest)	dest = new Seq(0);
	else	dest->refresh();
	int	destid = dest->sid;
	memcpy(dest, this, sizeof(Seq));
	dest->jxt = 0;
	if (is_eldest()) --vrtl;
	if (this_is_alias) {
	    vrtl = sid;
#if USE_WEIGHT
	    inex.vtwt = !weight;
	} else {
	    dest->vrtl = destid;
	    dest->inex.vtwt = !dest->weight;
#else
	} else {
	    dest->vrtl = destid;
#endif
	}
	return (dest);
}

Seq::~Seq()
{
	if (noseq > 0) --noseq;
	delete mnhash;
	refresh();
}

void initseq(Seq** seqs, int n)
{
	for ( ; n--; ++seqs) *seqs = new Seq();
}

// set eldest to the last member of aliases

void rectiseq(Seq** seqs, Seq** cur)
{
	for ( ; cur > seqs; --cur) {
	    if ((*cur)->is_eldest()) continue;
	    int	old_sid = (*cur)->sid;
	    Seq**	wrk = cur;
	    while (--wrk >= seqs) {
		if ((*wrk)->sid == old_sid && (*wrk)->is_eldest()) {
		    (*cur)->sid = (*cur)->vrtl;
		    (*cur)->vrtl = (*wrk)->vrtl;
		    (*wrk)->vrtl = (*wrk)->sid;
		    if ((*wrk)->lens != (*cur)->lens)
			(*cur)->vrtl = (*cur)->sid;
		    break;
		}
	    }
	    while (--wrk >= seqs) {
		if ((*wrk)->sid == old_sid)
		    (*wrk)->sid = (*cur)->sid;
	    }
	}
}

void clearseq(Seq** seqs, int n)
{
	rectiseq(seqs, seqs + n - 1);
	for ( ; n--; ++seqs) delete *seqs;
}

void cleanseq(Seq** seqs, int n, const int& num)
{
	rectiseq(seqs, seqs + n - 1);
	for ( ; n--; ++seqs) (*seqs)->refresh(num);
}

void reportseq(Seq** sqs, int n)
{
	for (int i = 0; i < n; ++i) {
	    Seq*&	sd = sqs[i];
	    printf("%d\t%d\t%d\t%d\t%d\t%d\t%p\t%p\n", i + 1,
	    sd->sid, sd->vrtl, sd->len, sd->area_, sd->did, sd->seq_, sd->lens);
	}
}

void Seq::fillpad()
{
	int	pad = code->amb_code;
	memset(seq_ - many, pad, many);
	memset(at(len), pad, many);
}

void swapseq(Seq** x, Seq** y)
{
	if (!x || !y) return;
	swap(*x, *y);
 }

void Seq::seqalloc(const int& num, const int& len, const bool& keep)
{
	int length = len? len: DEFSEQLEN;
	int	area = num * (length + 2);
	if (seq_ && area <= area_) {
	    if (num != many) {
		seq_ += num - many;
		end_ -= num - many;
		many = byte = num;
	    }
	    return;
	}

	CHAR*	ss = 0;
	try {
	    ss = new CHAR[area];
	} catch (std::bad_alloc ba) {
	    fatal(NoSeqSpace);
	}
	if (seq_) {
	    seq_ -= many;
	    if (keep) memcpy(ss, seq_, area_);
	    delete[] seq_;
	}
	many = byte = num;
	area_ = area;
	seq_ = ss + many;
	end_ = ss + area_ - many;
}

CHAR* Seq::seq_realloc()
{
	int	area = area_;
	area_ *= 2;
	CHAR*	ss = 0;

	try {
	    ss = new CHAR[area_];
	} catch (std::bad_alloc ba) {
	    fatal(NoSeqSpace);
	}
	memcpy(ss, seq_ -= many, area);
	delete[] seq_ ;
	int	sz = end_ - seq_;
	seq_ = ss + many;
	end_ = ss + area_ - many;
	return (ss + sz);
}

void Seq::initialize()
{
	sid = ++noseq; did = 0; vrtl = 0;
	seq_ = end_ = 0; area_ = 0;
	anti_ = 0; len = left = right = base_ = 0;
	CdsNo = tlen = wllvl = 0;
	lens = 0; nbr = 0; code = 0; spath = 0; msaname = 0; sname = 0; 
	descr = 0; cmps = 0; jxt = 0; exons = 0; exin = 0; sigII = 0; 
	inex = def_inex;
	mnhash = 0;
#if USE_WEIGHT
	weight = pairwt = 0;
#endif
}

Seq::Seq(const int& num, const int& length)
{
	initialize();
	if ((many = num)) {
	    byte = many = num;
	    lens = new int[many];
	    vclear(lens, many);
	    nbr = new int[many];
	    vclear(nbr, many);
	    sname = new Strlist(num);
	}
	if (num) seqalloc(num, length);
}

Seq::Seq(const char* fname)
{
	initialize();
	if (!fname) return;
	FILE*	fd = 0;
#if USE_ZLIB
	gzFile  gzfd = 0;
	fd = openseq(fname, &gzfd);
	if (gzfd) {
	    if (!fgetseq(gzfd, cdr(fname))) fatal("%s is empty !\n", fname);
	    fclose(gzfd);
	    return;
	}
#else
	fd = openseq(fname);
#endif
	if (!fd) fatal("%s not found !\n", fname);
	if (!fgetseq(fd, cdr(fname))) fatal("%s is empty !\n", fname);
	fclose(fd);
}

Seq::Seq(Seq& sd, const int* which, const int& snl)
{
	initialize();
	if (which) sd.extseq(this, which, snl);
	else	sd.copyseq(this, snl);
}

Seq::Seq(Seq* sd, const int* which, const int& snl)
{
	initialize();
	if (which) sd->extseq(this, which, snl);
	else	sd->copyseq(this, snl);
}

const char* Seq::path2fn(const char* pname) const
{
	const char*	ns = strrchr(pname, PATHDELM);
	ns = ns? ns + 1: pname;
	if (*ns == _CHEAD || *ns == _NHEAD) ++ns;
	return ns;
}

#if USE_ZLIB
FILE* Seq::openseq(const char* str, gzFile* gzfd)
#else
FILE* Seq::openseq(const char* str)
#endif
{
	char	qname[LINE_MAX];
	char	pname[LINE_MAX];

	car(qname, str);
	makefnam(qname, seqdfn, pname);
	if (is_gz(pname)) {
#if USE_ZLIB
	    *gzfd = gzopen(pname, "r");
	    if (*gzfd) spath = strrealloc(spath, pname);
	    return (0);
#else
	    fatal(gz_unsupport, pname);
#endif
	} 
	FILE*	fd = fopen(pname, "r");
	if (fd) spath = strrealloc(spath, pname);
	return (fd);
}

void Seq::comple()
{
	if (!seq_) return;
	delete exin; exin = 0;
	if (isprotein()) {
	    prompt("Can't complement Protein sequence !\n");
	    return;
	}
	if (inex.molc == TRON) {
	    tron2nuc(true);
	    nuc2tron();
	} else {
	    CHAR*	ss = at(0);
	    CHAR*	tt = at(len);

	    for ( ; ss < tt; ss++)
		*ss = complcod[*ss];
	}
	inex.sens ^= COMPLE;
}

void Seq::rev_attr()
{
static	const	int	a2t[4] = {0, 2, 1, 3};
	inex.sens ^= REVERS;
	inex.polA = a2t[inex.polA];
	INT t = inex.exgl;
	inex.exgl = inex.exgr;
	inex.exgr = t;
	swap(left, right);
	left = len - left;
	right = len - right;
}

void Seq::reverse()
{
	delete exin; exin = 0;
	if (!seq_) return;
	CHAR* seq = at(0);
	CHAR* rsq = at(len - 1);
	CHAR* temp = new CHAR[many];

	for ( ; seq < rsq; seq += many, rsq -= many) {
	    memcpy(temp, seq, many);
	    memcpy(seq, rsq, many);
	    memcpy(rsq, temp, many);
	}
	delete[] temp;
	rev_attr();
}

void Seq::comrev(Seq** dest)
{
	copyseq(*dest);
	(*dest)->comrev();
	anti_ = dest;
}

void antiseq(Seq** seqs)
{
	Seq**	as = (*seqs)->anti_;
	if (as) {
	    swapseq(seqs, as);
	    (*seqs)->anti_ = as;
	    (*as)->anti_ = seqs;
	} else (*seqs)->comrev();
}

JUXT* Seq::revjxt()
{
        JUXT*   wjx = jxt;
        JUXT*   lst = jxt + CdsNo;

        for ( ; wjx < lst; ++wjx) {
            wjx->jx = lst->jx - wjx->jx - wjx->jlen;
            wjx->jy = len - wjx->jy - wjx->jlen;
        }
        return (vreverse(jxt, CdsNo));
}

/* assume seq.many = 1 */

CHAR* spliceTron(CHAR* spliced, const CHAR* b5, const CHAR* b3, const int& n)
{
	int	nn = n + n;
	CHAR*	sp = spliced + nn;

	for (int i = 0; i++ < nn; )
	    *sp++ = aa2nuc[*b5++];
	for (int i = 0; i++ < nn; )
	    *sp++ = aa2nuc[*b3++];
	sp = spliced + nn;
	for (int i = 0; i < nn; sp += n)
	    spliced[i++] = nuc2tron3(sp, n);
	return (spliced);
}

void Seq::nuc2tron()
{
	if (inex.molc == TRON) return;

	CHAR*	que[2];
	int	parity = 0;

	CHAR*	buf = new CHAR[2 * many];
	que[0] = buf;
	que[1] = buf + many;
	memset(que[1], AMB, many);
	CHAR*	soc = at(0) - many;
	CHAR*	trm = at(len);
	for ( ; soc < trm; soc += many) {
	    for (int i = 0; i < many; ++i)
	    	que[parity][i] = nuc2tron3(soc + i, many);
	    parity = 1 - parity;
	    memcpy(soc, que[parity], many);
	}
	memset(soc, AMB, many);
	delete[] buf;
	setSeqCode(this, TRON);
	inex.trcv = 1;
	delete[] cmps; cmps = 0;
}

void Seq::tron2nuc(const bool& cmpl)
{
	if (!istron()) return;

	CHAR*	soc = at(0) - many;
	CHAR*	trm = at(len);
	memset(soc, NAMB, many);
	for (soc += many; soc < trm; ++soc) { 
	    *soc = aa2nuc[*soc];
	    if (cmpl) *soc = complcod[*soc];
	}
	memset(soc, NAMB, many);
	setSeqCode(this, DNA);
	inex.trcv = 0;
	delete[] cmps; cmps = 0;
}

RANGE* Seq::setrange(const char* attr, int* ncr) const
{
	if (!attr) return (0);
	int 	lr = 0;		// left boundary
	int	c;
const	char*	pa = attr;
	int	num = 0;

	while ((c = *pa++)) {
	    if (isdigit(c)) {
		int	n = atoi(pa - 1) - base_;
		if (lr % 2) ++num;
		if (lr++ == 0) {
		    left = n - 1;
		    if (left < 0 || (len && left >= len))
			left = 0;
		} else {
		    right = n;
		    if (right <= left || (len && right > len))
			right = len;
		}
		while (isdigit(*pa)) ++pa;
	    } else if (c == ',' && ++lr % 2) ++num;
	    else if (!isspace(c)) break;
	}
	if (lr % 2) ++num;
	if (ncr) *ncr = num;
	if (!ncr || !num) return (0);
	RANGE*	slices = new RANGE[num + 2];
	lr = 0;
	pa = attr;
	RANGE*	rng = slices;
	*rng = fullrng;		// default
	while ((c = *pa++)) {
	    if (isdigit(c)) {
		int	n = atoi(pa - 1) - base_;
		if (lr++ % 2 == 0) rng->left = n - 1;
		else	(rng++)->right = n;
		while (isdigit(*pa)) ++pa;
	    } else if (c == ',' && ++lr % 2) ++rng;
	    else if (!isspace(c)) break;
	}
	if (lr < 2) ++rng;
	*ncr = rng - slices;
	*rng = endrng;
	return (slices);
}

Seq* Seq::attrseq(const char* pa)
{
	if (!pa) return this;
	float	fdelg = 0.;
	while (int c = *pa++) {
	    if (isdigit(c)) {
		while (isdigit(*pa)) ++pa;
	    } else {
		switch (toupper(c)) {
		  case CM:	comple(); break;
		  case RV:	reverse(); break;
		  case CR:	comrev(); break;
		  case 'C':	inex.form = CIRCLE; break;
		  case 'L':	inex.form = LINEAR; break;
		  case 'G':	inex.intr = 1; break;
		  case _PROF:	inex.prof = 1; break;
		  case _TRM:	inex.trcv = 1; break;
		  case _DELG:
		    while (isspace(*pa)) ++pa;
		    if (isdigit(*pa)) {
			fdelg = atof(pa);
			while (isdigit(*++pa)) ;
		     } else
			fdelg = DelFrac;
		  default:  break;
		}
	    }
	}
	if (fdelg > 0.) elim_column(DEL_GRC, fdelg);
	if (inex.dela) elim_column(DEL_AMB, 0);
	if (inex.trcv) nuc2tron();
	return (this);
}

int Seq::isAmb(const CHAR& r) const
{
static	CHAR nucamb[] = {0,0,0,0,2,0,2,2,3,0,2,2,3,2,3,3,4};
static	CHAR aasamb[] = {0,0,20,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,2};

	switch (inex.molc) {
	    case PROTEIN: return aasamb[r];;
	    case DNA: case RNA:	case GENOME: return (nucamb[r]);
	    default:	return (0);
	}
}

//	reset terminal gap codes

int Seq::sname2memno(const char* memid)
{
	if (many == 1) return (0);
	if (!memid) return (-1);
	if (!mnhash) {
	    mnhash = new StrHash<int>(many);
	    for (int i = 0; i < many; ++i)
		mnhash->assign((*sname)[i], i);
	}
	KVpair<INT, int>* kv = mnhash->find(memid);
	return (kv? kv->val: -1);
}

void Seq::exg_seq(int gl, int gr) 
{
	CHAR*	ss = at(0);
	CHAR*	tt = at(len);
	inex.exgl = gl? 1: 0;
	inex.exgr = gr? 1: 0;

	gl = gl || alprm.tgapf < 1;
	gr = gr || alprm.tgapf < 1;
	CHAR	new_codel = (CHAR) (gl? nil_code: gap_code);
	CHAR	new_coder = (CHAR) (gr? nil_code: gap_code);
	if (gl || !algmode.qck) vset(ss - many, new_codel, many);
	if (gr || !algmode.qck) vset(tt, new_coder, many);
	inex.nils = gl || gr;
	for (int i = 0; i < many; ++i) {
	    CHAR*	s = ss + i;
	    for ( ; s < tt; s += many) {
		if (IsGap(*s)) *s = new_codel;
		else	break;
	    }
	    if (s >= tt && (gl || !gr)) continue;
	    for (s = tt - many + i; s > ss; s -= many) {
		if (IsGap(*s)) *s = new_coder;
		else	break;
	    }
	}
	if (anti_) {
	    (*anti_)->inex.exgl = inex.exgr;
	    (*anti_)->inex.exgr = inex.exgl;
	}
}

void Seq::test_gap_amb() 
{
	CHAR*	ss = at(left);
	CHAR*	tt = at(right);

	inex.dels = inex.ambs = 0;
	for ( ; ss < tt; ++ss) {
	    if (IsGap(*ss)) inex.dels = 1;
	    int	namb = isAmb(*ss);
	    if (!namb) continue;
	    if (isprotein()) {
		if (namb == 20) inex.ambs |= 1;
		else	inex.ambs |= 2;
	    } else	inex.ambs = 1;
	}
}

Seq* Seq::postseq(const CHAR* last)
{
	left = 0;
	len = tlen = right = many? (last - seq_ + many - 1) / many: 0;
	if (many && len) {
	    fillpad();
	    test_gap_amb();
	    if (many == 1 && !lens[0]) lens[0] = len;
	} else {
	    inex.dels = inex.ambs = 0;
	}
	return this;
}

static int randaa()
{
static const float pcmp[21] = {0.00000, 0.07686, 0.12792, 0.17047, 0.22173,
	0.24201, 0.28307, 0.34489, 0.41961, 0.44259, 0.49516, 0.58627, 0.64577,
	0.66918, 0.70971, 0.76025, 0.82847, 0.88699 ,0.90132, 0.93363, 1.00000};
	double	rn = drand48();
	int	a = int(rn * 20);
	while (1) {
	    if (rn <  pcmp[a]) --a; else
	    if (rn >= pcmp[a]) ++a;
	    else break;
	}
	return (a + ALA);
}

Seq* Seq::deamb(const int& bzx)
{
	if (!inex.ambs) return this;
	CHAR*	ss = at(left) - 1;
	CHAR*	tt = at(right);

	if (isdrna()) {
	    while (++ss < tt) {
		int	namb = isAmb(*ss);
		if (!namb) continue;
		int	c = *ss - _;
		int	l = 1;
		for (int n = rand() % namb; ; l <<= 1)
		    if ((c & l) && !n--) break;
		*ss = l + _;
	    }
	} else {
	    while (++ss < tt) {
		if ((bzx & 1) && *ss == AMB) *ss = randaa();
		if (bzx & 2)
		    if (*ss == ASX) *ss = rand() & 0x1000? ASN: ASP;
		if (bzx & 4)
		    if (*ss == GLX) *ss = rand() & 0x1000? GLN: GLU;
	    }
	}
	inex.ambs &= ~bzx;
	return this;
}

void Seq::copyattr(Seq* dest) const
{
	if (msaname) dest->msaname = strrealloc(dest->msaname, msaname);
	dest->code = code;
	dest->inex = inex;
	dest->did = did;
}

void Seq::fillnbr(Seq* dest) const
{
const	CHAR*	ss = at(0);
const	CHAR*	tt = at(left);

	if (!dest->lens) dest->lens = new int[many];
	vcopy(dest->lens, lens, many);
	if (!dest->nbr) dest->nbr = new int[many];
	for (int i = 0; i < many; ++i) {
	    dest->nbr[i] = nbr[i];
	    for (const CHAR* s = ss++; s < tt; s += many)
		if (IsntGap(*s)) dest->nbr[i]++;
	}
}

void Seq::copynbr(Seq* dest) const
{
	if (!dest->lens) dest->lens = new int[many];
	vcopy(dest->lens, lens, many);
	if (!dest->nbr) dest->nbr = new int[many];
	vcopy(dest->nbr, nbr, many);
}

void Seq::copylbl(Seq* dest) const
{
	if (dest->sname) dest->sname->assign(*sname);
	else dest->sname = new Strlist(*sname);
}

#if USE_WEIGHT
void Seq::copyweight(Seq* dest) const
{
	if (weight) {
	    if (!dest->weight) dest->weight = new FTYPE[many];
	    vcopy(dest->weight, weight, many);
	}
	if (pairwt) {
	    int nn = many * (many - 1) / 2;
	    dest->pairwt = new FTYPE[nn];
	    vcopy(dest->pairwt, pairwt, nn);
	}
}
#endif

// dest is a new seq sliced from this's segment

Seq* Seq::cutseq(Seq* dest, const int& snl) const
{
const	CHAR*	s = at(left);
	int	l = right - left;

	if (dest) dest->refresh(many, l);
	else	  dest = new Seq(many, l);
	if (snl & CPY_NBR) fillnbr(dest);
	if (snl & CPY_LBL) copylbl(dest);
	if (l <= 0)	return(dest);		// nothing to be copied
	copyattr(dest);
#if USE_WEIGHT
	copyweight(dest);
#endif
	memcpy(dest->at(0), s, l * many);
	dest->postseq(dest->at(l));
	if (snl & CPY_SBI) cutSigII(dest, this);
	return (dest);
}

// dest is a copy of this's segment
// seq numbering is preserved

Seq* Seq::copyseq(Seq* dest, const int& snl) const
{
	if (snl & CPY_SEQ) {
	    if (dest)	dest->refresh(many, len);
	    else	dest = new Seq(many, len);
	    dest->left = left;
	    dest->right = right;
	    dest->len = len;
	    dest->tlen = tlen;
	    dest->base_ = base_;
	    memcpy(dest->seq_ - many, seq_ - many, many * (len + 2));
	}
	if (snl & CPY_NBR) copynbr(dest);
	if (snl & CPY_LBL) copylbl(dest);
	if (snl & CPY_SBI) dest->sigII = copySigII(sigII);
	copyattr(dest);
#if USE_WEIGHT
	copyweight(dest);
#endif
	return (dest);
}

Seq* Seq::duplseq(Seq* dest) const
{
	int	l = 0;
	int	r = len;

	swap(l, left);
	swap(r, right);
	dest = copyseq(dest, CPY_ALL);
	left = l;
	right = r;
	return (dest);
}

// catenate this after head with a spacer of length cushion

Seq* Seq::catseq(Seq* head, const int& cushion)
{
	CHAR*	s = at(left);
	int	l = right - left + cushion;

	if (l <= 0)	return(head);		// nothing to be append
	if (!head || head->left == head->right)
	    return (cutseq(head, CPY_ALL));
	if (head->many != many) return (head);	// cannot catenate
	catSigII(head, this, cushion);
	head->seqalloc(head->many, l + head->len, true);
	if (cushion) memset(head->at(head->len), N, cushion * many);
	memcpy(head->at(head->len + cushion), s, (right - left) * many);
	head->right = head->len += l;
	head->fillpad();
	head->test_gap_amb();
	return (head);
}

Seq* Seq::extseq(Seq* dest, const int* which, const int& snl, const FTYPE& nfact)
{
	int	dmany = 0;
	for ( ; which[dmany] >= 0; ++dmany) ;
const	CHAR*	sseq = at(left);
	if (snl & CPY_SEQ) {
	    if (dest)	dest->refresh(dmany, right - left);
	    else	dest = new Seq(dmany, right - left);
	    CHAR*	dseq = dest->at(0);
	    copyattr(dest);
	    dest->inex.vect = dest->inex.dels = dest->inex.ambs = 
	    dest->inex.prof = dest->inex.algn = 0;
	    dest->inex.sngl = dmany == 1;
	    for (int i = left; i < right; ++i, sseq += many) {
		for (const int* wk = which; *wk >= 0; ) {
		    const CHAR	s = sseq[*wk++];
		    if (isGap(s)) dest->inex.dels = 1;
		    if (isAmb(s)) dest->inex.ambs = 1;
		    *dseq++ = s;
		}
	    }
	    dest->postseq(dseq);
	}
	if (!dest) return (0);
	if (msaname) dest->msaname = strrealloc(dest->msaname, msaname);
	if (snl & CPY_NBR) {
	    const int*	wk = which;
	    for (int i = 0; *wk >= 0; ++i, ++wk) {
		dest->lens[i] = lens[*wk];
		dest->nbr[i] = nbr[*wk];
	    }
	}
	if (snl & CPY_LBL) {
	    const int*	wk = which;
	    if (dest->sname) dest->sname->reset(dmany);
	    else dest->sname = new Strlist(dmany);
	    for (int i = 0; *wk >= 0; ++i, ++wk)
		dest->sname->push((*sname)[*wk]);
	}
	if (snl & CPY_SBI) {
	    SigII*	dst = dest->sigII = extSigII(this, which, nfact, snl & RLT_SBI);
	    if (dst) dst->resetend(dest->len);
	}
#if USE_WEIGHT
// weights are normalized to sum up to many
	if (weight) {
	    if (!dest->weight) dest->weight = new FTYPE[dmany];
	    if (dmany == 1) dest->weight[0] = 1;
	    else { 
		const int*	wk = which;
		for (int i = 0; *wk >= 0; ++i, ++wk) {
		    dest->weight[i] = weight[*wk] / nfact;
		}
	    }
	}
#endif
	return (dest);
}

Seq* Seq::splice(Seq* dest, RANGE* rng, const int bias)
{
	RANGE	svr;

	saverange(&svr);
	if (dest) dest->refresh(many, len);
	else	  dest = new Seq(many, len);
	rng = fistrng(rng);
	if (bias) {
	    for ( ; neorng(rng); ++rng) {
		left = rng->left - bias;
		right = rng->right - bias;
		dest = catseq(dest);
	    }
	} else {
	    while (neorng(rng)) {
		restrange(rng++);
		dest = catseq(dest);
	    }
	}
	restrange(&svr);
	return (dest);
}

void Seq::elim_amb()
{
	if (many > 1) return;
	int	n = left;
	GAPS	gbf = {0, n};
	CHAR*	ss = at(left);
	CHAR*	tt = at(right);
	CHAR*	dd = ss;
	Mfile*	mfd = 0;

	if (sigII) {
	    mfd = new Mfile(sizeof(GAPS));
	    mfd->write(&gbf);	/* reserve */
	    mfd->write(&gbf);	/* first record */
	}
	for ( ; ss < tt; ++n, ++ss) {
	    if (isAmb(*ss)) {
		if (sigII) {
		    if (!gbf.gln) gbf.gps = n;
		    gbf.gln++;
		}
	    } else {
		*dd++ = *ss;
		if (sigII && gbf.gln) {
		    mfd->write(&gbf);
		    gbf.gln = 0;
		}
	    }
	}
	if (sigII) {
	    if (gbf.gln) mfd->write(&gbf);
	    gbf.gps = right;
	    gbf.gln = gaps_end;
	    if (gbf.gln) mfd->write(&gbf);
	    n = mfd->size();
	    GAPS*	gaps = (GAPS*) mfd->flush();
	    delete mfd;
	    gaps->gps = left;
	    gaps->gln = n;
	    sigII->rmGapPfq(gaps);
	    delete[] gaps;
	}
	if (dd == tt) return;
	int k = at(len) - tt;
	if (k) memcpy(dd, ss, k);
	postseq(dd + k);
}

Seq* Seq::rndseq()
{
	int 	width = right - left;
	CHAR*	seq = at(left);
	CHAR*	as = seq;
	CHAR*	buf = new CHAR[many];

	for (int i = left; i < right; i++, as += many) {
	    CHAR*	tmp = seq + many * (rand() % width);
	    memcpy(buf, tmp, many);
	    memcpy(tmp, as, many);
	    memcpy(as, buf, many);
	}
	delete[] buf;
	return (this);
}

char* Seq::onecds(RANGE& wexon, char* ps, int& par)
{
	int	qar = par;
	wexon = zerorng;
	while (*ps && !isdigit(*ps)) ++ps;
	if (!*ps) return (0);
	wexon.left = atoi(ps) - 1;
	while (isdigit(*ps)) ++ps;
	while (*ps && !isdigit(*ps))
	    if (*ps++ == ')') --qar;
	if (!*ps) return (0);
	wexon.right = atoi(ps);
	while (isdigit(*ps)) ++ps;
	if (!*ps) return (0);
	while (*ps && !isdigit(*ps))
	    if (*ps++ == ')') --qar;
	par = qar;
	return (ps);
}

int Nprim_code(int c)
{
	int	rn = ncredctab[c];

	if (rn < 4) return (rn);
	if (c == 0) c = 15;
	do {
	    rn = (rand() >> 8) & 3;
	} while (!(c & red2nuc[rn]));
	return (rn);
}

int en_code(int c, const SEQ_CODE* code)
{
	if (isalpha(c)) {
	    c = toupper(c);
	    c = code->encode? code->encode[c - 'A']: aton(c);
	    return (c < code->max_code)? c: IGNORE;
	} else switch (c) {
	    case _UNP:	return (gap_code);
	    case _TRM:	return (code->amb_code);
	    case _EOS:	return (EOF);
	    default:	return (IGNORE);
	}
}

CHAR* tosqcode(CHAR* ns, const SEQ_CODE* code)
{
	CHAR*	ps = ns;
	CHAR*	bs = ns;

	while (*ns) {
	    int	c = en_code(*ns++, code);
	    if (c > 0) *ps++ = c;
	    else if (c == 0 || c == EOF) break;
	}
	*ps = 0;
	return (bs);
}

int infermolc(const char* fname)
{
	Seq	sd(1);
	char	str[MAXL];
	if (is_gz(fname)) {
#if USE_ZLIB
	    gzFile	gzfd = gzopen(fname, "r");
	    if (!gzfd) fatal("%s not found !\n", fname);
	    if (!fgets(str, MAXL, gzfd)) fatal("%s is empty !\n", fname);
	    sd.infermolc(gzfd, str);
	    fclose(gzfd);
#else
	    fatal(gz_unsupport, fname);
#endif
	} else {
	    FILE*	fd = fopen(fname, "r");
	    if (!fd) fatal("%s not found !\n", fname);
	    if (!fgets(str, MAXL, fd)) fatal("%s is empty !\n", fname);
	    sd.infermolc(fd, str);
	    fclose(fd);
	}
	return (sd.inex.molc);
}

// remove polyA sequences if present

SHORT PolyA::rmpolyA(Seq* sd, const int& q_mns) const
{
	if (sd->isprotein() || sd->inex.intr) return (1);
	int	rv = q_mns? q_mns: 3;
	if (thr <= 0) return (rv);	// has treated
	if (sd->inex.polA) return (sd->inex.polA == 3? rv: sd->inex.polA);
	int	polya = 0;
	int	polyt = 0;
	int	scr = 0;
	CHAR*	maxa = 0;
	CHAR*	maxt = 0;
static	const	int	mn[2] = {-5, 1};

	CHAR*	lend = sd->at(0);
	CHAR*	s = sd->at(sd->len);
	CHAR*	t = s;
	if (q_mns != 2) {
	    while (--s >= lend) {
		scr += mn[*s == A];
		if (scr > polya) {
		    polya = scr;
		    if (scr > thr) maxa = s;
		}
		if (scr < polya - thr) break;
	    }
	}
	if (q_mns != 1) {
	    for (scr = 0, s = lend; s < t; ++s) {
		scr += mn[*s == T];
		if (scr > polyt) {
		    polyt = scr;
		    if (scr > thr) maxt = s;
		}
		if (scr < polyt - thr) break;
	    }
	}
	if (maxa && maxt) {
	    if (polya >= polyt) maxt = 0;
	    else		maxa = 0;
	}
	sd->inex.polA = 3;
	if (maxa) {
	    sd->tlen = maxa - lend;
	    if (sd->right > sd->tlen) sd->right = sd->tlen;
	    sd->inex.polA = 1;
	} else if (maxt) {
	    int	polytlen = maxt - lend;
	    if (sd->left < polytlen) sd->left = polytlen;
	    sd->tlen = sd->len - polytlen;
	    sd->inex.polA = 2;
	    sd->comrev();
	}
	return (q_mns == 3? 3: 
	    ((sd->inex.polA == 1 || sd->inex.polA == 2)? sd->inex.polA: rv));
}

StrPhrases::StrPhrases(const char* fname)
{
	FILE*	fd = ftable.fopen(fname, "r");

	if (!fd) {
	    strphrase[0] = strphrase[1] = 0;
	    strprhase_buf = 0;
	    return;
	}
	fseek(fd, 0L, SEEK_END);
	long	filelen = ftell(fd);
	Mfile	mfd(sizeof(STRPHRASE));
	STRPHRASE	rec;
	char	str[MAXL];
	strprhase_buf = new char[filelen];
	char*	bf = strprhase_buf;
	for (int i = 1; i <= 2; ++i) {
	    rewind(fd);
	    while (char* ps = fgets(str, MAXL, fd)) {
		if (*ps == '#' || *ps == '\n') continue;
		rec.sens = atoi(ps);
		ps = cdr(ps);
		int	ioc = atoi(ps);
		if (ioc != i) continue;
		ps = cdr(ps);
		rec.phrase = bf;
		while(*ps && !isspace(*ps)) *bf++ = *ps++;
		*bf++ = '\0';
		mfd.write(&rec);
	    }
	    rec.sens = 0;
	    rec.phrase = 0;
	    mfd.write(&rec);
	    if (i == 1) filelen = mfd.size();
	}
	strphrase[0] = (STRPHRASE*) mfd.flush();
	strphrase[1] = strphrase[0] + filelen;
	fclose(fd);
}

static StrPhrases* strphrs = 0;

void Seq::setstrand(const int idorcom, const char* text)
{
	if (!strphrs) strphrs = new StrPhrases(StrandPhrase);
	STRPHRASE*	sp = strphrs->strphrase[idorcom];

	if (!sp) return;
	for ( ; sp->phrase; ++sp) {
	    if (strstr(text, sp->phrase)) {
		inex.sens = sp->sens;
		inex.ori  = 1;
		break;
	    }
	}
}

void eraStrPhrases() {delete strphrs;}

/**********************************************
*	Get multiple sequences in native format
***********************************************/

char* Seq::sqline(int i, char* ps)
{
	char	str[16];

	if (sname->unfilled()) {
	    if (sscanf(ps, "%d", nbr + i)) --nbr[i];
	    else	nbr[i] = 0;
	    char*	pl = strchr(ps, _LABL);
	    if (pl) while (*++pl && isspace(*pl));
	    if (pl && *pl) {
		sname->push(pl);
	    } else {
		sprintf(str, "m%d", i + 1);
		sname->push(str);
	    }
	}
	if (isnbr(ps)) ps++;
	while (isdigit(*ps) || isspace(*ps)) ++ps;
	return (ps);
}

void Seq::header_nat_aln(const int& n, const FTYPE& sumwt)
{
	sname->setfill();
#if USE_WEIGHT
	if (weight) {
	    if (n != many)
		fprintf(stderr, "Insufficient weights %d < %d!\n", n, many);
	    VTYPE	f = sumwt / many;
	    if (f < fepsilon) vset(weight, (VTYPE) 1, many);	// equal weight
	    else if (fabs(1. - f) > fepsilon) {		// sum up to many
		for (int i = 0; i < many; i++)
		    weight[i] /= f;
	    }
	    if (sigII) sigII->rescale_dns(f);
	}
#endif
}

int cmpPfqPos(const PFQ* a, const PFQ* b)
{
	int	i = a->pos - b->pos;

	if (i) return (i);
	return (a->num - b->num);
}

CHAR* Seq::ToInferred(CHAR* src, CHAR* lastseq, const int& step)
{
	FTYPE	total = 0;

	for (int i = 2; i < 27; ++i) total += cmps[i];
	FTYPE	ncode = cmps[aton('A')] + cmps[aton('C')] + cmps[aton('G')]
		+ cmps[aton('T')] + cmps[aton('U')];
	int	molc = PROTEIN;
	if (total == cmps[aton('N')]) molc = DNA;
	else if (100 * ncode/(total - cmps[aton('N')]) > MinPctNucChar)
	    molc = (cmps[aton('T')] >= cmps[aton('U')])? DNA: RNA;
	else if (100 * (cmps[aton('J')] + cmps[aton('O')]
		 + cmps[aton('U')]) / total > MinPctTronChar) {
	    molc = TRON;
	    prompt("Warning: %s is regarded as TRON sequence !\n", sqname());
	}
	setSeqCode(this, molc);
	CHAR*	tab = code->encode;
	CHAR*	dst = src;
	for ( ; src < lastseq; src += step) {
	    if (IsGap(*src)) {*dst = gap_code; dst += step;}
	    else if ((*dst = tab[*src - A]) < code->max_code) dst += step;
	}
	return (dst);
}

void Seq::estimate_len(FILE* fd, const int& nos, const SeqDb* dbf)
{
	long	fpos = ftell(fd);
	if (nos == 1) {		// FASTA single sequence
	    char	str[MAXL];
	    while (fgets(str, MAXL, fd)) {
		if (*str == _NHEAD || *str == _CHEAD || *str == _EOS) {
		    if ((strlen(str) + 1) == MAXL) flush_line(fd);
		    break;
		}
		if (*str == _COMM || *str == _LCOMM || *str == _WGHT) {
		    if ((strlen(str) + 1) == MAXL) flush_line(fd);
		    continue;
		}
		len += strlen(str) - 1;
	    }
	} else {
	    fseek(fd, 0L, SEEK_END);
	    area_ = (int) ftell(fd);
	    len = area_ / nos;
	}
	fseek(fd, fpos, SEEK_SET);
}

#if USE_ZLIB
void Seq::estimate_len(gzFile gzfd, const int& nos, const SeqDb* dbf)
{
	long	fpos = ftell(gzfd);
	char	str[MAXL];
	while (fgets(str, MAXL, gzfd)) {
	    if (*str == _NHEAD || *str == _CHEAD || *str == _EOS ||
		*str == _COMM || *str == _LCOMM || *str == _WGHT) {
		if ((strlen(str) + 1) == MAXL) flush_line(gzfd);
		continue;
	    }
	    if (nos == 1 || (dbf && dbf->FormID == FASTA)) {
		len += strlen(str) - 1;
	    } else {
		Strlist	stl(str, stddelim);
		if (stl.size() > 1) 
		    len += strlen(stl[1]) - 1;
	    }
	}
	area_ = len;
	len /= nos;
	fseek(gzfd, fpos, SEEK_SET);
}
#endif

Seq* Seq::getseq(const char* str, DbsDt* dbf)
{
	FILE*	fd = 0;
#if USE_ZLIB
	gzFile	gzfd = 0;
#endif
	char	input[MAXL];

	if (!str) {		//  Interactive	mode
	    for (;;) {
		fprintf(stderr, "seq %d := %s", sid, seqdfn);
		if (!fgets(input, MAXL, stdin)) return this;
		for (str = input; isspace(*str); ) ++str;
		if (isattrib(str)) return attrseq(str);
		if (*str == DBSID) {	// Get from Database File
		    if (getdbseq(dbf, str)) return attrseq(cdr(str));
		    else continue;
		}
#if USE_ZLIB
		fd = openseq(str, &gzfd);
		if (fd || gzfd) break;
#else
		fd = openseq(str);
		if (fd) break;
#endif
	    }
	} else {
	    while (isspace(*str)) ++str;
	    if (isattrib(str)) return attrseq(str);
	    if (*str == DBSID) 	/* Get from Database File */
		return getdbseq(dbf, str);
	    if (!strcmp(str, "-")) {
		fd = stdin;
		sprintf(input, "std%d", sid);
		spath = strrealloc(spath, input);
		if (sname) sname->assign(input);
		else	sname = new Strlist(input, 0);
	    } else {
#if USE_ZLIB
		fd = openseq(str, &gzfd);
#else
		fd = openseq(str);
#endif
	    }
	}
#if USE_ZLIB
	if (gzfd) {
	    Seq*	sd = fgetseq(gzfd, cdr(str));
	    fclose(gzfd);
	    return (sd);
	}
#endif
	if (!fd) return (0);
	Seq*	sd = fgetseq(fd, cdr(str));
	if (fd != stdin) fclose(fd);
	return (sd);
}

Seq* Seq::apndseq(const char* aname)
{
	Seq* temp = new Seq(many);
	temp->getseq(aname);
	if (temp->empty()) return this;
	catseq(temp);
	delete temp;
	return this;
}

Seq* inputseq(Seq* seqs[], char* str)
{
static	char aa[] = "aa";
static	char nt[] = "nt";
static	char tc[] = "tc";
	int 	i = 0;
	char*	res;
	Seq*	sd = 0;

	sscanf(str, "%d", &i);
	if (0 < i && --i < noseq) {
	    sd = seqs[i];
	    while (isdigit(*str)) str++;
	    if (*str == _APPN)	{
		sd->apndseq(0);
		str++;
	    } else
		sd->getseq(0);
	    if (sd->empty()) return(sd);
	    switch (sd->inex.molc) {
		case PROTEIN:	res = aa; break;
		case TRON:	res = tc; break;
		case DNA: case RNA: case GENOME:
		default:	res = nt; break;
	    }
	    fprintf(stderr, " %c%s", sd->senschar(), sd->sqname());
	    if (sd->inex.vect) putc('%', stderr);
	    fprintf(stderr, "[%d] ( %d %s\'s )  [ %d - %d ]\n", 
		sd->many, sd->len, res, sd->SiteNo(sd->left),
		sd->SiteNo(sd->right - 1));
	}
	return (sd);
}

int samerange(const Seq* a, const Seq* b)
{
	if (!a || !b || !a->sname || !b->sname) return (0);
	bool	half = !strcmp((*a->sname)[0], (*b->sname)[0]) && 
	    a->SiteNo(a->left) == b->SiteNo(b->left) &&
	    a->SiteNo(a->right) == b->SiteNo(b->right);
	if (!half) return(0);
	return (((a->inex.sens ^ b->inex.sens) & REVERS)? -1 : 1);
}

void save_range(const Seq* seqs[], RANGE* temp, int n) 
{
	while (n--) {
	    temp->left = (*seqs)->left;
	    (temp++)->right = (*seqs++)->right;
	}
}

void rest_range(const Seq* seqs[], const RANGE* temp, int n)
{
	while (n--) {
	    (*seqs)->left = temp->left;
	    (*seqs++)->right = (temp++)->right;
	}
}

void setthr(const double& thr)
{
	if ((int) thr == QUERY)
		promptin("Threshold (%.1f) : ", &alprm.thr);
	else
		alprm.thr = (float) thr;
	algmode.thr = alprm.thr != 0;
}

void setminus(int yon)
{
	switch (yon) {
	    case 0:	algmode.mns = 0; break;
	    case 1:	algmode.mns = 1; break;
	    default:	
		yon = tolower(progetc("Both direction ? [y/n] : "));
		if (yon == 'y')		algmode.mns = true;
		else if (yon == 'n')	algmode.mns = false;
	    break;
	}
}

int setoutmode(int nsa)
{
	if (nsa >= 0) algmode.nsa = nsa;
	else {
	    nsa = algmode.nsa;
	    promptin("OutputMode [0-3]\t(%d) : ", &nsa);
	    algmode.nsa = nsa;
	}
	return (nsa);
}

void setalgmode(int nsa, int reg)
{
	setoutmode(nsa);
	if (reg >= 0) algmode.alg = reg;
	else if (reg == QUERY) {
	    reg = algmode.alg;
	    promptin("Algorithm\t[0-4]\t(%d) : ", &reg);
	    algmode.alg = reg;
	}
}

void setstrip(const int& sh)
{
	if (sh >= 0)
	    alprm.sh = sh;
	else if (sh == QUERY)
	    promptin("BandWidth (%d) : ", &alprm.sh);
}

bool Seq::test_aligned() 
{
const 	CHAR*	ss = at(left);
const 	CHAR*	tt = at(right);

	if (inex.algn) return (true);
	for (int i = 0; i < many; ++i) {
	    int	alt = 0;
	    int	prv = 1;
	    for (const CHAR* s = ss + i; s < tt; s += many) {
		int	now = (IsGap(*s));
		if (now ^ prv) ++alt;
		if (alt > 2) {
		    inex.algn = 3;
		    return (true);
		}
		prv = now;
	    }
	}
	return (false);
}

FTYPE* Seq::composition()
{
	if (cmps && inex.cmps) return (cmps);
	CHAR*	ss = at(left);
	CHAR*	tt = at(right);
	if (!cmps) cmps = new FTYPE[code->max_code];
	vclear(cmps, code->max_code);
	inex.cmps = 1;
#if USE_WEIGHT
	if (weight) {
	    while (ss < tt) {
		FTYPE*	ww = weight;
		FTYPE*	wt = weight + many;
		while (ww < wt) cmps[*ss++] += *ww++;
	    }
	    return (cmps);
	}
#endif
	while (ss < tt) cmps[*ss++]++;
	return (cmps);
}

FTYPE* Seq::composition(FTYPE* cmp) const
{
	if (!cmp) return (0);
	if (cmps && inex.cmps)
	    return (vcopy(cmp, cmps, code->max_code));
	CHAR*	ss = at(left);
	CHAR*	tt = at(right);
	vclear(cmp, code->max_code);
#if USE_WEIGHT
	if (weight) {
	    while (ss < tt) {
		FTYPE*	ww = weight;
		FTYPE*	wt = weight + many;
		while (ww < wt) cmp[*ss++] += *ww++;
	    }
	    return (cmp);
	}
#endif
	while (ss < tt) cmp[*ss++]++;
	return (cmp);
}

int Seq::countgap(const int& mem, int from, int to) const
{
	if (from < left) from = left;
const	CHAR*	ss = at(from) + mem;
	if (to > right) to = right;
const	CHAR*	tt = at(to);
	int	n = 0;

	for ( ; ss < tt; ss += many)
	    if (IsGap(*ss)) ++n;
	return (n);
}

VTYPE Seq::countunps() const
{
const	CHAR*	ss = at(left);
const	CHAR*	tt = at(right);
	VTYPE	unps = 0;

	for (int i = 0; i < many; ++i) {
	    int	nu = 0;
	    for (const CHAR* wk = ss++; wk < tt; wk += many)
		if (IsTrueGap(*wk)) ++nu;
#if USE_WEIGHT
	    if (weight)
		unps += (VTYPE) (nu * weight[i]);
	    else
#endif
		unps += nu;
	}
	return (unps);
}

bool Seq::isgap(const CHAR* ps) const
{
const	CHAR*	ts = ps + many;
	for (; ps < ts; ++ps) 
	    if (!isGap(*ps)) return (false);
	return (true);
}

bool Seq::nogap(const CHAR* ps) const
{
const	CHAR*	ts = ps + many;
	for (; ps < ts; ++ps) 
	    if (isGap(*ps)) return (false);
	return (true);
}

void Seq::pregap(int* gl) const
{
const	CHAR*	ss = at(left);

	for (int i = 0; i < many; ++i) {
	    const	CHAR*	s = ss++;
	    int	n = left;
	    for ( ; n > 0; --n) {
		s -= many;
		if (IsntGap(*s)) break;
	    }
	    *gl++ = left - n;
	}
}

// return information on deleted columns when (which & RET_GAP)

GAPS* Seq::elim_column(const int& which, float frac)
{
	bool	store_gap = (which & RET_GAP) || sigII;
	int	cond = which & ~RET_GAP;
	FTYPE	w = 0;
#if USE_WEIGHT
	FTYPE*	wt = weight;
#else
	FTYPE*	wt = 0;
#endif

	if (frac > 1.) frac /= 100.;
	float	thr = (1. - frac) * many;
	if (!wt) {
	    if (cond == DEL_GRC) cond = DEL_GAP;
	    if (cond == DEL_ARC) cond = DEL_AMB;
	}
	Mfile*	mfd = store_gap? new Mfile(sizeof(GAPS)): 0;
	CHAR*	ss = at(left);
	CHAR*	tt = at(right);
	CHAR*	dd = ss;
	int	n = left;
	GAPS	gbf = {n, 0};
	if (mfd) {
	    mfd->write(&gbf);			// make space
	    mfd->write(&gbf);			// 
	}
	for ( ; ss < tt; ++n) {
	    int	k = 0;
	    int	ndel = 0;
	    int	namb = 0;
#if USE_WEIGHT
	    if (wt) {
		w = 0;
		wt = weight;
	    }
	    for (int i = 0; i < many; ++i) {
		if (isGap(*ss)) {
		    ++ndel;
		    if (cond == DEL_GAP) ++k; else
		    if (cond == DEL_GRC) w += *wt;
		} else if (isAmb(*ss)) {
		    ++namb;
		    if (cond == DEL_AMB) ++k; else
		    if (cond == DEL_ARC) w += *wt;
		}
		*dd++ = *ss++;
		if (wt) ++wt;
	    }
#else
	    for (int i = 0; i < many; ++i) {
		if (isGap(*ss)) {
		    ++ndel;
		    if (cond == DEL_GAP) ++k;
		} if (isAmb(*ss)) {
		    ++namb;
		    if (cond == DEL_AMB) ++k;
		}
		*dd++ = *ss++;
	    }
#endif
	    bool y = false;			// delete this column?
	    switch (which & ~RET_GAP) {
		case DEL_GAP:
		case DEL_AMB:
		     y = k == many; break;
		case DEL_GRC:
		case DEL_ARC:
		    if (cond == DEL_GAP || cond == DEL_AMB) w = k;
		    y = w >= thr; break;
	    }
	    if (y) {
		dd -= many;
		if (!gbf.gln) gbf.gps = n;	// new gap
		gbf.gln++;			// elong gap
	    } else if (gbf.gln) {
		if (ndel) inex.dels = 1;
		if (namb) inex.ambs = 1;
		if (store_gap) mfd->write(&gbf);
		gbf.gln = 0;			// close gap
	    }
	}
	GAPS*	gaps = 0;
	if (mfd) {
	    if (gbf.gln) mfd->write(&gbf);
	    gbf.gps = n; gbf.gln = gaps_end;
	    mfd->write(&gbf);
	    int	m = mfd->size();
	    gaps = (GAPS*) mfd->flush();
	    gaps->gps = n; gaps->gln = m;
	    delete mfd;
	}
	RANGE	rng;
	rng.left = left;
	rng.right = index(dd);
	int k = at(len) - at(right);
	if (k) memcpy(dd, ss, k);
	postseq(dd + k);
	restrange(&rng);
	if (sigII) sigII->rmGapPfq(gaps);
	if (which & RET_GAP) return (gaps);
	delete[] gaps;
	return (0);
}

/*****************************************************
	Generate a random sequence 
*****************************************************/

Seq* Seq::randseq(double* pcmp)
{
	char	name[MAXL];
	int	molc = inex.molc;

	refresh(1, len);
	inex.molc = molc;
	sprintf(name, "RAND%d", sid);
	if (sname) sname->assign(name);
	else	sname = new Strlist(name, 0);
	CHAR*	s = at(0);
	RandNumGen	rn(pcmp, isprotein()? 20: 4);
	for (int n = 0; n < len; ++n)
	    *s++ = r2s(rn.get());
	return (postseq(s));
}

/*****************************************************
	Random substitution of 'n' times
	if which >= 0 which's member is substituted
 	else substitutions may occur in any members
*****************************************************/

Seq* Seq::substseq(int n, const int& which)
{
	composition();
	INT	nelm = isprotein()? 20: 4;
	double*	pcmp = new double[nelm];
	double	sum = 0;
	for (INT i = 0; i < nelm; ++i) {
	    INT	k = r2s(i);
	    sum += cmps[k];
	    pcmp[i] = cmps[k];
	}
	if (sum == 0) {
	    delete[] pcmp;
	    return (this);
	}
	for (INT i = 0; i < nelm; ++i) pcmp[i] /= sum;

	CHAR*	ss = at(left);
	int	area = right - left;
	if (which >= 0)  ss += which;
	else	    area *= many;
	RandNumGen	rn(pcmp, nelm);
	while (n > 0) {
	    int	k = (int) (drand48() * area);
	    CHAR*	ps = ss + (which >= 0? k * many: k);
	    if (IsGap(*ps)) continue;
	    while ((k = r2s(rn.get())) == *ps) ;
	    *ps = k;
	    --n;
	}
	delete[] pcmp;
	return (this);
}

#if USE_WEIGHT

FTYPE* Seq::saveseqwt() const
{
	if (!weight) return (0);
	FTYPE*	tmpwt = new FTYPE[many];
	return (vcopy(tmpwt, weight, many));
}

void Seq::restseqwt(const FTYPE* tmpwt)
{
	if (!tmpwt) return;
	if (!weight) weight = new FTYPE[many];
	vcopy(weight, tmpwt, many);
}

#endif
