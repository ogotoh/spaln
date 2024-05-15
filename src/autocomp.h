/*****************************************************************************
*
*	Automatic performance of program "perform" for a given set
*	of protein or nucleotide sequences.
*
*	Subcommand: 'e' all; 		'f' single scan; 	'g' group
*	Subcommand: 'i' internal; 	'j' juxtapose		'k' addon
*	Subcommand: 'o' onesequence;
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-2023)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
*
*****************************************************************************/

#ifndef  _AUTOCOMP_
#define  _AUTOCOMP_

#include "calcserv.h"
#include "aln.h"
#include "phyl.h"

static	const	int	CONLINE = 80;
static	const	int	MAX_ARGS = 128;
static	const	int	MAX_ARGTEXT = 4096;
static	const	char	tmpchr[] = "_s";
static	const	char	logfn[] = "aln.log";

//	specific to mseq

struct InputSeqTest {Strlist* sname; int num, many, maxmany, space, bad, maxlen, minlen; 
			INT molc, dels;};

template <class seq_t> class SeqLoader;

template <class seq_t>
class AlnServer : public CalcServer<seq_t> {
protected:
	void	initialize();
	int	readargs();
	int	serialJob();
	void	getoptions(int& argc, const char**& argv);
	char*	getaseq(char* ps);
	void	as_initialize()
	  {
	    logfd = 0; logfile = 0;
	    sql1 = sql2 = 0;
	  }

public:
const	char*	logfile;
	FILE*	logfd;
	char	tmpdir[LINE_MAX];
	Subset*	ss;
	SeqLoader<seq_t>*	sql1;
	SeqLoader<seq_t>*	sql2;
	AlnServer(seq_t** sqs, int nos, int im, void* p,
	    int 	(*mj)(CalcServer<seq_t>*, seq_t**, ThQueue<seq_t>*),
	    void	(*uj)(CalcServer<seq_t>*) = 0,
	    void	(*cj)(CalcServer<seq_t>*) = 0,
	    int 	(*bj)(CalcServer<seq_t>*, seq_t**, ThQueue<seq_t>*) = 0) :
	    CalcServer<seq_t>(im, p, mj, uj, cj, sqs, nos, 0, true),
	    sub_job(bj)
	  {
	    as_initialize();
	    ss = 0;
	  }
	AlnServer(seq_t* hsd, int im, Subset* sst, INT asis, void* p,
	    int 	(*mj)(CalcServer<seq_t>*, seq_t**, ThQueue<seq_t>*),
	    void	(*uj)(CalcServer<seq_t>*) = 0,
	    void	(*cj)(CalcServer<seq_t>*) = 0,
	    int 	(*bj)(CalcServer<seq_t>*, seq_t**, ThQueue<seq_t>*) = 0) :
	    CalcServer<seq_t>(im, p, mj, uj, cj, 0, 0, 0, false), 
	    ss(sst), sub_job(bj)
	  {
	    as_initialize();
	    this->msd = hsd;
	    sql1 = new SeqLoader<seq_t>(this, 0, 0, this->getgrp22());
	    hsd->inex.est = asis;
	  }
	AlnServer(int& argc, const char**& argv, int im, int defim, void* p = 0,
	    int 	(*mj)(CalcServer<seq_t>*, seq_t**, ThQueue<seq_t>*) = 0,
	    void	(*uj)(CalcServer<seq_t>*) = 0,
	    void	(*cj)(CalcServer<seq_t>*) = 0,
	    int 	(*bj)(CalcServer<seq_t>*, seq_t**, ThQueue<seq_t>*) = 0) :
	    CalcServer<seq_t>(defim, p, mj, uj, cj, 0, 0, 0, false), 
	    sub_job(bj)
	  {	// 	This constructor requires localoption()
	    as_initialize();
	    ss = 0;
	    getoptions(argc, argv);
	    this->input_mode = im? im: this->calc_mode;
	    if (dbs_dt[0]) setdefdbf(dbs_dt[0]);
	    setup_output(algmode.nsa);
	    int	fin = 0;
	    switch (this->input_mode) {
	      case IM_PARA: fin = 1; 
	      case IM_ALTR: case IM_EVRY: case IM_GRUP: case IM_SCAN:
		this->input_ns = 2;
		break;
	      case IM_ADON: case IM_TREE:
		this->input_ns = 0;
		break;
	      default:
		this->input_ns = 1;
		break;
	    }
	    this->memb = new InFiles(this->catalog, argc, argv, 
		this->input_mode, this->calc_mode);
	    sql1 = new SeqLoader<seq_t>(this, fin, 0, this->getgrp22());
	  }
	AlnServer(int& argc, const char**& argv, int im, const char* clg = 0, void* p = 0,
	    int 	(*mj)(CalcServer<seq_t>*, seq_t** sqs, ThQueue<seq_t>*) = 0,
	    void	(*uj)(CalcServer<seq_t>*) = 0,
	    void	(*cj)(CalcServer<seq_t>*) = 0,
	    int 	(*bj)(CalcServer<seq_t>*, seq_t** sqs, ThQueue<seq_t>*) = 0) :
	    CalcServer<seq_t>(im, p, mj, uj, cj, 0, 0, 0, false), 
	    sub_job(bj)
	  {	// 	getoptions() should have been called elsewhere
	    as_initialize();
	    this->catalog = clg;
	    ss = 0;
	    int	fin = 0;
	    switch (this->input_mode) {
	      case IM_PARA: fin = 1; 
	      case IM_ALTR: case IM_EVRY: case IM_GRUP: case IM_SCAN:
		this->input_ns = 2; break;
	      case IM_ADON: case IM_TREE:
		this->input_ns = 0; return;
	      default:
		this->input_ns = 1; break;
	    }
	    this->memb = new InFiles(this->catalog, argc, argv, 
		this->input_mode, this->calc_mode);
	    sql1 = new SeqLoader<seq_t>(this, fin);
	  }
	~AlnServer()
	  {
	    delete sql1; delete sql2;
	  }
	char*	restsq(char* str, int n, int nseq);
	void	setshuffle(int in, int nj, int wh);
	int	nextvars();
	int	localoption(int& argc, const char**& argv);
	int	auto_comp(bool multhr = true);
	int	autocomp(bool multhr = true) {
	    return (auto_comp(multhr)? 0: (this->interactive? 1: 2));
	}
	int	makemsa();
	void	menucomp();	// general purpose
	void	common_usage();
// application specific functions
	void	setparam(int level);
	void	menu_job();
	int	(*sub_job)(CalcServer<seq_t>* svr, seq_t** sqs, ThQueue<seq_t>*);
};

template <class seq_t>
class ShuffleServer : public CalcServer<seq_t> {
public:
	int	nextvars() {
	    int	no_read = 0;
	    if (this->idx_a++ >= this->njumble) return (0);
	    if (this->which & 1) {this->in_face[0]->rndseq(); ++no_read;}
	    if (this->input_ns == 2 && this->which & 2)
		{this->in_face[1]->rndseq(); ++no_read;}
	    return (no_read);
	}
	ShuffleServer(int nj, int wh, void* p, int inn, seq_t** sqs,
	    int	(*mj)(CalcServer<seq_t>*, seq_t** sqs, ThQueue<seq_t>*)) :
	    CalcServer<seq_t>(inn == 1? IM_SNGL: IM_ALTR, p, mj)
	  {
	   	this->njumble = nj;
		this->which = wh;
		vclear(this->in_face, 3);
		this->vararray = 0;
		this->idx_a = 0;
		*this->in_face[0] = *sqs[0];
		if (inn == 2) *this->in_face[1] = *sqs[1];
		if (this->which & 1) this->in_face[0]->rndseq();
	 	if (this->which & 2 && this->in_face[1]) this->in_face[1]->rndseq();
	  }
	~ShuffleServer() {}
};

template <class seq_t>
class SeqLoader : public VarLoader<seq_t> {
	AlnServer<seq_t>* svr;
	DbsDt*	dbs;
	DbsGrp*	dbs_grp;
	int	dbs_sub;
	long	dbs_pos;
public:
	SeqLoader(AlnServer<seq_t>* as, int fin = 0, int bs = 0, int cl = 0, int dbs_no = 0);
	~SeqLoader() {}
	InSt	at(int sid, seq_t** inface= 0);
	InSt	nextseq(seq_t** sd, bool readin = true);
};

template <class seq_t>
class MakeMsa {
	char	wrkstr[LINE_MAX];
	char*	wrkvar;
	char*	seqnam(Tnode* root);
	AlnServer<seq_t>* svr;
public:
	MakeMsa(AlnServer<seq_t>* sr);
	MakeMsa(AlnServer<seq_t>* sv, Tnode* root);
	~MakeMsa() {}
};

template <class seq_t>
SeqLoader<seq_t>::SeqLoader(AlnServer<seq_t>* as, int fin, int bs, int cl, int dbs_no)
	: VarLoader<seq_t>((CalcServer<seq_t>*) as, fin, bs, cl), svr(as), 
	  dbs(0), dbs_grp(0), dbs_sub(0), dbs_pos(0)
{
	while (dbs_no >= 0) if ((dbs = dbs_dt[dbs_no--])) break;

	if (svr->calc_mode == IM_SCAN) {
	    dbs_grp = dbs->dbsgrp;
	}
	this->reset();
}

template <class seq_t>
InSt SeqLoader<seq_t>::nextseq(seq_t** sd, bool readin)
{
	if (svr->vararray) {
	    if (this->var_no >= this->ceilno) return (IS_END);
	    *sd = svr->vararray[this->var_no];	// 
	} else if (svr->msd && !svr->msd->empty()) {	// extract from msd msa
	    (*sd)->sid = this->var_no;
	    if (svr->ss) {
		if (this->var_no == svr->ss->num) return (IS_END);
		if (!(*sd)->inex.est)
		    svr->msd->extseq(*sd, svr->ss->group[this->var_no], CPY_SEQ | CPY_LBL);
	    } else {
		if (this->var_no == svr->msd->many) return (IS_END);
		if (!(*sd)->inex.est) {
		    int	grp[2] =  {this->var_no, -1};
		    svr->msd->extseq(*sd, grp, CPY_ALL);
		}
	    }
	} else if (dbs && dbs_grp) {	// read directly from database
	    if (dbs_pos == dbs_grp[1].recnbr) {++dbs_sub; ++dbs_grp;}
	    if (dbs_sub == dbs->numgrp) return (IS_END);
	    (*sd)->read_dbseq(dbs, dbs_pos++);
	} else {	// read from file
	    const char* attr = "";
	    const char*	attr2 = 0;
	    bool	first = false;
	    while (true) {
		while (!this->active_file()) {
		    if (!(this->mname && *this->mname)) return (IS_END);
		    if (!**this->mname) {
			if (svr->input_mode == IM_GRUP) return (IS_END);
			svr->setgrp2();
			++this->mname;
		    } else {
			attr = cdr(*this->mname);
			if (svr->input_mode == IM_MULT) attr2 = "M";
			first = true;
			if (**this->mname == DBSID) break;
#if USE_ZLIB
			this->fd = (*sd)->openseq(*this->mname, &this->gzfd);
#else
			this->fd = (*sd)->openseq(*this->mname);
#endif
			if (this->fixedin != 1) {
			    if (!this->active_file())
				prompt("%s not found !\n", *this->mname);
			    ++this->mname;
			} else if (!this->active_file()) return (IS_END);
		    }
		}
		if (*this->mname && **this->mname == DBSID) {	// Get from Database File */
		    Seq*	sdb = (*sd)->getdbseq(dbs, *this->mname, -1, readin);
		    if (this->fixedin != 1) ++this->mname;
		    else if (!sdb) return (IS_END);
		    break;
		} else {
		    if (this->fd && (*sd)->fgetseq(this->fd, attr, attr2)) break;
#if USE_ZLIB
		    if (this->gzfd && (*sd)->fgetseq(this->gzfd, attr, attr2)) break;
#endif
		    this->close_file();
		    if (this->fixedin == 1) return (IS_END);
		    if (first) return (IS_ERR);			// empty or absent
		}
	    }
	}
	++this->var_no;
	return ((*sd)->many? IS_OK: IS_ERR);
}

template <class seq_t>
InSt SeqLoader<seq_t>::at(int sid, seq_t** inface)
{
	if (dbs) {	// reset dbs
	    dbs_grp = dbs->dbsgrp;
	    dbs_pos = dbs_sub = 0;
	} if (svr->msd) {
	    this->var_no = sid;
	} else {
	    if (sid < this->var_no) {
		this->var_no = 0;
		this->mname = this->members;
		this->close_file();
	    }
	    while (this->var_no < sid)
		if (nextseq(inface) == IS_END) return (IS_END);
	}
	return (IS_OK);
}

template <class seq_t>
int AlnServer<seq_t>::readargs()
{
	FILE*	fd = ftable.fopen(AlnParam, "r");
	if (!fd) return (ERROR);

	char	argstr[MAX_ARGTEXT];
	char*	pa = argstr;
	const	char*	argv[MAX_ARGS];
	int	argc = 0;
	char	str[MAXL];
	int	quote = 0;
	int	dquote = 0;
	argv[argc++] = 0;
	while (char* ps = fgets(str, MAXL, fd)) {
	    for (int spc = 1; *ps; ++ps) {
		if (*ps == '#') break;
		if (*ps == '\'') quote = 1 - quote;
		if (*ps == '\"') dquote = 1 - dquote;
		if (isspace(*ps) && !quote && !dquote) {
		    if (!spc) {
			*pa++ = '\0';
			spc = 1;
		    }
		} else {
		    if (argc >= MAX_ARGS) goto readend;
		    if (spc) {
			argv[argc++] = pa;
			spc = 0;
		    }
		    *pa++ = *ps;
		}
	    }
	}
readend:
	fclose(fd);
const	char** av = argv;
	getoptions(argc, av);
	return (OK);
}

template <class seq_t>
void AlnServer<seq_t>::common_usage()
{
	fputs("Usage:\n", stderr);
	fprintf(stderr, "\t%s [-option_list] [seq_file1] ... \n", 
		progname? progname: "Prog");
	fputs("Common Options:\n", stderr);
	fputs("\t-A#\t;Algorithm\n", stderr);
	fputs("\t-d \"Database used\"\n", stderr);
	fputs("\t-D \"Database scanned\"\n", stderr);
	fputs("\t-F[5-9|C|L-O]\t; output format\n", stderr);
	fputs("\t-H#\t;Threshold\n", stderr);
	fputs("\t-i[a|e|f|g|i|l|p]:\tInput mode\n", stderr);
	fputs("\t-K [d|p|r|t];\t DNA|Protein|RNA|TRCODE\n", stderr);
	fputs("\t-l#\t;# of residues per line\n", stderr);
	fputs("\t-n \t;Output mode = 0\n", stderr);
	fputs("\t-O#\t;Output mode\n", stderr);
	fputs("\t-o \"Output file name\"\n", stderr);
	fputs("\t-p[a|b|..]\t;Modulate output form\n", stderr);
	fputs("\t-q#\"NSize of queue\"\n", stderr);
	fputs("\t-s \"Directory of sequence files\"\n", stderr);
	fputs("\t-t#\"Number of threads\"\n", stderr);
	fputs("\t-T \"Directory of tables\"\n", stderr);
	usage();
	exit (1);
}

template <class seq_t>
void AlnServer<seq_t>::getoptions(int& argc, const char**& argv)
{
	DbsDt**	dbs = dbs_dt;
static	int	visit = 0;

	if (visit++ == 0) {
	    topath(tmpdir, ".");
	    progname = *argv;
	}
	for ( ; --argc > 0 && **++argv == OPTCHAR; ) {
	  const	char*	opt= argv[0] + 1;
	  int	c = *opt;
	  if (!c) break;
	  const	char*	val = argv[0] + 2;
	  int		k = 0;
	  if (!localoption(argc, argv)) {
	    switch (c) {
	      case '?': case 'h': common_usage();
	      case 'A': 
		if ((val = getarg(argc, argv, true)))
		    algmode.alg = atoi(val);
		else	algmode.bnd = 0;
		break;
	      case 'C': 
		if ((val = getarg(argc, argv, true)))
		    initcodon(atoi(val));
		break;
	      case 'D':
		this->calc_mode = IM_SCAN;
	      case 'd':
		if ((val = getarg(argc, argv))) {
		    delete *dbs;
		    *dbs = new DbsDt(val);
		}
		break;
	      case 'F':
		if (('a' <= *val && *val <= 'z') ||
		    ('0' <= *val && *val <= '9'))
		    setprmode(*val, SILENT, SILENT); else
		if (*val == 'C' || ('K' <= *val && *val <= 'O'))
		    setprmode(SILENT, *val, SILENT); else
		if ('A' <= *val && *val <= 'Z')
		    setprmode(SILENT, SILENT, *val);
		break;
	      case 'H':
		if ((val = getarg(argc, argv, true)))
		     setthr(atof(val));
		else	setthr((double) (INT_MIN + 3));
		break;
	      case 'i': 		// input order
		switch (tolower(*val)) {
		  case 'a': case 'j':	// every juxt. pairs
		    this->calc_mode = IM_ALTR; break;
		  case 'e':		// lower-left combinaton
		    this->calc_mode = IM_EVRY; break;
		  case 'f':		// first vs others
		    this->calc_mode = IM_FvsO; break;
		  case 'g':		// group vs group, 1 file
		    this->calc_mode = IM_GRUP; break;
		  case 'i':		// self comparison
		    this->calc_mode = IM_SELF; break;
		  case 'l':		// last vs others
		    this->calc_mode = IM_OvsL; break;
		  case 'p':		// parallel, 2 files
		    this->calc_mode = IM_PARA; break;
		  case 'q':		// interactive
		    this->interactive = true; break;
		  case ':': this->catalog = val + 1;
		  default:
		    this->calc_mode = IM_SNGL; break;
		}
		if (*val && *val != ':' && val[1] == ':')
		    this->catalog = val[2]? val + 2: CATALOG;
		break;
	      case 'K':
		if ((val = getarg(argc, argv))) setdefmolc(*val);
		break;
	      case 'l':
		if ((val = getarg(argc, argv, true))) setlpw(atoi(val));
		break;
	      case 'm':
		if ((val = getarg(argc, argv, false))) {
		    if ((opt = strchr(val, ':')))
			k = atoi(opt + 1) - 1;
		    if (0 <= k && k < max_simmtxes) mdm_file[k] = val;
		}
		break;
	      case 'n':
		algmode.nsa = 0;
		algmode.rng = (*val && *val != '0');
		break;
	      case 'O':
		if ((val = getarg(argc, argv, true)))
		    algmode.nsa = atoi(val);
		break;
	      case 'o':
		OutPrm.out_file = ((val = getarg(argc, argv)))? val: "";
		break;
	      case 'p':
		switch (*val) {
		  case 'a': 
			if (isdigit(val[1])) ++val;
			else if (!val[1] && isdigit(argv[1][0])) {
			    val = *++argv; --argc;
			} else val = 0;
			if (val) polyA.setthr(val);
			else	polyA.setthr(def_polya_thr);
			break;
		  case 'd': OutPrm.descrp = 1; break;
		  case 'D': OutPrm.debug = 1; break;
		  case 'e': OutPrm.trimend = !OutPrm.trimend; break;
		  case 'f': OutPrm.deflbl = 1; break;
		  case 'g': OutPrm.taxoncode = 3; break;	// genus
		  case 'h': OutPrm.ColorEij = 2; break;
		  case 'i': OutPrm.ColorEij = 1; break;
		  case 'j': OutPrm.spjinf = 1 - OutPrm.spjinf; break;
		  case 'l': logfile = logfn; break;
		  case 'n': OutPrm.overwrite = 2; break;
		  case 'o': OutPrm.overwrite = 1; break;
		  case 'O': OutPrm.olrsum = 1; break;
		  case 'p': OutPrm.taxoncode = 5; break;	// phylum
		  case 'q': setprompt(false, false); break;
		  case 'Q': setprompt(false, true);  break;
		  case 's': OutPrm.sortodr = isdigit(val[1])? atoi(val+1): 1; break;
		  case 't': OutPrm.deflbl = 2; break;
		  case 'T': OutPrm.supTcodon = 1; break;
		  case 'w': algmode.thr = 0; break;
		  case 'W': OutPrm.printweight = 1; break;
		  case 'x': OutPrm.supself = 1; break;
		  case 'J': case 'm': case 'u': case 'v':
			setexprm_z(argc, argv); break;
		  default: break;
		}
		break;
#if M_THREAD
	      case 'q':
		if ((val = getarg(argc, argv, true)))
		    max_queue_num = atoi(val);
		else	max_queue_num = 0;
		break;
#endif
	      case 'Q':
		if ((val = getarg(argc, argv, true)))
		    algmode.qck = atoi(val);
		else	algmode.qck = 1;
		break;
	      case 'R':
		if ((val = getarg(argc, argv, true))) {
		    const char*	s2 = strchr(val, '.');
		    this->njumble = atoi(val);
		    if (s2) this->which = atoi(++s2);
		}
		break;
	      case 's':
		if ((val = getarg(argc, argv, false)))
		    setdfn(val);
		break;
	      case 'S': 
		if (*val == '-')	algmode.mns = 2;
		else if (isdigit(*val)) algmode.mns = atoi(val);
		else if (*val)	algmode.mns = 3;
		else if (argc > 1 && isdigit(argv[1][0]))
		    {algmode.mns = atoi(*++argv); --argc;}
		else algmode.mns = 0;
		break;
#if M_THREAD
	      case 't':
	        thread_num = (val = getarg(argc, argv, true))?
		    atoi(val): -1;
		break;
#endif
	      case 'T':
		if ((val = getarg(argc, argv, false)))
		    ftable.setpath(val, gnm2tab);
		readargs();
		break;
	      case 'W':
		if ((val = getarg(argc, argv, false)))
		    topath(tmpdir, val);
		break;
	      case 'u': case 'v': case 'w': 
		readalprm(argc, argv); break;
	      case 'y': readalprm(argc, argv, 2); break;
	      case 'z': setexprm_z(argc, argv); break;
	      default: break;
	    }
	  }
	  if (*dbs && (dbs + 1 < dbs_dt + MAX_DBS)) ++dbs;
	}
}

template <class seq_t>
void AlnServer<seq_t>::setshuffle(int ns, int jmbl, int wch)
{
	this->input_ns = ns;
again:
	if (jmbl == QUERY) promptin("# of trials (%d) : ", &this->njumble);
	else if (jmbl >= 0) this->njumble = jmbl;
	if (this->njumble) {
	    int	mask = (2 << (ns - 1)) - 1;
	    if (wch == 0) this->which = mask;
	    else if (wch > 0) this->which = wch & mask;
	    else if (wch == QUERY) {
		promptin("Which sequence is randomized ? [1-%d] : ", &this->which);
		if (!(this->which &= mask)) goto again;
		if (tolower(progetc("Use new random series ? [n/y] : ")) == 'y')
		    srand(rand());
	    }
	}
}

template <class seq_t>
int AlnServer<seq_t>::serialJob()
{
	int	nn = 0;
	do {
	    if (this->njumble && sub_job) {
		ShuffleServer<seq_t> ssvr(this->njumble, this->which, this->prm, 
		    this->input_ns, this->in_face, sub_job);
		ssvr.auto_comp();
		sum2avsd(ssvr.avsd, ssvr.njumble);
		vcopy(this->avsd, ssvr.avsd, 2);
	    }
	    if (this->main_job(this, this->in_face, 0) == OK) ++nn;
	} while (nextvars());
	return (nn);
}

template <class seq_t>
int AlnServer<seq_t>::auto_comp(bool multhr)
{
	if (this->vararray) return CalcServer<seq_t>::auto_comp(multhr);

	seq_t**	ins_a = this->in_face;
	seq_t**	ins_b = this->in_face + 1;
	InSt	inst_a = IS_ERR;
	InSt	inst_b = IS_OK;
	bool	readin = true;

	do {
	 switch (this->input_mode) {
	  case IM_SNGL: readin = this->jobcode != 'n';
	  case IM_MULT:
	    if ((inst_a = sql1->nextseq(ins_a, readin)) == IS_OK)
		(*ins_a)->sid = this->idx_a = 0;
	    setdefmolc((*ins_a)->inex.molc);
	    break;
	  case IM_EVRY:				// all combinations
	    if (inst_a == IS_ERR && (inst_a = sql1->nextseq(ins_a)) == IS_OK)
		(*ins_a)->sid = this->idx_a = 0;
	    if (inst_a == IS_OK) {
		if (!sql2) {
		    sql2 = new SeqLoader<seq_t>(this);
		    sql2->at(sql1->var_no, ins_b);
		}
		if ((inst_b = sql2->nextseq(ins_b)) == IS_OK)
		   (*ins_b)->sid = this->idx_b = 1;
	    }
	    setdefmolc((*ins_a)->inex.molc);
	    break;
	  case IM_FvsO:				// 1st vs others
	    if (inst_a == IS_ERR && (inst_a = sql1->nextseq(ins_b)) == IS_OK)
		(*ins_b)->sid = this->idx_b = 0;
	    if (inst_a == IS_OK) {
		if ((inst_b = sql1->nextseq(ins_a)) == IS_OK)
		    (*ins_a)->sid = this->idx_a = 1;
	    }
	    break;
	  case IM_OvsL:
	    this->idx_b = 0;
	    while ((inst_a = sql1->nextseq(ins_b)) != IS_END) ++this->idx_b;
	    if (!this->idx_b) return (0);
	    sql1->at((*ins_b)->sid = --this->idx_b, ins_b);
	    sql1->nextseq(ins_b);
	    sql1->reset();
	    if ((inst_a = sql1->nextseq(ins_a)) == IS_OK)
		(*ins_a)->sid = this->idx_a = 0;
	    break;
	  case IM_SCAN:
	    if (inst_a == IS_ERR && (inst_a = sql1->nextseq(ins_a)) == IS_OK)
		(*ins_a)->sid = this->idx_a = 0;
	    if (inst_a == IS_OK) {
		if (*ins_b) {
		    if (!sql2) sql2 = new SeqLoader<seq_t>(this, 0, 0, 1);
		    if ((inst_b = sql2->nextseq(ins_b)) == IS_OK)
			(*ins_b)->sid = this->idx_b = 0;
		}
	    }
	    setdefmolc((*ins_a)->inex.molc);
	    break;
	  case IM_GRUP:				// inter-groups
	    if (inst_a == IS_ERR && (inst_a = sql1->nextseq(ins_a)) == IS_OK)
		(*ins_a)->sid = this->idx_a = 0;
	    if (inst_a == IS_OK) {
		if (!sql2) sql2 = new SeqLoader<seq_t>(this, 0, this->getgrp22(), 1);
		if ((inst_b = sql2->nextseq(ins_b)) == IS_OK)
		    (*ins_b)->sid = this->idx_b = sql2->var_no;
	    }
	    setdefmolc((*ins_a)->inex.molc);
	    break;
	  case IM_ALTR:
	    if ((inst_a = sql1->nextseq(ins_a)) == IS_OK)
		(*ins_a)->sid = this->idx_a = 0;
	    if ((inst_b = sql1->nextseq(ins_b)) == IS_OK)
		(*ins_b)->sid = this->idx_b = 1;
	    break;
	  case IM_PARA:				// parallel input
	    if ((inst_a = sql1->nextseq(ins_a)) == IS_OK)
		(*ins_a)->sid = this->idx_a = 0;
	    if (!sql2) sql2 = new SeqLoader<seq_t>(this, 2, 0, 1);
	    if ((inst_b = sql2->nextseq(ins_b)) == IS_OK)
		(*ins_b)->sid = this->idx_b = 0;
	    break;
	  default:
	    fatal("input mode %d not supported !\n", this->input_mode);
	 }
	 if (inst_a == IS_END || inst_b == IS_END) return(0);
	} while (inst_a == IS_ERR || inst_b == IS_ERR);
	int	nprocessed = 0;
	if (this->setup_job) this->setup_job(this);
#if M_THREAD
	if (multhr && thread_num && !this->njumble) 
	    nprocessed  = this->MasterWorker(); else
#endif
	nprocessed = serialJob();
	if (this->cleanup_job) this->cleanup_job(this);
	return (nprocessed);
}

template <class seq_t>
int AlnServer<seq_t>::nextvars()
{
	if (this->vararray) return CalcServer<seq_t>::nextvars();

	seq_t**	ins_a = this->in_face;
	seq_t**	ins_b = this->in_face + 1;
	InSt	inst_a = IS_OK;
	InSt	inst_b = IS_OK;
	int	no_read = 0;

	do {
	 switch (this->input_mode) {
	  case IM_SNGL:			// single sequence
	  case IM_MULT:			// multiple sequence
	    inst_a = sql1->nextseq(ins_a, this->jobcode != 'n');
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++this->idx_a;
	    }
	    break;
	  case IM_ALTR:			// alternating pairs
	    inst_a = sql1->nextseq(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = this->idx_a += 2;
	    } else if (inst_a == IS_END) break;
	    inst_b = sql1->nextseq(ins_b);
	    if (inst_b == IS_OK) {
		++no_read; (*ins_b)->sid = this->idx_b += 2;
	    }
	    break;
	  case IM_PARA:			// parallel inputs
	    inst_a = sql1->nextseq(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++this->idx_a;
	    } else if (inst_a == IS_END) break;
	    inst_b = sql2->nextseq(ins_b);
	    if (inst_b == IS_OK) {
		++no_read; (*ins_b)->sid = ++this->idx_b;
	    }
	    break;
	  case IM_EVRY:			// all combinations
	    if (sql1->var_no + 1 == sql2->var_no) {
		inst_b = sql2->nextseq(ins_b);
		if (inst_b == IS_OK) {
		    ++no_read; ++this->idx_b;
		} else if (inst_b == IS_END) break;
		sql1->reset(); this->idx_a = -1;
		if (inst_b == IS_ERR) break;
	    }
	    inst_a = sql1->nextseq(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++this->idx_a;
		(*ins_b)->sid = this->idx_b;	// reset sid each time
	    }
	    break;
	  case IM_FvsO:			// 1st vs others
	    inst_a = sql1->nextseq(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++this->idx_a;
	    }
	    break;
	  case IM_OvsL:
	    if ((this->idx_a + 1) >= this->idx_b) break;
	    inst_a = sql1->nextseq(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++this->idx_a;
	    }
	    break;
	  case IM_SCAN:			// a = dbs, b = mem
	    inst_a = sql1->nextseq(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++this->idx_a;
	    } else if (inst_a == IS_END && sql2) {
		sql1->reset(); this->idx_a = -1;// rewind
		inst_a = sql1->nextseq(ins_a);
		if (inst_a == IS_OK) {
		    ++no_read; (*ins_a)->sid = ++this->idx_a;
		}
		inst_b = sql2->nextseq(ins_b);
		if (inst_b == IS_OK) {
		    ++no_read; (*ins_b)->sid = ++this->idx_b;
		}
	    }
	    break;
	  case IM_GRUP:			// inter-groups
	    inst_b = sql2->nextseq(ins_b);
	    if (inst_b == IS_OK) {
		++no_read; (*ins_b)->sid = ++this->idx_b;
		(*ins_a)->sid = this->idx_a;	// reset sid each time
	    } else if (inst_b == IS_END) {
		inst_a = sql1->nextseq(ins_a);
		if (inst_a == IS_OK) {
		    ++no_read; (*ins_a)->sid = ++this->idx_a;
		} else if (inst_a == IS_END) break;
		sql2->reset();	// rewind
		inst_b = sql2->nextseq(ins_b);
		if (inst_b == IS_OK) {
		    ++no_read;
		    (*ins_b)->sid = this->idx_b = sql2->var_no;
		}
	    }
	    break;
	  case IM_SELF:			// self compariosn
	    inst_a = sql1->nextseq(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++this->idx_a;
		this->alias_of(*ins_b, *ins_a);
	    }
	    break;
	 }
	 if (inst_a == IS_END || inst_b == IS_END) return(0);
	} while (inst_a == IS_ERR || inst_b == IS_ERR);
	return (no_read);
}

/********************************************************************
*
*	Interactive interface
*
********************************************************************/

template <class seq_t>
char* AlnServer<seq_t>::restsq(char* str, int n, int nseq)
{
	if (n >= nseq) return (strcpy(str, "gq"));
	str[0] = '\0';
	int	slen = CONLINE;
	for (char* s = str; slen > 0 && ++n <= nseq; s = str + strlen(str))
	    slen -= snprintf(s, slen, "%d ", n);
	return (str);
}

template <class seq_t>
char* AlnServer<seq_t>::getaseq(char* ps) 
{
	while (*ps) {
	    while (*ps && isspace(*ps)) ++ps;
	    if (!isdigit(*ps)) break;
	    char*	qs = ps;
	    for ( ; *qs && !isspace(*qs); ++qs) ;
	    if (*qs) *qs++ = '\0';
	    inputseq(this->in_face, ps);
	    ps = qs;
	}
	return (ps);
}

template <class seq_t>
void AlnServer<seq_t>::menucomp()
{
	char    cmd[CONLINE];
	int     nreadseq = sql1->var_no? 1: 0;
	if (sql2 && sql2->var_no) ++nreadseq;
	restsq(cmd, nreadseq, this->input_ns? this->input_ns: 1);
	setprompt(true, false);
	int	ngo = 0;
	while (true) {
	    char*       ps = getaseq(cmd);
	    while (*ps) {
		switch (*ps++) {
		    case 'p': // Set Parameter Grp. 1
			setparam(1);
			break;
		    case 'P': // Set Parameter Grp. 2
			setparam(2);
			break;
		    case 'g': // Calculate
			if (!ngo++ && this->setup_job)
			    this->setup_job(this);
			this->main_job(this, this->in_face, 0); break;
		    case 'q': // quit
			if (this->cleanup_job) this->cleanup_job(this);
			return;
		    default:  break;
		}
	    }
	    ps = progets(cmd, "\n#(+)[read/app]/p[aram]/g[o]/q[uit] : ");
	}
}

extern	int	testInputSeq(CalcServer<Seq>* svr, Seq* sqs[], ThQueue<Seq>* q);
extern	InputSeqTest	zero_InputSeqTest;

#endif	// _AUTOCOMP_
