/*****************************************************************************
*
*	Execution of a program with a combination of a series of data
*
*	Subcommand: 'e' all; 		'f' single scan; 	'g' group
*	Subcommand: 'i' internal; 	'j' juxtapose		'k' addon
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

#ifndef  _CALCSERV_
#define  _CALCSERV_

#define QDEBUG	0

#include "cmn.h"

extern	void	usage();

static	const	char*	CATALOG = "catalog";
static	const	char	HARD_DELIM = '*';
static	const	char	erropen[] = "%s can't open!";
static	const	char*	progname = 0;
static	const	char	no_group_2[] = "Specify 2nd group !\n";

class	DbsDt;

#if M_THREAD

#include <pthread.h>
#include <unistd.h>
#include <sched.h>

template <class var_t>
class ThQueue {
	int	rp, wp;
	int	remain;
	int	done;
	pthread_cond_t	not_full;
	pthread_cond_t	not_empty;
public:
	int	qsize;
	var_t**	varque;
	pthread_mutex_t	mutex;
	ThQueue(int qs);
	~ThQueue() {
	    for (int i = 0; i < qsize; ++i) delete varque[i];
	    delete[] varque;
	}
	void	fill();
	void	enqueue(var_t** fsd, int n = 1);
	void	dequeue(var_t** fsd, int n = 1);
};

template <class var_t>
void ThQueue<var_t>::enqueue(var_t** vars, int n)
{
	pthread_mutex_lock(&mutex);
	while (remain == qsize)
	    pthread_cond_wait(&not_full, &mutex);
	while (--n >= 0) {
#if QDEBUG
	    if (vars[n]) fprintf(stderr, "e%d: %d\n", n, vars[n]->sid);
	    else	 fprintf(stderr, "e%d: z\n", n);
#endif // QDEBUG
	    swap(varque[wp], vars[n]);
	    ++wp; ++remain;
	    if (wp == qsize) wp = 0;
	}
	pthread_cond_signal(&not_empty);
	pthread_mutex_unlock(&mutex);
}

template <class var_t>
void ThQueue<var_t>::dequeue(var_t** vars, int n)
{
	pthread_mutex_lock(&mutex);
	while (remain == 0)
	     pthread_cond_wait(&not_empty, &mutex);
	while (--n >= 0) {
	    swap(vars[n], varque[rp]);
#if QDEBUG
	    if (vars[n]) fprintf(stderr, "d%d: %d\n", n, vars[n]->sid);
	    else	 fprintf(stderr, "d%d: z\n", n);
#endif // QDEBUG
	    ++rp; --remain;
	    if (rp == qsize) rp = 0;
	}
	pthread_cond_signal(&not_full);
	pthread_mutex_unlock(&mutex);
}

template <class var_t>
void m_thread_Lock(ThQueue<var_t>* q)
{
	if (q) pthread_mutex_lock(&q->mutex);
}

template <class var_t>
void m_thread_Unlock(ThQueue<var_t>* q)
{
	if (q) pthread_mutex_unlock(&q->mutex);
}

template <class var_t>
ThQueue<var_t>::ThQueue(int qs) : qsize(qs)
{
	varque = new var_t*[qsize];
	vclear(varque, qsize);
	rp = wp = remain = 0;
	pthread_mutex_init(&mutex, 0);
	pthread_cond_init(&not_full, 0);
	pthread_cond_init(&not_empty, 0);
}

template <class var_t>
void ThQueue<var_t>::fill()
{
	var_t**	vars = varque;
	var_t**	vart = varque + qsize;
	while (vars < vart) *vars++ = new var_t;
}

#else	// M_THREAD

template <class var_t>
class ThQueue {
	var_t dummy;
};

template <class var_t>
void m_thread_Lock(ThQueue<var_t>* q) {}

template <class var_t>
void m_thread_Unlock(ThQueue<var_t>* q) {}

#endif	// M_THREAD

struct InFiles {
	int	membersize = 0;
	char**	name = 0;
	char**	grp2 = 0;
	int	idx_g2 = 0;
	InFiles(const char* fname, int ac, const char** arv,
	    int imode, int cmode = 0) {
	    int	ctgnmbr = 0;	// members in catalog
	    int	argsize = 0;	// total lengths of argx
	    int	ctgsize = 0;	// size of catalog file
	    const	char**	av = arv;
	    for (int i = 0; i++ < ac; ++av) argsize += strlen(*av) + 1;
	    FILE*	fd = fname? fopen(fname, "r"): 0;
	    if (fd) {
		char	str[MAXL];
		while (fgets(str, MAXL, fd)) {
		    if (*str == _LCOMM) flush_line(fd);
		    else	++ctgnmbr;
		}
		ctgsize = ftell(fd);
	    }
	    membersize = ac + ctgnmbr;
	    if (!membersize) return;	// no input seq
	    char**	ptr = name = new char*[membersize + 2];
	    char*	wk = *name = new char[argsize + ctgsize + 2];
	    if (fd) {			// read catalog
		rewind(fd);
		if (fread(*name, ctgsize, 1, fd) != 1)
		    fatal("Can't read %s !\n", fname);
		if (wk[ctgsize - 1] != '\n') wk[ctgsize] = '\n';
		for (int i = 0; i++ < ctgnmbr; ++ptr) {
		    while (*wk == _LCOMM)
			while (*wk++ != '\n') ; // skip comment
		    for (*ptr = wk; *wk != '\n'; ++wk) ;
		    *wk++ = '\0';
       		    if (!**ptr) grp2 = ptr + 1;	// blank line
		    else if (**ptr == HARD_DELIM) {
			**ptr = '\0';
			grp2 = ptr + 1;
			break;		// hard delimiter
		    }
		}
		fclose(fd);
	    }
	    av = arv;
	    for (int i = 0; i++ < ac; ) {	// argments
		if (!grp2 && ptr > name && imode != IM_SNGL) {
		    *ptr = wk - 1;		// delimiter
		    grp2 = ++ptr;		// 2nd group members
		}
		strcpy(wk, *av);
		*ptr++ = wk;
		wk += strlen(*av++) + 1;
	    }
	    *ptr = 0;
	    if (cmode == IM_GRUP) {
		if (!(grp2 && *grp2)) fatal("2nd group is missing !\n");
		idx_g2 = grp2 - name - 1;
	    } else	grp2 = 0;
	}
	~InFiles() {
	    if (name) {
		delete[] *name;
		delete[] name;
	    }
	}
};

template <class var_t> class VarLoader;

struct Minvar {
	int	sid;
	bool	vrtl;
	int	fget(FILE* fd, const char* fn = 0) {return (0);}
};

template <class var_t>
void clearvars(var_t** vars, int n, bool del = true)
{
	while (--n >= 0) {
	    if (del && *vars && !(*vars)->vrtl) delete *vars;
	    *vars++ = 0;
	}
}

// run 'main_job' for each member of vrars

template <class var_t>
class CalcServer {
protected:
	int	idx_a, idx_b, idx_g;
	bool	interactive;
#if M_THREAD
	int	MasterWorker();
#endif
	int	serialJob();
	void	getoptions(int& argc, const char**& argv); 
	void	initvars(var_t** vars, int n);
	void	common_usage();
public:
const	char*	catalog;
	InFiles*	memb;
	VarLoader<var_t>*	ldr1;
	VarLoader<var_t>*	ldr2;
	int	jobcode;
	int	njumble;
	int	which;
	double	avsd[2];
	var_t*	msd;
	var_t*	in_face[3];		// input interface
	int	input_ns;
	int	input_mode;
	int	calc_mode;
	void*	prm;
	var_t**	vararray;		// 
	int	vararray_no;
	int	(*main_job)(CalcServer<var_t>* svr, var_t** vars, 
			ThQueue<var_t>* q);
	int	(*optn_job)(int& ac, const char**& av);
	void	(*setup_job)(CalcServer<var_t>* svr);
	void	(*cleanup_job)(CalcServer<var_t>* svr);
	int	auto_comp(bool multhr = thread_num);
	int	autocomp(bool multhr = thread_num) {
	    return (auto_comp(multhr)? 0: (interactive? 1: 2));
	}
	void	alias_of(var_t* dst, var_t* src) {
	    *dst = *src;
	    dst->vrtl = true;
	}
	void	initialize() 
	  {
	    catalog = 0; msd = 0;
	    njumble = which = 0;
	    vclear(avsd, 2);
	    idx_a = idx_b = jobcode = 0;
	    memb = 0; ldr1 = ldr2 = 0;
	    interactive = false;
	    initvars(in_face, 3);
	  }
virtual	int	nextvars();
	void	setgrp2(int num2 = 0) {idx_g = num2? num2: idx_a + 1;}
	int	getgrp2() {return (idx_g);}
	int	getgrp22() {return (memb? memb->idx_g2: idx_g);}
	void	fpavsd(double val);
	void	set_catalog(const char* cg) {catalog = (cg && *cg)? cg: CATALOG;}
	int	memsize(VarLoader<var_t>* ldr = 0)
	  {
	    if (!ldr) return (memb? memb->membersize: vararray_no);
	    InSt	inst;
	    do {
		var_t	tmp;
		var_t*	ptmp = &tmp;
		inst = ldr->nextvar(&ptmp);
	    } while (inst != IS_END);
	    int	rv = ldr->var_no;
	    ldr->reset();
	    return (rv);
	  }
	int	calcsize(int num, int cm = IM_NONE) 
	  {
	    if (cm == IM_NONE) cm = calc_mode;
	    switch (cm) {
		case IM_ALTR: case IM_PARA: num /= 2; break;
		case IM_EVRY: num = elem(num, 0); break;
		case IM_GRUP: num = idx_g * (num - idx_g); break;
		default: break;
	    }
	    return (num);
	  }
	int	calcnbr(int a, int b)
	  {
	    int	nn = 0;
	    switch (calc_mode) {
		case IM_ALTR: nn = a / 2; break;
		case IM_EVRY: nn = elem(a, b); break;
		case IM_FvsO: nn = a - 1; break;
		case IM_GRUP: nn = a + idx_g * (b - idx_g); break;
		default: nn = a; break;
	    }
	    return (nn);
	  }
	void	change_job(int (*mj)(CalcServer<var_t>*, var_t**, ThQueue<var_t>*))
	  {
		if (ldr1) ldr1->reset();
		if (ldr2) ldr2->reset();
		main_job = mj;
	  }
	CalcServer(int im, void* p,
	    int		(*mj)(CalcServer<var_t>*, var_t**, ThQueue<var_t>*),
	    void	(*uj)(CalcServer<var_t>*) = 0,
	    void	(*cj)(CalcServer<var_t>*) = 0,
	    var_t** vars = 0, int nos = 0, int nos2 = 0, bool ldr = true) :
	    idx_g(nos2), input_mode(im), calc_mode(im), prm(p), 
	    vararray(vars), vararray_no(nos), 
	    main_job(mj), setup_job(uj), cleanup_job(cj)
	  {
	    initialize();
	    switch (input_mode) {
		case IM_GRUP: if (!idx_g) fatal(no_group_2);
		case IM_ALTR: case IM_EVRY: case IM_SCAN:
		case IM_PARA: 
		    input_ns = 2; break;
		case IM_ADON: case IM_TREE:
		    input_ns = 0; return;
		default:
		    input_ns = 1; break;
	    }
	    if (ldr) ldr1 = new VarLoader<var_t>(this, 0, 0, idx_g);
	  }
	CalcServer(int& argc, const char**& argv, int im, int defim, void* p = 0,
	    int		(*mj)(CalcServer<var_t>*, var_t**, ThQueue<var_t>*) = 0,
	    int 	(*oj)(int& ac, const char**& av) = 0,
	    void	(*uj)(CalcServer<var_t>*) = 0,
	    void	(*cj)(CalcServer<var_t>*) = 0,
	    int nos = 0, int nos2 = 0) :
	    idx_g(nos2), calc_mode(defim), prm(p), 
	    vararray(0), vararray_no(nos), 
	    main_job(mj), optn_job(oj), setup_job(uj), cleanup_job(cj)
	  {
	    initialize();
	    getoptions(argc, argv);
	    if (im) input_mode = im;
	    int	fin = 0;
	    switch (input_mode) {
		case IM_PARA: fin = 1; input_ns = 2; break;
		case IM_GRUP: if (!idx_g) fatal(no_group_2);
		case IM_ALTR: case IM_EVRY: input_ns = 2; break;
		default: input_ns = 1; break;
	    }
	    memb = new InFiles(catalog, argc, argv, input_mode, calc_mode);
	    ldr1 = new VarLoader<var_t>(this, fin, 0, idx_g);
	  }
virtual	~CalcServer()
	  {
	    delete memb; delete ldr1; delete ldr2;
	    clearvars(in_face, 3, !vararray);
	  }
};

template <class var_t>
class VarLoader {
protected:
	CalcServer<var_t>* svr;
	char**	members;	// list of file names
	char**	mname;		// working members
	int	fixedin;	// input from a fixed file
	int	baseno;		// > 0 for group2
	int	ceilno;
	const	char*	fname;
	FILE*	fd;
#if USE_ZLIB
	gzFile	gzfd;
#endif
	void	close_file() {
	    fname = 0;
	    if (fd) {fclose(fd); fd = 0;}
#if USE_ZLIB
	    if (gzfd) {fclose(gzfd); gzfd = 0;}
#endif
	}
public:
	int	var_no;
	void	reset();
	VarLoader(CalcServer<var_t>* cs, int fin = 0, int bs = 0, int cl = 0) 
	    : svr(cs), fixedin(fin), baseno(bs), ceilno(cl), fname(0), fd(0)
#if USE_ZLIB
		, gzfd(0)
#endif
	{
	    if (!svr->memb) mname = members = 0;
	    if (!ceilno) ceilno = svr->memsize();
	    reset();
	}
	~VarLoader() {
	    if (fd) fclose(fd);
#if USE_ZLIB
	    if (gzfd) fclose(gzfd);
#endif
	}
	bool	active_file() {
#if USE_ZLIB
	    return (fd || gzfd);
#else
	    return (fd);
#endif
	}
	InSt	at(int sid, var_t** inface = 0);
	InSt	nextvar(var_t** var);
};

#if M_THREAD

template <class var_t>
struct thread_arg_t {
	ThQueue<var_t>*	q;
	int	cpuid;
	var_t**	vars;
	DbsDt*	dbf;
	void*	svr;
};

template <class var_t>
struct mast_arg_t {
	ThQueue<var_t>*	q;
	int	tnum;
	void*	svr;
};

template <class svr_t, class var_t>
void* master_func(void* arg)
{
	mast_arg_t<var_t>*	targ = (mast_arg_t<var_t>*) arg;
	ThQueue<var_t>*	q = targ->q;
	svr_t*	svr = (svr_t*) targ->svr;
	var_t**	intfc = svr->in_face;
	int	insn = 2;
	int	save = svr->input_mode == IM_EVRY? 1: 0;

	do {
	    if (svr->input_mode == IM_EVRY || svr->input_mode == IM_GRUP ||
		(svr->input_mode == IM_SCAN && svr->input_ns == 2)) {
		if (insn == 2) {
		    if (svr->vararray) intfc[2] = intfc[save];
		    else	*intfc[2] = *intfc[save];
		} else {
		    if (svr->vararray) intfc[save] = intfc[2];
		    else	*intfc[save] = *intfc[2];
		}
	    }
	    q->enqueue(intfc, svr->input_ns);
	} while ((insn = svr->nextvars()));

// fill 0's to indicate end of inputs

	var_t*	nullvar[2] = {0, 0};
	for (int n = 0; n < targ->tnum; ++n) {
	    q->enqueue(nullvar, svr->input_ns);
	    clearvars(nullvar, svr->input_ns, !svr->vararray);
	}
	return (void*) 0;
}

template <class svr_t, class var_t>
static void* worker_func(void* arg)
{
	thread_arg_t<var_t>* targ = (thread_arg_t<var_t>*) arg;
	svr_t*	svr = (svr_t*) targ->svr;

#ifdef __CPU_SET
	cpu_set_t	mask;
	__CPU_ZERO(&mask);
	__CPU_SET(targ->cpuid, &mask);
	if (sched_setaffinity(0, sizeof(mask), &mask) == -1)
	    prompt("Warning: faild to set CPU affinity !\n");
#endif // __CPU_SET

	while (true) {
	    targ->q->dequeue(targ->vars, svr->input_ns);
	    if (!targ->vars[0]) break;
	    if (svr->input_ns == 2 && !targ->vars[1]) break;
	    svr->main_job(svr, targ->vars, targ->q);
	}
	return (void*) 0;
}

template <class var_t>
int CalcServer<var_t>::MasterWorker()
{
	cpu_num = sysconf(_SC_NPROCESSORS_CONF);

	if (thread_num < 0)	thread_num = cpu_num;
	if (max_queue_num <= 0)	max_queue_num = int(FACT_QUEUE * thread_num);
	max_queue_num = (max_queue_num + input_ns - 1) / input_ns * input_ns;
	pthread_t	master;
	pthread_t*	worker = new pthread_t[thread_num];
	thread_arg_t<var_t>*	targ = new thread_arg_t<var_t>[thread_num];
	int	nsvar = (input_mode == IM_FvsO || input_mode == IM_OvsL)? 2: input_ns;

	ThQueue<var_t>	q(max_queue_num);
	if (!vararray) q.fill();
	mast_arg_t<var_t>	maarg = {&q, thread_num, (void*) this};
	var_t**	vars = new var_t*[nsvar * thread_num];
	for (int n = 0; n < thread_num; ++n) {
	    targ[n].q = &q;
	    targ[n].cpuid = n % cpu_num;
	    targ[n].svr = (void*) this;
	    targ[n].vars = vars + nsvar * n;
	    initvars(targ[n].vars, input_ns);
	    if (input_mode == IM_FvsO || input_mode == IM_OvsL)
		targ[n].vars[1] = in_face[1];
	}
	pthread_create(&master, 0, master_func<CalcServer<var_t>, var_t>, (void*) &maarg);
	for (int n = 0; n < thread_num; ++n)
	    pthread_create(worker + n, 0, worker_func<CalcServer<var_t>, var_t>, 
		(void*) (targ + n));
	for (int n = 0; n < thread_num; ++n)
	    pthread_join(worker[n], 0);
	if (input_mode == IM_FvsO || input_mode == IM_OvsL)
	    for (int n = 1; n < 2 * thread_num; n += 2) vars[n] = 0;
	clearvars(vars, nsvar * thread_num, !vararray);
	clearvars(q.varque, q.qsize, !vararray);
	delete[] vars;
	delete[] targ;
	delete[] worker;
	return (thread_num);
}

#endif	//	M_THREAD

template <class var_t>
void VarLoader<var_t>::reset()
{
	if (fd) {fclose(fd); fd = 0;} var_no = 0;
#if USE_ZLIB
	if (gzfd) {fclose(gzfd); gzfd = 0;} var_no = 0;
#endif
	if (svr->memb) {
	    if (baseno && svr->memb->grp2) {
		members = svr->memb->grp2;
	    } else {
		 members = svr->memb->name;
		if (fixedin) members += fixedin - 1;
	    }
	    mname = members;
	} else {
	    var_no = baseno;
	}
}

template <class var_t>
InSt VarLoader<var_t>::nextvar(var_t** var)
{
	bool	first = false;
	if (svr->vararray) {
	    if (var_no >= ceilno) return (IS_END);
	    *var = svr->vararray[var_no];
	} else if (fd == stdin) {	// read from file
	    if ((*var)->fget(fd, 0) == EOF) return (IS_END);
	} else {
	    int	rv = EOF;
	    while (rv == EOF) {
		while (!fd) {
		    if (!(mname && *mname)) return (IS_END);
		    fname = *mname;
		    if (!**mname) {	// blank line
			if (svr->input_mode == IM_GRUP) return (IS_END);
			svr->setgrp2();
			++mname;
		    } else {
			fd = fopen(fname, "r");
			if (fixedin != 1) {
			    ++mname;
			    first = true;
			    if (!fd) prompt("%s not found !\n", fname);
			} else if (!fd) return (IS_END);
		    }
		}
		if (fd && fname) {
		    rv = (*var)->fget(fd, fname);
		    if (feof(fd)) close_file();
		}
		if (rv == EOF) {
		    close_file();
		    if (fixedin == 1) return (IS_END);
		    if (first) return (IS_ERR);		// empty or absent
		} else if (rv == IGNORE) {
		    return (IS_ERR);
		}
	    }
	}
	++var_no;
	return (IS_OK);
}

template <class var_t>
InSt VarLoader<var_t>::at(int sid, var_t** inface)
{
	if (inface) {
	   if (sid < var_no) {
		var_no = 0;
		mname = members;
		close_file();
	   }
	   while (var_no < sid)
		if (nextvar(inface) == IS_END) return (IS_END);
	} else
	    var_no = sid + baseno;
	return (IS_OK);
}

template <class var_t>
void CalcServer<var_t>::common_usage()
{
	fputs("Usage:\n", stderr);
	fprintf(stderr, "\t%s [-option_list] [seq_file1] ... \n", 
		progname? progname: "Prog");
	fputs("Common Options:\n", stderr);
	fputs("\t-A#\t;Algorithm\n", stderr);
	fputs("\t-O#\t;Output mode\n", stderr);
	fputs("\t-i[a|e|f|g|i|l|p]:\tInput mode\n", stderr);
	fputs("\t-q#\"NSize of queue\"\n", stderr);
	fputs("\t-t#\"Number of threads\"\n", stderr);
	usage();
	exit (1);
}

template <class var_t>
void CalcServer<var_t>::getoptions(int& argc, const char**& argv) 
{
static	int	visit = 0;

	if (visit++ == 0) progname = *argv;
	for ( ; --argc > 0 && **++argv == OPTCHAR; ) {
	  const	char*	opt= argv[0] + 1;
	  int	c = *opt;
	  if (!c) break;
	  if (!(optn_job && optn_job(argc, argv))) {
	    const	char*	val = argv[0] + 2;
	    switch (c) {
	      case '?': case 'h': common_usage();
	      case 'A': 
		if ((val = getarg(argc, argv, true)))
		    algmode.alg = atoi(val);
		else	algmode.bnd = 0;
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
	      case 'O':
		if ((val = getarg(argc, argv, true)))
		    algmode.nsa = atoi(val);
		break;
	      case 'q':
		if ((val = getarg(argc, argv, true)))
		    max_queue_num = atoi(val);
		else	max_queue_num = 0;
		break;
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
	      case 't':
	        thread_num = (val = getarg(argc, argv, true))?
		    atoi(val): -1;
		break;
	      default: break;
	    }
	  }
	}
}

template <class var_t>
void  CalcServer<var_t>::initvars(var_t** vars, int n)
{
	while (--n >= 0)
	    *vars++ = vararray? 0: new var_t();
}

template <typename X>
void sum2avsd(X* avsd, int num)
{
	if (num < 1) return;
	if (num < 2) avsd[1] = 0.;
	else {
	    avsd[0] /= num;
	    avsd[1] = sqrt((avsd[1] - avsd[0] * avsd[0]
		 * num) / (num - 1));
	}
}

template <class var_t>
void CalcServer<var_t>::fpavsd(double val)
{
	if (!njumble) return;

	double	dev = (avsd[1] > 0.)? (val - avsd[0]) / avsd[1]: 0;
	if (minmax == MINIMUM) dev = -dev;
	fprintf(out_fd, 
	    "Dev = %6.2lf  AV = %7.2lf  SD = %7.2lf   (%d jumbles)\n",
	    dev, avsd[0], avsd[1], njumble);
}

template <class var_t>
int CalcServer<var_t>::serialJob()
{
	int	nn = 0;
	do {
	    if (main_job(this, this->in_face, 0) == OK) ++nn;
	} while (nextvars());
	return (nn);
}

template <class var_t>
int CalcServer<var_t>::auto_comp(bool multhr)
{
	var_t**	ins_a = in_face;
	var_t**	ins_b = in_face + 1;
	InSt	inst_a = IS_ERR;
	InSt	inst_b = IS_OK;

	if (ldr2) ldr2->reset();
	do {
	 switch (input_mode) {
	  case IM_SNGL:
	    if ((inst_a = ldr1->nextvar(ins_a)) == IS_OK)
		(*ins_a)->sid = idx_a = 0;
	    break;
	  case IM_EVRY:				// all combinations
	    if (inst_a == IS_ERR && (inst_a = ldr1->nextvar(ins_a)) == IS_OK)
		(*ins_a)->sid = idx_a = 0;
	    if (inst_a == IS_OK) {
		if (!ldr2) ldr2 = new VarLoader<var_t>(this);
		ldr2->at(ldr1->var_no, memb? ins_b: 0);
		if ((inst_b = ldr2->nextvar(ins_b)) == IS_OK)
		   (*ins_b)->sid = idx_b = 1;
	    }
	    break;
	  case IM_FvsO:				// 1st vs others
	    if (inst_a == IS_ERR && (inst_a = ldr1->nextvar(ins_b)) == IS_OK)
		(*ins_b)->sid = idx_b = 0;
	    if (inst_a == IS_OK) {
		if ((inst_b = ldr1->nextvar(ins_a)) == IS_OK)
		    (*ins_a)->sid = idx_a = 1;
	    }
	    break;
	  case IM_OvsL:
	    idx_b = memsize() - 1;
	    ldr1->at(idx_b, memb? ins_b: 0);
	    ldr1->nextvar(ins_b);
	    (*ins_b)->sid = idx_b;
	    ldr1->reset();
	    if ((inst_a = ldr1->nextvar(ins_a)) == IS_OK)
		(*ins_a)->sid = idx_a = 0;
	    break;
	  case IM_SCAN:
	    if (inst_a == IS_ERR && (inst_a = ldr1->nextvar(ins_a)) == IS_OK)
		(*ins_a)->sid = idx_a = 0;
	    if (inst_a == IS_OK) {
		if (*ins_b) {
		    if (!ldr2) ldr2 = new VarLoader<var_t>(this);
		    if ((inst_b = ldr2->nextvar(ins_b)) == IS_OK)
			(*ins_b)->sid = idx_b = 0;
		}
	    }
	    break;
	  case IM_GRUP:				// inter-groups
	    if (inst_a == IS_ERR && (inst_a = ldr1->nextvar(ins_a)) == IS_OK)
		(*ins_a)->sid = idx_a = 0;
	    if (inst_a == IS_OK) {
		if (!ldr2) ldr2 = new VarLoader<var_t>(this, 0, getgrp22());
		idx_b = ldr2->var_no;
		if ((inst_b = ldr2->nextvar(ins_b)) == IS_OK)
		    (*ins_b)->sid = idx_b;
	    }
	    break;
	  case IM_ALTR:
	    if ((inst_a = ldr1->nextvar(ins_a)) == IS_OK)
		(*ins_a)->sid = idx_a = 0;
	    if ((inst_b = ldr1->nextvar(ins_b)) == IS_OK)
		(*ins_b)->sid = idx_b = 1;
	    break;
	  case IM_PARA:				// parallel input
	    if ((inst_a = ldr1->nextvar(ins_a)) == IS_OK)
		(*ins_a)->sid = idx_a = 0;
	    if (!ldr2) ldr2 = new VarLoader<var_t>(this, 2);
	    if ((inst_b = ldr2->nextvar(ins_b)) == IS_OK)
		(*ins_b)->sid = idx_b = 0;
	    break;
	  default: break;
	 }
	 if (inst_a == IS_END || inst_b == IS_END) return(0);
	} while (inst_a == IS_ERR || inst_b == IS_ERR);
	int	nprocessed = 0;
	if (setup_job) setup_job(this);
#if M_THREAD
	if (multhr && thread_num)
	    nprocessed = MasterWorker(); else
#endif
	nprocessed = serialJob();
	if (cleanup_job) cleanup_job(this);
	return (nprocessed);
}

// return the number of successfully read sequences
// assume idx_a < idx_b except for IM_FvsO

template <class var_t>
int CalcServer<var_t>::nextvars()
{
	var_t**	ins_a = in_face;
	var_t**	ins_b = in_face + 1;
	InSt	inst_a = IS_OK;
	InSt	inst_b = IS_OK;
	int	no_read = 0;

	do {
	 switch (input_mode) {
	  case IM_SNGL:			// single input
	    inst_a = ldr1->nextvar(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++idx_a;
	    }
	    break;
	  case IM_ALTR:			// alternating pairs
	    inst_a = ldr1->nextvar(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = idx_a += 2;
	    } else if (inst_a == IS_END) break;
	    inst_b = ldr1->nextvar(ins_b);
	    if (inst_b == IS_OK) {
		++no_read; (*ins_b)->sid = idx_b += 2;
	    } 
	    break;
	  case IM_PARA:			// parallel inputs
	    inst_a = ldr1->nextvar(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++idx_a;
	    } else if (inst_a == IS_END) break;
	    inst_b = ldr2->nextvar(ins_b);
	    if (inst_b == IS_OK) {
		++no_read; (*ins_b)->sid = ++idx_b;
	    }
	    break;
	  case IM_EVRY:			// all combinations
	    if (ldr1->var_no + 1 == ldr2->var_no) {
		inst_b = ldr2->nextvar(ins_b);
		if (inst_b == IS_OK) {
		    ++no_read; ++idx_b;
		} else if (inst_b == IS_END) break;
		ldr1->reset(); idx_a = -1;
	    }
	    inst_a = ldr1->nextvar(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++idx_a;
	    }
	    break;
	  case IM_FvsO:			// 1st vs others
	    inst_a = ldr1->nextvar(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++idx_a;
	    }
	    break;
	  case IM_OvsL:
	    if ((idx_a + 1) >= idx_b) break;
	    inst_a = ldr1->nextvar(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++idx_a;
	    }
	    break;
	  case IM_SCAN:
	    inst_a = ldr1->nextvar(ins_a);
	    if (inst_a == IS_OK) {
		++no_read; (*ins_a)->sid = ++idx_a;
	    } else if (inst_a == IS_END && ldr2) {
		ldr1->reset(); idx_a = -1;// rewind
		inst_a = ldr1->nextvar(ins_a);
		if (inst_a == IS_OK) {
		    ++no_read; (*ins_a)->sid = ++idx_a;
		}
		inst_b = ldr2->nextvar(ins_b);
		if (inst_b == IS_OK) {
		    ++no_read; (*ins_b)->sid = ++idx_b;
		}
	    }
	    break;
	  case IM_GRUP:			// inter-groups
	    inst_b = ldr2->nextvar(ins_b);
	    if (inst_b == IS_OK) {
		++no_read; (*ins_b)->sid = ++idx_b;
		(*ins_a)->sid = idx_a;	// reset sid each time
	    } else if (inst_b == IS_END) {
		inst_a = ldr1->nextvar(ins_a);
		if (inst_a == IS_OK) {
		    ++no_read; (*ins_a)->sid = ++idx_a;
		} else if (inst_a == IS_END) break;
		ldr2->reset();	// rewind
		idx_b = ldr2->var_no;
		inst_b = ldr2->nextvar(ins_b);
		if (inst_b == IS_OK) {
		    ++no_read;
		    (*ins_b)->sid = idx_b;
		}
	    }
	    break;
	 }
	 if (inst_a == IS_END || inst_b == IS_END) return (0);
	} while (inst_a == IS_ERR || inst_b == IS_ERR);
	return (no_read);
}

#endif	// _CALCSERV_
