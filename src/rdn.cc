/*****************************************************************************
*
*	Pick up some members from input multiple sequence alignment
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

#include "aln.h"
#if USE_WEIGHT
#include "mseq.h"
#endif
#include "autocomp.h"

struct	CARRY {int id, num; FTYPE frq;};
struct	STAT3 {long val, amc, tmc, mmc, gap, unp;};
struct	HOMSET {int* sib; int* last;};
struct	Submembers
{
	int	num;
	int*	which;
	bool	check(int manys);
	Submembers(Seq* sd, int argc, const char** argv, FILE* fd);
	Submembers(int many, FILE* fi);
	~Submembers() {delete[] which;}
};

static	int	pcmpf(const CARRY* a, const CARRY* b);
static	int	ccmpf(const CARRY* a, const CARRY* b);
static	void	compress(CHAR* dseq, int dmany, CHAR* sseq, FTYPE* wt, int smany, int max_code);
static	void	prntout(Seq* dst);
static	Seq*	choseseq(Seq* dst, Seq* src);
static	int	icmpf(const int* a, const int* b);
static	Seq*	eliminseq(Seq* dst, Seq* src);
static	void	breakup(Seq* sd, Seq* tmp, int delgap);
static	Seq*	rdcseq(Seq* dst, Seq* src, int many);
static	Seq*	justseq(Seq* dst, int ichi);
#if FVAL
static	Seq*	uniqseq(Seq* dst, Seq* src);
#endif

static	Submembers*	submem = 0;
static	FILE*	fromfile = 0;

#if REPORT_STAT
static	void	statseq(STAT3* stat, Seq* a);
static	void	eval(Seq* sd, STAT3* stat, CHAR* x, CHAR* y, int nos);
static	STAT3	matches;
#endif

void usage()
{
	fputs("rdn -cS MSA [list_of_member_nos]\n", stderr);
	fputs("\tex: rdn -csd MSA_10 2 4 6\n", stderr);
	fputs("\t-cS S (- [bcdejlprsu] \n", stderr);
	fputs("\t  b: decomposition\n", stderr);
	fputs("\t  c: centering\n", stderr);
	fputs("\t  d: delete common gaps\n", stderr);
	fputs("\t  e: eliminate spcified member\n", stderr);
	fputs("\t  j: justify at both ends\n", stderr);
	fputs("\t  l: left justify\n", stderr);
	fputs("\t  p: representatives\n", stderr);
	fputs("\t  r: right justify\n", stderr);
	fputs("\t  s: seletct spcified member\n", stderr);
#if FVAL
	fputs("\t  u: make unique\n", stderr);
#endif
}

template <class seq_t>
void AlnServer<seq_t>::setparam(int lvl)
{
	setalprm();
	setlpw(QUERY);
	setform(QUERY);
	setprmode(QUERY, QUERY, QUERY);
}

static void prntout(Seq* dst)
{
	if (!dst) return;

#if REPORT_STAT
	int 	nos = dst->many;
	dst->fphseq(out_fd);
	nos = nos * (nos - 1) / 2;
	fparam(out_fd);
	fprintf(out_fd, "Dist = %3ld (%3ld), t-id. = %3ld (%3ld), ",
	    matches.val, matches.val/nos, matches.tmc, matches.tmc/nos);
	fprintf(out_fd, "t-mm = %3ld (%3ld)\n", matches.mmc, matches.mmc/nos);
	fprintf(out_fd, "a-id. = %3ld, Gaps = %3ld, Unpairs = %3ld",
	    matches.amc, matches.gap, matches.unp);
#endif
	dst->typeseq();
}

static int pcmpf(const CARRY* a, const CARRY* b)
{
	float	df = b->frq - a->frq;
	if (df == 0) return (0);
	return (df > 0? 1: -1);
}

static int ccmpf(const CARRY* a, const CARRY* b)
{
	return (b->num - a->num);
}

static 	CARRY*	count;

static void compress(CHAR* dseq, int dmany, CHAR* sseq, FTYPE* wt, int smany, int max_code)
{
	CARRY*	pc = count;

	for (int i = 0; i < max_code; i++, pc++) {
	    pc->id = i;
	    pc->frq = 0.;
	}
	for (int i = 0; i < smany; i++) 
	    count[*sseq++].frq += *wt++;
	qsort(count, max_code, sizeof(CARRY), (CMPF) pcmpf);
	int	n = 0;
	pc = count;
	for (int i = 0; i++ < max_code; pc++) {
	    pc->num = (int) pc->frq;
	    pc->frq -= pc->num;
	    n += pc->num;
	}
	qsort(count, max_code, sizeof(CARRY), (CMPF) pcmpf);
	for (int i = 0; n++ < dmany; ) 
	    ++count[i++].num;
	qsort(count, max_code, sizeof(CARRY), (CMPF) ccmpf);
	pc = count;
	for (n = 0; n < dmany; pc++) {
	    n += pc->num;
	    while (pc->num-- > 0)
		*dseq++ = pc->id;
	}
}

bool Submembers::check(int manys)
{
	int*	wch = which;
	for ( ; *wch > 0; ++wch) 
	    if (0 > --(*wch) || *wch >= manys) return (false);
	*wch = -1;
	return (true);
}

Submembers::Submembers(Seq* sd, int argc, const char** argv, FILE* fd)
	: num(0), which(0)
{
	StrHash<int>	strsubmems(sd->many);
	Dhash<int, int>	numsubmems(sd->many, 0);
	for (int c = 0; c < argc; ++c) {
	    if (isnumber(argv[c])) {
		int	from = atoi(argv[c]);
		const char*	hyphen = strchr(argv[c], '-');
		if (hyphen) {
		    if (isnumber(++hyphen)) {
			int	to = atoi(hyphen);
			if (from > to) swap(from, to);
			if (to > sd->many) to = sd->many;
			while (from <= to) numsubmems.incr(from++);
		    } else	numsubmems.incr(from);
		} else	numsubmems.incr(from);
	    } else strsubmems.incr(argv[c]);
	}
	if (fd) {
	    char	term[MAXL];
	    char*	ps = term;
	    int 	c = 0;
	    while ((c = fgetc(fd)) != EOF) {
		if (isspace(c)) {
		    if (ps > term) {
			*ps = '\0';
			strsubmems.incr(term);
			ps = term;
		    }
		} else	*ps++ = c;
	    }
	}
	int	nn = numsubmems.count();
	int	ns = strsubmems.count();
	if ((nn + ns) == 0) return;
	which = new int [nn + ns + 1];
	int*	ml = which;
	if (nn)
	    for (int m = 0; m < sd->many; ++m) 
		if (numsubmems.find(m + 1)) *ml++ = m;
	if (ns)
	    for (int m = 0; m < sd->many; ++m) 
		if (strsubmems.find((*sd->sname)[m])) *ml++ = m;
	num = ml - which;
	*ml = -1;
	if (num) return;
	delete[] which;
	num = 0; which = 0;
}
 
Submembers::Submembers(int many, FILE* fi)
{
	which = new int[many + MAXL];

	do {
	    num = fgetiarray(which, many, fi);
	    if (!num) {
		delete[] which;
		return;
	    }
	} while (!check(many));
}

static Seq* choseseq(Seq* dst, Seq* src)
{
	if (!submem) submem = new Submembers(src->many, stdin);
	if (!submem->num) return(0);
	dst = src->extseq(dst, submem->which, CPY_ALL + RLT_SBI);
	return (dst);
}

static int icmpf(const int* a, const int* b)
{
	return (*a - *b);
}

static Seq* eliminseq(Seq* dst, Seq* src)
{
	if (!submem) submem = new Submembers(src->many, stdin);
	if (!submem->num || submem->num == src->many) return(0);
	int*	delet = submem->which;
	qsort((UPTR) delet, submem->num, sizeof(int), (CMPF) icmpf);
	int*	which = new int[src->many + 2];
	int	k = 0;
	for (int i = 0, j = 0; i < src->many; ++i) {
	    if (i == delet[j])	++j;
	    else	which[k++] = i;
	}
	which[k] = -1;
	dst = src->extseq(dst, which, CPY_ALL + RLT_SBI);
	delete[] which;
	return (dst);
}

#if FVAL

static Seq* uniqseq(Seq* dst, Seq* src)
{
	FTYPE*  dist = calcdist(src);
	int*	which = new int[src->many + 2];
	int*	parent = new int[src->many];
	int*	linkm = new int[src->many];
	int*	lastm = new int[src->many];
	double*	dstsum = new double[src->many];
	double	uniq_thr = alprm.thr;

	for (int i = 0; i < src->many; ++i) {
	    parent[i] = linkm[i] = lastm[i] = -1;
	    which[i] = i;
	    dstsum[i] = 0.;
	}
	int k = 0;
	for (int j = 1; j < src->many; ++j) {
	    for (int i = 0; i < j; ++i, ++k) {
		if (dist[k] <= uniq_thr) {
		    int	p = parent[i];
		    if (p == -1) p = i;
		    parent[j] = p;
		    linkm[j] = i;
		    lastm[p] = j;
		    dstsum[i] += dist[k];
		    dstsum[j] += dist[k];
		}
	    }
	}
	for (int p = 0; p < src->many; ++p) {
	    if (lastm[p] >= 0) {
		int i = lastm[p];
		for (int j = i; j >= 0; j = linkm[j])
		    if (dstsum[j] < dstsum[i]) i = j;
		for (int j = lastm[p]; j >= 0; j = linkm[j])
		    if (i != j) which[j] = -1;
	    }
	}
	for (int i = k = 0; i < src->many; ++i)
	    if (which[i] >= 0) which[k++] = which[i];
	which[k] = -1;
	dst = src->extseq(dst, which, CPY_ALL + RLT_SBI);
	delete[] dist; delete[] parent; delete[] linkm; delete[] lastm;
	delete[] which; delete[] dstsum;
	return (dst);
}

#endif	// FVAL

static void breakup(Seq* sd, Seq* tmp, int delgap)
{
	int	buf[2] = {0, -1};
	char	str[MAXL];
	bool	ul = sd->sname && (*sd->sname)[0];

	topath(str, OutPrm.out_file);
	if (!ul) strcat(str, sd->sqname());
	char*	ps = str + strlen(str);
	int	slen = MAXL - strlen(str);

	for (int i = 0; i < sd->many; ++i) {
	    if (ul)	strncpy(ps, (*sd->sname)[i], slen);
	    else	snprintf(ps, slen, ".%d", i+1);
	    FILE*	fd = fopen(str, "w");
	    if (!fd) fatal(no_file);
	    buf[0] = i;
	    tmp = sd->extseq(tmp, buf, CPY_ALL);
	    if (delgap) tmp->elim_column(DEL_GAP);
	    tmp->typeseq(fd);
	    fclose(fd);
	}
}

static Seq* rdcseq(Seq* dst, Seq* src, int many)
{
	char	str[MAXL];
	FTYPE	fac = (FTYPE) many / (FTYPE) src->many;
	int 	max_code = src->code->max_code;
	FTYPE*	wt = new FTYPE[src->many];

	if (!count) count = new CARRY[max_code];
#if USE_WEIGHT
	if (src->weight) 
	    for (int i = 0; i < src->many; i++)
		wt[i] = src->weight[i] * many;
	else
#endif
	    for (int i = 0; i < src->many; i++)
		wt[i] = fac;
	dst->refresh(many, (src->right - src->left + 2));
	CHAR*	sseq = src->at(src->left);
	CHAR*	dseq = dst->at(0);
	for (int i = src->left; i < src->right; 
	    i++, sseq += src->many, dseq += many)
		compress(dseq, many, sseq, wt, src->many, max_code);
	for (int i = 0; i < many; i++) {
	    snprintf(str, MAXL, "RDC%d", i);
	    dst->sname->push(str);
	}
	dst->postseq(dseq);
	delete[] wt;
	return (dst);
}

static Seq* justseq(Seq* dst, int ichi)
{
	Seq*	src = dst->copyseq(0, CPY_ALL);
	int 	m = 0;
	CHAR*	sb = src->at(src->left);
	CHAR*	st = src->at(src->right);
	CHAR*	db = dst->at(dst->left);
	CHAR*	dt = dst->at(dst->right);
	if (ichi == 1) OutPrm.trimendgap = 1;
	for (int k = 0; k < src->many; k++) {
	    CHAR*	ds = db;
	    int	n = 0, l = 0;
	    if (ichi > 1) {
		CHAR*	ss = sb;
		for ( ; ss < st; ss += src->many)
		    if (IsGap(*ss)) n++;
		    else	l++;
		if (ichi == 2 || ichi == 4) m = n / 2;
		if (ichi < 4) {
		    while (n-- > m) {
			*ds = gap_code;
			ds += dst->many;
		    }
		}
	    }
	    if (ichi == 4) st = src->at(l / 2);
	    CHAR* ss = sb;
	    for ( ; ss < st; ss += src->many) {
		if (IsntGap(*ss)) {
		    *ds = *ss;
		    ds += dst->many;
		}
	    }
	    if (ichi == 4) {
		while (n--) {
		    *ds = gap_code;
		    ds += dst->many;
		}
		st = src->at(src->right);
		for ( ; ss < st; ss += src->many) {
		    if (IsntGap(*ss)) {
			*ds = *ss;
			ds += dst->many;
		    }
		}
	    }
	    while (ds < dt) {
		*ds = gap_code;
		ds += dst->many;
	    }
	    sb++; st++; db++; dt++;
	}
	delete (src);
	return (dst);
}

#if REPORT_STAT
static void eval(Seq* sd, STAT3* stat, CHAR* x, CHAR* y, int nos)
{
	int	k= 0;

	stat->val -= (*sim1)(x);
	for (int j = 0; j < nos; j++) {
	    if (x[j] == gap_code) {
		stat->unp++;
		if (IsntGap(y[j])) stat->gap++;
	    }
	    for (int i = 0; i < j; i++) { 
		if (x[j] == gap_code) continue;
		if (x[i] == x[j]) {
		    stat->mch++ ;
		    k++;
		} else if (IsntGap(x[i]))
		    stat->mmc++;
	    }
	}
	if (k == nos * (nos-1) / 2) stat->mch++;
}

static void statseq(STAT3* stat, Seq* a)
{
	CHAR*	as = a->at(a->left);
	CHAR*	ps = as - a->many;

	setsim1(a);
	clearmem(stat, sizeof(STAT3));
	for (int c = a->left; c < a->right; c++, ps = as, as += a->many) {
	    eval(a, stat, as, ps, a->many);
	}
}
#endif

static	char	ans[MAXL] = "";
static	char*	cmdl = 0;

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
	const	char*	opt = argv[0] + 1;
	const	char*	val = argv[0] + 2;

	switch (*opt) {
	  case 'c': cmdl = strcpy(ans, *argv + 2); break;
	  case 'f': 
	    if ((val = getarg(argc, argv)))
		fromfile = fopen(val, "r");
	     break;
	  default:	return (0);
	}
	return (1);
}

int rdn_main(CalcServer<Seq>* svr, Seq* seqs[], ThQueue<Seq>* q)
{
	Seq*&	src = seqs[0];
	Seq*	dstsq = 0;
	int	manys = 0;
	int	ichi = 0;
	double	gfrac = 0.;
static	const char mssg[] = 
	    "p[rop]/s[el]/l[eft]/c[ent]/r[ight]/B[reak]/d[el] : ";

	if (!cmdl) progets(ans, mssg);

	for (char* wk = ans; *wk; ++wk) {
	    switch (*wk) {
		case 'b':
		case 'B':
		    breakup(src, seqs[1], strchr(ans, 'd') || *wk == 'B');
		    return (OK);
		case 'c': ichi = 2; break;	/* centering	*/
		case 'd': gfrac = atof(wk+1);	/* delete gap	*/
		    if (gfrac > 0.) {
			if (manys == 0) manys = src->many;
			dstsq = src->copyseq(seqs[1], CPY_ALL);
			dstsq->elim_column(DEL_GRC, gfrac);
			while(*++wk && isspace(*wk)) ;
			while(*++wk && !isspace(*wk)) ;
		    }
		    break;
		case 'e':
		    manys = src->many;
		    dstsq = eliminseq(seqs[1], src);
		    if (!dstsq) return (!OK);
		    break;
		case 'j': ichi = 4; break;	/* justify at both ends */
		case 'l': ichi = 1; break;	/* left justify	*/
		case 'p': 
		    promptin("How many ? (%d) : ", &manys);
		    if (manys > 0) {
			dstsq = rdcseq(seqs[1], src, manys);
			if (!dstsq) return (!OK);
		    }
		    break;
		case 'r': ichi = 3; break;	/* right justify */
		case 's':
		    manys = src->many;
		    dstsq = choseseq(seqs[1], src);
		    if (!dstsq) return (!OK);
		    break;
#if FVAL
		case 'u':
		    manys = src->many;
		    dstsq = uniqseq(seqs[1], src);
		    if (!dstsq) return (!OK);
		    break;
#endif
		default: break;
	    }
	}
	if (manys == 0) dstsq = src->copyseq(seqs[1], CPY_ALL);
	if (ichi) dstsq = justseq(dstsq, ichi);
	if (ichi || strchr(ans, 'd'))  dstsq->elim_column(DEL_GAP);
	snprintf(ans, MAXL, "%s%d", src->sqname(), dstsq->many);
	strrealloc(dstsq->spath, ans);
	prntout(dstsq);
	*ans = '\0';
	return (OK);
}

int main(int argc, const char** argv)
{
	alprm.thr = 0;
	alprm2.spb = 1;
	setprmode(Row_Last, 'L', SILENT);
	optimize(GLOBAL, MINIMUM);
	AlnServer<Seq>	svr(argc, argv, IM_MULT, IM_SNGL, (void*) 0, &rdn_main);
	if (argc) {
	    svr.sql1->nextseq(svr.in_face);
	    --argc; ++argv;
	    if (argc || fromfile)
		submem = new Submembers(svr.in_face[0], argc, argv, fromfile);
	}
	svr.menucomp();
	delete submem;
	if (fromfile) fclose(fromfile);
	return (0);
}
