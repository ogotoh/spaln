/*****************************************************************************
*
*	Estimation of nucleotide divergence between aligned sequences
*
*	Ref. 1:	  Jukes,T.H. and Cantor,C.R. (1969)
*	    "Evolution of protein molecules." in Mammalian Protein Metabolism", 
*	    (Munro,H.N., ed.) pp. 21-132. Academic Press, New York.
*	Ref. 2:   Kimura,M. (1981)
*	    "Estimation of evolutionary distances between homologous
*	    nucleotide sequences." Proc. Natl. Acad. Sci. USA 78, 5773-5777.
*	Ref. 3:   Tajima,F. and Nei,M. (1984)
*	    "Estimation of evolutionary distance between nucleotide sequences."
*	    Mol. Biol. Evol. 1, 269-286.
*
*	Output format:
*	uncorr'd_div.  Jukes-Cantor  Kimura-2  Tajima-Nei    Seq(i)  Seq(j)  
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

#include <math.h>
#include "aln.h"
#include "autocomp.h"
#include "divseq.h"

#if M_THREAD
#include <pthread.h>
#include <unistd.h>
#include <sched.h>
#endif

enum method {jc, p2, tn, end};
enum base {Ad, Ci, Gu, Th};

static	const	int	NELEM = 5;
static	const	int	GAP = NELEM - 1;
static	const	INT	BasePerShort = 4 * sizeof(SHORT);

static	bool	putmat = false;
static	const	char*	groups = 0;

static	int	skip_terms = 2;
static	float*	distance = 0;

class Dvn
{
	double	frq[NELEM][NELEM];
	int	base;
public:
	void	printmat(FILE* fd);
	void	countup(int a, int b);
	void	subfrq(Seq* sd, int* group1, int* group2);
	double	fracnid();
	double	kimura2();
	double	tajnei();
	Dvn(Seq* sd, int* group1, int* group2);
	~Dvn() {}
};

class Bitsqtab
{
	int	width;
	int	bitln;
	SHORT** bitsq;
	CHAR*	lutab;
public:	
	Bitsqtab(Seq* sd);
	~Bitsqtab()	{
	    if (bitsq) {delete[] *bitsq; delete[] bitsq;}
	    delete[] lutab;
	}
	void	lookuptab();
	double	frac_mismatch(int i, int j);
};

struct Member
{
	int	sid;
	int	pid;
	bool	vrtl;
	Member(int id = 0) : sid(id), vrtl(false) {}
	void	setids(int id, int pd) {sid = id; pid = pd;}
	int	fget(FILE* fd, const char* fn = 0) {return (0);}
};

struct Species : public Member
{
	char*	sname;
	int*	paralogs;
	int	many;
	Species(int id, char* str, Seq*& msd, int sz);
	Species() : sname(0), paralogs(0), many(0) {}
	~Species()	{delete[] sname; delete[] paralogs;}
};

void usage()
{
	fputs("Usage:\tdvn [options] msa\n", stderr);
	fputs(" or\tdvn [options] -ix:catalog msa\n", stderr);
	fputs(" or\tdvn [options] [-ix] -g spec_list msa\n", stderr);
 	fputs("\nOptions:\n", stderr);
	fputs("\t\t-An\tn=0: simple; n=1: quick;\n", stderr);
	fputs("\t\t-On\tn=0: %%div; n=1: %%UC;\n", stderr);
	fputs("\t\t-k skip_terms (def = 2)\n", stderr);
	fputs("\t\t-m \tprint matrix\n", stderr);
	exit(1);
}

template <class seq_t>
void AlnServer<seq_t>::setparam(int lvl)
{
	switch (progetc("Print freq. matrix ? [y/n] : ")) {
	    case 'y':
	    case 'Y': putmat = true; break;
	    case 'n':
	    case 'N': putmat = false; break;
	    default : break;
	}
}

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
	int	n = 1;
	switch (argv[0][1]) {
	  case 'g': groups = getarg(argc, argv); break;
	  case 'h': usage();
	  case 'k': skip_terms = atoi(getarg(argc, argv, true)); break;
	  case 'm': putmat = true; break;
	  default: n = 0; break;
	}
	return (n);
}

Species::Species(int id, char* str, Seq*& msd, int sz) : Member(id)
{
	char*	qs = str;
	char*	ps = car(qs);
	sname = new char[strlen(ps) + 1];
	strcpy(sname, ps);
	StrHash<int> sh(sz);
	many = 0;
	for (++qs; *(ps = car(qs)); ++qs) {
	    sh.incr(ps);
	    ++many;
	}
	paralogs = new int[many];
	int*	pi = paralogs;
	for (int i = 0; i < msd->many; ++i) {
	    if (sh.find((*msd->sname)[i]))
		*pi++ = i;
	}
}

void Dvn::printmat(FILE* fd)
{
	putc('\n', fd);
	for (int i = 0; i < NELEM; i++) {
	    for (int j = 0; j < NELEM; j++) 
		fprintf(fd, "%6.1lf ", frq[i][j]);
	    putc('\n', fd);
	}
}

void Dvn::countup(int a, int b)
{
static	const	int	frac[] = 
		{1, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
	int	aa = a - base;
	int	bb = b - base;
	double	f = 1./((double) frac[aa]*frac[bb]);

	if (IsGap(a)) {
	    for (int j = 0; j < GAP; ++j, bb >>= 1)
		if (bb & 1) frq[GAP][j] += f;
	} else if (IsGap(b)) {
	    for (int i = 0; i < GAP; ++i, aa >>= 1)
		if (aa & 1) frq[i][GAP] += f;
	} else {
	    for (int i = 0; i < GAP; ++i, aa >>= 1)
		if (aa & 1)
		    for (int j = 0, cc = bb; j < GAP; ++j, cc >>= 1)
			if (cc & 1) frq[i][j] += f;
	}
}

Dvn::Dvn(Seq* sd, int* group1, int* group2)
{
	base = sd->code->base_code - 1;
	for (int i = 0; i < NELEM; ++i)
	    for (int j = 0; j < NELEM; ++j) frq[i][j] = 0.;
	CHAR*	cs = sd->at(sd->left);
	for (int i = sd->left; i < sd->right; ++i, cs += sd->many) {
	    for (int* g1 = group1; *g1 >= 0; ++g1)
		for (int* g2 = group2; *g2 >= 0; ++g2)
		    countup(cs[*g1], cs[*g2]);
	}
}

double Dvn::fracnid()
{
	double	id = 0.;
	double	nd = 0.;

	for (int i = 0; i < GAP; ++i)
	    for (int j = 0; j < GAP; ++j) {
		if (i == j)	id += frq[i][j];
		else		nd += frq[i][j];
	    }
	return (nd / (id + nd));
}

double Dvn::kimura2()
{
	double	r = frq[Ad][Ad] + frq[Ci][Ci] + frq[Gu][Gu] + frq[Th][Th];
	double	p = frq[Ad][Gu] + frq[Gu][Ad] + frq[Ci][Th] + frq[Th][Ci];
	double	q = frq[Ad][Ci] + frq[Ad][Th] + frq[Gu][Ci] + frq[Gu][Th]
	  + frq[Ci][Ad] + frq[Ci][Gu] + frq[Th][Ad] + frq[Th][Gu];
	r += p + q;
	if (r != 0.) {
	    p /= r;
	    q /= r;
	}
	return (-.5 * log((1. - 2. * p - q) * sqrt(1. - 2. * q)));
}

double Dvn::tajnei()
{
	double	g[NELEM];
	double  n = 0;
	double	p = fracnid();
	for (int i = 0; i < GAP; ++i) {
	    g[i] = 0.;
	    for (int j = 0; j < GAP; j++)   // g[i]: 2 * Num of i-th base
		g[i] += frq[j][i] + frq[i][j];
	    if (g[i] == 0.) return (999.999);
	    n += g[i];			// n: 2 * Num of sites
	}
	for (int i = 0; i < GAP; i++) g[i] /= n;
					// Now g[i]: average base comp.
	n /= 2;				// Now n: Num of sites
	double	c = 0;
	for (int i = 1; i < GAP; ++i) {
	    for (int j = 0; j < i; ++j) {
		double	b = frq[i][j] + frq[j][i];
		c += b * b / g[i] / g[j];
	    }
	}
	if (c == 0.) return (0.);
	c /= 2 * n * n;			// c has been normalized
	double	b = 0;
	for (int i = 0; i < GAP; ++i) b += g[i]*g[i];
	b = (1. - b + p * p / c) / 2.;
	return (-b * log(1. - p / b));
}

/***********************************************************
*	Convert a nucleotide MSA into a 2-bit matrix
*	no amb nor gap is allowed
***********************************************************/

void Bitsqtab::lookuptab()
{
static	const	INT	tsize = USHRT_MAX + 1;
	CHAR*	w = lutab = new CHAR[tsize];

	for (INT i = 0; i < tsize; ++i) {
	    CHAR	n = 0;
	    for (INT j = i; j; j >>= 2) {
		if (j & 3) ++n;
	    }
	    *w++ = n;
	}
}

Bitsqtab::Bitsqtab(Seq* sd)
	: width(sd->right - sd->left), bitln((width + BasePerShort - 1) / BasePerShort)
{
	bitsq = new SHORT*[sd->many];
	SHORT*	bitwk = *bitsq = new SHORT[sd->many * bitln];
	CHAR*	s = sd->at(sd->left);

	for (int i = 0; ; ++s) {
	    SHORT	bit = 0;
	    CHAR*	ws = s;
	    for (int j = 0; j++ < width; ws += sd->many) {
		CHAR	rnc = ncredctab[*ws];
		if (rnc < 4) bit += rnc;
		else bit += rand() % 4;
		if (j % BasePerShort == 0) {
		    *bitwk++ = bit;
		    bit = 0;
		} else {
		    bit <<= 2;
		}
	    }
	    if (bit) *bitwk = bit;
	    if (++i == sd->many) break;
	    bitwk = bitsq[i] = bitsq[i-1] + bitln;
	}
	lookuptab();
}

double Bitsqtab::frac_mismatch(int i, int j)
{
	SHORT*	si = bitsq[i];
	SHORT*	sj = bitsq[j];
	int	mmc = 0;
	for (int k = 0; k < bitln; ++k) {
	    mmc += lutab[*si ^ *sj];
	    ++si;
	    ++sj;
	}
	return ((double) mmc / (double) width);
}

/**************************************************************
 * 	main part
**************************************************************/

static	void	stndrddvn(Seq* sd);
static	void	lookupdvn(Seq* sd);
static	void	simpledvn(Seq* sd);
static	double	simple_pair(const Seq*  sd, const int& i, const int& j);
static	Bitsqtab*	bst = 0;
static	double	DistThr = alprm.thr / 100.;

static	const	char	frmt0[] = "%8.4lf %-12s %-12s\n";
static	const	char	frmt1[] = "%8.4lf %8.4lf %-12s %-12s\n";

static void print_dvn(double nid, const char* as, const char* bs)
{
	if (algmode.nsa == 0)
	    fprintf(out_fd, frmt0, 100. * nid, as, bs);
	else
	    fprintf(out_fd, frmt1, 100. * nid, 100. * jukcan(nid), as, bs);
}

static void stndrddvn(Seq* sd)
{
	for (int j = 1; j < sd->many; ++j) {
	    for (int i = 0; i < j; ++i) {
		double nid = pairdvn(sd, i, j);
		if (nid <= DistThr) 
		    print_dvn(nid, (*sd->sname)[i], (*sd->sname)[j]);
	    }
	}
}

static double simple_pair(const Seq* sd, const int& i, const int& j)
{
const	CHAR*	as = sd->at(sd->left);
const	CHAR*	ts = sd->at(sd->right);
const	CHAR*	bs = as + j;
	int	nid = 0;
	for (as += i; as < ts; as += sd->many, bs += sd->many)
	    if (*as != *bs) ++nid;
	return ((double) nid / ((double) (sd->right - sd->left)));
}

static void simpledvn(Seq* sd)
{
	for (int j = 1; j < sd->many; ++j) {
	    for (int i = 0; i < j; ++i) {
		double nid = simple_pair(sd, i, j);
		if (nid <= DistThr) 
		    print_dvn(nid, (*sd->sname)[i], (*sd->sname)[j]);
	    }
	}
}

static void lookupdvn(Seq* sd)
{
	for (int j = 1; j < sd->many; ++j) {
	    for (int i = 0; i < j; ++i) {
		double	nid = bst->frac_mismatch(i, j);
		if (nid <= DistThr)
		    print_dvn(nid, (*sd->sname)[i], (*sd->sname)[j]);
	    }
	}
}

static int dvn_main(CalcServer<Member>* svr, Member** vars, ThQueue<Member>* q)
{
	Seq*	sd = (Seq*) svr->prm;
	Member*&	a = vars[0];
	Member*&	b = vars[1];
	int&	i = a->pid;
	int&	j = b->pid;
	double nid = bst? bst->frac_mismatch(i, j): pairdvn(sd, i, j);
	if (distance)
	    distance[svr->calcnbr(a->sid, b->sid)] = nid;
	else if (nid <= DistThr) {
	    m_thread_Lock(q);
	    print_dvn(nid, (*sd->sname)[i], (*sd->sname)[j]);
	    m_thread_Unlock(q);
	}
	return (OK);
}

static int dvn_output(CalcServer<Member>* svr, Member** vars, ThQueue<Member>* q)
{
	Seq*	sd = (Seq*) svr->prm;
	Member*&	a = vars[0];
	Member*&	b = vars[1];
	double	nid = distance[svr->calcnbr(a->sid, b->sid)];
	if (nid <= DistThr)
	    print_dvn(nid, (*sd->sname)[a->pid], (*sd->sname)[b->pid]);
	return (OK);
}

static void calcdvn(AlnServer<Seq>& svr)
{
	Seq*&   msd = svr.in_face[0];
	size_t	ng1 = 0;
	size_t	ng2 = 0;
	size_t	many = 0;
	Member*	memb = 0;
	if (svr.catalog) {
	    StrHash<int> sh(msd->many);
	    FILE* fc = fopen(svr.catalog, "r");
	    if (!fc) fatal("%s not found !\n", svr.catalog);
	    char	str[MAXL];
	    while (fgets(str, MAXL, fc)) {
		if (isBlankLine(str)) {
		    if (!ng1) ng1 = ng2;
		    continue;
		}
		char*	qs = str;
		char*	ps = car(qs);
		sh.incr(ps);
		++ng2;
	    }
	    fclose(fc);
	    many = ng2;
	    ng2 -= ng1;
	    if (!ng1) swap(ng1, ng2);
	    Member*	mwk = memb = new Member[many + 1];
	    for (int i = many = 0; i < msd->many; ++i) {
		if (sh.find((*msd->sname)[i]))
		    (mwk++)->setids(many++, i);
	    }
	} else {
	    Member*     mwk = memb = new Member[msd->many + 1];
	    for (int i = 0; i < msd->many; ++i) (mwk++)->setids(many++, i);
	}
	Member**	mmbr = new Member*[many + 1];
	for (INT i = 0; i <= many; ++i) mmbr[i] = &memb[i];
	size_t	nn = svr.calc_mode == IM_GRUP? ng1 * ng2: svr.calcsize(many);
	if (!nn) fatal("set '-ig:catalog' option!\n");
	distance = (thread_num > 1)? new float[nn]: 0;
	if (svr.calc_mode != IM_GRUP) ng1 = many;
	CalcServer<Member> csv(svr.calc_mode, (void*) msd, dvn_main, 0, 0, mmbr, ng1, ng2);
	csv.autocomp();
	if (distance) {
	    csv.change_job(dvn_output);
	    csv.autocomp(false);
	    delete[] distance;
	}
	delete[] memb; delete[] mmbr;
}

static int spcdvn_main(CalcServer<Species>* svr, Species** spcs, ThQueue<Species>* q)
{
	Seq*	msd = (Seq*) svr->prm;
	Species*&	a = spcs[0];
	Species*&	b = spcs[1];
	double	nid = algmode.alg > 1? 0: DBL_MAX;
	for (int i = 0; i < a->many; ++i) {
	    for (int j = 0; j < b->many; ++j) {
		double	d = pairdvn(msd, a->paralogs[i], b->paralogs[j]);
		if (algmode.alg) nid += d;
		else if (d < nid) nid = d;
	    }
	}
	if (algmode.alg > 1) nid /= (a->many * b->many);
	if (distance) 
	    distance[svr->calcnbr(a->sid, b->sid)] = nid;
	else if (nid <= DistThr) {
	    m_thread_Lock(q);
	    print_dvn(nid, a->sname, b->sname);
	    m_thread_Unlock(q);
	}
	return (OK);
}

static int spcdvn_output(CalcServer<Species>* svr, Species** spcs, ThQueue<Species>* q)
{
	Species*&	a = spcs[0];
	Species*&	b = spcs[1];
	double	nid = distance[svr->calcnbr(a->sid, b->sid)];
	if (nid <= DistThr)
	    print_dvn(nid, a->sname, b->sname);
	return (OK);
}

static void comp_group(AlnServer<Seq>& svr)
{
	Seq*&   msd = svr.in_face[0];
	StrHash<int>*	sh = 0;
	char	str[MAXL];
	size_t	ng1 = 0;
	size_t	ng2 = 0;
	if (svr.catalog) {
	    sh = new StrHash<int>(msd->many);
	    FILE* fc = fopen(svr.catalog, "r");
	    if (!fc) fatal("%s not found !\n", svr.catalog);
	    while (fgets(str, MAXL, fc)) {
		if (isBlankLine(str)) {
		    if (!ng1) ng1 = ng2;
		    continue;
		}
		char*	qs = str;
		char*	ps = car(qs);
		sh->incr(ps);
		++ng2;
	    }
	    fclose(fc);
	    ng2 -= ng1;
	    if (!ng1) swap(ng1, ng2);
	}
	FILE*	fg = fopen(groups, "r");
	if (!fg) fatal("%s not found !\n", groups);
	size_t	ng = 0;
	size_t	np = 0;
	size_t	slen = 0;
	size_t	maxl = 0;
	size_t	maxn = 0;
	while (fgets(str, MAXL, fg)) {
	    char*	ps = str;
	    if (!np) {		// top of line
		for (int i = 0; i < skip_terms; ++i)
		    if (!*(ps = cdr(ps))) break;
		if (sh) {
		    char*	qs = ps;
		    ps = car(qs);
		    if (sh->find(ps)) {
			++np;
			ps = ++qs;
		    } else {
			while (str[strlen(str) - 1] != '\n')
			    if (!fgets(str, MAXL, fg)) break;
			continue;
		    }
		}
	    }
	    while (*(ps = cdr(ps))) ++np;
	    slen += ps - str;
	    if (ps[-1] == '\n') {
		++ng;
		if (slen > maxl) maxl = slen;
		if (np > maxn) maxn = np;
		slen = np = 0;
	    }
	}
	rewind(fg);
	Species**	spec = new Species*[ng];
	Species**	spwk = spec;
	char*	buf = new char[++maxl];
	int	sid = 0;
	while (fgets(buf, maxl, fg)) {
	    char*	ps = buf;
	    for (int i = 0; i < skip_terms; ++i)
		if (!*(ps = cdr(ps))) break;
	    if (sh) {
		char*	qs = ps;
		ps = cdr(qs);
		if (sh->find(ps)) *qs = ' ';
		else continue;
	    }
	    *spwk++ = new Species(sid++, ps, msd, maxn);
	}
	delete sh;
	fclose(fg);
	size_t	nn = svr.calc_mode == IM_GRUP? ng1 * ng2: svr.calcsize(ng);
	if (!nn) fatal("set '-ig:catalog' option!\n");
	distance = (thread_num > 1)? new float[nn]: 0;
	if (svr.calc_mode != IM_GRUP) ng1 = ng;
	CalcServer<Species> csv(svr.calc_mode, (void*) msd, spcdvn_main, 0, 0, spec, ng1, ng2);
	csv.autocomp();
	if (distance) {
	    csv.change_job(spcdvn_output);
	    csv.autocomp(false);
	    delete[] distance;
	}
	for (INT i = 0; i < ng; ++i) delete spec[i];
	delete[] spec;
}

int main(int argc, const char** argv)
{
const	char*	tmpcat = 0;
	alprm.thr = 100;	// no threshold
	AlnServer<Seq>	svr(argc, argv, IM_MULT, IM_EVRY);
	swap(svr.catalog, tmpcat);
	Seq*&   msd = svr.in_face[0];
	if (svr.sql1->nextseq(&msd) != IS_OK) usage();
	swap(svr.catalog, tmpcat);
	msd->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	if  (algmode.alg == 1) bst = new Bitsqtab(msd);
	DistThr = alprm.thr / 100.;
	if (groups) comp_group(svr);
	else if (svr.catalog || thread_num > 1) calcdvn(svr);
	else {
          switch (algmode.alg) {
            case 1:  lookupdvn(msd); break;
	    case 2:  simpledvn(msd); break;
            default: stndrddvn(msd); break;
          }
	}
	EraDbsDt();
	delete bst;
	return (0);
}
