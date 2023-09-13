/*****************************************************************************
*
*	fitild.c; Fit intron length distribution to a composite 
*	Frechet distributions with BFGS algorithm
*
*	Requires GSL library 
*	https://www.gnu.org/software/gsl/
*
*	Requires libLBFGS library 
*	http://www.chokkan.org/software/liblbfgs/
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

#include "calcserv.h"
#include "ildpdf.h"
#include <gsl/gsl_multimin.h>
#if USE_BFGS
#include <lbfgs.h>
#endif

enum	Optim {NOOP, BFGS, GSL, MULT};
enum	GslAlg {ASIS, CONJ_FR, CONJ_PR, BFGS1, BFGS2, NMSIMP1, NMSIMP2, NMRAND};

typedef	double (*Dbl2DblFunc)(double);

static	const	double Frechet_mn_step[4] = {10., 10., 1., 0.1};
static	const	double LogNormal_mn_step[3] = {1., 0.1, 0.1};
static	const	double Geometric_mn_step[3] = {0.01, 10, 0.1};
static	double	Initial_mn_step[4] = {1., 1., 1., 1.}; 
static	const	GslAlg	gslalg[] = {NMSIMP1, CONJ_FR, ASIS};
static	const	int	REDUCE_DIM = 33;
static	const	double	Min_fraction = 0.02;
static	const	double	long_quantile = 0.995;

static	int	verbs = 0;
static	char*	specid = 0;
static	int	nDistParam = NoDistParam[0];
static	int	minsamples = 1000;
static	double	max_kappa = 30;
static	int	given_n_modes = 0;

class GslOptim {
	GslAlg	fdf_alg;
	gsl_vector*	vf;
	gsl_vector*	vg;
	gsl_vector*	ss;
	gsl_multimin_function_fdf	fdf_func;
	gsl_multimin_function	nm_func;
	void*	instance;
	double*	gparam;
	double	fx;
	int	status;
	int	niter;
	int	nparam;
	int	mode_size;
	void	fdf_iterate(const gsl_multimin_fdfminimizer_type* T, int max_iter);
	void	nm_iterate(const gsl_multimin_fminimizer_type* T, int max_iter);
public:
	GslOptim(void* inst, double* p, int n, GslAlg alg);
	~GslOptim() {
	    gsl_vector_free(vf); gsl_vector_free(vg); gsl_vector_free(ss);
	}
	double	optimize(int max_iter, double* p);
	int	niteration() {return (niter);}
	int	status_code() {return (status);}
};

static double evaluate(void* instance, const double* p, double* g,
	const int  n, const double step = 1.);
static	double	gsl_f(const gsl_vector * x, void * inst);
static	void	gsl_df(const gsl_vector * x, void * inst, gsl_vector * g);
static	void	gsl_fdf(const gsl_vector * x, void * inst, double* fx, gsl_vector * g);

Verbout	vbs;

void Verbout::printout(int n, const double* p, double fx, const char* spc)
{
	if (verbs & pro_phase) {
	    fprintf(stderr, "%4d", ++iter_no);
	    for (int k = 0; k < n; ++k) fprintf(stderr, " %7.2f", *p++);
	    fprintf(stderr, ": %15.7e", fx);
	    if (spc) fprintf(stderr, "\t%s", spc);
	    if (specid) fprintf(stderr, "\t%s", specid);
	    fputc('\n', stderr);
	}
}

/*************************************************************************************
*
*	maximize log likelihood
*	ll = Sigma_i f_i * log(Freche(x_i, theta))
*	p[] = {m1, t1, k1, a1, m2, t2, k2, a2, m3, t3, k3}
*
*************************************************************************************/

static double evaluate(void* instance, const double* p, double* g,
	const int  n, const double step)
{ 
	double	a = 1.;			// amplitude
	for (int j = 0; j < n; j += nDistParam) {
	    if ((j + nDistParam) < n) {	// not_lastj
		const double&	aj = p[j + nDistParam - 1];
		if (aj < 0.) return (OutOfRange);
		a -= aj;
	    }
	    int	q = j + (nDistParam == 4? 1: 0);
	    if (p[q] <= 0.) return (OutOfRange);
	    if (defpdf == GEOMETRIC) {
	        if (p[q] >= 1.) return (OutOfRange);
	    } else if (p[q + 1] < 0.) return (OutOfRange);
	}
	if (a < 0.) return (OutOfRange);
	Ild*	inst = (Ild*) instance;
	double	fx = 0.0;
	LenFrq* ptr = inst->intlf;
	if (g) vclear(g, n);
	double*	numer = g? new double[n + 1]: 0;
	for (int i = 0; i < inst->ntotal; ++i, ++ptr) {
	    double	len = ptr->len;
	    double	frq = ptr->frq;
	    double	denom = 0;	// Sigma_j p_j
	    if (numer) vclear(numer, n + 1);	// Sigma_j dp_j/dt
	    for (int j = 0; j < n; j += nDistParam) {	// mu, th, kp, al
		bool	not_lastj = (j + nDistParam) < n;
		const double&	mu = p[j];
		if (defpdf == FRECHET) {
		    const double&	al = not_lastj? p[j + 3]: a;
		    const double&	th = p[j + 1];
		    const double&	kp = p[j + 2];
		    if (len - mu < minl_mu || kp > max_kappa) {
			denom += al * SmallProb;
			continue;
		    }
		    double	z = th / (len - mu);
		    double	w = pow(z, kp);
		    double	f = errno? SmallProb: frechet(z, w, th, kp);
		    errno = 0;
		    denom += al * f;
		    if (numer) {
			numer[j + 0] += al * f * lfrechet_m(z, w, th, kp);
			numer[j + 1] += al * f * lfrechet_t(z, w, th, kp);
			numer[j + 2] += al * f * lfrechet_k(z, w, th, kp);
			numer[j + 3] += f;
		   }
		} else if (defpdf == LOGNORMAL) {
		    const double&	al = not_lastj? p[j + 2]: a;
		    const double&	sg = p[j + 1];
		    if (sg <= epsilon) {
			denom += al * SmallProb;
			continue;
		    }
		    double	z = (log(len) - mu) / sg;
		    double	f = lognormal(len, z, mu, sg);
		    denom += al * f;
		    if (numer) {
			numer[j + 0] += al * f * llognormal_m(len, z,  mu, sg);
			numer[j + 1] += al * f * llognormal_s(len, z,  mu, sg);
			numer[j + 2] += f;
		    }
		} else if (defpdf == GEOMETRIC) {
		    const double&	al = not_lastj? p[j + 2]: a;
		    const double&	q = p[j];
		    const double&	d = p[j + 1];
		    if (len < d || q < 0. || q > 1.) {
			denom += al * SmallProb;
			continue;
		    }
		    double	f = geometric(len, q, d);
		    denom += al * f;
		    if (numer) {
			numer[j + 0] += al * f * lgeometric_q(len, q, d);
			numer[j + 1] += al * f * lgeometric_d(len, q, d);
			numer[j + 2] += f;
		    }
		} else if (defpdf == GAMMA) {
		    const double&	al = not_lastj? p[j + 3]: a;
		    const double&	th = p[j + 1];
		    const double&	kp = p[j + 2];
		    if (len - mu < minl_mu || kp > max_kappa) {
			denom += al * SmallProb;
			continue;
		    }
		    double	z = (len - mu) / th;
		    double	f = errno? SmallProb: gamma(z, th, kp);
		    errno = 0;
		    denom += al * f;
		    if (numer) {
			numer[j + 0] += al * f * lgamma_m(z, th, kp);
			numer[j + 1] += al * f * lgamma_t(z, th, kp);
			numer[j + 2] += al * f * lgamma_k(z, th, kp);
			numer[j + 3] += f;
		   }
		}
	    }
	    if (denom > 0) {
		fx -= frq * log(denom);
		if (numer) {
		    for (int j = 0; j < n; ++j) {
			if (j % nDistParam == 3) numer[j] -= numer[n];
			g[j] -= frq * numer[j] / denom;
		    }
		}
	    }
	}
	delete[] numer;
	vbs.printout(n, p, fx);
	return fx;
}

IldPrm::IldPrm(const char* ip_stat, Ild* inst, const char* choose)
	: ildpdf(FRECHET), m_size(NoDistParam[ildpdf]), vrtl(false), 
	  sname(0), n_param(0), sid(0), n_sample(0), k_th_qtil(0), cdf_table(0)
{
	FILE*	fd = ftable.fopen(ip_stat, "r");
	if (!fd) fatal("%s not found !\n", ip_stat);
	double	min_aic = FLT_MAX;
	char	str[MAXL];
	char	genspc[MAXL] = "";
	double	tparam[paramsize];
	StatDist	savepdf = defpdf;
	defpdf = ildpdf;
	while (fgets(str, MAXL, fd)) {
	    if (*str == '#') {
		char*	ps = str + 1;
		while (*ps && isspace(*ps)) ++ps;
		if (!wordcmp(ps, "LOGNORMAL")) {
		    defpdf = ildpdf = LOGNORMAL;
		    nDistParam = m_size = NoDistParam[ildpdf];
		} else if (!wordcmp(ps, "GEOMETRIC")) {
		    defpdf = ildpdf = GEOMETRIC;
		    nDistParam = m_size = NoDistParam[ildpdf];
		}
		continue;
	    } else if (choose && wordcmp(str, choose)) continue;
	    int	n = get_IldPrm(str, minsamples, tparam);
	    if (!n) continue;
	    if (given_n_modes && given_n_modes != n_modes) continue;
	    double	halfaic = n + evaluate((void*) inst, tparam, 0, n);
	    if (0 < halfaic && halfaic < min_aic) {
		vcopy(dparam, tparam, n);
		car(genspc, str);
		n_param = n;
		min_aic = halfaic;
	    }
	}
	fclose(fd);
	if (!n_param) fatal("no entry in %s !\n", ip_stat);
	sname = new char[strlen(genspc) + 1];
	strcpy(sname, genspc);
	complete();
	defpdf = savepdf;
}

/*************************************************************************************
*
*	Use Gsl functions for optimization
*
*************************************************************************************/

static double gsl_f(const gsl_vector *v, void *inst)
{
	double*	x = new double[v->size];
	for (INT i = 0; i < v->size; ++i)
	     x[i] = gsl_vector_get(v, i);
	double fx = (double) evaluate(inst, x, 0, v->size);
	delete[] x;
	return (fx);
}

static void gsl_df(const gsl_vector *vx, void *inst, gsl_vector *df)
{
	double*	x = new double[vx->size];
	double*	g = new double[vx->size];
	for (INT i = 0; i < vx->size; ++i)
	     x[i] = gsl_vector_get(vx, i);
	evaluate(inst, x, g, vx->size);
	for (INT i = 0; i < vx->size; ++i)
	     gsl_vector_set(df, i, g[i]);
	delete[] x;
	delete[] g;
}

static void gsl_fdf(const gsl_vector *vx, void *inst, double *fx, gsl_vector *df)
{
	double*	x = new double[vx->size];
	double*	g = new double[vx->size];
	for (INT i = 0; i < vx->size; ++i)
	     x[i] = gsl_vector_get(vx, i);
	*fx = evaluate(inst, x, g, vx->size);
	for (INT i = 0; i < vx->size; ++i)
	     gsl_vector_set(df, i, g[i]);
	delete[] x;
	delete[] g;
}

GslOptim::GslOptim(void* inst, double* p, int n, GslAlg alg) 
	: fdf_alg(alg), instance(inst), gparam(p), nparam(n), 
	  mode_size(NoDistParam[defpdf])
{
	vf = gsl_vector_alloc(n);
	vg = gsl_vector_alloc(n);
	ss = gsl_vector_alloc(n);
	for (int i = 0; i < n; ++i) gsl_vector_set(vf, i, p[i]);
	fdf_func.n = nm_func.n = n;
	fdf_func.f = nm_func.f = &gsl_f;
	fdf_func.df = &gsl_df;
	fdf_func.fdf = &gsl_fdf;
	fdf_func.params = nm_func.params = inst;
	gsl_set_error_handler_off();
}

void GslOptim::fdf_iterate(const gsl_multimin_fdfminimizer_type* T, int max_iter)
{
	gsl_multimin_fdfminimizer*
	    s = gsl_multimin_fdfminimizer_alloc(T, nparam);
	gsl_multimin_fdfminimizer_set(s, &fdf_func, vf, gslprm.mn_step, gslprm.mn_conv);
	for (niter = 0; niter < max_iter; ++niter) {
	    status = gsl_multimin_fdfminimizer_iterate(s);
	    if (status) break;
	    status = gsl_multimin_test_gradient(s->gradient, gslprm.mn_grad);
	    if (status != GSL_CONTINUE) break;
	}
	if (status != GSL_FAILURE) {
	    for (int i = 0; i < nparam; ++i)
		gparam[i] = gsl_vector_get(s->x, i);
	}
	fx = gsl_multimin_fdfminimizer_minimum(s);
	gsl_multimin_fdfminimizer_free(s);
}

void GslOptim::nm_iterate(const gsl_multimin_fminimizer_type* T, int max_iter)
{
	vbs.set_phase(4);
	for (int i = 0; i < nparam; ++i)
	    gsl_vector_set(ss, i, Initial_mn_step[i % mode_size]);
	gsl_multimin_fminimizer*
	    s = gsl_multimin_fminimizer_alloc(T, nparam);
	gsl_multimin_fminimizer_set(s, &nm_func, vf, ss);
	for (niter = 0; niter < max_iter; ++niter) {
	    status = gsl_multimin_fminimizer_iterate(s);
	    if (status) break;
	    double	size = gsl_multimin_fminimizer_size(s);
	    status = gsl_multimin_test_size(size, gslprm.mn_step);
	    if (status != GSL_CONTINUE) break;
	}
	double	a = (nparam > nDistParam)? gsl_vector_get(s->x, nDistParam - 1): 1;
	double	t = gsl_vector_get(s->x, 1);
	if (a < 0. || 1. < a || t < 0.) status = GSL_FAILURE;
	else if (defpdf == FRECHET && gsl_vector_get(s->x, 2) < 0.)
	    status = GSL_FAILURE;
	else {
	    for (int i = 0; i < nparam; ++i)
		gparam[i] = gsl_vector_get(s->x, i);
	}
	fx = gsl_multimin_fminimizer_minimum(s);
	gsl_multimin_fminimizer_free(s);
}

double GslOptim::optimize(int max_iter, double* p)
{
	vbs.set_phase(4);
	switch (fdf_alg) {
	    case CONJ_FR:
		fdf_iterate(gsl_multimin_fdfminimizer_conjugate_fr, max_iter); break;
	    case CONJ_PR:
		fdf_iterate(gsl_multimin_fdfminimizer_conjugate_pr, max_iter); break;
	    case BFGS1:
		fdf_iterate(gsl_multimin_fdfminimizer_vector_bfgs, max_iter); break;
	    case BFGS2:
		fdf_iterate(gsl_multimin_fdfminimizer_vector_bfgs2, max_iter); break;
	    case NMSIMP1:
		nm_iterate(gsl_multimin_fminimizer_nmsimplex, max_iter); break;
	    case NMSIMP2:
		nm_iterate(gsl_multimin_fminimizer_nmsimplex2, max_iter); break;
	    default:
		nm_iterate(gsl_multimin_fminimizer_nmsimplex2rand, max_iter); break;
	}
	return (fx);
}

int Ild::fitild(IldPrm* dp, double* pfx)
{
	IldPrm	reserve(*dp);
	int	status = -1;
	double	lfx = 0.;
	if (!pfx) pfx = &lfx;
	for (const GslAlg* palg = gslalg; *palg; ++palg) {
	    GslOptim	gslopt(this, dp->dparam, dp->n_param, *palg);
	    *pfx = gslopt.optimize(gslprm.opt_max_iter, dp->dparam);
	    status = gslopt.status_code();
	    if (status == GSL_SUCCESS || status == GSL_ENOPROG) break;
	    if (status == REDUCE_DIM && dp->reduce_dim()) {
		GslOptim gslopt1(this, dp->dparam, dp->n_param, *palg);
		gslopt1.optimize(gslprm.opt_max_iter, dp->dparam);
		status = gslopt1.status_code();
		if (status == GSL_SUCCESS || status == GSL_ENOPROG) break;
	    }
	    dp->copy_from(reserve);
	}
	dp->complete();
	if (!dp->proper()) status = -1;
	return (status);
}

#ifdef MAIN

/**************************************************************************
*
*	Main functions
*
**************************************************************************/

void usage()
{
	fputs("Usage:\n\tfitild [options] -d IpTable xxx.ild\n", stderr);
	fputs(" or\tfitild [options] xxx.ild Init_params\n", stderr);
	fputs("\nOptions:\n", stderr);
	fputs("\t-h\t: display this help\n", stderr);
	fputs("\t-q\t: supress progress message\n", stderr);
	fputs("\t-a\t: as is: no optimization\n", stderr);
	fputs("\t-m\t: multiple methods (default)\n", stderr);
#if USE_BFGS
	fputs("\t-b\t: BFGS (Broyden-Fletcher-Golodfarb-Shanno)\n", stderr);
#endif
	fputs("\t-b[1|2]\t: BFGS optimization (GSL)\n", stderr);
	fputs("\t-f\t: conjugate_fr (Fletcher-Reeves)\n", stderr);
	fputs("\t-p\t: conjugate_pr (Polak-Ribi`ere)\n", stderr);
	fputs("\t-n[1|2]\t: nmsimplex (Nelder-Mead)\n\n", stderr);
	fputs("\t-c#\t: number of components (0: not specified)\n", stderr);
	fputs("\t-e#\t: degree of convergence (1e-4)\n", stderr);
	fputs("\t-i\t: improve params of the same species id\n", stderr);
	fputs("\t-k#\t: maxmum allowd kapper value (30)\n", stderr);
	fputs("\t-r#\t: max number of iterations (1000)\n", stderr);
	fputs("\t-s#\t: step size in GSL minimization function (1e-2)\n", stderr);
	fputs("\t-H#\t: skip template with less than # samples\n", stderr);
	fputs("\t-L#\t: lower bound of data (0)\n", stderr);
	fputs("\t-U#\t: upper bound of data (inf)\n", stderr);
	fputs("\t-X\t: X-axis in linear scale (log10 scale)\n\n", stderr);
	fputs("\t-g[out]\t: graphic output\n", stderr);
	fputs("\t-l#\t: lower bound of plot (in log_10 (1) unless -X)\n", stderr);
	fputs("\t-u#\t: upper bound of plot (in log_10 (6) unless -X)\n", stderr);
	fputs("\t-x#\t: number of bins (100)\n\n", stderr);
	fputs("\t-ME\t: exponential (geometoric) model\n", stderr);
	fputs("\t-MF\t: Frechet model\n", stderr);
	fputs("\t-MG\t: Gamma model\n", stderr);
	fputs("\t-MN\t: log normal model\n", stderr);
	fputs("\nInit_Params:\n", stderr);
	fputs("\t1. m_1 t_1 k_1\n", stderr);
	fputs(" or\ta_1 m_1 t_1 k_1 m_2 t_2 k_2\n", stderr);
	fputs(" or\ta_1 m_1 t_1 k_1 m_2 t_2 k_2 a_2 m_3 t_3 k_3\n", stderr);
	fputs("\nExamples:\n", stderr);
	fputs("\tfitild -b2 -d IldModel.txt xxx.ild\n", stderr);
	fputs("\tfitild -n2 xxx.ild 0.8 -200 300 50 -150 300 3.5\n", stderr);
	fputs("\tfitild -g -a xxx.ild 0.8 -200 300 50 -150 300 3.5\n", stderr);
	fputs("\tfitild -gxxx.eps -d IldModel.txt xxx.ild\n", stderr);
	exit(0);
}

int main(int argc, const char *argv[])
{
	char	pdfchar = 'F';
	const	char*	gsl_alg = 0;
static	const	char*	ip_stat = 0;
static	Optim	optim = MULT;
static	const	char*	graph_out = 0;
	const	char*	pn = 0;
	GslAlg	alg = NMSIMP1;		// default
	bool	identifier = false;	// choose template by id
	bool	x_axis_in_log_scale = true;

	while (--argc && **++argv == '-') {
	  gsl_alg = *argv + 1;
	  switch (*gsl_alg) {
	    case 'a': optim = NOOP; break;
	    case 'b': 
		if (gsl_alg[1] == '1') {
		    optim = GSL; alg = BFGS1;
		} else if (gsl_alg[1] == '2') {
		    optim = GSL; alg = BFGS2;
		}
#if USE_BFGS
		else	optim = BFGS;
#endif
		break;
	    case 'c':
		if ((pn = getarg(argc, argv, true)))
		    gslprm.mn_conv = atof(pn);
		break;
	    case 'd': ip_stat = getarg(argc, argv); break;
	    case 'f': optim = GSL; alg = CONJ_FR; break;
	    case 'g': graph_out = getarg(argc, argv);
		if (!graph_out) graph_out = "";
		else if (!gnuplot_terminal(graph_out)) {
		    ++argc; --argv; graph_out = "";
		}
		break;
	    case 'h': case '?': usage();
	    case 'i': identifier = true; break;
	    case 'j':
		if ((pn = getarg(argc, argv, true)))
		    given_n_modes = atoi(pn);
		break;
	    case 'k':
		if ((pn = getarg(argc, argv, true)))
		    max_kappa = atof(pn);
		break;
	    case 'l':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.llmt = atof(pn);
		break;
	    case 'm': optim = MULT; break;
	    case 'n': optim = GSL; 
		alg = gsl_alg[1] == '1'? NMSIMP1: NMSIMP2; break;
	    case 'p': optim = GSL; alg = CONJ_PR; break;
	    case 'q': setprompt(false, false); break;
	    case 'r':
		if ((pn = getarg(argc, argv, true)))
		    gslprm.opt_max_iter = atoi(pn);
		break;
	    case 's':
		if ((pn = getarg(argc, argv, true)))
		    gslprm.mn_step = atof(pn);
		break;
	    case 'u':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.ulmt = atof(pn);
		break;
	    case 'x':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.ndiv = atoi(pn);
		break;
	    case 'H':
		if ((pn = getarg(argc, argv, true)))
		    minsamples = atoi(pn);
		break;
	    case 'L':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.minx = atoi(pn);
		break;
	    case 'M':	// statistical model
	     switch (gsl_alg[1]) {
	      case 'E': defpdf = GEOMETRIC; pdfchar = 'E'; break;
	      case 'F': defpdf = FRECHET; pdfchar = 'F'; break;
	      case 'G': defpdf = GAMMA; pdfchar = 'G'; break;
	      case 'N': defpdf = LOGNORMAL; pdfchar = 'N'; break;
	      case 'W': defpdf = WEIBULL; pdfchar = 'W'; break;
	     }
	     break;
	    case 'P':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.psdcnt = atof(pn);
		break;
	    case 'V':
		verbs = (pn = getarg(argc, argv, true))? atoi(pn): 3;
		break;
	    case 'U':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.maxx = atoi(pn);
		 break;
	    case 'X': x_axis_in_log_scale = false; break;
	    default: break;
	  }
	}

// read date frrm file
const	char*	flename = *argv++;
	if (argc <= 0 || !flename) usage();
	const	char*	dot = strchr(flename, '.');
	int	slen = dot? dot - flename: strlen(flename);
	specid = new char[slen + 1];
	strncpy(specid, flename, slen);
	specid[slen] = '\0';
	Ild	intdist(flename);
	if (!intdist.ntotal) return 1;
	double	scale_position_step = 1.98 * log10(intdist.maxx()) - 1.88;
	if (defpdf == FRECHET || defpdf == GAMMA || defpdf == WEIBULL) {
	    vcopy(Initial_mn_step, Frechet_mn_step, 4);
	    Initial_mn_step[0] = Initial_mn_step[1] = scale_position_step;
	} else if (defpdf == GEOMETRIC) {
	    vcopy(Initial_mn_step, Geometric_mn_step, 3);
	    Initial_mn_step[1] = scale_position_step;
	} else if (defpdf == LOGNORMAL) {
	    vcopy(Initial_mn_step, LogNormal_mn_step, 3);
	}

// set initial values
	vbs.set_phase(1);
	IldPrm*	dp = ip_stat? new IldPrm(ip_stat, &intdist, identifier? specid: 0):
		(--argc? new IldPrm(argc, argv, defpdf):// from command line
		new IldPrm(&intdist, defpdf));		// single component
	nDistParam = NoDistParam[defpdf];

	if (!dp->convert(defpdf)) fatal("Fail to convert !\n");
	double*    resrv = new double[dp->n_param];
	vcopy(resrv, dp->dparam, dp->n_param);	// reserve initial params
	double	fx = 0;
	int	status = GSL_FAILURE;
	int	niter = 0;
	if (optim == MULT) {
	    for (const GslAlg* palg = gslalg; *palg; ++palg) {
		GslOptim	gslopt(&intdist, dp->dparam, dp->n_param, *palg);
		fx = gslopt.optimize(gslprm.opt_max_iter, dp->dparam);
		status = gslopt.status_code();
		niter = gslopt.niteration();
		dp->complete();
		if ((status == GSL_SUCCESS || status == GSL_ENOPROG) && dp->proper()) break;
		vcopy(dp->dparam, resrv, dp->n_param);
	    }
	    prompt("%s: GSL optimization terminated with status code = %d after %d iterations\n",
		specid, status, niter);
	}
#if USE_BFGS
	else if (optim == BFGS) {
// Initialize the parameters for the L-BFGS optimization
	    lbfgs_parameter_t param;
	    lbfgs_parameter_init(&param);

//    Start the L-BFGS optimization
	    status = lbfgs(dp->n_param, dp->dparam, &fx, evaluate, NULL, &intdist, &param);
	    prompt("%s: L-BFGS optimization terminated with status code = %d\n", specid, status);
	    if (status == GSL_FAILURE) vcopy(dp->dparam, resrv, dp->n_param);
	}
#endif
	else if (optim == GSL) {
	    GslOptim	gslopt(&intdist, dp->dparam, dp->n_param, alg);
	    fx = gslopt.optimize(gslprm.opt_max_iter, dp->dparam);
	    status = gslopt.status_code();
	    prompt("%s: GSL optimization terminated with status code = %d after %d iterations\n", 
		specid, status, gslopt.niteration());
	    if (status == GSL_FAILURE) vcopy(dp->dparam, resrv, dp->n_param);
	}
	delete[] resrv;
	dp->complete();
	if (dp->proper()) {
	  double	q999 = intdist.quantile(long_quantile);
	  double	f999 = dp->quantile(long_quantile, q999);
	  fx = evaluate(&intdist, dp->dparam, 0, dp->n_param, 0);
	  Dbl2DblFunc	enf = 0;
	  Dbl2DblFunc	def = 0;
	  if (x_axis_in_log_scale) {
	    enf = log10;
	    def = exp10;
	  }
	  printf("%-10s %c\t%7d\t%7d\t%7.1f\t%7.0f\t", specid, pdfchar, int(intdist.ftotal),
	   intdist.intlf->len, intdist.quantile(0.5), (0 < f999 && f999 < 1.e6)? f999: q999);
	  printf("%7.4f\t", dp->negentropy_10());
	  dp->printprm(stdout);
//	  residual root mean suare deviation, AIC, BIC
	  printf(" %11.4e %11.4e %11.4e\n", intdist.rmsd(dp),
	     2 * (dp->n_param + fx), dp->n_param * log(intdist.ftotal) + 2 * fx);
	  if (graph_out) {
	    Lild	lild(intdist, enf, def);
	    GnuPlotLild	gp(&lild, dp);
	    gp.plot(*graph_out? graph_out: 0);
	  }
	} else {
	  printf("%-10s %c\t%7d: Error !\n", specid, pdfchar, int(intdist.ftotal));
	}
	delete dp;
	delete[] specid;
	return 0;
}

#endif	// MAIN

