/*****************************************************************************
*
*	<<< ildpdf.c >>>
*
*	probabilistic distribution functions
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

#include "ildpdf.h"

LildPrm	lildprm = {0, 100, 0, INT_MAX - 1, 0, 4, 12, 1., 6., 1e-5, 0.999};
//LildPrm	rildprm = {0, 100, INT_MAX - 1, 0, 4, 12, 0., 0.04, 1e-5, 0.999};
GslPrm	gslprm = {1000, 10000, 100, 0., 1.e-3, 0.01, 1.e-4, 1.e-3};
IldOutMode	ildoutmode = IldOutMode::PDF;

struct  CdfParam {IldPrm* ip; double y;};
static	SpecList*	speclist = 0;

StatDist	defpdf = FRECHET;
static	const	int	param_order3[] = {2, 0, 1, 3, 4, 5, 6, 7};
static	const	int	param_order4[] = {3, 0, 1, 2, 4, 5, 6, 7, 8, 9, 10};
static	const	int*	param_order[2] = {param_order3, param_order4};

static	const	char*	def_form = "EPS";
static	const	int	n_colors = 110;
static	const	char*	colors[] = {
"white", "black", "dark-grey", "red", "light-green",
"light-blue", "dark-magenta", "dark-cyan", "dark-orange", "dark-yellow", "royalblue",
"goldenrod", "dark-green", "purple", "steelblue", "dark-red", "dark-chartreuse", 
"orchid", "aquamarine", "brown", "yellow", "turquoise", "grey0", "grey10", "grey20",
"grey30", "grey40", "grey50", "grey60", "grey70", "grey", "grey80", "grey90", "grey100",
"light-red", "light-green", "light-blue", "light-magenta", "light-cyan",
"light-goldenrod", "light-pink", "light-turquoise", "gold", "green", "dark-green",
"spring-green", "forest-green", "sea-green", "blue", "dark-blue", "midnight-blue",
"navy", "medium-blue", "skyblue", "cyan", "magenta", "dark-turquoise", "dark-pink",
"coral", "light-coral", "orange-red", "salmon", "dark-salmon", "khaki", "dark-khaki",
"dark-goldenrod", "beige", "olive", "orange", "violet", "dark-violet", "plum",
"dark-plum", "dark-olivegreen", "orangered4", "brown4", "sienna4", "orchid4"
"mediumpurple3", "slateblue1", "yellow4", "sienna1", "tan1", "sandybrown",
"light-salmon", "pink", "khaki1", "lemonchiffon", "bisque", "honeydew", "slategrey",
"seagreen", "antiquewhite", "chartreuse", "greenyellow", "gray", "light-gray",
"light-grey", "dark-gray", "slategray", "gray0", "gray10", "gray20", "gray30",
"gray40", "gray50", "gray60", "gray70", "gray80", "gray90", "gray100"
};

static  double  cdf_func(double x, void* param);
static  double  dcdf_func(double x, void* param);
static  void    fdf_func(double x, void* param, double *y, double *dy);
static	double	rqrl_func(double x, void* param);
static	double	drqrl_func(double x, void* param);
static  void    frqrl_func(double x, void* param, double *y, double *dy);
static	int	mark_color = 9;

SpecList::SpecList(int ac, const char** av) 
	: StrHash<int>(std::max(ac, min_hash))
{
	char	str[MAXL];
	const char**        av0 = av;
	while (ac-- && **av != '-') {
	    strcpy(str, *av++);
	    char* dot = strrchr(str, '.');
	    char* sls = strrchr(str, '/');
	    if (dot) *dot = '\0';
	    incr(sls? sls + 1: str);
	}
	n_spec = av - av0;
}

SpecList* make_speclist(int nn, const char** mm) {
	delete speclist;
	return (speclist = new SpecList(nn, mm));
}

void delete_speclist() {
	delete speclist; speclist = 0;
}

/***********************************************************************************
*
*	Intron length distribution represented by 1-3 componets of Frechet distributions
*
***********************************************************************************/

IldPrm::IldPrm(int np, const char** argv, StatDist pdf)
	: ildpdf(pdf), m_size(NoDistParam[ildpdf]), n_param(np)
{
const	int*	order = param_order[m_size - 3];
	for (int i = 0; i < n_param; ++i)
	    dparam[order[i]] = atof(argv[i]);
	if (n_param == m_size) --n_param;
	complete();
}

IldPrm::IldPrm(Ild* inst, StatDist pdf)
	: ildpdf(pdf), n_modes(1), m_size(NoDistParam[ildpdf]), 
	  n_param(m_size - 1), n_sample(int(inst->ftotal))
{
	dparam[n_param] = 1.;
	if (ildpdf == GEOMETRIC) {
	    dparam[1] = inst->minx();
	    double	m = inst->mean();
	    double	p = m - dparam[1] + 1;
	    dparam[0] = 2. / (p + sqrt(p * p - 4 * m));
	} else if (ildpdf == LOGNORMAL) {
	    inst->logmean(dparam);
	} else if (ildpdf == GAMMA) {
	    double	v = 1.;
	    double	m = inst->mean(&v);
	    dparam[0] = 0;		// gamma, mu
	    dparam[1] = v / m;		// beta, th
	    dparam[2] = m * m / v;	// alpha, ki
	} else {
	    double	xx[4] = {
		inst->quantile(0.25),
		inst->quantile(0.5),
		inst->quantile(0.75),
		0.
	    };
	    single_mode(xx);
	}
	if (lildprm.n_qtl) calc_qtl(lildprm.n_qtl);
}

int IldPrm::get_IldPrm(char* ps, int minsamples, double* tparam)
{
	if (!tparam) tparam = dparam;
	Strlist	terms(ps, stddelim);
	int	n_terms = terms.size() - n_skip_term - 3;	// 3 indicators
	if (n_terms < m_size) return (0);	// no param
	n_sample = atoi(terms[2]);
	if (n_sample < minsamples) return (0);  // unreliable
	min_x = atoi(terms[3]);
	max_x = atoi(terms[5]);
const	int*	order = param_order[m_size - 3];
	int     n = 0;
	for ( ; n < n_terms; ++n) {
	    char*	tm = terms[n + n_skip_term];
	    if (strchr(tm, 'n') || strchr(tm, 'i')) return (0);
	    tparam[order[n]] = atof(tm);
	}
	n_modes = (n + 1) / m_size;
	if (tparam == dparam) n_param = n_modes * m_size - 1;
	return (n_modes * m_size - 1);
}

IldPrm::IldPrm(const char* ip_stat, const char* gs)
	: ildpdf(FRECHET), m_size(NoDistParam[ildpdf])
{
	FILE*   fd = fopen(ip_stat, "r");
	if (!fd) fatal("%s not found !\n", ip_stat);
	char	str[MAXL];
	int	n = 0;
	while (fgets(str, MAXL, fd)) {
	    if (*str == '#') {
		char*	ps = str + 1;
		while (*ps && isspace(*ps)) ++ps;
		if (!wordcmp(ps, "LOGNORMAL")) {
		    ildpdf = LOGNORMAL;
		    m_size = NoDistParam[ildpdf];
		} else if (!wordcmp(ps, "GEOMETRIC")) {
		    ildpdf = GEOMETRIC;
		    m_size = NoDistParam[ildpdf];
		}  else if (!wordcmp(ps, "GAMMA")) {
		    ildpdf = GAMMA;
		    m_size = NoDistParam[ildpdf];
		}  else if (!wordcmp(ps, "WEIBULL")) {
		    ildpdf = WEIBULL;
		    m_size = NoDistParam[ildpdf];
		}
		continue;
	    } else if (wordcmp(str, gs)) continue;
	    if ((n = get_IldPrm(str))) break;
	}
	fclose(fd);
	if (n) {
	    complete();
	    sname = new char[strlen(gs) + 1];
	    strcpy(sname, gs);
	} else {
	    prompt("%s not found in %s\n", gs, ip_stat);
	}
}

double IldPrm::pdf_function(double len, int ildent)
{
	double	f = 0;	// Sigma_j p_j
	for (int j = 0, m = 0; j < n_param; j += m_size) {
	    if (ildent && ++m != ildent) continue;
	    const double&	mu = dparam[j];
	    if (ildpdf == FRECHET) {
		if (len - mu < minl_mu) continue;
		const double&	th = dparam[j + 1];
		const double&	kp = dparam[j + 2];
		const double&	al = dparam[j + 3];
		if (al <= 0.) continue;
		double	z = th / (len - mu);
		double	w = pow(z, kp);
		f += al * frechet(z, w, th, kp);
	    } else if (ildpdf == LOGNORMAL) {
		const double&	sg = dparam[j + 1];
		const double&	al = dparam[j + 2];
		if (al <= 0.) continue;
		double	z = (log(len) - mu) / sg;
		f += al * lognormal(len, z, mu, sg);
	    } else if (ildpdf == GEOMETRIC) {
		const double&   q = dparam[j];
		const double&	d = dparam[j + 1];
		const double&	al = dparam[j + 2];
		if (al <= 0.) continue;
		f += al * geometric(len, q, d);
	    } else if (ildpdf == GAMMA) {
		const double&	th = dparam[j + 1];
		const double&	kp = dparam[j + 2];
		const double&	al = dparam[j + 3];
		if (al <= 0.) continue;
		double	z = (len - mu) / th;
		f += al * gamma(z, th, kp);
	    } else if (ildpdf == WEIBULL) {
		const double&	th = dparam[j + 1];
		const double&	kp = dparam[j + 2];
		const double&	al = dparam[j + 3];
		if (al <= 0.) continue;
		double	z = (len - mu) / th;
		double	w = pow(z, kp);
		f += al * weibull(z, w, th, kp);
	    }
	}
	return f;
}

double IldPrm::cumulative(double len, int ildent)
{
	double	f = 0;
	for (int j = 0, m = 0; j < n_param; j += m_size) {
	    if (ildent && ++m != ildent) continue;
	    const double&	mu = dparam[j];
	    if (ildpdf == FRECHET) {
		if (len - mu < minl_mu) continue;
		const double&	th = dparam[j + 1];
		const double&	kp = dparam[j + 2];
		const double&	al = dparam[j + 3];
		if (al <= 0.) continue;
		double	z = pow(th / (len - mu), kp);
		f += al * exp(-z);
	    } else if (ildpdf == LOGNORMAL) {
		const double&	sg = dparam[j + 1];
		const double&	al = dparam[j + 2];
		if (al <= 0.) continue;
		double	z = (log(len) - mu) / sg;
		f += al * 0.5 * (1. + erf(z / sqrt(2.)));
	    } else if (ildpdf == GEOMETRIC) {
		const double&   q = dparam[j];
		const double&	d = dparam[j + 1];
		const double&	al = dparam[j + 2];
		if (al <= 0.) continue;
		f += al * (1. - pow(1. - q, len - d));
	    } else if (ildpdf == GAMMA) {
		const double&	th = dparam[j + 1];
		const double&	kp = dparam[j + 2];
		const double&	al = dparam[j + 2];
		if (al <= 0.) continue;
		double	z = (len - mu) / th;
		f += al * gsl_sf_gamma_inc_P(kp, z);
	    } else if (ildpdf == WEIBULL) {
		const double&	th = dparam[j + 1];
		const double&	kp = dparam[j + 2];
		const double&	al = dparam[j + 2];
		if (al <= 0.) continue;
		double	z = (len - mu) / th;
		f += al * (1 - exp(-pow(z, kp)));
	    }
	}
	return (f);
}

double IldPrm::mean()
{
	double	f = 0;	// Sigma_j p_j
	for (int j = 0; j < n_param; j += m_size) {	// mu, th, kp, al
	    if (ildpdf == FRECHET) {
		f += dparam[j + 3] * frechet_mean(dparam + j);
	    } else if (ildpdf == LOGNORMAL) {
		f += dparam[j + 2] * lognormal_mean(dparam + j);
	    } else if (ildpdf == GEOMETRIC) {
		f += dparam[j + 2] * geometric_mean(dparam + j);
	    } else if (ildpdf == GAMMA) {
		f += dparam[j + 3] * gamma_mean(dparam + j);
	    } else if (ildpdf == WEIBULL) {
		f += dparam[j + 3] * weibull_mean(dparam + j);
	    }
	}
	return f;
}

double IldPrm::quantile(double y, double x)
{
	if (x <= 0) x = mean();
	double	initial_x = x;
	CdfParam	cp = {this, y};
	gsl_function_fdf FDF;
	FDF.f = &cdf_func;
	FDF.df = &dcdf_func;
	FDF.fdf = &fdf_func;
	FDF.params = (void*) &cp;
	gsl_set_error_handler_off();
const   gsl_root_fdfsolver_type* T = gsl_root_fdfsolver_newton;
	gsl_root_fdfsolver*     s = gsl_root_fdfsolver_alloc(T);
	gsl_root_fdfsolver_set(s, &FDF, x);
	int     status = GSL_CONTINUE;
	for (int i = 0; status == GSL_CONTINUE && i < gslprm.rsl_max_iter; ++i) {
	    status = gsl_root_fdfsolver_iterate(s);
	    double      x0 = x;
	    x = gsl_root_fdfsolver_root(s);
	    status = gsl_root_test_delta(x, x0, gslprm.epsabs, gslprm.epsrel);
	}
	gsl_root_fdfsolver_free(s);
	return (x > 0? x: initial_x + 1);
}

// estimate single-mode parameters from (1/4, 1/2, 3/4) quantile values

int IldPrm::single_mode(double q123_4[])
{
	double	x = -0.5;	// -1/kappa
	double	rqrl_param[4];
	for (int i = 0; i < 3; ++i)
	    rqrl_param[i] = -log(double(i + 1) / 4);
	rqrl_param[3] = (q123_4[1] - q123_4[0]) / (q123_4[2] - q123_4[1]);
	gsl_function_fdf FDF;
	FDF.f = &rqrl_func;
	FDF.df = &drqrl_func;
	FDF.fdf = &frqrl_func;
	FDF.params = (void*) rqrl_param;
	gsl_set_error_handler_off();
const   gsl_root_fdfsolver_type* T = gsl_root_fdfsolver_newton;
	gsl_root_fdfsolver*     s = gsl_root_fdfsolver_alloc(T);
	gsl_root_fdfsolver_set(s, &FDF, x);
	int     status = GSL_CONTINUE;
	for (int i = 0; status == GSL_CONTINUE && i < gslprm.rsl_max_iter; ++i) {
	    status = gsl_root_fdfsolver_iterate(s);
	    double      x0 = x;
	    x = gsl_root_fdfsolver_root(s);
	    status = gsl_root_test_delta(x, x0, gslprm.epsabs, gslprm.epsrel);
	}
	gsl_root_fdfsolver_free(s);
	dparam[2] = ildpdf == WEIBULL? -x: -1 / x;
	dparam[1] = (q123_4[2] - q123_4[0]) / 
	    (pow(rqrl_param[2], x) - pow(rqrl_param[0], x));
	dparam[0] = q123_4[1] - dparam[1] * pow(rqrl_param[1], x);
	return (status);
}

bool IldPrm::proper()
{
	bool	rv = true;
	for (int j = 0; j < n_param; j += m_size) {
	    double	a = dparam[j + m_size - 1];
	    rv = rv && (0. <= a && a <= 1.);
	    if (ildpdf == FRECHET || ildpdf == GAMMA || ildpdf == WEIBULL) {
		rv = rv && dparam[j + 1] >= 0 && dparam[j + 2] >= 0;
	    } else if (ildpdf == LOGNORMAL) {
		rv = rv && dparam[j + 1] >= 0;
	    } else if (ildpdf == GEOMETRIC) {
		rv = rv && dparam[j] >= 0;
	    }
	}
	return rv;
}

char* IldPrm::strget(char* str, char** terms)
{
	clear();
	char*	ps = str;
	char*   qs = car(ps); if (!qs || !*++ps) return (0);
	if (speclist && !speclist->find(qs)) return (0);
	sname = new char[strlen(qs) + 1]; strcpy(sname, qs);
	if (!(qs = car(ps)) || !*++ps) return (0);
	if (!(qs = car(ps)) || !*++ps) return (0);
	if ((n_sample = atoi(qs)) < lildprm.minfreq) return (0);
	if (!(qs = car(ps)) || !*++ps) return (0);
	else {min_x = atoi(qs);}
	if (!(qs = car(ps)) || !*++ps) return (0);
	if (!(qs = car(ps)) || !*++ps) return (0);
	else {max_x = atoi(qs);}
	if (!(qs = car(ps)) || !*++ps) return (0);
	int	n = 0;
	int	nmax = max_modes * m_size - 1;
	char*	rs = ps;
const	int*	order = param_order[m_size - 3];
	while ((qs = car(ps))) {
	    if (n < nmax) dparam[order[n]] = atof(qs);
	    if (terms) terms[n] = qs;
	    ++n;
	    if (!*++ps) break;
	}
	n_modes = (n - 3 + 1) / m_size;
	n_param = n_modes * m_size - 1;
	complete();
	return (rs);
}

int IldPrm::fget(FILE* fd, const char* fn)
{
	char    str[MAXL];
	char*   ps = 0;
	while ((ps = fgets(str, MAXL, fd)))
	    if (*ps != '#') break;
	if (!ps) return (EOF);
	return (strget(str)? OK: IGNORE);
}

void IldPrm::complete()
{
	dparam[n_param] = 0;
	double	a = 0;
	double	medians[3] = {DBL_MAX, DBL_MAX, DBL_MAX};
	double	fraction[3] = {0., 0., 0.};
	int	n_comp = 0;
	for (int j = 0; j < n_param; ++n_comp, j += m_size) {
	    a += dparam[j + m_size - 1];
	    switch (ildpdf) {
		case FRECHET:
		    medians[n_comp] = frechet_median(dparam + j); break;
		case LOGNORMAL:
		    medians[n_comp] = lognormal_median(dparam + j); break;
		case GEOMETRIC:
		    medians[n_comp] = geometric_median(dparam + j); break;
		case GAMMA:
		    medians[n_comp] = gamma_mean(dparam + j); break;
		case WEIBULL:
		    medians[n_comp] = weibull_mean(dparam + j); break;
	    }
	    fraction[n_comp] = dparam[j + m_size - 1];
	}
	n_modes = n_comp;
	dparam[n_param] = fraction[n_comp - 1] = 1. - a;
// sort components on their median in ascending order
	while (--n_comp > 0) {
	    int	nn = m_size * n_comp; 
	    for (int i = 0; i < n_comp; ++i) {
		bool	swp = false;
		double	dif = 2 * (medians[i] - medians[n_comp]) 
			/ (medians[i] + medians[n_comp]);
		if (dif > 0.1) swp = true;
		else if (dif > -0.1)
		    swp = fraction[i] < fraction[n_comp];
		if (swp) {
		    int	ii = m_size * i;
		    for (int k = 0; k < m_size; ++k)
			swap(dparam[ii + k], dparam[nn + k]);
		}
	    }
	}
	if (lildprm.n_qtl) calc_qtl(lildprm.n_qtl);
}

int IldPrm::printprm(FILE* fd)
{
const	int*	order = param_order[m_size - 3];
	int	jj = n_modes == 1? m_size: n_param;
	for (int j = 0; j < jj; ++j)
	    fprintf(fd, " %11.4e", dparam[order[j]]);
	return (jj);
}

bool IldPrm::convert(StatDist tipdf)
{
	if (tipdf == ildpdf) return (true);
	if (ildpdf != FRECHET) return (false);
	if (tipdf != LOGNORMAL && tipdf != GEOMETRIC) return (false);
	double	tparam[paramsize];
	vcopy(tparam, dparam, paramsize);
	int	n_size = NoDistParam[tipdf];
	double*	dp = dparam;
	double*	tp = tparam;
	for (int j = 0; j < n_param; j += m_size, tp += m_size, dp += n_size) {
	    double	linear_mu = frechet_median(tp);
	    if (tipdf == LOGNORMAL) {
		dp[0] = log(linear_mu);
		dp[1] = (log(frechet_quantile(tp, 0.75)) - log(frechet_quantile(tp, 0.25))) 
			/ (2 * l14coef);
	    } else if (tipdf == GEOMETRIC) {
		double	qhalf = frechet_quantile(tp, 0.75) - frechet_quantile(tp, 0.25);
		double	p = exp(-log(3)/qhalf);
		dp[0] = 1 - p;
		dp[1] = linear_mu + log(2) / log(p);
	    }
	    dp[2] = tp[3];	// a
	}
	m_size = n_size;
	n_param = n_modes * m_size - 1;
	ildpdf = tipdf;
	return (true);
}

static double plogp(double len, void* p)
{
	IldPrm* dp = (IldPrm*) p;
	double  pv = dp->pdf_function(len);
	return pv > 0.? pv * log(pv): 0.;
}

double IldPrm::negentropy_10()
{
	gsl_function    gslfunc;
	gslfunc.function = plogp;
	gslfunc.params = (void*) this;
	double  result = 0;
	double  abserr= 0;
	gsl_set_error_handler_off();
	gsl_integration_workspace* ws = gsl_integration_workspace_alloc(gslprm.n_worksp);
	gsl_integration_qagiu(&gslfunc, 0., gslprm.epsabs, gslprm.epsrel, gslprm.n_worksp,
	ws, &result, &abserr);
	gsl_integration_workspace_free(ws);
	return (result / ln_10);
}

void IldPrm::calc_qtl(int nq)
{
	k_th_qtil = new float[nq];
	double	delta = 1. / (double) nq;
	double	y = delta;
	double	z = min_x;
	for (int k = 0; k < nq; ++k, y += delta) {
	    z = quantile(y, z);
	    k_th_qtil[k] = log(z);
	}
}

void IldPrm::calc_cdf()
{
	cdf_table = new double[max_x - min_x] - min_x;
	for (int x = min_x; x < max_x; ++x)
	    cdf_table[x] = cumulative(x);
}
	
/*************************************************************
*
*	Fit binned data to B-spline functions
*
*************************************************************/

Bspline::Bspline(PutIntoBins& bins, double* dt, int ko, int nc, int n_th_d)
	: k_order(ko), n_breaks(nc - ko + 2), n_coeffs(nc)
{
	bsws = gsl_bspline_alloc(ko, n_breaks);
#if GSL_MAJOR_VERSION == 1
	if (n_breaks) bsdws = gsl_bspline_deriv_alloc(ko);
#endif
	Bk = gsl_vector_alloc(nc);
	int	n_dt = bins.size();
	gsl_vector*	x = gsl_vector_alloc(n_dt);
	gsl_vector*	y = gsl_vector_alloc(n_dt);
	dB = gsl_matrix_alloc(n_dt, nc);
	coeff = gsl_vector_alloc(nc);
	cov = gsl_matrix_alloc(nc, nc);
	gsl_multifit_linear_workspace*	
	    mw = gsl_multifit_linear_alloc(n_dt, nc);
	for (int n = 0; n < n_dt; ++n) {
	    gsl_vector_set(x, n, bins.xval[n]);
	    gsl_vector_set(y, n, bins.freq[n]);
	}
	gsl_bspline_knots_uniform(bins.xval[0], bins.xval[n_dt - 1], bsws);
	for (int n = 0; n < n_dt; ++n) {
	    gsl_bspline_eval((double(bins.xval[n])), Bk, bsws);
	    for (int j = 0; j < nc; ++j) {
		double Bj = gsl_vector_get(Bk, j);
		gsl_matrix_set(dB, n, j, Bj);
	    }
	}
	double	chisq;
	gsl_multifit_linear(dB, y, coeff, cov, &chisq, mw);
	gsl_vector_free(x);
	gsl_vector_free(y);
	gsl_multifit_linear_free(mw);
}

double Bspline::fitted(double x)
{
	double	yval, yerr;
	gsl_bspline_eval(x, Bk, bsws);
	gsl_multifit_linear_est(Bk, coeff, cov, &yval, &yerr);
	return (yval);
}

double Bspline::derivative(double x)
{
	double	der1, yerr;
#if GSL_MAJOR_VERSION == 1
	gsl_bspline_deriv_eval(x, 1, dB, bsws, bsdws);
#else
	gsl_bspline_deriv_eval(x, 1, dB, bsws);
#endif
	for (int i = 0; i < n_coeffs; ++i)
	    gsl_vector_set(Bk, i, gsl_matrix_get(dB, i, 1));
	gsl_multifit_linear_est(Bk, coeff, cov, &der1, &yerr);
	return (der1);
}

/***********************************************************************************
*
*	Observed intron length distribution with x-axis in linear scale
*
***********************************************************************************/

static int lcompf(const LenFrq* a, const LenFrq* b)
{
	return (a->len - b->len);
}

void Ild::calc_qtl(int nq)
{
	k_th_qtil = new float[nq];
	double	sum = 0.;
	double	delta = 1. / (double) nq;
	double	y = delta;
	int	k = 0;
	for (int i = 0; i < ntotal; ++i) {
	    sum += intlf[i].frq;
	    if (sum >= y) {
		double	x1 = i? log((double) intlf[i - 1].len): 0.;
		double	x2 = log((double) intlf[i].len);
		k_th_qtil[k++] = x2 - (x2 - x1) * (sum - y) / intlf[i].frq;
		y += delta;
	    }
	}
}

Ild::Ild(const char* fn, int id)
	: nfact(0), sid(0), vrtl(false), ntotal(0), 
	  ftotal(0), fname(fn), intlf(0), k_th_qtil(0)
{
	FILE*	fd = fopen(fn,"r");
	char	str[MAXL];
	if (!fd) {
	    strcpy(str, fname);
	    strcat(str, ".ild");
	    fd = fopen(str, "r");
	    if (!fd) fatal("Can't open %s\n", str);
	}
	int	len;
	double	frq;
	while (fgets(str, MAXL, fd)) {
	    int	nr = sscanf(str,"%d %lf", &len, &frq);
	    if (nr == 2 && len > lildprm.minx && frq > 0 && len <= lildprm.maxx) {
		ftotal += frq; 
		++ntotal;
	    }
	}
	if (ntotal) {
	    intlf = new LenFrq[ntotal + 1];
	    rewind(fd);
	    LenFrq*	ptr = intlf;
	    while (fgets(str, MAXL, fd)) {
		int	frq;
		int	nr = sscanf(str, "%d %d", &ptr->len, &frq);
		if (nr == 2 && ptr->len > lildprm.minx && frq > 0 &&  ptr->len <= lildprm.maxx)
		    (ptr++)->frq = frq;
	    }
	}
	qsort((UPTR) intlf, ntotal, sizeof(LenFrq), (CMPF) lcompf);
	fclose(fd);
	if (lildprm.n_qtl) calc_qtl(lildprm.n_qtl);
}

/***********************************************************************************
*
*	Random ild distribution generator	
*
***********************************************************************************/

static	int cmpd(const void* a, const void* b)
{
	const	double*	da = (const double*) a;
	const	double*	db = (const double*) b;
	if (*da < *db) return (-1);
	if (*da > *db) return (1);
	return (0);
}

static	int cmpi(const int* a, const int* b)
{
	return (*a - *b);
}

#if M_THREAD
Ild::Ild(IldPrm* ildprm, int n_sample, Drand48_data* drnad_buff, const char* snm)
#else
Ild::Ild(IldPrm* ildprm, int n_sample, const char* snm)
#endif
	: nfact(0), sid(0), vrtl(false), ntotal(0), 
	  ftotal(0), fname(snm), intlf(0), k_th_qtil(0)
{
	double*	uniform = new double[n_sample];

// generage a uniform distribution [0, 1)

	for (int i = 0; i < n_sample; ++i) {
#if M_THREAD
	    drand48_r(drnad_buff, &uniform[i]);
#else
	    uniform[i] = drand48();
#endif
	}
	qsort((UPTR) uniform, (INT) n_sample, sizeof(double), (CMPF) cmpd);
	if (!ildprm->cdf_table) ildprm->calc_cdf();
	int	x = ildprm->min_x;
	int	n = 0;
	double	p = 0.;
	Dhash<int, double>	lf(ildprm->max_x, 0.);
	while (uniform[n] < ildprm->cdf_table[ildprm->min_x]) ++n;
	while (n < n_sample && x < ildprm->max_x) {
	    if (uniform[n] < ildprm->cdf_table[x]) {
		bool	lower_half = 2 * uniform[n] < p + ildprm->cdf_table[x];
		lf.incr(lower_half? x - 1: x);
		++n; ++ftotal;
	    } else {
		p = ildprm->cdf_table[x++];
	    }
	}
	double	dx = x;
	while (n < n_sample) {
	    dx = ildprm->quantile(uniform[n], dx);
	    x = int (dx + 0.5);
	    if (x > 2 * ildprm->max_x) break;
	    lf.incr(x);
	    ++n; ++ftotal;
	}
	intlf = (LenFrq*) lf.squeeze(&ntotal);
	qsort((UPTR) intlf, ntotal, sizeof(LenFrq), (CMPF) lcompf);
	if (lildprm.n_qtl) calc_qtl(lildprm.n_qtl);
	delete[] uniform;
}

// bootstrap resampling

#if M_THREAD
Ild::Ild(Ild* src, Drand48_data* drnad_buff)
#else
Ild::Ild(Ild* src)
#endif
	: nfact(0), sid(0), vrtl(false), ntotal(0), 
	  ftotal(src->ftotal), fname(0), intlf(0), k_th_qtil(0)
{
	INT	n_sample = (INT) src->ftotal;
	int*	ordinal = new int[n_sample];

// generage a uniform distribution [0, 1)
	for (INT i = 0; i < n_sample; ++i) {
	    double	u;
#if M_THREAD
	    drand48_r(drnad_buff, &u);
#else
	    u = drand48();
#endif
	    ordinal[i] = int(src->ftotal * u);
	}
	qsort((UPTR) ordinal, n_sample, sizeof(int), (CMPF) cmpi);
	INT	n = 0;
	double	cumulative = 0.;
	LenFrq*	tlf = src->intlf + src->ntotal;
	intlf = new LenFrq[src->ntotal];
	LenFrq*	nlf = intlf;
	for (LenFrq* slf = src->intlf; slf < tlf; ++slf) {
	    nlf->len = slf->len;
	    int	nn = 0;
	    cumulative += slf->frq;
	    while (n < n_sample && ordinal[n++] < cumulative) ++nn;
	    if (nn) (nlf++)->frq = nn;
	}
	ntotal = nlf - intlf;
	if (lildprm.n_qtl) calc_qtl(lildprm.n_qtl);
	delete[] ordinal;
}

Ild::~Ild()
{
	if (!vrtl) {delete[] intlf; delete[] k_th_qtil;}
}

int Ild::fget(FILE* fd, const char* fn)
{
	fname = fn;
	Mfile	mfd(sizeof(LenFrq));
	LenFrq	wrk;
	char	str[MAXL];
	for (ftotal = 0; fgets(str, MAXL, fd); ) {
	    char*	ps = str;
	    while (*ps && isspace(*ps)) ++ps;
	    if (!isdigit(*ps)) continue;
	    if (sscanf(ps, "%d %lf", &wrk.len, &wrk.frq) != 2) break;
	    if (wrk.len <= lildprm.minx || wrk.len > lildprm.maxx) continue;
	    ftotal += wrk.frq;
	    mfd.write(&wrk);
	}
	ntotal = mfd.size();
	wrk.len = INT_MAX; wrk.frq = 0; 
	mfd.write(&wrk);
	intlf = (LenFrq*) mfd.flush();
	if (ntotal < lildprm.minfreq && ftotal < lildprm.minfreq) {
	    delete[] intlf; intlf = 0;
	    return (IGNORE);
	}
	qsort((UPTR) intlf, ntotal, sizeof(LenFrq), (CMPF) lcompf);
	for (int n = 0; n < ntotal; ++n) intlf[n].frq /= ftotal;
	nfact = 1. / ftotal;
	if (lildprm.n_qtl) calc_qtl(lildprm.n_qtl);
	return (OK);
}

double Ild::rmsd(IldPrm* dp)
{
	LenFrq*	ptr = intlf;
	double	var = 0.;
	for (int i = 0; i < ntotal; ++i, ++ptr) {
	    double	len = ptr->len;
	    double	frq = ptr->frq;
	    double	d = dp->pdf_function(len) - frq / ftotal;
	    var += d * d;
	}
	return (sqrt(var / ntotal));
}

double Ild::mean(double* v)
{
	if (!ftotal) return (0);
	LenFrq*	ptr = intlf;
	double	av = 0;
	double	vr = 0;
	for (int i = 0; i < ntotal; ++i, ++ptr) {
	    av += ptr->len * ptr->frq;
	    if (v) vr += ptr->len * ptr->len * ptr->frq;
	}
	if (nfact == 0) {
	    av /= ftotal;
	    vr /= ftotal;
	}
	if (v) *v = vr - av * av;
	return (av);
}

double Ild::logmean(double* avvr)	// geometric mean
{
	if (!ftotal) return (0);
	LenFrq*	ptr = intlf;
	double	av = 0;
	double	vr = 0;
	for (int i = 0; i < ntotal; ++i, ++ptr) {
	    double	lx = log(ptr->len);
	    av += lx * ptr->frq;
	    if (avvr) vr += lx * lx * ptr->frq;
	}
	if (nfact == 0) {
	    av /= ftotal;
	    vr /= ftotal;
	}
	if (avvr) {		// retrun log mean and variance
	    avvr[0] = av;
	    avvr[1] = vr - av * av;
	}
	return (exp(av));
}

double Ild::quantile(double y)
{
	if (!ftotal || y <= 0.) return (0);
	LenFrq*	ptr = intlf;
	double	f = 1;
	double	sum = 0.;
	for (int i = 0; sum < y && i < ntotal; ++i, ++ptr) {
	    f = ptr->frq;
	    if (nfact == 0) f /= ftotal;
	    sum += f;
	}
	LenFrq*	prv = ptr - 1;
	return (ptr->len - (ptr->len - prv->len) * (sum - y) / f);
}

void Ild::normalize(double to)
{
	if (!ftotal || nfact == (to /= ftotal)) return;
	nfact = to;
	LenFrq* ptr = intlf;
	for (int i = 0; i < ntotal; ++i, ++ptr)
	    ptr->frq *= nfact;
}

double Ild::kolmo_smir(IldPrm* dp)
{
	if (!ftotal) return (0);
	LenFrq*	ptr = intlf;
	double	ks = 0.;
	double	sum = 0;
	double	f;
	for (int i = 0; i < ntotal; ++i, ++ptr) {
	    f = ptr->frq;
	    if (nfact == 0) f /= ftotal;
	    sum += f;
	    double	d = fabs(sum - dp->cumulative(ptr->len));
	    if (d > ks) ks = d;
	}
	return (sqrt(double(ntotal)) * ks);
}

double Ild::kolmo_smir(Ild* b)
{
	if (!ftotal || !b->ftotal) return (0);
	LenFrq* lfa = intlf;
	LenFrq* lfb = b->intlf;
	LenFrq* tfa = lfa + ntotal;
	LenFrq* tfb = lfb + b->ntotal;
	double  ks = 0.;
	double	asum = 0;
	double	bsum = 0;
	while (lfa < tfa || lfb < tfb) {
	    int	alen = lfa < tfa? lfa->len: INT_MAX;
	    int	blen = lfb < tfb? lfb->len: INT_MAX;
	    if (alen <= blen) {
		double	af = lfa->frq;
		if (nfact == 0) af /= ftotal;
		asum += af;
		++lfa;
	    }
	    if (alen >= blen) {
		double	bf = lfb->frq;
		if (b->nfact == 0) bf /= b->ftotal;
		bsum += bf;
		++lfb;
	    }
	    double	d = fabs(asum - bsum);
	    if (d > ks) ks = d;
	}
	return (sqrt(double(ntotal)) * ks);
}

void Ild::print_lf(const char* fn)
{
	FILE*    fd = fopen(fn,"w");
	if (!fd) fd = stdout;
	fputs("len\tfreq\n", fd);
	LenFrq*	lf = intlf;
	LenFrq*	tf = lf + ntotal;
	for ( ; lf < tf; ++lf) {
	    if (nfact == 0)
		fprintf(fd, "%d\t%d\n", lf->len, (int) lf->frq);
	    else
		fprintf(fd, "%d\t%15.7e\n", lf->len, lf->frq);
	}
	if (fd != stdout) fclose(fd);
}

/***********************************************************************************
*
*	Observed intron length distribution with x-axis in log_10 scale
*
***********************************************************************************/

int Lild::fget(FILE* fd, const char* fn)
{
	fname = fn;
	char    str[MAXL];
	reset();
	while (fgets(str, MAXL, fd)) {
	    char*       ps = str;
	    while (*ps && isspace(*ps)) ++ps;
	    if (!isdigit(*ps)) continue;
	    int len; double frq;
	    if (sscanf(ps, "%d %lf", &len, &frq) != 2) break;
	    accumulate(len, frq);
	}
	normalize(1., lildprm.psdcnt);
	ftotal = samples();
	if (ftotal < lildprm.minfreq) {
	    erase();
	    return (IGNORE);
	}
	return (OK);
}

Lild::Lild(Ild& ild, DblDbl t, DblDbl r) 
	: PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, t, r),
	  bsp(0), sid(0), ftotal(0), fname(ild.fname)
{
	LenFrq*	tlf = ild.end();
	reset();
	for (LenFrq* ilf = ild.begin(); ilf < tlf; ++ilf) 
	    accumulate(ilf->len, ilf->frq);
	normalize();
	ftotal = samples();
}

void Lild::smooth()
{
	if (!bsp) fit2bspline(lildprm.k_order, lildprm.n_coeff);
	double*	dt = begin();
	for (int i = 0; i < n_data; ++i) 
	    *dt++ = bsp->fitted(xval[i]);
}

int Lild::estimate_n_mode()
{
	if (!bsp) fit2bspline(lildprm.k_order, lildprm.n_coeff, 1);
	int	n_mode = 0;
	double	prv = 0;
	for (int i = 0; i < n_data; ++i) {
	    double	der1 = bsp->derivative(xval[i]);
	    if (prv > 0 && der1 <= 0) ++n_mode;
	    prv = der1;
	}
	return (n_mode);
}

static double cdf_func(double x, void* param)
{
	CdfParam* cp = (CdfParam*) param;
	return (cp->ip->cumulative(x) - cp->y);
}

static double dcdf_func(double x, void* param)
{
	CdfParam* cp = (CdfParam*) param;
	return (cp->ip->pdf_function(x));
}

static void fdf_func(double x, void* param, double *y, double *dy)
{
	CdfParam* cp = (CdfParam*) param;
	*y = cp->ip->cumulative(x) - cp->y;
	*dy = cp->ip->pdf_function(x);
}

static double rqrl_func(double x, void* param)
{
	double*	mlogq = (double*) param;
	return ((pow(mlogq[0], x) + mlogq[3] * pow(mlogq[2], x)
	    - (1 +  mlogq[3]) * pow(mlogq[1], x)));
}

static double drqrl_func(double x, void* param)
{
	double*	mlogq = (double*) param;
	double	lmd[3];
	for (int i = 0; i < 3; ++i) lmd[i] = pow(mlogq[i], x);
	return (log(mlogq[0]) * lmd[0] + mlogq[3] * log(mlogq[2]) * lmd[2]
	    - (1. + mlogq[3]) * log(mlogq[1]) * pow(mlogq[1], x));
}

static void frqrl_func(double x, void* param, double* y, double *dy)
{
	double*	mlogq = (double*) param;
	double	lmd[3];
	for (int i = 0; i < 3; ++i) lmd[i] = pow(mlogq[i], x);
	*y = pow(mlogq[0], x) + mlogq[3] * pow(mlogq[2], x)
	    - (1 +  mlogq[3]) * pow(mlogq[1], x);
	*dy = log(mlogq[0]) * lmd[0] + mlogq[3] * log(mlogq[2]) * lmd[2]
	    - (1. + mlogq[3]) * log(mlogq[1]) * pow(mlogq[1], x);
}

/***********************************************************************************
*
*	Compare two analytical Ilds 
*
***********************************************************************************/

/******* distance measure functions *******/

static double compspd_L2(double x, void* prm)
{
	IldPrm*	a = (IldPrm*) prm;
	double	af = a->pdf_function(x);
	return (af * af);
}

// unique area

static double compspd_UA(double x, void* prm)
{
	IldPrm**	ab = (IldPrm**) prm;
	double  af = ab[0]->pdf_function(x);
	double  bf = ab[1]->pdf_function(x);
	return min(af, bf);
}

static double compspd_CS(double x, void* prm)
{
	IldPrm**	ab = (IldPrm**) prm;
	double  af = ab[0]->pdf_function(x);
	double  bf = ab[1]->pdf_function(x);
	return (af * bf);
}

// Manhattan distance
static double compspd_MH(double x, void* prm)
{
	IldPrm**	ab = (IldPrm**) prm;
	return (fabs(ab[0]->pdf_function(x) - ab[1]->pdf_function(x)));
}

// Euclid distance
static double compspd_EC(double x, void* prm)
{
	IldPrm**	ab = (IldPrm**) prm;
	double  df = ab[0]->pdf_function(x) - ab[1]->pdf_function(x);
	return (df * df);
}

// Jensen-Shannon distance

static double compspd_JS(double x, void* prm)
{
	IldPrm**	ab = (IldPrm**) prm;
	double  af = ab[0]->pdf_function(x);
	double  bf = ab[1]->pdf_function(x);
	if (af < fepsilon || bf < fepsilon) return (0);
	double  m = (af + bf) / 2;
	return (af * log(af / m) +
		bf * log(bf / m)) / 2.;
}

// Kullback-Leibler distance

static double compspd_KL(double x, void* prm)
{
	IldPrm**	ab = (IldPrm**) prm;
	double  af = ab[0]->pdf_function(x);
	double  bf = ab[1]->pdf_function(x);
	if (af < fepsilon || bf < fepsilon) return (0);
	return ((af * log(af / bf) +
		bf * log(bf / af)) / 2.);
}

// Distance between two analytical ilds

double l2norm(IldPrm* a)
{
	gsl_set_error_handler_off();
	gsl_integration_workspace* w =
	    gsl_integration_workspace_alloc(gslprm.n_worksp);
	double  result, abserr;
	gsl_function    fnc;
	fnc.function = &compspd_L2;
	fnc.params = (void*) a;
	int	status = gsl_integration_qagiu(&fnc, a->min_x,
	    gslprm.epsabs, gslprm.epsrel, gslprm.n_worksp, w, &result, &abserr);
	gsl_integration_workspace_free(w);
	if (status) {
	    prompt("%s : ERROR %d: %s\n",
		a->sname, status, gsl_strerror(status));
	}
	return (result);
}
	
double dist_ilds(IldPrm* a, IldPrm* b, DistMethod mthd)
{
	gsl_set_error_handler_off();
	gsl_integration_workspace* w =
	    gsl_integration_workspace_alloc(gslprm.n_worksp);
	double  result, abserr;
	gsl_function    fnc;
	double  minx = min(a->min_x, b->min_x) + 1;
	IldPrm*	ab[2] = {a, b};

	fnc.params = (void*) ab;
	if (lildprm.n_qtl) {
	  switch (mthd) {
	    case EC: return compqtl_EC<IldPrm>(a, b);
	    case E2: return compqtl_EC<IldPrm>(a, b, false);
	    case MH: return compqtl_MH<IldPrm>(a, b);
	    default: fatal("Distance mode %d not supported !\n", int(mthd));
	  }
	} else {
	  switch (mthd) {
	    case UA:
	    case JA: fnc.function = &compspd_UA; break;
	    case CS: fnc.function = &compspd_CS; break;
	    case EC: fnc.function = &compspd_EC; break;
	    case E2: fnc.function = &compspd_EC; break;
	    case MH: fnc.function = &compspd_MH; break;
	    case JS: fnc.function = &compspd_JS; break;
	    case KL: fnc.function = &compspd_KL; break;
	    default: fatal("Distance mode %d not supported !\n", int(mthd));
	  }
	}
	int status = gsl_integration_qagiu(&fnc, minx,
	    gslprm.epsabs, gslprm.epsrel, gslprm.n_worksp, w, &result, &abserr);
	gsl_integration_workspace_free(w);
	if (status) {
	    prompt("%s %s : ERROR %d: %s\n",
		a->sname, b->sname, status, gsl_strerror(status));
	}
	double	aa, bb;
	switch (mthd) {
	    case UA: result = 1. - result; break;
	    case JA: result = (1. - result) / ( 1 - result / 2); break;
	    case EC: result = sqrt(result); break;
	    case CS:
		aa = l2norm(a);
		bb = l2norm(b);
		if (aa > 0 && bb > 0) result /= sqrt(aa * bb);
		else	result = 0;
		if (result < 0 || result > 1) result = 0;
		result = acos(result);
		break;
	    default: break;
	}
	return (result);
}

/***********************************************************************************
*
*	Compare two observed Ilds with x-axis in linear scalse
*
***********************************************************************************/

/******* distance measure functions *******/

// unique area
double compild_UA(Ild* a, Ild* b)
{
	LenFrq* lfa = a->intlf;
	LenFrq* lfb = b->intlf;
	LenFrq* tfa = lfa + a->ntotal;
	LenFrq* tfb = lfb + b->ntotal;
	double  abc = 0.;
	while (lfa < tfa && lfb < tfb) {
	    if (lfa->len == lfb->len) {
		abc += min(lfa->frq, lfb->frq);
		++lfa; ++lfb;
	    } else if (lfa->len < lfb->len) {
		++lfa;
	    } else {
		++lfb;
	    }
	}
	return (1. - abc);
}

// Jaccard distance
double compild_JA(Ild* a, Ild* b)
{
	LenFrq* lfa = a->intlf;
	LenFrq* lfb = b->intlf;
	LenFrq* tfa = lfa + a->ntotal;
	LenFrq* tfb = lfb + b->ntotal;
	double  abc = 0.;
	while (lfa < tfa && lfb < tfb) {
	    if (lfa->len == lfb->len) {
		abc += min(lfa->frq, lfb->frq);
		++lfa; ++lfb;
	    } else if (lfa->len < lfb->len) {
		++lfa;
	    } else {
		++lfb;
	    }
	}
	return (1. - abc) / ( 1 - abc / 2);
}

//  Euclid distance
double compild_EC(Ild* a, Ild* b, bool sqr = true)
{
	LenFrq* lfa = a->intlf;
	LenFrq* lfb = b->intlf;
	LenFrq* tfa = lfa + a->ntotal;
	LenFrq* tfb = lfb + b->ntotal;
	double  abc = 0.;
	while (lfa < tfa && lfb < tfb) {
	    if (lfa->len == lfb->len) {
		double  d = lfa->frq - lfb->frq;
		abc += d * d;
		++lfa; ++lfb;
	    } else if (lfa->len < lfb->len) {
		abc += lfa->frq * lfa->frq;
		++lfa;
	    } else {
		abc += lfb->frq * lfb->frq;
		++lfb;
	    }
	}
	while (lfa < tfa) {
	    abc += lfa->frq * lfa->frq;
	    ++lfa;
	}
	while (lfb < tfb) {
	    abc += lfb->frq * lfb->frq;
	    ++lfb;
	}
	return (sqr? sqrt(abc): abc);
}

//  Mantattan distance
double compild_MH(Ild* a, Ild* b)
{
	LenFrq* lfa = a->intlf;
	LenFrq* lfb = b->intlf;
	LenFrq* tfa = lfa + a->ntotal;
	LenFrq* tfb = lfb + b->ntotal;
	double  abc = 0.;
	while (lfa < tfa && lfb < tfb) {
	    if (lfa->len == lfb->len) {
		abc += fabs(lfa->frq - lfb->frq);
		++lfa; ++lfb;
	    } else if (lfa->len < lfb->len) {
		abc += lfa->frq;
		++lfa;
	    } else {
		abc += lfb->frq;
		++lfb;
	    }
	}
	while (lfa < tfa) {
	    abc += lfa->frq;
	    ++lfa;
	}
	while (lfb < tfb) {
	    abc += lfb->frq;
	    ++lfb;
	}
	return (abc);
}

// Cosine distance
double compild_CS(Ild* a, Ild* b)
{
	LenFrq* lfa = a->intlf;
	LenFrq* lfb = b->intlf;
	LenFrq* tfa = lfa + a->ntotal;
	LenFrq* tfb = lfb + b->ntotal;
	double  abc = 0.;
	double  aa = 0;
	double  bb = 0;
	while (lfa < tfa && lfb < tfb) {
	    if (lfa->len == lfb->len) {
		abc += lfa->frq * lfb->frq;
		aa += lfa->frq * lfa->frq;
		bb += lfb->frq * lfb->frq;
		++lfa; ++lfb;
	    } else if (lfa->len < lfb->len) {
		aa += lfa->frq * lfa->frq;
		++lfa;
	    } else{
		bb += lfb->frq * lfb->frq;
		++lfb;
	    }
	}
	while (lfa < tfa) {
	    aa += lfa->frq * lfa->frq;;
	    ++lfa;
	}
	while (lfb < tfb) {
	     bb += lfb->frq * lfb->frq;
	    ++lfb;
	}
	if (aa > 0 && bb > 0) abc /= sqrt(aa * bb);
	else	abc = 0;
	if (abc < 0 || abc > 1) abc = 0;
	return (acos(abc));
}

// Jensen-Shannon distance

double compild_JS(Ild* a, Ild* b)
{
	LenFrq* lfa = a->intlf;
	LenFrq* lfb = b->intlf;
	LenFrq* tfa = lfa + a->ntotal;
	LenFrq* tfb = lfb + b->ntotal;
	double  js = 0.;
	while (lfa < tfa && lfb < tfb) {
	    if (lfa->len == lfb->len) {
		double  m = (lfa->frq + lfb->frq) / 2;
		double  d = (lfa->frq * log(lfa->frq / m) +
			     lfb->frq * log(lfb->frq / m)) / 2.;
		js += d;
		++lfa; ++lfb;
	    } else if (lfa->len < lfb->len) {
		double  m = lfa->frq / 2;
		double  d = lfa->frq * log(lfa->frq / m) / 2;
		js += d;
		++lfa;
	    } else {
		double  m = lfb->frq / 2;
		double  d = lfb->frq * log(lfb->frq / m) / 2;
		js += d;
		++lfb;
	    }
	}
	while (lfa < tfa) {
	    double  m = lfa->frq / 2;
	    double  d = lfa->frq * log(lfa->frq / m) / 2;
	    js += d;
	    ++lfa;
	}
	while (lfb < tfb) {
	    double  m = lfb->frq / 2;
	    double  d = lfb->frq * log(lfb->frq / m) / 2;
	    js += d;
	    ++lfb;
	}
	return (js);
}

// Kullback-Leibler distance
double compild_KL(Ild* a, Ild* b)
{
	LenFrq* lfa = a->intlf;
	LenFrq* lfb = b->intlf;
	LenFrq* tfa = lfa + a->ntotal;
	LenFrq* tfb = lfb + b->ntotal;
	double  kl = 0.;
	while (lfa < tfa && lfb < tfb) {
	    if (lfa->len == lfb->len) {
		double  d = (lfa->frq * log(lfa->frq / lfb->frq) +
			     lfb->frq * log(lfb->frq / lfa->frq)) / 2.;
		kl += d;
		++lfa; ++lfb;
	    } else if (lfa->len < lfb->len) ++lfa;
	    else	++lfb;
	}
	return (kl);
}

double dist_ilds(Ild* a, Ild* b, DistMethod mthd)
{
	if (lildprm.n_qtl) {
	  switch (mthd) {
	    case EC: return compqtl_EC<Ild>(a, b);
	    case E2: return compqtl_EC<Ild>(a, b, false);
	    case MH: return compqtl_MH<Ild>(a, b);
	    default: fatal("Distance mode %d not supported !\n", int(mthd));
	  }
	} else {
// Distance between two observed ilds with x-axis in linear scale
	  switch (mthd) {
	    case CS: return compild_CS(a, b);
	    case EC: return compild_EC(a, b);
	    case E2: return compild_EC(a, b, false);
	    case JS: return compild_JS(a, b);
	    case JA: return compild_JA(a, b);
	    case KL: return compild_KL(a, b);
	    case MH: return compild_MH(a, b);
	    case UA: return compild_UA(a, b);
	    default: fatal("Distance mode %d not supported !\n", int(mthd));
	  }
	}
	return (DistIldErr);
}

/***********************************************************************************
*
*	Compare two observed Ilds with x-axis in linear scalse
*
***********************************************************************************/

/******* distance measure functions *******/

// Unique area distance
double complild_UA(Lild* a, Lild* b)
{
	double* lfa = a->begin();
	double* lfb = b->begin();
	double  abc = 0.;
	for (int n = 0; n < a->size(); ++n)
	    abc += min(lfa[n], lfb[n]);
	return (1. - abc);
}

// Jaccard distance
double complild_JA(Lild* a, Lild* b)
{
	double* lfa = a->begin();
	double* lfb = b->begin();
	double  abc = 0.;
	for (int n = 0; n < a->size(); ++n)
	    abc += min(lfa[n], lfb[n]);
	return (1. - abc) / ( 1 - abc / 2);
}

//  Euclid distance
double complild_EC(Lild* a, Lild* b, bool sqr = true)
{
	double* lfa = a->begin();
	double* lfb = b->begin();
	double  abc = 0.;
	for (int n = 0; n < a->size(); ++n) {
	    double	ns = a->step(n);
	    double      d = lfa[n] - lfb[n];
	    abc += d * d / ns;
	}
	return (sqr? sqrt(abc): abc);
}

//  Mantattan distance
double complild_MH(Lild* a, Lild* b)
{
	double* lfa = a->begin();
	double* lfb = b->begin();
	double  abc = 0.;
	for (int n = 0; n < a->size(); ++n)
	    abc += fabs(lfa[n] - lfb[n]);
	return (abc);
}

// Cosine distance
double complild_CS(Lild* a, Lild* b)
{
	double* lfa = a->begin();
	double* lfb = b->begin();
	double  abc = 0.;
	double  aa = 0;
	double  bb = 0;
	for (int n = 0; n < a->size(); ++n) {
	    double	ns = a->step(n);
	    abc += lfa[n] * lfb[n] / ns;
	    aa += lfa[n] * lfa[n] / ns;
	    bb += lfb[n] * lfb[n] / ns;
	}
	if (aa > 0 && bb > 0) abc /= sqrt(aa * bb);
	else	abc = 0;
	if (abc < 0 || abc > 1) abc = 0;
	return (acos(abc));
}

// Jensen-Shannon distance
double complild_JS(Lild* a, Lild* b)
{
	double* lfa = a->begin();
	double* lfb = b->begin();
	double  js = 0.;
	for (int n = 0; n < a->size(); ++n) {
	    double	ns = a->step(n);
	    double	la = lfa[n] / ns;
	    double	lb = lfb[n] / ns;
	    double      m = (la + lb) / 2;
	    if (la < fepsilon || lb < fepsilon) continue;
	    double      d = (la * log(la / m) +
		     lb * log(lb / m)) / 2.;
	    js += d * ns;
	}
	return (js);
}

// Kullback-Leibler distance
double complild_KL(Lild* a, Lild* b)
{
	double* lfa = a->begin();
	double* lfb = b->begin();
	double  kl = 0.;
	for (int n = 0; n < a->size(); ++n) {
	    double      ns = a->step(n);
	    double	la = lfa[n] / ns;
	    double	lb = lfb[n] / ns;
	    if (la < fepsilon || lb < fepsilon) continue;
	    double      d = (la * log(la / lb) +
		     lb * log(lb / la)) / 2.;
	    kl += d * ns;
	}
	return (kl);
}

// Distance between two observed ilds with x-axis in log_10 scale

double dist_ilds(Lild* a, Lild* b, DistMethod mthd)
{
	switch (mthd) {
	    case CS: return complild_CS(a, b);
	    case EC: return complild_EC(a, b);
	    case E2: return complild_EC(a, b, false);
	    case JS: return complild_JS(a, b);
	    case JA: return complild_JA(a, b);
	    case KL: return complild_KL(a, b);
	    case MH: return complild_MH(a, b);
	    case UA: return complild_UA(a, b);
	    default: fatal("Distance mode %d not supported !\n", int(mthd));
	}
	return (DistIldErr);
}

/***********************************************************************************
*
*	Graphical display of ilds
*
***********************************************************************************/

void GnuPlotLild::initialize(IldPrm* dprm, bool ildents)
{
	n_clmns = n_lild + n_fitf + 2;
	if (ildents) {
	    for (int n = 0; n < n_fitf; ++n, ++dprm)
		n_clmns += dprm->n_modes;
	}
	n_rows = size();
	normal_factor = invtransform? (invtransform(width) - 1.) / invtransform(width / 2): 1.;
	data = new double*[n_clmns];
	*data = new double[n_clmns * (n_rows + 2)] + 1;
	for (int i = 1; i < n_clmns; ++i)
	    data[i] = data[i - 1] + (n_rows + 2);
	if (invtransform) {
	    int 	px = xval[-1];
	    data[1][-1] = px;
	    data[0][-1] = nlitransform(double(px));
	    for (int n = 0; n <= n_rows; ++n) {
		int	ix = xval[n];
		if ((ix - px) == 1) {
		    data[1][n] = ix;
		    data[0][n] = nlitransform(data[1][n]);
		} else {
		    ++px;
		    data[0][n] = (nlitransform(double(px)) + nlitransform(double(ix))) / 2;
		    data[1][n] = invtransform(data[0][n]);
		}
		px = ix;
	    }
	} else {
	    for (int n = -1; n <= n_rows; ++n)
		data[0][n] = data[1][n] = xval[n];
	}
	clmn = 2;
}

void GnuPlotLild::add(IldPrm* dprm, bool ildents)
{
	int	c = clmn;
	for (int n = 0; n <= n_rows; ++n) {
	    c = clmn;
	    double	dx = data[1][n];
	    if (ildoutmode == IldOutMode::CDF) {
		data[c++][n] = dprm->cumulative(dx);
		if (!ildents) continue;
		for (int m = 1; m <= dprm->n_modes; ++m)
		    data[c++][n] = dprm->cumulative(dx, m);
	    } else if (ildoutmode == IldOutMode::PDF) {
		double	y = normal_factor * dprm->pdf_function(dx);
		if (nlitransform) y *= dx;
		data[c++][n] = y;
		if (!ildents) continue;
		for (int m = 1; m <= dprm->n_modes; ++m) {
		    double y = normal_factor * dprm->pdf_function(dx, m);
		    if (nlitransform) y *= dx;
		    data[c++][n] = y;
		}
	    } else if (ildoutmode == IldOutMode::Penalty) {
		data[c++][n] = log10(dprm->pdf_function(dx));
	    }
	}
	sname.push(dprm->sname);
	clmn = c;
}

void GnuPlotLild::add(Lild* lild)
{
const	char*	sl = strrchr(lild->fname, '/');
	if (!sl) sl = lild->fname;
	char	str[MAXL];
	char*	ps = str;
	while (*sl && *sl != '.') *ps++ = *sl++;
	*ps = '\0';
	sname.push(str);
	if (ildoutmode == IldOutMode::CDF) {
	    double  cdf = 0;
	    int	n = -1;
	    for (double* lf = lild->begin(); lf < lild->end(); ++lf, ++n)
		data[clmn][n] = cdf += *lf; 
	} else
	    vcopy(data[clmn] - 1, lild->begin() - 1, n_rows + 2);
	++clmn;
}

void GnuPlotLild::add(Ild* ild)
{
const	char*	sl = strrchr(ild->fname, '/');
	if (!sl) sl = ild->fname;
	char	str[MAXL];
	char*	ps = str;
	while (*sl && *sl != '.') *ps++ = *sl++;
	*ps = '\0';
	sname.push(str);
	double	cdf = 0;
	for (LenFrq* lf = ild->begin(); lf < ild->end(); ++lf) {
	    int	n = (lf->len - int(llmt)) / ndiv;
	    if (n < 0) n = -1;
	    if (n > ulmt) n = int(ulmt);
	    data[clmn][n] = (ildoutmode == IldOutMode::CDF)? 
		cdf += lf->frq: lf->frq;
	}
	++clmn;
}

GnuPlotLild::GnuPlotLild(Lild* lild, IldPrm* dprm)
	: PutIntoBins(*lild), n_lild(lild? 1: 0), n_fitf(dprm? 1: 0), data(0)
{
	initialize();
	if (dprm) add(dprm);
	if (lild) add(lild);
}

GnuPlotLild::GnuPlotLild(Ild* ild, IldPrm* dprm)
	: PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, 0, 0, false),
	  n_lild(1), n_fitf(dprm? 1: 0), data(0), normal_factor(1.)
{
	initialize();
	if (dprm) add(dprm);
	add(ild);
}

GnuPlotLild::GnuPlotLild(IldPrm* dprm, int np, Lild** lild, int ni, bool ildents)
	: PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, log10, exp10, false),
	  n_lild(ni), n_fitf(np), data(0)
{
	initialize(dprm, ildents);
	if (dprm) {
	    for (int i = 0; i < np; ++i)
		add(dprm + i, ildents);
	}
	if (lild) {
	    for (int i = 0; i < ni; ++i)
		add(lild[i]);
	}
}

GnuPlotLild::GnuPlotLild(IldPrm* dprm, int np, Ild** ilds, int ni, bool ildents)
	: PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, 0, 0, false),
	  n_lild(ni), n_fitf(np), data(0)
{
	initialize(dprm, ildents);
	if (dprm) {
	    for (int i = 0; i < np; ++i)
		add(dprm + i, ildents);
	}
	for (int i = 0; i < ni; ++i)
		add(ilds[i]);
}

GnuPlotLild::GnuPlotLild(Lild** lilds, int num)
	: PutIntoBins(*lilds[0]), n_lild(num), n_fitf(0), data(0)
{
	initialize();
	for (int i = 0; i < num; ++i)
	    vcopy(data[clmn++] - 1, lilds[i]->begin() - 1, n_rows);
}

GnuPlotLild::GnuPlotLild(Ild** ilds, int num)
	: PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, 0, 0, false),
	  n_lild(num), n_fitf(0), data(0)
{
	initialize();
	for (int i = 0; i < num; ++i, ++clmn) {
	    for (LenFrq* lf = ilds[i]->begin(); lf < ilds[i]->end(); ++lf) {
		int	n = (lf->len - int(llmt)) / ndiv;
		if (n < 0) n = -1;
		if (n > ulmt) n = int(ulmt);
		data[clmn][n] = lf->frq;
	    }
	}
}

GnuPlotLild::GnuPlotLild(IldPrm** dprms, int num, DblDbl t, DblDbl r)
	: PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, t, r, false),
	  n_lild(0), n_fitf(num), data(0)
{
	initialize();
	for (int i = 0; i < num; ++i)
	    add(dprms[i]);
}

void GnuPlotLild::plot(const char* oform, Strlist* iname, StrHash<int>* phyl_code, IldPrm* dprm)
{
static	const	char	pdffmt[] = "lines lw %d lc rgb '%s'";
static	const	char	smlfmt[] = "%s ps 2 lc rgb '%s'";
const	char*	smllt = n_fitf? "points": "linespoints";
	FILE*	fd = fopen(gpdata, "w");
	if (!fd) fatal("Can't write to %s !\n", gpdata);
const	char**	header = new const char*[n_clmns];
	fputs("loglen\tlen", fd);
	int	c = 0;
	for (int i = 0; i < n_fitf; ++i) {
	    fprintf(fd, "\t%s", header[c++] = sname[i]);
	    if (!dprm) continue;
	    for (int m = 0; m < dprm->n_modes; ++m)
		fprintf(fd, "\t%s", header[c++] = sname[i]);
	}
	if (iname) {
	    for (int i = 0; i < n_lild; ++i)
		fprintf(fd, "\t%s", header[c++] = (*iname)[i]);
	} else {
	    for (int i = 0; i < n_lild; ++i) 
		fprintf(fd, "\t%s", header[c++] = sname[n_fitf + i]);
	}
	fputc('\n', fd);
	for (int n = 0; n < n_rows; ++n) {
	    fprintf(fd, "%7.5f\t%7.1f", data[0][n], data[1][n]); 
	    for (int c = 2; c < clmn; ++c)
		fprintf(fd, "\t%7.5f", data[c][n]);
	    fputc('\n', fd);
	}
	fclose(fd);

	fd = fopen(gpcmd, "w");
	if (!fd) fatal("Can't write to %s !\n", gpcmd);
	if (oform) {
	    const char*	trm = strrchr(oform, '.');
	    if (trm) ++trm;
	    else	trm = def_form;
	    fprintf(fd, "set terminal %s\n", trm);
	    fprintf(fd, "set output '%s'\n", oform);
	}
	c = 0;
	fprintf(fd, "set xrange [%.2f:%.2f]\n", lildprm.llmt, lildprm.ulmt);
	fprintf(fd, "set key autotitle columnheader\n");
	fprintf(fd, "plot '%s' using 1:3 ", gpdata);	// first
	fprintf(fd, "title '%s' with ", header[c++]);
	for (int i = 0; i < n_fitf; ++i) {
	    if (i == 0) fprintf(fd, pdffmt, 4, "black");
	    else {
		int	cc = c++;
		fprintf(fd, ",\\\n'' using 1:($%d) ", cc + 3);
		fprintf(fd, "title '%s' with ", header[cc]);
		if (phyl_code) {
		    KVpair<INT, int>* kv = phyl_code->find(header[cc]);
		    if (kv) cc = kv->val + 2;
		}
		fprintf(fd, pdffmt, 4, colors[cc % n_colors + 1]);
	    }
	    if (dprm) {
		for (int m = 0; m < dprm->n_modes; ++m) {
		    fprintf(fd, ",\\\n'' using 1:($%d) ", c + 3);
		    fprintf(fd, "title '%s' with ", header[c++]);
		    fprintf(fd, pdffmt, 2, colors[c % n_colors + 1]);
		}
		++dprm;
	    }
	}
	if (!n_fitf)    fprintf(fd, smlfmt, smllt, "black");
	for (int i = mark_color; c < clmn - 2; ++i, ++c) {
	    fprintf(fd, ",\\\n'' using 1:($%d) ", c + 3);
	    fprintf(fd, "title '%s' with ", header[c]);
	    int	cc = i;
	    if (phyl_code) {
		KVpair<INT, int>* kv = phyl_code->find(header[c]);
		if (kv) cc = kv->val + 2;
	    }
	    fprintf(fd, smlfmt, smllt, colors[cc % n_colors + 1]);
	}
	fputs("\npause -1 'push RETURN key to quit !'\n", fd);
	fclose(fd);
	char	cmd[MAXL];
	snprintf(cmd, MAXL, "gnuplot %s", gpcmd);
const	int	srv = system(cmd);
	if (srv) fatal("%s failed with code = %d !\n", cmd, srv);
	delete[] header;
}

static const char* terminal_exts[] = 
{"emf", "ps", "eps", "latex", "svg", "jpeg", "png", "gif", "tab", 0};
static const char* terminal_types[] = 
{"", "postscript", "postscript eps", "epslatex", "", "", "", "", "table", 0};

const char* gnuplot_terminal(const char* ext)
{
const	char*	dot = strrchr(ext, '.');
	if (dot) ext = dot + 1;
	char     str[MAXL];
	strcpy(str, ext);
	for (char* ps = str; *ps; ++ps)
	    *ps = tolower(*ps);
	for (int i = 0; terminal_exts[i]; ++i)
	    if (!strcmp(str, terminal_exts[i]))
		return (terminal_types[i]);
	return (0);
}

