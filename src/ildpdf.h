/*****************************************************************************
*
*	<<< ildpdf.h >>>
*
*	header to probabilistic distribution functions
*
*	Requires GSL library 
*	https://www.gnu.org/software/gsl/
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

#ifndef	_ILDPDF_H_
#define	_ILDPDF_H_

#include "stdtype.h"
#include "ild_data.h"
#include "dist2.h"
#include <math.h>
#include <vector>
#include <gsl/gsl_version.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>

enum	StatDist {FRECHET, LOGNORMAL, GEOMETRIC, GAMMA, WEIBULL};
enum	GplotMode {INITIAL, CONT, LAST};
enum	IldOutMode {PDF, CDF, Penalty};

typedef	double	(*DblDbl)(double x);

struct	dLfp {int len; double frq;};
struct	LildPrm {int minfreq, ndiv, minx, maxx, n_qtl, k_order, n_coeff; 
	double llmt, ulmt, psdcnt, max_qtl;};
struct	GslPrm {int n_worksp, opt_max_iter, rsl_max_iter; 
	double epsabs, epsrel, mn_step, mn_conv, mn_grad;};
struct	FitIndex {double val[3];};
typedef	std::vector<FitIndex>::iterator	FiiItr;

	typedef struct drand48_data     Drand48_data;

static	const	int	NoDistParam[] = {4, 3, 3, 4, 4, 4};
static	const	double	sqrt2pi = 2.50662827463;// sqrt(2 * pi)
static	const	double	ln_10 = 2.30258509299;	// ln(10)
static	const	double	ln_2 = 0.69314718056;	// ln(2)
						// erf^(-1)(0.5) * sqrt(2)
static	const	double	Eu_gamma = 0.5772156649;
static	const	double	l14coef = 0.47693627 * 1.41421356;
static	const	int	max_modes = 3;
static	const	int	paramsize = 12;
static	const	int	n_skip_term = 7;
static	const	char	gpcmd[] = "ild2lild.gp";
static	const	char	gpdata[] = "ild2lild.dat";
static	const	int	def_cdf_tab_size = 10000;
static	const	double	OutOfRange = DBL_MAX / 2;
static	const	double	DistIldErr = -OutOfRange;
static	const	double	SmallProb = 1.e-20;
static	const	int	min_hash = 1000;
static	const	double	minl_mu = 0.1;

extern	LildPrm	lildprm;
extern	GslPrm	gslprm;
extern	StatDist	defpdf;
extern	IldOutMode	ildoutmode;
extern	const	char*	gnuplot_terminal(const char* ext);

class Verbout {
	int	pro_phase = 1;
	int	iter_no = 0;
public:
	Verbout() {}
	void	set_phase(int phs) {pro_phase = phs; iter_no = 0;}
	void	printout(int n, const double* p, double fx, const char* spc = 0);
};

class SpecList : public StrHash<int> {
	int	n_spec = 0;
public:
	SpecList(int ac, const char** av);
	~SpecList() {}
	int	size() const {return (n_spec);}
};

// dparam[] = {m1, t1, k1, a1, m2, t2, k2, a2, m3, t3, k3, a3}

class	Ild;

class	IldPrm {
	StatDist	ildpdf;
public:
	int	n_modes = 0;
	int	m_size = 0;
	bool	sibs = false;
	char*	genspc = 0;
	int	n_param = 0;
	double	dparam[paramsize];
	double	min_aic = DBL_MAX;
	int	sid = 0;
	INT	n_sample = 0;
	float*	k_th_qtil = 0;
	double*	cdf_table = 0;
	int	min_x = 0;
	int	max_x = 0;
	StatDist	statdist() {return (ildpdf);}
	double	pdf_function(double len, int ildent = 0);
	double	cumulative(double len, int ildent = 0);
	double	mean();
	double	quantile(double y, double iv = 0);
	double	negentropy_10();
	bool	proper();
	bool	empty() const {return (!n_param);}
	char*	strget(char* str, char** terms = 0);
template <typename file_t>
	int	fget(file_t& fd, const char* fn = 0) {
	    char    str[MAXL];
	    char*   ps = 0;
	    while ((ps = fgets(str, MAXL, fd)))
		if (*ps != '#') break;
	    if (!ps) return (EOF);
	    return (strget(str)? OK: IGNORE);
	}
	void	clear() {
	    if (!sibs) {
		delete[] genspc; delete[] k_th_qtil; delete[] cdf_table;
	    }
	    n_modes = n_param = sid = n_sample = min_x = max_x = 0;
	    sibs = false;
	    genspc = 0; k_th_qtil = 0; cdf_table = 0;
	}
	void	complete();
	int	single_mode(double q123_4[]);
	void	copy_from(const IldPrm& src) {
	    char*	sn = genspc;
	    *this = src;
	    genspc = sn;
	}
	bool	reduce_dim() {
	    if (n_param > m_size) {
		n_param -= m_size;
		--n_modes;
		return (true);
	    }
	    return (false);
	}
	int	get_IldPrm(char* str, const INT& db_min_samples = 0, 
		double* tparam = 0);
	bool	good() const {return (n_sample >= min_samples);}
	bool	convert(StatDist topdf);
	int	printprm(FILE* fd);
	void	calc_qtl(int nq);
	void	calc_cdf();
	char*	entry() const {return (genspc);}
	double	aic(double fx) {return(2 * (n_param + fx));}
	double	bic(double fx) {
	    return ((n_param * log((double) n_sample)) + 2 * fx);
	}
	IldPrm() {
	    m_size = NoDistParam[ildpdf = defpdf];
	}
	IldPrm(const IldPrm& src) {*this = src; sibs = true;}
	IldPrm(int np, const char** argv, StatDist pdf = FRECHET);
	IldPrm(const char* ip_stat, Ild* inst, const char* choose = 0);
	IldPrm(const char* ip_stat, const char* gs);
	IldPrm(Ild* inst, StatDist pdf = FRECHET);
	~IldPrm(){
	    if (!sibs) {
		delete[] genspc; genspc = 0;
		delete[] k_th_qtil; k_th_qtil = 0;
		if (cdf_table) {delete[] (cdf_table + min_x); cdf_table = 0;}
	    }
	}
};

class Ild {
	double	nfact = 0;
public:
	int	sid = 0;
	bool	sibs = false;
	int	ntotal = 0;
	double	ftotal = 0;
const	char*	entry = 0;
	dLfp*	intlf = 0;
	float*	k_th_qtil = 0;
	bool	empty() const {return ftotal == 0;}
	Ild() {}
	Ild(const char* fn, int id = 0) : sid(id) {
	    FILE*	fd = 0;
	    fget(fd, fn);
	}
template <typename file_t>
	Ild(file_t fd, const char* fn, int id) : sid(id) {
	    fget(fd, fn);
	}
	Ild(Ild& src) {
	    *this = src;
	    sibs = true;
	}
	Ild(Ild* src, Drand48_data* drand_buff);
	Ild(IldPrm* ildprm, const INT& nsample, 
		Drand48_data* drand_buff, const char* snm = 0);
	~Ild() {
	    if (!sibs) {delete[] entry; delete[] intlf; delete[] k_th_qtil;}
	}
	dLfp*	begin() {return intlf;}
	dLfp*	end() {return intlf + ntotal;}
template <typename file_t>
	int	fget(file_t& fd, const char* fn = 0);
	void	print_lf(const char* fn = 0);
	double	rmsd(IldPrm* dp);
	double	mean(double* v = 0);
	double	logmean(double* avvr = 0);	// exp(mean(log(x))
	double	quantile(double y);
	double	median() {return quantile(0.5);}
	double	kolmo_smir(Ild* ild2);
	double	kolmo_smir(IldPrm* dp);
	void	normalize(double to = 1.);
	int	fitild(IldPrm* dp, double* pfx = 0);
	void	calc_qtl(int nq);
	INT	n_samples() {return (INT(ftotal));}
	int	minx() const {return (intlf->len);}
	int	maxx() const {return (intlf[ntotal - 1].len);}
	bool	good() const {return (ftotal >= min_samples);}
};

// compatibile with CalcServer::VarLoader

template <typename file_t>
int Ild::fget(file_t& dmyfd, const char* fn)
{
	if (dmyfd) {
	    fclose(dmyfd);
	    dmyfd = 0;
	}
	Ild_data	ild_dt(fn);
	char	str[LINE_MAX];
	strcpy(str, fn);
	char*	sl = strrchr(str, '/');
	char*	dot = strchr(str, '.');
	sl = sl? sl + 1: str;
	if (dot) *dot = '\0';
	entry = strrealloc(0, sl);
	ntotal = ild_dt.size();
	intlf = new dLfp[ntotal + 1];
	dLfp*	wlf = intlf;
	sLfp*	ilf = ild_dt.lenfrq;
	sLfp*	elf = ilf + ntotal;
	for (ntotal = 0; ilf < elf; ++ilf) {
	    if (ilf->len > lildprm.minx && ilf->len <= lildprm.maxx) {
		wlf->len = ilf->len;
		(wlf++)->frq = ilf->frq;
		++ntotal;
		ftotal += ilf->frq;
	    }
	}
	if (ntotal < lildprm.minfreq && ftotal < lildprm.minfreq) {
	    delete[] intlf; intlf = 0;
	    return (IGNORE);
	}
	if (lildprm.n_qtl) calc_qtl(lildprm.n_qtl);
	return (OK);
}

// fit data to B-spline

class Bspline {
protected:
	int	k_order = 0;
	int	n_breaks = 0;
	int	n_coeffs = 0;
	gsl_bspline_workspace*	bsws;
#if GSL_MAJOR_VERSION == 1
	gsl_bspline_deriv_workspace*	bsdws;
#endif
	gsl_vector*	Bk = 0;
	gsl_matrix*	dB = 0;
	gsl_vector*	coeff = 0;
	gsl_matrix*	cov = 0;
public:
	Bspline(PutIntoBins& bins, double* dt, 
		int ko = 4, int nc = 12, int n_th_d = 0);
	~Bspline() {
	    gsl_bspline_free(bsws);
#if GSL_MAJOR_VERSION == 1
	    if (bsdws) gsl_bspline_deriv_free(bsdws);
#endif
	    gsl_vector_free(Bk);
	    gsl_matrix_free(dB);
	    gsl_vector_free(coeff);
	    gsl_matrix_free(cov);
	}
	double	fitted(double x);
	double	derivative(double x);
};

//	distribute into bins

class Lild : public PutIntoBins {
friend	class	GnuPlotLild;
	Bspline*	bsp = 0;
public:
	int	sid = 0;
	bool	sibs = false;
	double	ftotal = 0;
const	char*	entry = 0;
	Lild(DblDbl t = log10, DblDbl r = exp10)
	    : PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, t, r)
	{}
	Lild(Lild& src) : PutIntoBins(src) { *this = src; }
	Lild(Ild& ild, DblDbl t = log10, DblDbl r = exp10);
	Lild(const char* fn, const int& sd, DblDbl t = log10, DblDbl r = exp10);
template <typename file_t>
	Lild(file_t fd, const char* fn, DblDbl t = log10, DblDbl r = exp10)
	    :  PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, t, r),
		bsp(0), sid(0), ftotal(0) {
	    fget(fd, fn);
	}
	~Lild() {delete bsp; delete[] entry;}
template <typename file_t>
	int	fget(file_t& fd, const char* fn = 0);
	void	fit2bspline(int ko, int n_coef, int n_th_d = 0) {
	    bsp = new Bspline(*this, begin(), n_coef, n_th_d);
	}
	void	smooth();
	int	estimate_n_mode();
	bool	good() const {return (true);}
};

// compatibile with CalcServer::VarLoader

template <typename file_t>
int Lild::fget(file_t& dmyfd, const char* fn)
{
	if (dmyfd) {
	    fclose(dmyfd);
	    dmyfd = 0;
	}
	Ild_data	ild_dt(fn);
	reset();
	char	str[LINE_MAX];
	strcpy(str, fn);
	char*	sl = strrchr(str, '/');
	char*	dot = strchr(str, '.');
	sl = sl? sl + 1: str;
	if (dot) *dot = '\0';
	entry = strrealloc(0, sl);
	int	ntotal = ild_dt.size();
	sLfp*	ilf = ild_dt.lenfrq;
	sLfp*	elf = ilf + ntotal;
	for ( ; ilf < elf; ++ilf)
	    accumulate(ilf->len, (double) ilf->frq);
	normalize(1., lildprm.psdcnt);
	ftotal = samples();
	if (ftotal < lildprm.minfreq) {
	    ftotal = 0;
	    erase();
	    return (IGNORE);
	}
	return (OK);
}

class GnuPlotLild : public PutIntoBins {
	int	n_lild = 0;
	int	n_fitf = 0;
	int	n_clmns = 0;
	int	n_rows = 0;
	int	clmn = 0;
	Strlist	sname;
	double**	data = 0;
	double	normal_factor = 1.;
	void	initialize(IldPrm* dprm = 0, bool ildents = false); 
	void	add(IldPrm* dprm, bool ildents = false);
	void	add(Ild* lild);
	void	add(Lild* lild);
public:
	void	plot(const char* oform = 0, Strlist* iname = 0, 
		StrHash<int>* phyl_code = 0, IldPrm* dprm = 0);
	GnuPlotLild(int c = 1, int f = 0, DblDbl t = log10, DblDbl r = exp10) 
	    : PutIntoBins(lildprm.ndiv, lildprm.llmt, lildprm.ulmt, t, r, false),
		n_lild(c), n_fitf(f), data(0), normal_factor(1.)
		{initialize();}
	GnuPlotLild(Lild* lild, IldPrm* dprm = 0);
	GnuPlotLild(Ild* ild, IldPrm* dprm = 0);
	GnuPlotLild(IldPrm* dprm, int nprm, Lild** lild, int nild = 1, bool compo = false);
	GnuPlotLild(IldPrm* dprm, int nprm, Ild** ilds, int nild = 1, bool compo = false);
	GnuPlotLild(Lild** lilds, int num);
	GnuPlotLild(Ild** ilds, int num);
	GnuPlotLild(IldPrm** lilds, int num, DblDbl t = log10, DblDbl r = exp10);
	~GnuPlotLild() {
	    if (data) {
		if (*data) delete[] (*data - 1);
		delete[] data;
	    }
	} 
};

inline	double frechet(double z, double w, double t, double k) {
	return (k / t * z * w * exp(-w));
}

// dlogF/dx

inline	double lfrechet_x(double z, double w, double t, double k) {
	return (z * (k * (w - 1) - 1) / t);
}

// dlogF/dm = -dlogF/dx

inline	double lfrechet_m(double z, double w, double t, double k) {
	return (z * (1 + k * (1 - w)) / t);
}

// dlogF/dt

inline	double lfrechet_t(double z, double w, double t, double k) {
	return (k * ( 1 - w) / t);
}

// dlogF/dk

inline	double lfrechet_k(double z, double w, double t, double k) {
	return (1/k + log(z) * ( 1 - w));
}

inline double lognormal(double x, double z, double m, double s) {
//	double	z = (log(x) - m) / s;
	return (exp(-z * z / 2) / (x * s * sqrt2pi));
}

inline double llognormal_m(double x, double z, double m, double s) {
	return (z / s);
}

inline double llognormal_s(double x, double z, double m, double s) {
	return ((z * z - 1) / s);
}

inline double frechet_median(double* p)
{
	return (p[0] + p[1] / pow(0.69314718056, 1 / p[2]));
}				//  ln(2)

inline double frechet_mean(double* p)
{
	return (p[2] > 1? (p[0] + p[1] * gsl_sf_gamma(1 - 1 / p[2])): frechet_median(p));
}

inline double frechet_mode(double* p)
{
        return (p[0] + p[1] * pow(p[2] / (1 + p[2]), 1 / p[2]));
}

inline double frechet_entropy(double* p)
{
	return (1 - Eu_gamma / p[2] + Eu_gamma + log(p[1] / p[2]));
}

inline double freche3_4quantile(double* p)
{
	return (p[0] + p[1] / pow(0.28768207245, 1. / p[2]));
}				//  ln(4/3)

inline double frechet_quantile(double* p, double q)
{
	return (p[0] + (q > 0? p[1] / pow(-log(q), 1. / p[2]): 0));
}

inline double lognormal_mean(double* p)
{
	return (exp(p[0] + p[1] * p[1] / 2));
}

inline double lognormal_median(double* p)
{
	return (exp(p[0]));
}

inline double lognormal_mode(double* p)
{
	return (exp(p[0] - p[1] * p[1]));
}

inline double lognormal_entropy(double* p)
{
	return (log(p[1] * sqrt2pi) + p[0] + 0.5);
}

inline double lognormal1_4quantile(double* p)
{
	return (exp(p[0] - p[1] * l14coef));
}

inline double lognormal3_4quantile(double* p)
{
	return (exp(p[0] + p[1] * l14coef));
}

extern double lognormal_quantile(double* p);

inline double geometric(double x, double q, double d)
{
	return ((x > d + 1)? q * pow(1. - q, x - d - 1): 0);
}

inline double lgeometric_q(double x, double q, double d)
{
	return ((x > d + 1)? 1 / q - (x - d - 1) / (1 - q): 0);
}

inline double lgeometric_d(double x, double q, double d)
{
        return (-log(1 - q));
}

inline double geometric_mean(double* p)
{
	return (p[1] + 1. / p[0]);
}

inline double geometric_median(double* p)
{
	return (p[1] - log(2.) / log(1 - p[0]));
}

inline double geometric_mode(double* p)
{
	return (p[1]);
}

inline double geometric_quantile(double* p, double y)
{
	return (p[1] + log(1 - y) / log(1 - p[0]));
}

inline	double gamma(double z, double t, double k) {
	return (pow(z, k - 1) * exp(-z) / t / gsl_sf_gamma(k));
}

// dlogF/dx

inline	double lgamma_x(double z, double t, double k) {
	return (((k - 1) / z - 1) / t);
}

// dlogF/dm = -dlogF/dx

inline	double lgamma_m(double z, double t, double k) {
	return (((1 - k) / z + 1) / t);
}

// dlogF/dt

inline	double lgamma_t(double z, double t, double k) {
	return ((z - k) / t);
}

// dlogF/dk

inline	double lgamma_k(double z, double t, double k) {
	return (log(z) - gsl_sf_psi(k));
}

inline double gamma_mean(double* p)
{
	return (p[1] * p[2] + p[0]);
}

inline double gamma_mode(double* p)
{
	return (p[2] > 1)? ((p[2] - 1) * p[1] + p[0]): 0;
}

inline	double weibull(double z, double w, double t, double k) {
	return (z > 0? k / t * w * exp(-w) / z: 0);
}

// dlogF/dx

inline	double lweibull_x(double z, double w, double t, double k) {
	return ((k * (1 - w) - 1) / t / z);
}

// dlogF/dm = -dlogF/dx

inline	double lweibull_m(double z, double w, double t, double k) {
	return ((k * (w - 1) + 1) / t / z);
}

// dlogF/dt

inline	double lweibull_t(double z, double w, double t, double k) {
	return (-k * ( 1 + w) / t);
}

// dlogF/dk

inline	double lweibull_k(double z, double w, double t, double k) {
	return (1/k + log(z) * ( 1 - w));
}

inline double weibull_mean(double* p)
{
	return (p[1] + gsl_sf_gamma((p[2] + 1) / p[2]));
}

inline double weibull_median(double* p)
{
	return (p[1] * pow(log(2), 1 / p[2]));
}

/********************************************************
*	Compare horizontally at N - 1 quantile positions
*	x-axis in log scalse
*********************************************************/

//  Euclid distance

template <class ild_t>
double compqtl_EC(ild_t* a, ild_t* b, bool sqr = true)
{
	double  abc = 0.;
	for (int k = 0; k < lildprm.n_qtl - 1; ++k) {
	   float	d = a->k_th_qtil[k] - b->k_th_qtil[k];
	   abc += d * d;
	}
	abc /= (lildprm.n_qtl - 1);
	return (sqr? sqrt(abc): abc);
}

//  Mantattan distance

template <class ild_t>
double compqtl_MH(ild_t* a, ild_t* b)
{
	double  abc = 0.;
	for (int k = 0; k < lildprm.n_qtl - 1; ++k) {
	   abc += fabs(a->k_th_qtil[k] - b->k_th_qtil[k]);
	}
	abc /= (lildprm.n_qtl - 1);
	return (abc);
}

// Distance between two ilds

extern double dist_ilds(IldPrm* a, IldPrm* b, DistMethod mthd);
extern double dist_ilds(Ild* a, Ild* b, DistMethod mthd);
extern double dist_ilds(Lild* a, Lild* b, DistMethod mthd);

// list of species ananlyzed

extern SpecList* make_speclist(int nn, const char** mm);
extern void delete_speclist();

#endif	// _ILDPDF_H_
