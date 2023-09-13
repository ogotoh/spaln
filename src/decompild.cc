/*****************************************************************************
*
*	decompild.c
*	characteristics of individual components of ILDs
*
*	Usage: decompild [-N] IldModel.txt
*
*	Requires GSL library 
*	https://www.gnu.org/software/gsl/
*
*	Osamu Gotoh, ph.D.      (-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.      (2001-2023)
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

#include"ildpdf.h"

void usage()
{
	fputs("Usage: decompild [-N] IldModel.txt\n", stdout);
	fputs("Options:\n", stdout);
	fputs("-N: lognormal (default: Frechet)\n", stdout);
	exit (1);
}

int main(int argc, const char* argv[])
{
	while (--argc > 0 && **++argv == '-') {
	    if (argv[0][1] == 'N') defpdf = LOGNORMAL; 
	    if (argv[0][1] == 'G') defpdf = GEOMETRIC; 
	}
	FILE* fd = argc? fopen(*argv, "r"): stdin;
	if (!fd) usage();
	char	str[MAXL];
	IldPrm	ildp;
	IldPrm	term;
	while (fgets(str, MAXL, fd)) {
	    if (*str == '#') continue;
	    ildp.get_IldPrm(str);
	    ildp.complete();
	    Strlist	stl(str, " ");
	    int	nterm = stl.size();
	    if (nterm - 11 < ildp.n_param) continue;
	    double*	dp = ildp.dparam;
	    double	qt = 0.;
	    int	n_sample = atoi(stl[2]);
	    printf("%s\t%s\t%7d\t%14s\t%3d", stl[0], stl[1], n_sample, stl[nterm - 4], ildp.n_modes);
	    for (int m = 0; m < ildp.n_modes; ++m, dp += ildp.m_size) {
		double	a = dp[ildp.m_size - 1];
		double	q14 = 0., q34 = 0., mode = 0., median = 0.;
		qt += a;
		if (defpdf == FRECHET) {
		    q14 = frechet_quantile(dp, 0.25);
		    q34 = frechet_quantile(dp, 0.75);
		    mode = frechet_mode(dp);
		    median = frechet_median(dp);
		} else if  (defpdf == LOGNORMAL) {
		    q14 = lognormal1_4quantile(dp);
		    q34 = lognormal3_4quantile(dp);
		    mode = lognormal_mode(dp);
		    median = lognormal_median(dp);
		} else if (defpdf == GEOMETRIC) {
		    q14 = geometric_quantile(dp, 0.25);
		    q34 = geometric_quantile(dp, 0/75);
		    mode = geometric_mode(dp);
		    median = geometric_median(dp);
		}
//		double	mean = frechet_mean(dp);
//		printf("\t%7.4f\t%7.2f\t%7.2f\t%7.2f\t%7.2f", 
//		    a, mode, median, mean, q34 - q14);
		printf("\t%7.4f\t%7.2f\t%7.2f\t%7.2f", a, mode, median, q34 - q14);
		if (m < ildp.n_modes - 1) {
		    double	boundary = ildp.quantile(qt, q34);
		    if (boundary < 0) boundary = ildp.quantile(qt, median);
		    if (boundary < 0) boundary = ildp.quantile(qt, q14);
		    if (boundary < 0) boundary = ildp.quantile(qt, 0.);
		    printf("\t%7.2f", boundary);
		}
	    }
	    putchar('\n');
	}
	if (fd != stdin) fclose(fd);
	return (0);
}
