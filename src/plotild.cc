/*****************************************************************************
*
*	plotild.c; Display intron length distribution with gnuplot
*
*	Requires gnuplot
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
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
*
*****************************************************************************/

#include "ildpdf.h"

static  const   char*   graph_out = 0;
static	const	char*	ip_stat = 0;
static	const	char*	ild_dir = "";

static	int	phylc = 0;
static	double	alpha = 0.;
static	DblDbl	xtransform = log10;
static	DblDbl	xinvtransf = exp10;
static	bool	xlogscale = true;
static	bool	ildents = false;
static	char	pdfchar = 'F';

void usage()
{
	fputs("Usage:\n\tplotild [options] [y1.ild y2.ild ...] [-d IldModel x1 x2 ...]\n", stderr);
	fputs("\nOptions:\n", stderr);
	fputs("\n\t-g[out]\t: graphic output\n", stderr);
	fputs("\t-i\t: ild components (false)\n", stderr);
	fputs("\t-l#\t: lower bound of plot in log_10 (1)\n", stderr);
	fputs("\t-s$\t: directory of ild files\n", stderr);
	fputs("\t-u#\t: upper bound of plot in log_10 (6)\n", stderr);
	fputs("\t-x#\t: number of bins (100)\n", stderr);
	fputs("\t-C\t: plot cdf (plot pdf)\n", stderr);
	fputs("\t-L#\t: shortest intron length (0)\n", stderr);
        fputs("\t-ME\t: exponential (geometoric) model\n", stderr);
        fputs("\t-MF\t: Frechet model\n", stderr);
        fputs("\t-MG\t: Gamma model\n", stderr);
        fputs("\t-MN\t: log normal model\n", stderr);
	fputs("\t-P\t: plot penalty (plot pdf)\n", stderr);
	fputs("\t-U#\t: longest intron length(inf)\n", stderr);
	fputs("\t-X\t: X-axix in linear scale (false)\n", stderr);
	exit(1);
}

void getoptions(int& argc, const char**& argv)
{
	const	char*	pn = 0;
	while (--argc && **++argv == '-') {
	  switch (argv[0][1]) {
	    case 'a':
		if ((pn = getarg(argc, argv, true)))
		    alpha = atof(pn);
		break;
	    case 'c':
		if ((pn = getarg(argc, argv, true)))
		    phylc = atoi(pn);
		break;
	    case 'd': ip_stat = getarg(argc, argv); break;
	    case 'g': graph_out = getarg(argc, argv);
		if (graph_out && !gnuplot_terminal(graph_out)) {
		    ++argc; --argv; graph_out = 0;
		}
		break;
	    case 'h': case '?': usage();
	    case 'i': ildents = true; break;
	    case 'l':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.llmt = atof(pn);
		break;
	    case 'q': setprompt(false, false); break;
	    case 's': ild_dir  = getarg(argc, argv); break;
	    case 'u':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.ulmt = atof(pn);
		break;
	    case 'x':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.ndiv = atoi(pn);
		break;
	    case 'C': ildoutmode = IldOutMode::CDF; break;
	    case 'L':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.minx = atoi(pn);
		break;
            case 'M':   // statistical model
             switch (argv[0][2]) {
              case 'E': defpdf = GEOMETRIC; pdfchar = 'E'; break;
              case 'F': defpdf = FRECHET; pdfchar = 'F'; break;
              case 'G': defpdf = GAMMA; pdfchar = 'G'; break;
              case 'N': defpdf = LOGNORMAL; pdfchar = 'N'; break;
              case 'W': defpdf = WEIBULL; pdfchar = 'W'; break;
             }
             break;
	    case 'P': ildoutmode = IldOutMode::Penalty; break;
	    case 'U':
		if ((pn = getarg(argc, argv, true)))
		    lildprm.maxx = atoi(pn);
		break;
	    case 'X': xlogscale = false; break;
	    default: break;
	  }
	}
}

static double trans(double x)
{
	return (1 - pow(x, -1. / alpha));
}

static double invtr(double y)
{
	return (pow(1 - y, -alpha));
}

static void getiprms(IldPrm*& iprm, int& nprm, int& argc, const char**& argv)
{
	FILE*	fd = ftable.fopen(ip_stat, "r");
	if (!fd) fatal("%s not found !\n", ip_stat);
	SpecList*	speclist = make_speclist(argc, argv);
	int	nn = speclist->size();
	argc -= nn; argv += nn;
	IldPrm*	ipprm = new IldPrm[nn + nprm + 1];
	if (iprm) {
	    vcopy(ipprm, iprm, nprm);
	    delete[] iprm;
	}
	iprm = ipprm;
	int	rv;
	while ((rv = iprm[nprm].fget(fd)) != EOF)
	    if (rv == OK) ++nprm;
	fclose(fd);
	delete_speclist();
}

static StrHash<int>* read_gnm2tab(int mems)
{
	FILE*	fd = ftable.fopen("gnm2tab", "r");
	if (!fd) fatal("Can't open gnm2tab !\n");
	char	str[MAXL];
	StrHash<int>	domphy(mems);
	StrHash<int>*	phyl_code = new StrHash<int>(mems);
	while (fgets(str, MAXL, fd)) {
	    char*	ps = str;
	    if (*ps == '#' || *ps == '\n') continue;
	    char*	id = car(ps);
	    char*	qs = 0;
	    int	c = 1;
	    for ( ; c < phylc; ++c) {
		if (!*++ps || !(qs = car(ps))) break;
	    }
	    if (c < phylc) continue;
	    KVpair<INT, int>*	kv = domphy.pile(qs);
	    phyl_code->assign(id, kv->val);
	}
	fclose(fd);
	return (phyl_code);
}

int main(int argc, const char *argv[])
{
	IldPrm*	iprm = 0;
	Lild**	lild = 0;
	Ild**	ilds = 0;
	int	nprm = 0;
	int	nild = 0;
	getoptions(argc, argv);
	if (alpha > 0) {	// 1/x
	    xtransform = trans;
	    xinvtransf = invtr;
	}
	int	mems = argc;
	if (ip_stat) {
	    getiprms(iprm, nprm, argc, argv);
	    ++argc; --argv;
	    getoptions(argc, argv);
	    ip_stat = 0;
	}
	if (argc) {
	    if (xlogscale) lild = new Lild*[argc];
	    else	ilds = new Ild*[argc];
	    char	str[MAXL];
	    strcpy(str, ild_dir);
	    int	ild_dir_len = strlen(str);
	    char*	ps = str + ild_dir_len;
	    if (ild_dir_len && ps[-1] != '/')
		*ps++ = '/';
	    for ( ; argc && **argv != '-'; --argc, ++argv, ++nild) {
		strcpy(ps, *argv);
		FILE*	fd = fopen(str, "r");
		if (!fd) fatal("%s not found !\n", str);
		if (xlogscale) {
		    lild[nild] = new Lild(fd, *argv, xtransform, xinvtransf);
		} else {
		    ilds[nild] = new Ild(fd, *argv, nild);
		    ilds[nild]->normalize();
		}
		fclose(fd);
	    }
	}
	if (argc > 0) {
	    ++argc; --argv;
	    getoptions(argc, argv);
	    if (ip_stat) getiprms(iprm, nprm, argc, argv);
	}
	StrHash<int>*	phyl_code = phylc? read_gnm2tab(mems): 0;
	if (nprm || nild) {
	    if (xlogscale) {
		GnuPlotLild	gp(iprm, nprm, lild, nild, ildents);
		gp.plot(graph_out, 0, phyl_code, ildents? iprm: 0);
	    } else {
		GnuPlotLild	gp(iprm, nprm, ilds, nild, ildents);
		gp.plot(graph_out, 0, phyl_code, ildents? iprm: 0);
	    }
	}
	delete[] iprm;
	if (lild) {
	    for (int i = 0; i < nild; ++i) delete lild[i];
	    delete[] lild;
	}
	if (ilds) {
	    for (int i = 0; i < nild; ++i) delete ilds[i];
	    delete[] ilds;
	}
	delete phyl_code;
	return 0;
}

