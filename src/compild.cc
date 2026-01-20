/*****************************************************************************
*
*	Compare intron length distributions 
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
#include "dist2.h"
#include "ildpdf.h"

static	DistMethod metrics = KL;
static	OutMode	omode = DMX;
static	int	calc_mode = IM_EVRY;
static	bool	logscale = false;
static	bool	smooth = false;
static	const	char*	ip_stat = 0;
static	const	char*	graph_out = 0;
static	int	text_out = 1;
static	bool	gzip = false;
const	char*	oname = 0;

void usage()
{
	fputs("Usage:\tcompild [options] Ild1[.dgz]..IldN[.dgz]\n", stderr);
	fputs("  or\tcompild [options] -d IldModel.txt\n", stderr);
	fputs("Options:\n", stderr);
	fputs("\t-i[a|e|f|g|l|p] (input mode)\n", stderr);
	fputs("\t\ta: alternative\n",  stderr);
	fputs("\t\te: every pair\n",  stderr);
	fputs("\t\tf: first and others\n",  stderr);
	fputs("\t\tg: between groups\n",  stderr);
	fputs("\t\tl: last and others\n",  stderr);
	fputs("\t\tp: alternative\n",  stderr);
	fputs("\t-m[a|c|e|f|j|k|m|s] (metrics)\n", stderr);
	fputs("\t\ta: disjoint area\n",  stderr);
	fputs("\t\tc: cosine\n",  stderr);
	fputs("\t\te: Euclid\n",  stderr);
	fputs("\t\tf: Euclid^2\n",  stderr);
	fputs("\t\tj: Jaccard\n",  stderr);
	fputs("\t\tk: Kullback-Leibler\n",  stderr);
	fputs("\t\tm: Manhattan\n",  stderr);
	fputs("\t\ts: Jensen-Shannon\n",  stderr);
	fputs("\t-q:  suppress warnig messages\n", stderr);
	fputs("\t-xN  maximum intron length\n", stderr);
	fputs("\t-HN  minimum frequency\n", stderr);
	fputs("\t-L   log scalse x-axis\n", stderr);
	fputs("\t-O   each pair in one line\n", stderr);
	fputs("\t-QN  Horizontal at N quantile points\n", stderr);
	fputs("Notes:\n", stderr);
	fputs("\t1) Not all distance measures are supported depending on\n", stderr);
	fputs("\t   w/wo -d IldModel.txt or -L options\n", stderr);
	fputs("\t2) With -d IldModel.txt, only -ie input mode is currently supported\n", stderr);
	exit (1);
}

static void getoption(int& argc, const char**& argv)
{
	const	char*	val;
	INT	nthr = 0;
	while (--argc && **++argv == '-') {
	  switch (argv[0][1]) {
	    case 'b':
		gzip = 1;
		oname = getarg(argc, argv);
		break;
	    case 'd': ip_stat = getarg(argc, argv); break;
	    case 'g':
		val = getarg(argc, argv);
		if (!graph_out) graph_out = "";
		else if (!gnuplot_terminal(graph_out)) {
		    ++argc; --argv; graph_out = "";
		}
		break;
	    case 'h': case '?': usage();
	    case 'i': 		// input order
		switch (tolower(argv[0][2])) {
		  case 'a': case 'j':	// every juxt. pairs
		    calc_mode = IM_ALTR; break;
		  case 'e':		// lower-left combinaton
		    calc_mode = IM_EVRY; break;
		  case 'f':		// first vs others
		    calc_mode = IM_FvsO; break;
		  case 'g':		// group vs group, 1 file
		    calc_mode = IM_GRUP; break;
		  case 'i':		// self comparison
		    calc_mode = IM_SELF; break;
		  case 'l':		// last vs others
		    calc_mode = IM_OvsL; break;
		  case 'p':		// parallel, 2 files
		    calc_mode = IM_PARA; break;
		  default:
		    calc_mode = IM_SNGL; break;
		}
		break;
	    case 'm':
		switch (argv[0][2]) {
		    case 'a': metrics = UA; break;
		    case 'c': metrics = CS; break;
		    case 'e': metrics = EC; break;
		    case 'f': metrics = E2; break;
		    case 'j': metrics = JA; break;
		    case 'k': metrics = KL; break;
		    case 'm': metrics = MH; break;
		    case 's': metrics = JS; break;
		}
		break;
	    case 'o':
		text_out = 1;
		oname = getarg(argc, argv);
		break;
	    case 'q': setprompt(0, 0); break;
	    case 't': val = getarg(argc, argv, true);
	        thread_num = val? atoi(val): -1;
		break;
	    case 'x': lildprm.maxx = atoi(getarg(argc, argv, true)); break;
	    case 'E': gslprm.epsrel = atof(getarg(argc, argv, true)); break;
	    case 'H':
		nthr = atoi(getarg(argc, argv, true));
		lildprm.minfreq = nthr;
		set_min_samples(nthr);
		break;
	    case 'L': logscale = true; break;
	    case 'O': omode = PAIR; break;
	    case 'P': lildprm.psdcnt = atof(getarg(argc, argv, true)); break;
	    case 'Q': 
		lildprm.n_qtl = atoi(getarg(argc, argv, true));
		if (lildprm.n_qtl < 2) lildprm.n_qtl = 0;
		break;
	    case 'S': smooth = true; break;
	    default: break;
	  }
	}
}

/*************************************************************************
*
*	Experimental ILDs in linear scale x-axis
*
*************************************************************************/

int compild_main(CalcServer<Ild>* svr, Ild* ild[], ThQueue<Ild>* q = 0)
{
	Ild*&	a = ild[0];
	Ild*&	b = ild[1];
	Dist2<Ild>*	icp = (Dist2<Ild>*) svr->prm;
	FTYPE	dst = def_scale * dist_ilds(a, b, metrics);
	int	nn = svr->calcnbr(a->sid, b->sid);
	m_thread_Lock(q);
	icp->dist[nn] = dst;
	m_thread_Unlock(q);
	return (OK);
}

/*************************************************************************
*
*	Experimental ILDs in log scale x-axis
*
*************************************************************************/

int complild_main(CalcServer<Lild>* svr, Lild* ild[], ThQueue<Lild>* q = 0)
{
	Lild*&	a = ild[0];
	Lild*&	b = ild[1];
	Dist2<Lild>*	icp = (Dist2<Lild>*) svr->prm;
	FTYPE	dst = def_scale * dist_ilds(a, b, metrics);
	int	nn = svr->calcnbr(a->sid, b->sid);
	m_thread_Lock(q);
	icp->dist[nn] = dst;
	m_thread_Unlock(q);
	return (OK);
}

/*************************************************************************
*
*	Statistically modeled ILDs
*
*************************************************************************/

static int compspd_main(CalcServer<IldPrm>* svr, IldPrm* ab[], ThQueue<IldPrm>* q = 0)
{
	IldPrm*&	a = ab[0];
	IldPrm*&	b = ab[1];
	FTYPE	dst = def_scale * dist_ilds(a, b, metrics);
	int	nn = svr->calcnbr(a->sid, b->sid);
	Dist2<IldPrm>*	icp = (Dist2<IldPrm>*) svr->prm;
	m_thread_Lock(q);
	icp->dist[nn] = dst;
	m_thread_Unlock(q);
	return (OK);
}

static int spd_out2(CalcServer<IldPrm>* svr, IldPrm* spd[], ThQueue<IldPrm>* q = 0)
{
	IldPrm*&	a = spd[0];
	IldPrm*&	b = spd[1];
	int&	i = a->sid;
	int&	j = b->sid;
	Dist2<IldPrm>*	icp = (Dist2<IldPrm>*) svr->prm;
	FTYPE	dst = icp->dist[svr->calcnbr(i, j)];
	fprintf(out_fd, "%14.6e\t%d\t%d\t%s\t%s\n", dst, 
	    a->n_sample, b->n_sample,
	    (*icp->sname)[i], (*icp->sname)[j]);
	return (OK);
}

/*************************************************************************
*
*	main function of compild
*
*************************************************************************/

int ild_out(CalcServer<Ild>* svr, Ild* ild[], ThQueue<Ild>* q = 0)
{
	Ild*&	a = ild[0];
	Ild*&	b = ild[1];
	int&	i = a->sid;
	int&	j = b->sid;
	Dist2<Ild>*	icp = (Dist2<Ild>*) svr->prm;
	FTYPE	dst = icp->dist[svr->calcnbr(i, j)];
	fprintf(out_fd, "%14.6e\t%.1f\t%.1f\t%s\t%s\n", dst, 
	    icp->vars[i]->ftotal, icp->vars[j]->ftotal,
	    (*icp->sname)[i], (*icp->sname)[j]);
	return (OK);
}

int lild_out(CalcServer<Lild>* svr, Lild* ild[], ThQueue<Lild>* q = 0)
{
	Lild*&	a = ild[0];
	Lild*&	b = ild[1];
	Dist2<Lild>*	icp = (Dist2<Lild>*) svr->prm;
	int&	i = a->sid;
	int&	j = b->sid;
	double	dst = icp->dist[svr->calcnbr(i, j)];
	fprintf(out_fd, "%14.6e\t%.1f\t%.1f\t%s\t%s\n", dst, 
	    icp->vars[i]->ftotal, icp->vars[j]->ftotal,
	    (*icp->sname)[i], (*icp->sname)[j]);
	return (OK);
}

int main(int argc, const char** argv)
{
	int	ac = argc;
const	char**	av = argv;
	getoption(ac, av);

	if (ip_stat) {		// statistical models
	    Dist2<IldPrm> dist2(ip_stat, calc_mode);
	    CalcServer<IldPrm> ssvr(calc_mode, &dist2, &compspd_main, 
		0, 0, dist2.vars, dist2.mem_num, dist2.grp2, true, false);
	    int	rv = ssvr.autocomp();
	    if (rv) return (rv);
	    if (omode == PAIR || calc_mode != IM_EVRY) {
		ssvr.change_job(&spd_out2);
		rv = ssvr.autocomp(false);
	    } else {
		dist2.to_file(oname, text_out, gzip);
	    }
	    if (graph_out) {
		GnuPlotLild gp(dist2.vars, dist2.mem_num);
		gp.plot(*graph_out? graph_out: 0, dist2.sname);
	    }
	    return (rv);
	}
	if (logscale) {
	    Dist2<Lild>	dist2(argc, argv, calc_mode);
	    CalcServer<Lild> ssvr(calc_mode, &dist2, &complild_main, 
		0, 0, dist2.vars, dist2.mem_num, dist2.grp2);
	    int	rv = ssvr.autocomp();
	    if (rv) return (rv);
	    if (omode == PAIR || calc_mode != IM_EVRY) {
		ssvr.change_job(&lild_out);
		rv = ssvr.autocomp(false);
	    } else dist2.to_file(oname, text_out, gzip);
	    if (graph_out) {
		GnuPlotLild	gp(dist2.vars, dist2.mem_num);
		gp.plot(*graph_out? graph_out: 0, dist2.sname);
	    }
	    return (rv);
	} else {
	    Dist2<Ild>	dist2(argc, argv, calc_mode);
	    CalcServer<Ild> ssvr(calc_mode, &dist2, &compild_main, 
		0, 0, dist2.vars, dist2.mem_num, dist2.grp2);
	    int	rv = ssvr.autocomp();
	    if (rv) return (rv);
	    if (omode == PAIR || calc_mode != IM_EVRY) {
		ssvr.change_job(&ild_out);
		rv = ssvr.autocomp(false);
	    } else dist2.to_file(oname, text_out, gzip);
	    return (rv);
	}
}
