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

#include "dist2.h"
#include "ildpdf.h"

static	DistMethod distmeth = EC;
static	OutMode	omode = DMX;
static	bool	logscale = false;
static	bool	smooth = false;
static	const	char*	ip_stat = 0;
static	const	char*	graph_out = 0;

void usage()
{
	fputs("Usage:\tcompild [options] [Ild1..IldN]\n", stderr);
	fputs("  or\tcompild [options] -d IldModel.txt\n", stderr);
	fputs("Options:\n", stderr);
	fputs("\t-i[a|e|f|g|l|p] (input mode)\n", stderr);
	fputs("\t\ta: alternative\n",  stderr);
	fputs("\t\te: every pair\n",  stderr);
	fputs("\t\tf: first and others\n",  stderr);
	fputs("\t\tg: between groups\n",  stderr);
	fputs("\t\tl: last and others\n",  stderr);
	fputs("\t\tp: alternative\n",  stderr);
	fputs("\t-[a|c|e|f|j|k|m|s] (distance measure)\n", stderr);
	fputs("\t\ta: disjoint area\n",  stderr);
	fputs("\t\tc: cosine\n",  stderr);
	fputs("\t\te: Euclid\n",  stderr);
	fputs("\t\tf: Euclid^2\n",  stderr);
	fputs("\t\tj: Jaccard\n",  stderr);
	fputs("\t\tk: Kullback-Leibler\n",  stderr);
	fputs("\t\tm: Manhattan\n",  stderr);
	fputs("\t\ts: Jensen-Shannon\n",  stderr);
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
	while (--argc && **++argv == '-') {
	  switch (argv[0][1]) {
	    case 'd': ip_stat = getarg(argc, argv); break;
	    case 'L': logscale = true; break;
	    default: break;
	  }
	}
}

int localoption(int& argc, const char**& argv)
{
	int	n = 1;

	switch (argv[0][1]) {
	    case 'a': distmeth = UA; break;
	    case 'c': distmeth = CS; break;
	    case 'e': distmeth = EC; break;
	    case 'f': distmeth = E2; break;
	    case 'g': graph_out = getarg(argc, argv);
		if (!graph_out) graph_out = "";
		else if (!gnuplot_terminal(graph_out)) {
		    ++argc; --argv; graph_out = "";
		}
		break;
	    case 'h': case '?': usage();
	    case 'j': distmeth = JA; break;
	    case 'k': distmeth = KL; break;
	    case 'm': distmeth = MH; break;
	    case 'q': setprompt(0, 0); break;
	    case 's': distmeth = JS; break;
	    case 'x': lildprm.maxx = atoi(getarg(argc, argv, true)); break;
	    case 'E': gslprm.epsrel = atof(getarg(argc, argv, true)); break;
	    case 'H': lildprm.minfreq = atoi(getarg(argc, argv, true)); break;
	    case 'O': omode = PAIR; break;
	    case 'P': lildprm.psdcnt = atof(getarg(argc, argv, true)); break;
	    case 'Q': 
		lildprm.n_qtl = atoi(getarg(argc, argv, true));
		if (lildprm.n_qtl < 2) lildprm.n_qtl = 0;
		break;
	    case 'S': smooth = true; break;
	    default: n = 0; break;
	}
	return (n);
}

/*************************************************************************
*
*	Experimental ILDs in linear scale x-axis
*
*************************************************************************/

int ild_main(CalcServer<Ild>* svr, Ild** ild, ThQueue<Ild>* q = 0)
{
	Dist2<Ild>*	icp = (Dist2<Ild>*) svr->prm;
	Ild*&	ic = ild[0];
	const	char*	sl = strrchr(ic->fname, '/');
	m_thread_Lock(q);
	icp->sname.push(sl? sl + 1: ic->fname);
	m_thread_Unlock(q);
	icp->vars[ic->sid] = new Ild;
	icp->nmidx[ic->sid] = icp->inodr++;
	swap(icp->vars[ic->sid], ic);
	return (OK);
}

int compild_main(CalcServer<Ild>* svr, Ild* ild[], ThQueue<Ild>* q = 0)
{
	Ild*&	a = ild[0];
	Ild*&	b = ild[1];
	Dist2<Ild>*	icp = (Dist2<Ild>*) svr->prm;
	double	dst = dist_ilds(a, b, distmeth);
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

int lild_main(CalcServer<Lild>* svr, Lild** ild, ThQueue<Lild>* q = 0)
{
	Dist2<Lild>*	icp = (Dist2<Lild>*) svr->prm;
	Lild*&	ic = ild[0];
	const	char*	sl = strrchr(ic->fname, '/');
	icp->sname.push(sl? sl + 1: ic->fname);
	if (smooth) ic->smooth();
	icp->vars[icp->inodr] = new Lild();
	swap(icp->vars[icp->inodr], ild[0]);
	icp->nmidx[ic->sid] = icp->inodr++;
	return (OK);
}

int complild_main(CalcServer<Lild>* svr, Lild* ild[], ThQueue<Lild>* q = 0)
{
	Lild*&	a = ild[0];
	Lild*&	b = ild[1];
	Dist2<Lild>*	icp = (Dist2<Lild>*) svr->prm;
	double	dst = dist_ilds(a, b, distmeth);
	int	nn = svr->calcnbr(a->sid, b->sid);
	m_thread_Lock(q);
	icp->dist[nn] = (FTYPE) dst;
	m_thread_Unlock(q);
	return (OK);
}

/*************************************************************************
*
*	Statistically modeled ILDs
*
*************************************************************************/

static int spd_main(CalcServer<IldPrm>* svr, IldPrm** spd, ThQueue<IldPrm>* q = 0)
{
	Dist2<IldPrm>*	icp = (Dist2<IldPrm>*) svr->prm;
	IldPrm*&	sp = spd[0];
	icp->sname.push(sp->sname);
	icp->nmidx[sp->sid] = icp->inodr++;
	icp->vars[sp->sid] = new IldPrm;
	swap(icp->vars[sp->sid], sp);
	return (OK);
}

static int compspd_main(CalcServer<IldPrm>* svr, IldPrm* ab[], ThQueue<IldPrm>* q = 0)
{
	IldPrm*&	a = ab[0];
	IldPrm*&	b = ab[1];
	double	result = dist_ilds(a, b, distmeth);
	int	nn = svr->calcnbr(a->sid, b->sid);
	Dist2<IldPrm>*	icp = (Dist2<IldPrm>*) svr->prm;
	m_thread_Lock(q);
	icp->dist[nn] = result;
	m_thread_Unlock(q);
	return (OK);
}

static int spd_out2(CalcServer<IldPrm>* svr, IldPrm* spd[], ThQueue<IldPrm>* q = 0)
{
	IldPrm*&	a = spd[0];
	IldPrm*&	b = spd[1];
	Dist2<IldPrm>*	icp = (Dist2<IldPrm>*) svr->prm;
	double	dst = icp->dist[svr->calcnbr(a->sid, b->sid)];
	int	i = icp->nmidx[a->sid];
	int	j = icp->nmidx[b->sid];
	fprintf(out_fd, "%14.6e\t%d\t%d\t%s\t%s\n", (float) dst, 
	    a->n_sample, b->n_sample,
	    icp->sname[i], icp->sname[j]);
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
	Dist2<Ild>*	icp = (Dist2<Ild>*) svr->prm;
	double	dst = icp->dist[svr->calcnbr(a->sid, b->sid)];
	int	i = icp->nmidx[a->sid];
	int	j = icp->nmidx[b->sid];
	fprintf(out_fd, "%14.6e\t%.1f\t%.1f\t%s\t%s\n", 100. * dst, 
	    icp->vars[a->sid]->ftotal, icp->vars[b->sid]->ftotal,
	    icp->sname[i], icp->sname[j]);
	return (OK);
}

int lild_out(CalcServer<Lild>* svr, Lild* ild[], ThQueue<Lild>* q = 0)
{
	Lild*&	a = ild[0];
	Lild*&	b = ild[1];
	Dist2<Lild>*	icp = (Dist2<Lild>*) svr->prm;
	double	dst = icp->dist[svr->calcnbr(a->sid, b->sid)];
	int	i = icp->nmidx[a->sid];
	int	j = icp->nmidx[b->sid];
	fprintf(out_fd, "%14.6e\t%.1f\t%.1f\t%s\t%s\n", 100. * dst, 
	    icp->vars[a->sid]->ftotal, icp->vars[b->sid]->ftotal,
	    icp->sname[i], icp->sname[j]);
	return (OK);
}

int main(int argc, const char** argv)
{
	int	ac = argc;
const	char**	av = argv;
	getoption(ac, av);

	if (ip_stat) {		// statistical models
	    if (ac) {
		make_speclist(ac, av);
		*av++ = ip_stat;
		argc = av - argv;
	    }
	    Dist2<IldPrm> dist2;
	    CalcServer<IldPrm> csvr(argc, argv, IM_SNGL, IM_EVRY, &dist2,
		&spd_main, &localoption);
	    dist2.prepare(csvr);
	    CalcServer<IldPrm> ssvr(dist2.calc_mode, &dist2, &compspd_main, 
		0, 0, dist2.vars, dist2.num, dist2.grp2);
	    int	rv = ssvr.autocomp();
	    if (rv) return (rv);
	    if (omode == PAIR || dist2.calc_mode != IM_EVRY) {
		ssvr.change_job(&spd_out2);
		rv = ssvr.autocomp(false);
	    } else dist2.outdmx();
	    if (graph_out) {
		GnuPlotLild gp(dist2.vars, dist2.num);
		gp.plot(*graph_out? graph_out: 0);
	    }
	    return (rv);
	}
	AddExt	adext(argc, argv, ".ild");
	const	char**	argild = adext.add_ext();
	if (logscale) {
	    Dist2<Lild>	dist2;
	    CalcServer<Lild> csvr(argc, argild, IM_SNGL, IM_EVRY,
		 (void*) &dist2, &lild_main, &localoption);
	    dist2.prepare(csvr);
	    CalcServer<Lild> ssvr(dist2.calc_mode, &dist2, &complild_main, 
		0, 0, dist2.vars, dist2.num, dist2.grp2);
	    int	rv = ssvr.autocomp();
	    if (rv) return (rv);
	    if (omode == PAIR || dist2.calc_mode != IM_EVRY) {
		ssvr.change_job(&lild_out);
		rv = ssvr.autocomp(false);
	    } else dist2.outdmx();
	    if (graph_out) {
		GnuPlotLild	gp(dist2.vars, dist2.num);
		gp.plot(*graph_out? graph_out: 0, &dist2.sname);
	    }
	    return (rv);
	} else {
	    Dist2<Ild>	dist2;
	    CalcServer<Ild> csvr(argc, argild, IM_SNGL, IM_EVRY,
		 &dist2, &ild_main, &localoption);
	    dist2.prepare(csvr);
	    CalcServer<Ild> ssvr(dist2.calc_mode, &dist2, &compild_main, 
		0, 0, dist2.vars, dist2.num, dist2.grp2);
	    int	rv = ssvr.autocomp();
	    if (rv) return (rv);
	    if (omode == PAIR || dist2.calc_mode != IM_EVRY) {
		ssvr.change_job(&ild_out);
		rv = ssvr.autocomp(false);
	    } else dist2.outdmx();
	    return (rv);
	}
}
