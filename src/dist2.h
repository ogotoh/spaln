/*****************************************************************************
*
*	header for pairwise calculation of distance
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

#ifndef  _DIST2_
#define  _DIST2_

#include "cmn.h"
#include "calcserv.h"

enum DistMethod {EC, JS, JA, KL, MH, UA, E2, CS, XS};
enum OutMode {DMX, PAIR};

template <class var_t>
struct Dist2 {
	int	calc_mode;
	int	num;
	int	grp2;;
	var_t**	vars;
	Strlist	sname;
	int*	nmidx;
	int	inodr;
	FTYPE*	dist;
	FTYPE*	getdist() {FTYPE* rv = 0; swap(rv, dist); return (rv);}
	void	prepare(CalcServer<var_t>& svr);
	void	outdmx();
	Dist2<var_t>() : calc_mode(0), num(0), grp2(0), vars(0),
		nmidx(0), inodr(0), dist(0) {}
	~Dist2<var_t>() {;
	    for (int n = 0; n < num; ++n) delete vars[n];
	    delete[] vars; delete[] nmidx; delete[] dist;
	}
};

template <class var_t>
void Dist2<var_t>::prepare(CalcServer<var_t>& csvr)
{
	num = csvr.memsize(csvr.ldr1);
	if (num < 2) fatal("Too small # of data");
	vars = new var_t*[num];
	nmidx = new int[num];
	inodr = 0;
	calc_mode = csvr.calc_mode;
	grp2 = csvr.getgrp22();
	int	nproc = csvr.auto_comp();
	swap(num, nproc);
	nproc -= num;
//	if (nproc && !ip_stat) prompt("%d inputs has been abondanded!\n", nproc);
	int     dist_size = csvr.calcsize(num, csvr.calc_mode);
	if (dist_size <= 0) fatal("Too small number of inputs !\n");
	else	dist = new FTYPE[dist_size];
}

template <class var_t>
void Dist2<var_t>::outdmx()
{
	for (int i = 0; i < num; ++i) printf("%s\n", sname[i]);
	putchar('\n');
	for (int j = 1, k = 0; j < num; ++j) {
	    int	n = 0;
	    for (int i = 0; i < j; ++i, ++k) {
		printf(" %15.7e", 100. * dist[k]);
		if (++n == 5) {
		    putchar('\n');
		    n = 0;
		}
	    }
	    if (n) putchar('\n');
	}
}

#endif	// _DIST2_

