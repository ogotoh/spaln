/*******************************************************************************
*
*	subset 
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "stdtype.h"

void Subset::rsubset(int n, const char* ps, FILE* fd)
{
	pool = new int[2 * n];
	group = new int*[n + 1];
	elms = 0;
	if (ps && *ps == '0') ps = 0;
	int*	pl = pool;
	int**	pg = group;

	if (n > 2 && (fd || ps)) {
	    if (fd) fgetiarray(0, 0, fd);	/* initialize buffer	*/
	    for (int j = n; j; *pl++ = EOTAB) {
		int	i;
		if (ps) {
		    i = sgetiarray(pl, j, &ps);
		    if (ps && *ps) ++ps;
		} else if (fd) {
		    i = fgetiarray(pl, j, fd);
		} else	break;
		if (!i) break;
		if (i > 0) {
		    *pg++ = pl;
		    j -= i;
		    while (i--) {
			(*pl++)--;
			++elms;
		    }
		} else {
		    j += i;
		    while (i++) {
			(*pl)--;
			*pg++ = pl++;
			*pl++ = EOTAB;
			++elms;
		    }
		}
	    }
	} else {
	    for (int j = 0; j < n; ++j) {
		*pg++ = pl;
		*pl++ = j;
		*pl++ = EOTAB;
	    }
	    elms = n;
	}
	*pg = 0;
	num = pg - group;
}

Subset::Subset(Subset& src)	// copy
{
	elms = src.elms;
	num = src.num;
	pool = new int[2 * elms];
	group = new int*[num + 1];
	int*	pl = pool;
	int**	pg = group;
	int*	sp = src.pool;
	for (int i = 0; i < num; ++i) {
	    *pg++ = pl;
	    while ((*pl++ = *sp++) != EOTAB) ;
	}
	*pg = 0;
}

Subset::Subset(int nn, const char* ps)
{
	FILE*	fg = 0;

	if (ps) {
	    if (*ps) {
		while (*ps && isspace(*ps)) ++ps;
		if (isdigit(*ps)) {
		    rsubset(nn, ps, 0);
		    return;
		}
		fg = fopen(ps, "r");
		if (!fg) fprintf(stderr, "%s can't open!\n", ps);
	    } else	fg = stdin;
	}
	rsubset(nn, 0, fg);
	if (fg && fg != stdin) fclose(fg);
}

Subset::~Subset()
{
	delete[] group;
	delete[] pool;
}

