/**********************************************************************
*
*	Generate and mutate random nucleotide/amino acid
*	sequences for Monte Carlo Test
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
*	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
*
*****************************************************************************/

#include <math.h>
#include "seq.h"
#include "montseq.h"

static	int	randres();
extern double	drand48();
static double	cumu[21];
static int	Nelm;
static int	rno = 0;

/*
*	Knuth, pg 117, 3.4.1, vol. 2 2nd edition
*/

double rgauss()
{
    double v1, v2, s, x1;

    do {
        v1 = 2.*drand48()-1.;
        v2 = 2.*drand48()-1.;
        s = v1*v1+v2*v2;
    } while (s >= 1);
    s = sqrt((-2.*log(s))/s);
    x1 = v1*s;
    return x1;
}

/*	Random numbers that follow Poisson distribution	*/

int rpoisson(double mu)
{
	int	i = 0;
	double	p = drand48();

	mu = exp(-mu);
	for ( ; p >= mu; ++i) p *= drand48();
	return (i);
}

void initrandseq(double acgt[], int molc)
{
	double	acum = 0.;

	Nelm = molc == PROTEIN? 20: 4;
	double	t = 0;
	for (int i = 0; i < Nelm; i++) t += acgt[i];
	for (int i = 0; i < Nelm; i++) {
	    cumu[i] = acum;
	    acgt[i] /= t;
	    acum += acgt[i];
	}
	cumu[Nelm] = 1.;
	rno = 0;
}

static int randres()
{
	double	x = drand48();
	int	k = (int) (Nelm * x);

	while (cumu[k+1] <= x) ++k;
	while (cumu[k] > x) --k;
	return (Nelm == 4? (1 << k) + _: k + ALA);
}

Seq* randseq(Seq* sd, int len)
{
	char	name[MAXL];

	sd->refresh(1, len);
	snprintf(name, MAXL, "RAND%d", ++rno);
	if (sd->sname) sd->sname->assign(name);
	else	sd->sname = new Strlist(&name[0], "");
	CHAR*	s = sd->at(0);
	for (int n = 0; n < len; n++)
	    *s++ = randres();
	sd->postseq(s);
	return (sd->rndseq());
}

/*
	Random substitution of 'n' times
	if which >= 0 which's member is substituted
 	else substitutions may occur in any members
*/

Seq* substseq(Seq* sd, int n, int which)
{
	CHAR*	ss = sd->at(sd->left);
	int	area = sd->right - sd->left;

	if (which >= 0)  ss += which;
	else	    area *= sd->byte;
	while (n > 0) {
	    int	k = (int) (drand48() * area);
	    CHAR*	ps = ss + (which >= 0? k * sd->byte: k);
	    if (IsGap(*ps)) continue;
	    while ((k = randres()) == *ps) ;
	    *ps = k;
	    --n;
	}
	return (sd);
}

