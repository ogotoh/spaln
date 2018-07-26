/****************************************************************************
*
*	Subroutines for nucleotide/amino-acid sequence divergence
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
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japang
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "cmn.h"
#include "divseq.h"

static	DistCal	ReAlignB = ThisAln;
static	int	GapEval = 3;
static	int	GapEvalB = 3;
static	int	PickUpOneB = 0;
static	FTYPE	gapfrac = 0.8;

DISTPRM	distPrm = {0.5, 0.5, 0.0, 1.0, 0, ThisAln, 0};

FTYPE pairdvn(Seq* sd, int i, int j)
{
	CHAR*   ps = sd->at(sd->left);
	CHAR*   tt = sd->at(sd->right);

	int     mch = 0, mmc = 0, unp = 0, gap = 0;
	int     gsi = 0, gsj = 0;
	for ( ; ps < tt; ps += sd->many) {
	    if (IsGap(ps[i])) {
		if (!IsGap(ps[j])) {
		    ++unp;
		    if (gsi <= gsj) ++gap;
		    gsj = 0;
		}
		++gsi;
	    } else {
		if (IsGap(ps[j])) {
		    ++unp;
		    if (gsi >= gsj) ++gap;
		    ++gsj;
		} else {
		    if (ps[i] == ps[j]) ++mch;
		    else	++mmc;
		    gsj = 0;
		}
		gsi = 0;
	    }
	}
	FTYPE	gapunp = gapfrac * gap + (1 - gapfrac) * unp;
	return (1. - (FTYPE) mch / (gapunp + mch + mmc));
}

FTYPE pairdvn(Seq* sd, int* gi, int* gj)
{
	CHAR*   ps = sd->at(sd->left);
	CHAR*   tt = sd->at(sd->right);

	int	ni = 0, nj = 0;
	for (int* g = gi; *g >= 0; ++g) ++ni;
	for (int* g = gj; *g >= 0; ++g) ++nj;
	int*	gsi = new int[ni]; vclear(gsi, ni);
	int*	gsj = new int[nj]; vclear(gsj, nj);
	int	mch = 0, mmc = 0, unp = 0, gap = 0;
	for ( ; ps < tt; ps += sd->many) {
	  for (int i = 0; gi[i] >= 0; ++i) {
	    for (int j = 0; gj[j] >= 0; ++j) {
	      if (IsGap(ps[gi[i]])) {
		if (!IsGap(ps[gj[j]])) {
		    ++unp;
		    if (gsi[i] <= gsj[j]) ++gap;
		    gsj[j] = 0;
		}
		++gsi[i];
	      } else {
		if (IsGap(ps[gj[j]])) {
		    ++unp;
		    if (gsi[i] >= gsj[j]) ++gap;
		    ++gsj[j];
		} else {
		    if (ps[gi[i]] == ps[gj[j]]) ++mch;
		    else	++mmc;
		    gsj[j] = 0;
		}
		gsi[i] = 0;
	      }
	    }
	  }
	}
	delete[] gsi;
	delete[] gsj;
	FTYPE	gapunp = gapfrac * gap + (1 - gapfrac) * unp;
	return (1. - (FTYPE) mch / (gapunp + mch + mmc));
}

double jukcan(double nid)
{
	if (nid <= 0.) return (0.);
	nid = 1. - 4./3. * nid;
	if (nid <= 0.) return (OUTOFRANGE);
	return (-3./4. * log(nid));
}

double jukcanpro(double nid)
{
	if (nid <= 0.) return (0.);
	nid = 1. - 20./19. * nid;
	if (nid <= 0.) return (OUTOFRANGE);
	return (-19./20. * log(nid));
}

double pamcorrect(double x, int molc)
{
	if (x <= 0) return (0);
	if (distPrm.corr_mhits == 0) return (100. * x);
	double	y;
	double	alpha = 1.;
	if (molc == PROTEIN) {
	    if (distPrm.corr_mhits == 1)
		y = 1. - x / (alpha = 0.95);
	    else
		y = (x <= 0.7)? 1. - (0.987151 + 0.220560*x)*x:
		    -1.260444 + (8.603930  - (13.869219 - 6.521836*x)*x)*x; 
	} else
	    y = 1. - x / (alpha = 0.75);

	return (y > 0.)? -100. * alpha * log(y): OUTOFRANGE;
}

void setdivseq(int realn, int geval, int pick)
{
	setintval("Pre_alignment[0]/Rough_est[1]/Re_calculate[2] [%d]: ", 
	    (int*) &distPrm.realign,(int *) &ReAlignB, realn);
	setintval("Treat gap as missmatch [0-3] [%d]: ",
	    &GapEval, &GapEvalB, geval);
	switch (GapEval) {
	    case 0: distPrm.gap_v_rate = distPrm.gap_u_rate = 0; break;
	    case 1: distPrm.gap_v_rate = 1; distPrm.gap_u_rate = 0; break;
	    case 2: distPrm.gap_v_rate = 0; distPrm.gap_u_rate = 1; break;
	    case 3: distPrm.gap_v_rate = 0.5; distPrm.gap_u_rate = 0.5; break;
	    default: break;
	}
	setintval("Pick up one[1]/All[0] [%d]: ", 
	    &distPrm.pickup_one, &PickUpOneB, pick);
}

void setexprm_z(const char* ps)
{
const	char*	val = ps + 1;

	switch (*ps) {
	  case 'a':	// Realign
	    distPrm.realign = (*val && isnumber(val))? 
		(DistCal) atoi(val): DynAln;
	    break;
	  case 'j': case 'J':	// spj divergence equiv.
	    distPrm.contbt_spj = (*val && isnumber(val))?
		atof(val): alprm2.spb;
	    if (distPrm.contbt_spj > 1) distPrm.contbt_spj /= 100;
	    distPrm.contbt_seq = 1. - distPrm.contbt_spj;
	    break;
	  case 'm':	// Correct multi-hit
	    distPrm.corr_mhits = (*val && isdigit(*val))? atoi(val): 0;
	    break;
	  case 'u':	// Gap extension equiv.
	    if (*val && isnumber(val)) 
		distPrm.gap_u_rate = atof(val);
	    else {
		distPrm.gap_u_rate = 1.;
		distPrm.gap_v_rate = 0.;
	    }
	    break;
	  case 'v':	// Gap open equiv.
	    if (*val && isnumber(val))
		distPrm.gap_v_rate = atof(val);
	    else {
		distPrm.gap_v_rate = 1.;
		distPrm.gap_u_rate = 0.;
	    }
	    break;
	  default: break;
	}
}

double degdiv(FSTAT* stt)
{
	double	fd = stt->mmc;
	fd += distPrm.gap_v_rate * stt->gap;
	fd += distPrm.gap_u_rate * stt->unp;
	double	fn = fd + stt->mch;
	return (fn > 0.)? (fd / fn): 0.;
}

void put_stat(FILE* fd, FSTAT* fstt, double* ppc)
{
	static char frmt1[] = "%7.2lf %7.2lf %7.2lf ";
	static char frmt2[] = "%6.2f %5.1f %5.1f %5.1f %5.1f ";
	double	pc = ppc? *ppc: degdiv(fstt);

	if (distPrm.corr_mhits) {
	    double  sd = (fstt->mch > 0.)? sqrt(pc * (1. - pc) / fstt->mch): 0.;
	    fprintf(fd, frmt1, 100. * pc, 100. * sd, pamcorrect(pc));
	} else {
	    fprintf(fd, "%7.2f ", 100. * pc);
	}
	fprintf(fd, frmt2, fstt->val,
	    fstt->mch, fstt->mmc, fstt->gap, fstt->unp);
}

double degdiv(Gsinfo* gsi)
{
	if (!gsi) return (0);
	double  divseq = degdiv(&gsi->fstat);
	SigII*& sigii = gsi->sigII;
	if (distPrm.contbt_spj <= 0 || !sigii) return (divseq);
	double	divspj = sigii->n_common();
	divspj = 1. - 2 * divspj / sigii->lstnum;
	return (distPrm.contbt_seq * divseq + distPrm.contbt_spj * divspj);
}

void put_stat(FILE* fd, Gsinfo* gsi)
{
	double  pc = degdiv(gsi);
	put_stat(fd, &gsi->fstat, &pc);
}

