/*****************************************************************************
*
*	Collection of headers for calculation of sequence divergence
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

#ifndef  _DIVSEQ_H_
#define  _DIVSEQ_H_

#include <math.h>
#include "aln.h"

static	const	INT	OUTOFRANGE = 1024;

enum DistCal {ThisAln, Composition, DynAln, DynScr, GeneOrg};

struct DISTPRM {
	float	gap_v_rate;
	float	gap_u_rate;
	float	contbt_spj;
	float	contbt_seq;
	int	corr_mhits;
	DistCal	realign;
	int	pickup_one;
};

extern	double	jukcan(double nid);
extern	double	jukcanpro(double nid);
extern	double	pamcorrect(double x, int molc = PROTEIN);

extern	DISTPRM	distPrm;
extern	void	setdivseq(int realn, int geval, int pick);
extern	void	setexprm_z(const char* s);
extern	double	degdiv(FSTAT* stt);
extern	FTYPE	pairdvn(Seq* sd, int i, int j);
extern	FTYPE	pairdvn(Seq* sd, int* g1, int* g2);
extern	void	put_stat(FILE* fd, FSTAT* fstt, double* ppc = 0);
extern	double	degdiv(Gsinfo* gsi);
extern	void	put_stat(FILE* fd, Gsinfo* gsi);

#endif	// _DIVSEQ_H_
