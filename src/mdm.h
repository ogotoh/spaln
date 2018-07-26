/*****************************************************************************
*
*	Header file for use of Mutation-Data-Matrix
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

#ifndef _MDM_H_
#define _MDM_H_

static	const	int	PAMSTEP = 10;
static	const	int	MAXPAM = 300;
static	const	INT	PAMLEVELS = MAXPAM / PAMSTEP;
static	const	INT	AASCMB = (AAS*(AAS+1)/2);
static	const	int	MAXPWR = 10;
static	const	int	GivenMat = -1;
static	const	int	EmptyMat = -2;
static	const	char	mdm_cmp[] = "mdm_cmp";
static	const	char	mdm_tab[] = "mdm_mtx";
static	const	char	mdmld_tab[] = "mdmld_mtx";
static	const	char	mdmpwr_tab[] = "mampwr_mtx";

#endif
