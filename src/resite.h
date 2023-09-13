/*****************************************************************************
*
*	Header file for restriction site
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
*****************************************************************************/

#ifndef _RESITE_H_
#define _RESITE_H_

#define RENSIZE	10
#define RERSIZE	15
#define RECSIZE 4
#define RESSIZE 30

typedef struct {
	char	name[RENSIZE+1];
	char	chr[RERSIZE+1];
	char	seq[2 * RERSIZE + 2];
	int	cut;
	int 	rct;
} RESITE;

extern	RESITE*	nextresite();
extern	RESITE*	firstresite();
extern	RESITE*	recogseq(char* src);
extern	int	respos(int loc[], int maxloc, Seq* ns, RESITE* res);
extern	void	resites(FILE* fd, Seq* ns, char* pat);

#endif
