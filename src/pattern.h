/*****************************************************************************
*
*	Header file for pattern search
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

#ifndef _PATTERN_
#define _PATTERN_

#define	    MAXPATLEN	80
#define	    MAXPATSIZ	MAXL
#define	    MAXLOC	2048
#define	    HITEI	0x80

typedef struct {CHAR* pat; short from, till;} VARPAT;

extern	char*	parspat(char* pat);
extern	int	findseq(int loc[], int maxloc, Seq* sd, char* pat);
extern	void	putloc(FILE* fd, Seq* sd, int* pl, int num);
extern	void	findpat(FILE* fd, Seq* sd, char* pat);

#endif
