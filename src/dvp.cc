/*****************************************************************************
*
*	Estimation of amino-acide sequence divergence
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

#include "aln.h"
#include "phyl.h"

static	void	setparam();
static	int	dvp(FILE* fd, Seq* seq, int* group1, int* group2);
extern	int	main(int argc, char** argv);

static void setparam()
{
	setdivseq(QUERY, SILENT, SILENT);
}

static int dvp(FILE* fd, Seq* seq, int* group1, int* group2)
{
	VSTAT	stat;
	static char* frmt = "  %-12s %-12s\n";

	divseq(&stat, seq, group1, group2);
	put_stat(fd, &stat);
	fprintf(fd, frmt, seq->lbl[*group1], seq->lbl[*group2]);
	return (OK);
}

int main(int argc, char** argv)
{
	int 	n;
	Seq* 	seq;

	initseq(&seq, 1);
	n = getoption(argc, argv);
	getargseq(argc - n, argv + n, &seq, getseq, 1);
	menu1(seq, inputseq, setparam, dvp);
	eraseq(seq);
	return (0);
}
