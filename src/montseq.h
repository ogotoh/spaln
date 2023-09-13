/**********************************************************************
*
*	Generate and mutate random nucleotide/amino acid
*	sequences for Monte Carlo Test
*
**	Osamu Gotoh, ph.D.	(-2001)
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

extern	double	rgauss();
extern	int	rpoisson(double mu);
extern	void	initrandseq(double acgt[], int molc);
extern	Seq*	randseq(Seq* sd, int len);
extern	Seq*	substseq(Seq* sd, int n, int which);
