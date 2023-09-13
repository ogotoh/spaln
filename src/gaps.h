/*****************************************************************************
*
*	Header for abstract structures of alingmnet
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

#ifndef _GAPS_H_
#define _GAPS_H_

static	const	int	gaps_end = EOS;

inline	int	gaps_span(const GAPS* gg) {return gg->gps;}
inline	int	gaps_size(const GAPS* gg) {return gg->gln;}
inline	GAPS*	gaps_begin(GAPS* gg) {return gg + 1;}
inline	const	GAPS*	gaps_begin(const GAPS* gg) {return gg + 1;}
inline	bool	gaps_free(const GAPS* gg) {return gg->gln == 3;}
inline	bool	gaps_intr(const GAPS* g) {return (g->gln != gaps_end);}

extern	void	putskl(const SKL* skl);
extern	void	toimage(GAPS* gaps[], int numseq);
extern	void	unfoldgap(GAPS* gp, int step, bool hl = false);
extern	void	skl2gaps(GAPS* gaps[], const SKL* skl, bool hl = false);
extern	void	skl2gaps3(GAPS* gaps[], const SKL* skl, int pro);
extern	void	swapskl(SKL* skl);
extern	bool	badskl(const SKL* skl, const Seq** sqs = 0);
extern	SKL*	stdskl(SKL** pskl);
extern	SKL*	stdskl3(SKL** pskl);
extern	SKL*	gap2skl(const GAPS* ga, const GAPS* gb);
extern	int	sklpartner(const SKL* skl, int m, int given);
extern	SKL*	trimskl(const Seq* seqs[], SKL* skl);
extern	bool	sameskl(const SKL* a, const SKL* b);

#endif
