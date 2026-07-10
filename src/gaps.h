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

class Seq;

struct Gaps {
	int	capa = 0;	// max number of GAPS elements
	int	elms = 0;	// current number of GAPS elements
	int	alen = 0;	// alignment length
	GAPS*	gaps = 0;	// 
	Gaps(const int& cap, GAPS* gps = 0)
	    : capa(cap), elms(gps? cap: 0) {
	    if (capa) gaps = gps? gps: new GAPS[capa];
	}
	Gaps(const int& left, const int& right)
	    : capa(2), elms(2), alen(right - left) {
	    gaps = new GAPS[capa];
	    gaps[0].gps = left; 
	    gaps[1].gps = right;
	    gaps[0].gln = gaps[1].gln = 0;
	}
	Gaps(const Gaps& src) 
	    : capa(src.capa), elms(src.elms), alen(src.alen) {
	    if (capa) {
		gaps = new GAPS[capa];
		vcopy(gaps, src.gaps, elms);
	    }
	}
	~Gaps() {delete[] gaps;}
	Gaps& operator=(const Gaps& src) {
	    elms = src.elms;
	    alen = src.alen;
	    if (capa < src.capa) {
		capa = src.capa;
		delete[] gaps;
		gaps = new GAPS[capa];
	    }
	    vcopy(gaps, src.gaps, elms);
	    return (*this);
	}
	int	capacity() const {return (capa);}
	int	size() const {return (elms);}
	GAPS*	begin() const {return (gaps);}
	GAPS*	end() const {return (gaps? (gaps + elms): 0);}
	GAPS*	last() const {return (gaps? (gaps + elms - 1): 0);}
	bool	gapless() const {
	    return (!elms || (elms == 1 && gaps->gln == 0) ||
		(elms == 2 && gaps->gln == 0 && gaps[1].gln == 0));
	}
	void	reset(const int& left = -1, const int& right = 0) {
	    if (capa < 2) resize(2, false);
	    if (left < 0) {
		elms = alen = 0;
		return;
	    }
	    gaps[0].gps = left; 
	    gaps[0].gln = 0;
	    alen = (left < right)? right - left: 0;
	    if (alen == 0) {
		elms = 1;
		return;
	    }
	    elms = 2;
	    gaps[1].gps = right;
	    gaps[1].gln = 0;
	}
	void	resize(const int& newcapa, const bool& keep = false) {
	    if (newcapa <= capa) return;
	    capa = newcapa;
	    GAPS*	ngaps = new GAPS[newcapa];
	    if (keep) vcopy(ngaps, gaps, elms);
	    delete[] gaps;
	    gaps = ngaps;
	}
	void	push_back(const GAPS& gp) {
	    if (elms >= capa) resize(2 * capa, true);
	    gaps[elms++] = gp;
	}
};

extern	void	putskl(const SKL* skl);
extern	void	swapskl(SKL* skl);
extern	bool	badskl(const SKL* skl, const Seq** sqs = 0);
extern	SKL*	stdskl(SKL** pskl);
extern	SKL*	stdskl3(SKL** pskl);
extern	int	sklpartner(const SKL* skl, int m, int given);
extern	SKL*	trimskl(const Seq* seqs[], SKL* skl);
extern	bool	sameskl(const SKL* a, const SKL* b);
extern	void	skl2gaps(Gaps** gaps, const SKL* skl);
extern	void	skl2gaps3(Gaps** gaps, const SKL* skl, const int& pro);
extern	void	toimage(Gaps** gaps, int numseq);
extern	void	unfoldgap(Gaps* gp, int step);
extern	SKL*	gap2skl(const Gaps* ga, const Gaps* gb, const Seq** sqs = 0);

#endif
