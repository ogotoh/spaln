/*******************************************************************************
*
*	subsets
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

#ifndef  _SETS_H_
#define  _SETS_H_

static	const	int	EOTAB = -1;

struct Subset {
	int**	group;
	int*	pool;
	int	num;
	int	elms;
	Subset(Subset& ss);
	Subset(int n) {rsubset(n);}
	Subset(int n, int e) : num(n), elms(e) {
	    pool = new int[n + e];
	    group = new int*[n + 1];
	}
	Subset(int n, FILE* fd) {
	    rsubset(n, 0, fd);
	}
	Subset(int n, const char* str);
	~Subset();
	void	rsubset(int n, const char* ps = 0, FILE* fd = 0);
};

#endif
