#ifndef _BOYER_MOORE_
#define _BOYER_MOORE_

/**************************************************************************
*
*	Boyer-Moore algorithm for text search 
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

#include "seq.h"

class BoyerMoore {
const	CHAR*	text;
const	int	tlen;
	CHAR*	patt;
const	int	plen;
const	int	pl_1;
const	int	ll;
const	int	nalph;
const	int	step;
	int	sah;
	int*	dlt1;
	int*	dlt2;
	int	idx[3];
	int*	found;
const	int	nfound;
	PrQueue<int>*	pq;
	bool	eq(CHAR a, CHAR b)  {return (a == b);}
	bool	eqa(CHAR a, CHAR b) {return (a == b || a == AMB || b == AMB);}
	bool	eqn(CHAR a, CHAR b) {return ((a && b)? (--a & --b): false);}
	bool	eqta(CHAR a, CHAR b) {	// tron vs aa
	    return (a == b || (a == SER2 && b == SER) || b == AMB);
	}
	bool	eqtn(CHAR a, CHAR b) {	// tron vs nucl
	    a = aa2nuc[a];
	    return ((a && b)? (--a & --b): false);
	}
	bool	(BoyerMoore::*eqsrc)(CHAR a, CHAR b);
	bool	(BoyerMoore::*eqtab)(CHAR a, CHAR b);
	void	bmdelta1();
	void	bmdelta2();
public:
	int	nexthit(int l = -1, int r = -1);
	int	nexthit3(int l = -1, int r = -1);
	int	shift_after_hit() {return (sah);}
	int	size() {return (pq? pq->size(): 0);}
	bool	scanned(int n) {
	    if (size()) return (false);
	    n -= ll;
	    return ((step > 0)?
		min3(idx[0], idx[1], idx[2]) >= n:
		max3(idx[0], idx[1], idx[2]) <= n);
	}
	bool	finished() {
	    return ((step > 0)?
		min3(idx[0], idx[1], idx[2]) >= tlen:
		max3(idx[0], idx[1], idx[2]) <= 0);
	}
	BoyerMoore(const Seq* sd, const Seq* pd, int step_ = 1);
	BoyerMoore(const char* t, const char* p, int dir = 1);
	~BoyerMoore() {
	    delete[] patt; delete[] dlt1; delete[] dlt2; 
	    delete[] found; delete pq;
	}
};

#endif	// _BOYER_MOORE_
