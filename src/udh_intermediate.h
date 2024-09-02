/*****************************************************************************
*
*	Define 'Intermediate' structure used in
*	unidirectional Hirschberg linear-space algorithm
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
*	Osamu Gotoh, Ph.D.	(2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
*
*****************************************************************************/

#ifndef	_Udh_intermediate_H_
#define	_Udh_intermediate_H_

#include "aln.h"

class UdhIntermediate {
	int*	buf = 0;
public:
const	int	mi;		// intermediate m-coordinate
	int**	hlnk = 0;	// horizontal link
	int**	vlnk = 0;	// vertical link
	int**	lwrb = 0;	// lower bound diagoonal
	int**	uprb = 0;	// upper bound diagoonal
	UdhIntermediate(const int& m, const WINDOW& wdw, 
		const int& nol_, const bool& lub = false) 
	    : mi (m) {
const	    int	nol = std::min(NOL, nol_);
	    hlnk = new int*[(lub? 4: 2) * nol];
	    vlnk = hlnk + nol;
const	    size_t	u_size = nol * wdw.width;
const	    size_t	buf_siz = (lub? 4: 2) * u_size;
	    buf = new int[buf_siz];
	    vset(buf, end_of_ulk, 2 * u_size);
	    hlnk[0] = buf - wdw.lw + 1;
	    vlnk[0] = hlnk[0] + u_size;
	    if (lub) {
		lwrb = vlnk + nol;
		uprb = lwrb + nol;
		lwrb[0] = vlnk[0] + u_size;
		uprb[0] = lwrb[0] + u_size;
		vset(lwrb[0] + wdw.lw - 1, INT_MAX, u_size);
		vset(uprb[0] + wdw.lw - 1, INT_MIN, u_size);
	    }
	    for (int k = 1; k < nol; ++k) {
		hlnk[k] = hlnk[k - 1] + wdw.width;
		vlnk[k] = vlnk[k - 1] + wdw.width;
		if (lub) {
		    lwrb[k] = lwrb[k - 1] + wdw.width;
		    uprb[k] = uprb[k - 1] + wdw.width;
		}
	    }
	}
	~UdhIntermediate() {delete[] hlnk; delete[] buf;}
};

class Udh_Imds {
const	int	n_imd;
	UdhIntermediate**	imds = 0;
public:
	Udh_Imds(const int& n, int mi, const int& mm, const WINDOW& wdw, 
	    const int& nol = NOL, const bool& lub = false)
	 : n_imd(n) {
	    imds = new UdhIntermediate*[n_imd];
	    for (int i = 0; i < n_imd; ++i)
		imds[i] = new UdhIntermediate(mi += mm, wdw, nol, lub);
	}
	~Udh_Imds() {
	    for (int i = 0; i < n_imd; )
		delete imds[i++];
	    delete[] imds;
	}
	UdhIntermediate*	operator[](const int& i) const {
	    return ((0 <= i && i < n_imd)? imds[i]: 0);
	}
};

using	Dim10 = int[10];

#endif	// _Udh_intermediate_H_
