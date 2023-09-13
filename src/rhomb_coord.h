/*****************************************************************************
*
*	rhomboidic bitmap coordinate system
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

#ifndef _RHOMB_COORD_H_
#define _RHOMB_COORD_H_

/*--------------------------------------------
		********************
	       ********************

		............

	********************
--------------------------------------------*/

#include "cmn.h"

enum class TraceBackCode {
	STOP = 0,	// stop
	DIAG = 1,	// diagonal
	HORI = 2,	// hori
	HORL = 3,	// long hori
	HOR1 = 4,	// 1nt ins
	HOR2 = 5,	// 2nt ins
	VERT = 8,	// vert
	VERL = 9,	// long vert
	VER1 = 10,	// 1nt del
	VER2 = 11,	// 2nt del
	ACCM = 13,	// 3' end of exon phase -1
	ACCR = 14,	// 3' end of exon
	ACCZ = 14,	// 3' end of exon phase 0
	ACCP = 15,	// 3' end of exon phase +1

	NHOR = 16,	// new hori
	NVER = 32,	// new vert
	NHOL = 64,	// new long hori
	DONM = 64,	// 5' end of intron phase -1
	NVEL = 128,	// new long vert
	DONR = 128,	// 5' end of intron
	DONZ = 128,	// 5' end of intron phase 0
	DONP = 256,	// 5' end of intron phase +1
	NIND = 240	// NHOR + NVER + NHOL + NVEL
};

static const char*	bad_dir = "Unexpected dir %d\n";

template <typename var_t>
class Anti_rhomb_coord {
const	int	m_base;
const	int	n_base;
const	int	step;
const	int	m_width;
const	int	n_width;
	var_t*	bbuf;
	int	cur_m;
	int	cur_n;
	var_t*	cur_p;
public:
	Anti_rhomb_coord(const int mmax, const int nmax, 
	const int mbase = 0, const int nbase = 0, const int st = 1);
	~Anti_rhomb_coord()	{delete[] bbuf;}
	void	initialize_m0(const var_t dir);
	var_t*	set_point(const int& m, const int& n) {
	    cur_m = m - m_base; 
	    cur_n = n - n_base;
	    cur_p = bbuf + (step * cur_m + cur_n) * m_width + cur_m;
	    return (cur_p);
	}
	void	to_right(var_t*& p) {p += m_width;}
	INT	to_left(int& m, int& n, const int s = 1) {
		m = cur_m;
		n = cur_n -= s;
		if (n < 0) {
		    cur_n = n = 0;
		    return (0);
		}
		cur_p -= s * m_width;
		return (*cur_p);
	}
	INT	to_upper(int& m, int& n, const int s) {
		m = --cur_m;
		n = cur_n -= s;
		if (m < 0) {
		    cur_m = m = 0;
		    cur_n = n += s;
		    return (0);
		} else if (n < 0) {
		    cur_n = n = 0;
		    return (0);
		}
		cur_p -= ((step + s) * m_width + 1);
		return (*cur_p);
	}
	INT	go_back(INT dir, int& m, int& n);
	void	traceback(int m, int n, Mfile* mfd);
};

template <typename var_t>
Anti_rhomb_coord<var_t>::Anti_rhomb_coord(const int mmax, const int nmax, 
const int mbase, const int nbase, const int st) :
	m_base(mbase), n_base(nbase), step(st), 
	m_width(mmax - mbase + 1),
	n_width(nmax - nbase + 1 + step * m_width)
{
	size_t	bufsize = m_width * n_width + 32;	// play 
	try {
	    bbuf = new var_t[bufsize];
	} catch (std::bad_alloc& ba) {
	    fatal("no more memory for Anti_rhomb_coord");
	}
	vclear(bbuf, bufsize);
}

template <typename var_t>
void Anti_rhomb_coord<var_t>::initialize_m0(const var_t dir)
{
	var_t*	p = bbuf;
	for (int n = 1; n < n_width; ++n)
	    *(p += m_width) = dir;
}

template <typename var_t>
INT Anti_rhomb_coord<var_t>::go_back(INT code, int& m, int& n)
{
	TraceBackCode	dir = static_cast<TraceBackCode>(code & 15);

	switch (dir) {
	  case TraceBackCode::STOP:
	    break;
	  case TraceBackCode::DIAG:
	    do {
		if (!(code = to_upper(m, n, step))) return (0);
	    } while (static_cast<TraceBackCode>(code & 15) == TraceBackCode::DIAG);
	    break;
	  case TraceBackCode::HORI:
	    while (!(code & static_cast<int>(TraceBackCode::NHOR))) {
		if (!(code = to_left(m, n, step))) return (0);
	    }
	    dir = static_cast<TraceBackCode>(code & 15);
	    if (dir != TraceBackCode::HOR1 && dir != TraceBackCode::HOR2)
		code = to_left(m, n, step);
	    break;
	  case TraceBackCode::HORL:
	    while (!(code & static_cast<int>(TraceBackCode::NHOL))) {
		if (!(code = to_left(m, n, step))) return (0);
	    }
	    dir = static_cast<TraceBackCode>(code & 15);
	    if (dir != TraceBackCode::HOR1 && dir != TraceBackCode::HOR2)
		code = to_left(m, n, step);
	    break;
	  case TraceBackCode::VERT:
	    while (!(code & static_cast<int>(TraceBackCode::NVER))) {
		if (!(code = to_upper(m, n, 0))) return (0);
	    }
	    dir = static_cast<TraceBackCode>(code & 15);
	    if (dir != TraceBackCode::VER1 && dir != TraceBackCode::VER2)
		code = to_upper(m, n, 0);
	    break;
	  case TraceBackCode::VERL: 
	    while (!(code & static_cast<int>(TraceBackCode::NVEL))) {
		if (!(code = to_upper(m, n, 0))) return (0);
	    }
	    dir = static_cast<TraceBackCode>(code & 15);
	    if (dir != TraceBackCode::VER1 && dir != TraceBackCode::VER2)
		code = to_upper(m, n, 0);
	    break;
	  case TraceBackCode::ACCR:
	    do {
		if (!(code = to_left(m, n, 1))) return (0);
	    } while (!(code & static_cast<int>(TraceBackCode::DONR)));
	    break;
	  case TraceBackCode::ACCM:
	    do {
		if (!(code = to_left(m, n, 1))) return (0);
	    } while (!(code & static_cast<int>(TraceBackCode::DONM)));
	    break;
	  case TraceBackCode::ACCP:
	    do {
		if (!(code = to_left(m, n, 1))) return (0);
	    } while (!(code & static_cast<int>(TraceBackCode::DONP)));
	    code = to_upper(m, n, step);
	    ++m; n += step;
	    break;
	  case TraceBackCode::HOR1:
	    code = to_left(m, n, 1);
	    break;
	  case TraceBackCode::HOR2:
	    code = to_left(m, n, 2);
	    break;
	  case TraceBackCode::VER1:
	    code = to_upper(m, n, 1);
	    break;
	  case TraceBackCode::VER2:
	    code = to_upper(m, n, 2);
	    break;
	  default:
	    fatal(bad_dir, code);
	    break;
	}
	return (code);
}

template <typename var_t>
void	Anti_rhomb_coord<var_t>::traceback(int m, int n, Mfile* mfd)
{
	INT code = *set_point(m, n);
	m -= m_base;
	n -= n_base;
	while (code) {
	    SKL	wsk = {m + m_base, n + n_base};
	    mfd->write(&wsk);
	    code = go_back(code, m, n);
	}
	SKL	wsk = {m + m_base, n + n_base};
	mfd->write(&wsk);
}

#endif
