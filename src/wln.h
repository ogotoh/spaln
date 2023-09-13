/*****************************************************************************
*
*	Lookup table for fast substring matching ( Wilber-Lipman algorithm)
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

#ifndef _WLN_H_
#define _WLN_H_

#include "bitpat.h"
#include <fcntl.h>
#include <unistd.h>

enum {ON_SCORE, ON_RDIAG, ON_POSIT};

static	const	INT	MaxWlpLevel = 3;
static	const	int	WlnPamNo = 2;
static	int	shortquery = 50;

class	mSeq;

using	KV_r_nh = KVpair_lt<int, int>;

struct WLPRM {
const	char*	bitpat;
const	char*	redpat;
	INT	elem;
	INT	tpl;
	INT	mask;
	INT	width;
	INT	gain;
	INT	gain1;
	INT	thr;
	int	xdrp;
	int	cutoff;
	VTYPE	vthr;
	INT*	ConvTab;
};

/*****************************************************
	Wilip
*****************************************************/

struct WLUNIT {
	int	num;
	int	nid;
	int	tlen;
	int	llmt;
	int	ulmt;
	JUXT*	jxt;
	VTYPE	scr;
};

class Wlp;

class Wilip {
	JUXT*	top = 0;
	int	nwlu = 0;
	WLUNIT*	wlu = 0;
	bool	int_wlp = true;
public:
	Wilip(const Seq* seqs[], const PwdB* pwd, const int level);
	Wilip(const Seq* seqs[], Wlp* wln);
	Wilip(mSeq* seqs[], const PwdB* pwd, const int level);
	~Wilip() {delete[] top; if (int_wlp) delete[] wlu;}
	int	size() const {return nwlu;}
	WLUNIT*	begin() {return wlu;}
	WLUNIT*	end() {return (wlu? wlu + nwlu: 0);}
	void	sort_on_scr();
	void	shift_y(int bias, int rbias);
};

struct JXT_SCORE {
	int	num;
	VTYPE	max;
	VTYPE	sum;
	RANGE	rx;
	RANGE	ry;
	WINDOW	wdw;
};

class JxtQueue {
	JUXT**	jqueue;
	JUXT**	kqueue;
	int	nrep;
	int	qp, kq;
public:
	JxtQueue(int n) : nrep(n), qp(0) {
	    jqueue = new JUXT*[nrep];
	    kqueue = new JUXT*[nrep];
	    vclear(jqueue, nrep);
	}
	~JxtQueue()	{delete[] jqueue; delete[] kqueue;}
	int	push(JUXT* wjx);
	JUXT*	pop(bool last = false) {
	    if (last) {
		if (!qp) return (0);
		while (qp && !jqueue[--qp]);
		return (jqueue[qp]);
	    } else	return (kq? kqueue[--kq]: 0);
	}
};

struct HSP {
	int	lx;
	int	ly;
	int	rx;
	int	ry;
	int	ux;
	int	rr;
	int	nid;
	int	len;
	int	irno;
	VTYPE	jscr;
	VTYPE	sscr;
	VTYPE	sumh;
	HSP*	ulnk;
};	/* high-scoring pair */

struct HSPPRM {
	int	DirRep;	/* Max repeat length on a splicing boundary */
	float	RepPen;	/* Penalty for back-going direction */
};

struct JXTD {
union {
 	int	score;
	int	nhit;
};
union {
	int	mxscr;
	int	ml;
};
union {
	int	prevj;
	int	mr;
};
	int	lastj;
	int	maxj;
};

class Wlp {
protected:
const 	Seq*	a = 0;
const 	Seq*	b = 0;
const 	PwdB*	pwd = 0;
const	int	bbt = 0;	// 1 or 3
const	int	mm = 0;		// query length
const	int	sect_l = 0;	// section length
	Mfile*	mfd = 0;
const	WLPRM*	wlprm = 0;
	WLPRM	tmprm;		// temporaly
const	int	tplwt = 0;
const	int	awspan = 0;
const	int	bwspan = 0;
	int	precutoff = 0;
const	JXTD	ixtd;
const	int	min_lnkscr;
const	bool	dhit = false;
	INT*	header = 0;
	INT	ttl = 0;
	INT*	position = 0;
	JXTD*	jxtd = 0;
	Bitpat_wq*	bpp = 0;
	INT*	lookup(INT* s, int kk);
virtual	INT*	foldseq();
	INT*	kmercount();
	RANGE*	lowic(const INT* cnt);
	VTYPE	eval(JUXT* jxt);
	JUXT*	reeval(JUXT* jxt, int& num);
	void	enter(const JXTD* wxtd, const int& r, const bool& on_k = false);
	void	storedh(const int& r, int ml, int mr);
	void	scan_b(INT m, INT n);
virtual	void	dmsnno();
virtual	void	dmsnno31();
	VTYPE	LinkHspScr(HSP* mcl, HSP* ncl);
	HSP*	mkhsps(const JUXT* jxt, int& n);
	WLUNIT* jxtcore(int& num, JUXT** ptop);

public:
	Wlp(const Seq* seqs[], const PwdB* _pwd, const int level);
	Wlp(Seq* seqs[]);
virtual	~Wlp() {
	    delete bpp; delete mfd; delete[] header; 
	    delete[] position; delete[] jxtd;
	}
	bool	ng() const {return (!position);}
	void	reset(const Seq* g) {b = g; if (mfd) mfd->reset();}
	JUXT*	run_dmsnno(int& njxt);
	WLUNIT* willip(JUXT** ptop, int& nwlu, JUXT* jxt);
friend	class	Wilip;
};

struct Wlprms {
	int	DvsP;
	VTYPE	Vab;
	Simmtx*	simmtx;
	VTYPE	EndBonus;
	VTYPE	RepPen;
	Wlprms(int dvsp);
	~Wlprms();
	VTYPE	sim2(const CHAR* as, const CHAR* bs) {return simmtx->mtx[*as][*bs];}
	void	initilize(INT level);
};

inline	void	set_shortquery(const int s) {shortquery = s;}
extern	void	makeWlprms(int dvsp);
extern	void	eraWlprms();
extern	WLPRM*	setwlprm(INT level);
extern	void	setexprm_x(int& argc, const char**& argv);
extern	int	geneorient(Seq* seqs[], const PwdB* pwd);

#endif	// _WLN_H_
