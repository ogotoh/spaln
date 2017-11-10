/*****************************************************************************
*
*	Lookup table for fast substring matching ( Wilber-Lipman algorithm)
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
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#ifndef _WLN_H_
#define _WLN_H_

#include "bitpat.h"

enum {ON_SCORE, ON_RDIAG, ON_POSIT};

class	mSeq;

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

struct WLUNIT {
	int	num;
	int	nid;
	int	tlen;
	int	llmt;
	int	ulmt;
	JUXT*	jxt;
	VTYPE	scr;
};

class Wilip {
	JUXT*	top;
	int	nwlu;
	WLUNIT*	wlu;
public:
	Wilip(Seq* seqs[], PwdB* pwd, INT level);
	Wilip(mSeq* seqs[], PwdB* pwd, INT level);
	~Wilip() {delete[] top; delete[] wlu;}
	int	size() {return nwlu;}
	WLUNIT*	begin() {return wlu;}
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
 	int	score;
 	int	mxscr;
	int	prevj;
	int	lastj;
	int	maxj;
};

class Wlp {
	Seq*	a;
	Seq*	b;
	int	bbt;
	Mfile*	mfd;
	PwdB*	pwd;
	WLPRM*	wlprm;
	Bitpat_wq*	bpp;
public:
	Wlp();
	Wlp(Seq* seqs[], PwdB* _pwd, INT level);
	~Wlp() {delete bpp; delete mfd;}
	INT*	lookup(INT* s, int kk);
	INT*	foldseq();
	VTYPE	eval(JUXT* jxt);
	JUXT*	reeval(JUXT* jxt, int* num);
	void	enter(JXTD* jxtd, int r);
	void	dmsnno(INT* m, INT jj, INT* t);
	void	dmsnno31(INT* m, INT jj3, INT* t);
	VTYPE	LinkHspScr(HSP* mcl, HSP* ncl);
	HSP*	mkhsps(JUXT* jxt, int n);
	WLUNIT* jxtcore(int* num, JUXT** ptop);
	WLUNIT* willip(int* nwlu, JUXT** ptop);
};

struct Wlprms {
	int	DvsP;
	VTYPE	Vab;
	Simmtx*	simmtx;
	VTYPE	EndBonus;
	VTYPE	RepPen;
	Wlprms(int dvsp);
	~Wlprms();
	VTYPE	sim2(CHAR* as, CHAR* bs) {return simmtx->mtx[*as][*bs];}
	void	initilize(INT a_molc, INT b_molc, INT level);
};

extern	void	makeWlprms(int dvsp);
extern	void	eraWlprms();
extern	WLPRM*	setwlprm(INT level);
extern	void	setexprm_x(const char* ps);
extern	JUXT*	revjxt(JUXT* jxt, int n);
extern	int	geneorient(Seq* seqs[], PwdB* pwd);

#endif
