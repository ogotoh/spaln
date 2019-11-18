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
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef _WLN_H_
#define _WLN_H_

#include "bitpat.h"
#include <fcntl.h>
#include <unistd.h>

enum {ON_SCORE, ON_RDIAG, ON_POSIT};

static	const	char	strlut[] = "LookupTabs";
static	const	INT	lut_version = 10;
static	const	INT	MaxWlpLevel = 3;

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

/*****************************************************
	LookupTable
*****************************************************/

struct LookupTab {
	SHORT*	header;
	SHORT*	position;
};

struct LookupTabInfo {
	SHORT	version;
	short	dvsp;
	SHORT	elem;
	SHORT	tuple;
	SHORT	BitPat;
	SHORT	Nshift;
	INT	tabsize;
	INT	possize;
	INT	blklen;
	INT	n_lut;
	INT	reserve;
};

class MakeLookupTabs : public LookupTab, public WordTab {
	bool	master;
	LookupTabInfo	info;
	INT&	tab_size;
	INT&	pos_size;	// in table
	INT&	n_pos;		// = blklen
	INT&	n_lut;
	int	flut;		// save table
	int	flux;		// table index
	void lutreset() {
	    if (header) vclear(header, tab_size + 1);
	    if (position) vclear(position, pos_size);
	}
public:
	INT	prelude;
	MakeLookupTabs(const Seq* sq, const char* fname, INT bl, 
	    const WLPRM* wlp, INT afact, INT wq_size);
	MakeLookupTabs(const MakeLookupTabs& src);
	~MakeLookupTabs();
	void store();
};

class LookupTabs {
	char*	fname;
	LookupTabInfo	info;
	INT	unitsz;
	INT	unitno;
	size_t	lutno;
	INT&	tab_size;
	INT&	pos_size;
	INT&	n_pos;
	INT&	n_lut;
	bool	on_memory;
	INT	prvf;
	INT	prve;
	int	flut;   // stored table
	int	flux;	// table index
	SHORT*	lutbuf;
	LookupTab*	lutabs;
	long*	lutidx;
	long	basep;	// base file position
public:
	INT	blklen() {return(n_pos);}
	INT	nshift() {return(info.Nshift);}
	LookupTabs(const char* fn, int n_lut = 1);
	LookupTabs(LookupTabs& src);	// copy constructor
	~LookupTabs() {
	    if (flut >= 0) {
		close(flut);
		delete[] fname;
		delete[] lutbuf;
		delete[] lutabs;
	    }
	    if (flux >= 0) {
		close(flux);
		delete[] lutidx;
	    }
	}
	LookupTab*	readluts(INT from, INT end = 0);
	void	readheader() {
	    if (read(flut, &info, sizeof(LookupTabInfo)) != sizeof(LookupTabInfo))
		fatal(read_error, fname);
	    basep = lseek(flut, 0, SEEK_CUR);
	}
	void	readin(int lidx, int anlut);
	bool	error() {return (!lutabs);}
	void	reportinfo() {
	    printf("LUT:%s Version %u, elem %u, tuple %u, BitPat %u, Nshift %u\n", 
		fname, info.version, info.elem, info.tuple, info.BitPat, info.Nshift);
	    printf("TabSize %u, PosSize %u, UnitLen %u, No_Unit %u\n",
		tab_size, pos_size, n_pos, n_lut);
	}
	void	testlut(INT from, INT to);
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

class Wilip {
	JUXT*	top;
	int	nwlu;
	WLUNIT*	wlu;
public:
	Wilip(const Seq* seqs[], const PwdB* pwd, INT level);
	Wilip(Seq* seqs[], const PwdB* pwd, LookupTabs* lut, INT lb, INT rb);
	Wilip(mSeq* seqs[], const PwdB* pwd, INT level);
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
const 	Seq*	a;
const 	Seq*	b;
	int	bbt;
	Mfile*	mfd;
const 	PwdB*	pwd;
	WLPRM*	wlprm;
	Bitpat_wq*	bpp;
	int	awspan;	
	int	bwspan;
public:
	Wlp() : a(0), b(0), bbt(0), mfd(0), pwd(0), wlprm(0), bpp(0) {}
	Wlp(const Seq* seqs[], const PwdB* _pwd, INT level);
	~Wlp() {delete bpp; delete mfd;}
	INT*	lookup(INT* s, int kk);
	INT*	foldseq();
	VTYPE	eval(JUXT* jxt);
	JUXT*	reeval(JUXT* jxt, int& num);
	void	enter(JXTD* jxtd, int r, bool on_k = false);
	void	dmsnno(INT jj, INT* m, INT* t);
	void	dmsnno31(INT jj3, INT* m, INT* t);
	void	dmsnno(INT jj, LookupTabs* lut, INT lb, INT rb);
	void	revcoord(JUXT* jxt, const int& n);
	VTYPE	LinkHspScr(HSP* mcl, HSP* ncl);
	HSP*	mkhsps(const JUXT* jxt, int n);
	WLUNIT* jxtcore(int& num, JUXT** ptop);
	JUXT*	run_dmsnno(int& njxt, LookupTabs* lut = 0, INT lb = 0, INT rb = 0);
	WLUNIT* willip(JUXT** ptop, int& nwlu, JUXT* jxt);
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

extern	void	makeWlprms(int dvsp);
extern	void	eraWlprms();
extern	WLPRM*	setwlprm(INT level);
extern	WLPRM*	selectwlprm(INT sz, int dvsp, WLPRM* wlp = 0);
extern	void	setexprm_x(int& argc, const char**& argv);
extern	JUXT*	revjxt(JUXT* jxt, const int& n);
extern	int	geneorient(Seq* seqs[], const PwdB* pwd, int max_n = 0);

#endif
