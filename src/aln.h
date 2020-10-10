/*****************************************************************************
*
*	Collection of headers commonly used for sequence comparison
*
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
*	Osamu Gotoh, Ph.D. (2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef	_ALBH_
#define	_ALBH_

#include "seq.h"
#include "simmtx.h"

static	const	int	NOL = 3;
static	const	int	NOD = 2 * NOL - 1;
static	const	int	INTR = 2;
static	const	int	NCANDS = INTR * NOL;
static	const	int	MAXR = 4;
static	const	int	MAX_COLONY = 512;
static	const	int	Def_COLONY = 16;
static	const	int	MAX_PARALOG = 4;
static	const	int	LARGEN = (INT_MAX / 4 * 3);
static	const	VTYPE	SKIP = INT_MIN;
static	const	char	AlnParam[] ="AlnParam";
static	const	float	def_alprm2z = 2.;
static	const	int	end_of_ulk = INT_MIN;
static	const	int	dir2nod[16] = 
	{-1, -1, 0, 0, 2, 2, 2, 4, 1, 1, 1, 3, 2, 1, -1, -1};
static	const	float	defSss[2] = {0.3, 0.50};
static	const	int	N_Out_Modes = 16;

enum TraceBackDir {
	DEAD, RSRV, DIAG, NEWD, VERT, SLA1, SLA2, VERL, 
	HORI, HOR1, HOR2, HORL, NEWV, NEWH,
	SPIN=16, SPJC=32, SPF2=64, SPJCI=SPJC + SPIN,
	SPALL=SPIN+SPJC+SPF2
};

enum DistMes {QFdiv, QFidn, QCdiv, QCidn, QJuCa, Qpamd, NJuCa};
enum OutFm {GFF_FORM, ALN_FORM, PWA_FORM, BED_FORM, EXN_FORM, ITN_FORM, 
	CDS_FORM, AAS_FORM, CIG_FORM, VLG_FORM, SAM_FORM, BAM_FORM,
	BIN_FORM, MAP1_FORM, MAP2_FORM, PSJ_FORM};

//				  0,0,d,n,v,v,v,v,h,h,h,h,v,h,0,0
static const bool _is_diag[16] = {0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0};
static const bool _is_vert[16] = {0,0,0,0,1,1,1,1,0,0,0,0,1,0,0,0};
static const bool _is_hori[16] = {0,0,0,0,0,0,0,0,1,1,1,1,0,1,0,0};

#define isdiag(x) (_is_diag[(x)->dir & 15])
#define isvert(x) (_is_vert[(x)->dir & 15])
#define ishori(x) (_is_hori[(x)->dir & 15])
#define isntdiag(x) (!(_is_diag[(x)->dir & 15]))
#define isntvert(x) (!(_is_vert[(x)->dir & 15]))
#define isnthori(x) (!(_is_hori[(x)->dir & 15]))

enum ALN_MODE {BAD_ALN,
	NGP_ALN, HLF_ALN, RHF_ALN, GPF_ALN, NTV_ALN,
	NGP_ALB, HLF_ALB, RHF_ALB, GPF_ALB, NTV_ALB,
	NGP_ALS, HLF_ALS, RHF_ALS, GPF_ALS, NTV_ALS,
	NGP_ALH, HLF_ALH, RHF_ALH, GPF_ALH, NTV_ALH
};
enum Skl3M {AlgnTrb = 0x1, JuncTrb = 0x2, STEP3m = 0x4, STEP3n = 0x8, A_RevCom = 0x10};
enum {GLOBAL, LOCAL};

struct	ISLAND {VTYPE val; long upr, lwr;};
struct	BPPRM {float factor, g_alpha, g_beta; int maxb3d;};

struct BOUND {int la, lb, ua, ub;};

struct RVP {
	VTYPE	val;
	long	ptr;
};

struct RVPD {
	VTYPE	val;
	long	ptr;
	int	dir;
};

struct RVDWL {
	VTYPE	val;
	int	dir;
	int	lwr;
	int	upr;
	int	lst;
};

struct RVDJ {
	VTYPE	val;
	int	dir;
	int	jnc;
};

struct RVPDJ {
	VTYPE	val;
	long	ptr;
	int	dir;
	int	jnc;
	void reseth();
};

struct RVWU {
	VTYPE	val;
	int	upr;
	int	lwr;
	int	ulk;
};

struct RVWUJ {
	VTYPE	val;
	int	upr;
	int	lwr;
	int	ulk;
	int	jnc;
};

struct RVDWU {
	VTYPE	val;
	int	dir;
	int	upr;
	int	lwr;
	int	ulk;
};

struct RVDWUJ {
	VTYPE	val;
	int	dir;
	int	upr;
	int	lwr;
	int	ulk;
	int	jnc;
};

struct COLONY {
	VTYPE	val;
	int	lwr;
	int	upr;
	int	mlb;
	int	nlb;
	int	mrb; 
	int	nrb;
	short	clno;
	short	mark;
};

struct RVDWC {
	VTYPE	val;
	int	dir;
	int	lwr;
	int	upr;
	int	mlb;
	int	nlb;
	COLONY* clny;
};

struct RVDWJC {
	VTYPE	val;
	int	dir;
	int	jnc;
	int	upr;
	int	lwr;
	int	mlb;
	int	nlb;
	COLONY* clny;
};

static	const	RVP	black_vp = {NEVSEL, 0L};
static	const	RVP	white_vp = {0, 0L};
static	const	RVPD	black_vpd  = {NEVSEL, 0L, 0};
static	const	RVDWL	black_vdwl = {NEVSEL, 0, INT_MIN, INT_MAX, 0};
static	const	RVDWC	black_vdwc = {NEVSEL, 0, INT_MIN, INT_MAX, 0, 0};
static	const	RVDWC	white_vdwc = {0, 0, INT_MAX, INT_MIN, 0, 0, 0};
static	const	RVDJ	black_vdj   = {NEVSEL, 0, 0};
static	const	RVPDJ	black_vpdj = {NEVSEL, 0L, 0, 0};
static	const	RVWU	black_vwu  = {NEVSEL, INT_MIN, INT_MAX, end_of_ulk};
static	const	RVWUJ	black_vwuj = {NEVSEL, INT_MIN, INT_MAX, end_of_ulk, 0};
static	const	RVDWU	black_vdwu = {NEVSEL, 0, INT_MIN, INT_MAX, end_of_ulk};
static	const	RVDWUJ	black_vdwuj = {NEVSEL, 0, INT_MIN, INT_MAX, end_of_ulk, 0};
static	const	RVDWJC	black_vdwjc = {NEVSEL, 0, 0, INT_MIN, INT_MAX, 0, 0, 0};

inline void RVPDJ::reseth() {*this = black_vpdj;}

class Colonies {
protected:
	int	no_clny;
	COLONY*	clny;
public:
	Colonies(int n = 0);
	~Colonies() {delete[] clny;}
	int	size() const {return no_clny;}
	void	sortcolonies();
	void	detectoverlap(COLONY* cc);
	void	removeoverlap();
	void	removelowscore();
	COLONY*	at(int n = 0) const {return (clny + n);}
	COLONY*	end() const {return (clny + no_clny);}
	COLONY*	next() {
	    if (no_clny < (int) OutPrm.MaxOut) ++no_clny;
	    COLONY* cix = clny + no_clny;
	    vclear(cix);
	    cix->clno = no_clny;
	    return (cix);
	}
	int	index(COLONY* cc) const {return (cc - clny);}
};

/*********
	x: frameshift; y: splice signal; z: coding potential; o: termination codon; 
	bti: bias to trnsini; bsp: bias to splice sig; jneibr: neighborhood to ei junc.
**********/

struct PwdB {
const	Simmtx*	simmtx;
	int	DvsP;
	int	Noll;
	int	Nrow;
	int	Nrwb;
	VTYPE	Vab;	// a_thk * b_thk
	VTYPE	Vthr;	// Vab * thr
	VTYPE	BasicGOP;	// Vab * v
	VTYPE	BasicGEP;	// Vab * u
	VTYPE	LongGEP;	// Vab * u1
	VTYPE	diffu;	// Vab * (u - u1)
	int	codonk1;
	VTYPE	LongGOP;	// GasicGOP + diffu * k1
	VTYPE	GOP[NOL];
	int	MaxGapL;
	VTYPE	ExtraGOP;
	VTYPE	GapE1;
	VTYPE	GapE2;
	VTYPE	GapW1;
	VTYPE	GapW2;
	VTYPE	GapW3;
	VTYPE	GapW3L;
	Premat*	pmt;
	CodePot* codepot;
	ExinPot* exinpot;
	EijPat*	eijpat;

	PwdB(const Seq** seqs, const ALPRM* alp = 0);
	~PwdB();
	bool	same(const CHAR& ar, const CHAR& br) const {
	  switch (DvsP) {
	    case TxP: return(ar == br || (ar == SER2 && br == SER));
	    case PxT: return(ar == br || (br == SER2 && ar == SER));
	    default: return (ar == br);
	  }
	}
	VTYPE	sim2(const CHAR* as, const CHAR* bs) const {return simmtx->mtx[*as][*bs];}
	IntronPenalty*	IntPen;
	VTYPE	GapPenalty(int i) const {
	    if (i == 0) return 0;
	    return	(i > codonk1)? LongGOP + i * LongGEP:
		BasicGOP + i * BasicGEP;
	}
	VTYPE	GapExtPen(int i) const {
	    return (i > codonk1)? LongGEP: BasicGEP;
	}
	VTYPE	UnpPenalty(int d) const {
	    VTYPE	unp = d * BasicGEP;
	    if (d <= codonk1) return unp;
	    return (unp + diffu * (d - codonk1));
	}
	VTYPE	GapPenalty3(int i) const {return GapPenalty3(i, BasicGOP);}
	VTYPE	GapPenalty3(int i, VTYPE bgop) const;
	VTYPE	UnpPenalty3(int i) const {
	    int	d = i / 3;
	    VTYPE	unp = d * BasicGEP;
	    if (i <= codonk1) return unp;
	    return (unp - diffu * (d - alprm.k1) + (i % 3? ExtraGOP: 0));
	}
	VTYPE	GapExtPen3(int i) const {
	    return (i > codonk1)? LongGEP: BasicGEP;
	}
	VTYPE	TermGapExtPen3(int i) const {
	   return (i < alprm2.termk1)? 0: BasicGEP;
	}
};

/*	Select one or more output modes of aln and spaln */

class AlnOutModes {
	int	n_out_modes;
	int	out_mode[N_Out_Modes];
	char*	prefix;
public:
	FILE**	fds;
	void	getopt(const char* arg);
	void	setup(const char* prefix);
	void	alnoutput(Seq** sqs, Gsinfo* GsI);
	int	end() {return (n_out_modes);}
	AlnOutModes() :  n_out_modes(0), prefix(0), fds(0) {
	     vset(out_mode, -1, N_Out_Modes);
	}
	~AlnOutModes() {
	    delete[] prefix;
	    if (fds) {
		for (int n = 0; n < n_out_modes; ++n)
		    if (fds[n] && fds[n] != stdout) fclose(fds[n]);
		delete[] fds;
	    }
	}
};

/*	Headers to aln2.c	*/

extern	void	putvar(VTYPE x);
extern	VTYPE	selfAlnScr(const Seq* sd, const Simmtx* sm);
extern	int	prePwd(int molc, bool use_mdm = false);
extern	int	prePwd(const Seq* sd, bool use_mdm = false);
extern	int	prePwd(const Seq** seqs, bool use_mdm = false);
extern	void	stripe(const Seq* seqs[], WINDOW* wdw, int sh);
extern	void	stripe31(const Seq* seqs[], WINDOW* wdw, int shld);
extern	Seq*	synthseq(Seq* c, const Seq* a, const Seq* b, const SKL* skl);
extern	void	Intron53(Seq* sd, const PwdB* pwd, bool both_ori);
extern	void	Intron53N(Seq* sd, VTYPE f, const PwdB* pwd, bool both_ori);

extern	VTYPE	HomScoreB_ng(const Seq* seqs[], const PwdB* pwd);
extern	VTYPE	HomScoreH_ng(const Seq* seqs[], const PwdB* pwd);
extern	VTYPE	HomScoreS_ng(const Seq* seqs[], const PwdB* pwd);
extern	SKL*	alignB_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr);
extern	SKL*	alignH_ng(const Seq* seqs[], const PwdB* pwd, VTYPE* scr);
extern	SKL*	alignS_ng(Seq* seqs[], const PwdB* pwd, VTYPE* scr, int ori = 3);
extern	SKL*	nogap_skl(const Seq* a, const Seq* b = 0);
extern	VTYPE	skl_rngB_ng(const Seq* seqs[], Gsinfo* gsi, const PwdB* pwd);
extern	VTYPE	skl_rngH_ng(const Seq* seqs[], Gsinfo* gsi, const PwdB* pwd);
extern	VTYPE	skl_rngS_ng(const Seq* seqs[], Gsinfo* gsi, const PwdB* pwd);

extern	VTYPE	alnScoreD(const Seq* seqs[], const Simmtx* sm, int* ends = 0);
extern	FTYPE	alnscore2dist(Seq* sqs[], const PwdB* pwd, int* ends = 0, FTYPE denom = 0);

/*	Header to sqpr.c	*/

extern	int print2(Seq* seqs[], const GAPS** gps, double fscr, 
	Gsinfo* GsI, int nbr, int ttl, int skip, FILE* fd = 0);
extern	void	repalninf(Seq* seqs[], const Gsinfo* gsi, int mode = -1, FILE* fd = 0);
extern	void	printgene(Seq* seqs[], const Gsinfo* gsi, FILE* fd = 0);

extern	VTYPE	spSigII(const Seq* sd);
extern	BPPRM	bpprm;

#endif
