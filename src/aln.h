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
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#ifndef	_ALBH_
#define	_ALBH_

#include "seq.h"
#include "simmtx.h"

#define	NOL	3
#define	NOD	(2 * NOL + 1)
#define	MAX_COLONY	512
#define	Def_COLONY	16
#define	MAX_PARALOG	4
#define	MARGIN(m, A)	((m) == (A)->left || (m) == (A)->right)
#define	AlnParam	"AlnParam"

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

#define	LARGEN	(INT_MAX / 4 * 3)
#define EACHSETSIM	0

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
enum Skl3M {AlgnTrb = 0x1, JuncTrb = 0x2, STEP3m = 0x4, STEP3n = 0x8};
enum {GLOBAL, LOCAL};

struct	ISLAND {VTYPE val; long upr, lwr;};
struct	BPPRM {float factor; int maxb3d;};

/*********
	x: frameshift; y: splice signal; z: coding potential; o: termination codon; 
	bti: bias to trnsini; bsp: bias to splice sig; jneibr: neighborhood to ei junc.
**********/

struct PwdB {
	Simmtx*	simmtx;
	int	DvsP;
	int	Noll;
	int	Nrow;
	int	Nrwb;
	int	MaxGapL;
	VTYPE	Vab;	// a_thk * b_thk
	VTYPE	Vthr;	// Vab * thr
	VTYPE	BasicGOP;	// Vab * v
	VTYPE	BasicGEP;	// Vab * u
	VTYPE	LongGEP;	// Vab * u1
	VTYPE	diffu;	// Vab * (u - u1)
	int	codonk1;
	VTYPE	LongGOP;	// GasicGOP + diffu * k1
	VTYPE	GOP[NOL];
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

	PwdB(Seq** seqs, ALPRM* alp = 0);
	~PwdB();
	bool	same(const CHAR& ar, const CHAR& br) {
	  switch (DvsP) {
	    case TxP: return(ar == br || (ar == SER2 && br == SER));
	    case PxT: return(ar == br || (br == SER2 && ar == SER));
	    default: return (ar == br);
	  }
	}
	VTYPE	sim2(CHAR* as, CHAR* bs) {return simmtx->mtx[*as][*bs];}
	IntronPenalty*	IntPen;
	VTYPE	GapPenalty(int i) {
	    if (i == 0) return 0;
	    return	(i > codonk1)? LongGOP + i * LongGEP:
		BasicGOP + i * BasicGEP;
	}
	VTYPE	GapExtPen(int i) {
	    return (i > codonk1)? LongGEP: BasicGEP;
	}
	VTYPE	UnpPenalty(int d) {
	    VTYPE	unp = d * BasicGEP;
	    if (d <= codonk1) return unp;
	    return (unp + diffu * (d - codonk1));
	}
	VTYPE	prematT(const CHAR* ps) {return prematT(ps);}
	VTYPE	GapPenalty3(int i) {return GapPenalty3(i, BasicGOP);}
	VTYPE	GapPenalty3(int i, VTYPE bgop);
	VTYPE	UnpPenalty3(int i) {
	    int	d = i / 3;
	    VTYPE	unp = d * BasicGEP;
	    if (i <= codonk1) return unp;
	    return (unp - diffu * (d - alprm.k1) + (i % 3? ExtraGOP: 0));
	}
	VTYPE	GapExtPen3(int i) {
	    return (i > codonk1)? LongGEP: BasicGEP;
	}
	VTYPE	TermGapExtPen3(int i) {
	   return (i < alprm2.termk1)? 0: BasicGEP;
	}
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

class Colonies {
protected:
	int	no_clny;
	COLONY*	clny;
public:
	Colonies(int n = 0);
	~Colonies() {delete[] clny;}
	int	size() {return no_clny;}
	void	sortcolonies();
	void	detectoverlap(COLONY* cc);
	void	removeoverlap();
	void	removelowscore();
	COLONY*	at(int n = 0) {return (clny + n);}
	COLONY*	end() {return (clny + no_clny);}
	COLONY*	next() {
	    if (no_clny < (int) OutPrm.MaxOut) ++no_clny;
	    COLONY* cix = clny + no_clny;
	    vclear(cix);
	    cix->clno = no_clny;
	    return (cix);
	}
	int	index(COLONY* cc) {return (cc - clny);}
};

extern	BPPRM	bpprm;
static	const	float	def_alprm2z = 2.;

/*	Headers to aln2.c	*/

extern	void	putvar(VTYPE x);
extern	VTYPE	selfAlnScr(Seq* sd, Simmtx* sm);
extern	int	prePwd(int molc, bool use_mdm = false);
extern	int	prePwd(Seq* sd, bool use_mdm = false);
extern	int	prePwd(Seq** seqs, bool use_mdm = false);
extern	void	stripe(Seq* seqs[], WINDOW* wdw, int sh);
extern	void	stripe31(Seq* seqs[], WINDOW* wdw, int shld);
extern	Seq*	synthseq(Seq* c, Seq* a, Seq* b, SKL* skl);
extern	void	Intron53(Seq* sd, PwdB* pwd);
extern	void	Intron53N(Seq* sd, VTYPE f, PwdB* pwd);

extern	SKL*	alignB_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr);
extern	SKL*	alignH_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr);
extern	SKL*	alignS_ng(Seq* seqs[], PwdB* pwd, VTYPE* scr);
extern	VTYPE	skl_rngB_ng(Seq* seqs[], Gsinfo* gsi, PwdB* pwd);
extern	VTYPE	skl_rngH_ng(Seq* seqs[], Gsinfo* gsi, PwdB* pwd);
extern	VTYPE	skl_rngS_ng(Seq* seqs[], Gsinfo* gsi, PwdB* pwd);

extern	VTYPE	alnScoreD(Seq* seqs[], Simmtx* sm, int* ends = 0);
extern	FTYPE	alnscore2dist(Seq* sqs[], PwdB* pwd, int* ends = 0);

/*	Header to sqpr.c	*/

extern	int print2(Seq* seqs[], GAPS* gps[], double fscr, 
	Gsinfo* GsI, int nbr, int ttl, int skip);
extern	void	repalninf(Seq* seqs[], Gsinfo* gsi);
extern	void	printgene(Seq* seqs[], Gsinfo* gsi);

extern	VTYPE	spSigII(Seq* sd);

#endif
