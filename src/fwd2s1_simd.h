/*****************************************************************************
*
*	Vetorized version of DP-based spliced alignment of cDNA/EST vs DNA 
*	Accelerated by SIMD instructions of Intel intrinsics.
*	Curently accepted architectures: SSE4.1, AVX, AVX2, and AVX-512BW
*	Assumes capability of store/load in/from unaligned memory
*
*	5' and 3' splice site signals, intron-length distribution,
*	coding potential, and conservation of intron sites  are considered.
*	Assumes there is no internal gap in the reference cDNA/EST sequence.
*
*	Implement unidirectional Hirschberg linear-space algorithm
*	togegher with ordinary forward-traceback algorithm
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

#ifndef	_FWD2S1_SIMD_H_
#define _FWD2S1_SIMD_H_

#include "aln.h"
#include "vmf.h"

#include "simd_functions.h"
#include "udh_intermediate.h"

using	var_t = short;

static	const	int	check_scr = int(0.9 * SHRT_MAX);
static	const	int	PS_A = 5;
static	const	int	PV_A = 6;
static	const	int	nevsel_play = 1024;

struct Rvulmn {		// used in Local mode
	VTYPE	val;	// score
	int	ulk;	// lower end diagonal
	int	ml;	// left end m-coord
	int	mr;	// right end m-coord
	int	nr;	// right end n-coord
};

struct Rvlujd {
	var_t	val;	// score
	int	ulk;	// lower end diagonal
	int	jnc;	// junction
	int	ml;	// left end m-coord
	short	dir;	// direction
};

class SimdAln2s1 : public Simd_functions<var_t> {
#if __ARM_NEON
Simd_functions<SHORT>	usf;
#endif
protected:
#if __AVX512BW__
using	regist_v = __m512i;
using	regist_m = __mmask32;
#elif __AVX2__
using	regist_v = __m256i;
using	regist_m = var_v;
#elif __SSE4_1__
using	regist_v = __m128i;
using	regist_m = var_v;
#elif __ARM_NEON
using	regist_v = int16x8_t;
using	regist_m = uint16x8_t;
#endif
const	int	nelem = Nelem;
const	Seq**	seqs;
const	Seq*&	a;	// query
const	Seq*&	b;	// genomic
const	PwdB*	pwd;
const	WINDOW&	wdw;
const	bool	Local;
	SpJunc*	spjcs;
const	Cip_score* cip;	// bonus to conserved splicing junction
const	int	mode;	// 
const	bool	used;	// 2 or 4 bytes
const	bool	spj;	// splicing
const	bool	dagp;	// double affine
const	int	nod;
const	int	Ncand;
	Vmf*	vmf;
	int	mm = 0;
	UdhIntermediate*	imd = 0;
const	int	Np1;
const	int	avmch;
const	var_t	nevsel;
const	Rvulmn	black_Rvulmn;
	int	rlst;
	size_t	buf_size;
	var_t*	abuf = 0;
	var_t*	ps_a = 0;
	var_t*	pv_a = 0;
	var_t*	hv_a[2] = {0, 0};
	var_t*	fv_a = 0, *fv2_a = 0;
	var_t*	ev_a = 0, *ev2_a = 0;
	var_t*	s5_a = 0, *s3_a = 0;
	var_t*	hb_a[2] = {0, 0};
	var_t*	fb_a = 0, *fb2_a = 0;
	var_t*	eb_a = 0, *eb2_a = 0;
	var_t*	pb_a = 0;
	var_t*	qb_a = 0;
	var_t*	hc_a[2] = {0, 0};
	var_t*	fc_a = 0, *fc2_a = 0;
	var_t*	ec_a = 0, *ec2_a = 0;
	var_t*	hd_a[2] = {0, 0};
	var_t*	fd_a = 0, *fd2_a = 0;
	var_t*	ed_a = 0, *ed2_a = 0;
 	var_t*	vbuf = 0;
 	var_t*	bbuf = 0;
 	var_t*	cbuf = 0;
 	var_t*	dbuf = 0;
	var_t*	hv = 0, *fv = 0, *fv2 = 0;	// score
	var_t*	hb = 0, *fb = 0, *fb2 = 0;	// boundary
	var_t*	hc = 0, *fc = 0, *fc2 = 0;	// coordinate
	var_t*	hd = 0, *fd = 0, *fd2 = 0;	// coordinate2
	var_t*	hfesv[2][Nelem][7];
	var_t*	hfesb[2][Nelem][5];
	var_t*	hfesc[2][Nelem][5];
	var_t*	hfesd[2][Nelem][5];

class Sjsites {			// nested class
	SimdAln2s1& hs1;
	Rvlujd*	rcd;
	short*	idx;
	short	ncand;
	Rvlujd	black_r;
public:
	void	reset() {
	    ncand = -1;
	    for (int i = 0; i <= hs1.Ncand; ++i) {
		rcd[i] = black_r;
		idx[i] = i;
	    }
	}
	Sjsites(SimdAln2s1& _hs1) : hs1(_hs1) {
	    black_r.val = hs1.nevsel;
	    black_r.ulk = black_r.jnc = black_r.ml = black_r.dir = 0;
	    rcd = new Rvlujd[hs1.Ncand + 1];
	    idx = new short[hs1.Ncand + 1];
	}
	~Sjsites() {delete[] rcd; delete[] idx;}
	void	put(const int& j, const int& m, const int& n, const int& p);
	int	get(const int& j, const int& m, const int& n, const int& p);
};	// end of Sjsites

	Sjsites**	fsjss = 0;
	Queue2<int>*	donor_q = 0;
	Queue2<int>*	accep_q = 0;

	void	fhinitS1();
	int	fhlastS1(Rvulmn& maxh);
	int	from_spj(Queue2<int>& l, const int& m, 
		const int& n, const int& q);
	void	to_spj(Queue2<int>& l, const int& m, 
		const int& n, const int& q);
	void	vadd(var_t* ar, const int& n, const var_t& c) {
	    for (var_t* er = ar + n; ar < er; ++ar) {
		VTYPE	x = *ar + c;
		*ar = (x < nevsel)? nevsel: static_cast<var_t>(x);
	    }
	}
	int	checkpoint(const var_t& pv) const {
	    return ((sizeof(var_t) > 2)? INT_MAX:
	    (int((check_scr - abs(pv)) / avmch / nelem * nelem)));
	}
public:
	VTYPE	scoreonlyS1();
	VTYPE	scoreonlyS1_wip();
	VTYPE	forwardS1(int pp[]);
	VTYPE	forwardS1_wip(Mfile* mfd);
	VTYPE	hirschbergS1(Dim10* cpos, const int& n_im);
	VTYPE	hirschbergS1_wip(Dim10* cpos, const int& n_im);
// Constructor
	SimdAln2s1(const Seq** sqs, const PwdB* _pwd, const WINDOW& _w, 
	    SpJunc* _spj, const Cip_score* _cip, int _mode = 1, Vmf* _vmf = 0) :
	    seqs(sqs), a(sqs[0]), b(sqs[1]), pwd(_pwd), wdw(_w),
	    Local(algmode.lcl & 16), spjcs(_spj), cip(_cip), 
	    mode(_mode), used(sizeof(var_t) == 2 && (_mode & 4)), 
	    spj(b->inex.intr), dagp(pwd->Noll == 3), 
	    nod(dagp? 5: 3), Ncand(NCAND + (dagp? 2: 0)), vmf(_vmf), 
	    Np1(nelem + 1), avmch(pwd->simmtx->AvTrc()), 
#if FVAL
	    nevsel(NEVSEL),
#else
	    nevsel((sizeof(var_t) == 2? SHRT_MIN: INT_MIN) + nevsel_play),
#endif
	    black_Rvulmn{nevsel, end_of_ulk, a->left, a->right, b->right},
	    rlst(INT_MAX), buf_size(wdw.width + 2 * nelem)
{

#define	Add(a, b)	this->add(a, b)
#define	Blend(a, b, m)	this->blend(a, b, m)
#define	Clear()		this->clear()
#define	Cmp_eq(a, b)	this->cmp_eq(a, b)
#define	Cmp_gt(a, b)	this->cmp_gt(a, b)
#define	Load(a)		this->load(a)
#define	Splat(c)	this->splat(c)
#define	Store(a, v)	this->store(a, v)
#define And(a, b)	this->bit_and(a, b)
#define Or(a, b)	this->bit_or(a, b)
#define AndNot(a, b)	this->bit_andnot(a, b)
#define Cast16to8(a)	this->cast16to8(a)

/*****************************************************************
*	mode	v	b	c	d	o	
*	  0	4 + 2	0 + 0	0 + 0	0 + 0	2 + 2	score only
*	  1	4 + 2	4 + 2	0 + 0	0 + 0	2 + 2	bimap
*	  2	4 + 2	4 + 2	4 + 2	0 + 0	2 + 2	udh1
*	  3	4 + 2	4 + 2	4 + 2	0 + 0	2 + 4	vmf1
*	  4	4 + 2	4 + 2	4 + 2	4 + 2	2 + 2	udh2
*	  5	4 + 2	4 + 2	4 + 2	4 + 2	2 + 4	vmf2
******************************************************************/

	    size_t		bufsiz =  10 * Np1 + 6 * nelem;
	    if (mode > 1)	bufsiz += (8 * Np1 + 6 * nelem);
	    abuf = new var_t[bufsiz];
	    ps_a = abuf;
	    pv_a = ps_a + nelem;
	    hv_a[0] = pv_a + nelem;
	    hv_a[1] = hv_a[0] + Np1;
	    ev_a = hv_a[1] + Np1;
	    fv_a = ev_a + nelem;
	    ev2_a = fv_a + Np1;
	    fv2_a = ev2_a + nelem;
	    s5_a = fv2_a + Np1;
	    s3_a = s5_a + Np1;
	    hb_a[0] = s3_a + Np1;
	    hb_a[1] = hb_a[0] + Np1;
	    eb_a = hb_a[1] + Np1;
	    fb_a = eb_a + nelem;
	    eb2_a = fb_a + Np1;
	    fb2_a = eb2_a + nelem;
	    if (mode > 1) {
		pb_a = fb2_a + Np1;
		qb_a = pb_a + nelem;
		hc_a[0] = qb_a + nelem;
		hc_a[1] = hc_a[0] + Np1;
		ec_a = hc_a[1] + Np1;
		fc_a = ec_a + nelem;
		ec2_a = fc_a + Np1;
		fc2_a = ec2_a + nelem;
		hd_a[0] = fc2_a + Np1;
		hd_a[1] = hd_a[0] + Np1;
		ed_a = hd_a[1] + Np1;
		fd_a = ed_a + nelem;
		ed2_a = fd_a + Np1;
		fd2_a = ed2_a + nelem;
	    }

	    for (int p = 0; p < 2; ++p) {
		for (int j = 0; j < nelem; ++j) {
		    hfesv[p][j][PS_A] = ps_a + j;
		    hfesv[p][j][PV_A] = pv_a + j;
		    hfesv[p][j][0] = hv_a[p] + j + 1;
		    hfesv[p][j][1] = ev_a + j;
		    hfesv[p][j][2] = fv_a + j + 1;
		    hfesv[p][j][3] = ev2_a + j;
		    hfesv[p][j][4] = fv2_a + j + 1;
		    hfesb[p][j][0] = hb_a[p] + j + 1;
		    hfesb[p][j][1] = eb_a + j;
		    hfesb[p][j][2] = fb_a + j + 1;
		    hfesb[p][j][3] = eb2_a + j;
		    hfesb[p][j][4] = fb2_a + j + 1;
		    if (mode <= 1) continue;
		    hfesc[p][j][0] = hc_a[p] + j + 1;
		    hfesc[p][j][1] = ec_a + j;
		    hfesc[p][j][2] = fc_a + j + 1;
		    hfesc[p][j][3] = ec2_a + j;
		    hfesc[p][j][4] = fc2_a + j + 1;
		    hfesd[p][j][0] = hd_a[p] + j + 1;
		    hfesd[p][j][1] = ed_a + j;
		    hfesd[p][j][2] = fd_a + j + 1;
		    hfesd[p][j][3] = ed2_a + j;
		    hfesd[p][j][4] = fd2_a + j + 1;
		}
	    }

	    vbuf = new var_t[((mode > 1)? (used? 4: 3): 2) * 
		pwd->Noll * buf_size];
	    vec_set(vbuf, nevsel, pwd->Noll * buf_size);

	    hv = vbuf - wdw.lw + 1;
	    fv = hv + buf_size;
	    if (dagp) fv2 = fv + buf_size;
	    bbuf = vbuf + pwd->Noll * buf_size;
	    hb = bbuf - wdw.lw + 1;
	    fb = hb + buf_size;
	    if (dagp) fb2 = fb + buf_size;
	    if (mode > 1) {
		cbuf = bbuf + pwd->Noll * buf_size;
		hc = cbuf - wdw.lw + 1;
		fc = hc + buf_size;
		if (dagp) fc2 = fc + buf_size;
		dbuf = cbuf + pwd->Noll * buf_size;
		hd = dbuf - wdw.lw + 1;
		fd = hd + buf_size;
		if (dagp) fd2 = fd + buf_size;
	    }

	    if (spj) {
		fsjss = new Sjsites*[nelem];
		for (int k = 0; k < nelem; ++k)
		    fsjss[k] = new Sjsites(*this);
		donor_q = new Queue2<int>(nelem);
		accep_q = new Queue2<int>(nelem);
	    }
	}
	~SimdAln2s1() {
	    delete[] abuf; delete[] vbuf;
	    if (fsjss) {
		for (int k = 0; k < nelem; ++k) delete fsjss[k];
		delete[] fsjss;
		delete donor_q;
		delete accep_q;
	    }
	}
 	int	s2_i(const var_t& lb, const var_t& ub = 0) {
	    if (used) {
		int	i = ub << 16;
		return (i | SHORT(lb));
	    } else
		return (INT(lb));
	}
	var_t	i_s2(var_t* ub, const int& i) {
	    if (used) {
		*ub = i >> 16;
		return (i & 0xffff);
	    } else 
		return (INT(i));
	}
friend	class	Sjsites;
};

#endif	// _FWD2S1_SIMD_H_
