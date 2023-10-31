/*****************************************************************************
*
*	Vetorized version of DP-based spliced alignment of protein vs DNA 
*	Accelerated by SIMD instructions of Intel intrinsics.
*	Curently accepted architectures: SSE4.1, AVX, AVX2, and AVX-512BW
*	Assumes capability of store/load in/from unaligned memory
*
*	5' and 3' splice site signals, intron-length distribution,
*	coding potential, and conservation of intron sites  are considered.
*	Assumes there is no internal gap in the reference protein sequence.
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

#ifndef	_FWD2H1_SIMD_H_
#define _FWD2H1_SIMD_H_

#include "aln.h"
#include "vmf.h"

#if _SIMD_PP_
#include <simdpp/simd.h>
#endif	// _SIMD_PP_

#include "simd_functions.h"
#include "udh_intermediate.h"
#include "rhomb_coord.h"

using	TBU_t = SHORT;

static	const	int	check_scr = int(0.9 * SHRT_MAX);
static	const	int	PS_A = 4;
static	const	int	PV_A = 5;
static	const	int	nevsel_play = 1024;
static	const	int	max_ncand = 4;

struct Rvulmn {		// used in Local mode
	VTYPE	val;	// score
	int	ulk;	// lower end diagonal
	SHORT	ml;	// left end m-coord
	SHORT	mr;	// right end m-coord
	int	nr;	// right end n-coord
};

template <typename var_t>
struct Rvlujd {
	var_t	val;	// score
	int	ulk;	// lower end diagonal
	int	jnc;	// junction
	SHORT	ml;	// left end m-coord
	short	dir;	// direction
	short	phs;	// phase
};

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
#if _SIMD_PP_
class SimdAln2h1 {
#else
class SimdAln2h1 : public Simd_functions<var_t, Nelem, regist_v> {
#endif

protected:
const	Seq**	seqs;
const	Seq*&	a;
const	Seq*&	b;
const	PwdB*	pwd;
const	WINDOW&	wdw;
const	bool	Local;
	SpJunc*	spjcs;
const	Cip_score* cip;
const	int	mode;	// 0: score only, 1: traceback, 2: linear space
const	bool	used;
const	bool	spj;
	Vmf*	vmf;
	int	mm = 0;
	int	mm3 = 0;
	UdhIntermediate*	imd = 0;
const	int	Np1;
const	int	avmch;
const	var_t	nevsel;
const	Rvulmn	black_Rvulmn;
	int	rlst[3];
	size_t	buf_size;
	var_t*	abuf;
	var_t*	sm_a;		// aa similarity
	var_t*	cp_a[3];	// coding potential
	var_t*	s5_a[6];	// 5' splice signal
	var_t*	s3_a[6];	// 3' splice signal
	var_t*	p5_a[6];	// 5' splice phase
	var_t*	p3_a[6];	// 3' splice phase
	var_t*	ps_a[3];
	var_t*	pv_a[3];
	var_t*	qv_a[3];
	var_t*	hv_a[6];
	var_t*	fv_a[6];
	var_t*	ev_a[3];
	var_t*	hb_a[6];
	var_t*	fb_a[6];
	var_t*	eb_a[3];
	var_t*	qb_a[3];
	var_t*	hc_a[6];
	var_t*	fc_a[6];
	var_t*	fc2_a[6];
	var_t*	ec_a[3];
	var_t*	ec2_a[3];
	var_t*	qc_a[3];
	var_t*	hd_a[6];
	var_t*	fd_a[6];
	var_t*	fd2_a[6];
	var_t*	ed_a[3];
	var_t*	ed2_a[3];
	var_t*	qd_a[3];
 	var_t*	vbuf = 0;
 	var_t*	bbuf = 0;
 	var_t*	cbuf = 0;
 	var_t*	dbuf = 0;
	var_t*	hv, *fv, *fv2;	// score
	var_t*	hb, *fb, *fb2;	// boundary
	var_t*	hc, *fc, *fc2;	// coordinate
	var_t*	hd, *fd, *fd2;	// coordinate2
	var_t*	hfesv[6][Nelem][6];
	var_t*	hfesb[6][Nelem][4];
	var_t*	hfesc[6][Nelem][4];
	var_t*	hfesd[6][Nelem][4];
#if _SIMD_PP_
#if __AVX2__
regist_v simdpp_cast16to8(regist_v v) {
	SHORT	buf[16];
	simdpp::store(buf, v);
	__m256i v_v = _mm256_loadu_si256((__m256i*) buf);
	__m256i	b_v = _mm256_loadu_si256((__m256i*) b32_a);
	v_v = _mm256_shuffle_epi8(v_v, b_v);
	v_v = _mm256_permute4x64_epi64(v_v, 216);
	_mm256_storeu_si256((__m256i*) buf, v_v);
	return (simdpp::load(buf));
}
#elif __SSE4_1__
regist_v simdpp_cast16to8(regist_v v) {
	SHORT	buf[8];
	simdpp::store(buf, v);
	__m128i v_v = _mm_loadu_si128((__m128i const*) buf);
	__m128i b_v = _mm_loadu_si128((__m128i const*) b32_a);
	v_v = _mm_shuffle_epi8(v_v, b_v);
	_mm_storeu_si128((__m128i*) buf, v_v);
	return (simdpp::load(buf));
}
#endif	// __AVX2__
#endif	// _SIMD_PP_

class Sjsites {			// nested class
	SimdAln2h1& hb1;
	Rvlujd<var_t>	rcd[max_ncand + 1];
	short	idx[max_ncand + 1];
	short	ncand;
	Rvlujd<var_t>	black_r;
public:
	void	reset() {
	    ncand = -1;
	    for (int i = 0; i <= max_ncand; ++i) {
		rcd[i] = black_r;
		idx[i] = i;
	    }
	}
	Sjsites(SimdAln2h1& _hb1) : hb1(_hb1) {
	    black_r.val = hb1.nevsel;
	    black_r.ulk = black_r.jnc = black_r.ml = black_r.dir = 0;
	    black_r.phs = -2;
	}
	~Sjsites() {}
	void	put(const int& j, const int& m, int n, int q);
	void	get(const int& j, const int& m, const int& n, const int& q);
};	// end of Sjsites

	Sjsites**	fsjss = 0;
	Queue2<int>*	donor_q[3];
	Queue2<int>*	accep_q[3];

	void	fhinitH1(Anti_rhomb_coord<TBU_t>* trb = 0);
	int	fhlastH1(Rvulmn& maxh, Anti_rhomb_coord<TBU_t>* trb = 0);
	void	from_spj(Queue2<int>& l, const int& m, 
		const int& n, const int& q);
	void	to_spj(Queue2<int>& l, const int& m, 
		const int& n, const int& q);
	void	vadd(var_t* ar, const int& n, const var_t& c) {
	    for (var_t* er = ar + n; ar < er; ++ar) {
		VTYPE	x = *ar + c;
		*ar = (x < nevsel)? nevsel: static_cast<var_t>(x);
	    }
	}
	void	vec_add(var_t* ar, const int& n, const var_t& c);
	var_t	vec_max(const var_t* ar, const int& n);
	int	checkpoint(const var_t& pv) const {
	    return ((sizeof(var_t) > 2)? INT_MAX:
	    (int((check_scr - pv) / avmch / Nelem * Nelem)));
	}
public:
	VTYPE	forwardH1(int pp[]);
	VTYPE	forwardH1_wip(Mfile* mfd);
	VTYPE	hirschbergH1(Dim10* cpos, const int& n_im);
	VTYPE	hirschbergH1_wip(Dim10* cpos, const int& n_im);
// Constructor
	SimdAln2h1(const Seq** sqs, const PwdB* _pwd, const WINDOW& _w, 
	    SpJunc* _spj, const Cip_score* _cip, int _mode = 1, Vmf* _vmf = 0) :
	    seqs(sqs), a(sqs[0]), b(sqs[1]), pwd(_pwd), wdw(_w), 
	    Local(algmode.lcl & 16), spjcs(_spj), cip(_cip), 
	    mode(_mode), used(sizeof(var_t) == 2 && (_mode & 4)), 
	    spj(b->inex.intr), vmf(_vmf), Np1(Nelem + 1), 
	    avmch(pwd->simmtx->AvTrc()), 
#if FVAL
	    nevsel(NEVSEL),
#else
	    nevsel((sizeof(var_t) == 2? SHRT_MIN: INT_MIN) + nevsel_play),
#endif
	    black_Rvulmn{nevsel, end_of_ulk, 
		SHORT(a->left), SHORT(a->right), b->right},
	    buf_size(wdw.width + 6 * Nelem)
{

#if _SIMD_PP_
#define	Add(a, b)	simdpp::add(a, b)
#define	Blend(a, b, m)	simdpp::blend(a, b, m)
#define	Clear()		simdpp::splat(0)
#define	Cmp_eq(a, b)	simdpp::cmp_eq(a, b)
#define	Cmp_gt(a, b)	simdpp::cmp_gt(a, b)
#define	Load(a)		simdpp::load(a)
#define	Splat(c)	simdpp::splat(c)
#define	Store(a, v)	simdpp::store(a, v)
#define And(a, b)	simdpp::bit_and(a, b)
#define Or(a, b)	simdpp::bit_or(a, b)
#define AndNot(a, b)	simdpp::bit_andnot(a, b)
#define To_mask(a)	simdpp::to_mask((regist_v) simdpp::load(a))
#define Cast16to8(a)	simdpp_cast16to8(a) 
#define AllZero(a)	(simdpp::reduce_max(a) == 0)
#else	// Intel Intrinsics
#define	Add(a, b)	this->add(a, b)
#define	Blend(a, b, m)	this->blend(a, b, m)
#define	Clear()		this->clear()
#define	Cmp_eq(a, b)	this->cmp_eq(a, b)
#define	Cmp_gt(a, b)	this->cmp_gt(a, b)
#define	Load(a)		this->load(a)
#define	Splat(c)	this->splat(c)
#define	Store(a, v)	this->store(a, v)
#define And(a, b)	this->bit_and(a, b)
#define AndNot(a, b)	this->bit_andnot(a, b)
#define Or(a, b)	this->bit_or(a, b)
#define To_mask(a)	this->load(a)
#define Cast16to8(a)	this->cast16to8(a)
#define AllZero(a)	this->all_zero(a)
#endif	// _SIMD_PP_

/*****************************************************************
*	mode	v	b	c	d	o	
*	  1	12 + 3	12 + 3	 0 + 0	 0 + 0 28 + 6	bimap
*	  alg2	 0 + 0	 0 + 0	 0 + 0	 0 + 0	0 + 0	bimap2
*	  2	 0 + 0	 0 + 0	12 + 3	 0 + 0	0 + 0	udh1
*	  3	 0 + 0	 0 + 0	12 + 3	 0 + 0	0 + 0	vmf1
*	  4	 0 + 0	 0 + 0	12 + 3	12 + 3	0 + 0	udh2
*	  5	 0 + 0	 0 + 0	12 + 3	12 + 3	0 + 0	vmf2
******************************************************************/

	    size_t		abufsiz  = 52 * Np1 + 12 * Nelem;
	    if (mode > 1)	abufsiz += 24 * Np1 + 18 * Nelem;

	    abuf = new var_t[abufsiz];
	    sm_a = abuf;
	    cp_a[0] = sm_a + Np1;
	    s5_a[0] = cp_a[0] + 3 * Np1;
	    s3_a[0] = s5_a[0] + 6 * Np1;
	    p5_a[0] = s3_a[0] + 6 * Np1;
	    p3_a[0] = p5_a[0] + 6 * Np1;
	    ps_a[0] = p3_a[0] + 6 * Np1;
	    pv_a[0] = ps_a[0] + 3 * Nelem;
	    hv_a[0] = pv_a[0] + 3 * Nelem;
	    fv_a[0] = hv_a[0] + 6 * Np1;
	    ev_a[0] = fv_a[0] + 6 * Np1;
	    hb_a[0] = ev_a[0] + 3 * Nelem;
	    fb_a[0] = hb_a[0] + 6 * Np1;
	    eb_a[0] = fb_a[0] + 6 * Np1;
	    if (mode > 1) {
		qv_a[0] = eb_a[0] + 3 * Nelem;
		hc_a[0] = qv_a[0] + 3 * Nelem;
		fc_a[0] = hc_a[0] + 6 * Np1;
		ec_a[0] = fc_a[0] + 6 * Np1;
		hd_a[0] = ec_a[0] + 3 * Nelem;
		fd_a[0] = hd_a[0] + 6 * Np1;
		ed_a[0] = fd_a[0] + 6 * Np1;
		qb_a[0] = ed_a[0] + 3 * Nelem;
		qc_a[0] = qb_a[0] + 3 * Nelem;
		qd_a[0] = qc_a[0] + 3 * Nelem;
	    }

	    for (int q = 1; q < 6; ++q) {
		s5_a[q] = s5_a[q - 1] + Np1;
		s3_a[q] = s3_a[q - 1] + Np1;
		p5_a[q] = p5_a[q - 1] + Np1;
		p3_a[q] = p3_a[q - 1] + Np1;
		hv_a[q] = hv_a[q - 1] + Np1;
		fv_a[q] = fv_a[q - 1] + Np1;
		hb_a[q] = hb_a[q - 1] + Np1;
		fb_a[q] = fb_a[q - 1] + Np1;
		if (mode > 1) {
		    hc_a[q] = hc_a[q - 1] + Np1;
		    fc_a[q] = fc_a[q - 1] + Np1;
		    hd_a[q] = hd_a[q - 1] + Np1;
		    fd_a[q] = fd_a[q - 1] + Np1;
		}
	    }
	    for (int p = 1; p < 3; ++p) {
		cp_a[p] = cp_a[p - 1] + Np1;
		ps_a[p] = ps_a[p - 1] + Nelem;
		pv_a[p] = pv_a[p - 1] + Nelem;
		ev_a[p] = ev_a[p - 1] + Nelem;
		eb_a[p] = eb_a[p - 1] + Nelem;
		if (mode > 1) {
		    qv_a[p] = qv_a[p - 1] + Nelem;
		    qb_a[p] = qb_a[p - 1] + Nelem;
		    ec_a[p] = ec_a[p - 1] + Nelem;
		    ed_a[p] = ed_a[p - 1] + Nelem;
		    qc_a[p] = qc_a[p - 1] + Nelem;
		    qd_a[p] = qd_a[p - 1] + Nelem;
		}
	    }
	    for (int q = 0; q < 6; ++q) {
		int	p = q % 3;
		for (int j = 0; j < Nelem; ++j) {
		    hfesv[q][j][0] = hv_a[q] + j + 1;
		    hfesv[q][j][1] = ev_a[p] + j;
		    hfesv[q][j][2] = fv_a[q] + j + 1;
		    hfesv[q][j][PS_A] = ps_a[p] + j;
		    hfesv[q][j][PV_A] = pv_a[p] + j;
		    hfesb[q][j][0] = hb_a[q] + j + 1;
		    hfesb[q][j][1] = eb_a[p] + j;
		    hfesb[q][j][2] = fb_a[q] + j + 1;
		    if (mode > 1) {
			hfesb[q][j][3] = qb_a[p] + j;
			hfesv[q][j][3] = qv_a[p] + j;
			hfesc[q][j][0] = hc_a[q] + j + 1;
			hfesc[q][j][1] = ec_a[p] + j;
			hfesc[q][j][2] = fc_a[q] + j + 1;
			hfesc[q][j][3] = qc_a[p] + j;
			hfesd[q][j][0] = hd_a[q] + j + 1;
			hfesd[q][j][1] = ed_a[p] + j;
			hfesd[q][j][2] = fd_a[q] + j + 1;
			hfesd[q][j][3] = qd_a[p] + j;
		    }
		}
	    }

	    vset(rlst, INT_MAX, 3);
	    vbuf = new var_t[((mode > 1)? (used? 8: 6): 4) * buf_size];
	    bbuf = vbuf + 2 * buf_size;
	    hv = vbuf - wdw.lw + 3;
	    fv = hv + buf_size;
	    hb = bbuf - wdw.lw + 3;
	    fb = hb + buf_size;
	    if (mode > 1) {
		cbuf = bbuf + 2 * buf_size;
		dbuf = used? (cbuf + 2 * buf_size): 0;
		hc = cbuf - wdw.lw + 3;
		fc = hc + buf_size;
		hd = used? (dbuf - wdw.lw + 3): 0;
		fd = used? (hd + buf_size): 0;
	    }

	    if (spj) {
		fsjss = new Sjsites*[Nelem];
		for (int k = 0; k < Nelem; ++k)
		    fsjss[k] = new Sjsites(*this);
		for (int p = 0; p < 3; ++p) {
		    donor_q[p] = new Queue2<int>(Nelem);
		    accep_q[p] = new Queue2<int>(Nelem);
		}
	    } else {
		vclear(donor_q, 3);
		vclear(accep_q, 3);
	    }
	}
	~SimdAln2h1() {
	    delete[] abuf; delete[] vbuf;
	    if (fsjss) {
		for (int k = 0; k < Nelem; ++k) delete fsjss[k];
		delete[] fsjss;
		for (int p = 0; p < 3; ++p) {
		    delete donor_q[p];
		    delete accep_q[p];
		}
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

/*************************************************************************
	nested class: Sjsites
*************************************************************************/

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2h1<var_t, Nelem, regist_v, regist_m>::Sjsites::
get(const int& j, const int& m, const int& n, const int& q)
{
const	int	acc = n - 1;
const	int	r = acc - 3 * (m + 1);
const	Rvlujd<var_t>*	maxprd[3];
const	Rvlujd<var_t>*	brd = 0;
const	CHAR*	as = hb1.a->at(m);
const	bool	is_imd = hb1.imd && (m + 1) == hb1.imd->mi;
	vclear(maxprd, 3);
	for (int l = 0; l <= ncand; ++l) {
const	    Rvlujd<var_t>*	prd = rcd + idx[l];
const	    int	rr = r + prd->phs;
	    if (rr < hb1.wdw.lw || rr >= hb1.wdw.up) continue;
const	    Rvlujd<var_t>*&	mrd = maxprd[prd->dir];
const	    int&	don = prd->jnc;
const	    int&	d = prd->dir;
	    if (d == 2 && prd->phs == 1) continue;
	    if ((acc - don) < IntronPrm.minl) continue;
	    VTYPE	x = hb1.cip? hb1.cip->cip_score(3 * (m + 1) - prd->phs): 0;
	    x += prd->val + hb1.spjcs->spjscr(don, acc);
	    if (d == 0 && prd->phs) {
const	CHAR*	cs = hb1.spjcs->spjseq(don, acc);
		if (prd->phs == -1) ++cs;
		x += hb1.pwd->simmtx->mtx[*as][*cs];
	    }
	    if (prd->phs == -1) x -= hb1.b->exin->data_p[acc].sigE;
const	    int	qq = modN<6>(q + prd->phs);
	    var_t**	hfesmv = hb1.hfesv[qq][j];
	    if (x <= *hfesmv[d]) continue;
	    if (!mrd || x > mrd->val) {
		mrd = prd;
		if (!brd || x > brd->val) brd = prd;
	    }
	    *hfesmv[d] = static_cast<var_t>(x);
	    *hfesmv[PS_A] |= psp_bit[d];	// mark spliced
	    var_t**	hfesmb = hb1.hfesb[qq][j];
	    var_t**	hfesmc = hb1.hfesc[qq][j];
	    var_t**	hfesmd = hb1.hfesd[qq][j];
	    if (hb1.vmf) {
		*hfesmb[d] = prd->ml;
const		int ptr = hb1.vmf->add(m + 1, acc + prd->phs, 
		    hb1.vmf->add(m + 1, don + prd->phs, prd->ulk));
		*hfesmc[d] = hb1.i_s2(hfesmd[d], ptr);
	    } else {
		*hfesmb[d] = prd->ml;
		*hfesmc[d] = hb1.i_s2(hfesmd[d], prd->ulk);
	    }
	    if (d && *hfesmv[d] > *hfesmv[0]) {
		*hfesmv[0] = *hfesmv[d];
		*hfesmb[0] = *hfesmb[d];
		*hfesmc[0] = *hfesmc[d];
		if (hb1.used) *hfesmd[0] = *hfesmd[d];
	    }
	    if (j + 1 == Nelem) {
		hb1.hv[rr] = *hfesmv[0];
		if (is_imd) {
		    hb1.imd->hlnk[0][rr] = prd->ulk;
		} else {
		    hb1.hb[rr] = *hfesmb[0];
		    hb1.hc[rr] = *hfesmc[0];
		    if (hb1.used) hb1.hd[rr] = *hfesmd[0];
		}
		if (d == 2) {
		    hb1.fv[rr] = *hfesmv[d];
		    if (is_imd) {
			hb1.imd->hlnk[1][rr] = prd->ulk;
		    } else {
			hb1.fb[rr] = *hfesmb[d];
			hb1.fc[rr] = *hfesmc[d];
			if (hb1.used) hb1.fd[rr] = *hfesmd[d];
		    }
		}
	    }
	}

	if (is_imd && brd) {
const	    int&	maxd = brd->dir;
const	    Rvlujd<var_t>*	prd = maxprd[maxd];
const	    int	qq = modN<6>(q + prd->phs);
const	    int&	lstr = hb1.rlst[qq % 3] = acc + prd->phs - hb1.mm3;
	    hb1.imd->hlnk[0][lstr] = prd->ulk;
	    var_t**	hfesmv = hb1.hfesv[qq][j];
	    var_t**	hfesmc = hb1.hfesc[qq][j];
	    var_t**	hfesmd = hb1.hfesd[qq][j];
	    *hfesmc[maxd] = hb1.i_s2(hfesmd[maxd], lstr);
	    *hfesmv[PV_A] = maxd;
	    if (maxd) {
		*hfesmc[0] = *hfesmc[maxd];
		*hfesmd[0] = *hfesmd[maxd];
		return;
	    }
	    if ((prd = maxprd[1]) && 
		*hfesmv[1] > (*hfesmv[0] + hb1.pwd->BasicGOP)) {
		hb1.imd->hlnk[1][lstr] = prd->ulk;
		*hfesmc[1] = hb1.i_s2(hfesmd[1], lstr + hb1.wdw.width);
	    }
	    if (maxprd[2] && 
		*hfesmv[2] > (*hfesmv[0] + hb1.pwd->BasicGOP)) {
		*hfesmc[2] = hb1.i_s2(hfesmd[2], lstr + hb1.wdw.width);
	    }
	}
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2h1<var_t, Nelem, regist_v, regist_m>::Sjsites::
put(const int& j, const int& m, int n, int q)
{
const	int	don = n - 1;
const	VTYPE	sigJ = hb1.b->exin->sig53(don, don, IE5);
const	bool	is_imd = hb1.imd && (m + 1) == hb1.imd->mi;
	for (int phs = 1; phs > -2; --n, --phs, q = modN<6>(q - 1)) {
const	    int	rr = don - 3 * (m + 1) + phs;
	    if (rr < hb1.wdw.lw || rr >= hb1.wdw.up) continue;
	    var_t**	hfesmv = hb1.hfesv[q][j];
	    var_t**	hfesmb = hb1.hfesb[q][j];
	    var_t**	hfesmc = hb1.hfesc[q][j];
	    var_t**	hfesmd = hb1.hfesd[q][j];
const	    int&	h = *hfesmv[PV_A];
const	    VTYPE	thrscr = *hfesmv[0] + hb1.pwd->BasicGOP;
	    for (int k = (h && phs < 1)? 1: 0; k < 3; ++k) {
		if (*hfesmv[PS_A] & psp_bit[k]) continue;	// prevent orphan exon
const		int	crossspj = (phs == 1 && k == 0)? 3: k;	// qv_a: xv_a[q]
const		var_t&	from = *hfesmv[crossspj];
		if (k && from <= thrscr) continue;		// prune
const		VTYPE	x = from + sigJ;
		if (x <= hb1.nevsel) continue;
		short	l = ncand < max_ncand? ++ncand: max_ncand;
		while (--l >= 0) {		// sort on score
		    if (x >= rcd[idx[l]].val)
			std::swap(idx[l], idx[l + 1]);
		    else
			break;
		}
		if (++l < max_ncand) {
		    Rvlujd<var_t>*	prd = rcd + idx[l];
		    prd->val = static_cast<var_t>(x);
		    prd->ml = *hfesmb[k];
const		    int	r = n - hb1.mm3;
const		    int	rl = hb1.s2_i(*hfesmc[crossspj], *hfesmd[crossspj]);
		    if (is_imd) {
			if (k == 1) hb1.imd->hlnk[0][r] = hb1.rlst[q % 3];
			else if (r != rl) hb1.imd->vlnk[k / 2][r] = rl;
			prd->ulk = r;
		    } else {
			prd->ulk = rl;
		    }
		    prd->jnc = don;
		    prd->dir = k;
		    prd->phs = phs;
		} else --ncand;
	    }
	}
}

/*************************************************************************
	normal member functions of SimdAln2h1
*************************************************************************/

// add const c to each element of array ar of size n
template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2h1<var_t, Nelem, regist_v, regist_m>::
vec_add(var_t* ar, const int& n, const var_t& c) {
#if _SIMD_PP_
	vadd(ar, n, c);
#else
const	int	nn = n / Nelem * Nelem;
	if (nn) {
regist_v    c_v = Splat(c);
regist_v    n_v = Splat(nevsel);
	    for (var_t* er = ar + nn; ar < er; ar += Nelem) {
		regist_v	a_v = Load(ar);
		a_v = Add(a_v, c_v);
		regist_m	msk_m = Cmp_gt(a_v, n_v);
		a_v = Blend(a_v, n_v, msk_m);
		Store(ar, a_v);
	    }
	}
	if (n > nn) vadd(ar, n - nn, c);
#endif
}

// max in array ar of size n
template <typename var_t, int Nelem, typename regist_v, typename regist_m>
var_t SimdAln2h1<var_t, Nelem, regist_v, regist_m>::
vec_max(const var_t* ar, const int& n) {
const	int	nn = n / Nelem * Nelem;
	var_t	buf[Nelem + Nelem / 2];
regist_v    a_v = Splat(ar[0]);
	Store(buf, a_v); Store(buf + Nelem / 2, a_v); 
	if (n > nn) vcopy(buf, ar + nn, n - nn);
	a_v = Load(buf);

	if (nn) {
const	    var_t*	er = ar + nn;
	    for ( ; (ar + Nelem) < er; ar += Nelem) {
regist_v	b_v = Load(ar);
regist_m	msk_m = Cmp_gt(a_v, b_v);
		a_v = Blend(a_v, b_v, msk_m);
	    }
	    Store(buf, a_v);
	}
// O(log_2(Nelem)) algorithm
	for (int d = Nelem; (d >>= 1); ) {
regist_v    b_v = Load(buf + d);
regist_m    msk_m = Cmp_gt(a_v, b_v);
	    a_v = Blend(a_v, b_v, msk_m);
	    Store(buf, a_v);
	}
	return (buf[0]);
}

// prepare variables and initialize DP matrix
template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2h1<var_t, Nelem, regist_v, regist_m>::
fhinitH1(Anti_rhomb_coord<TBU_t>* trb)
{
	vec_set(vbuf, nevsel, 2 * buf_size);
const	int	rl = b->left - 3 * a->left;
const	bool	usec = mode > 1;
	TBU_t*	row0 = trb? trb->set_point(a->left, b->left): 0;

	if (vmf) {	// ordinary traceback
	    vec_clear(bbuf, 2 * buf_size);	// non_diagonal
	    int	ptr = vmf->add(0, 0, 0);	// dummy
	    if (!(a->inex.exgl && b->inex.exgl))
		ptr = vmf->add(a->left, b->left, ptr);
	    if (used) {
		var_t	upb;
		var_t	lwb = i_s2(&upb, ptr);
	 	if (a->inex.exgl) {
		    vec_clear(hc + rl, wdw.up - rl);
		    vec_clear(hd + rl, wdw.up - rl);
		} else {
		    vec_set(hc + rl, lwb, wdw.up - rl);
		    vec_set(hd + rl, upb, wdw.up - rl);
		}
		if (b->inex.exgl) {
		    vec_clear(hc + wdw.lw, rl - wdw.lw);
		    vec_clear(hd + wdw.lw, rl - wdw.lw);
		} else {
		    vec_set(hc + wdw.lw, lwb, rl - wdw.lw);
		    vec_set(hd + wdw.lw, upb, rl - wdw.lw);
		}
	    } else {
		if (a->inex.exgl) 
		    vec_clear(hc + rl, wdw.up - rl);
		else
		    vec_set(hc + rl, var_t(ptr), wdw.up - rl);
		if (b->inex.exgl)
		    vec_clear(hc + wdw.lw, rl - wdw.lw);
		else
		    vec_set(hc + wdw.lw, var_t(ptr), rl - wdw.lw);
	    }
	    if (b->inex.exgl == 2)
		fc[rl] = i_s2(fd? fd + rl: 0, ptr);
	} else if (usec) {			// Hirschberg
	    vec_set(bbuf, var_t(a->left), 2 * buf_size);
	    int	r = wdw.lw;
	    int	rr = a->inex.exgl? rl: wdw.up;
	    var_t*	c = usec? hc + r: 0;
	    var_t*	d = used? hd + r: 0;
	    while (r < rr)
		*c++ = i_s2(d? d++: 0, r++);
	    for (int i = 0, r = rl; r >= wdw.lw; --r)
		hb[r] = a->left + (i++ / 3);	// ml
	}

	if (b->inex.exgl == 1)	// vertical semi-global
	    vec_clear(hv + wdw.lw, rl - wdw.lw);
	else if (b->inex.exgl == 2) {
	    fv[rl] = 0;
	    if (usec) fc[rl] = i_s2(fd? fd + rl: 0, rl);
	}

	int	rr = b->right - 3 * a->left;
	if (wdw.up < rr) rr = wdw.up;
	int	r = rl;
// global
	if (!a->inex.exgl) {
	    if (b->inex.exgl) {
		fv[r] = 0;	// no gop for terminal gap
		if (usec) fc[r] = hc[r];
		if (used) fd[r] = hd[r];
	    }
	    hv[r++] = 0;
	    hv[r++] = pwd->GapW1;
	    hv[r++] = pwd->GapW2;
	    hv[r++] = pwd->GapW3;
	    if (pwd->BasicGEP) {
		int	x = (nevsel - pwd->GapW3) / pwd->BasicGEP + r;
		if (x < rr) rr = x;
		for ( ; r < rr; ++r)
		    hv[r] = hv[r - 3] + pwd->BasicGEP;
	    } else if (rr > r)
		vec_set(hv + r, hv[r - 1], rr - r);
	    return;
	}

// semi-global
	var_t*	h = hv + r;
	var_t*	e = vmf? hb + r: 0;
	var_t*	c = usec? (hc + r): 0;
	var_t*	d = used? (hd + r): 0;
	int	n = b->left;
	int	lend[3] = {r, r + 1, r + 2};
const	SGPT6*	bb = b->exin->score_p(n + 1);
	for (int p = 0; p < 3; ++r, ++n, ++bb, ++p) {
	    *h++ = (bb->sigS > 0)? bb->sigS: 0;
	    if (vmf) {
		int	ptr = vmf->add(a->left, n, 0);
		*c++ = i_s2(d? d++: 0, ptr);
		*e++ = 1;	// HORI
	    } else if (usec)
		*c++ = i_s2(d? d++: 0, r);
	    if (trb) trb->to_right(row0);
	}
	for (int p = 0; r < rr; ++r, ++h, ++n, ++bb, p = next_p[p]) {
	    *h = h[-3];
	    if (c) *c = c[-3];
	    if (d) *d = d[-3];
const	    int	gl = r - lend[p];
	    if (!(a->inex.exgl & 1) && gl == 3) *h += pwd->BasicGOP;
	    if (!(a->inex.exgl & 2)) *h += pwd->GapExtPen3(gl);
	    *h += bb[-3].sigE;
	    if (*h < nevsel) break;
	    var_t	x = h[-1] + pwd->GapW1;
	    if (x > *h) {
		*h = x;
		if (c) *c = c[-1];
		if (d) *d = d[-1];
		if (row0) *row0 = static_cast<TBU_t>(TraceBackCode::HOR1);
	    }
	    x = h[-2] + pwd->GapW2;
	    if (x > *h) {
		*h = x;
		if (c) *c = c[-2];
		if (d) *d = d[-2];
		if (row0) *row0 = static_cast<TBU_t>(TraceBackCode::HOR2);
	    }
	    x = (bb->sigS > 0)? bb->sigS: 0;
	    if (x > *h) {
		*h = x;
		lend[p] = r;
		if (vmf) {
		    int	ptr = vmf->add(a->left, n, 0);
		    *c = i_s2(d, ptr);
		    *e = 1;
		} else if (usec)
		    *c = i_s2(d, r);
	    } else if (row0) *row0 = static_cast<TBU_t>(TraceBackCode::HORI);
;	// HORI
	    if (c) ++c;
	    if (d) ++d;
	    if (e) ++e;
	    if (trb) trb->to_right(row0);
	} // end of r-loop
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
int SimdAln2h1<var_t, Nelem, regist_v, regist_m>::
fhlastH1(Rvulmn& maxh, Anti_rhomb_coord<TBU_t>* trb)
{
	int	glen[3] = {0, 0, 0};
	bool	tcdn[3] = {false, false, false};
const	int	m3 = 3 * a->right;
	int	rw = wdw.lw;
	int	rf = b->left - m3;
	if (rf > rw) rw = rf;
	else	rf = rw;
const	int	rr = b->right - m3;
	int	maxr = rr;
	int	maxt = rr;
	var_t*	h = hv + rw;
	var_t*	h9 = hv + maxr;
	VTYPE&	maxv = maxh.val;
const	SGPT6*	bb = b->exin->score_p(rw + m3 - 2);
	TBU_t*	rowM = trb? trb->set_point(a->right, rw + m3): 0;

	maxv = *h9;
	if (a->inex.exgr) {
	    for (int p = 0; h <= h9; ++h, ++rf, ++bb, p = next_p[p]) {
		VTYPE	x = NEVSEL;
		VTYPE	y = NEVSEL;
		if (rf - rw >= 3 && !tcdn[p]) {
		    x = h[-3] + bb->sigE;
		    if (!(a->inex.exgr & 2)) x += pwd->GapExtPen3(glen[p]);
		    if (!(a->inex.exgr & 1) && glen[p] == 3) x += pwd->BasicGOP;
		    if (bb->sigT > 0) y = h[-3] + bb->sigT;
		}
		if (rf - rw >= 3) tcdn[p] = bb->sigT > 0;
		if (*h > x && *h > y) {		/* \ */
		    glen[p] = 0;
		    if (*h > maxv) {
			maxv = *h;
			maxr = maxt = rf;
		    }
		    if (rowM) *rowM = static_cast<TBU_t>(TraceBackCode::DIAG);
		} else if (x >= y) {		// _
		    glen[p] += 3;
		    *h = static_cast<var_t>(x);
		    if (x > maxv) {
			maxv = x;
			maxt = rf;
		    }
		    if (rowM) *rowM = static_cast<TBU_t>(TraceBackCode::HORI);
		} else if (y > maxv) {		// termination codon
		    glen[p] += 3;
		    maxv = y;
		    maxt = rf;
		    if (rowM) *rowM = static_cast<TBU_t>(TraceBackCode::HORI);
		}
		if (rowM) {
		    if (glen[p] == 3)		// \_
			*rowM |= static_cast<TBU_t>(TraceBackCode::NHOR);
		    trb->to_right(rowM);
		}
	    }
	} else {
	    bb += h9 - h;
const	    VTYPE	y = h9[-3] + bb->sigT;
	    if (y > *h9) {
		*h9 = y;
		maxt = b->right - m3;
		maxr = maxt - 3;
	    }
	}
	if (b->inex.exgr) {
	    rw = std::min(wdw.up - 1, b->right - 3 * a->left);
	    for (h = hv + rw; h > h9; --h, --rw) {
		VTYPE	x = *h + (rw % 3? pwd->ExtraGOP: 0);
		if (x > maxv) {
		    maxv = x;
		    maxt = maxr = rw;
		}
	    }
	}
	if (mode == 2 || mode == 4) hb[maxt] = hb[maxr];
	maxh.ulk = (mode > 1)? s2_i(hc[maxr], used? hd[maxr]: 0): maxr;
	int	p = maxr - rr;
	if (vmf) {
	    int	m9 = a->right;
	    int n9 = b->right;
	    if (p > 0) {
		m9 -= (p + 2) / 3;
		if (p %= 3) n9 -= (3 - p); 
	    } else if (p < 0) n9 += p;
	    maxh.ulk = vmf->add(m9, n9, maxh.ulk);
	    if (maxr != maxt)	// terminal gap
		maxh.ulk = vmf->add(a->right, maxt + m3, maxh.ulk);
	} else {
	    if (p > 0)	maxh.mr = (b->right - maxr) / 3;	    
	    else	maxh.nr = maxt + m3;
	}
	return (maxt);
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2h1<var_t, Nelem, regist_v, regist_m>::
from_spj(Queue2<int>& list, const int& m, const int& n, const int& q)
{
	for (int k = list.begin_i(); list.isnt_end(); k = list.next_i()) {
const	    int&	nj = list[k];
const	    int	j = (n - nj) / 3;
const	    int	mj = m + j;
	    if ((mj + 1) < a->right)
		fsjss[j]->get(j, mj, nj, q - 1);
	}
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2h1<var_t, Nelem, regist_v, regist_m>::
to_spj(Queue2<int>& list, const int& m, const int& n, const int& q)
{
	for (int k = list.begin_i(); list.isnt_end(); k = list.next_i()) {
const	    int&	nj = list[k];
const	    int	j = (n - nj) / 3;
const	    int	mj = m + j;
	    if ((mj + 1) < a->right)
		fsjss[j]->put(j, mj, nj, q);
	}
}

// ordinary DP forward algorithm
// assume number of stored traceback records < INT_MAX

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
VTYPE SimdAln2h1<var_t, Nelem, regist_v, regist_m>::
forwardH1(int* pp)
{
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);
const	regist_v	two_v = Add(one_v, one_v);
const	regist_v	ninf_v = Splat(nevsel);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	g1_v = Splat(pwd->GapW1);
const	regist_v	g2_v = Splat(pwd->GapW2);
const	regist_v	g3_v = Splat(pwd->GapW3);
regist_v	hd_v, fd_v, ed_v, qd_v;	// used only when used == true

	fhinitH1();

	VTYPE	accscr = 0;
const	int	md = checkpoint(0);
	int	mc = md + a->left;
	for (int ml = a->left; ml < a->right; ml += Nelem) {
const	    int	j9 = std::min(Nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + 3 * ml);
const	    int	n9 = std::min(b->right, wdw.up + 3 * (ml + j9) + 1) + 3 * j9;
	    int	n0 = n - 3 * j8;
const	    int	mp1 = ml + 1;
	    int	q = (n + 3 * mp1) % 6;
	    vec_set(hv_a[0], nevsel, 12 * Np1 + 3 * Nelem);	// [hv_a..fv_a]
	    vec_clear(hb_a[0], 12 * Np1 + 3 * Nelem);
	    vec_clear(ps_a[0], 6 * Nelem);	// ps_a, pv_a
	    vclear(sm_a, 4 * Np1);
	    if (spj) {
		for (int j = 0; j < j9; ++j)
		    fsjss[j]->reset();
		for (int p = 0; p < 3; ++p) {
		    donor_q[p]->clear();
		    accep_q[p]->clear();
		}
	    }
	    int	r = n - 3 * mp1;

	    for ( ; n < n9; ++n, ++n0, ++r, q = modN<6>(q + 1)) {
const		int	p = q % 3;
const		int	nb = std::max(0, n - b->right + 1);
const		int	kb = (nb - 1) / 3;
const		int	ke = std::min(j9, (n - b->left) / 3);
const		SGPT6*	bb = b->exin->score_p(n);

		if (spj && !nb) {
		    if (isPhs1(bb->phs5)) donor_q[p]->push(n);
		    if (isPhs1(bb->phs3)) accep_q[p]->push(n);
		    cp_a[p][0] = b->exin->good(bb - 2)? bb[-2].sigE: 0;
		}

regist_v	cv_v = Load(cp_a[p]);	// coding potential
		Store(cp_a[p] + 1, cv_v);

//	ordinary insertion
		int	qq = modN<6>(q - 1);
regist_v	hv_v = Load(hv_a[qq] + 1);
regist_v	hc_v = Load(hc_a[qq] + 1);
		if (used) hd_v = Load(hd_a[qq] + 1);
		hv_v = Add(hv_v, g1_v);	// frame shift
		qq = modN<6>(q - 2);
regist_v	qv_v = Load(hv_a[qq] + 1);
regist_v	qc_v = Load(hc_a[qq] + 1);
		if (used) qd_v = Load(hd_a[qq] + 1);
		qv_v = Add(qv_v, g2_v);	// frame shift
regist_m	msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		hc_v = Blend(hc_v, qc_v, msk_m);
		if (used) hd_v = Blend(hd_v, qd_v, msk_m);
		qq = modN<6>(q - 3);
		qv_v = Load(hv_a[qq] + 1);
		qc_v = Load(hc_a[qq] + 1);
		if (used) qd_v = Load(hd_a[qq] + 1);
		qv_v = Add(qv_v, g3_v);	// new insertion
		qv_v = Add(qv_v, cv_v);
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		hc_v = Blend(hc_v, qc_v, msk_m);
		if (used) hd_v = Blend(hd_v, qd_v, msk_m);
		qv_v = Load(ev_a[p]);
		qc_v = Load(ec_a[p]);
		if (used) qd_v = Load(ed_a[p]);
		qv_v = Add(qv_v, ge_v);	// extension
		qv_v = Add(qv_v, cv_v);
		msk_m = Cmp_gt(qv_v, hv_v);
regist_v	ev_v = Blend(qv_v, hv_v, msk_m);
regist_v	ec_v = Blend(qc_v, hc_v, msk_m);
		if (used) ed_v = Blend(qd_v, hd_v, msk_m);
		Store(ev_a[p], ev_v);
		Store(ec_a[p], ec_v);
		if (used) Store(ed_a[p], ed_v);

//	deletion
		fv_a[qq][0] = fv[r + 3];
		fc_a[qq][0] = fc[r + 3];
		if (used) fd_a[qq][0] = fd[r + 3];
		hv_v = Load(fv_a[qq]);
		hc_v = Load(fc_a[qq]);
		if (used) hd_v = Load(fd_a[qq]);
		hv_v = Add(hv_v, ge_v);	// gap extension

		hv_a[qq][0] = hv[r + 3];
		hc_a[qq][0] = hc[r + 3];
		if (used) hd_a[qq][0] = hd[r + 3];
		qv_v = Load(hv_a[qq]);
		qc_v = Load(hc_a[qq]);
		if (used) qd_v = Load(hd_a[qq]);
		qv_v = Add(qv_v, g3_v);	// gap open
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		hc_v = Blend(hc_v, qc_v, msk_m);
		if (used) hd_v = Blend(hd_v, qd_v, msk_m);

		qq = modN<6>(q - 4);		// frame shift
		hv_a[qq][0] = hv[r + 2];
		hc_a[qq][0] = hc[r + 2];
		if (used) hd_a[qq][0] = hd[r + 2];
		qv_v = Load(hv_a[qq]);
		qc_v = Load(hc_a[qq]);
		if (used) qd_v = Load(hd_a[qq]);
		qv_v = Add(qv_v, g2_v);
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		hc_v = Blend(hc_v, qc_v, msk_m);
		if (used) hd_v = Blend(hd_v, qd_v, msk_m);

		qq = modN<6>(q - 5);		// frame shift
		hv_a[qq][0] = hv[r + 1];
		hc_a[qq][0] = hc[r + 1];
		if (used) hd_a[qq][0] = hd[r + 1];
		qv_v = Load(hv_a[qq]);
		qc_v = Load(hc_a[qq]);
		if (used) qd_v = Load(hd_a[qq]);
		qv_v = Add(qv_v, g1_v);
		msk_m = Cmp_gt(hv_v, qv_v);
regist_v	fv_v = Blend(hv_v, qv_v, msk_m);
regist_v	fc_v = Blend(hc_v, qc_v, msk_m);
		if (used) fd_v = Blend(hd_v, qd_v, msk_m);
		Store(fv_a[q] + 1, fv_v);
		Store(fc_a[q] + 1, fc_v);
		if (used) Store(fd_a[q] + 1, fd_v);

//	diagonal match
const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 3 * kb - 2);
		if (nb) vclear(sm_a, Nelem);
		for (int k = kb; k < ke; ++k, ++as, bs -= 3)
		    sm_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(sm_a);

		hv_a[q][0] = hv[r];
		hb_a[q][0] = hb[r];
		hc_a[q][0] = hc[r];
		if (used) hd_a[q][0] = hd[r];
		qv_v = Load(hv_a[q]);
		qc_v = Load(hc_a[q]);
		if (used) qd_v = Load(hd_a[q]);
		hv_v = Add(hv_v, qv_v);
		hv_v = Add(hv_v, cv_v);
		Store(qv_a[p], qv_v);
		Store(qc_a[p], qc_v);
		if (used) Store(qd_a[p], qd_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
		hc_v = Blend(fc_v, qc_v, msk_m);
		if (used) hd_v = Blend(fd_v, qd_v, msk_m);
regist_v	hb_v = Blend(two_v, zero_v, msk_m);

//	best of three
		msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !(e > t)
		hv_v = Blend(ev_v, hv_v, msk_m);
		hc_v = Blend(ec_v, hc_v, msk_m);
		if (used) hd_v = Blend(ed_v, hd_v, msk_m);
		hb_v = Blend(one_v, hb_v, msk_m);
		Store(pv_a[p], hb_v);	// diag: 0, hori: 1, vert: 2

regist_v	qb_v = Load(ps_a[p]);		// post splicing
		qb_v = And(hb_v, qb_v);
		Store(ps_a[p], qb_v);

//	Find optimal path

		if (!Local) {
		    msk_m = Cmp_gt(hv_v, ninf_v);	// clamp underflow
		    hv_v = Blend(hv_v, ninf_v, msk_m);
		} else if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);
		    hb_v = Blend(one_v, hb_v, msk_m);
		    hc_v = Blend(zero_v, hc_v, msk_m);
		    if (used) hd_v = Blend(zero_v, hd_v, msk_m);
		}
		msk_m = Cmp_eq(hb_v, zero_v);	// diag: 1, others: 0
		hb_v = Blend(one_v, zero_v, msk_m);
		qb_v = Load(hb_a[q]);
		qb_v = AndNot(hb_v, qb_v);		// t && !q = nd && di
		Store(qb_a[p], qb_v);
		Store(hv_a[q] + 1, hv_v);
		Store(hb_a[q] + 1, hb_v);
		Store(hc_a[q] + 1, hc_v);
		if (used) Store(hd_a[q] + 1, hd_v);
		if (LocalR) {	// max(H)
		    var_t*	mx = vmax(hv_a[q] + 1, j9);
		    if ((*mx + accscr) > maxh.val) {
			int	k = mx - hv_a[q];
			maxh.val = *mx + accscr;
			maxh.ulk = s2_i(hc_a[q][k], hd_a[q][k]);
			maxh.mr = ml + k;
			maxh.nr = n - 3 * k + 3;
		    }
		}
		int	mk = ml + kb;
		int	nk = n - 3 * (kb + 1);
		for (int k = kb; k < ke; ++k, ++mk, nk -= 3) {
const		    int	kp1 = k + 1;
		    if (qb_a[p][k]) {
			int ptr = s2_i(hc_a[q][kp1], hd_a[q][kp1]);
			ptr = vmf->add(mk, nk, ptr);
			hc_a[q][kp1] = i_s2(hd_a[q] + kp1, ptr);
		    }
		}

//	intron 3' boundary fprward phase
		if (spj && accep_q[p]->remain()) {
		    if (accep_q[p]->head() < n0) accep_q[p]->pull();
		    if (accep_q[p]->remain())
			from_spj(*accep_q[p], ml, n, q);
		}

//	intron 5' boundary forward phase
		if (spj && donor_q[p]->remain()) {
		    if (donor_q[p]->head() < n0) donor_q[p]->pull();
		    if (donor_q[p]->remain())
			to_spj(*donor_q[p], ml, n, q);
		}

//	prepare for next cycle
		int	r0 = r - 6 * j8;
		if (j9 == ke && wdw.lw <= r0 && r0 <= wdw.up) {
		    hv[r0] = hv_a[q][j9];
		    hb[r0] = hb_a[q][j9];
		    hc[r0] = hc_a[q][j9];
		    if (used) hd[r0] = hd_a[q][j9];
		    fv[r0] = fv_a[q][j9];
		    fc[r0] = fc_a[q][j9];
		    if (used) fd[r0] = fd_a[q][j9];
		}
	    }	// end of n loop
	    if (ml == mc) {
		var_t	c = vec_max(hv + wdw.lw - 3, wdw.width);
		int	d = checkpoint(c);
		if (d < md / 2) {	// all-round down score
		    vec_add(hv + wdw.lw - 3, wdw.width, -c);
		    vec_add(fv + wdw.lw - 3, wdw.width, -c);
		    accscr += c;
		    mc += md;
		} else			// postpone
		    mc += d;
	    }
	}	// end of ml loop

	if (!LocalR || maxh.mr == a->right) {
	    fhlastH1(maxh);
	    maxh.val = maxh.val + accscr;
	    if (pp) *pp = maxh.ulk;
	} else if (pp) {
	    *pp = vmf->add(maxh.mr, maxh.nr, maxh.ulk);
	}
	return (maxh.val);
}

// unidirectional Hirschberg method, multi-intermediate version

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
VTYPE SimdAln2h1<var_t, Nelem, regist_v, regist_m>::
hirschbergH1(Dim10* cpos, const int& n_im)
{
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);
const	regist_v	two_v = Add(one_v, one_v);
const	regist_v	ninf_v = Splat(nevsel);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	g1_v = Splat(pwd->GapW1);
const	regist_v	g2_v = Splat(pwd->GapW2);
const	regist_v	g3_v = Splat(pwd->GapW3);

	fhinitH1();

	mm = (a->right - a->left + n_im) / (n_im + 1);
	Udh_Imds	udhimds(n_im, a->left, mm, wdw, pwd->Noll);
	imd = udhimds[0];
	mm = a->left + (imd->mi - a->left - 1) / Nelem * Nelem;
	mm3 = 3 * imd->mi;
	int	k9 = imd->mi - mm;
	int	k8 = k9 - 1;

const	int	md = checkpoint(0);
	int	mc = md + a->left;
	VTYPE	accscr = 0;

regist_v	hb_v, fb_v, eb_v, qb_v;	// used only when LocalL == true
regist_v	hd_v, fd_v, ed_v, qd_v;	// used only when used == true

	for (int ml = a->left, i = 0; ml < a->right; ml += Nelem) {
const	    int	j9 = std::min(Nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + 3 * ml);
	    int	n9 = std::min(b->right, wdw.up + 3 * (ml + j9) + 1) + 3 * j9;
	    int	n0 = n - 3 * j8;
	    int	mp1 = ml + 1;
	    int	q = modN<6>(n + 3 * mp1);
	    int	r = n - 3 * mp1;
	    vec_set(hv_a[0], nevsel, 12 * Np1 + 3 * Nelem);	// [hv_a..fv_a]
	    vec_clear(hb_a[0], 12 * Np1 + 3 * Nelem);
	    vec_clear(ps_a[0], 6 * Nelem);	// ps_a, pv_a
	    vclear(sm_a, 4 * Np1);
const	    bool	is_imd_ = ml == mm;

	    if (spj) {
		for (int j = 0; j < j9; ++j)
		    fsjss[j]->reset();
		for (int p = 0; p < 3; ++p) {
		    donor_q[p]->clear();
		    accep_q[p]->clear();
		}
	    }

	    for ( ; n < n9; ++n, ++n0, ++r, q = modN<6>(q + 1)) {
const		int	rj = r - 6 * k8;
const		int	p = q % 3;
const		int	nb = std::max(0, n - b->right + 1);
const		int	kb = (nb - 1) / 3;
const		int	ke = std::min(j9, (n - b->left) / 3);
const		SGPT6*	bb = b->exin->score_p(n);
const		bool	is_imd = is_imd_ && rj >= wdw.lw && rj <= wdw.up;

		if (spj && !nb) {
		    if (isPhs1(bb->phs5)) donor_q[p]->push(n);
		    if (isPhs1(bb->phs3)) accep_q[p]->push(n);
		    cp_a[p][0] = b->exin->good(bb - 2)? bb[-2].sigE: 0;
		}

regist_v	cv_v = Load(cp_a[p]);
		Store(cp_a[p] + 1, cv_v);

//	ordinary insertion
		int	qq = modN<6>(q - 1);
regist_v	hv_v = Load(hv_a[qq] + 1);
		if (LocalL) hb_v = Load(hb_a[qq] + 1);
regist_v	hc_v = Load(hc_a[qq] + 1);
		if (used) hd_v = Load(hd_a[qq] + 1);
		hv_v = Add(hv_v, g1_v);	// frame shift
		qq = modN<6>(q - 2);
regist_v	qv_v = Load(hv_a[qq] + 1);
		if (LocalL) qb_v = Load(hb_a[qq] + 1);
regist_v	qc_v = Load(hc_a[qq] + 1);
		if (used) qd_v = Load(hd_a[qq] + 1);
		qv_v = Add(qv_v, g2_v);	// frame shift
regist_m	msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		if (LocalL) hb_v = Blend(hb_v, qb_v, msk_m);
		hc_v = Blend(hc_v, qc_v, msk_m);
		if (used) hd_v = Blend(hd_v, qd_v, msk_m);
		qq = modN<6>(q - 3);
		qv_v = Load(hv_a[qq] + 1);
		if (LocalL) qb_v = Load(hb_a[qq] + 1);
		qc_v = Load(hc_a[qq] + 1);
		if (used) qd_v = Load(hd_a[qq] + 1);
		qv_v = Add(qv_v, g3_v);	// new insertion
		qv_v = Add(qv_v, cv_v);
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		if (LocalL) hb_v = Blend(hb_v, qb_v, msk_m);
		hc_v = Blend(hc_v, qc_v, msk_m);
		if (used) hd_v = Blend(hd_v, qd_v, msk_m);
		qv_v = Load(ev_a[p]);
		if (LocalL) qb_v = Load(eb_a[p]);
		qc_v = Load(ec_a[p]);
		if (used) qd_v = Load(ed_a[p]);
		qv_v = Add(qv_v, ge_v);	// extension
		qv_v = Add(qv_v, cv_v);
		msk_m = Cmp_gt(qv_v, hv_v);
regist_v	ev_v = Blend(qv_v, hv_v, msk_m);
		if (LocalL) eb_v = Blend(qb_v, hb_v, msk_m);
regist_v	ec_v = Blend(qc_v, hc_v, msk_m);
		if (used) ed_v = Blend(qd_v, hd_v, msk_m);
		Store(ev_a[p], ev_v);
		if (LocalL) Store(eb_a[p], eb_v);
		Store(ec_a[p], ec_v);
		if (used) Store(ed_a[p], ed_v);

//	deletion
		fv_a[qq][0] = fv[r + 3];
		if (LocalL) fb_a[qq][0] = fb[r + 3];
		fc_a[qq][0] = fc[r + 3];
		if (used) fd_a[qq][0] = fd[r + 3];
		hv_v = Load(fv_a[qq]);
		if (LocalL) hb_v = Load(fb_a[qq]);
		hc_v = Load(fc_a[qq]);
		if (used) hd_v = Load(fd_a[qq]);
		hv_v = Add(hv_v, ge_v);	// gap extension

		hv_a[qq][0] = hv[r + 3];
		if (LocalL) hb_a[qq][0] = hb[r + 3];
		hc_a[qq][0] = hc[r + 3];
		if (used) hd_a[qq][0] = hd[r + 3];
		qv_v = Load(hv_a[qq]);
		if (LocalL) qb_v = Load(hb_a[qq]);
		qc_v = Load(hc_a[qq]);
		if (used) qd_v = Load(hd_a[qq]);
		qv_v = Add(qv_v, g3_v);	// gap open
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		if (LocalL) hb_v = Blend(hb_v, qb_v, msk_m);
		hc_v = Blend(hc_v, qc_v, msk_m);
		if (used) hd_v = Blend(hd_v, qd_v, msk_m);

		qq = modN<6>(q - 4);		// frame shift
		hv_a[qq][0] = hv[r + 2];
		if (LocalL) hb_a[qq][0] = hb[r + 2];
		hc_a[qq][0] = hc[r + 2];
		if (used) hd_a[qq][0] = hd[r + 2];
		qv_v = Load(hv_a[qq]);
		if (LocalL) qb_v = Load(hb_a[qq]);
		qc_v = Load(hc_a[qq]);
		if (used) qd_v = Load(hd_a[qq]);
		qv_v = Add(qv_v, g2_v);
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		if (LocalL) hb_v = Blend(hb_v, qb_v, msk_m);
		hc_v = Blend(hc_v, qc_v, msk_m);
		if (used) hd_v = Blend(hd_v, qd_v, msk_m);

		qq = modN<6>(q - 5);		// frame shift
		hv_a[qq][0] = hv[r + 1];
		if (LocalL) hb_a[qq][0] = hb[r + 1];
		hc_a[qq][0] = hc[r + 1];
		if (used) hd_a[qq][0] = hd[r + 1];
		qv_v = Load(hv_a[qq]);
		if (LocalL) qb_v = Load(hb_a[qq]);
		qc_v = Load(hc_a[qq]);
		if (used) qd_v = Load(hd_a[qq]);
		qv_v = Add(qv_v, g1_v);
		msk_m = Cmp_gt(hv_v, qv_v);
regist_v	fv_v = Blend(hv_v, qv_v, msk_m);
		if (LocalL) fb_v = Blend(hb_v, qb_v, msk_m);
regist_v	fc_v = Blend(hc_v, qc_v, msk_m);
		if (used) fd_v = Blend(hd_v, qd_v, msk_m);
		Store(fv_a[q] + 1, fv_v);
		if (LocalL) Store(fb_a[q] + 1, fb_v);
		Store(fc_a[q] + 1, fc_v);
		if (used) Store(fd_a[q] + 1, fd_v);

//	diagonal match
const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 3 * kb - 2);
		if (nb) vclear(sm_a, Nelem);
		for (int k = kb; k < ke; ++k, ++as, bs -= 3)
		    sm_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(sm_a);

		hv_a[q][0] = hv[r];
		if (LocalL) hb_a[q][0] = hb[r];
		hc_a[q][0] = hc[r];
		if (used) hd_a[q][0] = hd[r];
		qv_v = Load(hv_a[q]);
		if (LocalL) qb_v = Load(hb_a[q]);
		qc_v = Load(hc_a[q]);
		if (used) qd_v = Load(hd_a[q]);
		hv_v = Add(hv_v, qv_v);
		hv_v = Add(hv_v, cv_v);
		Store(qv_a[p], qv_v);
		if (LocalL) Store(qb_a[p], qb_v);
		Store(qc_a[p], qc_v);
		if (used) Store(qd_a[p], qd_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
		if (LocalL) hb_v = Blend(fb_v, qb_v, msk_m);
		hc_v = Blend(fc_v, qc_v, msk_m);
		if (used) hd_v = Blend(fd_v, qd_v, msk_m);
		fb_v = Blend(two_v, zero_v, msk_m);

//	best of three
		msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !e > t)
		hv_v = Blend(ev_v, hv_v, msk_m);
		if (LocalL) hb_v = Blend(eb_v, hb_v, msk_m);
		hc_v = Blend(ec_v, hc_v, msk_m);
		if (used) hd_v = Blend(ed_v, hd_v, msk_m);
		fb_v = Blend(one_v, fb_v, msk_m);
		Store(pv_a[p], fb_v);		// diag: 0, hori: 1, vert: 2

		qb_v = Load(ps_a[p]);		// post splicing
		qb_v = And(fb_v, qb_v);
		Store(ps_a[p], qb_v);

//	Find optimal path

		if (!Local) {
		    msk_m = Cmp_gt(hv_v, ninf_v);	// clamp underflow
		    hv_v = Blend(hv_v, ninf_v, msk_m);
		} else if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);
		}
		Store(hv_a[q] + 1, hv_v);
		if (LocalL) Store(hb_a[q] + 1, hb_v);
		Store(hc_a[q] + 1, hc_v);
		if (used) Store(hd_a[q] + 1, hd_v);
		if (LocalL && !accscr) {		// left end
		    for (int k = kb; k < ke; ++k) {
const			int	kp1 = k + 1;
			if (hv_a[q][kp1] == 0) {
			    hb_a[q][kp1] = static_cast<var_t>(ml + k);
			    hc_a[q][kp1] = i_s2(hd_a[q] + kp1, r - 6 * k);
			}
		    }
		}
		if (LocalR) {			// right end
		    var_t*	mx = vmax(hv_a[q] + 1, j9);
		    if (*mx + accscr > maxh.val) {
			int	k = mx - hv_a[q];
			maxh.val = *mx + accscr;
			maxh.ml = hb_a[q][k];	// ml
			maxh.ulk = s2_i(hc_a[q][k], hd_a[q][k]);
			maxh.mr = ml + k + 1;
			maxh.nr = n - 3 * k;
		    }
		}

//	intron 3' boundary fprward phase
		if (spj && accep_q[p]->remain()) {
		    if (accep_q[p]->head() < n0) accep_q[p]->pull();
		    if (accep_q[p]->remain())
			from_spj(*accep_q[p], ml, n, q);
		}

//	intron 5' boundary forward phase
		if (spj && donor_q[p]->remain()) {
		    if (used) Store(qd_a[p], qd_v);
		    if (donor_q[p]->head() < n0) donor_q[p]->pull();
		    if (donor_q[p]->remain())
			to_spj(*donor_q[p], ml, n, q);
		}

//	intermediate row
		if (is_imd) {
		    if (pv_a[p][k8] == 0) rlst[p] = rj;		// diag
		    if (pv_a[p][k8] == 1) 			// hori
			imd->hlnk[0][rj] = rlst[p];
		    imd->vlnk[0][rj] = s2_i(hc_a[q][k9], hd_a[q][k9]);
		    hc_a[q][k9] = i_s2(hd_a[q] + k9, rj);
		    imd->vlnk[1][rj] = s2_i(fc_a[q][k9], fd_a[q][k9]);
		    fc_a[q][k9] = i_s2(fd_a[q] + k9, rj + wdw.width);
		}

//	prepare for next cycle
		int	r0 = r - 6 * j8;
		if (j9 == ke && wdw.lw <= r0 && r0 < wdw.up) {
		    hv[r0] = hv_a[q][j9];
		    if (LocalL) hb[r0] = hb_a[q][j9];
		    hc[r0] = hc_a[q][j9];
		    if (used) hd[r0] = hd_a[q][j9];
		    fv[r0] = fv_a[q][j9];
		    if (LocalL) fb[r0] = fb_a[q][j9];
		    fc[r0] = fc_a[q][j9];
		    if (used) fd[r0] = fd_a[q][j9];
		}
	    }	// end of n loop
	    if (ml == mc) {
		var_t	c = vec_max(hv + wdw.lw - 3, wdw.width);
		int	d = checkpoint(c);
		if (d < md / 2) {	// down score all
		    vec_add(hv + wdw.lw - 3, wdw.width, -c);
		    vec_add(fv + wdw.lw - 3, wdw.width, -c);
		    accscr += c;
		    mc += md;
		} else			// postpone
		    mc += d;
	    }
	    if (is_imd_ && ++i < n_im) {	// reset intermediate
		imd = udhimds[i];
		mm3 = 3 * imd->mi;
		mm = a->left + (imd->mi - a->left - 1) / Nelem * Nelem;
		k9 = imd->mi - mm;
		k8 = k9 - 1;
	    }
	}	// end of ml loop

	if (LocalR && maxh.mr < a->right) {
	    a->right = maxh.mr;
	    b->right = maxh.nr;
	} else {	// global or semi-global
const	    int	r = fhlastH1(maxh);
	    maxh.val += accscr;
	    maxh.ml = LocalL? hb[r]: a->left;
	    a->right = maxh.mr;
	    b->right = maxh.nr;
	}

	int	i = n_im;
	while (--i >= 0 && udhimds[i]->mi > a->right) ;
	if (i < 0 && udhimds[0]->mi > a->right) cpos[0][2] = b->right;
	int	r = maxh.ulk;
	for ( ; i >= 0 && (imd = udhimds[i])->mi > maxh.ml; --i) {
	    int	c = 0;
	    int	d = 0;
	    for ( ; r > wdw.up; r -= wdw.width) ++d;
	    if (imd->vlnk[d][r] < end_of_ulk) {	// cross intermediate
		cpos[i][c++] = imd->mi;
		cpos[i][c++] = (d > 0)? 1: 0;
		mm3 = 3 * imd->mi;
		for (int rp = imd->hlnk[d][r]; rp < end_of_ulk && r != rp;
		    rp = imd->hlnk[d][r = rp]) {
		    cpos[i][c++] = r + mm3;
		}
		cpos[i][c++] = r + mm3;
		cpos[i][c] = end_of_ulk;
		r = imd->vlnk[d][r];
		if (r == end_of_ulk) break;
	    } else
		cpos[i][0] = end_of_ulk;	// don't cross intermediate
	}
	for (int d = 0; r > wdw.up; r -= wdw.width) ++d;
	if (LocalL) {
	    a->left = maxh.ml;
	    b->left = r + 3 * a->left;
	} else {
const	    int	rl = b->left - 3 * a->left;
	    if (b->inex.exgl && rl > r) {
		a->left = (b->left - r) / 3;
		for (int j = 0; j < n_im && udhimds[j]->mi < a->left; ++j)
		    cpos[j][0] = end_of_ulk;
	    }
	    if (a->inex.exgl && rl < r) b->left = 3 * a->left + r;
	}
	if (udhimds[n_im - 1]->mi < a->left)
	    cpos[0][2] = b->left;
	return (maxh.val);
}

#include "fwd2h1_wip_simd.h"

#endif	// _FWD2H1_SIMD_H_
