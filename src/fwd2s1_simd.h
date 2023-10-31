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

#if _SIMD_PP_
#include <simdpp/simd.h>
#endif	// _SIMD_PP_

#include "simd_functions.h"
#include "udh_intermediate.h"

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

template <typename var_t>
struct Rvlujd {
	var_t	val;	// score
	int	ulk;	// lower end diagonal
	int	jnc;	// junction
	int	ml;	// left end m-coord
	short	dir;	// direction
};

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
#if _SIMD_PP_
class SimdAln2s1 {
#else
class SimdAln2s1 : public Simd_functions<var_t, Nelem, regist_v> {
#endif	// _SIMD_PP_

protected:
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
	SimdAln2s1& hs1;
	Rvlujd<var_t>*	rcd;
	short*	idx;
	short	ncand;
	Rvlujd<var_t>	black_r;
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
	    rcd = new Rvlujd<var_t>[hs1.Ncand + 1];
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
	void	vec_add(var_t* ar, const int& n, const var_t& c);
	var_t	vec_max(const var_t* ar, const int& n);
	int	checkpoint(const var_t& pv) const {
	    return ((sizeof(var_t) > 2)? INT_MAX:
	    (int((check_scr - abs(pv)) / avmch / Nelem * Nelem)));
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
	    Np1(Nelem + 1), avmch(pwd->simmtx->AvTrc()), 
#if FVAL
	    nevsel(NEVSEL),
#else
	    nevsel((sizeof(var_t) == 2? SHRT_MIN: INT_MIN) + nevsel_play),
#endif
	    black_Rvulmn{nevsel, end_of_ulk, a->left, a->right, b->right},
	    rlst(INT_MAX), buf_size(wdw.width + 2 * Nelem)
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
#define Or(a, b)	this->bit_or(a, b)
#define AndNot(a, b)	this->bit_andnot(a, b)
#define To_mask(a)	this->load(a)
#define Cast16to8(a)	this->cast16to8(a)
#endif	// _SIMD_PP_

/*****************************************************************
*	mode	v	b	c	d	o	
*	  0	4 + 2	0 + 0	0 + 0	0 + 0	2 + 2	score only
*	  1	4 + 2	4 + 2	0 + 0	0 + 0	2 + 2	bimap
*	  2	4 + 2	4 + 2	4 + 2	0 + 0	2 + 2	udh1
*	  3	4 + 2	4 + 2	4 + 2	0 + 0	2 + 4	vmf1
*	  4	4 + 2	4 + 2	4 + 2	4 + 2	2 + 2	udh2
*	  5	4 + 2	4 + 2	4 + 2	4 + 2	2 + 4	vmf2
******************************************************************/

	    size_t		bufsiz =  10 * Np1 + 6 * Nelem;
	    if (mode > 1)	bufsiz += (8 * Np1 + 6 * Nelem);
	    abuf = new var_t[bufsiz];
	    ps_a = abuf;
	    pv_a = ps_a + Nelem;
	    hv_a[0] = pv_a + Nelem;
	    hv_a[1] = hv_a[0] + Np1;
	    ev_a = hv_a[1] + Np1;
	    fv_a = ev_a + Nelem;
	    ev2_a = fv_a + Np1;
	    fv2_a = ev2_a + Nelem;
	    s5_a = fv2_a + Np1;
	    s3_a = s5_a + Np1;
	    hb_a[0] = s3_a + Np1;
	    hb_a[1] = hb_a[0] + Np1;
	    eb_a = hb_a[1] + Np1;
	    fb_a = eb_a + Nelem;
	    eb2_a = fb_a + Np1;
	    fb2_a = eb2_a + Nelem;
	    if (mode > 1) {
		pb_a = fb2_a + Np1;
		qb_a = pb_a + Nelem;
		hc_a[0] = qb_a + Nelem;
		hc_a[1] = hc_a[0] + Np1;
		ec_a = hc_a[1] + Np1;
		fc_a = ec_a + Nelem;
		ec2_a = fc_a + Np1;
		fc2_a = ec2_a + Nelem;
		hd_a[0] = fc2_a + Np1;
		hd_a[1] = hd_a[0] + Np1;
		ed_a = hd_a[1] + Np1;
		fd_a = ed_a + Nelem;
		ed2_a = fd_a + Np1;
		fd2_a = ed2_a + Nelem;
	    }

	    for (int p = 0; p < 2; ++p) {
		for (int j = 0; j < Nelem; ++j) {
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
		fsjss = new Sjsites*[Nelem];
		for (int k = 0; k < Nelem; ++k)
		    fsjss[k] = new Sjsites(*this);
		donor_q = new Queue2<int>(Nelem);
		accep_q = new Queue2<int>(Nelem);
	    }
	}
	~SimdAln2s1() {
	    delete[] abuf; delete[] vbuf;
	    if (fsjss) {
		for (int k = 0; k < Nelem; ++k) delete fsjss[k];
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

/*************************************************************************
	nested class: Sjsites
*************************************************************************/

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
int SimdAln2s1<var_t, Nelem, regist_v, regist_m>::Sjsites::
get(const int& j, const int& m, const int& n, const int& p)
{
const	int	r = n - m;
	if (r < hs1.wdw.lw || r >= hs1.wdw.up) return (0);
	var_t**	hfesmv = hs1.hfesv[p][j];
	var_t**	hfesmb = hs1.hfesb[p][j];
	var_t**	hfesmc = hs1.hfesc[p][j];
	var_t**	hfesmd = hs1.hfesd[p][j];
	int	rv = 0;
const	VTYPE	sigJ = hs1.cip? hs1.cip->cip_score(m): 0;
const	Rvlujd<var_t>*	brd = 0;
const	Rvlujd<var_t>*	maxprd[hs1.nod];
const	bool	is_imd = hs1.imd && m == hs1.imd->mi;
	vclear(maxprd, hs1.nod);
	for (int l = 0; l <= ncand; ++l) {
const	    Rvlujd<var_t>*	prd = rcd + idx[l];
const	    Rvlujd<var_t>*&	mrd = maxprd[prd->dir];
const	    int&	don = prd->jnc;
const	    int&	d = prd->dir;
	    if ((n - don) < IntronPrm.minl) continue;
const	    VTYPE	x = prd->val + sigJ + hs1.spjcs->spjscr(don, n);
	    if (x <= *hfesmv[d]) continue;
	    if (!mrd || x > mrd->val) {
		mrd = prd;
		if (!brd || x > brd->val) brd = prd;
	    }
	    *hfesmv[d] = static_cast<var_t>(x);
	    *hfesmv[PS_A] |= psp_bit[d];	// mark spliced
	    if (d & 1) rv |= ((d + 1) / 2);	// hori
	    if (hs1.vmf) {
		*hfesmb[d] = prd->ml;
		int ptr = hs1.vmf->add(m, n, 
		    hs1.vmf->add(m, prd->jnc, prd->ulk));
		*hfesmc[d] = hs1.i_s2(hfesmd[d], ptr);
	    } else {
		*hfesmb[d] = prd->ml;
		if (hs1.mode > 1) *hfesmc[d] = hs1.i_s2(hfesmd[d], prd->ulk);
	    }
	    if (d && *hfesmv[d] > *hfesmv[0]) {
		*hfesmv[0] = *hfesmv[d];
		*hfesmb[0] = *hfesmb[d];
		if (hs1.mode > 1) {
		    *hfesmc[0] = *hfesmc[d];
		    *hfesmd[0] = *hfesmd[d];
		}
	    }
	}
	if (is_imd && brd) {
const	    int&	maxd = brd->dir;
const	    Rvlujd<var_t>*    prd = maxprd[maxd];
	    hs1.imd->hlnk[0][hs1.rlst = r] = prd->ulk;
	    if (hs1.mode > 1) *hfesmc[maxd] = hs1.i_s2(hfesmd[maxd], r);
	    *hfesmv[PV_A] = maxd;
	    if (maxd != 0) {
		*hfesmc[0] = *hfesmc[maxd];
		*hfesmd[0] = *hfesmd[maxd];
		return (rv);
	    }
	    for (int k = 1; k < hs1.pwd->Noll; ++k) {
		int	d = 2 * k - 1;
		if (maxprd[d] && 
		    *hfesmv[d] > (*hfesmv[0] + hs1.pwd->GOP[k])) {
		    *hfesmc[d] = hs1.i_s2(hfesmd[d], r + k * hs1.wdw.width);
		}
		if (maxprd[++d] && 
		    *hfesmv[d] > (*hfesmv[0] + hs1.pwd->GOP[k])) {
		    hs1.imd->hlnk[k][r] = prd->ulk;
		    *hfesmc[d] = hs1.i_s2(hfesmd[d], r + k * hs1.wdw.width);
		}
	    }
	}
	return (rv);
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2s1<var_t, Nelem, regist_v, regist_m>::Sjsites::
put(const int& j, const int& m, const int& n, const int& p)
{
const	int	r = n - m;
	if (r < hs1.wdw.lw || r >= hs1.wdw.up) return;
	var_t**	hfesmv = hs1.hfesv[p][j];
	var_t**	hfesmb = hs1.hfesb[p][j];
	var_t**	hfesmc = hs1.hfesc[p][j];
	var_t**	hfesmd = hs1.hfesd[p][j];

const	int	h = *hfesmv[PV_A];
const	bool	is_imd = hs1.imd && m == hs1.imd->mi;
const	VTYPE	sigJ = hs1.b->exin->sig53(n, n, IE5);
	for (int k = h? 1: 0; k < hs1.nod; ++k) {
	    if (*hfesmv[PS_A] & psp_bit[k]) continue;	// prevent orphan exon
const	    VTYPE	from = *hfesmv[k];
	    if (k) {
const	        VTYPE	thrscr = *hfesmv[0] + hs1.pwd->GOP[(k + 1) / 2];
	        if (from <= thrscr) continue;		// prune
	    }
const	    VTYPE	x = from + sigJ;
	    if (x <= hs1.nevsel) continue;
	    short	l = ncand < hs1.Ncand? ++ncand: hs1.Ncand;
	    while (--l >= 0) {
		if (x >= rcd[idx[l]].val)
		    std::swap(idx[l], idx[l + 1]);
		else
		    break;
	    }
	    if (++l < hs1.Ncand) {
		Rvlujd<var_t>*	prd = rcd + idx[l];
		prd->val = static_cast<var_t>(x);
		prd->ml = *hfesmb[k];
const		int	rl = hs1.mode? hs1.s2_i(*hfesmc[k], *hfesmd[k]): 0;
		if (is_imd) {
		    if (k & 1)	hs1.imd->hlnk[0][r] = hs1.rlst;
		    else if (r != rl) hs1.imd->vlnk[k / 2][r] = rl;
		    prd->ulk = r;
		} else if (hs1.mode) {
		    prd->ulk = rl;
		}
		prd->jnc = n;
		prd->dir = k;
	    } else --ncand;
	}
}

/*************************************************************************
	normal member functions of SimdAln2s1
*************************************************************************/

// add const c to each element of array ar of size n
template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
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
var_t SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
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
void SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
fhinitS1()
{
const	int	ru = wdw.up + 2 * Nelem;
const	int	rl = b->left - a->left;
	if (b->inex.exgl)	// vertical semi-global
	    vec_clear(hv + wdw.lw, rl - wdw.lw);
	int	rr = b->right - a->left;
	if (wdw.up < rr) rr = wdw.up;
	if (a->inex.exgl)	// semi-global
	    vec_clear(hv + rl, rr - rl + 1);
	else {			// global
	    int	r = rl;
	    hv[r++] = 0;
	    hv[r] = pwd->GapPenalty(1);
	    if (pwd->BasicGEP) {
		int	x = (nevsel - pwd->BasicGOP) / pwd->BasicGEP + rl;
		if (x < rr) rr = x;
		while (++r < rr)
		    hv[r] = hv[r - 1] + pwd->BasicGEP;
	    } else if (rr > r)
		vec_set(hv + r, hv[r], rr - r);
	}

	if (vmf) {	// ordinary traceback
	    vec_clear(bbuf, pwd->Noll * buf_size);	// non_diagonal
	    int	ptr = vmf->add(0, 0, 0);		// dummy
	    ptr = vmf->add(a->left, b->left, 0);
	    hc[rl] = i_s2(used? hd + rl: 0, ptr);
	    if (a->inex.exgl) {
		vec_clear(hc + rl + 1, ru - rl);
		if (used) vec_clear(hd + rl + 1, ru - rl);
	    } else {
		vec_set(hc + rl + 1, hc[rl], ru - rl);
		if (used) vec_set(hd + rl + 1, hd[rl], ru - rl);
	    }
	    if (b->inex.exgl) {
		vec_clear(hc + wdw.lw - 1, rl - wdw.lw + 1);
		if (used) vec_clear(hd + wdw.lw - 1, rl - wdw.lw + 1);
	    } else {
		vec_set(hc + wdw.lw - 1, hc[rl], rl - wdw.lw + 1);
		if (used) vec_set(hd + wdw.lw -1, hd[rl], rl - wdw.lw + 1);
	    }
	} else if (mode > 1) {			// Hirschberg
	    vec_set(bbuf, var_t(a->left), pwd->Noll * buf_size);
	    hc[rl] = i_s2(used? hd + rl: 0, rl);
	    if (a->inex.exgl) {		// horizontal semi-global
		int	r = rl;
		var_t*	c = hc + r;
		var_t*	d = hd? hd + r: 0;
		while (r < ru)
		    *c++ = i_s2(d? d++: 0, r++);
	    } else {
		vec_set(hc + rl + 1, hc[rl], ru - rl);
		if (used) vec_set(hd + rl + 1, hd[rl], ru - rl);
	    }
	    if (b->inex.exgl) {		// vertical semi-global
		int	r = wdw.lw - 1;
		var_t	m = var_t(a->left);
		var_t*	e = hb + rl;
		var_t*	c = hc + r;
		var_t*	d = hd? hd + r: 0;
		while (r < rl) {
		    *c++ = i_s2(d? d++: 0, r++);
		    *e-- = m++;
		}
	    } else {
		vec_set(hc + wdw.lw - 1, hc[rl], rl - wdw.lw + 1);
		if (used) vec_set(hd + wdw.lw - 1, hd[rl], rl - wdw.lw + 1);
	    }
	}
	if (mode > 1) {
	    vec_copy(fc + wdw.lw - 1, hc + wdw.lw - 1, buf_size);
	    if (used) vec_copy(fd + wdw.lw - 1, hd + wdw.lw - 1, buf_size);
	    if (fc2) vec_copy(fc2 + wdw.lw - 1, hc + wdw.lw - 1, buf_size);
	    if (fd2) vec_copy(fd2 + wdw.lw - 1, hd + wdw.lw - 1, buf_size);
	}
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
int SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
fhlastS1(Rvulmn& maxh)
{
const	int	rr = b->right - a->right;
	int	maxr = rr;

	if (a->inex.exgr) {
	    int	r = std::max(wdw.lw, b->left - a->right);
	    maxr = vmax(hv + r, rr - r) - hv;
	}
	if (b->inex.exgr) {
	    int	r = std::min(wdw.up - 1, b->right - a->left);
	    int	max_vert = vmax(hv + rr, r - rr) - hv;
	    if (hv[max_vert] > hv[maxr]) maxr = max_vert;
	}
	maxh.val = hv[maxr];
	if (maxr > rr)	maxh.mr = b->right - maxr;
	else		maxh.nr = a->right + maxr;
	if (mode > 1) maxh.ulk = s2_i(hc[maxr], used? hd[maxr]: 0);
	if (vmf)
	    maxh.ulk = vmf->add(maxh.mr, maxh.nr, maxh.ulk);
	return (maxr);
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
int SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
from_spj(Queue2<int>& list, const int& m, const int& n, const int& p)
{
	int	rv = 0;
	for (int k = list.begin_i(); list.isnt_end(); k = list.next_i()) {
const	    int&	nj = list[k];
const	    int	j = n - nj;
	    rv |= fsjss[j]->get(j, m + j, nj, p);
	}
	return (rv);
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
void SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
to_spj(Queue2<int>& list, const int& m, const int& n, const int& p)
{
	for (int k = list.begin_i(); list.isnt_end(); k = list.next_i()) {
const	    int&	nj = list[k];
const	    int	j = n - nj;
	    fsjss[j]->put(j, m + j, nj, p);
	}
}

// score-only DP forward algorithm

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
VTYPE SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
scoreonlyS1()
{
	if (algmode.alg > 1) return (scoreonlyS1_wip());
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);
const	regist_v	two_v = Add(one_v, one_v);
const	regist_v	ninf_v = Splat(nevsel);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	gn_v = Splat(pwd->BasicGEP + pwd->BasicGOP);
const	regist_v	ge2_v = Splat(pwd->LongGEP);
const	regist_v	gn2_v = Splat(pwd->LongGEP + pwd->LongGOP);

	fhinitS1();

	VTYPE	accscr = 0;
const	int	md = checkpoint(0);
	int	mc = md + a->left;
	for (int ml = a->left; ml < a->right; ml += Nelem) {
const	    int	j9 = std::min(Nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + ml);
const	    int	n9 = std::min(b->right, wdw.up + (ml + j9) + 1) + j9;
	    int	n0 = n - j8;
	    int	r = n - (ml + 1);
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * Nelem);
	    vec_clear(ps_a, 2 * Nelem);

regist_v    ev_v = ninf_v;
regist_v    ev2_v = ninf_v;

	    if (spj) {
		for (int j = 0; j < j9; ++j)
		    fsjss[j]->reset();
		donor_q->clear();
		accep_q->clear();
	    }

	    for (int p = 0; n < n9; ++n, ++n0, ++r, p = 1 - p) {
regist_v	fv_v, hb_v;
const		int	q = 1 - p;
const		int	r0 = r - 2 * j8;
const		int	kb = std::max(0, n - b->right);
const		int	ke = std::min(j9, n - b->left);

		if (spj && !kb) {
		    if (b->exin->isDonor(n)) donor_q->push(n);
		    if (b->exin->isAccpt(n)) accep_q->push(n);
		}

//	ordinary insertion
regist_v	qv_v = Load(hv_a[q] + 1);
regist_v	hv_v = Add(qv_v, gn_v);	// new insertion
		ev_v = Add(ev_v, ge_v);	// extension
regist_m	msk_m = Cmp_gt(ev_v, hv_v);
		ev_v = Blend(ev_v, hv_v, msk_m);
		if (dagp) {
		    hv_v = Add(qv_v, gn2_v);	// new insertion
		    ev2_v = Add(ev2_v, ge2_v);	// extension
		    msk_m = Cmp_gt(ev2_v, hv_v);
		    ev2_v = Blend(ev2_v, hv_v, msk_m);
		}
		if (spj) {
		    Store(ev_a, ev_v);
		    if (dagp) Store(ev2_a, ev2_v);
		}

//	deletion
		hv_a[q][0] = hv[r + 1];
		qv_v = Load(hv_a[q]);
		hv_v = Add(qv_v, gn_v);	// gap open
		fv_a[0] = fv[r + 1];
		fv_v = Load(fv_a);
		fv_v = Add(fv_v, ge_v);	// gap extension
		msk_m = Cmp_gt(fv_v, hv_v);
		fv_v = Blend(fv_v, hv_v, msk_m);
		Store(fv_a + 1, fv_v);
		if (dagp) {
		    hv_v = Add(qv_v, gn2_v);	// gap open
		    fv2_a[0] = fv2[r + 1];
regist_v	    fv2_v = Load(fv2_a);
		    fv2_v = Add(fv2_v, ge2_v);	// gap extension
		    msk_m = Cmp_gt(fv2_v, hv_v);
		    fv2_v = Blend(fv2_v, hv_v, msk_m);
		    Store(fv2_a + 1, fv2_v);
		    msk_m = Cmp_gt(fv2_v, fv_v);
		    fv_v = Blend(fv2_v, fv_v, msk_m);
		}

//	diagonal match

const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 1 - kb);
		Store(pv_a, zero_v);
		for (int k = kb; k < ke; ++k, ++as, --bs)
		    pv_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(pv_a);

		hv_a[p][0] = hv[r];
		qv_v = Load(hv_a[p]);
		hv_v = Add(hv_v, qv_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
		hb_v = Blend(two_v, zero_v, msk_m);

//	best of three
		if (dagp) {
		    msk_m = Cmp_gt(ev2_v, ev_v);
		    qv_v = Blend(ev2_v, ev_v, msk_m);
		    msk_m = Cmp_gt(qv_v, hv_v);
		    hv_v = Blend(qv_v, hv_v, msk_m);
		} else {
		    msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !(e > t)
		    hv_v = Blend(ev_v, hv_v, msk_m);
		}
		hb_v = Blend(one_v, hb_v, msk_m);
		Store(pv_a, hb_v);	// diag: 0, hori: 1, vert: 2

		if (spj) {
regist_v	    qb_v = Load(ps_a);		// post splicing
		    qb_v = And(hb_v, qb_v);
		    Store(ps_a, qb_v);
		}

		if (!Local) {
		    msk_m = Cmp_gt(hv_v, ninf_v);	// clamp underflow
		    hv_v = Blend(hv_v, ninf_v, msk_m);
		} else if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);	// Kadane-Gries
		}

//	Find optimal path

		Store(hv_a[p] + 1, hv_v);
		if (LocalR) {	// max(H)
const		    var_t*	mx = vmax(hv_a[p] + 1, j9);
		    if ((*mx + accscr) > maxh.val)
			maxh.val = *mx + accscr;
		}

//	intron 3' boundary
		if (spj && accep_q->remain()) {
		    if (accep_q->head() < n0) accep_q->pull();
		    if (accep_q->remain()) {
const			int	rv = from_spj(*accep_q, ml + 1, n, p);
			if (rv & 1) ev_v = Load(ev_a);
			if (rv & 2) ev2_v = Load(ev2_a);
		    }
		}

//	intron 5' boundary
		if (spj && donor_q->remain()) {
		    if (donor_q->head() < n0) donor_q->pull();
		    if (donor_q->remain())
			to_spj(*donor_q, ml + 1, n, p);
		}

//	prepare for next cycle
		if (j9 == ke && wdw.lw <= r0 && r0 <= wdw.up) {
		    hv[r0] = hv_a[p][j9];
		    fv[r0] = fv_a[j9];
		    if (dagp) fv2[r0] = fv2_a[j9];
		}
	    }	// end of n loop
	    if (ml == mc) {
const		var_t	c = vec_max(hv + wdw.lw, wdw.up - wdw.lw);
const		int	d = checkpoint(c);
		if (d < md / 2) {	// all-round down score
		    vec_add(hv + wdw.lw - 1, wdw.width, -c);
		    vec_add(fv + wdw.lw - 1, wdw.width, -c);
		    if (dagp) vec_add(fv2 + wdw.lw - 1, wdw.width, -c);
		    accscr += c;
		    mc += md;
		} else			// postpone
		    mc += d;
	    }
	}	// end of ml loop

	if (!LocalR) {
	    fhlastS1(maxh);
	    maxh.val += accscr;
	}
	return (maxh.val);
}

// assume number of stored traceback records < INT_MAX

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
VTYPE SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
forwardS1(int* pp)
{
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);
const	regist_v	two_v = Add(one_v, one_v);
const	regist_v	ninf_v = Splat(nevsel);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	gn_v = Splat(pwd->BasicGEP + pwd->BasicGOP);
const	regist_v	ge2_v = Splat(pwd->LongGEP);
const	regist_v	gn2_v = Splat(pwd->LongGEP + pwd->LongGOP);
regist_v	hd_v, fd_v, fd2_v, qd_v; // used only when used == true

	fhinitS1();

	VTYPE	accscr = 0;
const	int	md = checkpoint(0);
	int	mc = md + a->left;

	for (int ml = a->left; ml < a->right; ml += Nelem) {
const	    int	j9 = std::min(Nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + ml);
const	    int	n9 = std::min(b->right, wdw.up + (ml + j9) + 1) + j9;
	    int	n0 = n - j8;
	    int	r = n - (ml + 1);
	    vec_clear(ps_a, 2 * Nelem);
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * Nelem);
	    vec_clear(hb_a[0], 4 * Np1 + 2 * Nelem);
	    vec_clear(pb_a, Nelem);

	    if (spj) {
		for (int j = 0; j < j9; ++j)
		    fsjss[j]->reset();
		donor_q->clear();
		accep_q->clear();
	    }

regist_v    ev_v = ninf_v;
regist_v    ec_v = zero_v;
regist_v    ed_v = zero_v;
regist_v    ev2_v = ninf_v;
regist_v    ec2_v = zero_v;
regist_v    ed2_v = zero_v;

	    for (int p = 0; n < n9; ++n, ++n0, ++r, p = 1 - p) {
const		int	q = 1 - p;
const		int	r0 = r - 2 * j8;
const		int	kb = std::max(0, n - b->right);
const		int	ke = std::min(j9, n - b->left);

		if (spj && !kb) {
		    if (b->exin->isDonor(n)) donor_q->push(n);
		    if (b->exin->isAccpt(n)) accep_q->push(n);
		}

//	ordinary insertion

regist_v	qv_v = Load(hv_a[q] + 1);
regist_v	hc_v = Load(hc_a[q] + 1);
		if (used) hd_v = Load(hd_a[q] + 1);
regist_v	hv_v = Add(qv_v, gn_v);	// new insertion
		ev_v = Add(ev_v, ge_v);	// extension
regist_m	msk_m = Cmp_gt(ev_v, hv_v);
		ev_v = Blend(ev_v, hv_v, msk_m);
		ec_v = Blend(ec_v, hc_v, msk_m);
		if (used) ed_v = Blend(ed_v, hd_v, msk_m);
		if (spj) {
		    Store(ev_a, ev_v);
		    Store(ec_a, ec_v);
		    if (used) Store(ed_a, ed_v);
		}
		if (dagp) {
		    hv_v = Add(qv_v, gn2_v);	// new insertion
		    ev2_v = Add(ev2_v, ge2_v);	// extension
		    msk_m = Cmp_gt(ev2_v, hv_v);
		    ev2_v = Blend(ev2_v, hv_v, msk_m);
		    ec2_v = Blend(ec2_v, hc_v, msk_m);
		    if (used) ed2_v = Blend(ed2_v, hd_v, msk_m);
		    if (spj) {
			Store(ev2_a, ev2_v);
			Store(ec2_a, ec2_v);
			if (used) Store(ed2_a, ed2_v);
		    }
		}

//	deletion
		hv_a[q][0] = hv[r + 1];
		hc_a[q][0] = hc[r + 1];
		if (used) hd_a[q][0] = hd[r + 1];
		qv_v = Load(hv_a[q]);
		hc_v = Load(hc_a[q]);
		if (used) hd_v = Load(hd_a[q]);
		fv_a[0] = fv[r + 1];
		fc_a[0] = fc[r + 1];
		if (used) fd_a[0] = fd[r + 1];
regist_v	fv_v = Load(fv_a);
regist_v	fc_v = Load(fc_a);
		if (used) fd_v = Load(fd_a);
		fv_v = Add(fv_v, ge_v);	// gap extension
		hv_v = Add(qv_v, gn_v);	// gap open
		msk_m = Cmp_gt(fv_v, hv_v);
		fv_v = Blend(fv_v, hv_v, msk_m);
		fc_v = Blend(fc_v, hc_v, msk_m);
		if (used) fd_v = Blend(fd_v, hd_v, msk_m);
		Store(fv_a + 1, fv_v);
		Store(fc_a + 1, fc_v);
		if (used) Store(fd_a + 1, fd_v);
		if (dagp) {
		    fv2_a[0] = fv2[r + 1];
		    fc2_a[0] = fc2[r + 1];
		    if (used) fd2_a[0] = fd2[r + 1];
regist_v	    fv2_v = Load(fv2_a);
regist_v	    fc2_v = Load(fc2_a);
		    if (used) fd2_v = Load(fd2_a);
		    fv2_v = Add(fv2_v, ge2_v);	// gap extension
		    hv_v = Add(qv_v, gn2_v);	// gap open
		    msk_m = Cmp_gt(fv2_v, hv_v);
		    fv2_v = Blend(fv2_v, hv_v, msk_m);
		    fc2_v = Blend(fc2_v, hc_v, msk_m);
		    if (used) fd2_v = Blend(fd2_v, hd_v, msk_m);
		    Store(fv2_a + 1, fv2_v);
		    Store(fc2_a + 1, fc2_v);
		    if (used) Store(fd2_a + 1, fd2_v);
		    msk_m = Cmp_gt(fv2_v, fv_v);
		    fv_v = Blend(fv2_v, fv_v, msk_m);
		    fc_v = Blend(fc2_v, fc_v, msk_m);
		    if (used) fd_v = Blend(fd2_v, fd_v, msk_m);
		}

//	diagonal match

const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 1 - kb);
		Store(pv_a, zero_v);
		for (int k = kb; k < ke; ++k, ++as, --bs)
		    pv_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(pv_a);

		hv_a[p][0] = hv[r];
		hb_a[p][0] = hb[r];
		hc_a[p][0] = hc[r];
		if (used) hd_a[p][0] = hd[r];
		qv_v = Load(hv_a[p]);
regist_v	hb_v = Load(hb_a[p]);
		hc_v = Load(hc_a[p]);
		if (used) hd_v = Load(hd_a[p]);
		hv_v = Add(hv_v, qv_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
		hc_v = Blend(fc_v, hc_v, msk_m);
		if (used) hd_v = Blend(fd_v, hd_v, msk_m);
		hb_v = Blend(two_v, zero_v, msk_m);

//	best of three
		if (dagp) {
		    msk_m = Cmp_gt(ev2_v, ev_v);
regist_v	    fv2_v = Blend(ev2_v, ev_v, msk_m);
regist_v	    fc2_v = Blend(ec2_v, ec_v, msk_m);
		    if (used) qd_v = Blend(ed2_v, ed_v, msk_m);
		    msk_m = Cmp_gt(fv2_v, hv_v);	// t >= e = !(e > t)
		    hv_v = Blend(fv2_v, hv_v, msk_m);
		    hc_v = Blend(fc2_v, hc_v, msk_m);
		    if (used) hd_v = Blend(qd_v, hd_v, msk_m);
		} else {
		    msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !(e > t)
		    hv_v = Blend(ev_v, hv_v, msk_m);
		    hc_v = Blend(ec_v, hc_v, msk_m);
		    if (used) hd_v = Blend(ed_v, hd_v, msk_m);
		}
		hb_v = Blend(one_v, hb_v, msk_m);
		Store(pv_a, hb_v);	// diag: 0, hori: 1, vert: 2

		if (spj) {
regist_v	    qb_v = Load(ps_a);		// post splicing
		    qb_v = And(hb_v, qb_v);
		    Store(ps_a, qb_v);
		}

//	Find optimal path

		if (!Local) {
		    msk_m = Cmp_gt(hv_v, ninf_v);	// clamp underflow
		    hv_v = Blend(hv_v, ninf_v, msk_m);
		} else if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);	// Kadane-Gries
		    hb_v = Blend(one_v, hb_v, msk_m);
		    hc_v = Blend(zero_v, hc_v, msk_m);
		    if (used) hd_v = Blend(zero_v, hd_v, msk_m);
		}
		if (ke < j9) fc_v = hb_v;		// reserve
		msk_m = Cmp_eq(hb_v, zero_v);
		hb_v = Blend(one_v, zero_v, msk_m);	// diag: 1, others: 0
regist_v	qb_v = Load(hb_a[p]);
		qb_v = AndNot(hb_v, qb_v);		// t && !q = nd && di
		if (ke < j9) {
		    msk_m = To_mask(pb_a);
		    hb_v = Blend(hb_v, fc_v, msk_m);
		    pb_a[ke] = -1;			// all 1
		}
		Store(qb_a, qb_v);
		Store(hv_a[p] + 1, hv_v);
		Store(hb_a[p] + 1, hb_v);
		Store(hc_a[p] + 1, hc_v);
		if (used) Store(hd_a[p] + 1, hd_v);
		int	mk = ml + kb;
		int	nk = n - 1 - kb;
		for (int k = kb; k < ke; ++k, ++mk, --nk) {
const		    int	kp1 = k + 1;
		    if (qb_a[k]) {
			int ptr = s2_i(hc_a[p][kp1], hd_a[p][kp1]); 
			ptr = vmf->add(mk, nk, ptr);
			hc_a[p][kp1] = i_s2(hd_a[p] + kp1, ptr);
		    }
		}
		if (LocalR) {	// max(H)
		    var_t*	mx = vmax(hv_a[p] + 1, j9);
		    if ((*mx + accscr) > maxh.val) {
			int	k = mx - hv_a[p];
			maxh.val = *mx + accscr;
			maxh.ulk = s2_i(hc_a[p][k], hd_a[p][k]);
			maxh.mr = ml + k;
			maxh.nr = n - k + 1;
		    }
		}

//	intron 3' boundary
		if (spj && accep_q->remain()) {
		    if (accep_q->head() < n0) accep_q->pull();
		    if (accep_q->remain()) {
const			int	rv = from_spj(*accep_q, ml + 1, n, p);
			if (rv & 1) {
			    ev_v = Load(ev_a);
			    ec_v = Load(ec_a);
			    if (used) ed_v = Load(ed_a);
			}
			if (rv & 2) {
			    ev2_v = Load(ev2_a);
			    ec2_v = Load(ec2_a);
			    if (used) ed2_v = Load(ed2_a);
			}
		    }
		}

//	intron 5' boundary
		if (spj && donor_q->remain()) {
		    if (donor_q->head() < n0) donor_q->pull();
		    if (donor_q->remain())
			to_spj(*donor_q, ml + 1, n, p);
		}

//	prepare for next cycle
		if (j9 == ke && wdw.lw <= r0 && r0 <= wdw.up) {
		    hv[r0] = hv_a[p][j9];
		    hb[r0] = hb_a[p][j9];
		    hc[r0] = hc_a[p][j9];
		    if (used) hd[r0] = hd_a[p][j9];
		    fv[r0] = fv_a[j9];
		    fc[r0] = fc_a[j9];
		    if (used) fd[r0] = fd_a[j9];
		    if (dagp) {
			fv2[r0] = fv2_a[j9];
			fc2[r0] = fc2_a[j9];
			if (used) fd2[r0] = fd2_a[j9];
		    }
		}
	    }	// end of n loop

	    if (ml == mc) {
const		var_t	c = vec_max(hv + wdw.lw, wdw.up - wdw.lw);
const		int	d = checkpoint(c);
		if (d < md / 2) {	// all-round down score
		    vec_add(hv + wdw.lw - 1, wdw.width, -c);
		    vec_add(fv + wdw.lw - 1, wdw.width, -c);
		    if (dagp) vec_add(fv2 + wdw.lw - 1, wdw.width, -c);
		    accscr += c;
		    mc += md;
		} else			// postpone
		    mc += d;
	    }
	}	// end of ml loop

	if (LocalR) {
	    if (pp) *pp = vmf->add(maxh.mr, maxh.nr, maxh.ulk);
	} else {
	    fhlastS1(maxh);
	    maxh.val += accscr;
	    if (pp) *pp = maxh.ulk;
	}
	return (maxh.val);
}

// unidirectional Hirschberg method, multi-intermediate version

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
VTYPE SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
hirschbergS1(Dim10* cpos, const int& n_im)
{
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);
const	regist_v	two_v = Add(one_v, one_v);
const	regist_v	ninf_v = Splat(nevsel);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	gn_v = Splat(pwd->BasicGEP + pwd->BasicGOP);
const	regist_v	ge2_v = Splat(pwd->LongGEP);
const	regist_v	gn2_v = Splat(pwd->LongGEP + pwd->LongGOP);

	fhinitS1();

	mm = (a->right - a->left + n_im) / (n_im + 1);
	Udh_Imds	udhimds(n_im, a->left, mm, wdw, pwd->Noll);
	imd = udhimds[0];
	mm = a->left + (imd->mi - a->left - 1) / Nelem * Nelem;
	int	k9 = imd->mi - mm;
	int	k8 = k9 - 1;

const	int	md = checkpoint(0);
	int	mc = md + a->left;
	VTYPE	accscr = 0;

regist_v	hb_v, fb_v, qb_v;	// used only when LocalL == true
regist_v	hd_v, fd_v, fd2_v, qd_v;// used only when used == true

	for (int ml = a->left, i = 0; ml < a->right; ml += Nelem) {
const	    int	j9 = std::min(Nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + ml);
	    int	n9 = std::min(b->right, wdw.up + (ml + j9) + 1) + j9;
	    int	n0 = n - j8;
	    int	r = n - (ml + 1);
	    vec_clear(ps_a, 2 * Nelem);
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * Nelem);
	    vec_clear(hb_a[0], 4 * Np1 + 2 * Nelem);
const	    bool	is_imd_ = ml == mm;

	    if (spj) {
		for (int j = 0; j < j9; ++j)
		    fsjss[j]->reset();
		donor_q->clear();
		accep_q->clear();
	    }

regist_v    ev_v = ninf_v;
regist_v    eb_v = zero_v;
regist_v    ec_v = zero_v;
regist_v    ed_v = zero_v;
regist_v    ev2_v = ninf_v;
regist_v    eb2_v = zero_v;
regist_v    ec2_v = zero_v;
regist_v    ed2_v = zero_v;

	    for (int p = 0; n < n9; ++n, ++n0, ++r, p = 1 - p) {
const		int	q = 1 - p;
const		int	r0 = r - 2 * j8;
const		int	rj = r - 2 * k8;
const		int	kb = std::max(0, n - b->right);
const		int	ke = std::min(j9, n - b->left);
const		bool	is_imd = is_imd_ && rj >= wdw.lw && rj <= wdw.up;

		if (spj && !kb) {
		    if (b->exin->isDonor(n)) donor_q->push(n);
		    if (b->exin->isAccpt(n)) accep_q->push(n);
		}

//	ordinary insertion
regist_v	qv_v = Load(hv_a[q] + 1);
		if (LocalL) hb_v = Load(hb_a[q] + 1);
regist_v	hc_v = Load(hc_a[q] + 1);
		if (used) hd_v = Load(hd_a[q] + 1);
regist_v	hv_v = Add(qv_v, gn_v);	// new insertion
		ev_v = Add(ev_v, ge_v);	// extension
regist_m	msk_m = Cmp_gt(ev_v, hv_v);
		ev_v = Blend(ev_v, hv_v, msk_m);
		if (LocalL) eb_v = Blend(eb_v, hb_v, msk_m);
		ec_v = Blend(ec_v, hc_v, msk_m);
		if (used) ed_v = Blend(ed_v, hd_v, msk_m);
		if (spj) {
		    Store(ev_a, ev_v);
		    if (LocalL) Store(eb_a, eb_v);
		    Store(ec_a, ec_v);
		    if (used) Store(ed_a, ed_v);
		}
		if (dagp) {
		    hv_v = Add(qv_v, gn2_v);	// new insertion
		    ev2_v = Add(ev2_v, ge2_v);	// extension
		    msk_m = Cmp_gt(ev2_v, hv_v);
		    ev2_v = Blend(ev2_v, hv_v, msk_m);
		    if (LocalL) eb2_v = Blend(eb2_v, hb_v, msk_m);
		    ec2_v = Blend(ec2_v, hc_v, msk_m);
		    if (used) ed2_v = Blend(ed2_v, hd_v, msk_m);
		    if (spj) {
			Store(ev2_a, ev2_v);
			if (LocalL) Store(eb2_a, eb2_v);
			Store(ec2_a, ec2_v);
			if (used) Store(ed2_a, ed2_v);
		    }
		}

//	deletion
		hv_a[q][0] = hv[r + 1];
		if (LocalL) hb_a[q][0] = hb[r + 1];
		hc_a[q][0] = hc[r + 1];
		if (used) hd_a[q][0] = hd[r + 1];
		qv_v = Load(hv_a[q]);
		if (LocalL) hb_v = Load(hb_a[q]);
		hc_v = Load(hc_a[q]);
		if (used) hd_v = Load(hd_a[q]);
		fv_a[0] = fv[r + 1];
		if (LocalL) fb_a[0] = fb[r + 1];
		fc_a[0] = fc[r + 1];
		if (used) fd_a[0] = fd[r + 1];
regist_v	fv_v = Load(fv_a);
		if (LocalL) fb_v = Load(fb_a);
regist_v	fc_v = Load(fc_a);
		if (used) fd_v = Load(fd_a);
		fv_v = Add(fv_v, ge_v);	// gap extension
		hv_v = Add(qv_v, gn_v);	// gap open
		msk_m = Cmp_gt(fv_v, hv_v);
		fv_v = Blend(fv_v, hv_v, msk_m);
		if (LocalL) fb_v = Blend(fb_v, hb_v, msk_m);
		fc_v = Blend(fc_v, hc_v, msk_m);
		if (used) fd_v = Blend(fd_v, hd_v, msk_m);
		Store(fv_a + 1, fv_v);
		if (LocalL) Store(fb_a + 1, fb_v);
		Store(fc_a + 1, fc_v);
		if (used) Store(fd_a + 1, fd_v);
		if (dagp) {
		    fv2_a[0] = fv2[r + 1];
		    if (LocalL) fb2_a[0] = fb2[r + 1];
		    fc2_a[0] = fc2[r + 1];
		    if (used) fd2_a[0] = fd2[r + 1];
regist_v	    fv2_v = Load(fv2_a);
regist_v	    fb2_v = Load(fb2_a);
regist_v	    fc2_v = Load(fc2_a);
		    if (used) fd2_v = Load(fd2_a);
		    fv2_v = Add(fv2_v, ge2_v);	// gap extension
		    hv_v = Add(qv_v, gn2_v);	// gap open
		    msk_m = Cmp_gt(fv2_v, hv_v);
		    fv2_v = Blend(fv2_v, hv_v, msk_m);
		    if (LocalL) fb2_v = Blend(fb2_v, hb_v, msk_m);
		    fc2_v = Blend(fc2_v, hc_v, msk_m);
		    if (used) fd2_v = Blend(fd2_v, hd_v, msk_m);
		    Store(fv2_a + 1, fv2_v);
		    Store(fb2_a + 1, fb2_v);
		    Store(fc2_a + 1, fc2_v);
		    if (used) Store(fd2_a + 1, fd2_v);
		    msk_m = Cmp_gt(fv2_v, fv_v);
		    fv_v = Blend(fv2_v, fv_v, msk_m);
		    if (LocalL) fb_v = Blend(fb2_v, fb_v, msk_m);
		    fc_v = Blend(fc2_v, fc_v, msk_m);
		    if (used) fd_v = Blend(fd2_v, fd_v, msk_m);
		}

//	diagonal match

const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 1 - kb);
		Store(pv_a, zero_v);
		for (int k = kb; k < ke; ++k, ++as, --bs)
		    pv_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(pv_a);

		hv_a[p][0] = hv[r];
		if (LocalL) hb_a[p][0] = hb[r];
		hc_a[p][0] = hc[r];
		if (used) hd_a[p][0] = hd[r];
		qv_v = Load(hv_a[p]);
		if (LocalL) hb_v = Load(hb_a[p]);
		hc_v = Load(hc_a[p]);
		if (used) hd_v = Load(hd_a[p]);
		hv_v = Add(hv_v, qv_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
		if (LocalL) hb_v = Blend(fb_v, hb_v, msk_m);
		hc_v = Blend(fc_v, hc_v, msk_m);
		if (used) hd_v = Blend(fd_v, hd_v, msk_m);
		fb_v = Blend(two_v, zero_v, msk_m);

//	best of three
		if (dagp) {
		    msk_m = Cmp_gt(ev2_v, ev_v);
		    qv_v = Blend(ev2_v, ev_v, msk_m);
		    if (LocalL) qb_v = Blend(eb2_v, eb_v, msk_m);
regist_v	    qc_v = Blend(ec2_v, ec_v, msk_m);
		    if (used) qd_v = Blend(ed2_v, ed_v, msk_m);
		    msk_m = Cmp_gt(qv_v, hv_v);
		    hv_v = Blend(qv_v, hv_v, msk_m);
		    if (LocalL) hb_v = Blend(qb_v, hb_v, msk_m);
		    hc_v = Blend(qc_v, hc_v, msk_m);
		    if (used) hd_v = Blend(qd_v, hd_v, msk_m);
		} else {
		    msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !e > t)
		    hv_v = Blend(ev_v, hv_v, msk_m);
		    if (LocalL) hb_v = Blend(eb_v, hb_v, msk_m);
		    hc_v = Blend(ec_v, hc_v, msk_m);
		    if (used) hd_v = Blend(ed_v, hd_v, msk_m);
		}
		fb_v = Blend(one_v, fb_v, msk_m);
		Store(pv_a, fb_v);	// diag: 0, hori: 1, vert: 2

		if (spj) {
regist_v	    qb_v = Load(ps_a);		// post splicing
		    qb_v = And(fb_v, qb_v);
		    Store(ps_a, qb_v);
		}

//	Find optimal path

		if (!Local) {
		    msk_m = Cmp_gt(hv_v, ninf_v);	// clamp underflow
		    hv_v = Blend(hv_v, ninf_v, msk_m);
		} else if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);
		}
		Store(hv_a[p] + 1, hv_v);
		if (LocalL) Store(hb_a[p] + 1, hb_v);
		Store(hc_a[p] + 1, hc_v);
		if (used) Store(hd_a[p] + 1, hd_v);
		if (LocalL && !accscr) {	// left end
		    for (int k = kb; k < ke; ++k) {
const			int	kp1 = k + 1;
			if (hv_a[p][kp1] == 0) {
			    hb_a[p][kp1] = static_cast<var_t>(ml + kp1);
			    hc_a[p][kp1] = i_s2(hd_a[p] + kp1, r - 2 * k);
			}
		    }
		}
		if (LocalR) {			// right end
		    var_t*	mx = vmax(hv_a[p] + 1, j9);
		    if (*mx + accscr > maxh.val) {
			int	k = mx - hv_a[p];
			maxh.val = *mx + accscr;
			maxh.ml = hb_a[p][k];	// ml
			maxh.ulk = s2_i(hc_a[p][k], hd_a[p][k]);
			maxh.mr = ml + k;
			maxh.nr = n - k + 1;
		    }
		}

//	intron 3' boundary
		if (spj && accep_q->remain()) {
		    if (accep_q->head() < n0) accep_q->pull();
		    if (accep_q->remain()) {
const			int rv = from_spj(*accep_q, ml + 1, n, p);
			if (rv & 1) {
			    ev_v = Load(ev_a);
			    if (LocalL) eb_v = Load(eb_a);
			    ec_v = Load(ec_a);
			    if (used) ed_v = Load(ed_a);
			}
			if (rv & 2) {
			    ev2_v = Load(ev2_a);
			    if (LocalL) eb2_v = Load(eb2_a);
			    ec2_v = Load(ec2_a);
			    if (used) ed2_v = Load(ed2_a);
			}
		    }
		}

//	intron 5' boundary
		if (spj && donor_q->remain()) {
		    if (donor_q->head() < n0) donor_q->pull();
		    if (donor_q->remain())
			to_spj(*donor_q, ml + 1, n, p);
		}

//	intermediate row
		if (is_imd) {
		    if (pv_a[k8] == 0) rlst = rj;	// diag
		    if (pv_a[k8] == 1) 			// hori
			imd->hlnk[0][rj] = rlst;
		    imd->vlnk[0][rj] = s2_i(hc_a[p][k9], hd_a[p][k9]);
		    hc_a[p][k9] = i_s2(hd_a[p] + k9, rj);
		    imd->vlnk[1][rj] = s2_i(fc_a[k9], fd_a[k9]);
		    fc_a[k9] = i_s2(fd_a + k9, rj + wdw.width);
		    if (dagp)
			fc2_a[k9] = i_s2(fd2_a + k9, rj + 2 * wdw.width);
		}

//	prepare for next cycle
		if (j9 == ke && wdw.lw <= r0 && r0 <= wdw.up) {
		    hv[r0] = hv_a[p][j9];
		    if (LocalL) hb[r0] = hb_a[p][j9];
		    hc[r0] = hc_a[p][j9];
		    if (used) hd[r0] = hd_a[p][j9];
		    fv[r0] = fv_a[j9];
		    if (LocalL) fb[r0] = fb_a[j9];
		    fc[r0] = fc_a[j9];
		    if (used) fd[r0] = fd_a[j9];
		    if (dagp) {
			fv2[r0] = fv2_a[j9];
			if (LocalL) fb2[r0] = fb2_a[j9];
			fc2[r0] = fc2_a[j9];
			if (used) fd2[r0] = fd2_a[j9];
		    }
		}
	    }	// end of n loop
	    if (ml == mc) {
const		var_t	c = vec_max(hv + wdw.lw, wdw.up - wdw.lw);
const		int	d = checkpoint(c);
		if (d < md / 2) {	// down score all
		    vec_add(hv + wdw.lw - 1, wdw.width, -c);
		    vec_add(fv + wdw.lw - 1, wdw.width, -c);
		    if (dagp) vec_add(fv2 + wdw.lw - 1, wdw.width, -c);
		    accscr += c;
		    mc += md;
		} else			// postpone
		    mc += d;
	    }
	    if (is_imd_ && ++i < n_im) {	// reset intermediate
		imd = udhimds[i];
		mm = a->left + (imd->mi - a->left - 1) / Nelem * Nelem;
		k9 = imd->mi - mm;
		k8 = k9 - 1;
	    }
	}	// end of ml loop

	if (LocalR) {
	    a->right = maxh.mr;
	    b->right = maxh.nr;
	} else {	// global or semi-global
	    int	r = fhlastS1(maxh);
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
	    if (wdw.lw < imd->vlnk[d][r] && imd->vlnk[d][r] < wdw.up) {
		cpos[i][c++] = imd->mi;		// cross intermediate
		cpos[i][c++] = (d > 0)? 1: 0;
		for (int rp = imd->hlnk[d][r]; rp < end_of_ulk && r != rp;
		    rp = imd->hlnk[d][r = rp]) {
		    cpos[i][c++] = r + imd->mi;
		}
		cpos[i][c++] = r + imd->mi;
		cpos[i][c] = end_of_ulk;
		r = imd->vlnk[d][r];
		if (wdw.lw <= r || r >= wdw.up) break;
	    } else
		cpos[i][0] = end_of_ulk;	// don't cross intermediate
	}
	for (int d = 0; r > wdw.up; r -= wdw.width) ++d;
	if (LocalL) {
	    a->left = maxh.ml;
	    b->left = r + a->left;
	} else {
const	    int	rl = b->left - a->left;
	    if (b->inex.exgl && rl > r) {
		a->left = b->left - r;
		for (int j = 0; j < n_im && udhimds[j]->mi < a->left; ++j)
		    cpos[j][0] = end_of_ulk;
	    }
	    if (a->inex.exgl && rl < r) b->left = a->left + r;
	}
	if (udhimds[n_im - 1]->mi < a->left)
	    cpos[0][2] = b->left;
	return (maxh.val);
}

#include "fwd2s1_wip_simd.h"

#endif	// _FWD2S1_SIMD_H_
