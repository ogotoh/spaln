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
*	with simplified length-independent (flat) intron penalty
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

#ifndef	_FWD2S1_WIP_SIMD_H_
#define _FWD2S1_WIP_SIMD_H_

#include "rhomb_coord.h"
#include "udh_intermediate.h"

// score-only DP forward algorithm

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
VTYPE SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
scoreonlyS1_wip()
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
const	regist_v	mil_v = Splat(IntronPrm.llmt);
const	int	ipen = spj? pwd->IntPen->Penalty(): nevsel;
	regist_v	quant_v[IntronPrm.nquant];
	regist_v	mean_v[IntronPrm.nquant];
	if (spj) {
	    for (int j = 0; j < IntronPrm.nquant; ++j) {
		quant_v[j] = Splat(pwd->IntPen->qm[j].len);
		mean_v[j] = Splat(pwd->IntPen->qm[j].pen);
	    }
	}

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
const	    SGPT2*	bb = spj? b->exin->score_n(n): 0;
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * Nelem);
	    vec_clear(s5_a, 2 * Np1);
	    vec_clear(ps_a, 2 * Nelem);

regist_v    ev_v = ninf_v;
regist_v    hv2_v = ninf_v;
regist_v    fv2_v = ninf_v;
regist_v    ev2_v = ninf_v;
regist_v    hil_v = zero_v;

	    for (int p = 0; n < n9; ++n, ++n0, ++r, p = 1 - p) {
const		int	q = 1 - p;
const		int	r0 = r - 2 * j8;
const		int	kb = std::max(0, n - b->right);
const		int	ke = std::min(j9, n - b->left);

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
		    msk_m = Cmp_gt(ev2_v, ev_v);
		    ev_v = Blend(ev2_v, ev_v, msk_m);
		}

//	deletion
		hv_a[q][0] = hv[r + 1];
		qv_v = Load(hv_a[q]);
		hv_v = Add(qv_v, gn_v);	// gap open
		fv_a[0] = fv[r + 1];
regist_v	fv_v = Load(fv_a);
		fv_v = Add(fv_v, ge_v);	// gap extension
		msk_m = Cmp_gt(fv_v, hv_v);
		fv_v = Blend(fv_v, hv_v, msk_m);
		if (dagp) {
		    hv_v = Add(qv_v, gn2_v);	// long gap open
		    fv2_a[0] = fv2[r + 1];
		    fv2_v = Load(fv2_a);
		    fv2_v = Add(fv2_v, ge2_v);	// long gap extension
		    msk_m = Cmp_gt(fv2_v, hv_v);
		    fv2_v = Blend(fv2_v, hv_v, msk_m);
		    Store(fv2_a + 1, fv2_v);
		    msk_m = Cmp_gt(fv2_v, fv_v);
		    fv_v = Blend(fv2_v, fv_v, msk_m);
		}
		Store(fv_a + 1, fv_v);

//	diagonal match

const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 1 - kb);
		if (kb) Store(pv_a, zero_v);
		for (int k = kb; k < ke; ++k, ++as, --bs)
		    pv_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(pv_a);

		hv_a[p][0] = hv[r];
		qv_v = Load(hv_a[p]);
		hv_v = Add(hv_v, qv_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
regist_v	pb_v = Blend(two_v, zero_v, msk_m);

//	best of three
		msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !(e > t)
		hv_v = Blend(ev_v, hv_v, msk_m);
		pb_v = Blend(one_v, pb_v, msk_m);	// 0: d, 1: h, 2: v

//	intron 3' boundary
		if (spj) {
		    s3_a[0] = kb? 0: bb->sig3;
regist_v	    ss_v = Load(s3_a);		// accepter signal
		    Store(s3_a + 1, ss_v);
		    qv_v = Add(hv2_v, ss_v);
regist_v	    pv_v = mean_v[0];
		    for (int j = 1; j < IntronPrm.nquant; ++j) {
			msk_m = Cmp_gt(hil_v, quant_v[j - 1]);
			pv_v = Blend(mean_v[j], pv_v, msk_m);
		    }
		    qv_v = Add(qv_v, pv_v);
		    msk_m = Cmp_gt(hil_v, mil_v);	// filter shorter than
		    qv_v = Blend(qv_v, ninf_v, msk_m);	// lower limit of intron
		    msk_m = Cmp_gt(qv_v, hv_v);
		    hv_v = Blend(qv_v, hv_v, msk_m);
		    qv_v = Blend(one_v, zero_v, msk_m);
		}

//	Find optimal path

#if _SIMD_PP_
		if (!Local) {
		    msk_m = Cmp_gt(hv_v, ninf_v);	// clamp underflow
		    hv_v = Blend(hv_v, ninf_v, msk_m);
		    if (spj || dagp) {
			msk_m = Cmp_gt(fv2_v, ninf_v);	// clamp underflow
			fv2_v = Blend(fv2_v, ninf_v, msk_m);
			msk_m = Cmp_gt(ev2_v, ninf_v);	// clamp underflow
			ev2_v = Blend(ev2_v, ninf_v, msk_m);
		    }
		}
#endif
		if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);	// Kadane-Gries
		}

		Store(hv_a[p] + 1, hv_v);
		Store(fv_a + 1, fv_v);
		if (dagp) Store(fv2_a + 1, fv2_v);
		if (LocalR) {	// max(H)
const		    var_t*	mx = vmax(hv_a[p] + 1, j9);
		    if ((*mx + accscr) > maxh.val)
			maxh.val = *mx + accscr;
		}

//	intron 5' boundary
		if (spj) {
		    s5_a[0] = kb? 0: bb->sig5 + ipen;
regist_v	    ss_v = Load(s5_a);		// accepter signal
		    Store(s5_a + 1, ss_v);
		    qv_v = Add(hv_v, ss_v);
		    msk_m = Cmp_gt(qv_v, hv2_v);
		    hv2_v = Blend(qv_v, hv2_v, msk_m);
		    hil_v = Blend(zero_v, hil_v, msk_m);
		    hil_v = Add(hil_v, one_v);
		}

//	prepare for next cycle
		if (j9 == ke && wdw.lw <= r0 && r0 <= wdw.up) {
		    hv[r0] = hv_a[p][j9];
		    fv[r0] = fv_a[j9];
		    if (dagp) fv2[r0] = fv2_a[j9];
		}
		if (bb) ++bb;
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

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
VTYPE SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
forwardS1_wip(Mfile* mfd)
{
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);		// DIAG
const	regist_v	two_v = Add(one_v, one_v);	// HORI
const	regist_v	c3_v = Add(one_v, two_v);	// HORL
const	regist_v	c4_v = Add(two_v, two_v);	// HORL
const	regist_v	c8_v = Add(c4_v, c4_v);		// VERT
const	regist_v	c9_v = Add(one_v, c8_v);	// VERL
const	regist_v	acc_v = Splat(static_cast<var_t>(TraceBackCode::ACCR));		// ACC
const	regist_v	c16_v = Add(c8_v, c8_v);	// NHOR
const	regist_v	c32_v = Add(c16_v, c16_v);	// NVER
const	regist_v	c64_v = Add(c32_v, c32_v);	// NHOL
const	regist_v	c128_v = Add(c64_v, c64_v);	// NVEL
const	regist_v	ninf_v = Splat(nevsel);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	gn_v = Splat(pwd->BasicGEP + pwd->BasicGOP);
const	regist_v	ge2_v = Splat(pwd->LongGEP);
const	regist_v	gn2_v = Splat(pwd->LongGEP + pwd->LongGOP);
const	regist_v	mil_v = Splat(IntronPrm.llmt);
const	int	ipen = spj? pwd->IntPen->Penalty(): nevsel;
	regist_v	quant_v[IntronPrm.nquant];
	regist_v	mean_v[IntronPrm.nquant];
	if (spj) {
	    for (int j = 0; j < IntronPrm.nquant; ++j) {
		quant_v[j] = Splat(pwd->IntPen->qm[j].len);
		mean_v[j] = Splat(pwd->IntPen->qm[j].pen);
	    }
	}

	fhinitS1();
	Anti_rhomb_coord<CHAR>	trb(a->right, b->right, a->left, b->left, 1);
	if (!a->inex.exgl) trb.initialize_m0(2);	// HORI

	VTYPE	accscr = 0;
const	int	md = checkpoint(0);
	int	mc = md + a->left;
const	int	mw = a->right - a->left;
const	int	mb = a->right - 2 * Nelem;
const	int	mt = a->left + mw / Nelem * Nelem;
	vclear(pv_a, Nelem);
	for (int k = 0; k < a->right - mt; ++k) pv_a[k] = 1;
regist_v    qv_v = Load(pv_a);	
const	regist_m	lmask = Cmp_eq(qv_v, one_v);

	for (int ml = a->left; ml < a->right; ml += Nelem) {
const	    int	j9 = std::min(Nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + ml);
	    int	n9 = std::min(b->right, wdw.up + (ml + j9) + 1) + j9;
	    int	n0 = n - j8;
	    int	r = n - (ml + 1);
const	    SGPT2*	bb = spj? b->exin->score_n(n): 0;
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * Nelem);
	    vec_clear(s5_a, 2 * Np1);
	    vec_clear(ps_a, 2 * Nelem);

regist_v    ev_v = ninf_v;
regist_v    hv2_v = ninf_v;
regist_v    fv2_v = ninf_v;
regist_v    ev2_v = fv2_v;
regist_v    hil_v = zero_v;	// intron length

	    for (int p = 0; n <= n9; ++n, ++n0, ++r, p = 1 - p) {
const		int	q = 1 - p;
const		int	r0 = r - 2 * j8;
const		int	kb = std::max(0, n - b->right);
const		int	ke = std::min(j9, n - b->left);
		CHAR*	trb_buf = trb.set_point(ml + 1, n);

//	horizontal ordinary insertion
regist_v	qv_v = Load(hv_a[q] + 1);
regist_v	hv_v = Add(qv_v, gn_v);		// new insertion
		ev_v = Add(ev_v, ge_v);		// extension
regist_m	msk_m = Cmp_gt(ev_v, hv_v);
		ev_v = Blend(ev_v, hv_v, msk_m);
regist_v	hb_v = Blend(zero_v, c16_v, msk_m);	// NHOR
		if (dagp) {
		    hv_v = Add(qv_v, gn2_v);	// new long insertion
		    ev2_v = Add(ev2_v, ge2_v);	// long extension
		    msk_m = Cmp_gt(ev2_v, hv_v);
		    ev2_v = Blend(ev2_v, hv_v, msk_m);
		    qv_v = Blend(zero_v, c64_v, msk_m);	// NHOL
		    hb_v = Or(hb_v, qv_v);
		}

//	vertical deletion
		hv_a[q][0] = hv[r + 1];
		qv_v = Load(hv_a[q]);
		fv_a[0] = fv[r + 1];
regist_v	fv_v = Load(fv_a);
		fv_v = Add(fv_v, ge_v);	// gap extension
		hv_v = Add(qv_v, gn_v);	// gap open
		msk_m = Cmp_gt(fv_v, hv_v);
		fv_v = Blend(fv_v, hv_v, msk_m);
		Store(fv_a + 1, fv_v);
		qv_v = Blend(zero_v, c32_v, msk_m);	// NVER
		hb_v = Or(hb_v, qv_v);
		if (dagp) {
		    fv2_a[0] = fv2[r + 1];
		    fv2_v = Load(fv2_a);
		    fv2_v = Add(fv2_v, ge2_v);	// long gap extension
		    hv_v = Add(qv_v, gn2_v);	// long gap open
		    msk_m = Cmp_gt(fv2_v, hv_v);
		    fv2_v = Blend(fv2_v, hv_v, msk_m);
		    Store(fv2_a + 1, fv2_v);
		    qv_v = Blend(zero_v, c128_v, msk_m);// NVEL
		    hb_v = Or(hb_v, qv_v);
		}

//	diagonal match

const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 1 - kb);
		if (kb) Store(pv_a, zero_v);
		for (int k = kb; k < ke; ++k, ++as, --bs)
		    pv_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(pv_a);

		hv_a[p][0] = hv[r];
		qv_v = Load(hv_a[p]);
		hv_v = Add(hv_v, qv_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
regist_v	pb_v = Blend(c8_v, one_v, msk_m);	// VERT
		if (dagp) {
		    msk_m = Cmp_gt(fv2_v, hv_v);
		    hv_v = Blend(fv2_v, hv_v, msk_m);
		    pb_v = Blend(c9_v, pb_v, msk_m);	// VERL
		}

//	best of three
		msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !e > t)
		hv_v = Blend(ev_v, hv_v, msk_m);
		pb_v = Blend(two_v, pb_v, msk_m);	// HORI
			// diag: 1, hori: 2, vert: 8
		if (dagp) {
		    msk_m = Cmp_gt(ev2_v, hv_v);
		    hv_v = Blend(ev2_v, hv_v, msk_m);
		    pb_v = Blend(c3_v, pb_v, msk_m);	// HORL
		}

//	intron 3' boundary

regist_v	ab_v = zero_v;
		if (spj) {
		    s3_a[0] = kb? 0: bb->sig3;
regist_v	    ss_v = Load(s3_a);		// accepter signal
		    Store(s3_a + 1, ss_v);
		    qv_v = Add(hv2_v, ss_v);
regist_v	    pv_v = mean_v[0];
		    for (int j = 1; j < IntronPrm.nquant; ++j) {
			msk_m = Cmp_gt(hil_v, quant_v[j - 1]);
			pv_v = Blend(mean_v[j], pv_v, msk_m);
		    }
		    qv_v = Add(qv_v, pv_v);
		    msk_m = Cmp_gt(hil_v, mil_v);
		    qv_v = Blend(qv_v, ninf_v, msk_m);	// filter short intron
		    msk_m = Cmp_gt(qv_v, hv_v);		// is acceptor ?
		    hv_v = Blend(qv_v, hv_v, msk_m);
		    pb_v = Blend(acc_v, pb_v, msk_m);
		    ab_v = Blend(two_v, zero_v, msk_m);
		}

//	Find optimal path

#if _SIMD_PP_
		if (!Local) {
		    msk_m = Cmp_gt(hv_v, ninf_v);	// clamp underflow
		    hv_v = Blend(hv_v, ninf_v, msk_m);
		    if (spj || dagp) {
			msk_m = Cmp_gt(fv2_v, ninf_v);	// clamp underflow
			fv2_v = Blend(fv2_v, ninf_v, msk_m);
			msk_m = Cmp_gt(ev2_v, ninf_v);	// clamp underflow
			ev2_v = Blend(ev2_v, ninf_v, msk_m);
		    }
		}
#endif
		if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);
		    hb_v = Blend(zero_v, hb_v, msk_m);
		}
		Store(hv_a[p] + 1, hv_v);
		if (LocalR) {	// max(H)
		    var_t*	mx = vmax(hv_a[p] + 1, j9);
		    if ((*mx + accscr) > maxh.val) {
			int	k = mx - hv_a[p];
			maxh.val = *mx + accscr;
			maxh.mr = ml + k;
			maxh.nr = n - k + 1;
		    }
		}

//	intron 5' boundary
		if (spj) {
		    s5_a[0] = kb? 0: bb->sig5 + ipen;
regist_v	    ss_v = Load(s5_a);			// donor signal
		    Store(s5_a + 1, ss_v);
		    qv_v = Add(hv_v, ss_v);
		    msk_m = Cmp_eq(ab_v, zero_v);
		    qv_v = Blend(qv_v, ninf_v, msk_m);	// prevent empty intron
		    msk_m = Cmp_gt(qv_v, hv2_v);	// is donor
		    hv2_v = Blend(qv_v, hv2_v, msk_m);
		    qv_v = Blend(c128_v, zero_v, msk_m);
		    hb_v = Or(hb_v, qv_v);	// set donor mark
		    hil_v = Blend(zero_v, hil_v, msk_m);
		    hil_v = Add(hil_v, one_v);
		} // end of spj

//	prepare for next cycle
		if (j9 == ke && wdw.lw <= r0 && r0 <= wdw.up) {
		    hv[r0] = hv_a[p][j9];
		    fv[r0] = fv_a[j9];
		    if (dagp) fv2[r0] = fv2_a[j9];
		}
		hb_v = Or(hb_v, pb_v);
		if (ml == mt) hb_v = Blend(hb_v, zero_v, lmask);
		hb_v = Cast16to8(hb_v);
		if (ml > mb) {
		    qv_v = Load((var_t*) trb_buf);
		    hb_v = Or(hb_v, qv_v);
		}
		Store((var_t*) trb_buf, hb_v);
		if (bb) ++bb;
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
	}	// end of ml loop

	if (!LocalR) {
	    fhlastS1(maxh);
	    maxh.val += accscr;
	}
	trb.traceback(maxh.mr, maxh.nr, mfd);
	return (maxh.val);
}

template <typename var_t, int Nelem, typename regist_v, typename regist_m>
VTYPE SimdAln2s1<var_t, Nelem, regist_v, regist_m>::
hirschbergS1_wip(Dim10* cpos, const int& n_im)
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
const	regist_v	mil_v = Splat(IntronPrm.llmt);
const	int	ipen = spj? pwd->IntPen->Penalty(): nevsel;
	regist_v	quant_v[IntronPrm.nquant];
	regist_v	mean_v[IntronPrm.nquant];
	if (spj) {
	    for (int j = 0; j < IntronPrm.nquant; ++j) {
		quant_v[j] = Splat(pwd->IntPen->qm[j].len);
		mean_v[j] = Splat(pwd->IntPen->qm[j].pen);
	    }
	}

	fhinitS1();

	mm = (a->right - a->left + n_im) / (n_im + 1);
	Udh_Imds	udhimds(n_im, a->left, mm, wdw, pwd->Noll);
	imd = udhimds[0];
	mm = a->left + (imd->mi - a->left - 1) / Nelem * Nelem;
	int	k9 = imd->mi - mm;
	int	k8 = k9 - 1;

const	int	md = checkpoint(0);
	int	mc = a->left + md;
	VTYPE	accscr = 0;

	for (int ml = a->left, i = 0; ml < a->right; ml += Nelem) {
const	    int	j9 = std::min(Nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + ml);
	    int	n9 = std::min(b->right, wdw.up + (ml + j9) + 1) + j9;
	    int	n0 = n - j8;
	    int	r = n - (ml + 1);
	    int	donor_r = r;
const	    SGPT2*	bb = spj? b->exin->score_n(n): 0;
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * Nelem);
	    vec_clear(s5_a, 2 * Np1);
	    vec_clear(ps_a, 2 * Nelem);
const	    bool	is_imd_ = ml == mm;

regist_v    ev_v = ninf_v;
regist_v    eb_v = zero_v;
regist_v    ec_v = zero_v;
regist_v    ed_v = zero_v;
regist_v    hv2_v = ninf_v;
regist_v    hb2_v = zero_v;
regist_v    hc2_v = zero_v;
regist_v    hd2_v = zero_v;
regist_v    fv2_v = ninf_v;
regist_v    fb2_v = zero_v;
regist_v    fc2_v = zero_v;
regist_v    fd2_v = zero_v;
regist_v    ev2_v = fv2_v;
regist_v    eb2_v = zero_v;
regist_v    ec2_v = zero_v;
regist_v    ed2_v = zero_v;
regist_v    hil_v = zero_v;

	    for (int p = 0; n < n9; ++n, ++n0, ++r, p = 1 - p) {
const		int	q = 1 - p;
const		int	r0 = r - 2 * j8;
const		int	rj = r - 2 * k8;
const		int	kb = std::max(0, n - b->right);
const		int	ke = std::min(j9, n - b->left);
const		bool	is_imd = (is_imd_ && rj >= wdw.lw && rj <= wdw.up);

//	ordinary insertion
regist_v	hb_v, fb_v, qb_v;
regist_v	hd_v, fd_v, qd_v;
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
		if (dagp) {
		    hv_v = Add(qv_v, gn2_v);	// new insertion
		    ev2_v = Add(ev2_v, ge2_v);	// extension
		    msk_m = Cmp_gt(ev2_v, hv_v);
		    ev2_v = Blend(ev2_v, hv_v, msk_m);
		    if (LocalL) eb2_v = Blend(eb2_v, hb_v, msk_m);
		    ec2_v = Blend(ec2_v, hc_v, msk_m);
		    if (used) ed2_v = Blend(ed2_v, hd_v, msk_m);
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
		    fv2_v = Load(fv2_a);
		    if (LocalL) fb2_v = Load(fb2_a);
		    fc2_v = Load(fc2_a);
		    if (used) fd2_v = Load(fd2_a);
		    fv2_v = Add(fv2_v, ge2_v);	// gap extension
		    hv_v = Add(qv_v, gn2_v);	// gap open
		    msk_m = Cmp_gt(fv2_v, hv_v);
		    fv2_v = Blend(fv2_v, hv_v, msk_m);
		    if (LocalL) fb2_v = Blend(fb2_v, hb_v, msk_m);
		    fc2_v = Blend(fc2_v, hc_v, msk_m);
		    if (used) fd2_v = Blend(fd2_v, hd_v, msk_m);
		    Store(fv2_a + 1, fv2_v);
		    if (LocalL) Store(fb2_a + 1, fb2_v);
		    Store(fc2_a + 1, fc2_v);
		    if (used) Store(fd2_a + 1, fd2_v);
		}

//	diagonal match

const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 1 - kb);
		if (kb) Store(pv_a, zero_v);
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
regist_v	pb_v = Blend(two_v, zero_v, msk_m);

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
		pb_v = Blend(one_v, pb_v, msk_m);
		if (is_imd)    Store(pv_a, pb_v);	// d: 0, h: 1, v: 2

//	intron 3' boundary

		if (spj) {
		    s3_a[0] = kb? 0: bb->sig3;
regist_v	    ss_v = Load(s3_a);			// accepter signal
		    Store(s3_a + 1, ss_v);
		    qv_v = Add(hv2_v, ss_v);
regist_v	    pv_v = mean_v[0];
		    for (int j = 1; j < IntronPrm.nquant; ++j) {
			msk_m = Cmp_gt(hil_v, quant_v[j - 1]);
			pv_v = Blend(mean_v[j], pv_v, msk_m);
		    }
		    qv_v = Add(qv_v, pv_v);
		    msk_m = Cmp_gt(hil_v, mil_v);
		    qv_v = Blend(qv_v, ninf_v, msk_m);	// filter short intron
		    msk_m = Cmp_gt(qv_v, hv_v);		// is acceptor
		    hv_v = Blend(qv_v, hv_v, msk_m);
		    if (LocalL) hb_v = Blend(hb2_v, hb_v, msk_m);
		    hc_v = Blend(hc2_v, hc_v, msk_m);
		    if (used) hd_v = Blend(hd2_v, hd_v, msk_m);
		    qv_v = Blend(one_v, zero_v, msk_m);
		    if (is_imd) {
			Store(ps_a, qv_v);
		        if (ps_a[k8]) {
			    imd->hlnk[0][rj] = donor_r;
			    imd->hlnk[1][rj] = (donor_r + wdw.width);
			    rlst = rj;
			}
		    }
		}

//	Find optimal path

#if _SIMD_PP_
		if (!Local) {
		    msk_m = Cmp_gt(hv_v, ninf_v);	// clamp underflow
		    hv_v = Blend(hv_v, ninf_v, msk_m);
		    if (spj || dagp) {
			msk_m = Cmp_gt(fv2_v, ninf_v);	// clamp underflow
			fv2_v = Blend(fv2_v, ninf_v, msk_m);
			msk_m = Cmp_gt(ev2_v, ninf_v);	// clamp underflow
			ev2_v = Blend(ev2_v, ninf_v, msk_m);
		    }
		}
#endif
		if (LocalL && !accscr) {
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

//	intron 5' boundary
		if (spj) {
		    s5_a[0] = kb? 0: bb->sig5 + ipen;
regist_v	    ss_v = Load(s5_a);		// accepter signal
		    Store(s5_a + 1, ss_v);
		    qv_v = Add(hv_v, ss_v);
		    msk_m = Cmp_gt(qv_v, hv2_v);
		    hv2_v = Blend(qv_v, hv2_v, msk_m);
		    if (LocalL) hb2_v = Blend(hb_v, hb2_v, msk_m);
		    hc2_v = Blend(hc_v, hc2_v, msk_m);
		    if (used) hd2_v = Blend(hd_v, hd2_v, msk_m);
		    hil_v = Blend(zero_v, hil_v, msk_m);
		    hil_v = Add(hil_v, one_v);
		    if (is_imd) {
			qv_v = Blend(one_v, zero_v, msk_m);
			Store(pb_a, qv_v);
			if (pb_a[k8]) {		// spliced
			    donor_r = rj;
const			    int	rl = s2_i(hc_a[p][k9], hd_a[p][k9]);
const			    int&	k = pv_a[k8];
			    if (!(k & 1) && rj != rl)
				imd->vlnk[k / 2][rj] = rl;
			}
		    }
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
		    if (dagp) {
			fc2_a[k9] = i_s2(fd2_a + k9, rj + 2 * wdw.width);
			imd->vlnk[2][rj] = s2_i(fc2_a[k9], fd2_a[k9]);
		    }
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
		if (bb) ++bb;
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
	    for ( ; r >= wdw.up; r -= wdw.width) ++d;
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
		cpos[i][0] =  end_of_ulk;	// don't cross intermediate
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

#endif	// _FWD2S1_WIP_SIMD_H_
