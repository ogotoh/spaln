/*****************************************************************************
*
*	Simplified version of vetorized DP-based spliced alignment of 
*	protein vs DNA 
*	Accelerated by SIMD instructions of Intel intrinsics.
*	Curently accepted architectures: SSE4.1, AVX, AVX2, and AVX-512BW
*	Assumes capability of store/load in/from unaligned memory
*
*	5' and 3' splice site signals, coding potential, and conservation
*	of intron sites are considered.
*	Assumes there is no internal gap in the reference cDNA/EST sequence.
*
*	Implement unidirectional Hirschberg linear-space algorithm
*	with simplified well-shaped intron penalty function
*	togegher with bit-map based forward-traceback algorithm
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
*	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
*
*****************************************************************************/

#ifndef	_FWD2H1_WIP_SIMD_H_
#define _FWD2H1_WIP_SIMD_H_

static	const	int	MaxNquant = 5;

static	const	short	donor_code[4] = {
	static_cast<short>(TraceBackCode::DONM), 
	static_cast<short>(TraceBackCode::DONZ), 
	static_cast<short>(TraceBackCode::DONP), 0};
static	const	short	accpr_code[4] = {
	static_cast<short>(TraceBackCode::ACCM), 
	static_cast<short>(TraceBackCode::ACCZ), 
	static_cast<short>(TraceBackCode::ACCP), 0};
static	const	int	min_ssv = -1000;

VTYPE SimdAln2h1::forwardH1_wip(Mfile* mfd)
{
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);		// DIAG
const	regist_v	two_v = Add(one_v, one_v);	// HORI
const	regist_v	c4_v = Add(two_v, two_v);	// HOR1
const	regist_v	c5_v = Add(c4_v, one_v);	// HOR2
const	regist_v	c8_v = Add(c4_v, c4_v);		// VERT
const	regist_v	c10_v = Add(c8_v, two_v);	// VER1
const	regist_v	c11_v = Add(c10_v, one_v);	// VER2
const	regist_v	c16_v = Add(c8_v, c8_v);	// NHOR
const	regist_v	c32_v = Add(c16_v, c16_v);	// NVER
const	regist_v	donor_v[3] =
	{Splat(donor_code[0]), Splat(donor_code[1]), Splat(donor_code[2])};
const	regist_v	accpr_v[3] =
	{Splat(accpr_code[0]), Splat(accpr_code[1]), Splat(accpr_code[2])};
const	regist_v	ninf_v = Splat(nevsel);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	g1_v = Splat(pwd->GapW1);
const	regist_v	g2_v = Splat(pwd->GapW2);
const	regist_v	g3_v = Splat(pwd->GapW3);
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

	Anti_rhomb_coord<TBU_t>	trb(a->right, b->right, a->left, b->left, 3);
	if (!a->inex.exgl) trb.initialize_m0(4);	// VERT
	fhinitH1(&trb);

	VTYPE	accscr = 0;
const	int	md = checkpoint(0);
	int	mc = md + a->left;
const	int	mw = a->right - a->left;
const	int	mb = a->right - Nelem;
const	int	mt = a->left + mw / Nelem * Nelem;
	vclear(sm_a, Nelem);
	for (int k = 0; k < a->right - mt; ++k) sm_a[k] = 1;
regist_v    qv_v = Load(sm_a);
const	regist_m	lmask = Cmp_eq(qv_v, one_v);

	for (int ml = a->left; ml < a->right; ml += Nelem) {
const	    int	j9 = std::min(nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + 3 * ml);
	    int	n9 = std::min(b->right, wdw.up + 3 * (ml + j9) + 1) + 3 * j9;
	    int	n0 = n - 3 * j8;
const	    int	mp1 = ml + 1;
	    int	q = (n + 3 * mp1) % 6;
	    int	r = n - 3 * mp1;
	    vec_set(hv_a[0], nevsel, 12 * Np1 + 3 * Nelem);	// [hv_a..fv_a]
	    vec_clear(sm_a, 28 * Np1);		// sm_a, cp_a, s5/3_a, s5/3_a

regist_v    hiv_v[3] = {ninf_v, ninf_v, ninf_v};
regist_v    hil_v[3] = {zero_v, zero_v, zero_v};	// intron length

	    for ( ; n <= n9; ++n, ++n0, ++r, q = modN<6>(q + 1)) {
const		int	p = q % 3;
const		int	nb = std::max(0, n - b->right + 1);
const		int	kb = (nb - 1) / 3;
const		int	ke = std::min(j9, (n - b->left) / 3);
const		SGPT6*	bb = b->exin->score_p(n);
		TBU_t*	trb_buf = trb.set_point(mp1, n);

		cp_a[p][0] = (b->exin->good(bb - 2))? bb[-2].sigE: 0;
regist_v	cv_v = Load(cp_a[p]);	// coding potential
		Store(cp_a[p] + 1, cv_v);

//	horizontal ordinary insertion
		int	qq = modN<6>(q - 1);
regist_v	hv_v = Load(hv_a[qq] + 1);
		hv_v = Add(hv_v, g1_v);		// frame shift
		qq = modN<6>(q - 2);
		qv_v = Load(hv_a[qq] + 1);
		qv_v = Add(qv_v, g2_v);		// frame shift
regist_m	msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
regist_v        eb_v = Blend(c4_v, c5_v, msk_m);	// 1 or 2-nt ins
		qq = modN<6>(q - 3);
		qv_v = Load(hv_a[qq] + 1);
		qv_v = Add(qv_v, g3_v);		// new insertion
		qv_v = Add(qv_v, cv_v);
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		eb_v = Blend(eb_v, two_v, msk_m);
regist_v	ev_v = Load(ev_a[p]);
		ev_v = Add(ev_v, ge_v);		// extension
		ev_v = Add(ev_v, cv_v);
		msk_m = Cmp_gt(ev_v, hv_v);
		ev_v = Blend(ev_v, hv_v, msk_m);
regist_v	hb_v = Blend(zero_v, c16_v, msk_m);	// NHOR
		eb_v = Blend(two_v, eb_v, msk_m);
		Store(ev_a[p], ev_v);

//	vertical deletion
		fv_a[qq][0] = fv[r + 3];
regist_v	fv_v = Load(fv_a[qq]);
		fv_v = Add(fv_v, ge_v);		// gap extension

		hv_a[qq][0] = hv[r + 3];
		hv_v = Load(hv_a[qq]);
		hv_v = Add(hv_v, g3_v);		// gap open

		qq = modN<6>(q - 4);		// frame shift
		hv_a[qq][0] = hv[r + 2];
		qv_v = Load(hv_a[qq]);
		qv_v = Add(qv_v, g2_v);
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
regist_v	pb_v = Blend(c8_v, c10_v, msk_m);

		qq = modN<6>(q - 5);		// frame shift
		hv_a[qq][0] = hv[r + 1];
		qv_v = Load(hv_a[qq]);
		qv_v = Add(qv_v, g1_v);
		msk_m = Cmp_gt(hv_v, qv_v);
		hv_v = Blend(hv_v, qv_v, msk_m);
		pb_v = Blend(pb_v, c11_v, msk_m);

		msk_m = Cmp_gt(fv_v, hv_v);
		fv_v = Blend(fv_v, hv_v, msk_m);
		qv_v = Blend(zero_v, c32_v, msk_m);	// NVER
		hb_v = Or(hb_v, qv_v);
		pb_v = Blend(c8_v, pb_v, msk_m);
		Store(fv_a[q] + 1, fv_v);

//	diagonal match
const		CHAR*	as = a->at(ml + kb);
const		CHAR*	bs = b->at(n - 3 * kb - 2);
		if (nb) Store(sm_a, zero_v);
		for (int k = kb; k < ke; ++k, ++as, bs -= 3)
		    sm_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(sm_a);

		hv_a[q][0] = hv[r];
regist_v	dv_v = Load(hv_a[q]);
		hv_v = Add(hv_v, dv_v);
		hv_v = Add(hv_v, cv_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
		pb_v = Blend(pb_v, one_v, msk_m);	// VERT

//	best of three
		msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !e > t)
		hv_v = Blend(ev_v, hv_v, msk_m);
		pb_v = Blend(eb_v, pb_v, msk_m);	// HORI
			// diag: 1, hori: 2, vert: 8

//	intron 3' boundary

regist_v	ab_v = zero_v;	// is acceptor
		if (spj) {
		  for (int k = 0; k < 2; ++k) {
const		    bool	leg = !nb && bb->phs3 > -2 
			&& (!k || bb->phs3 == 2);
const		    int phase = leg? (bb->phs3 == 2? 
			(k? 1: -1): (k? 2: bb->phs3)): 2;
const		    int	pk = 2 * p + k;
		    s3_a[pk][0] = phase < 2? bb[-phase].sig3: min_ssv;
		    p3_a[pk][0] = accpr_code[phase + 1];
regist_v	    ss_v = Load(s3_a[pk]);		// donor signal
regist_v	    ph_v = Load(p3_a[pk]);
		    Store(s3_a[pk] + 1, ss_v);
		    Store(p3_a[pk] + 1, ph_v);
		    if (AllZero(ph_v)) continue;

		    for (int f = k? 2: 0; f < 3; ++f) {
			qv_v = Add(hiv_v[f], ss_v);
regist_v		pv_v = mean_v[0];
			for (int j = 1; j < IntronPrm.nquant; ++j) {
			    msk_m = Cmp_gt(hil_v[f], quant_v[j - 1]);
			    pv_v = Blend(mean_v[j], pv_v, msk_m);
			}
			qv_v = Add(qv_v, pv_v);
			msk_m = Cmp_eq(ph_v, accpr_v[f]);
			qv_v = Blend(qv_v, ninf_v, msk_m);
			msk_m = Cmp_gt(hil_v[f], mil_v);	// filter 
			qv_v = Blend(qv_v, ninf_v, msk_m);	// short intron
			msk_m = Cmp_gt(qv_v, hv_v);		// is acceptor ?
			hv_v = Blend(qv_v, hv_v, msk_m);
			pb_v = Blend(accpr_v[f], pb_v, msk_m);
			qv_v = Blend(ph_v, zero_v, msk_m);
			ab_v = Or(ab_v, qv_v);			// is acceptor
		    }
		  }
		}

//	Find optimal path

		if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);
		    hb_v = Blend(zero_v, hb_v, msk_m);
		}
		Store(hv_a[q] + 1, hv_v);

		if (LocalR) {	// max(H)
		    var_t*	mx = vmax(hv_a[q] + 1, j9);
		    if ((*mx + accscr) > maxh.val) {
			int	k = mx - hv_a[q];
			maxh.val = *mx + accscr;
			maxh.mr = ml + k;
			maxh.nr = n - 3 * k + 3;
		    }
		}

//	intron 5' boundary

		if (spj) {
const	regist_m  nem_m = Cmp_eq(ab_v, zero_v);		// non-empty exon
		  for (int k = 0; k < 2; ++k) {
const		    bool	leg = !nb && bb->phs5 > -2
			&& (!k || bb->phs5 == 2);
const		    int phase = leg? (bb->phs5 == 2? 
			(k? 1: -1): (k? 2: bb->phs5)): 2;
const		    int	pk = 2 * p + k;
		    s5_a[pk][0] = phase < 2? (bb[-phase].sig5 + ipen): min_ssv;
		    p5_a[pk][0] = donor_code[phase + 1];
regist_v	    ss_v = Load(s5_a[pk]);		// donor signal
regist_v	    ph_v = Load(p5_a[pk]);
		    Store(s5_a[pk] + 1, ss_v);
		    Store(p5_a[pk] + 1, ph_v);
		    if (AllZero(ph_v)) continue;
		    qv_v = Add(hv_v, ss_v);

		    for (int f = k? 2: 0; f < 3; ++f) {
regist_v	        pv_v = (f == 2)? Add(dv_v, ss_v): qv_v;
			pv_v = Blend(pv_v, ninf_v, nem_m);	// non-empty exon
			msk_m = Cmp_eq(ph_v, donor_v[f]);
			pv_v = Blend(pv_v, ninf_v, msk_m);	// select phase
			msk_m = Cmp_gt(pv_v, hiv_v[f]);
			hiv_v[f] = Blend(pv_v, hiv_v[f], msk_m);
			pv_v = Blend(donor_v[f], zero_v, msk_m);// is donor
			hb_v = Or(hb_v, pv_v);	// set intron mark
			hil_v[f] = Blend(zero_v, hil_v[f], msk_m);
		    }
		  }
		  for (int f = 0; f < 3; ++f)
			hil_v[f] = Add(hil_v[f], one_v);
		} // end of spj

//	prepare for next cycle
		int	r0 = r - 6 * j8;
		if (j9 == ke && wdw.lw <= r0 && r0 <= wdw.up) {
		    hv[r0] = hv_a[q][j9];
		    fv[r0] = fv_a[q][j9];
		}
		hb_v = Or(hb_v, pb_v);
		if (ml == mt) hb_v = Blend(hb_v, zero_v, lmask);
		if (ml > mb) {
		    qv_v = Load((var_t*) trb_buf);
		    hb_v = Or(hb_v, qv_v);
		}
		Store((var_t*) trb_buf, hb_v);
	    }	// end of n loop

	    if (ml == mc) {
const		var_t	c = vec_max(hv + wdw.lw - 3, wdw.width);
const		int	d = checkpoint(c);
		if (d < md / 2) {	// down score all
		    vec_sub_c(hv + wdw.lw - 3, c, wdw.width);
		    vec_sub_c(fv + wdw.lw - 3, c, wdw.width);
		    accscr += c;
		    mc += md;
		} else			// postpone
		    mc += d;
	    }
	}	// end of ml loop
	if (!LocalR || maxh.mr == a->right) {
	    fhlastH1(maxh, &trb);
	    maxh.val += accscr;
	}
	if (mfd) trb.traceback(maxh.mr, maxh.nr, mfd);
	return (maxh.val);
}

VTYPE SimdAln2h1::hirschbergH1_wip(Dim10* cpos, const int& n_im)
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
const	regist_v	mil_v = Splat(IntronPrm.llmt);
const	regist_v	donor_v[3] =
	{Splat(donor_code[0]), Splat(donor_code[1]), Splat(donor_code[2])};
const	regist_v	accpr_v[3] =
	{Splat(accpr_code[0]), Splat(accpr_code[1]), Splat(accpr_code[2])};
const	int	ipen = spj? pwd->IntPen->Penalty(): nevsel;
	regist_v	quant_v[IntronPrm.nquant];
	regist_v	mean_v[IntronPrm.nquant];
	if (spj) {
	    for (int j = 0; j < IntronPrm.nquant; ++j) {
		quant_v[j] = Splat(pwd->IntPen->qm[j].len);
		mean_v[j] = Splat(pwd->IntPen->qm[j].pen);
	    }
	}
	fhinitH1();

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
const	    int	j9 = std::min(nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + 3 * ml);
	    int	n9 = std::min(b->right, wdw.up + 3 * (ml + j9) + 1) + 3 * j9;
	    int	n0 = n - j8;
	    int	mp1 = ml + 1;
	    int	q = modN<6>(n + 3 * mp1);
	    int	r = n - 3 * mp1;
	    int	donor_r[3] = {r, r, r};
	    vec_set(hv_a[0], nevsel, 12 * Np1 + 3 * Nelem);
	    vec_clear(hb_a[0], 12 * Np1 + 3 * Nelem);
	    vec_clear(sm_a, 28 * Np1);
const	    bool	is_imd_ = ml == mm;

regist_v    hiv_v[3] = {ninf_v, ninf_v, ninf_v};
regist_v    hib_v[3] = {zero_v, zero_v, zero_v};
regist_v    hic_v[3] = {zero_v, zero_v, zero_v};
regist_v    hid_v[3] = {zero_v, zero_v, zero_v};
regist_v    hil_v[3] = {zero_v, zero_v, zero_v};

	    for ( ; n < n9; ++n, ++n0, ++r, q = modN<6>(q + 1)) {
const		int	p = q % 3;
const		int	rj = r - 6 * k8;
const		int	nb = std::max(0, n - b->right + 1);
const		int	kb = (nb - 1) / 3;
const		int	ke = std::min(j9, (n - b->left) / 3);

const	        SGPT6*	bb = b->exin->score_p(n);
regist_v	hb_v, fb_v, eb_v, qb_v;
regist_v	hd_v, fd_v, ed_v, qd_v;
		bool	is_imd = is_imd_ && rj >= wdw.lw && rj <= wdw.up;

		cp_a[p][0] = (b->exin->good(bb - 2))? bb[-2].sigE: 0;
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
		hv_v = Add(hv_v, ge_v);		// gap extension

		hv_a[qq][0] = hv[r + 3];
		if (LocalL) hb_a[qq][0] = hb[r + 3];
		hc_a[qq][0] = hc[r + 3];
		if (used) hd_a[qq][0] = hd[r + 3];
		qv_v = Load(hv_a[qq]);
		if (LocalL) qb_v = Load(hb_a[qq]);
		qc_v = Load(hc_a[qq]);
		if (used) qd_v = Load(hd_a[qq]);
		qv_v = Add(qv_v, g3_v);		// gap open
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
		if (nb) Store(sm_a, zero_v);
		for (int k = kb; k < ke; ++k, ++as, bs -= 3)
		    sm_a[k] = pwd->simmtx->mtx[*as][*bs];
		hv_v = Load(sm_a);

		hv_a[q][0] = hv[r];
		if (LocalL) hb_a[q][0] = hb[r];
		hc_a[q][0] = hc[r];
		if (used) hd_a[q][0] = hd[r];
regist_v	dv_v = Load(hv_a[q]);
		if (LocalL) hb_v = Load(hb_a[q]);
		hc_v = Load(hc_a[q]);
		if (used) hd_v = Load(hd_a[q]);
		hv_v = Add(hv_v, dv_v);
		hv_v = Add(hv_v, cv_v);

//	best of two
		msk_m = Cmp_gt(fv_v, hv_v);	// t >= f = !(f > t)
		hv_v = Blend(fv_v, hv_v, msk_m);
		if (LocalL) hb_v = Blend(fb_v, hb_v, msk_m);
		hc_v = Blend(fc_v, hc_v, msk_m);
		if (used) hd_v = Blend(fd_v, hd_v, msk_m);
regist_v 	pb_v = Blend(two_v, zero_v, msk_m);

//	best of three
		msk_m = Cmp_gt(ev_v, hv_v);	// t >= e = !e > t)
		hv_v = Blend(ev_v, hv_v, msk_m);
		if (LocalL) hb_v = Blend(eb_v, hb_v, msk_m);
		hc_v = Blend(ec_v, hc_v, msk_m);
		if (used) hd_v = Blend(ed_v, hd_v, msk_m);
		pb_v = Blend(one_v, pb_v, msk_m);// diag: 0, hori: 1, vert: 2
		if (is_imd) Store(pv_a[p], pb_v);

//	intron 3' boundary

regist_v	ab_v = zero_v;	// is acceptor
		if (spj) {
		  for (int k = 0; k < 2; ++k) {
const		    bool	leg = !nb && bb->phs3 > -2
			 && (!k || bb->phs3 == 2);
const		    int phase = leg? (bb->phs3 == 2? 
			(k? 1: -1): (k? 2: bb->phs3)): 2;
const		    int	pk = 2 * p + k;
		    s3_a[pk][0] = phase < 2? bb[-phase].sig3: min_ssv;
		    p3_a[pk][0] = accpr_code[phase + 1];
regist_v	    ss_v = Load(s3_a[pk]);		// donor signal
regist_v	    ph_v = Load(p3_a[pk]);
		    Store(s3_a[pk] + 1, ss_v);
		    Store(p3_a[pk] + 1, ph_v);
		    if (AllZero(ph_v)) continue;

		    for (int f = 0; f < 3; ++f) {
			qv_v = Add(hiv_v[f], ss_v);
regist_v		pv_v = mean_v[0];
			for (int j = 1; j < IntronPrm.nquant; ++j) {
			    msk_m = Cmp_gt(hil_v[f], quant_v[j - 1]);
			    pv_v = Blend(mean_v[j], pv_v, msk_m);
			}
			qv_v = Add(qv_v, pv_v);
			msk_m = Cmp_eq(ph_v, accpr_v[f]);
			qv_v = Blend(qv_v, ninf_v, msk_m);
			msk_m = Cmp_gt(hil_v[f], mil_v);	// filter 
			qv_v = Blend(qv_v, ninf_v, msk_m);	// short intron
			msk_m = Cmp_gt(qv_v, hv_v);		// is acceptor ?
			hv_v = Blend(qv_v, hv_v, msk_m);
			if (LocalL) hb_v = Blend(hib_v[f], hb_v, msk_m);
			hc_v = Blend(hic_v[f], hc_v, msk_m);
			if (used) hd_v = Blend(hid_v[f], hd_v, msk_m);
			qv_v = Blend(one_v, zero_v, msk_m);
			ab_v = Or(ab_v, qv_v);
			if (is_imd) {
			    Store(sm_a, qv_v);
		            if (sm_a[k8]) {
				imd->hlnk[0][rj] = donor_r[f];
				imd->hlnk[1][rj] = donor_r[f] + wdw.width;
				rlst[p] = rj;
			    }
			}
		    }
		  }
		}

//	Find optimal path

		if (LocalL && !accscr) {
		    msk_m = Cmp_gt(zero_v, hv_v);
		    hv_v = Blend(zero_v, hv_v, msk_m);
		}
		Store(hv_a[q] + 1, hv_v);
		if (LocalL) Store(hb_a[q] + 1, hb_v);
		Store(hc_a[q] + 1, hc_v);
		if (used) Store(hd_a[q] + 1, hd_v);
		if (LocalL && !accscr) {	// left end
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

//	intron 5' boundary

		if (spj) {
const	regist_m  nem_m = Cmp_eq(ab_v, zero_v);		// non-empty exon
		  for (int k = 0; k < 2; ++k) {		// GTGT
const		    bool	leg = !nb && bb->phs5 > -2
			 && (!k || bb->phs5 == 2);
const		    int phase = leg? (bb->phs5 == 2? 
			(k? 1: -1): (k? 2: bb->phs5)): 2;
const		    int	pk = 2 * p + k;
		    s5_a[pk][0] = phase < 2? (bb[-phase].sig5 + ipen): min_ssv;
		    p5_a[pk][0] = donor_code[phase + 1];
regist_v	    ss_v = Load(s5_a[pk]);		// donor signal
regist_v	    ph_v = Load(p5_a[pk]);
		    Store(s5_a[pk] + 1, ss_v);
		    Store(p5_a[pk] + 1, ph_v);
		    qv_v = Add(hv_v, ss_v);
		    if (AllZero(ph_v)) continue;

		    for (int f = k? 2: 0; f < 3; ++f) {
regist_v	        pv_v = (f == 2)? Add(dv_v, ss_v): qv_v;
			pv_v = Blend(pv_v, ninf_v, nem_m);	// non-empty exon
			msk_m = Cmp_eq(ph_v, donor_v[f]);
			pv_v = Blend(pv_v, ninf_v, msk_m);	// select phase
			msk_m = Cmp_gt(pv_v, hiv_v[f]);
			hiv_v[f] = Blend(pv_v, hiv_v[f], msk_m);
			hil_v[f] = Blend(zero_v, hil_v[f], msk_m);
			hil_v[f] = Add(hil_v[f], one_v);
			if (LocalL) hib_v[f] = Blend(hb_v, hib_v[f], msk_m);
			hic_v[f] = Blend(hc_v, hic_v[f], msk_m);
			if (used) hid_v[f] = Blend(hd_v, hid_v[f], msk_m);
			if (is_imd) {
			    pv_v = Blend(one_v, zero_v, msk_m);// is donor
			    Store(sm_a, pv_v);
			    if (sm_a[k8]) donor_r[f] = rj;
			}
		    }	// end of f
		  }	// end of k
		} // end of spj

//	intermediate row
		if (is_imd) {
		    Store(sm_a, ab_v);
		    if (pv_a[p][k8] == 0) rlst[p] = rj;	// diag
		    if (!sm_a[k8] && pv_a[p][k8] == 1) 	// hori
			imd->hlnk[0][rj] = rlst[p];
		    imd->vlnk[0][rj] = s2_i(hc_a[q][k9], hd_a[q][k9]);
		    hc_a[q][k9] = i_s2(hd_a[q] + k9, rj);
		    imd->vlnk[1][rj] = s2_i(fc_a[q][k9], fd_a[q][k9]);
		    fc_a[q][k9] = i_s2(fd_a[q] + k9, rj + wdw.width);
		}

//	prepare for next cycle
		int	r0 = r - 6 * j8;
		if (j9 == ke && wdw.lw <= r0 && r0 <= wdw.up) {
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
const		int	d = checkpoint(c);
		if (d < md / 2) {	// down score all
		    vec_sub_c(hv + wdw.lw - 3, c, wdw.width);
		    vec_sub_c(fv + wdw.lw - 3, c, wdw.width);
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
		for (int rp = imd->hlnk[d][r]; 
		    wdw.lw <= rp && rp < wdw.up && r != rp;
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
	for ( ; r > wdw.up; r -= wdw.width) ;
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
	if ((udhimds[++i] && udhimds[i]->mi < a->left) || cpos[i][2] < b->left)
	    maxh.val = NEVSEL;
	return (maxh.val);
}

#endif	// _FWD2H1_WIP_SIMD_H_
