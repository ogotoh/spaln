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

#ifndef	_FWD2S1_SIMD_CC_
#define _FWD2S1_SIMD_CC_

/*************************************************************************
	nested class: Sjsites
*************************************************************************/

int SimdAln2s1::Sjsites::
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
const	Rvlujd*	brd = 0;
const	Rvlujd*	maxprd[hs1.nod];
const	bool	is_imd = hs1.imd && m == hs1.imd->mi;
	vclear(maxprd, hs1.nod);
	for (int l = 0; l <= ncand; ++l) {
const	    Rvlujd*	prd = rcd + idx[l];
const	    Rvlujd*&	mrd = maxprd[prd->dir];
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
const	    Rvlujd*    prd = maxprd[maxd];
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

void SimdAln2s1::Sjsites::
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
		Rvlujd*	prd = rcd + idx[l];
		prd->val = static_cast<var_t>(x);
		prd->ml = *hfesmb[k];
const		int	rl = hs1.mode? hs1.s2_i(*hfesmc[k], *hfesmd[k]): 0;
		if (is_imd) {
		    if (k & 1)	hs1.imd->hlnk[0][r] = hs1.rlst;
		    prd->ulk = r;
		} else if (hs1.mode) {
		    prd->ulk = rl;
		}
		prd->jnc = n;
		prd->dir = k;
	    } else --ncand;
	}
}


// prepare variables and initialize DP matrix
void SimdAln2s1::fhinitS1()
{
const	int	ru = wdw.up + 2 * nelem;
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

int SimdAln2s1::fhlastS1(Rvulmn& maxh)
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

int SimdAln2s1::
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

void SimdAln2s1::
to_spj(Queue2<int>& list, const int& m, const int& n, const int& p)
{
	for (int k = list.begin_i(); list.isnt_end(); k = list.next_i()) {
const	    int&	nj = list[k];
const	    int	j = n - nj;
	    fsjss[j]->put(j, m + j, nj, p);
	}
}

// score-only DP forward algorithm

VTYPE SimdAln2s1::
scoreonlyS1()
{
	if (algmode.alg > 1) return (scoreonlyS1_wip());
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	gn_v = Splat(pwd->BasicGEP + pwd->BasicGOP);
const	regist_v	ge2_v = Splat(pwd->LongGEP);
const	regist_v	gn2_v = Splat(pwd->LongGEP + pwd->LongGOP);
const	regist_v	two_v = Add(one_v, one_v);
const	regist_v	ninf_v = Splat(nevsel);

	fhinitS1();

	VTYPE	accscr = 0;
const	int	md = checkpoint(0);
	int	mc = md + a->left;
	for (int ml = a->left; ml < a->right; ml += nelem) {
const	    int	j9 = std::min(nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + ml);
const	    int	n9 = std::min(b->right, wdw.up + (ml + j9) + 1) + j9;
	    int	n0 = n - j8;
	    int	r = n - (ml + 1);
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * nelem);
	    vec_clear(ps_a, 2 * nelem);

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
		    vec_sub_c(hv + wdw.lw - 1, c, wdw.width);
		    vec_sub_c(fv + wdw.lw - 1, c, wdw.width);
		    if (dagp) vec_sub_c(fv2 + wdw.lw - 1, c, wdw.width);
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

VTYPE SimdAln2s1::
forwardS1(int* pp)
{
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	one_v = Splat(1);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	gn_v = Splat(pwd->BasicGEP + pwd->BasicGOP);
const	regist_v	ge2_v = Splat(pwd->LongGEP);
const	regist_v	gn2_v = Splat(pwd->LongGEP + pwd->LongGOP);
const	regist_v	two_v = Add(one_v, one_v);
const	regist_v	ninf_v = Splat(nevsel);
regist_v	hd_v, fd_v, fd2_v, qd_v; // used only when used == true

	fhinitS1();

	VTYPE	accscr = 0;
const	int	md = checkpoint(0);
	int	mc = md + a->left;

	for (int ml = a->left; ml < a->right; ml += nelem) {
const	    int	j9 = std::min(nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + ml);
const	    int	n9 = std::min(b->right, wdw.up + (ml + j9) + 1) + j9;
	    int	n0 = n - j8;
	    int	r = n - (ml + 1);
	    vec_clear(ps_a, 2 * nelem);
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * nelem);
	    vec_clear(hb_a[0], 4 * Np1 + 2 * nelem);
	    vec_clear(pb_a, nelem);

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
		    msk_m = Cmp_gt(zero_v, hv_v);	// if hv < 0
		    hv_v = Blend(zero_v, hv_v, msk_m);	// Kadane-Gries
		    hb_v = Blend(one_v, hb_v, msk_m);	// non diagonal
		    hc_v = Blend(zero_v, hc_v, msk_m);
		    if (used) hd_v = Blend(zero_v, hd_v, msk_m);
		}
		msk_m = Cmp_eq(hb_v, zero_v);
		hb_v = Blend(one_v, zero_v, msk_m);	// diag: 1, others: 0
regist_v	qb_v = Load(hb_a[p]);
		qb_v = AndNot(hb_v, qb_v);		// t && !q = nd && di
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
		    vec_sub_c(hv + wdw.lw - 1, c, wdw.width);
		    vec_sub_c(fv + wdw.lw - 1, c, wdw.width);
		    if (dagp) vec_sub_c(fv2 + wdw.lw - 1, c, wdw.width);
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

VTYPE SimdAln2s1::
hirschbergS1(Dim10* cpos, const int& n_im)
{
	Rvulmn	maxh = black_Rvulmn;
const	bool	LocalL = Local && a->inex.exgl && b->inex.exgl;
const	bool	LocalR = Local && a->inex.exgr && b->inex.exgr;
const	regist_v	zero_v = Clear();
const	regist_v	ninf_v = Splat(nevsel);
const	regist_v	one_v = Splat(1);
const	regist_v	ge_v = Splat(pwd->BasicGEP);
const	regist_v	gn_v = Splat(pwd->BasicGEP + pwd->BasicGOP);
const	regist_v	ge2_v = Splat(pwd->LongGEP);
const	regist_v	gn2_v = Splat(pwd->LongGEP + pwd->LongGOP);
const	regist_v	two_v = Add(one_v, one_v);

	fhinitS1();

	mm = (a->right - a->left + n_im) / (n_im + 1);
	Udh_Imds	udhimds(n_im, a->left, mm, wdw, pwd->Noll);
	imd = udhimds[0];
	mm = a->left + (imd->mi - a->left - 1) / nelem * nelem;
	int	k9 = imd->mi - mm;
	int	k8 = k9 - 1;

const	int	md = checkpoint(0);
	int	mc = md + a->left;
	VTYPE	accscr = 0;

regist_v	hb_v, fb_v, qb_v;	// used only when LocalL == true
regist_v	hd_v, fd_v, fd2_v, qd_v;// used only when used == true

	for (int ml = a->left, i = 0; ml < a->right; ml += nelem) {
const	    int	j9 = std::min(nelem, a->right - ml);
const	    int j8 = j9 - 1;
	    int	n  = std::max(b->left, wdw.lw + ml);
	    int	n9 = std::min(b->right, wdw.up + (ml + j9) + 1) + j9;
	    int	n0 = n - j8;
	    int	r = n - (ml + 1);
	    vec_clear(ps_a, 2 * nelem);
	    vec_set(hv_a[0], nevsel, 4 * Np1 + 2 * nelem);
	    vec_clear(hb_a[0], 4 * Np1 + 2 * nelem);
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
		    vec_sub_c(hv + wdw.lw - 1, c, wdw.width);
		    vec_sub_c(fv + wdw.lw - 1, c, wdw.width);
		    if (dagp) vec_sub_c(fv2 + wdw.lw - 1, c, wdw.width);
		    accscr += c;
		    mc += md;
		} else			// postpone
		    mc += d;
	    }
	    if (is_imd_ && ++i < n_im) {	// reset intermediate
		imd = udhimds[i];
		mm = a->left + (imd->mi - a->left - 1) / nelem * nelem;
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
	    if (wdw.lw <= imd->vlnk[d][r] && imd->vlnk[d][r] < wdw.up) {
		cpos[i][c++] = imd->mi;		// cross intermediate
		cpos[i][c++] = (d > 0)? 1: 0;
		for (int rp = imd->hlnk[d][r]; 
		    wdw.lw <= rp && rp < wdw.up && r != rp;
		    rp = imd->hlnk[d][r = rp]) {
		    cpos[i][c++] = r + imd->mi;
		}
		cpos[i][c++] = r + imd->mi;
		cpos[i][c] = end_of_ulk;
		r = imd->vlnk[d][r];
		if (r == end_of_ulk) break;
	    } else
		cpos[i][0] = end_of_ulk;	// don't cross intermediate
	}
	for ( ; r > wdw.up; r -= wdw.width) ;
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
	if ((udhimds[++i] && udhimds[i]->mi < a->left) || cpos[i][2] < b->left)
	    maxh.val = NEVSEL;
	return (maxh.val);
}

#include "fwd2s1_wip_simd.h"

#endif	// _FWD2S1_SIMD_CC_
