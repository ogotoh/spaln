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

#include "aln.h"
#include "wln.h"
#include "dbs.h"
#include "mfile.h"
#include <math.h>

#define	MaxLevel	3

static	int	cmpwlscr(const WLUNIT* a, const WLUNIT* b);
static	int	rdiagcmp(const JUXT* a, const JUXT* b);
static	int	yposicmp(const JUXT* a, const JUXT* b);
static	JUXT*	jxtsort(JUXT* jxt, int n, int key);

static	WLPRM	wlprms[MaxLevel] = {{0}, {0}, {0}};
static	HSPPRM	hspprm = {20, 10};
static	Wlprms*	wlparams = 0;
static  int     max_no_jxts = 128;
static  int     max_stuck = 2;
static	int	play = 10;
static	const	char*	WlnDefBitPat[MaxBitPat] = {"", "1", "101", "1101",
	"11011","1101011", "110011011", "1101101011", "110010110111",
	"11101100101011", "110110010110111", "1111011001011011"};

void	makeWlprms(int dvsp)
	{if (!wlparams) wlparams = new Wlprms(dvsp);}

void	eraWlprms()
	{delete wlparams; wlparams = 0;}

WLPRM*	setwlprm(INT level)
	{return (level < MaxLevel)? wlprms + level: 0;}

/* assume the following variables are constant after once set */

void Wlprms::initilize(INT a_molc, INT b_molc, INT level)
{
/*	bitpat, redpat, elem, tpl, mask, width, gain, gain1, thr, xdrp	*/
static	WLPRM	ncprm[] = 
	{{"", 0, 4, 8, 65536, 12, 1, 4, 50, -1, 0, 0, 0},
	 {"", 0, 4, 6,  4096,  9, 1, 3, 40, -1, 0, 0, 0},
	 {"", 0, 4, 4,   256,  5, 1, 2, 30, -1, 0, 0, 0}};
static	WLPRM	aaprm[] = 
	{{"", DefConvPat[DefConvPatNo[16]], 16, 5, 1048576, 7, 4, 12, 50, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[14]], 14, 4, 38416, 5, 4, 8, 40, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[12]], 14, 3, 1728, 4, 4, 8, 30, -1, 0, 0, 0}};
//	 {0,  20, 20, 2, 400, 2, 4, 4, 30, -1, 0, 0, 0}};
static	WLPRM	trprm[] = 
	{{"", DefConvPat[DefConvPatNo[16]], 23, 4, 279841, 5, 4, 8, 50, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[14]], 14, 4, 38416, 5, 4, 8, 35, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[12]], 12, 3, 1728, 4, 4, 8, 30, -1, 0, 0, 0}};

	WLPRM*	wlprm = wlprms + level;
	WLPRM*	wlp;
	int	Ncode = ZZZ + 1;
	SEQ_CODE*	mol_code = setSeqCode(0, b_molc);

	switch (a_molc) {
	    case PROTEIN: wlp = aaprm + level; break;
	    case TRON:	wlp = trprm + level; break;
	    default:	wlp = ncprm + level; break;
	}
	if (wlprm->elem == 0) {
	    wlprm->elem = wlp->elem;
	    wlprm->redpat = wlp->redpat;
	}
	if (wlprm->tpl == 0)	wlprm->tpl = wlp->tpl;
	if (wlprm->gain == 0)	wlprm->gain = wlp->gain;
	if (wlprm->thr == 0)	wlprm->thr = wlp->thr;
	if (wlprm->xdrp == 0)	wlprm->xdrp = wlp->xdrp;
	if (wlprm->bitpat == 0)	wlprm->bitpat = wlp->bitpat;
	if (wlprm->bitpat) {
	    if (!*wlprm->bitpat) {	/* default bit pattern */
		if (wlprm->tpl < MaxBitPat)
		    wlprm->bitpat = WlnDefBitPat[wlprm->tpl];
		else
                    fatal("Ktuple must be < %d\n", MaxBitPat);
            } else if (*wlprm->bitpat != '1')  /* continuous pattern */
                wlprm->bitpat = 0;          /* else given pattern */
        }
        if (wlprm->bitpat) {
	    const char* bp = wlprm->bitpat;
	    wlprm->gain1 = wlprm->gain;
            for (wlprm->tpl = 0; *bp; ++bp) {
                if (*bp == '1') ++wlprm->tpl;
		else if (bp[-1] == '1') wlprm->gain1 += wlprm->gain;
	    }
            wlprm->width = bp - wlprm->bitpat;
            if (wlprm->width == wlprm->tpl)
		wlprm->bitpat = 0; /* continuous */
        } else
            wlprm->width = wlprm->tpl;
	if (!wlprm->ConvTab) wlprm->ConvTab = new INT[Ncode];
	for (int r = 0; r < Ncode; ++r) (wlprm->ConvTab)[r] = 127;
const	char*	ps = wlprm->redpat;
	if (!ps) ps = DefConvPat[DefConvPatNo[DvsP? 20: FourN]];
	for (wlprm->elem = 0; *ps; ++ps) {
	    if (!isalpha(*ps)) ++wlprm->elem;
	    else {
		int r = en_code(*ps, mol_code);
		if (r >= mol_code->base_code && r < mol_code->max_code)
		    (wlprm->ConvTab)[r] = wlprm->elem;
	    }
	}
	wlprm->mask = (INT) ipower(wlprm->elem, wlprm->tpl);
	wlprm->vthr = (VTYPE) (Vab * wlprm->thr);
	double	nmlfact = simmtx->avrmatch(wlprm->ConvTab);
//	double	nmlfact = simmtx->AvTrc();
	wlprm->cutoff = (int) (wlprm->gain * wlprm->vthr / nmlfact);
}

Wlprms::Wlprms(int dvsp) : DvsP(dvsp)
{
static	const	INT	a_molc[5] = {DNA, PROTEIN, UNKNOWN, PROTEIN, TRON};
static	const	INT	b_molc[5] = {DNA, TRON, UNKNOWN, PROTEIN, TRON};
	if (DvsP == 2) fatal("Ilegal combination of seq types !i\n");
	Vab = (VTYPE) alprm.scale;
	simmtx = getSimmtx(1);
	EndBonus = (VTYPE) simmtx->AvTrc();
	RepPen = (VTYPE) (Vab * hspprm.RepPen);
	for (INT level = 0; level < MaxLevel; ++level)
	    initilize(a_molc[DvsP], b_molc[DvsP], level);
	if (MaxLevel > 2)
	    wlprms[1].thr = (wlprms[0].thr + wlprms[2].thr) / 2;
}

Wlprms::~Wlprms()
{
	WLPRM*	wlp = wlprms;

	for (int i = 0; i < MaxLevel; ++i, ++wlp) {
	    delete[] wlp->ConvTab; wlp->ConvTab = 0;
	}
}

inline void resetwlprm(WLPRM* wp, int k)
{
	if (k < 0 || 20 < k) k = 20;
	wp->redpat = DefConvPat[DefConvPatNo[k]];
}

void setexprm_x(const char* ps)
{
const	char*	vl = ps + 1;
	WLPRM*	wlp0 = wlprms;
	WLPRM*	wlp2 = wlprms + MaxLevel - 1;

	switch (*ps) {
	  case 'B': wlp0->bitpat = (*vl == '0' || *vl == '-')? 0: vl; break;
	  case 'D': hspprm.DirRep = atoi(vl); break;
	  case 'G': wlp0->gain = atoi(vl); break;
	  case 'H': wlp0->thr = atoi(vl); break;
	  case 'K': wlp0->tpl = atoi(vl); break;
	  case 'R': 
	    if (isdigit(*vl)) resetwlprm(wlp0, atoi(vl));
	    else	wlp0->redpat = vl;
	    break;
	  case 'X': wlp0->xdrp = *vl? atoi(vl): 0; break;
	  case 'b': wlp2->bitpat = (*vl == '0' || *vl == '-')? 0: vl; break;
	  case 'd': hspprm.RepPen = atof(vl); break;
	  case 'g': wlp2->gain = atoi(vl); break;
	  case 'h': wlp2->thr = atoi(vl); break;
	  case 'j': max_no_jxts = atoi(vl); break;
	  case 'k': wlp2->tpl = atoi(vl); break;
	  case 'r': 
	    if (isdigit(*vl)) resetwlprm(wlp2, atoi(vl));
	    else	wlp2->redpat = vl;
	    break;
	  case 's': max_stuck = atoi(vl); break;
	  case 'x': wlp2->xdrp = *vl? atoi(vl): 0; break;
	}
}

Wlp::Wlp()
{
	a = b = 0;
	mfd = 0;
	pwd = 0;
	wlprm = 0;
	bpp = 0;
}

Wlp::Wlp(Seq* seqs[], PwdB* _pwd, INT level) : pwd(_pwd)
{
static	const	int	byte[5] = {1, 3, 0, 1, 3};
	bbt = byte[pwd->DvsP];
	a = seqs[0];
	b = seqs[1];
	wlprm = setwlprm(level);
	mfd = new Mfile(sizeof(JUXT));
	bpp = 0;
}

INT* Wlp::lookup(INT* s, int kk)
{
	INT	m = wlprm->mask;
	INT*	t = new INT[m];

	vclear(t, m);
	kk -= wlprm->width - 1;
	for (int k = 0; k++ < kk; ) {
	    if (*s < wlprm->mask) {
		m = t[*s];
		t[*s] = k;
		*s++ = m;
	    } else
		*s++ = 0;
	}
	return (t);
}

INT* Wlp::foldseq()
{
	INT	n = a->right - a->left;
	INT*	kmer = new INT[n];
	CHAR*	ps = a->at(a->left);
	CHAR*	ts = a->at(a->right);

	for (INT* s = kmer; ps < ts; ++s) {
	    INT	c = wlprm->ConvTab[*ps++];
	    if (bpp->good(c)) {
		INT	w = bpp->word(c);
		*s = bpp->flawless()? w: wlprm->mask;
	    } else {
		*s = wlprm->mask;
		bpp->flaw();
	    }
	}
	return (kmer);
}

// extend maching region if possible and calculate matching score

VTYPE Wlp::eval(JUXT* jxt)
{
	VTYPE	scr = 0;
	CHAR*	as = a->at(jxt->jx);
	CHAR*	at = a->at(min(jxt->jx + jxt->jlen, a->right));
	CHAR*	bs = b->at(jxt->jy);
	int	x = a->inex.polA == 2? a->len - a->tlen: 0;
	CHAR*	ax = a->at(x);
	CHAR*	bx = b->at(b->left);
	EXIN*	bb = 0;

	if (bbt == 3) {
	    ++bs;
            if (b->exin) {
		bb = b->exin->score(jxt->jy + 1);
		if (jxt->jx == 0 && *as == MET && scr > bb->sigS)
		    scr = wlprm->vthr / 2;
	    }
	}
	if (scr <= 0) {
	    int	lend = wlprm->tpl - jxt->jx;
	    if (a->inex.exgl && lend > 0)
		scr += wlparams->EndBonus * lend;
	}
	while (--as >= ax && (bs -= bbt) >= bx &&
	    wlprm->ConvTab[*as] == wlprm->ConvTab[*bs]) {
		jxt->jx--; jxt->jy -= bbt;	// backward extension
		if (bb) bb -= bbt;
	}
	if (as < ax) bs -= bbt;
	ax = a->at(a->right);
	bx = b->at(b->right);
	jxt->jlen = jxt->nid = 0;
	while (++as < ax && (bs += bbt) < bx) {
	    if (as >= at && wlprm->ConvTab[*as] != wlprm->ConvTab[*bs])
		break;				// forward extension
	    jxt->jlen++;
	    scr += wlparams->sim2(as, bs);
	    if (*as == *bs || (*as == SER && *bs == SER2)) jxt->nid++;
	    if (bb) {
		scr += bb->sigE;
		bb += bbt;
	    }
	}
	if (as == ax && bb && bb->sigT > 0) scr += wlprm->vthr / 2;
	else {
	    int	rend = wlprm->tpl - a->len + jxt->jx + jxt->jlen;
	    if (a->inex.exgr && rend > 0)
		scr += wlparams->EndBonus * rend;
	}
	return (scr);
}

// restore coordinates, reevaluate hsp score, remove low-score hsp

JUXT* Wlp::reeval(JUXT* jxt, int* num)
{
	int 	half = samerange(a, b);
	JUXT*	jwk = jxt;
	int	j = 0;
	int	l2fa = wlprm->width - 1;
	int	l2fb = bbt * l2fa;

	for (int i = 0; i < *num; jwk++, i++) {
	    jwk->jx += a->left - l2fa;	// start position
	    jwk->jy += b->left - l2fb;	// of the first word
	    if (!half || jwk->jx < jwk->jy) {
		jwk->jscr = eval(jwk);
		if (jwk->jscr > wlprm->vthr) jxt[j++] = *jwk;
	    }
	}
	if (j < *num) {
	    jxt[j] = jxt[*num];
	    *num = j;
	}
	return (jxt);
}

// reverse

JUXT* revjxt(JUXT* jxt, int n)
{
	JUXT*	wjx = jxt;
	JUXT*	lst = jxt + n;

	for ( ; wjx < lst; ++wjx) {
	    wjx->jx = lst->jx - wjx->jx - wjx->jlen;
	    wjx->jy = lst->jy - wjx->jy - wjx->jlen;
	}
	for (wjx =jxt; wjx < --lst; ++wjx) {
	    JUXT	tmp = *wjx;
	    *wjx = *lst;
	    *lst = tmp;
	}
	return (jxt);
}

static int scorecmp(const JUXT* a, const JUXT* b)
{
	if (b->jscr == a->jscr) return (0);
	return (b->jscr > a->jscr? 1: -1);
}

static int rdiagcmp(const JUXT* a, const JUXT* b)
{
	return (a->jx + a->jy - b->jx - b->jy);
}

static int yposicmp(const JUXT* a, const JUXT* b)
{
	if (a->jy < b->jy) return (-1);
	if (a->jy > b->jy) return (1);
	return (a->jx - b->jx);
}

static JUXT* jxtsort(JUXT* jxt, int n, int key)
{
	CMPF	cmpf;

	if (n < 2) return (jxt);
	switch (key) {
	    case ON_SCORE:	cmpf = (CMPF) scorecmp; break;
	    case ON_RDIAG:	cmpf = (CMPF) rdiagcmp; break;
	    case ON_POSIT:	cmpf = (CMPF) yposicmp; break;
	    default: return (0);
	}
	qsort((UPTR) jxt, (INT) n, sizeof(JUXT), (CMPF) cmpf);
	return (jxt);
}

void Wlp::enter(JXTD* jxtd, int r)
{
	JUXT	jbuf;

	jbuf.jx = jxtd->lastj;
	jbuf.jy = ((bbt == 3)? 3 * jxtd->lastj - 1: jxtd->lastj) + r;
	jbuf.jlen = jxtd->maxj - jxtd->lastj + wlprm->width;
	jbuf.jscr = jxtd->score;
	mfd->write(&jbuf);
}

void Wlp::dmsnno(INT* m, INT jj, INT* t)
{
	int	intvl = -wlprm->width - 1;
	int	tplwt = wlprm->tpl * wlprm->gain;
	INT	kk = b->right - b->left - wlprm->width + 1;
	CHAR*	bs = b->at(b->left);
	JXTD*	jxtd = new JXTD[jj];
	JXTD*	wxtd = jxtd;
	JXTD	ixtd = {0, 0, intvl, 0, 0};

	for (INT q = 0; q < jj; ++q, ++wxtd) *wxtd = ixtd;
	for (INT k = 0; k < kk; ) {
	    INT	c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) {
		INT	w = bpp->word(c);
		INT	j = (bpp->flawless())? t[w]: 0;
		for ( ; j; j = m[j]) {
		    int	r = k - --j;
		    int	q = (r + jj) % jj;
		    wxtd = jxtd + q;
		    intvl = j - wxtd->prevj - wlprm->width;
		    if (intvl > 0) {
			int	land = wxtd->mxscr - wlprm->cutoff;
			wxtd->score -= intvl;		// xdrop-off
			if (land > wxtd->score || wxtd->score < 0) {
			    if (land > 0) enter(wxtd, r);
			    wxtd->score = tplwt;
			    if (j < wlprm->width)		// end bonus
				wxtd->score += wlprm->gain * (wlprm->width - j);
			    wxtd->mxscr = wxtd->score;
			    wxtd->lastj = wxtd->maxj = j;
			} else wxtd->score += tplwt;
		    } else if (j - wxtd->lastj == 1)
			wxtd->score += wlprm->gain1;
		    else	wxtd->score += wlprm->gain;
		    if (!m[j]) {				// end bonus
			r = j + wlprm->width + wlprm->width - jj;
			if (r > 0) wxtd->score += wlprm->gain * r;
		    }
		    if (wxtd->score > wxtd->mxscr) {
			wxtd->mxscr = wxtd->score;
			wxtd->maxj = j;
		    }
		    wxtd->prevj = j;
		}
	    } else bpp->flaw();
	    int	r = ++k - jj;
	    int	q = k % jj;
	    wxtd = jxtd + q;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	    *wxtd = ixtd;
	}
	for (int r = kk - jj; r < (int) kk; ++r) {
	    int	q = (r + jj) % jj;
	    wxtd = jxtd + q;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	}
	delete[] jxtd;
}

void Wlp::dmsnno31(INT* m, INT jj3, INT* t)
{
	INT	jj = jj3 / 3;
	int	intvl = -wlprm->width - 1;
	int	tplwt = wlprm->tpl * wlprm->gain;
	INT	kk = b->right - b->left - 3 * (wlprm->width - 1);
	CHAR*	bs = b->at(b->left);
	JXTD*	jxtd = new JXTD[jj3];
	JXTD*	wxtd = jxtd;
	JXTD	ixtd = {0, 0, intvl, 0, 0};

	for (INT q = 0; q < jj3; ++q, ++wxtd) *wxtd = ixtd;
	INT	k = 0;
	for (INT p = 0; k < kk; ++p) {
	    if (p == 3) p = 0;
	    INT c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) {
		INT	w = bpp->word(c, p);
		INT	j = bpp->flawless(p)? t[w]: 0;
		for ( ; j; j = m[j]) {
		    int	r = k - --j * 3;
		    int	q = (r + jj3) % jj3;
		    wxtd = jxtd + q;
		    intvl = j - wxtd->prevj - wlprm->width;
		    if (intvl > 0) {
			int	land = wxtd->mxscr - wlprm->cutoff;
			wxtd->score -= intvl;
			if (land > wxtd->score || wxtd->score < 0) {
			    if (land > 0) enter(wxtd, r);	/* end of HSR */
			    wxtd->score = tplwt;
			    if (j < wlprm->width)		/* end bonus */
				wxtd->score += wlprm->gain * (wlprm->width - j);
			    wxtd->mxscr = wxtd->score;
			    wxtd->maxj = j;
			    wxtd->lastj = j;
			} else wxtd->score += tplwt;
		    } else if (j - wxtd->lastj == 1)
			wxtd->score += wlprm->gain1;
		    else	wxtd->score += wlprm->gain;
		    if (!m[j]) {				/* end bonus */
			r = j + wlprm->width + wlprm->width - jj;
			if (r > 0) wxtd->score += wlprm->gain * r;
		    }
		    if (wxtd->score > wxtd->mxscr) {
			wxtd->mxscr = wxtd->score;
			wxtd->maxj = j;
		    }
		    wxtd->prevj = j;
		}
	    } else bpp->flaw(p);
	    int	r = ++k - jj3;
	    int	q = k % jj3;
	    wxtd = jxtd + q;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	    *wxtd = ixtd;
	}	// end of k-loop
	for (int r = kk - jj3; r < (int) kk; ++r) {
	    int	q = (r + jj3) % jj3;
	    wxtd = jxtd + q;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	}
	delete[] jxtd;
}

VTYPE Wlp::LinkHspScr(HSP* mcl, HSP* ncl)
{
	int	dr = ncl->rr - mcl->rr;
	int	dd = ncl->lx - mcl->rx;
	VTYPE	scr = NEVSEL;

	if (dr < 0) dr = -dr;
	else if (dr && algmode.lsg)
	    scr = pwd->IntPen->Penalty(dr, true);
	dr /= bbt;
	VTYPE	pen = pwd->GapPenalty(dr);
	if (scr > pen) pen = scr;
	if ((dd += hspprm.DirRep) < 0) {/* overlapping segment */
	    pen += (VTYPE) (dd * wlparams->RepPen);
	}
	return (pen);
}

                                                                                      
HSP* Wlp::mkhsps(JUXT* jxt, int n)
{
	HSP*	hsp = new HSP[n];

	for (HSP* wcl = hsp; n--; ++wcl, ++jxt) {
	    wcl->lx = jxt->jx;
	    wcl->ly = jxt->jy;
	    wcl->rx = jxt->jx + jxt->jlen;
	    wcl->ry = jxt->jy + bbt * jxt->jlen;
	    wcl->rr = jxt->jy - bbt * jxt->jx;
	    wcl->nid = jxt->nid;
	    wcl->len = jxt->jlen;
	    wcl->jscr = jxt->jscr;
	    wcl->sscr = 0;
	    wcl->ulnk = 0;
	    wcl->irno = 0;
	    wcl->sumh = 0;
	    wcl->ux = INT_MAX;
	}
	return (hsp);
}

static int cmpwlscr(const WLUNIT* a, const WLUNIT* b)
{
	if (a->scr == b->scr) return (b->nid - a->nid);
	return (b->scr > a->scr)? 1: -1;
}

static int cmpwlpos(const WLUNIT* a, const WLUNIT* b)
{
	int	d = a->llmt - b->llmt;
	if (d == 0) return (a->ulmt - b->ulmt);
	return (d);
}

int JxtQueue::push(JUXT* wjx)
{
	kq = 0;
	if (qp < nrep) ++qp;
	JUXT**	mjx = 0;
	int	q = 0;
	for ( ; q < qp; ++q) {
	    JUXT*&	sjx = jqueue[q];
	    if (!sjx) {sjx = wjx; break;}
	    else if (sjx->jy + sjx->jlen - wjx->jy < hspprm.DirRep) {	// proceed
		kqueue[kq++] = sjx;
		sjx = wjx; 
		break;
	    } else if (wjx->jscr > sjx->jscr &&
		    (!mjx || (*mjx)->jscr > sjx->jscr)) mjx = &jqueue[q];
	}
	if (q == qp && mjx) *mjx = wjx;
	return (kq);
}

static void prunejxt(JUXT* jxt, int* num)
{
	JUXT*	tjx = jxt + *num - 1;
	JUXT*	djx = jxt;
	JxtQueue	hq(max_stuck);
	for (JUXT* wjx = jxt; wjx < tjx; ++wjx) {
	    if (hq.push(wjx)) {
		while (JUXT* rjx = hq.pop())
		    *djx++ = *rjx;
	    }
	}
	while (JUXT* rjx = hq.pop(true))
	    *djx++ = *rjx;
	*num = djx - jxt;
	jxtsort(jxt, *num, ON_POSIT);
}
	    
/* quadratic sparse DP */

WLUNIT* Wlp::jxtcore(int* num, JUXT** jxt)
{
	int	irno = 0;
	Mfile	wmfd(sizeof(WLUNIT));
	jxtsort(*jxt, *num, ON_POSIT);
	if (*num > max_no_jxts) prunejxt(*jxt, num);
	HSP*	ccl = mkhsps(*jxt, *num);
	HSP*	lcl = ccl + *num;
	HSP*	ncl = ccl;
	VTYPE	sumh = 0;
	int	maxclny = *num + *num;
	int*	maxx = new int[maxclny];
	vclear(maxx, maxclny);
	for ( ; ncl < lcl; ++ncl) {
	    ncl->ulnk = 0;
	    VTYPE	sscr = 0;
	    HSP*	qcl = ncl;
	    for (HSP* mcl = ncl; --mcl >= ccl; ) {
		if (ncl->rx <= mcl->rx || ncl->ry <= mcl->ry ||
		    ncl->lx <= mcl->lx || mcl->ux <= ncl->lx ||
		    (mcl->rx - ncl->lx) * 2 > ncl->rx - mcl->lx ||
		    maxx[mcl->irno] > ncl->lx + play)
		    continue;
		VTYPE	h = mcl->sscr + LinkHspScr(mcl, ncl);
		if (h > sscr) {
		    sscr = h;
		    qcl = mcl;
		}
	    }
	    ncl->sscr = sscr += ncl->jscr;
	    if (qcl != ncl) {		// link to
		ncl->ulnk = qcl;
		ncl->irno = qcl->irno;
		sumh = sscr + (qcl->sumh - qcl->sscr);
		if (qcl->ux > ncl->rx) qcl->ux = ncl->rx;
	    } else {
		ncl->irno = ++irno;
		sumh += ncl->jscr;
	    }
	    ncl->sumh = sumh;
	    if (ncl->irno >= maxclny) {
		int	n = maxclny;
		maxclny += maxclny;
		int*	t = new int[maxclny];
		vclear(t + n, n);
		vcopy(t, maxx, n);
		delete[] maxx;
		maxx = t;
	    }
	    if (ncl->rx > maxx[ncl->irno]) maxx[ncl->irno] = ncl->rx;
	}
	HSP*	qcl = ccl;
	for (ncl = ccl; ncl < lcl; ++ncl)
	    if (ncl->sscr > qcl->sscr) qcl = ncl;
	VTYPE	maxh = (!algmode.lsg && algmode.mlt < 2)?
		qcl->sscr - wlprm->vthr: wlprm->vthr;
	WLUNIT	wbf;
	wbf.num = 0;
	delete[] *jxt;
	wbf.jxt = *jxt = new JUXT[*num + irno];
	JUXT*	wjx = wbf.jxt;
	while (qcl->sscr >= maxh) {
	    for (ncl = qcl; qcl && qcl->sscr > 0; qcl = qcl->ulnk)
		wjx++;		/* count */
	    if (qcl) {		/* partial overlap */
		for (qcl = ncl; qcl && qcl->sscr > 0; qcl = qcl->ulnk)
		    qcl->sscr = 0;
		wjx = wbf.jxt;
		goto nextchain;
	    }
	    wbf.num = wjx - wbf.jxt;	/* count legitimate segments */
	    wjx->jx = a->right;
	    wjx->jy = b->right;
	    wjx->jlen = 0;
	    wbf.scr = wbf.nid = wbf.tlen = 0;
	    for (qcl = ncl; qcl; qcl = qcl->ulnk) {
		--wjx;		/* traceback */
		wjx->jx = qcl->lx;
		wjx->jy = qcl->ly;
		wjx->jlen = qcl->len;
		wjx->nid = qcl->nid;
		wjx->jscr = qcl->jscr;
		wbf.scr += qcl->jscr;
		wbf.nid += qcl->nid;
		wbf.tlen += qcl->len;
		qcl->sscr = 0;	/* erase once read */
	    }
	    wmfd.write(&wbf);
	    wbf.jxt = wjx += wbf.num + 1;
nextchain:
	    for (ncl = qcl = ccl; ncl < lcl; ++ncl)
		if (ncl->sscr > qcl->sscr) qcl = ncl;
	}
	vclear(&wbf);
	wmfd.write(&wbf);
	*num = wmfd.size() - 1;
	delete[] ccl;
	delete[] maxx;
// set lower and upper bounds
	WLUNIT*	wlu = (WLUNIT*) wmfd.flush();
	WLUNIT*	wlul = wlu;
	WLUNIT*	wlur = wlu + *num;
	for ( ; wlul < wlur; ++wlul) {
	    JUXT*	jxtr = wlul->jxt + wlul->num - 1;
	    wlul->llmt = wlul->jxt->jy;
	    wlul->ulmt = jxtr->jy + bbt * jxtr->jlen;
	}
	qsort((UPTR) wlu, (INT) *num, sizeof(WLUNIT), (CMPF) cmpwlpos);
// reset lower and upper bounds
	int	llmt = b->left;
	for (wlul = wlu; wlul < wlur; ++wlul) {
	    if (!wlul->num) continue;
	    wlul->llmt = llmt;
	    for (WLUNIT* wluu = wlul + 1; wluu < wlur; ++wluu) {
		if (wlul->ulmt < wluu->llmt) {	// in order
	    	    llmt = wlul->ulmt;
		    wlul->ulmt = wluu->llmt;
		    break;
		} else {			// overlap
		    JUXT*	jxtl = wluu->jxt;
		    JUXT*	jxtr = wluu->jxt + wluu->num;
		    while (++jxtl < jxtr) {	// find ulmt
			if (wlul->ulmt < jxtl->jy) {
			    wlul->ulmt = jxtl->jy;
			    break;
			}
		    }
		    if (jxtl == jxtr) {		// included
			if (wlul->scr >= wluu->scr) {
			    wluu->num = 0;
			    continue;
			} else {
			    wlul->num = 0;
			} 
		    }
		    jxtl = wlul->jxt;
		    jxtr = wlul->jxt + wlul->num;
		    while (--jxtr >= jxtl) {
			int	u = jxtr->jx + bbt * jxtr->jlen;
			if (u < wluu->llmt) {
			    llmt = u;
			    break;
			}
		    }
		    if (jxtr < jxtl) wlul->num = 0;
		    break;
		}
	    }
	}
	for (WLUNIT* wlul = wlu; wlul < wlur; ) {
	    if (wlul->num) ++wlul;
	    else	gswap(*wlul, *--wlur);
	}
	*num = wlur - wlu;
	wlur[-1].ulmt = b->right;
	qsort((UPTR) wlu, (INT) *num, sizeof(WLUNIT), (CMPF) cmpwlscr);
	return (wlu);
}

WLUNIT* Wlp::willip(int* nwlu, JUXT** ptop)
{
	if (a->right - a->left <= (int) wlprm->width || 
	    b->right - b->left <= (int) (bbt * wlprm->width)) {
	    *ptop = 0;
	    return (0);
	}
	bpp = new Bitpat_wq(wlprm->elem, 1, false, bitmask(wlprm->width), wlprm->bitpat);
	int	jj = a->right - a->left;
	int	jj3 = bbt * jj;
	INT*	bo = foldseq();
	INT*	t = lookup(bo, jj);
	bpp->clear();
	if (bbt == 1)	dmsnno(bo, jj, t);
	else {
	    delete bpp;
	    bpp = new Bitpat_wq(wlprm->elem, 3, false,
		bitmask(wlprm->width), wlprm->bitpat);
	    dmsnno31(bo, jj3, t);
	}
	delete[] t; delete[] bo;
	WLUNIT*	wlu = 0;
	*nwlu = mfd->size();
	if (*nwlu) {
	    JUXT	jbuf = {a->right, b->right, 0};
	    mfd->write(&jbuf);
	}
	JUXT*	jxt = (JUXT*) mfd->flush();
	delete mfd; mfd = 0;
	if (*nwlu) *ptop = reeval(jxt, nwlu);
	if (*nwlu == 0) {
	    delete[] jxt;
	    *ptop = 0;
	} else if (*nwlu == 1) {
	    wlu = new WLUNIT[2];
	    wlu->num = 1;
	    wlu->jxt = jxt;
	    wlu->scr = jxt->jscr;
	    wlu->nid = jxt->nid;
	    wlu->tlen = jxt->jlen;
	    wlu->llmt = b->left;
	    wlu->ulmt = b->right;
	    memset(wlu + 1, '\0', sizeof(WLUNIT));
	} else {
	    wlu = jxtcore(nwlu, ptop);
	}
	return (wlu);
}

Wilip::Wilip(Seq* seqs[], PwdB* pwd, INT level)
{
	Wlp	wln(seqs, pwd, level);
	wlu = wln.willip(&nwlu, &top);
}

int geneorient(Seq* seqs[], PwdB* pwd)
{
	int	ori = -1;
	Wilip*	wl[2];
	INT	level = 0;

	for ( ; level < MaxLevel; ++level) {
	    wl[0] = new Wilip(seqs, pwd, level);
	    WLUNIT*	wlu0 = wl[0]->begin();
	    antiseq(seqs + 1);            /* revers genome */
	    wl[1] = new Wilip(seqs, pwd, level);
	    WLUNIT*	wlu1 = wl[1]->begin();
	    if (wlu0 && !wlu1) {
	        antiseq(seqs + 1);
		ori = 0;
		break;
	    } else if (wlu1 && !wlu0) {
		ori = 1;
		break;
	    } else if (!wlu0 && !wlu1) {
		delete wl[0]; delete wl[1];
	        antiseq(seqs + 1);
		continue;
	    }
	    if (wlu1->scr > wlu0->scr) {
		ori = 1;
		break;
	    }
	    antiseq(seqs + 1);
	    if (wlu0->scr > wlu1->scr || level + 1 == MaxLevel) {
		ori = 0;
		break;
	    }
	    delete wl[0]; delete wl[1];
	}
	if (ori < 0) {
	    delete[] seqs[1]->jxt;
	    seqs[1]->jxt = 0;
	    seqs[1]->CdsNo = 0;
	    return (0);
	}
	WLUNIT*	wlu0 = wl[0]->begin();
	WLUNIT*	wlu1 = wl[1]->begin();
	int	n = min((int) OutPrm.MaxOut, wl[0]->size() + wl[1]->size());
	for (int k = 1; k <= n + 1; ++k) {
	    delete[] seqs[k]->jxt;
	    WLUNIT*&	wlu = (!wlu1 || (wlu0 && wlu0->scr >= wlu1->scr))?
		wlu0: wlu1;
	    if (wlu == wlu1 && k > 1) seqs[k]->comrev();
	    seqs[k]->jxt = new JUXT[wlu->num + 1];
	    memcpy(seqs[k]->jxt, wlu->jxt, (wlu->num + 1) * sizeof(JUXT));
	    seqs[k]->CdsNo = wlu->num;
	    seqs[k]->wllvl = level;
	    if (k == 1) ++k;		// skip seqs[2]
	    ++wlu;
	}
	delete wl[0]; delete wl[1];
	return (n);
}
