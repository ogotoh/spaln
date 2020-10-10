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

static	int	cmpwlscr(const WLUNIT* a, const WLUNIT* b);
static	int	rdiagcmp(const JUXT* a, const JUXT* b);
static	int	yposicmp(const JUXT* a, const JUXT* b);
static	JUXT*	jxtsort(JUXT* jxt, int n, int key);

static  const	int     max_stuck = 2;
static	WLPRM	wlprms[MaxWlpLevel + 1] = {{0}, {0}, {0}, {0}};
static	HSPPRM	hspprm = {20, 10};
static	Wlprms*	wlparams = 0;
static  int     max_no_jxts = 128;
static	INT	lutwdshft = 0;
static	const	char*	WlnDefBitPat[MaxBitPat] = {"", "1", "101", "1101",
	"11011","1101101", "110011011", "1101101011", "110010110111",
	"11101100101011", "110110010110111", "1111011001011011"};
static	const	int	too_many = 100;


void	makeWlprms(int dvsp)
	{if (!wlparams) wlparams = new Wlprms(dvsp);}

void	eraWlprms()
	{delete wlparams; wlparams = 0;}

WLPRM*	setwlprm(INT level)
	{return (level <= MaxWlpLevel)? wlprms + level: 0;}

static void fillin_wlprm(WLPRM* wlprm, int dvsp)
{
	if (wlprm->bitpat) {
	    if (!*wlprm->bitpat) {	/* default bit pattern */
		if (wlprm->tpl < MaxBitPat)
		    wlprm->bitpat = WlnDefBitPat[wlprm->tpl];
		else
		    fatal("Ktuple must be < %d\n", MaxBitPat);
	    } else if (*wlprm->bitpat != '1')  /* continuous pattern */
		wlprm->bitpat = 0;	  /* else given pattern */
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
	int	Ncode = ZZZ + 1;
	if (!wlprm->ConvTab) wlprm->ConvTab = new INT[Ncode];
	vset(wlprm->ConvTab, 127U, Ncode);
	SEQ_CODE*	mol_code = setSeqCode(0, dvsp == 3? PROTEIN: (dvsp? TRON: DNA));
const	char*	ps = wlprm->redpat;
	if (!ps) ps = DefConvPat[DefConvPatNo[dvsp? 20: FourN]];
	for (wlprm->elem = 0; *ps; ++ps) {
	    if (!isalpha(*ps)) ++wlprm->elem;
	    else {
		int r = en_code(*ps, mol_code);
		if (r >= mol_code->base_code && r < mol_code->max_code)
		    (wlprm->ConvTab)[r] = wlprm->elem;
	    }
	}
	wlprm->mask = (INT) ipower(wlprm->elem, wlprm->tpl);
	VTYPE	Vab = (VTYPE) alprm.scale;
	wlprm->vthr = (VTYPE) (Vab * wlprm->thr);
	Simmtx*	simmtx = getSimmtx(WlnPamNo);
	double	nmlfact = simmtx->avrmatch(wlprm->ConvTab);
	wlprm->cutoff = (int) (wlprm->gain * wlprm->vthr / nmlfact);
}

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
	{{"", DefConvPat[DefConvPatNo[16]], 20, 4, 279841, 5, 4, 8, 50, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[14]], 14, 4, 38416, 5, 4, 8, 40, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[12]], 12, 3, 1728, 4, 4, 8, 30, -1, 0, 0, 0}};

WLPRM*	selectwlprm(INT sz, int dvsp, WLPRM* wlp)
{
	if (!wlp) wlp = wlprms + MaxWlpLevel;
	vclear(wlp);
	if (dvsp == 0) {		// DNA vs DNA
	    wlp->tpl = INT(log((double) sz) / log(4.) + 0.5);
	    if (wlp->tpl > 8) wlp->tpl = 8;
	    if (wlp->tpl < 4) wlp->tpl = 4;
	    wlp->bitpat = WlnDefBitPat[wlp->tpl];
	    wlp->gain = 1;
	    wlp->thr = wlp->tpl * 5 + 10;
	    fillin_wlprm(wlp, dvsp);
	} else {
	    int	diff = INT_MAX;
	    wlp->elem  = 16;
	    wlp->tpl = 0;
	    for (INT a = 16; a <= 20; ++a) {	// alphabets
		INT	t = INT(log((double) sz) / log(double(a)));
		for ( ; t < 6U; ++t) {
		    INT	mask = ipower(a, t);
		    if (mask > 65538) break;
		    int	d = mask - sz;
		    if (d >= 0) {
			if (d < diff) {		// sup
			    diff = d; wlp->elem = a; wlp->tpl = t;
			}
			break;
		    }
		}
	    }
	    wlp->redpat = DefConvPat[DefConvPatNo[wlp->elem]];
	    wlp->bitpat = WlnDefBitPat[wlp->tpl];
	    wlp->gain = wlp->elem / 4;
	    wlp->thr = wlp->gain * (wlp->tpl + wlp->elem) / 2;
	    fillin_wlprm(wlp, dvsp);
	}
	return (wlp);
}
    
/* assume the following variables are constant after once set */

void Wlprms::initilize(INT level)
{
	WLPRM*	wlprm = wlprms + level;
	WLPRM*	wlp;

	switch (DvsP) {
	    case 3:	wlp = aaprm + level; break;
	    case 1:	wlp = trprm + level; break;
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
	fillin_wlprm(wlprm, DvsP);
}

Wlprms::Wlprms(int dvsp) : DvsP(dvsp)
{
	if (DvsP == 2) fatal("Ilegal combination of seq types !i\n");
	Vab = (VTYPE) alprm.scale;
	simmtx = getSimmtx(WlnPamNo);
	EndBonus = (VTYPE) simmtx->AvTrc();
	RepPen = (VTYPE) (Vab * hspprm.RepPen);
	for (INT level = 0; level < MaxWlpLevel; ++level)
	    initilize(level);
	if (MaxWlpLevel > 2)
	    wlprms[1].thr = (wlprms[0].thr + wlprms[2].thr) / 2;
}

Wlprms::~Wlprms()
{
	WLPRM*	wlp = wlprms;

	for (INT i = 0; i <= MaxWlpLevel; ++i, ++wlp) {
	    delete[] wlp->ConvTab; wlp->ConvTab = 0;
	}
}

inline void resetwlprm(WLPRM* wp, int k)
{
	if (k < 0 || 20 < k) k = 20;
	wp->redpat = DefConvPat[DefConvPatNo[k]];
}

void setexprm_x(int& argc, const char**& argv)
{
const	char&	opt = argv[0][2];
	if (!opt) return;
	bool	num = !(opt == 'B' || opt == 'R' || opt == 'b' || opt == 'r');
const	char*	vl = getarg(argc, argv, num, 3);
	if (!vl) return;
	WLPRM*	wlp0 = wlprms;
	WLPRM*	wlp2 = wlprms + MaxWlpLevel - 1;

	switch (opt) {
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
	  case 's': lutwdshft = atoi(vl); break;
	  case 'x': wlp2->xdrp = *vl? atoi(vl): 0; break;
	}
}

static	const	int	byte[5] = {1, 3, 0, 1, 3};

Wlp::Wlp(const Seq* seqs[], const PwdB* _pwd, INT level)
	: a(seqs[0]), b(seqs[1]), pwd(_pwd), bpp(0)
{
	bbt = byte[pwd->DvsP];
	wlprm = setwlprm(level);
	awspan = wlprm->width - 1;
	bwspan = 3 * wlprm->width - 1;
	mfd = new Mfile(sizeof(JUXT));
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
	INT	n = a->right - a->left - awspan;
	if (n <= 0) return (0);
	INT*	kmer = new INT[n + 1];
const	CHAR*	ps = a->at(a->left);
const	CHAR*	ts = a->at(a->left + awspan);

	while (ps < ts) {	// prelude
	    INT	c = wlprm->ConvTab[*ps++];
	    if (bpp->good(c)) bpp->word(c);
	    else	bpp->flaw();
	}
	ts = a->at(a->right);
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
	kmer[n] = kmer[0];
	return (kmer);
}

// extend maching region if possible and calculate matching score

VTYPE Wlp::eval(JUXT* jxt)
{
	VTYPE	scr = 0;
const	CHAR*	as = a->at(jxt->jx);
const	CHAR*	at = a->at(min(jxt->jx + jxt->jlen, a->right));
const	CHAR*	bs = b->at(jxt->jy);
const	CHAR*	ax = a->at(a->left);
const	CHAR*	bx = b->at(b->left);
const	EXIN*	bb = 0;

	if (bbt == 3) {
	    ++bs;
	    if (b->exin) {
		bb = (const EXIN*) b->exin->score(jxt->jy + 1);
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
	if (bbt == 3) --bx;
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
	    int	rend = wlprm->tpl - a->tlen + jxt->jx + jxt->jlen;
	    if (a->inex.exgr && rend > 0)
		scr += wlparams->EndBonus * rend;
	}
	return (scr);
}

// restore coordinates, reevaluate hsp score, remove low-score hsp

JUXT* Wlp::reeval(JUXT* jxt, int& num)
{
	int 	half = samerange(a, b);
	JUXT*	jwk = jxt;
	int	j = 0;

	for (int i = 0; i < num; jwk++, i++) {
	    jwk->jx += a->left;		// start position
	    jwk->jy += b->left;		// of the first word
	    if (!half || jwk->jx < jwk->jy) {
		jwk->jscr = eval(jwk);
		if (jwk->jscr > wlprm->vthr) jxt[j++] = *jwk;
	    }
	}
	if (j < num) {
	    jxt[j] = jxt[num];
	    num = j;
	}
	return (jxt);
}

// reverse

JUXT* revjxt(JUXT* jxt, const int& n)
{
	JUXT*	wjx = jxt;
	JUXT*	lst = jxt + n;

	for ( ; wjx < lst; ++wjx) {
	    wjx->jx = lst->jx - wjx->jx - wjx->jlen;
	    wjx->jy = lst->jy - wjx->jy - wjx->jlen;
	}
	return (vreverse(jxt, n));
}

// reverse coordinates

void Wlp::revcoord(JUXT* jxt, const int& n)
{
	JUXT*	wjx = jxt;
	JUXT*	lst = jxt + n;

	for ( ; wjx < lst; ++wjx) {
	    wjx->jx = lst->jx - wjx->jx - wjx->jlen;
	    wjx->jy = lst->jy - wjx->jy - bbt * wjx->jlen;
	}
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

void Wlp::enter(JXTD* jxtd, int r, bool on_k)
{
	JUXT	jbuf;

	if (on_k) {
	    jbuf.jx = (jxtd->lastj - r) / bbt;
	    jbuf.jy = jxtd->lastj;
	    jbuf.jlen = (jxtd->maxj - jxtd->lastj) / bbt + wlprm->width;
	} else {
	    jbuf.jx = jxtd->lastj;
	    jbuf.jy = bbt * jxtd->lastj + r;
	    jbuf.jlen = jxtd->maxj - jxtd->lastj + wlprm->width;
	}
	jbuf.jscr = jxtd->mxscr;
	mfd->write(&jbuf);
}

void Wlp::dmsnno(INT jj, INT* m, INT* t)
{
	int	intvl = -wlprm->width - 1;
	int	tplwt = wlprm->tpl * wlprm->gain;
	INT	kk = b->right - b->left - awspan;
const	CHAR*	bs = b->at(b->left);
const	CHAR*	ts = b->at(b->left + awspan);
	JXTD*	jxtd = new JXTD[jj];
	JXTD	ixtd = {0, 0, intvl, 0, 0};

	vset(jxtd, ixtd, jj);
	while (bs < ts) {
	    INT	c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) bpp->word(c);
	    else	bpp->flaw();
	}
	for (INT k = 0; k < kk; ) {
	    INT	c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) {
		INT	w = bpp->word(c);
		INT	j = (bpp->flawless())? t[w]: 0;
		for ( ; j; j = m[j]) {
		    int	r = k - --j;
		    int	q = (r + jj) % jj;
		    JXTD*	wxtd = jxtd + q;
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
	    JXTD*	wxtd = jxtd + q;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	    *wxtd = ixtd;
	}
	for (int r = kk - jj; r < (int) kk; ++r) {
	    int	q = (r + jj) % jj;
	    JXTD*	wxtd = jxtd + q;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	}
	delete[] jxtd;
}

void Wlp::dmsnno31(INT jj3, INT* m, INT* t)
{
	INT	jj = jj3 / 3;
	int	intvl = -wlprm->width - 1;
	int	tplwt = wlprm->tpl * wlprm->gain;
	INT	kk = b->right - b->left - bwspan;
const	CHAR*	bs = b->at(b->left);
const	CHAR*	ts = b->at(b->left + bwspan) - 1;
	JXTD*	jxtd = new JXTD[jj3];
	JXTD	ixtd = {0, 0, intvl, 0, 0};
	INT	p = 0;

	vset(jxtd, ixtd, jj3);
	for ( ; bs < ts; p = (p + 1) % 3) {
	    INT	c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) {
		bpp->word(c, p);
	    } else	bpp->flaw(p);
	}
	INT	k = 0;
	for ( ; k < kk; p = (p + 1) % 3) {
	    INT c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) {
		INT	w = bpp->word(c, p);
		INT	j = bpp->flawless(p)? t[w]: 0;
		for ( ; j; j = m[j]) {
		    int	r = k - --j * 3;
		    int	q = (r + jj3) % jj3;
		    JXTD*	wxtd = jxtd + q;
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
	    JXTD*	wxtd = jxtd + q;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	    *wxtd = ixtd;
	}	// end of k-loop
	for (int r = kk - jj3; r < (int) kk; ++r) {
	    int	q = (r + jj3) % jj3;
	    JXTD*	wxtd = jxtd + q;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	}
	delete[] jxtd;
}

void Wlp::dmsnno(INT jj, LookupTabs* lut, INT lb, INT rb)
{
	int	intvl = -(wlprm->width + 1);
	int	tplwt = wlprm->tpl * wlprm->gain;
	INT	Nshift = lut->nshift();
	int	neighbor = wlprm->width * Nshift;
	int	adjacent = bbt * Nshift;
	int	gain = wlprm->gain * Nshift;
	INT	jjj = jj * bbt;
	INT	n_pos = lut->blklen();
	JXTD*	jxtd = new JXTD[n_pos + jjj];
	JXTD	ixtd = {0, 0, intvl, 0, 0};
	LookupTab* lt = lut->readluts(lb, rb);
	INT 	l = 0;
	Dhash<INT, int>	hsp(max_no_jxts, INT_MIN);

	vset(jxtd, ixtd, jjj + n_pos);		// initialize
	for ( ; lb++ < rb; ++lt) {
const	  CHAR*	as = a->at(a->left);
const	  CHAR*	ts = a->at(a->left + awspan);
	  while (as < ts ) {
	    INT	c = wlprm->ConvTab[*as++];
	    if (bpp->good(c)) bpp->word(c);
	    else	bpp->flaw();
	  }
	  ts = a->at(a->right);
	  for (int j = 0; as < ts; ++j) {
	    INT	c = wlprm->ConvTab[*as++];
	    if (bpp->good(c)) {
		INT	w = bpp->word(c);
		if (!bpp->flawless()) continue;
		for (INT p = lt->header[w]; p < lt->header[w + 1]; ++p) {
		    int	k = lt->position[p];
		    int	r = k - bbt * j;
		    INT	q = r + jjj;
		    k += l; r += l;
		    JXTD*	wxtd = jxtd + q;
		    intvl = k - wxtd->prevj - neighbor;
		    if (intvl > 0) {
			int	land = wxtd->mxscr - wlprm->cutoff;
			wxtd->score -= intvl;		// xdrop-off
			if (land > wxtd->score || wxtd->score < 0) {
			    if (land > 0) enter(wxtd, r, true);
			    wxtd->score = tplwt;
			    if (j < int(wlprm->width))		// end bonus
				wxtd->score += wlprm->gain * (wlprm->width - j);
			    wxtd->mxscr = wxtd->score;
			    wxtd->lastj = wxtd->maxj = k;
			} else wxtd->score += tplwt;
		    } else if (k - wxtd->lastj == adjacent)
			wxtd->score += wlprm->gain1;
		    else	wxtd->score += gain;
		    int	eb = j + wlprm->width + wlprm->width - jj;	// end bonus
		    if (eb > 0) wxtd->score += wlprm->gain * eb;
		    if (wxtd->score > wxtd->mxscr) {
			wxtd->mxscr = wxtd->score;
			wxtd->maxj = k;
			if (wxtd->mxscr > wlprm->cutoff)
			    hsp.assign(q, r);
		    }
		    wxtd->prevj = k;
		}
	    } else bpp->flaw();
	  }
	  l += n_pos;
	  for (KVpair<INT, int>* kv = hsp.begin(); kv != hsp.end(); ++kv) {
	    JXTD*	wxtd = jxtd + kv->key;
	    intvl = l - wxtd->prevj - neighbor;
	    if (intvl > 0 && wxtd->mxscr > wlprm->cutoff)
		enter(wxtd, kv->val, true);
	  }
	  vcopy(jxtd, jxtd + n_pos, jjj);
	  vset(jxtd + jjj, ixtd, n_pos);
	}
	for (KVpair<INT, int>* kv = hsp.begin(); kv != hsp.end(); ++kv) {
	    JXTD*	wxtd = jxtd + kv->key;
	    if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, kv->val, true);
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
	if (dd < 0)	// overlapping segment 
	    pen += (mcl->jscr + ncl->jscr) * dd / (mcl->len + ncl->len);
	return (pen);
}

										      
HSP* Wlp::mkhsps(const JUXT* jxt, int n)
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

/* quadratic sparse DP */

WLUNIT* Wlp::jxtcore(int& num, JUXT** jxt)
{
	int	irno = 0;
	Mfile	wmfd(sizeof(WLUNIT));
	jxtsort(*jxt, num, ON_POSIT);
	HSP*	ccl = mkhsps(*jxt, num);
	HSP*	lcl = ccl + num;
	HSP*	ncl = ccl;
	VTYPE	sumh = 0;
	int	maxclny = num + num;
	int*	maxx = new int[maxclny];
	vclear(maxx, maxclny);
	for ( ; ncl < lcl; ++ncl) {
	    ncl->ulnk = 0;
	    VTYPE	sscr = 0;
	    HSP*	qcl = ncl;
	    for (HSP* mcl = ncl; --mcl >= ccl; ) {
		if (ncl->rx <= mcl->rx || ncl->ry <= mcl->ry ||
		    ncl->lx <= mcl->lx || mcl->ux <= ncl->lx ||
		    (mcl->rx - ncl->lx) * 2 > ncl->rx - mcl->lx)
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
	wbf.jxt = *jxt = new JUXT[num + irno];
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
	num = wmfd.size() - 1;
	delete[] ccl;
	delete[] maxx;
// set lower and upper bounds
	WLUNIT*	wlu = (WLUNIT*) wmfd.flush();
	WLUNIT*	wlul = wlu;
	WLUNIT*	wlur = wlu + num;
	for ( ; wlul < wlur; ++wlul) {
	    JUXT*	jxtr = wlul->jxt + wlul->num - 1;
	    wlul->llmt = wlul->jxt->jy;
	    wlul->ulmt = jxtr->jy + bbt * jxtr->jlen;
	}
	qsort((UPTR) wlu, (INT) num, sizeof(WLUNIT), (CMPF) cmpwlpos);
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
	    else	swap(*wlul, *--wlur);
	}
	num = wlur - wlu;
	wlur[-1].ulmt = b->right;
	qsort((UPTR) wlu, (INT) num, sizeof(WLUNIT), (CMPF) cmpwlscr);
	return (wlu);
}

JUXT* Wlp::run_dmsnno(int& njxt, LookupTabs* lut, INT lb, INT rb)
{
	if (a->right - a->left <= (int) wlprm->width || 
	    b->right - b->left <= (int) (bbt * wlprm->width)) {
	    njxt = 0;
	    return (0);
	}
	int	jj = a->right - a->left;
	int	jj3 = bbt * jj;
	if (rb > lb) {
	    bpp = new Bitpat_wq(wlprm->elem, 1, false, 
		bitmask(wlprm->width), wlprm->bitpat);
	    dmsnno(jj, lut, lb, rb);
	} else {
	    bpp = new Bitpat_wq(wlprm->elem, 1, false, 
		bitmask(wlprm->width), wlprm->bitpat);
	    INT*	bo = foldseq();
	    INT*	t = lookup(bo, jj);
	    if (bbt == 1) {
		bpp->clear();
		dmsnno(jj, bo, t);
	    } else {
		delete bpp;
		bpp = new Bitpat_wq(wlprm->elem, 3, false,
		    bitmask(wlprm->width), wlprm->bitpat);
		dmsnno31(jj3, bo, t);
	    }
	    delete[] t; delete[] bo;
	}
	JUXT*	jxt = 0;
	njxt = mfd->size();
	if (njxt) {
	    JUXT	jbuf = {a->right, b->right, 0};
	    mfd->write(&jbuf);
	    jxt = (JUXT*) mfd->flush();
	}
	delete mfd; mfd = 0;
	if (lut && a->inex.sens) revcoord(jxt, njxt);
	return (jxt);
}

WLUNIT* Wlp::willip(JUXT** ptop, int& nwlu, JUXT* jxt)
{
	WLUNIT*	wlu = 0;
	*ptop = reeval(jxt, nwlu);
	if (!nwlu) {
	    delete[] jxt;
	    *ptop = 0;
	    return (0);
	}
	if (nwlu == 1) {
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

Wilip::Wilip(const Seq* seqs[], const PwdB* pwd, INT level)
	: top(0), nwlu(0), wlu(0)
{
	Wlp	wln(seqs, pwd, level);
	JUXT*	jxt = wln.run_dmsnno(nwlu);
	if (!jxt) return;
	wlu = wln.willip(&top, nwlu, jxt);
}

Wilip::Wilip(Seq* seqs[], const PwdB* pwd, LookupTabs* lut, INT lb, INT rb)
	: top(0), nwlu(0), wlu(0)
{
	Wlp	wln((const Seq**) seqs, pwd, MaxWlpLevel);
	JUXT*	jxt = wln.run_dmsnno(nwlu, lut, lb, rb);
	if (!jxt) return;
	if (lut && seqs[0]->inex.sens) {
	    seqs[0]->comrev();
	    seqs[1]->comrev();
	}
	wlu = wln.willip(&top, nwlu, jxt);
}

// max_n must be smaller than OutPrm.MaxOut2

int geneorient(Seq* seqs[], const PwdB* pwd, int max_n)
{
	int	ori = -1;
	Wilip*	wl[2];
	INT	level = 0;

	for ( ; level < MaxWlpLevel; ++level) {
	    wl[0] = new Wilip((const Seq**) seqs, pwd, level);
	    WLUNIT*	wlu0 = wl[0]->begin();
	    antiseq(seqs + 1);	    /* revers genome */
	    wl[1] = new Wilip((const Seq**) seqs, pwd, level);
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
	    if (wlu0->scr > wlu1->scr || level + 1 == MaxWlpLevel) {
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
	max_n = min(max_n, wl[0]->size() + wl[1]->size());
	for (int k = 1; k <= max_n + 1; ++k) {
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
	return (max_n);
}

/*****************************************************
	LookupTable
*****************************************************/

MakeLookupTabs::MakeLookupTabs(const Seq* sq, const char* fname, INT blklen, 
	const WLPRM* wlp, INT afact, INT wq_size)
	: WordTab(sq, wlp->tpl, lutwdshft = lutwdshft? lutwdshft: (sq->istron()? 2: 1),
	  wlp->elem, 0, bpcompress(wlp->bitpat), 0, 1, blklen, afact, wq_size),
	  master(true), tab_size(info.tabsize), pos_size(info.possize),
	  n_pos(info.blklen), n_lut(info.n_lut), flut(-1), flux(-1)	// reference
{
	char	str[LINE_MAX];
	strcpy(str, fname);
	char* dot = strchr(str, '.');
	if (dot) *dot = '\0';
	else	dot = str + strlen(str);
	strcat(str, sq->isdrna()? LUN_EXT: LUP_EXT);
	flut = open(str, O_RDWR | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if (flut < 0) fatal("%s: open error %d\n", str, flut);
	if (write(flut, &info, sizeof(LookupTabInfo)) != sizeof(LookupTabInfo))
	    fatal(write_error, "lookup table");	// make space
	*dot = '\0';
	strcat(str, sq->isdrna()? LXN_EXT: LXP_EXT);
	flux = open(str, O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IROTH);
	if (flux < 0) fatal("%s: open error %d\n", str, flux);
	long	fpos = lseek(flut, 0L, SEEK_CUR);
	if (write(flux, &fpos, sizeof(long)) != sizeof(long))
	    fatal(write_error, "lookup index");
	info.version = lut_version;
	info.dvsp = sq->istron()? 1: 0;
	info.elem = wlp->elem;
	info.tuple = wlp->tpl;
	info.BitPat = BitPat;
	info.Nshift = lutwdshft;
	info.n_lut = info.reserve = 0;
	tab_size = wlp->mask;
	n_pos = blklen;
	pos_size = (info.dvsp + 1) * n_pos + lutwdshft - 1;
	pos_size /= lutwdshft;
	header = new SHORT[tab_size + 1];
	position = new SHORT[pos_size];
	prelude = max_width() * (sq->istron()? 3: 1) - 1;
	lutreset();
}

MakeLookupTabs::MakeLookupTabs(const MakeLookupTabs& src) :
	WordTab(src), master(false),
	info(src.info), tab_size(src.tab_size), pos_size(src.pos_size),
	n_pos(src.n_pos), n_lut(src.n_lut), flut(src.flut), 
	flux(src.flux), prelude(src.prelude)
{
	header = new SHORT[tab_size + 1];
	position = new SHORT[pos_size];
}

MakeLookupTabs::~MakeLookupTabs()
{
	delete[] header;
	delete[] position;
	info.n_lut = n_lut;

	if (flut >= 0 && master) {
	    lseek(flut, 0L, SEEK_SET);
	    if (write(flut, &info, sizeof(LookupTabInfo)) != sizeof(LookupTabInfo))
		fatal(write_error, "lookup table");
	    close(flut);
	}
	if (flux >= 0 && master) close(flux);
}

void MakeLookupTabs::store()
{
	WordTab::list2lut<SHORT>(header, position);
	if (header[tab_size]) {
	    ssize_t	t_bytes = sizeof(SHORT) * (tab_size + 1);
	    ssize_t	p_bytes = sizeof(SHORT) * header[tab_size];
	    if (write(flut, header, t_bytes) != t_bytes ||
		write(flut, position, p_bytes) != p_bytes)
		    fatal(write_error, strlut);
	}
	long	fpos = lseek(flut, 0L, SEEK_CUR);
	if (write(flux, &fpos, sizeof(long)) != sizeof(long))
	    fatal(write_error, strlut);
	++n_lut;
	WordTab::reset(2);
}

LookupTabs::LookupTabs(const char* fn, int sz)
	: fname(0), tab_size(info.tabsize), pos_size(info.possize),
	  n_pos(info.blklen), n_lut(info.n_lut), on_memory(sz <= 0),
	  prvf(0), prve(0), flut(-1), flux(-1), lutbuf(0), lutabs(0)
{
/******	prepare space for lut ******/

	n_pos = 0;
	fname = strrealloc(0, fn);
	flut = open(fname, O_RDONLY);
	if (flut < 0) return;
	readheader();
	unitsz = pos_size + tab_size + 1;	// # vars
	unitno = on_memory? n_lut: sz;
	lutno = unitsz * unitno;
	lutabs = new LookupTab[unitno];
	lutbuf = new SHORT[lutno];
	LookupTab*	lut = lutabs;
	LookupTab*	lst = lutabs + unitno;
	SHORT*	bufp = lutbuf - pos_size;
	for ( ; lut < lst; ++lut) {
	    lut->header = bufp += pos_size;
	    lut->position = bufp += tab_size + 1;
	}

/******	read lut index file ******/

	lutidx = new long[n_lut + 1];
	char	str[LINE_MAX];
	strcpy(str, fname);
	char*	dot = strrchr(str, '.');
	*dot = 0;
	strcat(str, info.dvsp? LXP_EXT: LXN_EXT);
	flux = open(str, O_RDONLY);
	if (flux < 0) return;
	ssize_t	s_idx = sizeof(long) * (n_lut + 1);
	if (read(flux, lutidx, s_idx) != s_idx)
	    fatal(read_error, str);

/******	load entire table into memory ******/

	if (on_memory) {
	    readin(0, n_lut);
	    close(flut);
	    flut = -1;
	}
}

LookupTabs::LookupTabs(LookupTabs& src)		// copy constructor
	: fname(0), info(src.info), unitsz(src.unitsz), unitno(src.unitno),
	  lutno(src.lutno), tab_size(info.tabsize), pos_size(info.possize), n_pos(info.blklen), 
	  n_lut(info.n_lut), on_memory(false), prvf(0), prve(0), flut(-1), flux(-1),
	  lutbuf(src.lutbuf), lutabs(0), lutidx(src.lutidx), basep(src.basep)
{
	if (on_memory) return;
	flut = open(src.fname, O_RDONLY);
	if (flut < 0) return;
	lutabs = new LookupTab[unitno];
	lutbuf = new SHORT[lutno];
	LookupTab*	lut = lutabs;
	LookupTab*	lst = lutabs + unitno;
	SHORT*	bufp = lutbuf - n_pos;
	for ( ; lut < lst; ++lut) {
	    lut->header = bufp += n_pos;
	    lut->position = bufp += tab_size + 1;
	}
}

void LookupTabs::readin(int lidx, int anlut)
{
	lseek(flut, lutidx[lidx], SEEK_SET);
	SHORT*	buf = lutbuf;
	while (anlut--) {
	    ssize_t	len = lutidx[lidx + 1] - lutidx[lidx];
	    if (len) {
		if (read(flut, buf, len) != len)
		    fatal(read_error, strlut);
	    } else
		vclear(buf, unitsz);
	    buf += unitsz;
	    ++lidx;
	}
}

LookupTab* LookupTabs::readluts(INT from, INT end)
{
	if (!end) end = from + 1;
	int	anlut = end - from;	// active # of luts
	if (anlut <= 0) return (0);
	if (on_memory)
	    return (lutabs + from);
	if (prvf <= from && end <= prve)
	    return (lutabs + from - prvf);
	readin(from, anlut);
	prvf = from;
	prve = end;
	return (lutabs);
}

void LookupTabs::testlut(INT from, INT end)
{
	if (end > n_lut) end = n_lut;
	for ( ; from < end; ++from) {
	    readin(from, 1);
	    int	ttl = 0;
	    int	nnz = 0;
	    for	(INT w = 0; w < tab_size; ++w) {
		int	sz = lutabs->header[w+1] - lutabs->header[w];
		ttl += sz;
		if (!sz) ++nnz;
		if (sz < 0 || sz > too_many) {
		    printf("%u: %d\n", from, sz);
		    break;
		}
	    }
	    printf("%d\t%d\t%d\n", from, nnz, ttl);
	}
}
