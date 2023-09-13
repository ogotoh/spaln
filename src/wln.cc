/*****************************************************************************
*
*	Lookup table for fast substring matching ( Wilber-Lipman algorithm)
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

#include "aln.h"
#include "wln.h"
#include <math.h>

static	int	cmpwlscr(const WLUNIT* a, const WLUNIT* b);
static	int	scorecmp(const JUXT* a, const JUXT* b);
static	int	rdiagcmp(const JUXT* a, const JUXT* b);
static	int	yposicmp(const JUXT* a, const JUXT* b);
static	int	cmphcl(const HSP** a, const HSP** b);
static	JUXT*	jxtsort(JUXT* jxt, const int& n, const int& key);

static	const	int	min_hit = 3;
static	const	int	def_afact = 10;
static	WLPRM	wlprms[MaxWlpLevel + 1] = {{0}, {0}, {0}, {0}};
static	HSPPRM	hspprm = {20, 10};
static	Wlprms*	wlparams = 0;
static	int	afact = 0;
static	const	char*	WlnDefBitPat[MaxBitPat] = {"", "1", "101", "1101",
	"11011","1101101", "110011011", "1101101011", "110010110111",
	"11101100101011", "110110010110111", "1111011001011011"};

void	makeWlprms(int dvsp)
	{if (!wlparams) wlparams = new Wlprms(dvsp);}

void	eraWlprms()
	{delete wlparams; wlparams = 0;}

WLPRM*	setwlprm(INT level)
	{return (level <= MaxWlpLevel)? wlprms + level: 0;}

static void fillin_wlprm(WLPRM* wlprm, int dvsp)
{
	if (wlprm->bitpat) {
	    if (!*wlprm->bitpat) {		// default bit pattern
		if (wlprm->tpl < MaxBitPat)
		    wlprm->bitpat = algmode.crs? WlnDefBitPat[wlprm->tpl]: 0;
		else
		    fatal("Ktuple must be < %d\n", MaxBitPat);
	    } else if (*wlprm->bitpat != '1')	// continuous pattern
		wlprm->bitpat = 0;	  	// else given pattern
	}
	wlprm->gain1 = wlprm->gain;
	if (wlprm->bitpat) {
const	    char* bp = wlprm->bitpat;
	    for (wlprm->tpl = 0; *bp; ++bp) {
		if (*bp == '1') ++wlprm->tpl;
		else if (bp[-1] == '1') wlprm->gain1 += wlprm->gain;
	    }
	    wlprm->width = bp - wlprm->bitpat;
	    if (wlprm->width == wlprm->tpl)
		wlprm->bitpat = 0;		// continuous
	} else
	    wlprm->width = wlprm->tpl;
const	int	Ncode = ZZZ + 1;
	if (!wlprm->ConvTab) wlprm->ConvTab = new INT[Ncode];
	vset(wlprm->ConvTab, 127U, Ncode);
const	SEQ_CODE* mol_code = setSeqCode(0, dvsp == 3? PROTEIN: (dvsp? TRON: DNA));
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
const	VTYPE	Vab = (VTYPE) alprm.scale;
	wlprm->vthr = (VTYPE) (Vab * wlprm->thr);
const	Simmtx*	simmtx = getSimmtx(WlnPamNo);
const	double	nmlfact = simmtx->avrmatch(wlprm->ConvTab);
	wlprm->cutoff = algmode.crs? (int) (wlprm->gain * wlprm->vthr / nmlfact):
			(3 * wlprm->gain * wlprm->tpl - 1);
}

//	bitpat, redpat, elem, tpl, mask, width, gain, gain1, thr, xdrp
static	WLPRM	ncprm[] = 
	{{"", 0, 4, 8, 65536, 12, 1, 4, 50, -1, 0, 0, 0},
	 {"", 0, 4, 6,  4096,  9, 1, 3, 40, -1, 0, 0, 0},
	 {"", 0, 4, 4,   256,  5, 1, 2, 30, -1, 0, 0, 0}};
static	WLPRM	aaprm[] = 
	{{"", DefConvPat[DefConvPatNo[16]], 16, 5, 1048575, 7, 4, 12, 50, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[14]], 14, 4, 38415, 5, 4, 8, 40, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[20]], 20, 3, 7999, 4, 4, 8, 20, -1, 0, 0, 0}};
//	 {0,  20, 20, 2, 400, 2, 4, 4, 30, -1, 0, 0, 0}};
static	WLPRM	trprm[] = 
	{{"", DefConvPat[DefConvPatNo[20]], 20, 4, 279841, 5, 4, 8, 60, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[14]], 14, 4, 38416, 5, 4, 8, 40, -1, 0, 0, 0},
//	 {"", DefConvPat[DefConvPatNo[10]], 10, 4, 9999, 4, 4, 8, 32, 0, 0, 0, 0},
//	 {"", DefConvPat[DefConvPatNo[12]], 12, 5, 248831, 5, 4, 8, 36, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[20]], 20, 3, 7999, 3, 4, 8, 20, -1, 0, 0, 0},
	 {"", DefConvPat[DefConvPatNo[10]], 10, 4, 9999, 4, 4, 8, 32, 0, 0, 0, 0}};

// assume the following variables are constant after once set

void Wlprms::initilize(INT level)
{
	WLPRM*	wlprm = wlprms + level;
	WLPRM*	wlpr;

	switch (DvsP) {
	    case 1:
		wlpr = trprm + ((algmode.crs && level == 2)? 3: level);
		break;
	    case 3:	wlpr = aaprm + level; break;
	    default:	wlpr = ncprm + level; break;
	}
	if (wlprm->elem == 0) {
	    wlprm->elem = wlpr->elem;
	    wlprm->redpat = algmode.crs? wlpr->redpat: 0;
	}
	if (wlprm->tpl == 0)	wlprm->tpl = wlpr->tpl;
	if (wlprm->gain == 0)	wlprm->gain = wlpr->gain;
	if (wlprm->thr == 0)	wlprm->thr = wlpr->thr;
	if (wlprm->xdrp == 0)	wlprm->xdrp = wlpr->xdrp;
	if (wlprm->bitpat == 0)	wlprm->bitpat = wlpr->bitpat;
	fillin_wlprm(wlprm, DvsP);
}

Wlprms::Wlprms(int dvsp) : DvsP(dvsp), 
	Vab((VTYPE) alprm.scale),
	simmtx(getSimmtx(WlnPamNo)), 
	EndBonus((VTYPE) simmtx->AvTrc() / 2), 
	RepPen((VTYPE) (Vab * hspprm.RepPen))
{
	if (DvsP == 2) fatal("Ilegal combination of seq types !i\n");
	for (INT level = 0; level < MaxWlpLevel; ++level)
	    initilize(level);
	if (MaxWlpLevel > 2)
	    wlprms[1].thr = (wlprms[0].thr + wlprms[2].thr) / 2;
}

Wlprms::~Wlprms()
{
	WLPRM*	wlpr = wlprms;

	for (INT i = 0; i <= MaxWlpLevel; ++i, ++wlpr) {
	    delete[] wlpr->ConvTab; wlpr->ConvTab = 0;
	}
}

inline void resetwlprm(WLPRM* wp, int a)
{
	if (a < 0 || 20 < a) a = 20;
	wp->redpat = algmode.crs? DefConvPat[DefConvPatNo[a]]: 0;
	wp->elem = a;
}

void setexprm_x(int& argc, const char**& argv)
{
const	char&	opt = argv[0][2];
	if (!opt) return;
	bool	num = !(opt == 'B' || opt == 'R' || opt == 'b' || opt == 'r');
const	char*	vl = getarg(argc, argv, num, 3);
	if (opt == 'a' && !vl) afact = def_afact;
	if (!vl) return;
	WLPRM*	wlp0 = wlprms;
	WLPRM*	wlp1 = wlprms + 1;

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
	  case 'a': afact = atoi(vl); break;
	  case 'b': wlp1->bitpat = (*vl == '0' || *vl == '-')? 0: vl; break;
	  case 'd': hspprm.RepPen = atof(vl); break;
	  case 'g': wlp1->gain = atoi(vl); break;
	  case 'h': wlp1->thr = atoi(vl); break;
	  case 'k': wlp1->tpl = atoi(vl); break;
	  case 'r': 
	    if (isdigit(*vl)) resetwlprm(wlp1, atoi(vl));
	    else	wlp1->redpat = vl;
	    break;
	  case 'x': wlp1->xdrp = *vl? atoi(vl): 0; break;
	}
}

static	const	int	byte[5] = {1, 3, 0, 1, 3};

Wlp::Wlp(const Seq* seqs[], const PwdB* _pwd, const int level)
	: a(seqs[0]), b(seqs[1]), pwd(_pwd), bbt(byte[pwd->DvsP]), 
	  mm(a->right - a->left), sect_l(bbt * mm), 
	  wlprm(setwlprm(std::max(level, 0))), 
	  tplwt(wlprm->tpl * wlprm->gain), 
	  awspan(wlprm->width - 1), bwspan(3 * wlprm->width - 1), 
	  precutoff(wlprm->cutoff - wlprm->gain * wlprm->tpl), 
	  ixtd{0, 0, -int(wlprm->width + 1), 0, 0},
	  min_lnkscr(-wlprm->vthr / 2),
	  dhit(algmode.crs && level > 1)
{
	if (mm <= awspan) return;
	if (level < 0 && a->len < shortquery) {
	    tmprm = *wlprm;
	    tmprm.cutoff = wlprm->cutoff * a->len / shortquery;
	    tmprm.vthr = wlprm->vthr * a->len / shortquery;
	    precutoff = precutoff * a->len / shortquery;
	    wlprm = &tmprm;
	}
	mfd = new Mfile(sizeof(JUXT));
	position = foldseq();
	header = lookup(position, mm);
}

Wlp::Wlp(Seq* seqs[]) : 
	a(seqs[0]), b(seqs[1]), mm(a->right - a->left), wlprm(setwlprm(0)), 
	awspan(wlprm->width - 1), ixtd{0, 0, 0, 0,0}, min_lnkscr(0)
{
	if (!afact) return;
	position = foldseq();
	INT*	cnt = kmercount();
	RANGE*	masked = lowic(cnt);	// low information content
	if (masked) {
	    if (masked->left > a->right - masked->right) {
		a->right = masked->left;
		b->inex.sigt = 0;
	    } else {
		a->left = masked->right;
		b->inex.sigs = 0;
	    }
	}
	delete[] cnt;
	delete[] masked;
}

INT* Wlp::lookup(INT* s, int kk)
{
	if (!s) return (0);
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

INT* Wlp::kmercount()
{
	if (mm <= awspan) return (0);
const	int	nk = mm - awspan;	// number of kmers
	INT*	count = new INT[wlprm->mask + 1];

	vclear(count, wlprm->mask + 1);
const	INT*	ts = position + nk;

	ttl = 0;
	for (const INT* ps = position; ps < ts; ++ps) {
	    if (*ps <= wlprm->mask) {
		++count[*ps];
		++ttl;
	    }
	}
	return (count);
}

INT* Wlp::foldseq()
{
const	int	nk = mm - awspan;	// number of kmers
	if (nk <= 0) return (0);
	INT*	kmer = new INT[nk + 1];
const	CHAR*	ps = a->at(a->left);
const	CHAR*	ts = a->at(a->left + awspan);

	delete bpp;
	bpp = new Bitpat_wq(wlprm->elem, 1, false, 
	    bitmask(wlprm->width), wlprm->bitpat);
	while (ps < ts) {		// prelude
const	    INT	c = wlprm->ConvTab[*ps++];
	    if (bpp->good(c)) bpp->word(c);
	    else	bpp->flaw();
	}
	ts = a->at(a->right);
	for (INT* s = kmer; ps < ts; ++s) {
const	    INT	c = wlprm->ConvTab[*ps++];
	    if (bpp->good(c)) {
		INT	w = bpp->word(c);
		*s = bpp->flawless()? w: wlprm->mask + 1;
	    } else {
		*s = wlprm->mask + 1;
		bpp->flaw();
	    }
	}
	kmer[nk] = kmer[0];
	return (kmer);
}

RANGE* Wlp::lowic(const INT* cnt)
{
	if (afact <= 0 || mm <= awspan) return (0);
const	int	nk = mm - awspan;
const	INT*	ps = position;
const	INT*	ts = position + nk;
const	INT	cthr = INT(ttl / afact);

	int	scr = 0;
	int	maxscr = 0;
	int	ml = a->left;
	int	mr = a->left;
	int	ms = a->left;
	for (int m = a->left; ps < ts; ++ps, ++m) {
	    if (*ps <= wlprm->mask && cnt[*ps] > cthr) {
		if (++scr > maxscr) {
		    maxscr = scr;
		    mr = m;
		    ml = ms;
		}
	    } else if (--scr < 0) {
		scr = 0;
		ms = m;
	    }
	}
	RANGE*	rng = 0;
	if (maxscr > int(cthr)) {
	    rng = new RANGE[1];
	    rng->left = ml + 1;
	    rng->right = mr + wlprm->width;
	}
	return (rng);
}

// extend maching region if possible and calculate matching score

VTYPE Wlp::eval(JUXT* jxt)
{
	VTYPE	scr = 0;
const	CHAR*	as = a->at(jxt->jx);
const	CHAR*	bs = b->at(jxt->jy);
const	CHAR*	at = as;
const	CHAR*	bt = bs;
const	CHAR*	ax = a->at(0);
const	CHAR*	bx = b->at(0);
const	CHAR*	az = a->at(a->tlen);
const	SGPT6*	bb = 0;

	if (bbt == 3) {
	    ++bs;
	    if (b->exin) {
		bb = (const SGPT6*) b->exin->score_p(jxt->jy + 1);
		if (bb && jxt->jx == 0 && *as == MET && bb->sigS > 0)
		    scr = wlprm->vthr / 2;
	    }
	}
	if (scr <= 0 && a->inex.exgl) {
	    int	lend = wlprm->tpl - jxt->jx;
	    if (lend > 0)
		scr += wlparams->EndBonus * std::min(lend, int(wlprm->tpl));
	}
	while (--as >= ax && (bs -= bbt) >= bx) {
	    if ((as < at || bs < bt) && 	// backward extension
		wlprm->ConvTab[*as] != wlprm->ConvTab[*bs]) break;
	    jxt->jx--; jxt->jy -= bbt;
	    ++jxt->jlen;
	    if (*as == *bs || (*as == SER && *bs == SER2)) ++jxt->nid;
	    if (bb) bb -= bbt;
	}
	if (as < ax) bs -= bbt;
	at = a->at(std::min(jxt->jx + jxt->jlen, a->right));
	bt = b->at(std::min(jxt->jy + bbt * jxt->jlen, b->right));
	ax = a->at(a->tlen);
	bx = b->at(b->len);
	if (bbt == 3) --bx;
	jxt->jlen = jxt->nid = 0;
	VTYPE	maxscr = scr;
const	CHAR*	al = as;
	int	start = 0;
	int	restart = 0;
	int	end = 0;
	int	nid = 0;
	while (++as < ax && (bs += bbt) < bx) {
	    if ((as >= at || bs >= bt) && 	// forward extension
		wlprm->ConvTab[*as] != wlprm->ConvTab[*bs]) break;
	    ++jxt->jlen;
	    scr += wlparams->sim2(as, bs);
	    if (*as == *bs || (*as == SER && *bs == SER2)) ++jxt->nid;
	    if (bb) {
		scr += bb->sigE;
		bb += bbt;
	    }
// Kadane-Gries algorithm
	    if (scr < 0) {
		scr = 0;
		restart = as - al;
		jxt->jlen = jxt->nid = 0;
	    }
	    if (scr > maxscr) {
		maxscr = scr;
		start = restart;
		end = jxt->jlen;
		nid = jxt->nid;
	    }
	}
	jxt->jx += start;
	jxt->jy += bbt * start;
	jxt->jlen = end;
	jxt->nid = nid;
// end of Kadane-Gries algorithm

	int	nmmc = std::min(jxt->jlen - jxt->nid, 3);
	if (algmode.crs == 0 && bbt == 3 && nmmc)
	    scr -= nmmc * wlprm->vthr;
	if (as == az && bb && bb->sigT > 0) scr += wlprm->vthr / 2;
	else {
const	    int	rend = wlprm->tpl - a->right + jxt->jx + jxt->jlen;
	    if (a->inex.exgr && rend > 0) {
		scr += wlparams->EndBonus * std::min(rend, int(wlprm->tpl));
	    } else if (nid == end)
		scr += wlparams->EndBonus * 4;
	}
	return (scr);
}

// restore coordinates, reevaluate hsp score, remove low-score hsp

JUXT* Wlp::reeval(JUXT* jxt, int& num)
{
const	bool 	half = samerange(a, b);
	JUXT*	jwk = jxt;
	int	j = 0;

	for (int i = 0; i < num; ++jwk, ++i) {
	    jwk->jx += a->left;		// start position
	    jwk->jy += b->left;		// of the first word
	    if (!half || jwk->jx < jwk->jy) {
		jwk->jscr = eval(jwk);
		if (jwk->jscr > wlprm->vthr)
		    jxt[j++] = *jwk;
	    }
	}
	if (j < num) {
	    jxt[j] = jxt[num];
	    num = j;
	}
	return (jxt);
}

void	Wlp::storedh(const int& r, int ml, int mr)
{
	ml -= wlprm->width;
	mr += 2 * wlprm->width;
	int	x = (r < 0)? (bbt - r - 1) / bbt: ml;
	if (x < 0) x = 0;
const	int	y = r + bbt * x + (bbt == 3? 1: 0);
const	CHAR*	as = a->at(a->left + x);
const	CHAR*	bs = b->at(b->left + y);
const	CHAR*	at = a->at(std::min(a->right, a->left + mr));
const	CHAR*	bt = b->at(b->right);
	VTYPE	scr = 0;
	VTYPE	maxscr = 0;
	int	ms = ml = mr = x;
	for (int m = x; as < at && bs < bt; ++as, bs += bbt, ++m) {
	    scr += wlparams->sim2(as, bs);
	    if (scr <= 0) {
		scr = 0;
		ms = m;
	    } else if (scr > maxscr) {
		maxscr = scr;
		mr = m;
		ml = ms;
	    }
	}
	if (maxscr > wlprm->vthr) {
	    JUXT	jbuf = {ml + 1, r + bbt * (ml + 1), mr - ml};
	    mfd->write(&jbuf);
	}
}

static int scorecmp(const JUXT* a, const JUXT* b)
{
	if (b->jscr == a->jscr) return (0);
	return (b->jscr > a->jscr? 1: -1);
}

static int rdiagcmp(const JUXT* a, const JUXT* b)
{
const	int	dr = (a->jx + a->jy - b->jx - b->jy);
	if (dr) return (dr);
	return (a->jy - b->jy);
}

static int yposicmp(const JUXT* a, const JUXT* b)
{
	if (a->jy < b->jy) return (-1);
	if (a->jy > b->jy) return (1);
	return (a->jx - b->jx);
}

static JUXT* jxtsort(JUXT* jxt, const int& n, const int& key)
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

void Wlp::enter(const JXTD* wxtd, const int& r, const bool& on_k)
{
	JUXT	jbuf;

	if (on_k) {
	    jbuf.jx = (wxtd->lastj - r) / bbt;
	    jbuf.jy = wxtd->lastj;
	    jbuf.jlen = (wxtd->maxj - wxtd->lastj) / bbt + wlprm->width;
	} else {
	    jbuf.jx = wxtd->lastj;
	    jbuf.jy = bbt * wxtd->lastj + r;
	    jbuf.jlen = wxtd->maxj - wxtd->lastj + wlprm->width;
	}
	jbuf.jscr = wxtd->mxscr;
	mfd->write(&jbuf);
}

void Wlp::scan_b(INT m, INT n)
{
	for ( ; m; m = position[m]) {
	    int	r = n - --m * bbt;
	    JXTD*	wxtd = jxtd + (r + sect_l) % sect_l;
	    if (dhit) {
		if (!wxtd->nhit) wxtd->ml = m;
		wxtd->mr = m;
		++wxtd->nhit;
		continue;
	    }
const	    int	intvl = m - wxtd->prevj - wlprm->width;
	    if (intvl > 0) {
const		int	land = wxtd->mxscr - wlprm->cutoff;
		wxtd->score -= wlprm->gain * intvl;
		if (land > wxtd->score || wxtd->score < 0) {
		    if (land > 0) enter(wxtd, r);	// end of HSP
		    wxtd->score = tplwt;
		    if (m < wlprm->width)	// end bonus
			wxtd->score += wlprm->gain * (wlprm->width - m);
		    wxtd->mxscr = wxtd->score;
		    wxtd->maxj = wxtd->lastj = m;
		} else wxtd->score += tplwt;
	    } else if (m - wxtd->lastj == 1)
		wxtd->score += wlprm->gain1;
	    else	wxtd->score += wlprm->gain;
	    if (wxtd->score > wxtd->mxscr) {
		wxtd->mxscr = wxtd->score;
		wxtd->maxj = m;
	    }
	    wxtd->prevj = m;
	}
}

void Wlp::dmsnno()
{
const	INT	nn = b->right - b->left - awspan;
const	CHAR*	bs = b->at(b->left);
const	CHAR*	ts = b->at(b->left + awspan);
	jxtd = new JXTD[mm];

	vset(jxtd, ixtd, mm);
	while (bs < ts) {
const	    INT	c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) bpp->word(c);
	    else	bpp->flaw();
	}
	for (INT n = 0; n < nn; ) {
const	    INT	c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) {
		INT	w = bpp->word(c);
		INT	m = (bpp->flawless())? header[w]: 0;
		if (m) scan_b(m, n);
	    } else bpp->flaw();
const	    int	r = ++n - mm;
	    JXTD*	wxtd = jxtd + n % mm;
	    if (dhit) {
		if (wxtd->nhit >= min_hit)
		    storedh(r, wxtd->ml, wxtd->mr);
	    } else if (wxtd->mxscr > precutoff) {
		int	d = wxtd->maxj + 2 * wlprm->width - mm;
		if (d > 0) wxtd->mxscr += wlprm->gain * d;	// end bonus
		if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	    }
	    *wxtd = ixtd;
	}
	for (int r = nn - mm; r < (int) nn; ++r) {
	    JXTD*	wxtd = jxtd + (r + mm) % mm;
	    if (dhit) {
		if (wxtd->nhit >= min_hit)
		    storedh(r, wxtd->ml, wxtd->mr);
	    } else if (wxtd->mxscr > precutoff) {
		int	d = wxtd->maxj + 2 * wlprm->width - mm;
		if (d > 0) wxtd->mxscr += wlprm->gain * d;	// end bonus
		if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	    }
	}
}

void Wlp::dmsnno31()
{
const	INT	nn = b->right - b->left - bwspan;
const	CHAR*	bs = b->at(b->left);
const	CHAR*	ts = b->at(b->left + bwspan) - 1;
	jxtd = new JXTD[sect_l];
	INT	p = 0;

	vset(jxtd, ixtd, sect_l);
	for ( ; bs < ts; p = next_p[p]) {
const	    INT	c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) {
		bpp->word(c, p);
	    } else	bpp->flaw(p);
	}
	for (INT n = 0; n < nn; p = next_p[p]) {
const	    INT c = wlprm->ConvTab[*bs++];
	    if (bpp->good(c)) {
const		INT	w = bpp->word(c, p);
		INT	m = bpp->flawless(p)? header[w]: 0;
		if (m) scan_b(m, n);
	    } else bpp->flaw(p);
const	    int	r = ++n - sect_l;
	    JXTD*	wxtd = jxtd + n % sect_l;
	    if (dhit) {
		if (wxtd->nhit >= min_hit)
		    storedh(r, wxtd->ml, wxtd->mr);
	    } else if (wxtd->mxscr > precutoff) {
		int	d = wxtd->maxj + 2 * wlprm->width - mm;
		if (d > 0) wxtd->mxscr += wlprm->gain * d;	// end bonus
		if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	    }
	    *wxtd = ixtd;
	}	// end of n-loop
	for (int r = nn - sect_l; r < (int) nn; ++r) {
	    JXTD*	wxtd = jxtd + (r + sect_l) % sect_l;
	    if (dhit) {
		if (wxtd->nhit >= min_hit)
		    storedh(r, wxtd->ml, wxtd->mr);
	    } else if (wxtd->mxscr > precutoff) {
const		int	d = wxtd->maxj + 2 * wlprm->width - mm;
		if (d > 0) wxtd->mxscr += wlprm->gain * d;	// end bonus
		if (wxtd->mxscr > wlprm->cutoff) enter(wxtd, r);
	    }
	}
}

VTYPE Wlp::LinkHspScr(HSP* mcl, HSP* ncl)
{
	int	dr = ncl->rr - mcl->rr;
const	int	dd = std::min(ncl->lx - mcl->rx, ncl->ly - mcl->ry);
	VTYPE	scr = NEVSEL;

	if (dr < 0) dr = -dr;
	else if (dr && algmode.lsg)
	    scr = pwd->IntPen->PenaltyPlus(dr);	// spliced
	dr /= bbt;
	VTYPE	pen = pwd->GapPenalty(dr);	// ordinay gap
	if (pen > scr) scr = pen;
	if (dd < 0)	// overlapping segment 
	    scr += (mcl->jscr + ncl->jscr) * dd / (mcl->len + ncl->len);
	return (scr);
}

HSP* Wlp::mkhsps(const JUXT* jxt, int& n)
{
	HSP*	hsp = new HSP[n];
	HSP*	wcl = hsp;
	HSP*	pcl = wcl;

	for (const JUXT* txt = jxt + n; jxt < txt; ++jxt) {
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
// remove imcomaptible or embraced hsps
	    if (wcl == hsp || wcl->ly > pcl->ry || 
		wcl->rx > pcl->rx || wcl->rx < pcl->lx) {
		pcl = wcl++;
		continue;
	    }
const	    int	aovr = wcl->rx - std::max(wcl->lx, pcl->lx);
	    HSP*	mcl = wcl->jscr * pcl->len > pcl->jscr * wcl->len? wcl: pcl;
const	    VTYPE	ovrscr = mcl->jscr * aovr / mcl->len +
		pwd->GapPenalty(abs(pcl->rr - wcl->rr));
	    if (ovrscr > 0) pcl = wcl++;
	    else if (wcl->jscr > pcl->jscr) vcopy(pcl, wcl, 1);
	}
	n = wcl - hsp;
	return (hsp);
}

static int cmpwlscr(const WLUNIT* a, const WLUNIT* b)
{
	if (a->scr == b->scr) return (b->nid - a->nid);
	return (b->scr > a->scr)? 1: -1;
}

static int cmpwlpos(const WLUNIT* a, const WLUNIT* b)
{
const	int	d = a->llmt - b->llmt;
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

// quadratic sparse DP

static	int	cmphcl(const HSP** a, const HSP** b) {
	return ((*b)->sscr - (*a)->sscr);
}

WLUNIT* Wlp::jxtcore(int& num, JUXT** jxt)
{
	int	irno = 0;
	Mfile	wmfd(sizeof(WLUNIT));
	jxtsort(*jxt, num, ON_RDIAG);
	HSP*	ccl = mkhsps(*jxt, num);
	HSP*	lcl = ccl + num;
	HSP*	ncl = ccl;
	VTYPE	sumh = 0;
	HSP**	phcl = new HSP*[num + 1];
	for ( ; ncl < lcl; ++ncl) {
	    ncl->ulnk = 0;
	    VTYPE	sscr = 0;
	    HSP*	qcl = ncl;
	    for (HSP* mcl = ncl; --mcl >= ccl; ) {
		if (ncl->rx <= mcl->rx || ncl->ry < mcl->ry ||
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
		if (ncl->sscr > phcl[qcl->irno]->sscr)
		    phcl[qcl->irno] = ncl;
		sumh = sscr + (qcl->sumh - qcl->sscr);
		if (qcl->ux > ncl->rx) qcl->ux = ncl->rx;
	    } else {
		phcl[irno] = ncl;
		ncl->irno = irno++;
		sumh += ncl->jscr;
	    }
	    ncl->sumh = sumh;
	}
	qsort((UPTR) phcl, (INT) irno, sizeof(HSP*), (CMPF) cmphcl);
	phcl[irno] = 0;
	HSP**	whcl = phcl;
	HSP*	qcl = *whcl++;
const	VTYPE	maxh = (!algmode.lsg && algmode.mlt < 2)?
		qcl->sscr - wlprm->vthr: wlprm->vthr;
	WLUNIT	wbf;
	wbf.num = 0;
	delete[] *jxt;
	wbf.jxt = *jxt = new JUXT[num + irno];
	JUXT*	wjx = wbf.jxt;
	while (qcl && qcl->sscr >= maxh) {
	    for (ncl = qcl; qcl && qcl->sscr > 0; qcl = qcl->ulnk)
		wjx++;		// count
	    if (qcl) {		// partial overlap
		for (qcl = ncl; qcl && qcl->sscr > 0; qcl = qcl->ulnk)
		    qcl->sscr = 0;
		wjx = wbf.jxt;
		goto nextchain;
	    }
	    wbf.num = wjx - wbf.jxt;	// count legitimate segments
	    wjx->jx = a->right;
	    wjx->jy = b->right;
	    wjx->jlen = 0;
	    wjx->jscr = 0;
	    wbf.scr = ncl->sscr;
	    wbf.nid = wbf.tlen = 0;
	    for (qcl = ncl; qcl; qcl = qcl->ulnk) {
		--wjx;		// traceback
		wjx->jx = qcl->lx;
		wjx->jy = qcl->ly;
		wjx->jlen = qcl->len;
		wjx->nid = qcl->nid;
		wjx->jscr = qcl->jscr;
		wbf.nid += qcl->nid;
		wbf.tlen += qcl->len;
		qcl->sscr = 0;	// erase once read
	    }
	    wmfd.write(&wbf);
	    wbf.jxt = wjx += wbf.num + 1;
nextchain:
	    qcl = *whcl++;
	}
	vclear(&wbf);
	wmfd.write(&wbf);
	num = wmfd.size() - 1;
	delete[] phcl;
	delete[] ccl;
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
	WLUNIT* wlum = wlu;
	for (WLUNIT* wlul = wlu; wlul < wlur; ) {
	    if (wlul->num) {
		if (wlul->ulmt > wlum->ulmt) wlum = wlul;
		++wlul;
	    } else	swap(*wlul, *--wlur);
	}
	num = wlur - wlu;
	wlum->ulmt = b->right;
	qsort((UPTR) wlu, (INT) num, sizeof(WLUNIT), (CMPF) cmpwlscr);
	return (wlu);
}

JUXT* Wlp::run_dmsnno(int& njxt)
{
	if (!mfd || 
	    a->right - a->left < (int) wlprm->width ||
	    b->right - b->left < (int) (bbt * wlprm->width)) {
	    njxt = 0;
	    return (0);
	}
	if (bbt == 1) {
	    if (bpp) bpp->clear();
	    dmsnno();
	} else {
	    delete bpp;
	    bpp = new Bitpat_wq(wlprm->elem, 3, false,
		    bitmask(wlprm->width), wlprm->bitpat);
	    dmsnno31();
	}
	JUXT*	jxt = 0;
	njxt = mfd->size();
	if (njxt) {
const	    JUXT	jbuf = {a->right - a->left, b->right - b->left, 0};
	    mfd->write(&jbuf);
	    jxt = (JUXT*) mfd->flush();
	}
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
	    vclear(wlu + 1);
	} else {
	    wlu = jxtcore(nwlu, ptop);
	}
	return (wlu);
}

Wilip::Wilip(const Seq* seqs[], const PwdB* pwd, const int level)
	: int_wlp(true)
{
	Wlp	wln(seqs, pwd, level);
	if (wln.ng()) return;
	JUXT*	jxt = wln.run_dmsnno(nwlu);
	if (!jxt) return;
	wlu = wln.willip(&top, nwlu, jxt);
}

Wilip::Wilip(const Seq* seqs[], Wlp* wln)
	: int_wlp(false)
{
	JUXT*	jxt = wln->run_dmsnno(nwlu);
	if (!jxt) return;
	wlu = wln->willip(&top, nwlu, jxt);
}

void Wilip::shift_y(int bias, int ylen)
{
	WLUNIT*	llu = wlu;	// last Wlunit
	WLUNIT*	tlu = wlu + nwlu;
	for (WLUNIT* w = wlu; w < tlu; ++w) {
	    if (bias) {
		JUXT*	txt = w->jxt + w->num;
		for (JUXT* wjx = w->jxt; wjx <= txt; ++wjx)
		    wjx->jy += bias;
		if (w->llmt) w->llmt += bias;
		w->ulmt += bias;
	    }
	    if (w->ulmt > llu->ulmt) llu = w;
	}
	llu->ulmt = ylen;
	llu->jxt[llu->num].jy = ylen;
}

void Wilip::sort_on_scr()
{
	if (nwlu > 1)
	    qsort((UPTR) wlu, (INT) nwlu, sizeof(WLUNIT), (CMPF) cmpwlscr);
}

// max_n must be smaller than OutPrm.MaxOut2

int geneorient(Seq* seqs[], const PwdB* pwd)
{
	Wilip*	wl[2];
	INT	level = 0;
	Seq**	b = seqs + 1;
	Seq**	c = seqs + 2;
	(*b)->comrev(c);
	for ( ; level < MaxWlpLevel; ++level) {
	    wl[0] = new Wilip((const Seq**) seqs, pwd, level);
	    WLUNIT*	wlu0 = wl[0]->begin();
	    swapseq(b, c);			    // revers genome
	    wl[1] = new Wilip((const Seq**) seqs, pwd, level);
	    WLUNIT*	wlu1 = wl[1]->begin();
	    swapseq(b, c);			    // revers genome
	    if (wlu0 || wlu1) break;
	    delete wl[0]; delete wl[1];
	}
	if (level == MaxWlpLevel) return (0);

	WLUNIT*	wlu0 = wl[0]->begin();
	WLUNIT*	wlu1 = wl[1]->begin();
	WLUNIT*	end0 = wl[0]->end();
	WLUNIT*	end1 = wl[1]->end();
const	bool	master_cmpl = !wlu0 || (wlu1 && wlu1->scr >= wlu0->scr);
	if (master_cmpl) swapseq(b, c);
	INT	n = 0;
	for (INT k = 1; wlu0 < end0 || wlu1 < end1; ++k) {
const	    bool	cmpl = !wlu0 || (wlu1 && wlu1->scr >= wlu0->scr);
const	    bool	rvs = cmpl != master_cmpl;
	    WLUNIT*&	wlu = cmpl? wlu1: wlu0;
	    Seq*&	sd = seqs[k];
	    if (k > 1)	sd = (rvs? *c: *b)->aliaseq(sd);
	    delete[] sd->jxt;
	    sd->jxt = new JUXT[wlu->num + 1];
	    vcopy(sd->jxt, wlu->jxt, wlu->num + 1);
	    sd->CdsNo = wlu->num;
//	    sd->left = wlu->llmt;
//	    sd->right = wlu->ulmt;
	    sd->wllvl = level;
	    if (k == 1) ++k;
	    if (++n >= OutPrm.MaxOut) break;
	    ++wlu;
	}
	delete wl[0]; delete wl[1];
	return (n);
}
