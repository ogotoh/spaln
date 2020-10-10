/*****************************************************************************
*
*	Quickly calculate distance between sequences 
*	Using oligomer compositions
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-4-7 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "seq.h"
#include "bitpat.h"

const	char*	DefConvPat[] = {
	"A|C|G|TU|BDHKMNRSJVWXY",			// FourN
	"C|SJTPANDEQBZ|HRK|MILV|FYW|X|U",		// Dayh6
	"ASJT|CP|DEHKNQR|FWY|G|ILMV|X|U",		// SEB6
	"ASJT|CP|DHN|EKQR|FWY|G|ILMV|X|U",		// SEB7
	"ASJT|C|DHN|EKQR|FWY|G|ILMV|P|X|U",		// SEB8
	"ASJT|C|DEN|H|KQR|FWY|G|ILMV|P|X|U",		// SEB9
	"ASJT|C|DEN|FY|G|H|ILMV|KQR|P|W|X|U",		// SEB10
	"A|C|DEN|FY|G|H|ILMV|KQR|P|SJT|W|X|U",		// SEB11
	"A|C|DN|EQ|FY|G|H|ILMV|KR|P|SJT|W|X|U",		// SEB12
	"A|C|DN|EQ|FY|G|H|IV|KR|LM|P|SJT|W|X|U",	// SEB13
	"A|C|D|EQ|FY|G|H|IV|KR|LM|N|P|SJT|W|X|U",	// SEB14
	"A|C|D|E|FY|G|H|ILMV|KR|N|P|Q|SJ|T|W|X|U",	// SEB15
	"A|C|DE|Q|F|Y|G|H|IV|KR|L|M|N|P|SJT|W|X|U",	// SEB16
	"A|C|DE|Q|F|Y|G|H|IV|K|R|L|M|N|P|SJT|W|X|U",	// SEB17
	"A|C|DE|Q|F|Y|G|H|IV|K|R|L|M|N|P|SJ|T|W|X|U",	// SEB18
	"A|C|DE|Q|F|Y|G|H|I|V|K|R|L|M|N|P|SJ|T|W|X|U",	// SEB19
	"A|R|N|D|C|Q|E|G|H|I|L|K|M|F|P|SJ|T|W|Y|V|X|U"	// EachAa
};

const int DefConvPatNo[21] = 
	{0,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

const	char*	DefBitPat[MaxBitPat] = {
"", "1", "101", "1011,10011", 
"101011,1000111", "10100111,100101101", "1010011011,1010100111",
"1001110111,100011011011", "100110110111,1010010111011",
"1001110110111,10100101011111", "100111001101111,1010011010101111",
"1000111101111011,1001110101111011", "101001101011111011,100011100011111111",
"10101001110100111111,101100011011010011111", 
"100011110001111110111,1010110010101011111011",
"101011001001011101101111,1001100001100111101111011",
};
ReducWord::ReducWord(const Seq* sd, INT elms, const char* ap) : 
	master(true), Nalpha(elms), g2r(0)
{
	molc = sd->inex.molc;

	ConvTabSize = sd->code->max_code;
	if (!sd->isdrna()) {
	    ++ConvTabSize;
	    int	aaPat = 20 - 6 + SEB6;			// 20 letters
	    if (6 <= Nalpha && Nalpha <= 20) {
		aaPat = Nalpha - 6 + SEB6;
		ap = 0;
	    } else if (Nalpha)
		fatal("Alphabet size %d is not supported !\n", Nalpha);
	    if (!ap) ap = DefConvPat[aaPat];
 	} else if (!ap) ap = DefConvPat[FourN];
	aConvTab = new CHAR[26];
	iConvTab = new CHAR[ConvTabSize];
	vset(aConvTab, CHAR(UCHAR_MAX), 26);
	for (Nalpha = 0; *ap; ap++) {
	    if (!isalpha(*ap)) ++Nalpha;
	    else {
		aConvTab[toupper(*ap) - 'A'] = Nalpha;
		int	c = en_code(*ap, sd->code);
		if (c >= sd->code->base_code && c < sd->code->ceil_code) 
		    iConvTab[c] = Nalpha;
	    }
	}
	if (!sd->isdrna())
	    aConvTab['U' - 'A'] = iConvTab[sd->code->max_code] = --Nalpha;
	a2r = sd->isprotein()? aConvTab: ntconv;
	if (sd->istron()) {
            g2r = new CHAR[64];
            for (int g = 0; g < 64; ++g)
                g2r[g] = aConvTab[acodon[gencode[g]] - 'A'];
        }
}

// least significant bit of npat = left most position

Bitpat::Bitpat(INT npat, const char* spat) : master(true)
{
	if (spat) {
	    width = strlen(spat);
	    weight = 0;
	    for (int w = 0; w < width; ++w)
		if (spat[w] == '1') ++weight;
	    exam = new int[2 * weight];
	    int	wt = 0;
	    for (int w = 0; w < width; ++w)
		if (spat[w] == '1') exam[wt++] = w;
	    for (int w = 0; w < width; ++w)
		if (spat[width - w - 1] == '1') exam[wt++] = w;
	} else {
	    width = weight = 0;
	    for (INT x = npat; x; x >>= 1) {
		if (x & 1) ++weight;
		++width;
	    }
	    if (width == weight) exam = 0;	// continuouse seed
	    else {				// discontinuouse seed
		exam = new int[2 * weight];
		INT	x = 1;
		int	wt = 0;
		for (int w = 0; w < width; ++w, x <<= 1)
		    if (npat & x) exam[wt++] = w;
		x = 1 << (width - 1);
		for (int w = 0; w < width; ++w, x >>= 1)
		    if (npat & x) exam[wt++] = w;
	    }
	}
	wshift = 2 * (weight - 1);
}

void Bitpat_wq::clear()
{
	if (exam) vset(queue, BAD_RES, qsize);
	else	vclear(queue, qsize);
	vset(fstat, msb, nframe);
	vclear(qp, noq);
}

Bitpat_wq::Bitpat_wq(INT elms, int nf, bool rvs, INT npat, const char* spat)
	: Bitpat(npat, spat), nalpha(elms), nframe(nf), reverse(rvs)
{
	if (exam) {		// discontinuouse seed
	    noq = nframe;
	} else {		// continuouse seed
	    noq = nframe > 1? 2: 1; 
	}
	msb = 1 << (weight - 1);
	qsize = noq * width;
	queue = new INT[qsize];
	fstat = new INT[nframe];
	qp = new CHAR[noq];
	clear();
	TabSize = ipower(nalpha, weight);
}

Bitpat_wq::Bitpat_wq(const Bitpat_wq& src)
	: Bitpat(src)
{
	*this = src;
	queue = new INT[qsize];
	fstat = new INT[nframe];
	qp = new CHAR[noq];
	clear();
}
	
void Bitpat_wq::flaw(int f)
{
	if (!exam) {
	    queue[f] = 0;
	    fstat[f] = msb;
	} else {
	    queue[qp[f] + f * width] = BAD_RES;
	    if (++qp[f] == width) qp[f] = 0;
	}
}

INT Bitpat_wq::word(INT c, int f)
{
	if (!exam) {
	    fstat[f] >>= 1;
	    queue[f] = (queue[f] * nalpha + c) % TabSize;
	    return (fstat[f]? BadWord: queue[f]);
	}
	int	offset = f * width;
	queue[qp[f] + offset] = c;
	if (++qp[f] == width) qp[f] = 0;
	int*	exm = reverse? exam + weight: exam;
	INT	w = 0;
	for (int k = 0; k < weight; ++k) {
	    int	q = qp[f] + exm[k];
	    if (q >= width) q -= width;
	    if (queue[q + offset] == BAD_RES) {
		fstat[f] = 1;
		return (BadWord);
	    }
	    w = w * nalpha + queue[q + offset];
	}
	fstat[f] = 0;
	return (w);
}

INT bpcompress(const char* sp, INT* w)
{
	INT	c = 0;
	if (w) *w = 0;
	for (INT b = 1; *sp; b <<= 1, ++sp)
	    if (*sp == '1') {
		c |= b;
		if (w) ++(*w);
	    } else if (*sp != '0') break;
	return (c);
}

// Word count table

WordTab::WordTab(const Seq* sd, INT tpl, INT nsft, INT elms, const char* ap,
	INT bp, INT bp2, INT nbt, INT np, INT afact, INT minorf)
	: ReducWord(sd, elms, ap), Nbitpat(nbt), kmer(tpl),
	BitPat(bp), BitPat2(bp2), Nshift(nsft), npos(np), toomany(afact), p(0),
	heads(0), wposs(0), ss(ss6[0]), word_at_0s(0), head(0), wpos(0),
	w_qp(0), wq_size(minorf), w_queue(0)
{
	bpp = new Bitpat_wq*[Nbitpat];
	nfrm = sd->istron()? 6: 1;
	word_at_0s = new INT[2 * Nbitpat];
	INT	possize = npos;
	if (nfrm == 6) {
	    toomany *= 2;
	    possize *= 2;
	    vclear(ss6, 6);
	    if (wq_size) {
		w_queue = new INT*[Nbitpat];
		*w_queue = new INT[2 * Nbitpat * wq_size];
		for (INT k = 1; k < Nbitpat; ++k)
		    w_queue[k] = w_queue[k - 1] + 2 * wq_size;
	    }
	} else	ss = wq_size = 0;
	toomany /= Nshift;
	if (Nbitpat == 1) bpp[0] = new Bitpat_wq(Nalpha, nfrm, 0, BitPat);
	else {
	    bpp[0] = new Bitpat_wq(Nalpha, nfrm, 0, bitmask(kmer));
	    for (INT k = 1; k < Nbitpat; ++k)
		bpp[k] = new Bitpat_wq(Nalpha, nfrm, (k - 1) % 2, k < 3? BitPat: BitPat2);
	}
	if (npos) {
	    if (Nbitpat == 1) {
		head = new INT[bpp[0]->TabSize];
		wpos = new INT[possize];
	    } else {
		heads = new INT*[Nbitpat];
		*heads = new INT[bpp[0]->TabSize * Nbitpat];
		wposs = new INT*[Nbitpat];
		*wposs = new INT[possize * Nbitpat];
	    }
	    for (INT k = 1; k < Nbitpat; ++k) {
		heads[k] = heads[k - 1] + bpp[0]->TabSize;
		wposs[k] = wposs[k - 1] + possize;
	    }
	}
	reset(3);
}

WordTab::WordTab(const WordTab& src)
	: ReducWord(src), Nbitpat(src.Nbitpat), kmer(src.kmer), 
	  BitPat(src.BitPat), BitPat2(src.BitPat2), Nshift(src.Nshift), 
	  npos(src.npos), toomany(src.toomany), p(0), nfrm(src.nfrm), 
	  heads(0), wposs(0), ss(ss6[0]), word_at_0s(0), head(0), wpos(0),
	  w_qp(0), wq_size(src.wq_size), w_queue(0)
{
	bpp = new Bitpat_wq*[Nbitpat];
	word_at_0s = new INT[2 * Nbitpat];
	INT	possize = npos;
	if (nfrm == 6) {
	    possize *= 2;
	    vclear(ss6, 6);
	    if (wq_size) {
		w_queue = new INT*[Nbitpat];
		*w_queue = new INT[2 * Nbitpat * wq_size];
		for (INT k = 1; k < Nbitpat; ++k)
		    w_queue[k] = w_queue[k - 1] + 2 * wq_size;
	    }
	} else	ss = 0;
	if (Nbitpat == 1) bpp[0] = new Bitpat_wq(Nalpha, nfrm, 0, BitPat);
	else {
	    bpp[0] = new Bitpat_wq(Nalpha, nfrm, 0, bitmask(kmer));
	    for (INT k = 1; k < Nbitpat; ++k)
		bpp[k] = new Bitpat_wq(Nalpha, nfrm, (k - 1) % 2, k < 3? BitPat: BitPat2);
	}
	if (npos) {
	    if (Nbitpat == 1) {
		head = new INT[bpp[0]->TabSize];
		wpos = new INT[possize];
	    } else {
		heads = new INT*[Nbitpat];
		*heads = new INT[bpp[0]->TabSize * Nbitpat];
		wposs = new INT*[Nbitpat];
		*wposs = new INT[possize * Nbitpat];
	    }
	    for (INT k = 1; k < Nbitpat; ++k) {
		heads[k] = heads[k - 1] + bpp[0]->TabSize;
		wposs[k] = wposs[k - 1] + possize;
	    }
	}
	reset(3);
}

WordTab::~WordTab()
{
	if (Nbitpat == 1) {
	    delete[] head; delete[] wpos; delete bpp[0];
	} else {
	    for (INT k = 0; k < Nbitpat; ++k) delete bpp[k];
	    if (heads) {
		delete[] *heads; delete[] heads;
	    }
	    if (wposs) {
		delete[] *wposs; delete[] wposs;
	    }
	}
	delete[] bpp; delete[] word_at_0s;
	if (w_queue) {
	    delete[] *w_queue; delete[] w_queue;
	}
}

void WordTab::save(INT w, INT i, INT k)
{
	if (Nbitpat == 1) {
	    if (i == 0) word_at_0s[0] = w;
	    else if (i > head[w]) {	// avoid cycle
		wpos[i] = head[w];
		head[w] = i;
	    } else if (i < head[w])	// something bad
		fprintf(stderr, "WordTab::save: bad order %d < %d: %d\n", i, head[w], w);
	    return;
	}
	if (i == 0) word_at_0s[k] = w;
	else if (i != heads[k][w]) {
	    wposs[k][i] = heads[k][w];
	    heads[k][w] = i;
	}
}

void WordTab::c2w(INT uc, int i)	// letter to word
{					// i: first site of word
	if (uc == Nalpha) {		// ambiguous
	    ss = 0;
	    for (INT k = 0; k < Nbitpat; ++k) bpp[k]->flaw();
	    return;
	} else	++ss;
	if (wq_size) i -= wq_size;
	for (INT k = 0; k < Nbitpat; ++k) {
	    Bitpat_wq*&	bps = bpp[k];
	    INT	w = bps->word(uc);
	    int nw = ss - bps->width;
	    if (bps->flawless() && i >= 0 && nw >= 0 && !(nw % Nshift))
		save(w, i, k);
	}
}

void WordTab::c2w6(INT uc, int i)	// letter to word
{
	cc[p] = cc[p + 3] = 0;	// initialize codon
	if (uc == 4) {		// ambiguous
	    vset(xx, 4U, 3);
	} else {
	    for (int q = 0; q < 3; ++q) {
		cc[q] = ((cc[q] << 2) + uc) & 63;
		cc[q + 3] = (((3 - uc) << 4) + (cc[q + 3] >> 2)) & 63;
		xx[p] >>= 1;
	    }
	}
	if (++p == 3) p = 0;
	int	wq_base = 0;
	int	wqp = w_qp;
	if (wq_size) {
	    if (++w_qp == wq_size) w_qp = 0;
	    i -= wq_size;
	}
	for (int q = p; q < 6; q += 3, wqp += wq_size, wq_base += wq_size) {
	    SHORT&	sss = ss6[q];
	    int	orf_len = 3 * sss;
	    if (!xx[p]) uc = g2r[cc[q]];
	    if (!xx[p] && bpp[0]->good(uc)) ++sss;
	    else	sss = 0;
	    for (INT k = 0; k < Nbitpat; ++k) {
		Bitpat_wq*&	bps = bpp[k];
		INT	w = BadWord;
		if (sss) {
		    w = bps->word(uc, q);
		    if (bps->flawless(q)) {
			int	nw = sss - bps->width;
			if (nw < 0 || nw % Nshift) w = BadWord;
		    }
		} else {		// termination codon
		    bps->flaw(q);
		    int	nw = orf_len / 3 - bps->width;
		    if (0 <= nw && orf_len < wq_size) {
			int	sp = nw % Nshift;
			nw = (nw + sp) / Nshift;
			int	qp = wqp - (sp + 1) * 3;
			nw = (nw + sp) / Nshift;
			for ( ; nw-- >= 0; qp -= 3 * Nshift) {
			    if (qp < wq_base) qp += wq_size;
			    w_queue[k][qp] = BadWord;
			}
		    }
		}
		if (wq_size) swap(w, w_queue[k][wqp]);
		if (w != BadWord && i >= 0) save(w, i, k);
	    }
	}
}

void WordTab::c2w6_pp(int bpos, int n)
{
	for (bpos -= wq_size; n-- > 0; ++bpos) {
	    int	wqp = w_qp;
	    if (++w_qp == wq_size) w_qp = 0;
	    for (int q = 0; q < 2; ++q, wqp += wq_size) {
		for (INT k = 0; k < Nbitpat; ++k) {
		    INT	w = BadWord;
		    swap(w, w_queue[k][wqp]);
		    if (w != BadWord && bpos >= 0) save(w, bpos, k);
		}
	    }
	}
}

void WordTab::reset(int mode)
{
	if (mode & 1) {
	    vclear(ss6, 6);
	    if (nfrm == 6) {
		vclear(cc, 6);
		vset(xx, 4U, 3);
	    }
	    for (INT k = 0; k < Nbitpat; ++k)
		bpp[k]->clear();
	    if (wq_size) {
		w_qp = 0;
		vset(*w_queue, BadWord, 2 * Nbitpat * wq_size);
	    }
	}
	if (mode & 2) {
	    if (head) vclear(head, bpp[0]->TabSize);
	    if (wpos) vclear(wpos, npos);
	    if (heads) vclear(*heads, bpp[0]->TabSize * Nbitpat);
	    if (wposs) vclear(*wposs, npos * Nbitpat);
	    if (word_at_0s) vset(word_at_0s, SupSiteNo, 2 * Nbitpat);
	}
}
