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
	"C|SJTPAG|NDEQBZ|HRK|MILV|FYW|X",		// Dayh6
	"ASJT|CP|DEHKNQR|FWY|G|ILMV|X",			// SEB6
	"ASJT|CP|DHN|EKQR|FWY|G|ILMV|X",		// SEB7
	"ASJT|C|DHN|EKQR|FWY|G|ILMV|P|X",		// SEB8
	"ASJT|C|DEN|H|KQR|FWY|G|ILMV|P|X",		// SEB9
	"ASJT|C|DEN|FY|G|H|ILMV|KQR|P|W|X",		// SEB10
	"A|C|DEN|FY|G|H|ILMV|KQR|P|SJT|W|X",		// SEB11
	"A|C|DN|EQ|FY|G|H|ILMV|KR|P|SJT|W|X",		// SEB12
	"A|C|DN|EQ|FY|G|H|IV|KR|LM|P|SJT|W|X",		// SEB13
	"A|C|D|EQ|FY|G|H|IV|KR|LM|N|P|SJT|W|X",		// SEB14
	"A|C|D|E|FY|G|H|ILMV|KR|N|P|Q|SJ|T|W|X",	// SEB15
	"A|C|DE|Q|F|Y|G|H|IV|KR|L|M|N|P|SJT|W|X",	// SEB16
	"A|C|DE|Q|F|Y|G|H|IV|K|R|L|M|N|P|SJT|W|X",	// SEB17
	"A|C|DE|Q|F|Y|G|H|IV|K|R|L|M|N|P|SJ|T|W|X",	// SEB18
	"A|C|DE|Q|F|Y|G|H|I|V|K|R|L|M|N|P|SJ|T|W|X",	// SEB19
	"A|R|N|BD|C|Q|EZ|G|H|I|L|K|M|F|P|SJ|T|W|Y|V|X"	// EachAa
};

const int DefConvPatNo[21] = 
	{0,1,1,1,1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};

const	char*	DefBitPat[MaxBitPat] = {"", "1", "101", "1101",
	"110101","1110101", "110011101", "1101110101", "110010110111",
	"11101100101011", "110110010110111", "1111011001011011",
	"111101100101101011", "1110101110011101101", "110111010111001011011"};

ReducWord::ReducWord(Seq* sd, INT elms, const char* ap) : Nalpha(elms)
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
	Nalpha = 0;
	for (const char* sp = ap; *sp; ++sp)
	    if (!isalpha(*sp)) ++Nalpha;
	aConvTab = new INT[26];
	iConvTab = new INT[ConvTabSize];
	vset(aConvTab, INT(SHRT_MAX), 26);
	for (Nalpha = 0; *ap; ap++) {
	    if (!isalpha(*ap)) ++Nalpha;
	    else {
		aConvTab[toupper(*ap) - 'A'] = Nalpha;
		int	c = en_code(*ap, sd->code);
		if (c >= sd->code->base_code && c < sd->code->ceil_code) 
		    iConvTab[c] = Nalpha;
	    }
	}
	if (sd->isprotein()) iConvTab[sd->code->max_code] = Nalpha;
}

// least significant bit of npat = left most position

Bitpat::Bitpat(INT elms, INT npat, const char* spat) : Nalpha(elms)
{
	if (spat) {
	    width = strlen(spat);
	    weight = 0;
	    for (int w = 0; w < width; ++w)
		if (spat[w] == '1') weight++;
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
	TabSize = ipower(Nalpha, weight);
}

void Bitpat_wq::clear()
{
	for (int q = 0; q < qsize; ++q) queue[q] = exam? BAD_RES: 0;
	for (int f = 0; f < nframe; ++f)  fstat[f] = msb;
	for (int f = 0; f < noq; ++f) qp[f] = 0;
}

Bitpat_wq::Bitpat_wq(INT elms, int nf, bool rvs, INT npat, const char* spat)
	: Bitpat(elms, npat, spat), nframe(nf), reverse(rvs)
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
	qp = new int[noq];
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
	    return (queue[f] = (queue[f] * Nalpha + c) % TabSize);
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
		return (0);
	    }
	    w = w * Nalpha + queue[q + offset];
	}
	fstat[f] = 0;
	return (w);
}

INT bpcompress(const char* sp, INT& w)
{
	INT	c = w = 0;
	for (INT b = 1; *sp; b <<= 1, ++sp)
	    if (*sp == '1') {
		c |= b;
		++w;
	    } else if (*sp != '0') break;
	return (c);
}
