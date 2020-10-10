/**************************************************************************
*
*	Boyer-Moore algorithm for text search 
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

#include "boyer_moore.h"

void BoyerMoore::bmdelta1()
{
	vset(dlt1, plen, nalph);
	int	skip = plen;
	for (int j = 0; j < plen; ++j)
	    dlt1[patt[j]] = --skip;
	if (abs(step) == 3) dlt1[SER2] = dlt1[SER];
}

// modifed from Ishibata

void BoyerMoore::bmdelta2()
{
	int*	g = new int[plen];
	int	pm = 2 * plen;
	for (int j = 0; j < plen; )
	    dlt2[j++] = --pm;
	int	j = plen;
	for (int k = plen; --k >= 0; ) {
	    g[k] = j;
	    patt[plen] = patt[k];	// sentinel
	    while (!(this->*this->eqtab)(patt[j], patt[k])) {
		dlt2[j] = min(dlt2[j], pl_1 - k);
		j = g[j];
	    }
	    --j;
	}
	sah = max((j + 1), 2) * step;
	int	s = j;
	pm = plen;
	for (j = 0; j < plen; ++j) {
	    dlt2[j] = min(dlt2[j], s + pm--);
	    if (j >= s) s = g[s];
	}
	delete[] g;
}

int BoyerMoore::nexthit(int l, int r)
{
	if (l >= 0) l = max(l - ll, 0);
	if (r >= 0) r = max(r - ll, 0);
	if (step > 0) {
	    int i = (l >= 0? l: idx[0]) + pl_1;
	    idx[0] = r = r >= 0? r: tlen;
	    while (i < r) {
		int	j = pl_1;
		while (j >= 0 && (this->*this->eqsrc)(text[i], patt[j])) {
		    i -= step; --j;
		}
		if (j < 0) {
		    idx[0] = i + sah;
		    return (i + step + ll);
		}
		i += step * max(dlt1[text[i]], dlt2[j]);
	    }
	} else {
	    int i = (r >= 0? l: idx[0]) + step * pl_1;
	    idx[0] = l = l >= 0? l: 0;
	    while (i >= l) {
		int	j = 0;
		while (j < plen && (this->*this->eqsrc)(text[i], patt[j])) {
		    i -= step, ++j;
		}
		if (j >= plen) {
		    idx[0] = i + sah;
		    return (i + step * plen + ll);
		}
		i += step * max(dlt1[text[i]], dlt2[j]);
	    }
	}
	return (-1);
}

int BoyerMoore::nexthit3(int l, int r)
{
	if (pq->size()) return (pq->shift() + ll);
	if (l >= 0) l = max(l - ll, 0);
	if (r >= 0) r = max(r - ll, 0);
	for (int k = 0; k < 3; ++k) {
	    if (step > 0) {
		int i = l >= 0? l + k: idx[k];
		idx[k] = r = (r >= 0? r: tlen) + k;	// centinel
		for (i += step * pl_1; i < r; ) {
		    int	j = pl_1;
		    while (j >= 0 && (this->*this->eqsrc)(text[i], patt[j])) {
			i -= step; --j;
		    }
		    if (j < 0) {
			pq->put(i + step);
			idx[k] = i + sah;
			break;
		    }
		    i += step * max(dlt1[text[i]], dlt2[j]);
		}
	    } else {
		int i = r >= 0? r - k: idx[k];
		idx[k] = l = (l >= 0? l: 0) - k;		// centinel
		for (i += step * pl_1; i >= l; ) {
		    int	j = 0;
		    while (j < plen && (this->*this->eqsrc)(text[i], patt[j])) {
			i -= step, ++j;
		    }
		    if (j >= plen) {
			pq->put(i + step * plen);
			idx[k] = i + sah;
			break;
		    }
		    i += step * max(dlt1[text[i]], dlt2[j]);
		}
	    }
	}
	return (pq->size()? pq->shift() + ll: -1);
}

BoyerMoore::BoyerMoore(const char* t, const char* p, int dir) :
	text((const unsigned char*) t), tlen((int) strlen(t)), 
	plen((int) strlen(p)), pl_1(plen - 1), ll(0), nalph(128), 
	step(dir), sah(2 * step), found(0), nfound(1), pq(0)
{
	eqsrc = &BoyerMoore::eq;
	dlt1 = new int[nalph];
	dlt2 = new int[plen];
	patt = new CHAR[pl_1];
	memcpy(patt, p, plen + 1);
	if (step < 0) vreverse(patt, plen);
	bmdelta1();
	bmdelta2();
	if (step < 0) {
	    vreverse(patt, plen);
	    vreverse(dlt2, plen);
	}
}

BoyerMoore::BoyerMoore(const Seq* sd, const Seq* pd, int s) :
	text(sd->at(sd->left)), tlen(sd->right - sd->left),
	plen(pd->right - pd->left), pl_1(plen - 1), ll(sd->left),
	nalph(sd->code->max_code), step(s), sah(2 * s),
	found(0), nfound(abs(s)), pq(0)
{
	if (sd->istron()) {
	    eqsrc = pd->isprotein()? &BoyerMoore::eqta: &BoyerMoore::eqtn;
	    eqtab = pd->isprotein()? &BoyerMoore::eqa: &BoyerMoore::eqn;
	} else if (sd->isprotein() && pd->isprotein()) 
	    eqsrc = eqtab = &BoyerMoore::eqa;
	else if (sd->isdrna() && pd->isdrna())
	    eqsrc = eqtab = &BoyerMoore::eqn;
	else	fatal("compare distince molecular types %d vs %d!\n",
		sd->inex.molc, pd->inex.molc);
	dlt1 = new int[nalph];
	dlt2 = new int[plen];
	patt = new CHAR[plen + 1];
	memcpy(patt, pd->at(pd->left), plen);
	if (step < 0) vreverse(patt, plen);
	bmdelta1();
	bmdelta2();
	if (step < 0) {
	    vreverse(patt, plen);
	    vreverse(dlt2, plen);
	}
	if (nfound > 1) {
	    found = new int[nfound + 1];
	    pq = new PrQueue<int>(found, nfound, 0, step < 0, true);
	}
	idx[0] = step > 0? 0: tlen;
	if (step > 0) 
	     for (int k = 1; k < 3; ++k) idx[k] = nfound == 3? k: tlen;
	else 
	     for (int k = 1; k < 3; ++k) idx[k] = nfound == 3? tlen - k: 0;
}

#if MAIN

int main(int argc, const char** argv)
{
const	char	text[] = "abcdabcbcabcad";
const	char	patt[] = "abcbcabc";
	int	step = 1;
	if (argc > 1 && argv[1][0] == '-') {
	    --argc; ++argv;
	    step = -1;
	}
const	char*	pat = (argc > 1)? *++argv: patt;
	BoyerMoore	bm(text, pat, step);
	int	h = -1;
	int	l = 0;
	int	r = (int) strlen(text);
	if (step > 0) {
	    while ((h = bm.nexthit(l, r)) >= 0) {
		printf("%d ", h);
		l = h + bm.shift_after_hit();
	    }
	} else {
	    while ((h = bm.nexthit(l, r)) >= 0) {
		printf("%d ", h);
		r = h - bm.shift_after_hit() + 1;
	    }
	}
	putchar('\n');
	return (0);
}

#endif	// MAIN
