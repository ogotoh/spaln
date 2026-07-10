/*****************************************************************************
*
*	Word frequencies in database
*	which must have been formated by makdbs	
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

#include "seq.h"
#include "dbs.h"
#include "eijunc.h"
#include <math.h>

struct Header {
	INT	lnwords;
	INT	snwords;
	SHORT	letter;
	SHORT	maxk;
	SHORT	mink;
	SHORT	max_sk;
};

class Kmers {
protected:
	INT	supw = 0;
	int	molc = DNA;
	long*	scounts = 0;	// short kmers
	INT*	lcounts = 0;	// long kmers
	INT*	nok = 0;
	float	cvrate = 0;
const	char*	decoder;
template <typename file_t>
	int	read_binary(file_t fd, const char* fname);
template <typename file_t>
	int	read_text(file_t fd, const char* fname);
template <typename file_t>
	void	write_binary(file_t fd, const char* fname) {
	    if (fwrite(&hdr, sizeof(Header), 1, fd) != 1)
		fatal(write_error, fname);
	    if (scounts && 
		(fwrite(scounts, sizeof(long), hdr.snwords, fd) != hdr.snwords))
		fatal(write_error, fname);
	    if (fwrite(lcounts, sizeof(INT), hdr.lnwords, fd) != hdr.lnwords)
		fatal(write_error, fname);
	}
	void	make_nok();
template <typename file_t>
	void	write_text(file_t fd);
public:
	Header	hdr = {0, 0, 0, 0, 0, 0};
	long**	sresults = 0;
	INT**	lresults = 0;
	Kmers() {}
	Kmers(INT a, INT xkk, INT nkk = 0);
	Kmers(const char* fname);
	~Kmers() {
	    delete[] scounts; delete[] sresults;
	    delete[] lcounts; delete[] lresults;
	    delete[] nok;
	}
	void	outputCount(const char* bout = 0, const int& text = 2);
	void	get(INT* dest, const INT& j) const {
	    vcopy(dest, lresults[j], nok[j]);
	}
	void	sget(long* dest, INT j) const {
	    vcopy(dest, sresults[j], nok[j]);
	}
const	long* 	sat(INT j) const {return (sresults[j]);}
const	INT*	at(INT j) const {return (lresults[j]);}
	INT	size(INT j) const {return (nok[j]);}
	bool	islong(INT j) const {return (j < hdr.max_sk);}
	float	cover_rate() const {return (cvrate);}
	void	get(float* dest, const INT& j, bool normalize = false);
};

class KmersFe : public Kmers {
	INT	xx = 0;
	INT	ww = 0;
	INT	sp = 0;
	long	ok = 0;
	INT	amb = 0;
	INT	llmt = 0;
public:
	void	count(INT c);
	void	reset() {sp = xx = ww = 0;}
	void	fromText(FILE* fd);
	void	fromSeq(Seq* sd);
template <typename file_t>
	void	fromNucDbs(file_t fd);
template <typename file_t>
	void	fromAaDbs(file_t fd);
template <typename file_t>
	void	readCount(file_t fd, const int& mol);
	void	readCount(const char* eij, const char* dbs = 0);
	void	readCount(int argc, const char** argv);
	KmersFe(INT a, INT xkk, INT nkk = 0, INT l = 0) :
	    Kmers(a, xkk, nkk), llmt(l)
	{reset();}
};

template <typename file_t>
int Kmers::read_binary(file_t fd, const char* fname)
{
	if (fd && fread(&hdr, sizeof(Header), 1, fd) != 1)
	    fatal(read_error, fname);
	if (hdr.maxk == 0) {	// old style
	    INT	harray[4];
	    memcpy(&harray[0], &hdr, sizeof(Header));
	    hdr.snwords = 0;
	    hdr.letter = harray[1];
	    hdr.maxk = harray[2];
	    hdr.mink = harray[3];
	    hdr.max_sk = 0;
	}
	supw = ipower(hdr.letter, hdr.maxk);
	make_nok();
	if (hdr.max_sk) {
	    scounts = new long[hdr.snwords];
	    sresults = new long*[hdr.max_sk + 1];
	    if (fd && 
		fread(scounts, sizeof(long), hdr.snwords, fd) != hdr.snwords)
		fatal(read_error, fname);
	}
	lcounts = new INT[hdr.lnwords];
	lresults = new INT*[hdr.maxk - hdr.mink + 1];
	if (fd && fread(lcounts, sizeof(INT), hdr.lnwords, fd) != hdr.lnwords)
	    fatal(read_error, fname);
	for (INT j = 0; j < hdr.mink; ++j) {	// skip up to mink
	    lresults[j] = lcounts;
	    if (sresults) sresults[j] = scounts;
	}
	long*	sct = scounts;
	for (INT j = hdr.mink; j < hdr.max_sk; ++j) {
	    sresults[j] = sct;
	    sct += nok[j];
	}
	INT*	lct = lcounts;
	for (INT j = hdr.mink; j < hdr.maxk; ++j) {
	    lresults[j] = lct;
	    lct += nok[j];
	}
	return (lct - lcounts);
}

template <typename file_t>
int Kmers::read_text(file_t fd, const char* fname)
{
	char	str[MAXL];
	char*	ps;
	char*	qs;
	char	cc = 0;
	vclear(&hdr);
	hdr.mink = USHRT_MAX;
	long	maxcount[2] = {0L, 0L};
	while ((ps = fgets(str, MAXL, fd))) {
	    ++hdr.lnwords;
	    hdr.maxk = strlen(car(ps, qs, cc));
	    if (hdr.maxk < 3) {
		++hdr.snwords;
		long	ncount = atol(car(++qs));
		if (maxcount[hdr.maxk - 1] < ncount)
		    maxcount[hdr.maxk - 1] = ncount;
	    }
	    if (hdr.mink == USHRT_MAX) hdr.mink = hdr.maxk - 1;
	    if (hdr.maxk == 1) ++hdr.letter;
	}
	rewind(fd);
	if (maxcount[1] > UINT_MAX) hdr.max_sk = 2;
	else if (maxcount[0] > UINT_MAX) hdr.max_sk = 1;
	else	hdr.max_sk = 0;
	lcounts = new INT[hdr.lnwords];
	lresults = new INT*[hdr.maxk + 1];
	if (hdr.max_sk) {
	    if (hdr.max_sk == 1) hdr.snwords = hdr.letter;
	    scounts = new long[hdr.snwords];
	    sresults = new long*[hdr.max_sk + 1];
	}
	INT*	lct = lcounts;
	long*	sct = scounts;
	INT	j = 0;
	for ( ; j < hdr.mink; ++j) {
	    lresults[j] = lct;
	    if (hdr.max_sk) sresults[j] = sct;
	}
	int	pk = 0;
	while ((ps = fgets(str, MAXL, fd))) {
	    int	k = strlen(car(ps, qs, cc));
	    long	frq = atol(car(++qs));
	    if (pk != k) {
		if (k <= hdr.max_sk)
		    sresults[j] = sct;
		lresults[j++] = lct;
		pk = k;
	    }
	    *lct++ = (INT) frq;
	    if (k <= hdr.max_sk) *sct++ = frq;
	}
	lresults[j] = lct;
	if (sresults) sresults[hdr.max_sk] = sct;
	if (!hdr.letter) {
	    if (hdr.maxk - hdr.mink > 2) {
		hdr.letter = (lresults[hdr.mink + 2] - lresults[hdr.mink + 1]) /
	    	(lresults[hdr.mink + 1] - lresults[hdr.mink]);
	    } else {
		INT	n = lct - lcounts;
		double r = pow((double) n, 1. / ((double) hdr.maxk));
		hdr.letter = INT(r + 1.e-5);
	    }
	}
	return (lct - lcounts);
}

template <typename file_t>
void Kmers::write_text(file_t fd)
{
	INT	n = ipower(hdr.letter, hdr.mink);
	char	str[MAXL];
	for (INT j = hdr.mink; j < hdr.maxk; ++j) {
	    n *= hdr.letter;
	    for (INT i = 0; i <  n; ++i) {
		INT	c = i;
		for (INT m = n; m /= hdr.letter; ) {
	 	    fputc(decoder[c / m], fd);
		    c %= m;
		}
		if (j < hdr.max_sk) snprintf(str, MAXL, "\t%7ld\n", sresults[j][i]);
		else	snprintf(str, MAXL, "\t%7d\n", lresults[j][i]);
		fputs(str, fd);
	    }
	}
}

template <typename file_t>
void KmersFe::fromNucDbs(file_t fd)
{
			/*     - A C M G R S V T W Y H K D B N */
static	int	nc16to4[]  = {15,0,1,4,2,4,4,4,3,4,4,4,4,4,4,4};
	int	parity = 0;
	int	c;

	while (true) {
	    int	a;
	    if (!parity) {
		c = fgetc(fd);
		if (c == EOF) break;
		a = (c >> 4) & 15;
	    } else {
		a = c & 15;
	    }
	    if ((a = nc16to4[a]) <= 4 && ++sp > llmt) {
		count(a);
		parity = 1 - parity;
	    } else {
		reset();
		parity = 0;
	    }
	}
}

template <typename file_t>
void KmersFe::fromAaDbs(file_t fd)
{
	while (true) {
	    int	c = fgetc(fd);
	    if (c == EOF) break;
	    if (c == SEQ_DELIM) reset();
	    else if (++sp > llmt)
		count((c >= ALA && c <= VAL)? c - ALA: 20);
	}
}

template <typename file_t>
void KmersFe::readCount(file_t fd, const int& mol)
{
	if (mol == DNA)	fromNucDbs(fd);
	else		fromAaDbs(fd);
}

