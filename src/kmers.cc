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

#define	SKIP	0xff

static	CHAR	ntconv[26] =
	{0,4,1,4,7,7,2,4,7,7,4,7,4,4,7,4,7,4,4,3,3,4,4,7,4,7};
static	CHAR	aaconv[26];
static	INT	left_m = 0;
static	INT	right_m = 0;
static	bool	no_amb = false;
static	bool	homopoly = false;
static	bool	uptow = true;
static	const	char*	catalog;

static	void	Usage();

static void Usage()
{
	fputs("Usage:\n", stderr);
	fputs("\tkmers -K[D|P] [-options] [-d dbs.seq | fasta]\n", stderr);
	fputs("\tkmers -d dbs -K[D|P] [-options] -e xxx.eij\n", stderr);
	fputs("Options:\n", stderr);
	fputs("\t-h (print this message)\n", stderr);
	fputs("\t-H (homo oligomer only)\n", stderr);
	fputs("\t-K[P|D] (protein|DNA)\n", stderr);
	fputs("\t-l LeftMargin\n", stderr);
	fputs("\t-r RightMarging\n", stderr);
	fputs("\t-w MaxWordLen (upto w)\n", stderr);
	fputs("\t-W MaxWordLen (only w)\n", stderr);
	exit (1);
}

class Kmers {
	INT	letter;
	INT	width;
	INT	mask;
	INT	words;
	INT*	counts;
	INT**	results;
	INT	xx;
	INT	ww;
	INT	sp;
	INT	ok, amb;
public:
	Kmers(INT a, INT k);
	~Kmers() {
	    delete[] counts; delete[] results;
	}
	void	fromText(FILE* fd, CHAR* encoder);
	void	fromSeq(Seq* sd);
template <typename file_t>
	void	fromNucDbs(file_t fd);
template <typename file_t>
	void	fromAaDbs(file_t fd);
template <typename file_t>
	void	readCount(file_t fd, const int& molc);
	void	readCount(const char* eij, const char* dbs = 0);
	void	readCount(int argc, const char** argv, const int& molc);
	void	outputCount(char* decode);
	void	count(INT c);
	void	reset();
};

Kmers::Kmers(INT a, INT kk) : letter(a), width(kk)
{
	mask = ipower(letter, width);
	INT	words = homopoly? letter * width:
		(mask * letter - 1) / (letter - 1);
	counts = new INT[words];
	results = new INT*[width];
	INT	n = 1;
	INT	m = 0;
	for (INT k = 0; k < width; ++k) {
	    results[k] = counts + m;
	    m += (homopoly? letter: (n *= letter));
	}
	vclear(counts, words);
	ok = amb = 0;
	reset();
}

void Kmers::reset()
{
	sp = xx = ww = 0;
}

void Kmers::count(INT c)
{
	INT	w = 0;
	if (c >= letter) {	// ambcode
	    ++amb;
	    xx = 0;
	} else if (homopoly) {
	    ++ok;
	    w = c;
	    if (xx == 0 || w == ww) {
		if (xx < width - 1) ++xx;
	    } else xx = 0;
	} else {
	    ++ok;		// good letter
	    if (xx < width) ++xx;
	    w = (letter * ww + c) % mask;
	}
	ww = w;
	if (homopoly) {
	    for (INT k = 0; k <= xx; ++k)
		results[k][w]++;
	} else {
	    INT	n = letter;
	    for (INT k = 0; k < xx; ++k, n *= letter)
		results[k][w % n]++;
	}
}

void Kmers::fromText(FILE* fd, CHAR* encoder)
{
	char	str[MAXL];
	while (char* ps = fgets(str, MAXL, fd)) {
	    if (*str == '>') {
		reset();
		while (strlen(str) == MAXL - 1)
		    if (!(ps = fgets(str, MAXL, fd))) break;
		continue;
	    }
	    while (ps && *ps) {
		INT	c = *ps++;
		if (c == ';' || c == '#') break;// comments
		if (!isalpha(c)) continue;	// ignore
		c = encoder[tolower(c)-'a'];
		if (c <= letter && ++sp > left_m)
		    count(c);
	    }
	}
}

void Kmers::fromSeq(Seq* sd)
{
	reset();
	CHAR*	ps = sd->at(left_m);
	CHAR*	ts = sd->at(sd->right - right_m);
	CHAR*	encode = sd->isprotein()? aaredctab: ncredctab;

	for ( ; ps < ts; ++ps) {
	    INT	c = encode[*ps];
	    if (c <= letter) count(c);
	}
}


template <typename file_t>
void Kmers::fromNucDbs(file_t fd)
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
	    if ((a = nc16to4[a]) <= 4 && ++sp > left_m) {
		count(a);
		parity = 1 - parity;
	    } else {
		reset();
		parity = 0;
	    }
	}
}

template <typename file_t>
void Kmers::fromAaDbs(file_t fd)
{
	while (true) {
	    int	c = fgetc(fd);
	    if (c == EOF) break;
	    if (c == SEQ_DELIM) reset();
	    else if (++sp > left_m)
		count((c >= ALA && c <= VAL)? c - ALA: 20);
	}
}

template <typename file_t>
void Kmers::readCount(file_t fd, const int& molc)
{
	if (molc == DNA)	fromNucDbs(fd);
	else			fromAaDbs(fd);
}

void Kmers::readCount(int argc, const char** argv, const int& molc)
{
	SeqServer	svr(argc, argv, IM_SNGL, catalog, molc);
	Seq	sd(1);
	sd.setmolc(molc);
	InSt	ist;
	while ((ist = svr.nextseq(&sd)) != IS_END) 
	    if (ist == IS_OK) fromSeq(&sd);
}

void Kmers::readCount(const char* eij, const char* dbs)
{
	EiJuncSeq	eijseq(INTRON, eij, dbs);
	do {
	    Seq*	sd = eijseq.nextseq();
	    if (!sd || (no_amb && sd->inex.ambs)) continue;
	    fromSeq(sd);
	} while (eijseq.goahead());
}
	
void Kmers::outputCount(char* decoder)
{
	prompt("No. res = %d, amb = %d\n", ok, amb);
	INT	k = uptow? 0: (width -  1);
	if (homopoly) {
	    for ( ; k < width; ++k) {
		for (INT j = 0; j < letter; ++j) {
		    for (INT m = 0; m <= k; ++m)
			putchar(decoder[j]);
		    printf("\t%7d\n", results[k][j]);
		}
	    }
	} else {
	    INT	n = letter;
	    for ( ; k < width; ++k, n *= letter) {
		for (INT j = 0; j <  n; ++j) {
		    INT	c = j;
		    for (INT m = n; m /= letter; ) {
			putchar(decoder[c / m]);
			c %= m;
		    }
		    printf("\t%7d\n", results[k][j]);
		}
	    }
	}
}

int main(int argc, const char** argv)
{
	int	wd = -1;
	int	molc = DNA;
const	char*	dbs = 0;
const	char*	eij = 0;
	char	str[LINE_MAX] = {'\0'};

	while (--argc > 0 && (++argv)[0][0] == '-') { 
	    const	char*	pn = 0;
	    const	char*	cmd = argv[0] + 1;
	    if (!*cmd) break;	/* stdin */
	    switch (*cmd) {
		case '?': case 'h': Usage();
		case 'd':
		    if ((pn = getarg(argc, argv)))
			dbs = pn;
		    break;
		case 'e':
		    if ((pn = getarg(argc, argv)))
			eij = pn;
		    break;
		case 'i':
		    if ((pn = getarg(argc, argv)))
			catalog = *pn == ':'? pn + 1: pn;
		    break;
		case 'l':
		    if ((pn = getarg(argc, argv)))
			left_m = atoi(pn);
		    break;
		case 'r':
		    if ((pn = getarg(argc, argv)))
			right_m = atoi(pn);
		    break;
		case 'w':
		   if ((pn = getarg(argc, argv)))
			wd = atoi(pn);
		    break;
		case 'W':
		    if ((pn = getarg(argc, argv)))
			wd = atoi(pn);
		    uptow = false;
		    break;
		case 'K':
		    if ((pn = getarg(argc, argv)) &&
		    	(*pn=='a'||*pn=='A'||*pn=='p'||*pn=='P'))
			molc = PROTEIN;
		    break;
		case 'H': homopoly = true; break;
		default: break;
	    }
	}
	if (wd < 0) wd = molc == DNA? 8: 5;
	if (wd < 1) fatal("Wordlen must be positive!");
	if (!dbs && molc == PROTEIN)
	    for (int i = 0; i < 26; ++i) aaconv[i] = aacode[i] - ALA;
	char*	strdbs = dbs? path2dbf(str, dbs, ".seq"): 0;
	Kmers	kmers(molc == DNA? 4: 20, (INT) wd);
	if (eij) kmers.readCount(eij, dbs);
	else if (strdbs) {
	    if (is_gz(strdbs)) {
#if USE_ZLIB
		gzFile_s*	gzfd = gzopen(strdbs, "rb");
		if (!gzfd) fatal(not_found, strdbs);
		kmers.readCount(gzfd, molc);
#else
		fatal(gz_unsupport, strdbs);
#endif
	    } else {
		FILE*	fd = fopen(strdbs, "rb");
		if (!fd) fatal(not_found, strdbs);
		kmers.readCount(fd, molc);
	    }
	} else if (argc == 0) 
	    kmers.fromText(stdin, molc == DNA? ntconv: aaconv);
	else kmers.readCount(argc, argv, molc);
	kmers.outputCount(molc == DNA? Nucl: Amin);
	return(0);
}
