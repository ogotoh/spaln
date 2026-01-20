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

#include "kmers.h"

static	INT	right_m = 0;
static	bool	no_amb = false;
static	bool	homopoly = false;
static	const	char*	catalog;
static	bool	gzip = false;

void	Kmers::make_nok()
{
	if (!nok) nok = new INT[hdr.maxk];
	if (hdr.mink) vclear(nok, hdr.mink);
	INT	n = 1;
	for (INT j = hdr.mink; j < hdr.maxk; ++j)
	    nok[j] = (n *= hdr.letter);
}

Kmers::Kmers(INT a, INT xkk, INT nkk) : 
	supw(ipower(a, xkk)), molc(a == 4? DNA: PROTEIN),
	decoder(molc == DNA? Nucl: Amin)
{
	hdr.letter = a;
	hdr.maxk = xkk;
	hdr.mink = nkk;
	INT	n = ipower(hdr.letter, hdr.mink);
	hdr.lnwords = homopoly? hdr.letter * hdr.maxk:
		hdr.letter * (supw - n) / (hdr.letter - 1);
	lcounts = new INT[hdr.lnwords];
	lresults = new INT*[hdr.maxk + 1];
	make_nok();
	for (INT j = 0; j < hdr.mink; ++j)
	    lresults[j] = lcounts;
	INT	m = 0;
	for (INT j = hdr.mink; j < hdr.maxk; ++j) {
	    lresults[j] = lcounts + m;
	    m += (homopoly? hdr.letter: nok[j]);
	}
	lresults[hdr.maxk] = lcounts + m;
	vclear(lcounts, hdr.lnwords);

	hdr.snwords = nok[0] + nok[1];
	if (hdr.snwords) {
	    m = 0;
	    hdr.max_sk = 2;	// tentative
	    scounts = new long[hdr.snwords];
	    sresults = new long*[3];
	    for (INT j = hdr.mink; j < hdr.max_sk; ++j) {
		sresults[j] = scounts + m;
		m += nok[j];
	    }
	    sresults[hdr.maxk] = scounts + m;
	    vclear(scounts, hdr.snwords);
	}
}

Kmers::Kmers(const char* fname)
{
	ReadFile	fp(fname);
	int	nread = 0;
	if (fp.gzfd) {
	    if (fp.dtype == 1) nread = read_text(fp.gzfd, fname); else
	    if (fp.dtype == 2) nread = read_binary(fp.gzfd, fname);
	} else if (fp.fd) {
	    if (fp.dtype == 1) nread = read_text(fp.fd, fname); else
	    if (fp.dtype == 2) nread = read_binary(fp.fd, fname);
	} else
	    fatal(not_found, fname);
	if (!nread) fatal(read_error, fname);
	molc = hdr.letter == 4? DNA: PROTEIN;
	decoder = molc == DNA? Nucl: Amin;
	make_nok();
}

void KmersFe::count(INT c)
{
	INT	w = 0;
	if (c >= hdr.letter) {	// ambcode
	    ++amb;
	    xx = 0;
	} else if (homopoly) {
	    ++ok;
	    w = c;
	    if (xx == 0 || w == ww) {
		if (xx < INT(hdr.maxk - 1)) ++xx;
	    } else xx = 0;
	} else {
	    ++ok;		// good hdr.letter
	    if (xx < hdr.maxk) ++xx;
	    w = (hdr.letter * ww + c) % supw;
	}
	ww = w;
	if (homopoly) {
	    for (INT j = 0; j <= xx; ++j)
		lresults[j][w]++;
	} else if (xx > hdr.mink) {
	    for (INT j = hdr.mink; j < xx; ++j) {
		if (j < hdr.max_sk) sresults[j][w % nok[j]]++;
		lresults[j][w % nok[j]]++;
	    }
	}
}

void KmersFe::fromText(FILE* fd)
{
static	CHAR	ntconv[26] =
	{0,4,1,4,7,7,2,4,7,7,4,7,4,4,7,4,7,4,4,3,3,4,4,7,4,7};
static	CHAR	aaconv[26];
	CHAR*	encoder = ntconv;
	if (molc == PROTEIN) {
	    for (int i = 0; i < 26; ++i) aaconv[i] = aacode[i] - ALA;
	    encoder = aaconv;
	}
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
		if (c <= hdr.letter && ++sp > llmt)
		    count(c);
	    }
	}
}

void KmersFe::fromSeq(Seq* sd)
{
	reset();
	CHAR*	ps = sd->at(llmt);
	CHAR*	ts = sd->at(sd->right - right_m);
	CHAR*	encode = sd->isprotein()? aaredctab: ncredctab;

	for ( ; ps < ts; ++ps) {
	    INT	c = encode[*ps];
	    if (c <= hdr.letter) count(c);
	}
}

void KmersFe::readCount(int argc, const char** argv)
{
	SeqServer	svr(argc, argv, IM_SNGL, catalog, molc);
	Seq	sd(1);
	sd.setmolc(molc);
	InSt	ist;
	while ((ist = svr.nextseq(&sd)) != IS_END) 
	    if (ist == IS_OK) fromSeq(&sd);
	prompt("No. res = %ld, amb = %d\n", ok, amb);
}

void KmersFe::readCount(const char* eij, const char* dbs)
{
	EiJuncSeq	eijseq(INTRON, eij, dbs);
	do {
	    Seq*	sd = eijseq.nextseq();
	    if (!sd || (no_amb && sd->inex.ambs)) continue;
	    fromSeq(sd);
	} while (eijseq.goahead());
	prompt("No. res = %ld, amb = %d\n", ok, amb);
}

void Kmers::get(float* dest, const INT& j, bool normalize)
{				// get (k=j+1)-mer
	float*	dst = dest;
	float	total = 0;
	if (j < hdr.max_sk) {
	    long*	cnt = sresults[j];
	    long*	tcn = cnt + nok[j];
	    while (cnt < tcn) total += (*dst++ = (float) *cnt++);
	    if (!normalize) return;
	    while (dest < dst) *dest++ /= total; 
	} else {
	    INT*	cnt = lresults[j];
	    INT*	tcn = cnt + nok[j];
	    while (cnt < tcn) total += (*dst++ = (float) *cnt++);
	    if (!normalize) return;
	    while (dest < dst) *dest++ /= total; 
	}
	cvrate = total / nok[j];
}

void Kmers::outputCount(const char* outf, const int& text)
{
	if (homopoly) {
	    FILE*	afd = stdout;
	    if (outf) {
		afd = fopen(outf, "w");
		if (!afd) fatal(no_file, outf);
	    }
	    for (INT j = hdr.mink; j < hdr.maxk; ++j) {
		for (INT i = 0; i < hdr.letter; ++i) {
		    for (INT m = 0; m <= j; ++m) fputc(decoder[i], afd);
		    fprintf(afd, "\t%7d\n", lresults[j][i]);
		}
	    }
	    if (outf) fclose(afd);
	    return;
	}
	WriteFile	fp(outf, text, gzip);
	if (fp.gzfd) {
	    if (text)	write_text(fp.gzfd);
	    else	write_binary(fp.gzfd, outf);
	} else if (fp.fd) {
	    if (text)	write_text(fp.fd);
	    else	write_binary(fp.fd, outf);
	} else
	    fatal(no_file, outf);
}

#ifdef MAIN

static	void	Usage();
static	INT	left_m = 0;
static	bool	convert = false;
static	bool	binary = false;

static void Usage()
{
	fputs("Description:\n", stderr);
	fputs("\tCreate k-mer table from DNA/protein sequences\n", stderr);
	fputs("\toptionally with exon-intron junction information (.eij)\n", stderr);
	fputs("\t\tor\n", stderr);
	fputs("\tConvert between text and binary forms\n", stderr);
	fputs("Usage:\n", stderr);
	fputs("\tkmers -K[D|P] [-options] [-d dbs | fasta]\n", stderr);
	fputs("\tkmers -K[D|P] [-options] -d dbs -e X.eij\n", stderr);
	fputs("\tkmers -g -b X.wdfq -c X.wdfq (text -> binary)\n", stderr);
	fputs("\tkmers [-c X.wdfq] X.wdfq[.dat|.dgz] (binary -> text)\n", stderr);
	fputs("Options:\n", stderr);
	fputs("\t-b Binary_out (\".dat\" may be added)\n", stderr);
	fputs("\t-f [OutFile] (Convert binary <-> text)\n", stderr);
	fputs("\t-g: gzipped output (\".gz\" may be added)\n", stderr);
	fputs("\t  if -b option is also set \".dgz\" may be added\n", stderr);
	fputs("\t-h (print this message)\n", stderr);
	fputs("\t-H (homo oligomer only)\n", stderr);
	fputs("\t-K[P|D] (protein|DNA)\n", stderr);
	fputs("\t-l LeftMargin within intron\n", stderr);
	fputs("\t-r RightMargin within intron\n", stderr);
	fputs("\t-w WordLen (from 1 to w)\n", stderr);
	fputs("\t-W WordLen (only w)\n", stderr);
	exit (1);
}

int main(int argc, const char** argv)
{
	int	wd = -1;
	int	molc = DNA;
	bool	uptow = true;
const	char*	dbs = 0;
const	char*	outf = 0;
const	char*	eij = 0;
	char	str[LINE_MAX] = {'\0'};

	while (--argc > 0 && (++argv)[0][0] == '-') { 
	    const	char*	pn = 0;
	    const	char*	cmd = argv[0] + 1;
	    if (!*cmd) break;	/* stdin */
	    switch (*cmd) {
		case '?': case 'h': Usage();
		case 'b': binary = true;
		    if ((pn = getarg(argc, argv)))
			outf = pn;
		    break;
		case 'd':
		    if ((pn = getarg(argc, argv)))
			dbs = pn;
		    break;
		case 'e':
		    if ((pn = getarg(argc, argv)))
			eij = pn;
		    break;
		case 'f': convert = true;
		    if (argc > 2 && (pn = getarg(argc, argv)))
			outf = pn;
		    break;
		case 'g': gzip = true; break;
		case 'i':
		    if ((pn = getarg(argc, argv)))
			catalog = *pn == ':'? pn + 1: pn;
		    break;
		case 'l':
		    if ((pn = getarg(argc, argv)))
			left_m = atoi(pn);
		    break;
		case 'o':
		    if ((pn = getarg(argc, argv)))
			outf = pn;
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
	if (convert && argc) {	// convert binary <-> text
	    Kmers	kmers(argv[0]);
	    if (binary && !outf) outf = argv[0];
	    kmers.outputCount(outf, binary? 0: 2);
	    return (0);
	}
	if (wd < 0) wd = molc == DNA? 8: 5;
	if (wd < 1) fatal("Wordlen must be positive!");
	char*	strdbs = 0;
	if (dbs) {
	    strdbs = path2dbf(str, dbs, ".seq");
	    if (!strdbs) fatal(not_found, str);
	}
	KmersFe	kfe(molc == DNA? 4: 20, (INT) wd, uptow? 0: wd - 1, left_m);
	if (eij) kfe.readCount(eij, dbs);
	else if (strdbs) {
	    ReadFile	fp(strdbs, gz_ext);
	    if (fp.gzfd) kfe.readCount(fp.gzfd, molc);
	    else if (fp.fd) kfe.readCount(fp.fd, molc);
	    else	fatal(not_found, strdbs);
	} else if (argc == 0) 
	    kfe.fromText(stdin);
	else kfe.readCount(argc, argv);
	kfe.outputCount(outf, binary? 0: 2);
	return(0);
}

#endif	// MAIN
