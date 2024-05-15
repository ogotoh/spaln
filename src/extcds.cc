/*****************************************************************************
*
*	Extract coding regions from GenBank flat file
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
#include "extcds.h"

#define DEBUG	0
#define KITSUI	0

/*
static	const	int	MINORF = 300;
static	const	int	MAXAAS = 20000;
static	const	int	ddelim = ((SEQ_DELIM << 4) + SEQ_DELIM);
*/

#define strmatch(str, pat) (!strncmp(str, pat, sizeof(pat) - 1))
#define Fwrite(dt, sz, it, fd) \
if (fwrite(dt, sz, it, fd) != (it)) fatal("Write error")

static	int	countcds(List* wexon);
static	Seq*	extract(Seq* src, List* wexon, int c);
static	void	onecds(List* wexon, char* ps);
static	bool	findsegno(List* wexon, char* ps);
static	void	preonecds(List* wexon, char* ps);
static	void	fillexon(List* wexon, List* fexon, char* ps, long ptr);
static	int	cmpf(List* x, List* y);
static	void	rm_overlap(List* fexon, int nn);
static	int	collectcds();
static	int	parbalance(char* ps);
static	void	joinedlist(List* fexon);
static	Seq*	extcds(Seq* sd);
extern	void	extractcds(Seq* sd);

static FILE*	fsrc;
static List	exons[CDSTYPES][EXONMAX];
static SEGM	segments[SEGMAX];
static int	segidx;
static int	segnum;
static const	char*	InpRec = "Not found!\n";

static int countcds(List* wexon)
{
	int	n = 0;

	for ( ; wexon->from != EOE; ++wexon) {
	    int	l = wexon->to - wexon->from;
	    n += abs(l);
	}
	return (n);
}

static Seq* extract(Seq* src, List* wexon, int c)
{
	int	nn = countcds(wexon);

	if (nn <= 0) return (0);
	Seq*	dst = new Seq(1, nn + 2);
	CHAR*	dd = dst->at(0);
	for ( ; wexon->from != EOE; ++wexon) {
	    memcpy(dd, src->at(wexon->from), wexon->to - wexon->from);
	    dd += wexon->to - wexon->from;
	}
	dst->postseq(dd);
	if (c) dst->comrev();
	return (dst);
}

static void onecds(List* wexon, char* ps)
{
	wexon->from = wexon->to = 
	wexon->texon.lpart = wexon->texon.rpart = 0;
	if (*ps == '<') {
	    wexon->texon.lpart = 1;
	    ++ps;
	}
	char*	qs = ps;
	for ( ; isdigit(*qs); ++qs)
	    if (!*qs) return;
	*qs++ = '\0';
	if (!*qs++) return;
	wexon->from = atoi(ps) - 1;
	if (*qs == '>') {
	    wexon->texon.rpart = 1;
	    ++qs;
	}
	for (ps = qs; isdigit(*qs); ++qs) ;
	*qs = '\0';
	if (ps == qs)	wexon->from = 0;
	else		wexon->to = atoi(ps);
}

static bool findsegno(List* wexon, char* ps)
{
	for (int i = 1; i <= segnum; ++i) {
	    if (strstr(segments[i].acc, ps)) {
		wexon->texon.segno = i;
		return (true);
	    }
	}
	return (false);
}

static void preonecds(List* wexon, char* ps)
{
	if (*ps == 'j' && strmatch(ps, "join(")) ps += 5;
	if (*ps == 'c' && strmatch(ps, "complement(")) {
	    wexon->texon.cmpl = 1;
	    ps += 11;
	}
	if (*ps == 'j' && strmatch(ps, "join(")) ps += 5;
	if (char* qs = strchr(ps, ':')) {
	    *qs++ = '\0';
	    if (segidx == 0)
		fprintf(stderr, "%-10s : %s -- Unexpected SEG\n",
			InpRec, ps);
	    else if (!findsegno(wexon, ps))
		fprintf(stderr, "%-10s : %s -- Can't find Seg No. %d\n",
			InpRec, ps, wexon->texon.segno);
	    ps = qs;
	}
	onecds(wexon, ps);
}

static void fillexon(List* wexon, List* fexon, char* ps, long ptr)
{
	static	EXTYPE zeroex = {0, 0, 0, 0, 0, 0};

	if (wexon->texon.join) return;
	ps += datacolumn;
	wexon->texon = zeroex;
	if (*ps == 'c' && strmatch(ps, "complement(")) {
	    wexon->texon.cmpl = 1;
	    ps += 11;
	}
	if (*ps == 'j' && strmatch(ps, "join(")) {
	    fexon->texon.join = 1;
	    fexon->texon.cmpl = wexon->texon.cmpl;
	    fexon->from = ptr;
	} else if ((*ps == 'o' && strmatch(ps, "order")) 
		|| (*ps == 'r' && strmatch(ps, "replace"))) {
	    fexon->texon.undef = 1;
	} else if (fexon->texon.join == 0) {
	    onecds(wexon, ps);
	}
}

static int cmpf(List* x, List* y)
{
	return (x->from - y->from);
}

static void rm_overlap(List* fexon, int nn)
{
	qsort((UPTR) fexon, (INT) nn, sizeof(List), (CMPF) cmpf);
	for (List* wexon = fexon + 1; wexon->from != EOE; ++wexon) {
	    if (wexon->from > fexon->to)
		*++fexon = *wexon;
	    else if (fexon->to < wexon->to)
		fexon->to = wexon->to;
	}
	(++fexon)->from = EOE;
}

static int feature = 0;

static int collectcds()
{
	char	str[MAXL];
	char*	ps;
	long	ptr = ftell(fsrc);
	char	stracc[LLINE];
	List*	wexon[CDSTYPES];
	int	lbound[CDSTYPES];
	int	rbound[CDSTYPES];
	int	cds_overlap[CDSTYPES];
	int	endtype = ENTRY;
	int	codetype;
	int	mrna = false;

	for (INT i = 0; i < CDSTYPES; ++i) {
	    wexon[i] = exons[i];
	    wexon[i]->texon.join = lbound[i] = rbound[i] = 0;
	    cds_overlap[i] = false;
	}
	for ( ; fgets_wocr(str, MAXL, fsrc); ptr = ftell(fsrc)) {
	  if (*str == 'O' && strmatch(str, "ORIGIN")) break;
	  switch (*str) {
	    case 'L':
		if (strmatch(str, "LOCUS")) {
		    wexon[MRNA]->from = 0;
		    sscanf(str, "%*s %*s %ld %*s %s", &wexon[MRNA]->to, stracc);
		    if (!strcmp(stracc, "ss-mRNA")) mrna = true;
		}
		break;
	    case 'S':
		if (strmatch(str, "SEGMENT"))
		    fprintf(stderr, "Sorry! not supported\n");
		break;
	    case 'A':
		if (strmatch(str, "ACCESSION"))
		    strcpy(stracc, str + sizeof("ACCESSION"));
		break;
	    case 'F':
		if (strmatch(str, "FEATURES")) feature = 1;
		break;
#if DEBUG
	    case 'D':
		fputs(str, stdout);
		break;
#endif
	    case ' ':
		if (!feature) continue;
		ps = str + typecolumn;
		codetype = UNDEF;
		switch (*ps) {
		  case 'e':
		    if (strmatch(ps, "exon")) codetype = EXON;
		    break;
		  case 'C':
		    if (strmatch(ps, "CDS")) codetype = CDS;
		    break;
		  case 'm':
		    if (strmatch(ps, "mRNA")) codetype = MRNA;
		    break;
		  case ' ':
		    ps = str + datacolumn;
		    if (*ps++ != '/') continue;
		    if (*ps == 't' && strmatch(ps, "translation")) {
			endtype = TRANS;
			goto eoent;
		    } else if (*ps == 'p' && strmatch(ps, "pseudo")) {
			endtype = NOCDS;
			goto eoent;
		    }
		    break;
		  default:
		    break;
		}
		if (codetype != UNDEF) {
		    fillexon(wexon[codetype], exons[codetype], str, ptr);
		    if (!exons[codetype]->texon.join) {
			if (wexon[codetype]->to > lbound[codetype] ||
			    wexon[codetype]->from < rbound[codetype])
				cds_overlap[codetype] = true;
			if (lbound[codetype] > wexon[codetype]->from)
			    lbound[codetype] = wexon[codetype]->from;
			if (rbound[codetype] < wexon[codetype]->to)
			    rbound[codetype] = wexon[codetype]->to;
		    }
		    ++wexon[codetype];
		}
		break;
	    default:	break;
	  }
	}
eoent:
	if (mrna && wexon[MRNA] == exons[MRNA]) ++wexon[MRNA];
	int	n = 0;
	for (INT i = 0; i < CDSTYPES; ++i) {
	    wexon[i]->from = EOE;
	    if (wexon[i] > exons[i]) {
		++n;
		if (!exons[i]->texon.join && cds_overlap[i])
		    rm_overlap(exons[i], wexon[i] - exons[i]);
	    }
	}
	return (endtype);
}

static 	char	jstr[JOINMAX];
static	char*	js;

static int parbalance(char* ps)
{
	int	n = 0;

	while (*ps) {
	    if (js - jstr < JOINMAX - 1 && !isspace(*ps)) *js++ = *ps;
	    switch (*ps++) {
		case '(': ++n; break;
		case ')': --n; break;
		default:  break;
	    }
	}
	return (n);
}

static void joinedlist(List* fexon)
{
	char	str[MAXL];
	char*	ps;
	long	ptr = ftell(fsrc);
	List*	wexon = fexon;
	int	par = 0;

	fseek(fsrc, fexon->from, SEEK_SET);
	js = jstr;
	do {
	    if (!fgets_wocr(str, MAXL, fsrc))
		fatal("%-10s Unexpected EOF", InpRec);
	    ps =  str + datacolumn;
	    par += parbalance(ps);
	} while (par > 0);
	*js = 0;
	if (js - jstr >= JOINMAX)
	    fprintf(stderr, "%d > Buffer for joined segs\n", JOINMAX);
	for (ps = strtok(jstr, ","); ps; ps = strtok(0, ","))
	    preonecds(wexon++, ps);
	wexon->from = EOE;
	fseek(fsrc, ptr, SEEK_SET);
}

static Seq* extcds(Seq* sd)
{
	INT	i = 0;

	for ( ; i < CDSTYPES; ++i)	
	    if (exons[i]->from != EOE) break;
	if (i == CDSTYPES) return(0);	/* No CDS */
					/* Break on varid */
	int	c = exons[i]->texon.cmpl;
	if (exons[i]->texon.join) joinedlist(exons[i]);
	return extract(sd, exons[i], c);
}

void extractcds(Seq* sd)
{
	int	endtype;
	int	n = 0;
	char	str[MAXL];
	FILE*	fd = qout("");

	fsrc = fopen(sd->spath, "r");
	if (!fsrc) fatal("Bad source %s!", sd->spath);
	feature = 0;
	for (;;) {
	    endtype = collectcds();
	    if (endtype == NOCDS) continue;
	    if (n++ > 0 && endtype != TRANS) break;
	    Seq*	scds = extcds(sd);
	    if (scds) {
		snprintf(str, MAXL, "%s_%d", sd->sqname(), n);
		scds->sname = new Strlist(&str[0], "");
		scds->typeseq(fd);
		delete scds;
	    }
	}
	fclose(fsrc);
	qclose();
}
