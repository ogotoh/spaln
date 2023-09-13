/**************************************************************************
*
*	Find Patterns 
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
#include "resite.h"
#include "pattern.h"

#define equaln(a, b) (((a - base) & (b - base)) && ((a) != N || (b) == N))
#define equalp(a, b) ((a) == (b) || (b) == AMB)
#define equal(a, b) (Nuc? equaln(a, b): equalp(a, b))
#define Equaln(a, b) ((a - base) & (b - base))
#define Equalp(a, b) ((a) == (b) || (a) == AMB || (b) == AMB)
#define Equal(a, b) (Nuc? Equaln(a, b): Equalp(a, b))

static	int	match(int res, CHAR* pt, int base);
static	void	bmdelta1(int* dlt, CHAR* str, int slen, SEQ_CODE* defcode);
static	void	bmdelta2(int* dlt, CHAR* pat, int slen, int base);
static	int	simplepat(int* loc, int maxv, Seq* sd, char* pat, SEQ_CODE* defcode);
static	int	forematch(CHAR* sw, VARPAT* pat, int j, int base);
static	void	delta1(int* dlt, VARPAT* pat, int mlen, SEQ_CODE* defcode);
static	CHAR*	aagroup(CHAR* pt, char* mb);
static	VARPAT*	genpat(char* ps);
static	int	complexpat(int* loc, int maxv, Seq* sd, char* src, SEQ_CODE* defcode);
static	int	icmpf(int* a, int* b);
static	SEQ_CODE*	aas_code;

static int	Nuc = 0;
static char	lenerr[] = "Too long pattern !\n";
static char	sizerr[] = "Too complicated pattern !\n";
static char	adh[] = "+|";
static char	sep[] = " \t;+|";

static int match(int res, CHAR* pt, int base)
{
	int	unm = (res ^ *pt) & HITEI;

	res &= ~HITEI;
	while (int ptn = (*pt++ & ~HITEI))
	    if (equal(res, ptn)) return (!unm);
	return (unm);
}

static void bmdelta1(int* dlt, CHAR* str, int slen, SEQ_CODE* defcode)
{
	int	base = defcode->base_code - 1;

	for (int i = 0; i < defcode->max_code; i++)
	    dlt[i] = slen;
	for ( ; --slen; str++)
	    for (int i = base + 1; i < defcode->max_code; ++i)
		if (Equal(i, *str))
		    dlt[i] = slen;
}

static void bmdelta2(int* dlt, CHAR* pat, int slen, int base)
{
	for (int k = slen; k > 0; k--) {
	    int	i = slen;
	    for (int j = slen - k; j > 0; ) {
		--i; --j;
		if (!Equal(pat[i], pat[j])) {
		    dlt[i] = k;
		    i = 0;
		    break;
		}
	    }
	    while (--i >= 0) dlt[i] = k;
	}
}

static int simplepat(int* loc, int maxv, Seq* sd, char* pat, SEQ_CODE* defcode)
{
	int	base = defcode->base_code - 1;
	int*	pl = loc;
	int 	dlt2[MAXL];
	int 	patlen = strlen(pat);
	CHAR*	sl = sd->at(sd->left + patlen - 1);
	CHAR*	sr = sd->at(sd->right);

	if (sd->byte > 1) {
	    fputs("Can't find pattern in multiple sequence\n", stderr);
	    return (0);
	}
	int* 	dlt1 = new int[defcode->max_code + 1];
	if (*pat == '<') {
	    pat++;
	    patlen--;
	    sr = sl--;
	}
	if (CHAR* sw = (CHAR*) strchr(pat, '>')) {
	    *sw = '\0';
	    patlen--;
	    sl = sr - 1;
	}
	if (*pat >= defcode->max_code) tosqcode((CHAR*) pat, defcode);
	bmdelta1(dlt1, (CHAR*) pat, patlen, defcode);
	bmdelta2(dlt2, (CHAR*) pat, patlen, base);
/*	for ( ; sl < sr; sl += max(dlt1[*sl], dlt2[j]))	*/
	for ( ; sl < sr; sl += max(dlt1[*sl], int(1))) {
	    CHAR* sw = sl;
	    for (int j = patlen - 1; equal(*sw, pat[j]); j--, sw--)
		if (!j) {
		    *pl++ = sd->index(sw);
		    if (--maxv < 0) {
			fputs("Too large sites !\n", stderr);
			goto loopend;
		    }
		    break;
		}
	}
loopend:
	delete[] dlt1;
	return (pl - loc);
}

static CHAR* sleft;

static int forematch(CHAR* sw, VARPAT* pat, int j, int base)
{
	while (--j >= 0) {
	    int i = 0;
	    for ( ; i < pat[j].from; i++)
		if (!match(*sw--, pat[j].pat, base)) return (j);
	    for ( ; i < pat[j].till; i++) {
		if (forematch(sw, pat, j, base) < 0) return (-1);
		if (!match(*sw--, pat[j].pat, base)) return (j);
	    }
	}
	sleft = sw + 1;
	return (j);
}

static void delta1(int* dlt, VARPAT* pat, int mlen, SEQ_CODE* defcode)
{
	int	base = defcode->base_code - 1;

	for (int i = 0; i < defcode->max_code; i++)
	    dlt[i] = mlen;
	for ( ; pat->pat; pat++) {
	    mlen -= pat->from;
	    for (int i = 0; i < defcode->max_code; i++)
		if (match(i, pat->pat, base))
		    dlt[i] = mlen;
	}
}

static CHAR	pattern[MAXL+1];
static VARPAT	expat[MAXPATLEN+1];
static int  inflag = 0;
static int  outflag = 0;

static CHAR* aagroup(CHAR* pt, char* mb)
{
	while (*mb) {
	    if (pt - pattern >= MAXL) {
		fputs(sizerr, stderr);
		return (NULL);
	    }
	    *pt++ = (outflag)?
		(en_code(*mb++, aas_code) | HITEI): en_code(*mb++, aas_code);
	}
	return (pt);
}

static VARPAT*	genpat(char* ps)
{
static char hpho[] = "cfilmvwy";
static char hphi[] = "dehknqr";
static char plus[] = "hkr";
static char minus[] = "deqn";
static char alif[] = "ilvm";
static char arom[] = "fyw";
static char smll[] = "agpst";

	CHAR*	pat = pattern;
	VARPAT*	varpat = expat;

	varpat->pat = pat;
	varpat->from = varpat->till = 1;
	for ( ; *ps; ps++) {
	    switch (*ps) {
		case '[':   inflag = 1; break;
		case ']':   inflag = 0; break;
		case '{':   outflag = 1; break;
		case '}':   outflag = 0; break;
		case '-':
		    *pat++ = 0;
		    if (++varpat - expat >= MAXPATLEN) {
			fputs(lenerr, stderr);
			goto endloop;
		    }
		    varpat->pat = pat; 
		    varpat->from = varpat->till = 1;
		    break;
		case '(':
		    varpat->from = varpat->till = atoi(++ps);
		    while (*++ps != ')')
			if (*ps == ',')
 			    varpat->till = atoi(++ps);
		    break;
		case '<':
		case '>': break;
		case 'o':
		    if (!(pat = aagroup(pat, hpho))) goto endloop;
		    break;
		case 'j':
		    if (!(pat = aagroup(pat, hphi))) goto endloop;
		    break;
		case '+':
		    if (!(pat = aagroup(pat, plus))) goto endloop;
		    break;
		case '_':
		    if (!(pat = aagroup(pat, minus))) goto endloop;
		    break;
		case '@':
		    if (!(pat = aagroup(pat, alif))) goto endloop;
		    break;
		case '$':
		    if (!(pat = aagroup(pat, arom))) goto endloop;
		    break;
		case '.':
		    if (!ps[1] || isspace(ps[1])) goto endloop;
		    if (!(pat = aagroup(pat, smll))) goto endloop;
		    break;
		default:
		    if (pat - pattern >= MAXL) {
			fputs(sizerr, stderr);
			goto endloop;
		    }
		    *pat = en_code(*ps, aas_code);
		    if (*pat > 0) {
			if (outflag) *pat |= HITEI;
			pat++;
		    }
		    break;
	    }
	}
endloop:
	*pat = 0;
	(++varpat)->pat = NULL;
	return (expat);
}

static int complexpat(int* loc, int maxv, Seq* sd, char* src, SEQ_CODE* defcode)
{
	if (sd->byte > 1) {
	    prompt("Can't find pattern in multiple sequence\n");
	    return (0);
	}

	aas_code = sd->code;
	VARPAT*	pat = genpat(src);
	VARPAT*	pt = pat;
	int	minlen = 0;
	int	asblen = 0;
	while (pt->pat) {
	    asblen += pt->till - pt->from;
	    minlen += (pt++)->from;
	}
	int 	patlen = pt - pat;
	CHAR*	sl = sd->at(sd->left + minlen);
	CHAR*	sr = sd->at(sd->right);
	if (strchr(src, '<')) sr = sl + asblen;
	if (strchr(src, '>')) sl = sr - asblen;
	--sl;
	int* 	dlt1 = new int[defcode->max_code];
	delta1(dlt1, pat, minlen, defcode);
/*
	delta2(dlt2, pat, patlen);
	for ( ; sl < sr; sl += max(dlt1[*sl], dlt2[j])) {
*/
	int	base = defcode->base_code - 1;
	int*	pl = loc;
	for ( ; sl < sr; sl += max(dlt1[*sl], int(1))) {
	    int	j = forematch(sl, pat, patlen, base);
	    if (j < 0) {
		j = 0;
		*pl++ = sd->index(sleft);
		if (--maxv < 0) {
		    fputs("Too large sites !\n", stderr);
		    goto loopend;
		}
	    }
	}
loopend:
	delete[] dlt1;
	return (pl - loc);
}

char* parspat(char* pat)
{
	char*	parsed = pat;
	char*	wk = pat;

	while (strchr(sep, *pat)) pat++;    /* Skip delim   */
	while (*pat) {
	    while (*pat && !strchr(sep, *pat))
		*wk++ = *pat++;
	    int	sw = 0;
	    for ( ; *pat && strchr(sep, *pat); pat++)
		if (strchr(adh, *pat)) sw = 1;
	    if (*pat) *wk++ = sw? '+': ' ';
	}
	*wk ++ = '\0';
	return (parsed);
}

static int icmpf(int* a, int* b)
{
	return (*a - *b);
}

int findseq(int loc[], int maxloc, Seq* sd, char* pat)
{
static	char	cmppat[] = "-[]{}";
static	char*	buf = 0;
	char*	pb = buf = strrealloc(buf, pat);
	int 	num = 0;
	bool	cmplx = false;
	SEQ_CODE* defcode = sd->code;

	while (*pb)	    /* test whether pattern is complex or not	*/
	    if (strchr(cmppat, *pb++)) {
		cmplx = true;
		break;
	    }

	Nuc = sd->isdrna();
	for (pb = strtok(buf, adh); pb; pb = strtok(NULL, adh))
	    if (cmplx)
		num += complexpat(loc + num, maxloc - num, sd, pb, defcode);
	    else
		num += simplepat(loc + num, maxloc - num, sd, pb, defcode);
	qsort((UPTR) loc, (INT) num, sizeof(int), (CMPF) icmpf);
	return (num);
}

void putloc(FILE* fd, Seq* sd, int* pl, int num)
{
	int		i = 0;

	for ( ; i < num; ++pl) {
		if (!(i % 10)) putc('\t', fd);
		fprintf(fd, "%5d ", sd->SiteNo(*pl));
		if (!(++i % 10)) putc('\n', fd);
	}
	if (i % 10) putc('\n', fd);
}

void findpat(FILE* fd, Seq* sd, char* pat)
{
	int	num;
	char*	pt = parspat(pat);
	char*	pb;
	char*	ps;
	int	loc[MAXLOC];
	Seq*	single = NULL;
	int	grp[2];
	int	i;
	int	rsv;

	grp[1] = -1;
	if (sd->many > 0) {sd->fphseq(fd); fputc('\n', fd);}
	for (i = 0; i < sd->many; ++i) {
	    if (sd->many == 1) single = sd;
	    else {
		grp[0] = i;
		single = sd->extseq(single, grp, CPY_ALL);
	    }
	    for (ps = pt; *(pb = ps); ++ps) {
		while (*ps && !isspace(*ps)) ps++;
		rsv = *ps;
		*ps = '\0';
		num = findseq(loc, MAXLOC, single, pb);
		if (num > 0) {
		    if (sd->many == 1) single->fphseq(fd);
		    else {fprintf(fd, "%-10s %3d\t", (*sd->sname)[i], i+1);}
		    fprintf(fd, "  %s  %d\n", pb, num);
		    putloc(fd, single, loc, num);
		}
		if (!(*ps = rsv)) break;
	    }
	}
	if (sd->many > 1) delete single;
}

