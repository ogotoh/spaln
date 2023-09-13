/*****************************************************************************
*
*	Restriction enzyme cleavage sites
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

#define	ALL_ENZ	0
#define UNIQENZ 1

static	FILE*	openresite(int enzset);
static	int	cmpf(char* key, char* str);
static	void	cmplseq(CHAR* dst, CHAR* src);
static	RESITE*	getresite(char* buf);
static	int	icmpf(int* a, int* b);

static	RESITE	resite;
static	FILE*	resfd;
static	int	width;
static	long	nrec;
static	const	char*	refname[2] = {"renzyme", "resuniq"};

static FILE* openresite(int enzset)
{
	if (!resfd) {
	    resfd = ftable.fopen(refname[enzset], "r");
	    if (!resfd) return (0);
	    int	c;
	    for (width = 0; (c = getc(resfd)) != '\n'; width++) 
		if (c == EOF) return (0);
	    fseek(resfd, 0L, 2);
	    nrec = ftell(resfd) / ++width;
	}
	rewind(resfd);
	return (resfd);
}

static int cmpf(char* key, char* str)
{
	while (*key && !isspace(*key) && *str && !isspace(*str)) {
	    int	c = *key++ - toupper(*str++);
	    if (c) return (c);
	}
	return (0);
}

static void cmplseq(CHAR* dst, CHAR* src)
{
	extern CHAR complcod[];

	dst += strlen((char*) src);
	*dst-- = '\0';
	while (*src) *dst-- = complcod[*src++];
}

static RESITE* getresite(char* buf)
{
	resite.rct = 0;
	sscanf(buf, "%s %s %d %d", 
	    resite.name, resite.chr, &resite.cut, &resite.rct);
	tosqcode((CHAR*) strcpy(resite.seq, resite.chr), 
	    setSeqCode(0, DNA));
	*buf = '+';
	cmplseq((CHAR*) buf + 1, (CHAR*) resite.seq);
	if (strcmp(buf + 1, resite.seq)) strcat(resite.seq, buf);
	return (&resite);
}

RESITE* nextresite()
{
	char	buf[MAXL];

	if (!fgets(buf, MAXL, resfd)) return (0);
	return (getresite(buf));
}

RESITE*	firstresite()
{
	if (!openresite(UNIQENZ)) return (0);
	return (nextresite());
}

RESITE* recogseq(char* src)
{
	char	buf[MAXL];

	if (!openresite(ALL_ENZ)) return (0);
	if (fbisrch(buf, strupr(src), resfd, 0L, nrec, width, (CMPF) cmpf))
		return (getresite(buf));
	else	return (0);
}

static int icmpf(int* a, int* b)
{
	return (*a - *b);
}

int respos(int loc[], int maxloc, Seq* ns, RESITE* res)
{
	int 	num = findseq(loc, maxloc, ns, res->seq);
	qsort(loc, num, sizeof(int), (CMPF) icmpf);
	return (num);
}

void resites(FILE* fd, Seq* ns, char* pat)
{
	int	num = 0;
	char*	pt = parspat(pat);
	char*	pb;
	int	loc[MAXLOC];
	RESITE*	res;

	while (*(pb = pt)) {
	    while (*pt && *pt != ' ' && *pt != '+') pt++;
	    if (*pt) *pt++ = '\0';
	    res = recogseq(pb);
	    if (!res) {
		fprintf(stderr, "%s not found\n", pb);
		continue;
	    }
	    int	nn = respos(loc + num, MAXLOC - num, ns, res);
	    num += nn;
	    fprintf(fd, "%s  (%-10s %-10s %2d )  %d\n", ns->sqname(), 
		res->name, res->chr, res->cut, nn);
	}
	qsort(loc, num, sizeof(int), (CMPF) icmpf);
	putloc(fd, ns, loc, num);
}

#ifdef DEBUG

static	void	deco(CHAR* ps);

static void deco(CHAR* ps)
{
	extern char nucl[];

	while (*ps) {
	    *ps = nucl[*ps];
	    ps++;
	}
}


int main()
{
	char	pat[MAXL];

	while (true) {
	    *pat = '\0';
	    progets(pat, "Enzyme : ");
	    if (!*pat) return (0);
	    RESITE*	rs = recogseq(pat);
	    if (rs) {
		deco(rs->seq);
		printf("%s %s %d\n", rs->name, rs->seq, rs->cut);
	    } else
		printf("%s not found!\n", pat);
	}
}

#endif
