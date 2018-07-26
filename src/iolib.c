/*****************************************************************************
*
*	Collection of functions for general use
*
*	qopen		dopen		qout
*	makefnam	partfnam	topath
*	getarray	fgetarray	sgetarray
*	promptin	prompt		progetc		progets
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

#define  _IOLIBC_
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>
#include "adddef.h"
#include "clib.h"
#include "iolib.h"

#define isformchar(c) ((c)=='h'||(c)=='l'||(c)=='L')
#define isnumer(c) (isdigit(c)||c=='.'||c=='+'||c=='-'||c=='e'||c=='E')

FILE*	out_fd = stdout;
static	const	int	MAXARG = 4;
static	const	int	NEG = -1;
static	const	int	gzCompress_rate = 10;

const	char*	no_file = "\'%s\' cannot be created!\n";
const	char*	not_found = "\'%s\' is not found!\n";
const	char*	gz_unsupport = "compressed %s is not supported!\n";

static	const	char*	gnm2tab = "gnm2tab";
static  const   char*   font_tag[3] = {
	"<b><font color=\"%s\">",
	"<b><font style=\"background-color:%s\">",
	"<b><font color=\"%s\" style=\"background-color:%s\">"
};

static	void	make_fn(const char* fnam, const char* ss[]);
static	void	scanfrmt( char* p, const char* s);
static	int	frmtc(const char** ss, char* pb);
static	void	frmts(char* ps);
static	int	raw_scn(char* str, char* format, va_list args);
static	void	display(const char* s, va_list args);

#if LEAKTRACE
#include "MemoryTrace.hpp"
LeakTrace	lt;
#endif

char* fgets_wocr(char* str, INT maxl, FILE* fd)
{
	int	c = -1;
	char*	p = str;

	while (--maxl > 0 && (c = getc(fd)) != -1)
	    if (c != '\r' && (*p++ = c) == '\n')
		break;
	*p = '\0';
	return ((c == -1 && p == str)? 0: str);
}

static void make_fn(const char* fnam, const char* ss[])
{
	ss[0] = (char*) fnam;
	ss[1] = ss[2] = ss[3] = ss[4] = 0;
	if (!fnam) return;
	while (*fnam) {
	    switch (*fnam) {
		case ':':	ss[1] = fnam + 1; break;
		case PATHDELM:	ss[2] = fnam + 1; ss[3] = 0; break;
		case '.':	ss[3] = fnam; break;
	    }
	    fnam++;
	}
	if (!ss[1]) ss[1] = ss[0];
	if (!ss[2]) ss[2] = ss[1];
	if (!ss[3]) ss[3] = fnam;
	ss[4] = fnam;
}

char* topath( char* res, const char* org)
{
	char* trm;

	if (!org) *res = '\0';
	else	strcpy(res, org);
	trm = res + strlen(res) - 1;
	if (*res && *trm != PATHDELM && *trm != ':')
	    *++trm = PATHDELM;
	*++trm = '\0';
	return (res);
}
	
char* makefnam(const char* fnam, const char* defn, char* result)
{
	char const*	sg[5];
	char const*	sd[5];
	char	buf[MAXL];
	int 	i, l;

	if (fnam == result) fnam = strcpy(buf, fnam);
	make_fn(fnam, sg);
	make_fn(defn, sd);
	*result = '\0';
	for (i = 0; i < 4; i++) {
	    if ((l = sg[i+1] - sg[i]))
		strncat(result, sg[i], l);
	    else
		strncat(result, sd[i], sd[i+1] - sd[i]);
	}
	return (car(result, result));
}

/*	Extract parts of filename.  Charcters in where indicate
* 	'd': drive, 'p': path, 'b' body, 'e': extention.
*/

char* partfnam(char* part, const char* fname, const char* where)
{
	char *top = part;

	if (strchr(fname, ':')) {
	    if (strchr(where, 'd'))
		while ((*part++ = *fname++) != ':');
	    else
		while (*fname++ != ':');
	}
	while (strchr(fname, PATHDELM))
	    if (strchr(where, 'p'))
		while ((*part++ = *fname++) != PATHDELM);
	    else
		while (*fname++ != PATHDELM);
	if (strchr(fname, '.')) {
	    if (strchr(where, 'b'))
		while (*fname != '.') *part++ = *fname++;
	    else
		while (*fname != '.') fname++;
	    if (strchr(where, 'e')) {
		while ((*part++ = *fname++));
	    } else
		*part = '\0';
	}
	else {
	    if (strchr(where, 'b')) {
		while ((*part++ = *fname++));
	    } else *part = '\0';
	}
	return (top);
}

size_t file_size(FILE* fd)
{
	if (!fd) return (0);
	size_t	fpos = ftell(fd);
	fseek(fd, 0L, SEEK_END);
	size_t	fs = ftell(fd);
	fseek(fd, fpos, SEEK_SET);
	return (fs);
}

size_t file_size(const char* fn)
{
	FILE*	fd = fopen(fn, "r");
	if (!fd) return (0);
	fseek(fd, 0L, SEEK_END);
	size_t	fs = ftell(fd);
	fclose(fd);
	return (fs);
}
	
FILE* fopenpbe(const char* path, const char* name, 
	const char* extent, const char* opt, int lvl, char* str)
{
	char	buf[MAXL];
	if (!str) str = buf;

	FILE*	fd = 0;
	char*	ps;
const 	char*	pt = path? path: "";

	do {
	    for (ps = str; *pt; *ps++ = *pt++) {
		if (*pt == ';') {
		    pt++;
		    break;
		}
	    }
	    if (ps > str && ps[-1] != PATHDELM)
		*ps++ = PATHDELM;
	    *ps = '\0';
	    if (extent) {
		partfnam(str + strlen(str), name, "b");
		if (*extent && *extent != '.') strcat(str, ".");
		strcat(str, extent);
	    } else
		strcat(str, name);
	    if ((fd = fopen(str, opt))) return (fd);
	} while (*pt);

	if (lvl >= 0)
	    fprintf(stderr, "%s: cannot be open!\n", str);
	if (lvl > 0) exit (lvl);
	return (0);
}

FILE* qopen(const char* dfname, const char* mode)
{
	char	fname[MAXL];
	char	pname[MAXL];
	FILE*	fd;

	topath(pname, dfname);
	do {
	    progets(fname, "File name := %s", pname);
	    if (!*fname) return(0);
	    makefnam(fname, pname, fname);
	} while (!(fd = fopen(fname, mode))) ;
	return (fd);
}

static int	dvc = 0;
static FILE*	prn = 0;

FILE* qout(const char* dfn)
{
	dvc = 0;
askdev:	
	promptin("Print to -1[Non]/0[stderr]/1[stdout]/2[PRN]/3[File] (%d): ", 
	&dvc);
	if (dvc < 0) prn = 0;
	else if (dvc == 0) prn = stderr;
	else if (dvc == 1) prn = stdout;
	else if (dvc == 2) {
	} else if (dvc == 3) {
	    prn = qopen(dfn, "w");
	    if (!prn) goto askdev;
	} else {
	    prn = qopen(dfn, "a");
	    if (!prn) goto askdev;
	}
	return (prn);
}

void qclose()
{
	if (dvc == 2) pclose(prn); else 
	if (dvc > 2) fclose(prn);
	dvc = 0;
}

void setintval(const char* mssg, int* pval, int* back, int given)
{
	if (given == QUERY) {
	    *back = *pval;
	    if (mssg) promptin(mssg, pval);
	} else if (given == POPUP)	gswap(*pval, *back);
	else if (given != SILENT)	{*back = *pval; *pval = given;}
}

void setdblval(const char* mssg, double* pval, double* back, const double* given)
{
	int	gvn = (int) *given;

	if (gvn == QUERY) {
	    *back = *pval;
	    if (mssg) promptin(mssg, pval);
	} else if (gvn == POPUP) gswap(pval, back);
	else if (gvn != SILENT) {
	    *back = *pval;
	    *pval = *given;
	}
}

FILE*	Ftable::fopen(const char* fname, const char* mode)
{
	char	str[MAXL];
	FILE*	fd = 0;
	for (int i = 0; i < n_tabpath; ++i) {
	    if (tabpath[i]) {
		strcpy(str, tabpath[i]);
		if (subdir) {
		    strcat(str, "/");
		    strcat(str, subdir);
		    fd = fopenpbe(str, fname, 0, mode, -1);
		    if (fd) return (fd);
		}
		fd = fopenpbe(tabpath[i], fname, 0, mode, -1);
		if (fd) return (fd);
	    }
	}
	fatal("%s not found! Confirm whether ALN_TAB is correctly set!\n", fname);
	return (0);
}

// fullpath must have enough space to accmodate full path name

char* Ftable::getpath(char* fullpath, const char* fname)
{
	FILE*	fd = 0;
	for (int i = 0; i < n_tabpath; ++i) {
	    if (tabpath[i]) {
		if (subdir) {
		    char  str[MAXL];
		    strcat(str, tabpath[i]);
		    strcat(str, "/");
		    strcat(str, subdir);
		    fd = fopenpbe(str, fname, 0, "r", -1, fullpath);
		    if (fd) break;
		}
		fd = fopenpbe(tabpath[i], fname, 0, "r", -1, fullpath);
		if (fd) break;
	    }
	}
	if (fd) {
	    fclose(fd);
	    return (fullpath);
	}
	return (0);
}

FILE*	Ftable::fopen(const char* fname, const char* envpath, const char* defpath)
{
	FILE*	fd = this->fopen(fname, "r");
	if (fd) return (fd);

	char	str[MAXL];
	char*	path = envpath? getenv(envpath): 0;
	if (path) {
	    strcpy(str, path);
	    fd = fopenpbe(str, fname, 0, "r", -1);
	    if (fd) return (fd);
	}
	if (defpath) {
	    strcpy(str, defpath);
	    fd = fopenpbe(str, fname, 0, "r", -1);
	    if (fd) return (fd);
	}
	return (0);
}

DIR* Ftable::dopen(const char* dname, bool test)
{
	for (int i = 0; i < 3; ++i) {
	    if (tabpath[i]) {
		char	str[MAXL];
		strcpy(str, tabpath[i]);
		strcat(str, "/");
		strcat(str, dname);
		DIR*	dp = opendir(str);
		if (dp) {
		    setSubDir(dname);
		    if (test) closedir(dp);
		    return (dp);
		}
	    }
	}
	return (0);
}

void Ftable::setpath(const char* ps, const char* convtab)
{
	setSubDir(ps);
	DIR*    dp = dopen(ps, true);
	if (dp) return;
	if (!convtab) {
	    prompt(not_found, ps);
	    return;
	}
	FILE*   fd = this->fopen(convtab, "r");
	if (!fd) {
	    prompt(not_found, convtab);
	    return;
	}
	char    str[MAXL];
	while (fgets(str, MAXL, fd)) {
	    if (*str == '#') continue;
	    char*       qs = str;
	    char*       gs = car(qs);
	    if (!*++qs) continue;
	    char*	ts = car(qs);
	    if (strcmp(gs, ps) == 0) {
		dp = dopen(ts, true);
		break;
	    }
	    if (!*++qs) continue;
	    gs = car(qs);
	    if (strcmp(gs, ps)) continue;
	    dp = dopen(ts, true);
	    break;
	}
	fclose(fd);
	if (!dp) prompt(not_found, ps);
}

Ftable ftable;

/*********************************************************************
*
*	Non-portable codes for getting a variable number of data 
*
*********************************************************************/


static INTERACT crt = {1, 0};

void setprompt(int prom, int ech)
{
	crt.prompt = prom;
	crt.echo   = ech;
}

int getprompt()
{
	return (crt.prompt);
}

static void scanfrmt( char* p, const char* s)
{
	int  flag = 0;

	while (*s) {
	    if (*s == '%') {
		if (s[1] == '%') s++;
		else {
		    *p++ = *s;
		    flag = 1;
		}
	    } else if (isalpha(*s)) {
		if (flag == 1) {
		    if (isformchar(*s)) *p++ = *s++;
		    *p++ = *s;
		    *p++ = ' ';
		    flag = 0;
		}
	    }
	    s++;
	}
	*p = '\0';
}

static int frmtc(const char** ss, char* pb)
{
	int 	c = 0;
const	char*	ps = *ss;

cartop:
	while (*ps && *ps != '%') *pb ++ = *ps++;
	if (*ps && *(ps + 1) == '%') {
	    *pb++ = *ps++; ps++;
	    goto cartop;
	}
	while (*ps && !isalpha(c = *pb++ = *ps++)) ;
	if (c == 'l' || c == 'L')
	    c = toupper(*pb++ = *ps++);
	if (c == 'h' || c == 'H')
	    c = tolower(*pb++ = *ps++);
	*ss = ps;
	*pb = '\0';
	return (c);
}

static void frmts(char* ps)
{
cartop:
	while (*ps != '%')
	    if (!*ps++) return;
	if (*++ps == '%') {
	    ++ps;
	    goto cartop;
	}
	*ps++ = 's';
	*ps = '\0';
}

static int raw_scn(char* str, char* format, va_list args)
{
	char	csv;
	char	psv;
	int 	n = 0;
	char*	fmt = format;
	char*	ps = str;
	char*	ptr;

	for ( ; *(str = next_wd(&ps)) && *(format = next_wd(&fmt));
	     n++) {
		csv = *fmt;
		psv = *ps;
		*fmt = *ps = '\0';
		ptr = va_arg(args, char *);
		if (ps > str) sscanf(str, format, ptr);
		*fmt++ = csv;
		*ps++ = psv;
	}
	return (n);
}

static void display(const char* s, va_list args)
{
	char	str[MAXL];
	int 	c;
	char*	c_ptr;
	int*	i_ptr;
	INT*	u_ptr;
	long*	l_ptr;
	double*	d_ptr;

	while (*s) {
	    c = frmtc(&s, str);
	    switch (c) {
		case 'c':
		    c_ptr = va_arg(args, char *);
		    fprintf(stderr, str, (char) *c_ptr); break;
		case 's':
		    c_ptr = va_arg(args, char *);
		    fprintf(stderr, str, c_ptr); break;
		case 'd':
		case 'o':
		case 'x':
		    i_ptr = va_arg(args, int *);
		    if (*i_ptr == INT_MAX || *i_ptr == INT_MIN)
			frmts(str);
		    if (*i_ptr == INT_MAX)
			fprintf(stderr, str, "+++");
		    else if (*i_ptr == INT_MIN)
			fprintf(stderr, str, "---");
		    else
			fprintf(stderr, str, *i_ptr); break;
		case 'u':
		    u_ptr = va_arg(args, INT *);
		    fprintf(stderr, str, *u_ptr); break;
		case 'D':
		case 'O':
		case 'X':
		    l_ptr = va_arg(args, long *);
		    if (*l_ptr == LONG_MAX || *l_ptr == LONG_MIN)
			frmts(str);
		    if (*l_ptr == LONG_MAX)
			fprintf(stderr, str, "+++");
		    else if (*l_ptr == LONG_MIN)
			fprintf(stderr, str, "---");
		    else fprintf(stderr, str, *l_ptr); break;
		case 'e':
		case 'f':
		case 'g':
		case 'E':
		case 'F':
		case 'G':
		    d_ptr = va_arg(args, double *);
		    if ((int) *d_ptr == QUERY) {
			frmts(str);
			fprintf(stderr, str, "*");
		    } else
			fprintf(stderr, str, *d_ptr); break;
		default:  fputs(str, stderr); break;
	    }
	}
}

void fatal(const char* format,...)
{
	va_list	args;

	va_start(args, format);
	(void) vfprintf(stderr, format, args);
	putc('\n', stderr);
	va_end(args);
	exit(1);
}

void prompt(const char* format,...)
{
	va_list	args;

	if (crt.prompt) {
	    va_start(args, format);
	    vfprintf(stderr, format, args);
	    va_end(args);
	}
}

int progetc(const char* format,...)
{
	char str[MAXL];
	va_list	args;

	if (crt.prompt) {
	    va_start(args, format);
	    vfprintf(stderr, format, args);
	    va_end(args);
	}
	if (!fgets(str, MAXL, stdin)) return (EOF);
	if (crt.echo) fputs(str, stderr);
	return (*str);
}

char* progets(char* str, const char* format,...)
{
	va_list	ap;

	if (crt.prompt) {
	    va_start(ap, format);
	    vfprintf(stderr, format, ap);
	    va_end(ap);
	}
	if (!fgets(str, MAXL, stdin)) return (0);
	if (crt.echo) fputs(str, stderr);
	str[strlen(str) - 1] = '\0';
	return (str);
}

int promptin(const char* format,...)
{
	char	str[MAXL];
	char	frmt[MAXL];
	int 	n;
	va_list	args;

	if (crt.prompt) {
	    va_start(args, format);
	    display(format, args);
	    va_end(args);
	}
	scanfrmt(frmt, format);
	if (!*frmt || !fgets(str, MAXL, stdin)) return (0);
	va_start(args, format);
#ifdef SYSTEM_V
	n = vsscanf(str, frmt, args);
#else  
	n = raw_scn(str, frmt, args);
#endif
	va_end(args);
	if (crt.echo) fputs(str, stderr);
	return (n);
}

int sinputf(const char* str, const char* format,...)
{
	char	buf[MAXL];
	char	frmt[MAXL];
	int 	n;
	va_list	args;

	scanfrmt(frmt, format);
	if (!*frmt) return(0);
	va_start(args, format);
#ifdef SYSTEM_V
	n = vsscanf(str, frmt, args);
#else
	n = raw_scn(strcpy(buf, str), frmt, args);
#endif
	va_end(args);
	return (n);
}

int	sgetiarray(int* array, int size, const char** pps)
{
	const	char*	ns = 0;
	int	num;
	int 	n = 0;
	int 	i = 0;
	int	flg = 0;

	for ( ; **pps; ++(*pps)) {
	    if (isdigit(**pps)) {
		if (!ns) ns = *pps;
	    } else {
		if (ns) {
		    num = atoi(ns);
		    ns = 0;
		    if (flg == 1) {
			for (i = array[-1]; i < num && n < size; ++n)
			    *array++ = ++i;
			flg = 0;
		    } else if (flg == NEG) {
			for (i = array[-1]; i < num && n < size; ++n) {
			    *++array = ++i; ++array;
			}
			break;
		    } else if (n < size) {
			*array++ = num;
			++n;
		    }
		}
		if (**pps == '-' && n) flg = 1;
		if (**pps == '/') {
		    if (flg == 1) flg = NEG;
		    else 	break;
		}
	    }
	}
	return (flg == NEG? -n: n);
}

int	fgetiarray(int* array, int size, FILE* fd)
{
static	char	buf[MAXL] = {'\0'};
static	char*	ps = buf;
	int	num;
	int 	n = 0;

	if (size == 0 || !array) {
	    ps = buf;
	    *buf = '\0';
	    return (0);
	}
	while (n < size) {
	    if (!ps || !*ps || *ps == '\n') {
		if (fd != stdin) ps = fgets(buf, MAXL, fd);
		else {
		    prompt(": ");
		    ps = fgets(buf, MAXL, stdin);
		}
		if (!ps) return (n);
		if (*ps == ';') ps += 2;
	    }
	    num = sgetiarray(array + n, size - n, (const char**) &ps);
	    n += num;
	    if (*ps == '/' || num < 0) break;
	} 
	if (ps && *ps) ++ps;
	return (n);
}

int	getiarray(int* array, int size, const char* ps)
{
	FILE*	fd = 0;

	if (ps && *ps) {
	    while (*ps && isspace(*ps)) ++ps;
	    if (isdigit(*ps)) return (sgetiarray(array, size, &ps));
	    fd = fopen(ps, "r");
	    if (!fd) fprintf(stderr, "%s can't open!\n", ps);
	}
	return (fgetiarray(array, size, fd));
}

void EscCharCtl::putctl(int fg, int bg, int attr)
{
	int	nterm = 0;
	if (fg) {
	    fputs(esc_code, fd);
	    fputc('3', fd);
	    fputc(fg, fd);
	    ++nterm;
	}
	if (bg) {
	    if (nterm++) fputc(';', fd);
	    else fputs(esc_code, fd);
	    fputc('4', fd);
	    fputc(bg, fd);
	}
	if (attr) {
	    if (nterm++) fputc(';', fd);
	    else fputs(esc_code, fd);
	    fputc(attr, fd);
	}
	if (nterm) fputc('m', fd);
}

void HtmlCharCtl::putctl(const char* fg, const char* bg)
{
	if (fg && bg)	fprintf(fd, font_tag[2], fg, bg);
	else if (fg)	fprintf(fd, font_tag[0], fg);
	else if (bg)	fprintf(fd, font_tag[1], bg);
}

HtmlCharCtl::HtmlCharCtl(FILE* _fd, const char* title) : fd(_fd)
{
	fputs("<html>\n", fd);
	if (title) {
	    fputs("<head>\n", fd);
	    fprintf(fd, "<title>%s</title>\n", title);
	    fputs("</head>\n", fd);
	}
	fputs("<body>\n", fd);
	fputs("<p>\n<pre>\n", fd);
}

HtmlCharCtl::~HtmlCharCtl()
{
	fputs("</pre>\n</p>\n</body>\n", fd);
}

Gnm2tab::Gnm2tab(int field) : StrHash<int>()
{
	FILE*	fd = ftable.fopen(gnm2tab, "r");
	if (!fd) fatal("Can't open %s !\n", gnm2tab);
	char	str[MAXL];
	int	specs = 0;
	while (char* ps = fgets(str, MAXL, fd)) {
	    if (*ps == '#' || *ps == '\n') continue;
	    ++specs;
	}
	resize(specs);
	rewind(fd);
	domphy = new StrHash<int>(specs / 10);
	while (char* ps = fgets(str, MAXL, fd)) {
	    if (*ps == '#' || *ps == '\n') continue;
	    char*	id = car(ps);
	    char*	qs = 0;
	    int 	c = 1;
	    for ( ; c < field; ++c)
		if (!*++ps || !(qs = car(ps))) break;
	    if (c < field || !qs) continue;
	    KVpair<INT, int>*	kv = domphy->pile(qs);
	    assign(id, kv->val);
	}
	fclose(fd);
}

int Gnm2tab::taxon_code(const char* sqid, char** taxon)
{
	char	genspc[10];
	char*	ps = genspc;
	const char*	tg = sqid + 8;
	while (*sqid && sqid < tg) *ps++ = tolower(*sqid++);
	*ps = '\0';
	KVpair<INT, int>* kv = find(genspc);
	if (taxon) *taxon = kv? domphy->strkey(kv->val): 0;
	return (kv? kv->val: ERROR);
}

#if USE_ZLIB

gzFile gzopenpbe(const char* path, const char* name, 
	const char* extent, const char* opt, int lvl, char* str)
{
	char	buf[MAXL];
	if (!str) str = buf;

	gzFile	fd = 0;
	char*	ps;
const 	char*	pt = path? path: "";

	do {
	    for (ps = str; *pt; *ps++ = *pt++) {
		if (*pt == ';') {
		    pt++;
		    break;
		}
	    }
	    if (ps > str && ps[-1] != PATHDELM)
		*ps++ = PATHDELM;
	    *ps = '\0';
	    if (extent) {
		partfnam(str + strlen(str), name, "b");
		if (*extent && *extent != '.') strcat(str, ".");
		strcat(str, extent);
	    } else
		strcat(str, name);
	    if ((fd = gzopen(str, opt))) return (fd);
	} while (*pt);

	if (lvl >= 0)
	    fprintf(stderr, "%s: cannot be open!\n", str);
	if (lvl > 0) exit (lvl);
	return (0);
}

#endif

