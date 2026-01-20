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
static	const	int	NEG = -1;

	       //lpw blk Nout MaxOut eij ovl fnm rm lg spj 
OUTPRM	OutPrm = {60, 0, 16, 1, 4, 10, 5, 0, 1, 3, 1};

static  const   char*   font_tag[3] = {
	"<b><font color=\"%s\">",
	"<b><font style=\"background-color:%s\">",
	"<b><font color=\"%s\" style=\"background-color:%s\">"
};

static	void	make_fn(const char* fname, const char* ss[]);
static	void	scanfrmt(char* p, const char* s);
static	int	frmtc(const char** ss, char* pb);
static	void	frmts(char* ps);
static	int	raw_scn(char* str, char* format, va_list args);
static	void	display(const char* s, va_list args);
static	INTERACT	crt = {1, 0};

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
	ss[0] = fnam;
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

char* topath(char* res, const char* org)
{
	if (!org) *res = '\0';
	else	strcpy(res, org);
	char*	trm = res + strlen(res) - 1;
	if (*res && *trm != PATHDELM && *trm != ':')
	    *++trm = PATHDELM;
	*++trm = '\0';
	return (res);
}
	
char* makefnam(const char* fnam, const char* defn, char* result)
{
	char const*	sg[5];
	char const*	sd[5];
	char	buf[LINE_MAX];

	if (fnam == result) fnam = strcpy(buf, fnam);
	make_fn(fnam, sg);
	make_fn(defn, sd);
	*result = '\0';
	for (int i = 0; i < 4; ++i) {
const	    int	l = sg[i+1] - sg[i];
	    if (l)
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
	while (strchr(fname, PATHDELM)) {
	    if (strchr(where, 'p'))
		while ((*part++ = *fname++) != PATHDELM);
	    else
		while (*fname++ != PATHDELM);
	}
const	char*	dot = strrchr(fname, '.');
	if (dot) {
	    if (strchr(where, 'b'))
		while (fname != dot) *part++ = *fname++;
	    else
		while (fname != dot) fname++;
	    if (strchr(where, 'e')) {
		while ((*part++ = *fname++));
	    } else
		*part = '\0';
	} else {
	    if (strchr(where, 'b'))
		while ((*part++ = *fname++));
	    else *part = '\0';
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

// sizeof(buf) must be larger than strlen(fn) + strlen(ext)

char* add_ext(const char* fn, const char* ext, char* buf)
{
	strcpy(buf, fn);
	if (!ext || !*ext) return (buf);
	if (*ext == '.') ++ext;
	char*	sls = strrchr(buf, '/');
	if (!sls) sls = buf;
	char*	dot = strrchr(sls, '.');
	if (dot && !strcmp(dot + 1, ext)) return (buf);
	if (dot && !strcmp(dot + 1, "gz")) {	// .ext.gz
const	    char*	e = ext + strlen(ext);
	    bool	eq = true;
	    while (--dot > buf && --e > ext && *dot != '.') 
		if (!(eq = *dot == *e)) break;
	    if (eq) return(buf);
	}
	if (!dot || strlen(dot) > 1) strcat(buf, ".");
	return (strcat(buf, ext));
}

char* get_ext(char* txt, const char* ext)
{
const	char*	e = ext? ext + strlen(ext): 0;
	char*	sls = strrchr(txt, '/');
	if (!sls) sls = txt;
	char*	dot = strrchr(sls, '.');
	char*	t = (dot && !strcmp(dot + 1, "gz"))?
	    dot: (txt + strlen(txt));
	while (--t > sls && *t != '.')
	    if (ext && --e >= ext && *t != *e) break;
	return (((ext && e == ext) || *t == '.')? t: 0);
}

FILE* fopenpbe(const char* path, const char* name, const char* extent, 
	const char* opt, const int& lvl, char* str)
{
	char	buf[LINE_MAX];
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
	    if (extent && *extent) {
		partfnam(str + strlen(str), name, "b");
		if (*extent != '.') strcat(str, ".");
		strcat(str, extent);
	    } else
		partfnam(str + strlen(str), name, "be");
	    if ((fd = fopen(str, opt))) return (fd);
	} while (*pt);

	if (lvl >= 0)
	    fprintf(stderr, "%s: cannot be open!\n", str);
	if (lvl > 0) exit (lvl);
	return (0);
}

FILE* wfopen(const char* name, const char* mode)
{
	if (is_file(name)) {
	    int	c = 'y';
	    if (OutPrm.overwrite == 0) {
		INT	p = crt.prompt;
		crt.prompt = 1;
		c = progetc("Overwrite existing file \"%s\"? [y/n] ", name);
		crt.prompt = p;
	    } else if (OutPrm.overwrite > 1)
		c = 'n';
	    if (tolower(c) != 'y') return (0);
	}
	return (fopen(name, mode));
}

FILE* qopen(const char* dfname, const char* mode)
{
	char	fname[LINE_MAX];
	char	pname[LINE_MAX];
	FILE*	fd;

	topath(pname, dfname);
	do {
	    progets(fname, "File name := %s", pname);
	    if (!*fname) return(0);
	    makefnam(fname, pname, fname);
	} while (!(fd = fopen(fname, mode))) ;
	return (fd);
}

gzFile wgzopen(const char* name, const char* mode)
{
	char	str[LINE_MAX];
	strcpy(str, name);
	if (!is_gz(name)) strcat(str, gz_ext);
	if (is_file(str)) {
	    int	c = 'y';
	    if (OutPrm.overwrite == 0) {
		INT	p = crt.prompt;
		crt.prompt = 1;
		c = progetc("Overwrite existing file \"%s\"? [y/n] ", str);
		crt.prompt = p;
	    } else if (OutPrm.overwrite > 1)
		c = 'n';
	    if (tolower(c) != 'y') return (0);
	}
	return (gzopen(str, mode));
}

gzFile gzopenpbe(const char* path, const char* name, 
	const char* extent, const char* opt, const int& lvl, char* str)
{
	char	buf[LINE_MAX];
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
	    if (extent && *extent) {
		partfnam(str + strlen(str), name, "b");
		if (*extent != '.') strcat(str, ".");
		strcat(str, extent);
	    } else
		partfnam(str + strlen(str), name, "be");
	    if ((fd = gzopen(str, opt))) return (fd);
	} while (*pt);

	if (lvl >= 0)
	    fprintf(stderr, "%s: cannot be open!\n", str);
	if (lvl > 0) exit (lvl);
	return (0);
}

size_t lgzread(gzFile gzfd, char* buf, size_t sz)
{
	char*	tq = buf + sz;
	for (char* sq = buf; sq < tq; sq += INT_MAX) {
const	    int	n = std::min(long(INT_MAX), tq - sq);
	    int	actural_read = gzread(gzfd, sq, n);
	    if (actural_read <= 0) {
		int	z_errnum = 0;
		fatal("gz_read_error %s!\n", gzerror(gzfd, &z_errnum));
	    }
	}
	return (sz);
}

size_t lgzwrite(gzFile gzfd, const char* buf, size_t sz)
{
const	char*	tq = buf + sz;
	for (const char* sq = buf; sq < tq; sq += INT_MAX) {
const	    int	n = std::min(long(INT_MAX), tq - sq);
	    int	actural_write = gzwrite(gzfd, sq, n);
	    if (actural_write <= 0) {
		int	z_errnum = 0;
		fatal("gz_write_error %s!\n", gzerror(gzfd, &z_errnum));
	    }
	}
	return (sz);
}

//	read-only 

int ReadFile::open(const char* fname, const char* extnt)
{
	char	str[LINE_MAX];
	strcpy(str, fname);
// open fname depending on its extent
	char*	sls = strrchr(str, '/');
	if (!sls) sls = str;
	char*	dat = strstr(sls, dat_ext);
	char*	dot = strrchr(sls, '.');	// try .dgz
	if (dot && !strcmp(dot, dgz_ext) && is_filled_file(str) 
	    && (gzfd = gzopen(str, "r"))) return (dtype = 2);
	if (dot && !strcmp(dot, gz_ext) && is_filled_file(str) 
	    && (gzfd = gzopen(str, "r"))) return (dtype = dat? 2: 1);
	if (dot && !strcmp(dot, dat_ext) && is_filled_file(str) 
	    && (fd = fopen(str, "r"))) return (dtype = 2);
	if (is_filled_file(str) && (fd = fopen(str, "r")))
	    return (dtype = 1);

// add extensions
	if (dat) dot = dat;
	if (dot) {
	    *dot = '\0';
	    if (is_filled_file(str) && (fd = fopen(str, "r")))
		return (dtype = 1);
	} else {
	    dot = str + strlen(str);
	}
	if (extnt) {
	    char	ext[MAXL];
	    strcpy(ext, extnt);
	    dat = strstr(ext, dat_ext);
	    bool	isgz = is_gz(ext);
	    if (dat && isgz) {
		strcpy(dot, dgz_ext);		// try .dgz
		if (is_filled_file(str) && (gzfd = gzopen(str, "r")))
		    return (dtype = 2);
		strcpy(dot, datgz);		// try .dat.gz
		if (is_filled_file(str) && (gzfd = gzopen(str, "r")))
		    return (dtype = 2);
	    } else if (isgz) {
		strcpy(dot, gz_ext);		// try .gz
		if (is_filled_file(str) && (gzfd = gzopen(str, "r")))
		    return (dtype = 1);
	    } else if (dat) {
		strcpy(dot, dat_ext);		// try .gz
		if (is_filled_file(str) && (fd = fopen(str, "r")))
		    return (dtype = 2);
	    }
	}
	return (dtype = 0);
}

//	write-only 

void WriteFile::open(const char* fname, const int& text, const bool& gzip, 
	const char* in_name)
{
	char	str[LINE_MAX];
	strcpy(str, fname);
const	char*	sls = strrchr(str, '/');
	if (!sls) sls = str;
const	char*	dot = strrchr(sls, '.');
const	bool	is_dgz = dot && !strcmp(dot, dgz_ext);
const	bool	addat = !is_dgz && !text && (!dot || !strstr(dot, dat_ext));
const	bool	addgz = !is_dgz && gzip && (!dot || strcmp(dot, gz_ext));
const	bool	adtxt = text > 1 && (!dot || !strstr(dot, txt_ext));
	if (addgz && addat) strcat(str, dgz_ext);
	else if (addat) strcat(str, dat_ext);
	else if (addgz) strcat(str, gz_ext);
	else if (adtxt) strcat(str, txt_ext);
	if (in_name && !strcmp(in_name, str))
	    prompt("%s is the same name as input !\n", str);
	if (gzip || is_dgz) gzfd = gzopen(str, "wb");
	else	fd = fopen(str, "wb");
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
	} else if (given == POPUP)	swap(*pval, *back);
	else if (given != SILENT)	{*back = *pval; *pval = given;}
}

void setdblval(const char* mssg, double* pval, double* back, const double*& given)
{
	int	gvn = (int) *given;

	if (gvn == QUERY) {
	    *back = *pval;
	    if (mssg) promptin(mssg, pval);
	} else if (gvn == POPUP) swap(pval, back);
	else if (gvn != SILENT) {
	    *back = *pval;
	    *pval = *given;
	}
}

/*********************************************************************
*
*	class members of Ftable 
*
*********************************************************************/

FILE*	Ftable::fopen(const char* fname, const char* mode)
{
	char	str[LINE_MAX];
	FILE*	fd = 0;
	for (int i = 0; i < n_tabpath; ++i) {
	    if (!tabpath[i]) continue;
	    strcpy(str, tabpath[i]);
	    if (subdir) {
		strcat(str, "/");
		strcat(str, subdir);
		if (is_dir(str)) {
		    fd = fopenpbe(str, fname, 0, mode, -1);
		    if (fd) return (fd);
		}
	    }
	    fd = fopenpbe(tabpath[i], fname, 0, mode, -1);
	    if (fd) return (fd);
	}
	return (0);
}

// fullpath must have enough space to accmodate full path name

char* Ftable::getpath(char* fullpath, const char* fname)
{
	for (int i = 0; i < n_tabpath; ++i) {
	    if (!tabpath[i]) continue;
	    strcpy(fullpath, tabpath[i]);
	    strcat(fullpath, "/");
	    if (subdir) {
		strcat(fullpath, subdir);
		if (is_dir(fullpath)) {
		    strcat(fullpath, "/");
		    strcat(fullpath, fname);
		    if (is_file(fullpath)) return (fullpath);
		}
	    }
	    strcat(fullpath, fname);
	    if (is_file(fullpath)) return (fullpath);
	}
	return (0);
}

FILE*	Ftable::fopen(const char* fname, const char* envpath, const char* defpath)
{
	FILE*	fd = this->fopen(fname, "r");
	if (fd) return (fd);

	char	str[LINE_MAX];
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

DIR* Ftable::dopen(const char* dname, const bool& test)
{
	for (int i = 0; i < 3; ++i) {
	    if (!tabpath[i]) continue;
	    char	str[LINE_MAX];
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
	return (0);
}

void Ftable::setpath(const char* ps, const char* gnm2tab)
{
	setSubDir(ps);
	if (is_dir(ps)) return;
	char    str[LINE_MAX];
// examine dirs
	for (int i = 0; i < 3; ++i) {
	    if (!tabpath[i]) continue;
	    strcpy(str, tabpath[i]);
	    strcat(str, "/");
	    strcat(str, ps);
	    if (is_dir(str)) return;
	}
// consult with gnm2tab
	if (!gnm2tab) {
	    prompt(not_found, ps);
	    return;
	}
	FILE*   fd = this->fopen(gnm2tab, "r");
	if (!fd) {
	    prompt(not_found, gnm2tab);
	    return;
	}
	bool	found = false;
	char*	ts = 0;
	while (fgets(str, LINE_MAX, fd)) {
	    if (*str == '#') continue;
	    char*       qs = str;
	    char*       gs = car(qs);	// identifier
	    if (!*++qs) continue;
	    ts = car(qs);		// ssp group
	    if (strcmp(gs, ps) == 0) {
		found = true;
		break;
	    }
	}
	fclose(fd);
	if (!found) prompt(not_found, ps);
	setSubDir(ts);
	if (is_dir(ts)) return;
	char	sspgrp[MAXL];
	strcpy(sspgrp, ts);
	for (int i = 0; i < 3; ++i) {
	    if (!tabpath[i]) continue;
	    strcpy(str, tabpath[i]);
	    strcat(str, "/");
	    strcat(str, sspgrp);
	    if (is_dir(str)) return;
	}
	prompt(not_found, sspgrp);
}

gzFile Ftable::gzopen(const char* fname, const char* mode)
{
	char	str[LINE_MAX];
	gzFile	fd = 0;
const	char*	ext = is_gz(fname)? 0: gz_ext;
	for (int i = 0; i < n_tabpath; ++i) {
	    if (!tabpath[i]) continue;
	    strcpy(str, tabpath[i]);
	    if (subdir) {
		strcat(str, "/");
		strcat(str, subdir);
		if (is_dir(str)) {
		    fd = gzopenpbe(str, fname, ext, mode, -1);
		    if (fd) return (fd);
		}
	    }
	    fd = gzopenpbe(tabpath[i], fname, ext, mode, -1);
	    if (fd) return (fd);
	}
	return (0);
}

gzFile Ftable::gzopen(const char* fname, const char* envpath, const char* defpath)
{
	gzFile	fd = this->gzopen(fname, "r");
	if (fd) return (fd);
const	char*	ext = is_gz(fname)? 0: gz_ext;

	char	str[LINE_MAX];
	char*	path = envpath? getenv(envpath): 0;
	if (path) {
	    strcpy(str, path);
	    fd = gzopenpbe(str, fname, ext, "r", -1);
	    if (fd) return (fd);
	}
	if (defpath) {
	    strcpy(str, defpath);
	    fd = gzopenpbe(str, fname, ext, "r", -1);
	    if (fd) return (fd);
	}
	return (0);
}

// open for read

int Ftable::ropen(ReadFile& fp, const char* fname, const char* extnt)
{
	fp.open(fname, extnt);
	if (fp.dtype) return (fp.dtype);
	char	str[LINE_MAX];
	for (int i = 0; i < n_tabpath; ++i) {
	    if (!tabpath[i]) continue;
	    strcpy(str, tabpath[i]);
	    if (subdir) {
		strcat(str, "/");
		strcat(str, subdir);
		if (is_dir(str)) {
		    strcat(str, "/");
		    strcat(str, fname);
		    fp.open(str, extnt);
		    if (fp.dtype) return (fp.dtype);
		}
		strcpy(str, tabpath[i]);
	    }
	    strcat(str, "/");
	    strcat(str, fname);
	    fp.open(str, extnt);
	    if (fp.dtype) return (fp.dtype);
	}
	return (0);
}

Ftable ftable;

/*********************************************************************
*
*	Non-portable codes for getting a variable number of data 
*
*********************************************************************/

void setprompt(const int& prom, const int& ech)
{
	crt.prompt = prom;
	crt.echo   = ech;
}

int getprompt()
{
	return (crt.prompt);
}

static void scanfrmt(char* p, const char* s)
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
			fprintf(stderr, str, *i_ptr);
		    break;
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
		    else fprintf(stderr, str, *l_ptr);
		    break;
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
			fprintf(stderr, str, *d_ptr);
		    break;
		default:  fputs(str, stderr); break;
	    }
	}
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

int	fgetiarray(int* array, const int size, FILE* fd)
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

int	getiarray(int* array, const int size, const char* ps)
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

void EscCharCtl::putctl(const int& fg, const int& bg, const int& attr)
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

Gnm2tab::Gnm2tab(const int field) : StrHash<int>()
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

