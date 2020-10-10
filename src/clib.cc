/*****************************************************************************
*
*	Collection of functions for general use
*
*	vmax	vmin	vavsd
*	max3	min3	ax4	min4
*	ipower	lpower	
*	next_wd	wordcmp
*	isBlankLine   car	cdr
*	chop	replace
*	fbisrch
*	comb2	comb3	decomb
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

#include <stdio.h>
#include <ctype.h>
#include <limits.h>
#include "clib.h"
#include "iolib.h"

const	char*	no_space = "No more memory !\n";
const	char*	stddelim = " \t\n\r";

int	thread_num = 0;
int	cpu_num = 1;
int	max_queue_num = 0;

ALGMODE	algmode = {1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0};
//	{nsa, alg, bnd, mlt, aut lcl, lsg, mns, thr, rng, qck, blk, any, crs, slv, lvl. dim}

int ipower(int x, int n)
{
	if (n < 0) return(0);
	int v = 1;
	for ( ; n; n /= 2, x *= x)
		if (n % 2) v *= x;
	return (v);
}

long lpower(long x, int n)
{
	if (n < 0) return(0L);
	long v = 1L;
	for ( ; n; n /= 2, x *= x)
		if (n % 2) v *= x;
	return (v);
}

char* strrealloc(char* dst, const char* src)
{
	delete[] dst;
	if (!src) return (0);
	dst = new char[strlen(src) + 1];
	return (strcpy(dst, src));
}

bool isnumber(const char* ps)
{
	int	n = 0;
	int	d = 0;
	int	e = 0;
	if (*ps == '-') ++ps;
	while (int c = *ps++) {
	    if (isspace(c) || c == '-') break; else
	    if (isdigit(c)) ++n; else
	    if (c == '.') ++d; else
	    if (c == 'e' || c == 'E') {
		++e;
		if (*ps == '-') ++ps;
	    } else return (false);
	}
	return (n > 0 && d < 2 && e < 2);
}

bool isBlankLine(char* str)
{
	while (int c = *str++)
	    if (!isspace(c)) return (false);
	return (true);
}

const char* getarg(int& argc, const char**& argv, bool num, int pv)
{
	if (argv[0][pv]) return argv[0] + pv;
	if (!argv[1]) return (0);
	char	na = argv[1][0];
	if ((argc > 1 && (na != OPTCHAR && (!num || isdigit(na) || na == '.')))
	  || (na == '-' && num && (isdigit(argv[1][1]) || argv[1][1] == '.'))) {
	    --argc;
	    return (*++argv);
	}
	return (0);
}

char* replace(char* dest, const char* sorc, const char* to, const char* from)
{
	int 	slen = strlen(from);
	char	temp[MAXL];
	const	char*	wk = from;
	char*	pt = temp;

	while (char ch = *pt++ = *sorc++) {
	    if (*wk == ch) wk++;
	    else	wk = from;
	    if (!*wk) {
		pt -= slen;
		for (wk = to; *wk; ) *pt++ = *wk++;
		wk = from;
	    }
	}
	return (strcpy(dest, temp));
}

char* car(char*& ps)
{
	while (*ps && isspace(*ps)) ++ps;
	char*	ss = ps;
	while (*ps && !isspace(*ps)) ++ps;
	*ps = '\0';
	return (ss);
}

char* car(char* token, const char* ps)
{
	char*	pt = token;

	while (*ps && isspace(*ps)) ps++;
	while (*ps && !isspace(*ps)) *pt++ = *ps++;
	*pt = '\0';
	return (token);
}

char* car(char* ps, char*& pt, char& cc)
{
	while (*ps && isspace(*ps)) ++ps;
	char*	ss = ps;
	while (*ps && !isspace(*ps)) ++ps;
	pt = ps;
	cc = *ps;
	*ps = '\0';
	return (ss);
}

// copy up to n - 1 characters and \0 terminate

char* carn(char* token, const char* ps, int n)
{
	char*	pt = token;

	while (*ps && isspace(*ps)) ps++;
	while (*ps && !isspace(*ps) && --n > 0) *pt++ = *ps++;
	*pt = '\0';
	return (token);
}

char* cdr(const char* ps)
{
	while (*ps && isspace(*ps)) ++ps;
	while (*ps && !isspace(*ps)) ++ps;
	while (*ps && isspace(*ps)) ++ps;
	return ((char*) ps);
}

char* cdr(const char* ps, const char* delim)
{
	while (*ps && strchr(delim, *ps)) ++ps;
	while (*ps && !strchr(delim, *ps)) ++ps;
	while (*ps && strchr(delim, *ps)) ++ps;
	return ((char*) ps);
}

int chop(char* ps)
{
	if (!*ps) return (0);
	ps += strlen(ps) - 1;
	int	c = *ps;
	*ps = '\0';
	return (c);
}

char* next_wd(char** sp)
{
	char c;

	while ((c = **sp) && isspace(c)) ++*sp;
	char*	top = *sp;
	while ((c = **sp) && !isspace(c) && c != ',') ++*sp;
	return (top);
}

int wordcmp(const char* a, const char* b, int n)
{
 	if (!a) a = "";
 	if (!b) b = "";

	while (*a && *b && !isspace(*a) && !isspace(*b) && n--)
	    if (int i = *a++ - *b++) return (i);
	if (!n || ((!*a || isspace(*a)) && (!*b || isspace(*b))))
	    return (0);
	if (*a && !isspace(*a))
	    return (1);
	return (-1);
}

UPTR fbisrch( UPTR found, const UPTR key, FILE* fd, long left, long right, int width, CMPF cmpf)
{

	while (right - left >= 0) {
	    long	mid = (right + left) / 2;
	    if (fseek(fd, mid * width, 0) == -1L) return (0);
	    if (fread(found, width, 1, fd) == 0) return (0);
	    int	comp = (*cmpf)(key, found);
	    if (comp < 0)		right = mid - 1;
	    else if (comp > 0)	left = mid + 1;
	    else	return (found);
	}
	return (0);
}

int comb2(int i, int j)
{
	if (i == j) return(-1);
	if (i > j) swap(i, j);
	return (j*(j-1)/2 + i);
}

int comb3(int i, int j, int k)
{
	if (i == j || j == k || k == i) return(-1);
	if (i > j) swap(i, j);
	if (i > k) swap(i, k);
	if (j > k) swap(j, k);
	return (k*(k-1)*(k-2)/6 + j*(j-1)/2 + i);
}

// number of combinations

INT mCn(INT m, INT n)
{
	INT	k = m - n;
	if (k < n) n = k;
	INT	c = 1;
	m -= (n - 1);
	for (k = 1; k <= n; c = c * m++ / k++) ;
	return (c);
}

void decomb(int* t, int n, int i)
{
	while (i-- > 0) {
	    int	r = 0;
	    int k = i;
	    for (int s = 0; (s = mCn(k + 1, i + 1)) <= n; ++k)
		r = s;
	    t[i] = k;
	    n -= r;
	}
}

char* strlwr(char* str)
{
	for (char* ps = str; *ps; ps++) *ps = tolower(*ps);
	return (str);
}

char* strupr(char* str)
{
	for (char* ps = str; *ps; ps++) *ps = toupper(*ps);
	return (str);
}

// assume *str consists of writable 7-bit charactors

INT str2uint(const char* str)
{
	INT	i = 0;
	while (*str)
	    i = i * 17 + (*str++ - ' ');
	return i;
}

int countbit(unsigned long x)
{
	int	n = 0;

	for ( ; x; x >>= 1)
	    if (x & 1) ++n;
	return (n);
}

double ktof(const char* str)
{
	double	x = atof(str);
	const	char*	ps = str;
	while (*ps && (isdigit(*ps) || *ps == '.')) ++ps;
	switch (*ps) {
	    case 'k': case 'K': x *= KILO; break;
	    case 'm': case 'M': x *= MEGA; break;
	    case 'g': case 'G': x *= GIGA; break;
	    default: break;
	}
	return (x);
}

Strlist::Strlist(int m, int len) 
	: totallen(0), lastlen(0), maxlen(0), many(0), filled(false)
{
	if (m == 0) {
	    sunitsize = 0; strbuf = 0;
	} else {
	    sunitsize = (len + defsunit - 1) / defsunit * defsunit;
	    strbuf = new char[sunitsize];
	    *strbuf = '\0'; 
	}
	if (m > 1) {
	    punitsize = (m + defpunit - 1) / defpunit * defpunit;
	    idxlst = new INT[punitsize];
	    *idxlst = 0;
	} else {
	    punitsize = 0; idxlst = 0;
	}
}

void Strlist::reset(int m)
{
	if (m > 1 && (INT) m > punitsize) {
	    punitsize = (m + defpunit - 1) / defpunit * defpunit;
	    delete[] idxlst;
	    idxlst = new INT[punitsize];
	    *idxlst = 0;
	}
	*strbuf = '\0';
	many = totallen = lastlen = maxlen = 0;
	filled = false;
}

Strlist::Strlist(const char* str)
	: idxlst(0), punitsize(0), many(1), filled(true)
{
	totallen = lastlen = maxlen = strlen(str) + 1;
	sunitsize = (totallen + defsunit - 1) / defsunit * defsunit;
	strbuf = new char[sunitsize];
	strcpy(strbuf, str);
}

Strlist::Strlist(char* str, const char* delim)
	: lastlen(0), maxlen(0), many(0)
{
	bool	spc = true;
	for (const char* ps = delim; *ps; ++ps)
	    if (!(spc = isspace(*ps))) break;
	for (char* ps = str; ps && *ps; ++many)
	    ps = spc? cdr(ps): cdr(ps, delim);
	totallen = sunitsize = strlen(str) + 1;
	strbuf = new char[totallen];
	punitsize = (many + defpunit - 1) / defpunit * defpunit;
	idxlst = new INT[punitsize + 1];

	char*	ps = str;
	char*	pb = strbuf;
	for (INT ttl = 0, m = 0; *ps; ) {
	    idxlst[m] = ttl;
	    if (spc)
		while (*ps && isspace(*ps)) ++ps;
	    else
		while (*ps && strchr(delim, *ps)) ++ps;
	    char*	qb = pb;
	    if (spc)
		while (*ps && !isspace(*ps)) *pb++ = *ps++;
	    else
		while (*ps && !strchr(delim, *ps)) *pb++ = *ps++;
	    *pb++ = '\0';
	    if ((lastlen = pb - qb)) {
		++m;
		ttl += lastlen;
		if (lastlen > maxlen) maxlen = lastlen;
	    }
	}
	filled = true;
}

void Strlist::format()
{
	char*	ps = strbuf;
	if (!many) {
	    for (INT ttl = 0; ttl < totallen; ++many) {
		lastlen = strlen(ps) + 1;
		ttl += lastlen;
		ps += lastlen;
	    }
	    ps = strbuf;
	}
	punitsize = (many + defpunit - 1) / defpunit * defpunit;
	idxlst = new INT[punitsize];
	for (INT ttl = many = 0; ttl < totallen; ++many) {
	    idxlst[many] = ttl;
	    lastlen = strlen(ps) + 1;
	    if (lastlen > maxlen) maxlen = lastlen;
	    ttl += lastlen;
	    ps += lastlen;
	}
	filled = true;
}

Strlist::Strlist(FILE* fd, int m) 
	: lastlen(0), maxlen(0), many(m)
{
	fseek(fd, 0L, SEEK_END);
	totallen = ftell(fd);
	sunitsize = (totallen + defsunit - 1) / defsunit * defsunit;
	strbuf = new char[sunitsize];
	fseek(fd, 0L, SEEK_SET);
	if (fread(strbuf, sizeof(char), totallen, fd) != (INT) totallen)
	    fatal("Strlist file may be corupped !\n");
	format();
}

#if USE_ZLIB
Strlist::Strlist(gzFile gzfd, int m) 
	: totallen(0), lastlen(0), maxlen(0), many(m)
{
	int	c;
	while ((c = gzgetc(gzfd)) != EOF) ++totallen;
	sunitsize = (totallen + defsunit - 1) / defsunit * defsunit;
	strbuf = new char[sunitsize];
	fseek(gzfd, 0L, SEEK_SET);
	if (fread(strbuf, sizeof(char), totallen, gzfd) <= 0)
	    fatal("Strlist file may be corupped !\n");
	format();
}
#endif

char* Strlist::assign(const Strlist& src)
{
	if (sunitsize < src.sunitsize) {
	    sunitsize = src.sunitsize;
	    delete[] strbuf;
	    strbuf = new char[sunitsize];
	}
	totallen = src.totallen;
	lastlen = src.lastlen;
	maxlen = src.maxlen;
	memcpy(strbuf, src.strbuf, totallen);
	if (punitsize < src.punitsize) {
	    punitsize = src.punitsize;
	    delete[] idxlst;
	    idxlst = new INT[punitsize];
	    *idxlst = 0;
	}
	many = src.many;
	if (src.idxlst) vcopy(idxlst, src.idxlst, many);
	return (strbuf);
}

Strlist::Strlist(Strlist& src)
{
	*this = src;
	strbuf = new char[sunitsize];
	memcpy(strbuf, src.strbuf, totallen);
	idxlst = punitsize? new INT[punitsize]: 0;
	if (src.idxlst) vcopy(idxlst, src.idxlst, many);
}

char* Strlist::assign(const char* str)
{
	totallen = lastlen = maxlen = strlen(str) + 1;
	if (!strbuf) strbuf = new char[sunitsize];
	if (totallen >= sunitsize) {
	    sunitsize *= 2;
	    delete[] strbuf;
	    strbuf = new char[sunitsize];
	}
	return strcpy(strbuf, str);
}

char* Strlist::push(const char* str)
{
	if (!str) str = "";
	if (many >= punitsize) {
	    punitsize = (many + defpunit) / defpunit * defpunit;
	    INT*	tmp = new INT[punitsize];
	    if (idxlst) {
		vcopy(tmp, idxlst, many);
		delete[] idxlst;
	    } else *tmp = 0;
	    idxlst = tmp;
	}
	idxlst[many++] = totallen;
	lastlen = strlen(str) + 1;
	if (lastlen > maxlen) maxlen = lastlen;
	INT	sumlen = totallen + lastlen;
	char*	word = 0;
	if (sumlen >= sunitsize) {
	    if (lastlen > totallen)
		sunitsize = (sumlen + defsunit - 1) / defsunit * defsunit;
	    else	sunitsize *= 2;
	    char* tmp = new char[sunitsize];
	    word = new char[lastlen];
	    memcpy(tmp, strbuf, totallen);
	    strcpy(word, str);		// str may be a part of strbuf
	    delete[] strbuf;
	    strbuf = tmp;
	}
	char*	rv = strbuf + totallen;
	strcpy(rv, word? word: str);
	delete[] word;
	totallen = sumlen;
	return (rv);
}

AddExt::AddExt(int argc, const char** argv, const char* ext) 
	: argv_(argv), arg_ext(0)
{
	int     tl = strlen(argv[0]);
	int     woext = 0;
	for (int i = 1; i < argc; ++i) {
	    if (argv[i][0] == OPTCHAR) tl += strlen(argv[i]);
	    else {
		const char* dot = strrchr(argv[i], '.');
		if (!dot) {
		    tl += strlen(argv[i]) + 4;
		    ++woext;
		} else if (strcmp(dot, ext)) {
		    tl += dot - argv[i] + 4;
		    ++woext;
		} else      tl += strlen(argv[i]);
	    }
	}
	if (!woext) return;
	arg_ext = new char*[argc + 1];
	char*   ps = new char[tl + argc];
	arg_ext[0] = strcpy(ps, argv[0]);
	for (int i = 1; i < argc; ++i) {
	    ps += strlen(ps) + 1;
	    if (argv[i][0] == OPTCHAR) {
		arg_ext[i] = strcpy(ps, argv[i]);
	    } else {
		const char* dot = strrchr(argv[i], '.');
		if (!dot) {
		    arg_ext[i] = strcpy(ps, argv[i]);
		    strcat(ps, ext);
		} else if (strcmp(dot, ext)) {
		    arg_ext[i] = strncpy(ps, argv[i], dot - argv[i]);
		    strcat(ps, ext);
		} else {    
		    arg_ext[i] = strcpy(ps, argv[i]);
		}
	    }
	}
	arg_ext[argc] = 0;
}

//      distribute into bins
//      Input must be sorted on x

PutIntoBins::PutIntoBins(int n, double l, double u, double (*t)(double x), 
	double (*it)(double x), bool savedata)
	: ndiv(n), llmt(l), ulmt(u), width(0), nlitransform(t), invtransform(it), 
	  n_data(0), idx(0), sample_size(0), 
	  vrtl(false), xval(0), freq(0) {
	if (ndiv <= 0) fatal("number of bins <= 0 !\n");
	if (ulmt < llmt) swap(llmt, ulmt);
	width = (ulmt - llmt) / ndiv;
	xval = new int[ndiv + 2] + 1;
	size();
	if (savedata) {
	    freq = new double[(ndiv + 2)];
	    vclear(freq++, (ndiv + 2));
	}
}

int PutIntoBins::size()
{
	if (n_data) return (n_data);
	double	dx = llmt;
	int	px = int(invtransform? invtransform(dx + epsilon): dx) - 1;
	xval[-1] = px;
	for (int i = n_data = 0; i < ndiv; ++i, dx += width) {
	    int	x = int(invtransform? invtransform(dx + epsilon): dx);
	    if (x == px) continue;
	    xval[n_data++] = px = x;
	}
	return (n_data);
}

void PutIntoBins::accumulate(int x, double y)
{
	double	dx = x;
	if (nlitransform) dx = nlitransform(dx);
	while (idx < n_data && xval[idx] < x) ++idx;
	freq[idx] += y;
	sample_size += int(y);
}

void PutIntoBins::normalize(double to, double psc)
{
	if (sample_size == 0) return;
	to /= (1. + psc * (n_data + 2));
	for (int n = 0; n < n_data; ++n)
	    freq[n] = to * (freq[n] / sample_size + psc);
}

static float CriticalValues[98][7] = {
	{0.6836, 0.7808, 0.8850, 0.9411, 0.9763, 0.9881, 0.9940},
	{0.4704, 0.5603, 0.6789, 0.7651, 0.8457, 0.8886, 0.9201},
	{0.3730, 0.4508, 0.5578, 0.6423, 0.7291, 0.7819, 0.8234},
	{0.3173, 0.3868, 0.4840, 0.5624, 0.6458, 0.6987, 0.7437},
	{0.2811, 0.3444, 0.4340, 0.5077, 0.5864, 0.6371, 0.6809},
	{0.3177, 0.3858, 0.4793, 0.5539, 0.6321, 0.6808, 0.7226},
	{0.2883, 0.3519, 0.4404, 0.5114, 0.5876, 0.6346, 0.6757},
	{0.2660, 0.3260, 0.4102, 0.4778, 0.5509, 0.5972, 0.6375},
	{0.2836, 0.3461, 0.4329, 0.5018, 0.5752, 0.6207, 0.6603},
	{0.2646, 0.3242, 0.4076, 0.4740, 0.5455, 0.5903, 0.6290},
	{0.2496, 0.3065, 0.3868, 0.4519, 0.5215, 0.5651, 0.6031},
	{0.2370, 0.2921, 0.3697, 0.4327, 0.5006, 0.5431, 0.5803},
	{0.2263, 0.2793, 0.3547, 0.4162, 0.4827, 0.5243, 0.5613},
	{0.2170, 0.2685, 0.3418, 0.4015, 0.4665, 0.5076, 0.5435},
	{0.2090, 0.2590, 0.3304, 0.3889, 0.4523, 0.4929, 0.5285},
	{0.2019, 0.2507, 0.3203, 0.3779, 0.4406, 0.4807, 0.5161},
	{0.1958, 0.2433, 0.3116, 0.3682, 0.4298, 0.4692, 0.5041},
	{0.1901, 0.2365, 0.3035, 0.3590, 0.4199, 0.4588, 0.4931},
	{0.1851, 0.2306, 0.2962, 0.3511, 0.4106, 0.4493, 0.4833},
	{0.1804, 0.2252, 0.2897, 0.3436, 0.4025, 0.4404, 0.4742},
	{0.1765, 0.2202, 0.2839, 0.3369, 0.3952, 0.4329, 0.4664},
	{0.1725, 0.2156, 0.2781, 0.3302, 0.3874, 0.4244, 0.4573},
	{0.1690, 0.2114, 0.2731, 0.3245, 0.3814, 0.4179, 0.4507},
	{0.1659, 0.2076, 0.2683, 0.3191, 0.3752, 0.4117, 0.4437},
	{0.1628, 0.2039, 0.2641, 0.3143, 0.3702, 0.4062, 0.4383},
	{0.1599, 0.2005, 0.2597, 0.3095, 0.3645, 0.4007, 0.4325},
	{0.1573, 0.1973, 0.2560, 0.3052, 0.3600, 0.3952, 0.4266},
	{0.1550, 0.1946, 0.2524, 0.3012, 0.3552, 0.3905, 0.4219},
	{0.1526, 0.1917, 0.2491, 0.2972, 0.3507, 0.3857, 0.4169},
	{0.1503, 0.1890, 0.2458, 0.2937, 0.3472, 0.3818, 0.4125},
	{0.1485, 0.1868, 0.2430, 0.2904, 0.3432, 0.3776, 0.4081},
	{0.1465, 0.1844, 0.2400, 0.2870, 0.3394, 0.3736, 0.4039},
	{0.1447, 0.1821, 0.2373, 0.2838, 0.3360, 0.3701, 0.4001},
	{0.1430, 0.1800, 0.2347, 0.2810, 0.3331, 0.3669, 0.3967},
	{0.1413, 0.1780, 0.2323, 0.2783, 0.3296, 0.3632, 0.3928},
	{0.1398, 0.1762, 0.2300, 0.2757, 0.3267, 0.3599, 0.3895},
	{0.1382, 0.1743, 0.2277, 0.2728, 0.3235, 0.3564, 0.3867},
	{0.1367, 0.1725, 0.2254, 0.2705, 0.3210, 0.3542, 0.3839},
	{0.1354, 0.1709, 0.2235, 0.2682, 0.3186, 0.3514, 0.3810},
	{0.1342, 0.1695, 0.2215, 0.2660, 0.3160, 0.3485, 0.3779},
	{0.1329, 0.1679, 0.2196, 0.2638, 0.3134, 0.3457, 0.3750},
	{0.1318, 0.1665, 0.2181, 0.2619, 0.3110, 0.3431, 0.3723},
	{0.1305, 0.1651, 0.2162, 0.2599, 0.3093, 0.3415, 0.3704},
	{0.1295, 0.1637, 0.2146, 0.2582, 0.3070, 0.3390, 0.3674},
	{0.1284, 0.1623, 0.2129, 0.2560, 0.3045, 0.3364, 0.3649},
	{0.1273, 0.1612, 0.2114, 0.2542, 0.3026, 0.3345, 0.3634},
	{0.1263, 0.1598, 0.2096, 0.2525, 0.3006, 0.3324, 0.3606},
	{0.1253, 0.1587, 0.2083, 0.2508, 0.2989, 0.3305, 0.3584},
	{0.1245, 0.1575, 0.2070, 0.2494, 0.2972, 0.3287, 0.3565},
	{0.1235, 0.1564, 0.2056, 0.2476, 0.2951, 0.3266, 0.3545},
	{0.1226, 0.1555, 0.2044, 0.2464, 0.2938, 0.3250, 0.3525},
	{0.1219, 0.1545, 0.2030, 0.2446, 0.2918, 0.3229, 0.3511},
	{0.1209, 0.1533, 0.2016, 0.2433, 0.2901, 0.3208, 0.3485},
	{0.1202, 0.1525, 0.2006, 0.2418, 0.2887, 0.3198, 0.3474},
	{0.1194, 0.1515, 0.1994, 0.2405, 0.2872, 0.3180, 0.3452},
	{0.1186, 0.1505, 0.1982, 0.2392, 0.2856, 0.3165, 0.3436},
	{0.1179, 0.1495, 0.1970, 0.2379, 0.2842, 0.3148, 0.3420},
	{0.1172, 0.1487, 0.1959, 0.2365, 0.2823, 0.3126, 0.3398},
	{0.1165, 0.1478, 0.1949, 0.2354, 0.2814, 0.3119, 0.3389},
	{0.1159, 0.1470, 0.1938, 0.2341, 0.2800, 0.3104, 0.3377},
	{0.1152, 0.1463, 0.1928, 0.2329, 0.2786, 0.3091, 0.3361},
	{0.1146, 0.1455, 0.1921, 0.2320, 0.2776, 0.3077, 0.3346},
	{0.1140, 0.1448, 0.1911, 0.2308, 0.2763, 0.3065, 0.3332},
	{0.1134, 0.1442, 0.1902, 0.2299, 0.2752, 0.3051, 0.3316},
	{0.1128, 0.1433, 0.1893, 0.2290, 0.2742, 0.3041, 0.3307},
	{0.1122, 0.1426, 0.1883, 0.2279, 0.2730, 0.3028, 0.3297},
	{0.1116, 0.1418, 0.1873, 0.2268, 0.2718, 0.3014, 0.3284},
	{0.1112, 0.1413, 0.1868, 0.2260, 0.2706, 0.3004, 0.3267},
	{0.1106, 0.1406, 0.1858, 0.2249, 0.2701, 0.2997, 0.3255},
	{0.1101, 0.1400, 0.1850, 0.2238, 0.2684, 0.2979, 0.3241},
	{0.1096, 0.1394, 0.1843, 0.2231, 0.2673, 0.2967, 0.3231},
	{0.1091, 0.1387, 0.1834, 0.2221, 0.2667, 0.2961, 0.3227},
	{0.1088, 0.1383, 0.1829, 0.2215, 0.2657, 0.2947, 0.3213},
	{0.1082, 0.1378, 0.1822, 0.2206, 0.2647, 0.2938, 0.3199},
	{0.1078, 0.1372, 0.1815, 0.2199, 0.2638, 0.2928, 0.3190},
	{0.1072, 0.1365, 0.1807, 0.2189, 0.2626, 0.2918, 0.3181},
	{0.1068, 0.1360, 0.1799, 0.2180, 0.2618, 0.2909, 0.3169},
	{0.1064, 0.1354, 0.1793, 0.2174, 0.2609, 0.2900, 0.3159},
	{0.1060, 0.1349, 0.1786, 0.2164, 0.2602, 0.2889, 0.3150},
	{0.1056, 0.1344, 0.1781, 0.2160, 0.2593, 0.2881, 0.3140},
	{0.1051, 0.1338, 0.1773, 0.2151, 0.2583, 0.2873, 0.3132},
	{0.1047, 0.1334, 0.1769, 0.2143, 0.2574, 0.2860, 0.3120},
	{0.1043, 0.1328, 0.1762, 0.2137, 0.2565, 0.2852, 0.3113},
	{0.1039, 0.1324, 0.1755, 0.2131, 0.2561, 0.2844, 0.3100},
	{0.1035, 0.1319, 0.1749, 0.2123, 0.2550, 0.2838, 0.3094},
	{0.1031, 0.1315, 0.1745, 0.2117, 0.2543, 0.2826, 0.3085},
	{0.1028, 0.1309, 0.1739, 0.2110, 0.2539, 0.2824, 0.3078},
	{0.1025, 0.1306, 0.1733, 0.2102, 0.2528, 0.2809, 0.3063},
	{0.1020, 0.1302, 0.1726, 0.2097, 0.2522, 0.2806, 0.3064},
	{0.1016, 0.1297, 0.1721, 0.2091, 0.2514, 0.2796, 0.3050},
	{0.1014, 0.1293, 0.1715, 0.2084, 0.2508, 0.2788, 0.3043},
	{0.1011, 0.1289, 0.1712, 0.2079, 0.2503, 0.2784, 0.3038},
	{0.1006, 0.1285, 0.1706, 0.2072, 0.2492, 0.2774, 0.3027},
	{0.1004, 0.1280, 0.1701, 0.2068, 0.2487, 0.2768, 0.3020},
	{0.1001, 0.1276, 0.1696, 0.2061, 0.2479, 0.2762, 0.3013},
	{0.0998, 0.1274, 0.1691, 0.2056, 0.2475, 0.2755, 0.3005},
	{0.0994, 0.1268, 0.1685, 0.2049, 0.2467, 0.2744, 0.2994},
	{0.0992, 0.1266, 0.1682, 0.2044, 0.2461, 0.2740, 0.2995}
};

Dixon::Dixon(float a) : alpha(a)
{
static	float	pval[7] = {0.30, 0.20, 0.10, 0.05, 0.02, 0.01, 0.005};
	for (elt = 0; elt < 7; ++elt) {
	    if (alpha > pval[elt]) break;
	}
	--elt;
}

int Dixon::dixon(int* res, float* data, int* odr, int num, float min_deno)
{
	float	dtmax = (float) data[odr[num - 1]];
	float	dtmin = (float) data[odr[0]];
	float	deno = dtmax - dtmin;
	float	rs = 0.;
	float	rl = 0.;
	if (num > 2 && num <= 7) {
	    if (deno > min_deno) {
		rs = (data[odr[1]] - dtmin) / deno;
		rl = (dtmax - data[odr[num - 2]]) / deno;
	    }
	} else if (num >= 8 && num <= 10) {
	    deno = data[odr[num - 2]] - dtmin;
	    if (deno > min_deno) rs = (data[odr[1]] - dtmin) / deno;
	    deno = dtmax - data[odr[1]];
	    if (deno > min_deno) rl = (dtmax - data[odr[num - 1]]) / deno;
	} else if (num >= 11) {
	    deno = data[odr[num - 1]] - dtmin;
	    if (deno > min_deno) rs = (data[odr[2]] - dtmin) / deno;
	    deno = dtmax - data[odr[1]];
	    if (deno > min_deno) rl = (dtmax - data[odr[num - 2]]) / deno;
	} else {
	    return 0;
	}
	int	nn = num > 100? 100: num;
	float	thr = CriticalValues[nn - 3][elt];
	int	ol = 0;
	if (rl >= thr) res[ol++] = odr[num - 1];
	if (rs >= thr) res[ol++] = -*odr++ - 1;
	if (ol) ol += dixon(res + ol, data, odr, num - ol, min_deno);
	return ol;
}

RandNumGen::RandNumGen(double* pcmp, INT dim, INT rslv)
	: resol(rslv)
{
	buf = new INT[resol];
	INT	k = 0;
	double	acc = 0;
	for (INT i = 0; i < dim; ++i) {
	    acc += *pcmp++;
	    INT	r = min((INT) (acc * resol), resol);
	    while (k < r) buf[k++] = i;
	}
}

//	Knuth, pg 117, 3.4.1, vol. 2 2nd edition

double rgauss()
{
	double	s, v1;

	do {
	    v1 = 2.*drand48()-1.;
	    double	v2 = 2.*drand48()-1.;
	    s = v1*v1+v2*v2;
	} while (s >= 1);
	s = sqrt((-2.*log(s))/s);
	return (v1 * s);
}

//	Random numbers that follow Poisson distribution

int rpoisson(double mu)
{
	double	p = drand48();

	mu = exp(-mu);
	int	i = 0;
	for ( ; p >= mu; ++i) p *= drand48();
	return (i);
}

