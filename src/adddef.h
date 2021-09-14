/*****************************************************************************
*
*	Collection of commonly used headers
*
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

#ifndef  _ADDDEF_
#define  _ADDDEF_


#define UPTR	void*

#include <stdlib.h>
#include <string.h>

typedef unsigned char CHAR;
typedef	unsigned short SHORT;
typedef unsigned int  INT;
typedef unsigned long LONG;

typedef int (*CMPF)(const UPTR, const UPTR);

#define ERROR	(-1)
#define OK	0
#define ABORT	INT_MIN
#define MAXL	256
#ifndef LINE_MAX
#define LINE_MAX	2048
#endif
#define ON	1
#define OFF	0
#define YES	1
#define NO	0
#define ESC	0x1b
static	const	char	OPTCHAR = '-';
#define MAX(X, Y) (((X) > (Y))? (X): (Y))
#define MIN(X, Y) (((X) < (Y))? (X): (Y))
#define ABS(X)	  (((X) < 0)? (-(X)): (X))

/*	Added headers for conversion from Optimizing C to Turbo C */

#ifdef  BSD

#define cmpmem(x, y, l) bcmp(x, y, l)
#define memcmp(x, y, l) bcmp(x, y, l)
#define movmem(x, y, l) bcopy(x, y, l)
#define memcpy(t, f, l) bcopy(f, t, l)
#define clearmem(x, l)  bzero(x, l)
#define alloc(x)	calloc(x, 1)
#undef	tolower
#undef	toupper
#define srand(x)	srandom(x)
#define rand()		random()
#define remove(x)	unlink(x)

#else

#define cmpmem(x, y, l) memcmp(x, y, l)
#define movmem(x, y, l) memmove(y, x, l)
#define clearmem(x, l)  memset(x, '\0', l)
#define alloc(x)	calloc(x, 1)

#endif

#ifdef news800
#define strchr(s, c)	index(s, c)
#define strrchr(s, c)	rindex(s, c)
#endif

//	if no zlib.h #define USE_ZLIB 0

#ifndef USE_ZLIB
#define USE_ZLIB 1
#endif
#if USE_ZLIB
#include <zlib.h>
static	const	char	gz_ext[] = ".gz";
#if ZLIB_VERNUM < 0x12c0
inline	size_t  fread(void* buf, size_t s, size_t c, gzFile gzfd) {
	int	rv = gzread(gzfd, buf, s * c);
	return (rv < 0? 0: rv / s);
}
inline	size_t fwrite(const void* buf, size_t s, size_t c, gzFile gzfd) {
	return (gzwrite(gzfd, buf, s * c) / s);
}
inline	int fputs(const char* buf, gzFile gzfd) {
	return (gzputs(gzfd, buf));
}
inline	char* fgets(char* buf, int len, gzFile gzfd) {
	return (gzgets(gzfd, buf, len));
}
inline	int fputc(int c, gzFile gzfd) {
	return (gzputc(gzfd, c));
}
inline	int fgetc(gzFile gzfd) {
	return (gzgetc(gzfd));
}
inline	int ungetc(int c, gzFile gzfd) {
	return (gzungetc(c, gzfd));
}
inline	int rewind(gzFile gzfd) {
	return (gzrewind(gzfd));
}
inline	int fseek(gzFile gzfd, long offset, int whence) {
	return (gzseek(gzfd, offset, whence));
}
inline	long ftell(gzFile gzfd) {
	return (gztell(gzfd));
}
inline	int feof(gzFile gzfd) {
	return (gzeof(gzfd));
}
inline	int fclose(gzFile gzfd) {
	return (gzclose(gzfd));
}
inline	int fflush(gzFile gzfd) {
	return (gzflush(gzfd, Z_SYNC_FLUSH));
}
#endif	// ZLIB_VERNUM
#endif	// USE_ZLIB

/*	Others	*/

#ifdef VMS

#define PATHDELM	'/'
#define TEMPDIR		""
#define re_alloc(x) ((x)? realloc(x): malloc(x))

#else

#define PATHDELM	'/'
#define TEMPDIR		"/tmp/"
#define SHELL "csh"

#endif

/*	Error messages  */

extern	const	char*	no_space;
extern	const	char*	no_file;
extern	const	char*	not_found;
extern	const	char*	read_error;
extern	const	char*	write_error;

/*	etc.	*/

extern	const	char*	stddelim;

/*  End of adddef.h  */

#endif
