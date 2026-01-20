/*****************************************************************************
*
*	Collection of commonly used headers
*
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

#ifndef  _ADDDEF_
#define  _ADDDEF_

#define UPTR	void*

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <zlib.h>

typedef unsigned char CHAR;
typedef	unsigned short SHORT;
typedef unsigned int  INT;
typedef unsigned long LONG;

typedef int (*CMPF)(const UPTR, const UPTR);

static	const	int	ERROR = -1;
static	const	int	OK = 0;
static	const	int	ABORT = INT_MIN;
static	const	int	MAXL = 256;
static	const	int	ON = 1;
static	const	int	OFF = 0;
static	const	int	YES = 1;
static	const	int	NO = 0;
static	const	char	ESC = 0x1b;
static	const	char	OPTCHAR = '-';
static	const	char	read_error[] = "\'%s\' read error!\n";
static	const	char	write_error[] = "\'%s\' write error!\n";
#ifndef LINE_MAX
#define LINE_MAX	2048
#endif
#define MAX(X, Y) (((X) > (Y))? (X): (Y))
#define MIN(X, Y) (((X) < (Y))? (X): (Y))

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

/*	Others	*/

#define PATHDELM	'/'
#define TEMPDIR		"/tmp/"

inline	void fatal(const char* format,...)
{
	va_list	args;

	va_start(args, format);
	(void) vfprintf(stderr, format, args);
	putc('\n', stderr);
	va_end(args);
	exit(1);
}

#endif	// End of adddef.h

