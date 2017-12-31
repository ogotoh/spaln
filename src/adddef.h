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
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/

#ifndef  _ADDDEF_
#define  _ADDDEF_


#define UPTR	void*

#include <limits.h>
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
#define MAXL	_POSIX2_LINE_MAX
#define NAMSIZ	MAXL
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

extern const	char*	no_space;
extern const	char*	no_file;
extern const	char*	not_found;

/*	etc.	*/

extern const	char*	stddelim;   /* = " \t\n\r\f"    */

/*  End of adddef.h  */

#endif
