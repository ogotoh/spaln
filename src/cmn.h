/*****************************************************************************
*
*	Collection of headers commonly used for sequence analysis
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

#ifndef _CMN_H_
#define _CMN_H_

#include "stdtype.h"
#include <math.h>

#ifndef FVAL
#define	FVAL	0
#endif

#ifndef DVAL
#define	DVAL	0
#endif

#ifndef LVAL
#define	LVAL	0
#endif

#if DVAL
typedef	double	FTYPE;
#else
typedef	float	FTYPE;
#endif	// DVAL

#define	USE_WEIGHT	FVAL

#if FVAL

typedef	FTYPE	VTYPE;
typedef	FTYPE	STYPE;
typedef	float	PVTYPE;	// printf
static	const	FTYPE	VABORT = -1.e127;
static	const	FTYPE	fepsilon = 1.e-7;
static	const	VTYPE	fInfinit = FLT_MAX;
#if DVAL
static	const	VTYPE	NEVSEL = -(DBL_MAX / 16 * 7);
#else
static	const	VTYPE	NEVSEL = -(FLT_MAX / 16 * 7);
#endif
inline	bool gt(FTYPE a, FTYPE b) {return (a >  (b + fepsilon * MAX(1., fabs(b))));}
inline	bool ge(FTYPE a, FTYPE b) {return (a >= (b - fepsilon * MAX(1., fabs(b))));}
inline	bool lt(FTYPE a, FTYPE b) {return (a <  (b - fepsilon * MAX(1., fabs(b))));}
inline	bool le(FTYPE a, FTYPE b) {return (a <= (b + fepsilon * MAX(1., fabs(b))));}

#else	// !FVAL

#define	VABORT	ABORT
#if LVAL
typedef	long 	VTYPE;
typedef	int	STYPE;
static	const	VTYPE	NEVSEL = (LONG_MIN / 16 * 7);
static	const	VTYPE	fInfinit = LONG_MAX;
#else	// !LVAL
typedef	int	VTYPE;
typedef	short	STYPE;
static	const	VTYPE	NEVSEL = (INT_MIN / 16 * 7);
static	const	VTYPE	fInfinit = INT_MAX;
#endif	// LVAL
typedef	int	PVTYPE;
static	const	FTYPE	fepsilon = 0;
inline	bool gt(int a, int b) {return (a >  b);}
inline	bool ge(int a, int b) {return (a >= b);}
inline	bool lt(int a, int b) {return (a <  b);}
inline	bool le(int a, int b) {return (a <= b);}

#endif	// FVAL

#if LVAL
typedef	long	IVTYPE;
#else
typedef	int	IVTYPE;
#endif	// LVAL

static	const	int	EOS = INT_MIN;
static	const	int	POS_INT = (INT_MAX / 8 * 7);
static	const	int	NEG_INT = (INT_MIN / 8 * 7);
static	const	int	NSENT = 2;
static	const	int	IGNORE = EOF - 1;

enum InputMode {IM_NONE, IM_SNGL, IM_ALTR, IM_MULT, IM_EVRY, IM_FvsO, IM_OvsL,
	IM_GRUP, IM_SELF, IM_PARA, IM_ADON, IM_TREE, IM_SCAN, IM_RNDM};
enum InSt {IS_OK, IS_ERR, IS_END};
enum {NM='>', CM='^', RV='~', CR='<'};
enum {NORMAL, COMPLE, REVERS, COMREV}; 
enum {UNKNOWN, PROTEIN, DNA, RNA, TRON, GENOME};
enum {LINEAR, CIRCLE};
enum {RAWSEQ, VECTOR, PROVEC, VECPRO};
enum {Ade, Cyt, Gua, Thy};
enum {___,_,A,C,M,G,R,S,V,T,U=9,W,Y,H,K,D,B,N,NTS=16,Z};
//      0 1 2 3 4 5 6 7 8 9 9  10111213141516  16   17
enum {NIL,UNP,AMB,ALA,ARG,ASN,ASP,CYS,GLN,GLU,GLY,HIS,ILE,LEU,LYS,MET,PHE,
//      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16
PRO,SER,THR,TRP,TYR,VAL,ASX,SER2=23,GLX,SEC=24,TRM2=24,AAS=24,TRM,ZZZ};
//17 18  19  20  21  22  23   23     24          24     24     25  26
enum {MINIMUM, MAXIMUM};

inline	int	elem(int i, int j) {
	return (i < j? j*(j-1)/2+i: i*(i-1)/2+j);
}
inline	int	ncomb(int n) {
	return (n * (n - 1) / 2);
}

struct	RANGE	{int left, right;};
struct	ORF	{int pos, len, frm;};
struct	SKL	{int m, n;};
struct	SKLP	{int m, n, p;};
struct	VSKLP	{VTYPE val; int m, n, p;};
struct	WINDOW	{int lw, up, width;};
struct 	FSTAT	{FTYPE mch, mmc, gap, unp, val;};
struct	GAPS	{int gps, gln;};

static	const	char	_COMM = ';';
static	const	char	_LCOMM = '#';
static	const	RANGE	zerorng = {0, 0};
static	const	RANGE	endrng = {INT_MAX, INT_MAX};

inline	bool	neorng(const RANGE* r) {return(r->left != INT_MAX);}
inline	bool	emptyrng(const RANGE* r) {return(r->left == NSENT);}
inline	int	sizerng(const RANGE* r) {return (r->left);}
inline	RANGE*	fistrng(RANGE* r) {return (r + 1);}
inline	RANGE*	lastrng(RANGE* r) {return (r + r->left - NSENT);}
inline	const	RANGE*	fistrng(const RANGE* r) {return (r + 1);}
inline	const	RANGE*	lastrng(const RANGE* r) {return (r + r->left - NSENT);}

extern	FILE*	out_fd;
extern	int	minmax;

#endif	// _CMN_H_
