/*****************************************************************************
*
*	Quickly calculate distance between sequences 
*	Using oligomer compositions
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-4-7 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#ifndef	_BITPAT_
#define	_BITPAT_

static	const	INT	MaxBitPat = 15;
static	const	CHAR	BAD_RES = 0xff;
static	const	INT	DefAaKtuple = 5;
static	const	INT	DefNucKtuple = 11;
extern	const	char*	DefConvPat[];
extern	const	char*	DefBitPat[];
extern	const	int	DefConvPatNo[];

enum {FourN, Dayh6, SEB6};

class ReducWord {
protected:
	int	molc;
public:
	INT	ConvTabSize;
	INT	Nalpha;
	INT*	aConvTab;	// text -> reduced alphabet
	INT*	iConvTab;	// internal code -> reduced alphabet
	ReducWord(Seq* sd, INT elms = 0, const char* ap = 0);
	~ReducWord() {delete[] aConvTab; delete[] iConvTab;}
};

struct Bitpat {			// bit pattern for instant word generation
	INT	Nalpha;
	INT	TabSize;
	int	weight;
	int	width;
	int	wshift;
	int*	exam;
	bool	good(INT c) {return (c < Nalpha);}
	Bitpat(INT elms, INT npat, const char* spat = 0);
	~Bitpat() {delete[] exam;}
};

class Bitpat_wq : public Bitpat {	// bit pattern with queue 
	INT*	queue;
	INT*	fstat;
	int*	qp;
	INT	msb;
	int	noq;
	int	qsize;
	int	nframe;
	bool	reverse;
public:
	Bitpat_wq(INT elms, int nf, bool rvs, INT npat, const char* spat = 0);
	~Bitpat_wq() {delete[] queue; delete[] fstat; delete[] qp;}
	void	clear();
	void	flaw(int f = 0);
	bool	flawless(int f = 0) {return (fstat[f] == 0);}
	INT	word(INT c, int f = 0);
};

extern	INT	bpcompress(const char* sp, INT& w);
inline	INT	bitmask(INT k) {return ((1 << k) - 1);}

#endif
