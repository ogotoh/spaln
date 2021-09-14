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

static	const	INT	MaxBitPat = 16;
static	const	INT	BAD_RES = 0xff;
static	const	INT	DefAaKtuple = 5;
static	const	INT	DefNucKtuple = 11;
static	const	INT	SupSiteNo = UINT_MAX;
static	const	INT	BadWord = UINT_MAX;
static	const	CHAR	ntconv[26] =
	{0,4,1,4,7,7,2,4,7,7,4,7,4,4,7,4,7,4,4,3,3,4,4,7,4,7};
//       a b c d e f g h i j k l m n o p q r s t u v w x y z

extern	const	char*	DefConvPat[];
extern	const	char*	DefBitPat[];
extern	const	int	DefConvPatNo[];

enum {FourN, Dayh6, SEB6};

class ReducWord {
protected:
	int	molc;
	bool	master;
public:
	INT	ConvTabSize;
	INT	Nalpha;
	CHAR*	aConvTab;	// text -> reduced alphabet
const	CHAR*	a2r;		// aa -> reduced aa
	CHAR*	g2r = 0;	// codon -> raa
	CHAR*	g2r_c = 0;	// reverse complementary codon -> raa
	CHAR*	iConvTab;	// internal code -> reduced alphabet
	ReducWord(const Seq* sd, INT elms = 0, const char* ap = 0);
	ReducWord(const ReducWord& src) {
	    *this = src;
	    master = false;
	}
	~ReducWord() {
	    if (master) {
		delete[] aConvTab;
		delete[] iConvTab;
		delete[] g2r;
	    }
	}
};

class Bitpat {			// bit pattern for instant word generation
protected:
	bool	master;
public:
	int*	exam;
	int	weight;
	int	width;
	int	wshift;
	Bitpat(INT npat, const char* spat = 0);
	Bitpat(const Bitpat& src) {
	    *this = src;
	    master = false;
	}
	~Bitpat() {if (master) delete[] exam;}
	bool	spaced() const {return (width > weight);}
};

class Bitpat_wq : public Bitpat {	// bit pattern with queue 
const	INT	nalpha;
const	INT	msb;
const	int	nframe;
const	bool	reverse;
const	int	noq;
const	int	qsize;
	INT*	queue = 0;
	INT*	fstat = 0;
	CHAR*	qp = 0;
public:
const	INT	TabSize;
	Bitpat_wq(INT elms, int nf, bool rvs, INT npat, const char* spat = 0);
	Bitpat_wq(const Bitpat_wq& src);
	~Bitpat_wq() {delete[] queue; delete[] fstat; delete[] qp;}
	void	clear();
	bool	good(INT c) const {return (c < nalpha);}
	void	flaw(int f = 0);
	bool	flawless(int f = 0) const {return (fstat[f] == 0);}
	INT	word(INT c, int f = 0);
};

struct WordState {
	int	p;
	SHORT	ss;
	SHORT	ss6[6];
	WordState(int nss) : p(0), ss(0)
	    {vclear(ss6, 6);}
	~WordState() {}
};

class WordTab : public ReducWord {
protected:
const	INT     Nbitpat;
const	INT	kmer;
const	INT	BitPat;
const	INT	BitPat2;
const	INT	Nshift;
const	INT	npos;
	int	wq_size;
const	int	nfrm;
	INT	toomany;
	int	p = 0;
	Bitpat_wq**	bpp = 0;
	SHORT	cc[6];
	INT	xx[3];
	INT**	heads = 0;
	INT**	wposs = 0;
	SHORT	ss6[6];
	SHORT&	ss = ss6[0];
	INT*	word_at_0s = 0;
	INT*	head = 0;
	INT*	wpos = 0;
	int	w_qp = 0;
	INT**	w_queue = 0;
public:
	WordTab(const Seq* sd, INT tpl, INT nsft = 1, INT elms = 0, 
	    const char* ap = 0, INT bp = 0, INT bp2 = 0, INT nbt = 1, 
	    INT np = 0, INT afact = 0, INT minorf = 0);
	WordTab(const WordTab& src);
	~WordTab();
	void    reset(int mode = 1);
	void	save(INT w, INT i, INT k = 0);
	void	c2w(INT uc, int i = -11);
	void    c2w6(INT uc, int i = -1);
	void	c2w6_pp(int bpos, int n);
template <typename var_t>
	void	list2lut(var_t* header, var_t* position) const {
	    var_t   m = 0;
	    for (INT i = 0; i < bpp[0]->TabSize; ++i) {
		header[i] = m;
		var_t       n = 0;
		for (INT p = head[i]; p; p = wpos[p]) ++n;
		if (i == word_at_0s[0] || i == word_at_0s[1]) ++n;
		if (n == 0 || n > toomany) continue;	// no or too many
		m += n; n = m;
		for (INT p = head[i]; p; p = wpos[p])
		    position[--n] = (var_t) (p % npos);
		if (i == word_at_0s[0] || i == word_at_0s[1])
		    position[--n] = 0;
	    }
	    header[bpp[0]->TabSize] = m;
	}
	int	max_width() const {
	    int	w = 0;
	    for (INT k = 0; k < Nbitpat; ++k)
		if (bpp[k]->width > w) w = bpp[k]->width;
	    return (w);
	}
	void	snapshot(WordState* ws) const {
	    ws->p = p; ws->ss = ss; vcopy(ws->ss6, ss6, 6);
	}
};

extern	INT	bpcompress(const char* sp, INT* w = 0);
inline	INT	bitmask(INT k) {return ((1 << k) - 1);}

#endif
