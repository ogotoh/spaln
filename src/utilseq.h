/*****************************************************************************
*
*	Micelleious utilities subroutines
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

#ifndef _UTILSEQ_H_
#define _UTILSEQ_H_ 1

#include "seq.h"

extern	CHAR	gencode[];
class	EiJuncSeq;

enum class Iefp {IF, IP, IB, EF, EP, EB, CF, CP, CB, CU, GF, NG};

static	const	int	itn_lm = 6;	// intron 5' immune margin
static	const	int	itn_rm = 16;	// intron 3' immune margin
static	const	int	CP_NTERM = 4;
static	const	int	ID_SIZE = 10;
static	const	int	STOP = INT_MAX;
static	const	float	maxtonic = 5.;
static	const	CHAR	max_add_size = 24;	// capacity of additional bytes
static	const	char	text_ext[] = "txt";	// text file
static	const	char	data_ext[] = "dat";	// binary data file
static	const	char	iefp_ext[12][4] = 	// these are also binary
	{"ifq", "ipt", "ifp", "efq", "ept", "efp", 
	 "cfq", "cdp", "cfp", "cus", "gfq", ""};
static	const	char	iefp_tid[12][16] = 
	{"IntronFrqTab", "IntronPotTab", "IntronFpPTab", 
	 "ExonFrqTab", "ExonPotTab", "ExonFpPTab", 
	 "CodeFrqTab", "CodePotTab", "CodeFpPTab", 
	 "CodonUsage", "GenomeFrqTab", ""};

struct	PAT_MATRIX {
	int	rows, cols, offset;
	float*	mtx;
};

struct	ExpectMmm {float min, mean, max;};

struct PatMatHead {
	CHAR	vtype, vsize, nelm, add;	// 0/1:int/float, sizeof, 4/20/23
	int	rows, cols;
};

class CodonUse {
template <typename file_t>
	bool	read_binaryCu(file_t fd);
template <typename file_t>
	bool	read_textCu(file_t fd);
template <typename file_t>
	bool	write_binaryCu(file_t fd);
template <typename file_t>
	bool	write_textCu(file_t fd);
public:
	char	id[ID_SIZE];
	float	usage[64];
	long	ncodons = 0;
	CodonUse() {vclear(id, ID_SIZE);}
	CodonUse(const char* fname, bool from_file = true);
template <typename file_t>
	int	fget(file_t& fd, const char* fn);
	void	to_file(const char* oname, int text);
	void	normalize();
	void	setid(const char* src) {
	    char*	ps = id;
	    for (int i = 0; ++i < ID_SIZE && (*ps++ = *src++); ) ;
	}
};

template <typename file_t>
bool CodonUse::read_binaryCu(file_t fd)
{
	int	nrd = fread(this, sizeof(CodonUse), 1, fd);
	return (nrd == 1);
}

template <typename file_t>
bool CodonUse::read_textCu(file_t fd)
{
	char	str[MAXL];
	char*	ps = str;
	float*	u = usage;
	float*	ut = u + 64;
	do {
	    if (!fgets(str, MAXL, fd)) return (false);
	} while (isBlankLine(str));
	if (*str == '>') {
	    sscanf(str, "%*s [%ld:%*d]", &ncodons);
	    setid(str + 1);
	} else {
	    for ( ; !*ps; ps = cdr(ps))
		*u++ = (float) atof(ps) / 100.;
	}
	while (u < ut && *(ps = fgetw(str, MAXL, fd)))
	    *u++ = (float) atof(ps);
	return (u - usage == 64);
}

template <typename file_t>
int CodonUse::fget(file_t& dmyfd, const char* fn)
{
	if (dmyfd) {
	    fclose(dmyfd);
	    dmyfd = 0;
	}
	ReadFile	fp(fn);
	if (fp.fd) {
	    if (fp.dtype == 1) read_textCu(fp.fd); else
	    if (fp.dtype == 2) read_binaryCu(fp.fd); 
	} else if (fp.gzfd) {
	    if (fp.dtype == 1) read_textCu(fp.gzfd); else
	    if (fp.dtype == 2) read_binaryCu(fp.gzfd);
	} else
	    fatal(not_found, fn);
	return (fp.dtype);
}

template <typename file_t>
bool CodonUse::write_binaryCu(file_t fd)
{
	int	nrd = fwrite(this, sizeof(CodonUse), 1, fd);
	return (nrd == 1);
}

template <typename file_t>
bool CodonUse::write_textCu(file_t ofd)
{
static	const	char	acgt[] = {'A', 'C', 'G', 'T'};
	char	str[MAXL];
	snprintf(str, MAXL, ">%s [%ld:%d]\n", id, ncodons, 64);
	fputs(str, ofd);
	char*	ps = str;
	if (algmode.nsa & 4) {
	    *ps++ = '#';
	    for (int i = 0; i < 4; ++i) {
		snprintf(ps, MAXL - (ps - str), "\t%5c", acgt[i]);
		ps += strlen(ps);
	    }
	    *ps++ = '\n';
	    *ps = '\0';
	    fputs(str, ofd);
	}
	float*	u = usage;
	for (int f = 0; f < 4; ++f) {
	    for (int s = 0; s < 4; ++s) {
		char*	ps = str;
		if (algmode.nsa & 4) {
		    *ps++ = acgt[f];
		    *ps++ = acgt[s];
		}
		for (int t = 0; t < 4; ++t) {
		    snprintf(ps, MAXL - (ps - str), "\t%7.4f", 100. * *u++);
		    ps += strlen(ps);
		}
		*ps++ = '\n';
		*ps = '\0';
		fputs(str, ofd);
	    }
	}
	return (true);
}

class PatMat {
public:
	int	rows = 0, cols = 0, offset = 0;
	float	tonic = 0, min_elem = 0; 
	int	transvers = 0, skip = 0;	// used in reading
	int	nsupport = 0, nalpha = 0, morder = 0;
	ExpectMmm	mmm;		// min, mean, max
	float* mtx = 0;

	PatMat(const PatMat& src);	// copy constructor
	PatMat(FILE* fd, bool binary = false);
	PatMat(const char* fname = 0);
	PatMat(const int r, const int c, const int o = 0, float* m = 0);
	~PatMat() {delete[] mtx;}
template <typename file_t>
	void	read_text(file_t fd);
template <typename file_t>
	void	read_binary(file_t fd, const char* fname);
	float	pwm_score(const Seq* sd, const CHAR* ps, const CHAR* redctab = 0) const;
	float*	calcPatMat(const Seq* sd) const;
	int	order() const {return (morder);}
	int	columns() const {return cols;}
	void	clearmtx() {vclear(mtx, rows * cols);}
	CHAR*	setredctab(const Seq* sd) const;
	void	increment(const Seq* sd, int pos, const CHAR* redctab = 0);
	PatMat&	operator=(const PatMat& src);
};

class ExinPot {
protected:
	int	exin = static_cast<int>(Iefp::IP);
	int	morder = 0;
	int	nphase = 1;	// 1 or 3
	int	ndata = 0;
	int	nsupport = 0;
	int	lm = 0;
	int	rm = 0;
	float	avm = 0;	// mean ipt per 1000 bp
	float	avpot = 0;	// mean self scores per intron
	float	sdpot = 0;	// SD of self scores
	float	avlen = 0;	// average length of nsupport introns
	float	sdm = 0;	// SD of ipt per 100 bp
	float*	data = 0;
	void	count_kmers_1(const Seq* sd, float* fq);
	void	count_kmers_3(const Seq* sd, float* fq, CodonUse* cu);
	void	reform_1(const float& total, float* bkg = 0);
	void	reform_3(const float& total);
	float*	calcScr_1(const Seq* sd, float* scr = 0) const;
	float*	calcScr_3(const Seq* sd, float* scr = 0) const;
public:
	ExinPot(const int& ein, const int& mo, const int& nf = 1, const int it = 0) :
	    exin(ein), morder(mo), nphase(nf), ndata(ipower(4, mo + 1)), 
	    lm(ein / 3? 0: itn_lm), rm(ein / 3? 0: itn_rm)
	{
	    if (it && static_cast<Iefp>(exin) != Iefp::NG)
		from_file(iefp_tid[exin]);
	}
	ExinPot(const char*& fname) {
	    from_file(fname);
	    lm = exin / 3? 0: itn_lm;
	    rm = exin / 3? 0: itn_rm;
	}
	ExinPot(const int& exn) : exin(exn), 
	    lm(exin / 3? 0: itn_lm), rm(exin / 3? 0: itn_rm) {
	    if (static_cast<Iefp>(exin) != Iefp::NG)
		from_file(iefp_tid[exin]);
	    else
		fatal("Invalid exin code (%d) !\n", exin);
	}
	~ExinPot() {delete[] data;}
	bool	ispot() const {
const	    Iefp	iefp = static_cast<Iefp>(exin);
	    return (iefp == Iefp::IP || iefp == Iefp::IB ||
		iefp == Iefp::EP || iefp == Iefp::EB ||
		iefp == Iefp::CP || iefp == Iefp::CB);
	}
	bool	isfpp() const {
const	    Iefp	iefp = static_cast<Iefp>(exin);
	    return (iefp == Iefp::IB || iefp == Iefp::EB || 
		    iefp == Iefp::CB);
	}
	INT	size() const {return (ndata);}
	size_t	dsize() const {
	    return (nphase * ndata + (isfpp()? nphase * ndata: 0));
	}
	bool	from_file(const char* fname);
	float*	getKmers(const char* wdfq, const bool foregrd = true);
	float*	getKmers(EiJuncSeq* eijseq, CodonUse* cu = 0);
	float*	getKmers(int argc, const char** argv, CodonUse* cu = 0);
	void	reform(float* background = 0);
	bool	makeExinPot(const float* gfq);
template <typename file_t>
	bool	read_binary(file_t fd, const char* fname);
template <typename file_t>
	bool	read_text(file_t fd, const char* fname);
template <typename file_t>
	bool	write_binary(file_t fd);
	void	writeBinary(const char* oname, bool gzip = false);
	float*	begin() const {return (data);}
	float*	end() const {return (data + nphase * ndata);}
	float*	fbegin() const {
	    return (data + nphase * (isfpp()? ndata: 0));
	}
	float*	fend() const {
	    return (data + dsize());
	}
	float*	calcScr(const Seq* sd, float* scr = 0) const {
	    return ((nphase == 1)? calcScr_1(sd, scr): calcScr_3(sd, scr));
	}
	VTYPE	intpot(const SGPT6* b5,const  SGPT6* b3) const;
//	VTYPE	intpot(const SGPT2* b5,const  SGPT2* b3) const;
	VTYPE	avrpot(float f = 1.) const {return (VTYPE) (f * avpot);}
	VTYPE	self_score() const {return (avm);}
};

struct EijPat {
	PatMat* pattern5;
	PatMat* pattern3;
	PatMat* patternI;
	PatMat* patternT;
	PatMat*	patternB;
	float   tonic3;
	float   tonic5;
	float   tonicB;
	EijPat(int hh);
	~EijPat();
};

extern	void	setorf(int len, int ic = SILENT);
extern	int	setcodon(int c);
extern	int	codon_id(const CHAR* s, int byte);
extern	void	de_codon_4(CHAR* ncs, int n);
extern	int	toaa(const CHAR* ns);
extern	int	toaa3(const CHAR* ns, int inc);
extern	int	nuc2tron3(const CHAR* ns, int inc);
extern	int	initcodon(int code);
extern	int	initcodon(const char* genspc);
extern	void	mkinvtab();
extern	int	getCodonUsage(const char* fname);
extern	int	setCodonUsage(int gc);
extern	int	fname2exin(const char* fname, int& file_type);

#endif
