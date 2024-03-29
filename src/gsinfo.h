/*****************************************************************************
*
*	Header for scoring conserved exon-intron organizations
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

#ifndef _GSINFO_H_
#define _GSINFO_H_

class	Seq;
struct	Iiinfo;

#if USE_WEIGHT
struct	PFQ	{int pos, num, gps; VTYPE dns;};
#else
struct	PFQ	{int pos, num, gps;};
#endif

extern  VTYPE   SpbFact;
extern	PFQ	pfqend;

struct SigII {
	int	pfqnum = 0;
	int	lstnum = 0;
	int	step = 0;
	PFQ*	pfq = 0;
	int*	lst = 0;
	int**	eijtab = 0;
	int*	lone = 0;

	SigII(const Seq* sd = 0);		// 
	SigII(const SigII& src);	// copy constructor
	SigII(int p, int l, int s);	// known properties
	SigII(const Iiinfo& iif);
	SigII(const Seq** sq, const GAPS** gsrc, FTYPE* wtlst = 0);
template <typename file_t>
	SigII(file_t fd, char* str, FTYPE* wt = 0);	// read from seq file
	SigII(const int* poss, int num, int step);	// list of positions
	~SigII();
	void	rmGapPfq(const GAPS* gg);
	void	relist(int bias);
	void	swaplst(int an, int bn);
	void	renumlst(const int* lst_odr);
	void	putSigII(FILE* fd) const;
	void	mkeijtab(int many);	// bit table of exon-exon junctions
	void	printlones(FILE* fd, int i, int many) const;
	void	printmates(FILE* fd, int i, int many) const;
	float	eij_dist(int i, int j, int* abc = 0) const;
	void	pfqrepos(const RANGE* cr);
	void	locate(PFQ*& wfq, int*& wst, int pos);
	int	to_gene_end(int m, bool rend = false) const;
	void	resetend(int len) {if (pfq) pfq[pfqnum].pos = len * step;}
	int	n_common() const;
#if USE_WEIGHT
	void	rescale_dns(VTYPE f);
	void	reset_dns(FTYPE wt);
#endif
};

template <typename file_t>
SigII::SigII(file_t fd, char* str, FTYPE* wt)
{
	sscanf(str, "%*s %d %d", &pfqnum, &lstnum);
	if (pfqnum == 0) return;
	pfq = new PFQ[pfqnum + 1];
	PFQ*    wfq = pfq;
	int     i = 0;
	while (fgets(str, MAXL, fd)) {
	    if (wordcmp(str, ";b")) break;
	    for (char* ps = cdr(str); ps && *ps; ps = cdr(ps)) {
		wfq->pos = atoi(ps);
		ps = cdr(ps);
		wfq->num = atoi(ps);
		wfq->gps = 0;
#if USE_WEIGHT
		wfq->dns = (lstnum && wt)? 0: wfq->num;
#endif
		++wfq;
		if (++i >= pfqnum) goto readlst;
	    }
	}
	prompt("Insufficient SP boundaries: %d %d\n", i, pfqnum);
readlst:
	*wfq = pfqend;
	if (!lstnum) return;
	lst = new int[lstnum];
	int*    wst = lst;
	i = 0;
#if USE_WEIGHT
	int     n = 0;
#endif
	wfq = pfq;
	while (fgets(str, MAXL, fd)) {
	    if (wordcmp(str, ";m")) break;
	    for (char* ps = cdr(str); ps && *ps; ps = cdr(ps)) {
		int     m = atoi(ps) - 1;
		*wst++ = m;
#if USE_WEIGHT
		if (wt) {
		    wfq->dns += (VTYPE) wt[m];
		    if (++n == wfq->num) {
			++wfq; n = 0;
		    }
		}
#endif
		if (++i >= lstnum) return;
	    }
	}
	prompt("Insufficient SP boundary list: %d %d\n", i, lstnum);
}

class Cip_score {
	Dhash<int, VTYPE>*	cip_hash = 0;
public:
	Cip_score(const Seq* sd);
	~Cip_score() {
	    delete cip_hash;
	}
	VTYPE	cip_score(int c) const {
	    if (!cip_hash) return (0);
	    KVpair<int, VTYPE>* kv = cip_hash->find(c);
	    return (kv? kv->val: 0);
	}
};

class PfqItr {
friend	struct	SigII;
friend	struct	Iiinfo;
friend	class	Cip_score;
private:
	int	pfqnum;
	int	lstnum;
	int	step;
	PFQ*	pfq;
	int*	lst;
	PFQ*	wfq = 0;
	int*	wst = 0;
	PFQ*	tfq;
#if USE_WEIGHT
	FTYPE*	weight;
#endif
public:
	PfqItr& operator++() {
	    if (wfq < tfq) {
		if (wst) wst += wfq->num;
		++wfq;
	    }
	    return (*this);
	}
	PfqItr& operator--() {
	    if (wfq > pfq) {
		if (wst) wst -= wfq->num;
		 --wfq;
	    }
	    return (*this);
	}
	PfqItr& operator+=(int n) {
	    while (n-- > 0 && wfq < tfq) {
		if (wst) wst += wfq->num;
		++wfq;
	    }
	    return (*this);
	}
	PfqItr& operator-=(int n) {
	    while (n-- > 0 && wfq > pfq) {
		if (wst) wst -= wfq->num;
		 --wfq;
	    }
	    return (*this);
	}
	bool operator==(int n) const {	// exact match
	    return (wfq && (n == wfq->pos));
	}
	bool operator<(int n) const {
	    return (wfq && (wfq->pos < n * step));
	}
	bool operator<=(int n) const {
	    return (wfq && (wfq->pos < ++n * step));
	}
	bool operator&&(int m) const {	// codon match
	    if (SpbFact == 0 || !wfq) return (false);
	    if (step == 1) return (wfq->pos == m);
	    m *= step;
	    return ((m <= wfq->pos) && (wfq->pos < m + step));
	}
	bool eq(int m) const {		// codon +-1 match
	    if (SpbFact == 0 || !wfq) return (false);
	    if (step == 1) return (wfq->pos == m);
	    m = step * m - 1;
	    return ((m <= wfq->pos) && (wfq->pos < m + step));
	}
	bool lt(int m) const {
	    return (wfq && (wfq->pos + 1 < m * step));
	}
	bool le(int m) const {
	    return (wfq && (wfq->pos + 1 < ++m * step));
	}
	int size() const {return pfqnum;}
	bool end() const {return (wfq >= tfq);}
	PfqItr& reset(int n = 0) {
	    wfq = pfq; wst = lst; n *= step;
	    while (wfq < tfq && wfq->pos < n) {
		if (wst) wst += wfq->num;
		++wfq;
	    }
	    return (*this);
	}
#if USE_WEIGHT
	PfqItr(SigII& sgi, int n = 0, FTYPE* wt = 0);
	VTYPE match_score(int n) { return (n == wfq->pos? SpbFact * wfq->dns: 0); }
#else
	PfqItr(SigII& sgi, int n = 0);
	VTYPE match_score(int n) { return (n == wfq->pos? SpbFact * wfq->num: 0); }
#endif
	PfqItr(const Seq* sd, int n = 0);
	VTYPE match_score(PfqItr& bpi, bool all_phase = true) const {
	    if (step != 1) {
		if ((wfq->pos - bpi.wfq->pos) % step) return (0);
		if (!all_phase && wfq->pos % step) return (0);
	    }
#if USE_WEIGHT
	    return (SpbFact * wfq->dns * bpi.wfq->dns);
#else
	    return (SpbFact * wfq->num * bpi.wfq->num);
#endif
	}
};

struct	Iiinfo {
const	Seq*&	a;
const	Seq*&	b;
const	int&	amany;
const	int	step;
	PfqItr*	api = 0;
	PfqItr*	bpi = 0;
	SigII*	sgi = 0;
	PfqItr*	cpi = 0;
	int	agap = 0;
	int	bgap = 0;
	Iiinfo(const Seq* seqs[], int m, int n, bool save = false);
	~Iiinfo() {delete api; delete bpi; delete cpi; delete sgi;}
	VTYPE	StoreIIinfo(int m, int n);
	SigII*	finalize(int len);
};

struct EISCR {
	int	left;
	int	right;
	int	rleft;
	int	rright;
	int	mch;
	int	mmc;
	int	gap;
	int	unp;
	int	mch5;
	int	mmc5;
	int	gap5;
	int	unp5;
	int	mch3;
	int	mmc3;
	int	gap3;
	int	unp3;
	int	phs;
	VTYPE	escr;
	VTYPE	iscr;
	VTYPE	sig3;
	VTYPE	sig5;
};

struct CIGAR {
	char	ope;
	int	len;
};

struct VULGAR {
	char	ope;
	int	alen;
	int	blen;
};

class Cigar {
	int	num = 0;
	Mfile*	cmfd = 0;
public:
	CIGAR*	rec = 0;
	Cigar() {cmfd = new Mfile(sizeof(CIGAR));}
	~Cigar() {delete[] rec;}
	void	push(char mk, int len) {
	    CIGAR	buf = {mk, len};
	    cmfd->write(&buf);
	}
	void	flush() {
	    num = cmfd->size();
	    rec = num? (CIGAR*) cmfd->flush(): 0;
	    delete cmfd; cmfd = 0;
	}
	int	size() const {return num;}
};

class Vulgar {
	int	num = 0;
	Mfile*	vmfd = 0;
public:
	VULGAR*	rec;
	Vulgar() {vmfd = new Mfile(sizeof(VULGAR));}
	~Vulgar() {delete[] rec;}
	void	push(char mk, int alen, int blen) {
	    VULGAR	buf = {mk, alen, blen};
	    vmfd->write(&buf);
	}
	void	postproc();
	void	flush() {
	    num = vmfd->size();
	    rec = num? (VULGAR*) vmfd->flush(): 0;
	    delete vmfd; vmfd = 0;
	}
	int	size() const {return num;}
};

class Eijnc {
	int	num = 0;
	Mfile*	emfd = 0;
	EISCR*	rec = 0;
const	int	qsize = 0;
	FSTAT*	fstque = 0;
	int	q = 0;
public:
	Eijnc(bool que = false);
	~Eijnc() {delete[] rec; delete[] fstque;}
	void	push(EISCR* eis) {emfd->write(eis);}
	void	flush() {
	    num = emfd->size();
	    rec = num? (EISCR*) emfd->flush(): 0;
	    delete emfd; emfd = 0;
	}
	EISCR*	begin() const {return(rec);}
	int	size() const {return num;}
	int	genleft() const {return rec? rec->left: 0;}
	int	genright() const {return rec? rec[num - 2].right: 0;}
	int	refleft() const {return rec? rec->rleft: 0;}
	int	refright() const {return rec? rec[num - 2].rright: 0;}
	void	store(EISCR& rbuf, FSTAT& now, FSTAT& prv, bool nearjnc);
	void	shift(EISCR& rbuf, FSTAT& now, bool nearjnc);
	void	unshift() {if (q) --q; else q = alprm2.jneibr - 1;}
};

struct Samfmt : public Cigar {
	int	flag = 0;
	int	pos = 0;
	int	mapq = 0;
const	char*	rnex = 0;
	int	pnext = 0;
	int	tlen = 0;
	char*	qual = 0;
mutable	int	left = 0;
mutable	int	right = 0;
	Samfmt() : Cigar() {}
	~Samfmt() {}
};

class Gsinfo {
mutable	FILE*	fd;
	void	repalninf0(const SKL* skl, Seq* seqs[]) const;
	void	repalninf1(const SKL* skl, Seq* seqs[]) const;
	void	repalninf2(const SKL* skl, Seq* seqs[]) const;
	void	repalninf3(const SKL* skl, Seq* seqs[]) const;
	void	repalninf4(const SKL* skl, Seq* seqs[]) const;
	void	repalninf5(const SKL* skl, Seq* seqs[]) const;
	void	repalninf6(const SKL* skl, Seq* seqs[]) const;
	void	Gff3Form(const Seq* gene, const Seq* qry) const;
	void	Gff3PWA(const Seq* gene, const Seq* qry) const;
	void	BedForm(const Seq* gene, const Seq* qry) const;
	void	ExonForm(const Seq* gene, const Seq* qry, int mode) const;
	void	IntronForm(const Seq* gene, const Seq* qry) const;
	void	CigarForm(const Seq* gene, const Seq* qry) const;
	void	VulgarForm(const Seq* gene, const Seq* qry) const;
	void	SamForm(const Seq* gene, Seq* qry) const;
public:
	int	end_error_thr;
	VTYPE	scr = 0;
	VTYPE	rscr = 0;
	SKL*	skl =  0;
	FSTAT	fstat;
	int	noeij = 0;
	RANGE*	CDSrng = 0;
	Eijnc*	eijnc = 0;
	Cigar*	cigar = 0;
	Vulgar*	vlgar = 0;
	Samfmt*	samfm = 0;
	SigII*	sigII = 0;
const	char*	prefix = 0;
	bool	intronless() const;
	Gsinfo(SKL* s = 0) : end_error_thr(int(alprm2.jneibr * 0.8)), skl(s) {
	    vclear(&fstat);
	}
	~Gsinfo() {
	    delete[] skl;
	    delete[] CDSrng;
	    delete eijnc;
	    delete sigII;
	    delete cigar;
	    delete vlgar;
	    delete samfm;
	}
	void	clear();
	void	SaveGsInfo(Iiinfo* iif, int len) {
	    delete sigII;
	    sigII = iif? iif->finalize(len): 0;
	}
	RANGE*	eiscr2rng();
	RANGE*	eiscrunfold(GAPS* gp);
	RANGE*	querygs(const Seq* qry) const;
	void	BoundaryInf(const Seq* sd) const;	// write to out_fd
	void	BoundarySeq(const Seq* sd) const;	// write to out_fd
	int	center(int k) const;	// center position
	void	repalninf(Seq* seqs[], int mode, FILE* fd = 0) const;
	void	printgene(Seq** seqs, int mode, FILE* _fd = 0) const;
	int	print2(Seq* seqs[], const GAPS** gaps, 
		int nbr, int ttl, int skip, FILE* _fd = 0) const;
	void	setprefix(const char* px) {prefix = px;}
};

extern	SigII*  copySigII(const SigII* src);
extern	SigII*	extSigII(const Seq* sorc, const int* which, 
		const FTYPE& nfact = 1, const bool& renum_lst = false);
extern	void	cutSigII(Seq* dstseq, const Seq* srcseq);
extern	void	catSigII(Seq* dstseq, const Seq* srcseq, int bias);
extern	FTYPE*	eijdmx(Seq* sd);
extern	void	fouteijdmx(FILE* fd, const Seq* sd, bool dmx = true);
extern	void	fouteij(FILE* fd, Seq* sd);
extern	VTYPE	spb_fact();
extern	void	unfoldPfq(PFQ* pfq, int num, const GAPS* gg, int step);

inline	bool	neoeij(const EISCR* eij) {return (eij->left != endrng.left);}
inline	bool	neopfq(const PFQ* pfq) {return (pfq && pfq->num);}
inline	bool	use_spb() {return (SpbFact > 0);}
inline	bool	isEIJ(int phs) {return (algmode.lsg && phs > -2);}
inline	bool	isPhs1(int phs) {return (algmode.lsg && phs > 0);}

#endif
