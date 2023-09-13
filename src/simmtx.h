/*****************************************************************************
*
*	Header to  Similarity Matrix
*	for proteins:	Mutation Data Matrix
*	for DNA/RNA:	Simple match/mismatch
*	for protein/DNA	Special aa/codon(tron) matching score
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

#ifndef SIMMTX_H
#define	SIMMTX_H	1

enum ComPmt {DxD, PxT, TxP, PxP, TxT, CxC, AxA};

static	const int	max_simmtxes = 3;
extern	const char*	mdm_file[max_simmtxes];

struct	DefPrm {float u, v, n, b; int p;};

class Simmtx {
	DefPrm*	param;
	float	nrmlf, avtrc;
	VTYPE	drange, minscr;
	VTYPE** SquareMtx();
public:
	int	dim;
	int	rows;
	int	cols;
	ComPmt	AxB;
const	char*	simfile;
	VTYPE**	mtx;
	VTYPE*	simunp;
	Simmtx(ComPmt ab, DefPrm* dp);
	Simmtx(ComPmt ab, const char* fname, DefPrm* dp);
	~Simmtx() {delete[] *mtx; delete[] mtx;}
	int	simgrade(int aa, int bb) const;
	float	Nrmlf() const {return nrmlf;}
	float	AvTrc() const {return avtrc;}
	void	Nmtx();
	void	Nmtx(const char* fname);
	void	Pmtx();
	void	Pmtx(const char* fname);
	void	Hmtx(ComPmt compmt, Simmtx* pm);
	void	fparam(FILE* fd) const;
	double	avrmatch(INT* ConvTab) const;
	void	printmtx(FILE* fd) const;
};

class Simmtxes {
	double*	mdmcomp;
	int	current;
	Simmtx*	storedSimmtx[max_simmtxes];
public:
	double*	get_mdmcmp();
	Simmtxes() : mdmcomp(0), current(0) {
	    vclear(storedSimmtx, max_simmtxes);
	}
	~Simmtxes();
	Simmtx* getsimmtx(int c = -1) {
	    if (0 <= c && c < max_simmtxes) current = c;
	    return (storedSimmtx[current]);
	}
friend	void	setSimmtxes(ComPmt ab, bool mdm, DefPrm* dp);
friend	void	resetSimmtxes(bool rmcomps);
};

extern	Simmtxes	simmtxes;
extern	int	minmax;

extern	void	setlsegs(int ls);
extern	void	setalprm();
extern	void	readalprm(int& argc, const char**& argv, int oc = 1);
extern	int	dvp2pam(double x);
extern	void	optimize(int gl, int mnmx);
extern	void	setdefNprm(float n, float u, float v, float b = 0., int c = 0);
extern	void	setdefPprm(int p, float u, float v, float b = 0., int c = 0);
extern	void	setpam(int pam, int scnd);
extern	void	setNpam(int q, float v);
extern	int	(*simgrade)(int aa, int bb);
extern	int	getpam(int scnd = 0);
extern	float	getsmn(int q = 4);
extern	void	setSimmtxes(ComPmt ab, bool mdm = false, DefPrm* dp = 0);
extern	void	resetSimmtxes(bool rmcomps = false);

inline	double*	get_mdmcmp() {return simmtxes.get_mdmcmp();}
inline	Simmtx* getSimmtx(int c = -1) {return simmtxes.getsimmtx(c);}

#endif
