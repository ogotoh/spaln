/*****************************************************************************
*
*	Fetch  Similarity Matrix
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

#include "aln.h"
#include "mdm.h"
#include "simmtx.h"
#include "wln.h"
#include <math.h>

#define level(i, j) (4 - 9 * countbit((i) & (j))/countbit(i)/countbit(j)/2)

#define USE_GLX 0

/* arbitrarily assign termination code vs others */

int	minmax = MAXIMUM;

Simmtxes	simmtxes;

const	char*	mdm_file[max_simmtxes] = {0, 0, 0};
static	const	char*	Badpam = "Illegal pam = %d !\n";
static	const	char*	strscale = "scale";

//		   u,      v,   u0, u1, v0, tgapf, thr, sclase, maxsp, gamma, k1, ls, sh, ubh, mtx_no
ALPRM	alprm = {FQUERY, FQUERY, 0., 0.6, 0, 1.0, 35., 1., 8., 0.5, 7, 1, 100, 0, 0};
//ALPRM	alprm = {FQUERY, FQUERY, 0., 0.6, 0, 1.0, 35., 1., 8., 0.5, 7, 1, 100, 2, 0};
//		   x,   y,      z,   o, m, bti, spb, Z, sss, jneibr termk1 desert
ALPRM2	alprm2 = {30., FQUERY, FQUERY, 30., 9., 8., 20., 0., -1., 10, 45, 150};
//		scnd, hydr, hpmt, hpwing, no_angle
ALPRM3	alprm3 = {0., 0., 0., 3, 0};

BPPRM	bpprm = {0., 1., 1., 100};

static	int	glocal = GLOBAL;
static	float	smn[] = {2, 1., 0., -1., FQUERY};
static	DefPrm	defNprm[max_simmtxes+1] = 
{{3., 8., -6., 0., 1}, {2., 6., -4., 0., 1}, {2., 4., -2., 0., 1}, {0}};
static	DefPrm	defPprm[max_simmtxes+1] = 
{{4., 10., 0., 0., 100}, {2., 9., 0., 0., 150}, {2., 9., 0., 0., 250}, {0}};

void optimize(int gl, int mnmx)
{
	glocal = gl;
	minmax = mnmx;
}

//	Use of Mutation-Data-Matrix

int dvp2pam(double x)
{
	if (x >= 1.) return (MAXPAM);
	x = (x <= 0.7)? 1. - (0.987151 + 0.220560*x)*x:
	    -1.260444 + (8.603930  - (13.869219 - 6.521836*x)*x)*x;
	if (x <= 0.) return (MAXPAM);
	int	pam = (int) (-100. * log(x)) + PAMSTEP / 2;
	pam = (pam + PAMSTEP - 1) / PAMSTEP * PAMSTEP;
	return ((pam > MAXPAM)? MAXPAM: pam);
}

VTYPE** Simmtx::SquareMtx()
{
	mtx = new VTYPE*[dim];
	VTYPE*	wmtx = new VTYPE[dim * dim];
	for (int i = 0; i < dim; ++i, wmtx += dim)
	    mtx[i] = wmtx;
	return mtx;
}

int Simmtx::simgrade(int aa, int bb) const
{
	if (IsGap(aa) || IsGap(bb)) return (0);
	if (AxB == DxD)
	    return (int) (4 * (mtx[aa][bb] - minscr) / drange);
	int 	cc = (int) ((25 * mtx[aa][bb] + drange  - 1) / drange);
	cc = max(cc, 0);
	return (cc);
}

Simmtx::Simmtx(ComPmt ab, DefPrm* dp) 
	: param(dp), AxB(ab), simfile(0)
{
	Simmtx*	pm = 0;
	switch (ab) {
	    case DxD: Nmtx(); break;
	    case PxP: Pmtx(); break;
	    case PxT: case TxP: case TxT:
		pm = new Simmtx(PxP, dp);
		Hmtx(ab, pm);
		delete pm; break;
	    default:
		fatal("Simmtx %d is not supported currently!\n", (int) ab);
	}
	simunp = mtx[gap_code];
}

Simmtx::Simmtx(ComPmt ab, const char* fname, DefPrm* dp) 
	: param(dp), AxB(ab), simfile(fname)
{
	Simmtx*	pm = 0;
	param->p = GivenMat;
	switch (ab) {
	    case DxD: Nmtx(fname); break;
	    case PxP: Pmtx(fname); break;
	    case PxT: case TxP: case TxT:
		pm = new Simmtx(PxP, fname, dp);
		Hmtx(ab, pm);
		delete pm; break;
	    default:
		fatal("Simmtx %d is not supported currently!\n", (int) ab);
	}
	simunp = mtx[gap_code];
}

Simmtxes::~Simmtxes()
{
	for (int i = 0; i < max_simmtxes; ++i) {
	    delete storedSimmtx[i];
	    storedSimmtx[i] = 0;
	}
	delete[] mdmcomp; mdmcomp = 0;
}

void Simmtx::Nmtx()
{
	dim = rows = NSIMD;
	mtx = SquareMtx();
	VTYPE	ntsunp = (VTYPE) -(alprm.scale * param->u);
	setNpam(4, smn[4]);
	for (int i = 1; i < NTS; ++i) {
	    int	ii = i + gap_code;
	    for (int j = 1; j < i; ++j) {
		int	jj = j + gap_code;
		mtx[ii][jj] = mtx[jj][ii] = 
		    (VTYPE) (alprm.scale * smn[level(i, j)]);
	    }
	    mtx[ii][ii] = (VTYPE) (alprm.scale * smn[level(i, i)]);
	    mtx[_][ii] = mtx[ii][_] = ntsunp;
	    mtx[NIL][ii] = mtx[ii][NIL] = 0;
	}
	mtx[_][_] = mtx[NIL][NIL] = 
	mtx[_][NIL] = mtx[NIL][_] = 0;
	avtrc = mtx[A][A] + mtx[C][C] + mtx[G][G] + mtx[T][T];
	nrmlf = avtrc /= 4;
	minscr = mtx[A][C];
	drange = mtx[A][A] - minscr;
}

void Simmtx::Nmtx(const char* fname)
{
	FILE*	fd = ftable.fopen(fname, BLASTMAT, DEF_MAT_PATH);
	if (!fd) fatal("Matrix File %s was not found!\n", fname);
	dim = rows = NSIMD;
	mtx = SquareMtx();
	int	aref[26];
	int	xref[26];
	int	nrow[26];
	float	scale = 1.;
	VTYPE	maxs = INT_MIN;
	VTYPE	ntsunp = (VTYPE) -(alprm.scale * param->u);
static	const	int	nbit[16] = {0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4};

	char	str[MAXL];
	*str = _LCOMM;
	while (*str == _LCOMM) {
	    if (!fgets(str, MAXL, fd))
		fatal("Matrix File %s was incomplete!\n", fname);
	    char*	eql = strchr(str + 1, '=');
	    if (eql) {
		*eql = '\0';
		if (!strcmp(str + 1, strscale)) {
		    scale = atof(eql + 1);
		    if (scale <= 0.) scale = 1.;
		}
	    }
	}

	for (int i = 0; i < 26; ++i) {
	    aref[i] = nrow[i] = -1;
	    xref[i] = nccode[i];
	}
	for (int i = A; i < Z; ++i) aref[nucl[i] - 'A'] = i;
	char*	ps = str;
	for (int i = 0; *ps; ++ps) {
	    if ('A' <= *ps && *ps <= 'Z') xref[i++] = aref[*ps - 'A'];
	    else if (*ps == '*' || (*ps == '-' && isspace(ps[1])))
		xref[i++] = gap_code;
	}
	for (int i = 0; (ps = fgets(str, MAXL, fd)) && i < NTS; ++i) {
	    int	c = tolower(*ps);
	    int	k;
	    if (c == '\n' || c == '#' || c == ';') {
		continue;		/* comment */
	    } else if (c == '*' || (c == '-' && isspace(ps[1]))) {
		k = gap_code;		/* unp */
		ps = cdr(str);
	    } else if (ps[1] == '=') {	/* gap penalty */
		if (c == 'u' || c == 'e') {
		    alprm.u = (VTYPE) atof(ps+2);
		    k = gap_code;
		    ps = cdr(str);
		} else if (c == 'v' || c == 'o') {
		    alprm.v = (VTYPE) atof(ps+2);
		    ps = cdr(str);
		    continue;
		} else {
		    ps = cdr(str);
		    continue;
		}
	    } else if ('a' <= c && c <= 'z') {
		k = aref[c - 'a'];
		if (k == -1) continue;
		ps = cdr(str);
	    } else {
		k = xref[i];
		ps = str;
	    }
	    int	j = 0;
	    for ( ; *ps && j < NTS; ps = cdr(ps), ++j) {
		mtx[k][xref[j]] = (VTYPE) 
		    ((atof(ps) / scale + param->b) * alprm.scale);
	    }
	    nrow[i] = j;
	}
	fclose(fd);
	minscr = INT_MAX;
	for (int p = 1; p < 16; p <<= 1) {
	    int	i = p + 1;
	    for (int j = A; j < Z; ++j) {
		int q = j - 1;
		if (nbit[q] == 1) {
		    if (maxs < mtx[i][j]) maxs = mtx[i][j];
		    if (minscr > mtx[i][j]) 
			minscr = mtx[i][j];
		    continue;
		}
		VTYPE	s = 0;
		for (int k = 1; k < 16; k <<= 1) {
		    if (q & k) s += mtx[i][k+1];
		}
		mtx[i][j] = mtx[j][i] = s / nbit[q];
	    }
	}
	for (int i = A; i < Z; ++i) {
	    mtx[i][_] = mtx[_][i] = ntsunp;
	    int	p = i - 1;
	    if (nbit[p] == 1) continue;
	    for (int j = A; j < Z; ++j) {
		int	q = j - 1;
		if (nbit[q] == 1) continue;
		for (int k = 1; k < p; k <<= 1) {
		    if (p & k) mtx[i][j] += mtx[k+1][j];
		}
		mtx[i][j] /= nbit[p];
	    }
	}
	avtrc = mtx[A][A] + mtx[C][C] + mtx[G][G] + mtx[T][T];
	nrmlf = avtrc /= 4;
	simunp = mtx[gap_code];
	drange = maxs - minscr;
}

void Simmtx::Pmtx()
{
	double	buf[AASCMB];
	double	buf2[PAMLEVELS];
	double	fscl = alprm.scale / 10.;
	double	fbias = 10. * param->b;
	VTYPE	unp_aas = (VTYPE) -(alprm.scale * param->u);

	dim = rows = ASIMD;
	mtx = SquareMtx();
	nrmlf = 1.;
	FILE*	fd = ftable.fopen(mdm_tab, "rb");
	if (!fd) fatal(not_found, mdm_tab);
	int	pam = param->p;
	INT	ii = (pam + PAMSTEP - 1) / PAMSTEP;
	pam = ii * PAMSTEP;	// discrete allowd pam values
	fseek(fd, ii * AASCMB * sizeof(double), SEEK_SET);
	if (fread(buf, sizeof(double), AASCMB, fd) != AASCMB)
	   fatal(Badpam, pam);
	fseek(fd, (PAMLEVELS + 1) * AASCMB * sizeof(double), SEEK_SET);
	if (fread(buf2, sizeof(double), PAMLEVELS, fd) == PAMLEVELS) {
	    nrmlf = (buf2[ii] + fbias) * fscl;
	    if (fread(buf2, sizeof(double), PAMLEVELS, fd) == PAMLEVELS)
		avtrc = (buf2[ii] + fbias) * fscl;
	    else	fatal(Badpam, pam);
	} else	fatal(Badpam, pam);
	fclose(fd);
	param->p = pam;
	INT	k = 0;
	for (INT i = 0; i < AAS; i++) {
	    ii = i + USE_EXG;
	    for (INT j = 0; j <= i; j++, k++) {
		INT	jj = j + USE_EXG;
		buf[k] += fbias;
		buf[k] *= fscl;
#if !FVAL
		buf[k] += 0.5;
#endif
		mtx[ii][jj] = mtx[jj][ii] = (VTYPE) buf[k];
	    }
	}
	for (int i = AMB; i < ASIMD; ++i) {
	    mtx[i][UNP] = mtx[UNP][i] = unp_aas;
// modified on 2020-03-17
	    mtx[i][SEC] = mtx[SEC][i] = mtx[i][CYS];
	}
	mtx[UNP][UNP] = 0;
	for (int i = 0; i < ASIMD; ++i)
	    mtx[i][NIL] = mtx[NIL][i] = 0;
	simunp = mtx[gap_code];
	minscr = mtx[TRP][CYS];
	drange = mtx[TRP][TRP] - minscr;
}

void Simmtx::Pmtx(const char* fname)
{
	FILE*	fd = ftable.fopen(fname, BLASTMAT, DEF_MAT_PATH);
	if (!fd) fatal("Matrix File %s was not found!\n", fname);
	int	aref[26];
	int	xref[26];
	int	nrow[26];
	int	amb = 0;
	int	asx = 0;
#if USE_GLX
	int	glx = 0;
#endif
	float	scale = 1.;
	VTYPE	unp_aas = (VTYPE) -(alprm.scale * param->u);

	dim = rows = ASIMD;
	mtx = SquareMtx();
	nrmlf = log(4.) * alprm.scale;	/* half bit */;
	char	str[MAXL];
	*str = _LCOMM;
	while (*str == _LCOMM) {
	    if (!fgets(str, MAXL, fd))
		fatal("Matrix File %s was incomplete!\n", fname);
	    char*	eql = strchr(str + 1, '=');
	    if (eql) {
		*eql = '\0';
		if (!strcmp(str + 1, strscale)) {
		    scale = atof(eql + 1);
		    if (scale <= 0.) scale = 1.;
		}
	    }
	}

	for (int i = 0; i < 26; ++i) {
	    aref[i] = nrow[i] = -1;
	    xref[i] = aacode[i];
	}
	for (int i = AMB; i <= GLX; ++i)
	    aref[amino[i] - 'A'] = i;
	char*	ps = str;
	for (int i = 0; *ps; ++ps) {
	    if ('A' <= *ps && *ps <= 'Z') {
		xref[i] = aref[*ps - 'A'];
		switch (xref[i++]) {
		    case AMB: ++amb; break;
		    case ASX: ++asx; break;
#if USE_GLX
		    case GLX: ++glx; break;
#endif
		}
	    } else if (*ps == '*' || (*ps == '-' && isspace(ps[1])))
		xref[i++] = gap_code;
	}
	for (int i = 0; (ps = fgets(str, MAXL, fd)) && i < AAS; ++i) {
	    int	k;
	    if (*ps == '*' || (*ps == '-' && isspace(ps[1]))) {
		k = gap_code;
		ps = cdr(str);
	    } else if ('A' <= *ps && *ps <= 'Z') {
		k = aref[*ps - 'A'];
		if (k == -1) continue;
		ps = cdr(str);
	    } else {
		k = xref[i];
		ps = str;
	    }
	    int	j = 0;
	    for ( ; *ps && j < AAS; ps = cdr(ps), ++j) {
		mtx[k][xref[j]] = (VTYPE) 
		    ((atof(ps) / scale + param->b) * alprm.scale);
	    }
	    nrow[i] = j;
	}
//	trianglurer matrix
	avtrc = 0.;
	for (int k = AMB; k < ASX; ++k) {
	    for (int j = nrow[k]; j < ASX; ++j)
		if (j >= 0) mtx[k][j] = mtx[j][k];
	    mtx[k][UNP] = mtx[UNP][k] = unp_aas;
	    avtrc += mtx[k][k];
	}
	avtrc /= 20;

	mtx[UNP][UNP] = 0;
	for (int i = 0; i < ASIMD; ++i)
	    mtx[i][NIL] = mtx[NIL][i] = 0;

	if (!amb) {
	    for (int j = AMB; j < AAS; ++j) 
		mtx[AMB][j] = mtx[j][AMB] = 0;
	}
	if (!asx)
	    for (int j = AMB; j < AAS; ++j)
		mtx[ASX][j] = mtx[j][ASX] = 
		    (mtx[ASN][j] + mtx[ASP][j]) / 2;
//	 modified on 2020-03-17
#if USE_GLX
	if (!glx)
	    for (int j = AMB; j < AAS; ++j)
		mtx[GLX][j] = mtx[j][GLX] = 
		    (mtx[GLN][j] + mtx[GLU][j]) / 2;
#else
	for (int j = AMB; j < AAS; ++j)
	    mtx[SEC][j] = mtx[j][SEC] = mtx[j][CYS];
#endif	// USE_GLX
	fclose(fd);
	simunp = mtx[gap_code];
	minscr = mtx[TRP][CYS];
	drange = mtx[TRP][TRP] - minscr;
}

void Simmtx::Hmtx(ComPmt compmt, Simmtx* pm)
{
	dim = TSIMD;
	rows = compmt == TxT? dim: SER2;
	mtx = SquareMtx();
	VTYPE	unp_aas = (VTYPE) -(alprm.scale * param->u);
	VTYPE	trm_aas = (VTYPE) -(alprm.scale * alprm2.o);
	VTYPE	trm_trm = (VTYPE) (alprm.scale * pm->mtx[ALA][ALA]);
	nrmlf = pm->nrmlf;
	avtrc = pm->avtrc;
	for (int i = 0; i < SER2; ++i) {
	    for (int j = 0; j < i; ++j)
		mtx[i][j] = mtx[j][i] = pm->mtx[i][j];
	    mtx[i][i] = pm->mtx[i][i];
	}
	for (int i = 0; i < TSIMD; ++i)
	    mtx[i][SER2] = mtx[SER2][i] = mtx[SER][i];
	for (int i = AMB; i < TSIMD; ++i) {
	    mtx[UNP][i] = mtx[i][UNP] = unp_aas;
	    mtx[SEC][i] = mtx[i][SEC] = 
	    mtx[TRM][i] = mtx[i][TRM] = trm_aas;
	}
	mtx[UNP][UNP] = 0;
	mtx[SEC][SEC] = mtx[CYS][CYS];
	if (compmt == TxT)
	    mtx[TRM][TRM] = mtx[TRM][TRM2] = 
	    mtx[TRM2][TRM] = mtx[TRM2][TRM2] = trm_trm;
	for (int i = 0; i < TSIMD; ++i)
	    mtx[NIL][i] = mtx[i][NIL] = 0;
	simunp = mtx[gap_code];
	minscr = mtx[TRP][CYS];
	drange = mtx[TRP][TRP] - minscr;
}

double Simmtx::avrmatch(INT* ConvTab) const
{
	if (AxB == DxD) return (double) nrmlf;
	int	count[21];	// Radix sort
	vclear(count, 21);
	for (int i = ALA; i <= VAL; ++i) count[ConvTab[i]]++;
	int	m = 0;
	int	n = 0;
	for (int* c = count; *c; ++c) {
	    int	t = *c;
	    *c = m;
	    m += t;
	    ++n;
	}
	int	origin[20];
	for (int i = ALA; i <= VAL; ++i) {
	    m = ConvTab[i];
	    origin[count[m]++] = i;
	}
	m = 0;
	double	ev = 0.;
	double*	comp = get_mdmcmp() - ALA;
	for (int i = 0; i < n; ++i) {
	    double	f = 0;
	    double	v = 0;
	    for (int j = m; j < count[i]; ++j) {
		int	a = origin[j];
		f += comp[a];
		for (int k = m; k < count[i]; ++k) {
		    int	b = origin[k];
		    v += comp[a] * comp[b] * mtx[a][b];
		}
	    }
	    ev += v / f;
	    m = count[i];
	}
	return (ev);
}

void Simmtx::fparam(FILE* fd) const
{
	if (AxB == DxD)
	    fprintf(fd, "s[=] (%.1lf), s[#] (%.1lf), ", smn[0], param->n);
	else
	    fprintf(fd, "PAM = %d, BIAS = %.1lf, ", param->p, param->b);
	fprintf(fd, "u = %.1lf, v = %.1lf\n", param->u, alprm.v);
}

void Simmtx::printmtx(FILE* fd) const
{
	for (int i = 0; i < dim; ++i) {
	    for (int j = 0; j < dim; ++j) {
		fprintf(fd, " %2d", int(mtx[i][j]));
	    }
	    fputc('\n', fd);
	}
}

double* Simmtxes::get_mdmcmp()
{
	if (mdmcomp) return (mdmcomp);
	mdmcomp = new double[20];
	FILE*	fd = ftable.fopen(mdm_cmp, "rb");
	if (!fd) fatal(not_found, mdm_cmp);
	if (fread(mdmcomp, sizeof(double), 20, fd) != 20)
	    fatal("Fail to read mdm_cmp file !\n");
	fclose(fd);
	return (mdmcomp);
}

void setpam(int pam, int mtxno)
{
	if (mtxno < max_simmtxes) defPprm[mtxno].p = pam;
}

int getpam(int mtxno)
{
	return (defPprm[mtxno < max_simmtxes? mtxno: 0].p);
}

float getsmn(int q)
{
	return (smn[q]);
}

void setNpam(int q, float v)
{
	if (q == 0 || q == 2 || q == 4)	smn[q] = v;
	smn[1] = (smn[0] + smn[2]) / 2;
	smn[3] = (smn[2] + smn[4]) / 2;
}

void setdefNprm(float n, float u, float v, float b, int c)
{
	if (c < 0 || c > 1) return;
	if (n > FPOPUP) defNprm[c].n = n;
	if (u > FPOPUP) defNprm[c].u = u;
	if (v > FPOPUP) defNprm[c].v = v;
	if (b > FPOPUP) defNprm[c].b = b;
}

void setdefPprm(int p, float u, float v, float b, int c)
{
	if (c < 0 || c > 1) return;
	if (p > POPUP)  defPprm[c].p = p;
	if (u > FPOPUP) defPprm[c].u = u;
	if (v > FPOPUP) defPprm[c].v = v;
	if (b > FPOPUP) defNprm[c].b = b;
}

void setalprm()
{
	int	lcl = algmode.lcl;

	promptin("u  (%.1f), v  (%.1f) : ", &alprm.u, &alprm.v);
	promptin("s[=] (%.1f), s[:] (%.1f), s[#] (%.1f) : ", 
		smn, smn + 2, smn + 4);
	promptin("pam (%d), bias (%.1f) : ", &defPprm->p, &defPprm->b);
	if (glocal == GLOBAL) {
	    if (alprm.ls > 1) {
		promptin("u1 (%.1f), k1 (%d) : ", &alprm.u1, &alprm.k1);
		if (alprm.k1 < 1) alprm.k1 = 1;
	    }
	    promptin("Local? (0-15) %d : ", &lcl);
	    algmode.lcl = lcl;
	}
	promptin("x (%.3f) y (%.3f) z (%.3f) i (%.3f) o (%.3f) : ",
	    &alprm2.x, &alprm2.y, &alprm2.z, &alprm2.y, &alprm2.o);
}

void setlsegs(int ls)
{
	if (ls > 0) alprm.ls = ls;
	else if (ls == QUERY) promptin("lsegs  (%d) : ", &alprm.ls);
	alprm.ls = min(alprm.ls, NOL);
	alprm.ls = max(alprm.ls, 1);
}

void readalprm(int& argc, const char**& argv, int oc)
{
const	char&	opt = argv[0][oc];
	if (!opt) return;
	bool	num = !(opt == 'h' || opt == 'r' || opt == 'I' || opt == 'r');
const	char*	vl = getarg(argc, argv, num, ++oc);
	if (!vl && !(opt == 'H' || opt == 'S')) return;

	const	char*	ps;
	const	char*	cs = vl? strchr(vl, ':'): 0;
	int	k = cs? atoi(cs + 1): 0;
	if (k > max_simmtxes - 1) k = 0;
	if (vl) {
	  switch (opt) {
	    case 'a': algmode.any = atoi(vl); break;	// boundary stringency
	    case 'b': defPprm[k].b = atof(vl); break;	// bias added to mat elem
	    case 'c': alprm2.jneibr = atoi(vl); break;	// 
	    case 'd': alprm2.desert = atoi(vl); break;	// giveup alignment if nohit
	    case 'e': alprm.u0 = atof(vl); break;	// background gep
	    case 'f': alprm.v0 = atof(vl); break;	// background gop
	    case 'g': alprm.gamma = atof(vl); break;	// gap parameter gamma
	    case 'h': 
		if ((ps = strchr(vl, ','))) alprm3.hpwing = atoi(ps + 1);
		if (!ps || ps > vl) alprm3.hydr = atof(vl);
		break;					// hydrophobicity
	    case 'i': IntronPrm.ip = atof(vl); break;	// intron penalty
	    case 'j': alprm.u1 = atof(vl); break;	// 2nd gap extension
	    case 'k': alprm.k1 = atoi(vl); break;	// flex point
	    case 'l': setlsegs(atoi(vl)); break;	// # linear pieces
	    case 'm': smn[0] = atof(vl); break;		// match
	    case 'n':
		smn[4] = atof(vl);			// mismatch
		if (smn[4] > 0) smn[4] = -smn[4];
		break;
	    case 'o': alprm2.o = atof(vl); break;	// premature termination
	    case 'p': defPprm[algmode.crs].p = atoi(vl); break;	// 1st pam
	    case 'q': defPprm[WlnPamNo].p = atoi(vl); break;	// 2nd pam
	    case 'r': 
		if ((ps = strchr(vl, ','))) alprm3.no_angle = atoi(ps + 1);
		if (!ps || ps > vl) alprm3.hpmt = atof(vl);
		break;					// hydrophobicity moment
	    case 's': alprm3.scnd = atof(vl); break;	// secondary structure prop.
	    case 't': alprm.tgapf = atof(vl); break;	// reduced terminal gap factor
	    case 'u': alprm.u = atof(vl); break;	// gap extension
	    case 'v': alprm.v = atof(vl); break;	// gap open
	    case 'w': alprm.sh = atoi(vl); break;	// band width
	    case 'x': alprm2.x = atof(vl); break;	// frame shift
	    case 'y': alprm2.y = atof(vl); break;	// boundary signal
	    case 'z': alprm2.z = atof(vl); break;	// coding potential
	    case 'A': alprm2.bti = atof(vl); break;	// init/term codon
	    case 'B': bpprm.factor = atof(vl); break;	// branch point
	    case 'D': bpprm.maxb3d = atoi(vl); break;	// max distance bp to 3'ss
	    case 'E': IntronPrm.elmt = atoi(vl); break;	// min exon len
	    case 'G': bpprm.g_alpha = atof(vl); 
		vl = strchr(vl, ',');
		if (vl) bpprm.g_beta = atof(vl + 1); 
		break;	// Gamma distribution parameters of pranch position
	    case 'H':
		setthr(atof(vl)); break;
	    case 'I': 
		IntronPrm.a1 = 1.; IntronPrm.a2 = 0.;
		sscanf(vl + 1, "%d %d %f %f %f %f %f %f %f %f %f %f %f %f", 
		&IntronPrm.llmt, &IntronPrm.rlmt, 
		&IntronPrm.mean, &IntronPrm.a1,
		&IntronPrm.m1, &IntronPrm.t1, &IntronPrm.k1,
		&IntronPrm.m2, &IntronPrm.t2, &IntronPrm.k2,
		&IntronPrm.a2, &IntronPrm.m3, &IntronPrm.t3, &IntronPrm.k3);
		break;
	    case 'J': alprm2.spb = atof(vl); break;	// matching intron position
	    case 'K': alprm2.termk1 = atoi(vl); break;	// max terminal gap length without penalty
	    case 'L': IntronPrm.llmt = atoi(vl); break;	// lower limit of intron
	    case 'M': IntronPrm.maxl = int(ktof(vl)); break;	// maximum expected length of intron
//	    case 'N': alprm2.nrmlipot = 1; break;	// normalize intron potential
	    case 'Q': IntronPrm.nquant = atoi(vl); break;	// number of steps of rough ILD
	    case 'S': alprm2.sss = atof(vl); 
		if (alprm2.sss > 1.) alprm2.sss /= 100.;
		break;
	    case 'T': IntronPrm.tlmt = atoi(vl); break;	// 
	    case 'U': alprm.ubh = atoi(vl); break;	// min vects for uni-dir Hirschberg
	    case 'V': alprm.maxsp = atof(vl); break;	// max traceback volume
	    case 'W': alprm2.w = atof(vl); break;	// match factor in very short alignment
	    case 'X': algmode.crs = atoi(vl); break;	// cross-species				// set cross-species switch
	    case 'Y': IntronPrm.fact = atof(vl); break;	// amplitude of intron pen.
	    case 'Z': alprm2.Z = atof(vl); break;	// intron potential
	    default:  break;
	  }
	} else {
	  switch (opt) {
	    case 'a': algmode.any = 3; break;
	    case 'H': setthr((double) (INT_MIN + 3)); break;
	    case 'S': alprm2.sss = 1.; break;
	    case 'X': algmode.crs = 2; break;
	  }
	}
}

void setSimmtxes(ComPmt ab, bool mdm, DefPrm* dp)
{
	int	upto = dp? 1: max_simmtxes;
	if (!dp) dp = (ab == DxD? defNprm: defPprm);
	int	odr[3] = {0, 1, 2};
	if (algmode.crs) swap(odr[0], odr[1]);
	DefPrm*	ddp = upto == 1? dp: dp + odr[0];
	if (smn[4] == FQUERY) smn[4] = ddp->n;
	else	ddp->n = smn[4];
	if (alprm.v == FQUERY) alprm.v = ddp->v;
	else	ddp->v = alprm.v;
	if (alprm.u == FQUERY) alprm.u = ddp->u;
	else	ddp->u = alprm.u;
	Simmtx**	psm = simmtxes.storedSimmtx;
	for (int mid = 0; mid < upto; ++mid, ++psm) {
	    ddp = dp + odr[mid];
	    if (!*psm) {
		const char* fname = mdm_file[mid];
		if (fname && !mdm)	*psm = new Simmtx(ab, fname, ddp);
		else if (ddp->p) *psm = new Simmtx(ab, ddp);
	    }
	}
	simmtxes.current = 0;
}

void resetSimmtxes(bool rmcomps)
{
	for (int i = 0; i < max_simmtxes; ++i) {
	    delete simmtxes.storedSimmtx[i];
	    simmtxes.storedSimmtx[i] = 0;
	}
	if (rmcomps) {
	    delete[] simmtxes.mdmcomp;
	    simmtxes.mdmcomp = 0;
	}
}
