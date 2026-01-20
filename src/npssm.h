/*****************************************************************************
*
*	Headers to position specific score matrix based on m-th Markov model
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
*	Osamu Gotoh, Ph.D.	(2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<gotoh.osamu.67a@st.kyoto-u.ac.jp>>
*
*****************************************************************************/

#include "seq.h"
#include <math.h>

static	const	int	int_size = (int) sizeof(int);
static	const	int	nelm = 4;
static	const	int	max_cols = 32;
static	const	int	frq_size = nelm * max_cols;
static	const	int	nelmx2 = nelm * nelm;
static	const	int	nelmx3 = nelm * nelmx2;
static	const	int	bias[4] = {0, nelm, nelm + nelmx2, nelm + nelmx2 + nelmx3};

static	int	left = -1;
static	int	right = -1;
static	int	jnk = -1;
static	float	depsilon = -1;
static	bool	entropy = false;
static	const	char*	wdfq = 0;
static	const	char*	pssm_in = 0;
static	const	char*	trim_in = 0;

/*************************************************************************
*	TriFreq
*************************************************************************/

struct	PsFrqMat;

class TriFreq {
	struct	Buf	{int many; char sname[MAXL];};
	INT*	trimbuf = 0;		// trimer count
	float	sumi;
	float*	scomp = 0;
	float*	rcomp = 0;
	float**	sdifq = 0;
	float**	rdifq = 0;
	float***	strif = 0;
	float***	rtrif = 0;
	ExpectMmm	mmm = {FLT_MAX, 0, -FLT_MAX};
	float*	klds = 0;		// relative entropy (KLD)
	float*	lods = 0;		// log odds
	float	maxkld = 0;		// max KLD value
	int	nmaxkld = 0;		// max KLD site
	int	nmaxg = 0;		// G-richest site
	int	nmaxt = 0;		// T-richest site
public:
	char*	sname = 0;		// seq identifier
	int	morder = 0;
	int	sites = 0;
	int	lod_size = 0;
	int	nsupport = 0;
	TriFreq() {read_trimer();}
	TriFreq(const int mo, const Seq* sd = 0);
	TriFreq(const char* fname);
	~TriFreq();
template <typename file_t>
	void	read_trimer(file_t fd, const char* fname, int i);
	void	read_trimer();
	void	count_trimer(Seq& sd);
	void	background3(const char* wdfq);
	void	ref_normalize();
template <typename file_t>
	void	rel_frequency(file_t ofd, const int n);
	void	markovmodel(WriteFile& wfp, const int phase, const char* outf);
template <typename file_t>
	void	markovmodel(file_t fd, const int phase);
	void	setrange();
	void	calculate_mmm(Seq& sd);
	void	calculate_mean(TriFreq* ref_tfq = 0);
	void	output(WriteFile& wfp);
template <typename file_t>
	void	text_output(file_t fd);
template <typename file_t>
	bool	write_binary(file_t fd, Buf& buf, PatMatHead& header);
template <typename file_t>
	void	binary_output(file_t fd);
	friend	struct	PsFrqMat;
};

template <typename file_t>
void TriFreq::read_trimer(file_t fd, const char* fname, int i)
{
	PatMatHead	header;
	if (fread(&header, sizeof(PatMatHead), 1, fd) != 1 ||
	    header.vsize != int_size || header.nelm != nelm ||
	    header.cols != nelmx3) fatal(read_error, fname);
	if (!sites) sites = header.rows;
	else if (sites != header.rows) fatal(read_error, fname);
	INT	dsize = sites * nelmx3;
	union {
	    char	str[MAXL];
	    int	many;
	};
	if (fread(str, sizeof(char), header.add, fd) != (INT) header.add)
	    fatal(read_error, fname);
	str[header.add] = '\0';
	if (!nsupport) nsupport = many;
	if (!sname) sname = strrealloc(sname, str + int_size);
	if (i == 0) {
	    trimbuf = new INT[dsize];
	    if (fread(trimbuf, int_size, dsize, fd) != dsize)
		fatal(read_error, fname);
	} else {
	    INT* tmpbuf = new INT[dsize];
	    if (fread(tmpbuf, int_size, dsize, fd) != dsize)
		fatal(read_error, fname);
	    for (INT n = 0; n < dsize; ++n)
		trimbuf[n] += tmpbuf[n];
	    delete[] tmpbuf;
	}
}

template <typename file_t>
bool TriFreq::write_binary(file_t fd, Buf& buf, PatMatHead& header)
{
	INT	dsize = sites * nelmx3;
	if ((fwrite(&header, sizeof(PatMatHead), 1, fd) != 1) ||
	    (fwrite(&buf, sizeof(char), header.add, fd) != header.add) ||
	    (fwrite(trimbuf, sizeof(int), dsize, fd) != dsize))
		return (false);
	return (true);
}

template <typename file_t>
void TriFreq::text_output(file_t fd)
{
	char	str[MAXL];
	if (wdfq || pssm_in) {
	    sprintf(str, "%d %d %d 1 0 %7.3f %7.3f %7.3f %d\n",
		right - left, lod_size, jnk - left + 1, 
		mmm.min, mmm.mean, mmm.max, nsupport);
	    fputs(str, fd);
	} else if (!wdfq && trim_in) {
const	    char*	sl = strrchr(trim_in, '/');
	    sprintf(str, ">%s %d %d %d", sl? sl + 1: trim_in, 
		right - left, lod_size, jnk - left + 1);
	    fputs(str, fd);
	}
	float*	lod = lods + left * lod_size;
	for (int n = left; n < right; ++n) {
	    for (int i = 0; i < lod_size; ++i) {
		sprintf(str, " %9.5f", *lod++);
		fputs(str, fd);
	    }
	    fputc('\n', fd);
	}
}

template <typename file_t>
void TriFreq::binary_output(file_t fd)
{
const	float*	lod = lods + left * lod_size;
	if (wdfq) {	// pssm
	    PatMat	pm(lod_size, right - left, jnk - left + 1);
	    pm.min_elem = *vmin(lod, right - left);
	    pm.transvers = 0;
	    pm.nsupport = nsupport;
	    if (mmm.min == FLT_MAX) mmm.min = 0;
	    if (mmm.max == -FLT_MAX) mmm.max = 0;
	    pm.mmm = mmm;
	    if (fwrite(&pm, sizeof(PatMat), 1, fd) != 1)
		fatal(write_error);
	} else {	// relative frequency
	    PatMatHead	header = {1, sizeof(float), nelm, 0, right - left, lod_size};
	    if (fwrite(&header, sizeof(PatMatHead), 1, fd) != 1)
		fatal(write_error);
	}
const	size_t	dsize = (right - left) * lod_size;
	if (fwrite(lod, sizeof(float), dsize, fd) != dsize)
	    fatal(write_error);
}

template <typename file_t>
void TriFreq::markovmodel(file_t ofd, const int phase)
{
	char	str[MAXL];
	if (entropy) {
	    sprintf(str, ">%s [%d:%d]\n", sname, nsupport, sites);
	    fputs(str, ofd);
	}
	float	maxg = 0;
	float	maxt = 0;
	INT*	trimer = trimbuf;
	int	mo = 0;
const	int	mm = phase? morder: 0;
	float*	lod = lods;
	for (int n = 0; n < sites - 2 || (mm + mo++) < 2; ++n) {
const	    int	np1 = n + 1;
	    vclear(scomp, nelm);
	    if (morder > 0) vclear(*sdifq, nelmx2);
	    if (morder > 1) vclear(**strif, nelmx3);
	    for (int i = 0, m = 0; i < nelm; ++i) {
		for (int j = 0; j < nelm; ++j) {
		    for (int k = 0; k < nelm; ++k, ++m) {
			if (mo == 0) {
			    scomp[i] += trimer[m];
			    if (morder > 0) sdifq[i][j] += trimer[m];
			    if (morder > 1) strif[i][j][k] += trimer[m];
			} else if (mo == 1) {
			    scomp[j] += trimer[m];
			    if (morder > 0) sdifq[j][k] += trimer[m];
			} else if (mo == 2) {
			    scomp[k] += trimer[m];
			}
		    }
		}
	    }
	    sumi = 0;
	    for (int i = 0; i < nelm; ++i) sumi += scomp[i];
	    if (n < sites - 3) trimer += nelmx3;
	    if (!wdfq) {
		if (ofd)	rel_frequency(ofd, n);
		continue;
	    }
	    float	h = 0;
	    for (int i = 0; i < nelm; ++i) {
		float	p = scomp[i] / sumi;
		float  r = p / rcomp[i];
		h += (p > 0? p * log(r): 0);
		if (phase) {
const		    float	sig = log10((r + depsilon) / (1 + depsilon));
		    *lod++ = sig;
		}
	    }
	    if (entropy) {
		sprintf(str, "%d\t%15.7le", np1, h);
		fputs(str, ofd);
		for (int i = 0; i < nelm; ++i) {
		    sprintf(str, "\t%7d", int(scomp[i]));
		    fputs(str, ofd);
		}
	    } else if (phase == 0) {
		klds[n] = h;
		if (h > maxkld) {maxkld = h; nmaxkld = n;}
		if (scomp[2] > maxg) {maxg = scomp[2]; nmaxg = n;}
		if (scomp[3] > maxt) {maxt = scomp[3]; nmaxt = n;}
		continue;
	    }
	    if (morder > 0) {
		h = 0;
		for (int i = 0; i < nelm; ++i) {
		    float  ww = (scomp[i] > 0)? rcomp[i] / scomp[i]: 0.;
		    for (int j = 0; j < nelm; ++j) {
			float p = sdifq[i][j] / sumi;
			float r = ww? ww * sdifq[i][j] / rdifq[i][j]: 0;
			if (entropy) h += (p > 0? p * log(p / rdifq[i][j]): 0);
	    		else {
const			    float	sig = log10((r + depsilon) / (1 + depsilon));
			    *lod++ = sig;
			}
		    }
		}
	        if (entropy && morder == 1) {
		    sprintf(str, "%d\t%15.7le ", np1, h);
		    fputs(str, ofd);
		}
	    }
	    if (morder > 1) {
		h = 0;
		for (int i = 0; i < nelm; ++i) {
		  for (int j = 0; j < nelm; ++j) {
		    float	ww = (sdifq[i][j] > 0)? rdifq[i][j] / sdifq[i][j]: 0;
		    for (int k = 0; k < nelm; ++k) {
			float p =  strif[i][j][k] / sumi;
			float r = ww? ww * strif[i][j][k] / rtrif[i][j][k]: 0;
			if (entropy) h += (p > 0? p * log(p / rtrif[i][j][k]): 0);
			else {
const			    float	sig = log10((r + depsilon) / (1 + depsilon));
			    *lod++ = sig;
			}
		    }
		  }
		}
		if (entropy && morder == 2) {
		    sprintf(str, "%d\t%15.7le ", np1, h);
		    fputs(str, ofd);
		}
	    }
	    if (entropy) fputc('\n', ofd);
	}
}

/*************************************************************************
*	PsFrqMat: Position specific frequecy matrix
*************************************************************************/

struct PsFrqMat {
	char	id[ID_SIZE];
	SHORT	cols = 0;
	INT	many = 0;
	union {
	    float	pfm[max_cols][nelm];
	    float	frq[frq_size];
	};
template <typename file_t>
	void	read_binary(file_t fd, const char* fname);
template <typename file_t>
	void	read_text(file_t fd, const int n_sample = 0);
template <typename file_t>
	int	fget(file_t& dmyfd, const char* fn = 0);
	void	fuse(const char* fname);
	void	catenate(const char* fname);
	void	name2id(const char* name);
template <typename file_t>
	void	write_binary(file_t fd) {
	    if (fwrite(this, sizeof(PsFrqMat), 1, fd) != 1)
		fatal(write_error, id);
	}
template <typename file_t>
	void	write_text(file_t fd);
	void	output(WriteFile& wfp);
	PsFrqMat() {}
	PsFrqMat(const char* fname, const int n_sample = 0) {
	    vclear(id, ID_SIZE);
	    vclear(frq, frq_size);
const	    char*	la = strchr(fname, '<');
	    Strlist	sl(fname, "< ,+");
	    ReadFile	fp(sl[0]);
	    if (!fp.dtype) fatal(not_found, sl[0]);
	    if (fp.gzfd) {
		if (fp.dtype == 1) read_text(fp.gzfd, n_sample); else
		if (fp.dtype == 2) read_binary(fp.gzfd, sl[0]);
	    } else if (fp.fd) {
		if (fp.dtype == 1) read_text(fp.fd, n_sample); else
		if (fp.dtype == 2) read_binary(fp.fd, sl[0]);
	    }
	    name2id(sl[0]);
	    for (INT i = 1; i < sl.size(); ++i) {
		if (la) catenate(sl[i]);
		else	fuse(sl[i]);
	    }
	}
	PsFrqMat(const TriFreq* tfq) {
	    vclear(id, ID_SIZE);
	    name2id(tfq->sname);
	    many = tfq->nsupport;
	    if (left < 0) left = 0;
	    if (right < 0) right = std::min(tfq->sites, max_cols);
	    cols = right - left;
	    vcopy(frq, tfq->lods + nelm * left, nelm * cols);
	}
};

template <typename file_t>
void PsFrqMat::read_binary(file_t fd, const char* fname)
{
	if (fread(this, sizeof(PsFrqMat), 1, fd) != 1)
	    fatal(read_error, fname);
}

template <typename file_t>
void PsFrqMat::read_text(file_t fd, const int n_sample)
{
	char	str[MAXL];
	float*	f = frq;
	many = n_sample;
	while (fgets(str, MAXL, fd)) {
	    if (*str == '>') {
		sscanf(str, "%*s [%d:%hd]", &many, &cols);
		continue;
	    }
	    for (char* ps = str; *ps; ps = cdr(ps))
		*f++ = (float) atof(ps);
	}
	cols = (f - frq) / nelm;
}

template <typename file_t>
void PsFrqMat::write_text(file_t fd)
{
	char	str[MAXL];
	sprintf(str, ">%s [%d:%d]\n", id, many, cols);
	fputs(str, fd);
	float*	f = frq;
	for (int c = 0; c < cols; ++c) {
	    for (int r = 0; r < nelm; ++r) {
		sprintf(str, " %9.5f", *f++);
		fputs(str, fd);
	    }
	    fputc('\n', fd);
	}
}

template <typename file_t>
int PsFrqMat::fget(file_t& dmyfd, const char* fn)
{
	if (dmyfd) {
	    fclose(dmyfd);
	    dmyfd = 0;
	}
	ReadFile	fp(fn);
	if (fp.fd) {
	    if (fp.dtype == 1) read_text(fp.fd); else
	    if (fp.dtype == 2) read_binary(fp.fd, fn); 
	} else if (fp.gzfd) {
	    if (fp.dtype == 1) read_text(fp.gzfd); else
	    if (fp.dtype == 2) read_binary(fp.gzfd, fn);
	} else
	    fatal(not_found, fn);
const	char*	sl = strrchr(fn, '/');
	sl = sl? sl + 1: fn;
	strncpy(id, sl, ID_SIZE - 1);
	char*	dot = strchr(id, '.');
	if (dot) *dot = '\0';
	return (fp.dtype);
}
