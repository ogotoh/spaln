/*****************************************************************************
*
*	Calculate position specific score matrix based on m-th Markov model
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
static	const	int	nelmx2 = nelm * nelm;
static	const	int	nelmx3 = nelm * nelmx2;
static	const	int	bias[4] = {0, nelm, nelm + nelmx2, nelm + nelmx2 + nelmx3};

static	int	left = -1;
static	int	right = -1;
static	int	jnk = -1;
static	float	depsilon = -1;
static	const	char*	wdfq = 0;
static	const	char*	trim_out = 0;
static	const	char*	pssm_in = 0;
static	const	char*	trim_in = 0;
static	bool	entropy = false;
static	bool	binary = false;
static	float	sig_thr = 0.025;	// 
static	FILE*	ofd = stdout;

static void usage()
{
	fputs("Usage:\n", stdout);
	fputs("\tnpssm [-m[0-2]] [-r X.wdfq] [-f X.triN] [-b] [other options] X.SPNx[.gz]\n", stdout);
	fputs("Options:\n", stdout);
	fputs("\t-b [S]\t: binary output\n", stdout);
	fputs("\t-e [S]\t: output (relative) entropy\n", stdout);
	fputs("\t-f S\t: existing X.tri[3|5]\n", stdout);
	fputs("\t-h\t: display this\n", stdout);
	fputs("\t-H F\t: noise level (0.025)\n", stdout);
	fputs("\t-i S\t: existing pat_mat[.psm]\n", stdout);
	fputs("\t-j N\t: signal position\n", stdout);
	fputs("\t-l N\t: lower boundary\n", stdout);
	fputs("\t-m N\t: m(0|1|2)-th order MM\n", stdout);
	fputs("\t-o S\t: output file name\n", stdout);
	fputs("\t-q F\t: pseudocount\n", stdout);
	fputs("\t-r S\t: reference.wdfq\n", stdout);
	fputs("\t-t S\t: trimer output\n", stdout);
	fputs("\t-u N\t: upper boundary\n", stdout);
	exit(1);
}

class TriFreq {
	int	morder = 0;
	int	sites = 0;
	int	lod_size = 0;
	int	nsupport = 0;
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
	float	maxkld =0;		// max KLD value
	int	nmaxkld = 0;		// max KLD site
	int	nmaxg = 0;		// G-richest site
	int	nmaxt = 0;		// T-richest site
	char*	sname = 0;		// seq identifier
public:
	TriFreq(const int mo, const Seq* sd);
	TriFreq(const char* pssm_in);
	~TriFreq();
	void	read_trimer();
	void	count_trimer(Seq& sd);
	void	background3();
	void	ref_normalize();
	void	rel_frequency(const int n);
	void	markovmodel(const int phase);
	void	setrange();
	void	calculate_mmm(Seq& sd);
	void	output();
	void	ascii_output();
	void	binary_output();
};

TriFreq::TriFreq(const int mo, const Seq* sd)
	: morder(std::max(0, mo)), sites(sd? sd->len: 0)
{
	if (sd) {
const	    char*	sn = sd->sqname();
const	    char*	sl = strrchr(sn, '/');
	    if (sl) sn = sl + 1;
	    sname = strrealloc(sname, sn);
	}
	if (trim_in) read_trimer();
	else	trimbuf = new INT[sites * nelmx3];
	scomp = new float[2 * nelm];
	rcomp = scomp + nelm;
	vclear(scomp, 2 * nelm);
	if (jnk < 0 || left < 0 || right < 0)
	    klds = new float[sites];
	lod_size = nelm;
	if (morder > 0) lod_size += nelmx2;
	if (morder > 1) lod_size += nelmx3;
	lods = new float[lod_size * sites];
	if (morder < 1) return;

	sdifq = new float*[2 * nelm];
	rdifq = sdifq + nelm;
	*sdifq = new float[2 * nelmx2];
	*rdifq = *sdifq + nelmx2;
	for (int k = 1; k < nelm; ++k) {
	    sdifq[k] = sdifq[k - 1] + nelm;
	    rdifq[k] = rdifq[k - 1] + nelm;
	}
	vclear(*sdifq, 2 * nelmx2);
	if (morder < 2) return;

	strif = new float**[2 * nelm];
	rtrif = strif + nelm;
	*strif = new float*[2 * nelmx2];
	*rtrif = *strif + nelmx2;
	**strif = new float[2 * nelmx3];
	**rtrif = **strif + nelmx3;
	for (int k = 0; k < nelm; ++k) {
	    if (k) {
		strif[k] = strif[k - 1] + nelm;
		rtrif[k] = rtrif[k - 1] + nelm;
		strif[k][0] = strif[k - 1][0] + nelmx2;
		rtrif[k][0] = rtrif[k - 1][0] + nelmx2;
	    }
	    for (int i = 1; i < nelm; ++i) {
		strif[k][i] = strif[k][i - 1] + nelm;
		rtrif[k][i] = rtrif[k][i - 1] + nelm;
	    }
	}
	vclear(**strif, 2 * nelmx3);
}

TriFreq::TriFreq(const char* pssm_in)
{
	Strlist	sl(pssm_in, " ,");
	for (INT i = 0; i < sl.size(); ++i) {
	    PatMat	pm(sl[i]);
	    if (i) {
		mmm.mean = (nsupport * mmm.mean + pm.nsupport * pm.mmm.mean);
		mmm.mean /= (nsupport += pm.nsupport);
		if (pm.mmm.min < mmm.min) mmm.min = pm.mmm.min;
		if (pm.mmm.max > mmm.max) mmm.max = pm.mmm.max;
		continue;
	    }
	    left = 0;
	    jnk = pm.offset - 1;
	    sites = right = pm.cols;
	    lod_size = pm.rows;
	    morder = pm.order();
	    mmm = pm.mmm;
	    nsupport = pm.nsupport;
	    if (!wdfq) std::swap(lods, pm.mtx);
	}
}

TriFreq::~TriFreq() {
	delete[] trimbuf;
	delete[] scomp;
	delete[] klds;
	delete[] lods;
	delete[] sname;
	if (morder < 1 || !sdifq) return;
	delete[] *sdifq;
	delete[] sdifq;
	if (morder < 2 || !strif) return;
	delete[] **strif;
	delete[] *strif;
	delete[] strif;
}

void TriFreq::read_trimer()
{
static	const	char imcompatible[] = "%s is incompatible !\n";

	INT*	tmpbuf = 0;
	Strlist	sl(trim_in, " ,");
	for (INT i = 0; i < sl.size(); ++i) {
	    PatMatHead	header;
	    FILE* 	fd = fopen(sl[i], "rb");
	    if (!fd || fread(&header, sizeof(PatMatHead), 1, fd) != 1 ||
		header.vsize != int_size || header.nelm != nelm ||
		header.cols != nelmx3) fatal(imcompatible, sl[i]);
	    if (!sites) sites = header.rows;
	    else if (sites != header.rows) fatal(imcompatible, sl[i]);
	    INT	dsize = sites * nelmx3;
	    union {
		char	str[MAXL];
		int	many;
	    };
	    if (fread(str, sizeof(char), header.add, fd) != (INT) header.add)
		fatal(imcompatible, sl[i]);
	    str[header.add] = '\0';
	    if (!nsupport) nsupport = many;
	    if (!sname) sname = strrealloc(sname, str + int_size);
	    if (i == 0) {
		trimbuf = new INT[dsize];
		if (fread(trimbuf, int_size, dsize, fd) != dsize)
		    fatal(imcompatible, sl[i]);
	    } else {
	        if (!tmpbuf) tmpbuf = new INT[dsize];
		if (fread(tmpbuf, int_size, dsize, fd) != dsize)
		    fatal(imcompatible, sl[i]);
		for (INT n = 0; n < dsize; ++n)
		    trimbuf[n] += tmpbuf[n];
	    }
	}
	delete[] tmpbuf;
}

void TriFreq::count_trimer(Seq& sd)
{
	INT*	trimer = trimbuf;
	srand((INT) time(0));
const	CHAR*	ss = sd.at(0);
const	CHAR*	tt = sd.at(sd.len - 2);
	for ( ; ss < tt; ss += sd.many, trimer += nelmx3) {
	    vclear(trimer, nelmx3);
const	    CHAR*	t = ss + sd.many;
	    for (const CHAR* s = ss; s < t; ++s) {
		int	a = ncredctab[*s];
const		CHAR*	r = s + sd.many;
		int	b = ncredctab[*r];
		r += sd.many;
		int	c = ncredctab[*r];
		if (a < nelm && b < nelm && c < nelm)
		    ++trimer[a * nelmx2 + b * nelm + c];
	    }
	}
	if (!trim_out) return;

// write trimer frequecy table

	ofd = fopen(trim_out, "wb");
	if (!ofd) fatal(no_file, trim_out);
	PatMatHead	header = {0, int_size, nelm, 0, sd.len, nelmx3};
	struct	Buf	{int many; char sname[MAXL];} buf;
	buf.many = sd.many;
const	char*	sn = sd.sqname();
const	char*	sl = strrchr(sn, '/');
	if (sl) sn = sl + 1;
	strcpy(buf.sname, sn);
	buf.sname[max_add_size - int_size] = 0;	//	 truncate
	header.add = CHAR(int_size + (strlen(buf.sname) + 3) / 4 * 4);
	fwrite(&header, sizeof(PatMatHead), 1, ofd);
	fwrite(&buf, sizeof(char), header.add, ofd);
	fwrite(trimbuf, sizeof(int), sites * nelmx3, ofd);
	fclose(ofd);
}

void TriFreq::ref_normalize()
{
	float	rff = 0;
	if (morder == 0) {
	    for (int i = 0; i < nelm; ++i) rff += rcomp[i];
	    for (int i = 0; i < nelm; ++i) rcomp[i] /= rff;
	    return;
	} else if (morder == 1) {
	    float*	w = *rdifq;
	    for (int i = 0; i < nelmx2; ++i) rff += *w++;
	    for (int i = 0; i < nelm; ++i) {
		rcomp[i] = 0;
		for (int j = 0; j < nelm; ++j) {
		    rdifq[i][j] /= rff;
		    rcomp[i] += rdifq[i][j];
		}
	    }
	    return;
	}
	float*	w = **rtrif;
	for (int i = 0; i < nelmx3; ++i) rff += *w++;
	for (int i = 0; i < nelm; ++i) {
	    rcomp[i] = 0;
	    for (int j = 0; j < nelm; ++j) {
		rdifq[i][j] = 0;
		for (int k = 0; k < nelm; ++k) {
		    rtrif[i][j][k] /= rff;
		    rdifq[i][j] += rtrif[i][j][k];
		}
		rcomp[i] += rdifq[i][j];
	    }
	}
}

void TriFreq::background3()
{
	Strlist	sl(wdfq, " ,");
	for (INT i = 0; i < sl.size(); ++i) {
	    FILE*	fd = fopen(sl[i], "r");
	    if (!fd) fatal("%s can't open !\n", sl[i]);
	    char	kmer[16];
	    INT 	freq;
	    char	str[MAXL];
	    int 	code[3] = {0, 0, 0};
	    while (fgets(str, MAXL, fd)) {
		sscanf(str, "%s %u", kmer, &freq);
		int	k = strlen(kmer);
		if (--k > morder) break;
		switch (k) {
		    case 0: rcomp[code[0]++] += freq; break;
		    case 1: (*rdifq)[code[1]++] += freq; break;
		    case 2: (**rtrif)[code[2]++] += freq; break;
		}
	    }
	    fclose(fd);
	}
	ref_normalize();
}

void TriFreq::rel_frequency(const int n)
{
static	const	float	maxh = log(4.);
const	int	np1 = n + 1;
	float	h = 0;
	float*	lod = lods + n * lod_size;
	for (int i = 0; i < nelm; ++i) {
	    float	r = scomp[i] / sumi;
	    if (entropy) h += (r > 0? r * log(r): 0);
	    else	*lod++ = r;
	}
	if (entropy && morder == 0)
	    fprintf(ofd, "%d\t%9.5lf %9.5lf ", np1, -h, maxh + h);
	if (morder > 0) {
	    h = 0;
	    for (int i = 0; i < nelm; ++i)
		for (int j = 0; j < nelm; ++j) {
		    float	r = sdifq[i][j] / sumi;
		    if (entropy) h += (r > 0? r * log(r): 0);
		    else	*lod++ = r;
		}
	    if (entropy &&  morder == 1)
		fprintf(ofd, "%d\t%9.5lf %9.5lf ", np1, -h, 2 * maxh + h);
	}
	if (morder > 1) {
	  h = 0;
	  for (int i = 0; i < nelm; ++i)
	    for (int j = 0; j < nelm; ++j)
	      for (int k = 0; k < nelm; ++k) {
		float	r = strif[i][j][k] / sumi;
		if (entropy) h += (r > 0? r * log(r): 0);
		else	*lod++ = r;
	      }
	  if (entropy &&  morder == 2)
		fprintf(ofd, "%d\t%9.5lf %9.5lf ", np1, -h, 3 * maxh + h);
	}
	if (entropy) fputc('\n', ofd);
}

void TriFreq::markovmodel(const int phase)
{
	if (phase == 0 && !klds && !entropy) return;
	if (entropy)
	    fprintf(ofd, ">%s [%d:%d]\n", sname, nsupport, sites);
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
		rel_frequency(n);
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
		fprintf(ofd, "%d\t%15.7le", np1, h);
		for (int i = 0; i < nelm; ++i) {
		    fprintf(ofd, "\t%7d", int(scomp[i]));
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
	        if (entropy && morder == 1) fprintf(ofd, "%d\t%15.7le ", np1, h);
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
		if (entropy && morder == 2) fprintf(ofd, "%d\t%15.7le ", np1, h);
	    }
	    if (entropy) fputc('\n', ofd);
	}
}

void TriFreq::setrange()
{
	if (!klds) return;
	if (jnk < 0)
	    jnk = (nmaxg + 1 == nmaxt)? nmaxg - 1: nmaxg;	// ...|gt...
	sig_thr *= maxkld;
	if (left < 0) {
	    for (left = 0; left < nmaxkld; ++left)
		if (klds[left] >= sig_thr) break;
	}
	if (right < 0) {
	    for (right = sites; --right > nmaxkld; --right)
		if (klds[right] >= sig_thr) break;
	    ++right;
	}
}

// calculate signal score for each member in sd

void TriFreq::calculate_mmm(Seq& sd)
{
const	CHAR*	ss = sd.at(left);
const	CHAR*	tt = sd.at(right - 2);
	float	ng_res = INT_MIN;
	float*	lod = lods + left * lod_size;
	float* score = new float[sd.many];
	vclear(score, sd.many);
	for (int n = 0; ss < tt; ++n, ss += sd.many, lod += lod_size) {
	    const	CHAR*	t = ss + sd.many;
	    float*	scr = score;
	    for (const CHAR* s = ss; s < t; ++s, ++scr) {
		int	m = 0;
		int	c = 0;
		float	h = 0;
		for (const CHAR* r = s; m <= morder; ++m, r += sd.many) {
		    if (ncredctab[*r] >= nelm) break;
		    c = nelm * c + ncredctab[*r];
		    if ((n + m) < morder) h += lod[c + bias[m]];
		}
		if (--m == morder) *scr += (h + lod[c + bias[m]]);
		else	*scr = ng_res;
	    }
	}
const	float*	send = score + sd.many;
	ng_res += 100;
	nsupport = 0;
	for (const float* scr = score; scr < send; ++scr) {
	    if (*scr < ng_res) continue;
	    if (*scr < mmm.min) mmm.min = *scr;
	    if (*scr > mmm.max) mmm.max = *scr;
	    mmm.mean += *scr;
	    ++nsupport;
	}
	if (nsupport) mmm.mean /= nsupport;
	delete[] score;
}

void TriFreq::ascii_output()
{
	if (wdfq)
	    fprintf(ofd, "%d %d %d 1 0 %7.3f %7.3f %7.3f %d\n",
		right - left, lod_size, jnk - left + 1, 
		mmm.min, mmm.mean, mmm.max, nsupport);
	float*	lod = lods + left * lod_size;
	for (int n = left; n < right; ++n) {
	    for (int i = 0; i < lod_size; ++i)
		fprintf(ofd, "%9.5f ", *lod++);
	    putc('\n', ofd);
	}
}

void TriFreq::binary_output()
{
static	const	char	write_error[] = "fail to write binary data !\n";
const	float*	lod = lods + left * lod_size;
	if (wdfq) {	// pssm
	    PatMat	pm(lod_size, right - left, jnk - left + 1);
	    pm.min_elem = *vmin(lod, right - left);
	    pm.transvers = 0;
	    pm.nsupport = nsupport;
	    pm.mmm = mmm;
	    if (fwrite(&pm, sizeof(PatMat), 1, ofd) != 1)
		fatal(write_error);
	} else {	// relative frequency
	    PatMatHead	header = {1, sizeof(float), nelm, 0, right - left, lod_size};
	    if (fwrite(&header, sizeof(PatMatHead), 1, ofd) != 1)
		fatal(write_error);
	}
const	size_t	dsize = (right - left) * lod_size;
	if (fwrite(lod, sizeof(float), dsize, ofd) != dsize)
	    fatal(write_error);
}

void TriFreq::output()
{
	if (left < 0) left = 0;
	if (right < 0) right = sites;
	else if (right > sites) {
	    prompt("Warning! output range is shrinked %d <- %d !\n", 
		sites, right);
	    right = sites;
	}
	if (binary)	binary_output();
	else		ascii_output();
}

int main(int argc, const char** argv)
{
	int	mo = -1;
const	char*	outfn = 0;
	TriFreq*	tfq = 0;
	while (--argc > 0 && **++argv == OPTCHAR) {
	    const char* val = argv[0] + 2;
	    switch (argv[0][1]) {
		case 'b':
		    if ((val = getarg(argc, argv))) outfn = val;
		    binary = true;
		    break;
		case 'e':
		    if ((val = getarg(argc, argv))) outfn = val;
		    entropy = true;
		    break;
		case 'f':
		    if ((val = getarg(argc, argv))) trim_in = val;
		    break;
		case 'h': usage();
		case 'H':
		    if ((val = getarg(argc, argv, true))) sig_thr = atof(val);
		    if (sig_thr > 1.) sig_thr /= 100.;
		    break;
		case 'i':	// existing pssm
		    if ((val = getarg(argc, argv))) pssm_in = val;
		    break;
		case 'j':	// lower bound
		    if ((val = getarg(argc, argv, true))) jnk = atoi(val);
		    break;
		case 'l':	// lower bound
		    if ((val = getarg(argc, argv, true))) left = atoi(val);
		    break;
		case 'm':
		    if ((val = getarg(argc, argv, true)))
			mo = std::min(atoi(val), 2);
		    break;
		case 'o':
		    if ((val = getarg(argc, argv))) outfn = val;
		    break;
		case 'q':	// psuedocount
		    if ((val = getarg(argc, argv, true))) depsilon = atof(val);
		    break;
		case 'r':
		    if ((val = getarg(argc, argv))) wdfq = val;
		    break;
		case 't':
		    if ((val = getarg(argc, argv))) trim_out = val;
		    break;
		case 'u':	// upper bound
		    if ((val = getarg(argc, argv, true))) right = atoi(val);
		    break;
	    }
	}
	if (outfn) {
	    char	str[MAXL];
	    strcpy(str, outfn);
	    if (binary) {
const		char*	dot = strrchr(str, '.');
		if (!dot || strcmp(dot, patmat_ext))
		    strcat(str, patmat_ext);
	    }
	    if (pssm_in && !strcmp(pssm_in, str))
		fatal("Output %s is the same as input !\n", str);
	    ofd = fopen(str, "w");
	    if (!ofd) fatal(no_file, str);
	}
	if (pssm_in) {			// convert or synthesize existing pssm(s)
	    tfq = new TriFreq(pssm_in);
	    if (trim_in) tfq->read_trimer();
	    if (wdfq) {
		tfq->background3();	// Get backgroud data
		tfq->setrange();
		tfq->markovmodel(1);
	    } else
		wdfq = "";
	} else if (argc) {		// create pssm from MSA
	    Seq	sd(*argv);
	    tfq = new TriFreq(mo, &sd);
	    if (!trim_in) tfq->count_trimer(sd);
	    if (mo < 0 && trim_out) {
		delete tfq;
		return (0);
	    }
	    if (depsilon < 0) depsilon = 1./((float) sd.many);
	    if (wdfq) tfq->background3();	// Get backgroud data
	    tfq->markovmodel(0);
	    if (!entropy && wdfq) {
		tfq->setrange();
		tfq->markovmodel(1);
		tfq->calculate_mmm(sd);	// min, mean, max
	    }
	} else if (trim_in) {
	    tfq = new TriFreq(mo, 0);
	    if (wdfq) tfq->background3();	// Get backgroud data
	    tfq->markovmodel(0);
	} else	usage();

	if (!entropy) tfq->output();
	fclose(ofd);
	delete tfq;
	return (0);
}
