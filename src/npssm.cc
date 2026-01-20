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

#include "kmers.h"
#include "npssm.h"
#include <math.h>

static	bool	binary = false;
static	float	sig_thr = 0.025;

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
	lod_size = bias[morder + 1];
	scomp = new float[2 * lod_size];
	rcomp = scomp + lod_size;
	vclear(scomp, 2 * lod_size);
	if (jnk < 0 || left < 0 || right < 0)
	    klds = new float[sites];
	lods = new float[lod_size * sites];
	if (morder < 1) return;

	sdifq = new float*[2 * nelm];
	rdifq = sdifq + nelm;
	*sdifq = scomp + bias[1];
	*rdifq = rcomp + bias[1];
	for (int k = 1; k < nelm; ++k) {
	    sdifq[k] = sdifq[k - 1] + nelm;
	    rdifq[k] = rdifq[k - 1] + nelm;
	}
	if (morder < 2) return;

	strif = new float**[2 * nelm];
	rtrif = strif + nelm;
	*strif = new float*[2 * nelmx2];
	*rtrif = *strif + nelmx2;
	**strif = scomp + bias[2];
	**rtrif = rcomp + bias[2];
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
}

TriFreq::TriFreq(const char* fname)
{
	Strlist	sl(fname, " ,");
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
	if (depsilon < 0) depsilon = 1./((float) nsupport);
}

TriFreq::~TriFreq() {
	delete[] trimbuf;
	delete[] scomp;
	delete[] klds;
	delete[] lods;
	delete[] sname;
	if (morder < 1 || !sdifq) return;
	delete[] sdifq;
	if (morder < 2 || !strif) return;
	delete[] *strif;
	delete[] strif;
}

void TriFreq::read_trimer()
{
	Strlist	sl(trim_in, " ,");
	for (INT i = 0; i < sl.size(); ++i) {
	    ReadFile	fp(sl[i]);
	    if (fp.fd)	read_trimer(fp.fd, sl[i], i); else 
	    if (fp.gzfd) read_trimer(fp.gzfd, sl[i], i);
	    else	prompt(not_found, sl[i]);
	}
}

static	const	char*	rwdfq = 0;
static	const	char*	trim_out = 0;

void TriFreq::count_trimer(Seq& sd)
{
	Buf	buf;
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
	buf.many = sd.many;
const	char*	sn = sd.sqname();
const	char*	sl = strrchr(sn, '/');
	if (sl) sn = sl + 1;
	strcpy(buf.sname, sn);
	buf.sname[max_add_size - int_size] = 0;	//	 truncate
	PatMatHead	header = {0, int_size, nelm, 0, sd.len, nelmx3};
	header.add = CHAR(int_size + (strlen(buf.sname) + 3) / 4 * 4);

	WriteFile	tfp(trim_out, 0, true);
	bool	ok = false;
	if (tfp.gzfd)	ok = write_binary(tfp.gzfd, buf, header);
	else	fatal(no_file, trim_out);
	if (!ok) fatal(write_error, trim_out);
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

void TriFreq::background3(const char* wdf)
{
	Strlist	sl(wdf, " ,");
	for (INT i = 0; i < sl.size(); ++i) {
const	    char*	fname = sl[i];
	    Kmers	kmers(fname);
	    for (int j = 0; j <= morder; ++j) {
		switch (j) {
		    case 0: kmers.get(rcomp, j); break;
		    case 1: kmers.get(*rdifq, j); break;
		    case 2: kmers.get(**rtrif, j); break;
		}
	    }
	}
	ref_normalize();
}

template <typename file_t>
void TriFreq::rel_frequency(file_t ofd, const int n)
{
static	const	float	maxh = log(4.);
const	int	np1 = n + 1;
	float	h = 0;
	float*	lod = lods + n * lod_size;
	char	str[MAXL];
	for (int i = 0; i < nelm; ++i) {
	    float	r = scomp[i] / sumi;
	    if (entropy) h += (r > 0? r * log(r): 0);
	    else	*lod++ = r;
	}
	if (entropy && morder == 0) {
	    sprintf(str, "%d\t%9.5lf %9.5lf ", np1, -h, maxh + h);
	    fputs(str, ofd);
	}
	if (morder > 0) {
	    h = 0;
	    for (int i = 0; i < nelm; ++i)
		for (int j = 0; j < nelm; ++j) {
		    float	r = sdifq[i][j] / sumi;
		    if (entropy) h += (r > 0? r * log(r): 0);
		    else	*lod++ = r;
		}
	    if (entropy &&  morder == 1) {
		sprintf(str, "%d\t%9.5lf %9.5lf ", np1, -h, 2 * maxh + h);
		fputs(str, ofd);
	    }
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
	  if (entropy &&  morder == 2) {
		sprintf(str, "%d\t%9.5lf %9.5lf ", np1, -h, 3 * maxh + h);
		fputs(str, ofd);
	  }
	}
	if (entropy) fputc('\n', ofd);
}

void TriFreq::markovmodel(WriteFile& wfp, const int phase, const char* outf)
{
	if ((phase == 0 && !klds && !entropy)) return;
	if (wfp.gzfd) markovmodel(wfp.gzfd, phase);
	else if (wfp.fd) markovmodel(wfp.fd, phase);
	else	fatal(no_file, outf);
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
	    for (right = sites; --right > nmaxkld; )
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

void TriFreq::calculate_mean(TriFreq* ref)
{
const	bool	self = !ref;
	float*	rrcomp = 0;
	if (ref && rwdfq) {
	    ref->background3(rwdfq);
	    rrcomp = ref->rcomp;
	}
	if (self) ref = this;
	float*	lod = lods + left * lod_size;
	float*	led = lods + right * lod_size;
	float*	rod = self? lod: ref->lods;
	float*	red = self? led: ref->lods + ref->sites * ref->lod_size;
const	int	mm1 = std::min(morder, ref->morder) + 1;
	if (depsilon < 0) depsilon = 1./(float(nsupport));
	float	repsilon = 1./(float(ref->nsupport));

	float	scrf = 0;
	float	scrp = 0;
	for (; lod < led && rod < red; lod += lod_size, rod += ref->lod_size) {
	    for (int i = 0; i < bias[mm1]; ++i) {
		float	rf = rcomp[i] * (exp10((double) lod[i]) 
		    * (1. + depsilon) - depsilon);
		scrf += rf * rod[i];
		if (!rrcomp) continue;
		rf = rrcomp[i] * (exp10((double) rod[i]) 
		    * (1. + repsilon) - repsilon);
		scrp += rf * lod[i];
	    }
	}
	mmm.mean = scrf;
	if (self) return;

/***	output extimated score values	***/

	char	str[MAXL];
	strcpy(str, sname);
	char*	sl = strrchr(str, '/');
	sl = sl? sl + 1: str;
	char*	dot = strchr(sl, '.');
	if (dot)	{
		*dot = '\0';
		if (dot - sl > 8) sl[8] = '\0';
	}
	strcat(sl, ":");
	int	slen = strlen(sl);
	char*	ps = sl + slen;
	strcat(ps, ref->sname);
	sl = strrchr(ps, '/');
	sl = sl? sl + 1: ps;
	dot = strchr(sl, '.');
	if (dot)	{
	    *dot = '\0';
	    if (dot - sl > 8) sl[8] = '\0';
	}
	ps += strlen(ps);
	printf("%s\t%7.3f", str, scrf);
	if (rrcomp) printf("\t%7.3f", scrp);
	putchar('\n');
}

void TriFreq::output(WriteFile& wfp)
{
	if (left < 0) left = 0;
	if (right < 0) right = sites;
	else if (right > sites) {
	    prompt("Warning! output range is shrinked %d <- %d !\n", 
		sites, right);
	    right = sites;
	}
	if (binary) {
	    if (wfp.fd)	binary_output(wfp.fd);
	    else	binary_output(wfp.gzfd);
	} else {
	    if (wfp.fd)	text_output(wfp.fd);
	    else	text_output(wfp.gzfd);
	}
}

void PsFrqMat::name2id(const char* name)
{
	char	str[MAXL];
	strcpy(str, name);
	char*	sl = strrchr(str, '/');
	sl = sl? sl + 1: str;
	char*	dot = strchr(sl, '.');
	if (dot) *dot = '\0';
	strncpy(id, sl, ID_SIZE - 1);
}

void PsFrqMat::fuse(const char* fname)
{
	PsFrqMat	tail(fname);
	many += tail.many;
	float	sum = many + tail.many;
	for (int i = 0; i < nelm * cols; ++i) {
	    frq[i] = many * frq[i] + tail.many * tail.frq[i];
	    frq[i] /= sum;
	}
	many = (INT) (sum + 0.5);
}

void PsFrqMat::catenate(const char* fname)
{
	PsFrqMat	tail(fname);
	many = std::min(many, tail.many);
	vcopy(frq + nelm * cols, tail.frq, nelm * tail.cols);
	cols += tail.cols;
}

void PsFrqMat::output(WriteFile& wfp)
{
	if (binary) {
	    if (wfp.fd)	write_binary(wfp.fd);
	    else	write_binary(wfp.gzfd);
	} else {
	    if (wfp.fd)	write_text(wfp.fd);
	    else	write_text(wfp.gzfd);
	}
}

/****************************************************************************
*	Main
****************************************************************************/

#ifdef MAIN

static	bool	gzip = false;
static	const	char*	psfm_in = 0;

static void usage()
{
	fputs("Description:\n", stdout);
	fputs("\tGenerate N-th order position-specific score matrix ", stdout);
	fputs("(PSSM: Splice[3|5] or X.Nm[3|5])\n", stdout);
	fputs("\tfrom gap-less MSA around 3' (X.SP3) or 5' (X.SP5) splice junctions ", stdout);
	fputs("and backgroud\n", stdout);
	fputs("\tkmer-frequencies, ", stdout);
	fprintf(stdout, "where column size of MSA must be <= %d\n", max_cols);
	fputs("\t\tor\n", stdout);
	fputs("\tConvert between text and binary forms of PSSM files\n", stdout);
	fputs("Usage:\n", stdout);
	fputs("\tnpssm -m[0-2] -r X.wdfq [-g -b] PSSM X.SP[3|5]",stdout);
	fputs(" (MSA -> PSSM)\tor\n", stdout);
	fputs("\tnpssm -f PSSM[.dgz] [[-g -b| -o] PSSM] (text <-> binary)\n", stdout);
	fputs("Exapmle:\n", stdout);
	fputs("\tkmers -KD -w6 -d X_g -g -b X.wdfq\n", stdout);
	fputs("\tnpssm -m2 -r X.wdfq -g -b Splice3 X.SP3\n", stdout);
	fputs("\tnpssm -f X.2m3.dgz (read from gzipped binary, output ", stdout);
	fputs("text form to stdout)\n", stdout);
	fputs("\tnpssm -f X.2m3 -g -b X.02m[.dgz] (conversion from text to ", stdout);
	fputs("gzipped binary)\n", stdout);
	fputs("\tnpssm -p X.pfm.dat (convert binary to text format)\n", stdout);
	fputs("Options: (*: for private use only)\n", stdout);
	fputs("\t-b S\t: binary output file nam \n", stdout);
	fputs("\t-f S\t: existing PSSM w/wo extension\n", stdout);
	fputs("\t-g\t: gzipped output\n", stdout);
	fputs("\t-h\t: display this\n", stdout);
	fputs("\t-H F\t: minimum KL value (minH) for each column to be ", stdout);
	fputs("regarded as significant (0.025 nat)\n", stdout);
	fputs("\t-j N\t: signal position within MSA (inferred from MSA)\n", stdout);
	fputs("\t-l N\t: lower boundary within MSA (inferred from MSA and minH)\n", stdout);
	fputs("\t-m N\t: N-th order Markov model (N = [0|1|2])\n", stdout);
	fputs("\t-n* N\t: sample size\n", stdout);
	fputs("\t-o S\t: text output file name (to stdout if omitted)\n", stdout);
	fputs("\t-pq\t: suppress warning messages\n", stdout);
	fputs("\t-p* S\t: existing position-specific frequency matix\n", stdout);
	fputs("\t-q F\t: pseudocount (1./sample size)\n", stdout);
	fputs("\t-r S\t: background genomic k-mer frequencies w/wo extension\n", stdout);
	fputs("\t-t* S\t: trimer output\n", stdout);
	fputs("\t-u N\t: upper boundary within MSA (inferred from MSA and minH)\n", stdout);
	exit(1);
}

int main(int argc, const char** argv)
{
	int	mo = -1;
const	char*	outfn = 0;
const	char*	infn = 0;
const	char*	test_pssm = 0;
static	int	many_in = 0;
	bool	psfm_out = false;
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
		case 'f':	// existing pssm
		    if ((val = getarg(argc, argv))) {
			infn = pssm_in = val; psfm_in = 0;
		    }
		    break;
		case 'g': gzip = true; break;
		case 'h': usage();
		case 'H':
		    if ((val = getarg(argc, argv, true))) sig_thr = atof(val);
		    if (sig_thr > 1.) sig_thr /= 100.;
		    break;
		case 'i':
		    if ((val = getarg(argc, argv))) infn = trim_in = val;
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
		case 'n':
		    if ((val = getarg(argc, argv, true))) many_in = atoi(val);
		    break;
		case 'o':
		    if ((val = getarg(argc, argv))) outfn = val;
		    break;
		case 'p':	// existing position-specific freq mtx
		    if (!strcmp(argv[0], "-pq")) {
			setprompt(0, 0);
		    } else {
			if ((val = getarg(argc, argv))) {
			    infn = psfm_in = val; pssm_in = 0;
			}
			psfm_out = true;
		    }
		    break;
		case 'q':	// psuedocount
		    if ((val = getarg(argc, argv, true))) depsilon = atof(val);
		    break;
		case 'r':	// background k-mer frequency
		    if ((val = getarg(argc, argv))) wdfq = val;
		    break;
		case 'R':	// background k-mer frequency of 
		    if ((val = getarg(argc, argv))) rwdfq = val;
		    break;
		case 'F':
		    if ((val = getarg(argc, argv))) test_pssm = val;
		    break;
		case 't':
		    if ((val = getarg(argc, argv))) trim_out = val;
		    break;
		case 'u':	// upper bound
		    if ((val = getarg(argc, argv, true))) right = atoi(val);
		    break;
	    }
	}
	if (!outfn) {binary = 0; gzip = false;}
	WriteFile	wfp(outfn, binary? 0: 1, gzip, infn);
	if (pssm_in && test_pssm) {
	    tfq = new TriFreq(pssm_in);
	    TriFreq	tst_tfq(test_pssm);
	    tfq->calculate_mean(&tst_tfq);
	    return (0);
	} else if (pssm_in) {			// convert or synthesize existing pssm(s)
	    tfq = new TriFreq(pssm_in);
	    if (trim_in) tfq->read_trimer();
	    if (wdfq) {
		tfq->background3(wdfq);	// Get backgroud data
		tfq->setrange();
		tfq->markovmodel(wfp, 1, outfn);
	    }
	} else if (psfm_in) {		// Frequency matrix
	    PsFrqMat	psfm(psfm_in, many_in);
	    psfm.output(wfp);
	} else if (argc) {
	    Seq	sd(*argv);
	    if (trim_in) {	// read trq
		tfq = new TriFreq(mo, 0);
	    } else {		// calculate trq from MSA
		tfq = new TriFreq(mo, &sd);
		tfq->count_trimer(sd);
	    }
	    if (mo < 0 && trim_out) {
		delete tfq;
		return (0);
	    }
	    if (!wdfq) {
		delete tfq;
		fatal("Set -r wdfq option !\n");
	    }
	    if (depsilon < 0) depsilon = 1./((float) sd.many);
	    tfq->background3(wdfq);		// Get backgroud data
	    tfq->markovmodel(wfp, 0, outfn);
	    if (!entropy && wdfq) {
		tfq->setrange();
		tfq->markovmodel(wfp, 1, outfn);
		tfq->calculate_mmm(sd);	// min, mean, max
	    }
	} else if (trim_in) {
	    tfq = new TriFreq(mo, 0);
	    if (depsilon < 0) depsilon = 1./((float) tfq->nsupport);
	    if (wdfq) {
		tfq->background3(wdfq);		// Get backgroud data
		tfq->markovmodel(wfp, 0, outfn);
		if (!entropy && wdfq) {
		    tfq->setrange();
		    tfq->markovmodel(wfp, 1, outfn);
		    tfq->calculate_mean();
		}
	    } else {
		tfq->markovmodel(wfp, 0, outfn);
	    }
	} else	usage();

	if (!entropy && tfq) {
	    if (psfm_out) {
		PsFrqMat	psfm(tfq);
		psfm.output(wfp);
	    } else
		tfq->output(wfp);
	}
	delete tfq;
	return (0);
}
#endif	// MAIN
