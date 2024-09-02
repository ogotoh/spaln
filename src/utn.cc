/*****************************************************************************
*
*	Utilities for nucleotide/amino-acid sequence management
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
#if USE_WEIGHT
#include "mseq.h"
#endif
#include "autocomp.h"
#include "pattern.h"
#include "resite.h"
#include <math.h>
#include <time.h>

#define	lod(n, d, r) (log10(((d) < Epsilon? Epsilon: ((n)/(d) + Epsilon)) / ((r) + Epsilon)))

class Composit {
	double*	phase[3];
	int*	gccmp[3];
	int* 	codontab;
	int	TotalNseq;
	int	TotalNres;
	int	nsize;
	int	nsize2;
	int	molc;
	SEQ_CODE*	sqcode;
public:
	bool	each;
	Composit(Seq* sd);
	~Composit() {delete[] phase[0]; delete[] gccmp[0]; delete[] codontab;}
	void	calculate(Seq* sd);
	void	codonuse(Seq* sd);
	void	monocomp(Seq* sd);
	void	monotron(Seq* sd);
	void	ditronfrq(Seq* sd);
	void	csvtronfrq(Seq* sd);
	void	monocodon(Seq* sd);
	void	dicodonfrq(Seq* sd);
	void	csvcodonfrq(Seq* sd);
	void	dinucfrq(Seq* sd);
};

static	char	spat[] = "Pattern : ";
static	char	sezm[] = "Enzymes : ";
static	const	char*	reffile = 0;
static	char	pattern[MAXL];

extern	void	extractcds(Seq* sd);
static	void	aa_code();
static	void	poly_a(Seq* sd);
static	void	Proutseq(FILE* fd, Seq* sd);
static	void	listseq(Seq* ns, char* sub);
static	int	icmpf(int* a, int* b);
static	void	mute(Seq* seqs[]);
static	void	muteseq(FILE* fd, Seq* seqs[]);
static	void	mapsite(Seq* sd, int s2n);
static	double*	readwdfq(int kk, int nelm);
static	double*	DiNucFrq(double* frq, Seq* sd);
static	int	putorf(FILE* fd, Seq* ns);
static	void	openfrm(Seq* ns);
static	int	printorf(FILE* fd, Seq* ns);
static	int	transorf(FILE* fd, Seq* ns);
static	void	codseq(Seq* ns);
static	void	transl(Seq* ns);
static	void	compstn(Seq* sd);
static	void	printnbr(Seq* sd);
static	char*	pgetstr(char* pat, char* prmpt);
static	void	allezm(FILE* fd, Seq* sd, char* pat);
static	void	find(Seq* ns);
static	void	resezm(Seq* ns);
static	void	nameseq(FILE* fd, Seq* ns);
static	Seq*	forgeseq(FILE* fd, Seq* sd);
static	void	shufseq(Seq* sd);

static	int	trl_from = 0;
static	int	trl_to = 0;
static	int	longestorf = 0;
static	int	q_mns = 1;
static	double	psubst = 0.1;
static	double	pindel = 0.1;
static	double	pbase = 0.6;
static	long	randseed = 1984781;
static	double	Epsilon = 1.e-5;
static	char	bad_ref[] = "Improper wdfq file: %s";

void usage()
{
	fputs("Specific Options:\n", stderr);
	fputs("\t-c:\t Composition\n", stderr);
	fputs("\t-f[pattern]:\t Find pattern\n", stderr);
	fputs("\t-l:\t Print sequence\n", stderr);
	fputs("\t-L:\t Print pieces of sequence\n", stderr);
	fputs("\t-n:\t Name and range\n", stderr);
	fputs("\t-o:\t Open reading frames\n", stderr);
	fputs("\t-r[enzyme]:\t Restriction sites\n", stderr);
	fputs("\t-t:\t Translate the longest orf\n", stderr);
	fputs("\t-T:\t Translate the longest orf with comment\n", stderr);
	fputs("\t-p:\t Set parameters\n", stderr);
	fputs("\t-Gn:\t Genetic Code [0|1|2|3|4: univ|mt|Dme|yeast|cal]\n", stderr);
	fputs("\t-In:\t Min ORF length\n", stderr);
	fputs("\t-Jn:\t From IC\n", stderr);
	exit (0);
}

template <class seq_t>
int AlnServer<seq_t>::localoption(int& argc, const char**& argv)
{
static	const	char*   Localcmds = "QZacfgkilnortx";
	const	char*	opt = argv[0] + 1;
	const	char*	val = argv[0] + 2;
	int	c = *opt;
	int	n = 1;

	if (!c) return (0);
	const	char*	al = strchr(Localcmds, c);
	if (!al) {		// sub-options
	  switch (c) {
	    case 'E': algmode.mlt = ((val = getarg(argc, argv, true)))?
		atoi(val): 2; break;
	    case 'G': 
		if (this->jobcode == 'x') {
		    pbase = atof(val);
		    if (pindel >= 1.) pindel /= 100.;
		} else
		    setcodon(atoi(val));
		break;
	    case 'I':
		if (this->jobcode == 'x') {
		    pindel = atof(val);
		    if (pindel >= 1.) pindel /= 100.;
		} else
		    setorf(SILENT, atoi(val));
		break;
	    case 'J': val = getarg(argc, argv, true);
		if (val) setorf(atoi(val), SILENT);
		break;
	    case 'M': 
		if (this->jobcode == 'x') {
		    psubst = atof(val);
		    if (psubst >= 1.) psubst /= 100.;
		}
		break;
	    case 'R':
		if (this->jobcode == 'x') randseed = atoi(val);
		else if (!*val && argc - 1 && argv[1][0] != '-') {
		    reffile = *++argv; argc--;
		} else
		    reffile = val;
		break;
	    case 'p':
		if (*val == 'a') ++val;
		if (isdigit(val[1])) ++val;
		else if (!val[1] && isdigit(argv[1][0])) {
		    val = *++argv; --argc;
		} else val = 0;
		if (val) polyA.setthr(val);
		else	polyA.setthr(def_polya_thr);
		break;
	    default:  n = 0; break;
	  }
	} else if (!this->jobcode) {	// prime option
	    this->jobcode = c;
	    c = argv[1]? argv[1][0]: 0;
	    if (this->jobcode == 'f' || this->jobcode == 'r') {
		if (*val)	strcpy(pattern, val);
		else if (c && c != OPTCHAR) {
		    ++argv; --argc;
		    strcpy(pattern, *argv);
		} else if (!pgetstr(pattern, this->jobcode == 'f'? spat: sezm)) {
		    exit (1);
		}
	    } else if (this->jobcode == 'o' || this->jobcode == 'O') {
		if (isdigit(c))	longestorf = atoi(val);
	    }
	    if (const char* col = strchr(val, ':')) this->set_catalog(col + 1);
	    if (this->jobcode == 'l') n = 0;
	} else {
	  switch (c) {
	    case 'o':
		longestorf = 1; break;
	    default: n = 0; break;
	  }
	}
	return (n);
}

template <class seq_t>
void AlnServer<seq_t>::setparam(int lvl)
{
	if (lvl & 1) {
	    setcodon(QUERY);
	    setorf(QUERY, QUERY);
	}
	if (lvl & 2) {
	    setlpw(QUERY);
	    setform(QUERY);
	    setprmode(QUERY, QUERY, QUERY);
	}
}

static void aa_code()
{
	char	ns[MAXL];

	while (true) {
	    fputs(": ", stderr);
	    char*	ps = fgets(ns, MAXL, stdin);
	    if (!ps || !*ps) return;
	    int	i = toaa(tosqcode((CHAR*) ns, setSeqCode(0, DNA)));
	    fprintf(stderr, "%s [%c] %2d\n", amino3[i], amino[i], i);
	}
}

static void poly_a(Seq* sd)
{
	if (!sd) return;
	int	polya = 0;
	int	polyt = 0;
	int	isa = 0;
	int	ist = 0;
	int	scr = 0;
	int	cnt = 0;
	int	thr = polyA.getthr() * sd->many;
	CHAR*	maxa = 0;
	CHAR*	maxt = 0;
	CHAR*	s = sd->at(sd->right);
	CHAR*	b = sd->at(sd->left);
	CHAR*	t = s;

	while (--s >= b) {
	    scr += (*s == A? 3: -6);
	    if (*s == A) cnt++;
	    if (scr > polya) {
		polya = scr;
		if (scr > thr) {
		    maxa = s;
		    isa = cnt;
		}
	    }
	    if (scr < polya - thr) break;
	}
	s = sd->at(0);
	for (scr = cnt = 0; s < t; ++s) {
	    scr += (*s == T? 3: -6);
	    if (*s == T) cnt++;
	    if (scr > polyt) {
		polyt = scr;
		if (scr > thr) {
		    maxt = s;
		    ist = cnt;
		}
	    }
	    if (scr < polyt - thr) break;
	}
	if (maxa && maxt) {
	    if (polya >= polyt) maxt = 0;
	    else		maxa = 0;
	}
	if (maxa) {
	    thr = sd->index(maxa);
	    printf("PolyA: %s\t%d\t%d\t%d\t%d\n",
		sd->sqname(), sd->len, thr + 1, isa, int(t - maxa));
	} else if (maxt) {
	    thr = sd->index(maxt) + 1;
	    printf("PolyT: %s\t%d\t%d\t%d\t%d\n",
		 sd->sqname(), sd->len, thr, ist, sd->index(maxt + 1));
	} else {
	    printf("NoTail: %s\t%d\n", sd->sqname(), sd->len);
	}
}

static void Proutseq(FILE* fd, Seq* sd)
{
	if (OutPrm.lpw > 8)  sd->typeseq(fd);
	else {
	    algmode.nsa = OutPrm.lpw;
	    sd->fpmem_len(fd);
	}
}

static void listseq(Seq* ns, char* sub)
{
	RANGE	temp;
	int 	from = ns->left + 1;
	Seq*	sd = 0;

	ns->saverange(&temp);
	int	subcmd = atoi(sub);
	if (subcmd > 0) {
	    while (true) {
		promptin("from (%d) - to (%d): ", &from, &ns->right);
		if (from <= 0) break;
		ns->left = from - 1;
		if (subcmd == 1) break;
		sd = ns->catseq(sd);
		from = -1;
		ns->right = temp.right;
	    }
	}
	if (from <= -subcmd) return;
	FILE*	fd = qout("");
	if (fd) {
	    if (sd) {
		Proutseq(fd, sd);
		delete sd;
	    } else {
		Proutseq(fd, ns);
	    }
	    qclose();
	}
	ns->restrange(&temp);
}

static int icmpf(int* a, int* b)
{
	return (a - b);
}

typedef struct {int pos, len;} INDEL;

static int ildcmpf(INDEL*a, INDEL*b);
static int ildcmpf(INDEL*a, INDEL*b)
{
	return (a->pos - b->pos);
}

static void mute(Seq* seqs[])
{
	double	lpbase = log(pbase);
	Seq*	sd = seqs[0];
	CHAR*	buf = new CHAR[sd->len + 2] + 1;
static	const	CHAR	unredctab[4] = {A, C, G, T};

	memcpy(buf - 1, sd->at(-1), sd->len + 2);
	srand48(time(0));
	int	m = (int) (sd->len * psubst);
	INDEL*	idl = new INDEL[m];
	INDEL*	wid = idl;
	int	dif = 0;
	for (int n = 0; n < m; ++n) {
	    int	pos = (int) (sd->len * drand48());
	    CHAR*	s = buf + pos;
	    double	p = drand48();
	    if (p < pindel) {	/* indel */
		if (pos < 5 || sd->len - pos < 5) continue;
		wid->pos = pos;
		wid->len = (int) (log(p / pindel) / lpbase) + 1;
		int	mut = (int) (p * 10000) % 2;
		if (mut) wid->len = -wid->len;
		dif += wid->len;
		wid++;
	    } else {			/* substitution */
		int	mut = (int) ((p - pindel) * 9999) % 3 + 1;
		mut += ncredctab[*s];
		*s = unredctab[mut % 4];
	    }
	}
	wid->pos = sd->len;
	wid->len = 0;
	dif += sd->len;
	qsort((UPTR) idl, ++wid - idl, sizeof(INDEL), (CMPF) ildcmpf);
	seqs[1]->refresh(1, dif + 2);
	seqs[0]->copyattr(seqs[1]);
	seqs[0]->copylbl(seqs[1]);
	CHAR*	t = seqs[1]->at(0);
	CHAR*	s = buf;
	memcpy(t, s, idl->pos);
	t += idl->pos;
	s += idl->pos;
	for (wid = idl; wid->len; ) {
	    if (wid->len > 0) {	/* insert */
		for (int m = 0; m < wid->len; ++m) {
		    int	n = (int) (sd->len * drand48());
		    *t++ = *(sd->at(n));
		}
	    } else			/* delete */
		s -= wid->len;
	    CHAR*	u = buf + (++wid)->pos;
	    while (s < u) *t++ = *s++;
	}
	seqs[1]->postseq(t);
	delete[] idl;
	delete[] (buf - 1);
}

static void muteseq(FILE* fd, Seq* seqs[])
{
	mute(seqs);
	Proutseq(stdout, seqs[1]);
}

static int isreverse(Seq* sd, int mem)
{
	return (sd->inex.sens & RV);
}

static void mapsite(Seq* sd, int s2n)
{
	int	positions[MAXLOC];
	int*	pos = positions + 1;

	while (true) {
	    int	n = fgetiarray(positions, MAXLOC, stdin);
	    if (n == 0 || *pos == 0) return;
	    int	mem = positions[0] - 1;
	    int	rev = isreverse(sd, mem);
	    qsort((UPTR) pos, (INT) (n-1), sizeof(int), (CMPF) icmpf);
	    if (s2n) {
		sd->pos2num(mem, pos);
		if (sd->nbr && sd->nbr[mem]) {
		    for (int i = 1; i < n; ++i) {
			positions[i] = rev?
			sd->nbr[mem] - positions[i] + 2:
			sd->nbr[mem] + positions[i];
		    }
		}
	    } else
		sd->num2pos(mem, pos);
	    printf("%16s ", (*sd->sname)[mem]);
	    for (int i = 1; i < n; ++i) {
		printf("%d ", positions[i]);
	    }
	    putchar('\n');
	}
}

static double* readwdfq(int kk, int nelm)
{
	FILE*	fd;
	char	str[MAXL];
	int	nn = nelm;
	const	char*	rf = reffile;

	if (!rf) {
	    progets(str, ".wdfq file ? ");
	    rf = strrealloc(0, str);
	}
	if (!rf || !(fd = fopen(rf, "r"))) fatal(not_found, rf);
	for (int k = 1; k < kk; ++k, nn *= nelm) {	/* skip lines */
	    for (int i = 0; i < nn; ++i)
		if (!fgets(str, MAXL, fd)) fatal(bad_ref, rf);
	}
	double*	RefFreq = new double[nn];
	double*	dd = RefFreq;
	for (int i = 0; i < nn; ++i) {
	    if (!fgets(str, MAXL, fd)) fatal(bad_ref, rf);
	    *dd++ = atof(cdr(str));
	}
	fclose(fd);
	if (!reffile) delete[] rf;
	return (RefFreq);
}

static int putorf(FILE* fd, Seq* ns)
{
static	const	int	cmplphs[3] = {0, 2, 1};
static	const	char	frmt[] = "%-14s %5d %5d; %5d nt, %4d aa, frame [%d], 1stc %5d, phase [%d]\n";

	if (!ns->isdrna()) return(ERROR);
	ORF*	orfs = ns->getorf();
	if (!orfs) return (0);
	ORF*	wk = orfs;
	for (int n = 0; wk->len; ++wk) {
	    int	phs = (wk->pos - wk->frm + 3) % 3;
	    int	fstc = wk->pos + cmplphs[phs];
	    int	frm = (ns->inex.sens == COMREV)?
		(ns->len - wk->frm + 2) % 3: wk->frm;
	    fprintf(fd, frmt, ns->sqname(), ns->SiteNo(wk->pos), 
		ns->SiteNo(wk->pos + wk->len - 1), wk->len, wk->len / 3, 
		frm + 1, ns->SiteNo(fstc), phs);
	    if (longestorf && ++n >= longestorf) break;
	}
	int	norf = wk - orfs;
	if (norf > 0) {
	    trl_from = orfs->pos;
	    trl_to = orfs->pos + orfs->len + 3;
	}
	delete[] orfs;
	return (norf);
}
	
static void openfrm(Seq* ns)
{
	FILE*	fd = qout("");

	if (!fd) return;
	putorf(fd, ns);
	qclose();
}

static int printorf(FILE* fd, Seq* ns)
{
	ORF*	orfs = ns->getorf();

	if (!orfs) return (ERROR);
	RANGE	tmp;
	ns->saverange(&tmp);
	if (algmode.nsa <= 1) orfs[1].len = 0;
	int	orfn = orfs[1].len? 1: 0;
	ORF*	wk = orfs;
	for ( ; wk->len; ++wk) {
	     ns->left = wk->pos;
	     ns->right = wk->pos + wk->len + 3;
	     Proutseq(fd, ns);
	}
	orfn = wk - orfs;
	delete[] orfs;
	ns->restrange(&tmp);
	return (orfn);
}

static int transorf(FILE* fd, Seq* ns)
{
	ORF*	orfs = ns->getorf();
	if (!orfs) return (ERROR);
	RANGE	tmp;
	ns->saverange(&tmp);
	if (algmode.nsa <= 1) orfs[1].len = 0;
	int	orfn = orfs[1].len? 1: 0;
	ORF*	wk = orfs;
	for ( ; wk->len; ++wk) {
	    ns->left = wk->pos;
//	    if (wk->frm) ns->left += 3 - wk->frm;
	    ns->right = wk->pos + wk->len;
	    ns->transout(fd, STOP, orfn++);
	}
	orfn = wk - orfs;
	delete[] orfs;
	ns->restrange(&tmp);
	return (orfn);
}

static void codseq(Seq* ns)
{
	if (trl_from >= trl_to) {
	    FILE*	fd = qout("");
	    printorf(fd, ns);
	    qclose();
	    return;
	}
	if (trl_from < ns->left) trl_from = ns->left;
	++trl_from;
	int 	to = trl_to;
	promptin("from (%d) - to (%d) : ", &trl_from, &to);
	FILE*	fd = qout("");
	if (!fd) return;
	RANGE	tmp;
	ns->saverange(&tmp);
	ns->left  = --trl_from;
	ns->right = min(ns->right, to);
	Proutseq(fd, ns);
	qclose();
	ns->restrange(&tmp);
}

static void transl(Seq* ns)
{
	RANGE	tmp;
	FILE*	fd;
	int 	to = STOP;

	if (trl_from < ns->left) trl_from = ns->left;
	++trl_from;
	promptin("from (%d) - to (%d) : ", &trl_from, &to);
	fd = qout("");
	if (!fd) return;
	ns->saverange(&tmp);
	ns->left  = --trl_from;
	ns->right = min(ns->right, to);
	ns->transout(fd, to, 0);
	qclose();
	ns->restrange(&tmp);
}

void Composit::codonuse(Seq* sd)
{
static	const	int ctab[] = {T, C, A, G};	/* U, C, A, G */
	int 	gc;

	if (sd) {
	    if (!codontab) {
		codontab = new int[64];
		clearmem(codontab, 64 * sizeof(int));
	    } else if (each)
		clearmem(codontab, 64 * sizeof(int));
	    CHAR*	s = sd->at(sd->left);
	    for (int i = sd->left; i < sd->right; i += 3, s += sd->many * 2)
		for (int j = 0; j < sd->many; ++j, ++s) 
		    if ((gc = codon_id(s, sd->many)) != ERROR) codontab[gc]++;
	    if (!each) return;
	}

//	putc('\n', out_fd);
	CHAR	trip[3];
	if (algmode.nsa == 0 || algmode.nsa == 2) {
	    for (int i = 0; i < 4; i++) {
		trip[0] = red2nuc[i];
		for (int j = 0; j < 4; j++) {
		    trip[1] = red2nuc[j];
		    for (int k = 0; k < 4; k++) {
			trip[2] = red2nuc[k];
			gc = codon_id(trip, 1);
			if (algmode.nsa == 2) {
			    for (int n = 0; n < 3; )
				putc(nucl[trip[n++]], out_fd);
			    fprintf(out_fd, "\t%s\t%6d\t%7.4f\n",  amino3[toaa(trip)], 
				codontab[gc], 100. *  codontab[gc] / TotalNres);
			} else {
			    fprintf(out_fd, "\t%14.7e", double (codontab[gc]) / TotalNres);
			}
		    }
		    if (algmode.nsa == 0) fputc('\n', out_fd);
		}
	    }
	    return;
	}
	for (int i = 0; i < 4; i++) {
	    trip[0] = ctab[i];
	    for (int j = 0; j < 4; j++) {
		trip[2] = ctab[j];
		for (int k = 0; k < 4; k++) {
		    trip[1] = ctab[k];
		    gc = codon_id(trip, 1);
		    fputs("    ", out_fd);
		    for (int n = 0; n < 3; )
			putc(nucl[trip[n++]], out_fd);
		    fprintf(out_fd, "  %s %6d", amino3[toaa(trip)], codontab[gc]);
		}
		putc('\n', out_fd);
	    }
	    putc('\n', out_fd);
	}
}

void Composit::monocomp(Seq* sd)
{
	if (sd) {
	    if (!phase[0]) {
		sqcode = sd->code;
		nsize = sqcode->max_code;
		phase[0] = new double[3 * nsize];
		for (int i = 1; i < 3; ++i) phase[i] = phase[i-1] + nsize;
		vclear(phase[0], 3 * nsize);
	    } else if (each)
		vclear(phase[0], 3 * nsize);
	    CHAR*	ss = sd->at(sd->left);
	    CHAR*	tt = sd->at(sd->right);

	    FTYPE*	wt = 0;
#if USE_WEIGHT
	    wt = sd->weight;
#endif
	    for (int n = 0; ss < tt; ++n) {
		FTYPE*	w = wt;
		for (int i = 0; i < sd->many; ++i) {
		    if (wt)	phase[n % 3][*ss++] += *w++;
		    else	phase[n % 3][*ss++] += 1.;
		}
	    }
	    if (!each) return;
	}

	double	gc[4];
	double	ff = gc[3] = 0;
	double	ag = 0.;
	for (int n = 0; n < 3; ++n) {
	    double	f = phase[n][A] + phase[n][C] + phase[n][G] + phase[n][T];
	    ff += f;
	    gc[n] =  phase[n][C] + phase[n][G];
	    gc[3] += gc[n];
	    if (f > 0.) gc[n] /= f;
	    ag += phase[n][A] + phase[n][G];
	}
	if (ff != 0.) {
	    gc[3] /= ff;
	    ag /= ff;
	}
	if (!each)
	    fprintf(out_fd, " NoSeq = %d TotalNcl = %d %d\n", TotalNseq, int(ff), TotalNres);
	if (algmode.alg == 0) {
	    fprintf(out_fd, "%5.2f", 100 * gc[3]);
	    if (longestorf)
		fprintf(out_fd, "  %5.2f  %5.2f  %5.2f", 100 * gc[0], 100 * gc[1], 100 * gc[2]);
	    else
		fprintf(out_fd, "  %5.2f", 100. * ag);
	    fputc('\n', out_fd);
	} else {
	    for (int i = sqcode->base_code; i < sqcode->max_code; ++i) {
		double	ff = 0.;
		for (int n = 0; n < 3; ++n) ff += phase[n][i];
		if (ff == 0.) continue;
		if (algmode.nsa) {
		    putc(sqcode->decode[i], out_fd);
		    putc('\t', out_fd);
		    for (int n = 0; n < 3; ++n) 
			fprintf(out_fd, "%10.2f\t", phase[n][i]);
		    fprintf(out_fd, "%10.2f\n", ff);
		} else {
		    fprintf(out_fd, "\t%c:%d", sqcode->decode[i], int(ff));
		}
	    }
	    if (!algmode.nsa) putc('\n', out_fd);
	}
}

void Composit::monotron(Seq* sd)
{
	if (sd) {
	    if (!phase[0]) {
		sqcode = sd->code;
		nsize = sqcode->max_code - sqcode->base_code;
	        phase[0] = new double[3 * nsize];
	        for (int i = 1; i < 3; ++i) phase[i] = phase[i-1] + nsize;
		vclear(phase[0], 3 * nsize);
	    } else if (each)
		vclear(phase[0], 3 * nsize);

	    CHAR*	ss = sd->at(sd->left + 1);
	    CHAR*	tt = sd->at(sd->right - 1);

	    FTYPE*	wt = 0;
#if USE_WEIGHT
	    wt = sd->weight;
#endif
	    for (int n = 0; ss < tt; ++n) {
		FTYPE*	w = wt;
		for (int i = 0; i < sd->many; ++i) {
		    int	d = *ss++ - sqcode->base_code;
		    if (d >= 0 && d < sqcode->max_code) {
			if (wt) phase[n % 3][d] += *w++;
			else    phase[n % 3][d] += 1.;
		    } else if (wt) w++;
		}
	    }
	    if (!each) return;
	}

	for (int i = 0; i < nsize; ++i) {
	    double	ff = 0.;
	    for (int n = 0; n < 3; ++n) ff += phase[n][i];
	    if (ff == 0.) continue;
	    putc(sqcode->decode[i + sqcode->base_code], out_fd);
	    putc('\t', out_fd);
	    for (int n = 0; n < 3; ++n) 
		fprintf(out_fd, "%10.2f\t", phase[n][i]);
	    fprintf(out_fd, "%10.2f\n", ff);
	}
}

void Composit::ditronfrq(Seq* sd)
{
	if (sd) {
	    if (!phase[0]) {
		sqcode = sd->code;
		nsize = sqcode->max_code - sqcode->base_code;
		nsize2 = nsize * nsize;
		phase[0] = new double[3 * nsize2];
		for (int i = 1; i < 3; ++i) phase[i] = phase[i-1] + nsize2;
		vclear(phase[0], 3 * nsize2);
	    } else if (each)
		vclear(phase[0], 3 * nsize2);
	    CHAR*	ss = sd->at(sd->left + 4);
	    CHAR*	tt = sd->at(sd->right - 4);
	    CHAR*	uu = ss + sd->many * 3; 

	    FTYPE*	wt = 0;
#if USE_WEIGHT
	    wt = sd->weight;
#endif
	    for (int n = 0; ss < tt; ++n) {
		int	p = n % 3;
		FTYPE*	w = wt;
		for (int i = 0; i < sd->many; ++i) {
		    int	d  = *ss++ - sqcode->base_code;
		    int	dd = *uu++ - sqcode->base_code;
		    if (d >= 0 && d < nsize && dd >= 0 && dd < nsize) {
			d += dd;
			if (wt) phase[p][d] += *w++;
		        else	phase[p][d] += 1.;
		    } else if (wt) w++;
		}
	    }
	    if (!each) return;
	}

	double	denom[3];
	if (algmode.nsa > 0) {	/* ditron Model */
	    double*	RefFreq = readwdfq(2, nsize);
	    fprintf(out_fd, "CodePotTab 3 %d 1 1 %d %d\n", 
		nsize2, TotalNseq, TotalNres);
	    for (int i = 0, d = 0; i < nsize; ++i) {
	      for (int n = 0; n < 3; ++n) denom[n] = 0;
	      double	ff = 0.;
	      int	dd = d;
	      for (int j = 0; j < nsize; ++j, ++d) {
		ff += RefFreq[d];
		for (int n = 0; n < 3; ++n) denom[n] += phase[n][d];
	      }
	      d = dd;
	      for (int j = 0; j < nsize; ++j, ++d) {
		double	rr = RefFreq[d] / ff;
		fprintf(out_fd, "%10.5f\t%10.5f\t%10.5f\n", 
		  lod(phase[0][d], denom[0], rr),
		  lod(phase[1][d], denom[1], rr),
		  lod(phase[2][d], denom[2], rr));
	      }
	    }
	    delete[] RefFreq;
	} else {
	    int	d = 0;
	    for (int i = 0; i < nsize; ++i) {
	      for (int j = 0; j < nsize; ++j, ++d) {
		double	ff = 0.;
		for (int n = 0; n < 3; ++n) ff += phase[n][d];
		if (ff == 0.) continue;
		putc(sqcode->decode[sqcode->base_code + i], out_fd);
		putc(sqcode->decode[sqcode->base_code + j], out_fd);
		putc('\t', out_fd);
		for (int n = 0; n < 3; ++n) 
		    fprintf(out_fd, "%10.2f\t", phase[n][d]);
		fprintf(out_fd, "%10.2f\n", ff);
	      }
	    }
	}
}

void Composit::csvtronfrq(Seq* sd)
{
	if (sd) {
	    if (!phase[0]) {
		sqcode = sd->code;
		nsize = sqcode->max_code - sqcode->base_code;
		nsize2 = nsize * nsize;
		phase[0] = new double[3 * nsize2];
		for (int i = 1; i < 3; ++i) phase[i] = phase[i-1] + nsize2;
		vclear(phase[0], 3 * nsize2);
	    } else if (each)
		vclear(phase[0], 3 * nsize2);
	    CHAR*	ss = sd->at(sd->left + 1);
	    CHAR*	tt = sd->at(sd->right - 1);

	    FTYPE*	wt = 0;
#if USE_WEIGHT
	    wt = sd->weight;
#endif
	    for (int n = 0; ss < tt; ++n) {
		int	b = *ss - sqcode->base_code;
		for (int i = 0; i < sd->many; ++i, ++ss) {
		    FTYPE*	w = wt;
		    if (i == 0) {	/* base */
			if (b >= 0 && b < sqcode->max_code) b *= nsize;
			else {
			    ss += sd->many;
			    break;
			}
		    } else {	/* derived */
			int	d = *ss - sqcode->base_code;
			if (d >= 0 && d < sqcode->max_code) {
			    if (wt)	phase[n % 3][b + d] += *w++;
			    else	phase[n % 3][b + d] += 1.;
			} else if (wt) w++;
		    }
		}
	    }
	    if (!each) return;
	}

	for (int i = 0, d = 0; i < nsize; ++i) {
	    for (int j = 0; j < nsize; ++j, ++d) {
		double	ff = 0.;
		for (int n = 0; n < 3; ++n) ff += phase[n][d];
		if (ff == 0.) continue;
		putc(sqcode->decode[i + sqcode->base_code], out_fd);
		putc(sqcode->decode[j + sqcode->base_code], out_fd);
		putc('\t', out_fd);
		for (int n = 0; n < 3; ++n) 
		    fprintf(out_fd, "%10.2f\t", phase[n][d]);
		fprintf(out_fd, "%10.2f\n", ff);
	    }
	}
}

void Composit::monocodon(Seq* sd)
{
	if (sd) {
	    nsize = 64;
	    if (!phase[0]) {
		phase[0] = new double[3 * nsize];
		for (int i = 1; i < 3; ++i) phase[i] = phase[i - 1] + nsize;
		vclear(phase[0], 3 * nsize);
	    } else if (each)
		vclear(phase[0], 3 * nsize);
	    CHAR*	ss = sd->at(sd->left);
	    CHAR*	tt = sd->at(sd->right - 2);

	    FTYPE*	wt = 0;
#if USE_WEIGHT
	    wt = sd->weight;
#endif
	    for (int n = 0; ss < tt; ++n) {
		for (int i = 0; i < sd->many; ++i) {
		    FTYPE*	w = wt;
		    int	gc = codon_id(ss++, sd->many);
		    if (gc != ERROR) {
			if (wt) phase[n % 3][gc] += *w++;
			else    phase[n % 3][gc] += 1.;
		    } else if (wt) w++;
		}
	    }
	    if (!each) return;
	}

	for (int i = 0, gc = 0; i < 4; ++i) {
	    for (int j = 0; j < 4; ++j) {
	      for (int k = 0; k < 4; ++k, ++gc) {
		double	ff = 0.;
		for (int n = 0; n < 3; ++n) ff += phase[n][gc];
		if (ff == 0.) continue;
		fprintf(out_fd, "%c%c%c\t", Nucl[i], Nucl[j], Nucl[k]);
		for (int n = 0; n < 3; ++n) 
		    fprintf(out_fd, "%10.2f\t", phase[n][gc]);
		fprintf(out_fd, "%10.2f\n", ff);
	      }
	    }
	}
}

void Composit::dicodonfrq(Seq* sd)
{
	if (sd) {
	    if (!phase[0]) {
		nsize = 64;
		nsize2 = 4096;	/* 64 * 64 */
		phase[0] = new double[3 * nsize2];
		for (int i = 1; i < 3; ++i) 
		    phase[i] = phase[i-1] + nsize2;
		vclear(phase[0], 3 * nsize2);
	    } else if (each) {
		vclear(phase[0], 3 * nsize2);
	    }
	    delete[] gccmp[0];
	    gccmp[0] = new int[3 * sd->many];
	    for (int i = 1; i < 3; ++i)
		gccmp[i] = gccmp[i-1] + sd->many;
	    vclear(gccmp[0], 3 * sd->many);

	    CHAR*	ss = sd->at(sd->left);
	    CHAR*	tt = sd->at(sd->right - 5);
	    int*	wgc = gccmp[0];
	    for (int i = 0; i < 3 * sd->many; ++i)	// first codon
		*wgc++ = codon_id(ss++, sd->many);
	    FTYPE*	wt = 0;
#if USE_WEIGHT
	    wt = sd->weight;
#endif
	    for (int n = 0; ss < tt; ++n) {
		int	p = n % 3;
		FTYPE*	w = wt;
		for (int i = 0; i < sd->many; ++i) {
		    int	b = codon_id(ss++, sd->many);
		    if (b != ERROR && gccmp[p][i] != ERROR) {
			int	d = nsize * gccmp[p][i] + b;
			if (wt) phase[p][d] += *w++;
			else	phase[p][d] += 1.;
		    } else if (wt) w++;
		    gccmp[p][i] = b;
		}
	    }
	    if (!each) return;
	}

/*	algmode.nsa:: 0: count, 1: 5th MM, 2; 4th MM, 3: dicodon	*/

	if (algmode.nsa > 0) {
	    int	nsize0 = 4;
	    double	denom[3];
	    double*	RefFreq = readwdfq(6, 4);
	    int	i = 1;
	    if (algmode.nsa < 6) 
	      for ( ; i < (int) algmode.nsa; ++i) nsize0 *= 4;
	    int	nsize1 = nsize2 / nsize0;
	    fprintf(out_fd, "CodePotTab 4 %d %d %d %d %d\n", 
		nsize2, 6 - i, i, TotalNseq, TotalNres);
	    for (int d = i = 0; i < nsize1; ++i) {
	      for (int p = 0; p < 3; ++p) denom[p] = 0;
	      double	ff = 0.;
	      int	dd = d;
	      for (int j = 0; j < nsize0; ++j, ++d) {
		ff += RefFreq[d];
		for (int p = 0; p < 3; ++p) denom[p] += phase[p][d];
	      }
	      double	sumd = 0.;
	      for (int p = 0; p < 3; ++p) sumd += denom[p];
	      for (int j = 0, d = dd; j < nsize0; ++j, ++d) {
		double	rr = RefFreq[d] / ff;
		double	sumn = 0;
		for (int p = 0; p < 3; ++p) {
		  sumn += phase[p][d];
		  fprintf(out_fd, "%10.5f\t", lod(phase[p][d], denom[p], rr));
		}
		fprintf(out_fd, "%10.5f\n", lod(sumn, sumd, rr));
	      }
	    }
	    delete[] RefFreq; RefFreq = 0;
	} else {	/* Frequency */
	    int	d = 0;
	    for (int i = 0; i < nsize; ++i) {
	      for (int j = 0; j < nsize; ++j, ++d) {
		double	ff = 0.;
		for (int p = 0; p < 3; ++p) ff += phase[p][d];
		if (ff == 0.) continue;
		int	b = i;
		for (int p = 0; p < 3; ++p, b = (b << 2) & 63)
		    putc(Nucl[b >> 4], out_fd);
		b = j;
		for (int p = 0; p < 3; ++p, b = (b << 2) & 63)
		    putc(Nucl[b >> 4], out_fd);
		putc('\t', out_fd);
		for (int p = 0; p < 3; ++p) 
		    fprintf(out_fd, "%10.2f\t", phase[p][d]);
		fprintf(out_fd, "%10.2f\n", ff);
	      }
	    }
	}
}

void Composit::csvcodonfrq(Seq* sd)
{
	if (sd) {
	    if (!phase[0]) {
		nsize = 64;
		nsize2 = 64 * 64;
		phase[0] = new double[3 * nsize2];
		for (int i = 1; i < 3; ++i) {
		    phase[i] = phase[i-1] + nsize2;
		}
		vclear(phase[0], 3 * nsize2);
	    } else if (each) {
		vclear(phase[0], 3 * nsize2);
	    }
	    CHAR*	ss = sd->at(sd->left);
	    CHAR*	tt = sd->at(sd->right - 2);

	    FTYPE*	wt = 0;
#if USE_WEIGHT
	    wt = sd->weight;
#endif
	    for (int n = 0; ss < tt; ++n) {
		int p = n % 3;
		int	d = 0;
		FTYPE*	w = wt;
		for (int i = 0; i < sd->many; ++i) {
		    int	b = codon_id(ss++, sd->many);
		    if (b != ERROR) {
			if (i == 0)	d = nsize * b;
			else if (wt)	phase[p][d + b] += *w++;
			else		phase[p][d + b] += 1.;
		    } else if (i == 0) {
			ss += sd->many - 1;
			break;
		    } else if (wt) w++;
		}
	    }
	    if (!each) return;
	}

	for (int i = 0, d = 0; i < nsize; ++i) {
	    for (int j = 0; j < nsize; ++j, ++d) {
		double	ff = 0.;
		for (int p = 0; p < 3; ++p) ff += phase[p][d];
		if (ff == 0.) continue;
		int	b = i;
		for (int p = 0; p < 3; ++p, b = (b << 2) & 63)
		    putc(Nucl[b >> 4], out_fd);
		b = j;
		for (int p = 0; p < 3; ++p, b = (b << 2) & 63)
		    putc(Nucl[b >> 4], out_fd);
		putc('\t', out_fd);
		for (int p = 0; p < 3; ++p) 
		    fprintf(out_fd, "%10.2f\t", phase[p][d]);
		fprintf(out_fd, "%10.2f\n", ff);
	    }
	}
}

static double* DiNucFrq(double* frq, Seq* sd)
{
	int	nsize = sd->code->max_code;
	int	nsize2 = nsize * nsize;
	CHAR*	ss = sd->at(sd->left);
	CHAR*	tt = sd->at(sd->right - 1);
	CHAR*	uu = ss + sd->many; 

	FTYPE*	wt = 0;
#if USE_WEIGHT
	wt = sd->weight;
#endif
	if (!frq) frq = new double[nsize2];
	vclear(frq, nsize2);
	for (int n = 0; ss < tt; ++n) {
	    FTYPE*	w = wt;
	    for (int i = 0; i < sd->many; ++i, ++ss, ++uu) {
		int	d = nsize * *ss + *uu;
		if (wt) frq[d] += *w++;
		else	frq[d] += 1.;
	    }
	}
	return (frq);
}

void Composit::dinucfrq(Seq* sd)
{
/*	Thermodynamic parameters for RNA/DNA hybrid
	nearest neighbor stacking interactions
	by Sugimoto et al. Biochemistry 1995
*/

static	const	double	DeltaH[4][4] = {
		{-7.8, -5.9, -9.1, -8.3},
		{-9.0, -9.3, -16.3, -7.0},
		{-5.5, -8.0, -12.8, -7.8},
		{-7.8, -8.6, -10.4, -11.5}
	};
static	const	double	DeltaS[4][4] = {
		{-21.9, -12.3, -23.5, -23.9},
		{-26.1, -23.2, -47.1, -19.7},
		{-13.5, -17.1, -31.9, -21.6},
		{-23.2, -22.9, -28.4, -36.4}
	};
static	double	DHinit = 1.9;
static	double	DSinit = -3.9;
	if (!sd) return;

	double	DH = DHinit;
	double	DS = DSinit;
	double*	frq = DiNucFrq(0, sd);
	double	nml = sd->many * (sd->right - sd->left - 1);
	sqcode = sd->code;
	double*	wk = frq + sqcode->base_code * nsize;
	for (int i = sqcode->base_code; i < nsize; ++i) {
	    int	m = sqcode->redctab[i];
	    if (i == sqcode->amb_code) {
		wk += nsize;
		continue;
	    }
	    int	j = sqcode->base_code;
	    wk += j;
	    for ( ; j < nsize; ++j, ++wk) {
		int	n = sqcode->redctab[j];
		if (j != sqcode->amb_code && *wk > 0.)
		    fprintf(out_fd, "%c%c\t%10.2f\t%10.5f\n", sqcode->decode[i], 
			sqcode->decode[j], *wk, *wk / nml);
		if (m <= 3 && n <= 3) {
		    DH += DeltaH[m][n] * *wk;
		    DS += DeltaS[m][n] * *wk;
		}		    
	    }
	}
	fprintf(out_fd, "\nStacks = %.2lf DH = %6.1lf DS = %6.1lf DG37 = %6.1lf\n",
		nml, DH, DS, DH - .31015 * DS);
	delete[] frq;
}

Composit::Composit(Seq* sd)
	: codontab(0), TotalNseq(0), TotalNres(0), nsize2(0), 
	  molc(sd->inex.molc), sqcode(sd->code), each(algmode.mlt)
{
	gccmp[0] = 0;
	nsize = molc == TRON? sqcode->max_code - sqcode->base_code: sqcode->max_code;
	switch (algmode.alg) {
	    case 0: case 1: case 3: break;
	    case 2: nsize2 = nsize * nsize; break;
	    case 4: case 5: if (molc != TRON) nsize = 64; break;
	    case 6: case 7: case 8: case 9:
		if (molc != TRON) nsize = 64; 
		nsize2 = nsize * nsize;
		break;
	    default:
		fatal("-O# :: 1: mono; 2: di; 3-5: codon; 6: dicodon\n");
	}
	if (!nsize2) nsize2 = nsize;
	phase[0] = new double[3 * nsize2];
	for (int i = 1; i < 3; ++i) phase[i] = phase[i - 1] + nsize2;
	vclear(phase[0], 3 * nsize2);
}

void Composit::calculate(Seq* sd)
{
	if (sd) {
	    int	ln = sd->right - sd->left;
	    TotalNseq++;
	    TotalNres += ln;
	    molc = sd->inex.molc;
	    if (each) {
		if (algmode.nsa && algmode.alg) {
		    sd->fphseq(out_fd);
		    fprintf(out_fd, " %6d nt\'s, %3d codons\n", ln, ln/3);
		} else if (algmode.alg) {
		    fprintf(out_fd, "%s", sd->sqname());
		} else {
		    fprintf(out_fd, "%-31s\t%7d\t%7d\t%d\t", sd->sqname(), 
			sd->SiteNo(sd->left), sd->SiteNo(sd->right - 1), ln);
		}
	    }
	} else if (each) return;
	switch (algmode.alg) {
	    case 0: case 1:
		if (molc == TRON) monotron(sd);
		else	monocomp(sd);
		break;
	    case 2:
		if (molc == TRON) ditronfrq(sd);
		else	dinucfrq(sd);
		break;
	    case 3:
		codonuse(sd); break;
	    case 4: case 5:
		if (molc == TRON) monotron(sd);
		else	monocodon(sd);
		break;
	    case 6: case 7:
		if (molc == TRON) ditronfrq(sd);
		else	dicodonfrq(sd);
		break;
	    case 8: case 9:
		if (molc == TRON) csvtronfrq(sd);
		else	csvcodonfrq(sd);
		break;
	    default:
		fatal("-O# :: 1: mono; 2: di; 3-5: codon; 6: dicodon\n");
	}
}

static void compstn(Seq* sd)
{
	out_fd = qout("");

	if (!out_fd) return;
	Composit	cmpst(sd);
	cmpst.calculate(sd);
	qclose();
}

static void printnbr(Seq* sd)
{
	int*	nbr = new int[sd->many];
	CHAR*	ps = sd->at(0);
	CHAR*	ts = sd->at(sd->left);
	FILE*	fd = qout("");

	if (!fd) return;
	for (int i = 0; i < sd->many; ++i)
	    nbr[i] = sd->nbr[i];
	while (ps < ts)
	    for (int i = 0; i < sd->many; ++i, ++ps)
		if (IsntGap(*ps)) ++nbr[i];
	ts = sd->at(sd->right);
	while (ps < ts) {
	    for (int i = 0; i < sd->many; ++i, ++ps) {
		if (IsntGap(*ps)) ++nbr[i];
		fprintf(fd, "%3d\t", nbr[i]);
	    }
	    putc('\n', fd);
	}
	qclose();
	delete[] nbr;
}

static char* pgetstr(char* pat, char* prmpt)
{
	*pat = '\0';
	progets(pat, prmpt);
	return (*pat? pat: 0);
}

static void allezm(FILE* fd, Seq* sd, char* pat)
{
	int	loc[MAXLOC];
	RESITE* res = firstresite();
	char	rsq[RERSIZE+1];
	int	minsite = 1;
	int	maxsite = POS_INT;

	pat = cdr(pat);
	if (*pat) {
	    maxsite = atoi(pat);
	    pat = cdr(pat);
	    if (*pat) {
		minsite = maxsite;
		maxsite = atoi(pat);
	    }
	}
	if (maxsite == 0) minsite = 0;
	*rsq = '\0';
	for ( ; res; res = nextresite()) {
	    if (!strcmp(res->chr, rsq)) continue;
	    int	num = respos(loc, MAXLOC, sd, res);
	    if (minsite <= num && num <= maxsite) {
		fprintf(fd, "%-10s %-10s %2d   %d\n", 
		    res->name, res->chr, res->cut, num);
		putloc(fd, sd, loc, num);
		strcpy(rsq, res->chr);
	    }
	}
}

static void find(Seq* ns)
{
	FILE*	fd = qout("");

	if (!fd) return;
	while (pgetstr(pattern, spat)) {
	    findpat(fd, ns, pattern);
	}
	qclose();
}

static void resezm(Seq* ns)
{
	FILE*	fd = qout("");

	if (!fd) return;
	while (pgetstr(pattern, sezm)) {
	    if (!strncmp(pattern, "all", 3))
		allezm(fd, ns, pattern);
	    else
		resites(fd, ns, pattern);
	}
	qclose();
}

static void nameseq(FILE* fd, Seq* ns)
{
	fprintf(fd, "%-16s [ %3d ] %5d %5d\n", ns->sqname(),
	    ns->many? ns->many: 1, 
	    ns->SiteNo(ns->left), ns->SiteNo(ns->right - 1));
}

static Seq* forgeseq(FILE* fd, Seq* sd)
{
static	double	acgt[4] = {0.25, 0.25, 0.25, 0.25};

	promptin("Length [%d] : ", &sd->len);
	promptin("A[%6.3lf] C[%6.3lf] G[%6.3lf] T[%6.3lf] : ", 
	    acgt, acgt + 1, acgt + 2, acgt + 3);
	sd->randseq(acgt);
	sd->typeseq(fd);
	return (sd);
}

static void shufseq(Seq* sd)
{
	FILE*	fd = qout("");

	if (!fd) return;
	sd->rndseq();
	sd->typeseq(fd);
	qclose();
}

void utn_setup(CalcServer<Seq>* svr)
{
	if (svr->jobcode == 'c' || svr->jobcode == 'C')
	    svr->prm = (void*) (new Composit(svr->in_face[0]));
}

void utn_cleanup(CalcServer<Seq>* svr)
{
	if (svr->prm) {
	    Composit*	cmpjob = (Composit*) svr->prm;
	    cmpjob->calculate(0);
	    delete cmpjob;
	    svr->prm = 0;
	}
}

int utn_main(CalcServer<Seq>* svr, Seq* seqs[], ThQueue<Seq>* q)
{
	Composit*	cmpjob = (Composit*) svr->prm;
	ORF*	orfs = 0;
	RANGE	rng = {0, 0};
	Seq*&	sd = seqs[0];
	if (longestorf && !(svr->jobcode == 'o' && svr->jobcode == 'O')) {
	    orfs = sd->getorf();
	    if (orfs) {
		sd->saverange(&rng);
		sd->left  = orfs->pos;
	        sd->right = orfs->pos + orfs->len;
		if (svr->jobcode != 'c' && svr->jobcode != 'C') sd->right += 3;
	    } else	return (ERROR);
	}
	if (svr->jobcode != 'a' && svr->jobcode != 'n') polyA.rmpolyA(sd, q_mns);
	switch (svr->jobcode) {
	    case 'a': poly_a(sd); break;
	    case 'B': fouteij(stdout, sd); break;
	    case 'c':	
	    case 'C': cmpjob->calculate(sd); break;
	    case 'f':
	    case 'F': findpat(stdout, sd, pattern); break;
	    case 'g': forgeseq(stdout, sd); break;
	    case 'l':
	    case 'L': Proutseq(stdout, sd); break;
	    case 'k':
	    case 'K': printorf(stdout, sd); break;
	    case 'n':
	    case 'N': nameseq(stdout, sd); break;
	    case 'o':
	    case 'O': putorf(stdout, sd); break;
	    case 'Q': mapsite(sd, 0); break;
	    case 't':
	    case 'T': transorf(stdout, sd); break;
	    case 'x': muteseq(stdout, seqs); break;
	    case 'Z': mapsite(sd, 1); break;
	    default:  usage();
	}
	if (orfs) {
	    sd->restrange(&rng);
	    delete[] orfs;
	}
	return (OK);
}

template <class seq_t>
void AlnServer<seq_t>::menu_job()
{
	char	cmd[CONLINE];
	algmode.mlt = 2;
	int	c = *restsq(cmd, 0, 1);
	setprompt(1, 0);
	Seq*&	sd = this->in_face[0];
	if (c != 'a') polyA.rmpolyA(sd, q_mns);
	do {
	    switch (c) {
		case '1': inputseq(this->in_face, cmd); break;
		case 'a': aa_code(); break;
		case 'A': poly_a(sd); break;
		case 'B': fouteij(0, sd); break;
		case 'c': algmode.alg &= ~2; compstn(sd); break;
		case 'C': algmode.alg |= 2; compstn(sd); break;
		case 'e': extractcds(sd); break;
		case 'f': find(sd); break;
		case 'g': forgeseq(stdout, sd); break;
		case 'k': 
		case 'K': codseq(sd); break;
		case 'l': 
		case 'L': listseq(sd, cmd+1); break;
		case 'n': printnbr(sd); break;
		case 'o':
		case 'O': openfrm(sd); break;
		case 'Q': mapsite(sd, 0); break;
		case 'r':
		case 'R': resezm(sd); break;
		case 'S': shufseq(sd); break;
		case 't':
		case 'T': transl(sd); break;
		case 'p': setparam(1); break;
		case 'P': setparam(2); break;
		case 'x': mute(this->in_face); 
			  listseq(this->in_face[1], cmd+1); break;
		case 'Z': mapsite(sd, 1); break;
		default:  break;
	    }
	    prompt("\n#(+)/Aa/Comp/Ext/Find/List/Map/Nbr/Orf/");
	    c = *progets(cmd, "Quit/Rezm/Shufl/Trans/Weight/Xout : ");
	} while (c != 'q');
}

int main(int argc, const char** argv)
{
	setdefmolc(DNA);
	(void) setSimmtxes(DxD);
	alprm2.spb = 1;
	initcodon(0);
	setprmode(Row_Last, 'L', SILENT);
	AlnServer<Seq>	svr(argc, argv, IM_SNGL, IM_SNGL, (void*) 0, &utn_main,
	    &utn_setup, &utn_cleanup);
	if (svr.autocomp() == 1) svr.menu_job();
	EraDbsDt();
	eraStrPhrases();
	resetSimmtxes(true);
	return (0);
}
