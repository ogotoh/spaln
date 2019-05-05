/*****************************************************************************
*
*	spaln [-Qn] [-On] [-Ls] [-Rs|-Ws] [-Ts] [-ds|-cs|-as] <genome|aa_seq> cDNA/Aa
*
*	Mapping and alignment of protein/cDNA sequences onto
*	genomic sequence or aa database seq
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-41-6 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "aln.h"
#include "divseq.h"
#include "vmf.h"
#include "wln.h"
#include "blksrc.h"
#include <math.h>

#ifndef	M_THREAD
#define	M_THREAD	1
#endif

#if M_THREAD
#include <pthread.h>
#include <unistd.h>
#include <sched.h>
#endif

#define	USE_FULL_RANGE	0
#define	TESTRAN		0
#define	QDEBUG	0

static	const	int	MAX_ARGS = 128;
static	const	int	MAX_ARGTEXT = 4096;
static	const	int	DefPam = 150;
static	const	int	WlpPam = 50;
static	const	int	DefOrfLen = 75;
static	const	int	spacer = 32;
static	const	int	def_max_extend_gene_rng = 3;
static	const	long	MaxWordNoSpace = 1 << 30;	// 1GB

enum QvsD {UNDF, AvsA, AvsG, CvsC, CvsG, FvsG, GvsA, GvsC, GvsG};
const char* const QvsDstr[] = {"UNDF", "AvsA", "AvsG", "CvsC", "CvsG", "FvsG", "GvsA", "GvsC", "GvsG"};

static	QvsD	QRYvsDB = CvsG;
static	InputMode	input_form = IM_SNGL;

class HalfGene {
public:
	int	left;
	int	right;
	int	segno;
	char	*sname;
	HalfGene(int l, int r, int n) :
	    left(l), right(r), segno(n)
		{sname = 0;}
	~HalfGene() {delete[] sname;}
};

class PolyA {
	int	thr;
public:
	PolyA()	{thr = (int) alprm.thr;}	// remove polyA
	void	setthr(const char* s) {thr = *s? atoi(s): 0;}
	SHORT	rmpolyA(Seq* sd);
};

#if M_THREAD
class ThQueue {
	int	rp, wp;
	int	remain;
	int	done;
	pthread_cond_t	not_full;
	pthread_cond_t	not_empty;
public:
	Seq**	sinp;
	Seq**	sque;
	pthread_mutex_t	mutex;
	Mfile*	mfd;
	ThQueue(Seq** sqs);
	~ThQueue() {delete[] sque;}
	void	enqueue(Seq** fsd, int n = 1);
	void	dequeue(Seq** fsd, int n = 1);
	void	putqueue();
};

struct thread_arg_t {
	ThQueue*	q;
	int	cpuid;
	Seq**	seqs;
	SeqServer*	svr;
	void*	pwd;
};

struct mast_arg_t {
	ThQueue*	q;
	SeqServer*	svr;
};

struct mist_arg_t {
	ThQueue*	q;
	int	nhf;
	HalfGene*	hfg;
};

static	void*	master_func(void* targ);
static	void*	mistress_func(void* arg);
static	void*	worker_func(void* targ);
static	void 	MasterWorker(Seq** sqs, SeqServer* svr, void* prm);
#else
class	ThQueue;
#endif

static	void	usage(const char* messg);
static	int	getoption(int argc, const char** argv);
static	void	readargs();
static	RANGE*	skl2exrng(SKL* skl);
static	int	spalign2(Seq* sqs[], PwdB* pwd, Gsinfo* GsI);
static	void	pairends(Seq* sqs[], PwdB* pwd, Gsinfo* GsI, int len);
static	int	alnoutput(Seq* sqs[], Gsinfo* gsi);
static	SrchBlk*	getblkinf(Seq* sqs[], const char* dbs, MakeBlk* mb);
static	void	seg_job(Seq** sqs, SeqServer* svr, SrchBlk* sbk);
static	void	all_in_func(Seq** sqs, SeqServer* svr, void* prm);
static	void	setdefparam();
static	PwdB*	SetUpPwd(Seq* sqs[]);
static	int	blkaln(Seq* sqs[], SrchBlk* bks, RANGE* rng, ThQueue* q);
static	int	match_2(Seq* sqs[], PwdB* pwd, ThQueue* q = 0);
static	int	quick4(Seq* sqs[], SrchBlk* bks, ThQueue* q = 0);
static	void	spaln_job(Seq* sqs[], void* prm, ThQueue* q = 0);
static	void	genomicseq(Seq** sqs, PwdB* pwd, bool reverse);
static	int	matepair(Seq** sqs, Gsinfo* gsi, VTYPE thr);

static	const	char*	ReadBlock = 0;
static	const	char*	catalog = 0;
static	const	char*	aadbs = 0;
static	const	char*	genomedb = 0;
static	const	char*	cdnadb = 0;
static	bool	WriteMolc = false;
static	int	QryMolc = UNKNOWN;
static	int	TgtMolc = UNKNOWN;
static	int	MinQueryLen = 0;
static	int	g_segment = 2 * MEGA;
static	int	q_mns = 3;
static	PolyA	polyA;
static	int	no_seqs = 0;
static	bool	pairedends = false;
static	bool	gsquery = QRYvsDB == GvsA || QRYvsDB == GvsC;
static	const	char*	version = "2.3.3c";
static	const	int	date = 190502;

static void usage(const char* messg)
{
	if (messg && *messg) fputs(messg, stderr);
	fprintf(stderr, "\n*** SPALN version %s <%d> ***\n\n", version, date);
	fputs("Usage:\n", stderr);
	fputs("spaln -WGenome.bkn -KD [W_Options] Genome.mfa\t(to write block inf.)\n", stderr);
	fputs("spaln -WGenome.bkp -KP [W_Options] Genome.mfa\t(to write block inf.)\n", stderr);
	fputs("spaln -WAAdb.bka -KA [W_Options] AAdb.faa\t(to write aa db inf.)\n", stderr);
	fputs("spaln [R_options] genomic_segment cDNA.fa\t(to align)\n", stderr);
	fputs("spaln [R_options] genomic_segment protein.fa\t(to align)\n", stderr);
	fputs("spaln [R_options] -dGenome cDNA.fa\t(to map & align)\n", stderr);
	fputs("spaln [R_options] -dGenome protein.fa\t(to map & align)\n", stderr);
	fputs("spaln [R_options] -aAAdb genomic_segment.fa\t(to search aa database)\n", stderr);
	fputs("spaln [R_options] -aAAdb protein.fa\t(to search aa database)\n", stderr);
	fputs("\nin the following, # = integer or real number; $ = string; default in ()\n\n", stderr);
	fputs("W_Options:\n", stderr);
	fputs("\t-Xk#\tWord size (11)\n", stderr);
	fputs("\t-Xb#\tBlock size (4096)\n", stderr);
	fputs("\t-XG#\tMaximum expected gene size (262144)\n", stderr);
	fputs("\t-Xa#\tAbundance factor (10)\n", stderr);
	fputs("\t-g\tgzipped output\n\n", stderr);
	fputs("R_Options (representatives):\n", stderr);
	fputs("\t-H#\tMinimum score for report (35)\n", stderr);
	fputs("\t-L or -LS or -L#\tsemi-global or local alignment (-L)\n", stderr);
	fputs("\t-M#\tNumber of outputs per query (1) (max=4 for genome vs cDNA|protein)\n", stderr);
	fputs("\t-O# (GvsA|C)\t0:Gff3_gene; 1:alignment; 2:Gff3_match; 3:Bed; 4:exon-inf;\n", stderr);
	fputs("\t\t\t5:intron-inf; 6:cDNA; 7:translated; 8:block-only;\n", stderr); 
	fputs("\t\t\t10:SAM; 12:binary (4)\n", stderr);
	fputs("\t-O# (AvsA)\t0:statistics; 1:alignment; 2:Sugar; 3:Psl; 4:XYL;\n", stderr);
	fputs("\t\t\t5:srat+XYL; 8:Cigar; 9:Vulgar; 10:SAM (4)\n", stderr); 
	fputs("\t-Q#\t0:DP; 1-3:HSP-Search; 4-7; Block-Search (3)\n", stderr);
	fputs("\t-R$\tRead block information file *.bkn, *.bkp or *.bka\n", stderr);
	fputs("\t-S#\tOrientation. 0:annotation; 1:forward; 2:reverse; 3:both (0)\n", stderr);
	fputs("\t-T$\tSubdirectory where species-specific parameters reside\n", stderr);
	fputs("\t-a$\tSpecify AAdb. Must run `makeidx.pl -ia' breforehand\n", stderr);
	fputs("\t-A$\tSame as -a but db sequences are stored in memory\n", stderr);
	fputs("\t-d$\tSpecify genome. Must run `makeidx.pl -i[n|p]' breforehand\n", stderr);
	fputs("\t-D$\tSame as -d but db sequences are stored in memory\n", stderr);
	fputs("\t-g\tgzipped output in combination with -O12\n", stderr);
	fputs("\t-iC\tpaired inputs. C=a:alternative; C=p:parallel\n", stderr);
	fputs("\t-l#\tNumber of characters per line in alignment (60)\n", stderr);
	fputs("\t-o$\tFile where results are written (stdout)\n", stderr);
	fputs("\t-pa\tDon\'t remove poly A\n", stderr);
	fputs("\t-pw\tReport results even if the score is below the threshold\n", stderr);
	fputs("\t-pq\tQuiet mode\n", stderr);
	fputs("\t-r$\tReport information about block data file\n", stderr);
	fputs("\t-u#\tGap-extension penalty (3)\n", stderr);
	fputs("\t-v#\tGap-open penalty (8)\n", stderr);
	fputs("\t-w#\tBand width for DP matrix scan (100)\n", stderr);
	fputs("\t-t[#]\tMutli-thread operation with # CPUs\n", stderr);
	fputs("\t-ya#\tStringency of splice site. 0->3:strong->weak\n", stderr);
	fputs("\t-yl3\tDdouble affine gap penalty\n", stderr);
	fputs("\t-ym#\tNucleotide match score (2)\n", stderr);
	fputs("\t-yn#\tNucleotide mismatch score (-6)\n", stderr);
	fputs("\t-yy#\tWeight for splice site signal (8)\n", stderr);
	fputs("\t-yz#\tWeight for coding potential (2)\n", stderr);
	fputs("\t-yB#\tWeight for branch point signal (0)\n", stderr);
	fputs("\t-yI$\tIntron length distribution\n", stderr);
	fputs("\t-yL#\tMinimum expected length of intron (30)\n", stderr);
	fputs("\t-yS[#]\tUse species-specific parameter set\n", stderr);
	fputs("\t-yX\tUse parameter set for cross-species comparison\n", stderr);
	fputs("\t-yZ#\tWeight for intron potential (0)\n", stderr);
	fputs("\t-XG#\tReset maximum expected gene size, suffix k or M is effective\n", stderr);
	fputs("\nExamples:\n", stderr);
	fputs("\tspaln -Wdictdisc_g.bkn -KD -Xk10 -Xb8192 dictdisc_g.gf\n", stderr);
	fputs("\tspaln -WSwiss.bka -KA -Xk5 Swiss.faa\n", stderr);
	fputs("\tspaln -O1 -LS \'chr1.fa 10001 40000\' cdna.nfa\n", stderr);
	fputs("\tspaln -Q7 -yX -t10 -TTetrapod -XG1M -ommu.exon -dmus_musc_g hspcdna.nfa\n", stderr);
	fputs("\tspaln -Q7 -O5 -t10 -Tdictdics -ddictdisc_g \'dictdisc.faa (101 200)\' > ddi.intron\n", stderr);
	fputs("\tspaln -Q7 -O0 -t10 -Tdictdics -aSwiss \'chr1.nfa 200001 210000\' > Chr1_200-210K.gff\n", stderr);
	fputs("\tspaln -Q4 -O0 -t10 -M10 -aSwiss dictdisc.faa > dictdisc.alignment_score\n", stderr);
	exit(1);
}

static int getoption(int argc, const char** argv)
{
const	char**	argbs = argv;
	DbsDt**	dbs = dbs_dt;

	while (--argc > 0 && **++argv == OPTCHAR) {
	    const char*	opt = argv[0] + 1;
	    int	c = *opt;
	    if (!c) break;
	    const char*	val = argv[0] + 2;
	    int	k = 0;
	    int	rv = 0;
	    switch (c) {
		case '?': case 'h': usage(0);
		case 'a': case 'A':
		    algmode.dim = c == 'A';
		    if ((val = getarg(argc, argv)))
			{delete *dbs; *dbs = new DbsDt(aadbs = val);}
		    break;
		case 'd': case 'D': 
		    algmode.dim = c == 'D';
		    if ((val = getarg(argc, argv)))
			{delete *dbs; *dbs = new DbsDt(genomedb = val);}
		    break;
		case 'c':
		    if ((val = getarg(argc, argv)))
			{delete *dbs; *dbs = new DbsDt(cdnadb = val);}
		    break;
		case 'C':
		    if ((val = getarg(argc, argv, true)))
			initcodon(atoi(val));
		    break;
		case 'f':
		    if ((val = getarg(argc, argv))) catalog = val;
		    break;
		case 'g': OutPrm.gzipped = 1; break;
		case 'F':
		    if (('a' <= *val && *val <= 'z') ||
			('0' <= *val && *val <= '9'))
			setprmode(*val, SILENT, SILENT); else
		    if (*val == 'C' || ('K' <= *val && *val <= 'O'))
			setprmode(SILENT, *val, SILENT); else
		    if ('A' <= *val && *val <= 'Z')
			setprmode(SILENT, SILENT, *val);
		    break;
		case 'G':
		    if ((val = getarg(argc, argv)))
			g_segment = ktoi(val);
		    break;
		case 'H':
		    if ((val = getarg(argc, argv, true)))
			setthr(atof(val));
		    else
			setthr((double) (INT_MIN + 3));
		    break;
		case 'i': 
		    switch (*val) {
			case 'a': case 'A': case '2': 
			    input_form = IM_ALTR; break;
			case 'p': case 'P': case '3': 
			    input_form = IM_PARA; break;
			case ':': catalog = val + 1; break;
			default:  input_form = IM_SNGL; break;
		    }
		    break;
		case 'I':
		    if ((val = getarg(argc, argv, true)))
			setorf(SILENT, atoi(val));
		    break;
		case 'J':
		    if ((val = getarg(argc, argv, true)))
			setorf(atoi(val), SILENT);
		    break;
		case 'K':
		    switch (*val) {
			case 'A': case 'a':
			    QryMolc = TgtMolc = PROTEIN;
			    algmode.lsg = 0;
			    break;
			case 'C': case 'c': 
			    QryMolc = TgtMolc = DNA; break;
			    algmode.lsg = 0;
			    break;
			case 'T': case 't':
			    QryMolc = TgtMolc = TRON; break;
			case 'P': case 'p':
			    QryMolc = PROTEIN; 
			    TgtMolc = TRON; break;
			case 'D': case 'd':
			default:
			    QryMolc = TgtMolc = DNA; break;
		    }
		    break;
		case 'L': 
		    if (isdigit(*val)) algmode.lcl = atoi(val);
		    else {
			switch (*val) {
			    case 'l': 
			    case 'L': algmode.lcl = 5; break;
			    case 'r':
			    case 'R': algmode.lcl = 10; break;
			    case 'd':
			    case 'D': algmode.lcl = 3; break;
			    case 'u':
			    case 'U': algmode.lcl = 12; break;
			    case 's':	// Smith-Waterman
			    case 'S': algmode.lcl = 16; break;
			    default:  algmode.lcl = 15; break;
			}
		    }
		    break;
		case 'l':
		    if ((val = getarg(argc, argv, true)))
			setlpw(atoi(val));
		    else if (*val == '-') setlpw(0);
		    break;
		case 'M': 
		    if ((val = getarg(argc, argv, true))) {
			OutPrm.MaxOut = atoi(val);
			algmode.mlt = OutPrm.MaxOut == 1? 1: 2;
		    } else {
			OutPrm.MaxOut = MAX_PARALOG;
			algmode.mlt = 2;
		    }
		    break;
		case 'm':
		    if ((val = getarg(argc, argv, false))) {
			if ((opt = strchr(val, ':')))
			    k = atoi(opt + 1) - 1;
			if (0 <= k && k < max_simmtxes) mdm_file[k] = val;
		    }
		    break;
		case 'n': algmode.nsa = 0; break;
		case 'N':
		    if ((val = getarg(argc, argv, true)))
			 MinQueryLen = atoi(val);
		    break;
		case 'O':
		    if ((val = getarg(argc, argv, true)))
			algmode.nsa = atoi(val);
		    break;
		case 'o':
		    if ((val = getarg(argc, argv)))
			OutPrm.out_file = val;
		    break;
		case 'p':
		    switch (*val) {
			case 'a': polyA.setthr(val + 1); break;
			case 'd': OutPrm.descrp = 1; break;
			case 'e': OutPrm.trimend = !OutPrm.trimend; break;
			case 'f': OutPrm.deflbl = 1; break;
			case 'i': OutPrm.ColorEij = 1; break;
			case 'j': OutPrm.spjinf = 1 - OutPrm.spjinf; break;
			case 'q': setprompt(0, 0); break;
			case 'Q': setprompt(0, 1);  break;
			case 't': OutPrm.deflbl = 2; break;
			case 'w': algmode.thr = 0; break;
			case 'x': OutPrm.supself = (val[1] == '2')? 2: 1; break;
			case 'J': case 'm': case 'u': case 'v': setexprm_z(val); break;
			default: break;
		    }
		    break;
		case 'Q':
		    if ((val = getarg(argc, argv, true))) {
			int	q = atoi(val);
			algmode.qck = q & 3;
			algmode.blk = q >> 2;
		    } else	algmode.qck = 1;
		    break;
#if M_THREAD
		case 'q':
		    if ((val = getarg(argc, argv, true)))
			max_queue_num = atoi(val);
		    else max_queue_num = 0;
		    break;
#endif
		case 'R':
		    if ((val = getarg(argc, argv)))
			ReadBlock = val;
		    break;
		case 'r':
		    if ((val = getarg(argc, argv)))
			ReportBlkInfo(val); 	// exit
		case 's':
		    if ((val = getarg(argc, argv)))
			setdfn(val);
		    break;
		case 'S':
		    if (*val == 'F') algmode.mns = 1; else
		    if (*val == 'R') algmode.mns = 2; else
		    if (*val == 'B') algmode.mns = 3; else
		    if (*val == '+' || *val == 'f') q_mns = 1; else
		    if (*val == '-' || *val == 'r') q_mns = 2; else
		    if (isdigit(*val)) q_mns = atoi(val); else
		    if (argc > 1 && isdigit(argv[1][0]))
			{q_mns = atoi(*++argv); --argc;}
		    else algmode.slv = 1;	// salvage mode | read annotation
		    break;
#if M_THREAD
		case 't':
		    if ((val = getarg(argc, argv, true)))
			thread_num = atoi(val);
		    else thread_num = -1;
		    break;
#endif
		case 'T':
		    if ((val = getarg(argc, argv)))
			ftable.setpath(val, GNM2TAB);
		    readargs();
		    break;
		case 'u': case 'v': case 'w': readalprm(opt); break;
		case 'U': algmode.lsg = 0; break;
		case 'V':
		    if ((val = getarg(argc, argv)))
			setVmfSpace(ktol(val));
		    break;
		case 'W': setQ4prm(opt);
		    if (!*val) val = argv[1];
		    rv = setQ4prm(opt, val);
		    if (rv && argc > 1) {++argv; --argc;}
		    WriteMolc = true; break;
		case 'y': readalprm(val); break;
		case 'x': setexprm_x(val); break;
		case 'X': 
		    if (!*++val) val = argv[1];
		    rv = setQ4prm(opt + 1, val);
		    if (rv) {++argv; --argc;}
		    break;
		case 'z': setexprm_z(val); break;
		default: break;
	    }
	    if (*dbs && (dbs + 1 < dbs_dt + MAX_DBS)) ++dbs;
	}
	if (OutPrm.supself) ++OutPrm.MaxOut;
	return (argv - argbs);
}

static void readargs()
{
	FILE*	fd = ftable.fopen(AlnParam, "r");
	if (!fd) return;
	char	argstr[MAX_ARGTEXT];
	char*	pa = argstr;
const	char*	argv[MAX_ARGS];
	int	argc = 0;
	char	str[MAXL];
	int	quote = 0;
	int	dquote = 0;
	argv[argc++] = 0;
	while (char* ps = fgets(str, MAXL, fd)) {
	    for (int spc = 1; *ps; ++ps) {
		if (*ps == '#') break;
		if (*ps == '\'') quote = 1 - quote;
		if (*ps == '\"') dquote = 1 - dquote;
		if (isspace(*ps) && !quote && !dquote) {
		    if (!spc) {
			*pa++ = '\0';
			spc = 1;
		    }
		} else {
		    if (argc >= MAX_ARGS) goto readend;
		    if (spc) {
			argv[argc++] = pa;
			spc = 0;
		    }
		    *pa++ = *ps;
		}
	    }
	}
readend:
	fclose(fd);
	getoption(argc, argv);
}

// Convert alignment coordinates to gene organization
static	RANGE*	skl2exrng(SKL* skl)
{
	int	num = (skl++)->n;
	RANGE*	rng = new RANGE[(num + 1)/ 2 + 2];
	RANGE*	wrng = rng;

	(++wrng)->left = (skl++)->n;	// 5' end
	for ( ; --num; ++skl) {
	    if (skl->m - skl[-1].m == 0 &&	// intron
		skl->n - skl[-1].n > IntronPrm.llmt) {
		    (wrng++)->right = skl[-1].n;
		    wrng->left = skl->n;
	    }
	}
	(wrng++)->right = skl[-1].n;	// 3' end
	*wrng = endrng;
	rng->left = ++wrng - rng;
	return (rng);
}

static int spalign2(Seq* sqs[], PwdB* pwd, Gsinfo* GsI)
{
	switch (QRYvsDB) {
	  case AvsG: 
	  case GvsA: GsI->skl = alignH_ng(sqs, pwd, &GsI->scr); break;
	  case GvsC:
	  case CvsG: GsI->skl = alignS_ng(sqs, pwd, &GsI->scr); break;
	  case FvsG: 
	  case AvsA: 
	  case CvsC: GsI->skl = alignB_ng(sqs, pwd, &GsI->scr); break;
	  default:
		prompt("Improper combination of db and query!\n");
		return (ERROR);
	}
	if (GsI->skl) {
	    if (GsI->skl->n == 0) return (ERROR);
	    if (QRYvsDB == AvsA || QRYvsDB == CvsC || 
		QRYvsDB == GvsG || QRYvsDB == FvsG) {
		GsI->scr = skl_rngB_ng(sqs, GsI, pwd);
	    } else if (GsI->skl->m & AlgnTrb) {
		if (QRYvsDB == CvsG || QRYvsDB == GvsC)
			GsI->scr = skl_rngS_ng(sqs, GsI, pwd);
		else	GsI->scr = skl_rngH_ng(sqs, GsI, pwd);
		GsI->eiscr2rng();		// raise score for each intron
		GsI->scr -= pwd->GapPenalty(Ip_equ_k) * (GsI->noeij - 1);
	    } else {
		GsI->CDSrng = skl2exrng(GsI->skl);
	    }
	} else {
//	    prompt("%s %s failed to align!\n", 
//		sqs[0]->sname(), sqs[1]->sname());
	    return (ERROR);
	}
	return (OK);
}

static void pairends(Seq* sqs[], PwdB* pwd, Gsinfo* GsI, int k)
{
	JUXT*&	jxt1 = sqs[1]->jxt;
	JUXT*	jxt2 = 0;
	int&	CdsNo1 = sqs[1]->CdsNo;
	int	CdsNo2 = 0;
	if (jxt1) {
	    JUXT*	src = jxt1;
	    JUXT*	trm = src + CdsNo1;
	    while (src->jx < k) ++src;
	    CdsNo1 = src - jxt1;
	    CdsNo2 = trm - src;
	    if (CdsNo2) {
		jxt2 = new JUXT[CdsNo2 + 1];
		k += spacer;
		for (JUXT* wsk = jxt2; src <= trm; ++wsk) {
		    *wsk = *src++;
		    wsk->jx -= k;
		}
	    }
	    if (CdsNo1) {
		trm = jxt1 + CdsNo1;
		JUXT*	src = trm - 1;
		trm->jx = src->jx + src->jlen;
		trm->jy = src->jy + src->jlen;
		trm->jlen = 0;
	    } else {delete[] jxt1; jxt1 = 0;}
	}
	swapseq(sqs, sqs + no_seqs - 2);
	GsI->skl = alignS_ng(sqs, pwd, &GsI->scr);
	if (GsI->skl) {
	    GsI->scr = skl_rngS_ng(sqs, GsI, pwd);
	    if (GsI->eijnc) GsI->eiscr2rng();
	}
	swapseq(sqs, sqs + no_seqs - 2);
	swapseq(sqs, sqs + no_seqs - 1);
	delete[] jxt1;
	jxt1 = jxt2;
	CdsNo1 = CdsNo2;
	++GsI;
	GsI->skl = alignS_ng(sqs, pwd, &GsI->scr);
	if (GsI->skl) {
	    GsI->scr = skl_rngS_ng(sqs, GsI, pwd);
	    if (GsI->eijnc) GsI->eiscr2rng();
	}
	swapseq(sqs, sqs + no_seqs - 1);
	delete[] jxt1; jxt1 = 0; CdsNo1 = 0;
}

static int alnoutput(Seq* sqs[], Gsinfo* GsI)
{
	int	print2Skip = 0;
	int	swp = sqs[1]->inex.intr;
	Seq*	gene = 0;
	GAPS*	gaps[2] = {0, 0};

	if (swp) {
	    swapseq(sqs, sqs + 1);
	    swapskl(GsI->skl);
	    gene = sqs[0];
	    delete[] gene->exons;
	    gene->exons = 0;
	}
	if (QRYvsDB == AvsG || QRYvsDB == GvsA) { // DNA vs protein 
	    if (algmode.nsa == ALN_FORM) {
		if (gene) {
		    fphseqs(sqs, 2);
		    GBcdsForm(GsI->CDSrng, gene);
		    print2Skip = 1;
		}
		if (getlpw()) {		// print alignment
		    skl2gaps3(gaps, GsI->skl, 2);
		    if (gene) gene->exons = GsI->eiscrunfold(gaps[0]);
		    unfoldgap(gaps[0], 1);
		    unfoldgap(gaps[1], 3);
		    print2(sqs, gaps, (double) GsI->scr, GsI, 1, 1, print2Skip);
		} else
		    repalninf(sqs, GsI);
	    } else if (gene) {
		printgene(sqs, GsI);
	    }
	} else if (QRYvsDB == CvsG || QRYvsDB == GvsC) {	// genome vs cDNA
	    if (algmode.nsa == ALN_FORM || (algmode.nsa == ITN_FORM && !gene)) {
		if (getlpw()) {
		    skl2gaps(gaps, GsI->skl);
		    if (gene) gene->exons = GsI->eiscrunfold(gaps[0]);
		    toimage(gaps, 2);
		}
		if (gene) {
		    fphseqs(sqs, 2);
		    GBcdsForm(GsI->CDSrng, gene);
		    print2Skip = 1;
		}
		if (getlpw())
		    print2(sqs, gaps, (double) GsI->scr, GsI, 1, 1, print2Skip);
		else
		    repalninf(sqs, GsI);
	    } else if (gene) {
		printgene(sqs, GsI);
	    }
	} else {
	    if (algmode.nsa == ALN_FORM) {
		if (getlpw()) {
		    skl2gaps(gaps, GsI->skl);
		    toimage(gaps, 2);
		    print2(sqs, gaps, (double) GsI->scr, GsI, 1, 1, print2Skip);
		} else {
		    repalninf(sqs, GsI);
		}
	    } else {
		repalninf(sqs, GsI);
	    }
	}
	if (gene) {delete[] gene->exons; gene->exons = 0;}
	if (swp)	{
	    swapseq(sqs, sqs + 1);
	    swapskl(GsI->skl);
	}
	delete[] gaps[0];
	delete[] gaps[1];
	return (OK);
}

//	setup parameters. called once in a run

static	PwdB* SetUpPwd(Seq* sqs[])
{
	Seq*	a = sqs[0];
	Seq*	b = sqs[1];
	int	ap = a->isprotein();
	int	bp = b->isprotein();

	if (!ap && bp) {
	    prompt("Has exchanged the input order of %s and %s !\n",
		b->spath, a->spath);
	    gswap(a, b); gswap(ap, bp);
	}
	if (ap && bp)	{QRYvsDB = AvsA; algmode.lsg = 0;}
	else if (ap)	QRYvsDB = AvsG;
	else if (bp)	QRYvsDB = GvsA;
	else if (a->inex.intr && b->inex.intr) QRYvsDB = GvsG;
	else if (b->inex.intr) QRYvsDB = CvsG;
	else if (a->inex.intr) QRYvsDB = GvsC;
	else if ((b->inex.intr = algmode.lsg)) QRYvsDB = CvsG;
	else		QRYvsDB = CvsC;
	if (ap || bp) q_mns &= ~2;

	gsquery = QRYvsDB == GvsA || QRYvsDB == GvsC;
	if ((QRYvsDB == AvsG || QRYvsDB == GvsA) && algmode.lcl & 16)
	    algmode.lcl |= 3;
	if (gsquery) swapseq(sqs, sqs + 1);
	if (QRYvsDB == GvsA || QRYvsDB == AvsG) sqs[1]->inex.intr = algmode.lsg;
	makeWlprms(prePwd(sqs));
	PwdB* pwd = new PwdB(sqs);
	if (gsquery) swapseq(sqs, sqs + 1);
	return (pwd);
}

static int match_2(Seq* sqs[], PwdB* pwd, ThQueue* q)
{
	Seq*&	a = sqs[0];
	Seq*&	b = sqs[1];
	int	dir = a->isprotein() + 2 * b->isprotein();

	if (dir != pwd->DvsP)
	    fatal("Don't mix different combinations of seq. types %s %d %d\n",
		(*a->sname)[0], dir, pwd->DvsP);
	if (gsquery) swapseq(sqs, sqs + 1);
	if (QRYvsDB == GvsA || QRYvsDB == AvsG) b->inex.intr = algmode.lsg;
	dir = ERROR;
	int	ori = a->inex.ori;
	if (pwd->DvsP != 3) {
	    if (algmode.mns == 0 || algmode.mns == 3 || ori == 3) geneorient(sqs, pwd);
	    if (algmode.mns == 1 && ori == 3) ori = b->inex.sens? 2: 1;
	    if (algmode.mns == 2 && ori == 3) ori = b->inex.sens? 1: 2;
	    if (algmode.lsg == 0 && ori == 3) ori = 1;
	    if (dbs_dt[0] && b->inex.intr && extend_gene_rng(sqs, pwd, dbs_dt[0]))
		genomicseq(sqs + 1, pwd, ori == 3);
	}
	if (algmode.lcl & 16) {	// local
	    a->exg_seq(1, 1);
	    b->exg_seq(1, 1);
	} else {		// global
	    a->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	    b->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	}
	Gsinfo	GsI[2];
	if (ori & 1) dir = spalign2(sqs, pwd, GsI);
	if (ori & 2) {
	    a->comrev();
	    antiseq(sqs + 1);
	    if (b->jxt)
		revjxt(b->jxt, b->CdsNo);
	    if (spalign2(sqs, pwd, GsI + 1) == OK)
		dir = (dir == 0)? GsI[1].scr > GsI[0].scr: 1;
	    if (dir == 0) {
		a->comrev();
		antiseq(sqs + 1);
	    }
	}
	if (dir != ERROR && (!algmode.thr || GsI[dir].fstat.val >= pwd->Vthr)) {
	    if (algmode.nsa == BED_FORM)
		GsI[dir].rscr = selfAlnScr(a, pwd->simmtx);
#if M_THREAD
	    if (q) {
		pthread_mutex_lock(&q->mutex);
		alnoutput(sqs, GsI + dir);
		pthread_mutex_unlock(&q->mutex);
	    } else
#endif
		alnoutput(sqs, GsI + dir);
	}
	if (dir == 1) {
	    a->comrev();
	    antiseq(sqs + 1);
	}
	if (gsquery) swapseq(sqs, sqs + 1);
	return (dir);
}

// remove polyA sequences if present
SHORT PolyA::rmpolyA(Seq* sd)
{
	if (sd->isprotein() || sd->inex.intr) return (1);
	if (thr <= 0) return (3);
	int	polya = 0;
	int	polyt = 0;
	int	scr = 0;
	CHAR*	maxa = 0;
	CHAR*	maxt = 0;
static	int	mn[2] = {-15, 3};

	CHAR*	lend = sd->at(0);
	CHAR*	s = sd->at(sd->len);
	CHAR*	t = s;
	if (q_mns != 2) {
	    while (--s >= lend) {
		scr += mn[*s == A];
		if (scr > polya) {
		    polya = scr;
		    if (scr > thr) maxa = s;
		}
		if (scr < polya - thr) break;
	    }
	}
	if (q_mns != 1) {
	    for (scr = 0, s = lend; s < t; ++s) {
		scr += mn[*s == T];
		if (scr > polyt) {
		    polyt = scr;
		    if (scr > thr) maxt = s;
		}
		if (scr < polyt - thr) break;
	    }
	}
	if (maxa && maxt) {
	    if (polya >= polyt) maxt = 0;
	    else		maxa = 0;
	}
	sd->inex.polA = 0;
	if (maxa) {
	    sd->tlen = maxa - lend;
	    if (sd->right > sd->tlen) sd->right = sd->tlen;
	    sd->inex.polA = 1;
	} else if (maxt) {
	    int	polytlen = maxt - lend;
	    if (sd->left < polytlen) sd->left = polytlen;
	    sd->tlen = sd->len - polytlen;
	    sd->inex.polA = 2;
	}
	return (sd->inex.polA? sd->inex.polA: 3);
}

static int matepair(Seq** sqs, Gsinfo* gsi, VTYPE thr)
{
	bool	f_map = gsi[0].skl && (!algmode.thr || gsi[0].fstat.val >= thr);
	bool	s_map = gsi[1].skl && (!algmode.thr || gsi[1].fstat.val >= thr);
	int	m = f_map + 2 * s_map;
	if (algmode.nsa == SAM_FORM) {
	    Samfmt*	f_sfm = gsi[0].samfm;
	    Samfmt*	s_sfm = gsi[1].samfm;
	    f_map = f_map && f_sfm;
	    s_map = s_map && s_sfm;
	    f_sfm->flag |= 0x1; s_sfm->flag |= 0x1;
	    int	f_last = 0, s_last = 0;
	    if (f_map) {
		SKL*	fst = gsi[0].skl + 1;
		SKL*	lst = gsi[0].skl + gsi[0].skl->n;
		f_sfm->pos = fst->n;
		f_last = lst->n;
		if (sqs[0]->inex.sens) f_sfm->flag |= 0x10;
		if (sqs[1]->inex.sens) gswap(f_sfm->pos, f_last);
	    }
	    if (s_map) {
		SKL*	fst = gsi[1].skl + 1;
		SKL*	lst = gsi[1].skl + gsi[1].skl->n;
		s_sfm->pos = fst->n;
		s_last = lst->n;
		if (sqs[0]->inex.sens) s_sfm->flag |= 0x10;
		if (sqs[1]->inex.sens) gswap(s_sfm->pos, s_last);
	    }
	    if (f_map && s_map) {
		f_sfm->flag |= 0x2; s_sfm->flag |= 0x2;
		if (s_sfm->flag & 0x10) f_sfm->flag |= 0x20;
		if (f_sfm->flag & 0x10) s_sfm->flag |= 0x20;
		if ((f_sfm->pos < s_sfm->pos) ^ sqs[1]->inex.sens) {
		    f_sfm->flag |= 0x40;
		    s_sfm->flag |= 0x80;
		} else {
		    s_sfm->flag |= 0x40;
		    f_sfm->flag |= 0x80;
		}
		f_sfm->tlen = s_last - f_last;
		s_sfm->tlen = -f_sfm->tlen;
		f_sfm->pnext = s_sfm->pos;
		s_sfm->pnext = f_sfm->pos;
		f_sfm->rnext = s_sfm->rnext = "=";
	    }
	}
	return (m);
}

inline int singleton(Seq** sqs, Gsinfo* gsi, VTYPE thr)
{
	return (gsi->skl && (!algmode.thr || gsi->fstat.val >= thr));
}

static int blkaln(Seq* sqs[], SrchBlk* bks, RANGE* rng, ThQueue* q)
{
	int	nparalog = 0;
	RANGE	grng = {0, 0};
	RANGE	wrng;
	RANGE	frng = {INT_MAX, 0};

	if (gsquery) {
	    sqs[0]->saverange(&grng);
	    int	n = (grng.right - grng.left) / 10;
	    frng.left = grng.left + n;
	    frng.right = grng.right - n;
	}
	if (QRYvsDB == AvsA || QRYvsDB == CvsC || QRYvsDB == GvsC)
	    nparalog = bks->finds(sqs);
	else if (QRYvsDB == GvsA)
	    nparalog = bks->findh(sqs);
	else {
	    if (algmode.lcl & 16) sqs[0]->exg_seq(1, 1);	// SWG local alignment
	    else	sqs[0]->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	    nparalog = bks->findblock(sqs);
	}
	if (nparalog == ERROR || algmode.nsa == MAP1_FORM || algmode.nsa == MAP2_FORM)
	    return (0);		// no alignment
#if USE_FULL_RANGE
// skip the next line if only the specified range is used as the query
	if (rng) sqs[0]->right = sqs[0]->right == rng->right? sqs[0]->len: rng->right;
#endif
	int	b_intr = sqs[1]->inex.intr;
	HalfGene	hfg(0, INT_MAX, sqs[0]->CdsNo);
	int	n_out = nparalog > 1? OutPrm.supself: 0;
	for (int n = 1; n <= nparalog && n_out < nparalog; ++n) {
	    if (n > 1) swapseq(sqs + 1, sqs + n);
	    if ((QRYvsDB == AvsA || QRYvsDB == CvsC) && OutPrm.supself) {
		int	sc = strcmp((*sqs[0]->sname)[0], (*sqs[1]->sname)[0]);  
		if (!sc || (OutPrm.supself > 1 && sc > 0)) continue;
	    }
	    sqs[1]->inex.intr = b_intr;
	    if (gsquery) {
		swapseq(sqs, sqs + 1);
		gswap(sqs[0]->CdsNo, sqs[1]->CdsNo);
		if (sqs[0]->jxt) {
		    delete[] sqs[1]->jxt;
		    sqs[1]->jxt = sqs[0]->jxt;
		    sqs[0]->jxt = 0;
		    sqs[0]->saverange(&wrng);
		    sqs[1]->restrange(&wrng);
		    sqs[0]->left = 0;
		    sqs[0]->right = sqs[0]->len;
		    sqs[0]->inex.ori = QRYvsDB == GvsA? 1: 3;
		}
		if (algmode.lcl & 16) sqs[0]->exg_seq(1, 1);
		else	sqs[0]->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	    }				// incompatible gene orientation
	    Gsinfo	GsI[4];
	    if (algmode.lcl & 16) sqs[1]->exg_seq(1, 1);
	    else	sqs[1]->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	    int	dir = ERROR;
	    if (sqs[0]->inex.ori & 1) {
		if (bks->pwd->DvsP != 3) genomicseq(sqs + 1, bks->pwd, false);
		if (pairedends) pairends(sqs, bks->pwd, GsI, sqs[no_seqs - 2]->len);
		else	dir = spalign2(sqs, bks->pwd, GsI);	// direct directions
	    }
	    if (algmode.lsg == 0 && sqs[0]->inex.ori == 3) sqs[0]->inex.ori = 1;
	    if ((sqs[0]->inex.ori & 2) && !GsI->intronless() &&
		!(pairedends && GsI[1].intronless())) {
		antiseq(sqs + 1);			// complementary oriections
		if (sqs[1]->jxt) revjxt(sqs[1]->jxt, sqs[1]->CdsNo); 
		if (bks->pwd->DvsP != 3) genomicseq(sqs + 1, bks->pwd, false);
		if (pairedends) {
		    sqs[no_seqs - 2]->comrev();
		    sqs[no_seqs - 1]->comrev();
		    pairends(sqs, bks->pwd, GsI + 2, sqs[no_seqs - 1]->len);
		    if (GsI[0].skl && GsI[1].skl && GsI[2].skl && GsI[3].skl) 
			dir = (GsI[2].scr + GsI[3].scr > GsI[0].scr + GsI[1].scr)? 2: 0;
		    else if (GsI[2].skl && GsI[3].skl) dir = 2;
		} else {
		    sqs[0]->comrev();
		    int	dir2 = spalign2(sqs, bks->pwd, GsI + 1);
		    if (dir2 == OK && dir == OK) dir = GsI[1].scr > GsI[0].scr;
		    else if (dir2 == OK) dir = 1;
		}
		if (dir == 0) {
		    if (pairedends) {
			sqs[no_seqs - 2]->comrev();
			sqs[no_seqs - 1]->comrev();
		    } else	sqs[0]->comrev();
		    antiseq(sqs + 1);
		}
	    }
	    delete sqs[1]->exin; sqs[1]->exin = 0;	// suppress Bounary output
	    if (dir < 0) continue;
	    Gsinfo*	gsinf = GsI + dir;
	    int	mp = pairedends? matepair(sqs, gsinf, bks->pwd->Vthr): 
		singleton(sqs, gsinf, bks->pwd->Vthr);
	    for (int k = 0; k < (pairedends? 2: 1); ++k, ++gsinf) {
		if (!(mp & (1 << k))) continue;
		if (pairedends) swapseq(sqs, sqs + no_seqs - 2 + k);
		int	suppress = 0;
		if (gsquery) {
		    int	r = gsinf->skl[0].n - 1;
		    int	sa = gsinf->skl[1].m;
		    int	sb = gsinf->skl[1].n;
		    int	eb = gsinf->skl[r].n;
		    int	db = sb - sqs[1]->left;
		    if (QRYvsDB == GvsA) sa *= 3;
		    if ((IntronPrm.tlmt < sa) && (db < IntronPrm.rlmt) &&
			(sqs[1]->left > IntronPrm.rlmt)) {
			if (eb > hfg.left) hfg.left = eb;
			if (algmode.lcl < 16) suppress |= 1;
		    }
		    sa = sqs[0]->right - gsinf->skl[r].m;
		    db = sqs[1]->right - eb;
		    if (QRYvsDB == GvsA) sa *= 3;
		    if ((IntronPrm.tlmt < sa) && (db < IntronPrm.rlmt) &&
			(sqs[1]->len - sqs[1]->right > IntronPrm.rlmt)) {
			if (sb < hfg.right) hfg.right = sb;
			if (algmode.lcl < 16) suppress |= 2;
		    }
		} else {
		    hfg.right = dir + 1;
		}
		if (!suppress && (!algmode.thr || gsinf->fstat.val >= bks->pwd->Vthr)) {
		    if (sqs[1]->left < frng.left) frng.left = sqs[1]->left;
		    if (sqs[1]->right > frng.right) frng.right = sqs[1]->right;
		    if (algmode.nsa == BED_FORM)
			gsinf->rscr = selfAlnScr(sqs[0], bks->pwd->simmtx);
#if M_THREAD
		    if (q) {
			pthread_mutex_lock(&q->mutex);
			alnoutput(sqs, gsinf);
			pthread_mutex_unlock(&q->mutex);
		    } else
#endif
			alnoutput(sqs, gsinf);
		    if (k == 0) ++n_out;
		}
		if (algmode.mlt == 1 && gsinf->eijnc) {
		    EISCR*	fst = gsinf->eijnc->begin();
		    sqs[0]->left = fst[0].rleft;
		    sqs[0]->right = fst[gsinf->noeij - 1].rright;
		}
		if (pairedends) swapseq(sqs, sqs + no_seqs - 2 + k);
	    } 
	    if (gsquery) {
		swapseq(sqs, sqs + 1);
		sqs[0]->restrange(&grng);
		gswap(sqs[1]->CdsNo, sqs[0]->CdsNo);
	    }
	    if (dir) {
		if (!pairedends) sqs[0]->comrev();
		antiseq(sqs + 1);
	    }
	}
	if (nparalog > 1) gswap(sqs[1], sqs[2]);
#if M_THREAD
	if (q && q->mfd) {
	    pthread_mutex_lock(&q->mutex);
	    hfg.sname = strrealloc(0, sqs[0]->sqname());
	    if (hfg.left > 0) hfg.left += IntronPrm.llmt;
	    else	hfg.left = frng.left;
	    if (hfg.right < INT_MAX) hfg.right -= IntronPrm.llmt;
	    else	hfg.right = frng.right;
	    hfg.left = sqs[0]->SiteNo(hfg.left);
	    hfg.right = sqs[0]->SiteNo(hfg.right);
	    q->mfd->write(&hfg);
	    pthread_mutex_unlock(&q->mutex);
	}
#endif
	return (hfg.right);
}

static	int	MinSegLen = INT_MAX;

static SrchBlk* getblkinf(Seq* sqs[], const char* dbs, MakeBlk* mb)
{
	Seq*	a = sqs[0];	// query
	Seq*	b = sqs[1];	// database
	int	ap = a->isprotein();
	int	bp = mb? mb->dbsmolc() == PROTEIN: ((dbs && (dbs == aadbs))? 1: b->isprotein());
	char	str[LINE_MAX];

	if (genomedb) {
	    b->inex.molc = DNA;
	    b->inex.intr = algmode.lsg;
	} else if (aadbs) {
	    b->inex.molc = PROTEIN;
	    if (a->isprotein()) algmode.lsg = 0;
	    a->inex.intr = algmode.lsg;
	} else if (cdnadb) {
	    b->inex.molc = DNA;
	    b->inex.intr = 0;
	    a->inex.intr = algmode.lsg;
	} else if (mb) {
	    b->inex.molc = mb->dbsmolc() == PROTEIN? PROTEIN: DNA;
	    b->inex.intr = bp? 0: algmode.lsg;
	    a->inex.intr = (b->inex.intr || ap)? 0: algmode.lsg;
	} else if (b->isprotein()) {
	    if (a->isprotein()) algmode.lsg = 0;
	    a->inex.intr = algmode.lsg;
	} else {
	    b->inex.intr = algmode.lsg;
	    a->inex.intr = 0;
	} 
	b->byte = b->many = 1;
	if (ap && bp)	{QRYvsDB = AvsA; algmode.lsg = 0;}
	else if (ap)	QRYvsDB = AvsG;
	else if (bp)	QRYvsDB = GvsA;
	else if (a->inex.intr && b->inex.intr) QRYvsDB = GvsG;
	else if (b->inex.intr) QRYvsDB = CvsG;
	else if (genomedb)	QRYvsDB = FvsG;
	else if (algmode.lsg) QRYvsDB = GvsC;
	else		QRYvsDB = CvsC;
	gsquery = QRYvsDB == GvsA || QRYvsDB == GvsC;
	if ((QRYvsDB == AvsG || QRYvsDB == GvsA) && algmode.lcl & 16)
	    algmode.lcl |= 3;
	if (gsquery) {
	    swapseq(sqs, sqs + 1);
	    gswap(ap, bp);
	}
	makeWlprms(prePwd(sqs));
	if (!ReadBlock && dbs) {
	    const char*	ext = 0;
	    switch (QRYvsDB) {
		case AvsA: case GvsA: ext = BKA_EXT; break;
		case GvsC: case CvsC: case GvsG: case FvsG: case CvsG:
		    ext = BKN_EXT; break;
		case AvsG: ext = BKP_EXT; break;
		case UNDF: fatal("Ilegal seq combination !\n"); break;
	    }
	    ReadBlock = path2dbf(str, dbs, ext);
	    if (!ReadBlock) fatal("Specify database (expected %s%s - mode %s) !\n", dbs, ext, QvsDstr[QRYvsDB]);
	}
	if ((QRYvsDB == CvsC || QRYvsDB == FvsG) && algmode.mns != 2) algmode.mns = 1;
	SrchBlk* bks = dbs? new SrchBlk(sqs, ReadBlock, genomedb):
		new SrchBlk(sqs, mb, b->inex.intr);
	if (MinQueryLen == 0) MinQueryLen = bks->MinQuery();
	if (algmode.mlt == 1) MinSegLen = MinQueryLen;
	if (ap) q_mns &= ~2;
	if (gsquery){
	    swapseq(sqs, sqs + 1);
	    gswap(ap, bp);
	}
	return bks;
}

static int quick4(Seq* sqs[], SrchBlk* bks, ThQueue* q)
{
	Seq*&	a = sqs[0];	// query
	if (bks->incompatible(a)) {
	    const char*	qn = a->sqname();
	    fatal("\nBad molecular type of %s !\n", *qn? qn: "query");
	}
	RANGE	orgrng;
	RANGE	covrng;

	a->inex.intr = (gsquery)? algmode.lsg: 0;
	RANGE*	kptrng = (gsquery)? 0: &orgrng;
	if (pairedends) {
	    a->inex.ori = q_mns;
	    sqs[0]->copyseq(sqs[no_seqs - 2], CPY_SEQ);
	    swapseq(sqs, sqs + no_seqs - 2);
	    sqs[1]->comrev();
	    sqs[1]->copyseq(sqs[no_seqs - 1], CPY_SEQ);
	    swapseq(sqs + 1, sqs + no_seqs - 1);
	    sqs[1]->catseq(sqs[0], spacer);		// 'N's for spacer
	} else {
	    if (q_mns == 1 || q_mns == 2) a->inex.ori = q_mns;
	    else if (a->isdrna() && !a->inex.intr)		// cDNA
		a->inex.ori = polyA.rmpolyA(a);
	}
	if (a->right - a->left < MinQueryLen) {
	    prompt("%s (%d - %d) %d is too short!\n", 
		a->sqname(), a->left, a->right, a->right - a->left);
	    return (0);
	}
#if TESTRAN
	a->rndseq();
#endif
	if (kptrng) sqs[0]->saverange(&orgrng);
	int	rbdry = blkaln(sqs, bks, kptrng, q);
	if (gsquery) return (rbdry);

// dispersed locations?

	sqs[0]->saverange(&covrng);
	if (rbdry && covrng.left - orgrng.left > MinSegLen) {
	    a->left = orgrng.left;
	    a->right = covrng.left;
	    blkaln(sqs, bks, 0, q);
	}
	if (rbdry && orgrng.right - covrng.right > MinSegLen)  {
	    a->left = covrng.right;
	    a->right = orgrng.right;
	    blkaln(sqs, bks, 0, q);
	}
	sqs[0]->restrange(&orgrng);
	return (0);
}

static void spaln_job(Seq* sqs[], void* prm, ThQueue* q)
{
	if (algmode.blk) (void) quick4(sqs, (SrchBlk*) prm, q);
	else {
	    sqs[0]->inex.ori = polyA.rmpolyA(sqs[0]);
	    if (q_mns == 1 || q_mns == 2) sqs[0]->inex.ori = q_mns;
	    match_2(sqs, (PwdB*) prm, q);
	}
}

static void seg_job(Seq** sqs, SeqServer* svr, SrchBlk* sbk)
{
	if (!sbk) fatal("No block inf !\n");
	if (QRYvsDB == GvsA) genomicseq(sqs, sbk->pwd, false);
	quick4(sqs, sbk);
	if (QRYvsDB == GvsC) return;
	antiseq(sqs);
	quick4(sqs, sbk);
	antiseq(sqs);
}

static void genomicseq(Seq** sqs, PwdB* pwd, bool reverse = false)
{
	sqs[0]->inex.intr = algmode.lsg;
	if (pwd->DvsP == 1) sqs[0]->nuc2tron();
	if (algmode.lsg) Intron53(sqs[0], pwd);
	if (!reverse || pwd->DvsP == 3) return;
	sqs[0]->comrev(sqs + 1);
	sqs[1]->setanti(sqs);
	if (pwd->DvsP == 1) sqs[1]->nuc2tron();
	if (algmode.lsg) Intron53(sqs[1], pwd);
}

static void all_in_func(Seq** sqs, SeqServer* svr, void* prm)
{
	int	nf = svr->input_form == IM_PARA? 1: 0;
	SrchBlk*	sbk = 0;
	PwdB*	pwd = (PwdB*) prm;
	if (algmode.blk) {
	    sbk = (SrchBlk*) prm;
	    if (!sbk->dbf) sbk->dbf = svr->target_dbf;
	    pwd = sbk->pwd;
	}
	do {
// perform the job
	    if (gsquery) seg_job(sqs, svr, sbk);
	    else	spaln_job(sqs, prm);
// read pair of sequences
	    if (svr->input_ns == 2) {
		switch (svr->nextseq(sqs[1], nf)) {
		    case IS_END: case IS_ERR: return;
		    default: 
// read genomic sequence
			genomicseq(sqs + 1, pwd, q_mns & 2);
			break;
		}
	    }
// read query sequence
	    switch (svr->nextseq(sqs[0], 0)) {
		case IS_END: return;
		case IS_ERR: continue;
		default: break;
	    }
	    if (pairedends) swapseq(sqs, sqs + 1);
	} while (true);
}

static void put_genome_entries()
{
	DbsDt*&	dbf = dbs_dt[0];
	fprintf(out_fd, "@HD\tVN:1.3\tSO:unsorted\n");
	for (INT i = 0; i < dbf->numidx; ++i) {
	    DbsRec*	rec = dbf->dbsrec(i);
	    fprintf(out_fd, "@SQ\tSN:%s\tLN:%ld\n", 
		dbf->entname(i), (long) rec->seqlen);
	}
}

#if M_THREAD
ThQueue::ThQueue(Seq** sqs) : sinp(sqs)
{
	sque = new Seq*[max_queue_num];
	initseq(sque, max_queue_num);
	mfd = (gsquery)?
	    new Mfile(sizeof(HalfGene)): 0;
	rp = wp = remain = 0;
	pthread_mutex_init(&mutex, 0);
	pthread_cond_init(&not_full, 0);
	pthread_cond_init(&not_empty, 0);
}

void ThQueue::enqueue(Seq** sqs, int n)
{
	pthread_mutex_lock(&mutex);
	while (remain == max_queue_num)
	    pthread_cond_wait(&not_full, &mutex);
	while (--n >= 0) {
	    if (sqs) {
#if QDEBUG
		fprintf(stderr, "e%d: %s %d %d %d\n", n, sqs[n]->sqname(),
		sqs[n]->sid, sqs[n]->len, sqs[n]->many);
#endif
		if (sqs[n]->sid > 0) swapseq(sque + wp, sqs + n);
	    } else	sque[wp]->refresh(0);
	    ++wp; ++remain;
	    if (wp == max_queue_num) wp = 0;
	}
	pthread_cond_signal(&not_empty);
	pthread_mutex_unlock(&mutex);
}

void ThQueue::dequeue(Seq** sqs, int n)
{
	pthread_mutex_lock(&mutex);
	while (remain == 0)
	     pthread_cond_wait(&not_empty, &mutex);
	while (--n >= 0) {
#if QDEBUG
	    if (sque[rp]->many)	{swapseq(sqs + n, sque + rp);
	    fprintf(stderr, "d%d: %s %d %d %d\n", n, sqs[n]->sqname(),
	    sqs[n]->sid, sqs[n]->len, sqs[n]->many);}
#else
	    if (sque[rp]->many)	swapseq(sqs + n, sque + rp);
#endif
	    else	sqs[n]->refresh(0);
	    ++rp; --remain;
	    if (rp == max_queue_num) rp = 0;
	}
	pthread_cond_signal(&not_full);
	pthread_mutex_unlock(&mutex);
}

// enqueue the (if necessary segmented) query

void ThQueue::putqueue()
{
	int	slen = sinp[0]->len;

	if (slen <= g_segment) {
	    enqueue(sinp);
	    return;
	}
	float	r = 0.9;	// overlap fraction
	int	nseg = int((slen - r * g_segment) / ((1 - r) * g_segment)) + 1;
	int	lseg = int(slen / (nseg - (nseg - 1) * r));
	int	step = int(lseg * r + 1);
	RANGE	rng = {0};
	for (int segno = 0; rng.left < slen; rng.left += step) {
	    rng.right = rng.left + lseg;
	    if (slen < rng.right) rng.right = slen;
	    (*sinp)->restrange(&rng);
	    (*sinp)->CdsNo = ++segno;
	    enqueue(sinp);
	}
}

static void* master_func(void* arg)
{
	mast_arg_t*	targ = (mast_arg_t*) arg;
	ThQueue*	q = targ->q;
	Seq**	fsd = q->sinp;
	int	nf = targ->svr->input_form == IM_PARA? 1: 0;

	 while (true) {
// read genomic sequence
	    if (targ->svr->input_ns == 2 && 
		(targ->svr->nextseq(fsd[1], nf) != IS_OK)) break;
// read query and push it in the queue
	    InSt	inst = targ->svr->nextseq(fsd[0], 0);
	    if (inst == IS_END) break;
	    if (inst == IS_OK) {
		if (QRYvsDB == GvsA) q->putqueue();
		else {
		    if (pairedends) swapseq(fsd, fsd + 1);
		    q->enqueue(fsd, targ->svr->input_ns);
		}
	    }
	}
	for (int n = 0; n < thread_num; ++n)
	    q->enqueue(0, targ->svr->input_ns);
	return (void*) 0;
}

static int hcmp(HalfGene* a, HalfGene* b)
{
	int	c = strcmp(a->sname, b->sname);

	if (c) return(c);
	if (a->segno < 0) {
	    if (b->segno > 0) return (1);
	    return (b->segno - a->segno);
	} else {
	    if (b->segno < 0) return (-1);
	    return (a->segno - b->segno);
	}
}

static void* mistress_func(void* arg)
{
	mist_arg_t*	targ = (mist_arg_t*) arg;
	ThQueue*		q = targ->q;
	Seq**		fsd = q->sinp;
	HalfGene*	whf = targ->hfg;
	HalfGene*	thf = whf + targ->nhf - 1;
	char		str[LINE_MAX];

	qsort((UPTR) whf, (INT) targ->nhf, sizeof(HalfGene), (CMPF) hcmp);
	for ( ; whf < thf; ++whf) {
	    if (strcmp(whf->sname, whf[1].sname)) continue;
	    int	l = whf[1].segno - whf->segno;
	    if (l == 1) { 
		sprintf(str, "%s %d %d", whf->sname, whf->right, whf[1].left);
		if ((*fsd)->getseq(str)) q->enqueue(fsd);
	    } else if (l == -1) {
		sprintf(str, "%s %d %d <", whf->sname, whf[1].left, whf->right);
		if ((*fsd)->getseq(str)) q->enqueue(fsd);
	    }
	}
	for (int l = 0; l < thread_num; ++l) q->enqueue(0);
	return (void*) 0;
}

static void* worker_func(void* arg)
{
	thread_arg_t*	targ = (thread_arg_t*) arg;
	SrchBlk*	sbk = algmode.blk? (SrchBlk*) targ->pwd: 0;
	PwdB*	pwd = sbk? sbk->pwd: (PwdB*) targ->pwd;

#ifdef __CPU_SET
	cpu_set_t	mask;
	__CPU_ZERO(&mask);
	__CPU_SET(targ->cpuid, &mask);
	if (sched_setaffinity(0, sizeof(mask), &mask) == -1)
	    prompt("Warning: faild to set CPU affinity !\n");
#endif

	while (true) {
	    targ->q->dequeue(targ->seqs, targ->svr->input_ns);
	    if (targ->seqs[0]->many == 0) break;
	    if (targ->svr->input_ns == 2 && targ->seqs[1]->many == 0) break;
	    if (gsquery) {
		seg_job(targ->seqs, targ->svr, sbk);
	    } else {
		if (targ->svr->input_ns == 2)
		    genomicseq(targ->seqs + 1, pwd, q_mns & 2);
		(void) spaln_job(targ->seqs, targ->pwd, targ->q);
	    }
	}
	return (void*) 0;
}

static void MasterWorker(Seq** sqs, SeqServer* svr, void* prm)
{
	mast_arg_t	maarg;
	mist_arg_t	miarg;
	pthread_t	master;

	max_queue_num = (max_queue_num + svr->input_ns - 1) / svr->input_ns * svr->input_ns;
	thread_arg_t*	targ = new thread_arg_t[thread_num];
	pthread_t*	worker = new pthread_t[thread_num];
	SrchBlk*	primaty = svr->target_dbf? (SrchBlk*) prm: 0;

	ThQueue	q(sqs);
	maarg.q = &q;
	maarg.svr = svr;
	targ[0].seqs = new Seq*[no_seqs * thread_num];
	initseq(targ[0].seqs, no_seqs * thread_num);
	for (int n = 0; n < thread_num; ++n) {
	    targ[n].q = &q;
	    targ[n].cpuid = n % cpu_num;
	    if (n > 0) targ[n].seqs = targ[n - 1].seqs + no_seqs;
	    targ[n].seqs[0]->inex = sqs[0]->inex;
	    if (svr->target_dbf) {
		targ[n].seqs[1]->inex = sqs[1]->inex;
		if (n) {
		    DbsDt*	dbf = new DbsDt(*svr->target_dbf);
		    dbf->fseq = svr->target_dbf->dbsfopen();
		    targ[n].pwd = (void*) new SrchBlk(primaty, dbf);
		} else {
		    primaty->reset(svr->target_dbf);
		    targ[n].pwd = (void*) primaty;
		}
	    } else {
		if (svr->input_form == IM_SNGL) {
		    sqs[1]->aliaseq(targ[n].seqs[1]);
		    sqs[2]->aliaseq(targ[n].seqs[2]);
		    targ[n].seqs[1]->setanti(targ[n].seqs + 2);
		    targ[n].seqs[2]->setanti(targ[n].seqs + 1);
		}
		targ[n].pwd = (void*) prm;
	    }
	    targ[n].svr = svr;
	}
	if (QRYvsDB == GvsA) q.putqueue();
	else	q.enqueue(sqs, svr->input_ns);
	pthread_create(&master, 0, master_func, (void*) &maarg);
	for (int n = 0; n < thread_num; ++n)
	    pthread_create(worker + n, 0, worker_func, (void*) (targ + n));
	for (int n = 0; n < thread_num; ++n)
	    pthread_join(worker[n], 0);
	if (q.mfd) {
	    cleanseq(targ[0].seqs, no_seqs * thread_num);
	    cleanseq(q.sque, max_queue_num);
	    miarg.nhf = q.mfd->size();
	    miarg.hfg = (HalfGene*) q.mfd->flush();
	    delete q.mfd; q.mfd = 0;
	    if (miarg.nhf > 1) {
		miarg.q = &q;
		pthread_create(&master, 0, mistress_func, (void*) &miarg);
		for (int n = 0; n < thread_num; ++n)
		    pthread_create(worker + n, 0, worker_func, (void*) (targ + n));
		for (int n = 0; n < thread_num; ++n)
		    pthread_join(worker[n], 0);
	    }
	    delete[] miarg.hfg;
	}
	clearseq(targ[0].seqs, no_seqs * thread_num);
	clearseq(q.sque, max_queue_num);
	if (svr->target_dbf) {
	    for (int n = 1; n < thread_num; ++n) {
		SrchBlk*	sbk = (SrchBlk*) targ[n].pwd;
		sbk->dbf->clean();
		delete sbk;
	    }
	}
	delete[] targ[0].seqs;
	delete[] targ;
	delete[] worker;
}
#endif

static void setdefparam()
{
	algmode.lcl = 15;	// default semi-global
	alprm.ls = 2;		// affine gap penalty
	algmode.lsg = 1;	// spliced alignment
	algmode.qck = 3;	// three level HSP search
	algmode.mlt = 0;	// sigle alignment for each query
	algmode.mns = 3;	// both strand and both direction
	algmode.thr = 1;	// filter out weak matches
	setalgmode(4, 0);	// rich exon info
	setNpam(4, -6);		// default int mismatch score
	setpam(DefPam, 0);	// change to this pam value
	setpam(WlpPam, 1);	// pam of secondary sim matrix
	setorf(DefOrfLen, 2);	// orf length and AG< .. >GT ends
	OutPrm.MaxOut = 1;
	OutPrm.SkipLongGap = 1;	// suppress display of long gaps
	OutPrm.fastanno = 1;	// add annotation in fasta output
#if !FVAL
	alprm.scale = 10;
#endif
}

int main(int argc, const char** argv)
{
const	char*	messg = 0;
const	char*	insuf = "No input seq file !\n";

	if (argc <= 1) usage(insuf);	// no argments
	setdefparam();
	optimize(GLOBAL, MAXIMUM);
	int	n = getoption(argc, argv);
	spb_fact();
	argv += n;
	argc -= n;
	if (WriteMolc) {
	    if (!TgtMolc) fatal("Specify molecular type by -K[A|D|P] option !\n");
	    MakeBlk*	mb = 0;
	    if (TgtMolc == PROTEIN && alprm2.spb > 0.) {
		SeqServer	svr(argc, argv, IM_SNGL, 0, TgtMolc);
		mb = makeblock(&svr);
	    } else {
		mb = makeblock(argc, argv, TgtMolc);
	    }
	    mb->WriteBlkInfo();
	    mb->delete_dbf();
	    delete mb;
	    resetSimmtxes(true);
	    return (0);
	}
	if (!catalog && argc <= 0) usage(insuf);
#if M_THREAD
	cpu_num = sysconf(_SC_NPROCESSORS_CONF);
	if (thread_num < 0) thread_num = cpu_num;
	if (max_queue_num == 0) max_queue_num = int(FACT_QUEUE * thread_num);
#endif	// M_THREAD
	pairedends = genomedb && input_form != IM_SNGL;
	if (QryMolc == TRON) QryMolc = PROTEIN;
	if (!algmode.mlt) OutPrm.MaxOut = 1;
	no_seqs = OutPrm.MaxOut + (pairedends? 4: 2);
	setup_output(algmode.nsa);	// output format
	Seq**	seqs = new Seq*[no_seqs];
	initseq(seqs, no_seqs);	// 0: query, 1: genomic, 2: reverse
	MakeBlk*	mb = 0;
const	char*	dbs = genomedb? genomedb: (aadbs? aadbs: cdnadb);
	if (algmode.blk && !dbs) {	// make db from 1st file
	    if (!argv[1]) usage(insuf);
	    if (!TgtMolc) {
		TgtMolc = infermolc(argv[0]);
	    	QryMolc = infermolc(argv[1]);
		if (TgtMolc != PROTEIN && QryMolc == PROTEIN)
		    TgtMolc = TRON;
	    }
	    if (TgtMolc == PROTEIN && alprm2.spb > 0.) {
		SeqServer	svr(1, argv++, IM_SNGL, 0, TgtMolc);
		mb = makeblock(&svr);
	    } else {
		mb = makeblock(1, argv++, TgtMolc);
	    }
	    if (--argc <= 0 && !catalog) {
		clearseq(seqs, no_seqs);
		delete[] seqs;
		usage(insuf);
	    }
	}
	SeqServer	svr(argc, argv, input_form, catalog, QryMolc);
	if (algmode.blk) {
	    if (svr.nextseq(seqs[0], 0) > IS_ERR) {
		messg = "Can't open query !\n";
		goto postproc;
	    }
	    if (!OutPrm.out_file && algmode.nsa == BIN_FORM)
		OutPrm.out_file = seqs[0]->spath;
	    SrchBlk*	bprm = getblkinf(seqs, dbs, mb);
	    delete mb;
	    if (seqs[0]->inex.intr || seqs[1]->inex.intr) makeStdSig53();
	    set_max_extend_gene_rng(def_max_extend_gene_rng);
	    if (algmode.nsa == SAM_FORM) put_genome_entries();
	    if (pairedends && svr.nextseq(seqs[1], svr.input_form == IM_PARA) != IS_OK)
		usage("No mate pair !\n");
#if M_THREAD
	    if (thread_num) MasterWorker(seqs, &svr, (void*) bprm); else
#endif
	    all_in_func(seqs, &svr, (void*) bprm);
	    delete bprm;
	} else {
	    if (svr.nextseq(seqs[1], 1) == IS_END) {
		messg = "Can't open genomic sequence !\n";
		goto postproc;
	    }
	    if (svr.nextseq(seqs[0], 0) != IS_OK) {
		messg = "Can't open query !\n";
		goto postproc;
	    }
	    PwdB*	pwd = SetUpPwd(seqs);
	    if (seqs[0]->inex.intr || seqs[1]->inex.intr) makeStdSig53();
	    set_max_extend_gene_rng(0);
	    genomicseq(seqs + 1, pwd, true);
#if M_THREAD
	    if (thread_num) MasterWorker(seqs, &svr, (void*) pwd); else
#endif
	    all_in_func(seqs, &svr, (void*) pwd);
	    delete pwd;
	}

	eraWlprms();
	eraStrPhrases();
postproc:
	EraStdSig53();
	closeGeneRecord();
	clearseq(seqs, no_seqs);
	delete[] seqs;
	EraDbsDt();
	close_output();
	resetSimmtxes(true);
	if (messg) usage(messg);
	return (0);
}
