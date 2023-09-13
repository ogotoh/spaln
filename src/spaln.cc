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
static	const	int	InsPam = 100;	// intra species
static	const	int	CrsPam = 150;	// cross species
static	const	int	WlpPam = 50;	// HSP search
static	const	int	DefOrfLen = 75;
static	const	int	def_max_extend_gene_rng = 3;
static	const	long	MaxWordNoSpace = 1 << 30;	// 1GB

enum QvsD {UNDF, AvsA, AvsG, CvsC, CvsG, FvsG, GvsA, GvsC, GvsG};

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

#if M_THREAD

class ThQueue {
	int	rp, wp;
	int	remain;
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
#endif	// M_THREAD

static	void	usage(const char* messg);
static	int	getoption(int argc, const char** argv);
static	void	readargs();
static	RANGE*	skl2exrng(SKL* skl);
static	int	spalign2(Seq* sqs[], PwdB* pwd, Gsinfo* GsI, int ori = 1);
static	SrchBlk*	getblkinf(Seq* sqs[], const char* dbs, MakeBlk* mb);
static	void	seg_job(Seq** sqs, SeqServer* svr, SrchBlk* sbk);
static	void	all_in_func(Seq** sqs, SeqServer* svr, void* prm);
static	void	setdefparam();
static	PwdB*	SetUpPwd(Seq* sqs[]);
static	int	blkaln(Seq* sqs[], SrchBlk* bks, RANGE* rng, ThQueue* q);
static	int	match_2(Seq* sqs[], PwdB* pwd, ThQueue* q = 0);
static	int	quick4(Seq* sqs[], SrchBlk* bks, ThQueue* q = 0);
static	void	spaln_job(Seq* sqs[], void* prm, ThQueue* q = 0);
static	void	genomicseq(Seq** sqs, PwdB* pwd, int ori = 1);

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
static	int	no_seqs = 3;
static	bool	gsquery = QRYvsDB == GvsA || QRYvsDB == GvsC;
static	const	char*	version = "3.0.0";
static	const	int	date = 230912;
static	AlnOutModes	outputs;

static void usage(const char* messg)
{
#if __AVX512BW__
const	char*	arch = "_AVX512";
#elif __AVX2__
const	char*	arch = "_AVX2";
#elif __SSE4_1__
const	char*	arch = "_SSE4.1";
#else
const	char*	arch = "";
#endif
	if (messg && *messg) fputs(messg, stderr);
	fprintf(stderr, "\n*** SPALN%s version %s <%d> ***\n\n", 
	    arch, version, date);
	fputs("Usage:\n", stderr);
	fputs("spaln -W[Genome.bkn] -KD [W_Options] Genome.mfa\t(to write block inf.)\n", stderr);
	fputs("spaln -W[Genome.bkp] -KP [W_Options] Genome.mfa\t(to write block inf.)\n", stderr);
	fputs("spaln -W[AAdb.bka] -KA [W_Options] AAdb.faa\t(to write aa db inf.)\n", stderr);
	fputs("spaln -W [Genome.mfa|AAdb.faa]\t(alternative to makdbs.)\n", stderr);
	fputs("spaln [R_options] genomic_segment cDNA.fa\t(to align)\n", stderr);
	fputs("spaln [R_options] genomic_segment protein.fa\t(to align)\n", stderr);
	fputs("spaln [R_options] -dGenome cDNA.fa\t(to map & align)\n", stderr);
	fputs("spaln [R_options] -dGenome protein.fa\t(to map & align)\n", stderr);
	fputs("spaln [R_options] -aAAdb genomic_segment.fa\t(to search aa database & align)\n", stderr);
	fputs("spaln [R_options] -aAAdb protein.fa\t(to search aa database)\n", stderr);
	fputs("\nin the following, # = integer or real number; $ = string; default in ()\n\n", stderr);
	fputs("W_Options:\n", stderr);
	fputs("\t-XC#\tnumber of bit patterns < 6 (1)\n", stderr);
	fputs("\t-XG#\tMaximum expected gene size (inferred from genome|db size)\n", stderr);
	fputs("\t-Xk#\tWord size (inferred from genome|db size)\n", stderr);
	fputs("\t-Xb#\tBlock size (inferred from genome|db size)\n", stderr);
	fputs("\t-Xa#\tAbundance factor (10)\n", stderr);
	fputs("\t-Xr#\tMinimum ORF length with -KP (30))\n", stderr);
	fputs("\t-g\tgzipped output\n", stderr);
	fputs("\t-t#\tMutli-thread operation with # threads\n\n", stderr);
	fputs("R_Options (representatives):\n", stderr);
	fputs("\t-A[0-3]\t0: scalar, 1..3: simd; 1: rigorous, 2: intermediate, 3: fast\n", stderr);
	fputs("\t-H#\tMinimum score for report (35)\n", stderr);
	fputs("\t-L or -LS or -L#\tsemi-global or local alignment (-L)\n", stderr);
	fputs("\t-M#[,#2]\tNumber of outputs per query (1) (4 if # is omitted)\n", stderr);
	fputs("\t\t#2 (4) specifies the max number of candidate loci\n", stderr);
	fputs("\t\tThis option is effective only for map-and-align modes\n", stderr);
	fputs("\t-O#[,#2,..] (GvsA|C)\t0:Gff3_gene; 1:alignment; 2:Gff3_match; 3:Bed; 4:exon-inf;\n", stderr);
	fputs("\t\t\t5:intron-inf; 6:cDNA; 7:translated; 8:block-only;\n", stderr); 
	fputs("\t\t\t10:SAM; 12:binary; 15:query+GS (4)\n", stderr);
	fputs("\t-O#[,#2,..] (AvsA)\t0:statistics; 1:alignment; 2:Sugar; 3:Psl; 4:XYL;\n", stderr);
	fputs("\t\t\t5:srat+XYL; 8:Cigar; 9:Vulgar; 10:SAM; (4)\n", stderr); 
	fputs("\t-Q#\t0:DP; 1-3:HSP-Search; 4-7; Block-Search (3)\n", stderr);
	fputs("\t-R$\tRead block information file *.bkn, *.bkp or *.bka\n", stderr);
	fputs("\t-S#\tOrientation. 0:annotation; 1:forward; 2:reverse; 3:both (3)\n", stderr);
	fputs("\t-T$\tSubdirectory where species-specific parameters reside\n", stderr);
	fputs("\t-a$\tSpecify AAdb. Must run `makeidx.pl -ia' breforehand\n", stderr);
	fputs("\t-A$\tSame as -a but db sequences are stored in memory\n", stderr);
	fputs("\t-d$\tSpecify genome. Must run `makeidx.pl -i[n|p]' breforehand\n", stderr);
	fputs("\t-D$\tSame as -d but db sequences are stored in memory\n", stderr);
	fputs("\t-g\tgzipped output in combination with -O12\n", stderr);
	fputs("\t-l#\tNumber of characters per line in alignment (60)\n", stderr);
	fputs("\t-o$\tFile/directory/prefix where results are written (stdout)\n", stderr);
	fputs("\t-pa#\tRemove 3' poly A >= # (0: don't remove)\n", stderr);
	fputs("\t-pF\tOutput full Fasta entry name\n", stderr);
	fputs("\t-pj\tSuppress splice junction information with -O[6|7]\n", stderr);
	fputs("\t-pn\tRetain existing output file\n", stderr);
	fputs("\t-po\tOverwrite existing output file\n", stderr);
	fputs("\t-pw\tReport results even if the score is below the threshold\n", stderr);
	fputs("\t-pT\tExclude termination codon from CDS\n", stderr);
	fputs("\t-r$\tReport information about block data file\n", stderr);
	fputs("\t-u#\tGap-extension penalty (3)\n", stderr);
	fputs("\t-v#\tGap-open penalty (8)\n", stderr);
	fputs("\t-w#\tBand width for DP matrix scan (100)\n", stderr);
	fputs("\t-t[#]\tMutli-thread operation with # threads\n", stderr);
	fputs("\t-ya#\tStringency of splice site. 0->3:strong->weak\n", stderr);
	fputs("\t-yl3\tDdouble affine gap penalty\n", stderr);
	fputs("\t-ym#\tNucleotide match score (2)\n", stderr);
	fputs("\t-yn#\tNucleotide mismatch score (-6)\n", stderr);
	fputs("\t-yo#\tPenalty for a premature termination codon (100)\n", stderr);
	fputs("\t-yx#\tPenalty for a frame shift error (100)\n", stderr);
	fputs("\t-yy#\tWeight for splice site signal (8)\n", stderr);
	fputs("\t-yz#\tWeight for coding potential (2)\n", stderr);
	fputs("\t-yB#\tWeight for branch point signal (0)\n", stderr);
	fputs("\t-yI$\tIntron length distribution\n", stderr);
	fputs("\t-yL#\tMinimum expected length of intron (30)\n", stderr);
	fputs("\t-yS[#]\tUse species-specific parameter set (0.0/0.5)\n", stderr);
	fputs("\t-yX0\tDon't use parameter set for cross-species comparison\n", stderr);
	fputs("\t-yZ#\tWeight for intron potential (0)\n", stderr);
	fputs("\t-XG#\tReset maximum expected gene size, suffix k or M is effective\n", stderr);
	fputs("\nExamples:\n", stderr);
	fputs("\tspaln -W -KP -E -t4 dictdisc_g.gf\n", stderr);
	fputs("\tspaln -W -KA -Xk5 Swiss.faa\n", stderr);
	fputs("\tspaln -O -LS \'chr1.fa 10001 40000\' cdna.nfa\n", stderr);
	fputs("\tspaln -Q0,1,7 -t10 -TTetrapod -XG2M -ommu/ -dmus_musc_g hspcdna.nfa\n", stderr);
	fputs("\tspaln -Q7 -O5 -t10 -Tdictdics -ddictdisc_g [-E] \'dictdisc.faa (101 200)\' > ddi.intron\n", stderr);
	fputs("\tspaln -Q7 -O0 -t10 -Tdictdics -aSwiss \'chr1.nfa 200001 210000\' > Chr1_200-210K.gff\n", stderr);
	fputs("\tspaln -Q4 -O0 -t10 -M10 -aSwiss dictdisc.faa > dictdisc.alignment_score\n", stderr);
	exit(1);
}

static int getoption(int argc, const char** argv)
{
const	char**	argbs = argv;
	DbsDt**	dbs = dbs_dt;

	while (--argc > 0 && **++argv == OPTCHAR) {
const	    char*	opt = argv[0] + 1;
	    int	c = *opt;
	    if (!c) break;
const	    char*	val = argv[0] + 2;
	    int	k = 0;
	    int	rv = 0;
	    switch (c) {
		case '?': case 'h': usage(0);
		case 'a':
		    algmode.dim = 0;
		    if ((val = getarg(argc, argv))) {
			delete *dbs;
			*dbs = new DbsDt(aadbs = val);
		    }
		    break;
		case 'A':
		    if ((val = getarg(argc, argv))) {
			if (isdigit(*val))
			    algmode.alg = atoi(val);
			else {
			    algmode.dim = 1;	// on memory dbs
			    delete *dbs;
			    *dbs = new DbsDt(aadbs = val);
			}
		    }
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
		case 'g': OutPrm.gzipped = 1; break;
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
			default:  input_form = IM_SNGL; break;
		    }
		    if ((val = strchr(val, ':'))) catalog = val + 1;
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
			case 'X': case 'x':
			    ignoreamb = true;
			case 'P': case 'p':
			    QryMolc = PROTEIN; 
			    TgtMolc = TRON; break;
			case 'N': case 'n': 
			    ignoreamb = true;
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
			if (OutPrm.MaxOut > 1) algmode.mlt = 2;
			if ((val = strchr(val, '.')))
			    OutPrm.MaxOut2 = atoi(++val);
		    } else {
			OutPrm.MaxOut = OutPrm.MaxOut2 = MAX_PARALOG;
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
			outputs.getopt(val);
		    break;
		case 'o':
		    if ((val = getarg(argc, argv)))
			OutPrm.out_file = val;
		    break;
		case 'p':
		    switch (*val) {
			case 'a': polyA.setthr(val + 1); break;
			case 'd': OutPrm.descrp = 1; break;
			case 'D': OutPrm.debug = 1; break;
			case 'e': OutPrm.trimend = !OutPrm.trimend; break;
			case 'f': OutPrm.deflbl = 1; break;
			case 'F': OutPrm.full_name = 1; break;
			case 'i': OutPrm.ColorEij = 1; break;
			case 'j': OutPrm.spjinf = 1 - OutPrm.spjinf; break;
			case 'n': OutPrm.overwrite = 2; break;
			case 'o': OutPrm.overwrite = 1; break;
			case 'q': setprompt(0, 0); break;
			case 'Q': setprompt(0, 1);  break;
			case 'r': algmode.mlt = 1; break;
			case 't': OutPrm.deflbl = 2; break;
			case 'w': algmode.thr = 0; break;
			case 'x': OutPrm.supself = (val[1] == '2')? 2: 1; break;
			case 'J': case 'm': case 'u': case 'v': 
			    setexprm_z(argc, argv); break;
			case 'T': OutPrm.supTcodon = 1; break;
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
			ftable.setpath(val, gnm2tab);
		    readargs();
		    break;
		case 'u': case 'v': case 'w': 
		    readalprm(argc, argv); break;
		case 'U': algmode.lsg = 0; break;
		case 'V':
		    if ((val = getarg(argc, argv)))
			setVmfSpace(ktol(val));
		    break;
		case 'W': 
		    if (!*val) val = argv[1];
		    rv = setQ4prm(opt, val);
		    if (rv && argc > 1) {++argv; --argc;}
		    WriteMolc = true;
		    alprm2.spb = 0;		// 
		    break;
		case 'x': setexprm_x(argc, argv); break;
		case 'X': 
		    if (!*++val) val = argv[1];
		    rv = setQ4prm(opt + 1, val);
		    if (rv) {++argv; --argc;}
		    break;
		case 'y': readalprm(argc, argv, 2); break;
		case 'z': setexprm_z(argc, argv); break;
		default: break;
	    }
	    if (*dbs && (dbs + 1 < dbs_dt + MAX_DBS)) ++dbs;
	}
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

void AlnOutModes::alnoutput(Seq** sqs, Gsinfo* GsI)
{
	int	print2Skip = 0;
const	bool	swp = sqs[1]->inex.intr || 
		(sqs[0]->isprotein() && sqs[1]->istron());
	Seq*	gene = 0;
	GAPS*	gaps[2] = {0, 0};

	if (swp) {
	    std::swap(sqs[0], sqs[1]);
	    swapskl(GsI->skl);
	    gene = sqs[0];
	    delete[] gene->exons;
	    gene->exons = 0;
	}
	GsI->setprefix(prefix);
	if (gene && GsI->eijnc) {	// 
	    int	gl = GsI->eijnc->genleft();
	    int	gr = GsI->eijnc->genright();
	    if (gl > gr) std::swap(gl, gr);
	    if (gl < gene->left || gr > gene->right)
		prompt("Bad gene coord: %s %d %d %d %d %c %s\n",
		    gene->sqname(), gene->left, gl, gr, gene->right, 
		    gene->inex.sens? '<': ' ', sqs[1]->sqname());
	}
	for (int n = 0; n < n_out_modes; ++n) {
	  if (!fds[n]) continue;
	  if (QRYvsDB == AvsG || QRYvsDB == GvsA) { // DNA vs protein 
	    if (out_mode[n] == ALN_FORM) {
		if (gene) {
		    fphseqs((const Seq**) sqs, 2, fds[n]);
		    GBcdsForm(GsI->CDSrng, gene, fds[n]);
		    print2Skip = 1;
		}
		if (getlpw()) {		// print alignment
		    skl2gaps3(gaps, GsI->skl, 2);
		    if (gene) gene->exons = GsI->eiscrunfold(gaps[0]);
		    unfoldgap(gaps[0], 1);
		    unfoldgap(gaps[1], 3);
		    GsI->print2(sqs, (const GAPS**) gaps, 
			 1, 1, print2Skip, fds[n]);
		} else
		    GsI->repalninf(sqs, out_mode[n], fds[n]);
	    } else if (gene) {
		GsI->printgene(sqs, out_mode[n], fds[n]);
	    }
	  } else if (QRYvsDB == CvsG || QRYvsDB == GvsC) {	// genome vs cDNA
	    if (out_mode[n] == ALN_FORM || (out_mode[n] == ITN_FORM && !gene)) {
		if (getlpw()) {
		    skl2gaps(gaps, GsI->skl);
		    if (gene) gene->exons = GsI->eiscrunfold(gaps[0]);
		    toimage(gaps, 2);
		}
		if (gene) {
		    fphseqs((const Seq**) sqs, 2, fds[n]);
		    GBcdsForm(GsI->CDSrng, gene, fds[n]);
		    print2Skip = 1;
		}
		if (getlpw())
		    GsI->print2(sqs, (const GAPS**) gaps, 
			1, 1, print2Skip, fds[n]);
		else
		    GsI->repalninf(sqs, out_mode[n], fds[n]);
	    } else if (gene) {
		GsI->printgene(sqs, out_mode[n], fds[n]);
	    }
	  } else {
	    if (out_mode[n] == ALN_FORM) {
		if (getlpw()) {
		    skl2gaps(gaps, GsI->skl);
		    toimage(gaps, 2);
		    GsI->print2(sqs, (const GAPS**) gaps, 
			1, 1, print2Skip, fds[n]);
		} else {
		    GsI->repalninf(sqs, out_mode[n], fds[n]);
		}
	    } else {
		GsI->repalninf(sqs, out_mode[n], fds[n]);
	    }
	  }
	}
	if (gene) {delete[] gene->exons; gene->exons = 0;}
	if (swp)	{
	    std::swap(sqs[0], sqs[1]);
	    swapskl(GsI->skl);
	}
	delete[] gaps[0];
	delete[] gaps[1];
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
		skl->n - skl[-1].n > IntronPrm.minl) {
		    (wrng++)->right = skl[-1].n;
		    wrng->left = skl->n;
	    }
	}
	(wrng++)->right = skl[-1].n;	// 3' end
	*wrng = endrng;
	rng->left = ++wrng - rng;
	return (rng);
}

static int spalign2(Seq* sqs[], PwdB* pwd, Gsinfo* GsI, int ori)
{
	if (pwd->DvsP != 3)
	    genomicseq(sqs + 1, pwd, ori);
	switch (QRYvsDB) {
	  case AvsG: 
	  case GvsA: GsI->skl = alignH_ng((const Seq**) sqs, pwd, GsI); break;
	  case GvsC:
	  case CvsG: GsI->skl = alignS_ng(sqs, pwd, GsI, ori); break;
	  case FvsG: 
	  case AvsA: 
	  case CvsC: GsI->skl = alignB_ng((const Seq**) sqs, pwd, GsI); break;
	  default:
		prompt("Improper combination of db and query!\n");
		return (ERROR);
	}
	if (GsI->skl) {
	    if (GsI->skl->n == 0) return (ERROR);
	    if (QRYvsDB == AvsA || QRYvsDB == CvsC || 
		QRYvsDB == GvsG || QRYvsDB == FvsG) {
		GsI->scr = skl_rngB_ng((const Seq**) sqs, GsI, pwd);
	    } else if (GsI->skl->m & AlgnTrb) {
		if (QRYvsDB == CvsG || QRYvsDB == GvsC)
			GsI->scr = skl_rngS_ng((const Seq**) sqs, GsI, pwd);
		else	GsI->scr = skl_rngH_ng((const Seq**) sqs, GsI, pwd);
		GsI->eiscr2rng();		// raise score for each intron
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

//	setup parameters. called once in a run

static	PwdB* SetUpPwd(Seq* sqs[])
{
	Seq*&	a = sqs[0];
	Seq*&	b = sqs[1];
	int	ap = a->isprotein();
	int	bp = b->isprotein();

	if (!ap && bp) {
	    prompt("Has exchanged the input order of %s and %s !\n",
		b->spath, a->spath);
	    std::swap(a, b); std::swap(ap, bp);
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
	if (gsquery) std::swap(sqs[0], sqs[1]);
	if (QRYvsDB == GvsA || QRYvsDB == AvsG) sqs[1]->inex.intr = algmode.lsg;
	makeWlprms(prePwd((const Seq**) sqs));
	PwdB* pwd = new PwdB((const Seq**) sqs);
	if (gsquery) std::swap(sqs[0], sqs[1]);
	return (pwd);
}

static int match_2(Seq* sqs[], PwdB* pwd, ThQueue* q)
{
	if (gsquery) std::swap(sqs[0], sqs[1]);
	Seq*&	a = sqs[0];
	Seq*&	b = sqs[1];

	int	dir = a->isprotein() + 2 * b->isprotein();
	if (dir != pwd->DvsP)
	    fatal("Don't mix different combinations of seq. types %s %d %d\n",
		(*a->sname)[0], dir, pwd->DvsP);
	if (algmode.lcl & 16) {	// local
	    a->exg_seq(1, 1);
	    b->exg_seq(1, 1);
	} else {		// global
	    a->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	    b->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	}
	if (QRYvsDB == GvsA || QRYvsDB == AvsG) b->inex.intr = algmode.lsg;
	INT	ori = a->inex.ori;
	int	nparalog = 0;
	if (pwd->DvsP == 3) nparalog = 1;
	else {
	    if (algmode.mns == 0 || algmode.mns == 3 || ori == 3) {
		if (pwd->DvsP == 1) b->nuc2tron();
		nparalog = geneorient(sqs, pwd);
	    }
	    if (algmode.mns == 1 && ori == 3) ori = b->inex.sens? 2: 1;
	    if (algmode.mns == 2 && ori == 3) ori = b->inex.sens? 1: 2;
	    if (algmode.lsg == 0 && ori == 3) ori = 1;
	    if (b->inex.intr && dbs_dt[0]) extend_gene_rng(sqs, pwd, dbs_dt[0]);
	}
	Gsinfo*	GsI = new Gsinfo[nparalog + 1];
	Gsinfo*	gsinf = GsI;
	INT	n_out = 0;
	Seq**	gener = sqs + 3;
	for (int n = 0; n < nparalog; ++n, ++gsinf) {
	    if (n) std::swap(b, gener[n - 1]);
	    dir = spalign2(sqs, pwd, gsinf, ori);
	    if (n) std::swap(b, gener[n - 1]);
	    if (dir != ERROR && (!algmode.thr || GsI->fstat.val >= pwd->Vthr)) {
		++n_out;
		if (algmode.nsa == BED_FORM)
		    GsI->rscr = selfAlnScr(a, pwd->simmtx);
	    } else  GsI->scr = NEVSEL;
	}
	delete b->exin; b->exin = 0;	// suppress Boundary output

//	insert sort GsI as nearly sorted
	int*	odr = 0;
	if (nparalog > 1) {
	    odr = new int[nparalog];
	    for (int n = 0; n < nparalog; ++n) odr[n] = n;
	    for (int n = 1; n < nparalog; ++n) {
		int	l = odr[n];
		VTYPE	v = GsI[l].fstat.val;
		int	m = n;
		while (--m >= 0 && v > GsI[odr[m]].fstat.val)
		    odr[m + 1] = odr[m];
		odr[m + 1] = l;
	    } 
	}
	if (OutPrm.MaxOut < n_out) n_out = OutPrm.MaxOut;
	for (INT n = 0; n < n_out; ++n) {
	    int	k = odr? odr[n]: 0;
	    if (k) std::swap(b, gener[k - 1]);
	    Gsinfo*	gsinf = GsI + k;
	    if (bool(gsinf->skl->m & A_RevCom) ^ bool(a->inex.sens))
		a->comrev();
	    if (gsquery) std::swap(sqs[0], sqs[1]);
	    if (algmode.nsa == BED_FORM)
		gsinf->rscr = selfAlnScr(a, pwd->simmtx);
#if M_THREAD
	    if (q) {
		pthread_mutex_lock(&q->mutex);
		outputs.alnoutput(sqs, gsinf);
		pthread_mutex_unlock(&q->mutex);
	    } else
#endif
		outputs.alnoutput(sqs, gsinf);
	    if (gsquery) std::swap(sqs[0], sqs[1]);
	    if (k) std::swap(b, gener[k - 1]);
	}
	if (dir == 1) {
	    a->comrev();
	    antiseq(sqs + 1);
	}
	if (gsquery) std::swap(sqs[0], sqs[1]);
	cleanseq(sqs + 2, nparalog + 2, 1);
	delete[] GsI; delete[] odr;
	return (dir);
}

static int blkaln(Seq* sqs[], SrchBlk* bks, RANGE* rng, ThQueue* q)
{
	int	nparalog = 0;
	RANGE	grng = {0, 0};
	RANGE	wrng;
	RANGE	frng = {INT_MAX, 0};
	Seq*&	a = sqs[0];
	Seq*&	b = sqs[1];
	int	b_intr = b->inex.intr;

	if (gsquery) {
	    a->saverange(&grng);
	    int	n = (grng.right - grng.left) / 10;
	    frng.left = grng.left + n;
	    frng.right = grng.right - n;
	}
	bks->setseqs(sqs);

	if (algmode.lcl & 16) a->exg_seq(1, 1);	// SWG local alignment
	else	a->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	if (QRYvsDB == AvsA || QRYvsDB == CvsC || QRYvsDB == GvsC)
	    nparalog = bks->finds(sqs);
	else if (QRYvsDB == GvsA)
	    nparalog = bks->findh(sqs);
	else
	    nparalog = bks->findblock(sqs);
	if (nparalog == ERROR || algmode.nsa == MAP1_FORM || algmode.nsa == MAP2_FORM)
	    return (0);		// no alignment
#if USE_FULL_RANGE
// skip the next line if only the specified range is used as the query
	if (rng) a->right = a->right == rng->right? a->len: rng->right;
#endif
	HalfGene	hfg(0, INT_MAX, a->CdsNo);
	Gsinfo*	GsI = new Gsinfo[nparalog + 1];
	Gsinfo*	gsinf = GsI;
	INT	n_out = 0;
	Seq**	gener = bks->get_gener();
const	int	basis = gener - sqs;
	for (int n = 0; n < nparalog; std::swap(b, gener[n++]), ++gsinf) {
	    std::swap(b, gener[n]);
	    if ((QRYvsDB == AvsA || QRYvsDB == CvsC) && OutPrm.supself) {
		int	sc = strcmp((*a->sname)[0], (*b->sname)[0]);  
		if (!sc || (OutPrm.supself > 1 && sc > 0)) {
		    gsinf->scr = NEVSEL;
		    continue;
		}
	    }
	    b->inex.intr = b_intr;
	    if (gsquery) {
		std::swap(sqs[0], sqs[1]);
		std::swap(a->CdsNo, b->CdsNo);
		if (a->jxt) {
		    delete[] b->jxt;
		    b->jxt = a->jxt;
		    a->jxt = 0;
		    a->saverange(&wrng);
		    b->restrange(&wrng);
		    a->left = 0;
		    a->right = a->len;
		    a->inex.ori = QRYvsDB == GvsA? 1: 3;
		}
		if (algmode.lcl & 16) a->exg_seq(1, 1);
		else	a->exg_seq(algmode.lcl & 4, algmode.lcl & 8);
	    }				// incompatible gene orientation
	    if (algmode.lcl & 16) b->exg_seq(1, 1);
	    else	b->exg_seq(algmode.lcl & 1, algmode.lcl & 2);
	    int	dir = spalign2(sqs, bks->pwd, gsinf, a->inex.ori);
	    delete b->exin; b->exin = 0;	// suppress Boundary output
	    if (dir < 0 || !gsinf->skl ||
		(algmode.thr && gsinf->scr < bks->pwd->Vthr)) {
		gsinf->scr = NEVSEL;
		continue;
	    }
	    if (gsquery) {
		int	suppress = 0;
		int	r = gsinf->skl[0].n - 1;
		int	sa = gsinf->skl[1].m;
		int	sb = gsinf->skl[1].n;
		int	eb = gsinf->skl[r].n;
		int	db = sb - b->left;
		if (QRYvsDB == GvsA) sa *= 3;
		if ((IntronPrm.tlmt < sa) && (db < IntronPrm.rlmt) &&
			(b->left > IntronPrm.rlmt)) {
		    if (eb > hfg.left) hfg.left = eb;
		    if (algmode.lcl < 16) suppress |= 1;
		}
		sa = a->right - gsinf->skl[r].m;
		db = b->right - eb;
		if (QRYvsDB == GvsA) sa *= 3;
		if ((IntronPrm.tlmt < sa) && (db < IntronPrm.rlmt) &&
			(b->len - b->right > IntronPrm.rlmt)) {
		    if (sb < hfg.right) hfg.right = sb;
		    if (algmode.lcl < 16) suppress |= 2;
		}
		std::swap(sqs[0], sqs[1]);
		a->restrange(&grng);
		std::swap(b->CdsNo, a->CdsNo);
		if (suppress) {
		    gsinf->scr = NEVSEL;
		    continue;
		}
	    } else {
		hfg.right = dir + 1;
	    }
	    if (algmode.mlt == 1 && gsinf->eijnc) {
		EISCR*	fst = gsinf->eijnc->begin();
		a->left = fst[0].rleft;
		a->right = fst[gsinf->noeij - 1].rright;
	    }
	    ++n_out;
	}	// end of nparalog loop
//	insert sort GsI as nearly sorted
	int*	odr = 0;
	if (nparalog > 1) {
	    odr = new int[nparalog];
	    for (int n = 0; n < nparalog; ++n) odr[n] = n;
	    for (int n = 1; n < nparalog; ++n) {
		int	l = odr[n];
		VTYPE	v = GsI[l].fstat.val;
		int	m = n;
		while (--m >= 0 && v > GsI[odr[m]].fstat.val)
		    odr[m + 1] = odr[m];
		odr[m + 1] = l;
	    } 
	}
	if (OutPrm.MaxOut < n_out) n_out = OutPrm.MaxOut;
	for (INT n = 0; n < n_out; ++n) {
	    int	k = odr? odr[n]: 0;
	    std::swap(b, gener[k]);
	    Gsinfo*	gsinf = GsI + k;
	    if (bool(gsinf->skl->m & A_RevCom) ^ bool(a->inex.sens))
		a->comrev();
	    if (gsquery) std::swap(sqs[0], sqs[1]);
	    else {
		if (b->left < frng.left) frng.left = b->left;
		if (b->right > frng.right) frng.right = b->right;
	    }
	    if (algmode.nsa == BED_FORM)
		gsinf->rscr = selfAlnScr(a, bks->pwd->simmtx);
#if M_THREAD
	    if (q) {
		pthread_mutex_lock(&q->mutex);
		outputs.alnoutput(sqs, gsinf);
		pthread_mutex_unlock(&q->mutex);
	    } else
#endif
		outputs.alnoutput(sqs, gsinf);
	    if (gsquery) std::swap(sqs[0], sqs[1]);
	    std::swap(b, gener[k]);
	}
#if M_THREAD
	if (q && q->mfd) {
	    pthread_mutex_lock(&q->mutex);
	    hfg.sname = strrealloc(0, a->sqname());
	    if (hfg.left > 0) hfg.left += IntronPrm.minl;
	    else	hfg.left = frng.left;
	    if (hfg.right < INT_MAX) hfg.right -= IntronPrm.minl;
	    else	hfg.right = frng.right;
	    hfg.left = a->SiteNo(hfg.left);
	    hfg.right = a->SiteNo(hfg.right);
	    q->mfd->write(&hfg);
	    pthread_mutex_unlock(&q->mutex);
	}
#endif
	cleanseq(sqs + 2, basis + nparalog - 1, 1);
	delete[] GsI; delete[] odr;
	return (hfg.right);
}

static	int	MinSegLen = INT_MAX;

static SrchBlk* getblkinf(Seq* sqs[], const char* dbs, MakeBlk* mb)
{
	Seq*&	a = sqs[0];	// query
	Seq*&	b = sqs[1];	// database
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
	    std::swap(sqs[0], sqs[1]);
	    std::swap(ap, bp);
	}
	makeWlprms(prePwd((const Seq**) sqs));
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
	    if (!ReadBlock) {
		char*	dot = strrchr(str, '.');
		if (dot) *dot = 0;
		fatal("%s not found. Specify database!\n", str);
	    }
	}
	if ((QRYvsDB == CvsC || QRYvsDB == FvsG) && algmode.mns != 2) algmode.mns = 1;
	SrchBlk* bks = dbs? new SrchBlk(sqs, ReadBlock, genomedb):
		new SrchBlk(sqs, mb, b->inex.intr);
	if (MinQueryLen == 0) MinQueryLen = bks->MinQuery();
	if (algmode.mlt == 1) MinSegLen = MinQueryLen;
	if (ap) q_mns &= ~2;
	if (gsquery){
	    std::swap(sqs[0], sqs[1]);
	    std::swap(ap, bp);
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
	if (a->right - a->left < MinQueryLen) {
	    prompt("%s (%d - %d) %d is too short!\n", 
		a->sqname(), a->left, a->right, a->right - a->left);
	    return (0);
	}
#if TESTRAN
	a->rndseq();
#endif
	if (gsquery) return (blkaln(sqs, bks, 0, q));

	RANGE	orgrng;
	RANGE	covrng;

	a->saverange(&orgrng);
	int	rbdry = blkaln(sqs, bks, &orgrng, q);

// dispersed locations?

	a->saverange(&covrng);
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
	a->restrange(&orgrng);
	return (0);
}

static void genomicseq(Seq** sqs, PwdB* pwd, int ori)
{
	sqs[0]->inex.intr = algmode.lsg;
	if (ori & 1) {
	    if (pwd->DvsP == 1) sqs[0]->nuc2tron();
	    if (!sqs[0]->exin)
		sqs[0]->exin = new Exinon(sqs[0], pwd, ori == 3);
	}
	if (ori & 2 && pwd->DvsP != 3) {
	    sqs[0]->comrev(sqs + 1);
	    sqs[1]->setanti(sqs);
	    if (pwd->DvsP == 1) sqs[1]->nuc2tron();
	    if (!sqs[1]->exin)
		sqs[1]->exin = new Exinon(sqs[1], pwd, ori == 3);
	}
}

static void spaln_job(Seq* sqs[], void* prm, ThQueue* q)
{
	Seq*&	a = sqs[0];
	Seq*&	b = sqs[1];
	PwdB*	pwd = (PwdB*) prm;
	a->inex.intr = (gsquery)? algmode.lsg: 0;
	if (a->isdrna() && !gsquery && q_mns)	 	// cDNA
	    a->inex.ori = polyA.rmpolyA(a, a->sigII? 1: q_mns);
	else	a->inex.ori = q_mns? q_mns: 3;
	if (pwd->DvsP != 3)	b->inex.sigs = b->inex.sigt = 1;
	if (algmode.crs) Wlp	wlp(sqs);		// mask low ic region
	if (algmode.blk) (void) quick4(sqs, (SrchBlk*) prm, q);
	else	match_2(sqs, pwd, q);
}

static void seg_job(Seq** sqs, SeqServer* svr, SrchBlk* sbk)
{
	if (!sbk) fatal("No block inf !\n");
	if (QRYvsDB == GvsA) genomicseq(sqs, sbk->pwd, 1);
	quick4(sqs, sbk);
	if (QRYvsDB == GvsC) return;
	antiseq(sqs);
	quick4(sqs, sbk);
	antiseq(sqs);
}

static void all_in_func(Seq** sqs, SeqServer* svr, void* prm)
{
	int	nf = svr->input_form == IM_PARA? 1: 0;
	SrchBlk*	sbk = 0;
	if (algmode.blk) {
	    sbk = (SrchBlk*) prm;
	    if (!sbk->dbf) sbk->dbf = svr->target_dbf;
	}
	do {
// perform the job
	    if (gsquery) seg_job(sqs, svr, sbk);
	    else	spaln_job(sqs, prm);
// read pair of sequences
	    if (svr->input_ns == 2) {
		switch (svr->nextseq(sqs[1], nf)) {
		    case IS_END: case IS_ERR: return;
		    default: break;
		}
	    }
// read query sequence
	    switch (svr->nextseq(sqs[0], 0)) {
		case IS_END:
		    closeGeneRecord();
		    return;
		case IS_ERR: continue;
		default: break;
	    }
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
		if (sqs[n]->sid > 0) std::swap(sque[wp], sqs[n]);
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
	    if (sque[rp]->many)	{std::swap(sqs[n], sque[rp]);
	    fprintf(stderr, "d%d: %s %d %d %d\n", n, sqs[n]->sqname(),
	    sqs[n]->sid, sqs[n]->len, sqs[n]->many);}
#else
	    if (sque[rp]->many)	std::swap(sqs[n], sque[rp]);
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
		else q->enqueue(fsd, targ->svr->input_ns);
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
	closeGeneRecord();
//	reportseq(targ[0].seqs, no_seqs * thread_num);
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
	setpam(InsPam, 0);	// intra-species PAM value
	setpam(CrsPam, 1);	// cross-species PAM value
	setpam(WlpPam, WlnPamNo);	// pam for HSP search
	setorf(DefOrfLen, 2);	// orf length and AG< .. >GT ends
	OutPrm.MaxOut = 1;
	OutPrm.SkipLongGap = 1;	// suppress display of long gaps
	OutPrm.fastanno = 1;	// add annotation in fasta output
	polyA.setthr(def_polya_thr);

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
	    if (MaxVmfSpace == DefMaxVMF) MaxVmfSpace = MaxWordNoSpace;
	    Seq*	sqs[2];
	    initseq(sqs, 2);
	    int	molc = TgtMolc;
	    if (!molc) molc = infermolc(argv[0]);
	    if (molc == PROTEIN) setSeqCode(sqs[0], PROTEIN);
	    else	setSeqCode(sqs[0], DNA);
	    if (molc == DNA)	setSeqCode(sqs[1], DNA);
	    else setSeqCode(sqs[1], PROTEIN);
	    makeWlprms(prePwd((const Seq**) sqs));
	    clearseq(sqs, 2);
	    MakeBlk*	mb = 0;
	    if (TgtMolc == PROTEIN && alprm2.spb > 0.) {
		SeqServer	svr(argc, argv, IM_SNGL, 0, TgtMolc);
		mb = makeblock(&svr);
	    } else
		mb = makeblock(argc, argv, molc, TgtMolc);
	    if (mb) {
		mb->WriteBlkInfo();
		mb->delete_dbf();
		delete mb;
	    }
	    resetSimmtxes(true);
	    eraWlprms();
	    return (0);
	}
	if (!catalog && argc <= 0) usage(insuf);
#if M_THREAD
	cpu_num = sysconf(_SC_NPROCESSORS_CONF);
	if (thread_num < 0) thread_num = cpu_num;
	if (max_queue_num == 0) max_queue_num = int(FACT_QUEUE * thread_num);
#endif	// M_THREAD
	if (QryMolc == TRON) QryMolc = PROTEIN;
	if (!algmode.mlt) OutPrm.MaxOut = 1;
	if (OutPrm.MaxOut > OutPrm.MaxOut2) OutPrm.MaxOut2 = OutPrm.MaxOut;
	setup_output(CDS_FORM, 0, false);	// output format
	int	nseqs = no_seqs += OutPrm.MaxOut2 + 1;
	Seq**	seqs = new Seq*[nseqs];
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
	initseq(seqs, nseqs);		// 0: query, 1: genomic, 2: reverse
	MakeBlk*	mb = 0;
const	char*	dbs = genomedb? genomedb: (aadbs? aadbs: cdnadb);
	if (algmode.blk && !dbs) {	// make db from 1st file
	    if (--argc <= 0) {
		clearseq(seqs, nseqs);
		delete[] seqs;
		usage(insuf);
	    }
	    if (!TgtMolc) {
		TgtMolc = infermolc(argv[0]);
	    	QryMolc = infermolc(argv[1]);
		if (TgtMolc != PROTEIN && QryMolc == PROTEIN)
		    TgtMolc = TRON;
	    }
	    SeqServer	svr(1, argv++, IM_SNGL, 0, TgtMolc);
	    mb = makeblock(&svr);
	    dbs_dt[0] = mb->get_dbf();
	}
	SeqServer	svr(argc, argv, input_form, catalog, QryMolc, TgtMolc);
	if (algmode.blk) {
	    setdefmolc(QryMolc);
	    if (svr.nextseq(a, 0) > IS_ERR) {
		messg = "Can't open query !\n";
		goto postproc;
	    }
	    if (outputs.setup(a->spath)) {
		SrchBlk*	bprm = getblkinf(seqs, dbs, mb);
		delete mb;
		no_seqs = OutPrm.MaxOut2 + bprm->NoWorkSeq + 2;	// query + last
		if (a->inex.intr || b->inex.intr) makeStdSig53();
		set_max_extend_gene_rng(def_max_extend_gene_rng);
		if (bprm->pwd->DvsP == 3) OutPrm.SkipLongGap = 0;
		if (algmode.nsa == SAM_FORM) put_genome_entries();
#if M_THREAD
		if (thread_num) MasterWorker(seqs, &svr, (void*) bprm); else
#endif
		all_in_func(seqs, &svr, (void*) bprm);
		delete bprm;
	    }
	} else {
	    no_seqs = OutPrm.MaxOut2 + 2;
	    if (svr.nextseq(b, 1) == IS_END) {
		messg = "Can't open genomic sequence !\n";
		goto postproc;
	    }
	    if (svr.nextseq(a, 0) != IS_OK) {
		messg = "Can't open query !\n";
		goto postproc;
	    }
	    PwdB*	pwd = SetUpPwd(seqs);
	    if (pwd->DvsP == 3) OutPrm.SkipLongGap = 0;
	    if (a->inex.intr || b->inex.intr) 
		makeStdSig53();
	    set_max_extend_gene_rng(0);
	    if (outputs.setup(a->spath)) {
#if M_THREAD
		if (thread_num) MasterWorker(seqs, &svr, (void*) pwd); else
#endif
		all_in_func(seqs, &svr, (void*) pwd);
	    }
	    delete pwd;
	}

	eraWlprms();
	eraStrPhrases();
postproc:
	EraStdSig53();
	clearseq(seqs, nseqs);
	delete[] seqs;
	EraDbsDt();
	close_output();
	resetSimmtxes(true);
	if (messg) usage(messg);
	return (0);
}
