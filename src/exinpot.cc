/*****************************************************************************
*
*	exinpot [-lN] [-rN] [-TS] input_seqs
*	exinpot [-lN] [-rN] [-TS] -e xxx.eij [-d gnm_g]
*
*	Calculate intron potential based on k-mer frequences
*
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

#include "seq.h"
#include "eijunc.h"

#define	OPTCHAR	'-'

static	int	min_orf = 120;
static	int	text = 1;
static	bool	gzip = false;
static	const	char*	fname = 0;
static	const	char*	ewdfq = 0;
static	const	char*	gwdfq = 0;
static	const	char*	iwdfq = 0;
static	const	char*	iname = 0;
static	const	char*	oname = 0;
static	const	char*	uname = 0;
static	const	char*	eij = 0;
static	const	char*	gnm = 0;
static	const	char*	cus = 0;
static	const	int	ng = static_cast<int>(Iefp::NG);

static	const	char* hfmt = "%7d\t%7.3f %7.3f %7.1f %9.3f %7.3f %d %d\n";

static void usage()
{
	fputs("Description:\n", stdout);
	fputs("\tStationary Markov model (MM) for calculating intron potential (.ipt)\n",
		 stdout);
	fputs("\texon potential (.ept) or phase-sensitive coding potential (.cdp)\n",
		 stdout);
	fputs("\tcodon usage table (.cus) is also generated together with cdp\n",
		 stdout);
	fputs("Usage:\n", stdout);
	fputs("\texinpot [-mN] -d gnm -r X.wdfq [-g] [-o|-b X.[ipt|ifp]] -e X.eij or\n", stdout);
	fputs("\texinpot [-mN] -r X.wdfq [-g] [-o|-b X.[cdp|cfp]] [-JMinOrf] -c CDS or\n", stdout);
	fputs("\texinpot -O[1|2|4] -f X.ipt.dgz [-e X.eij] or\n", stdout);
	fputs("\texinpot -O[1|2|4] -f X.cdp.dgz [CDS] or\n", stdout);
	fputs("\texinpot -O[1|2|4] -f X.ept.dgz [cDNA] or\n", stdout);
	fputs("\texinpot -u X.cud.dat [-o X.cud]\n", stdout);
	fputs("Examples:\n", stdout);
	fputs("\texinpot -m4 -d X_g -r X.wdfq -g -b X.ipt -e X.eij\n", stdout);
	fputs("\texinpot -m5 -r X.wdfq -g -b X.cdp -J300 -c cDNA.fna\n", stdout);
	fputs("Options:\n", stdout);
	fputs("\t-b S: binary output file\n", stdout);
	fputs("\t-c: CDS input\n", stdout);
	fputs("\t-d [S]: genome db\n", stdout);
	fputs("\t-e S: X.eij\n", stdout);
	fputs("\t-f S: existing ipt, cdp or ept [.dat|.dgz]\n", stdout);
	fputs("\t-g S: gzipped output (.gz or .dgz)\n", stdout);
	fputs("\t-i S: X.iwdfq\n", stdout);
	fputs("\t-m N: m-th order MM (4 for ipt and ept, 5 for cdp)\n", stdout);
	fputs("\t-o S: readable output file\n", stdout);
	fputs("\t-r S: X.wdfq (background kmer frequency\n", stdout);
	fputs("\t-u S: X.cus.dat (codon usage)\n", stdout);
	fputs("\t-C N: NCBI genetic code (0: universal)\n", stdout);
	fputs("\t-E S: X.ewdfq\n", stdout);
	fputs("\t-J N: minimum ORF (120: -c CDS is given)\n", stdout);
	fputs("\t-O N: output mode (4)\n", stdout);
	fputs("\t\tN & 1: summary\n", stdout);
	fputs("\t\tN & 2: individual seqs\n", stdout);
	fputs("\t\tN & 4: write X.[ipt|ept|cdp]\n", stdout);
	fputs("Terms in headline:\n", stdout);
	fputs("\t1\tidentifier\n", stdout);
	fputs("\t2\t# of phases (1 or 3)\n", stdout);
	fputs("\t3\t# of rows\n", stdout);
	fputs("\t4\torder of Markov model\n", stdout);
	fputs("\t5\t# of sample objects (intron or CDS) used for calculation\n", stdout);
	fputs("\t6\tmean score (KL divergence) per 1000 bp\n", stdout);
	fputs("\t7\tSD of per-1000-bp score\n", stdout);
	fputs("\t8\tmean length of object\n", stdout);
	fputs("\t9\tmean score per object\n", stdout);
	fputs("\t10\tSD of per-object score\n", stdout);
	fputs("\t11\t5' side margin of object\n", stdout);
	fputs("\t12\t3' side margin of object\n", stdout);
	exit(1);
}

class Ipt : public ExinPot {
public:
	bool	nrml = false;
	CodonUse*	codons = 0;
	void	reset() {
	    avm = sdm = avpot = sdpot = 0; avlen = 0; nsupport = 0;
	    nrml = false;
	}
	void	make_ipt(Seq* sd, const int n);
	void	normalize();
	void	recover();
	void	finish(const float* gfq);
template <typename file_t>
	void	write_text(file_t fd, const float* gfq);
	Ipt(const int x, const int m, const int p, const int it) 
		: ExinPot(x, m, p, it) {
	    if (it) {
		avm = sdm = 0;
		nrml = true;
	    }
	    if (p == 3 && !fname && algmode.nsa & 4)
		codons = new CodonUse(uname? uname: iname, false);
	}
	~Ipt() {delete codons;}
};

void Ipt::make_ipt(Seq* sd, const int n)
{
	int	l = sd->right - sd->left - lm - rm;
	if (l <= 0) return;
	float	v = 0;
	float*	ipt = calcScr(sd, &v);
	float	m = v / l;
	++nsupport;
	avpot += v; sdpot += v * v; avm += m; sdm += m * m; avlen += l;
	if (algmode.nsa & 2)
	    fprintf(out_fd, "%5d\t%12.4e\t%12.4e %7d\n", n, v, m, l);
	delete[] ipt;
}

void Ipt::normalize()
{
	if (!nsupport) return;	// empty
	nrml = true;
	avlen /= nsupport;
	avpot /= nsupport;
	sdpot = sqrt(double(sdpot - nsupport * avpot * avpot) / (nsupport - 1));
	avm /= nsupport;
	sdm = sqrt(double(sdm - nsupport * avm * avm) / (nsupport - 1));
	avm *= 1000.;	 	// per 1 kbp
	sdm *= 1000.;
}

void Ipt::recover()
{
	nrml = true;
const	float	approximate_avm = 1000. * avpot / avlen;
	if (avm > 20. * approximate_avm) {
	    avm = approximate_avm;	// older version
	    sdm = 0.;			// not available
	}
}

template <typename file_t>
void Ipt::write_text(file_t fd, const float* gfq) {
	char	str[LINE_MAX];
	if (!algmode.nsa || algmode.nsa & 4) {	// readable output
	    snprintf(str, LINE_MAX, "%s %d %d %d ", 
		iefp_tid[exin], nphase, size(), morder);
	    if (ispot() && nrml) {
		char*	ps = str + strlen(str);
		snprintf(ps, LINE_MAX - (ps - str), hfmt, 
		nsupport, avm, sdm, avlen, avpot, sdpot, lm, rm);
	    } else	strcat(str, "\n");
	    fputs(str, fd);
const	    float*	pot = begin();
const	    float*	dend = end();
	    if (!pot) {
		if (!gfq) return;	// nothing to do
		pot = gfq;
		dend = pot + ndata;
	    }
	    char*	ps = str;
	    for (int p = 0, n = 0; pot < dend; ++pot) {
		if (nphase == 3 && p == 0) {
		    int	c = n++;
		    for (int m = ndata; m /= 4; ) {
			*ps++ = Nucl[c / m];
			c %= m;
		    }
		    *ps++ = '\t';
		}
		if (algmode.nsa) {
		    snprintf(ps, LINE_MAX - (ps - str), "%15.7e", *pot);
		    ps += 15;
		} else {
		    snprintf(ps, LINE_MAX - (ps - str), "%7d", (int) *pot);
		    ps += 7;
		}
		if (nphase == 3) {
		    if ((p = next_p[p]) == 0) {
			*ps++ = '\n';
			*ps = '\0';
			fputs(str, fd);
			ps = str;
		    } else	*ps++ = '\t';
		} else {
		    *ps++ = '\n';
		    *ps = '\0';
		    fputs(str, fd);
		    ps = str;
		}
	    }
	}
}

void Ipt::finish(const float* gfq)
{
	if (ispot() && nsupport == 0)
	    fatal("No sequence data !\n");
	if (codons) {
	    codons->normalize();
	    codons->to_file(uname, text);
	}
	if (iname && algmode.nsa & 1) {
	    char	str[MAXL];
	    char*	ps = str;
const	    char*	sls = strrchr(iname, '/');
	    sls = sls? sls + 1: iname;
	    strcpy(str, sls);
	    if (str[8] == '_') str[8] = '\0';
	    else {
		char*	dot = strchr(str, '.');
		if (dot) *dot = '\0';
	    }
	    ps = str + strlen(str);
	    *ps ++ = ':';
	    *ps = '\0';
const	    char*	rname = fname? fname: iname;
	    sls = strrchr(rname, '/');
	    sls = sls? sls + 1: rname;
	    strcpy(ps, sls);
	    if (ps[8] == '_') ps[8] = '\0';
	    else {
		char*	dot = strchr(ps, '.');
		if (dot) *dot = '\0';
	    }
	    printf("%s\t", str);
	    if (nrml)
		printf(hfmt, nsupport, avm, sdm, avlen, avpot, sdpot, lm, rm);
	    else
		printf("%7d\t%7.2f %7.2f %7.2f\n", nsupport, avlen, avpot, avm);
	}
	if (text) {
	    WriteFile	fp(oname, 1, gzip);
	    if (fp.gzfd)	write_text(fp.gzfd, gfq);
	    else if (fp.fd)	write_text(fp.fd, gfq);
	    else	fatal(no_file, uname);
	} else {
	    writeBinary(oname, gzip);
	}
}

int main(int argc, const char** argv)
{
	algmode.nsa = 4;	// default output mode
	int	mo = 4;		// default MM order
	bool	cds = false;
	int	exin = ng;
	int	file_type = 0;
	while (--argc > 0 && **++argv == OPTCHAR) {
const	    char*	opt = argv[0] + 1;
	    int	c = *opt;

	    if (!c) break;
const	    char*	val = argv[0] + 2;
	    switch (c) {
		case 'b':
		    if ((val = getarg(argc, argv)))
			oname = uname = val;
		    text = 0;
		    break;
		case 'c': cds = true; break;
		case 'd':
		    if ((val = getarg(argc, argv)))
			gnm = val;
		    break;
		case 'e':
		    if ((val = getarg(argc, argv)))
			eij = val;
		    break;
		case 'f':
		    if ((val = getarg(argc, argv)))
		        fname = val;
		    break;
		case 'g':	gzip = true; break;
		case 'h':	usage();
		case 'i':
		    if ((val = getarg(argc, argv)))
		        iwdfq = val;
		    break;
		case 'm':
		    if ((val = getarg(argc, argv, true)))
		    	mo = atoi(val);
		    break;
		case 'o':
		    if ((val = getarg(argc, argv)))
		        oname = uname = val;
		    text = 2;
		    break;
		case 'p':
		    if (!strcmp(opt, "pq")) setprompt(0, 0);
		    break;
		case 'r':
		    if ((val = getarg(argc, argv)))
		        gwdfq = val;
		    break;
		case 'u':
		    if ((val = getarg(argc, argv)))
		        cus = val;
		    break;
		case 'C': 
		    if ((val = getarg(argc, argv, true)))
			initcodon(atoi(val));
		    break;
		case 'E':
		    if ((val = getarg(argc, argv)))
		        ewdfq = val;
		    break;
		case 'J': 
		    if ((val = getarg(argc, argv, true)))
			min_orf = atoi(val);
		    break;
		case 'O':
		    if ((val = getarg(argc, argv, true)))
		        algmode.nsa = atoi(val);
		    break;
		default: break;
	    }
	}

	if (fname) {
	    exin = fname2exin(fname, file_type);
	} else if (gwdfq) {
	    if (cds && argc)	exin = static_cast<int>(Iefp::CP); else
	    if (iwdfq || eij)	exin = static_cast<int>(Iefp::IP); else
	    if (ewdfq || argc)	exin = static_cast<int>(Iefp::EP); else
				exin = static_cast<int>(Iefp::GF);
	} else {
	    if (cds && argc)	exin = static_cast<int>(Iefp::CF); else
	    if (iwdfq || eij)	exin = static_cast<int>(Iefp::IF); else
	    if (ewdfq || argc)	exin = static_cast<int>(Iefp::EF);
	}

	Iefp	iefp = static_cast<Iefp>(exin);
	if (!gwdfq && !fname && (iefp == Iefp::GF ||
	    iefp == Iefp::IP || iefp == Iefp::IB ||
	    iefp == Iefp::EP || iefp == Iefp::EB ||
	    iefp == Iefp::CP || iefp == Iefp::CB)) usage();
	if (!(iefp == Iefp::CF || iefp == Iefp::CP || iefp == Iefp::CB))
	    min_orf = 0;
	setorf(min_orf);

	char	str[MAXL];
	bool	get_kmers = iwdfq || ewdfq || gwdfq;
	iname = eij? eij: (argc? *argv: 0);
	if (oname) {		// test compatibility
	    int	ft = 1;
const	    int	exn = fname2exin(oname, ft);
	    if (exn == 0)	// undefined file type
		oname = add_ext(oname, iefp_ext[exin], str);
	    else if ((fname && exn > exin) || (exin - exn) / 3) usage();
	    else exin = exn;
	    if (text) oname = add_ext(oname, text_ext, str);
	}
	int	nphs = ((exin / 3) == 2)? 3: 1;
	Ipt	ipt(exin, mo, nphs, fname? 0: file_type);
	float*	gfq = 0;
	if (gwdfq) {
	    gfq = ipt.getKmers(gwdfq, false);
	    ipt.reform(gfq);
	}
	EiJuncSeq*	eijseq = eij? new EiJuncSeq(INTRON, eij, gnm): 0;
	if (ewdfq) iwdfq = ewdfq;
	if (fname && ipt.from_file(fname)) {
	    ipt.recover();
	} else if (cus) {
	    CodonUse	codons(cus, true);
	    codons.to_file(uname, text);
	    exit(0);
	} else if (get_kmers) {
	    if (iwdfq)	ipt.getKmers(iwdfq); else
	    if (eij)	ipt.getKmers(eijseq, ipt.codons); else
	    if (argc)	ipt.getKmers(argc, argv, ipt.codons);
	    if (algmode.nsa)	ipt.reform();
	    if (gfq && !ipt.makeExinPot(gfq)) usage();
	} else {
	    usage();
	}

	if (!ipt.ispot()) {
// nothing to do
	} else if (eij) {
	    eijseq->reset();
	    ipt.reset();
	    int	n = 0;
	    do {
		Seq*	sd = eijseq->nextseq();
		if (!sd || sd->inex.ambs) continue;
		ipt.make_ipt(sd, n++);
	    } while (eijseq->goahead());
	    ipt.normalize();
	} else if (argc > 0) {
	    Seq	seq(1);
	    SeqServer	sqsvr(argc, argv, IM_SNGL);
	    ipt.reset();
	    int	n = 0;
	    while (sqsvr.nextseq(&seq, 0) != IS_END) {
		if (min_orf) {
		    ORF*	orfs = seq.getorf();
		    if (!orfs) continue;
		    seq.left = orfs->pos;
		    seq.right = orfs->pos + orfs->len;
		    delete[] orfs;
		}
		ipt.make_ipt(&seq, n++);
	    }
	    ipt.normalize();
	}
	ipt.finish(gfq);
	delete eijseq;
	delete[] gfq;
	fclose(out_fd);
	return (0);
}
