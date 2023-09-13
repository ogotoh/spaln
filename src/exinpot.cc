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

static	int	omode = 4;
static	int	min_orf = 0;
static	bool	binary = false;
static	const	char*	fname = 0;
static	const	char*	ewdfq = 0;
static	const	char*	gwdfq = 0;
static	const	char*	iwdfq = 0;
static	const	char*	oname = 0;
static	const	char*	eij = 0;
static	const	char*	gnm = 0;

static	const	char* hfmt = "%7.5f %7.5f %7.5f %d %d %d %7.2f %7.5f %7.5f %7.5f\n";
static	const	char* hfmt7 = "%7.5f %7.5f %7.5f %d %d %d %7.2f\n";

static void usage()
{
	fputs("Usage:\n", stdout);
	fputs("\texinpot -mN -g X.wdfq -b X.[ipt|ifp] [-d gnm] -e X.eij\n", stdout);
	fputs("\texinpot -mN -g X.wdfq -b X.[ept|efp] cdna.fna\n", stdout);
	fputs("\texinpot -m5 -c -g X.wdfq -b X.[cdp|cfp] cds.fna\n", stdout);
	fputs("\texinpot -m5 -c -g X.wdfq -b X.[cdp|cfp] -JMinOrf cdna.fna\n", stdout);
	fputs("\texinpot -f xxx.ipt -e X.eij\n", stdout);
	fputs("\texinpot -f xxx.cdp cds\n", stdout);
	fputs("\texinpot -f xxx.ept cdna\n", stdout);
	fputs("Options:\n", stdout);
	fputs("\t-b S: binary output file\n", stdout);
	fputs("\t-c [S]: xxx.ewdfq or CDS\n", stdout);
	fputs("\t-d [S]: genome db\n", stdout);
	fputs("\t-e S: xxx.eij\n", stdout);
	fputs("\t-f S: xxx.ipt\n", stdout);
	fputs("\t-g S: xxx.wdfq\n", stdout);
	fputs("\t-i S: xxx.iwdfq\n", stdout);
	fputs("\t-m N: m-th order MM\n", stdout);
	fputs("\t-o S: readable output file\n", stdout);
	fputs("\t-C N: NCBI genetic code (0: universal)\n", stdout);
	fputs("\t-J N: minimum ORF (0: CDS is given)\n", stdout);
	fputs("\t-O N: output mode (4)\n", stdout);
	fputs("\t\tN & 1: summary\n", stdout);
	fputs("\t\tN & 2: individual seqs\n", stdout);
	fputs("\t\tN & 4: write xxx.ipt\n", stdout);
	fputs("\t-T S: table directory\n", stdout);
	exit(1);
}

class Ipt : public ExinPot {
	float	avm = 0;
	float	minv = LONG_MAX, minm = LONG_MAX;
	float	maxv = LONG_MIN, maxm = LONG_MIN;
	bool	nrml = false;
public:
	void	make_ipt(Seq* sd, const int n);
	void	normalize();
	void	finish(const float* gfq);
	Ipt(int x, int m, int p) : ExinPot(x, m, p) {}
	~Ipt() {}
};

void Ipt::make_ipt(Seq* sd, const int n)
{
	int	l = sd->right - sd->left - lm - rm;
	if (l <= 0) return;
	float	v = 0;
	float*	ipt = calcScr(sd, &v);
	float	m = v / l;
	++nsupport;
	avpot += v; avm += m; avlen += l;
	if (v < minv) minv = v;
	if (m < minm) minm = m;
	if (v > maxv) maxv = v;
	if (m > maxm) maxm = m;
	if (omode & 2)
	    fprintf(out_fd, "%5d\t%12.4e\t%12.4e %7d\n", n, v, m, l);
	delete[] ipt;
}

void Ipt::normalize()
{
	if (!nsupport) return;	// empty
	nrml = true;
	avlen /= nsupport;
	avpot  /= nsupport;
	minm *= avlen;
	maxm *= avlen;
	avm *= avlen / nsupport;
}

void Ipt::finish(const float* gfq)
{
	if (ispot() && nsupport == 0)
	    fatal("No sequence data !\n");
	if (binary) writeBinary(oname);
	else if (omode & 1 && nrml)
	    fprintf(out_fd, hfmt,
		minm, avm, maxm, nsupport, lm, rm, 
		avlen, minv, avpot, maxv);
	else if (omode & 1)
	    fprintf(out_fd, "%7d %7.2f %7.2f %7.2f\n",
		nsupport, avlen, avpot, ess);
	else if (!omode || omode & 4) {	// readable output
	    fprintf(out_fd, "%s %d %d %7.1f ", 
		iefp_tid[exin], nphase, size(), total / ndata / nphase);
	    if (ispot() && nrml) fprintf(out_fd, hfmt7, 
		minm, avm, maxm, nsupport, lm, rm, avlen);
	    else	fputc('\n', out_fd);
const	    float*	pot = begin();
const	    float*	dend = end();
	    if (!pot) {
		if (!gfq) return;	// nothing to do
		pot = gfq;
		dend = pot + ndata;
	    }
	    for (int p = 0, n = 0; pot < dend; ++pot) {
		if (nphase == 3 && p == 0) {
		    int	c = n++;
		    for (int m = ndata; m /= 4; ) {
			fputc(Nucl[c / m], out_fd);
			c %= m;
		    }
		    fputc('\t', out_fd);
		}
		if (omode) fprintf(out_fd, "%15.7e", *pot);
		else	fprintf(out_fd, "%7d", (int) *pot);
		if (nphase == 3) {
		    if ((p = next_p[p]) == 0) fputc('\n', out_fd);
		    else	fputc('\t', out_fd);
		} else		fputc('\n', out_fd);
	    }
	}
}

int main(int argc, const char** argv)
{
	int	mo = 4;		// default MM order
	bool	cds = false;
	int	exin = static_cast<int>(Iefp::NG);
	while (--argc > 0 && **++argv == OPTCHAR) {
const	    char*	opt = argv[0] + 1;
	    int	c = *opt;

	    if (!c) break;
const	    char*	val = argv[0] + 2;
	    switch (c) {
		case 'b':
		    if ((val = getarg(argc, argv)))
			oname = val;
		    binary = true;
		    break;
		case 'c':
		    if ((val = getarg(argc, argv)))
		        ewdfq = val;
		    else	cds = true;
		    break;
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
		case 'g':
		    if ((val = getarg(argc, argv)))
		        gwdfq = val;
		    break;
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
		        oname = val;
		    break;
		case 'C': 
		    if ((val = getarg(argc, argv, true)))
			initcodon(atoi(val));
		    break;
		case 'J': 
		    if ((val = getarg(argc, argv, true)))
			min_orf = atoi(val);
		    break;
		case 'O':
		    if ((val = getarg(argc, argv, true)))
		        omode = atoi(val);
		    break;
		case 'T':
		    if ((val = getarg(argc, argv)))
			ftable.setpath(val);
		    break;
		default: break;
	    }
	}

	if (fname) {
	    int	file_type = 0;
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
	bool	get_kmers = iwdfq || ewdfq || eij || argc;
	if (oname) {		// test compatibility
	    int	file_type = 1;
const	    int	exn = fname2exin(oname, file_type);
	    if (file_type == 0)	// undefined file type
		oname = add_ext(oname, iefp_ext[exin], str);
	    else if ((fname && exn > exin) || (exin - exn) / 3) usage();
	    else exin = exn;
	    if (!binary) oname = add_ext(oname, text_ext, str);
	}
	outfd(oname);		// setup out_fd
	Ipt	ipt(exin, mo, ((exin / 3) == 2)? 3: 1);
	float*	gfq = 0;
	if (gwdfq) {
	    gfq = ipt.getKmers(gwdfq, false);
	    ipt.reform(gfq);
	}
	EiJuncSeq*	eijseq = eij? new EiJuncSeq(INTRON, eij, gnm): 0;
	if (ewdfq) iwdfq = ewdfq;
	if (fname && !ipt.readFile(fname)) usage();
	else if (get_kmers) {
	    if (iwdfq)	ipt.getKmers(iwdfq); else
	    if (eij)	ipt.getKmers(eijseq); else
	    if (argc)	ipt.getKmers(argc, argv);
	    if (omode)	ipt.reform();
	    if (gfq && !ipt.makeExinPot(gfq)) usage();
	} else if (!gwdfq && !fname)
	    usage();

	if (!ipt.ispot()) {
// nothing to do
	} else if (eij) {
	    eijseq->reset();
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
