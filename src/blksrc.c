/*****************************************************************************
*
*	Subroutines for block search
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
*	Osamu Gotoh, Ph.D.      (2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "aln.h"
#include "utilseq.h"
#include "wln.h"
#include "blksrc.h"

#define	Undef	0

static	const	int	PRUNE = 2;
static	const	int	TFACTOR = 100;
static	const	int	NCAND = 10;
static	const	INT	ddelim = SEQ_DELIM + (SEQ_DELIM << 4);
static	const	int	genlencoeff = 36;
static	const	INT	MinMaxGene = 16384;

extern	CHAR	gencode[];

struct WCPRM {
	INT	Nalpha;
	INT	Ktuple;
	INT	Bitpat2;
	INT	TabSize;
	INT	BitPat;
	INT	Nshift;
	INT	blklen;
	INT	MaxGene;
	short	Nbitpat;
	SHORT	afact;		// abundant / even
};

static	int	cmpf(BLKTYPE* a, BLKTYPE* b);
static	int	scmpf(KVpair<INT, int>* a, KVpair<INT, int>* b);
static	int	gcmpf(GeneRng* a, GeneRng* b);

//		  elem tuple bitpat2 tabsize bitpat shift blklen maxgene nbitpat afact
static	WCPRM	wcp = {0, 0, 1, 0, 1, 0, 0, 0, 0, 0};
static	WCPRM	wcp_af = {20, 5, 0, 3200000, 0, 5, 4096, 0, 1, 10};
static	WCPRM	wcp_ax = {18, 4, 27, 104976, 29, 4, 4096, 0, 4, 10};
static	WCPRM	wcp_cf = {4, 8, 0, 65536, 0, 8, 4096, 0, 1, 10};
static	WCPRM	wcp_cx = {4, 8, 3255, 65536, 7573, 8, 4096, 0, 5, 10};
static	char	dir[2] = {'>', '<'};
static	char	emssg[] = "%s is incompatible !\n";
static	char	mfmt[] = "Gd:%u No:%u My:%u MS:%d Mb:%u Tw:%lu Tl:%u %6.2f %6.2f %6.2f\n";

static	int	gene_rng_max_extend = -1;
static	const	INT	blkmergin = 1024;
static	float	hfact = 1.25;
static	int	gratio = 10;	// ratio of average gene and cDNA
static	const	char*	WriteFile = 0;
static	SHORT	Ncand = NCAND;	// max number of candidates to pass to 2nd phase
static	const	char*	ConvPat = 0;
static	const	char*	sBitPat = 0;
static	float	aaafact = 1;
static	SHORT	MinSigpr = 3;
static	INT	MaxNref = 16;
static	int	AllRef = 0;
static	int	blockpam = 20;
static	const	char*	blockmdm = 0;
static	float	AllowedOverlapFact = 0.2;
static	INT	MaxMmc = 15;
static	float	RbsBias = 3.; 
static	const	float	RbsFactSqrX = 0.0;
static	float	RbsBase = 3.;
static	float	RbsBaseX = 0.5;
static	float	RbsFact = DQUERY;
static	int	Nascr = 2;
const	char*	writeerrmssg = "Fail to write to %s !\n";

void setQ4prm(const char*  ps)
{
const	char*	val = ps + 1;

	switch (*ps) {
	    case 'A':	// aa classification pattern
		if (isdigit(*val)) wcp.Nalpha = atoi(val);
		else if (*val)	ConvPat = val;
		break;
	    case 'B':	// matching bit pattern
		if (*val) sBitPat = val;
		else sBitPat = "";
		break;
	    case 'C':	// number of bit patterns
		wcp.Nbitpat = atoi(val);
		if (wcp.Nbitpat <= 0) wcp.Nbitpat = 1;
		break;
	    case 'E':	// max number of extension
		gene_rng_max_extend = atoi(val);
		if (gene_rng_max_extend <= 0) gene_rng_max_extend = 1;
		break;
	    case 'G':	// Maximum gene length
		if (*val) wcp.MaxGene = ktoi(val);
		break;
	    case 'W':	// Write block information to the file
		if (*val) WriteFile = val;
		break;
	    case 'a':	// abundant / average ratio
		if (*val) wcp.afact = (SHORT) atoi(val);
		break;
	    case 'b':	// block length
		if (*val) wcp.blklen = ktoi(val);
		break;
	    case 'c':	// factor for hash mergin
		if (*val) hfact = atof(val);
		break;
	    case 'd':	// N-best all hit block score
		if (*val) Nascr = atoi(val);
		break;
	    case 'e':	// base of extreme value distribution
		if (*val) RbsBase = atof(val);
		break;
	    case 'f':	// factor of extreme value distribution
		if (*val) RbsFact = atof(val);
		break;
	    case 'g':	// ratio between average gene and mRNA lengths
		if (*val) gratio = atoi(val);
		break;
	    case 'h':	// bias of phase1 threshold
		if (*val) RbsBias = atof(val);
		break;
	    case 'k':	// tuple size
		if (*val) wcp.Ktuple = atoi(val);
		break;
	    case 'm':	// max allowed mishits
		if (*val) MaxMmc = atoi(val);
		break;
	    case 'n':	// number of additional recurrences
		MinSigpr = (*val)? atoi(val): SHRT_MAX;
		break;
	    case 'o':	// max number of ORFs
		if (*val) MaxNref = atoi(val);
		break;
	    case 'p':	// pam level
		if (*val) {
		    if (isdigit(*val)) blockpam = atoi(val);
		    else	blockmdm = val;
		}
		break;
	    case 'q':	// use average aa composition
		if (*val) aaafact = atof(val);
		break;
	    case 's':	// word shift
		if (*val) wcp.Nshift = atoi(val);
		break;
	    case 'u':	// toggle paralog output mode
		AllRef = 1 - AllRef;
		break;
	    case 'w':	// wln level used in TestOutput
		if (*val) algmode.lvl = atoi(val);
		if (algmode.lvl > 2) algmode.lvl = 2;
		break;
	    default: break;
	}
}

ContBlk::ContBlk()
{
	VerNo = BsVersion;
	clean();
}

ContBlk::~ContBlk()
{
	delete[] Nblk; delete[] blkp; delete[] blkb;
	delete[] wscr; delete[] ChrID;
}

// find words in each block

void Chash::countBlk()
{
	KVpair<INT, int>* h = hash;
	for ( ; h < hz; ++h) {
	    if (h->val) {
		pwc->Nblk[h->key]++;
		h->val = 0;	// reset hash
	    }
	}
}

void Chash::registBlk(INT m)
{
	for (KVpair<INT, int>* h = hash; h < hz; ++h) {
	    if (h->val) {
		INT*	wrk = pwc->blkp[h->key];
		if (wrk) {
		    INT	Nblk = pwc->Nblk[h->key]++;
		    wrk[Nblk] = m;
		}
		h->val = 0;	// reset hash
	    }
	}
}

// Word count table

WordTab::WordTab(Seq* sd, INT kmer, INT nsft, INT elms, const char* ap, 
	INT bp, INT bp2, INT nbt)
	: ReducWord(sd, elms, ap), Nbitpat(nbt), 
	BitPat(bp), BitPat2(bp2), Nshift(nsft), total(0),
	ch(0), ss6(0), ww(0), idx(0), count(0)
{
	bpp = new Bitpat_wq*[Nbitpat];
	int	nfrm = sd->istron()? 6: 1;
	Nss = nfrm * Nbitpat;
	ss = new SHORT[Nss];
	vclear(ss, Nss);
	if (nfrm == 6) {
	    ss6 = new SHORT*[6];
	    ww = new INT[Nbitpat];
	    vclear(ww, Nbitpat);
	    ss6[0] = ss;
	    for (int k = 1; k < 6; ++k)
		ss6[k] = ss6[k - 1] + Nbitpat;
	}
	if (Nbitpat == 1) bpp[0] = new Bitpat_wq(Nalpha, nfrm, 0, BitPat);
	else {
	    bpp[0] = new Bitpat_wq(Nalpha, nfrm, 0, bitmask(kmer));
	    for (INT k = 1; k < Nbitpat; ++k)
		bpp[k] = new Bitpat_wq(Nalpha, nfrm, (k - 1) % 2, k < 3? BitPat: BitPat2);
	}
	INT	TabSiz = ipower(Nalpha, kmer);
	tcount = new INT[TabSiz];
	vclear(tcount, TabSiz);
}

WordTab::~WordTab()
{
	for (INT k = 0; k < Nbitpat; ++k) delete bpp[k];
	delete[] bpp; delete[] ww; delete[] ss; delete[] ss6;
	delete[] count; delete[] tcount; delete[] idx; delete ch;
}

bool WordTab::c2w(int c, int  mode, INT i)	// letter to word
{
	INT	uc = aConvTab[toupper(c) - 'A'];
	if (uc > Nalpha) return (false);	// ignore
	if (uc == Nalpha) {			// ambigous
	    for (INT k = 0; k < Nbitpat; ++k) {
		ss[k] = 0;
		bpp[k]->flaw();
	    }
	} else {				// ok
	    for (INT k = 0; k < Nbitpat; ++k) {
		INT	w = bpp[k]->word(uc);
		if (bpp[k]->flawless()) {
		    if (tcount && (mode == 0 || mode == 2)) ++tcount[w];
		    if (!ss[k]) {
			switch (mode) {
			  case 1: idx[k][count[w]] = i;
			  case 0: ++count[w]; break;
			  case 2: 
			  case 3: ch->incr(w); break;
			}
		    }
		    if (++ss[k] == Nshift) ss[k] = 0;
		}
	    }
	}
	return (true);
}

bool WordTab::c2w6(int c, int& p,int  mode, INT i)	// letter to word
{
	INT	uc = ntconv[toupper(c) - 'A'];
	if (uc > 4) return (false);	// ignore
	cc[p] = cc[p + 3] = 0;		// initialize codon
	if (uc == 4) {			// ambiguous
	    for (int q = 0; q < 3; ++q) xx[q] = 4;
	} else {
	    for (int q = 0; q < 3; ++q) {
		cc[q] = ((cc[q] << 2) + uc) & 63;
		cc[q + 3] = (((3 - uc) << 4) + (cc[q + 3] >> 2)) & 63;
		xx[p] >>= 1;
	    }
	}
	if (++p == 3) p = 0;
	for (int q = p; q < 6; q += 3) {
	    if (xx[p]) {		// omit
		for (INT k = 0; k < Nbitpat; ++k) {
		    ss6[q][k] = 0;
		    bpp[k]->flaw(q);
		}
	    } else {
		uc = aConvTab[acodon[gencode[cc[q]]] - 'A'];
		if (bpp[0]->good(uc)) {
		    for (INT k = 0; k < Nbitpat; ++k)
			ww[k] = bpp[k]->word(uc, q);
		} else {		// termination codon
		    for (INT k = 0; k < Nbitpat; ++k) {
			ss6[q][k] = 0;
			bpp[k]->flaw(q);
		    }
		}
	    }
	    for (INT k = 0; k < Nbitpat; ++k) {
		if (bpp[k]->flawless(q)) {
		    INT	w = ww[k];
		    if (tcount && (mode == 0 || mode == 2)) ++tcount[w];
		    if (!ss6[q][k]) {
			switch (mode) {
			  case 1: idx[k][count[w]] = i;
			  case 0: ++count[w]; break;
			  case 2: 
			  case 3: ch->incr(w); break;
			}
		    }
		    if (++ss6[q][k] == Nshift) ss6[q][k] = 0;
		}
	    }
	}
	return (true);
}

void WordTab::reset()
{
	vclear(ss, Nss);
	if (ss6) {
	    vclear(cc, 6);
	    vset(xx, INT(4), 3);
	}
	for (INT k = 0; k < Nbitpat; ++k)
	    bpp[k]->clear();
}

void WordTab::countwd(const char* argv[])
{
	total = 0;
	int	c;
	for (const char** genome = argv; *genome; ++genome) {
	    INT	i = 0;
	    FILE*	fd = fopen(*genome, "r");
	    if (!fd) fatal("%s not found!\n", *genome);
	    reset();
	    while ((c = fgetc(fd)) != EOF) {
		if (isalpha(c)) {
		    if (c2w(c)) ++i;
		} else switch (c) {
		  case '>': reset();	// FASTA Header
		  case ';': case '#':	// comments
		    while ((c = fgetc(fd)) != EOF && c != '\n');
		    break;
		  case '/':
		    if ((c = fgetc(fd)) == '/')
			while ((c = fgetc(fd)) != EOF && c != '>');
		    ungetc(c, fd);
		  default: break;
		}
	    }
	    fclose(fd);
	}
}

void MakeBlk::findChrBbound()
{
	double	B = (double) (chrblk(pwc->ChrNo) - 1);

	pbc->BClw = pbc->BCup = 0.;
	for (size_t k = 0; k <= pwc->ChrNo; ++k) {
	    double	offdiag = k * B - pwc->ChrNo * (chrblk(k) - 1);
	    if (offdiag < pbc->BClw) pbc->BClw = offdiag;
	    if (offdiag > pbc->BCup) pbc->BCup = offdiag;
	}
	pbc->BClw /= B;
	pbc->BCup /= B;
}

template <typename file_t>
void MakeBlk::writeBlkInfo(file_t fd, const char* fn)
{
	if (!fd) fatal(writeerrmssg, fn);
	size_t	chrn = pwc->ChrNo + 1;
	if (!pwc->ChrID) ++pwc->VerNo;
	if (fwrite(&wcp, sizeof(WCPRM), 1, fd) != 1 ||
	    fwrite(pwc, sizeof(ContBlk), 1, fd) != 1 ||
	    fwrite(pbc, sizeof(Block2Chr), 1, fd) != 1 ||
	    (pwc->ChrID && fwrite(pwc->ChrID, sizeof(CHROMO), chrn, fd) != chrn) ||
	    fwrite(pwc->Nblk, sizeof(SHORT), wcp.TabSize, fd) != wcp.TabSize)
		fatal(writeerrmssg, fn);
	INT*	blkidx = new INT[wcp.TabSize];
	for (INT w = 0; w < wcp.TabSize; ++w)
	    blkidx[w] = pwc->blkp[w]? pwc->blkp[w] - pwc->blkb + 1: 0;
	if (fwrite(blkidx, sizeof(INT), wcp.TabSize, fd) != wcp.TabSize)
		fatal(writeerrmssg, fn);
	delete[] blkidx;
	SHORT*	sblkb = (SHORT*) pwc->blkb;
	if (fwrite(sblkb, sizeof(SHORT), pwc->WordSz, fd) != pwc->WordSz ||
	    fwrite(pwc->wscr, sizeof(short), wcp.TabSize, fd) != wcp.TabSize ||
	    fwrite(iConvTab, sizeof(int), pwc->ConvTS, fd) != pwc->ConvTS)
		fatal(writeerrmssg, fn);
	fclose(fd);
}

void MakeBlk::WriteBlkInfo()
{
static	const char* wfmt =
	"#Segs %u, TabSize %u, Words: %lu, GenomeSize %lu, GIDs %d\n";

	const char*	fn = WriteFile;
	if (!fn) fatal("Specify write file !\n");
	prompt(wfmt, chrblk(pwc->ChrNo) - 1,
	    wcp.TabSize, pwc->WordNo, pwc->glen, pwc->ChrNo);

	SHORT*	sblkb = (SHORT*) pwc->blkb;
	if (pwc->BytBlk == 2) {
	    SHORT*	sblk = sblkb;
	    INT*	iblk = pwc->blkb;
	    for (size_t i = 0; i < pwc->WordNo; ++i)
		*sblk++ = (SHORT) *iblk++;
	}
	findChrBbound();
#if USE_ZLIB
	if (is_gz(fn)) {
	    gzFile	gzfd = gzopen(fn, "wb");
	    if (!gzfd) fatal(no_file, fn);
	    writeBlkInfo(gzfd, fn);
	    return;
	}
#endif
	FILE*	fd = fopen(fn, "wb");
	if (!fd) fatal(no_file, fn);
	writeBlkInfo(fd, fn);
}

MakeBlk::MakeBlk(Seq* sd, DbsDt* dd) :
	WordTab(sd, wcp.Ktuple, wcp.Nshift, wcp.Nalpha,
	    ConvPat, wcp.BitPat, wcp.Bitpat2, wcp.Nbitpat), wdbf(dd)
{
	pwc = new ContBlk;
	pbc = new Block2Chr;
	vclear(pbc);
	wcp.Nalpha = Nalpha;
	wcp.TabSize = bpp[0]->TabSize;
	pwc->ConvTS = ConvTabSize;
	pwc->Nblk = new SHORT[wcp.TabSize];
	vclear(pwc->Nblk, wcp.TabSize);
	deltaa = 0.;
	acomp = 0;
	ch = new Chash(int(hfact * wcp.blklen), pwc);
	int	molc = sd->inex.molc;
	if (!molc && dd) molc = dd->curdb->defmolc;
	isaa = molc == PROTEIN;
	istron = molc == TRON;
	defcode = setSeqCode(0, isaa? PROTEIN: DNA);
}

static void setupbitpat(int molc, size_t gnmsz)
{
	if (molc == PROTEIN || molc == TRON) {
	    WCPRM&	wcp_t = algmode.crs? wcp_ax: wcp_af;
	    if (wcp.Nalpha == 0) wcp.Nalpha = wcp_t.Nalpha;
	    if (wcp.Ktuple == 0 && !gnmsz) wcp.Ktuple = wcp_t.Ktuple;
	    if (wcp.Bitpat2 == 1) wcp.Bitpat2 = wcp_t.Bitpat2;
	    if (wcp.BitPat == 1) wcp.BitPat = wcp_t.BitPat;
	    if (wcp.blklen == 0 && !gnmsz) wcp.blklen = wcp_t.blklen;
	    if (wcp.Nbitpat == 0) wcp.Nbitpat = wcp_t.Nbitpat;
	    if (wcp.afact == 0) wcp.afact = wcp_t.afact;
	} else {
	    WCPRM&	wcp_t = algmode.crs? wcp_cx: wcp_cf;
	    if (wcp.Nalpha == 0) wcp.Nalpha = wcp_t.Nalpha;
	    if (wcp.Ktuple == 0 && !gnmsz) wcp.Ktuple = wcp_t.Ktuple;
	    if (wcp.Bitpat2 == 1) wcp.Bitpat2 = wcp_t.Bitpat2;
	    if (wcp.BitPat == 1) wcp.BitPat = wcp_t.BitPat;
	    if (wcp.blklen == 0 && !gnmsz) wcp.blklen = wcp_t.blklen;
	    if (wcp.Nbitpat == 0) wcp.Nbitpat = wcp_t.Nbitpat;
	    if (wcp.afact == 0) wcp.afact = wcp_t.afact;
	}
	if (gnmsz && !(wcp.Ktuple && wcp.blklen && wcp.MaxGene)) {
	    if (!wcp.blklen) {
		wcp.blklen = int(sqrt(double(gnmsz)));
		wcp.blklen = int(wcp.blklen / 1024 + 1) * 1024;
	    }
	    if (!wcp.Ktuple) {
		double	loggs = log(double(gnmsz));
		switch (molc) {
		  case PROTEIN:	loggs *= 0.30; break;
		  case TRON:	loggs *= 0.27; break;
		  default:	loggs *= 0.59; break;
		}
		wcp.Ktuple = int(loggs);
		if (molc == PROTEIN) {
		    if (wcp.Ktuple < 3) wcp.Ktuple = 3;
		    if (wcp.Ktuple > 6) wcp.Ktuple = 6;
		}
	    }
	    if (!wcp.MaxGene) {
		wcp.MaxGene = int(genlencoeff * sqrt(double (gnmsz)) / 1024 + 1) * 1024;
		if (wcp.MaxGene < MinMaxGene) wcp.MaxGene = MinMaxGene;
	    }
	}
	wcp.TabSize = ipower(wcp.Nalpha, wcp.Ktuple);
	if (wcp.Nshift == 0) wcp.Nshift = wcp.Ktuple;
	if (sBitPat && *sBitPat != '1') sBitPat = DefBitPat[wcp.Ktuple];
	if (sBitPat) {
	    INT	w = 0;
	    wcp.BitPat = bpcompress(sBitPat, w);
	    if (w > wcp.Ktuple) wcp.Ktuple = w;
	} else {
	    wcp.BitPat = bitmask(wcp.Ktuple);
	    wcp.Nbitpat = 1;
	}
	if (wcp.Nbitpat > 3) {
	    const char*	bp2 = strchr(sBitPat, ',');
	    if (bp2) {
		INT	w = 0;
		wcp.Bitpat2 = bpcompress(bp2 + 1, w);
		if (w > wcp.Ktuple) wcp.Ktuple = w;
	    } else {
		wcp.Nbitpat = 3;
	    }
	}
}

MakeBlk* makeblock(int argc, const char** argv, int molc)
{
	Seq	sd(1);
	setSeqCode(&sd, molc);
	DbsDt*	wdbf = new DbsDt('f', molc);
	setupbitpat(molc, file_size(*argv));
	MakeBlk*	mb = new MakeBlk(&sd, wdbf);
	mb->idxblk(argc, argv, molc);
	return (mb);
}

MakeBlk* makeblock(SeqServer* svr)
{
	int	molc = svr->getmolc();
	Seq	sd(1);
	setSeqCode(&sd, molc);
	DbsDt*	wdbf = new DbsDt('f', molc);
	setupbitpat(molc, svr->total_seq_len(&sd));
	MakeBlk*	mb = new MakeBlk(&sd, wdbf);
	mb->idxblk(&sd, svr);
	return (mb);
}

MakeBlk* makeblock(Seq* sd)
{
	if (!sd || sd->many < 2) fatal("Irmnproper input seq !\n");
	DbsDt*	wdbf = new DbsDt('f', sd->inex.molc);
	setupbitpat(sd->inex.molc, sd->len * sd->many);
	MakeBlk*	mb = new MakeBlk(sd, wdbf);
	mb->idxblk(sd);
	return (mb);
}

void MakeBlk::prepacomp()
{
	DefPrm	dp = {0., 0., 0., 0., blockpam};
	Simmtx* pm = blockmdm?	new Simmtx(PxP, blockmdm, &dp): 
				new Simmtx(PxP, &dp);
	if (pm && pm->Nrmlf() == 0) fatal("Run makmdm !\n");
	double	mdmfact = pm? pm->Nrmlf(): 1;
	if (aaafact > 0) {
	    double*	pcomp = get_mdmcmp();
	    double	mdmaa[20];
	    acomp = new double[wcp.Nalpha];
	    vclear(acomp, wcp.Nalpha);
	    if (pm) vclear(mdmaa, 20);
	    for (int a = 0; a < 20; ++a) {
		int	p = a + ALA;
		int	q = iConvTab[p];
		acomp[q] += pcomp[a];
		if (pm) mdmaa[q] += pcomp[a] * (pm->mtx[p][p] / mdmfact + log(pcomp[a]));
	    }
	    double	avr = 0.;
	    for (INT q = 0; q < wcp.Nalpha; ++q) {
		if (pm) 
		    acomp[q] = mdmaa[q] / acomp[q];
		else
		    acomp[q] = -log(acomp[q]);
		avr += acomp[q];
	    }
	    avr /= wcp.Nalpha;
	    for (INT q = 0; q < wcp.Nalpha; ++q)
		acomp[q] = TFACTOR * aaafact * (acomp[q] - avr);
	    deltaa = acomp[0] - acomp[wcp.Nalpha - 1];
	    if (pm) delete pm;
	}
}

void MakeBlk::blkscrtab(size_t segn)
{
	INT	m = 0;
	INT	c = 0;
	INT	i = 0;
	double	avr = 0;
	double	alc = (aaafact > 0)? bpp[0]->weight * acomp[0]: 0;
	double	basescr = log((double) segn);
	for (INT w = m = 0; w < wcp.TabSize; ++w) {
	    if (tcount[w]) {
		++m;
		short	wscr = (short) (TFACTOR * (basescr - log((double) tcount[w] / Nbitpat)));
		if (aaafact > 0) wscr += (short) alc;
		pwc->wscr[w] = wscr;
		avr += wscr;
	    } else pwc->wscr[w] = 0;
	    if (aaafact > 0) {
		int	p = 0, q = 0;
		for (INT x = w + 1; (q = x % wcp.Nalpha) == 0; x /= wcp.Nalpha)
		    ++p;
		if (p) alc += p * deltaa;	// (q-1)AAA
		alc += acomp[q] - acomp[q-1];	//     q000
	    }
	}
	pwc->AvrScr = (SHORT) (avr /= m);
	short	MinScr = (short) (avr - TFACTOR * (1 + aaafact) * log((double) wcp.afact));
	if (MinScr < 0) MinScr = 0;

	for (INT w = i = c = m = pwc->WordNo = 0; w < wcp.TabSize; ++w) {
	    if (pwc->Nblk[w] == 0) {
		++c;
		pwc->wscr[w] = -1;
	    } else if (pwc->wscr[w] > MinScr) {
		++m;
		pwc->WordNo += pwc->Nblk[w];
		if (pwc->Nblk[w] > pwc->MaxBlk) pwc->MaxBlk = pwc->Nblk[w];
	    } else {
		++i;
		pwc->wscr[w] = 0;
	    }
	}
	if (WriteFile && segn <= USHRT_MAX) {
	    pwc->BytBlk = 2;
	    pwc->WordSz = pwc->WordNo;
	} else {
	    pwc->BytBlk = 4;
	    pwc->WordSz = 2 * pwc->WordNo;
	}
	INT*	wrk = pwc->blkb = new INT[pwc->WordNo];
	size_t	x = 0;
	for (INT w = 0; w < wcp.TabSize; ++w) {
	    x += tcount[w];
	    if (pwc->Nblk[w]) {
		if (pwc->wscr[w] > MinScr) {
		    pwc->blkp[w] = wrk;
		    wrk += pwc->Nblk[w];
		} else
		    pwc->blkp[w] = 0;
		pwc->Nblk[w] = 0;	// reset for next round
	    } else  pwc->blkp[w] = 0;
	}
	prompt(mfmt, m, c, i, MinScr, pwc->MaxBlk, pwc->WordNo, x,
	    100. * m / wcp.TabSize, 100. * i / wcp.TabSize, 100. * c / wcp.TabSize);
}

void MakeBlk::blkscrtab(size_t segn, INT blksz)
{
	INT	m = 0;
	INT	c = 0;
	INT	i = 0;
	double	avr = 0.;
	for (INT w = 0; w < wcp.TabSize; ++w)
	    if (tcount[w]) ++m;
	double	basescr = log((double) segn);
	short	MinScr = (short) -(TFACTOR * log((double) wcp.afact * blksz / m));
	if (MinScr < 0) MinScr = 0;
	for (INT w = i = c = m = pwc->WordNo = 0; w < wcp.TabSize; ++w) {
	    if (pwc->Nblk[w]) {
		short	wscr = (short) (TFACTOR * (basescr - log((double) tcount[w] / Nbitpat)));
		if (wscr > MinScr) {
		    ++m;
		    pwc->WordNo += pwc->Nblk[w];
		    pwc->wscr[w] = wscr;
		    avr += wscr;
		    if (pwc->Nblk[w] > pwc->MaxBlk) pwc->MaxBlk = pwc->Nblk[w];
		} else {
		    ++i;
		    pwc->wscr[w] = MinScr;
		}
	    } else {
		++c;
		pwc->wscr[w] = -1;
	    }
	}
	if (WriteFile && segn <= USHRT_MAX) {
	    pwc->BytBlk = 2;
	    pwc->WordSz = pwc->WordNo;
	} else {
	    pwc->BytBlk = 4;
	    pwc->WordSz = 2 * pwc->WordNo;
	}
	pwc->AvrScr = (SHORT) (avr / m);
	size_t	x = 0;
	INT*	wrk = pwc->blkb = new INT[pwc->WordNo];
	for (INT w = 0; w < wcp.TabSize; ++w) {
	    x += tcount[w];
	    if (pwc->Nblk[w]) {
		if (pwc->wscr[w] > MinScr) {
		    pwc->blkp[w] = wrk;
		    wrk += pwc->Nblk[w];
		} else
		    pwc->blkp[w] = 0;
		pwc->Nblk[w] = 0;	// reset for next round
	    } else  pwc->blkp[w] = 0;
	}
	prompt( mfmt, m, c, i, MinScr, pwc->MaxBlk, pwc->WordNo, x,
	    100. * m / wcp.TabSize, 100. * i / wcp.TabSize, 100. * c / wcp.TabSize);
}

// Read from amino or nucleotide sequence files

template <typename file_t>
void MakeBlk::first_phase(file_t fd, int& entlen, CHROMO& chrbuf)
{
	INT	i = 0; int p = 0; reset();
	int	c;
	while ((c = fgetc(fd)) != EOF) {
	    if (isalpha(c)) {
		bool ok = istron? c2w6(c, p, 2): c2w(c, 2);
		if (ok && ++i == wcp.blklen) {
		    ch->countBlk();	// block boundary
		    ++chrbuf.segn;
		    pwc->glen += i;
		    i = 0;
		}
	    } else switch (c) {
		case '>':		// FASTA Header
		    ++pwc->ChrNo;
		    if (i) {
			ch->countBlk();
			++chrbuf.segn;
			pwc->glen += i;
			i = 0;
		    }
		    i = 0; reset();
		    if (wdbf) {
			while ((c = fgetc(fd)) != EOF && !isspace(c))
			    ++entlen;
			if (c == '\n') break;
		    }
		case ';': case '#':	// comments
		    while ((c = fgetc(fd)) != EOF && c != '\n');
		    break;
		case '/': 
		    if ((c = fgetc(fd)) == '/')
			while ((c = fgetc(fd)) != EOF && c != '>');
		    ungetc(c, fd);
		    break;
		default: break;
	    }
	}
// last block
	if (i) {
	    ch->countBlk();
	    ++chrbuf.segn;
	    pwc->glen += i;
	}
	fclose(fd);
}

template <typename file_t>
int MakeBlk::second_phase(file_t fd, CHROMO& chrbuf, int m,
		DbsRec*& pr, CHAR*& ps, char*& pe,
		DbsRec& prvrec, CHROMO*& pchrid)
{
	int	bias = A - 1;
	if (isaa) bias = ALA - 1;
	size_t	n = 0, i = 0; int b = 0, p = 0;
	int	c;
	while ((c = fgetc(fd)) != EOF) {
	    if (isalpha(c)) {
		if (istron? c2w6(c, p, 3): c2w(c, 3)) {
		    ++n;
		    if (++i == wcp.blklen) {
			ch->registBlk(++m);	// block boundary
			++chrbuf.segn;
			chrbuf.spos += i;
			i = 0;
		    }
		    if (wdbf) {
			c = defcode->encode[toupper(c) - 'A'] - bias;
			if (isaa)	*ps++ = c;
			else if (n & 1)	b = c << 4;
			else 	*ps++ = b + c;
		    }
		}
	    } else switch (c) {
		case '>':
		    if (i) {
			ch->registBlk(++m);
			++chrbuf.segn;
			chrbuf.spos += i;
			i = 0;
		    }
		    chrbuf.fpos = ftell(fd) - 1;
		    *pchrid++ = chrbuf;
		    if (wdbf) {
			if (n) {
			    if (isaa)	*ps++ = SEQ_DELIM;
			    else if (n & 1) *ps++ = b + SEQ_DELIM;
			    else	*ps++ = ddelim;
			    prvrec.seqlen = n;
			    *pr++ = prvrec;
			}
			prvrec.entptr = pe - wdbf->entry;
			prvrec.seqptr = ps - wdbf->dbsseq;
			while ((c = fgetc(fd)) != EOF && !isspace(c))
			    *pe++ = c;
			*pe++ = '\0';
		    }
		    n = b = 0; reset();
		    if (c == '\n') break;
		case ';': case '#':	// comments
		    while ((c = fgetc(fd)) != EOF && c != '\n');
		    break;
		case '/': 
		    if ((c = fgetc(fd)) == '/')
			while ((c = fgetc(fd)) != EOF && c != '>');
		    ungetc(c, fd);
		    break;
		default: break;
	    }
	}
// last block
	if (i) {
	    ch->registBlk(++m);
	    ++chrbuf.segn;
	    chrbuf.spos += i;
	}
	chrbuf.fpos = ftell(fd);
	if (wdbf && n) {
	    if (isaa)	*ps++ = SEQ_DELIM;
	    else if (n & 1) *ps++ = b + SEQ_DELIM;
	    else	*ps++ = ddelim;
	    prvrec.seqlen = n;
	    *pr++ = prvrec;
	}
	fclose(fd);
	return (m);
}

void MakeBlk::idxblk(int argc, const char** argv, int molc)
{
// first phase

	if (istron) prepacomp();
	int	entlen = 0;
	pwc->ChrNo = pwc->glen = 0;
	CHROMO	chrbuf = {0L, 0, 0};
	const char**	seqdb = argv;
	for (int ac = 0; ac++ < argc && *seqdb; ++seqdb) {
	    bool	gz = is_gz(*seqdb);
	    if (gz) {
#if USE_ZLIB
		gzFile	gzfd = gzopen(*seqdb, "rb");
		if (gzfd)
	            first_phase(gzfd, entlen, chrbuf);
		else	fatal("%s not found!\n", *seqdb);
#else
		fatal(gz_unsupport);
#endif
	    } else {
		FILE*	fd = fopen(*seqdb, "r");
		if (!fd) fatal("%s not found!\n", *seqdb);
	        first_phase(fd, entlen, chrbuf);
	    }
	}
	pwc->ChrID = new CHROMO[pwc->ChrNo + 1];
	pwc->blkp = new INT*[wcp.TabSize];
	pwc->wscr = new short[wcp.TabSize];
	if (wdbf) {
	    size_t	sz = isaa? pwc->glen: (pwc->glen / 2 + pwc->ChrNo);
	    wdbf->prepare(entlen + pwc->ChrNo, pwc->ChrNo, sz + pwc->ChrNo);
	}
	if (istron)	blkscrtab(chrbuf.segn);
	else	blkscrtab(chrbuf.segn, pwc->glen / chrbuf.segn);

// second phase

	INT	m = 0;
	chrbuf.segn = 1;	// skip 0
	DbsRec*	pr = wdbf? wdbf->recidx: 0;
	CHAR*	ps = wdbf? wdbf->dbsseq: 0;
	char*	pe = wdbf? wdbf->entry: 0;
	DbsRec	prvrec = {0L, 0, 0};
	if (wdbf) pr->seqptr = 0;
	CHROMO*	pchrid = pwc->ChrID;
	seqdb = argv;
	for (int ac = 0; ac++ < argc && *seqdb; ++seqdb) {
	    bool	gz = is_gz(*seqdb);
	    if (gz) {
#if USE_ZLIB
		gzFile	gzfd = gzopen(*seqdb, "rb");
	        m = second_phase(gzfd, chrbuf, m, pr, ps, pe, prvrec, pchrid);
#endif
	    } else {
		FILE*	fd = fopen(*seqdb, "r");
		if (!fd) fatal("%s not found!\n", *seqdb);
	        m = second_phase(fd, chrbuf, m, pr, ps, pe, prvrec, pchrid);
	    }
	}
	*pchrid++ = chrbuf;
}

//	read from SeqServer

void MakeBlk::idxblk(Seq* sd, SeqServer* svr)
{
	bool	isaa = sd->isprotein();
	int	bias = (isaa? int(ALA): int(A)) - 1;
	char*	cvt =  isaa? acodon: nucl;
	int	maxlen = 0;

// first phase

	int	entlen = 0;
	int	pfqnum = 0;
	pwc->ChrNo = pwc->glen = 0;
	InSt	ist;
	while ((ist = svr->nextseq(sd, 0)) != IS_END) {
	    if (ist == IS_ERR) continue;
	    if (sd->len > maxlen) maxlen = sd->len;
	    CHAR*	ts = sd->at(sd->len);
	    for (CHAR* ps = sd->at(0); ps < ts; ++ps)
		if (c2w(cvt[*ps], 2)) ++pwc->glen;
	    ++pwc->ChrNo;
	    ch->countBlk();
	    entlen += strlen((*sd->sname)[0]);
	    if (sd->sigII) pfqnum += sd->sigII->pfqnum;
	    reset();
	}
	pwc->ChrID = 0;
	wcp.blklen = maxlen;
	pwc->blkp = new INT*[wcp.TabSize];
	pwc->wscr = new short[wcp.TabSize];
	if (wdbf) {
	    size_t	sz = isaa? pwc->glen: (pwc->glen / 2 + pwc->ChrNo);
	    wdbf->prepare(entlen + pwc->ChrNo, pwc->ChrNo, sz + pwc->ChrNo + 1, pfqnum);
	}
	blkscrtab(pwc->ChrNo, pwc->glen / pwc->ChrNo);

// second phase

	DbsRec*	pr = wdbf? wdbf->recidx: 0;
	CHAR*	ps = wdbf? wdbf->dbsseq: 0;
	char*	pe = wdbf? wdbf->entry: 0;
	int*	pg = wdbf? wdbf->gsiidx: 0;
	int*	pp = wdbf? wdbf->gsipool: 0;
	if (wdbf) vclear(pr);
	INT	m = 0;
	svr->reset();
	while ((ist = svr->nextseq(sd, 0)) != IS_END) {
	    if (ist == IS_ERR) continue;
	    if (wdbf) {
		pr->entptr = pe - wdbf->entry;
		pr->seqptr = ps - wdbf->dbsseq;
	    }
	    CHAR*	ts = sd->at(sd->len);
	    INT		n = 0;
	    int		b = 0;
	    for (CHAR* pq = sd->at(0); pq < ts; ++pq) {
		if (c2w(cvt[*pq], 3) && wdbf) {
		    ++n;
		    int	c = *pq - bias;
		    if (isaa)	*ps++ = c;
		    else if (n & 1)	b = c << 4;
		    else 	*ps++ = b + c;
		}
	    }
	    if (wdbf) {
		if (isaa)	*ps++ = SEQ_DELIM;
		else if (n & 1) *ps++ = b + SEQ_DELIM;
		else	*ps++ = ddelim;
		(pr++)->seqlen = n;
		for (char* ss = (*sd->sname)[0]; (*pe++ = *ss++); );
		if (wdbf->gsipool) {
		    *pg++ = pp - wdbf->gsipool;
		    if (sd->sigII) {
			for (int i = 0; i < sd->sigII->pfqnum; ++i) 
			    *pp++ = sd->sigII->pfq[i].pos;
		    }
		}
	    }
	    ch->registBlk(++m);
	    reset();
	}
	if (wdbf && pg) *pg = pp - wdbf->gsipool;
}

// read from MSA

void MakeBlk::idxblk(Seq* sd)
{
	bool	isaa = sd->isprotein();
	int	bias = (isaa? int(ALA): int(A)) - 1;
	char*	cvt =  isaa? acodon: nucl;

// first phase

	pwc->glen = 0;
	CHAR*	ss = sd->at(0);
	CHAR*	ts = sd->at(sd->len);
	for (int m = 0; m < sd->many; ++m) {
	    for (CHAR* ps = ss + m; ps < ts; ps += sd->many) {
		if (IsGap(*ps)) continue;
		if (c2w(cvt[*ps], 2)) ++pwc->glen;
	    }
	    ch->countBlk();
	    reset();
	}
	pwc->ChrID = 0;
	pwc->ChrNo = sd->many;
	pwc->blkp = new INT*[wcp.TabSize];
	pwc->wscr = new short[wcp.TabSize];
	if (wdbf) {
	    INT	sz = isaa? pwc->glen: (pwc->glen / 2 + pwc->ChrNo);
	    INT	si = sd->sigII? sd->sigII->lstnum: 0;
	    wdbf->prepare(*sd->sname, pwc->ChrNo, sz + pwc->ChrNo, si);
	}
	blkscrtab(sd->many, pwc->glen / sd->many);

// second phase

	DbsRec*	pr = wdbf? wdbf->recidx: 0;
	CHAR*	ps = wdbf? wdbf->dbsseq: 0;
	char*	pe = wdbf? wdbf->entry: 0;
	for (int m = 0; m < sd->many; ) {
	    size_t n = 0;
	    int	b = 0;
	    if (wdbf) {
		pr->seqptr = ps - wdbf->dbsseq;
		pr->entptr = pe - wdbf->entry;
	    }
	    for (CHAR* pq = ss + m; pq < ts; pq += sd->many) {
		if (IsGap(*pq)) continue;
		if (c2w(cvt[*pq], 3) && wdbf) {
		    ++n;
		    int	c = *pq - bias;
		    if (isaa)	*ps++ = c;
		    else if (n & 1)	b = c << 4;
		    else 	*ps++ = b + c;
		}
	    }
	    if (wdbf) {
		if (isaa)	*ps++ = SEQ_DELIM;
		else if (n & 1) *ps++ = b + SEQ_DELIM;
		else	*ps++ = ddelim;
		while(*pe++) ;
		(pr++)->seqlen = n;
	    }
	    ch->registBlk(++m);
	    reset();
	}
	if (wdbf && sd->sigII) {
	    int*	npfq = new int[sd->many + 1];
	    vclear(npfq, sd->many + 1);
	    int*	lst = sd->sigII->lst;
	    PFQ*	lfq = sd->sigII->pfq + sd->sigII->pfqnum;
	    for (PFQ* pfq = sd->sigII->pfq; pfq < lfq; ++pfq) 
		for (int j = 0; j < pfq->num; ++j)
		    npfq[*lst++]++;
	    for (int i = 0, n = 0; i <= sd->many; ++i) {
		wdbf->gsiidx[i] = n;
		n += npfq[i];
		npfq[i] = wdbf->gsiidx[i];
	    }
	    lst = sd->sigII->lst;
	    for (PFQ* pfq = sd->sigII->pfq; pfq < lfq; ++pfq) 
		for (int j = 0; j < pfq->num; ++j) 
		    wdbf->gsipool[npfq[*lst++]++] = pfq->pos;
	    delete[] npfq;
	}
}

template <typename file_t>
int SrchBlk::read_blk_dt(file_t fd)
{
	if (fread(pwc->Nblk, sizeof(SHORT), wcp.TabSize, fd) != wcp.TabSize)
	    return (ERROR);
	SHORT** sblkp = (SHORT**) pwc->blkp;
	SHORT*  sblkb = pwc->BytBlk == 3?
	    new SHORT[pwc->WordSz]: (SHORT*) pwc->blkb;
	if (pwc->VerNo >= 25) {
	    INT*	blkidx = new INT[wcp.TabSize];
	    if (fread(blkidx, sizeof(INT), wcp.TabSize, fd) != wcp.TabSize) {
		delete[] blkidx;
		return (ERROR);
	    }
	    // Assume pwc->WordNo < UINT_MAX
	    for (INT w = 0; w < wcp.TabSize; ++w) 
		pwc->blkp[w] = blkidx[w]? pwc->blkb + blkidx[w] - 1: 0;
	    delete[] blkidx;
	} else {
	    if (fread(sblkp, sizeof(SHORT*), wcp.TabSize, fd) != wcp.TabSize)
		return (ERROR);
	    long	offset = LONG_MAX;
	    for (INT w = 0; w < wcp.TabSize; ++w) {
	        if (pwc->Nblk[w] && sblkp[w]) {
		    if (pwc->BytBlk == 3) {
			if (offset == LONG_MAX) offset = sblkp[w] - sblkb;
			sblkp[w] -= offset;
		    } else if (pwc->VerNo < 23) {
			long    idx = sblkp[w] - sblkb;
			if (offset == LONG_MAX) offset = idx;
			pwc->blkp[w] = pwc->blkb + idx - offset;
		    } else {
			if (offset == LONG_MAX) offset = pwc->blkp[w] - pwc->blkb;
			pwc->blkp[w] -= offset;
		    }
		} else {
		    if (pwc->BytBlk == 3) sblkp[w] = 0;
		    else pwc->blkp[w] = 0;
		}
	    }
	}
	if (fread(sblkb, sizeof(SHORT), pwc->WordSz, fd) != pwc->WordSz ||
	    fread(pwc->wscr, sizeof(short), wcp.TabSize, fd) != wcp.TabSize ||
 	    fread(ConvTab, sizeof(int), pwc->ConvTS, fd) != pwc->ConvTS)
		return (ERROR);
	if (pwc->BytBlk == 4) return (OK);
	if (pwc->BytBlk == 2) {
	    SHORT*	sblk = sblkb + pwc->WordSz;
	    INT*	iblk = pwc->blkb + pwc->WordNo;
	    while (sblk > sblkb) *--iblk = *--sblk;
	    return (OK);
	}
	BYTE4	b4 = {0};
	BYTE2	b2 = {0};
	INT*	iblk = pwc->blkb;
	for (INT w = 0; w < wcp.TabSize; ++w) {
	    SHORT*	sblk = sblkp[w];
	    if (!sblk) {
		pwc->blkp[w] = 0;
		continue;
	    }
	    pwc->blkp[w] = iblk;
	    for (int i = 0; i < pwc->Nblk[w]; ++i) {
		if (i % 2) {
		    b2.s = *sblk++;
		    b4.c[0] = b2.c[1];
		    b2.s = *sblk++;
		    b4.c[1] = b2.c[0];
		    b4.c[2] = b2.c[1];
		} else {
		    b2.s = *sblk++;
		    b4.c[0] = b2.c[0];
		    b4.c[1] = b2.c[1];
		    b2.s = *sblk;
		    b4.c[2] = b2.c[0];
		}
		*iblk++ = b4.i;
	    }
	}
	delete[] sblkb;
	return (OK);
}

template <typename file_t>
static void read_pwc(ContBlk* cb, file_t fd, const char* fn)
{
	size_t	fpos = ftell(fd);
	if (fread(cb, sizeof(ContBlk), 1, fd) != 1 || cb->VerNo < 24) {
	    fseek(fd, fpos, SEEK_SET);
	    ContBlk22	tmp;
	    if (fread(&tmp, sizeof(ContBlk22), 1, fd) == 1) {
		if (tmp.VerNo < 21) {
		    wcp.Nbitpat = 1;
		    wcp.Bitpat2 = 0;
		    gswap(wcp.blklen, tmp.ConvTS);
		}
		cb->ConvTS = tmp.ConvTS;
		cb->WordNo = tmp.WordNo;
		cb->WordSz = tmp.WordSz;
		cb->ChrNo = tmp.ChrNo;
		cb->glen = tmp.glen;
		cb->AvrScr = tmp.AvrScr;
		cb->MaxBlk = tmp.MaxBlk;
		cb->BytBlk = tmp.BytBlk;
		cb->VerNo = tmp.VerNo;
	    } else	fatal(emssg, fn);
	}
}

template <typename file_t>
void SrchBlk::ReadBlkInfo(file_t fd, const char* fn)
{
	INT	GivenMaxGene = wcp.MaxGene;
	pwc = new ContBlk;
	pbc = new Block2Chr;
	if (fread(&wcp, sizeof(WCPRM), 1, fd) != 1) fatal(emssg, fn);
	read_pwc(pwc, fd, fn);
	if (pwc->VerNo < 20) 
	    fatal("%s: Version %d no longer supported !\n", fn, pwc->VerNo);
	if (GivenMaxGene) wcp.MaxGene = GivenMaxGene;
	if (fread(pbc, sizeof(Block2Chr), 1, fd) != 1)
	    fatal(emssg, fn);
	if (wcp.Ktuple == wcp.BitPat) wcp.BitPat = bitmask(wcp.BitPat);
	size_t	i = pwc->ChrNo + 1;
	pwc->ChrID = (pwc->VerNo <= BsVersion)? new CHROMO[i]: 0;
	pwc->Nblk = new SHORT[wcp.TabSize];
	pwc->blkp = new INT*[wcp.TabSize];
	pwc->blkb = new INT[pwc->WordNo];
	pwc->wscr = new short[wcp.TabSize];
	ConvTab = new INT[pwc->ConvTS];
	if (pwc->VerNo < 22) {	// for back compatibility
	    CHROMO21*	chr21 = new CHROMO21[i];
	    if (fread(chr21, sizeof(CHROMO21), i, fd) != i)
		fatal(emssg, fn);
	    for (INT j = 0; j < i; ++j) {
		pwc->ChrID[j].fpos = chr21[j].fpos;
		pwc->ChrID[j].spos = chr21[j].spos;
		pwc->ChrID[j].segn = chr21[j].segn;
	    }
	    delete[] chr21;
	} else if (pwc->ChrID && (fread(pwc->ChrID, sizeof(CHROMO), i, fd) != i))
		fatal(emssg, fn);
	if (read_blk_dt(fd) == ERROR) fatal(emssg, fn);
	pwc->VerNo = BsVersion;
	fclose(fd);
}

template <typename file_t>
static void readBlkInfo (file_t fd, const char* fn, ContBlk* pwc, Block2Chr* pbc)
{
	if (fread(&wcp, sizeof(WCPRM), 1, fd) != 1) fatal(emssg, fn);
	read_pwc(pwc, fd, fn);
	if (pwc->VerNo >= 14 && fread(pbc, sizeof(Block2Chr), 1, fd) != 1)
	    fatal(emssg, fn);
	if (pwc->VerNo < 14) --pwc->ChrNo;
	size_t	chrn = pwc->ChrNo + 1;
	pwc->ChrID = (pwc->VerNo <= BsVersion)? new CHROMO[chrn]: 0;
	if (pwc->ChrID && (fread(pwc->ChrID, sizeof(CHROMO), chrn, fd) != chrn))
	    fatal(emssg, fn);
	fclose(fd);
}

void ReportBlkInfo(const char* fn)
{
	ContBlk	wc;
	Block2Chr  bc;
const	char*	path = getenv(ALN_DBS);
	if (!path) path = DBS_DIR;
	if (is_gz(fn)) {
#if USE_ZLIB
	    gzFile	gzfd = gzopen(fn, "rb");
	    if (!gzfd &&
	        !(gzfd = gzopenpbe(path, fn, NGZ_EXT, "rb", -1)) &&
		!(gzfd = gzopenpbe(path, fn, PGZ_EXT, "rb", -1)) &&
		!(gzfd = gzopenpbe(path, fn, AGZ_EXT, "rb", -1))) 
		    fatal("Can't read %s !\n", fn);
	    readBlkInfo(gzfd, fn, &wc, &bc);
#else
	    fatal(gz_unsupport, fn);
#endif
	} else {
	    FILE*	fd = fopen(fn, "rb");
	    if (!fd &&
	  	!(fd = fopenpbe(path, fn, BKN_EXT, "rb", -1)) &&
		!(fd = fopenpbe(path, fn, BKP_EXT, "rb", -1)) &&
		!(fd = fopenpbe(path, fn, BKA_EXT, "rb", -1))) 
		    fatal("Can't read %s !\n", fn);
	    readBlkInfo(fd, fn, &wc, &bc);
	}
	size_t	chrn = wc.ChrNo + 1;
	int	MaxBlock = wc.VerNo < 21? wcp.MaxGene / wcp.blklen: wcp.Nbitpat;
	const	char*	ms = wc.VerNo < 21? "MaxBlock": "BitpatNo";
	printf("Versioin %d, Elem %u, Tuple %u, TabSize %u, BlockSize %u, MaxGene %u, %s %u, afact %u\n",
	    wc.VerNo, wcp.Nalpha, wcp.Ktuple, wcp.TabSize, wcp.blklen, wcp.MaxGene, ms, MaxBlock, wcp.afact);
	printf("GenomeSize %lu, ChrNo %lu, BlockNo %lu, WordNo %lu, WordSz %lu B, AvrScr %u\n",
	    (LONG) wc.glen, (LONG) wc.ChrNo, (LONG) (wc.ChrID? wc.ChrID[wc.ChrNo].segn: chrn), 
	    (LONG) wc.WordNo, (LONG) wc.WordSz, wc.AvrScr);
	delete[] wc.ChrID;
	exit (0);
}

void SrchBlk::init2(Seq* a)
{
	vclear(bh2->rscr, nseg);
	bh2->prqueue_b->reset();
	CHAR*	ss = a->at(a->left);
	CHAR*	ts = a->at(a->right - (bpp[0]->width + wcp.Nshift) + 1);
	INT	q = (ts - ss) % wcp.Nshift;
	MinGeneLen = poslmt = 0;
	for (INT p = 0; p < wcp.Nshift; ++p) {
	    bh2->as[0][p] = ss++;
	    bh2->as[1][q] = ts++;
	    if (++q == wcp.Nshift) q = 0;
	}
}

void SrchBlk::init4(Seq* a)
{
	MinGeneLen = bbt * (a->right - a->left) * 3 / 4;
	poslmt = (short) (gratio * a->len / wcp.blklen + 2);
	vclear(bh4->sigw, Ncand);
	vclear(bh4->testword, 4);
	vclear(bh4->nhit, 4);
	vclear(bh4->maxs, 4);
	vclear(bh4->mmct, 4);
	vclear(bh4->sign, 4);
	vclear(bh4->bscr[0], 4 * nseg);
	vclear(bh4->ascr[0], 4 * nseg);
	for (int d = 0; d < 4; ++d) {
	    bh4->prqueue_a[d]->reset();
	    bh4->prqueue_b[d]->reset();
	}
	CHAR*	ss = a->at(a->left);
	CHAR*	ts = a->at(a->right - (bpp[0]->width + wcp.Nshift));	// last cycle
	INT	q = (ts-- - ss) % wcp.Nshift;	// phase
	for (INT p = 0; p < wcp.Nshift; ++p) {
	    bh4->as[0][p] = bh4->as[2][p] = ss++;
	    bh4->as[1][q] = bh4->as[3][q] = ts++;
	    if (++q == wcp.Nshift) q = 0;
	}
}

void SrchBlk::setSegLen()
{
	INT	blk = 0;

	SegLen = new INT[nseg + 1];
	SegLen[0] = 0;	// dummy
	for (INT chrn = 0; chrn < pwc->ChrNo; ++chrn) {
	    int	len = chrsize(chrn);
	    for (; len > (int) wcp.blklen; len -= wcp.blklen)
		SegLen[++blk] = wcp.blklen;
	    SegLen[++blk] = len;
	}
}

int SrchBlk::findChrNo(INT blk)
{
	int	lw = (int) (pbc->BClw + pbc->BCce * (blk - 1)) - 1;
	int	up = (int) (pbc->BCup + pbc->BCce * (blk - 1)) + 1;

	if (lw < 0) lw = 0;
	if (up > (int) pwc->ChrNo) up = pwc->ChrNo;
	if (chrblk(lw) > blk) lw = 0;
	if (chrblk(up) < blk) up = pwc->ChrNo;
	while (up - lw > 1) {
	    int	md = (lw + up) / 2;
	    if (chrblk(md) > blk) up = md;
	    else if (chrblk(md + 1) > blk) return (md);
	    else 	lw = md;
	}
	if (chrblk(up) > blk) return (lw);
	else	return (up);
}

Seq* SrchBlk::setgnmrng(BPAIR* wrkbp, SeqList* sl)
{
	Seq**	gener = sl->seqs + 1;
	Seq**	tgtgr = sl->curgr;
	Seq**	wrkgr = tgtgr;
	INEX	g_inex = (*gener)->inex;
	INT	lb = wrkbp->lb;
	int	c = findChrNo(lb);
	INT	u = chrblk(c + 1) - 1;
	INT	d = wrkbp->d? 3: 0;
	INT	e = wrkbp->d? 2: 1;
	char	str[MAXL];

	INT	rb = wrkbp->rb;
	INT	x = (chrblk(c))? lb - chrblk(c): 0;
	INT	y = (chrblk(c))? rb - chrblk(c): 0;
	while (x && bh4->bscr[d][--lb] && !bh4->bscr[e][lb]) --x;
	while (rb < u && bh4->bscr[e][++rb] && !bh4->bscr[d][rb]) ++y;
	if (rb <= u) ++y;
	x *= wcp.blklen;
	y *= wcp.blklen;
	if (x > blkmergin) x -= blkmergin;
	else	x = 0;
	y += blkmergin;
	INT	l = y;
	INT	r = x;
	while (--wrkgr >= gener) {
	    if ((*wrkgr)->did == c && (*wrkgr)->inex.sens == d) {
		INT	wl = (*wrkgr)->SiteNo((*wrkgr)->left);
		INT	wr = (*wrkgr)->SiteNo((*wrkgr)->right);
		if (wl > wr) {u = wl; wl = wr; wr = u;}
		if (wl > x && wl < l) l = wl;	// left most
		if (wr < y && wr > r) r = wr;	// right most
	    }
	}
	if (l > r)	{x = r; y = l;}		// x=r < l < y || x < r < l=y
	else if (x + y > l + r)	x = r;		// x < l < r << y
	else			y = l;		// x << l < r < y
	if (x + MinGeneLen > y) return (0);
	sprintf(str, "Dbs%d %d %d %c", c, x + 1, y, dir[wrkbp->d / 2]);
	Seq*	sd = (*tgtgr)->getdbseq(dbf, str, c);
	if (sd) {
	    if (g_inex.molc == TRON && sd->isdrna())
		sd->nuc2tron();
	    g_inex.sens = sd->inex.sens;
	    sd->inex = g_inex;
	}
	return (sd);
}

int SrchBlk::MinQuery() {
	return wcp.Ktuple * (pwd->DvsP? 3: 5) + wcp.Nshift;
}

int SrchBlk::MaxGene() {
	return wcp.MaxGene;
}

void SrchBlk::setaaseq(Seq* a, int chn)
{
	a->getdbseq(dbf, "", chn);
}

static int cmpf(BLKTYPE* a, BLKTYPE* b)
{
	return (*a - *b);
}

static int scmpf(KVpair<INT, int>* a, KVpair<INT, int>* b)
{
	return (b->val - a->val);
}

static int gcmpf(GeneRng* a, GeneRng* b)
{
	VTYPE	d = b->scr - a->scr;

	if (d > 0) return 1;
	if (d < 0) return -1;
	return 0;
}

Randbs::Randbs(double avr, bool gdb) : trans_form(gdb? log: sqrt)
{
	if (gdb) {
	    if (RbsFact == DQUERY) RbsFact = RbsFactLog;
	} else {
	    if (RbsFact == DQUERY) RbsFact = algmode.crs? RbsFactSqrX: RbsFactSqr;
	}
	RbsCoef = RbsFact * avr;
	RbsCons = ((!gdb && algmode.crs)? RbsBaseX: RbsBase) * avr;
	RlnCoef = RlnFact * avr * wcp.Nbitpat;
	Phase1T = (int) (RbsBias * avr);
	for (INT i = 0; i < NRTAB; ++i) 
	    rscrtab[i] = (int) (RbsCoef * (*trans_form)((double) (i + 1)) + RbsCons);
}

int Randbs::randbs(INT mmc)
{
	if (mmc < NRTAB) return (rscrtab[mmc]);
	if (RbsCoef == 0) return ((int) RbsCons);
	return ((int) (RbsCoef * (*trans_form)((double) (mmc + 1)) + RbsCons));
}

int genemergin(int apos, int mingap, Seq* sd, bool rend)
{
	int	agap = rend? (sd->right - apos): (apos - sd->left);
	if (!agap) return (0);
	int	play = IntronPrm.llmt;
	if (sd->sigII) {
	    int	tge = sd->sigII->to_gene_end(apos, rend);
	    return (tge + play);
	}
	if (agap < mingap && algmode.crs) return (2 * agap + play);
	return (IntronPrm.maxl + play);
}

bool SrchBlk::grngoverlap(GeneRng* a, GeneRng* b)
{
	int	al = 0;
	int	bl = 0;
	int	ol = 0;
	JUXT*	aj = a->jxt;
	JUXT*	bj = b->jxt;
	JUXT*	at = aj + a->num;
	JUXT*	bt = bj + b->num;
	int	az = aj->jy + aj->jlen * bbt;
	int	bz = bj->jy + bj->jlen * bbt;

	while (aj < at && bj < bt) {
	    int	d = min(az, bz) - max(aj->jy, bj->jy);
	    if (d > 0) ol += d;
	    if (az <= bz) {
		al += az - aj->jy;
		++aj;
		az = aj->jy + aj->jlen * bbt;
	    }
	    if (bz <= az) {
		bl += bz - bj->jy;
		++bj;
		bz = bj->jy + bj->jlen * bbt;
	    }
	}
	return (ol > min(al, bl) * AllowedOverlapFact);
}

SrchBlk::SrchBlk(Seq* seqs[], const char* fn, bool gdb)
	: gnmdb(gdb), pwc(0), pbc(0), bh2(0),bh4(0), master(0), 
	  pwd(0), dbf(dbs_dt[0]), rdbt(0), bpp(0)
{
	if (is_gz(fn)) {
#if USE_ZLIB
	    gzFile	gzfd = gzopen(fn, "rb");
	    if (!gzfd) fatal(not_found, fn);
	    ReadBlkInfo(gzfd, fn);
#else
	    fatal(gz_unsupport, fn);
#endif
	} else {
	    FILE*	fd = fopen(fn, "r");
	    if (!fd) fatal(not_found, fn);
	    ReadBlkInfo(fd, fn);
	}
	initialize(seqs, fn);
}

SrchBlk::SrchBlk(Seq* seqs[], MakeBlk* mb, bool gdb) 
	: gnmdb(gdb), master(0), pwd(0), rdbt(0), bpp(0)
{
	pwc = new ContBlk;
	pbc = new Block2Chr;
	*pwc = *mb->pwc;
	*pbc = *mb->pbc;
	dbf = mb->wdbf; mb->wdbf = 0;
	if (dbf->curdb->defmolc != PROTEIN) dbf->curdb->defmolc = DNA;
	ConvTab = mb->iConvTab;
	mb->iConvTab = 0;
	mb->pwc->clean();
	initialize(seqs);
}

static INT max_intron_len(const char* fn)
{
	FILE*	fd = ftable.fopen(ipstat, "r");
	if (!fd) return (0);
	char	str[MAXL];
	const	char*	sl = strrchr(fn, '/');
	if (sl) fn = sl + 1;
	while (fgets(str, MAXL, fd)) {
	    if (strncmp(str, fn, 8)) continue;
	    char*	ps = cdr(str);
	    for (int i = 1; i < 5; ++i)
		ps = cdr(ps);
	    return ((INT) atoi(ps));
	}
	return (0);
}

void SrchBlk::initialize(Seq* seqs[], const char* fn)
{
	if (!pwc->blkp)
	    fatal("%s: Make and specify database !\n", fn);
	pwd = new PwdB(seqs);
	DRNA = !pwd->DvsP;
	if (DRNA ^ (wcp.Nalpha == 4))
	    fatal("%s: Block table is incompatible with query !\n", fn);
	if (!pwc->AvrScr)
	    fatal("%s: Block table may be destroyed !\n", fn);
	kk = wcp.Nbitpat / 2 + 1;
	bpp = new Bitpat*[kk];
	if (wcp.Nbitpat == 1) 
	    bpp[0] = new Bitpat(wcp.Nalpha, wcp.BitPat);
	else {
	    bpp[0] = new Bitpat(wcp.Nalpha, bitmask(wcp.Ktuple));
	    bpp[1] = new Bitpat(wcp.Nalpha, wcp.BitPat);
	}
	if (wcp.Nbitpat > 3)
	    bpp[2] = new Bitpat(wcp.Nalpha, wcp.Bitpat2);
	rdbt = new Randbs((double) pwc->AvrScr * bpp[0]->weight / wcp.Nshift, gnmdb);
#if TESTRAN
	tstrn = 0;
	MaxMmc = 0;
#endif
	maxmmc = (MaxMmc == 0 || (int) MaxMmc > (INT_MAX / bpp[0]->weight)
		 || algmode.lcl & 16)?
	    INT_MAX: bpp[0]->weight * MaxMmc / wcp.Nshift;
	vthr = (VTYPE) (alprm.scale * alprm.thr);
	ptpl = rdbt->base() / pwc->AvrScr;
	DeltaPhase2 = (int) (ptpl * pwc->AvrScr);
	MaxBlock = wcp.MaxGene / wcp.blklen;
	ExtBlock = 0;
	if (fn) ExtBlock = max_intron_len(fn) / wcp.blklen;
	if (!ExtBlock) ExtBlock = MaxBlock / 10;
	nseg = chrblk(pwc->ChrNo);
	bbt = pwd->DvsP == 1? 3: 1;
	Ncand = OutPrm.MaxOut + NCAND;
	setSegLen();
	if (gnmdb) {
	    bh4 = new Bhit4(nseg); bh2 = 0;
	} else {
	    bh2 = new Bhit2(nseg); bh4 = 0;
	}
}

SrchBlk::SrchBlk(SrchBlk* sbk, DbsDt* df)
{
	*this = *sbk;
	master = sbk;
	if (sbk->bh2) bh2 = new Bhit2(sbk->nseg);
	if (sbk->bh4) bh4 = new Bhit4(sbk->nseg);
	reset(df);
}

#if TESTRAN
Testran::Testran()
{
	trnbr = 0;
	vclear(trcnt, TESTRAN);
	vclear(trmax, TESTRAN);
	vclear(trmin, TESTRAN);
	vclear(trwrd, TESTRAN);
	vclear(travr, TESTRAN);
	vclear(trvar, TESTRAN);
}

void Testran::out(int avrscr)
{
	if (!trnbr) return;
	FILE*	fd = fopen("RandBscr", "w");
	if (!fd) fatal("Can't write RandBscr");
	fprintf(fd, "# Nmmc  Count   BSmax     BSmin   BSavr    BSsd  UAvr=%d\n", avrscr);
	for (int i = 0; i < trnbr; ++i) {
	    travr[i] /= trcnt[i];
	    trvar[i] -= trcnt[i] * travr[i] * travr[i];
	    trwrd[i] /= trcnt[i];
	    fprintf(fd, "%5d\t%5d\t%5d\t%7d\t%7.1f\t%7.1f\t%7d\n", i+1, trcnt[i], trmax[i],
		trmin[i], travr[i], sqrt(trvar[i] / (trcnt[i] - 1)), trwrd[i]);
	}
	fclose(fd);
}

void Testran::add(const Testran* tr)
{
	if (tr->trnbr > trnbr) trnbr = tr->trnbr;
	for (int i = 0; i < tr->trnbr; ++i) {
	    trcnt[i] += tr->trcnt[i];
	    trmax[i] += tr->trmax[i];
	    trmin[i] += tr->trmin[i];
	    trwrd[i] += tr->trwrd[i];
	    travr[i] += tr->travr[i];
	    trvar[i] += tr->trvar[i];
	}
}
#endif

SrchBlk::~SrchBlk()
{
	delete bh2; delete bh4;
	if (dbf != dbs_dt[0] && dbf != dbs_dt[1]) delete dbf;
	if (master) {
#if TESTRAN
	    master->add(tstrn);
	    delete tstrn;
#endif
	} else {
#if TESTRAN
	    tstrn->out(pwc->AvrScr);
	    delete tstrn;
#endif
	    delete pwc; delete pbc;
	    delete pwd; delete rdbt;
	    delete[] ConvTab; delete[] SegLen;
	    for (int k = 0; k < kk; ++k) delete bpp[k];
	    delete[] bpp;
	}
}

void set_max_extend_gene_rng(int n, bool forced)
{
	if (forced || gene_rng_max_extend < 0)
	    gene_rng_max_extend = n;
}

int get_max_extend_gene_rng()
{
	return (gene_rng_max_extend);
}

bool extend_gene_rng(Seq* sqs[], PwdB* pwd, DbsDt* dbf)
{
	if (gene_rng_max_extend <= 0) return (false);
	Seq*&   a = sqs[0];
	Seq*&   b = sqs[1];
	RANGE   prv = {a->right, a->left};
	RANGE   grng = {b->left, b->right};
	int     bbt = a->isprotein()? 3: 1;
	WLPRM*  wlprm = setwlprm(0);
	bool    rvs = b->inex.sens & 1;
	bool	extended = false;

	for (int ntry = 0; ntry < gene_rng_max_extend; ++ntry) {
	    JUXT*       lend = b->jxt;
	    JUXT*       wjxt = b->jxt + b->CdsNo - 1;
	    JUXT	rend = {wjxt->jx + wjxt->jlen, wjxt->jy + bbt * wjxt->jlen};
	    int lrextend = 0;
	    if (lend->jx < prv.left) {
		prv.left = lend->jx;
		int extend = genemergin(lend->jx, wlprm->width, a, false);
		if (lend->jy < extend) {
		    lrextend = 1;
		    if (rvs)    grng.right += extend;
		    else	grng.left -= extend;
		}
	    }
	    if (rend.jx > prv.right) {
		prv.right = rend.jx;
		int extend = genemergin(rend.jx, wlprm->width, a, true);
		if (rend.jy + extend > b->len) {
		    lrextend |= 2;
		    if (rvs)    grng.left -= extend;
		    else	grng.right += extend;
		}
	    }
	    if (lrextend) extended = true;
	    else	break;

	    char	str[MAXL];
	    if (rvs)
		sprintf(str, "$%s %d %d <", (*b->sname)[0],     
		    b->SiteNo(grng.right), b->SiteNo(grng.left));
	    else
		sprintf(str, "$%s %d %d", (*b->sname)[0],    
		    b->SiteNo(grng.left), b->SiteNo(grng.right));
	    Seq*    tmp = new Seq(1, abs(grng.right - grng.left));
	    gswap(b, tmp);
	    b->getdbseq(dbf, str);
	    if (pwd->DvsP == 1) b->nuc2tron();
	    Wilip       wl(sqs, pwd, 0);
	    WLUNIT* wlu = wl.begin();
	    if (wlu) {
		b->jxt = new JUXT[wlu->num + 1];
		vcopy(b->jxt, wlu->jxt, wlu->num + 1);
		b->CdsNo = wlu->num;
		delete tmp;
	    } else {
		gswap(b, tmp);
		delete tmp;
		break;
	    }
	}
	return (extended);
}

VTYPE SrchBlk::FindHsp(BPAIR* wrkbp, VTYPE maxjscr, SeqList* sl)
{
	Seq*	a = sl->seqs[0];
	Seq**	gener = sl->seqs + 1;
	int	c1 = findChrNo(wrkbp->lb);
	INT	zl = chrblk(c1);
	INT	zr = chrblk(c1+1) - 1;
	int	sr = chrsize(c1);
	RANGE	rng;
	RANGE	prv = {a->right, a->left};
	int	ntry = 0;
	WLPRM*	wlprm = setwlprm(algmode.crs);

	wrkbp->jscr = 0;
	Seq*   cursd = 0;
retry:
	cursd = setgnmrng(wrkbp, sl);
	if (!cursd) return (0);
	cursd->saverange(&rng);
	int	rvs = (*sl->curgr)->inex.sens & REVERS;
	if (sl->curgr > gener) swapseq(sl->curgr, gener);	// save
	if (a->isprotein()) sl->seqs[1]->nuc2tron();
	Wilip*	wl = new Wilip(sl->seqs, pwd, algmode.crs);
	WLUNIT*	wlu = wl->begin();
	if (sl->curgr > gener) swapseq(sl->curgr, gener);	// restore
	if (!wlu) {
	    delete wl;
	    return (0);
	}
	int	n = min((int) OutPrm.MaxOut, wl->size());
	JUXT	lend = {a->right, 0};
	JUXT	rend = {a->left, 0};
	for ( ; n--; ++wlu) {
	    JUXT*	wjxt = wlu->jxt;
	    if (lend.jx > wjxt->jx) lend = *wjxt;
	    wjxt += wlu->num - 1;
	    int	l = wjxt->jx + wjxt->jlen;
	    if (rend.jx < l) {
		rend.jx = l;
		rend.jy = wjxt->jy + bbt * wjxt->jlen;
	    }
	}
	int	lrextend = 0;				// try extended blocks
	if (lend.jx < prv.left && ((rvs && wrkbp->rb < zr) || (!rvs && wrkbp->lb > zl))) {
	    prv.left = lend.jx;
	    if (lend.jx > (int) wlprm->width) {
		lrextend = 1;
		if (rvs) {
		    INT	p = min(wrkbp->rb + ExtBlock, zr);
		    while (++wrkbp->rb < p)
			if (bh4->bscr[wrkbp->d][wrkbp->rb]) break;
		} else {
		    INT	p = max(wrkbp->lb - ExtBlock, zl);
		    while (--wrkbp->lb > p)
			if (bh4->bscr[wrkbp->e][wrkbp->lb]) break;
		}
	    }
	}
	if (rend.jx > prv.right && ((rvs && wrkbp->lb > zl) || (!rvs && wrkbp->rb < zr))) {
	    prv.right = rend.jx;
	    if (a->right - rend.jx > (int) wlprm->width) {
		lrextend |= 2;
		if (rvs) {
		    INT	p = max(wrkbp->lb - ExtBlock, zl);
		    while (--wrkbp->lb > p)
			if (bh4->bscr[wrkbp->e][wrkbp->lb]) break;
		} else {
		    INT p = min(wrkbp->rb + ExtBlock, zr);
		    while (++wrkbp->rb < p)
			if (bh4->bscr[wrkbp->d][wrkbp->rb]) break;
		}
	    }
	}
	if (lrextend && ntry++ < 2) {
	    delete wl;
	    goto retry;
	}
	wlu  = wl->begin();
	wrkbp->jscr = wlu->scr;
	Seq**	wrkgr;
	for ( ; wlu->num; ++wlu) {
	    JUXT*	wjxt = wlu->jxt;
	    if (wlu->scr > maxjscr) maxjscr = wlu->scr;
	    else if (wrkbp - bh4->bpair >= (int) OutPrm.MaxOut && 
		wlu->scr + vthr < maxjscr)
		break;
	    if (wlu > wl->begin()) {
		if (cursd != *sl->curgr) cursd->aliaseq(*sl->curgr);
		(*sl->curgr)->restrange(&rng);
	    }
	    (*sl->curgr)->jscr = (VTYPE) wlu->scr;
	    (*sl->curgr)->CdsNo = 0;
	    (*sl->curgr)->left = wlu->llmt;
	    (*sl->curgr)->right = wlu->ulmt;
	    int	y = a->right - a->left;
	    if (wlu->tlen > y) {
		double	x = (*sl->curgr)->jscr;
		(*sl->curgr)->jscr = int(x * y / wlu->tlen);
	    }
	    int	cl = (*sl->curgr)->SiteNo(wjxt->jy);
	    wjxt += wlu->num - 1;
	    int	cr = (*sl->curgr)->SiteNo(wjxt->jy + wjxt->jlen);
	    if (rvs) gswap(cl, cr);
	    int	tl = 0;
	    int	tr = sr;
	    for (wrkgr = sl->curgr; --wrkgr >= gener; ) {	// find nearest
		if (!strcmp((*(*wrkgr)->sname)[0], (*(*sl->curgr)->sname)[0]) &&
		    (*wrkgr)->inex.sens == (*sl->curgr)->inex.sens) {
		    wjxt = (*wrkgr)->jxt;
		    int	wl = (*wrkgr)->SiteNo(wjxt->jy);
		    wjxt += (*wrkgr)->CdsNo - 1;
		    int	wr = (*wrkgr)->SiteNo(wjxt->jy + wjxt->jlen);
		    if (rvs) gswap(wl, wr);
		    if (cr > wl && cl < wr) break;		// overlap
		    if (cr < wl && wl < tr) tr = wl;
		    if (wr < cl && wr > tl) tl = wr;
		}
	    }
	    if (wrkgr >= gener) continue;			// overlap
	    for (wrkgr = sl->curgr; --wrkgr >= gener; ) {	// sort on score
		if (wrkgr[1]->jscr > (*wrkgr)->jscr) {
		    swapseq(wrkgr, wrkgr + 1);
		    if (wrkgr[0]->vrtl && wrkgr[0]->sid == wrkgr[1]->sid) {
			gswap(wrkgr[0]->vrtl, wrkgr[1]->vrtl);
		    }
		}
		else break;
	    }
	    if (++wrkgr >= sl->lstgr) break;			// no more locus
	    if (sl->curgr < sl->lstgr) ++sl->curgr;		// last locus
	    delete[] (*wrkgr)->jxt;
	    (*wrkgr)->jxt = new JUXT[wlu->num + 1];
	    vcopy((*wrkgr)->jxt, wlu->jxt, wlu->num + 1);
	    (*wrkgr)->CdsNo = wlu->num;
	}
	delete wl;
	if (sl->curgr == sl->lstgr) (*sl->curgr)->refresh();
	return (maxjscr);
}

int SrchBlk::TestOutput(Seq* seqs[], int force)
{
	Seq*	a = seqs[0];
	Seq**	gener = seqs + 1;
	SHORT	d, e, i, j, k, phase1;
	INT	x, y, z;
	int	c1, c2;
	BLKTYPE	p, q;
	int	ReportAln = algmode.nsa != MAP1_FORM && algmode.nsa != MAP2_FORM;
	BPAIR*	wrkbp = bh4->bpair;
	BPAIR*	sigbp = bh4->bpair;
	BPAIR*	curbp = bh4->bpair;
	BPAIR*	lstbp = bh4->bpair + Ncand;
	Seq**	wrkgr;
	SHORT	sigm[4];
	SeqList sl;
static	char	ofmt[] = "%-7s %c %5d %5d %-7s %5d %5d %6.2f %6.2f %3d %3d %2d %2d %2d %2d %d %d\n";
static	char	ofmt2[] = "%-7s %c %5d %5d %-7s %5d %5d %6.2f %6.2f %3d %3d %2d %2d %2d %2d\n";

	sl.seqs = seqs;
	sl.lstgr = gener + OutPrm.MaxOut;
	sl.curgr = gener;
	curbp->bscr = 0;
	for (wrkgr = gener; wrkgr < sl.lstgr; ++wrkgr) (*wrkgr)->len = 0;
	for (d = 0; d < 4; ++d)		// sort on position
	    sigm[d] = bh4->extract_to_work(d);
	for (d = 0; d < 4; d += 2) {	// pair end
	    e = d + 1;
	    i = j = 0;
	    int	bscrThr = rdbt->randbs(bh4->mmct[d] + bh4->mmct[e]) + rdbt->Phase1T;
	    while (i < sigm[d] && j < sigm[e]) {
		p = bh4->sigb[d][i];
		q = bh4->sigb[e][j];
		c1 = findChrNo(p);
		c2 = findChrNo(q);
		int	s = q - p;
		if (c1 != c2) s *= MaxBlock;
		if (d) {
		    s = -s;
		    if (s > MaxBlock) {
			++j;
			continue;
		    } else if (s < 0) {
			++i;
			continue;
		    }
		} else {
		    if (s > MaxBlock) {
			++i;
			continue;
		    } else if (s < 0) {
			++j;
			continue;
		    }
		}
		curbp->bscr = bh4->bscr[d][p] + bh4->bscr[e][q];
		curbp->jscr = 0;
		curbp->chr = c1;
		if (s > poslmt) {
		    s -= poslmt;
		    curbp->bscr -= (int) (s * s);
		}
		if (curbp->bscr < bscrThr) {
		    if (p < q)	++i;
		    else	++j;
		    continue;
		}
		curbp->d = d;
		curbp->e = e;
		curbp->lb = d? q: p;
		curbp->rb = d? p: q;
		for (wrkbp = curbp; --wrkbp >= bh4->bpair; ) {
		    if (wrkbp[1].bscr > wrkbp->bscr) {
			*lstbp = wrkbp[1];
			wrkbp[1] = *wrkbp;
			*wrkbp = *lstbp;
		    } else break;	// sort in decending order
		}
		if (curbp < lstbp) {
		    ++curbp;
		    bh4->sigb[d][i] = bh4->sigb[e][j] = Undef;
		}
		++i; ++j;
	    }
	}
	for (d = 0; d < 4; ++d) {	// single end
	    if (d % 2)	{j = e = d - 1; k = d;}
	    else	{k = e = d + 1; j = d;}
	    for (i = 0; i < sigm[d]; ++i) {
		p = bh4->sigb[d][i];
		if (p == Undef) continue;
		c1 = findChrNo(p);
		curbp->chr = c1;
		curbp->bscr = bh4->bscr[d][p];
		curbp->jscr = 0;
		curbp->d = curbp->e = d;
		curbp->lb = curbp->rb = p;
		curbp->bscr += bh4->bscr[e][p];
		for (q = p; q && --q >= chrblk(c1); ) {
		    if (bh4->bscr[j][q]) {
			curbp->bscr += bh4->bscr[j][q];
			curbp->lb = q;
		    } else if (j != d && bh4->bscr[d][q]) {
			curbp->bscr += bh4->bscr[d][q];
			curbp->lb = q;
		    } else break;
		}
		for (q = p; ++q < chrblk(c1+1); ) {
		    if (bh4->bscr[k][q]) {
			curbp->bscr += bh4->bscr[k][q];
			curbp->rb = q;
		    } else if (k != d && bh4->bscr[d][q]) {
			curbp->bscr += bh4->bscr[d][q];
			curbp->rb = q;
		    } else break;
		}
		for (wrkbp = curbp; --wrkbp >= bh4->bpair; ) {
		    if (wrkbp[1].bscr > wrkbp->bscr) {
			gswap(wrkbp[0], wrkbp[1]);
		    } else break;	// sort in decending order
		}
		if (curbp < lstbp) ++curbp;
	    }
	}
	if (curbp == bh4->bpair) return (force? ERROR: 0);
	++force;
	lstbp = curbp;		// merge overlaps
	for (wrkbp = curbp = bh4->bpair; ++wrkbp < lstbp; ) {
	    d = wrkbp->d / 2;
	    for (sigbp = bh4->bpair; sigbp < wrkbp; ++sigbp) {
		if (sigbp->d / 2 == d &&
		    ((sigbp->lb <= wrkbp->rb) && (wrkbp->lb <= sigbp->rb))) {
			if (sigbp->lb > wrkbp->lb) sigbp->lb = wrkbp->lb;
			if (sigbp->rb < wrkbp->rb) sigbp->rb = wrkbp->rb;
			wrkbp->bscr = 0;	// overlap
		}
	    }
	}
	phase1 = 0;	// look for best block pair
	for (wrkbp = curbp = bh4->bpair; wrkbp < lstbp; ++wrkbp) {
	    if (wrkbp->bscr == 0) continue;
	    d = wrkbp->d; e = wrkbp->e;
	    if (force != 2 && (wrkbp->bscr < 
		rdbt->randbs(bh4->mmct[d] + bh4->mmct[e]) + rdbt->Phase1T))
		continue;
	    if (FindHsp(wrkbp, curbp->jscr, &sl)) phase1++;
	    if (wrkbp->jscr > curbp->jscr) curbp = wrkbp;
	}
	if (phase1) {			// phase1 passed
#if PRUNE == 2
	    VTYPE	scr = (*gener)->jscr - DeltaPhase2;
	    for (wrkgr = gener; wrkgr < sl.curgr; ++wrkgr) {  // close relatives
		if ((*wrkgr)->jscr < scr) break;
	    }
	    for (sl.curgr = wrkgr; wrkgr <= sl.lstgr && (*wrkgr)->many; ++wrkgr)
		(*wrkgr)->refresh((*wrkgr)->many);    // clean other candidates
#endif
	    if (ReportAln) return (sl.curgr - gener);
	} else if (force < 2) return (0);	// next reccurence
	else if (ReportAln && !DNA && algmode.slv) {	// examine all positive blocks
	    sl.curgr = gener;
	    curbp = bh4->bpair + 1;
	    curbp->bscr = 0;
	    for (d = 0; d < 4; d += 2) {
		bh4->bpair->d = d;
		bh4->bpair->e = e = d + 1;
		bh4->bpair->bscr = 0;
		bh4->bpair->chr = 0;
		SHORT	s = pwc->AvrScr;
		for (y = z = x = c2 = 0; ++x < nseg; ) {
		    if (bh4->bscr[d][x] || bh4->bscr[e][x]) {
			c1 = findChrNo(x);
			if (((!y && z) || c1 != c2) && bh4->bpair->bscr >= s) {
			    FindHsp(bh4->bpair, 0, &sl);
			    if (bh4->bpair->bscr > curbp->bscr) *curbp = *bh4->bpair;
			}
			if (!y ||  c1 != c2) {
			    bh4->bpair->lb  = x;
			    bh4->bpair->bscr = 0;
			    bh4->bpair->chr = c2 = c1;
			}
			if (bh4->bscr[d][x]) ++y; else y = 0;
			if (bh4->bscr[e][x]) ++z; else z = 0;
			bh4->bpair->bscr += bh4->bscr[d][x] + bh4->bscr[e][x];
			bh4->bpair->rb = x;
		    } else if (y || z) {
			if (bh4->bpair->bscr >= s) {
			    FindHsp(bh4->bpair, 0, &sl);
			    if (bh4->bpair->bscr > curbp->bscr) *curbp = *bh4->bpair;
			}
			y = z = bh4->bpair->bscr = 0;
		    }
		}
	    }
	    if (sl.curgr > gener) return (sl.curgr - gener);
	}
	d = curbp->d;				// map only or phase1 failed
	e = curbp->e;
	p = d? curbp->rb: curbp->lb;
	q = d? curbp->lb: curbp->rb;
	char*	cid = dbf->entname(curbp->chr);
	if (algmode.nsa == MAP1_FORM) {
	    printf(ofmt, a->sqname(), dir[d / 2], a->left, a->right, cid, p, q,
		(double) bh4->bscr[d][p] / TFACTOR,
		(double) bh4->bscr[e][q] / TFACTOR,
		bh4->as[d][bh4->maxs[d]] - a->at(0),
		bh4->as[e][bh4->maxs[e]] - a->at(0), 
		bh4->mmct[d], bh4->mmct[e], bh4->nhit[d], bh4->nhit[e], bh4->testword[d], bh4->testword[e]);
	} else if (algmode.nsa == MAP2_FORM) {
	    x = p - chrblk(curbp->chr);
	    y = q - chrblk(curbp->chr);
	    if (x > 0) x = (x - 1) * wcp.blklen;
	    y *= wcp.blklen;
	    printf("%s\t%s\t%7d\t%7d %c\t%7d\t%7d\n",
		a->sqname(), cid, x + 1, y, dir[d / 2], y - x, a->len);
	} else { 
	    prompt(ofmt2, a->sqname(), dir[d / 2], a->left, a->right,
		cid, p, q,
		(double) bh4->bscr[d][p] / TFACTOR,
		(double) bh4->bscr[e][q] / TFACTOR,
		bh4->as[d][bh4->maxs[d]] - a->at(0),
		bh4->as[e][bh4->maxs[e]] - a->at(0),
		bh4->mmct[d], bh4->mmct[e], bh4->nhit[d], bh4->nhit[e]);
	}
	return (ERROR);
}

Bhit4::Bhit4(int nseg)
{
	if (Nascr > Ncand) Nascr = Ncand;
	if (Nascr < 1) Nascr = 1;
	bscr[0] = new int[4 * nseg];
	ascr[0] = new int[4 * nseg];
	as[0] = new CHAR*[4 * wcp.Nshift];
	sigw = new BLKTYPE[Ncand];
	sigb[0] = new BLKTYPE[4 * Ncand];
	for (int d = 1; d < 4; ++d) {
	    bscr[d] = bscr[d-1] + nseg;
	    ascr[d] = ascr[d-1] + nseg;
	    as[d] = as[d-1] + wcp.Nshift;
	    sigb[d] = sigb[d-1] + Ncand;
	}
	bpair = new BPAIR[Ncand + 1];
	INT	na = Nascr + 1;
	INT	nc = Ncand + 1;
	blkscr = new BlkScr[4 * (nc + na)];
	vclear(blkscr, 4 * (nc + na));
	BlkScr*	bs = blkscr;
	for (int d = 0; d < 4; ++d, bs += nc) {
	    prqueue_a[d] = new PrQueue<BlkScr>(bs, Nascr, 0, false, true);
	    prqueue_b[d] = new PrQueue<BlkScr>(bs += na, Ncand, 0, false, true);
	}
}

Bhit4::~Bhit4()
{
	delete[] bscr[0];
	delete[] ascr[0];
	delete[] sigw;
	delete[] sigb[0];
	delete[] bpair;
	delete[] as[0];
	delete[] blkscr;
	for (int d = 0; d < 4; ++d) {
	    delete prqueue_a[d];
	    delete prqueue_b[d];
	}
}

void Bhit4::update_a(int d, BLKTYPE blk)
{
	BlkScr	sb = {blk, ascr[d][blk]};
	int	p = prqueue_a[d]->find(sb);
	prqueue_a[d]->put(sb, p);
}

int Bhit4::update_b(int d, BLKTYPE blk)
{
	BlkScr	sb = {blk, bscr[d][blk]};
	int	p = prqueue_b[d]->find(sb);
	prqueue_b[d]->put(sb, p);
	return (sign[d] = prqueue_b[d]->size());
}

SHORT Bhit4::extract_to_work(int d)
{
	for (SHORT i = 0; i < sign[d]; ++i) {
	    BlkScr&	bs = (*prqueue_b[d])[i];
	    sigw[i] = bs.key;
	}
	if (sign[d] < 2) {
	    sigb[d][0] = sigw[0];
	    return (sign[d]);
	}
	qsort((UPTR) sigw, (INT) sign[d], sizeof(BLKTYPE), (CMPF) cmpf);
	SHORT	j = 0;
	for (SHORT i = 0; i < sign[d]; ++i) {
	    bool	s = i && (sigw[i] > sigw[i-1] + 1);
	    if (s) ++j;					// discontinuous
	    if (d == 1 || d == 2) {
		sigb[d][j] = sigw[i];		// last block
	    } else {
		if (!i || s) sigb[d][j] = sigw[i];	// first block
	    }
	}
	return (++j);
}

Qwords::Qwords(int k, int nc, INT* ct, ContBlk* pwc, Bitpat** bp, Seq* a) :
	kk(k), DRNA(nc), ConvTab(ct), wc(pwc), bpp(bp)
{
	ww = new INT[kk];
	xx = new int[kk];
	front = new BLKTYPE[kk];
	endss = new CHAR*[kk];
	if (a) reset(a);
}

void Qwords::reset(Seq* a)
{
	for (int k = 0; k < kk; ++k)
	    endss[k] = a->at(a->right - bpp[k]->width);
}

int Qwords::querywords(CHAR* ss)
{
	if (kk == 1) {
	    int	i = *ww = *xx = 0;
	    for ( ; i < bpp[0]->weight; ++i) {
		INT	c = ConvTab[ss[bpp[0]->exam? bpp[0]->exam[i]: i]];
		if (c >= wcp.Nalpha) break;
		*ww = *ww * wcp.Nalpha + c;
	    }
	    if (!wc->blkp[*ww]) return (0);		// ubiquitous
	    if (wc->wscr[*ww] < 0) *xx = -1;		// no hit
	    if (i == bpp[0]->weight) return (wc->wscr[*ww]);
	} else {
	    INT*	w = ww;
	    int*	x = xx;
	    for (int k = 0; k < kk; ++k, ++w, ++x) {	// multiple bit patterns
		if (ss >= endss[k]) break;
		*w = *x = 0;
		int*	exam = bpp[k]->exam;
		int	i = 0;
		for ( ; i < bpp[k]->weight; ++i) {
		    INT	c = ConvTab[ss[exam? exam[i]: i]];
		    if (c >= wcp.Nalpha) break;		// ambiguous
		    *w = *w * wcp.Nalpha + c;
		}
		if (i < bpp[k]->weight) *x = -1;	// contain amb
	    }
	    int	c = 0;
	    int	wdscr = 0;
	    w = ww;
	    x = xx;
	    for (int k = 0; k < kk; ++k, ++w, ++x) {
		if (wc->wscr[*w] < 0) *x = -1;		// absent
		else if (!*x && wc->blkp[*w]) {		// not amb or ubiquitous
		    ++c;
		    wdscr += wc->wscr[*w];
		}
	    }
	    if (c) return (wdscr / c);			// find hit
	}
	return (-1);					// no hit or bad query
}

int Qwords::querywords(CHAR* ss, int d, bool rvs)
{
	if (kk == 1) {
	    int*	exam = bpp[0]->exam;
	    if (exam && rvs) exam += bpp[0]->weight;
	    int	i = *ww = *xx = 0;
	    for ( ; i < bpp[0]->weight; ++i) {
		INT	c = ConvTab[ss[exam? exam[i]: i]];
		if (c >= wcp.Nalpha) break;
		if (DRNA) {
		    if (d >= 2)	*ww = (*ww >> 2) + ((3 - c) << bpp[0]->wshift);
		    else	*ww = (*ww << 2) + c;
		} else {	
		    if (d >= 2)	*ww = (wcp.TabSize * c + *ww) / wcp.Nalpha;
		    else	*ww = *ww * wcp.Nalpha + c;
		}
	    }
	    if (!wc->blkp[*ww]) return (0);		// ubiquitous
	    if (wc->wscr[*ww] < 0) *xx = -1;		// no hit
	    if (i == bpp[0]->weight) return (wc->wscr[*ww]);
	} else {
	    int	k = 0;
	    INT*	w = ww;
	    int*	x = xx;
	    for ( ; k < kk; ++k, ++w, ++x) {		// multiple bit patterns
		if (ss >= endss[k]) break;
		*w = *x = 0;
		int	i = 0;
		int*	exam = bpp[k]->exam;
		if (exam && rvs) exam += bpp[k]->weight;
		for ( ; i < bpp[k]->weight; ++i) {
		    INT	c = ConvTab[ss[exam? exam[i]: i]];
		    if (c >= wcp.Nalpha) break;		// ambiguous
		    if (DRNA) {
			if (d >= 2)	*w = (*w >> 2) + ((3 - c) << bpp[k]->wshift);
			else	*w = (*w << 2) + c;
		    } else {
			if (d >= 2)	*w = (wcp.TabSize * c + *w) / wcp.Nalpha;
			else	*w = *w * wcp.Nalpha + c;
		    }
		}
		if (i < bpp[k]->weight) *x = -1;	// contain amb
	    }
	    int	c = 0;
	    int	wdscr = 0;
	    for (k = 0, w = ww, x = xx; k < kk; ++k, ++w, ++x) {
		if (wc->wscr[*w] < 0) *x = -1;		// absent
		else if (!*x && wc->blkp[*w]) {		// amb or ubiquitous
		    ++c;
		    wdscr += wc->wscr[*w];
		}
	    }
	    if (c) return (wdscr / c);			// find hit
	}
	return (-1);					// no hit or bad query
}

void Qwords::init_mrglist()
{
	for (int k = 0; k < kk; ++k)
	    front[k] = (xx[k] >= 0 && wc->blkp[ww[k]])? 
		wc->blkp[ww[k]][xx[k]]: 0;
}

BLKTYPE	Qwords::next_mrglist()
{
	BLKTYPE	blk = *front;
	if (kk == 1)
	    *front = (++(*xx) < wc->Nblk[*ww])? wc->blkp[*ww][*xx]: 0;
	else {
	    int*	x = xx;
	    BLKTYPE*	f = front;
	    for (int k = 1; k < kk; ++k, ++x, ++f) {
		if (*x < 0) continue;
		if (blk == 0) blk = *f;
		if (*f && *f < blk) blk = *f;
	    }
	    if (blk == 0) return (0);
	    x = xx;
	    f = front;
	    INT*	w = ww;
	    for (int k = 1; k < kk; ++k, ++x, ++f, ++w) {
		if (*x < 0) continue;
		if (blk == *f) 
		    *f = (++(*x) < wc->Nblk[*w])? wc->blkp[*w][*x]: 0;
	    }
	}
	return (blk);
}

int SrchBlk::findblock(Seq* seqs[])
{
	Seq*	a = seqs[0];	// query
	init4(a);
	int	nohit = 0;
	int	sigpr = 0;
	int	c = (a->right - a->left) / (wcp.Nshift + wcp.Nshift) - 1;
	CHAR*	rms = a->at(a->right) - wcp.Nshift * ((c < ptpl)? c: ptpl);
	int	meet[2] = {0, 0};
	Chash	hh(2 * pwc->MaxBlk, pwc);
	INT	nmmc = 0;
	int	notry = 0;
	int	maxbscr[4];
	vclear(maxbscr, 4);
	Qwords	qwd(kk, DRNA, ConvTab, pwc, bpp, a);
	INT	totalsign = 0;
	while (meet[0] + meet[1] == 0)  {
	  totalsign = 0;
	  for (SHORT d = 0; d < 4; ++d) {
	    if (meet[d / 2]) continue;
	    SHORT	prty = d % 2;		// direct:reverse
	    SHORT	e = prty? d - 1: d + 1;	// tally
	    bool	rvs = d >= 2;
	    int*&	rscr = bh4->bscr[d];
	    SHORT	maxp = 0;
	    for (INT s = 0; s < wcp.Nshift; s++) {
		CHAR**	ws = bh4->as[d] + s;
		CHAR*	ms = bh4->as[e][s];
		if (!prty && ms > rms) ms = rms;
		int	cscr = 0, q = 0, p = 0;
		hh.clear();
		do {				// continueous hits
		    CHAR*	ss = *ws;
		    if (prty)	*ws -= wcp.Nshift;	// right side
		    else	*ws += wcp.Nshift;	// left side
		    if (prty ^ (ss >= ms)) {		// has scanned
			meet[d / 2]++;
			break;
		    }
		    int wdscr = qwd.querywords(ss, d, rvs);
		    if (wdscr < 0) break;		// no hits
		    bh4->testword[d] += kk;
		    if (wdscr == 0) {q = 1; continue;}	// ubiquitous
		    qwd.init_mrglist();
		    ++p; q = 0;
		    cscr += wdscr;
		    while (BLKTYPE blk = qwd.next_mrglist()) {
			KVpair<INT, int>*	h = hh.incr(blk);
			bh4->ascr[d][blk] += wdscr;
			bh4->update_a(d, blk);
			if (p != h->val) {	// fail to conitue
			    h->val = 0;
			    if (prty) h = hh.incr(++blk);
			    else if (blk) h = hh.incr(--blk);
			}			// previous block?
			if (p == h->val) {	// yes let's continue
			    ++q;
			    rscr[blk] += wdscr;
			    if (rscr[blk] > maxbscr[d]) {
				maxbscr[d] = rscr[blk];
				bh4->maxs[d] = s;
			    } 
			    if (rscr[blk] >= rdbt->randbs(nmmc))
				bh4->update_b(d, blk);
			} else	h->val = 0;
		    }
		} while (q && cscr < rdbt->base());
		if (p > maxp) maxp = p;
		if (bh4->maxs[d] == s) nohit = !q;
	    }	// end of s loop
	    bh4->mmct[d] += nohit;
	    bh4->nhit[d] += maxp;
	    totalsign += bh4->sign[d];
	  }	// end of d loop
#if TESTRAN
	  if (nmmc < TESTRAN) {
	    if (nmmc > tstrn->trnbr) tstrn->trnbr = nmmc;
	    tstrn->trcnt[nmmc]++;
	    w = maxbscr[0] + maxbscr[1];
	    INT	i = maxbscr[2] + maxbscr[3];
	    if (i > w) {
		w = i;
		tstrn->trwrd[nmmc] += bh4->testword[2] + bh4->testword[3];
	    } else 
		tstrn->trwrd[nmmc] += bh4->testword[0] + bh4->testword[1];
	    if (w > tstrn->trmax[nmmc]) tstrn->trmax[nmmc] = w;
	    if (tstrn->trmin[nmmc] == 0 || w < tstrn->trmin[nmmc]) tstrn->trmin[nmmc] = w;
	    tstrn->travr[nmmc] += w;
	    tstrn->trvar[nmmc] += w * w;
	  }
#endif
	  if ((bh4->sign[0] && bh4->sign[1]) || (bh4->sign[2] && bh4->sign[3]))
		++sigpr;
	  if ((++nmmc % maxmmc == 0 && totalsign) || (sigpr > MinSigpr)) {
		if  ((c = TestOutput(seqs, 0))) return (c);		// found
		else if (++notry > MinSigpr) return (0);
	  }
	}	// end of mmc loop
// significant pair was not found
	if (!((bh4->sign[0] && bh4->sign[1]) || (bh4->sign[2] && bh4->sign[3]))) {
	    c = ERROR;
	    for (SHORT d = 0; d < 4; ++d) {
		for (int i = 0; i < Nascr; ++i) {
		    BlkScr&	bs = (*bh4->prqueue_a[d])[i];
		    if (bs.key)
			c = bh4->update_b(d, bs.key);
		}
	    }
	}
	if (c != ERROR) c = TestOutput(seqs, 1);
	return (c);
}

INT SrchBlk::bestref(Seq* seqs[], KVpair<INT, int>* sh, int n)
{
	Seq*	a = seqs[1];
	Wilip**	wl = new Wilip*[n];
	GeneRng	gr;
	Mfile	mfd(sizeof(GeneRng));

	swapseq(seqs, seqs + 1);
	for (int i = 0; i < n; ++i, ++sh) {
	    if (!sh->val) continue;
	    setaaseq(a, sh->key);
	    wl[i] = new Wilip(seqs, pwd, 0);
	    for (WLUNIT* wlu = wl[i]->begin(); wlu && wlu->num; ++wlu) {
		JUXT*	wjxt = wlu->jxt;
		gr.jxt = wjxt;
		gr.num = wlu->num;
		gr.lend = wlu->llmt;
		gr.rend = wlu->ulmt;
		gr.lgap = wjxt->jx - a->left;
		wjxt += wlu->num - 1;
		gr.rgap = a->right - wjxt->jx - wjxt->jlen;
		gr.len = gr.rend - gr.lend;
		gr.sid = sh->key;
		gr.scr = wlu->scr;
		mfd.write(&gr);
	    }
	}
	swapseq(seqs, seqs + 1);
	INT	k = (INT) mfd.size();
	INT	j = 0;
	GeneRng	*grs = (GeneRng*) mfd.flush();
	if (k) {
	    if (k > 1)					// by score
		qsort((UPTR) grs, k, sizeof(GeneRng), (CMPF) gcmpf);
	    GeneRng	*igr = grs;
	    for (INT i = 0; ++i <= k; ++igr) {
		GeneRng* jgr = grs;
		for ( ; jgr < igr; ++jgr)
		    if (grngoverlap(igr, jgr)) break;
		if (jgr < igr) continue;		// overlap
		setaaseq(seqs[++j], jgr->sid);
		seqs[j]->left = jgr->lend;
		seqs[j]->right = jgr->rend;
		delete[] seqs[j]->jxt;
		seqs[j]->jxt = new JUXT[jgr->num + 1];
		memcpy(seqs[j]->jxt, jgr->jxt, (jgr->num + 1) * sizeof(JUXT));
		seqs[j]->CdsNo = jgr->num;
		if (j == OutPrm.MaxOut) break;
	    }
	}
	delete[] grs;
	for (int i = 0; i < n; ++i) delete wl[i];
	delete[] wl;
	return (j);
}

// query: genomic sequence, database: protein sequence

Bhit2::Bhit2(int nseg)
{
	sigm = 0;
	rscr = new int[nseg];		// block score continuous hits
	blkscr = new BlkScr[Ncand + 1];
	prqueue_b = new PrQueue<BlkScr>(blkscr, Ncand, 0, false, true);
	as[0] = new CHAR*[2 * wcp.Nshift];
	as[1] = as[0] + wcp.Nshift;
}

Bhit2::~Bhit2()
{
	delete[] rscr;
	delete[] as[0];
	delete[] blkscr;
	delete	prqueue_b;
}

int Bhit2::update(BLKTYPE blk)
{
	BlkScr	sb = {blk, rscr[blk]};
	int	p = prqueue_b->find(sb);
	prqueue_b->put(sb, p);
	return (sigm = prqueue_b->size());
}

INT Bhit2::hsort()
{
	INT	w = OutPrm.MaxOut;
	if (sigm < OutPrm.MaxOut) w = sigm;
	prqueue_b->hsort();
	return (w);
}

int SrchBlk::findh(Seq* seqs[])
{
	Seq*	b = seqs[0];	// query
	Seq*	a = seqs[1];	// translated
	Chash	hh(2 * pwc->MaxBlk, pwc);
	ORF*	orf = b->getorf();
	if (!orf) return (ERROR);
	Chash	shash(2 * MaxNref, pwc);
	Qwords	qwd(kk,  DRNA, ConvTab, pwc, bpp);
	INT	norf = 0;
	int	c;
	for (ORF* wrf = orf; wrf->len; ++norf, ++wrf) {
	  INT	meet = 0;
	  b->translate(a, *wrf);
	  qwd.reset(a);
	  c = (a->right - a->left) / (wcp.Nshift + wcp.Nshift) - 1;
	  init2(a);
	  CHAR*	rms = a->at(a->right) - wcp.Nshift * ((c < ptpl)? c: ptpl);
	  INT	nmmc = 0;
	  int	sigpr = 0;
	  for ( ; nmmc < maxmmc && meet == 0; ++nmmc) {
	    for (int d = 0; d < 2; ++d) {
	      int	e = 1 - d;		// tally
	      int	maxp = 0;
	      for (INT s = 0; s < wcp.Nshift; s++) {
		CHAR**	ws = bh2->as[d] + s;
		CHAR*	ms = bh2->as[e][s];
		if (ms > rms) ms = rms;
		int	cscr = 0;
		int	p = 0;
		int	q = 0;
		hh.clear();
		do {		// continueous hits
		    CHAR*	ss = *ws;
		    if (d)	*ws -= wcp.Nshift;	// right side
		    else	*ws += wcp.Nshift;	// left side
		    if (d ^ (ss >= ms)) {		// has scanned
			meet++;
			break;
		    }
		    int wdscr = qwd.querywords(ss);
		    if (wdscr < 0) break;		// no hits
		    if (wdscr == 0) {q = 1; continue;}	// uniquitous
		    qwd.init_mrglist();
		    ++p; q = 0;
		    cscr += wdscr;
		    while (BLKTYPE blk = qwd.next_mrglist()) {
			KVpair<INT, int>*	h = hh.incr(blk);
			int	lenadjst = rdbt->lencor(SegLen[blk]);
			if (p == h->val) {	// let's continue
			    ++q;
			    bh2->rscr[blk] += wdscr;
			    if (bh2->rscr[blk] >= rdbt->randbs(nmmc) + lenadjst)
				bh2->update(blk);
			} else	h->val = 0;	// fail to conitue
		    }
		} while (q && cscr < rdbt->base());
		if (p > maxp) maxp = p;
	      }	// end of s loop
	    }	// end of d loop
	    if (bh2->sigm && ++sigpr > MinSigpr) break;
	  }	// end of mmc loop
	  hh.clear();
	  if (bh2->sigm) {		// heap sort in the order of block score
	    INT w = bh2->hsort();
	    for (INT i = 0; i < w; ++i) {
		c = findChrNo((*bh2->prqueue_b)[i].key);
		if (!shash.incr(c, bh2->rscr[(*bh2->prqueue_b)[i].key]))
		    fatal("Too many references are found!\n");
	    }
	  }
	} // end of wrf loop
	KVpair<INT, int>*	sh = shash.press(&c);
	if (c) {	// sort by scr in descending order
	    if (AllRef) {
		if (c > 1) qsort((UPTR) sh, (INT) c, sizeof(KVpair<INT, int>), (CMPF) scmpf);
		if ((INT) c > OutPrm.MaxOut) c = OutPrm.MaxOut;
		KVpair<INT, int>*	th = sh + c;
		for (c = 0; sh < th; ++sh)
		    setaaseq(seqs[++c], sh->key);
	    } else
		c = bestref(seqs, sh, c);
	} else
	    c = ERROR;
	delete[] orf;
	return (c);
}

// query and database: proteins or DNAs

int SrchBlk::finds(Seq* seqs[])
{
	Seq*	a = seqs[0];	// query
	init2(a);
	INT	testword[2] = {0, 0};
	INT	meet = 0;
	BLKTYPE	maxb = 0;
	Chash	hh(int(2 * pwc->MaxBlk), pwc);
	Qwords	qwd(kk,  DRNA, ConvTab, pwc, bpp, a);
	int	c = (a->right - a->left) / (wcp.Nshift + wcp.Nshift) - 1;
	CHAR*	rms = a->at(a->right) - wcp.Nshift * ((c < ptpl)? c: ptpl);
#if TESTRAN
	int	sigpr = 0;
#endif
	INT	nmmc = 0;
	for ( ; nmmc < maxmmc && meet == 0; ++nmmc) {
	    for (SHORT d = 0; d < 2; ++d) {
	      int	e = 1 - d;		// tally
	      int	maxp = 0;
	      for (INT s = 0; s < wcp.Nshift; s++) {
		CHAR**	ws = bh2->as[d] + s;
		CHAR*	ms = bh2->as[e][s];
		if (ms > rms) ms =  rms;
		int	cscr = 0;
		int	p = 0;
		int	q = 0;
		hh.clear();
		do {		// continueous hits
		    CHAR*	ss = *ws;
		    if (d)	*ws -= wcp.Nshift;	// right side
		    else	*ws += wcp.Nshift;	// left side
		    if (d ^ (ss >= ms)) {		// has scanned
			meet++;
			break;
		    }
		    int wdscr = qwd.querywords(ss);
		    if (wdscr < 0) break;		// no hits
		    ++testword[d];
		    if (wdscr == 0) {q = 1; continue;}	// uniquitous
		    qwd.init_mrglist();
		    ++p; q = 0;
		    cscr += wdscr;
		    while (BLKTYPE blk = qwd.next_mrglist()) {
			int	lenadjst = rdbt->lencor(SegLen[blk]);
			KVpair<INT, int>*	h = hh.incr(blk);
			if (p == h->val) {	// let's continue
			    ++q;
			    bh2->rscr[blk] += wdscr;
			    if (bh2->rscr[blk] > bh2->rscr[maxb]) maxb = blk;
			    if (bh2->rscr[blk] >= rdbt->randbs(nmmc) + lenadjst)
				bh2->update(blk);
			} else	h->val = 0;	// fail to conitue
		    }
		} while (q && cscr < rdbt->base());
		if (p > maxp) maxp = p;
	      }	// end of s loop
	   }	// end of d loop
#if TESTRAN
	   if (nmmc < TESTRAN) {
	    if (nmmc > tstrn->trnbr) tstrn->trnbr = nmmc;
	    tstrn->trcnt[nmmc]++;
	    INT	w = bh2->rscr[maxb];
	    tstrn->trwrd[nmmc] += testword[0] + testword[1];
	    if (w > tstrn->trmax[nmmc]) tstrn->trmax[nmmc] = w;
	    if (tstrn->trmin[nmmc] == 0 || w < tstrn->trmin[nmmc]) tstrn->trmin[nmmc] = w;
	    tstrn->travr[nmmc] += w;
	    tstrn->trvar[nmmc] += w * w;
	   }
	   if (bh2->sigm && ++sigpr > MinSigpr) break;
#endif
	}	// end of mmc loop
#if TESTRAN
	printf("%d\t%d\t%d\t%d\t%d\t%d\n", nmmc, testword[0] + testword[1], a->len,
		chrsize(maxb), bh2->rscr[maxb], rdbt->randbs(nmmc));
#endif
	if (bh2->sigm) {	// heap sort in the order of block score
	    INT	w = bh2->hsort();
	    Chash	shash(2 * w, pwc);
	    int	j = 1;
	    for (INT i = 0; i < w; ++i) {
		c = findChrNo((*bh2->prqueue_b)[i].key);
		KVpair<INT, int>*	h = shash.incr(c);
		if (h->val < 2) setaaseq(seqs[j++], c);
	    }
	    c = j - 1;
	} else
	    c = ERROR;
	return (c);
}

