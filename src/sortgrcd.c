/*******************************************************************************
*
*	sortgrcd version 2
*
*	Recover the output of spaln -O12
*	Filter the output by several criteria
*	Sort the records on the genomic positions
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

#include "sortgrcd.h"

#if M_THREAD

static	void*	worker_func(void* targ);
static	void	MasterWorker(GRFn* grcd, int* uppr, INT nchr);

#endif	// M_THREAD

DbsDt*&	gdbs = dbs_dt[0];
DbsDt*&	qdbs = dbs_dt[1];
extern	INT	supprime(INT n);
static	void	usage();
static	int	compf(GRFn* a, GRFn* b);
inline	int	cmpcif_lex(ChrInf* a, ChrInf* b) {return strcmp(gdbs->entname(a->key), gdbs->entname(b->key));}
inline	int	cmpcif_hit(ChrInf* a, ChrInf* b) {return (b->gcount - a->gcount);}
inline	int	cmpcif_did(ChrInf* a, ChrInf* b) {return (a->key - b->key);}
static	char*	fname(char* str, const char* head, const char* ext, int nth = 1);
static	GRFn*	findGeneEnd(int* nexon, GRFn* elocus, GRFn* llocus, int* GeneEnd);

static	INT	MaxeRcd = MAXRECD;
static	const	char	grext[] = ".grd";
static	const	char	erext[] = ".erd";
static	const	char	qrext[] = ".qrd";		// query seq name
static	const	char	GeneDelim[] = "!!!!! %d\n";
static	const	char	refmt[] = "Can't read %s\n";
static	int	OutMode = 4;
static	INT	Gff3EID = 0;
static	int	Csense = 0;
static	FiltParam	Filters[4] = {
	{INT_MAX, INT_MAX, 3, INT_MAX, INT_MAX, INT_MAX, INT_MIN, 0., 0.},
	{5, 3, 2, 10, 6, 3, 35., 75., 75.},
	{3, 2, 1, 6, 4, 2, 35., 93., 93.},
	{1, 1, 0, 2, 2, 1, 35., 97., 97.}
};
static	FiltParam	filter = 
	{INT_MAX, INT_MAX, 3, INT_MAX, INT_MAX, INT_MAX, INT_MIN, 0., 0.};
static	const	char*	canonical[3] = {"GTAG", "GCAG", "ATAC"};
static	int	printgrcd = 0;
static	InOrder	sort_by = ALPHABETIC;
static	bool	reverse = false;

static void usage()
{
	fputs("sortgrcd version 2.2: read binary grds and erds and sort them\n", stderr);
	fputs("\tin the order of chromosomal location in each direction\n", stderr);
	fputs("Usage: sortgrcd [options] *.grd\n", stderr);
	fputs("Note: version 2 supports outputs from spaln 2.1.0 or later\n", stderr);
	fputs("Note: version 2.1 supports -O3, 6, 7, and 8 options\n", stderr);
	fputs("Options:\n", stderr);
	fputs("\t-CN:\tMinimum % of coverage (0-100)\n", stderr);
	fputs("\t-EN:\tReport only the best (N=1) or all (N=2) results per gene locus (1)\n", stderr);
	fputs("\t-FN:\tFilter Level (0 -> 3: no -> stringent)\n", stderr);
	fputs("\t-GN:\tGenetic code (0: universal)\n", stderr);
	fputs("\t-HN:\tMinimum spaln score\n", stderr);
	fputs("\t-JN:\tMinimum ORF length (300)\n", stderr);
	fputs("\t-MN:\tMaximum total number of missmatches\n", stderr);
	fputs("\t-NN:\tMaximum total number of non-canonical boundaries\n", stderr);
	fputs("\t-ON:\tOutput format. 0:Gff3, 3:BED, 4:Native, 5:Intron, \n", stderr);
	fputs("\t\t\t6:cDNA, 7:translate, 8:CDS, 15:unique intron\n", stderr);
	fputs("\t-PN:\tMinimum Overall % identity (0-100)\n", stderr);
	fputs("\t-UN:\tMaximum total number of unpaired bases in gaps\n", stderr);
	fputs("\t-lN:\tNumber of residues per line for -O6 or -O7 (60)\n", stderr);
	fputs("\t-mN:\tMaximum allowed missmatches at both exon boundaries\n", stderr);
	fputs("\t-nN:\tallow non-canonical boundary? [0: no; 1: AT-AN; 2: 1bp mismatch; 3: any]\n", stderr);
	fputs("\t-uN:\tMaximum allowed unpaired bases in gaps at both exon boundaries\n", stderr);
	fputs("\t-gS:\tSpecify the .grp file name\n", stderr);
	fputs("\t-Sa:\tsort chromosomes in the alphabetical order of identifier (default)\n", stderr);
	fputs("\t-Sb:\tsort chromosomes in the order of abundance mapped on them\n", stderr);
	fputs("\t-Sc:\tsort chromosomes in the order of apparence in the genome db\n", stderr);
	fputs("\t-Sr:\tsort records mapped on minus strand in the reverse order of genomic positions\n", stderr);
	exit(1);
}

Chash::Chash(INT n)
{
	size = supprime(n);
	hash = new ChrInf[size];
	clear();
	nelm = 0;
}

ChrInf* Chash::resize(int wkey)
{
	ChrInf* hold = hash;
	ChrInf* ht = hold + size;

	size = supprime(2 * size);
	hash = new ChrInf[size];
	clear();
	nelm = 0;
	for (ChrInf* hw = hold; hw < ht; ++hw) {
	    if (hw->gcount) {
		ChrInf*	hnew = chashv(hw->key, false);
		hnew->gcount = hw->gcount;
		hnew->ecount = hw->ecount;
	    }
	}
	delete[] hold;
	return (chashv(wkey, true));
}

ChrInf* Chash::chashv(int wkey, bool incr = false)
{
	INT	v = wkey % size;
	INT	u = SecondHS - wkey % SecondHS;
	INT	v0 = v;
	ChrInf*	hv = hash + v;

	while (hv->gcount && hv->key != wkey) {
	    v = (v + u) % size;
	    if (v == v0) return(resize(wkey));
	    hv = hash + v;
	}
	if (!hv->gcount) {
	    hv->key = wkey;
	    nelm++;
	}
	if (!incr) return (hv);
	hv->gcount++;
	return (hv);
}

ChrInf* Chash::scream()
{
	ChrInf*	chrlist = new ChrInf[nelm];
	ChrInf*	hc = chrlist;
	ChrInf*	hz = hash + size;
	for (ChrInf* hv = hash; hv < hz; ++hv)
	    if (hv->gcount) *hc++ = *hv;
	return chrlist;
}

Ehash::Ehash(INT n)
{
	size = supprime(n);
	hash = new ExnInf[size];
	clear();
}

ExnInf* Ehash::ehashv(int left, int right)
{
	int	key = (left << 1) ^ right;
	int	v = key % size;
	int	u = SecondHS - key % SecondHS;
	int	v0 = v;
	ExnInf*	hv = hash + v;

	while (hv->exonID && (left != hv->left || right != hv->right)) {
	    v = (v + u) % size;
	    if (v == v0) return(resize(left, right));
	    hv = hash + v;
	}
	if (!hv->exonID) {
	    hv->left = left;
	    hv->right = right;
	    hv->exonID = ++Gff3EID;
	}
	hv->count++;
	return (hv);
}

ExnInf* Ehash::resize(int left, int right)
{
	ExnInf* hold = hash;
	ExnInf* ht = hold + size;

	size = supprime(2 * size);
	hash = new ExnInf[size];
	clear();
	for (ExnInf* hw = hold; hw < ht; ++hw) {
	    if (hw->count) {
		ExnInf*	hnew = ehashv(left, right);
		hnew->exonID = hw->exonID;
		hnew->count = hw->count;
	    }
	}
	delete[] hold;
	return (ehashv(left, right));
}

/* Write sorted records */
static	const	char	tfmt[] = "@ %s %c ( %d %d ) %s %d ( %d %d ) S: %.1f =: %.1f C: %.1f "
			"T#: %d T-: %d B#: %d B-: %d X: %d\n";

void Sortgrcd::Exonform(ExonRecord* ercd, GRFn* grfn, INT gcn)
{
static  char    intends[] = "  .  ";
static	const	char	sfmt[] = "%s\t%s\t%7.2f\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t"
			"%7d\t%7.1f\t%7d\t%7.1f\t%7.2f\t%7.2f %2d %d %d %d %s\n";
	bool	bundle = !grfn;
	if (bundle) {
	    grfn = grcd;
	    gcn = ngrcd;
	}
	INT	niso = 0;
	GeneRecord* pwrk = &grfn->gr;
	int	GeneRB = pwrk->Gend;	// gene right end
	GRFn*	glst = grfn + gcn;

	for ( ; grfn < glst; ++grfn) {
	    GeneRecord* gwrk = &grfn->gr;
	    if (gwrk->ng < 0) gwrk->ng = 0;
	    if (gwrk->Pmatch < filter.Pmatch || gwrk->Pcover < filter.Pcover || 
		gwrk->Gscore < filter.Gscore || gwrk->bmmc > filter.Bmmc ||
		gwrk->bunp > filter.Bunp || gwrk->ng > filter.ng) continue;
	    GERecN*	nwrk = nrcd + grfn->fn;
	    ExonRecord* ewrk = ercd + gwrk->Nrecord;
	    char*	Rname = (*nwrk->sname)[gwrk->Rid];
	    if (bundle) ewrk += nwrk->ecn;
	    bool chrchange = gwrk->Csense != pwrk->Csense || gwrk->Cid != pwrk->Cid;
	    if (chrchange || gwrk->Gstart > GeneRB) {
		if (OutMode < 8 && niso) fprintf(out_fd, GeneDelim, niso);
		niso = 0;
		GeneRB = gwrk->Gend;
	    } else if (gwrk->Gend > GeneRB) GeneRB = gwrk->Gend;
	    niso++;
	    pwrk = gwrk;
	    for (INT m = 0; m < gwrk->nexn; ++m, ++ewrk) {
		if (m == 0 && gwrk->nexn > 1) {
		    if (filter.ncan < 3) {
			int	i = 0;
			for ( ; i < 3; ++i)
			    if (!strncmp(ewrk[1].Iends, canonical[i], 2)) break;
			if (i == 3) continue;
		    }
		    if (ewrk->Bmmc > filter.bmmc || ewrk->Bunp > filter.bunp)
			continue;
		}
		if (m > 1 && m == gwrk->nexn - 1) {
		    if (filter.ncan < 3) {
			int	i = 0;
			for ( ; i < 3; ++i)
			    if (!strncmp(ewrk->Iends + 2, canonical[i] + 2, 2)) break;
			if (i == 3) continue;
		    }
		    if (ewrk->Bmmc > filter.bmmc || ewrk->Bunp > filter.bunp)
			continue;
		}
		intends[0] = ewrk->Iends[0]; intends[1] = ewrk->Iends[1];
		intends[3] = ewrk->Iends[2]; intends[4] = ewrk->Iends[3];
		fprintf(out_fd, sfmt, Rname, gdbs->entname(gwrk->Cid), 
		    ewrk->Pmatch, ewrk->Elen,
		    ewrk->Nmmc, ewrk->Nunp, ewrk->Rleft, ewrk->Rright, ewrk->Gleft, ewrk->Gright,
		    ewrk->Escore, ewrk->Ilen, ewrk->Iscore, ewrk->Sig3, ewrk->Sig5,
		    ewrk->Bmmc, ewrk->Bunp, ewrk->miss, ewrk->phase, intends);
	    }
	    if (OutMode >= 8) continue;
	    fprintf(out_fd, tfmt, gdbs->entname(gwrk->Cid), gwrk->Csense? '-': '+',
		gwrk->Gstart, gwrk->Gend, Rname, gwrk->Rlen, gwrk->Rstart, gwrk->Rend,
		gwrk->Gscore, gwrk->Pmatch, gwrk->Pcover, 
		gwrk->mmc, gwrk->unp, gwrk->bmmc, gwrk->bunp, gwrk->ng);
	}
	if (OutMode < 8 && niso) fprintf(out_fd, GeneDelim, niso);
}

static void cdsrng(RANGE* exon, Seq* sd)
{
	int&	n_cds = exon->left;
	RANGE*	cds = exon + 1;
	int	exonlen = 0;
	while (neorng(++exon)) {
	    exonlen += exon->right - exon->left;
	    if (exonlen >= sd->left) break;
	}
	cds->left = exon->right - exonlen + sd->left;
	cds->right = exon->right;
	n_cds = 0;
	while (neorng(++exon)) {
	    ++n_cds;
	    exonlen += exon->right - exon->left;
	    *++cds = *exon;
	    if (exonlen >= sd->right) break;
	}
	cds->right = exon->right - exonlen + sd->right;
	*++cds = endrng;
}

void Sortgrcd::print_cds(GeneRecord* gwrk, RANGE* exon, const char* Rname)
{
	Seq	gene(1);
	char	str[MAXL];
	int	l = gwrk->Gstart;
	int	r = gwrk->Gend;
	int	Glen = r - l;

	sprintf(str, "Dbs%d %d %d", gwrk->Cid, l, r);
	gene.getdbseq(gdbs, str, gwrk->Cid);
	if (gwrk->Csense) gene.comrev();
	for (RANGE* rng = fistrng(exon); neorng(rng); ++rng) {
	    rng->left = gene.SiteOn(rng->left);
	    rng->right = gene.SiteOn(rng->right) + 1;
	}
	Seq*	cds = gene.splice(0, exon);
	ORF*	orfs = (OutMode != 6)? cds->getorf(): 0;
	if (OutMode == 6 || orfs) {
	  const char*	cid = gdbs->entname(gwrk->Cid);
	  fprintf(out_fd, ">%s.%d %s %c [1:%d] ( %d - %d ) ",
	    cid, (l + r) / 2000, cid, gwrk->Csense? '-': '+', 
	    Glen, gwrk->Gstart, gwrk->Gend);
	}
	if (OutMode == 6) {
	    fprintf(out_fd, "%s %c 1:%d ( %d - %d ) N %7.2f\n",
		Rname, gwrk->Csense? '-': '+', gwrk->Rlen,
		gwrk->Rstart, gwrk->Rend, gwrk->Gscore);
	    GBcdsForm(exon, &gene, out_fd);
	    cds->typeseq(out_fd);
	} else if (orfs) {		// translated
	    cds->left = orfs->pos;
	    cds->right = orfs->pos + orfs->len;
	    fprintf(out_fd, "%s %c 1:%d ( %d - %d ) N %7.2f\n",
		Rname, gwrk->Csense? '-': '+', orfs->len,
		gwrk->Rstart, gwrk->Rend, gwrk->Gscore);
	    cdsrng(exon, cds);
	    GBcdsForm(exon, &gene, out_fd);
	    if (OutMode == 7)
		cds->transout(out_fd, STOP, 0);
	    else {
		cds->left = orfs->pos;
		cds->right = orfs->pos + orfs->len;
		cds->typeseq(out_fd);
	    }
	    delete[] orfs;
	}
	delete cds;
}

static int gposition(Seq* gene, RANGE* exon, int cdspos)
{
	int	cdslen = 0;
	while (neorng(++exon)) {
	    cdslen += exon->right - exon->left;
	    if (cdslen >= cdspos) break;
	}
	cdspos = exon->right - cdslen + cdspos;
	return gene->SiteNz(cdspos);
}

void Sortgrcd::print_bed(GeneRecord* gwrk, RANGE* exon, const char* Rname)
{
static	const	char*	BedFrom = "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%-11s\t%d\t";

	Seq	gene(1);
	char	str[MAXL];
	int	l = gwrk->Gstart;
	int	r = gwrk->Gend;
	int&	n_exon = exon->left;

	sprintf(str, "Dbs%d %d %d", gwrk->Cid, l, r);
	gene.getdbseq(gdbs, str, gwrk->Cid);
	if (gwrk->Csense) gene.comrev();
	for (RANGE* rng = fistrng(exon); neorng(rng); ++rng) {
	    rng->left = gene.SiteOn(rng->left);
	    rng->right = gene.SiteOn(rng->right) + 1;
	}
	Seq*	cds = gene.splice(0, exon);
	ORF*	orfs = cds->getorf();
	const	char*	cid = gdbs->entname(gwrk->Cid);
	const	char*	rgb = (gwrk->Csense)? RGB_BLUE: RGB_RED;
	RANGE*	exn = fistrng(exon);
	RANGE*	txn = exn + n_exon - 1;
	int	rend = exon[n_exon].right;
	int	cl = orfs? gposition(&gene, exon, orfs->pos): l;
	int	cr = orfs? gposition(&gene, exon, orfs->pos + orfs->len): r;
	if (gwrk->Csense && orfs) swap(cl, cr);
	if (gwrk->Csense) vreverse(exon + 1, n_exon);
	fprintf(out_fd, BedFrom, cid, l - 1, r, Rname,
	    int(10 * gwrk->Pcover),
	    gwrk->Csense? '-': '+', cl, cr, rgb, n_exon);
	for (exn = fistrng(exon); exn <= txn; ++exn)
	    fprintf(out_fd, "%d,", abs(exn->right - exn->left));
	putc('\t', out_fd);
	for (exn = fistrng(exon); exn <= txn; ++exn)
	    fprintf(out_fd, "%d%c", gwrk->Csense? rend - exn->right: exn->left,
		(exn < txn)? ',': '\n' );
	delete[] orfs;
	delete cds;
}

void Sortgrcd::Cdsform(ExonRecord* ercd, GRFn* grfn, INT gcn)
{
static	int	visit = 0;
static	const	char*	BedHeader = 
	"track name=Spaln description=\"%s\" useScore=1\n";

	if (OutMode == 3 && !visit++) // UCSC BED format
	    fprintf(out_fd, BedHeader, grdname);
	bool	bundle = !grfn;
	if (bundle) {
	    grfn = grcd;
	    gcn = ngrcd;
	}
	INT	niso = 0;
	GeneRecord* pwrk = &grfn->gr;
	char*	Rname = 0;
	int	GeneRB = pwrk->Gend;	// gene right end
	GRFn*	glst = grfn + gcn;
	RANGE	cds[max_no_exon] = {{0, 0}, {0, 0}};
	int&	n_cds = cds->left;
	int&	cdslen = cds->right;

	for ( ; grfn < glst; ++grfn) {
	    GeneRecord* gwrk = &grfn->gr;
	    if (gwrk->ng < 0) gwrk->ng = 0;
	    if (gwrk->Pmatch < filter.Pmatch || gwrk->Pcover < filter.Pcover || 
		gwrk->Gscore < filter.Gscore || gwrk->bmmc > filter.Bmmc ||
		gwrk->bunp > filter.Bunp || gwrk->ng > filter.ng) continue;
	    GERecN*	nwrk = nrcd + grfn->fn;
	    char*	rname = (*nwrk->sname)[gwrk->Rid];
	    ExonRecord* ewrk = ercd + gwrk->Nrecord;
	    if (bundle) ewrk += nwrk->ecn;
	    bool chrchange = gwrk->Csense != pwrk->Csense || gwrk->Cid != pwrk->Cid;
	    if (chrchange || gwrk->Gstart > GeneRB) {
		if (niso && algmode.mlt != 2) {
		    if (OutMode == 3) print_bed(pwrk, cds, Rname);
		    else	print_cds(pwrk, cds, Rname);
		}
		niso = n_cds = cdslen = cds[1].left = cds[1].right = 0;
		Rname = 0; pwrk = gwrk;
		GeneRB = gwrk->Gend;
	    } else if (gwrk->Gend > GeneRB) GeneRB = gwrk->Gend;
	    RANGE	exon[max_no_exon] = {{0, 0}, {0, 0}};
	    int&	n_exon = exon->left;
	    int&	exnlen = exon->right;;
	    INT	m = 0;
	    for ( ; m < gwrk->nexn; ++m, ++ewrk) {
		if (m == 0 && gwrk->nexn > 1) {
		    if (filter.ncan < 3) {
			int	i = 0;
			for ( ; i < 3; ++i)
			    if (!strncmp(ewrk[1].Iends, canonical[i], 2)) break;
			if (i == 3) continue;
		    }
		    if (ewrk->Bmmc > filter.bmmc || ewrk->Bunp > filter.bunp)
			continue;
		}
		if (m > 1 && m == gwrk->nexn - 1) {
		    if (filter.ncan < 3) {
			int	i = 0;
			for ( ; i < 3; ++i)
			    if (!strncmp(ewrk->Iends + 2, canonical[i] + 2, 2)) break;
			if (i == 3) continue;
		    }
		    if (ewrk->Bmmc > filter.bmmc || ewrk->Bunp > filter.bunp)
			continue;
		}
		if (++n_exon == max_no_exon)
		    fatal("Too many exons %d !\n", n_exon);
		exon[n_exon].left = ewrk->Gleft;
		exon[n_exon].right = ewrk->Gright;
		exnlen += abs(ewrk->Gright - ewrk->Gleft) + 1;
	    }
	    if (n_exon < (int) gwrk->nexn) continue;		// has defects
	    exon[n_exon + 1] = endrng; 
	    niso++;
	    if (algmode.mlt == 2) {
		if (OutMode == 3) print_bed(gwrk, exon, rname);
		else print_cds(gwrk, exon, rname);
	    } else if (exnlen > cdslen) {
		vcopy(cds, exon, n_exon + 2);
		Rname = rname;
		pwrk = gwrk;
	    }
	}
	if (niso && algmode.mlt != 2) {
	    if (OutMode == 3) print_bed(pwrk, cds, Rname);
	    else	print_cds(pwrk, cds, Rname);
	}
}

Ihash::Ihash(INT n)
{
	size = supprime(n);
	hash = new IntronInf[size];
	clear();
}

IntronInf* Ihash::ihashv(int left, int right, ExonRecord* ewrk, 
ExonRecord* fwrk, int mch, int mmc, int unp, char* intends, int rid)
{
	int	key = (left << 1) ^ right;
	int	v = key % size;
	int	u = SecondHS - key % SecondHS;
	int	v0 = v;
	IntronInf*	hv = hash + v;

	while (hv->count && (left != hv->left || right != hv->right)) {
	    v = (v + u) % size;
	    if (v == v0) return (0);
	    hv = hash + v;
	}
	if (!hv->count || mch > hv->mch) {
	    hv->left = left;
	    hv->right = right;
	    hv->Gleft = ewrk->Gleft;
	    hv->Gright = fwrk->Gright;
	    hv->Rleft = ewrk->Rleft;
	    hv->Rright = fwrk->Rright;
	    hv->mch = mch;
	    hv->mmc = mmc;
	    hv->unp = unp;
	    hv->ilen = fwrk->Ilen;
	    hv->Rid = rid;
	    strcpy(hv->intends, intends);
	}
	hv->count++;
	return (hv);
}

IntronInf* Ihash::resize(int left, int right, ExonRecord* ewrk,
ExonRecord* fwrk, int mch, int mmc, int unp, char* intends, int rid)
{
	IntronInf* hold = hash;
	IntronInf* ht = hold + size;

	size = supprime(2 * size);
	hash = new IntronInf[size];
	clear();
	for (IntronInf* hw = hold; hw < ht; ++hw) {
	    if (hw->count) {
		IntronInf*	hnew = ihashv(left, right, ewrk, fwrk,
			mch, mmc, unp, intends, rid);
		hnew->count = hw->count;
	    }
	}
	delete[] hold;
	return (ihashv(left, right, ewrk, fwrk, mch, mmc, unp, intends, rid));
}

static int icompf(IntronInf* a, IntronInf* b)
{
	return (Csense && reverse? b->left - a->left: a->left - b->left);
}

INT Ihash::sortihash()
{
	INT	i = 0;
	IntronInf*	bf = hash;
	IntronInf*      af = hash;

	for ( ; i < size; ++i, ++bf) 
	    if (bf->count) *af++ = *bf;
	i = af - hash;
	qsort((UPTR) hash, i, sizeof(IntronInf), (CMPF) icompf);
	return (i);
}

void Sortgrcd::Intronform(ExonRecord* ercd, GRFn* grfn, INT gcn)
{
static	int	visit = 0;
static	char	intends[] = "  ..  ";
static	const	char	fmt[] = "%s\t%c %9d %9d %3d %9d %9d\t%s\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t %s\n";
static	const	char	hdr[] = "# gID\tdir   Donor  Acceptor Phs     tgt_5     tgt_3\trefID\t  ref_l\t  ref_r\t  Match\tMisMach\t Unpair\tIntronL\tIntronEnd\n";
static	const	char	hdr15[] = "# gID\tdir   Donor  Acceptor Count   tgt_5     tgt_3\trefID\t  ref_l\t  ref_r\t  Match\tMisMach\t Unpair\tIntronL\tIntronEnd\n";

	bool	bundle = !grfn;
	if (bundle) {
	    grfn = grcd;
	    gcn = ngrcd;
	}
	INT	niso = 0;
	GRFn* 	lrfn = 0;
	GRFn* 	glst = grfn + gcn;
	GeneRecord* pwrk = &grfn->gr;
	int	GeneRB = pwrk->Gend;
	int	Nexon = 1;
	Ihash*	ihh = 0;

	for (int n = 0; grfn < glst; ++grfn, ++n) {
	    GeneRecord* gwrk = &grfn->gr;
	    bool	drna = gwrk->ng == 0;
	    if (gwrk->ng < 0) gwrk->ng = 0;
	    if (gwrk->Pmatch < filter.Pmatch || gwrk->Pcover < filter.Pcover || 
		gwrk->Gscore < filter.Gscore || gwrk->bmmc > filter.Bmmc ||
		gwrk->bunp > filter.Bunp || gwrk->ng > filter.ng) continue;
	    GERecN*	nwrk = nrcd + grfn->fn;
	    ExonRecord* fwrk = ercd + gwrk->Nrecord;
	    char*	Rname = (*nwrk->sname)[gwrk->Rid];
	    if (bundle) fwrk += nwrk->ecn;
	    if (gwrk->Csense != pwrk->Csense || 
		gwrk->Cid != pwrk->Cid || gwrk->Gstart > GeneRB) {
		if (OutMode < 8 && niso) fprintf(out_fd, GeneDelim, niso);
		niso = 0;
		GeneRB = gwrk->Gend;
	    } else if (gwrk->Gend > GeneRB) GeneRB = gwrk->Gend;
	    niso++;
	    if (OutMode == 15 && niso == 1) {	/* unique intron */
		lrfn = findGeneEnd(&Nexon, grfn, grfn + gcn - n, &GeneRB);
		ihh = new Ihash((INT) (Nexon * 1.2) + SecondHS);
	    }
	    pwrk = gwrk;
	    Csense = gwrk->Csense;
	    char	dir = gwrk->Csense? '-': '+';
	    for (INT m = 1; m < gwrk->nexn; ++m) {
		ExonRecord* ewrk = fwrk++;
		int	l, r;
		if (filter.ncan < 2) {	/* stringent */
		    for (l = 0; l < 3; ++l)
			if (!strncmp(fwrk->Iends, canonical[l], 4)) break;
		    if (l == 3) {	/* non-canonical */
			if (filter.ncan == 0 || strncmp(fwrk->Iends, canonical[2], 3))
			   continue;		/* AT-AN ? */
		    }
		} else if (filter.ncan == 2) {	/* allows 1 bp mismatch */
		    for (l = 0; l < 2; ++l)
			if (!strncmp(fwrk->Iends, canonical[l], 4)) break;
		    if (l == 2 && strncmp(fwrk->Iends, canonical[2], 3)) {
			for (l = r = 0; l < 4; ++l)
			    if (fwrk->Iends[l] != canonical[0][l]) ++r;
			if (r > 1) continue;
		    }
		}
		intends[0] = fwrk->Iends[0]; intends[1] = fwrk->Iends[1];
		intends[4] = fwrk->Iends[2]; intends[5] = fwrk->Iends[3];
		int	mch, mmc, unp;
		if (drna) {
		    mmc = ewrk->Bmmc + fwrk->Bmmc;
		    if (mmc > filter.bmmc) continue;
		    unp = ewrk->Bunp + fwrk->Bunp;
		    if (unp > filter.bunp) continue;
		    mch = 2 * alprm2.jneibr - mmc - unp;
		} else {
		    mmc = ewrk->Nmmc;
		    unp = ewrk->Nunp;
		    mch = ewrk->Elen - mmc - unp;
		}
		l = ewrk->Gright + (gwrk->Csense? -1: 1);
		r = fwrk->Gleft + (gwrk->Csense? 1: -1);
		if (OutMode == 15) {
		    ihh->ihashv(l, r, ewrk, fwrk, mch, mmc, unp, intends, gwrk->Rid);
		} else {
		    if (visit++ == 0) fputs(hdr, out_fd);
		    fprintf(out_fd, fmt, gdbs->entname(gwrk->Cid), dir, l, r, fwrk->phase,
			ewrk->Gleft, fwrk->Gright, Rname, ewrk->Rleft, 
			fwrk->Rright, mch, mmc, unp, fwrk->Ilen, intends);
		}
	    }
	    intends[0] = intends[1] = intends[4] = intends[5] = ' ';
	    if (OutMode == 15 && grfn == lrfn) {
		if (visit++ == 0) fputs(hdr15, out_fd);
		INT	l = ihh->sortihash();
		IntronInf*	ii = ihh->begin();
		for (INT r = 0; r < l; ++r, ++ii) {
		    fprintf(out_fd, fmt, gdbs->entname(gwrk->Cid), dir, ii->left, ii->right, ii->count,
			ii->Gleft, ii->Gright, (*nrcd->sname)[ii->Rid], ii->Rleft, ii->Rright,
			ii->mch, ii->mmc, ii->unp, ii->ilen, ii->intends);
		}
	    }
	    if (OutMode >= 8) continue;
	    fprintf(out_fd, tfmt, gdbs->entname(gwrk->Cid), dir,
		gwrk->Gstart, gwrk->Gend, Rname, gwrk->Rlen, gwrk->Rstart, gwrk->Rend,
		gwrk->Gscore, gwrk->Pmatch, gwrk->Pcover, 
		gwrk->mmc, gwrk->unp, gwrk->bmmc, gwrk->bunp, gwrk->ng);
	}
	if (OutMode < 8 && niso) fprintf(out_fd, GeneDelim, niso);
	delete ihh;
}
	
static GRFn* findGeneEnd(int* nexon, GRFn* elocus, GRFn* llocus, int* GeneEnd)
{
	GeneRecord* gwrk;
	GeneRecord* pwrk = &elocus->gr;

	*GeneEnd = pwrk->Gend;
	*nexon = pwrk->nexn;
	while (++elocus < llocus) {
	    gwrk = &elocus->gr;
	    if (gwrk->Csense != pwrk->Csense ||
		gwrk->Cid != pwrk->Cid ||
		gwrk->Gstart > *GeneEnd) break;
	    if (gwrk->Gend > *GeneEnd) *GeneEnd = gwrk->Gend;
	    *nexon += gwrk->nexn;
	    pwrk = gwrk;
	}
	return (--elocus);
}

/* Write sorted records in GFF3 format */
void Sortgrcd::Gff3form(ExonRecord* ercd, GRFn* grfn, INT gcn)
{
static	INT	Gff3GID = 0;
static	INT	Gff3MID = 0;
static	const	char*	fmt3r = "##sequence-region\t%s %d %d\n";
static	const	char*	fmt3g = "ID=gene%05d;Name=%s_%d;isoforms=%d\n";
static	const	char*	fmt3p = "%s\tALN\t%s\t%d\t%d\t%d\t%c\t%d\t";
static	const	char*	fmt3n = "%s\tALN\t%s\t%d\t%d\t%d\t%c\t.\t";
static	const	char*	fmt3m = "ID=mRNA%05d;Parent=gene%05d;Name=%s\n";
static	const	char*	fmt3c = "ID=%s%05d;Parent=mRNA%05d;Name=%s;";
static	const	char*	fmt3t = "Target=%s %d %d %c\n";
	bool	bundle = !grfn;
	if (bundle) {
	    grfn = grcd;
	    gcn = ngrcd;
	}
	int	niso = 0;	/* number of isoforms */
	GRFn*	slocus = grfn;	/* start of locus */
	GRFn*	elocus = 0;	/* end of locus */
	GRFn*	glst = grfn + gcn;
const	char*	ltype = "";
	char	mname[MAXL];
	int	rvg = 0, rvc = 0, hetero = 0, Nexon = 0;
	Ehash*	ehh = 0;

	if (!Gff3GID) fputs("##gff-version\t3\n", out_fd);
	for (int n = 0; grfn < glst; ++grfn, ++n) {
	    GeneRecord* gwrk = &grfn->gr;
	    if (gwrk->Pmatch < filter.Pmatch || gwrk->Pcover < filter.Pcover || 
		gwrk->Gscore < filter.Gscore || gwrk->bmmc > filter.Bmmc ||
		gwrk->bunp > filter.Bunp || gwrk->ng > filter.ng) continue;
	    GERecN*	nwrk = nrcd + grfn->fn;
	    ExonRecord* ewrk = ercd + gwrk->Nrecord;
	    char*	Rname = (*nwrk->sname)[gwrk->Rid];
	    if (bundle) ewrk += nwrk->ecn;
	    if (++niso == 1) {		/* first transcript */
		int	r;
		elocus = findGeneEnd(&Nexon, slocus = grfn, grfn + gcn - n, &r);
		ehh = new Ehash((INT) (Nexon * 2 + SecondHS));
		hetero = gwrk->ng;
		if (hetero)	ltype = "cds";
		else		ltype = "exon";
		rvg = gwrk->Csense? '-': '+';
		rvc = gwrk->Rsense? '-': '+';
		int	l = slocus->gr.Gstart;
		fprintf(out_fd, fmt3r, gdbs->entname(gwrk->Cid), l, r);
		fprintf(out_fd, fmt3n, gdbs->entname(gwrk->Cid), "gene", l, r, Nexon, rvg); 
		fprintf(out_fd, fmt3g, ++Gff3GID, gdbs->entname(gwrk->Cid), (l + r) / 2000,
		    elocus - slocus + 1);
	    }
	    if (grfn >= elocus) niso = 0;
	    fprintf(out_fd, fmt3n, gdbs->entname(gwrk->Cid), "mRNA", 
		gwrk->Gstart, gwrk->Gend,(int) gwrk->Gscore, rvg);
	    sprintf(mname, "%s_%d", gdbs->entname(gwrk->Cid), 
		(gwrk->Gstart + gwrk->Gend) / 2000);
	    fprintf(out_fd, fmt3m, ++Gff3MID, Gff3GID, mname);
	    for (INT m = 0; m < gwrk->nexn; ++m, ++ewrk) {
		int	l = ewrk->Gleft;
		int	r = ewrk->Gright;
		if (l > r) swap(l, r);
		if (hetero) {
		    fprintf(out_fd, fmt3p, gdbs->entname(gwrk->Cid), ltype, l, r,
			(int) ewrk->Escore, rvg, (3 - ewrk->phase) % 3);
		} else {
		    fprintf(out_fd, fmt3n, gdbs->entname(gwrk->Cid), ltype, l, r, 
			(int) ewrk->Escore, rvg);
		}
		ExnInf*	eif = ehh->ehashv(ewrk->Gleft, ewrk->Gright);
		l = ewrk->Rleft;
		r = ewrk->Rright;
		if (l > r) swap(l, r);
		fprintf(out_fd, fmt3c, ltype, eif->exonID, Gff3MID, mname);
		fprintf(out_fd, fmt3t, Rname, l, r, rvc);
	    }
	}
	delete ehh;
}

static int compf(GRFn* a, GRFn* b)
{
	int	gel = a->gr.Csense - b->gr.Csense;
	if (gel) return (gel);
	gel = a->gr.Gstart - b->gr.Gstart;
	if (gel) return (gel);
	gel = a->gr.Gend - b->gr.Gend;
	if (gel) return (gel);
	return (a->gr.nexn - b->gr.nexn);
}

static char* fname(char* str, const char* head, const char* ext, int nth)
{
	strcpy(str, head);
	char* pdot = str + strlen(str);
	while (nth && --pdot >= str)
	    if (*pdot == '.') --nth;
	if (nth) strcat(str, ext);
	else	strcpy(pdot, ext);
	return (str);
}

// read all E-records at once
ExonRecord* Sortgrcd::ReadRcd(int ac, const char** av)
{
	ercd = new ExonRecord[nercd];
	GERecN*	nwrk = nrcd;
	char	str[LINE_MAX];

	for (ExonRecord* ewrk = ercd; ac--; ++nwrk) {
	    int	nth = 1;
#if USE_ZLIB
	    if (is_gz(*av)) nth = 2;
#endif
	    FILE*	fe = fopen(fname(str, *av++, erext, nth), "r");
	    if (fe) {
		if (fread(ewrk, sizeof(ExonRecord), nwrk->ern, fe) != nwrk->ern) 
		    fatal(refmt, str);
		fclose(fe);
	    } else {
#if USE_ZLIB
		strcat(str, gz_ext);
		gzFile	gzfe = gzopen(str, "rb");
		if (!gzfe) fatal(not_found, str);
		if (fread(ewrk, sizeof(ExonRecord), nwrk->ern, gzfe) <= 0)
		    fatal(refmt, str);
		fclose(gzfe);
#else
		fatal(not_found, str);
#endif
	    }
	    ewrk += nwrk->ern;
	}
	return (ercd);
}

template <typename file_t>
int Sortgrcd::read_chr_rec(file_t fe, ExonRecord*& ewrk, GRFn* frcd, INT grn, INT fn)
{
	GRFn*	ftrm = frcd + grn;
	for (GRFn* grfn = frcd; grfn < ftrm; ++grfn) {
	    if (grfn->fn == fn) {
		GeneRecord*	gwrk = &grfn->gr;
		fseek(fe, gwrk->Nrecord * sizeof(ExonRecord), SEEK_SET);
		gwrk->Nrecord = ewrk - ercd;
		if (fread(ewrk, sizeof(ExonRecord), gwrk->nexn, fe) <= 0)
		    return (ERROR);
		ewrk += gwrk->nexn;
	    }
	}
	return (OK);
}

// read E-records with a specific Cid
ExonRecord* Sortgrcd::ReadChrRcd(int ac, const char** av, 
	INT nercd, GRFn* frcd, INT grn)
{
	ExonRecord* ewrk = ercd = new ExonRecord[nercd];
	char	str[LINE_MAX];

	for (INT fn = 0; ac--; ++fn) {
	    int	nth = 1;
#if USE_ZLIB
	    if (is_gz(av[fn])) nth = 2;
#endif
	    FILE* fe = fopen(fname(str, av[fn], erext, nth), "r");
	    if (fe) {
		read_chr_rec(fe, ewrk, frcd, grn, fn);
		fclose(fe);
	    } else {
#if USE_ZLIB
		strcat(str, gz_ext);
		gzFile	gzfe = gzopen(str, "rb");
		if (gzfe) read_chr_rec(gzfe, ewrk, frcd, grn, fn);
		else	fatal(refmt, str);
		fclose(gzfe);
#else
		fatal(refmt, str);
#endif
	    }
	}
	return (ercd);
}

#if M_THREAD

static void* worker_func(void* arg)
{
	thread_arg_t* targ = (thread_arg_t*) arg;
	int	incr[2] = {2 * (cpu_num - targ->cpuid) - 1, 2 * targ->cpuid + 1};

#ifdef __CPU_SET
	cpu_set_t	mask;
	__CPU_ZERO(&mask);
	__CPU_SET(targ->cpuid, &mask);
	if (sched_setaffinity(0, sizeof(mask), &mask) == -1)
	    prompt("Warning: faild to set CPU affinity !\n");
#endif
	int	parity = 1;
	for (INT cn = targ->cpuid; cn < targ->nchr; cn += incr[parity = 1 - parity]) {
	    qsort((UPTR) (targ->grcd + targ->occr[cn]), 
	    targ->occr[cn + 1] - targ->occr[cn], sizeof(GRFn), (CMPF) compf);
	}
	return (void*) 0;
}

static void MasterWorker(GRFn* grcd, int* occr, INT nchr)
{
	cpu_num = sysconf(_SC_NPROCESSORS_CONF);

	if (thread_num <= 0) thread_num = cpu_num;
	thread_arg_t*	targ = new thread_arg_t[thread_num];
	pthread_t*	worker = new pthread_t[thread_num];

	for (int n = 0; n < thread_num; ++n) {
	    targ[n].cpuid = n % cpu_num;
	    targ[n].grcd = grcd;
	    targ[n].occr = occr;
	    targ[n].nchr = nchr;
	    pthread_create(worker + n, 0, worker_func, (void*) (targ + n));
	}
	for (int n = 0; n < thread_num; ++n)
	    pthread_join(worker[n], 0);
	delete[] targ;
	delete[] worker;
}
#endif

// assort G-records by chromosome
void Sortgrcd::assort_by_chr(Chash* hh)
{
	int*	occr = new int[2 * nchr];
	int*	uppr =  occr + nchr;
	int	ngrd = 0;
	for (INT i = 0; i < nchr; ++i) {
	    occr[i] = ngrd;
	    uppr[i] = ngrd += chrlist[i].gcount;
	    ChrInf*	hv = hh->chashv(chrlist[i].key);
	    hv->ecount = i;
	}
	while (ngrd) {
	    for (INT i = 0; i < nchr; ++i) {
		if (occr[i] == uppr[i]) continue;
		ChrInf*	hv = hh->chashv(grcd[occr[i]].gr.Cid);
		INT	j = hv->ecount;
		if (i == j) occr[i]++;
		else {
		    GRFn	tmp = grcd[occr[i]];
		    grcd[occr[i]] = grcd[occr[j]];
		    grcd[occr[j]] = tmp;
		    occr[j]++;
		}
		if (--ngrd == 0) break;
	    }
	}

// sort G-records within each chromosome
#if M_THREAD
	uppr[-1] = 0;
	if (thread_num) MasterWorker(grcd, uppr - 1, nchr);
	else
#endif
	{
	    ngrd = 0;
	    for (INT i = 0; i < nchr; ngrd = uppr[i++])
		qsort((UPTR) (grcd + ngrd), uppr[i] - ngrd, 
			sizeof(GRFn), (CMPF) compf);
	}
	delete[] occr;
}

static int count_record(const char* bdy, const char* ext, size_t rcdsize)
{
const	char*	errmsg = "%s may be obolete or corrupted!\n";
	char	str[LINE_MAX];
	int	nth = 1;
#if USE_ZLIB
	if (is_gz(bdy)) nth = 2;
#endif	    
	FILE*	fd = fopen(fname(str, bdy, ext, nth), "r");
	if (fd) {
	    fseek(fd, 0L, SEEK_END);
	    long	fsz = ftell(fd);
	    if (fsz % rcdsize) fatal(errmsg, str);
	    fclose(fd);
	    return (fsz / rcdsize);
	}
#if USE_ZLIB
	strcat(str, gz_ext);
	gzFile	gzfd = gzopen(str, "rb");
	if (gzfd) {
	    INT	nr = 0;
	    char	buf[rcdsize];
	    while (fread(buf, rcdsize, 1, gzfd) > 0) ++nr;
	    fclose(gzfd);
	    return (nr);
	}
#endif
	fatal(errmsg, str);
	return (0);
}

// Count the numbers of records
Sortgrcd::Sortgrcd(int ac, const char** av) : argc(ac), grdname(*av)
{
	GERecN	gerNo = {0, 0, 0, 0, 0};
	Mfile	mfd(sizeof(GERecN));
	char	str[LINE_MAX];

	for ( ; ac--; ++av) {
	    gerNo.grn = count_record(*av, grext, sizeof(GeneRecord));
	    gerNo.ern = count_record(*av, erext, sizeof(ExonRecord));
	    int	nth = 1;
#if USE_ZLIB
	    if (is_gz(*av)) nth = 2;
#endif
	    FILE*	fd = fopen(fname(str, *av, qrext, nth), "r");
	    if (fd) {
		gerNo.sname = new Strlist(fd, gerNo.ern);
		fclose(fd);
	    } else {
#if USE_ZLIB
		strcat(str, gz_ext);
		gzFile	gzfd = gzopen(str, "rb");
		if (!gzfd) fatal(refmt, str);
		gerNo.sname = new Strlist(gzfd, gerNo.ern);
		fclose(gzfd);
#else
		fatal(refmt, str);
#endif
	    }
	    mfd.write(&gerNo);
	    gerNo.gcn += gerNo.grn;
	    gerNo.ecn += gerNo.ern;
	}
	nrcd = (GERecN*) mfd.flush();
	ngrcd = gerNo.gcn;
	nercd = gerNo.ecn;
	gdbs = new DbsDt((*gerNo.sname)[0]);
}

Sortgrcd::~Sortgrcd()
{
	for (int i = 0; i < argc; ++i) delete nrcd[i].sname;
	delete[] chrlist; delete[] grcd; delete[] nrcd;
}

void Sortgrcd::readGrcd(int ac, const char** av, INT hashsize)
{
// readn G-records
	Chash*	hh = new Chash(hashsize);
	grcd = new GRFn[ngrcd];
	GRFn*	grfn = grcd;
	char	str[LINE_MAX];

	for (INT fn = 0; ac--; ++fn, ++av) {
	    int	nth = 1;
#if USE_ZLIB
	    if (is_gz(*av)) nth = 2;
#endif
	    FILE* fd = fopen(fname(str, *av, grext, nth), "r");
	    if (fd) {
		while (fread(&grfn->gr, sizeof(GeneRecord), 1, fd) == 1) {
		    GeneRecord* gwrk = &grfn->gr;
		    if (gwrk->Csense) swap(gwrk->Gstart, gwrk->Gend);
		    ChrInf*	hv = hh->chashv(gwrk->Cid, true);
		    hv->ecount += gwrk->nexn;
		    (grfn++)->fn = fn;
		}
		fclose(fd);
	    } else {
#if USE_ZLIB
		strcat(str, gz_ext);
		gzFile	gzfd = gzopen(str, "rb");
		while (fread(&grfn->gr, sizeof(GeneRecord), 1, gzfd) > 0) {
		    GeneRecord* gwrk = &grfn->gr;
		    if (gwrk->Csense) swap(gwrk->Gstart, gwrk->Gend);
		    ChrInf*	hv = hh->chashv(gwrk->Cid, true);
		    hv->ecount += gwrk->nexn;
		    (grfn++)->fn = fn;
		}
		fclose(gzfd);
#else
		fatal(not_found, str);
#endif
	    }
	}

// sort G-records by chromosomal position
	chrlist = hh->scream();
	nchr = hh->nelm;
	CMPCIF	cmpcif = 0;
	switch (sort_by) {
	  case INPUT_ODR:	cmpcif = cmpcif_did; break;
	  case ALPHABETIC:	cmpcif = cmpcif_lex; break;
	  case ABUNDANCE:	cmpcif = cmpcif_hit; break;
	}
	qsort((UPTR) chrlist, nchr, sizeof(ChrInf), (CMPF) cmpcif);
	assort_by_chr(hh);
	delete hh;
}

void Sortgrcd::printGrcd()
{
	GRFn* glst = grcd + ngrcd;
	for (GRFn* wrfn = grcd; wrfn < glst; ++wrfn) {
	    GeneRecord* gwrk = &wrfn->gr;
	    GERecN*	nrec = nrcd + wrfn->fn;
	    char*	Rname = (*nrec->sname)[gwrk->Rid];
	    printf(tfmt, gdbs->entname(gwrk->Cid), gwrk->Csense? '-': '+',
		gwrk->Gstart, gwrk->Gend, Rname,
		gwrk->Rlen, gwrk->Rstart, gwrk->Rend,
		gwrk->Gscore, gwrk->Pmatch, gwrk->Pcover, 
		gwrk->mmc, gwrk->unp, gwrk->bmmc, gwrk->bunp, gwrk->ng);
	}
}

int main(int argc, const char** argv)
{
#if MONITOR
static	const	int	no_cp = 10;
	struct  tms*	tt = new tms[no_cp];
	struct  timeb*	tb = new timeb[no_cp];
	int	cp = 1;
	times(tt);
	ftime(tb);
#endif

const	char*	grpfn = 0;
	char	str[MAXL];
	INT	HashSize = 1033;

// Get options
	while (--argc && (++argv)[0][0] == OPTCHAR) {
	  const char*	val = argv[0] + 2;
	  int	flevel = 0;
	  switch (argv[0][1]) {
	    case 'F': 
		if ((val = getarg(argc, argv, true)))	// filter level
		    {flevel = atoi(val);}
		if (0 > flevel || flevel > 3) flevel = 3;
		filter = Filters[flevel];
		break;
	    case 'C': 
		if ((val = getarg(argc, argv, true)))	// minimum coverage
		    {filter.Pcover = atof(val);}
		if (filter.Pcover < 1.) filter.Pcover *= 100.;
		break;
	    case 'E':
		if ((val = getarg(argc, argv, true)))	// isoforms
		    {algmode.mlt = atoi(val);}
		break;
	    case 'G':
		if ((val = getarg(argc, argv, true)))	// genetic code
		    {initcodon(atoi(val));}
		else	{printgrcd = 1;}		// print grcd records only
		break;
	    case 'H': 
		if ((val = getarg(argc, argv, true)))	// minimum spaln score
		    {filter.Gscore = atof(val);}
		break;
	    case 'J': 
		if ((val = getarg(argc, argv, true)))	// minimum ORF length
		    {setorf(atoi(val));}
		break;
	    case 'M':
		if ((val = getarg(argc, argv, true)))	// max total number of
		    {filter.Bmmc = atoi(val);}
		break;		// mismatches near junction
	    case 'N':
		if ((val = getarg(argc, argv, true)))	// max number of non-canonical ends
		    {filter.ng = atoi(val);}
		break;
	    case 'O': 
		if ((val = getarg(argc, argv, true)))	// output format
		    {OutMode = atoi(val);}
		break;
	    case 'P':
		if ((val = getarg(argc, argv, true)))	// Minimum total % sequnce identity
		    {filter.Pmatch = atof(val);}
		if (filter.Pmatch < 1.) filter.Pmatch *= 100.;
		break;
	    case 'S': 
		if ((val = getarg(argc, argv))) {
		  switch (*val) {
		    case 'a': sort_by = ALPHABETIC; break;
		    case 'b': sort_by = ABUNDANCE; break;
		    case 'c': sort_by = INPUT_ODR; break;
		    case 'r': reverse = true; break;
		    default: break;	// no action
		  }
		}
		break;
	    case 'U':
		if ((val = getarg(argc, argv, true)))	// max total number of gaps
		    {filter.Bunp = atoi(val);}		// near junction break;
		break;
	    case 'V':
		if ((val = getarg(argc, argv))) {	// Maximum internal memory
		    MaxeRcd = atoi(val);		// size used for core sort
		    switch (val[strlen(val) - 1]) {
			case 'k': case 'K': MaxeRcd *= KILO; break;
			case 'm': case 'M': MaxeRcd *= MEGA; break;
			default: break;
		    }
		  }
		  break;
	    case 'g':
		if ((val = getarg(argc, argv)))		// .grp file name
		    {grpfn = val;}
		break;
	    case 'h':
		if ((val = getarg(argc, argv, true)))	// hash size
		    {HashSize = atoi(val);}
		break;
	    case 'l':
		if ((val = getarg(argc, argv)))		// max number of mismatches
		    {setlpw(atoi(val));}		// near each junction
		break;
	    case 'm':
		if ((val = getarg(argc, argv)))		// max number of mismatches
		    {filter.bmmc = atoi(val);}		// near each junction
		break;
	    case 'n':
		if ((val = getarg(argc, argv)))		// 0: disallow; 1: allow
		    {filter.ncan = atoi(val);}		// non-canonical ends
		else	{filter.ncan = 0;}
		break;
#if M_THREAD
	    case 't':
		if ((val = getarg(argc, argv)))		// number of thread
		    {thread_num = atoi(val);}
		else {thread_num = -1;}
		break;
#endif
	    case 'u':
		if ((val = getarg(argc, argv)))		// max total number of gaps
		    {filter.bunp = atoi(val);}		// near each junction
		break;
	    default: break;
	  }
	}
	if (!argc) usage();
	setup_output(OutMode);
	if (grpfn) {	// get seq size from *.grp 
	    FILE*	fd = fopen(grpfn, "r");
	    if (!fd) fatal("%s not found !\n", grpfn);
	    while (fgets(str, MAXL, fd)) ;
	    fclose(fd);
	    HashSize = atoi(cdr(str));	// read from the last line
	}
#if MONITOR
	times(tt + cp);
	ftime(tb + cp++);
#endif

// constructor that counts number of records
	Sortgrcd	srtgrcd(argc, argv);

#if MONITOR
	times(tt + cp);
	ftime(tb + cp++);
#endif
// read and sort G-records
	srtgrcd.readGrcd(argc, argv, HashSize);

#if MONITOR
	times(tt + cp);
	ftime(tb + cp++);
#endif
// output G-records if indicated
	if (printgrcd) {
	    srtgrcd.printGrcd();
#if MONITOR
	    goto report;
#else
	    exit(0);
#endif
	}

// output sorted records
	if (srtgrcd.nErcd() < MaxeRcd) {	// sort in core memory
#if MONITOR
	times(tt + cp);
	ftime(tb + cp++);
#endif
	    ExonRecord* ercd = srtgrcd.ReadRcd(argc, argv);
#if MONITOR
	times(tt + cp);
	ftime(tb + cp++);
#endif
	    switch (OutMode) {
	    	case 0: srtgrcd.Gff3form(ercd); break;
	    	case 5: case 13: case 15:
			srtgrcd.Intronform(ercd); break;
		case 3: case 6:	case 7: case 8:
			srtgrcd.Cdsform(ercd); break;
		case 4:	default:
			srtgrcd.Exonform(ercd); break;
	    }
	    delete[] ercd;
	} else {		// 
	    ChrInf*	hz = srtgrcd.chrend();
	    GRFn*	grfn = srtgrcd.begin();
	    INT ne = 0, ng = 0, nn = 0;
	    for (ChrInf* hv = srtgrcd.chrbegin(); hv < hz; ++hv) {
		if (ne + hv->ecount > MaxeRcd) {
		    if (!ng) {
			ne = hv->ecount;
			ng = hv->gcount;
			hv->ecount = hv->gcount = 0;
		    }
#if MONITOR
	times(tt + cp);
	ftime(tb + cp++);
#endif
		    ExonRecord* ercd = srtgrcd.ReadChrRcd(argc, argv, ne, grfn, ng);
#if MONITOR
	times(tt + cp);
	ftime(tb + cp++);
#endif
	    	    switch (OutMode) {
			case 0:
			    srtgrcd.Gff3form(ercd, grfn, ng); break;
			case 5: case 13: case 15:
			    srtgrcd.Intronform(ercd, grfn, ng); break;
			case 3: case 6: case 7: case 8:
			    srtgrcd.Cdsform(ercd, grfn, ng); break;
			case 4:	default:
			   srtgrcd. Exonform(ercd, grfn, ng); break;
		    }
		    delete[] ercd;
		    grfn += ng;
		    ne = ng = nn = 0;
		}
		ne += hv->ecount;
		ng += hv->gcount;
		++nn;
	    }
	    if (ne) {
		ExonRecord* ercd = srtgrcd.ReadChrRcd(argc, argv, ne, grfn, ng);
	    	switch (OutMode) {
		    case 0:
			srtgrcd.Gff3form(ercd, grfn, ng); break;
		    case 5: case 13: case 15:
			srtgrcd.Intronform(ercd, grfn, ng); break;
		    case 3: case 6: case 7: case 8:
			srtgrcd.Cdsform(ercd, grfn, ng); break;
		    case 4:	default:
			srtgrcd.Exonform(ercd, grfn, ng); break;
		}
		delete[] ercd;
	    }
	}
	EraDbsDt();

#if MONITOR
report:
	times(tt + cp);
	ftime(tb + cp);
	for (int i = 1; i < cp; ++i) {
	    int p1 = tb[i+1].time - tb[i].time;
	    int p2 = tb[i+1].millitm - tb[i].millitm;
	    if (p2 < 0) {
		--p1; p2 += 1000;
	    }
	    fprintf(stderr, "CPU = %ld/60s, elt = %d.%ds\n",
	        tt[i+1].tms_utime + tt[i+1].tms_stime -
		tt[i].tms_utime - tt[i].tms_stime, p1, p2);
	}
	delete[] tt;
	delete[] tb;
#endif
	return (0);
}
