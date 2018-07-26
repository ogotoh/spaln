/*****************************************************************************
*
*	Subroutines for printing aligned sequences
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

#include "divseq.h"
#include "utilseq.h"
#include "simmtx.h"

static	SeqDb*	out_form = 0;
static	SeqDb*	fst_form = 0;

static	const	char	BLANK = 0;
static	const	char	DEF_GAP_CHAR = '-';
static	const	char	GCG_GAP_CHAR = '.';
static	const	int	NO_WGHT = 5;
static	const	int	DEF_LINELENGTH = 60;
static	const	int	GCG_LINELENGTH = 50;
static	const	int	GCG_SeqBlkSz = 10;
static	const	INT	INTRONBIT = 0x80;
static	const	INT	RES_BITS = 0x1f;
static	const	char*	RGB_RED = "255,0,0";
static	const	char*	RGB_BLUE = "0,255,255";
static	const	char*	unknown_taxon = "K_Pp";

static	void	setgapchar(Seq* sd, int c);
static	void	repalninf0(SKL* skl, Seq* seqs[], Gsinfo* gsi);
static	void	repalninf1(SKL* skl, Seq* seqs[], Gsinfo* gsi);
static	void	repalninf2(SKL* skl, Seq* seqs[], Gsinfo* gsi);
static	void	repalninf3(SKL* skl, Seq* seqs[], Gsinfo* gsi);
static	void	repalninf4(SKL* skl, Seq* seqs[], Gsinfo* gsi);
static	void	repalninf6(SKL* skl, Seq* seqs[], Gsinfo* gsi);
static	int	cdsForm(RANGE* rng, Seq* sd, int mode);
static	void	put_SigII(PFQ* pfq, PFQ* tfq, int* lst, int numlst);
static	int	chempro(int* cmp, int ii);
static	int	logonuc(int* cmp, int ii, const SEQ_CODE* defcode);
static	int	csym(int cmp[], int rows, const SEQ_CODE* defcode);
static	int	SeqGCGCheckSum(char* seq, int len);
static	int	checksum(int* checks, GAPS* gaps, Seq* sd);
static	void	gcg_out(GAPS* gaps[], Seq* seqs[], int seqnum);
static	void	Gff3Form(Gsinfo* gsinf, Seq* gene, Seq* qry);
static	void	Gff3PWA(Gsinfo* gsinf, Seq* gene, Seq* qry);
static	void	BedForm(Gsinfo* gsinf, Seq* gene, Seq* qry);
static	void	ExonForm(Gsinfo* gsinf, Seq* gene, Seq* qry);
static	void	IntronForm(Gsinfo* gsinf, Seq* gene, Seq* qry);
static	void	CigarForm(Gsinfo* gsinf, Seq* gene, Seq* qry);
static	void	VulgarForm(Gsinfo* gsinf, Seq* gene, Seq* qry);
static	void	SamForm(Gsinfo* gsinf, Seq* gene, Seq* qry);
		//lpw blk Nout Ncolony eij ovl fnm rm trim lg lbl dsc odr spj olr color self
OUTPRM	OutPrm = {60, 0, 1, MAX_COLONY, 10, 5, 0, 1, 0, 3, 0, 0, 0, 0, 1, 0, 0, 0};
static	char	ditto = _IBID;
static	Row_Mode	prmode = Row_Last;
static	int	prtrc = TRON;
static	int	lLblMode = LABEL_Numb;
static	int	rLblMode = LABEL_Name;
static	int	clmark = 0;
static	char	mark[] = " .:=*";
static	int	emphsim = 0;
static	int	nmk = 4;
static	double	NGP = 0.25;
static	const	char*	moltype[] = {"TEXT", "PROTEIN", "DNA", "RNA", "TRON", "Genomic"};
static	const	char*	attrfrmt = " %c [1:%d]  ( %d - %d )";
static	const	char*	attrfrmt2 = " [%d:%d]  ( %d - %d )";
static	const	char*	attrfrmt3 = " %d:%d  ( %d - %d )";
static	const	char	grext[] = ".grd";		// gene
static	const	char	erext[] = ".erd";		// exon
static	const	char	qrext[] = ".qrd";		// query seq name
static	const	char	arext[] = ".ard";		// alignment
//	const	char	cigar[] = "MIDNSHP=XJA";	// J=5, A=3
//	const	char	vulgar[] = "MCGN53ISFJA";	// J=5, A=3

int setlpw(int lpwd)
{
	if (lpwd >= 0)
	    OutPrm.lpw = lpwd;
	else if (lpwd == QUERY)
	    promptin("Line Width (%d) : ", &OutPrm.lpw);
	return (OutPrm.lpw);
}

INT setdeflbl(int msf)
{
	if (msf == QUERY) {
	    msf = (int) OutPrm.deflbl;
	    promptin("mem label (%d) : ", &msf);
	} else if (msf >= 0)
	    OutPrm.deflbl = (INT) msf;
	return (OutPrm.deflbl);
}

FILE* setup_output(int omode, const char* def_fn)
{
	if (!(OutPrm.out_file && *OutPrm.out_file) && def_fn)
	    OutPrm.out_file = def_fn;
	if (OutPrm.out_file && *OutPrm.out_file) {
	    if (!(out_fd = fopen(OutPrm.out_file, "w")))
		fatal("Can't write to %s\n", OutPrm.out_file);
	} else
	    out_fd = stdout;
	switch (omode) {
	  case EXN_FORM: case CIG_FORM: case VLG_FORM: case SAM_FORM:
	    setform('b'); break;
	  case CDS_FORM: case AAS_FORM:
	    if (OutPrm.deflbl == 0) OutPrm.deflbl = 3;
	  case PSJ_FORM:
	    fst_form = setform('f');
	    setform('b'); break;
	}
	out_form = setform(0);
	return out_fd;
}

void close_output()
{
	if (out_fd && out_fd != stdout) {
	    fclose(out_fd);
	    out_fd = 0;
	}
}

/* 1 line seq. difference form */

static void repalninf0(SKL* skl, Seq* seqs[], Gsinfo* gsi)
{
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
	FSTAT*	fst = &gsi->fstat;
	int	rsv_mhits = 0;

	gswap(rsv_mhits, distPrm.corr_mhits);
	if (OutPrm.ColorEij) {
	    put_stat(out_fd, fst);
	    int 	apfq = a->sigII? a->sigII->pfqnum: 0;
	    int 	bpfq = b->sigII? b->sigII->pfqnum: 0;
	    int		dnm = apfq + bpfq;
	    float	ncmn = gsi->sigII? gsi->sigII->n_common(): 0;
	    float	spbdist = dnm? (1. - 2 * ncmn / dnm): 0;
	    fprintf(out_fd, "%6.2f\t%6d\t%6d\t", 100. * spbdist, apfq, bpfq);
	} else	put_stat(out_fd, gsi);
	fprintf(out_fd, "%6.1f\t%d %d %c\t%d %d %c\t%s\t%s\n",
	    (double) fst->val / (alprm.scale * a->many * b->many),
	    a->SiteLe() + 1, a->SiteRe(), a->Strand(),
	    b->SiteLe() + 1, b->SiteRe(), b->Strand(),
	    a->sqname(), b->sqname());
	gswap(rsv_mhits, distPrm.corr_mhits);
}

/* 2 lines classic SKL format */

static void repalninf1(SKL* skl, Seq* seqs[], Gsinfo* gsi)
{
	int	nn = skl->n;
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
	FSTAT&	fst = gsi->fstat;

	a->fphseq();
	putc(' ', out_fd);
	b->fphseq();
	fprintf(out_fd, "  %d  %.2f\n", nn, 
	    (double) fst.val / (alprm.scale * a->many * b->many));
	for (++skl; nn--; skl++)
	    fprintf(out_fd, "%d %d ", a->SiteNo(skl->m), b->SiteNo(skl->n));
	putc('\n', out_fd);
}

/* sugar-like format */

static void repalninf2(SKL* skl, Seq* seqs[], Gsinfo* gsi)
{
	int	nn = skl->n;
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
const	char*	anm = a->sqname();
const	char*	bnm = b->sqname();
	char	asn = a->Strand();
	char	bsn = b->Strand();
	Simmtx*	sm = getSimmtx(0);

	for (++skl; nn-- > 1; skl++) {
	    int	dm = skl[1].m - skl->m;
	    int	dn = skl[1].n - skl->n;
	    if (dm && dn) {
		CHAR*	as = a->at(skl->m);
		CHAR*	bs = b->at(skl->n);
		int	i = min(dm, dn);
		double	val = 0;
		for ( ; i--; ++as, ++bs)
		    val += sm->mtx[*as][*bs];
		val /= alprm.scale;
		fprintf(out_fd, "sugar: %s %d %d %c %s %d %d %c %d M %d %d\n", 
		    bnm, b->SiteNo(skl->n), b->SiteNo(skl[1].n), bsn,
		    anm, a->SiteNo(skl->m), a->SiteNo(skl[1].m), asn,
		    (int) val, dn, dm);
	    }
	}
}

/* psl-like format */

static void repalninf3(SKL* skl, Seq* seqs[], Gsinfo* gsi)
{
	int	nn = skl->n;
	int	mch = 0;
	int	mmc = 0;
	int	rep = 0;
	int	noN = 0;
	int	qgap = 0;
	int	qunp = 0;
	int	tgap = 0;
	int	tunp = 0;
	int	diag = 0;
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
	CHAR*	as;
	CHAR*	bs;
	Mfile	mfd(sizeof(SKLP));
	SKLP	tmp;
	SKLP*	xyl;

	for (++skl; nn-- > 1; skl++) {
	    int	dm = skl[1].m - skl->m;
	    int	dn = skl[1].n - skl->n;
	    if (!dm) {
		qgap++;
		qunp += dn;
		bs = b->at(skl->n);
		for (int i = 0; i < dn; ++i, ++bs)
		    if (*bs == AMB) noN++;
	    } else if (!dn) {
		tgap++;
		tunp += dm;
	    } else {
		++diag;
		tmp.m = a->SiteNz(skl->m);
		tmp.n = b->SiteNz(skl->n);
		tmp.p = dm;
		mfd.write(&tmp);
		as = a->at(skl->m);
		bs = b->at(skl->n);
		for (int i = 0; i < dm; ++i, ++as, ++bs) {
		    if (*as == *bs || (*as == SER && *bs == SER2)) ++mch;
		    else	++mmc;
		}
	    }
	}
	xyl = (SKLP*) mfd.flush();
	fprintf(out_fd, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c%c\t",
	    mch, mmc, rep, noN, qgap, qunp, tgap, tunp, a->Strand(), b->Strand());
	fprintf(out_fd, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t",
	    a->sqname(), a->len, a->SiteNz(a->left), a->SiteNz(a->right), 
	    b->sqname(), b->len, b->SiteNz(b->left), b->SiteNz(b->right), diag);
	for (int i = 0; i < diag; ++i)
	    fprintf(out_fd, "%ld,", xyl[i].p);
	putc('\t', out_fd);
	for (int i = 0; i < diag; ++i)
	    fprintf(out_fd, "%d,", xyl[i].m);
	putc('\t', out_fd);
	for (int i = 0; i < diag; ++i)
	    fprintf(out_fd, "%d,", xyl[i].n);
	putc('\n', out_fd);
	delete[] xyl;
}

/* 1 line compact XYL format (x, y, len) * n */

static void repalninf4(SKL* skl, Seq* seqs[], Gsinfo* gsi)
{
	int	nn = skl->n;
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
	FSTAT&	fst = gsi->fstat;

	fprintf(out_fd, "XYL: %s %d %d %c %s %d %d %c %.1lf : ",
	    a->sqname(), a->SiteLe() + 1, a->SiteRe(), a->Strand(),
	    b->sqname(), b->SiteLe() + 1, b->SiteRe(), b->Strand(),
	    (double) fst.val/ (alprm.scale) * a->many * b->many);
	for (++skl; nn-- > 1; skl++) {
	    int	dm = skl[1].m - skl->m;
	    int	dn = skl[1].n - skl->n;
	    if (dm && dn) fprintf(out_fd, " %d %d %d",
		a->SiteNo(skl->m), b->SiteNo(skl->n), dm);
	}
	fprintf(out_fd, "\n");
}

/* 2 line compact XYL format (x, y, len) * n */

static void repalninf6(SKL* skl, Seq* seqs[], Gsinfo* gsi)
{
	int	nn = skl->n;
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
	FSTAT&	fst = gsi->fstat;

	fprintf(out_fd, "XYL2: %s %d %d %c %s %d %d %c %7.1lf ",
	    a->sqname(), a->SiteLe() + 1, a->SiteRe(), a->Strand(),
	    b->sqname(), b->SiteLe() + 1, b->SiteRe(), b->Strand(),
	    (double) fst.val/ (alprm.scale) * a->many * b->many);
	fprintf(out_fd, "%6.2f %d %d %d %d %d\n",
	    100. * fst.mch / (fst.mch + fst.mmc + fst.gap),
	    (int) fst.mch, (int) fst.mmc, (int) fst.gap, (int) fst.unp,
	    nn);
	for (++skl; nn-- > 1; skl++) {
	    int	dm = skl[1].m - skl->m;
	    int	dn = skl[1].n - skl->n;
	    if (dm && dn) fprintf(out_fd, " %d %d %d",
		a->SiteNo(skl->m), b->SiteNo(skl->n), dm);
	}
	fprintf(out_fd, "\n");
}

void repalninf(Seq* seqs[], Gsinfo* gsi)
{
	if (!out_fd) return;
	switch (algmode.nsa) {
	  case 0:  repalninf0(gsi->skl, seqs, gsi); break;
	  case 1:  repalninf1(gsi->skl, seqs, gsi); break;
	  case 2:  repalninf2(gsi->skl, seqs, gsi); break;
	  case 3:  repalninf3(gsi->skl, seqs, gsi); break;
	  case 4:  repalninf4(gsi->skl, seqs, gsi); break;
	  case CIG_FORM: CigarForm(gsi, seqs[1], seqs[0]); return;
	  case VLG_FORM: VulgarForm(gsi, seqs[1], seqs[0]); return;
	  case SAM_FORM: 
	  case BAM_FORM: SamForm(gsi, seqs[1], seqs[0]); return;
	  default: repalninf6(gsi->skl, seqs, gsi); break;
	}
}

static void setgapchar(Seq* sd, int c)
{
	if (sd->isprotein())
	    amino[gap_code] = amino[nil_code] = c;
	else
	    nucl[gap_code] = nucl[nil_code] = c;
}

/*	CDS format	*/

static int cdsForm(RANGE* rng, Seq* sd, int mode)
{
	int	cmpl = sd->inex.sens == COMREV;
	char	str[MAXL];
	int	i, l, r;
	int	intv;
	int	indel = 0;
static	char	gftcds[] = "     CDS	     ";
static	char	eftcds[] = "FT   CDS	     ";
static	char	blank[]  = "		     ";
static	char	eblank[] = "FT		   ";
static	char	emblxx[] = "XX\n";
static	char	cdsmsg[] = ";C ";
static	char	modmsg[] = ";M ";
const	char*	head = "";
	int	c = 0;
	int	lb = 0;
const	char*	ModifyMssg = modmsg;

	RANGE*	wrng = rng = fistrng(rng);
	if (mode == 1) {
	  switch(out_form->FormID) {
	    case GenBank:
		fprintf(out_fd, "%s\n", out_form->ComLabel);
		if (out_form->ContSpc < strlen(blank))
		    blank[out_form->ContSpc] = '\0';
		ModifyMssg = blank;
		break;
	    case EMBL:
		fputs(emblxx, out_fd);
		ModifyMssg = "CC   ";
		break;
	  }
	} else if (mode == 2) {
	    switch(out_form->FormID) {
	      case GenBank:
		head = blank;
		fputs(out_form->FeaLabel, out_fd);
		putc('\n', out_fd);
		strcpy(str, gftcds);
		break;
	      case EMBL:
	    	fputs(emblxx, out_fd);
		head = eblank;
		strcpy(str, eftcds);
		break;
	      default:
		head = cdsmsg;
		strcpy(str, cdsmsg);
		break;
	    }
	    strcat(str, cmpl? "complement(join(": "join(");
	    c = strlen(str);
	    lb = strlen(head);
	}

	if (cmpl) {
	    while (neorng(wrng + 1)) ++wrng;
	    l = sd->SiteNo(wrng->right - 1);
	} else
	    l = sd->SiteNo(wrng->left);
	while (neorng(wrng) && wrng >= rng) {
	    if (cmpl) {
		r = sd->SiteNo(wrng->left);
		if (wrng > rng) {
		    i = sd->SiteNo(wrng[-1].right - 1);
		    intv = i - r - 1;
		    --wrng;
		} else {
		    if (mode == 2) sprintf(str + c, "%d..%d", l, r);
		    break;
		}
	    } else {
		r = sd->SiteNo(wrng->right - 1);
		if (neorng(wrng + 1)) {
		    i = sd->SiteNo(wrng[1].left);
		    intv = i - r - 1;
		    ++wrng;
		} else {
		    if (mode == 2) sprintf(str + c, "%d..%d", l, r);
		    break;
		}
	    }
	    if (intv == 0 || intv >= IntronPrm.llmt) {
		if (mode != 2) continue;
		sprintf(str + c, "%d..%d", l, r);
		l = i;
		strcat(str, ",");
		if ((c = strlen(str)) > 56) {
		    fputs(str, out_fd);
		    putc('\n', out_fd);
		    strcpy(str, head);
		    c = lb;
		}
	    } else ++indel;
	    if (mode != 1) continue;
	    if (intv <= 2 && intv > 0) {
		fputs(ModifyMssg, out_fd);
		fprintf(out_fd, "Deleted %d chars at %d\n", intv, r);
	    } else if (intv < 0) {
		fputs(ModifyMssg, out_fd);
		fprintf(out_fd, "Insert %d chars at %d\n", -intv, r);
	    }
	}
	if (mode == 2) {
	    fputs(str, out_fd);
	    putc(')', out_fd);
	    if (cmpl) putc(')', out_fd);
	    putc('\n', out_fd);
	}
	if (out_form->FormID == GenBank && mode == 1)
	    blank[out_form->ContSpc] = ' ';
	return (indel);
}

void GBcdsForm(RANGE* rng, Seq* sd)
{
	RANGE	svr;

	sd->saverange(&svr);
	if (out_form->FormID != GenBank && out_form->FormID != EMBL) {
	    sd->left = 0; sd->right = sd->len;
	}
	if (cdsForm(rng, sd, 0))
	    cdsForm(rng, sd, 1);
	cdsForm(rng, sd, 2);
	sd->restrange(&svr);
}

static	const	char*	fmt3p = "%s\tALN\t%s\t%d\t%d\t%d\t%c\t%d\t";
static	const	char*	fmt3n = "%s\tALN\t%s\t%d\t%d\t%d\t%c\t.\t";

/* Phase 1 and 2 of Gff3 are reversed from those of native format */

static void Gff3Form(Gsinfo* gsinf, Seq* gene, Seq* qry)
{
	int	rvg = gene->Strand();
	int	rvr = qry->Strand();
	int	cds = 0;
	int	hetero = qry->isprotein();
	int	len3 = (gene->len + 2) / 3 * 3;
	EISCR*	rng = gsinf->eijnc->begin();
	EISCR*	wkr = rng;
	EISCR*	skp = rng;
	SKL*	skl = gsinf->skl;
	SKL*	lst = skl + (hetero? skl->n: 0);
	SKL*	prv;
	char	mname[MAXL];
const	char*	ltype;
const	char*	fmt3;
static	int	Gff3MID = 0;
static	int	Gff3CID = 0;
static	const	char*	fmt3g = "ID=gene%05d;Name=%s\n";
static	const	char*	fmt3m = "ID=mRNA%05d;Parent=gene%05d;Name=%s\n";
static	const	char*	fmt3c = "ID=%s%05d;Parent=mRNA%05d;Name=%s;";
static	const	char*	fmt3t = "Target=%s %d %d %c\n";

	if (hetero) {
	    ltype = "cds";
	    fmt3 = fmt3p;
	} else {
	    ltype = "exon";
	    fmt3 = fmt3n;
	}
	while (neoeij(wkr + 1)) ++wkr;
	int	l = gene->SiteNo(rng->left);
	int	r = gene->SiteNo(wkr->right - 1);
	sprintf(mname, "%s_%d", (*gene->sname)[0], (l + r) / 2000);
	if (!Gff3MID++) fputs("##gff-version\t3\n", out_fd);
	l = gene->SiteNo(0);
	r = gene->SiteNo(gene->len - 1);
	if (rvg == '-') gswap(l, r);
	fprintf(out_fd, "##sequence-region\t%s %d %d\n", gene->sqname(), l, r);
	l = gene->SiteNo(rng->left);
	r = gene->SiteNo(wkr->right - 1);
	if (rvg == '-') gswap(l, r);
	fprintf(out_fd, fmt3n, (*gene->sname)[0], "gene", l, r, 
	    (int) (gsinf->scr / alprm.scale), rvg);
	fprintf(out_fd, fmt3g, Gff3MID, mname);
	fprintf(out_fd, fmt3n, (*gene->sname)[0], "mRNA", l, r, 
	    (int) (gsinf->scr / alprm.scale), rvg);
	fprintf(out_fd, fmt3m, Gff3MID, Gff3MID, mname);

	wkr = rng;
	prv = ++skl;
	while (neoeij(wkr)) {
	    if (wkr->iscr > NEVSEL) {
		l = gene->SiteNo(skp->left);
		r = gene->SiteNo(wkr->right - 1);
		if (rvg == '-') gswap(l, r);
		if (hetero) {
		    fprintf(out_fd, fmt3, (*gene->sname)[0], ltype, l, r,
			(int) (wkr->escr / alprm.scale), rvg, (len3 - cds) % 3);
		} else {
		    fprintf(out_fd, fmt3, (*gene->sname)[0], ltype, l, r, 
			(int) (wkr->escr / alprm.scale), rvg);
		}
		cds += wkr->right - skp->left;
		fprintf(out_fd, fmt3c, ltype, ++Gff3CID, Gff3MID, mname);
		fprintf(out_fd, fmt3t, (*qry->sname)[0], qry->SiteNo(wkr->rleft),
		    qry->SiteNo(wkr->rright - 1), rvr);
		while (skl <= lst && skl->m <= wkr->right + 1) {
		    l = skl->m - prv->m;	/* genome */
		    r = skl->n - prv->n;	/* protein */
		    if (l && (r || l <= IntronPrm.llmt)) cds -= l % 3;
		    prv = skl++;
		}
		skp = ++wkr;
	    } else {
		skp = wkr++;
	    }
	}
}

static void Gff3PWA(Gsinfo* gsinf, Seq* gene, Seq* qry)
{
	EISCR*	rng = gsinf->eijnc->begin();
	EISCR*	wkr = rng;
	EISCR*	skp = rng;
	SKL*	skl = gsinf->skl;
	SKL*	lst = skl + skl->n;
	SKL*	prv = ++skl;
	int	rvg = gene->Strand();
	int	rvr = qry->Strand();
	int	hetero = qry->isprotein();
	char	mname[MAXL];
	int	donor = 0;
const	char*	ltype;
static	int	Gff3MID = 0;
static	const	char*	fmt3c = "ID=match%05d;Name=%s;";
static	const	char*	fmt3t = "Target=%s %d %d %c;Gap=";

	swapskl(gsinf->skl);
	if (hetero)	ltype = "nucleotide_to_protein_match";
	else		ltype = "cDNA_match";
	while (neoeij(wkr + 1)) ++wkr;
	int	l = gene->SiteNo(rng->left);
	int	r = gene->SiteNo(wkr->right - 1);
	sprintf(mname, "%s_%d", (*gene->sname)[0], (l + r) / 2000);
	if (!Gff3MID++) fputs("##gff-version\t3\n", out_fd);
	l = gene->SiteNo(0);
	r = gene->SiteNo(gene->len - 1);
	if (rvg == '-') gswap(l, r);
	fprintf(out_fd, "##sequence-region\t%s %d %d\n", gene->sqname(), l, r);

	wkr = rng;
	while (neoeij(wkr)) {
	    if (wkr->iscr > NEVSEL) {
		l = gene->SiteNo(skp->left);
		r = gene->SiteNo(wkr->right - 1);
		if (rvg == '-') gswap(l, r);
		fprintf(out_fd, fmt3n, (*gene->sname)[0], ltype, l, r, 
		    (int) (wkr->escr / alprm.scale), rvg);
		fprintf(out_fd, fmt3c, Gff3MID, mname);
		fprintf(out_fd, fmt3t, (*qry->sname)[0], qry->SiteNo(wkr->rleft),
		    qry->SiteNo(wkr->rright - 1), rvr);
		for ( ; skl <= lst && skl->n <= wkr->right + 1; prv = skl++) {
		    int	dm = skl->m - prv->m;
		    int	dn = skl->n - prv->n;
		    r = skp->left - skl->n;	/* acceptor */
		    if (!dm && donor && -1 <= r && r <= 1)
			continue;		 /* intron */
		    l = wkr->right - skl->n;
		    donor = -1 <= l && l <= 1;
		    if (hetero) {
			if (!dn) {
			    if (dm) fprintf(out_fd, "I%d ", dm);
			    continue;
			}
			l = dn - dm * 3;
			if (l < 0) {
			    if (l > -3) dn += 2;
			    dm = dn / 3;
			}
			if (dm)		fprintf(out_fd, "M%d ", dm);
			if (l < 0) {
			    l = -l;
			    if (l % 3)	fprintf(out_fd, "R%d ", l);
			    else	fprintf(out_fd, "I%d ", l / 3);
			} else if (l)	{
			    if (l % 3)	fprintf(out_fd, "F%d ", l);
			    else	fprintf(out_fd, "D%d ", l / 3);
			}
		    } else {
			if (!dm && !dn)	continue;
			else if (!dn)	fprintf(out_fd, "I%d ", dm);
			else if (!dm )	fprintf(out_fd, "D%d ", dn);
			else		fprintf(out_fd, "M%d ", dm);
		    }
		}
		putc('\n', out_fd);
		skp = ++wkr;
	    } else {
		skp = wkr++;
	    }
	}
	swapskl(gsinf->skl);
}

#if PSL_FORM

/* psl-like format */

static void PslForm(Gsinfo* gsinf, Seq* t, Seq* q)
{
	int	mch = 0;
	int	mmc = 0;
	int	rep = 0;
	int	noN = 0;
	int	qgap = 0;
	int	qunp = 0;
	int	tgap = 0;
	int	tunp = 0;
	int	diag = 0;
	int	hetero = q->isprotein();
	EISCR*	fst = gsinf->eijnc->bigin();
	EISCR*	wkr = fst;
	SKL*	skl = gsinf->skl;
	SKL*	lst = skl + skl->n;
	SKL*	psk = ++skl;
	Mfile	mfd(sizeof(SKLP));
	SKLP	tmp;
	SKLP*	xyl;
static	int	visit = 0;
static	const	char*	PslHeader[3] = {
"match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q		Q   	Q    	Q  	T		T   	T    	T  	block	blockSizes 	qStarts	 tStarts\n",
"     	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	start	end	count\n",
"---------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
};

	if (!visit++)
	    for (int i = 0; i < 3; ++i) fputs(PslHeader[i], out_fd);
	swapskl(gsinf->skl);
	SKL	prv = *skl;
	for ( ; ++skl <= lst; psk = skl) {
	    int	dm = skl->m - psk->m;
	    int	dn = skl->n - psk->n;
	    if (dm && !dn)	{++qgap; qunp += dm;}
	    else if (!dm && dn)	{++tgap; tunp += dn;}
	    else if (dm && dn) {
		prv = *psk;
		for ( ; neorng(wkr) && wkr->rleft < psk->m; ++wkr) {
		    if (wkr->iscr > NEVSEL) {
			mch += wkr->mch;
			mmc += wkr->mmc;
		    }
		}
		if (psk->m == wkr->rleft) prv.n = wkr->left;
		++diag;
		tmp.m = q->SiteNz(prv.m);
		tmp.n = t->SiteNz(prv.n);
		tmp.p = dm;
		mfd.write(&tmp);
	    }
	}
	for ( ; neorng(wkr); ++wkr) {
	    if (wkr->iscr > NEVSEL) {
		mch += wkr->mch;
		mmc += wkr->mmc;
	    }
	}
	tmp.m = q->SiteNz(prv.m);
	tmp.n = t->SiteNz(prv.n);
	tmp.p = 0;
	mfd.write(&tmp);
	xyl = (SKLP*) mfd.flush();
	fprintf(out_fd, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c%c\t",
	    mch, mmc, rep, noN, qgap, qunp, tgap, tunp,
	    t->Strand(), q->Strand());
	--wkr;
	tunp = t->SiteNz(wkr->right);
	qunp = q->SiteNz(wkr->rright);
	wkr = fst;
	tgap = t->SiteNz(wkr->left);
	qgap = q->SiteNz(wkr->rleft);
	if (q->inex.sens) gswap(qgap, qunp);
	if (t->inex.sens) gswap(tgap, tunp);
	fprintf(out_fd, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t",
	    q->sqname(), q->len, qgap, qunp, 
	    t->sqname(), t->len, tgap, tunp, diag);
	if (!hetero && t->inex.sens & 2) {
	    tmp = xyl[diag];
	    for (int i = 0, rep = diag; i < rep; ++i, --rep)
		gswap(xyl[i], xyl[rep]);
	    for (int i = 0; i < diag; ++i)
		xyl[i].p = xyl[i+1].p;
	}
	for (int i = 0; i < diag; ++i)
	    fprintf(out_fd, "%ld,", xyl[i].p);
	putc('\t', out_fd);
	for (int i = 0; i < diag; ++i)
	    fprintf(out_fd, "%d,", xyl[i].m);
	putc('\t', out_fd);
	for (int i = 0; i < diag; ++i)
	    fprintf(out_fd, "%d,", xyl[i].n);
	putc('\n', out_fd);
	delete[] xyl;
	swapskl(gsinf->skl);
}

#endif

/* UCSC BED format */

static void BedForm(Gsinfo* gsinf, Seq* gene, Seq* qry)
{
static	int	visit = 0;
static	const	char*	BedHeader = 
	"track name=Spaln description=\"%s\" useScore=1\n";
static	const	char*	BedFrom = "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%s\t%d\t";

	if (!visit++) fprintf(out_fd, BedHeader, qry->sqname());
	EISCR*	fst = gsinf->eijnc->begin();
	EISCR*	wkr = fst;
	int	gstart = wkr->left;
	EISCR*	prv = wkr;
	while (neoeij(wkr)) prv = wkr++;
	int	gend = prv->right;
	int	nexon = wkr - fst;
	int*	epos = new int[nexon];
	int*	elen = new int[nexon];
	wkr = fst;
	prv = wkr;
	EISCR*	skp = wkr;
	int*	eps = epos;
	int*	eln = elen;
	bool	rvs = gene->inex.sens & REVERS;

	while (neoeij(wkr)) {
	    if (wkr->iscr > NEVSEL) {
		*eln++ = wkr->right - skp->left;
		if (rvs) {
		    *eps++ = gend - skp->right;
		} else {
		    *eps++ = skp->left - gstart;
		}
		prv = wkr++;
		if (!neoeij(wkr)) break;
		skp = wkr;
	    } else {
		skp = wkr++;
		if (!neoeij(wkr)) break;
	    }
	}
	nexon = eps - epos;
	gstart = gene->SiteNz(gstart);
	gend = gene->SiteNz(gend);
	if (rvs) {
	    gswap(gstart, gend);
	    vreverse(epos, nexon);
	    vreverse(elen, nexon);
	}
	float	bedscr = gsinf->rscr > 0? 
	    1000. * gsinf->scr / gsinf->rscr: gsinf->scr;
	const	char*	rgb = (gene->inex.sens & REVERS)? RGB_BLUE: RGB_RED;
	fprintf(out_fd, BedFrom, (*gene->sname)[0], gstart, gend, (*qry->sname)[0], 
	    min(1000, int(bedscr)), gene->Strand(), gstart, gend, rgb, nexon);
	eln = elen;
	for (int i = 0; i < nexon; ++i) fprintf(out_fd, "%d,", *eln++);
	putc('\t', out_fd);
	eps = epos;
	for (int i = 0; i < nexon; ++i) {
	    if (i < nexon - 1) fprintf(out_fd, "%d,", *eps++);
	    else	fprintf(out_fd, "%d\n", *eps++);
	}
	delete[] epos; delete[] elen;
}

/* similar to megablast -D 3 output format without E-value */
static	FILE*	fg = 0;
static	FILE*	fe = 0;
static	FILE*	fq = 0;
static void ExonForm(Gsinfo* gsinf, Seq* gene, Seq* qry)
{
	int	cds = 0;
	int	exon = 0;
	int	clen = 0;
	int	mch = 0;
	bool	binary = algmode.nsa == BIN_FORM;
	float	fscr = gsinf->scr / alprm.scale;
	VTYPE	iscr = 0;
	EISCR*	fst = gsinf->eijnc->begin();
	EISCR*	wkr = fst;
	EISCR*	prv = wkr;
	EISCR*	skp = wkr;
	char	str[MAXL];
	char*	convt = gene->isdrna()? nucl: ncodon;
	ExonRecord	er;
static	GeneRecord	gr;
static	int	visit = 0;
static	int	qryidx = 0;
static	char	intends[] = "  .  ";
static	char	fmt[] = "%s\t%s\t%7.2f\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t"
			"%7d\t%7.1f\t%7d\t%7.1f\t%7.2f\t%7.2f %2d %d %2d %d %s\n";
static	char	tfmt[] = "@ %s %c ( %d %d ) %s %d ( %d %d ) S: %.1f =: %.1f C: %.1f "
			"T#: %d T-: %d B#: %d B-: %d X: %d Pam: %d\n";
static	char	efmt[] = "Cant't write to gene record %s: # %ld\n";

	gr.nexn = gr.mmc = gr.unp = gr.bmmc = gr.bunp = gr.ng = 0;
	memset(&er, '\0', sizeof(ExonRecord));
	if (visit++ == 0) {
	    if (binary) {
		if (!OutPrm.out_file || !*OutPrm.out_file)
		    OutPrm.out_file = qry->spath;
		partfnam(str, OutPrm.out_file, "b");
		strcat(str, grext);
		if (!(fg = fopen(str, "wb"))) fatal(efmt, str, 0);
		partfnam(str, OutPrm.out_file, "b");
		strcat(str, erext);
		if (!(fe = fopen(str, "wb"))) fatal(efmt, str, 0);
		partfnam(str, OutPrm.out_file, "b");
		strcat(str, qrext);
		if (!(fq = fopen(str, "wb"))) fatal(efmt, str, 0);
		if ((fputs(dbs_dt[0]->dbsid, fq) == EOF) || (putc('\0', fq) == EOF))
		    fatal(efmt, qrext, gr.Nrecord);
		++qryidx;	// o-th q-recode contains the database name
	    } else {
		if (!out_fd) out_fd = stdout;
		fputs("# rID\t  gID\t   %id\t  ExonL\t MisMch\t Unpair\t "
		"ref_l\t  ref_r\t  tgt_l\t  tgt_r\t eScore\t IntrnL\t "
		"iScore\t Sig3/I\t Sig5/T  # -  X P DiNuc\n", out_fd);
	    }
	}
	while (neoeij(wkr)) {
	    if (qry->isdrna() && exon) {
		gr.bmmc += prv->mmc3 + wkr->mmc5;
		gr.bunp += prv->unp3 + wkr->unp5;
	    }
	    exon = wkr->mch + wkr->mmc + wkr->unp;
	    cds += exon;
	    if (wkr->iscr > NEVSEL) {
		clen += wkr->right - skp->left;
		mch += wkr->mch;
		gr.mmc += wkr->mmc;
		gr.unp += wkr->unp;
		++gr.nexn;
		er.Pmatch = exon? 100. * wkr->mch / exon: 0;
		er.Elen   = wkr->right - skp->left;
		er.Rleft  = qry->SiteNo(skp->rleft);
		er.Rright = qry->SiteNo(wkr->rright - 1);
		er.Gleft  = gene->SiteNo(skp->left);
		er.Gright = gene->SiteNo(wkr->right - 1);
		er.Bmmc   = prv->mmc3 + wkr->mmc5;
		er.Bunp   = prv->unp3 + wkr->unp5;
		er.Escore = (float) wkr->escr / alprm.scale;
		er.Iscore = (float) iscr / alprm.scale;
		er.Sig3   = (float) wkr->sig3 / alprm.scale;
		er.Sig5   = (float) wkr->sig5 / alprm.scale;
		if (binary) {
		    er.Nmmc = wkr->mmc;
		    er.Nunp = wkr->unp;
		    if (fwrite(&er, sizeof(ExonRecord), 1, fe) != 1)
			fatal(efmt, erext, gr.Nrecord);
		} else {
		    fprintf(out_fd, fmt, 
		    (*qry->sname)[0], (*gene->sname)[0], er.Pmatch, er.Elen,
		    wkr->mmc, wkr->unp, er.Rleft, er.Rright, er.Gleft, er.Gright,
		    er.Escore, er.Ilen, er.Iscore, er.Sig3, er.Sig5,
		    er.Bmmc, er.Bunp, er.miss, er.phase, intends);
		}
		iscr = wkr->iscr;
		prv = wkr++;
		if (!neoeij(wkr)) break;
		skp = wkr;
		er.Ilen = wkr->left - prv->right;
		er.phase = qry->isdrna()? cds % 3: (3 - prv->phs) % 3;
		er.miss = 0;
		intends[0] = er.Iends[0] = convt[*gene->at(prv->right)];
		intends[1] = er.Iends[1] = convt[*gene->at(prv->right+1)];
		intends[3] = er.Iends[2] = convt[*gene->at(wkr->left-2)];
		intends[4] = er.Iends[3] = convt[*gene->at(wkr->left-1)];
	    } else {
		gr.ng++;
		skp = wkr++;
		if (!neoeij(wkr)) break;
		er.miss = wkr->left - skp->right;
	    }
	}
	er.Ilen = qry->right - qry->left;
	gr.Gstart = gene->SiteNo(fst->left);
	gr.Gend   = gene->SiteNo(prv->right - 1);
	gr.Rstart = qry->SiteNo(qry->left);
	gr.Rend   = qry->SiteNo(qry->right - 1);
	gr.Gscore = fscr;
	gr.Pmatch = 100. * mch / cds;
	gr.Pcover = 100. * clen / er.Ilen;
	if (qry->isprotein()) gr.Pcover /= 3.;
	gr.Csense = gr.Gstart > gr.Gend;
	if (binary) {
	    gr.Cid = gene->did;
	    gr.Rid = qryidx++;
	    gr.Rsense = qry->inex.sens;
	    gr.Rlen = qry->len;
	    if (qry->isprotein() && gr.ng == 0) gr.ng = -1;
	    if (fwrite(&gr, sizeof(GeneRecord), 1, fg) != 1)
		fatal(efmt, grext, gr.Nrecord);
	    if ((fputs((*qry->sname)[0], fq) == EOF) || (putc('\0', fq) == EOF))
		fatal(efmt, qrext, gr.Nrecord);
	} else {
	    fprintf(out_fd, tfmt, (*gene->sname)[0], gr.Csense? '-': '+', 
	    gr.Gstart, gr.Gend,
	    (*qry->sname)[0], qry->len, gr.Rstart, gr.Rend,
	    gr.Gscore, gr.Pmatch, gr.Pcover, 
	    gr.mmc, gr.unp, gr.bmmc, gr.bunp, gr.ng, 0);
	}
	gr.Nrecord += gr.nexn;
	intends[0] = intends[1] = intends[3] = intends[4] = ' ';
}

void closeGeneRecord()
{
	if (fg) {fclose(fg); fg = 0;}
	if (fe) {fclose(fe); fe = 0;}
	if (fq) {fclose(fq); fq = 0;}
}

static void IntronForm(Gsinfo* gsinf, Seq* gene, Seq* qry)
{
	char	rvc = gene->Strand();
	EISCR*	fst = gsinf->eijnc->begin();
	EISCR*	wkr = fst;
	EISCR*	prv = wkr;
	char*	convt = gene->isdrna()? nucl: ncodon;
	int	cds = wkr->right - wkr->left;
	bool	isnuc = qry->isdrna();
static	int	visit = 0;
static	char	intends[] = "   ..   ";
static	char	tfmt[] = "@ %s %c ( %d %d ) %s %d ( %d %d )\n";
static	char	fmt[] = "%s\t%c %9d %9d  %d  %9d %9d\t%s\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t %s\n";

	if (visit++ == 0)
	    fputs("# gID\tdir   Donor  Acceptor Phs     tgt_5     tgt_3\trefID\t  ref_l\t  ref_r\t  Match\tMisMach\t Unpair\tIntronL\tIntronEnd\n", out_fd);
	while (neoeij(++wkr)) {
	    intends[0] = tolower(convt[*gene->at(prv->right-1)]);
	    intends[1] = convt[*gene->at(prv->right)];
	    intends[2] = convt[*gene->at(prv->right+1)];
	    intends[5] = convt[*gene->at(wkr->left-2)];
	    intends[6] = convt[*gene->at(wkr->left-1)];
	    intends[7] = tolower(convt[*gene->at(wkr->left)]);
	    int	intv = wkr->left - prv->right;
	    int	mch = isnuc? prv->mch3 + wkr->mch5: prv->mch;
	    int	mmc = isnuc? prv->mmc3 + wkr->mmc5: prv->mmc;
	    int	unp = isnuc? prv->unp3 + wkr->unp5: prv->unp;
	    if (prv->iscr > NEVSEL) {
		fprintf(out_fd, fmt, (*gene->sname)[0], rvc, gene->SiteNo(prv->right),
		gene->SiteNo(wkr->left - 1), qry->isdrna()? cds % 3: (3 - prv->phs) % 3,
		gene->SiteNo(prv->left), gene->SiteNo(wkr->right - 1),
		(*qry->sname)[0], qry->SiteNo(prv->rleft), 
		qry->SiteNo(wkr->rright - 1), mch, mmc, unp, intv, intends);
	    }
	    prv = wkr;
	    cds += wkr->right - wkr->left;
	}
	intends[0] = intends[1] = intends[2] = intends[5] = intends[6] = intends[7] = ' ';
	GeneRecord      gr;
	gr.Gstart = gene->SiteNo(fst->left);
	gr.Gend   = gene->SiteNo(prv->right) - 1;
	gr.Rstart = qry->SiteNo(qry->left);
	gr.Rend   = qry->SiteNo(qry->right) - 1;
	gr.Csense = gr.Gstart > gr.Gend;
	fprintf(out_fd, tfmt, (*gene->sname)[0], gr.Csense? '-': '+',
	    gr.Gstart, gr.Gend, (*qry->sname)[0], qry->len, gr.Rstart, gr.Rend);
}

/*	print a single (j-th) seqnuence		*/ 

void Seq::listseq(int j, bool sub)
{
	CHAR*	p = at(left) + j;
	INT	width = out_form->SeqBlkNo * out_form->SeqBlkSz;
	INT	block = out_form->SeqBlkSz;
const	char*	fname = (*sname)[j];

	if (out_form->DbName)
	    fprintf(out_fd, "%-11s %-10s %10d bp %10s\n",
		out_form->EntLabel, fname, right - left, moltype[inex.molc]);
	else {
	    if (out_form->EntLabel) {
		if (*out_form->EntLabel) { 
		    fprintf(out_fd, "%s%s", out_form->EntLabel, fname);
		    if (OutPrm.fastanno) {
			if (sub) fprintf(out_fd, " %d/%d", j + 1, many);
			fprintf(out_fd, attrfrmt2, sub? 1: many, len, 
			    SiteNo(left), SiteNo(right - 1));
		    }
		    putc('\n', out_fd);
		    width = block = OutPrm.lpw;
		} else { 
		    fprintf(out_fd, "%-15s\t", fname);
		}
	    } else	width = block = OutPrm.lpw;
	}
	if (descr && OutPrm.descrp)
	    fprintf(out_fd, "%s%c%s\n", out_form->DefLabel, 
		out_form->DbName? '\t': ' ', (*descr)[j]);
	if (out_form->SeqLabel) fprintf(out_fd, "%s\n", out_form->SeqLabel);
	if (out_form->SeqHead)  fputs(out_form->SeqHead, out_fd);
	int	i = 0;
	for (int m = left; m < right; ++m, p += many) {
	    if (OutPrm.trimendgap && IsGap(*p)) continue;
	    if (width && i % width == 0 && out_form->SeqForm)
		fprintf(out_fd, out_form->SeqForm, m + 1);
	    putc(code->decode[*p], out_fd);
	    if (!width) continue;
	    if (++i % width == 0) putc('\n', out_fd);
	    else if (i % block == 0) putc(' ', out_fd);
	}
	if (OutPrm.asterisk) putc('*', out_fd);
	if (!width || i % width) putc('\n', out_fd);
	if (out_form->EndLabel) fprintf(out_fd, "%s\n", out_form->EndLabel);
}

static void CigarForm(Gsinfo* gsi, Seq* gen, Seq* qry)
{
	if (!gsi->cigar || !gsi->cigar->rec || !gsi->skl) {
	    prompt("No cigar records !\n");
	    return;
	}
	SKL*	fst = gsi->skl + 1;
	SKL*	lst = gsi->skl + gsi->skl->n;
	fprintf(out_fd, "cigar: %s %d %d %c %s %d %d %c %d",
	    qry->sqname(), qry->SiteNm(fst->n), 
	    qry->SiteNm(lst->n), qry->inex.sens? '-': '+',
	    gen->sqname(), gen->SiteNm(fst->m), 
	    gen->SiteNm(lst->m), gen->inex.sens? '-': '+',
	    (int) gsi->fstat.val);
	CIGAR*	tm = gsi->cigar->rec + gsi->cigar->size();
	for (CIGAR* wk = gsi->cigar->rec; wk < tm; ++wk)
	    fprintf(out_fd, " %c %d", wk->ope, wk->len);
	putc('\n', out_fd);
}

static void VulgarForm(Gsinfo* gsi, Seq* gen, Seq* qry)
{
	if (!gsi->vlgar || !gsi->vlgar->rec || !gsi->skl) {
	    prompt("No vlgar records !\n");
	    return;
	}
	SKL*	fst = gsi->skl + 1;
	SKL*	lst = gsi->skl + gsi->skl->n;
	fprintf(out_fd, "vulgar: %s %d %d %c %s %d %d %c %d",
	    qry->sqname(), qry->SiteNm(fst->n), 
	    qry->SiteNm(lst->n), qry->inex.sens? '-': '+',
	    gen->sqname(), gen->SiteNm(fst->m), 
	    gen->SiteNm(lst->m), gen->inex.sens? '-': '+',
	    (int) gsi->fstat.val);
	VULGAR*	tm = gsi->vlgar->rec + gsi->vlgar->size();
	for (VULGAR* wk = gsi->vlgar->rec; wk < tm; ++wk)
	    fprintf(out_fd, " %c %d %d", wk->ope, wk->alen, wk->blen);
	putc('\n', out_fd);
}

static void SamForm(Gsinfo* gsi, Seq* gen, Seq* qry)
{
	Samfmt* samfm = gsi->samfm;
	if (!samfm || !samfm->rec || !gsi->skl) {
	    prompt("No sam records !\n");
	    return;
	}
	fprintf(out_fd, "%s\t%d\t%s\t%d\t%d\t", qry->sqname(),
	    samfm->flag, gen->sqname(), gen->SiteNo(samfm->pos), samfm->mapq);
	CIGAR*	tm = samfm->rec + samfm->size();
	for (CIGAR* wk = samfm->rec; wk < tm; ++wk)
	    fprintf(out_fd, "%d%c", wk->len, wk->ope);
	putc('\t', out_fd);
	fprintf(out_fd, "*\t%d\t%d\t", 
	   samfm->tlen? gen->SiteNo(samfm->pnext): 0, samfm->tlen);
	gswap(qry->left, samfm->left);
	gswap(qry->right, samfm->right);
	if (gen->inex.sens) qry->comrev();
	qry->typeseq(0, true);
	if (gen->inex.sens) qry->comrev();
	gswap(qry->left, samfm->left);
	gswap(qry->right, samfm->right);
	fprintf(out_fd, "\t*\tNM:i:%d\tAS:i:%d\n", 
	    int(gsi->fstat.mmc + gsi->fstat.unp), int(gsi->scr / alprm.scale));
}

/*************************************************************
*	exon-intron structure
*	omode = 0:	Gff3 format
*	omode = 1:	alignment
*	omode = 2:	gff-pwa format
*	omode = 3:	Bed format
*	omode = 4:	Exon-oriented
*	omode = 5:	Intron-oriented
*	omode = 6:	catinated exon seqs.
*	omode = 7:	translated seq.
**************************************************************/

void printgene(Seq* seqs[], Gsinfo* gsi)
{
	int	m = seqs[0]->inex.intr? 0: 1;
	Seq*&	gene = seqs[m];
	Seq*&	qry = seqs[1-m];
	Seq*	tsd = 0;
	INT	width = out_form->SeqBlkNo * out_form->SeqBlkSz;
	int	block = OutPrm.lpw;
	int	i = gene->right - gene->left;
	float	fscr = (float) gsi->scr / alprm.scale;
	RANGE*	rng = gsi->CDSrng;
	bool	swp = false;

	switch (algmode.nsa) {		// without sequence
	    case GFF_FORM: Gff3Form(gsi, gene, qry); return;
	    case PWA_FORM: Gff3PWA(gsi, gene, qry); return;
//	    case PSL_FORM: PslForm(gsi, gene, qry); return;
	    case BED_FORM: BedForm(gsi, gene, qry); return;
	    case EXN_FORM:
	    case BIN_FORM: ExonForm(gsi, gene, qry); return;
	    case ITN_FORM: IntronForm(gsi, gene, qry); return;
	    case CDS_FORM: case AAS_FORM:	// 
		tsd = gene->splice(0, rng, 0);
		tsd->tron2nuc(0);
		swp = OutPrm.deflbl == 2;
		break;
	    case CIG_FORM: CigarForm(gsi, gene, qry); return;
	    case VLG_FORM: VulgarForm(gsi, gene, qry); return;
	    case SAM_FORM: 
	    case BAM_FORM: SamForm(gsi, gene, qry); return;
	    case PSJ_FORM:	// predicted gene strucure of the reference
		rng = gsi->querygs(qry);
		swp = true;
		break;
	}
	if (swp) swapseq(seqs, seqs + 1);
	const	char*	ofn = OutPrm.out_file? strrchr(OutPrm.out_file, '/'): 0;
	if (OutPrm.out_file) ofn = ofn? ofn + 1: OutPrm.out_file;
	else if (OutPrm.deflbl != 3) ofn = gene->sqname();
	if (out_form->DbName) {
	    fprintf(out_fd, "%-11s %-10s %10d bp %10s %s %9.2f\n", 
		out_form->EntLabel, gene->sqname(), i, moltype[qry->inex.molc], 
		qry->sqname(), fscr);
	    block = out_form->SeqBlkSz;
	} else if (out_form->EntLabel || fst_form) {
	    int	gl = gene->left;
	    int	gr = gene->right;
	    int	ql = qry->left;
	    int	qr = qry->right;
	    if (gsi->eijnc) {
		gl = gsi->eijnc->genleft();
		gr = gsi->eijnc->genright();
		ql = gsi->eijnc->refleft();
		qr = gsi->eijnc->refright();
	    }
	    fprintf(out_fd, "%s", (fst_form? fst_form: out_form)->EntLabel);
	    if (ofn) {
		fprintf(out_fd, "%s", ofn);
		if (OutPrm.out_file) fprintf(out_fd, " %s", gene->sqname());
	    } else fprintf(out_fd, "%s.%d %s", gene->sqname(), 
		gene->SiteNo(gsi->center(1)) / 1000, gene->sqname());
	    fprintf(out_fd, attrfrmt, gene->Strand(), gene->len, 
		gene->SiteNo(gl), gene->SiteNo(gr - 1));
	    fprintf(out_fd, " %s %c", qry->sqname(), qry->Strand());
	    fprintf(out_fd, attrfrmt3, qry->many, qry->len,
		qry->SiteNo(ql), qry->SiteNo(qr - 1));
	    fputs(swp? (algmode.nsa == PSJ_FORM? " Y": " T"): " N", out_fd);
	    fprintf(out_fd, " %9.2f\n", fscr);
	}
	if (algmode.nsa != PSJ_FORM && swp) swapseq(seqs, seqs + 1);
	GBcdsForm(rng, gene);
	if (algmode.nsa == AAS_FORM) {		// translate
	    tsd->transout(out_fd, 0, 0);
	} else if (algmode.nsa == PSJ_FORM) {
	    swapseq(seqs, seqs + 1);
	    qry->typeseq();
	} else if (tsd) {			// without translate
	    m = tsd->left;
	    CHAR*	p = tsd->at(tsd->left);
	    if (out_form->SeqLabel) fprintf(out_fd, "%s\n", out_form->SeqLabel);
	    if (out_form->SeqHead)  fputs(out_form->SeqHead, out_fd);
	    for (i = 0; m < tsd->right; ++m, p += tsd->many) {
		if (width && i % width == 0 && out_form->SeqForm)
		    fprintf(out_fd, out_form->SeqForm, m + 1);
		putc(tsd->code->decode[*p], out_fd);
		if (!width) continue;
		if (++i % width == 0) putc('\n', out_fd);
		else if (i % block == 0) putc(' ', out_fd);
	    }
	    if (width && i % width) putc('\n', out_fd);
	    if (out_form->EndLabel) fprintf(out_fd, "%s\n", out_form->EndLabel);
	}
	if (OutPrm.spjinf && gene->exin && 
		(algmode.nsa == CDS_FORM || algmode.nsa == AAS_FORM))
	    gsi->BoundaryInf(gene);
	if (tsd != qry) delete tsd;
	if (rng != gsi->CDSrng) delete[] rng;
}

enum AaProp {GAP, SML, PLS, MNS, HPO, ARO = 6, NEU, PHI, PHO};

static int chempro(int* cmp, int ii)
{
	static int proch[2][ASIMD] = {
	    {GAP, GAP, GAP, SML, PLS, MNS, MNS, CYS, MNS, MNS, SML, PLS, HPO, 
	     HPO, PLS, HPO, ARO, SML, SML, SML, ARO, ARO, HPO, MNS, MNS}, 
	    {GAP, GAP, GAP, NEU, PHI, PHI, PHI, PHO, PHI, PHI, NEU, PHI, PHO, 
	     PHO, PHI, PHO, PHO, NEU, NEU, NEU, PHO, PHO, PHO, PHI, PHI}
	};
	static char chemcode[] = " .+_@C$.jo";
	int	p = proch[0][ii];
	int	s = proch[1][ii];
	bool	pc = true;
	bool	sc = true;

	for (int i = AMB; ++i < ASIMD; ) {
	    if (cmp[i]) {
		if (p != proch[0][i]) pc = false;
		if (s != proch[1][i]) sc = false;
	    }
	}
	if (pc) return chemcode[proch[0][ii]];
	if (sc) return chemcode[proch[1][ii]];
	return (' ');
}

static int logonuc(int* cmp, int ii, const SEQ_CODE* defcode)
{
	int	c = 0;
	int	n = cmp[gap_code]? 1: 0;
	int	b = A - 1;

	if (IsGap(ii)) return (' ');
	for (int i = 1; i < NTS; i++)
	    if (cmp[i+b]) c |= i;
	for (int i = 1; i < NTS; i <<= 1)
	    if (c & i) n++;
	if (n == 1) return(defcode->decode[c]);
	return (n == 2? tolower(defcode->decode[c+b]): ' ');
}

static int csym(int cmp[], int rows, const SEQ_CODE* defcode)
{
	int 	ii = 0;
	int 	cc = 0;
	int 	gg = rows;

	for (int i = 0; i < defcode->max_code; ++i) {
		gg -= cmp[i];
		if (cmp[i]) ++cc;
		if (cmp[i] > cmp[ii]) ii = i;
	}
	if (cmp[nil_code]) return (' ');	/*  blank	*/
	else {
	    cc += gg - 1;
	    if (cc == 0)		/*  Conserved	*/
		return (defcode->decode[ii]);
	    else if (emphsim > 0) {
		int	m = 10 * (rows - cmp[ii]) / rows;
		if (m < emphsim)
			return (tolower(defcode->decode[ii]));
		else
			return (' ');
	    } else if (emphsim < 0) {
		cc = nmk - cc  + 1;
		if (cc < 0) cc = 0;
		if (cc > nmk) cc = nmk;
		return (mark[cc]);
	    } else if (defcode->ceil_code == NTS)
		return (logonuc(cmp, ii, defcode));
	    else
		return (chempro(cmp, ii));
	}
}

void setprmode(int pmd, int lorn, int trc)
{
	static char msg1[] = 
	"Native Form 0[None]/1[Last]/2[Every]/3[First]/4[Ditto]\n";
	static char msg2[] = 
	"Foreign Form 5[Bare]/6[Phylp]/7[GCG]/8[CLW]/9[Ser]/10[GDE] (%d) : ";
	static char msg3[] = "None[0]/Number[2]/Name[3] (%d) ";
	static char msg4[] =
	"Codon Code P[rotein]/D[NA]/R[NA]/T[RC] ";

	SeqDb*	sdb = 0; 
	if (pmd == 'c') prmode = Form_CLW;
	else if (pmd == 'm') prmode = Form_GCG;
	else if (pmd == 'h') prmode = Form_Phylp;
	else if ('a' <= pmd && pmd <= 'z') {
	    sdb = setform(pmd);
	    if (sdb->FormID == FASTA || sdb->FormID == PIR || sdb->FormID == NEXUS)
		prmode = Form_Serial;
	} else if ('0' <= pmd && pmd <= '9')  prmode = (Row_Mode) (pmd - '0');
	else if (pmd >= 0) prmode = (Row_Mode) pmd;
	else if (pmd == QUERY) {
	    prompt(msg1);
	    promptin(msg2, &prmode);
	    switch (prmode) {
		case Row_Last:
		case Row_Every:
		    promptin("Consensus level [-1 - 9] (%d) ? ", &emphsim);
		    break;
		case Row_First: ditto = _SAME; break;
		case Row_Ditto: ditto = _IBID; break;
		default: break;
	    }
	}
	switch (prmode) {
	    case Form_Phylp:
		lLblMode = LABEL_Once;
		rLblMode = LABEL_None; break;
	    case Form_GCG:
	    case Form_CLW:
		lLblMode = LABEL_Name;
		rLblMode = LABEL_None; break;
	    case Form_GDE:	OutPrm.noseqline = 1;
	    case Form_Serial:	if (!sdb) setform('f'); break;
	    default:
		lLblMode = LABEL_Numb;
		switch (lorn) {
		    case 'L': rLblMode = LABEL_Name; break;
		    case 'N': rLblMode = LABEL_Numb; break;
		    case 'O': clmark = 0; break;
		    case 'M': clmark = 1; break;
		    case 'C': clmark = 2; break;
		    case QUERY:
			promptin(msg3, &rLblMode);
			promptin("Mark/Number column [0-2] (%d) ? ", &clmark);
			break;
		    default: break;
		}
	}
	if (trc == QUERY) trc = progetc(msg4);
	switch (trc) {
	    case 'a':
	    case 'A':
	    case 'p':
	    case 'P':	prtrc = PROTEIN; break;
	    case 'd':
	    case 'D':	prtrc = DNA; break;
	    case 'g':
	    case 'G':	prtrc = GENOME; break;
	    case 'r':
	    case 'R':	prtrc = RNA; break;
	    case 't':
	    case 'T':	prtrc = TRON; break;
	    default:	break;
	}
}

/*
	Following two functions were ported from CLUSTAW package 
	and modified.
*/

static int SeqGCGCheckSum(char* seq, int len)
{
	long	check = 0;
	
	for (int i = 0; i < len; ++i, ++seq)
	    check += ((i % 57) + 1) * toupper(*seq);
	return (check % 10000);
}

static int checksum(int* checks, GAPS* gaps, Seq* sd)
{
	GAPS*	gp;
	int	len = 0;
	char	gapchar = sd->code->decode[gap_code];

	for (gp = gaps; gaps_intr(++gp); ) len += gp->gln;
	len += gp->gps - gaps->gps;
	CHAR*	ps = sd->at(gaps->gps);
	CHAR*	ts = sd->at(gp->gps);

	char*	seq = new char[len];
	for (int i = 0; i < sd->many; ++i, ++ps) {
	    int c = gaps->gps;
	    gp = gaps;
	    CHAR*	ws = ps;
	    for (char* psq = seq; ws < ts; ) {
		if (c == gp->gps) {
		    if (gp->gln > 0) {
			memset(psq, gapchar, gp->gln);
			c += gp->gln;
			psq += gp->gln;
		    }
		    if (gaps_intr(gp)) ++gp;
		    else	break;
		} else {
		    *psq++ = sd->isprotein()? amino[*ws]: nucl[*ws];
		    ++c;
		    ws += sd->many;
		}
	    }
	    checks[i] = SeqGCGCheckSum(seq, len);
	}
	delete[] seq;
	return (len);
}

static void gcg_out(GAPS* gaps[], Seq* seqs[], int seqnum)
{
	setlpw(GCG_LINELENGTH);
	setgapchar(seqs[0], GCG_GAP_CHAR);
	OutPrm.BlkSz = GCG_SeqBlkSz;
	int	rows = 0;
	for (int i = 0; i < seqnum; ++i) rows += seqs[i]->many;
	int*	checks = new int[rows];
	int	len = rows = 0;
	for (int i = 0; i < seqnum; rows += seqs[i++]->many)
	    len = checksum(checks + rows, gaps[i], seqs[i]);
	long	grand_checksum = 0;
	for (int i = 0; i < rows; ++i) grand_checksum += checks[i];
	grand_checksum %= 10000;

	fputs("PileUp\n\n", out_fd);
	fprintf(out_fd, "   MSF:%5d  Type: ", len);
	if (seqs[0]->isprotein())	fprintf(out_fd,"P");
	else	fprintf(out_fd,"N");
	fprintf(out_fd,"    Check:%6ld   .. \n\n", grand_checksum);
	rows = 0;
	for (int i = 0; i < seqnum; ++i) {
	    for (int j = 0; j < seqs[i]->many; ++j, ++rows)  {
		fprintf(out_fd, " Name: %-24s oo  Len:%5d  Check:%6ld  Weight:",
		(*seqs[i]->sname)[j], seqs[i]->len, (long) checks[rows]);
#if USE_WEIGHT
		if (seqs[i]->weight)
		    fprintf(out_fd, "%8.4f\n", (double) seqs[i]->weight[j]);
		else
#endif
		    fprintf(out_fd, "  1.00\n");
	    }
	}
	fprintf(out_fd, "\n//\n");  
	delete[] checks;
}

int Seq::calcnbr(int gp, int i)
{
	int 	n = nbr[i];
	CHAR*	s = at(0) + i;

	for ( ; gp-- > 0; s += many) if (IsntGap(*s)) ++n;
	return (n);
}

void Seq::fphseq(int n, FILE* fd)
{
	if (!fd) fd = out_fd;
	char	str[MAXL];
static	const char* frm[] = {"%c%s [%d:%d] ", "%s ( %d - %d )",  "%-16s ( %3d - %3d )"};

	if (n <= 0) {
	    fputs(sqname(), fd);
	} else if (--n == 0) {
	    fprintf(fd, frm[0], senschar(), sqname(), many, len);
	} else {
	    if (n > 2) n = 1;
	    sprintf(str, frm[0], senschar(), sqname(), many, len);
	    fprintf(fd, frm[n], str, SiteNo(left), SiteNo(right - 1));
	}
}

void fphseqs(Seq* seqs[], int n)
{
	putc('\n', out_fd);
	while (n--) {
	    (*seqs++)->fphseq();
	    if (n) fputs(" - ", out_fd);
	}
	putc('\n', out_fd);
}

#if USE_WEIGHT
void Seq::fpweight()
{
	if (!weight) return;
	int 	i = 0;
	while (i < many) {
	    if (i % NO_WGHT == 0) putc(_WGHT, out_fd);
	    fprintf(out_fd, " %14.7e", weight[i]);
	    if (++i % NO_WGHT == 0) putc('\n', out_fd);
	}
	if (i % NO_WGHT) putc('\n', out_fd);
}
#endif

void fprint_seq_mem(Seq* seqs[], int n)
{
	for (int i = 0; i < n; ++i) {
	    seqs[i]->fphseq(algmode.rng? 2: 0);
	    putc(' ', out_fd);
	}
	putc('\n', out_fd);
}

/*	Map the sequence number to column position
	Array must be sorted and terminated with a minus value
	Both number and position start with 1	*/

void Seq::num2pos(int which, int* array)
{
	if (which >= many || inex.vect) {
	    prompt("Improper input or seq form!\n");
	    return;
	}

	CHAR*	ss = at(0) + which;
	CHAR*	tt = at(len);
	int	num = 0;
	int	pos = 0;

	while (*array > 0) {
	    while (num < *array) {
		if (IsntGap(*ss)) ++num;
		++pos;
		if ((ss += many) >= tt) {
		    *array = pos;
		    return;
		}
	    }
	    *array++ = pos;
	}
}

void Seq::pos2num(int which, int* array)
{
	if (which >= many || inex.vect) {
	    prompt("Improper input or seq form!\n");
	    return;
	}

	CHAR*	ss = at(0) + which;
	CHAR*	tt = at(len);
	int	num = 0;
	int	pos = 0;

	while (*array > 0) {
	    while (pos < *array) {
		if (IsntGap(*ss)) ++num;
		pos++;
		if ((ss += many) >= tt) {
		    *array = pos;
		    return;
		}
	    }
	    *array++ = num;
	}
}

bool Seq::findGate(RANGE* gate)
{
	double	upper = NGP;
	if (upper > 1.) upper /= 100;
	double	lower = many * upper;
	upper *= many * alprm.thr;
	saverange(gate);
	int	cg = 0;
	int	i = left;
	CHAR*	s = at(left);
	while (i++ < right) {		// forward
	    CHAR*	t = s + many;
	    int ng = 0;
	    for ( ; s < t; ++s)
		if (IsntGap(*s)) ++ng;
	    if (ng <= lower) {cg = 0; gate->left = i;}
	    else if ((cg += ng) > upper) break;
	}
	if (i == right) return (false);	// No conserved region
	s = at(right) - 1;
	cg = 0;
	i = right;
	while (i-- > gate->left) {	// backward
	    int ng = 0;
	    for (CHAR* t = s - many ; s > t; --s)
		if (IsntGap(*s)) ++ng;
	    if (ng <= lower) {cg = 0; gate->right = i;}
	    else if ((cg += ng) >= upper) break;
	}
	return (true);
}

int Seq::calcResNum(int i)
{
	int	n = 0;
	CHAR*	s = at(left) + i;
	CHAR*	t = at(right);

	for ( ; s < t; s += many)
	    if (IsntGap(*s)) ++n;
	return (n);
}

typedef	int	 DIM3[3];

void Seq::fpmem_len(FILE* fd)
{
	RANGE	gate;
	bool	border = (algmode.nsa & 4)
		&& many > 1 && findGate(&gate);
	DIM3*	insrt = 0;
	DIM3*	delet = 0;
	
	if (border) {
	    insrt = new DIM3[many];
	    delet = new DIM3[many];
	    vclear(insrt, many);
	    vclear(delet, many);
	    CHAR*	s = at(left);
	    for (int i = 0, k = 0; i < right; ) {
		int	unps = 0;
		for (int j = 0; j < many; ++s, ++j)
		    if (IsGap(*s)) ++unps;
		if (unps) {
		    s -= many;
		    for (int j = 0; j < many; ++s, ++j) {
			if (IsGap(*s))  delet[j][k] += many - unps;
			else		insrt[j][k] += unps;
		    }
		}
		if (++i == gate.left || i == gate.right) ++k;
	    }
	}
	PrintMember	prm(sname, true, " ");
	for (int j = 0; j < many; ++j) {
	    if (algmode.nsa & 1) fprintf(fd, "%5d ", j+1);
	    prm.put_member(fd, j);
	    if (algmode.nsa & 2) fprintf(fd, "%6d", calcResNum(j));
	    if (border) {
		for (int k = 0; k < 3; ++k) {
		    int	unps = (insrt[j][k] > delet[j][k])?
			insrt[j][k]: -delet[j][k];
		    fprintf(fd, " %6d", unps / (many - 1));
		}
	    }
	    fputc('\n', fd);
	}
	delete[] insrt; delete[] delet;
}

PrintMember::PrintMember(Strlist* sn, bool pad_space, const char* tl) 
	: sname(sn), g2t(0)
{
	if (OutPrm.taxoncode) {
	    g2t = ftable.read_gnm2tab(OutPrm.taxoncode);
	    if (pad_space)
		sprintf(fmt, "%%s:%%-%ds", sname->longest());
	    else
		strcpy(fmt, "%s:%s");
	} else {
	    if (pad_space)
		sprintf(fmt, "%%-%ds", sname->longest());
	    else
		strcpy(fmt, "%s");
	}
	if (tl) strcat(fmt, tl);
}

void PrintMember::put_member(FILE* fd, int i)
{
	if (g2t) {
	    char*	taxon = 0;
	    g2t->taxon_code((*sname)[i], &taxon);
	    fprintf(fd, fmt, (taxon? taxon: unknown_taxon), (*sname)[i]);
	} else
	    fprintf(fd, fmt, (*sname)[i]);
}

/*************************************************************
*	Output alignment
*************************************************************/

static	const	int iis_color[] = {CMC_BLACK, CMC_RED, CMC_GREEN, CMC_BLUE};
static	const	char* hml_color[] = {"black", "red", "green", "blue"};
static	const	char spc9[] = "	 ";
static	const	char spc12[] = "	    ";

PrintAln::PrintAln(GAPS* _gaps[], Seq* _seqs[], int _seqnum)
	: gaps(_gaps), seqs(_seqs), seqnum(_seqnum)
{
	globpt = 1;
	cpm = prmode;
	markeij = 0;
	if (OutPrm.ColorEij) {
	    for (int j = 0; j < seqnum; ++j) {
		if (seqs[j]->sigII && seqs[j]->sigII->pfq) {
		    markeij = OutPrm.ColorEij;
		    break;
		}
	    }
	}
	if (markeij) {
	    char	str[MAXL];
	    sprintf(str, "%s: %s", seqnum == 1? "Prrn": "Aln", seqs[0]->sqname());
	    for (int n = 1; n < seqnum; ++n) {
		strcat(str, " + ");
		strcat(str, seqs[n]->sqname());
	    }
	    ecc = markeij == 1? new EscCharCtl(out_fd): 0;
	    hcc = markeij == 2? new HtmlCharCtl(out_fd, str): 0;
	    pfqs = new PFQ*[seqnum];
	    tfqs = new PFQ*[seqnum];
	    lsts = new int*[seqnum];
	} else {
	    ecc = 0; hcc = 0; pfqs = tfqs = 0; lsts = 0;
	}
	for (int i = rows = 0; i < seqnum; ++i) rows += seqs[i]->many;
	gp  = new GAPS*[seqnum];
	seq = new CHAR*[seqnum];
	wkr = new CHAR*[seqnum];
	agap = new int[seqnum];
	nbr = new int[rows];
	left = new int[rows];
	image = new CHAR*[rows];
	wbuf = new CHAR[rows * OutPrm.lpw];
}

PrintAln::~PrintAln()
{
	prmode = cpm;
	delete[] image[0]; delete[] image; delete[] left; delete[] nbr;
	delete[] gp; delete[] seq; delete[] wkr; delete[] agap;
	delete ecc; delete hcc;
	delete[] pfqs; delete[] tfqs; delete[] lsts;
}

void PrintAln::seqline(const CHAR* qry, const CHAR* brc, 
	const CHAR* cur, int rpl, const char decode[])
{
	for (int i = 0; i < rpl; ++i, ++cur) {
	    int	res = *cur & RES_BITS;
	    int	ccd = (*cur >> 5) & 3;
	    if (OutPrm.BlkSz && (i % OutPrm.BlkSz) == 0)
		putc(' ', out_fd);
	    int	ch;
	    if (res == BLANK || (brc && *brc == BLANK))
		ch =' ';
	    else if (*cur & INTRONBIT)
		ch = tolower(decode[res]);
	    else if (qry)
		ch = *qry == res? ditto: tolower(decode[res]);
	    else
		ch = decode[res];
	    if (ccd) {
		if (ecc) ecc->putchr(ch, CMC_WHITE, iis_color[ccd], CMC_BOLD);
		if (hcc) hcc->putchr(ch, "white", hml_color[ccd]);
	    } else	putc(ch, out_fd);
	    if (qry) ++qry;
	    if (brc) ++brc;
	}
}

void PrintAln::putclmark()
{
	int	i;

	if (clmark == 1) {
	    for (i = globpt; i % 10; i++) putc(' ', out_fd);
	    for ( ; i < globpt + OutPrm.lpw; i += 10)
		fputs("       .  ", out_fd);
	} else {
	    fprintf(out_fd, ";%03d", globpt / OutPrm.lpw + 1);
	    for (i = globpt; i % 10; i++) putc(' ', out_fd);
	    fprintf(out_fd, "%6d", i);
	    for (i += 10; i < globpt + OutPrm.lpw; i += 10)
		fprintf(out_fd, "%10d", i);
	}
	putc('\n', out_fd);
}

void PrintAln::pair_mrk(const CHAR* p1, const CHAR* p2, int clms)
{
	Simmtx*	sm = getSimmtx(0);
	fputs(spc9, out_fd);
	for (int i = 0; i < clms; ++i) {
	    int	cc = (*p1 == BLANK || *p2 == BLANK)? 0:
	    	sm->simgrade(*p1++ & RES_BITS, *p2++ & RES_BITS);
	    if (emphsim < 0) cc = nmk - cc;
	    if (cc > nmk) cc = nmk;
	    if (cc < 0) cc = 0;
	    putc(mark[cc], out_fd);
	}
	putc('\n', out_fd);
}

void PrintAln::calc_mrk(int rows, int clms, const SEQ_CODE* defcode)
{
	int*	cmp = new int[defcode->max_code];

	fputs(spc9, out_fd);
	for (int l = 0; l < clms; ++l) {
	    vclear(cmp, defcode->max_code);
	    for (int i = 0; i < rows; ++i) {
		int	cc = image[i][l] & RES_BITS;
		if (cc == BLANK) ++cmp[nil_code];
		else		++cmp[cc];
	    }
	    putc(csym(cmp, rows, defcode), out_fd);
	}
	putc('\n', out_fd);
	delete[] cmp;
}

void PrintAln::prnt_aln(const int* lft, const int* right)
{
	CHAR**	img = image;
	CHAR*	prv = 0;
	CHAR*	brc = (htl == 3)? image[pro]: 0;
	CHAR*	ns;
	char*	decode;
	SEQ_CODE*	defcode = (*seqs)->code;

	putc('\n', out_fd);
	if (clmark) putclmark();
	int	k = 0;
	for (int j = 0; j < seqnum; ++j) {
	  Seq*&	sd = seqs[j];
	  decode = defcode->decode;
	  if (sd->inex.molc == TRON) {
	    switch (prtrc) {
	      case DNA: case RNA: case GENOME:	decode = ncodon; break;
	      case PROTEIN:	decode = acodon; break;
	      case TRON:
	      default:	break;
	    }
	  }
	  PrintMember	prm(sd->sname, false, "\n");
	  for (int i = 0; i < sd->many; ++i, ++k, ++img) {
	    if (htl == 3 && sd->inex.molc == TRON) {
	      decode = defcode->decode;
	      int n = 0;
	      switch (lLblMode) {
		case LABEL_Name: n = 15; break;
		case LABEL_Numb: n = 9; break;
		case LABEL_Once: n = globpt == 1? 14: 0; break;
		case LABEL_None:
	    	default:	break;
	      }
	      while (n--) putc(' ', out_fd);
	      seqline(prv, brc, *img, OutPrm.lpw, decode);
	      putc('\n', out_fd);
	      for (n = 0, ns = *img; n < OutPrm.lpw; ++n, ++ns) {
		if (*ns != BLANK) {
		    int	intbit = *ns & INTRONBIT;
		    *ns = aa2nuc[*ns & RES_BITS] | intbit;
		}
	      }
	      decode = nucl;
	    }
	    switch (lLblMode) {
		case LABEL_Once:
		    if (globpt == 1) fprintf(out_fd, "%-12s", (*sd->sname)[i]);
		    else	fputs(spc12, out_fd);
		    break;
		case LABEL_Name:
		    fprintf(out_fd, "%-24s", (*sd->sname)[i]); break;
		case LABEL_Numb:
		    fprintf(out_fd, "%8d ", sd->SiteNo(*lft) % 100000000);
		    ++lft; break;
		case LABEL_None:
	    	default:	break;
	    }
	    seqline(prv, 0, *img, OutPrm.lpw, decode);
	    if (j == gene) {
		ns = *img;
		for (int n = 0; n < OutPrm.lpw; ++n, ++ns)
		    *ns = *ns & ~INTRONBIT;
	    }
	    switch (rLblMode) {
		case LABEL_None: putc('\n', out_fd); break;
		case LABEL_Numb:
		    fprintf(out_fd, "%c%5d\n", _LABL, sd->SiteNo(*right - 1)); 
		    ++right; break;
		case LABEL_Name:
	    	default:
		    fprintf(out_fd, "%c ", _LABL);
		    prm.put_member(out_fd, i); break;
	    }
	    if (k < rows - 1 && prmode == Row_Every)
		pair_mrk(*img, img[1], OutPrm.lpw);
	    if (prmode == Row_First) prv = *image;
	    if (prmode == Row_Ditto) prv = *img;
	  }
	}
	if (prmode == Row_Last)
	    calc_mrk(rows, OutPrm.lpw, defcode);
	else if (prmode == Form_GCG)
	    putc('\n', out_fd);
}

void PrintAln::markiis(int k, int j, int clm)
{
	int phs = (pfqs[j]->pos % 3 + 1) << 5;
	for (int n = 0; n < pfqs[j]->num; ++n) {
	    int	i = lsts[j]? *lsts[j]++: 0;
	    image[k + i][clm] |= phs;
	}
	++pfqs[j];
}

void PrintAln::printaln()
{
	int	k = 0;
	int	gpos = 0;
	int 	active = seqnum;
	RANGE*	exon = 0;
	pro = gene = -1;
	int	maxleft = 0;
	for (int j = htl = 0; j < seqnum; ++j) {
	    Seq*&	sd = seqs[j];
	    if (sd->left > maxleft) maxleft = sd->left;
	    if (sd->isprotein()) {
		htl |= 1;
		pro = k;
	    } else if (seqs[j]->inex.molc == TRON)
		htl |= 2;
	    if (sd->inex.intr && sd->exons) {
		exon = sd->exons + 1;
		gene = j;
		gpos = gaps[j]->gps;
	    }
	    for (int i = 0; i < sd->many; ++i, ++k, wbuf += OutPrm.lpw) {
		nbr[k] = sd->calcnbr(gaps[j]->gps, i);
		image[k] = wbuf;
	    }
	    wkr[j] = seq[j] = sd->at(gaps[j]->gps);
	    gp[j] = gaps[j];
	    if (!gp[j]->gln) ++gp[j];
	    agap[j] = 0;
	}
	int	c_step = (htl == 1)? 3: 1;
	if (markeij) {
	    for (int j = 0; j < seqnum; ++j) {
		Seq*&	sd = seqs[j];
		if (sd->sigII) {
		    sd->sigII->locate(pfqs[j], lsts[j], sd->left);
		    tfqs[j] = sd->sigII->pfq + sd->sigII->pfqnum;
		    agap[j] = -gaps[j]->gps * c_step;
		} else {
		    pfqs[j] = tfqs[j] = 0; lsts[j] = 0;
		}
	    }
	}

	int	z = 0, gph = 0, phs = 0, intlen = 0;
	if (htl == 3) prmode = Row_None;
	do {
	    if (algmode.nsa == ALN_FORM && OutPrm.SkipLongGap) {
		int	jj = 0;			//	skip long insert
		int	gap = INT_MAX;
		for (int j = 0; j < seqnum; ++j) {
		    int	upr = gp[j]->gps - gaps[j]->gps;
		    if (gp[j]->gln > 0) upr += gp[j]->gln;
		    if (gap > upr) {
			gap = upr;
			jj = j;
		    }
		}
		int	pos = gp[jj]->gps - gaps[jj]->gps; 	// print one line
		int	upr = (gap - z - OutPrm.EijMergin) / OutPrm.lpw * OutPrm.lpw;
		if (gap < INT_MAX && z - pos > OutPrm.EijMergin && upr > 0) {
		    for (int j = k = 0; j < seqnum; k += seqs[j++]->many) {
			if (j == jj) {
			    if (seqs[j]->isprotein())
				gph = phs = (phs + upr) % htl;
			} else {
			    int	inc = upr;
			    if (seqs[j]->isprotein() && htl == 3) {
				int	r = upr % 3;
				inc /= 3;
				if (r && phs != 2 && (r + phs) > 1) ++inc;
			    }
			    wkr[j] += inc * seqs[j]->many;
			    for (int i = 0; i < seqs[j]->many; ++i)
				nbr[k + i] += inc;
			}
		    }
		    z += upr;
		    globpt += upr;
		    gpos += upr;
		    fprintf(out_fd, "\n;; skip %d nt\'s\n", upr);
		    continue;
		}
	    }
	    vcopy(left, nbr, rows);
	    int	cpos = (z + 1) * c_step;
	    int	pphs = phs;
	    for (int clm = 0; clm < OutPrm.lpw; ++clm, ++z, ++gpos, cpos += c_step, pphs = phs) {
		int	reij = 0;
		if (exon) {
		    while (exon->right <= gpos) {intlen = exon[1].left - exon->right; ++exon;}
		    if ((exon->right == gpos + 1) && phs == 1) reij = 3;
		}
		if (active) active = seqnum;
		for (int j = k = 0; j < seqnum; k += seqs[j++]->many) {
		    Seq*&	sd = seqs[j];
		    bool	prot = sd->isprotein() && htl == 3;
	loop:	    int	pos = gp[j]->gps - gaps[j]->gps;
		    bool	neog = gaps_intr(gp[j]);
		    int	gap = neog? gp[j]->gln: 0;
		    if (!active) {
			for (int i = 0; i < sd->many; ++i)
			    image[k+i][clm] = BLANK;
		    } else if (pos > z) {
			if (prot) {
			    for (int i = 0; i < sd->many; ++i) {
				if (phs == 1) {
				    image[k+i][clm] = *wkr[j] == nil_code? gap_code: *wkr[j];
				    if (IsntGap(*wkr[j])) ++nbr[k+i];
				    ++wkr[j];
				} else
				    image[k+i][clm] = BLANK;
			    }
			    gph = phs = (phs + 1) % htl;
			} else {
			    for (int i = 0; i < sd->many; ++i) {
				image[k+i][clm] = *wkr[j] == nil_code? gap_code: *wkr[j];
				if (j == gene) {
				    int	df = gpos - exon->left;
				    if (df < 0) image[k+i][clm] |= INTRONBIT;
				    else if (intlen > 2 && df <= 1 && pphs == 1) reij = 2 - df;
				    if (markeij && reij) image[k+i][clm] |= (reij << 5);
				}
				if (IsntGap(*wkr[j])) ++nbr[k+i];
				++wkr[j];
			    }
			}
			if (markeij && pfqs[j]) {
			    int	niis = 0;
			    while ((pfqs[j] + niis < tfqs[j]) &&
				pfqs[j][niis].pos + agap[j] + (htl == 3? 1 - pfqs[j][niis].pos % 3: 0) < cpos) 
				    ++niis;
			    while (niis--) markiis(k, j, clm);
			}
		    } else if (pos + gap > z || !neog) {
			if (prot) {
			    for (int i = 0; i < sd->many; i++) {
				if (gph == 1 && (!exon || gpos > exon->left))
				    image[k+i][clm] = gap_code;
				else
				    image[k+i][clm] = BLANK;
			    }
			    gph = (gph + 1) % htl;
			} else {
			    for (int i = 0; i < sd->many; i++)
				image[k+i][clm] = (!exon || gpos >= exon->left)?
				    gap_code: BLANK;
			}
			if (!neog && !--active && !clm--) return;
		    } else {
			agap[j] += gp[j]->gln * c_step;
			gp[j]++;
			goto loop;
		    }
		}
	    }
	    if (out_fd)  prnt_aln(left, nbr);
	    globpt += OutPrm.lpw;
	} while(active);
}

void pralnseq(GAPS* gaps[], Seq* seqs[], int seqnum)
{
	if (!OutPrm.lpw) return;
	PrintAln	praln(gaps, seqs, seqnum);
	praln.printaln();
}

static void put_SigII(PFQ* pfq, PFQ* tfq, int* lst, int numlst)
{
	int	lwd = OutPrm.lpw >= 10? OutPrm.lpw - 4: 56;
	char	str[MAXL];

	fprintf(out_fd, ";B %ld %d\n", (long) (tfq - pfq), numlst);
	int	c = 0;
	for ( ; pfq < tfq - 1; ++pfq) {
	    if (c == 0) fputs(";b", out_fd);
	    sprintf(str + c, " %d %d,", pfq->pos, pfq->num);
	    if ((c = strlen(str)) > lwd) {
		fputs(str, out_fd);
		fputs("\n", out_fd);
		c = 0;
	    }
	}
	if (c)	fputs(str, out_fd);
	else	fputs(";b", out_fd);
	fprintf(out_fd, " %d %d\n", pfq->pos, pfq->num);
	if (!lst || !numlst) return;
	int*	tst = lst + numlst;
	for (c = 0; lst < tst - 1; ++lst) {
	    if (c == 0) fputs(";m", out_fd);
	    sprintf(str + c, " %d", *lst + 1);
	    if ((c = strlen(str)) > lwd) {
		fputs(str, out_fd);
		fputs("\n", out_fd);
		c = 0;
	    }
	}
	if (c)	fputs(str, out_fd);
	else	fputs(";m", out_fd);
	if (lst < tst) fprintf(out_fd, " %d\n", *lst + 1);
}

void Seq::putSigII()
{
	PFQ*	pfq = sigII? sigII->pfq: 0;
	int*	lst = sigII? sigII->lst: 0;

	if (!pfq) {
	    fputs(";B 0 0\n", out_fd);
	    return;
	}
	int	step = isprotein()? 3: 1;
	int	upto = step * (left + base_);
	for ( ; pfq->pos < upto; ++pfq)
	    if (lst) lst += pfq->num;
	upto = step * (right + base_);
	PFQ*	wfq = pfq;
	int	n = 0;
	for ( ; wfq->pos < upto; ++wfq)
	    if (lst) n += wfq->num;
	put_SigII(pfq, wfq, lst, n);
}

void SigII::putSigII()
{
	if (!step || !pfqnum) return;
	put_SigII(pfq, pfq + pfqnum, lst, lstnum);
}

int print2(Seq* seqs[], GAPS* gaps[], double fscr, 
	Gsinfo* GsI, int nbr, int ttl, int skip)
{
	GAPS*	gp = gaps[0];
	int	lpw = OutPrm.lpw;
	INT	BlkSz = OutPrm.BlkSz;
	Simmtx*	sm = getSimmtx(0);

	if (!out_fd || !gp) return (ERROR);
	while (gaps_intr(gp)) ++gp;
	int	span = gp->gps - gaps[0]->gps;
	if (!span) return (ERROR);
	fscr /= alprm.scale;
	PrintAln	praln(gaps, seqs, 2);
	int	many;
	switch (prmode) {
	    case Form_Phylp:
		many = seqs[0]->many + seqs[1]->many;
		fprintf(out_fd, "%5d %5d\n", many, span);
		break;
	    case Form_GDE:
		many = seqs[0]->many + seqs[1]->many;
		fprintf(out_fd, "%5d %5d  %s vs %s\n", many, span, 
			seqs[0]->sqname(), seqs[1]->sqname());
		for (int i = 0; i < seqs[0]->many; ++i)
		    if (seqs[0]->sname) fprintf(out_fd, "%-10s\n", (*seqs[0]->sname)[i]);
		    else	fprintf(out_fd, "Seq.%-6d\n", i);
		for (int i = 0; i < seqs[1]->many; ++i)
		    if (seqs[1]->sname) fprintf(out_fd, "%-10s\n", (*seqs[1]->sname)[i]);
		    else	fprintf(out_fd, "Seq.%-6d\n", seqs[0]->many + i);
		break;
	    case Form_GCG:
		gcg_out(gaps, seqs, 2);
		break;
	    case Form_CLW: break;
	    default:
		if (skip < 1)	fphseqs(seqs, 2);
		if (OutPrm.ColorEij) break;
		sm->fparam(out_fd);
		if (skip == -1)	return (OK);
//		if (skip < 2)	fpavsd(GsI->fstat.val);
		if (skip == -2)	return (OK);
		if (skip < 3) {
		    FTYPE	percent = GsI->fstat.mch + 
				GsI->fstat.mmc + GsI->fstat.unp; 
		    percent = 100. * GsI->fstat.mch / percent;
		    fprintf(out_fd, "Score = %5.1lf (%5.1f), %.1f (=), %.1f (#),",
			GsI->fstat.val, fscr, GsI->fstat.mch, GsI->fstat.mmc);
		    fprintf(out_fd, " %.1f (g), %.1f (u), (%5.2f %%)\n",
			GsI->fstat.gap, GsI->fstat.unp, percent);
		}
		if (skip == -3)	return (OK);
#if USE_WEIGHT
		if (skip < 4 && seqs[0]->weight && seqs[1]->weight) {
		    seqs[0]->fpweight();
		    seqs[1]->fpweight();
		}
#endif
		if (skip == -4)	return (OK);
		if (GsI->sigII) GsI->sigII->putSigII();
		if (skip < 5) fprintf(out_fd, "ALIGNMENT   %d / %d\n", nbr, ttl);
		break;
	}
	praln.printaln();
	fputs("\n\n", out_fd);
	if (prmode == Form_GCG) {
	    setgapchar(seqs[0], DEF_GAP_CHAR);
	    OutPrm.lpw = lpw;
	    OutPrm.BlkSz = BlkSz;
	}
	return (OK);
}

static void nexus_head(FILE* fd, Seq* sd)
{
	fputs("#NEXUS\n\n", fd);
	fputs("begin data;\n", fd);
	fprintf(fd, "\tdimensions ntax=%d nchar=%d;\n", sd->many, sd->len);
	fputs("\tformat datatype=protein missing=? gap=- matchchar=.;\n\n", fd);
	fputs("\tmatrix\n", fd);
}

void Seq::typeseq(FILE* fd, bool in_line)
{
	if (!many) return;
	if (!fd) fd = out_fd;
	else gswap(fd, out_fd);
	if (in_line) {
	    CHAR*	t = at(right);
	    for (CHAR*	p = at(left); p < t; ++p) 
		putc(code->decode[*p], out_fd);
	} else if (many == 1) {
	    listseq(0);
	} else if (prmode >= Form_Serial) {
	    if (OutPrm.noseqline)
		fprintf(out_fd, "%5d %5d\t%s\n", many, len, sqname());
	    else if (out_form->FormID == NEXUS)
		nexus_head(fd, this);
	    for (int j = 0; j < many; ++j) listseq(j, true);
	    if (out_form->FormID == NEXUS) fputs(";\nend;\n", fd);
	} else {
	    GAPS	gap[2] = {{left, 0}, {right, gaps_end}};
	    GAPS*	gpp[2] = {gap, gap + 1};
	    Seq*	tmp[1] = {this};
	    PrintAln	praln(gpp, tmp, 1);
	    if (prmode <= Form_Native) {
		if (OutPrm.ColorEij) {
		    printf(">%s\n", sqname());
		} else {
		    fphseqs(tmp, 1);
#if USE_WEIGHT
		    if (weight && OutPrm.printweight) fpweight();
#endif
		    if (sigII) putSigII();
		}
	    } else if (prmode == Form_Phylp) {
		fprintf(out_fd, "%5d %5d\n", many, len);
	    } else if (prmode == Form_GCG) {
		gcg_out(gpp, tmp, 1);
	    }
	    if (OutPrm.lpw) praln.printaln();
	}
}
