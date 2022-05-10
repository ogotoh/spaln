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
#include "gsinfo.h"

static	SeqDb*	out_form = 0;
static	SeqDb*	fst_form = 0;

static	const	char	BLANK = 0;
static	const	char	DEF_GAP_CHAR = '-';
static	const	char	GCG_GAP_CHAR = '.';
#if USE_WEIGHT
static	const	int	NO_WGHT = 5;
#endif
static	const	int	GCG_LINELENGTH = 50;
static	const	int	GCG_SeqBlkSz = 10;
static	const	INT	INTRONBIT = 0x80;
static	const	INT	RES_BITS = 0x1f;
static	const	char*	RGB_RED = "255,0,0";
static	const	char*	RGB_BLUE = "0,255,255";
static	const	char*	unknown_taxon = "K_Pp";

static	void	setgapchar(const Seq* sd, int c);
static	int	cdsForm(const RANGE* rng, const Seq* sd, int mode, FILE* fd);
static	void	put_SigII(const PFQ* pfq, const PFQ* tfq, const int* lst, 
			int numlst, FILE* fd);
static	int	chempro(const int* cmp, int ii);
static	int	logonuc(const int* cmp, int ii, const SEQ_CODE* defcode);
static	int	csym(const int cmp[], int rows, const SEQ_CODE* defcode);
static	int	SeqGCGCheckSum(const char* seq, int len);
static	int	checksum(int* checks, const GAPS* gaps, const Seq* sd);
static	void	gcg_out(const GAPS** gaps, Seq* seqs[], int seqnum, FILE* fd);

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

int setlpw(const int& lpwd)
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

FILE* setup_output(int omode, const char* def_fn, bool setup_out_fd)
{
	out_fd = stdout;
	if (setup_out_fd) {
	    if (!(OutPrm.out_file && *OutPrm.out_file) && def_fn)
		OutPrm.out_file = def_fn;
	    if (OutPrm.out_file && !is_dir(OutPrm.out_file)) {
		if (!(out_fd = wfopen(OutPrm.out_file, "w")))
		    fatal("Can't write to %s\n", OutPrm.out_file);
	    }
	}
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

void Gsinfo::repalninf0(const SKL* skl, Seq* seqs[]) const
{
	Seq*&	a = seqs[0];
	Seq*&	b = seqs[1];
const	FSTAT*	fst = &fstat;
	int	rsv_mhits = 0;

	swap(rsv_mhits, distPrm.corr_mhits);
	if (OutPrm.ColorEij) {
	    put_stat(fd, fst);
	    int 	apfq = a->sigII? a->sigII->pfqnum: 0;
	    int 	bpfq = b->sigII? b->sigII->pfqnum: 0;
	    int		dnm = apfq + bpfq;
	    float	ncmn = sigII? sigII->n_common(): 0;
	    float	spbdist = dnm? (1. - 2 * ncmn / dnm): 0;
	    fprintf(fd, "%6.2f\t%6d\t%6d\t", 100. * spbdist, apfq, bpfq);
	} else	put_stat(fd, this);
	fprintf(fd, "%6.1f\t%d %d %c\t%d %d %c\t%s\t%s\n",
	    (double) fst->val / (alprm.scale * a->many * b->many),
	    a->SiteLe() + 1, a->SiteRe(), a->Strand(),
	    b->SiteLe() + 1, b->SiteRe(), b->Strand(),
	    a->sqname(), b->sqname());
	swap(rsv_mhits, distPrm.corr_mhits);
}

/* 2 lines classic SKL format */

void Gsinfo::repalninf1(const SKL* skl, Seq* seqs[]) const
{
	int	nn = skl->n;
 	Seq*&	a = seqs[0];
 	Seq*&	b = seqs[1];
const 	FSTAT&	fst = fstat;

	a->fphseq(fd);
	putc(' ', fd);
	b->fphseq(fd);
	fprintf(fd, "  %d  %.2f\n", nn, 
	    (double) fst.val / (alprm.scale * a->many * b->many));
	for (++skl; nn--; skl++)
	    fprintf(fd, "%d %d ", a->SiteNo(skl->m), b->SiteNo(skl->n));
	putc('\n', fd);
}

/* sugar-like format */

void Gsinfo::repalninf2(const SKL* skl, Seq* seqs[]) const
{
	int	nn = skl->n;
 	Seq*&	a = seqs[0];
 	Seq*&	b = seqs[1];
const	char*	anm = a->sqname();
const	char*	bnm = b->sqname();
	char	asn = a->Strand();
	char	bsn = b->Strand();
const 	Simmtx*	sm = getSimmtx(0);

	for (++skl; nn-- > 1; skl++) {
	    int	dm = skl[1].m - skl->m;
	    int	dn = skl[1].n - skl->n;
	    if (dm && dn) {
		const	CHAR*	as = a->at(skl->m);
		const	CHAR*	bs = b->at(skl->n);
		int	i = min(dm, dn);
		double	val = 0;
		for ( ; i--; ++as, ++bs)
		    val += sm->mtx[*as][*bs];
		val /= alprm.scale;
		fprintf(fd, "sugar: %s %d %d %c %s %d %d %c %d M %d %d\n", 
		    bnm, b->SiteNo(skl->n), b->SiteNo(skl[1].n), bsn,
		    anm, a->SiteNo(skl->m), a->SiteNo(skl[1].m), asn,
		    (int) val, dn, dm);
	    }
	}
}

/* psl-like format */

void Gsinfo::repalninf3(const SKL* skl, Seq* seqs[]) const
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
const 	CHAR*	as;
const 	CHAR*	bs;
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
	fprintf(fd, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c%c\t",
	    mch, mmc, rep, noN, qgap, qunp, tgap, tunp, a->Strand(), b->Strand());
	fprintf(fd, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t",
	    a->sqname(), a->len, a->SiteNz(a->left), a->SiteNz(a->right), 
	    b->sqname(), b->len, b->SiteNz(b->left), b->SiteNz(b->right), diag);
	for (int i = 0; i < diag; ++i)
	    fprintf(fd, "%d,", xyl[i].p);
	putc('\t', fd);
	for (int i = 0; i < diag; ++i)
	    fprintf(fd, "%d,", xyl[i].m);
	putc('\t', fd);
	for (int i = 0; i < diag; ++i)
	    fprintf(fd, "%d,", xyl[i].n);
	putc('\n', fd);
	delete[] xyl;
}

/* 1 line compact XYL format (x, y, len) * n */

void Gsinfo::repalninf4(const SKL* skl, Seq* seqs[]) const
{
	int	nn = skl->n;
 	Seq*&	a = seqs[0];
 	Seq*&	b = seqs[1];
const 	FSTAT&	fst = fstat;

	fprintf(fd, "XYL: %s %d %d %c %s %d %d %c %.1lf : ",
	    a->sqname(), a->SiteLe() + 1, a->SiteRe(), a->Strand(),
	    b->sqname(), b->SiteLe() + 1, b->SiteRe(), b->Strand(),
	    (double) fst.val/ (alprm.scale) * a->many * b->many);
	for (++skl; nn-- > 1; skl++) {
	    int	dm = skl[1].m - skl->m;
	    int	dn = skl[1].n - skl->n;
	    if (dm && dn) fprintf(fd, " %d %d %d",
		a->SiteNo(skl->m), b->SiteNo(skl->n), dm);
	}
	fprintf(fd, "\n");
}

/* 1 line exon boundary format */

void Gsinfo::repalninf5(const SKL* skl, Seq* seqs[]) const
{
	int	nn = skl->n;
 	Seq*&	b = seqs[1];

	b->fphseq(fd, 0);
	for (++skl; nn--; skl++)
	    fprintf(fd, " %d", b->SiteNo(skl->n));
	putc('\n', fd);
}

/* 2 line compact XYL format (x, y, len) * n */

void Gsinfo::repalninf6(const SKL* skl, Seq* seqs[]) const
{
	int	nn = skl->n;
 	Seq*&	a = seqs[0];
 	Seq*&	b = seqs[1];
const 	FSTAT&	fst = fstat;

	fprintf(fd, "XYL2: %s %d %d %c %s %d %d %c %7.1lf ",
	    a->sqname(), a->SiteLe() + 1, a->SiteRe(), a->Strand(),
	    b->sqname(), b->SiteLe() + 1, b->SiteRe(), b->Strand(),
	    (double) fst.val/ (alprm.scale) * a->many * b->many);
	fprintf(fd, "%6.2f %d %d %d %d %d\n",
	    100. * fst.mch / (fst.mch + fst.mmc + fst.gap),
	    (int) fst.mch, (int) fst.mmc, (int) fst.gap, (int) fst.unp,
	    nn);
	for (++skl; nn-- > 1; skl++) {
	    int	dm = skl[1].m - skl->m;
	    int	dn = skl[1].n - skl->n;
	    if (dm && dn) fprintf(fd, " %d %d %d",
		a->SiteNo(skl->m), b->SiteNo(skl->n), dm);
	}
	fprintf(fd, "\n");
}

static void setgapchar(const Seq* sd, int c)
{
	if (sd->isprotein())
	    amino[gap_code] = amino[nil_code] = c;
	else
	    nucl[gap_code] = nucl[nil_code] = c;
}

/*	CDS format	*/

static int cdsForm(const RANGE* rng, const Seq* sd, int mode, FILE* fd)
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

const	RANGE*	wrng = rng = fistrng(rng);
	if (mode == 1) {
	  switch(out_form->FormID) {
	    case GenBank:
		fprintf(fd, "%s\n", out_form->ComLabel);
		if (out_form->ContSpc < strlen(blank))
		    blank[out_form->ContSpc] = '\0';
		ModifyMssg = blank;
		break;
	    case EMBL:
		fputs(emblxx, fd);
		ModifyMssg = "CC   ";
		break;
	  }
	} else if (mode == 2) {
	    switch(out_form->FormID) {
	      case GenBank:
		head = blank;
		fputs(out_form->FeaLabel, fd);
		putc('\n', fd);
		strcpy(str, gftcds);
		break;
	      case EMBL:
	    	fputs(emblxx, fd);
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
		    fputs(str, fd);
		    putc('\n', fd);
		    strcpy(str, head);
		    c = lb;
		}
	    } else ++indel;
	    if (mode != 1) continue;
	    if (intv <= 2 && intv > 0) {
		fputs(ModifyMssg, fd);
		fprintf(fd, "Deleted %d chars at %d\n", intv, r);
	    } else if (intv < 0) {
		fputs(ModifyMssg, fd);
		fprintf(fd, "Insert %d chars at %d\n", -intv, r);
	    }
	}
	if (mode == 2) {
	    fputs(str, fd);
	    putc(')', fd);
	    if (cmpl) putc(')', fd);
	    putc('\n', fd);
	}
	if (out_form->FormID == GenBank && mode == 1)
	    blank[out_form->ContSpc] = ' ';
	return (indel);
}

void GBcdsForm(const RANGE* rng, const Seq* sd, FILE* _fd)
{
	RANGE	svr;

	FILE*	fd = _fd? _fd: out_fd;
	sd->saverange(&svr);
	if (out_form->FormID != GenBank && out_form->FormID != EMBL) {
	    sd->left = 0; sd->right = sd->len;
	}
	if (cdsForm(rng, sd, 0, fd))
	    cdsForm(rng, sd, 1, fd);
	cdsForm(rng, sd, 2, fd);
	sd->restrange(&svr);
}

static	const	char*	fmt3p = "%s\tALN\t%s\t%d\t%d\t%d\t%c\t%d\t";
static	const	char*	fmt3n = "%s\tALN\t%s\t%d\t%d\t%d\t%c\t.\t";

/* Phase 1 and 2 of Gff3 are reversed from those of native format */

void Gsinfo::Gff3Form(const Seq* gene, const Seq* qry) const
{
	int	rvg = gene->Strand();
	int	rvr = qry->Strand();
	int	cds = 0;
	int	hetero = qry->isprotein();
	int	len3 = (gene->len + 2) / 3 * 3;
const 	EISCR*	rng = eijnc->begin();
const 	EISCR*	wkr = rng;
const 	EISCR*	skp = rng;
const 	SKL*	wsk = skl;
const 	SKL*	lst = wsk + (hetero? wsk->n: 0);
const 	SKL*	prv;
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
	if (!Gff3MID++) fputs("##gff-version\t3\n", fd);
	l = gene->SiteNo(0);
	r = gene->SiteNo(gene->len - 1);
	if (rvg == '-') swap(l, r);
	fprintf(fd, "##sequence-region\t%s %d %d\n", gene->sqname(), l, r);
	l = gene->SiteNo(rng->left);
	r = gene->SiteNo(wkr->right - 1);
	if (rvg == '-') swap(l, r);
	fprintf(fd, fmt3n, (*gene->sname)[0], "gene", l, r, 
	    (int) (scr / alprm.scale), rvg);
	fprintf(fd, fmt3g, Gff3MID, mname);
	fprintf(fd, fmt3n, (*gene->sname)[0], "mRNA", l, r, 
	    (int) (scr / alprm.scale), rvg);
	fprintf(fd, fmt3m, Gff3MID, Gff3MID, mname);

	wkr = rng;
	prv = ++wsk;
	while (neoeij(wkr)) {
	    if (wkr->iscr > NEVSEL) {
		l = gene->SiteNo(skp->left);
		r = gene->SiteNo(wkr->right - 1);
		if (rvg == '-') swap(l, r);
		if (hetero) {
		    fprintf(fd, fmt3, (*gene->sname)[0], ltype, l, r,
			(int) (wkr->escr / alprm.scale), rvg, (len3 - cds) % 3);
		} else {
		    fprintf(fd, fmt3, (*gene->sname)[0], ltype, l, r, 
			(int) (wkr->escr / alprm.scale), rvg);
		}
		cds += wkr->right - skp->left;
		fprintf(fd, fmt3c, ltype, ++Gff3CID, Gff3MID, mname);
		fprintf(fd, fmt3t, (*qry->sname)[0], qry->SiteNo(wkr->rleft),
		    qry->SiteNo(wkr->rright - 1), rvr);
		while (wsk <= lst && wsk->m <= wkr->right + 1) {
		    l = wsk->m - prv->m;	/* genome */
		    r = wsk->n - prv->n;	/* protein */
		    if (l && (r || l <= IntronPrm.llmt)) cds -= l % 3;
		    prv = wsk++;
		}
		skp = ++wkr;
	    } else if ((wkr++)->iscr > NEVSEL) {
		skp = wkr;
	    }
	}
}

void Gsinfo::Gff3PWA(const Seq* gene, const Seq* qry) const
{
const 	EISCR*	rng = eijnc->begin();
const 	EISCR*	wkr = rng;
const 	EISCR*	skp = rng;
const 	SKL*	wsk = skl;
const 	SKL*	lst = wsk + wsk->n;
const 	SKL*	prv = ++wsk;
	int	rvg = gene->Strand();
	int	rvr = qry->Strand();
	int	hetero = qry->isprotein();
	char	mname[MAXL];
	int	donor = 0;
const	char*	ltype;
static	int	Gff3MID = 0;
static	const	char*	fmt3c = "ID=match%05d;Name=%s;";
static	const	char*	fmt3t = "Target=%s %d %d %c;Gap=";

	swapskl(skl);
	if (hetero)	ltype = "nucleotide_to_protein_match";
	else		ltype = "cDNA_match";
	while (neoeij(wkr + 1)) ++wkr;
	int	l = gene->SiteNo(rng->left);
	int	r = gene->SiteNo(wkr->right - 1);
	sprintf(mname, "%s_%d", (*gene->sname)[0], (l + r) / 2000);
	if (!Gff3MID++) fputs("##gff-version\t3\n", fd);
	l = gene->SiteNo(0);
	r = gene->SiteNo(gene->len - 1);
	if (rvg == '-') swap(l, r);
	fprintf(fd, "##sequence-region\t%s %d %d\n", gene->sqname(), l, r);

	wkr = rng;
	while (neoeij(wkr)) {
	    if (wkr->iscr > NEVSEL) {
		l = gene->SiteNo(skp->left);
		r = gene->SiteNo(wkr->right - 1);
		if (rvg == '-') swap(l, r);
		fprintf(fd, fmt3n, (*gene->sname)[0], ltype, l, r, 
		    (int) (wkr->escr / alprm.scale), rvg);
		fprintf(fd, fmt3c, Gff3MID, mname);
		fprintf(fd, fmt3t, (*qry->sname)[0], qry->SiteNo(wkr->rleft),
		    qry->SiteNo(wkr->rright - 1), rvr);
		for ( ; wsk <= lst && wsk->n <= wkr->right + 1; prv = wsk++) {
		    int	dm = wsk->m - prv->m;
		    int	dn = wsk->n - prv->n;
		    r = skp->left - wsk->n;	/* acceptor */
		    if (!dm && donor && -1 <= r && r <= 1)
			continue;		 /* intron */
		    l = wkr->right - wsk->n;
		    donor = -1 <= l && l <= 1;
		    if (hetero) {
			if (!dn) {
			    if (dm) fprintf(fd, "I%d ", dm);
			    continue;
			}
			l = dn - dm * 3;
			if (l < 0) {
			    if (l > -3) dn += 2;
			    dm = dn / 3;
			}
			if (dm)		fprintf(fd, "M%d ", dm);
			if (l < 0) {
			    l = -l;
			    if (l % 3)	fprintf(fd, "R%d ", l);
			    else	fprintf(fd, "I%d ", l / 3);
			} else if (l)	{
			    if (l % 3)	fprintf(fd, "F%d ", l);
			    else	fprintf(fd, "D%d ", l / 3);
			}
		    } else {
			if (!dm && !dn)	continue;
			else if (!dn)	fprintf(fd, "I%d ", dm);
			else if (!dm )	fprintf(fd, "D%d ", dn);
			else		fprintf(fd, "M%d ", dm);
		    }
		}
		putc('\n', fd);
		skp = ++wkr;
	    } else if ((wkr++)->iscr > NEVSEL) {
		skp = wkr;
	    }
	}
	swapskl(skl);
}

#if PSL_FORM

/* psl-like format */

void Gsinfo::PslForm(const Seq* t, const Seq* q) const
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
const 	EISCR*	fst = eijnc->bigin();
const 	EISCR*	wkr = fst;
const 	SKL*	wsk = skl;
const 	SKL*	lst = wsk + wsk->n;
const 	SKL*	psk = ++wsk;
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
	    for (int i = 0; i < 3; ++i) fputs(PslHeader[i], fd);
	swapskl(skl);
	SKL	prv = *wsk;
	for ( ; ++wsk <= lst; psk = wsk) {
	    int	dm = wsk->m - psk->m;
	    int	dn = wsk->n - psk->n;
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
	fprintf(fd, "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%c%c\t",
	    mch, mmc, rep, noN, qgap, qunp, tgap, tunp,
	    t->Strand(), q->Strand());
	--wkr;
	tunp = t->SiteNz(wkr->right);
	qunp = q->SiteNz(wkr->rright);
	wkr = fst;
	tgap = t->SiteNz(wkr->left);
	qgap = q->SiteNz(wkr->rleft);
	if (q->inex.sens) swap(qgap, qunp);
	if (t->inex.sens) swap(tgap, tunp);
	fprintf(fd, "%s\t%d\t%d\t%d\t%s\t%d\t%d\t%d\t%d\t",
	    q->sqname(), q->len, qgap, qunp, 
	    t->sqname(), t->len, tgap, tunp, diag);
	if (!hetero && t->inex.sens & 2) {
	    tmp = xyl[diag];
	    for (int i = 0, rep = diag; i < rep; ++i, --rep)
		swap(xyl[i], xyl[rep]);
	    for (int i = 0; i < diag; ++i)
		xyl[i].p = xyl[i+1].p;
	}
	for (int i = 0; i < diag; ++i)
	    fprintf(fd, "%ld,", xyl[i].p);
	putc('\t', fd);
	for (int i = 0; i < diag; ++i)
	    fprintf(fd, "%d,", xyl[i].m);
	putc('\t', fd);
	for (int i = 0; i < diag; ++i)
	    fprintf(fd, "%d,", xyl[i].n);
	putc('\n', fd);
	delete[] xyl;
	swapskl(skl);
}

#endif

/* UCSC BED format */

void Gsinfo::BedForm(const Seq* gene, const Seq* qry) const
{
static	int	visit = 0;
static	const	char*	BedHeader = 
	"track name=Spaln description=\"%s\" useScore=1\n";
static	const	char*	BedFrom = "%s\t%d\t%d\t%s\t%d\t%c\t%d\t%d\t%s\t%d\t";

	if (!visit++) fprintf(fd, BedHeader, qry->sqname());
const 	EISCR*	fst = eijnc->begin();
const 	EISCR*	wkr = fst;
	int	gstart = wkr->left;
const 	EISCR*	prv = wkr;
	while (neoeij(wkr)) prv = wkr++;
	int	gend = prv->right;
	int	nexon = wkr - fst;
	int*	epos = new int[nexon];
	int*	elen = new int[nexon];
	wkr = fst;
	prv = wkr;
const 	EISCR*	skp = wkr;
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
	    swap(gstart, gend);
	    vreverse(epos, nexon);
	    vreverse(elen, nexon);
	}
	float	bedscr = rscr > 0? 1000. * scr / rscr: scr;
	const	char*	rgb = (gene->inex.sens & REVERS)? RGB_BLUE: RGB_RED;
	fprintf(fd, BedFrom, (*gene->sname)[0], gstart, gend, (*qry->sname)[0], 
	    min(1000, int(bedscr)), gene->Strand(), gstart, gend, rgb, nexon);
	eln = elen;
	for (int i = 0; i < nexon; ++i) fprintf(fd, "%d,", *eln++);
	putc('\t', fd);
	eps = epos;
	for (int i = 0; i < nexon; ++i) {
	    if (i < nexon - 1) fprintf(fd, "%d,", *eps++);
	    else	fprintf(fd, "%d\n", *eps++);
	}
	delete[] epos; delete[] elen;
}

/* similar to megablast -D 3 output format without E-value */
static	FILE*	fg = 0;
static	FILE*	fe = 0;
static	FILE*	fq = 0;
#if USE_ZLIB
static	gzFile	gzfg = 0;
static	gzFile	gzfe = 0;
static	gzFile	gzfq = 0;
#endif
void Gsinfo::ExonForm(const Seq* gene, const Seq* qry, int mode) const
{
	int	cds = 0;
	int	mch = 0;
	bool	binary = mode == BIN_FORM;
	float	scale = alprm.scale;
	if (gene->exin && gene->exin->fact) scale *= gene->exin->fact;
	float	fscr = scr / alprm.scale;
	VTYPE	iscr = 0;
const	int	bbt = qry->isprotein()? 3: 1;
const 	EISCR*	fst = eijnc->begin();
const 	EISCR*	wkr = fst;
const 	EISCR*	prv = wkr;
const 	EISCR*	skp = wkr;
	char	str[LINE_MAX];
const 	char*	convt = gene->isdrna()? nucl: ncodon;
	ExonRecord	er;
static	GeneRecord	gr;
static	int	visit = 0;
static	int	qryidx = 0;
static	const	char	tmp_intends[] = "  .  ";
static	const	char	fmt[] = "%s\t%s\t%7.2f\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t"
			"%7d\t%7.1f\t%7d\t%7.1f\t%7.2f\t%7.2f %2d %d %2d %d %s\n";
static	const	char	tfmt[] = "@ %s %c ( %d %d ) %s [%d:%d] ( %d %d ) S: %.1f =: %.1f C: %.1f "
//			"T#: %d T-: %d B#: %d B-: %d X: %d Nexn: %d\n";
			"T#: %d T-: %d B#: %d B-: %d X: %d Nexn: %d HC: %4.1f\n";
static	const	char	efmt[] = "Cant't write to gene record %s: # %ld\n";

	char	intends[8];
	strcpy(intends, tmp_intends);
	gr.nexn = gr.mmc = gr.unp = gr.bmmc = gr.bunp = gr.ng = 0;
	memset(&er, '\0', sizeof(ExonRecord));
	if (visit++ == 0) {
	    if (binary) {
		strcpy(str, prefix);
		char*	bdy = str + strlen(str);
		if (is_gz(bdy)) {
		    char*	dot = strrchr(bdy, '.');
		    *dot = '\n';
		}
		strcpy(bdy, grext);
#if USE_ZLIB
		if (OutPrm.gzipped) {
		    strcat(str, gz_ext);
		    if (!(gzfg = wgzopen(str, "wb"))) fatal(efmt, str, 0);
		} else 
#endif
		    if (!(fg = wfopen(str, "wb"))) fatal(efmt, str, 0);
		strcpy(bdy, erext);
#if USE_ZLIB
		if (OutPrm.gzipped) {
		    strcat(str, gz_ext);
		    if (!(gzfe = wgzopen(str, "wb"))) fatal(efmt, str, 0);
		} else 
#endif
		    if (!(fe = wfopen(str, "wb"))) fatal(efmt, str, 0);
		strcpy(bdy, qrext);
#if USE_ZLIB
		if (OutPrm.gzipped) {
		    strcat(str, gz_ext);
		    if (!(gzfq = wgzopen(str, "wb"))) fatal(efmt, str, 0);
		    if ((fputs(dbs_dt[0]->dbsid, gzfq) == EOF) || (fputc('\0', gzfq) == EOF))
		    fatal(efmt, qrext, gr.Nrecord);
		} else 
#endif
		    if (!(fq = wfopen(str, "wb"))) fatal(efmt, str, 0);
		if (fq && ((fputs(dbs_dt[0]->dbsid, fq) == EOF) || (fputc('\0', fq) == EOF)))
		    fatal(efmt, qrext, gr.Nrecord);
		++qryidx;	// o-th q-recode contains the database name
	    } else {
		fputs("# rID\t  gID\t   %id\t  ExonL\t MisMch\t Unpair\t "
		"ref_l\t  ref_r\t  tgt_l\t  tgt_r\t eScore\t IntrnL\t "
		"iScore\t Sig3/I\t Sig5/T  # -  X P DiNuc\n", fd);
	    }
	}
	while (neoeij(wkr)) {
	    if (wkr->iscr > NEVSEL) {	// normal
		if (gr.nexn) {
		    gr.bmmc += prv->mmc3 + wkr->mmc5;
		    gr.bunp += prv->unp3 + wkr->unp5;
		}
		int	exon = wkr->right - skp->left;
		int	rlen = (wkr->rright - skp->rleft) * qry->many + wkr->unp;
		cds += exon;
		mch += wkr->mch;
		gr.mmc += wkr->mmc;
		gr.unp += wkr->unp;
		++gr.nexn;
		er.Pmatch = rlen? 100. * wkr->mch / rlen: 0;
		er.Elen   = exon;
		er.Rleft  = qry->SiteNo(skp->rleft);
		er.Rright = qry->SiteNo(wkr->rright - 1);
		er.Gleft  = gene->SiteNo(skp->left);
		er.Gright = gene->SiteNo(wkr->right - 1);
		er.Bmmc   = prv->mmc3 + wkr->mmc5;
		if (prv->unp3 % bbt || wkr->unp5 % bbt) er.Bunp = 9;
		else	er.Bunp = (prv->unp3 + wkr->unp5) / bbt;
		er.Escore = (float) wkr->escr / scale;
		er.Iscore = (float) iscr / scale;
		er.Sig3   = (float) wkr->sig3 / scale;
		er.Sig5   = (float) wkr->sig5 / scale;
		if (binary) {
		    er.Nmmc = wkr->mmc;
		    er.Nunp = wkr->unp;
		    if (fe && fwrite(&er, sizeof(ExonRecord), 1, fe) != 1)
			fatal(efmt, erext, gr.Nrecord);
#if USE_ZLIB
		    if (gzfe && fwrite(&er, sizeof(ExonRecord), 1, gzfe) != 1)
			fatal(efmt, erext, gr.Nrecord);
#endif
		} else {
		    fprintf(fd, fmt, 
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
	    } else {		// frame shift
		gr.ng++;
		if ((wkr++)->iscr > NEVSEL) skp = wkr;
		if (!neoeij(wkr)) break;
		er.miss = wkr->left - skp->right;
	    }
	}
	cds /= bbt;
	er.Ilen = qry->right - qry->left;
	gr.Gstart = gene->SiteNo(fst->left);
	gr.Gend   = gene->SiteNo(prv->right - 1);
	gr.Rstart = qry->SiteNo(qry->left);
	gr.Rend   = qry->SiteNo(qry->right - 1);
	gr.Gscore = fscr;
	gr.Pmatch = 100. * mch / er.Ilen;
	gr.Pcover = 100. * (mch + gr.mmc) / er.Ilen;
	gr.Csense = gr.Gstart > gr.Gend;
	if (binary) {
	    gr.Cid = gene->did;
	    gr.Rid = qryidx++;
	    gr.Rsense = qry->inex.sens;
	    gr.Rlen = qry->len;
	    if (qry->isprotein() && gr.ng == 0) gr.ng = -1;
	    if (fg && fwrite(&gr, sizeof(GeneRecord), 1, fg) != 1)
		fatal(efmt, grext, gr.Nrecord);
#if USE_ZLIB
	    if (gzfg && fwrite(&gr, sizeof(GeneRecord), 1, gzfg) != 1)
		fatal(efmt, grext, gr.Nrecord);
#endif
	    if (fq && (fputs((*qry->sname)[0], fq) == EOF || fputc('\0', fq) == EOF))
		fatal(efmt, qrext, gr.Nrecord);
#if USE_ZLIB
	    if (gzfq && (fputs((*qry->sname)[0], gzfq) == EOF || fputc('\0', gzfq) == EOF))
		fatal(efmt, qrext, gr.Nrecord);
#endif
	} else {
	    int	n = 0;
	    for (int j = 0; j < gene->CdsNo; ++j)
		n += gene->jxt[j].jlen;
	    float hc = 100. * n / (qry->right - qry->left);
	    fprintf(fd, tfmt, (*gene->sname)[0], gr.Csense? '-': '+', 
	    gr.Gstart, gr.Gend,
	    (*qry->sname)[0], qry->many, qry->len, gr.Rstart, gr.Rend,
	    gr.Gscore, gr.Pmatch, gr.Pcover, 
	    gr.mmc, gr.unp, gr.bmmc, gr.bunp, gr.ng, gr.nexn, hc);
	}
	gr.Nrecord += gr.nexn;
}

void closeGeneRecord()
{
	if (fg) {fclose(fg); fg = 0;}
	if (fe) {fclose(fe); fe = 0;}
	if (fq) {fclose(fq); fq = 0;}
#if USE_ZLIB
	if (gzfg) {fclose(gzfg); gzfg = 0;}
	if (gzfe) {fclose(gzfe); gzfe = 0;}
	if (gzfq) {fclose(gzfq); gzfq = 0;}
#endif
}

void Gsinfo::IntronForm(const Seq* gene, const Seq* qry) const
{
	char	rvc = gene->Strand();
const 	EISCR*	fst = eijnc->begin();
const 	EISCR*	wkr = fst;
const 	EISCR*	prv = wkr;
const 	char*	convt = gene->isdrna()? nucl: ncodon;
	int	cds = wkr->right - wkr->left;
	bool	isnuc = qry->isdrna();
static	int	visit = 0;
static	const char	tmp_intends[] = "   ..   ";
static	const char	tfmt[] = "@ %s %c ( %d %d ) %s [%d:%d] ( %d %d )\n";
static	const char	fmt[] = "%s\t%c %9d %9d  %d  %9d %9d\t%s\t%7d\t%7d\t%7d\t%7d\t%7d\t%7d\t %s\n";

	char	intends[10];
	strcpy(intends, tmp_intends);
	if (visit++ == 0)
	    fputs("# gID\tdir   Donor  Acceptor Phs     tgt_5     tgt_3\trefID\t  ref_l\t  ref_r\t  Match\tMisMach\t Unpair\tIntronL\tIntronEnd\n", fd);
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
		fprintf(fd, fmt, (*gene->sname)[0], rvc, gene->SiteNo(prv->right),
		gene->SiteNo(wkr->left - 1), qry->isdrna()? cds % 3: (3 - prv->phs) % 3,
		gene->SiteNo(prv->left), gene->SiteNo(wkr->right - 1),
		(*qry->sname)[0], qry->SiteNo(prv->rleft), 
		qry->SiteNo(wkr->rright - 1), mch, mmc, unp, intv, intends);
	    }
	    prv = wkr;
	    cds += wkr->right - wkr->left;
	}
	GeneRecord      gr;
	gr.Gstart = gene->SiteNo(fst->left);
	gr.Gend   = gene->SiteNo(prv->right) - 1;
	gr.Rstart = qry->SiteNo(qry->left);
	gr.Rend   = qry->SiteNo(qry->right) - 1;
	gr.Csense = gr.Gstart > gr.Gend;
	fprintf(fd, tfmt, (*gene->sname)[0], gr.Csense? '-': '+',
	    gr.Gstart, gr.Gend, (*qry->sname)[0], qry->many, qry->len, gr.Rstart, gr.Rend);
}

void Gsinfo::CigarForm(const Seq* gen, const Seq* qry) const
{
	if (!cigar || !cigar->rec || !skl) {
	    prompt("No cigar records !\n");
	    return;
	}
const 	SKL*	fst = skl + 1;
const 	SKL*	lst = skl + skl->n;
	fprintf(fd, "cigar: %s %d %d %c %s %d %d %c %d",
	    qry->sqname(), qry->SiteNm(fst->n), 
	    qry->SiteNm(lst->n), qry->inex.sens? '-': '+',
	    gen->sqname(), gen->SiteNm(fst->m), 
	    gen->SiteNm(lst->m), gen->inex.sens? '-': '+',
	    (int) fstat.val);
const 	CIGAR*	tm = cigar->rec + cigar->size();
	for (const CIGAR* wk = cigar->rec; wk < tm; ++wk)
	    fprintf(fd, " %c %d", wk->ope, wk->len);
	putc('\n', fd);
}

void Gsinfo::VulgarForm(const Seq* gen, const Seq* qry) const
{
	if (!vlgar || !vlgar->rec || !skl) {
	    prompt("No vlgar records !\n");
	    return;
	}
	const	SKL*	fst = skl + 1;
	const	SKL*	lst = skl + skl->n;
	fprintf(fd, "vulgar: %s %d %d %c %s %d %d %c %d",
	    qry->sqname(), qry->SiteNm(fst->n), 
	    qry->SiteNm(lst->n), qry->inex.sens? '-': '+',
	    gen->sqname(), gen->SiteNm(fst->m), 
	    gen->SiteNm(lst->m), gen->inex.sens? '-': '+',
	    (int) fstat.val);
	const VULGAR*	tm = vlgar->rec + vlgar->size();
	for (const VULGAR* wk = vlgar->rec; wk < tm; ++wk)
	    fprintf(fd, " %c %d %d", wk->ope, wk->alen, wk->blen);
	putc('\n', fd);
}

void Gsinfo::SamForm(const Seq* gen, Seq* qry) const
{
const	Samfmt*	samf = samfm;
	if (!samf || !samf->rec || !skl) {
	    prompt("No sam records !\n");
	    return;
	}
	fprintf(fd, "%s\t%d\t%s\t%d\t%d\t", qry->sqname(),
	    samf->flag, gen->sqname(), gen->SiteNo(samf->pos), samf->mapq);
	const CIGAR*	tm = samf->rec + samf->size();
	for (const CIGAR* wk = samf->rec; wk < tm; ++wk)
	    fprintf(fd, "%d%c", wk->len, wk->ope);
	putc('\t', fd);
	fprintf(fd, "*\t%d\t%d\t", 
	   samf->tlen? gen->SiteNo(samf->pnext): 0, samf->tlen);
	swap(qry->left, samf->left);
	swap(qry->right, samf->right);
	if (gen->inex.sens) qry->comrev();
	qry->typeseq(0, true);
	if (gen->inex.sens) qry->comrev();
	swap(qry->left, samf->left);
	swap(qry->right, samf->right);
	fprintf(fd, "\t*\tNM:i:%d\tAS:i:%d\n", 
	    int(fstat.mmc + fstat.unp), int(scr / alprm.scale));
}

void Gsinfo::repalninf(Seq* seqs[], int mode, FILE* _fd) const
{
	fd = _fd? _fd: out_fd;
	switch (mode) {
	  case 0:  repalninf0(skl, seqs); break;
	  case 1:  repalninf1(skl, seqs); break;
	  case 2:  repalninf2(skl, seqs); break;
	  case 3:  repalninf3(skl, seqs); break;
	  case 4:  repalninf4(skl, seqs); break;
	  case 5:  repalninf5(skl, seqs); break;
	  case CIG_FORM: CigarForm(seqs[1], seqs[0]); return;
	  case VLG_FORM: VulgarForm(seqs[1], seqs[0]); return;
	  case SAM_FORM: 
	  case BAM_FORM: SamForm(seqs[1], seqs[0]); return;
	  default: repalninf6(skl, seqs); break;
	}
}

int Gsinfo::print2(Seq* seqs[], const GAPS** gaps, 
	int nbr, int ttl, int skip, FILE* _fd) const
{
const	GAPS*	gp = gaps[0];
	if (!gp) return (ERROR);
	fd = _fd? _fd: out_fd;
	int	lpw = OutPrm.lpw;
	INT	BlkSz = OutPrm.BlkSz;
	Simmtx*	sm = getSimmtx(0);

	while (gaps_intr(gp)) ++gp;
	int	span = gp->gps - gaps[0]->gps;
	if (!span) return (ERROR);
	double	sscr = double(scr) / alprm.scale;	// spliced alignment score
	double	tscr = double(fstat.val) / alprm.scale;	// transcript alignment score
	PrintAln	praln(gaps, seqs, 2, fd);
	int	many;
	switch (prmode) {
	    case Form_Phylp:
		many = seqs[0]->many + seqs[1]->many;
		fprintf(fd, "%5d %5d\n", many, span);
		break;
	    case Form_GDE:
		many = seqs[0]->many + seqs[1]->many;
		fprintf(fd, "%5d %5d  %s vs %s\n", many, span, 
			seqs[0]->sqname(), seqs[1]->sqname());
		for (int i = 0; i < seqs[0]->many; ++i)
		    if (seqs[0]->sname) fprintf(fd, "%-10s\n", (*seqs[0]->sname)[i]);
		    else	fprintf(fd, "Seq.%-6d\n", i);
		for (int i = 0; i < seqs[1]->many; ++i)
		    if (seqs[1]->sname) fprintf(fd, "%-10s\n", (*seqs[1]->sname)[i]);
		    else	fprintf(fd, "Seq.%-6d\n", seqs[0]->many + i);
		break;
	    case Form_GCG:
		gcg_out(gaps, seqs, 2, fd);
		break;
	    case Form_CLW: break;
	    default:
		if (skip < 1)	fphseqs((const Seq**) seqs, 2, fd);
		if (OutPrm.ColorEij) break;
		sm->fparam(fd);
		if (skip == -1)	return (OK);
		if (skip == -2)	return (OK);
		if (skip < 3) {
		    FTYPE	percent = fstat.mch + 
				fstat.mmc + fstat.unp; 
		    percent = 100. * fstat.mch / percent;
		    fprintf(fd, "Score = %5.1lf (%5.1f), %.1f (=), %.1f (#),",
			sscr, tscr, fstat.mch, fstat.mmc);
		    fprintf(fd, " %.1f (g), %.1f (u), (%5.2f %%)\n",
			fstat.gap, fstat.unp, percent);
		}
		if (skip == -3)	return (OK);
#if USE_WEIGHT
		if (skip < 4 && seqs[0]->weight && seqs[1]->weight) {
		    seqs[0]->fpweight();
		    seqs[1]->fpweight();
		}
#endif
		if (skip == -4)	return (OK);
		if (sigII) sigII->putSigII(fd);
		if (skip < 5) fprintf(fd, "ALIGNMENT   %d / %d\n", nbr, ttl);
		break;
	}
	praln.printaln();
	fputs("\n\n", fd);
	if (prmode == Form_GCG) {
	    setgapchar(seqs[0], DEF_GAP_CHAR);
	    OutPrm.lpw = lpw;
	    OutPrm.BlkSz = BlkSz;
	}
	return (OK);
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

void Gsinfo::printgene(Seq** seqs, int mode, FILE* _fd) const
{
	fd = _fd? _fd: out_fd;
	if (mode < 0) mode = algmode.nsa;
	int	m = seqs[0]->inex.intr? 0: 1;
	Seq*&	gene = seqs[m];
	Seq*&	qry = seqs[1 - m];
	Seq*	tsd = 0;
	INT	width = out_form->SeqBlkNo * out_form->SeqBlkSz;
	int	block = OutPrm.lpw;
	int	i = gene->right - gene->left;
	float	fscr = (float) scr / alprm.scale;
	RANGE*	rng = CDSrng;
	bool	swp = false;

	switch (mode) {		// without sequence
	    case GFF_FORM: Gff3Form(gene, qry); return;
	    case PWA_FORM: Gff3PWA(gene, qry); return;
//	    case PSL_FORM: PslForm(gene, qry); return;
	    case BED_FORM: BedForm(gene, qry); return;
	    case EXN_FORM:
	    case BIN_FORM: ExonForm(gene, qry, mode); return;
	    case ITN_FORM: IntronForm(gene, qry); return;
	    case CDS_FORM: case AAS_FORM:	// 
		tsd = gene->splice(0, rng, 0);
		tsd->tron2nuc(0);
		swp = OutPrm.deflbl == 2;
		break;
	    case CIG_FORM: CigarForm(gene, qry); return;
	    case VLG_FORM: VulgarForm(gene, qry); return;
	    case SAM_FORM: 
	    case BAM_FORM: SamForm(gene, qry); return;
	    case PSJ_FORM:	// predicted gene strucure of the reference
		rng = querygs(qry);
		swp = true;
		break;
	}
	if (swp) swap(seqs[0], seqs[1]);
	const	char*	ofn = OutPrm.out_file? strrchr(OutPrm.out_file, '/'): 0;
	if (OutPrm.out_file) ofn = ofn? ofn + 1: OutPrm.out_file;
	else if (OutPrm.deflbl != 3) ofn = gene->sqname();
	if (out_form->DbName) {
	    fprintf(fd, "%-11s %-10s %10d bp %10s %s %9.2f\n", 
		out_form->EntLabel, gene->sqname(), i, moltype[qry->inex.molc], 
		qry->sqname(), fscr);
	    block = out_form->SeqBlkSz;
	} else if (out_form->EntLabel || fst_form) {
	    int	gl = gene->left;
	    int	gr = gene->right;
	    int	ql = qry->left;
	    int	qr = qry->right;
	    if (eijnc) {
		gl = eijnc->genleft();
		gr = eijnc->genright();
		ql = eijnc->refleft();
		qr = eijnc->refright();
	    }
	    fprintf(fd, "%s", (fst_form? fst_form: out_form)->EntLabel);
	    if (ofn) {
		fprintf(fd, "%s", ofn);
		if (OutPrm.out_file) fprintf(fd, " %s", gene->sqname());
	    } else fprintf(fd, "%s.%d %s", gene->sqname(), 
		gene->SiteNo(center(1)) / 1000, gene->sqname());
	    fprintf(fd, attrfrmt, gene->Strand(), gene->len, 
		gene->SiteNo(gl), gene->SiteNo(gr - 1));
	    fprintf(fd, " %s %c", qry->sqname(), qry->Strand());
	    fprintf(fd, attrfrmt3, qry->many, qry->len,
		qry->SiteNo(ql), qry->SiteNo(qr - 1));
	    fputs(swp? (mode == PSJ_FORM? " Y": " T"): " N", fd);
	    fprintf(fd, " %9.2f\n", fscr);
	}
	if (mode != PSJ_FORM && swp) swap(seqs[0], seqs[1]);
	GBcdsForm(rng, gene, fd);
	if (mode == AAS_FORM) {		// translate
	    tsd->transout(fd, 0, 0);
	} else if (mode == PSJ_FORM) {
	    swap(seqs[0], seqs[1]);
	    qry->typeseq(fd);
	} else if (tsd) {			// without translate
	    m = tsd->left;
	    CHAR*	p = tsd->at(tsd->left);
	    if (out_form->SeqLabel) fprintf(fd, "%s\n", out_form->SeqLabel);
	    if (out_form->SeqHead)  fputs(out_form->SeqHead, fd);
	    for (i = 0; m < tsd->right; ++m, p += tsd->many) {
		if (width && i % width == 0 && out_form->SeqForm)
		    fprintf(fd, out_form->SeqForm, m + 1);
		putc(tsd->code->decode[*p], fd);
		if (!width) continue;
		if (++i % width == 0) putc('\n', fd);
		else if (i % block == 0) putc(' ', fd);
	    }
	    if (width && i % width) putc('\n', fd);
	    if (out_form->EndLabel) fprintf(fd, "%s\n", out_form->EndLabel);
	}
	if (OutPrm.spjinf && gene->exin && 
		(mode == CDS_FORM || mode == AAS_FORM))
	    BoundaryInf(gene);
	if (tsd != qry) delete tsd;
	if (rng != CDSrng) delete[] rng;
}

/*	Select one or more output modes of aln and spaln */

void AlnOutModes::getopt(const char* val)
{
	char*	arg = strrealloc(0, val);
	Strlist	sl(arg, " \t\n,:");
	for (INT i = 0; i < sl.size(); ++i) {
	    algmode.nsa = atoi(sl[i]) & 15;
	    if (out_mode[algmode.nsa] < 0) {
		out_mode[algmode.nsa] = algmode.nsa;
		++n_out_modes;
	    }
	}
	delete[] arg;
}

int AlnOutModes::setup(const char* spath)
{
	if (!n_out_modes) ++n_out_modes;
	fds = new FILE*[n_out_modes];	// directory
	vclear(fds, n_out_modes);
const	char*	sl = strrchr(spath, '/');
	if (sl) ++sl;
	else	sl = spath;
	int	slen = strlen(sl);
	bool	isdir = OutPrm.out_file && is_dir(OutPrm.out_file);
	if (isdir) {
	    slen += strlen(OutPrm.out_file);
	    prefix = new char[slen + 2];
	    strcpy(prefix, OutPrm.out_file);
	    if (prefix[slen - 1] != '/') strcat(prefix, "/");
	    strcat(prefix, sl);
	} else if (OutPrm.out_file) {	// file
	    prefix = strrealloc(0, OutPrm.out_file);
	} else {			// from query
	    prefix = strrealloc(0, sl);
	}
	char*	dt = is_file(prefix)? strrchr(prefix, '.'): 0;
	if (dt) *dt = '\0';

	char	str[LINE_MAX];
	strcpy(str, prefix);
	char*	fn = str + strlen(str);
	if (n_out_modes == 1) {
	    out_mode[0] = algmode.nsa;
	    if (OutPrm.out_file) {
		if (isdir) {
		    OutPrm.out_file = 0;
		    if (algmode.nsa == BIN_FORM) {
		        fds[0] = stdout;
			return (1);
		    }
		    sprintf(fn, ".O%d", algmode.nsa);
		} else if (algmode.nsa == BIN_FORM) {
		    fds[0] = stdout;
		    return (1);
		}
		fds[0] = wfopen(str, "w");
		if (!fds[0]) return (0);
	    } else
		fds[0] = stdout;
	    return (1);
	}
	OutPrm.out_file = 0;
	int	n = 0;
	for (int i = 0; i < N_Out_Modes; ++i) {
	    if (out_mode[i] < 0) continue;
	    if (out_mode[i] == BIN_FORM) fds[n] = 0;
	    else {
		sprintf(fn, ".O%d", i);
		fds[n] = wfopen(str, "w");
	    }
	    if (fds[n]) out_mode[n++] = out_mode[i];
	}
	return (n);
}

enum AaProp {GAP, SML, PLS, MNS, HPO, ARO = 6, NEU, PHI, PHO};

static int chempro(const int* cmp, int ii)
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

static int logonuc(const int* cmp, int ii, const SEQ_CODE* defcode)
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

static int csym(const int cmp[], int rows, const SEQ_CODE* defcode)
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

void setprmode(const int& pmd, const int& lorn, int trc)
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

static int SeqGCGCheckSum(const char* seq, int len)
{
	long	check = 0;
	
	for (int i = 0; i < len; ++i, ++seq)
	    check += ((i % 57) + 1) * toupper(*seq);
	return (check % 10000);
}

static int checksum(int* checks, const GAPS* gaps, const Seq* sd)
{
	int	len = 0;
	char	gapchar = sd->code->decode[gap_code];

const	GAPS*	gp = gaps;
	for ( ; gaps_intr(++gp); ) len += gp->gln;
	len += gp->gps - gaps->gps;
const 	CHAR*	ps = sd->at(gaps->gps);
const 	CHAR*	ts = sd->at(gp->gps);

	char*	seq = new char[len];
	for (int i = 0; i < sd->many; ++i, ++ps) {
	    int c = gaps->gps;
	    gp = gaps;
	    const CHAR*	ws = ps;
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
	    checks[i] = SeqGCGCheckSum((const char*) seq, len);
	}
	delete[] seq;
	return (len);
}

static void gcg_out(const GAPS** gaps, Seq* seqs[], int seqnum, FILE* fd)
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

	fputs("PileUp\n\n", fd);
	fprintf(fd, "   MSF:%5d  Type: ", len);
	if (seqs[0]->isprotein())	fprintf(fd,"P");
	else	fprintf(fd,"N");
	fprintf(fd,"    Check:%6ld   .. \n\n", grand_checksum);
	rows = 0;
	for (int i = 0; i < seqnum; ++i) {
	    for (int j = 0; j < seqs[i]->many; ++j, ++rows)  {
		fprintf(fd, " Name: %-24s oo  Len:%5d  Check:%6ld  Weight:",
		(*seqs[i]->sname)[j], seqs[i]->len, (long) checks[rows]);
#if USE_WEIGHT
		if (seqs[i]->weight)
		    fprintf(fd, "%8.4f\n", (double) seqs[i]->weight[j]);
		else
#endif
		    fprintf(fd, "  1.00\n");
	    }
	}
	fprintf(fd, "\n//\n");  
	delete[] checks;
}

/*	print a single (j-th) seqnuence		*/ 

void Seq::listseq(FILE* fd, int j, bool sub)
{
	if (!fd) fd = out_fd;
	CHAR*	p = at(left) + j;
	INT	width = out_form->SeqBlkNo * out_form->SeqBlkSz;
	INT	block = out_form->SeqBlkSz;
const	char*	fname = (*sname)[j];

	if (out_form->DbName)
	    fprintf(fd, "%-11s %-10s %10d bp %10s\n",
		out_form->EntLabel, fname, right - left, moltype[inex.molc]);
	else {
	    if (out_form->EntLabel) {
		if (*out_form->EntLabel) { 
		    fprintf(fd, "%s%s", out_form->EntLabel, fname);
		    if (OutPrm.fastanno) {
			if (sub) fprintf(fd, " %d/%d", j + 1, many);
			fprintf(fd, attrfrmt2, sub? 1: many, len, 
			    SiteNo(left), SiteNo(right - 1));
		    }
		    putc('\n', fd);
		    width = block = OutPrm.lpw;
		} else { 
		    fprintf(fd, "%-15s\t", fname);
		}
	    } else	width = block = OutPrm.lpw;
	}
	if (descr && OutPrm.descrp)
	    fprintf(fd, "%s%c%s\n", out_form->DefLabel, 
		out_form->DbName? '\t': ' ', (*descr)[j]);
	if (out_form->SeqLabel) fprintf(fd, "%s\n", out_form->SeqLabel);
	if (out_form->SeqHead)  fputs(out_form->SeqHead, fd);
	int	i = 0;
	for (int m = left; m < right; ++m, p += many) {
	    if (OutPrm.trimendgap && IsGap(*p)) continue;
	    if (width && i % width == 0 && out_form->SeqForm)
		fprintf(fd, out_form->SeqForm, m + 1);
	    putc(code->decode[*p], fd);
	    if (!width) continue;
	    if (++i % width == 0) putc('\n', fd);
	    else if (i % block == 0) putc(' ', fd);
	}
	if (OutPrm.asterisk) putc('*', fd);
	if (!width || i % width) putc('\n', fd);
	if (out_form->EndLabel) fprintf(fd, "%s\n", out_form->EndLabel);
}

int Seq::calcnbr(int gp, int i)
{
	int 	n = nbr[i];
	CHAR*	s = at(0) + i;

	for ( ; gp-- > 0; s += many) if (IsntGap(*s)) ++n;
	return (n);
}

void Seq::fphseq(FILE* fd, int n) const
{
	if (!fd) fd = out_fd;
	char	str[MAXL];
static	const char* frm[] = {"%c%s [%d:%d] ", "%s ( %d - %d )",  "%-16s ( %3d - %3d )"};

	if (n <= 0) {
	    fputs(sqname(), fd);
	} else if (n == 1) {
	    fprintf(fd, frm[0], senschar(), sqname(), many, len);
	} else if (n == 2) {
	    fprintf(fd, frm[1], sqname(), SiteNo(left), SiteNo(right - 1));
	} else {
	    n = n == 3? 1: 2;
	    sprintf(str, frm[0], senschar(), sqname(), many, len);
	    fprintf(fd, frm[n], str, SiteNo(left), SiteNo(right - 1));
	}
}

void fphseqs(const Seq* seqs[], int n, FILE* fd)
{
	if (!fd) fd = out_fd;
	putc('\n', fd);
	while (n--) {
	    (*seqs++)->fphseq(fd);
	    if (n) fputs(" - ", fd);
	}
	putc('\n', fd);
}

#if USE_WEIGHT
void Seq::fpweight(FILE* fd) const
{
	if (!weight) return;
	if (!fd) fd = out_fd;
	int 	i = 0;
	while (i < many) {
	    if (i % NO_WGHT == 0) putc(_WGHT, fd);
	    fprintf(fd, " %14.7e", weight[i]);
	    if (++i % NO_WGHT == 0) putc('\n', fd);
	}
	if (i % NO_WGHT) putc('\n', fd);
}
#endif

void fprint_seq_mem(const Seq* seqs[], const int& n, FILE* fd)
{
	if (!fd) fd = out_fd;
	for (int i = 0; i < n; ++i) {
	    seqs[i]->fphseq(fd, algmode.rng? 3: 0);
	    putc(' ', fd);
	}
	putc('\n', fd);
}

/*	Map the sequence number to column position
	Array must be sorted and terminated with a minus value
	Both number and position start with 1	*/

void Seq::num2pos(int which, int* array) const
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

void Seq::pos2num(int which, int* array) const
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

int Seq::calcResNum(const int& i)
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

PrintMember::PrintMember(const Strlist* sn, bool pad_space, const char* tl) 
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

void PrintMember::put_member(FILE* fd, int i) const
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

PrintAln::PrintAln(const GAPS**  _gaps, Seq* _seqs[], const int& _seqnum, FILE* _fd)
	: gaps(_gaps), seqs(_seqs), seqnum(_seqnum), 
	  globpt(1), cpm(prmode), markeij(0), fd(_fd? _fd: out_fd)
{
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
	    ecc = markeij == 1? new EscCharCtl(fd): 0;
	    hcc = markeij == 2? new HtmlCharCtl(fd, str): 0;
	    pfqs = new PFQ*[seqnum];
	    tfqs = new PFQ*[seqnum];
	    lsts = new int*[seqnum];
	}
	for (int i = rows = 0; i < seqnum; ++i) rows += seqs[i]->many;
	gp  = new const GAPS*[seqnum];
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
		putc(' ', fd);
	    int	ch;
	    if (res == BLANK || (brc && *brc == BLANK))
		ch = ' ';
	    else if (*cur & INTRONBIT)
		ch = tolower(decode[res]);
	    else if (qry)
		ch = *qry == res? ditto: tolower(decode[res]);
	    else
		ch = decode[res];
	    if (ccd) {
		if (ecc) ecc->putchr(ch, CMC_WHITE, iis_color[ccd], CMC_BOLD);
		if (hcc) hcc->putchr(ch, "white", hml_color[ccd]);
	    } else	putc(ch, fd);
	    if (qry) ++qry;
	    if (brc) ++brc;
	}
}

void PrintAln::putclmark()
{
	int	i;

	if (clmark == 1) {
	    for (i = globpt; i % 10; i++) putc(' ', fd);
	    for ( ; i < globpt + OutPrm.lpw; i += 10)
		fputs("       .  ", fd);
	} else {
	    fprintf(fd, ";%03d", globpt / OutPrm.lpw + 1);
	    for (i = globpt; i % 10; i++) putc(' ', fd);
	    fprintf(fd, "%6d", i);
	    for (i += 10; i < globpt + OutPrm.lpw; i += 10)
		fprintf(fd, "%10d", i);
	}
	putc('\n', fd);
}

void PrintAln::pair_mrk(const CHAR* p1, const CHAR* p2, int clms)
{
	Simmtx*	sm = getSimmtx(0);
	fputs(spc9, fd);
	for (int i = 0; i < clms; ++i) {
	    int	cc = (*p1 == BLANK || *p2 == BLANK)? 0:
	    	sm->simgrade(*p1++ & RES_BITS, *p2++ & RES_BITS);
	    if (emphsim < 0) cc = nmk - cc;
	    if (cc > nmk) cc = nmk;
	    if (cc < 0) cc = 0;
	    putc(mark[cc], fd);
	}
	putc('\n', fd);
}

void PrintAln::calc_mrk(int rows, int clms, const SEQ_CODE* defcode)
{
	int*	cmp = new int[defcode->max_code];

	fputs(spc9, fd);
	for (int l = 0; l < clms; ++l) {
	    vclear(cmp, defcode->max_code);
	    for (int i = 0; i < rows; ++i) {
		int	cc = image[i][l] & RES_BITS;
		if (cc == BLANK) ++cmp[nil_code];
		else		++cmp[cc];
	    }
	    putc(csym(cmp, rows, defcode), fd);
	}
	putc('\n', fd);
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

	putc('\n', fd);
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
	      while (n--) putc(' ', fd);
	      seqline(prv, brc, *img, OutPrm.lpw, decode);
	      putc('\n', fd);
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
		    if (globpt == 1) fprintf(fd, "%-12s", (*sd->sname)[i]);
		    else	fputs(spc12, fd);
		    break;
		case LABEL_Name:
		    fprintf(fd, "%-24s", (*sd->sname)[i]); break;
		case LABEL_Numb:
		    fprintf(fd, "%8d ", sd->SiteNo(*lft, i) % 100000000);
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
		case LABEL_None: putc('\n', fd); break;
		case LABEL_Numb:
		    fprintf(fd, "%c%5d\n", _LABL, sd->SiteNo(*right - 1, i) % 100000000); 
		    ++right; break;
		case LABEL_Name:
	    	default:
		    fprintf(fd, "%c ", _LABL);
		    prm.put_member(fd, i); break;
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
	    putc('\n', fd);
}

void PrintAln::markiis(int k, int j, int clm)
{
	int phs = pfqs[j]->pos;
	if (htl == 0) --phs;
	phs = (phs % 3 + 1) << 5;
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
const	RANGE*	exon = 0;
	int	maxleft = 0;
	gene = -1;
	for (int j = htl = 0; j < seqnum; ++j) {
	    Seq*&	sd = seqs[j];
	    if (sd->left > maxleft) maxleft = sd->left;
	    if (sd->isprotein()) htl |= 1;
	    if (sd->inex.intr && sd->exons) {
		htl |= 2;		// 0: CvsC
		exon = sd->exons + 1;	// 1: AvsA
		gene = j;		// 2: GvsC
		gpos = gaps[j]->gps;	// 3: GvsA
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
	pro = (gene >= 0)? gene + 1: -1;
const	int	c_step = (htl == 1)? 3: 1;
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
	if (htl & 2) prmode = Row_None;
	do {
	    if (OutPrm.SkipLongGap) {
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
const		int	pos = gp[jj]->gps - gaps[jj]->gps; 	// print one line
const		int	upr = (gap - z - OutPrm.EijMergin) / OutPrm.lpw * OutPrm.lpw;
		if (gap < INT_MAX && z - pos > OutPrm.EijMergin && upr > 0) {
		    for (int j = k = 0; j < seqnum; k += seqs[j++]->many) {
			Seq*&	sd = seqs[j];
const			bool	step3 = htl == 3 && !sd->inex.intr;
			if (j == jj) {
			    if (step3) gph = phs = (phs + upr) % 3;
			} else {
			    int	inc = upr;
			    if (step3) {
				int	r = upr % 3;
				inc /= 3;
				if (r && phs != 2 && (r + phs) > 1) ++inc;
			    }
			    wkr[j] += inc * sd->many;
			    for (int i = 0; i < sd->many; ++i)
				nbr[k + i] += inc;
			}
		    }
		    z += upr;
		    globpt += upr;
		    gpos += upr;
		    fprintf(fd, "\n;; skip %d nt\'s\n", upr);
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
		    if (exon->right == gpos + 1 && phs == 1) reij = 3;
		}
		if (active) active = seqnum;
		for (int j = k = 0; j < seqnum; k += seqs[j++]->many) {
		    Seq*&	sd = seqs[j];
		    bool	prot = sd->isprotein() && htl == 3;
	loop:	    int		pos = gp[j]->gps - gaps[j]->gps;
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
			} else {
			    for (int i = 0; i < sd->many; ++i) {
				image[k+i][clm] = *wkr[j] == nil_code? gap_code: *wkr[j];
				if (j == gene) {
				    int	df = gpos - exon->left;
				    if (df < 0) image[k+i][clm] |= INTRONBIT;
				    else if (intlen > 2 && df < 2 && (htl & 2) && pphs == 1)
					reij = 2 - df;
				    if (markeij && reij) 
					image[k+i][clm] |= (reij << 5);
				}
				if (IsntGap(*wkr[j])) ++nbr[k+i];
				++wkr[j];
			    }
			}
			if ((htl & 2) && j == pro) gph = phs = (phs + 1) % 3;
			if (markeij && pfqs[j]) {
			    int	niis = 0;
			    while ((pfqs[j] + niis < tfqs[j]) &&
				pfqs[j][niis].pos + agap[j] + ((htl & 2)?
				    1 - pfqs[j][niis].pos % 3: 0) < cpos)
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
			} else {
			    for (int i = 0; i < sd->many; i++)
				image[k+i][clm] = (!exon || gpos >= exon->left)?
				    gap_code: BLANK;
			}
			if (htl & 2) gph = (gph + 1) % htl;
			if (!neog && !--active && !clm--) return;
		    } else {
			agap[j] += gp[j]->gln * c_step;
			gp[j]++;
			goto loop;
		    }
		}
	    }
	    if (fd)  prnt_aln(left, nbr);
	    globpt += OutPrm.lpw;
	} while(active);
}

void pralnseq(const GAPS** gaps, Seq* seqs[], const int& seqnum, FILE* fd)
{
	if (!OutPrm.lpw) return;
	PrintAln	praln(gaps, seqs, seqnum, fd);
	praln.printaln();
}

static void put_SigII(const PFQ* pfq, const PFQ* tfq, const int* lst, 
	int numlst, FILE* fd)
{
	int	lwd = OutPrm.lpw >= 10? OutPrm.lpw - 4: 56;
	char	str[MAXL];

	fprintf(fd, ";B %ld %d\n", (long) (tfq - pfq), numlst);
	int	c = 0;
	for ( ; pfq < tfq - 1; ++pfq) {
	    if (c == 0) fputs(";b", fd);
	    sprintf(str + c, " %d %d,", pfq->pos, pfq->num);
	    if ((c = strlen(str)) > lwd) {
		fputs(str, fd);
		fputs("\n", fd);
		c = 0;
	    }
	}
	if (c)	fputs(str, fd);
	else	fputs(";b", fd);
	fprintf(fd, " %d %d\n", pfq->pos, pfq->num);
	if (!lst || !numlst) return;
	const int*	tst = lst + numlst;
	for (c = 0; lst < tst - 1; ++lst) {
	    if (c == 0) fputs(";m", fd);
	    sprintf(str + c, " %d", *lst + 1);
	    if ((c = strlen(str)) > lwd) {
		fputs(str, fd);
		fputs("\n", fd);
		c = 0;
	    }
	}
	if (c)	fputs(str, fd);
	else	fputs(";m", fd);
	if (lst < tst) fprintf(fd, " %d\n", *lst + 1);
}

void Seq::putSigII(FILE* fd) const
{
const	PFQ*	pfq = sigII? sigII->pfq: 0;
const	int*	lst = sigII? sigII->lst: 0;

	if (!pfq) {
	    fputs(";B 0 0\n", fd);
	    return;
	}
	int	step = isprotein()? 3: 1;
	int	upto = step * (left + base_);
	for ( ; pfq->pos < upto; ++pfq)
	    if (lst) lst += pfq->num;
	upto = step * (right + base_);
const	PFQ*	wfq = pfq;
	int	n = 0;
	for ( ; wfq->pos < upto; ++wfq)
	    if (lst) n += wfq->num;
	put_SigII(pfq, wfq, lst, n, fd);
}

void SigII::putSigII(FILE* fd) const
{
	if (!step || !pfqnum) return;
	put_SigII(pfq, pfq + pfqnum, lst, lstnum, fd);
}

static void nexus_head(FILE* fd, const Seq* sd)
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
	if (in_line) {
	    CHAR*	t = at(right);
	    for (CHAR*	p = at(left); p < t; ++p) 
		putc(code->decode[*p], fd);
	} else if (many == 1) {
	    listseq(fd);
	} else if (prmode >= Form_Serial) {
	    if (OutPrm.noseqline)
		fprintf(fd, "%5d %5d\t%s\n", many, len, sqname());
	    else if (out_form->FormID == NEXUS)
		nexus_head(fd, this);
	    for (int j = 0; j < many; ++j) listseq(fd, j, true);
	    if (out_form->FormID == NEXUS) fputs(";\nend;\n", fd);
	} else {
const	    GAPS	gap[2] = {{left, 0}, {right, gaps_end}};
const	    GAPS*	gpp[2] = {gap, gap + 1};
	    Seq*	tmp[1] = {this};
	    PrintAln	praln(gpp, tmp, 1, fd);
	    if (prmode <= Form_Native) {
		if (OutPrm.ColorEij) {
		    printf(">%s\n", sqname());
		} else {
		    fphseqs((const Seq**) tmp, 1, fd);
#if USE_WEIGHT
		    if (weight && OutPrm.printweight) fpweight(fd);
#endif
		    if (sigII) putSigII(fd);
		}
	    } else if (prmode == Form_Phylp) {
		fprintf(fd, "%5d %5d\n", many, len);
	    } else if (prmode == Form_GCG) {
		gcg_out(gpp, tmp, 1, fd);
	    }
	    if (OutPrm.lpw) praln.printaln();
	}
}
