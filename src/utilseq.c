/*****************************************************************************
*
*	Micelleious utilities subroutines
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

enum Triplet {AAA = 0, AGA = 8, AGG = 10, AUA = 12, 
	CUA = 28, CUC, CUG, CUU,
	UAA = 48, UAG = 50, UCA = 52, UGA = 56, UUA = 60};

static	const	char	DEF_CDI_PATH[] = "/data2/seqdb/blast/cdi";
static	const	char	DEF_CDI[] = "BLASTCDI";
static	const	int	PassCom = 0x08;
static	const	int	MAXICD = 4;
static	const	int	MAXTCD = 16;
static	const	float	ambLimit = 0.2;

inline	bool	isNucAmb(int c) {return (c != A && c != C && c != G && c != T); }

CHAR gencode[64] = {
	LYS,ASN,LYS,ASN,THR,THR,THR,THR,ARG,SER,ARG,SER,ILE,ILE,MET,ILE,
	GLN,HIS,GLN,HIS,PRO,PRO,PRO,PRO,ARG,ARG,ARG,ARG,LEU,LEU,LEU,LEU,
	GLU,ASP,GLU,ASP,ALA,ALA,ALA,ALA,GLY,GLY,GLY,GLY,VAL,VAL,VAL,VAL,
	TRM,TYR,TRM,TYR,SER,SER,SER,SER,TRM,CYS,TRP,CYS,LEU,PHE,LEU,PHE
};

float CodonUsage[64] = {
	0.4164,0.5527,0.5836,0.4473,0.2709,0.4101,0.0810,0.2380,
	0.1941,0.2329,0.2337,0.1440,0.1326,0.5458,1.0000,0.3216,
	0.2625,0.5648,0.7375,0.4352,0.2782,0.3364,0.0793,0.3060,
	0.1090,0.1744,0.1655,0.1233,0.0760,0.2177,0.3970,0.1192,
	0.4375,0.4994,0.5625,0.5006,0.2669,0.3758,0.0587,0.2986,
	0.2884,0.3134,0.2319,0.1663,0.1093,0.2715,0.4310,0.1883,
	0.3780,0.5638,0.1829,0.4362,0.1475,0.2427,0.0386,0.1943,
	0.4390,0.4695,1.0000,0.5305,0.0569,0.5343,0.1332,0.4657
};

static float Human50CDU[64] = {
	0.4164,0.5527,0.5836,0.4473,0.2709,0.4101,0.0810,0.2380,
	0.1941,0.2329,0.2337,0.1440,0.1326,0.5458,1.0000,0.3216,
	0.2625,0.5648,0.7375,0.4352,0.2782,0.3364,0.0793,0.3060,
	0.1090,0.1744,0.1655,0.1233,0.0760,0.2177,0.3970,0.1192,
	0.4375,0.4994,0.5625,0.5006,0.2669,0.3758,0.0587,0.2986,
	0.2884,0.3134,0.2319,0.1663,0.1093,0.2715,0.4310,0.1883,
	0.3780,0.5638,0.1829,0.4362,0.1475,0.2427,0.0386,0.1943,
	0.4390,0.4695,1.0000,0.5305,0.0569,0.5343,0.1332,0.4657
};

static SeqDb TransGB = {
	GenBank,
	PROTEIN,
	"TransGB",
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	"        	     /translation=\"",
	"	             ",
	1, 50, 0, 0
};

static SeqDb TransEB = {
	EMBL,
	PROTEIN,
	"TransEB",
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	NULL,
	"FT                   /translation=\"",
	"FT                   ",
	1, 50, 0, 0
};

static	int	fromic = 1;
static	int	minorf = 300;

static	CHAR*	de_codon(CHAR* ncs, int n);
static	int	lcomp(ORF* a, ORF* b);
static	int	evenusage();
static	int	kmer(const CHAR* s, int m, const CHAR* elements, int k, int* x = 0);

static	int	curcode = 0;
static	int	ChangeCodon = 1;
static	CHAR	initcode[4 * MAXICD];
static	CHAR	termcode[4 * MAXTCD];

int setcodon(int c)
{
	while ((c = initcodon(c)) == ERROR) {
	    c = curcode;
	    promptin("EukUniv[0]/NCBI_Geneci_Code (%d) : ", &c);
	};
	return (curcode = c);
}

void setorf(int len, int ic)
{
	if (len == QUERY || ic == QUERY)
	    promptin("min-orf (%d), from-init-codon (%d) : ", &minorf, &fromic);
	else {
	    if (len >= 0) minorf = len;
	    if (ic >= 0) fromic = ic;
	}
}

int codon_id(const CHAR* s, int many)
{
	int	c = 0;

	for (int i = 0; i < 3; ++i, s += many) {
	    int	n = ncredctab[*s];
	    if (n > 3) return (ERROR);
	    c = 4 * c + n;
	}
	return (c);
}

static CHAR* de_codon(CHAR* ncs, int n)
{
	int	i = 3;

	ncs[3] = '+';
	while (i-- > 0) {
	    ncs[i] = red2nuc[n % 4];
	    n /= 4;
	}
	return (ncs + 4);
}

void de_codon_4(CHAR* ncs, int n)
{
	int	i = 3;

	while (i-- > 0) {
	    ncs[i] = n % 4;
	    n /= 4;
	}
}

static	const CHAR most_abund[4] = {LYS, ALA, GLY, LEU};

int toaa(const CHAR* ns)
{
	int	c1;
	int	c2 = ns[1];

	if (IsGap(c2))	return (UNP);		/* 2nd codon is del */
	c2 = ncredctab[c2];
	if (c2 >= 4) return (AMB);		/* 2nd codon is amb */
	c1 = ncredctab[ns[0]];
	if (c1 >= 4) return (most_abund[c2]);	/* 1st codon is amb */
	return (gencode[16 * c1 + 4 * c2 + ncelements[ns[2]]]);
}

int toaa3(const CHAR* ns, int inc)
{
	int	c1;
	int	c2 = ns[inc];

	if (IsGap(c2))	return (UNP);		/* 2nd codon is del */
	c2 = ncredctab[c2];
	if (c2 >= 4) return (AMB);		/* 2nd codon is amb */
	c1 = ncredctab[ns[0]];
	if (c1 >= 4) return (most_abund[c2]);	/* 1st codon is amb */
	return (gencode[16 * c1 + 4 * c2 + ncelements[ns[inc+inc]]]);
}

int nuc2tron3(const CHAR* ns, int inc)
{
	int	c2 = ns[inc];

	if (IsGap(c2))	return (UNP);		/* 2nd codon is del */
	c2 = ncredctab[c2];
	if (c2 >= 4) return (AMB);		/* 2nd codon is amb */
	int	c1 = ncredctab[ns[0]];
	int	aa = (c1 >= 4)? most_abund[c2]: /* 1st codon is amb */
		gencode[16 * c1 + 4 * c2 + ncelements[ns[inc+inc]]];

	switch (aa) {
	    case SER:
		if (ns[inc] == G) aa = SER2;
		break;
	    case TRM:
		if (ns[inc] == G) aa = TRM2;
		break;
	}
	return (aa);
}

int initcodon(int code)
{
	switch(code) {
	  case 0: 		// Reset to Universal codon
		gencode[AAA] = LYS; gencode[AGA] = gencode[AGG] = ARG; gencode[AUA] = ILE;
		gencode[CUA] = gencode[CUC] = gencode[CUG] = gencode[CUU] = LEU;
		gencode[UAA] = gencode[UAG] = gencode[UGA] = TRM;
		break;
	  case 1: case 11:	// Universal
		break;
	  case 2:		// Vertebrate Mitochondrial
		gencode[AGA] = gencode[AGG] = TRM; gencode[AUA] = MET;
		gencode[UGA] = TRP;
		break;
	  case 3:		// Yeast Mitochondrial,**incompatible with TRON code**
		gencode[AUA] = MET;
		gencode[CUA] = gencode[CUC] = gencode[CUG] = gencode[CUU] = THR;
		gencode[UGA] = TRP;
		break;
	  case 4:		// Mold, Protozoan, and Coelenterate Mitochondrial
		gencode[UGA] = TRP;
		break;
	  case 5:		// Invertebrate Mitochondrial
		gencode[AGA] = gencode[AGG] = SER; gencode[AUA] = MET;
		gencode[UGA] = TRP;
		break;
	  case 6:		// Ciliate, Dasycladacean and Hexamita Nuclear 
		gencode[UAA] = gencode[UAG] = GLN;
		break;
	  case 9:		// Echinoderm and Flatworm Mitochondrial
		gencode[AAA] = ASN; gencode[AGA] = gencode[AGG] = SER;
		gencode[UGA] = TRP;
		break;
	  case 10:	// Euplotid Nuclear
		gencode[UGA] = CYS;
		break;
	  case 12:	// Alternative Yeast Nuclear (Candida)
			// C. albicans, cylindracea, melibiosica, parapsilosis, rugos
		gencode[CUG] = SER;	//**incompatible with TRON code**
		break;
	  case 13:		// Ascidian Mitochondrial
		gencode[AGA] = gencode[AGG] = GLY; gencode[AUA] = MET;
		gencode[UGA] = TRP;
		break;
	  case 14:		// Alternative Flatworm Mitochondrial
		gencode[AAA] = ASN; gencode[AGA] = gencode[AGG] = SER;
		gencode[UAA] = TYR; gencode[UGA] = TRP;
		break;
	  case 16:	// Chlorophycean Mitochondrial
		gencode[UAG] = LEU;	//**incompatible with TRON code**
		break;
	  case 21:		// Trematode Mitochondrial
		gencode[AAA] = ASN; gencode[AGA] = gencode[AGG] = SER; gencode[AUA] = MET;
		gencode[UGA] = TRP;
		break;
	  case 22:	// Scenedesmus obliquus Mitochondrial
		gencode[UAG] = LEU; gencode[UCA] = TRM;	//**incompatible with TRON code**
		break;
	  case 23:	// Thraustochytrium Mitochondrial
		gencode[UUA] = TRM;	//**incompatible with TRON code**
		break;
	  case 24:		// Pterobranchia Mitochondrial
		gencode[AGA] = SER; gencode[AGG] = LYS; //**incompatible with TRON code**
		gencode[UGA] = TRP;
		break;
	  case 25:	// Division SR1 and Gracilibacteria
		gencode[UGA] = GLY;
		break;
	  case 26:	// Pachysolen tannophilus Nuclear
		gencode[CUG] = ALA;	//**incompatible with TRON code*
		break;
	  case 27:	// Karyorelict Nuclear
		gencode[UAG] = gencode[UAA] = GLN;
		break;
	  case 28:	// Condylostoma Nuclear
		gencode[UAA] = gencode[UAG] = GLN;
		break;
	  case 29:	// Mesodinium Nuclear
		gencode[UAA] = gencode[UAG] = TYR;
		break;
	  case 30:	// Peritrich Nuclear
		gencode[UAA] = gencode[UAG] = GLU;
		break;
	  case 31:	// Blastocrithidia Nuclear
		gencode[UGA] = TRP;
		break;
	  default: 
		return (ERROR);
	}
	int	i, j, k;
	CHAR*	tc = termcode;
	CHAR*	ic = initcode;

	for (i = j = k = 0; i < 64; i++) {
	    if (gencode[i] == MET && j++ < MAXICD)
		ic = de_codon(ic, i);
	    if (gencode[i] == TRM && k++ < MAXTCD)
		tc = de_codon(tc, i);
	}
	ic[-1] = tc[-1] = 0;
	ChangeCodon = code != curcode;
	return (code);
}

int initcodon(const char* genspc)
{
	int	code = 0;
	if (genspc)  {
	    if (!strncmp(genspc, "tetrther", 8)) code = 6; else
	    if (!strncmp(genspc, "paratetr", 8)) code = 6; else
	    if (!strncmp(genspc, "candalbi", 8)) code = 12; else
	    if (!strncmp(genspc, "candcyli", 8)) code = 12; else
	    if (!strncmp(genspc, "candmeli", 8)) code = 12; else
	    if (!strncmp(genspc, "candpara", 8)) code = 12; else
	    if (!strncmp(genspc, "candrugo", 8)) code = 12;
	}
	return (initcodon(code));
}

static int lcomp(ORF* a, ORF* b)
{
    return (b->len - a->len);
}

/* DFA to find ORF 

fromic == 0: Term< --- >Term
fromic == 1: <ATG  --- >Term
fromic == 2: ATG - AG< --- >GT - Term
fromic == 3: <ATG - AG --- >GT - Term
fromic == 4: ATG - AG< --- GT - >Term
fromic == 5: <ATG - AG --- GT - >Term

*/

ORF* Seq::getorf() const
{
static	const	int	tm[14][4] = {
/*       A  C  G  T */
	{1, 0, 2, 3},	/*  0: C */
	{1, 0, 4, 5},	/*  1: A */
	{1, 6, 2, 7},	/*  2: G */
	{8, 0, 9, 3},	/*  3: T */
	{1, 6, 2, 7},	/*  4: AG */
	{8, 0,10, 3},	/*  5: AT */ 
	{1, 0, 2, 3},	/*  6: GC */
	{8, 0, 9, 3},	/*  7: GT */
	{11,0,12, 5},	/*  8: TA */
	{13,6, 2, 7},	/*  9: TG */
	{13,6, 2, 7},	/* 10: ATG */
	{1, 0, 4, 5},	/* 11: TAA */
	{1, 6, 2, 7},	/* 12: TAG */
	{1, 0, 4, 5}};	/* 13: TGA */

	ORF	cur[3];
	ORF*	orfs;
	int	f;			// frame
	int	n = left;
	for (f = 0; f < 3; ++f) {
	    cur[f].pos = n;
	    cur[f].frm = fromic? 0: 1;	/* close/open */
	    cur[f].len = 0;
	}
const	CHAR*	ts = at(right);
	Mfile*	mfd = new Mfile(sizeof(ORF));
const	CHAR*	redctab = isdrna()? ncredctab: tnredctab;
	int	state = 0;
	for (const CHAR* ps = at(left); ps < ts; ++n) {
	    int	c = redctab[*ps++];
	    if (c < 4) {
		state = tm[state][c];
		switch (state) {
		  case 4:	/* AG */
		    if (fromic > 1) {	/* exon start */
			for (f = 0; f < 3; ++f) {
			    if (!cur[f].frm || (fromic % 2 == 0 && cur[f].frm == 1)) {
				cur[f].pos = n + 1;
				cur[f].frm = 2;
			    }
			}
		    }
		    break;
		  case 6:	/* GC */
		  case 7:	/* GT */
		    if (fromic > 1) {	/* exon end */
			for (f = 0; f < 3; ++f) {
			    if (cur[f].frm) {
				cur[f].len = n - 1 - cur[f].pos;
				if (fromic == 2 || fromic == 3)
				    cur[f].frm = 4;	/* mark */
			    }
			}
		    }
		    break;
		  case 10:	/* ATG */
		    f = (n - 2) % 3;
		    if (!cur[f].frm) {
			cur[f].pos = n - 2;
			cur[f].frm = 1;
		    }
		    break;
		  case 11:	/* TAA */
		    if (gencode[UAA] != TRM) break;
		  case 12:	/* TAG */
		    if (gencode[UAG] != TRM) break;
		  case 13:	/* TGA */
		    if (gencode[UGA] != TRM) break;
		    f = (n - 2) % 3;
		    if (cur[f].frm == 1 || cur[f].frm == 2)
			cur[f].len = n - 2 - cur[f].pos;
		    if (cur[f].len >= minorf) { 
			cur[f].frm = f;
			mfd->write(cur + f);
		    }
		    cur[f].pos = n + 1;
		    cur[f].frm = fromic? 0: 1;
		    cur[f].len = 0;
		    if (state == 12 && fromic > 1) {	/* AG */
			for (f = 0; f < 3; ++f) {
			    if (!cur[f].frm || (fromic % 2 == 0 && cur[f].frm == 1)) {
				cur[f].pos = n + 1;
				cur[f].frm = 2;
			    }
			}
		    }
		    break;
		  default: break;
		}
	    } else {		/* ambiguous base */
		for (f = 0; f < 3; ++f) {
		    cur[f].pos = n + 1;
		    cur[f].frm = fromic? 0: 1;
		    cur[f].len = 0;
		}
	    }
	}
	for (f = 0; f < 3; ++f) {
	    if (cur[f].frm) cur[f].len = n - cur[f].pos;
	    if (cur[f].len >= minorf) {
		cur[f].frm = f;
		mfd->write(cur + f);
	    }
	}
	cur->pos = INT_MAX;
	cur->len = cur->frm = 0;
	mfd->write(cur);
	n = mfd->size();
	orfs = (ORF*) mfd->flush();
	delete mfd;
	if (n < 2) {
	    delete[] orfs; orfs = 0;
	} else {
	    qsort(orfs, n, sizeof(ORF), (CMPF) lcomp);
	}
	return (orfs);
}

void Seq::passcom(FILE* fo) const
{
	FILE*	fi = fopen(spath, "r");
	if (!fi) return;

	int	c;
	while ((c = fgetc(fi)) != EOF) {
	    if (c != _COMM) continue;
	    for (;;) {
		fputc(c, fo);
		c = fgetc(fi);
		if (c == EOF || c == '\n') break;
	    }
	}
	fclose(fi);
}

Seq* Seq::translate(Seq* aas, ORF& orf) const
{
const	CHAR*	ns = at(orf.pos + 3 - (orf.pos - orf.frm) % 3);
const	CHAR*	ts = ns + orf.len;
	int	ln = orf.len / 3 + 2;

	if (aas)	aas->refresh(1, ln);
	else		aas = new Seq(1, ln);
	setSeqCode(aas, PROTEIN);
	CHAR*	as = aas->at(0);
	if (isdrna()) {
	    for (; ns < ts; ns += 3, ++as) {
		*as = toaa(ns);
		if (aas->isAmb(*as)) aas->inex.ambs = 1;
		if (aas->isGap(*as)) aas->inex.dels = 1;
	    }
	} else if (istron()) {
	    for (++ns; ns < ts; ns += 3, ++as) {
		*as = *ns == SER2? SER: *ns;
		if (aas->isAmb(*as)) aas->inex.ambs = 1;
		if (aas->isGap(*as)) aas->inex.dels = 1;
	    }
	} else {
	    prompt("Can't translate %s !\n", sqname());
	}
	return aas->postseq(as);
}

void Seq::ftranslate(FILE* fd, int at_term, SeqDb* form, int orfn) const
{
	int 	i = 0;
const	CHAR*	seq = at(left);
	int 	m = left;
	int 	width = form->SeqBlkNo * form->SeqBlkSz;
	int 	block = form->SeqBlkSz;

	switch (form->FormID) {
	  case	GenBank: form = &TransGB; break;
	  case	EMBL:	form = &TransEB; break;
	  case	Bare:	break;
	  default:
	    if (form->DbName)
		fprintf(fd, "%-9s %s", form->EntLabel, sqname());
	    else {
		if (form->EntLabel)
		    fprintf(fd, "%s%s",  form->EntLabel, sqname());
		width = block = getlpw();
	    }
	    if (inex.sens || orfn) putc('_', fd);
	    if (inex.sens) putc('C', fd);
	    if (orfn) fprintf(fd, "%d", orfn);
	    putc('\n', fd);
	    if (algmode.nsa) {
		if (inex.sens)
		    fprintf(fd, ";c CDS_C %d %d\n", SiteNo(left),
		    SiteNo(right - 1));
		else
		    fprintf(fd, ";c CDS %d %d\n", SiteNo(left),
			SiteNo(right - 1));
	    }
	}
	if (form->FormID == GenBank || form->FormID == EMBL) {
	    width = form->SeqBlkNo * form->SeqBlkSz;
	    block = form->SeqBlkSz;
	    i = width - strlen(form->SeqHead);
	}
	if (algmode.nsa & PassCom) passcom(fd);
	if (form->SeqLabel) fprintf(fd, "%s\n", form->SeqLabel);
	if (form->SeqHead)  fputs(form->SeqHead, fd);
	for ( ; m < right; seq += 3) {
	    if (i % width == 0 && form->SeqForm)
		fprintf(fd, form->SeqForm, i + 1);
	    int	aa = toaa(seq);
	    m += 3;
	    if (aa == TRM && (at_term == STOP ||
		(at_term == 0 && m >= right))) break;
	    putc(amino[aa], fd);
	    if (++i % width == 0) putc('\n', fd);
	    else if (i % block == 0) putc(' ', fd);
	}
	if (form == &TransGB || form == &TransEB) {
	    if (i % width == 0 && form->SeqForm) {
		fputs(form->SeqForm, fd);
		++i;
	    }
	    putc('\"', fd);
	}
	if (i % width) putc('\n', fd);
	if (form->EndLabel) fprintf(fd, "%s\n", form->EndLabel);
}

void Seq::fmtranslate(FILE* fd, int at_term, int orfn) const
{
	int*	label = new int[many];

	fprintf(fd, "%s%s", SeqDBs[FASTA].EntLabel, sqname());
	if (orfn)   fprintf(fd, "_%d [%d]\n\n", orfn, many);
	else	    fprintf(fd, " [%d]\n\n", many);
	if (algmode.nsa & PassCom) passcom(fd);
	for (int i = 0; i < many; i++) label[i] = 1;
	for (int m = left; m < right; m += 3 * OutPrm.lpw) {
	    const CHAR*	sq = at(m);
	    for (int i = 0; i < many; i++, sq++) {
		fprintf(fd, "%6d ", label[i]);
		int	c = (right - m) / 3;
		if (c > OutPrm.lpw) c = OutPrm.lpw;
		const CHAR*	b = sq;
		int j = 0;
		for ( ; j++ < c; b += 3 * many) {
		    int	aa = toaa3(b, many);
		    if (aa == TRM && (at_term == STOP ||
			(at_term == 0 && j == c && c < OutPrm.lpw))) break;
		    putc(amino[aa], fd);
		    if (aa != UNP) label[i]++;
		}
		while (j++ < OutPrm.lpw) putc(' ', fd);
		fprintf(fd, "| %s\n", (*sname)[i]);
	    }
	    putc('\n', fd);
	}
	delete[] label;
}

int Seq::transout(FILE* fd, int at_term, int orfn) const
{
	if (!isdrna()) {
	    prompt("Not a nucleotide sequence !\n");
	    return (ERROR);
	}
	if (many == 1) 
	    ftranslate(fd, at_term, setform(0), orfn);
	else
	    fmtranslate(fd, at_term, orfn);
	return (OK);
}

int*	invtranslate[ZZZ];
static	int	ndeg[ZZZ];

void mkinvtab()
{
static	int	invtab[64 + ZZZ];

	if (!ChangeCodon && invtranslate[0]) return;
	for (int a = 0; a < ZZZ; ++a) ndeg[a] = 0;
	for (int i = 0; i < 64; ++i) ndeg[gencode[i]]++;
	int	n = 0;
	for (int a = 0; a <= AAS; ) {
	    invtranslate[a] = invtab + n;
	    n += ndeg[a];
	    invtab[n++] = -1;
	    ndeg[a++] = 0;
	}
	for (int i = 0; i < 64; ++i) {
	    int	a = gencode[i];
	    invtranslate[a][ndeg[a]++] = i;
	}
	ChangeCodon = 0;
}

static int evenusage()
{
	mkinvtab();
	for (int i = 0; i < 64; ++i)
	    CodonUsage[i] = 1. / (float) ndeg[gencode[i]];
	return (OK);
}

int getCodonUsage(const char* fname)
{
	char	str[LINE_MAX];
	char*	ps;
	FILE*	fd;
static	char	errmsg[] = "Codon Usage File";

	if (!fname || !*fname) return (evenusage());
	fd = fopen(fname, "r");
	if (!fd) fd = ftable.fopen(fname, "r");
	if (!fd) {
	    ps = getenv(DEF_CDI);
	    if (ps) {
		topath(str, ps);
		strcat(str, fname);
		fd = fopen(str, "r");
	    }
	}
	if (!fd) {
	    topath(str, DEF_CDI_PATH);
	    strcat(str, fname);
	    fd = fopen(str, "r");
	}
	if (!fd) {
	    fputs(errmsg, stderr);
	    fprintf(stderr, "%s not found!\n", fname);
	    return (ERROR);
	}
	for (;;) {
	    if (!fgets(str, LINE_MAX, fd)) {
		fputs(errmsg, stderr);
		fprintf(stderr, "%s incomplete!\n", fname);
		return (ERROR);
	    }
	    if (*str == _LCOMM) flush_line(fd);
	    else	break;
	}
	int	i = 0;
	do {
	    sscanf(str, "%*s %*s %*s %f", CodonUsage + i++);
	    if (i == 64) return (OK);
	} while (fgets(str, LINE_MAX, fd));
/* readerr */
	fputs(errmsg, stderr);
	fprintf(stderr, "%s was incomplete!\n", fname);
	for (i = 0; i < 64; ++i) CodonUsage[i] = Human50CDU[i];
	return (ERROR);
}

int setCodonUsage(int gc)
{
	char	str[MAXL];

	if (gc == 40 || gc == 50 || gc == 60 || gc == 70) {
	    sprintf(str, "human%d.cdi", gc);
	    return getCodonUsage(str);
	}
	if (!progets(str,"Codon Usage file [0/40/50/60/70/filename] : ") ||
	    isBlankLine(str)) return (OK);
	if (isdigit(*str)) {
	    gc = atoi(str);
	    if (gc == 0) return evenusage();
	    sprintf(str, "human%d.cdi", gc);
	}
	return getCodonUsage(str);
}

void PatMat::readPatMat(FILE* fd)
{
	int	t = 0;
	int	skip = 0;
	char	str[MAXL];
	double	maxval = -FLT_MAX;
	double	minval = FLT_MAX;

	vclear(&mmm);
	do {
	    if (!fgets(str, MAXL, fd)) return;
	} while (isBlankLine(str));
	if (sscanf(str, "%d %d %d %d %d %f %f %f %d",
	    &rows, &cols, &offset, &t, &skip, 
	    &mmm.min, &mmm.mean, &mmm.max, &nsupport) < 3 ||
	    rows <= 0 || cols <= 0) return;
	tonic = mmm.min;
	if (-tonic > maxtonic) tonic = -maxtonic; 
	mtx = new double[rows * cols];
	while (skip-- > 0) {
	    int rc;
	    while ((rc = fgetc(fd)) != EOF && rc != '\n') ;
	}
	double*	wk = mtx;
	for (int rc = 0; rc < rows * cols; ++rc, ++wk) {
	    if (fscanf(fd, "%lf", wk) <= 0) {
		fputs("Insufficient data!\n", stderr);
		delete[] mtx; mtx = 0;
		return;
	    }
	    if (*wk > maxval) {maxval = *wk; maxidx = rc;}
	    if (*wk < minval) {minval = *wk; minidx = rc;}
	}
	if (t) swap(rows, cols);
	if (rows % 23 == 0) nalpha = 23;
	else if (rows % 4 == 0) nalpha = 4;
	else	nalpha = rows;
	while ((skip = getc(fd)) != EOF && skip != '\n');
}

PatMat& PatMat::operator=(const PatMat& src)
{
	rows = src.rows; cols = src.cols; 
	offset = src.offset; nalpha = src.nalpha;
	maxidx = src.maxidx; minidx = src.minidx; nsupport = src.nsupport;
	tonic = src.tonic; mmm = src.mmm;
	int	mtxsize = rows * cols;
	if (mtxsize != (src.rows * src.cols)) {
	    delete[] mtx;
	    mtx = new double[mtxsize];
	}
	vcopy(mtx, src.mtx, rows * cols);
	return (*this);
}

PatMat::PatMat(const PatMat& src) :
	maxidx(src.maxidx), minidx(src.minidx), nsupport(src.nsupport), 
	nalpha(src.nalpha), rows(src.rows), cols(src.cols), offset(src.offset), 
	tonic(src.tonic), mmm(src.mmm)
{
	int	mtxsize = rows * cols;
	mtx = new double[mtxsize];
	vcopy(mtx, src.mtx, mtxsize);
}

PatMat::PatMat(const int r, const int c, const int o, float* m)
	: maxidx(0), minidx(0), nsupport(0),
	  rows(r), cols(c), offset(o), tonic(0), mtx(0)
{
	if (r <= 0 || c <= 0) return;
	mmm.min = mmm.mean = mmm.max = 0;
	double	maxval = -FLT_MAX;
	double	minval = FLT_MAX;
	mtx = new double[r * c];
	if (!m) return;
	for (int i = 0; i < r * c; ++i) {
	    mtx[i] = m[i];
	    if (mtx[i] > maxval) {maxval = mtx[i]; maxidx = i;}
	    if (mtx[i] < minval) {minval = mtx[i]; minidx = i;}
	}
	if (rows % 4 == 0) nalpha = 4;
	else if (rows == 23) nalpha = 23;
	else	nalpha = cols;
}

PatMat::PatMat(FILE* fd) :
	maxidx(0), minidx(0), nsupport(0),
	rows(0), cols(0), offset(0), tonic(0), mtx(0)
{
	readPatMat(fd);
}

PatMat::PatMat(const char* fname) :
	maxidx(0), minidx(0), nsupport(0),
	rows(0), cols(0), offset(0), tonic(0), mtx(0)
{
	char	str[LINE_MAX];
	FILE*	fd = 0;
retry:
	if (fname && *fname) {
	    strcpy(str, fname);
	    fname = 0;
	} else
	    progets(str, "pattern file name? ");
	if (!*str || *str == '\n') return;
	fd = ftable.fopen(str, "r");
	if (fd) readPatMat(fd);
	else goto retry;
	fclose(fd);
}

CHAR* PatMat::setredctab(const Seq* sd) const
{
	CHAR*	redctab = 0;
	switch (sd->inex.molc) {
	    case DNA: case RNA: case GENOME:
		if (!sd->inex.cmpc) redctab = ncredctab;
		break;
 	    case TRON:
		if (rows % 23) redctab = tnredctab;
		break;
	}
	return (redctab);
}

void PatMat::increment(const Seq* sd, int pos, const CHAR* redctab)
{
const	CHAR*	ps = sd->at(pos);
const	CHAR*	ts = sd->at(pos + cols);
	double*	ptn = mtx;
	for ( ; ps < ts; ++ps, ptn += rows) {
	    int	k = redctab? redctab[*ps]: (*ps - sd->code->base_code);
	    if (k >= 0 && k < rows) ++ptn[k];
	}
}

float PatMat::pwm_score(const Seq* sd, const CHAR* ps, const CHAR* redctab) const
{
const	double*	ptn = mtx;
const	CHAR*	ts = min((const CHAR*) sd->at(sd->right), ps + cols);
	double	fit = 0;
	for ( ; ps < ts; ++ps, ptn += rows) {
	    int	k = redctab? redctab[*ps]: (*ps - sd->code->base_code);
	    if (k >= 0 && k < rows) fit += ptn[k];
	}
	return (float) fit;
}
   
float* PatMat::calcPatMat(const Seq* sd) const
{
	int	k = sd->right - sd->left;
	int	Mrkv = 0;

const 	CHAR*	redctab = setredctab(sd);
	if (rows == 20) Mrkv = 1;	// 1st-order Markov/
	else if (rows == 84) Mrkv = 2;	// 2nd-order Markov/
	else if (rows == 67) Mrkv = 3;	// Should be MODIFIED !!!
	float*	result = new float[k];
	float*	rest = result;
	float*	last = result + k;
	double	minval = cols * mtx[minidx];
const	CHAR*	aa = sd->at(0);			// left limit
const	CHAR*	zz = sd->at(sd->len - Mrkv);	// right limit
	int	n = sd->left - offset;		// start point
	if (Mrkv <= 1) {
	    for ( ; rest < last; ++n) {
		const	CHAR*	ss = sd->at(n);
		const	CHAR*	tt = sd->at(n + cols);
		if (tt > zz) tt = zz;
		double	fit = 0;
const		double*	ptn = mtx;
		if (n < 0) {ptn -= n * rows; ss = aa;}
		int	q = n + cols >= sd->len;	// number of bad chars
		for (int m = 0; ss < tt; ptn += rows, ++m) {
		    const	CHAR*	rr = ss + sd->many;
		    for ( ; ss < rr; ++ss) {
			k = redctab? redctab[*ss]: (*ss - sd->code->base_code);
			if (k < 0 || k >= nalpha) ++q;
			if (Mrkv && !q) {
			    if (m == 0) fit += ptn[k];
			    int j = redctab? redctab[ss[sd->many]]:
				(ss[sd->many] - sd->code->base_code);
			    if (j < 0 || j >= nalpha) ++q;
			    k = nalpha * k + j + nalpha;
			}
			double	tabval = q? 0: ptn[k];
//			double	tabval = q? minval: ptn[k];
			fit += tabval;
		    }
		}
		*rest++ = fit + tonic;
	    }
	} else if (Mrkv == 2) {
	    for ( ; rest < last; ++n) {
		const	CHAR*	ss = sd->at(n);
		const	CHAR*	tt = sd->at(n + cols);
		if (tt > zz) tt = zz;
		double	fit = 0;
const		double*	ptn = mtx;
		if (n < 0) {ptn -= n * rows; ss = aa;}
		int	q = n + cols >= sd->len;	// number of bad chars
		for (int m = 0; ss < tt; ptn += rows, ++m) {
		    if (sd->many == 1) {	/* in most cases */
			int i = redctab? redctab[*ss]: (*ss - sd->code->base_code);
			k = i;
			const	CHAR*	s1 = ++ss;
			if (i > 3) ++q;				// bad char
			if (m == 0 && q == 0) fit += ptn[k];	// 1st pos
			i = redctab? redctab[*s1]: (*s1 - sd->code->base_code);
			if (i > 3) ++q;
			else if (q == 0) {
			    k = nalpha * k + i;
			    if (m == 0) fit += ptn[k + nalpha];	// 1st order
			}
			++s1;
			i = redctab? redctab[*s1]: (*s1 - sd->code->base_code);
			if (i > 3) ++q;
			else if (q == 0) {
			    k = nalpha * k + i;
			    fit += ptn[k + 20];			// 2nd oder
			}
		    } else {
		      const	CHAR*	rr = ss + sd->many;
		      const	CHAR*	s2 = rr + sd->many;
		      for (const CHAR* s1 = rr; ss < rr; ++ss, ++s1, ++s2) {
			int	i = redctab? redctab[*ss]: (*ss - sd->code->base_code);
			k = i;
			if (i > 3) ++q;
			if (m == 0 && q == 0) fit += ptn[k];	// 1st pos
			i = redctab? redctab[*s1]: (*s1 - sd->code->base_code);
			if (i > 3) ++q;
			else if (q == 0) {
			    k = nalpha * k + i;
			    if (m == 0) fit += ptn[k + nalpha];	// 1st order
			}
			i = redctab? redctab[*s2]: (*s2 - sd->code->base_code);
			if (i > 3) ++q;
			else if (q == 0) {
			    k = nalpha * k + i;
			    fit += ptn[k + 20];			// 2nd oder
			}
		      }
		    }
		}
		if (q) fit = minval;
		*rest++ = fit + tonic;
	    }
	} else if (Mrkv == 3) {
	    for ( ; rest < last; ++n) {
		const	CHAR*	ss = sd->at(n);
		const	CHAR*	tt = sd->at(n + cols);
		if (tt > zz) tt = zz;
const		double*	ptn = mtx;
		if (n < 0) {ptn -= n * rows; ss = aa;}
		const	CHAR*	sss = ss;

		int	flag = 0;
		double	fit = 0;
		for (int m = 0; ss < tt; ptn += rows, ++m) {
		    const	CHAR*	rr = ss + sd->many;
		    for ( ; ss < rr; ++ss) {
			INT	x = (INT) ptn[0];
			INT	y = (INT) ptn[1];
			INT	z = (INT) ptn[2];
			k = redctab? redctab[sss[x]]: (sss[x] - sd->code->base_code);
			int	i = redctab? redctab[sss[y]]: (sss[y] - sd->code->base_code);
			int	j = redctab? redctab[sss[z]]: (sss[z] - sd->code->base_code);
			if (k > 10 || i > 10 || j > 10){
			    flag = 1;
			    break;
			}
			fit += ptn[16 * k + nalpha * i + j + 3];
		    }
		    if (flag == 1) {
			fit = 3 * ptn[minidx];
			break;
		    } 
		}
		*rest++ = fit + tonic;
	    }
	}
	return (result);
}

// read coding potential from file

CodePot::CodePot()
{
	char	str[MAXL+1];
	int	d = 0;
	int	trows = TSIMD * TSIMD;
	double	buf[CP_NTERM];
	FILE*	fd = ftable.fopen(CODEPOT, "r");

	if (!fd) {
	    fatal("CodePotTab can't open!\n", stderr);
	    return;
	}
	if (!fgets(str, MAXL, fd)) goto fail_to_read;
	if (!wordcmp(str, "CodePotTab")) {	/* header line */
	    sscanf(str, "%*s %d %d", &CodePotClmn, &d);
	    if (!fgets(str, MAXL, fd)) goto fail_to_read;
	    if (d % 23) {
		CodePotType = DNA;
		trows = d;
	    } else
		CodePotType = TRON;
	}
	if (CodePotClmn > CP_NTERM) CodePotClmn = CP_NTERM;
	CodePotBuf = new FTYPE[CodePotClmn * trows];
	d = 0;
	for (int i = 0; i < CodePotClmn; d += trows)
	    CodePotTab[i++] = CodePotBuf + d;
	d = 0;
	do {
	    char*	ps = str;
	    if (isalpha(*ps)) ps = cdr(str); 
	    for (int i = 0; ps && i < CodePotClmn; ps = cdr(ps))
		buf[i++] = atof(ps);
	    if (CodePotType == TRON && isalpha(*str) && isalpha(str[1])) {
		d = TSIMD * trccode[*str - 'A'] + trccode[str[1] - 'A'];
		for (int i = 0; i < CodePotClmn; ++i) CodePotTab[i][d] = buf[i];
	    } else {
		for (int i = 0; i < CodePotClmn; ++i) CodePotTab[i][d] = buf[i];
		++d;
	    }
	} while (fgets(str, MAXL, fd));
	fclose(fd);
	return;
fail_to_read:
	delete[] CodePotBuf; CodePotBuf = 0;
	fputs("CodePotTab can't be read in !\n", stderr);
	fclose(fd);
}

static int kmer(const CHAR* s, int m, const CHAR* elements, int k, int* x)
{
	int	c = elements[*s];
	int	d = c;
	if (c >= 4) {
	    d = 0;
	    if (x) *x = 0;
	}
const	CHAR*	t = s + k * m;

	while ((s += m) < t) {
	    c = elements[*s];
	    if (c < 4) {
		d = CP_NTERM * d + c;
		++*x;
	    } else {
		d = 0;
		if (x) *x = 0;
	    }
	}
	return (d);
}

/*	Calculate coding potential from DNA seq based on the 5th-order Markov model */

float* CodePot::calc5MMCodePot(const Seq* sd, int phase) const
{
	CHAR*   redctab = sd->inex.molc == TRON? tnredctab: ncredctab;
	CHAR*   elements = sd->inex.molc == TRON? trelements: ncelements;
	int	k = sd->right - sd->left;
	FTYPE	prf;
const	CHAR*	ss = sd->at(sd->left);
const	CHAR*	tt = sd->at(sd->right);

	float*	result = new float[k];
	float*	rest = result;
	int	x = 0;
	if (sd->left < 0)	*rest++ = 0;	/* fill pat */
	else	ss -= sd->many;
	if (sd->many == 1) {			/* in most cases */
	    int	d = kmer(ss, 1, redctab, 6, &x);	/* dummy when phase != 0 */
	    int	dm = d / CP_NTERM;
	    ss += 6;
	    while (ss < tt) {
		int	c = redctab[*ss++];
		FTYPE	val = 0;
		if (phase == 0) {
		    if (c < 4) {
			val = CodePotTab[2][dm] + CodePotTab[0][d]; /* -1 + 0 */
			d = (CP_NTERM * (dm = d) + c) % 4096;
			++x;
		    } else d = dm = x = 0;
		    if (x >= 6) val += CodePotTab[1][d];			/* +1 */
		    else	val = 0;
		} else {
		    if (x >= 6) {
			if (phase <= CodePotClmn) {
			    val = CodePotTab[phase - 1][d];
			} else if (phase == CP_NTERM + 1 && CodePotClmn >= CP_NTERM) {
			    val = CodePotTab[0][d] - CodePotTab[1][d] 
				- CodePotTab[2][d] - CodePotTab[3][d];
			} else {
			    val = CodePotTab[0][d] - CodePotTab[1][d] - CodePotTab[2][d];
			}
		    }
		    if (c < 4) {
			d = (CP_NTERM * d + c) % 4096;
			++x;
		    } else d = x = 0;
		}
		*rest++ = val;
	    }
	} else {
	    tt -= 6 * sd->many;
	    while (ss < tt) {
		FTYPE	val = 0;
		for (int i = 0; i++ < sd->many; ++ss) {
		    int	d = kmer(ss, sd->many, elements, 6);
		    if (phase == 0) {
			prf = CodePotTab[0][d];
			d = (ss < sd->at(0) - sd->many)? d / CP_NTERM: 
			    kmer(ss - sd->many, sd->many, elements, 6);
			prf += CodePotTab[2][d];
			d = kmer(ss + sd->many, sd->many, elements, 6);
			prf += CodePotTab[1][d];
		    } else if (phase <= CodePotClmn) {
			prf = CodePotTab[phase - 1][d];
		    } else if (phase == CP_NTERM + 1 && CodePotClmn >= CP_NTERM) {
			prf = CodePotTab[0][d] - CodePotTab[1][d] 
			- CodePotTab[2][d] - CodePotTab[3][d];
		    } else {
			prf = CodePotTab[0][d] - CodePotTab[1][d] - CodePotTab[2][d];
		    }
		    val += prf;
		}
		*rest++ = val;
	    }
	}
	while (rest < result + k) *rest++ = 0.;
	return (result);
}

/*	Calculate coding potential from DiTRON (1st-order Markov) model */

float* CodePot::calcDitCodePot(Seq* sd, int phase) const
{
	int	i = sd->right;
	int	k = sd->right - sd->left;
	int	d = sd->len - 3;
	int	n2t = 0;
const	CHAR*	ss = sd->at(sd->left);
const	CHAR*	uu = ss + 3 * sd->many;
const	CHAR*	tt = sd->at(min(i, d));

	if (sd->isdrna()) {sd->nuc2tron(); n2t = 1;}
	if (sd->inex.molc != TRON) {
	    prompt("%s must be DNA or TRON code!\n", sd->sqname());
	    return (0);
	}
	float*	result = new float[k];
	float*	rest = result;
	if (phase == 0 && sd->left < 0) {
	    *rest++ = 0;
	    ss += sd->many;
	    uu += sd->many;
	}
	FTYPE	val = 0;
	while (ss < tt) {
	    for (int i = 0; i++ < sd->many; ++ss, ++uu) {
		FTYPE	prf;
		d = TSIMD * *ss + *uu;
		if (phase == 0) {
		    prf = CodePotTab[0][d];
		    d = TSIMD * ss[-sd->many] + uu[-sd->many];
		    prf += CodePotTab[2][d];
		    d = TSIMD * ss[sd->many] + uu[sd->many];
		    prf += CodePotTab[1][d];
		} else if (phase <= CP_NTERM) {
		    prf = CodePotTab[phase - 1][d];
		} else if (phase == CP_NTERM + 1) {
		    prf = CodePotTab[0][d] - CodePotTab[1][d] 
			- CodePotTab[2][d] - CodePotTab[3][d];
		} else {
		    prf = CodePotTab[0][d] - CodePotTab[1][d] - CodePotTab[2][d];
		}
		val += prf;
	    }
	    *rest++ = val;
	}
	for (tt = sd->at(sd->right); ss < tt; ss += sd->many)
	    *rest++ = 0.;
	if (n2t) sd->tron2nuc(0);
	return (result);
}

float* CodePot::calcPrefCodePot(Seq* sd, int phase) const
{
	if (!CodePotBuf) return (0);
	switch (CodePotType) {
	    case DNA:	return calc5MMCodePot(sd, phase); 
	    case TRON:	return calcDitCodePot(sd, phase);
	    default:	return (0);
	}
}

bool ExinPot::readExinPot(const char* fname)
{
	FILE*	fd = ftable.fopen(fname, "r");
	if (!fd) {
	    prompt("%s can't open!\n", fname);
	    return false;
	}
	char	str[MAXL+1];
	int	exin = -1;
	int	sz = 1;
	FTYPE*  pot = 0;
	if (!fgets(str, MAXL, fd)) goto fail_to_read;
	if (!wordcmp(str, EXONPOT)) exin = 1;	/* header line */
	else if (!wordcmp(str, INTRONPOT)) exin = 0;
	else goto fail_to_read;
	if (sscanf(str, "%*s %*d %d %*f %f %*f %*d %d %d %f", 
	     &size, &avpot, &lm, &rm, &avlen) < 1) goto fail_to_read;
	avlen -= (lm + rm);
	pot = new FTYPE[size];
	for (morder = 0; sz < size; sz *= CP_NTERM) ++morder;
	if (sz != size) goto fail_to_read;
	if (exin) ExonPot = pot;
	else	IntronPot = pot;
	while (fgets(str, MAXL, fd)) {
	    char*	ps = str;
	    if (isalpha(*ps)) ps = cdr(str); 
	    *pot++ = atof(ps);
	}
	fclose(fd);
	return true;
fail_to_read:
	delete[] ExonPot; ExonPot = 0;
	delete[] IntronPot; IntronPot = 0;
	prompt("%s has wrong exinpot format !\n", fname);
	fclose(fd);
	return false;
}

ExinPot::ExinPot(int zZ, const char* fname)
	: morder(0), size(0), lm(0), rm(0), 
	  avpot(0), avlen(0), ExonPot(0), IntronPot(0)
{
	if (zZ & 1 && !readExinPot(fname? fname: EXONPOT)) return;
	if (zZ & 2) readExinPot(fname? fname: INTRONPOT);
}

// calculate "instaneous" exonic or "cumulative" intronic potential

float* ExinPot::calcExinPot(const Seq* sd, bool exon) const
{
	FTYPE*	pot = exon? ExonPot: IntronPot;
	if (!pot) return 0;
	double	acc = 0.;
const	CHAR*	ss = sd->at(sd->left);
const	CHAR*	tt = sd->at(sd->right);
const	CHAR*   redctab = sd->inex.molc == TRON? tnredctab: ncredctab;

	float*	result = new float[sd->right - sd->left];
	float*	rest = result;
	if (sd->left < 0)	*rest++ = 0;	/* fill pat */
	else	--ss;
	INT	c = redctab[*ss];
	int	x = c < 4;
	int	d = x? c: 0;
	while (++ss < tt) {
	    if ((c = redctab[*ss]) < 4) {
		++x;
	        d = (CP_NTERM * d + c) % size;
	    } else {
		d = x = 0;
	    }
	    float	pt = (x < morder)? 0: pot[d];
	    *rest++ = exon? pt: acc += pt;
	}
	return (result);
}

VTYPE ExinPot::intpot(const EXIN* b5, const EXIN* b3) const
{
	b5 += lm;
	b3 -= rm;
	if (b5 >= b3) return (0);
	return  (VTYPE) (b3->sigI - b5->sigI);	// unnormalize
//	return  (VTYPE) (avlen * (b3->sigI - b5->sigI) / (b3 - b5)):
}

