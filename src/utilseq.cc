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

#include "seq.h"
#include "eijunc.h"

enum Triplet {AAA = 0, AGA = 8, AGG = 10, AUA = 12, 
	CUA = 28, CUC, CUG, CUU,
	UAA = 48, UAG = 50, UCA = 52, UGA = 56, UUA = 60};

static	const	char	DEF_CDI_PATH[] = "/data2/seqdb/blast/cdi";
static	const	char	DEF_CDI[] = "BLASTCDI";
static	const	int	PassCom = 0x08;
static	const	int	MAXICD = 4;
static	const	int	MAXTCD = 16;

inline	bool	isNucAmb(int c) {return (c != A && c != C && c != G && c != T); }

CHAR gencode[64] = {
	LYS,ASN,LYS,ASN,THR,THR,THR,THR,ARG,SER,ARG,SER,ILE,ILE,MET,ILE,
	GLN,HIS,GLN,HIS,PRO,PRO,PRO,PRO,ARG,ARG,ARG,ARG,LEU,LEU,LEU,LEU,
	GLU,ASP,GLU,ASP,ALA,ALA,ALA,ALA,GLY,GLY,GLY,GLY,VAL,VAL,VAL,VAL,
	TRM,TYR,TRM,TYR,SER,SER,SER,SER,TRM2,CYS,TRP,CYS,LEU,PHE,LEU,PHE
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
		gencode[UAA] = gencode[UAG] = TRM; gencode[UGA] = TRM2;
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
	    if ((gencode[i] == TRM || gencode[i] == TRM2) && k++ < MAXTCD)
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

/*************************************

* DFA to find ORF 

fromic == 0: Term< --- >Term
fromic == 1: <ATG  --- >Term
fromic == 2: ATG - AG< --- >GT - Term
fromic == 3: <ATG - AG --- >GT - Term
fromic == 4: ATG - AG< --- GT - >Term
fromic == 5: <ATG - AG --- GT - >Term

*************************************/

ORF* Seq::getorf() const
{
static	const	int	tm[14][4] = {
//       A  C  G  T
	{1, 0, 2, 3},	//  0: C
	{1, 0, 4, 5},	//  1: A
	{1, 6, 2, 7},	//  2: G
	{8, 0, 9, 3},	//  3: T
	{1, 6, 2, 7},	//  4: AG
	{8, 0,10, 3},	//  5: AT
	{1, 0, 2, 3},	//  6: GC
	{8, 0, 9, 3},	//  7: GT
	{11,0,12, 5},	//  8: TA
	{13,6, 2, 7},	//  9: TG
	{13,6, 2, 7},	// 10: ATG
	{1, 0, 4, 5},	// 11: TAA
	{1, 6, 2, 7},	// 12: TAG
	{1, 0, 4, 5}};	// 13: TGA

	ORF	cur[3];
	ORF*	orfs;
	int	f;			// frame
	int	n = left;
	for (f = 0; f < 3; ++f) {
	    cur[f].pos = n;
	    cur[f].frm = fromic? 0: 1;	// close/open
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
		  case 4:		// AG
		    if (fromic > 1) {	// exon start
			for (f = 0; f < 3; ++f) {
			    if (!cur[f].frm || (fromic % 2 == 0 && cur[f].frm == 1)) {
				cur[f].pos = n + 1;
				cur[f].frm = 2;
			    }
			}
		    }
		    break;
		  case 6:		// GC
		  case 7:		// GT
		    if (fromic > 1) {	// exon end
			for (f = 0; f < 3; ++f) {
			    if (cur[f].frm) {
				cur[f].len = n - 1 - cur[f].pos;
				if (fromic == 2 || fromic == 3)
				    cur[f].frm = 4;	// mark
			    }
			}
		    }
		    break;
		  case 10:		// ATG
		    f = (n - 2) % 3;
		    if (!cur[f].frm) {
			cur[f].pos = n - 2;
			cur[f].frm = 1;
		    }
		    break;
		  case 11:		// TAA
		    if (gencode[UAA] != TRM) break;
		  case 12:		// TAG
		    if (gencode[UAG] != TRM) break;
		  case 13:		// TGA
		    if (gencode[UGA] != TRM2) break;
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
		    if (state == 12 && fromic > 1) {	// AG
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
	    } else {			// ambiguous base
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
	    if ((aa == TRM || aa == TRM2) && (at_term == STOP ||
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
		    if ((aa == TRM || aa == TRM2) && (at_term == STOP ||
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
// readerr
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
	mtx = new float[rows * cols];
	while (skip-- > 0) {
	    int rc;
	    while ((rc = fgetc(fd)) != EOF && rc != '\n') ;
	}
	float*	wk = mtx;
	for (int rc = 0; rc < rows * cols; ++rc, ++wk) {
	    if (fscanf(fd, "%f", wk) <= 0) {
		fputs("Insufficient data!\n", stderr);
		delete[] mtx; mtx = 0;
		return;
	    }
	    if (*wk < min_elem) min_elem = *wk;
	}
	if (t) std::swap(rows, cols);
	if (rows % 23 == 0) nalpha = 23;
	else if (rows % 4 == 0) nalpha = 4;
	else	nalpha = rows;
	morder = 0;
	for (int d = nalpha; d < rows; d = d * (d + 1))
	    ++morder;
	while ((skip = getc(fd)) != EOF && skip != '\n');
}

void PatMat::readBinPatMat(FILE* fd)
{
	if (fread(this, sizeof(PatMat), 1, fd) != 1)
	    fatal("incompatible !\n");
	if (transvers) std::swap(rows, cols);
	tonic = mmm.min;
	if (-tonic > maxtonic) tonic = -maxtonic;
	INT	dszie = rows * cols;
	mtx = new float[dszie];
	if (fread(mtx, sizeof(float), dszie, fd) != dszie)
	    fatal("incompatible !\n");
}

PatMat& PatMat::operator=(const PatMat& src)
{
	rows = src.rows; cols = src.cols; 
	offset = src.offset; nalpha = src.nalpha;
	nsupport = src.nsupport;
	tonic = src.tonic; mmm = src.mmm;
	int	mtxsize = rows * cols;
	if (mtxsize != (src.rows * src.cols)) {
	    delete[] mtx;
	    mtx = new float[mtxsize];
	}
	vcopy(mtx, src.mtx, rows * cols);
	return (*this);
}

PatMat::PatMat(const PatMat& src) :
	rows(src.rows), cols(src.cols), offset(src.offset), 
	tonic(src.tonic), min_elem(src.min_elem), 
	nsupport(src.nsupport), nalpha(src.nalpha), morder(src.morder), 
	mmm(src.mmm)
{
	int	mtxsize = rows * cols;
	mtx = new float[mtxsize];
	vcopy(mtx, src.mtx, mtxsize);
}

PatMat::PatMat(const int r, const int c, const int o, float* m)
	: rows(r), cols(c), offset(o)
{
	if (r <= 0 || c <= 0) return;
	mmm.min = mmm.mean = mmm.max = 0;
	if (rows % 4 == 0) nalpha = 4;
	else if (rows == 23) nalpha = 23;
	else	nalpha = rows;
	for (int d = nalpha; d < rows; d = d * (d + 1))
	    ++morder;
	mtx = new float[r * c];
	if (!m) return;
	vcopy(mtx, m, r * c);
	min_elem = *vmin(mtx, r * c);
}

PatMat::PatMat(FILE* fd, bool binary)
{
	if (binary) readBinPatMat(fd);
	else	readPatMat(fd);
}

PatMat::PatMat(const char* fname)
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
const	char*	dot = strrchr(str, '.');
const	bool	binary = dot && !strcmp(dot, patmat_ext);
	fd = ftable.fopen(str, "rb");
	if (fd) {
	    if (binary) readBinPatMat(fd);
	    else	readPatMat(fd);
	} else goto retry;
	if (fd) fclose(fd);
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
	float*	ptn = mtx;
	for ( ; ps < ts; ++ps, ptn += rows) {
	    int	k = redctab? redctab[*ps]: (*ps - sd->code->base_code);
	    if (k >= 0 && k < rows) ++ptn[k];
	}
}

float PatMat::pwm_score(const Seq* sd, const CHAR* ps, const CHAR* redctab) const
{
const	float*	ptn = mtx;
const	CHAR*	ts = min((const CHAR*) sd->at(sd->right), ps + cols);
	float	fit = 0;
	for ( ; ps < ts; ++ps, ptn += rows) {
	    int	k = redctab? redctab[*ps]: (*ps - sd->code->base_code);
	    if (k >= 0 && k < rows) fit += ptn[k];
	}
	return (float) fit;
}
   
float* PatMat::calcPatMat(const Seq* sd) const
{
	int	k = sd->right - sd->left;
	int	Mrkv = order();

const 	CHAR*	redctab = setredctab(sd);
	float*	result = new float[k];
	float*	rest = result;
	float*	last = result + k;
	float	minval = cols * min_elem;
const	CHAR*	aa = sd->at(0);			// left limit
const	CHAR*	zz = sd->at(sd->len - Mrkv);	// right limit
	int	n = sd->left - offset;		// start point
	if (Mrkv <= 1) {
	    for ( ; rest < last; ++n) {
		const	CHAR*	ss = sd->at(n);
		const	CHAR*	tt = sd->at(n + cols);
		if (tt > zz) tt = zz;
		float	fit = 0;
const		float*	ptn = mtx;
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
			float	tabval = q? 0: ptn[k];
//			float	tabval = q? minval: ptn[k];
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
		float	fit = 0;
const		float*	ptn = mtx;
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
const		float*	ptn = mtx;
		if (n < 0) {ptn -= n * rows; ss = aa;}
		const	CHAR*	sss = ss;

		int	flag = 0;
		float	fit = 0;
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
			fit = 3 * min_elem;
			break;
		    } 
		}
		*rest++ = fit + tonic;
	    }
	}
	return (result);
}

int fname2exin(const char* fname, int& file_type)
{
const	int	ng = static_cast<int>(Iefp::NG);
const	int	write_mode = file_type;
	char	str[MAXL];
	strcpy(str, fname);
	file_type = 0;
	char*	dot = strrchr(str, '.');
	if (dot) {
	    if (!strcmp(dot, text_ext)) {
		*dot = '\0';
		file_type = 1;		// text
		dot = strrchr(str, '.');
	    }
	    for (int exn = 0; exn < ng; ++exn) {
		if (!strcmp(dot, iefp_ext[exn])) {
		    if (!file_type) file_type = 2;
		    return (exn);	// binary
		}
	    }
	    if (write_mode) {
		fputs("Error: extension must be one of:\n\t", stderr);
		for (int exn = 0; exn < ng; ++exn)
		    fprintf(stderr, "%s ", iefp_ext[exn]);
		fputc('\n', stderr);
		exit (1);
	    }
	}

	file_type = 0;
	FILE*	fd = ftable.fopen(fname, "r");
	if (!fd) {
	    prompt(not_found, fname);
	    return (ng);
	}
	if (!fgets(str, MAXL, fd)) return (ng);
	int	exn = 0;
	for ( ; exn < static_cast<int>(Iefp::NG); ++exn) {
	    if (!wordcmp(str, iefp_tid[exn])) {
		file_type = 1;
		break;			// text
	    }
	}
	fclose(fd);
	return (exn);
}

bool ExinPot::readFile(const char* fname)
{
	int	file_type = 0;
	exin = fname2exin(fname, file_type);
	if (file_type == 2)	// binary
	    return readBinary(fname);
	else if (file_type == 0) {
	    prompt(incompatible, fname);
	    return (false);
	}
	FILE*	fd = ftable.fopen(fname, "r");
	if (!fd) {
	    prompt(not_found, fname);
	    return false;
	}
	char	str[MAXL];
	int	sz = 1;
	float*	pot = 0;
	if (!fgets(str, MAXL, fd)) goto fail_to_read;
	if (sscanf(str, "%*s %d %d %*f %f %*f %*d %d %d %f", 
	     &nphase, &ndata, &avpot, &lm, &rm, &avlen) < 1) goto fail_to_read;
	nphase = (nphase >= 3)? 3: 1;
	avlen -= (lm + rm);
	pot = data = new float[dsize()];
	for (morder = -1; sz < ndata; sz *= CP_NTERM) ++morder;
	if (sz != ndata) goto fail_to_read;
	while (fgets(str, MAXL, fd)) {
	    int	i = 0;
	    char* ps = isalpha(*str)? cdr(str): str;
	    for ( ; i++ < 3 && ps; ps = cdr(ps))
		*pot++ = atof(ps);
	    if (i < 3) break;
	}
	fclose(fd);
	if ((pot - data) == dsize()) return true;
fail_to_read:
	fclose(fd);
	prompt(incompatible, fname);
	return false;
}

// from wdfq file specific to nphase == 1

float* ExinPot::getKmers(const char* wdfq, const bool foregrd)
{
	FILE*	fd = fopen(wdfq, "r");
	if (!fd) {
	    prompt(not_found, wdfq);
	    return (0);
	}
const	int	plus = (foregrd && isfpp())? ndata: 0;
	float*	frq = new float[ndata + plus];
	float*	fq = frq + plus;
	if (foregrd) data = frq;
	char	str[MAXL];
	char	kmer[16];
	vclear(fq, ndata);
	int	n = 0;
	while (fgets(str, MAXL, fd)) {
	    int	freq;
	    if (sscanf(str, "%s %d", kmer, &freq) < 2) break;
	    int	k = strlen(kmer);
	    if (--k > morder) break;
	    if (k < morder) continue;
	    *fq++ = freq;
	    if (n++ == ndata) break;
	}
	fclose(fd);
	return (frq);
}

void ExinPot::count_kmers_1(const Seq* sd, float* fq)
{
const	CHAR*	ss = sd->at(sd->left + lm);
const	CHAR*	tt = sd->at(sd->right - rm - morder);
const	CHAR*   redctab = sd->istron()? tnredctab: ncredctab;

	int	x = morder + 1;
	int	w = 0;

	while (ss < tt) {
	    int	c = redctab[*ss++];
	    if (c < 4) {
		if (x) --x;
	        w = (4 * w + c) % ndata;
	    } else {
		w = 0;
		x = morder + 1;
	    }
	    if (!x) ++fq[w];
	}
}

void ExinPot::count_kmers_3(const Seq* sd, float* fq)
{
const	CHAR*	ss = sd->at(sd->left);
const	CHAR*	tt = sd->at(sd->right - 5);
const	CHAR*   redctab = sd->istron()? tnredctab: ncredctab;

	int	x = 6;
	int	w = 0;

	for (int p = 1; ss < tt; p = next_p[p]) {
	    int	c = redctab[*ss++];
	    if (c < 4) {
		if (x) --x;
	        w = (4 * w + c) % ndata;
	    } else {
		w = 0;
		x = 6;
	    }
	    if (!x) ++fq[3 * w + p];
	}
}

float* ExinPot::getKmers(EiJuncSeq* eijseq)
{
const	int	vsize = ndata * nphase;
const	int	plus = isfpp()? vsize: 0;
	data = new float[vsize + plus];
	float*	fq = data + plus;
	vclear(fq, vsize);
	do {
	    Seq*	sd = eijseq->nextseq();
	    if (!sd) continue;
	    if (nphase == 1)	count_kmers_1(sd, fq);
	    else		count_kmers_3(sd, fq);
	} while (eijseq->goahead());
	return (data);
}

float* ExinPot::getKmers(int argc, const char** argv)
{
const	int	vsize = ndata * nphase;
const	int	plus = isfpp()? vsize: 0;
	data = new float[vsize + plus];
	float*	fq = data + plus;
	vclear(fq, vsize);
	Seq	sd(1);
	SeqServer	sqsvr(argc, argv, IM_SNGL);
	InSt	inst;
	while ((inst = sqsvr.nextseq(&sd, 0)) != IS_END) {
	    if (inst == IS_ERR) continue;
	    if (nphase == 1)	count_kmers_1(&sd, fq);
	    else		count_kmers_3(&sd, fq);
	}
	return (data);
}

void ExinPot::reform_1(float* bkg)
{
	int	i = 0;
	float	s = 0;
	float*	frq = bkg? bkg: fbegin();
	float*	fre = frq + ndata;
	float*	pot = bkg? bkg: begin();
const	bool	cpb = ispot();		// conditional probability
const	bool	fpp = bkg? false: isfpp();
	while (frq < fre) {
	    if (cpb) {
		s += ++frq[i++];	// psuedo count 1
		if (i == 4) {
		    for (int j = 0; j++ < 4; ++frq) {
			*pot++ = *frq / s;
			if (fpp) *frq /= total;
		    }
		    s = i = 0;
		}
	    } else    *frq++ /= total;
	}
}

void ExinPot::reform_3()
{
	int	i = 0;
	float	s[3];
	vclear(s, 3);
	float*	frq = fbegin();
	float*	fre = fend();
	float*	pot = begin();
const	bool	cpb = ispot();		// conditional probability
const	bool	fpp = isfpp();
	for (int p = 1; frq < fre; p = next_p[p]) {
	    if (cpb) {
		s[p] += frq[i++] + 1;	// psuedo count 1
		if (i == 12) {
		    for (int j = 0, q = 1; j++ < 12; ++frq, q = next_p[q]) {
			*pot++ = (*frq + 1)/ s[q];
			if (fpp) *frq /= total;
		    }
		    i = 0;
		    vclear(s, 3);
		}
	    } else    *frq++ /= total;
	}
}

void ExinPot::reform(float* bkg)
{
	total = 0;
const	float*	fw = bkg? bkg: fbegin();
const	float*	fe = bkg? (fw + ndata): fend();
	while (fw < fe) total += *fw++;
	if (bkg || nphase == 1) reform_1(bkg);
	else	reform_3();
}

bool ExinPot::makeExinPot(const float* bkg)
{
	if (!bkg) return (false);
	float*	pot = begin();
	float*	frq = fbegin();
	float*	fed = fend();
	int	p = 0;
	while (frq < fed) {
const	    float	freq = *frq++;
	    *pot = log10(*pot / *bkg);
	    if (isfpp()) ess += freq * *pot++;
	    if (++p == nphase) {
		++bkg;
		p = 0;
	    }
	}
	if (isfpp()) ess *= 1000;	// mean score per 1000 bp
	return (true);
}

bool ExinPot::readBinary(const char* fname)
{
	FILE*	fd = fopen(fname, "rb");
	if (!fd) {
	    prompt(not_found, fname);
	    return (false);
	}
	if (fread(this, sizeof(ExinPot), 1, fd) != 1) {
	    prompt(incompatible, fname);
	    return (false);
	}
	data = new float[dsize()];
	if (fread(data, sizeof(float), dsize(), fd) != size_t(dsize())) {
	    prompt(incompatible, fname);
	    return (false);
	}
	fclose(fd);
	return (true);
}

bool ExinPot::writeBinary(const char* oname)
{
	int	file_type = 0;
	int	exn = fname2exin(oname, file_type);
	if (file_type == 0) exn = exin;
	else	std::swap(exn, exin);
	char	str[MAXL];
const	char*	fn = file_type? oname: add_ext(oname, iefp_ext[exin], str);
	FILE*	fd = fopen(fn, "wb");
	if (!fd) {
	    prompt(no_file, str);
	    return (false);
	}
	if (fwrite(this, sizeof(ExinPot), 1, fd) != 1) {
	    prompt(no_file, str);
	    return (false);
	}
	float*	pot = begin();
	if (fwrite(pot, sizeof(float), dsize(), fd) != size_t(dsize())) {
	    prompt(no_file, str);
	    return (false);
	}
	fclose(fd);
	std::swap(exn, exin);
	return (true);
}

// calculate "instaneous" exonic or "cumulative" intronic potential

float* ExinPot::calcScr_1(const Seq* sd, float* scr) const
{
const	float*	pot = begin();
const	int	len = sd->right - sd->left;
const	CHAR*	ss = sd->at(sd->left);
const	CHAR*	tt = sd->at(sd->right);
const	CHAR*   redctab = sd->istron()? tnredctab: ncredctab;
const	int	kk = morder + 1;

	if (sd->left >= 0) --ss;
	float*	result = new float[len];
	float*	rest = result - morder;
	int	x = kk;
	int	w = 0;
const	bool	epot = (exin == static_cast<int>(Iefp::EP));
	double	acc = 0.;			// cumulative
	for ( ; ss < tt; ++rest) {
	    INT	c = redctab[*ss++];
	    if (c < 4) {
		if (x) --x;
	        w = (4 * w + c) % ndata;
	    } else {
		w = 0;
		x = kk;
	    }
const	    float	pt = x? 0: pot[w];
	    if (rest >= result) {
		acc += pt;
		*rest = epot? pt: float(acc);
	    }
	}
	vclear(rest, result + len - rest);
	if (scr)
	    *scr = epot? float(acc): (result[len - 1 - rm] - result[lm]);
	return (result);
}

// Calculate coding potential from DNA seq based on the 5th-order Markov model

float* ExinPot::calcScr_3(const Seq* sd, float* scr) const
{
	CHAR*   redctab = sd->inex.molc == TRON? tnredctab: ncredctab;
const	int	len = sd->right - sd->left;
const	CHAR*	ss = sd->at(sd->left);
const	CHAR*	tt = sd->at(sd->right);

const	float*	cdp = begin();
	float*	result = new float[len];
	float*	rest = result - 5;
const	int	kk = morder + 1;
	int	x = kk;
	int	buf[3] = {0, 0, 0};
	int	w = 0;
double	acc = 0.;
	for (int p = 1; ss < tt; p = next_p[p], ++rest) {
	    int	c = redctab[*ss++];
	    if (c < 4) {
		buf[p] = 3 * (w = (4 * w + c) % ndata);
		if (x) --x;
	    } else {
		w = 0;
		x = kk;
	    }
	    float	val = 0;
	    if (!x) {
	 	val += cdp[buf[next_p[p]] + 2];	// -1
		val += cdp[buf[prev_p[p]]];	// 0
		val += cdp[buf[p] + 1];		// +1
	    }
	    if (rest >= result) {
		*rest = val;
		if (scr && p == 1) acc += val;
	    }
	}
	vclear(rest, result + len - rest);
	if (scr) *scr = float(acc);
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

