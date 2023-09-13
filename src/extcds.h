/*****************************************************************************
*
*	Extract coding regions from GenBank flat file
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

#ifndef  _EXTCDS_H_
#define  _EXTCDS_H_

static	const	INT	SEGMAX = 100;
static	const	INT	EXONMAX = 100;
static	const	int	JOINMAX = 1024;
static	const	INT	LLINE = 80;
static	const	INT	TLEN = 2;
static	const	INT	CDSTYPES = 3;
static	const	int	EOE = -1;
static	const	INT	DB_TYPE_PRO = 0x78857a4f;// Magic # for a protein sequence database
static	const	INT	DB_TYPE_NUC = 0x788325f8;// Magic # for a nt. sequence database
static	const	INT	AAFORMAT = 3; // Latest a.a. database format ID number
static	const	INT	NTFORMAT = 6; // Latest nt. database format ID number

enum {CDS, EXON, MRNA, UNDEF};
enum {NOCDS, ENTRY, TRANS};

struct EXTYPE {
	INT	segno:	8;
	INT	join:	1;
	INT	undef:	1;
	INT	cmpl:	1;
	INT	lpart:	1;
	INT	rpart:	1;
};

struct List {
	EXTYPE	texon;
	long	from, to;
} ;

struct SEGM {
	long	seqptr;
	char	acc[LLINE];
};

struct BLAST_DB_H {
	long	dbtype;
	long	dbformat;
	long	title_len;
	char	title[TLEN];
	long	n_seqs;
	long	total_len;
	long	max_len;
};

#endif
