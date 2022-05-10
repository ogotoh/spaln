/*****************************************************************************
*
*	Read exon/intron sequences according to coordinates in xxx.eij
*
*	Osamu Gotoh, ph.D.      (-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.      (2001-)
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

#include "eijunc.h"

#ifdef MAIN

void usage()
{
	fputs("Usage:\n", stdout);
	fputs("\teijunc -[3|5|b|c|e|i|j] [other otptions] xxx.eij\n", stdout);
	fputs("Options:\n", stdout);
	fputs("\t-3:\tacceptor\n", stdout);
	fputs("\t-5:\tdonor\n", stdout);
	fputs("\t-a:\tinclude sequences with ambiguous nucleotide\n", stdout);
	fputs("\t-b:\tbranch point\n", stdout);
	fputs("\t-c:\tcDNA\n", stdout);
	fputs("\t-dS:\tgenome dbs\n", stdout);
	fputs("\t-e:\texon\n", stdout);
	fputs("\t-h:\tprint help\n", stdout);
	fputs("\t-i:\tintron\n", stdout);
	fputs("\t-j:\tjoin ends\n", stdout);
	fputs("\t-oS:\toutput (stdout)\n", stdout);
	fputs("\t-wN:\twing size (bp)\n", stdout);
	exit (1);
}

int main(int argc, const char* argv[])
{
	const	char*	pn;
	const	char*	gnm = 0;
	const	char*	outf = 0;
	int	wing = 0;
	Eijmode	mode = DONOR;
	FILE*	fo = stdout;
	const	char*	ext = 0;
	bool	cdna = false;
	bool	no_amb = true;

	while (--argc && **++argv == '-') {
	  switch (argv[0][1]) {
	    case 'd': pn = getarg(argc, argv);
		gnm = pn? pn: ""; break;
	    case '3': mode = ACCPTR; ext = ".SP3"; break;
	    case '5': mode = DONOR; ext = ".SP5"; break;
	    case 'a': no_amb = false; break;
	    case 'b': mode = BRANCH; break;
	    case 'c': cdna = true;
	    case 'e': mode = EXON; break;
	    case 'h': usage();
	    case 'i': mode = INTRON; break;
	    case 'j': mode = JOIN; ext = ".SPb"; break;
	    case 'o': outf  = getarg(argc, argv);
		if (outf) fo = fopen(outf, "w");
		if (!fo) fatal("Can't write to %s\n", outf);
		break;
	    case 's': pn = getarg(argc, argv);
		if (pn) setdfn(pn); break;
	    case 'w': pn = getarg(argc, argv);
		if (pn) wing = atoi(pn); break;
	    default: break;
	  }
	}

	EiJuncSeq eijseq(mode, *argv, gnm, wing);
	int	nrcd = 0;
	setup_output(4);
	do {
	    Seq*	sd = eijseq.nextseq();
	    if (!sd || (no_amb && sd->inex.ambs)) continue;
	    if (!ext) {
		if (!cdna)
		    fprintf(fo, ">%s.%d\n", sd->sqname(), ++nrcd);
		else if (eijseq.new_entry)
		    fprintf(fo, ">%s %d\n", eijseq.transcript,++nrcd);
		sd->typeseq(fo);
	    } else if (sd->len == eijseq.width) ++nrcd;
	} while (eijseq.goahead());

	if (ext) {
	    if (!outf) {
		strcpy(eijseq.genspc + 8, ext);
		outf = eijseq.genspc;
	    }
	    fprintf(fo, ">%s [%d]\n\n", outf, nrcd);
	    eijseq.reset();
	    nrcd = 0;
	    while (eijseq.goahead()) {
	        Seq*	sd = eijseq.nextseq();
		if (!sd || sd->inex.ambs || sd->len != eijseq.width) continue;
		fputs("      1 ", fo);
		sd->typeseq(fo, true);
		fprintf(fo, "| %s.%d\n", sd->sqname(), ++nrcd);
	    }
	}
	return (0);
}

#endif	// MAIN

EiJuncSeq::EiJuncSeq(Eijmode mode, const char* eij, const char* gnmdb, int w)
	: eori(mode), wing(w)
{
	if (eij) {
	    eijfd = fopen(eij, "r");
	    if (!eijfd) fatal(not_found, eij);
	}

	switch (eori) {
	  case INTRON: case EXON: setlpw(60); break;
	  case BRANCH: setlpw(50); 
		wms1 = (wing? wing: def_brch) - 1; wing = 0; break;
	  case JOIN: 
		if (!wing) wing = 30;
		wms1 = wing - 1; 
		setlpw(width = 4 * wing);
		break;
	  default: 
		if (!wing) wing = def_wing; 
		wms1 = wing - 1; 
		setlpw(width = 2 * wing);
		break;
	}

// prepare dbs and background

	*str = *transcript = '\0';
	while (fgets(str, MAXL, eijfd) && (*str == '#' || *str == 0));
	char	genome[MAXL];
	if (gnmdb) strcpy(genome, gnmdb);
	if (!gnmdb || !*genome) {
	    strncpy(genome, str, 8);
	    strcpy(genome + 8, "_g");
	}
	*genome = tolower(*genome);
	strncpy(genspc, genome, 8);
	dbs = new DbsDt(genome);
	if (gnmdb && !dbs) fatal(no_file, genome);
	eijseq = new Seq(1);
}

Seq* EiJuncSeq::nextseq()
{
static	const	char	seqargfrm[] = "%s %d %d %c";
static	const	char	seqargfrm2[] = "%s %d %d %d %d %c";
	Strlist	stl(str, stddelim);
	new_entry = strcmp(transcript, stl[Ref_Column]);
	if (new_entry) strcpy(transcript, stl[Ref_Column]);
	char*	chr = stl[0];
	int	l = atoi(stl[2]);
	int	r = atoi(stl[3]);
	bool	rvs = *stl[1] == '-';
	char	sgn = rvs? '<': ' ';
	if (rvs) std::swap(l, r);
	eijseq->tlen = l;
	int	ll, lr, rl, rr;
	if (rvs && eori != INTRON) {
	    ll = r - wms1;
	    lr = r + wing;
	    rl = l - wing;
	    rr = l + wms1;
	    if (eori == BRANCH && rr > r) rr = r;
	    if (rl < 0) return (0);
	} else {
	    ll = l - wing;
	    lr = l + wms1;
	    rl = r - wms1;
	    if (eori == BRANCH && rl < l) rl = l;
	    rr = r + wing;
	    if (lr < 0) return (0);
	}
	if (eori == EXON) {
	    if (new_entry) {
		if (rvs) {
		    ll = atoi(stl[6]);
		    rr = atoi(stl[3]) - 1;
		} else {
		    ll = atoi(stl[5]);
		    rr = atoi(stl[2]) - 1;
		}
	    } else {
		if (rvs) {
		    ll = atoi(stl[2]) + 1;
		    rr = atoi(stl[5]);
		} else {
		    ll = atoi(stl[3]) + 1;
		    rr = atoi(stl[6]);
		}
	    }
	}

	char	seqarg[MAXL] = "$";
	char	sid[MAXL];
	char*	seqid = dbs? seqarg + 1: seqarg;
	switch (eori) {
	  case DONOR:
	    sprintf(seqid, seqargfrm, chr, ll, lr, sgn);
	    sprintf(sid, "%s%d", chr, l);
	    break;
	  case ACCPTR: case BRANCH:
	    sprintf(seqid, seqargfrm, chr, rl, rr, sgn);
	    sprintf(sid, "%s%d", chr, r);
	    break;
	  case JOIN:
	    if (rvs)
		sprintf(seqid, seqargfrm2, chr, rl, rr, ll, lr, sgn);
	    else
		sprintf(seqid, seqargfrm2, chr, ll, lr, rl, rr, sgn);
	    sprintf(sid, "%s%d.%d", chr, l, r);
	    break;
	  case INTRON: case EXON:
	    sprintf(seqid, seqargfrm, chr, ll, rr, sgn);
	    sprintf(sid, "%s%d.%d", chr, l, r);
	    break;
	}
	Seq*	tmp = eijseq->getseq(seqarg, dbs);
	if (tmp) tmp->tlen = rl - l;	// from 5' end to left
	return (tmp);
}

