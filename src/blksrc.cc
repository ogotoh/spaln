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
*	Osamu Gotoh, Ph.D.	(2003-)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/

#include "aln.h"
#include "utilseq.h"
#include "blksrc.h"

bool	ignoreamb = false;

static	const	int	TFACTOR = 100;
static	const	int	NCAND2PHS = 10;
static	const	INT	ddelim = SEQ_DELIM + (SEQ_DELIM << 4);
static	const	int	genlencoeff = 38;
static	const	INT	MinMaxGene = 16384;
static	const	INT	NoRetry = 2;

static	const	char	dir[2] = {'>', '<'};
static	const	char	emssg[] = "%s is incompatible !\n";
static	const	char	mfmt[] = "Gd:%u No:%u My:%u MS:%d Mb:%u Tw:%u Tl:%u %6.2f %6.2f %6.2f\n";
static	const	char	writeerrmssg[] = "Fail to write to %s !\n";

//		  elem tuple bitpat2 tabsize bitpat shift blklen maxgene nbitpat afact
static	BlkWcPrm	wcp = {0, 0, 1, 0, 1, 0, 0, 0, 0, 0};
static	BlkWcPrm	wcp_af = {20, 5, 0, 3200000, 0, 5, 4096, 0, 1, 10};
static	BlkWcPrm	wcp_ax = {18, 5, 0, 1889568, 0, 5, 4096, 0, 4, 10};
static	BlkWcPrm	wcp_cf = {4, 8, 0, 65536, 0, 8, 4096, 0, 1, 10};
static	BlkWcPrm	wcp_cx = {4, 8, 3255, 65536, 7573, 8, 4096, 0, 5, 10};

static	int	gene_rng_max_extend = -1;
static	float	hfact = 1.25;
static	int	gratio = 10;	// ratio of average gene and cDNA
static	const	char*	WriteFile = 0;
static	SHORT	Ncand = NCAND2PHS;	// max number of candidates to pass to 2nd phase
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
static	float	RbsFact = FQUERY;
static	int	Nascr = 2;
static	INT	MinOrf = 30;	// nt
static	float	ild_up_quantile = 0.996;
static	BLKTYPE	MaxBlock = 0;
static	BLKTYPE	ExtBlock = 0;
static	BLKTYPE	ExtBlockL = 0;

static	void	setupbitpat(int molc, size_t gnmsz);
static	bool	newer(char* str, const char** argv);
static	bool	testblk(char* str, int molc, const char** argv);
static	void	idx2SeqLenStat(const char* fn, SeqLenStat& sls);
static	bool	testlut(char* str, int molc, const char** argv);

int setQ4prm(const char* ps, const char* ss)
{
const	char*	val = ps + 1;

	int	rv = 0;
	if (ss && *ss == '-') ss = 0;
	if (ss && !*val) {val = ss; rv = 1;}
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
		ExtBlock = gene_rng_max_extend;
		break;
	    case 'G':	// Maximum gene length
		if (*val) wcp.MaxGene = ktoi(val);
		break;
	    case 'M':	// Maximum gene length
		if (*val) OutPrm.MaxOut2 = atoi(val);
		break;
	    case 'Q':	// Maximum gene length
		if (*val) ild_up_quantile = atof(val);
		break;
	    case 'S':	// activate salvage mode
		algmode.slv = 1; break;
	    case 'W':	// Write block information to the file
		WriteFile = val;
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
	    case 'r':	// Minum size of orf
		if (*val) MinOrf = atoi(val);
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
	return (rv);
}

/*****************************************************
	MakeDbs
*****************************************************/

MakeDbs::MakeDbs(const char* dbn, int molc)
	: isaa(molc == PROTEIN)
{
	entryprv[0] = '\0';
	dbname = strrealloc(0, dbn);
#if USE_ZLIB
	if (OutPrm.gzipped) {
	    gzseq = gzopenpbe(OutPrm.out_file, dbname, SGZ_EXT, "w", 2);
	    gzidx = gzopenpbe(OutPrm.out_file, dbname, IDZ_EXT, "w", 2);
	    gzent = gzopenpbe(OutPrm.out_file, dbname, ENZ_EXT, "w", 2);
	} else
#endif
	{
	    fseq = fopenpbe(OutPrm.out_file, dbname, SEQ_EXT, "w", 2);
	    fidx = fopenpbe(OutPrm.out_file, dbname, IDX_EXT, "w+", 2);
	    fent = fopenpbe(OutPrm.out_file, dbname, ENT_EXT, "w+", 2);
	}
	vclear(&rec);
	fgrp = fopenpbe(OutPrm.out_file, dbname, GRP_EXT, "w", 2);
}

static	DbsRec*	rbuf;
static	char*	cbuf;

static int cmpkey(const INT* a, const INT* b)
{
	return(strcmp(cbuf + rbuf[*a].entptr, cbuf + rbuf[*b].entptr));
}

template <typename file_t>
void MakeDbs::read_ent(file_t fd)
{
	size_t	flen = ftell(fd);
	cbuf = new char[flen];
	if (gzent) {
	    fclose(fd);
	    gzent = gzopenpbe(OutPrm.out_file, dbname, ENZ_EXT, "r", 2);
	    if (fread(cbuf, sizeof(char), flen, gzent) != flen)
		fatal("Corrupted entry file !\n");
	} else {
	    rewind(fd);
	    if (fread(cbuf, sizeof(char), flen, fd) != flen)
		fatal("Corrupted entry file !\n");
	}
}

template <typename file_t>
void MakeDbs::read_idx(file_t fd)
{
	if (gzidx) {
	    fclose(fd);
	    gzidx = gzopenpbe(OutPrm.out_file, dbname, IDZ_EXT, "r", 2);
	    if (fread(rbuf, sizeof(DbsRec), recnbr, gzidx) != recnbr)
		fatal("Corrupted index file !\n");
	} else {
	    rewind(fd);
	    if (fread(rbuf, sizeof(DbsRec), recnbr, fd) != recnbr)
		fatal("Corrupted index file !\n");
	}
}

template <typename file_t>
void MakeDbs::write_idx(file_t fd)
{
	if (fwrite(&rec, sizeof(DbsRec), 1, fd) != 1)
	    fatal(write_error, "dbs.rec");
	++recnbr;
}

template <typename file_t>
void MakeDbs::write_odr(file_t fd, const INT* order, const char* str)
{
	if (fwrite(order, sizeof(INT), recnbr, fd) != recnbr)
	    fatal(write_error, str);
	fclose(fd);
}

void MakeDbs::mkidx()
{
	if (!cridxf) return;		// has been sorted
	if (gzent)	read_ent(gzent);
	else 		read_ent(fent);
	rbuf = new DbsRec[recnbr];
	if (gzidx)	read_idx(gzidx);
	else		read_idx(fidx);
	INT*	order = new INT[recnbr];
	for (INT i = 0; i < recnbr; ++i) order[i] = i;
	qsort((UPTR) order, recnbr, sizeof(INT), (CMPF) cmpkey);
	char	str[LINE_MAX];
	FILE*	fodr = 0;
#if USE_ZLIB
	gzFile	gzodr = 0;
	if (OutPrm.gzipped)
	    gzodr = gzopenpbe(OutPrm.out_file, dbname, ODZ_EXT, "w", 2, str);
	else
	    fodr = fopenpbe(OutPrm.out_file, dbname, ODR_EXT, "w", 2, str);
	if (fodr)	write_odr(fodr, order, str);
	else		write_odr(gzodr, order, str);
#else
	fodr = fopenpbe(OutPrm.out_file, dbname, ODR_EXT, "w", 2, str);
	write_odr(fodr, str);
#endif
	delete[] cbuf;
	delete[] rbuf;
	delete[] order;
}

template <typename file_t>
int MakeDbs::write_recrd(file_t fd, int c)
{ 
	putseq(SEQ_DELIM);
	if (rec.seqlen) {
	    if (gzidx)	write_idx(gzidx);
	    else	write_idx(fidx);
	}
	if (c != EOF) {
	    if (gzent)	rec.entptr = ftell(gzent);
	    else	rec.entptr = ftell(fent);
	    if (gzseq)	rec.seqptr = ftell(gzseq);
	    else	rec.seqptr = ftell(fseq);
	    char*	ps = entrystr;
	    char*	pe = ps;
	    char*	pb = ps;
	    while ((c = fgetc(fd)) != EOF && !isspace(c)) {
		*ps++ = c;
		if (c == '|') {
		    pb = pe;
		    ps[-1] = '\0';
		    pe = ps = entrystr;
		}
	    }
	    if (pe == ps) pe = pb;	// last was '|'
	    else	*ps = '\0';
	    if (wordcmp(pe, entryprv) < 0) cridxf = true;
	    strcpy(entryprv, pe);
	    if (gzent) {
		fputs(pe, gzent); fputc('\0', gzent);
	    } else {
		fputs(pe, fent); fputc('\0', fent);
	    }
	}
	rec.seqlen = b = 0;
	return (c);
}

/*****************************************************
	MakeBlk
*****************************************************/

ContBlk::ContBlk() :
	ConvTS(0), WordNo(0), WordSz(0), ChrNo(0), glen(0)
{
	clean();
	VerNo = BsVersion;
}

ContBlk::~ContBlk()
{
	delete[] Nblk; 
	delete[] blkp; 
	delete[] blkb;
	delete[] wscr; 
	delete[] ChrID;
}

// find words in each block

void Chash::countBlk(ContBlk& cntblk)
{
	KVpair<INT, int>* h = hash;
	for ( ; h < hz; ++h) {
	    if (h->val) {
		cntblk.Nblk[h->key]++;
		h->val = 0;	// reset hash
	    }
	}
}

void Chash::registBlk(INT m, ContBlk& cntblk)
{
	for (KVpair<INT, int>* h = hash; h < hz; ++h) {
	    if (h->val) {
		INT*	wrk = cntblk.blkp[h->key];
		if (wrk) {
		    INT	nblk = cntblk.Nblk[h->key]++;
		    wrk[nblk] = m;
		}
		h->val = 0;	// reset hash
	    }
	}
}

// Block class

Block::Block(const Seq* sq, INT blklen, bool mk_blk) 
	: WordTab(sq, wcp.Ktuple, wcp.Nshift, wcp.Nalpha, ConvPat, 
	  wcp.BitPat, wcp.Bitpat2, wcp.Nbitpat, 0, wcp.afact, MinOrf),
	  sd(sq), isaa(sd->isprotein()), istron(sd->istron()),
	  mkblk(mk_blk), ch(0)
#if M_THREAD
	  , cps(0), prelude(0), margin(0), seq(0), tsq(0), bid(0)
#endif
{
	if (!istron) MinOrf = 0;
	prelude = max_width() * (istron? 3: 1) - 1;
	wcp.TabSize = bpp[0]->TabSize;
//	if (thread_num < 1 && mkblk) 
	    ch = new Chash(int(hfact * wcp.blklen));
	margin = prelude;
	margin += MinOrf;
}

void Block::c2w(INT uc, INT* tcc)	// letter to word
{
	if (uc == Nalpha) {		// ambigous
	    ss = 0;
	    for (INT k = 0; k < Nbitpat; ++k) bpp[k]->flaw();
	    return;
	} else	++ss;
	for (INT k = 0; k < Nbitpat; ++k) {
	    Bitpat_wq*&	bps = bpp[k];
	    INT	w = bps->word(uc);
	    if (bps->flawless()) {
		if (tcc) ++tcc[w];
		int	nw = ss - bps->width;
		if (!(nw % Nshift)) ch->incr(w);
	    }
	}
}

void Block::c2w6(INT uc, INT* tcc)	// letter to word
{
	cc[p] = cc[p + 3] = 0;  	// initialize codon
	if (uc == 4) {			// ambiguous
	    vset(xx, 4U, 3);
	} else {
	    for (int q = 0; q < 3; ++q) {
		cc[q] = ((cc[q] << 2) + uc) & 63;
		cc[q + 3] = (((3 - uc) << 4) + (cc[q + 3] >> 2)) & 63;
		xx[p] >>= 1;
	    }
	}
	if (++p == 3) p = 0;
	int	wqp = w_qp;
	int	wq_base = 0;
	if (++w_qp == wq_size) w_qp = 0;
	for (int q = p; q < 6; q += 3, wqp += wq_size, wq_base += wq_size) {
	    SHORT&	sss = ss6[q];		// orf length in codon
	    int	orf_len = 3 * sss;
	    if (!xx[p]) uc = g2r[cc[q]];
	    if (!xx[p] && bpp[0]->good(uc)) ++sss;
	    else	sss = 0;
	    for (INT k = 0; k < Nbitpat; ++k) {
		Bitpat_wq*&	bps = bpp[k];
		INT	w = BadWord;
		if (sss) {	// coding codon
		    w = bps->word(uc, q);
		    if (bps->flawless(q)) {
			if (tcc) ++tcc[w];
			int	nw = sss - bps->width;
			if (nw < 0 || nw % Nshift) w = BadWord;
		    }
		} else {	// termination codon
		    bps->flaw(q);
		    int	nw = orf_len / 3 - bps->width;
		    if (0 <= nw && orf_len < wq_size) {
			int	sp = nw % Nshift;
			nw = (nw + sp) / Nshift;
			int	qp = wqp - (sp + 1) * 3;
			for ( ; nw-- >= 0; qp -= 3 * Nshift) {
			    if (qp < wq_base) qp += wq_size;
			    w_queue[k][qp] = BadWord;
			}
		    }
		}
		if (wq_size) swap(w, w_queue[k][wqp]);
		if (w != BadWord) ch->incr(w);
	    }
	}
}

void Block::c2w6_pp(int n)	// process words remaining in queue
{
	while (n-- > 0) {
	    int	wqp = w_qp;
	    if (++w_qp == wq_size) w_qp = 0;
	    for (int q = 0; q < 2; ++q, wqp += wq_size) {
		for (INT k = 0; k < Nbitpat; ++k) {
		    INT	w = BadWord;
		    swap(w, w_queue[k][wqp]);
		    if (w != BadWord) ch->incr(w);
		}
	    }
	}
}

// MakeBlk class

MakeBlk::MakeBlk(const Seq* sq, DbsDt* dd, MakeDbs* mdbs, bool mk_blk) :
	Block(sq, wcp.blklen, mk_blk), sd(sq), 
	mkdbs(mdbs), wdbf(dd), pchrid(0), s2r(0), 
	bias(A - 1), s_size(margin + wcp.blklen), deltaa(0), acomp(0)
#if M_THREAD
	, did(0), c_qsize(0), bq(0), blks(0), seqbuf(0)
#endif	// M_THREAD
{
	cntblk.ConvTS = ConvTabSize;
	tcount = new INT[wcp.TabSize];
	vclear(tcount, wcp.TabSize);
	cntblk.Nblk = new SHORT[wcp.TabSize];
	vclear(cntblk.Nblk, wcp.TabSize);
	vclear(&chrbuf);
	int	molc = sd->inex.molc;
	if (!molc && dd) molc = dd->curdb->defmolc;
	defcode = setSeqCode(0, isaa? PROTEIN: DNA);
	if (isaa) bias = ALA - 1;
	if (istron) prepacomp();

	if (dd) {
const	    char*	cvt = isaa? acodon: nucl;
	    s2r = new CHAR[defcode->max_code];
	    s2r[nil_code] = s2r[gap_code] = Nalpha + 1;
	    s2r[defcode->amb_code] = Nalpha;
	    for (int i = defcode->base_code; i < defcode->max_code; ++i)
		s2r[i] = a2r[cvt[i] - 'A'];
	}
}

MakeBlk::~MakeBlk()
{
#if M_THREAD
	delete[] seqbuf;
	if (blks) {
	    for (int i = 0; i < c_qsize; ++i)
		delete blks[i];
	    delete[] blks;
	}
	if (bq) {
	    for (int i = 0; i < c_qsize; ++i)
		delete bq[i];
	    delete[] bq;
	}
#endif	// M_THREAD
	delete[] acomp; delete mkdbs; delete[] s2r;
	delete[] tcount;
}

void MakeBlk::findChrBbound(Block2Chr* pb2c)
{
	double	B = (double) (chrblk(cntblk.ChrNo) - 1);

	vclear(pb2c);
	for (size_t k = 0; k <= cntblk.ChrNo; ++k) {
	    double	offdiag = k * B - cntblk.ChrNo * (chrblk(k) - 1);
	    if (offdiag < pb2c->BClw) pb2c->BClw = offdiag;
	    if (offdiag > pb2c->BCup) pb2c->BCup = offdiag;
	}
	pb2c->BClw /= B;
	pb2c->BCup /= B;
}

template <typename file_t>
void MakeBlk::writeBlkInfo(file_t fd, const char* block_fn)
{
	findChrBbound(&b2c);
	if (!fd) fatal(writeerrmssg, block_fn);
	size_t	chrn = cntblk.ChrNo + 1;
	if (!cntblk.ChrID) ++cntblk.VerNo;
	if (fwrite(&wcp, sizeof(BlkWcPrm), 1, fd) != 1 ||
	    fwrite(&cntblk, sizeof(ContBlk), 1, fd) != 1 ||
	    fwrite(&b2c, sizeof(Block2Chr), 1, fd) != 1 ||
	    (cntblk.ChrID && fwrite(cntblk.ChrID, sizeof(CHROMO), chrn, fd) != chrn) ||
	    fwrite(cntblk.Nblk, sizeof(SHORT), wcp.TabSize, fd) != wcp.TabSize)
		fatal(writeerrmssg, block_fn);
	INT*	blkidx = new INT[wcp.TabSize];
	for (INT w = 0; w < wcp.TabSize; ++w)
	    blkidx[w] = cntblk.blkp[w]? cntblk.blkp[w] - cntblk.blkb + 1: 0;
	if (fwrite(blkidx, sizeof(INT), wcp.TabSize, fd) != wcp.TabSize)
		fatal(writeerrmssg, block_fn);
	delete[] blkidx;
	SHORT*	sblkb = (SHORT*) cntblk.blkb;
	if (fwrite(sblkb, sizeof(SHORT), cntblk.WordSz, fd) != cntblk.WordSz ||
	    fwrite(cntblk.wscr, sizeof(short), wcp.TabSize, fd) != wcp.TabSize ||
	    fwrite(iConvTab, sizeof(CHAR), cntblk.ConvTS, fd) != cntblk.ConvTS)
		fatal(writeerrmssg, block_fn);
	fclose(fd);
}

void MakeBlk::WriteBlkInfo()
{
static  const char* wfmt =
	"#Segs %u, TabSize %u, Words: %lu, GenomeSize %lu, GIDs %d\n";

	if (!WriteFile) fatal("Specify write file !\n");
	prompt(wfmt, chrblk(cntblk.ChrNo) - 1,
	    wcp.TabSize, cntblk.WordNo, cntblk.glen, cntblk.ChrNo);

	SHORT*  sblkb = (SHORT*) cntblk.blkb;
	if (cntblk.BytBlk == 2) {
	    SHORT*	sblk = sblkb;
	    INT*	iblk = cntblk.blkb;
	    for (size_t i = 0; i < cntblk.WordNo; ++i)
		*sblk++ = (SHORT) *iblk++;
	}

	char    block_fn[LINE_MAX];
	strcpy(block_fn, WriteFile);
	char*   dot = strrchr(block_fn, '.');
	if (dot && is_gz(block_fn)) {
	    *dot = '\0';
	    dot = strrchr(block_fn, '.');
	}
	if (dot) {
	    ++dot;
	    if (strcmp(dot, BKA_EXT) && strcmp(dot, BKP_EXT) && strcmp(dot, BKN_EXT)) {
		*--dot = '\0';
		dot = 0;
	    }
	}
	if (!dot) {
	    if (isaa) strcat(block_fn, BKA_EXT);
	    else if (istron) strcat(block_fn, BKP_EXT);
	    else strcat(block_fn, BKN_EXT);
	}

#if USE_ZLIB
	if (!is_gz(block_fn) && OutPrm.gzipped)
	    strcat(block_fn, gz_ext);
	if (is_gz(block_fn)) {
	    gzFile	gzfd = gzopen(block_fn, "wb");
	    if (!gzfd) fatal(no_file, block_fn);
	    writeBlkInfo(gzfd, block_fn);
	    return;
	}
#else
	if (is_gz(block_fn)) fatal(gz_unsupport, block_fn);
#endif
	FILE*   fd = fopen(block_fn, "wb");
	if (!fd) fatal(no_file, block_fn);
	writeBlkInfo(fd, block_fn);
}

static void setupbitpat(int molc, size_t gnmsz)
{
	if (molc == PROTEIN || molc == TRON) {
	    BlkWcPrm&	wcp_t = algmode.crs == 2? wcp_ax: wcp_af;
	    if (wcp.Nalpha == 0) wcp.Nalpha = wcp_t.Nalpha;
	    if (wcp.Ktuple == 0 && !gnmsz) wcp.Ktuple = wcp_t.Ktuple;
	    if (wcp.Bitpat2 == 1) wcp.Bitpat2 = wcp_t.Bitpat2;
	    if (wcp.BitPat == 1) wcp.BitPat = wcp_t.BitPat;
	    if (wcp.blklen == 0 && !gnmsz) wcp.blklen = wcp_t.blklen;
	    if (wcp.Nbitpat == 0) wcp.Nbitpat = wcp_t.Nbitpat;
	    if (wcp.afact == 0) wcp.afact = wcp_t.afact;
	} else {
	    BlkWcPrm&	wcp_t = algmode.crs == 2? wcp_cx: wcp_cf;
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
		if (wcp.blklen > 65536) wcp.blklen = 65536;	// 2**16
	    }
	    if (!wcp.Ktuple) {
		INT	max_kmer = 6;
		double	loggs = log(double(gnmsz));
		switch (molc) {
		  case PROTEIN:	loggs *= 0.30; break;
		  case TRON:	loggs *= 0.36; break;
		  default:	loggs *= 0.59; max_kmer = 16; break;
		}
		wcp.Ktuple = int(loggs);
		if (wcp.Ktuple < 3) wcp.Ktuple = 3;
		if (wcp.Ktuple > max_kmer) {wcp.Ktuple = max_kmer;}
	    }
	    if (!wcp.MaxGene) {
		wcp.MaxGene = int(genlencoeff * 
		sqrt(double (gnmsz)) / 1024 + 1) * 1024;
		if (wcp.MaxGene < MinMaxGene) wcp.MaxGene = MinMaxGene;
	    }
	}
	wcp.TabSize = ipower(wcp.Nalpha, wcp.Ktuple);
	if (wcp.Nshift == 0) wcp.Nshift = wcp.Ktuple;
	if (wcp.Nbitpat > 1 && (!sBitPat || *sBitPat != '1'))
	    sBitPat = DefBitPat[wcp.Ktuple];
	if (sBitPat) {
	    INT	w = 0;
	    wcp.BitPat = bpcompress(sBitPat, &w);
	    if (w > wcp.Ktuple) wcp.Ktuple = w;
	} else {
	    wcp.BitPat = bitmask(wcp.Ktuple);
	    wcp.Nbitpat = 1;
	}
	if (wcp.Nbitpat > 3) {
	    const char*	bp2 = strchr(sBitPat, ',');
	    if (bp2) {
		INT	w = 0;
		wcp.Bitpat2 = bpcompress(bp2 + 1, &w);
		if (w > wcp.Ktuple) wcp.Ktuple = w;
	    } else {
		wcp.Nbitpat = 3;
	    }
	}
}

static bool newer(char* str, const char** argv)
{
	bool    update = !file_size(str);
#if USE_ZLIB
	if (update) {
	    strcat(str, gz_ext);
	    update = !file_size(str);
	}
#endif
	if (update) return (true);	// absent
	time_t  exist_time = modify_time(str);
	for (const char** av = argv; *av; ++av)
	    if ((update = modify_time(*av) > exist_time))
		break;
	return (update);
}

static bool testblk(char* str, int molc, const char** argv)
{
	switch (molc) {
	    case PROTEIN: strcat(str, BKA_EXT); break;
	    case DNA:   strcat(str, BKN_EXT); break;
	    case TRON:  strcat(str, BKP_EXT); break;
	}
	return (newer(str, argv));
}

MakeBlk* makeblock(int argc, const char** argv, int molc, bool mk_blk)
{
	Seq	sd(1);
	setSeqCode(&sd, molc);
	MakeDbs*	makedbs = 0;
	char    str[LINE_MAX];
	if (!*WriteFile) WriteFile = *argv;
	strcpy(str, WriteFile);
	char*   dot = strrchr(str, '.');
	if (is_gz(str)) {
	    *dot = '\0';
	    dot = strrchr(str, '.');
	}
	if (dot)	*dot = '\0';
	else    dot = str + strlen(str);
	strcat(str, ".idx");
	bool    update = newer(str, argv);
	if (update) {   // idx file absent or older
	    *dot = '\0';
	    makedbs = new MakeDbs(str, molc);
	} else {
	    *dot = '\0';
	    update = testblk(str, molc, argv);
	    if (!update && algmode.alg) {
		*dot = '\0';
		update = testlut(str, molc, argv);
	    }
	}
	if (!update) return (0);	// up to date
	if (molc == PROTEIN) {
	    SeqLenStat	sls = {0, 0, 0};
	    if (!makedbs) {		// .idx exists
		*dot = '\0';
		strcat(str, ".idx");
		idx2SeqLenStat(str, sls);
	    } else {			// .idx to be made
		PreScan	ps(&sd, wcp.Nalpha);
		ps.lenStat(argc, argv, sls);
	    }
	    setupbitpat(molc, (size_t) sls.total);
	    wcp.blklen = (INT) sls.maxv;
	} else {
	    size_t	gnmsz = file_size(*argv);	// temporary genome size
	    if (is_gz(*argv)) gnmsz *= 2;
	    setupbitpat(molc, gnmsz);
	}
	MakeBlk*	mb = new MakeBlk(&sd, 0, makedbs, mk_blk);
#if M_THREAD
	if (thread_num >= 1 && molc != PROTEIN) 
	    mb->m_idxblk(argc, argv);
	else
#endif	//	M_THREAD
	mb->idxblk(argc, argv);
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
		cntblk.wscr[w] = wscr;
		avr += wscr;
	    } else cntblk.wscr[w] = 0;
	    if (aaafact > 0) {
		int	p = 0, q = 0;
		for (INT x = w + 1; (q = x % wcp.Nalpha) == 0; x /= wcp.Nalpha)
		    ++p;
		if (p) alc += p * deltaa;	// (q-1)AAA
		alc += acomp[q] - acomp[q-1];	//	q000
	    }
	}
	cntblk.AvrScr = (SHORT) (avr /= m);
	short	MinScr = (short) (avr - TFACTOR * (1 + aaafact) * log((double) wcp.afact));
	if (MinScr < 0) MinScr = 0;

	for (INT w = i = c = m = cntblk.WordNo = 0; w < wcp.TabSize; ++w) {
	    if (cntblk.Nblk[w] == 0) {
		++c;
		cntblk.wscr[w] = -1;
	    } else if (cntblk.wscr[w] > MinScr) {
		++m;
		cntblk.WordNo += cntblk.Nblk[w];
		if (cntblk.Nblk[w] > cntblk.MaxBlk) cntblk.MaxBlk = cntblk.Nblk[w];
	    } else {
		++i;
		cntblk.wscr[w] = 0;
	    }
	}
	if (WriteFile && segn <= USHRT_MAX) {
	    cntblk.BytBlk = 2;
	    cntblk.WordSz = cntblk.WordNo;
	} else {
	    cntblk.BytBlk = 4;
	    cntblk.WordSz = 2 * cntblk.WordNo;
	}
	INT*	wrk = cntblk.blkb = new INT[cntblk.WordNo];
	size_t	x = 0;
	for (INT w = 0; w < wcp.TabSize; ++w) {
	    x += tcount[w];
	    if (cntblk.Nblk[w]) {
		if (cntblk.wscr[w] > MinScr) {
		    cntblk.blkp[w] = wrk;
		    wrk += cntblk.Nblk[w];
		} else
		    cntblk.blkp[w] = 0;
		cntblk.Nblk[w] = 0;	// reset for next round
	    } else  cntblk.blkp[w] = 0;
	}
	prompt(mfmt, m, c, i, MinScr, cntblk.MaxBlk, cntblk.WordNo, x,
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
	for (INT w = i = c = m = cntblk.WordNo = 0; w < wcp.TabSize; ++w) {
	    if (cntblk.Nblk[w]) {
		short	wscr = (short) (TFACTOR * (basescr - log((double) tcount[w] / Nbitpat)));
		if (wscr > MinScr) {
		    ++m;
		    cntblk.WordNo += cntblk.Nblk[w];
		    cntblk.wscr[w] = wscr;
		    avr += wscr;
		    if (cntblk.Nblk[w] > cntblk.MaxBlk) cntblk.MaxBlk = cntblk.Nblk[w];
		} else {
		    ++i;
		    cntblk.wscr[w] = MinScr;
		}
	    } else {
		++c;
		cntblk.wscr[w] = -1;
	    }
	}
	if (WriteFile && segn <= USHRT_MAX) {
	    cntblk.BytBlk = 2;
	    cntblk.WordSz = cntblk.WordNo;
	} else {
	    cntblk.BytBlk = 4;
	    cntblk.WordSz = 2 * cntblk.WordNo;
	}

	cntblk.AvrScr = (SHORT) (avr / m);
	size_t	x = 0;
	INT*	wrk = cntblk.blkb = new INT[cntblk.WordNo];
	for (INT w = 0; w < wcp.TabSize; ++w) {
	    x += tcount[w];
	    if (cntblk.Nblk[w]) {
		if (cntblk.wscr[w] > MinScr) {
		    cntblk.blkp[w] = wrk;
		    wrk += cntblk.Nblk[w];
		} else
		    cntblk.blkp[w] = 0;
		cntblk.Nblk[w] = 0;	// reset for next round
	    } else  cntblk.blkp[w] = 0;
	}
	prompt( mfmt, m, c, i, MinScr, cntblk.MaxBlk, cntblk.WordNo, x,
	    100. * m / wcp.TabSize, 100. * i / wcp.TabSize, 100. * c / wcp.TabSize);
}

void MakeBlk::store_blk(bool first, Block* blk)
{
	if (!blk) blk = this;
	if (first) {
	    blk->ch->countBlk(cntblk);
	} else
	    blk->ch->registBlk(chrbuf.segn, cntblk);
	++chrbuf.segn;
	blk->reset();
}

/*****************************************************
	Pre Scan
*****************************************************/

template <typename file_t>
void PreScan::scan_genome(file_t fd, SeqLenStat& sls)
{
	int	c = 0;
	INT	posinentry = 0;

	while ((c = fgetc(fd)) != EOF) {
	    if (isalpha(c)) {
		INT	uc = a2r[toupper(c) - 'A'];		// ignore
		if (uc > Nalpha || (ignoreamb && uc == Nalpha)) continue;
		++posinentry;
	    } else switch (c) {
		case '>':		// FASTA Header
		    if (posinentry > 0) {
			++sls.num;
			sls.total += posinentry;
			if (sls.maxv < posinentry) sls.maxv = posinentry;
			posinentry = 0;
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
	if (posinentry > 0) {
	    ++sls.num;
	    sls.total += posinentry;
	    if (sls.maxv < posinentry) sls.maxv = posinentry;
	}
}

void PreScan::lenStat(int argc, const char** seqdb, SeqLenStat& sls)
{
	for (int ac = 0; ac++ < argc && *seqdb; ++seqdb) {
	    bool	gz = is_gz(*seqdb);
	    if (gz) {
#if USE_ZLIB
		gzFile	gzfd = gzopen(*seqdb, "rb");
		if (gzfd) {
		    scan_genome(gzfd, sls);
		    fclose(gzfd);
		} else	fatal("%s not found!\n", *seqdb);
#else
		fatal(gz_unsupport);
#endif
	    } else {
		FILE*	fd = fopen(*seqdb, "r");
		if (!fd) fatal("%s not found!\n", *seqdb);
		scan_genome(fd, sls);
		fclose(fd);
	    }
	}
}

static void idx2SeqLenStat(const char* fn, SeqLenStat& sls)
{
	FILE*	fidx = fopen(fn, "r");
	if (!fidx) fatal(not_found, fn);
	fseek(fidx, 0L, SEEK_END);
	long	fp = ftell(fidx);
	if (fp == 0 || fp % sizeof(DbsRec))
	    fatal("%s: Index file may be corrupted!\n", fn);
	INT	numidx = fp / sizeof(DbsRec);
	DbsRec*	recidx = new DbsRec[numidx];
	rewind(fidx);
	if (fread(recidx, sizeof(DbsRec), numidx, fidx) != numidx)
	    fatal("%s: Index file may be corrupted!\n", fn);
	sls.num = numidx;
	sls.maxv = 0;
	DbsRec*	widx = recidx;
	DbsRec*	tidx = recidx + numidx;
	for ( ; widx < tidx; ++widx) {
	    if (sls.maxv < widx->seqlen) sls.maxv = widx->seqlen;
	    sls.total += widx->seqlen;
	}
	delete[] recidx;
}

static bool testlut(char* str, int molc, const char** argv)
{
	switch (molc) {
	    case PROTEIN: return (false);
	    case DNA:   strcat(str, LUN_EXT); break;
	    case TRON:  strcat(str, LUP_EXT); break;
	}
	return (newer(str, argv));
}

// Read from amino acid or nucleotide sequence files

template <typename file_t>
void MakeBlk::scan_genome(file_t fd, INT* tc)
{
	int	c = 0;
	bool	mdbs = tc && mkdbs;
	int	posinblk = 0;
	int	rest = 0;
	while ((c = fgetc(fd)) != EOF) {
	    if (isalpha(c)) {
		INT	uc = a2r[toupper(c) - 'A'];		// ignore
		if (uc > Nalpha || (ignoreamb && uc == Nalpha)) continue;
		if (mkblk) {
		    INT* tcc = (posinblk < int(wcp.blklen))? tc: 0;
		    if (istron) c2w6(uc, tcc); else c2w(uc, tcc);
		    if (++posinblk == s_size) {		// block boundary
			store_blk(tc);
			posinblk = prelude + MinOrf;
		    }
		    ++chrbuf.spos;
		}
		if (mdbs) mkdbs->putseq(encode(c));
	    } else switch (c) {
		case '>':		// FASTA Header
		    rest = max(posinblk - int(wcp.blklen), 0);
		    if (mkblk && posinblk) {
			if (MinOrf) {
			    c2w6_pp(MinOrf - rest);
			    if (rest > 0) {
				store_blk(tc);
				c2w6_pp(rest);
			    }
			}
			store_blk(tc);
		     }
		     if (mkblk) {
			if (pchrid) *pchrid++ = chrbuf;
			else	++cntblk.ChrNo;
			posinblk = 0;
		    }
		    if (mdbs) {
			c = mkdbs->write_recrd(fd);
			if (c == '\n' || c == EOF) break;
		    }			// no break
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
	    if (c == EOF) break;
	}
// last block
	if (mkblk && posinblk) {
	    if (MinOrf) {
		rest = max(posinblk - int(wcp.blklen), 0);
		c2w6_pp(MinOrf - rest);
		if (rest > 0) {
		    store_blk(tc);
		    c2w6_pp(rest);
		}
	    }
	    store_blk(tc);
	}
	if (mdbs) {
	    mkdbs->write_recrd(fd, EOF);
	    mkdbs->stamp21();
	    mkdbs->mkidx();
	}
	if (fd) fclose(fd);
}

void MakeBlk::idxblk(int argc, const char** argv)
{
// first phase

const	char**	seqdb = argv;

	for (int ac = 0; ac++ < argc && *seqdb; ++seqdb) {
	    if (mkdbs) mkdbs->wrtgrp(*argv);
	    bool	gz = is_gz(*seqdb);
	    if (gz) {
#if USE_ZLIB
		gzFile	gzfd = gzopen(*seqdb, "rb");
		if (gzfd)
		    scan_genome(gzfd, tcount);
		else	fatal("%s not found!\n", *seqdb);
#else
		fatal(gz_unsupport);
#endif
	    } else {
		FILE*	fd = fopen(*seqdb, "r");
		if (!fd) fatal("%s not found!\n", *seqdb);
		scan_genome(fd, tcount);
	    }
	}
	if (mkdbs) mkdbs->wrtgrp("E_O_F");
	if (!mkblk) return;
	cntblk.glen = chrbuf.spos;
	cntblk.ChrID = pchrid = new CHROMO[cntblk.ChrNo + 1];
	cntblk.blkp = new INT*[wcp.TabSize];
	cntblk.wscr = new short[wcp.TabSize];
	if (istron)	blkscrtab(chrbuf.segn);
	else	blkscrtab(chrbuf.segn, chrbuf.spos / chrbuf.segn);

// second phase

	vclear(&chrbuf);
	++chrbuf.segn;		// skip 0-th block
	seqdb = argv;
	for (int ac = 0; ac++ < argc && *seqdb; ++seqdb) {
	    bool	gz = is_gz(*seqdb);
	    if (gz) {
#if USE_ZLIB
		gzFile	gzfd = gzopen(*seqdb, "rb");
		if (gzfd)
		    scan_genome(gzfd);
		else	fatal("%s not found!\n", *seqdb);
#else
		fatal(gz_unsupport);
#endif
	    } else {
		FILE*	fd = fopen(*seqdb, "r");
		if (!fd) fatal("%s not found!\n", *seqdb);
		scan_genome(fd);
	    }
	}
	*pchrid++ = chrbuf;
}

//	read from SeqServer

void MakeBlk::idxblk(Seq* sd, SeqServer* svr)
{
	int	maxlen = 0;

// first phase

	int	entlen = 0;
	int	pfqnum = 0;
	InSt	ist;
	wcp.blklen = UINT_MAX;
	while ((ist = svr->nextseq(sd, 0)) != IS_END) {
	    if (ist == IS_ERR) continue;
	    if (sd->len > maxlen) maxlen = sd->len;
	    CHAR*	ts = sd->at(sd->len);
	    for (CHAR* pq = sd->at(0); pq < ts; ++pq) {
		CHAR	uc = s2r[*pq];
		if (uc > Nalpha) continue;
		c2w(uc, tcount);
		++cntblk.glen;
	    }
	    ch->countBlk(cntblk);
	    ++cntblk.ChrNo;
	    entlen += strlen((*sd->sname)[0]);
	    if (sd->sigII) pfqnum += sd->sigII->pfqnum;
	    reset();
	}
	wcp.blklen = maxlen;
	cntblk.blkp = new INT*[wcp.TabSize];
	cntblk.wscr = new short[wcp.TabSize];
	size_t	sz = isaa? cntblk.glen: (cntblk.glen / 2 + cntblk.ChrNo);
	wdbf->prepare(entlen + cntblk.ChrNo, cntblk.
		ChrNo, sz + cntblk.ChrNo + 1, pfqnum);
	blkscrtab(cntblk.ChrNo, cntblk.glen / cntblk.ChrNo);

// second phase

	DbsRec*	pr = wdbf->recidx;
	CHAR*	ps = wdbf->dbsseq;
	char*	pe = wdbf->entry;
	int*	pg = wdbf->gsiidx;
	int*	pp = wdbf->gsipool;

	INT	m = 0;
	svr->reset();
	while ((ist = svr->nextseq(sd, 0)) != IS_END) {
	    if (ist == IS_ERR) continue;
	    pr->entptr = pe - wdbf->entry;
	    pr->seqptr = ps - wdbf->dbsseq;
	    INT		n = 0;
	    int		b = 0;
	    CHAR*	ts = sd->at(sd->len);
	    for (CHAR* pq = sd->at(0); pq < ts; ++pq) {
		CHAR    uc = s2r[*pq];
		if (uc > Nalpha) continue;
		c2w(uc);
		++n;
		int	c = *pq - bias;
		if (isaa)	*ps++ = c;
		else if (n & 1)	b = c << 4;
		else 	*ps++ = b + c;
	    }
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
	    ch->registBlk(++m, cntblk);
	    reset();
	}
	if (pg) *pg = pp - wdbf->gsipool;
}

// read from MSA

void MakeBlk::idxblk(Seq* sd)
{
// first phase

	CHAR*	ss = sd->at(0);
	CHAR*	ts = sd->at(sd->len);
	for (int m = 0; m < sd->many; ++m) {
	    for (CHAR* pq = ss + m; pq < ts; pq += sd->many) {
		if (IsGap(*pq)) continue;
		CHAR    uc = s2r[*pq];
		if (uc > Nalpha) continue;
		c2w(uc, tcount);
		++cntblk.glen;
	    }
	    ch->countBlk(cntblk);
	    reset();
	}
	cntblk.ChrNo = sd->many;
	cntblk.wscr = new short[wcp.TabSize];
	cntblk.blkp = new INT*[wcp.TabSize];
	INT	sz = isaa? cntblk.glen: (cntblk.glen / 2 + cntblk.ChrNo);
	INT	si = sd->sigII? sd->sigII->lstnum: 0;
	wdbf->prepare(*sd->sname, cntblk.ChrNo, sz + cntblk.ChrNo, si);
	blkscrtab(sd->many, cntblk.glen / sd->many);

// second phase

	DbsRec*	pr = wdbf->recidx;
	CHAR*	ps = wdbf->dbsseq;
	char*	pe = wdbf->entry;
	for (int m = 0; m < sd->many; ) {
	    INT n = 0;
	    int	b = 0;
	    pr->seqptr = ps - wdbf->dbsseq;
	    pr->entptr = pe - wdbf->entry;
	    for (CHAR* pq = ss + m; pq < ts; pq += sd->many) {
		if (IsGap(*pq)) continue;
		CHAR    uc = s2r[*pq];
		if (uc > Nalpha) continue;
		c2w(uc);
		++n;
		int	c = *pq - bias;
		if (isaa)	*ps++ = c;
		else if (n & 1)	b = c << 4;
		else 	*ps++ = b + c;
	    }
	    if (isaa)	*ps++ = SEQ_DELIM;
	    else if (n & 1) *ps++ = b + SEQ_DELIM;
	    else	*ps++ = ddelim;
	    while(*pe++) ;
	    (pr++)->seqlen = n;
	    ch->registBlk(++m, cntblk);
	    reset();
	}
	if (sd->sigII) {
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

#if M_THREAD

/*****************************************************
	multi-thread version
*****************************************************/

BlkQueue::BlkQueue(int qs) : 
	rp(0), wp(0), remain(0), qsize(qs), blkque(0)
{
	if (qs <= 0) return;
	blkque = new Block*[qsize];
	vclear(blkque, qsize);
	pthread_mutex_init(&mutex, 0);
	pthread_cond_init(&not_full, 0);
	pthread_cond_init(&not_empty, 0);
	pthread_cond_init(&not_busy, 0);
}

void BlkQueue::enqueue(Block* blk)
{
	pthread_mutex_lock(&mutex);
	while (remain == qsize)
	    pthread_cond_wait(&not_full, &mutex);
	blkque[wp] = blk;
	if (++wp == qsize) wp = 0;
	++remain;
//if (blk) fprintf(stderr, "e: %d %4d %d\n", qid, blk->bid, remain);
	pthread_cond_signal(&not_empty);
	pthread_mutex_unlock(&mutex);
}

Block* BlkQueue::dequeue()
{
	pthread_mutex_lock(&mutex);
	while (remain == 0)
	    pthread_cond_wait(&not_empty, &mutex);
	Block*	blk = blkque[rp];
//if (blk) fprintf(stderr, "d: %d %4d %d\n", qid, blk->bid, remain);
	if (++rp == qsize) rp = 0;
	--remain;
	pthread_cond_signal(&not_full);
	pthread_mutex_unlock(&mutex);
	return (blk);
}

Block::Block(Block& src, char* s)
	: WordTab(src), sd(src.sd), 
	  isaa(src.isaa), istron(src.istron), 
	  mkblk(src.mkblk), ch(0), cps(0), 
	  prelude(src.prelude), margin(src.margin), seq(s), 
	  esq(s + wcp.blklen + margin), tsq(esq), bid(0)
{
	if (mkblk) ch = new Chash(int(hfact * wcp.blklen));
}

void* worker(void* arg)
{
	Targ*	targ = (Targ*) arg;

	while (true) {
	    Block*	blk = targ->prev_q->dequeue();
	    if (!blk) 	break;
	    if (blk->bid % thread_num != targ->tid) {
		targ->next_q->enqueue(blk);
		continue;
	    }
const	    char*	ps = blk->begin();
const	    char*	ts = blk->end();
	    int	n = 0;
	    int	i = 0;
	    for ( ; ps < ts; ++ps, ++i, ++n) {
		INT uc = blk->a2r[toupper(*ps) - 'A'];
		if (blk->mkblk) {
		    INT* tcc = (n < int(wcp.blklen))? targ->tc: 0;
		    if (blk->istron) blk->c2w6(uc, tcc);
		    else	blk->c2w(uc, tcc);
		}
	    }
	    targ->next_q->enqueue(blk);
	}
	targ->next_q->enqueue(0);
	return (void*) 0;
}

void MakeBlk::harvest(Block* blk, bool first)
{
	if (mkdbs && first) {
	    const	char*	ts = blk->right();
	    for (const char* ss = blk->begin(); ss < ts; ++ss)
		mkdbs->putseq(encode(*ss));
	}
	int	n = blk->end() - blk->begin();
	int	rest = max(n - int(wcp.blklen), 0);
	bool	tail = blk->end() == blk->right();
	if (mkblk) {
	    if (MinOrf && tail) {
		blk->c2w6_pp(MinOrf - rest);
		if (rest) {
		    store_blk(first, blk);
		    blk->c2w6_pp(rest);
		}
	    }
	    store_blk(first, blk);
	}
	if (tail) blk->restore();
}

template <typename file_t>
void MakeBlk::m_scan_genome(file_t fd, bool first)
{
	int	c = 0;		// read letter
	int	a = 0;		// chromosomal block no
	int	b = 0;		// global block no
	int	n = 0;		// position within block
	int	cps = 0;	// chromosomal position
	char*	ps = 0;
	Block*	blk = 0;
	BlkQueue*&	prev_q = bq[c_qsize - 1];
	BlkQueue*&	next_q = bq[0];
	bool	mdbs = first && mkdbs;

	while ((c = fgetc(fd)) != EOF) {
	    if (isalpha(c)) {
		INT	uc = a2r[toupper(c) - 'A'];
		if (uc > Nalpha || (ignoreamb && uc == Nalpha)) continue;
		*ps++ = c;
		if (++n == s_size) {
		    blk->setbid(b++);
		    next_q->enqueue(blk);
		    if (++a < c_qsize) blk = blks[a];
		    else {
			blk = prev_q->dequeue();
			harvest(blk, first);
		    }
		    memcpy(blk->begin(), ps - margin, margin);
		    ps = blk->begin() + (n = margin);
		    blk->cpos(cps += wcp.blklen);
		}
	    } else switch (c) {
		case '>':		// FASTA Header
		    if (n) {
			blk->end(ps);
			blk->setbid(b++);
			next_q->enqueue(blk);
			++a;
			chrbuf.spos += cps + n;
		    }
		    if (a > c_qsize) a = c_qsize;
		    while (a--) {
			blk = prev_q->dequeue();
			if (blk) harvest(blk, first);
		    }
		    blk = blks[a = 0];
		    ps = blk->begin();
		    blk->cpos(cps = n = 0);
		    if (pchrid)	*pchrid++ = chrbuf;
		    else	++cntblk.ChrNo;
		    if (mdbs) {
			c = mkdbs->write_recrd(fd);
			if (c == '\n' || c == EOF) break;
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
	    if (c == EOF) break;
	}
// last block
	if (n) {
	    blk->end(ps);
	    blk->setbid(b++);
	    next_q->enqueue(blk);
	    ++a;
	    chrbuf.spos += cps + n;
	}
	if (a > c_qsize) a = c_qsize;
	while (a--) {
	    blk = prev_q->dequeue();
	    if (blk) harvest(blk, first);
	}
	if (mdbs) {
	    mkdbs->write_recrd(fd, EOF);
	    mkdbs->stamp21();
	    mkdbs->mkidx();
	}
	fclose(fd);
}

void MakeBlk::m_idxblk(int argc, const char** argv)
{
	c_qsize = thread_num + 1;
	int	qsize = max_queue_num > 1? max_queue_num: c_qsize;
	bq = new BlkQueue*[c_qsize];
	for (int n = 0; n < c_qsize; ++n) {
	    bq[n] = new BlkQueue(qsize);
bq[n]->qid = n;
	}
	blks = new Block*[c_qsize];
	seqbuf = new char[s_size * c_qsize];
	Targ*	targs = new Targ[thread_num];
	pthread_t*	handle = new pthread_t[thread_num];
	INT*	tcbuf = new INT[wcp.TabSize * thread_num];
	vclear(tcbuf, wcp.TabSize * thread_num);
	char*	wseq = seqbuf - s_size;
	for (int b = 0; b < c_qsize; ++b)
	    blks[b] = new Block(*this, wseq += s_size);
	INT*	tc = tcbuf;
	for (int n = 0; n < thread_num; ++n, tc += wcp.TabSize) {
	    targs[n].tid = thread_num - 1 - n;
	    targs[n].tc = tc;
	    targs[n].prev_q = bq[n];
	    targs[n].next_q = bq[n + 1];
	    pthread_create(handle + n, 0, worker, (void*) (targs + n));
	}

	const char**    seqdb = argv;
	for (int ac = 0; ac++ < argc && *seqdb; ++seqdb) {
	    if (mkdbs) mkdbs->wrtgrp(*argv);
	    bool	gz = is_gz(*seqdb);
	    if (gz) {
#if USE_ZLIB
		gzFile  gzfd = gzopen(*seqdb, "rb");
		if (gzfd)
		    m_scan_genome(gzfd, true);
		else    fatal("%s not found!\n", *seqdb);
#else
		fatal(gz_unsupport);
#endif
	    } else {
		FILE*   fd = fopen(*seqdb, "r");
		if (!fd) fatal("%s not found!\n", *seqdb);
		m_scan_genome(fd, true);
	    }
	}
	if (mkdbs) mkdbs->wrtgrp("E_O_F");
	cntblk.glen = chrbuf.spos;
	cntblk.blkp = new INT*[wcp.TabSize];
	cntblk.ChrID = pchrid = new CHROMO[cntblk.ChrNo + 1];
	cntblk.wscr = new short[wcp.TabSize];
	for (int n = 0; n < thread_num; ++n) {
	    for (INT w = 0; w < wcp.TabSize; ++w)
		tcount[w] += targs[n].tc[w];
	    targs[n].tc = 0;
	}
	delete[] tcbuf;
	if (istron)	blkscrtab(chrbuf.segn);
	else    blkscrtab(chrbuf.segn, cntblk.glen / chrbuf.segn);

// second phase

	vclear(&chrbuf);
	++chrbuf.segn;		// skip 0-th block
	seqdb = argv;
	for (int ac = 0; ac++ < argc && *seqdb; ++seqdb) {
	    bool	gz = is_gz(*seqdb);
	    if (gz) {
#if USE_ZLIB
		gzFile  gzfd = gzopen(*seqdb, "rb");
		if (gzfd)
		    m_scan_genome(gzfd, false);
		else    fatal("%s not found!\n", *seqdb);
#else
		fatal(gz_unsupport);
#endif
	    } else {
		FILE*   fd = fopen(*seqdb, "r");
		if (!fd) fatal("%s not found!\n", *seqdb);
		m_scan_genome(fd, false);
	    }
	}
	*pchrid++ = chrbuf;

	bq[0]->enqueue(0);			// mark end of input
	for (int n = 0; n < thread_num; ++n)
	    pthread_join(handle[n], 0);
	delete[] targs;
	delete[] handle;
}

#endif	// M_THREAD

/*****************************************************
	SrchBlk
*****************************************************/

static	int	cmpf(const BLKTYPE* a, const BLKTYPE* b);
static	int	scmpf(const KVpair<INT, int>* a, const KVpair<INT, int>* b);
static	int	gcmpf(const GeneRng* a, const GeneRng* b);

template <typename file_t>
int SrchBlk::read_blk_dt(file_t fd)
{
	if (fread(pbwc->Nblk, sizeof(SHORT), wcp.TabSize, fd) != wcp.TabSize)
	    return (ERROR);
	SHORT** sblkp = (SHORT**) pbwc->blkp;
	SHORT*  sblkb = pbwc->BytBlk == 3?
	    new SHORT[pbwc->WordSz]: (SHORT*) pbwc->blkb;
	if (pbwc->VerNo >= 25) {
	    INT*	blkidx = new INT[wcp.TabSize];
	    if (fread(blkidx, sizeof(INT), wcp.TabSize, fd) != wcp.TabSize) {
		delete[] blkidx;
		return (ERROR);
	    }
	    // Assume pbwc->WordNo < UINT_MAX
	    for (INT w = 0; w < wcp.TabSize; ++w) 
		pbwc->blkp[w] = blkidx[w]? pbwc->blkb + blkidx[w] - 1: 0;
	    delete[] blkidx;
	} else {
	    if (fread(sblkp, sizeof(SHORT*), wcp.TabSize, fd) != wcp.TabSize)
		return (ERROR);
	    long	offset = LONG_MAX;
	    for (INT w = 0; w < wcp.TabSize; ++w) {
		if (pbwc->Nblk[w] && sblkp[w]) {
		    if (pbwc->BytBlk == 3) {
			if (offset == LONG_MAX) offset = sblkp[w] - sblkb;
			sblkp[w] -= offset;
		    } else if (pbwc->VerNo < 23) {
			long    idx = sblkp[w] - sblkb;
			if (offset == LONG_MAX) offset = idx;
			pbwc->blkp[w] = pbwc->blkb + idx - offset;
		    } else {
			if (offset == LONG_MAX) offset = pbwc->blkp[w] - pbwc->blkb;
			pbwc->blkp[w] -= offset;
		    }
		} else {
		    if (pbwc->BytBlk == 3) sblkp[w] = 0;
		    else pbwc->blkp[w] = 0;
		}
	    }
	}
	if (fread(sblkb, sizeof(SHORT), pbwc->WordSz, fd) != pbwc->WordSz ||
	    fread(pbwc->wscr, sizeof(short), wcp.TabSize, fd) != wcp.TabSize)
 		return (ERROR);
	if (pbwc->VerNo <= 25) {
	    INT*	cvt = new INT[pbwc->ConvTS];
	    if (fread(cvt, sizeof(INT), pbwc->ConvTS, fd) != pbwc->ConvTS) 
		return (ERROR);
	    for (INT i = 0; i < pbwc->ConvTS; ++i)
		ConvTab[i] = (CHAR) cvt[i];
	    delete[] cvt;
	} else if (fread(ConvTab, sizeof(CHAR), pbwc->ConvTS, fd) != pbwc->ConvTS)
		return (ERROR);
	if (pbwc->BytBlk == 4) return (OK);
	if (pbwc->BytBlk == 2) {
	    SHORT*	sblk = sblkb + pbwc->WordSz;
	    INT*	iblk = pbwc->blkb + pbwc->WordNo;
	    while (sblk > sblkb) *--iblk = *--sblk;
	    return (OK);
	}
	BYTE4	b4 = {0};
	BYTE2	b2 = {0};
	INT*	iblk = pbwc->blkb;
	for (INT w = 0; w < wcp.TabSize; ++w) {
	    SHORT*	sblk = sblkp[w];
	    if (!sblk) {
		pbwc->blkp[w] = 0;
		continue;
	    }
	    pbwc->blkp[w] = iblk;
	    for (int i = 0; i < pbwc->Nblk[w]; ++i) {
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
static void read_pwc(ContBlk* pwc, file_t fd, const char* fn)
{
	size_t	fpos = ftell(fd);
	if (fread(pwc, sizeof(ContBlk), 1, fd) != 1 || pwc->VerNo < 24) {
	    fseek(fd, fpos, SEEK_SET);
	    ContBlk22	tmp;
	    if (fread(&tmp, sizeof(ContBlk22), 1, fd) == 1) {
		if (tmp.VerNo < 21) {
		    wcp.Nbitpat = 1;
		    wcp.Bitpat2 = 0;
		    swap(wcp.blklen, tmp.ConvTS);
		}
		pwc->ConvTS = tmp.ConvTS;
		pwc->WordNo = tmp.WordNo;
		pwc->WordSz = tmp.WordSz;
		pwc->ChrNo = tmp.ChrNo;
		pwc->glen = tmp.glen;
		pwc->AvrScr = tmp.AvrScr;
		pwc->MaxBlk = tmp.MaxBlk;
		pwc->BytBlk = tmp.BytBlk;
		pwc->VerNo = tmp.VerNo;
	    } else	fatal(emssg, fn);
	}
}

template <typename file_t>
void SrchBlk::ReadBlkInfo(file_t fd, const char* fn)
{
	INT	GivenMaxGene = wcp.MaxGene;
	pbwc = new ContBlk;
	pb2c = new Block2Chr;
	if (fread(&wcp, sizeof(BlkWcPrm), 1, fd) != 1) fatal(emssg, fn);
	read_pwc(pbwc, fd, fn);
	if (pbwc->VerNo < 20) 
	    fatal("%s: Version %d no longer supported !\n", fn, pbwc->VerNo);
	if (GivenMaxGene) wcp.MaxGene = GivenMaxGene;
	if (fread(pb2c, sizeof(Block2Chr), 1, fd) != 1)
	    fatal(emssg, fn);
	if (wcp.Ktuple == wcp.BitPat) wcp.BitPat = bitmask(wcp.BitPat);
	size_t	i = pbwc->ChrNo + 1;
	pbwc->ChrID = (pbwc->VerNo <= BsVersion)? new CHROMO[i]: 0;
	pbwc->Nblk = new SHORT[wcp.TabSize];
	pbwc->blkp = new INT*[wcp.TabSize];
	pbwc->blkb = new INT[pbwc->WordNo];
	pbwc->wscr = new short[wcp.TabSize];
	ConvTab = new CHAR[pbwc->ConvTS];
	if (pbwc->VerNo < 22) {	// for back compatibility
	    CHROMO21*	chr21 = new CHROMO21[i];
	    if (fread(chr21, sizeof(CHROMO21), i, fd) != i)
		fatal(emssg, fn);
	    for (INT j = 0; j < i; ++j) {
		pbwc->ChrID[j].spos = chr21[j].spos;
		pbwc->ChrID[j].segn = chr21[j].segn;
	    }
	    delete[] chr21;
	} else if (pbwc->VerNo < BsVersion) {
	    CHROMO25*	chr25 = new CHROMO25[i];
	    if (fread(chr25, sizeof(CHROMO25), i, fd) != i)
		fatal(emssg, fn);
	    for (INT j = 0; j < i; ++j) {
		pbwc->ChrID[j].spos = chr25[j].spos;
		pbwc->ChrID[j].segn = chr25[j].segn;
	    }
	    delete[] chr25;
	} else if (pbwc->ChrID && (fread(pbwc->ChrID, sizeof(CHROMO), i, fd) != i))
		fatal(emssg, fn);
	if (read_blk_dt(fd) == ERROR) fatal(emssg, fn);
	pbwc->VerNo = BsVersion;
	fclose(fd);
}

template <typename file_t>
static void readBlkInfo(file_t fd, const char* fn, ContBlk* pwc, Block2Chr* pbc)
{
	if (fread(&wcp, sizeof(BlkWcPrm), 1, fd) != 1) fatal(emssg, fn);
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
	char	str[LINE_MAX] = "";
	if (!path) path = DBS_DIR;
	if (is_gz(fn)) {
#if USE_ZLIB
	    gzFile	gzfd = gzopen(fn, "rb");
	    if (!gzfd &&
		!(gzfd = gzopenpbe(path, fn, NGZ_EXT, "rb", -1, str)) &&
		!(gzfd = gzopenpbe(path, fn, PGZ_EXT, "rb", -1, str)) &&
		!(gzfd = gzopenpbe(path, fn, AGZ_EXT, "rb", -1, str))) 
		    fatal("Can't read %s !\n", fn);
	    readBlkInfo(gzfd, fn, &wc, &bc);
#else
	    fatal(gz_unsupport, fn);
#endif
	} else {
	    FILE*	fd = fopen(fn, "rb");
	    if (!fd &&
	  	!(fd = fopenpbe(path, fn, BKN_EXT, "rb", -1, str)) &&
		!(fd = fopenpbe(path, fn, BKP_EXT, "rb", -1, str)) &&
		!(fd = fopenpbe(path, fn, BKA_EXT, "rb", -1, str))) 
		    fatal("Can't read %s !\n", fn);
	    readBlkInfo(fd, fn, &wc, &bc);
	}
	size_t	chrn = wc.ChrNo + 1;
	MaxBlock = wc.VerNo < 21? wcp.MaxGene / wcp.blklen: wcp.Nbitpat;
	const	char*	ms = wc.VerNo < 21? "MaxBlock": "BitpatNo";
	printf("Versioin %d, Elem %u, Tuple %u, TabSize %u, BlockSize %u, MaxGene %u, %s %u, afact %u\n",
	    wc.VerNo, wcp.Nalpha, wcp.Ktuple, wcp.TabSize, wcp.blklen, wcp.MaxGene, ms, MaxBlock, wcp.afact);
	printf("GenomeSize %lu, ChrNo %lu, BlockNo %lu, WordNo %lu, WordSz %lu B, AvrScr %u\n",
	    (LONG) wc.glen, (LONG) wc.ChrNo, (LONG) (wc.ChrID? wc.ChrID[wc.ChrNo].segn: chrn), 
	    (LONG) wc.WordNo, (LONG) wc.WordSz, wc.AvrScr);
	delete[] wc.ChrID;

// read lun/p information
	if (!*str) strcpy(str, fn);
	char*	dot = strrchr(str, '.');
        if (is_gz(str)) {
	    *dot = 0;
	    dot = strrchr(str, '.');
	}
	if (!strcmp(dot, BKN_EXT)) strcpy(dot, LUN_EXT);
	else if (!strcmp(dot, BKP_EXT)) strcpy(dot, LUP_EXT);
	else exit(0);
	exit (0);
}

void SrchBlk::init2(const Seq* sd)
{
	vclear(bh2->rscr, nseg);
	bh2->prqueue_b->reset();
const	CHAR*	ss = sd->at(sd->left);
const	CHAR*	ts = sd->at(sd->right - (bpp[0]->width + wcp.Nshift) + 1);
	INT	q = (ts - ss) % wcp.Nshift;
	MinGeneLen = 0;
	for (INT p = 0; p < wcp.Nshift; ++p) {
	    bh2->as[0][p] = ss++;
	    bh2->as[1][q] = ts++;
	    if (++q == wcp.Nshift) q = 0;
	}
}

void SrchBlk::init4(const Seq* sd)
{
	MinGeneLen = bbt * (sd->right - sd->left) * 3 / 4;
	vclear(bh4->sigw[0], Ncand * 2);
	vclear(bh4->sigw[1], Ncand * 2);
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
const	CHAR*	ss = sd->at(sd->left);
const	CHAR*	ts = sd->at(sd->right - (bpp[0]->width + wcp.Nshift));	// last cycle
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
	for (INT chrn = 0; chrn < pbwc->ChrNo; ++chrn) {
	    int	len = chrsize(chrn);
	    for (; len > (int) wcp.blklen; len -= wcp.blklen)
		SegLen[++blk] = wcp.blklen;
	    if (len > 0) {
		if (++blk == nseg) fatal("SrchBlk::setSegLen block # mismatch (%d)\n", blk);
		SegLen[blk] = len;
	    }
	}
	if (blk != (nseg - 1))
	    prompt("Warn: SegLen may be incorrect: %d != %d\n", blk, nseg);
}

int SrchBlk::findChrNo(const INT& blk)
{
	int	lw = (int) (pb2c->BClw + pb2c->BCce * (blk - 1)) - 1;
	int	up = (int) (pb2c->BCup + pb2c->BCce * (blk - 1)) + 1;

	if (lw < 0) lw = 0;
	if (up > (int) pbwc->ChrNo) up = pbwc->ChrNo;
	if (chrblk(lw) > blk) lw = 0;
	if (chrblk(up) < blk) up = pbwc->ChrNo;
	while (up - lw > 1) {
	    int	md = (lw + up) / 2;
	    if (chrblk(md) > blk) up = md;
	    else if (chrblk(md + 1) > blk) return (md);
	    else 	lw = md;
	}
	if (chrblk(up) > blk) return (lw);
	else	return (up);
}

Seq* SrchBlk::setgnmrng(const BPAIR* wrkbp)
{
	INT	x = (wrkbp->zl? wrkbp->lb - wrkbp->zl: 0) * wcp.blklen;
	INT	y = ((wrkbp->zl? wrkbp->rb - wrkbp->zl: 0) + 1) * wcp.blklen;
	char	str[LINE_MAX];
	sprintf(str, "Dbs%d %d %d", wrkbp->chr, x + 1, y);
	Seq*	sd = (*curgr)->getdbseq(dbf, str, wrkbp->chr);
	if (wrkbp->rvs) sd->comrev();
	return (sd);
}

int SrchBlk::MinQuery() const {
	return (2 * wcp.Ktuple + wcp.Nshift);
}

int SrchBlk::MaxGene() const {
	return wcp.MaxGene;
}

void SrchBlk::setaaseq(Seq* sd, int chn)
{
	sd->getdbseq(dbf, "", chn);
}

static int cmpf(const BLKTYPE* a, const BLKTYPE* b)
{
	return (*a - *b);
}

static int scmpf(const KVpair<INT, int>* a,const  KVpair<INT, int>* b)
{
	return (b->val - a->val);
}

static int gcmpf(const GeneRng* a, const GeneRng* b)
{
	VTYPE	d = b->scr - a->scr;

	if (d > 0) return 1;
	if (d < 0) return -1;
	return 0;
}

Randbs::Randbs(double avr, bool gdb)
{
	if (gdb) {
	    if (RbsFact == FQUERY) RbsFact = RbsFactLog;
	    trans_form = log;
	} else {
	    if (RbsFact == FQUERY) RbsFact = algmode.crs? RbsFactSqrX: RbsFactSqr;
	    trans_form = sqrt;
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
	int	play = IntronPrm.minl;
	if (sd->sigII) {
	    int	tge = sd->sigII->to_gene_end(apos, rend);
	    return (tge + play);
	}
	if (agap < mingap && algmode.crs) return (2 * agap + play);
	return (IntronPrm.maxl + play);
}

bool SrchBlk::grngoverlap(const GeneRng* a, const GeneRng* b)
{
	int	al = 0;
	int	bl = 0;
	int	ol = 0;
const 	JUXT*	aj = a->jxt;
const 	JUXT*	bj = b->jxt;
const 	JUXT*	at = aj + a->num;
const 	JUXT*	bt = bj + b->num;
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

SrchBlk::SrchBlk(Seq* sqs[], const char* fn, bool gdb) :
	gnmdb(gdb), dbf(dbs_dt[0])
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
	initialize(sqs, fn);
}

SrchBlk::SrchBlk(Seq* sqs[], MakeBlk* mb, bool gdb) :
	gnmdb(gdb), dbf(mb->wdbf)
{
	pbwc = new ContBlk;
	pb2c = new Block2Chr;
	*pbwc = mb->cntblk;
	*pb2c = mb->b2c;
	mb->wdbf = 0;
	if (dbf->curdb->defmolc != PROTEIN) dbf->curdb->defmolc = DNA;
	ConvTab = mb->iConvTab;
	mb->iConvTab = 0;
	mb->cntblk.clean();
	initialize(sqs);
}

SrchBlk::SrchBlk(SrchBlk* sbk, DbsDt* df)
{
	*this = *sbk;
	master = sbk;
	if (sbk->bh2) bh2 = new Bhit2(sbk->nseg);
	if (sbk->bh4) bh4 = new Bhit4(sbk->nseg);
	reset(df);
}

SrchBlk::~SrchBlk()
{
	delete bh2; delete bh4;
	if (dbf != dbs_dt[0] && dbf != dbs_dt[1]) delete dbf;
	if (master) {		// secondary
#if TESTRAN
	    master->add(tstrn);
	    delete tstrn;
#endif
	} else {		// primary
#if TESTRAN
	    tstrn->out(pbwc->AvrScr);
	    delete tstrn;
#endif
	    delete pbwc; delete pb2c;
	    delete pwd; delete rdbt;
	    delete[] ConvTab;
	    delete[] SegLen;
	    for (int k = 0; k < kk; ++k) delete bpp[k];
	    delete[] bpp;
	}
}

void SrchBlk::initialize(Seq* sqs[], const char* fn)
{
	if (!pbwc->blkp)
	    fatal("%s: Make and specify database !\n", fn);
	pwd = new PwdB((const Seq**) sqs);
	DRNA = !pwd->DvsP;
	NoWorkSeq = (DRNA && algmode.lsg)? 2: 0;
	if (DRNA ^ (wcp.Nalpha == 4))
	    fatal("%s: Block table is incompatible with query !\n", fn);
	if (!pbwc->AvrScr)
	    fatal("%s: Block table may be destroyed !\n", fn);
	kk = wcp.Nbitpat / 2 + 1;
	bpp = new Bitpat*[kk];
	if (wcp.Nbitpat == 1) {
	    bpp[0] = new Bitpat(wcp.BitPat);
	} else {
	    bpp[0] = new Bitpat(bitmask(wcp.Ktuple));
	    bpp[1] = new Bitpat(wcp.BitPat);
	}
	if (wcp.Nbitpat > 3) {
	    bpp[2] = new Bitpat(wcp.Bitpat2);
	}
	min_agap = bpp[0]->width;
	rdbt = new Randbs((double) pbwc->AvrScr * bpp[0]->weight / wcp.Nshift, gnmdb);
#if TESTRAN
	tstrn = 0;
	MaxMmc = 0;
#endif
	maxmmc = (MaxMmc == 0 || (int) MaxMmc > (INT_MAX / bpp[0]->weight)
		 || algmode.lcl & 16)?
	    INT_MAX: bpp[0]->weight * MaxMmc / wcp.Nshift;
	vthr = (VTYPE) (alprm.scale * 2 * alprm.thr);
	ptpl = rdbt->base() / pbwc->AvrScr;
	DeltaPhase2 = (int) (ptpl * pbwc->AvrScr);
	MaxBlock = wcp.MaxGene / wcp.blklen;
	if (ExtBlock == 0) 
	    ExtBlock = max_intron_len(ild_up_quantile, fn) / wcp.blklen + 1;
	ExtBlockL = MaxBlock / 2 + 1;
	nseg = chrblk(pbwc->ChrNo);
	bbt = pwd->DvsP == 1? 3: 1;
	set_shortquery(8 * wcp.Ktuple);
	Ncand = OutPrm.MaxOut + NCAND2PHS;
	if (gnmdb)	bh4 = new Bhit4(nseg);
	else {
	    bh2 = new Bhit2(nseg);
	    setSegLen();
	}
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

void set_max_extend_gene_rng(int n, bool forced)
{
	if (forced || gene_rng_max_extend < 0)
	    gene_rng_max_extend = n;
}

int get_max_extend_gene_rng()
{
	return (gene_rng_max_extend);
}

bool extend_gene_rng(Seq* sqs[], const PwdB* pwd, DbsDt* dbf)
{
	if (gene_rng_max_extend <= 0) return (false);
	Seq*&   a = sqs[0];
	Seq*&   b = sqs[1];
	RANGE   prv = {a->right, a->left};
	RANGE   grng = {b->left, b->right};
	int	bbt = a->isprotein()? 3: 1;
	WLPRM*  wlprm = setwlprm(0);
	bool    rvs = b->inex.sens & 1;
	bool	extended = false;

	for (int ntry = 0; ntry < gene_rng_max_extend; ++ntry) {
	    JUXT*	lend = b->jxt;
	    JUXT*	wjxt = b->jxt + b->CdsNo - 1;
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

	    char	str[LINE_MAX];
	    if (rvs)
		sprintf(str, "$%s %d %d <", (*b->sname)[0],
		    b->SiteNo(grng.right), b->SiteNo(grng.left));
	    else
		sprintf(str, "$%s %d %d", (*b->sname)[0],    
		    b->SiteNo(grng.left), b->SiteNo(grng.right));
	    Seq*    tmp = new Seq(1, abs(grng.right - grng.left));
	    swap(b, tmp);
	    b->getdbseq(dbf, str);
	    if (pwd->DvsP == 1) b->nuc2tron();
	    Wilip	wl((const Seq**) sqs, pwd, -1);
	    WLUNIT* wlu = wl.begin();
	    if (wlu) {
		b->jxt = new JUXT[wlu->num + 1];
		vcopy(b->jxt, wlu->jxt, wlu->num + 1);
		b->CdsNo = wlu->num;
		delete tmp;
	    } else {
		swap(b, tmp);
		delete tmp;
		break;
	    }
	}
	return (extended);
}

int SrchBlk::FindHsp(BPAIR* wrkbp)
{
const	INT	intr = (*gener)->inex.intr;
const	int	sr = chrsize(wrkbp->chr);
const	bool	rvs = wrkbp->rvs;
const	float	log_2[4] = {0.f, 1.0f, 1.585f, 2.0f};

	wrkbp->jscr = 0;
	Seq*	cursd = 0;
	RANGE	qrng, grng;
	query->saverange(&qrng);
	INT	retry_no = 0;
	RANGE	prv = {query->right, query->left};
	Wilip*	wilip = 0;
	BPAIR	orgbp = *wrkbp;
const	int	qlen = query->right - query->left;

retry:
	cursd = setgnmrng(wrkbp);
	if (!cursd) return (2);
	cursd->inex.intr = intr;
	if (query->isprotein()) cursd->nuc2tron();
	swap(*curgr, seqs[1]);		// temporally save
	delete wilip;
	wilip = new Wilip((const Seq**) seqs, pwd, -1);
	WLUNIT*	wlu = wilip->begin();
	swap(*curgr, seqs[1]);		// restore
	if (!wlu) {delete wilip; return (2);}
	cursd->saverange(&grng);
	int	n = min((int) OutPrm.MaxOut2, wilip->size());
	JUXT	lend = {query->right, 0};
	JUXT	rend = {query->left, 0};
	int	nbetter = 0;
	int	multi = 0;
	VTYPE	lcritjscr = critjscr;
	for ( ; n--; ++wlu) {
	    if (wlu->num > 1) ++multi;
	    if (wlu->tlen > qlen) {
		float	x = wlu->scr;
		wlu->scr = int(x * qlen / wlu->tlen);
	    }
	    wlu->scr += IntronPrm.sip * log_2[std::min(wlu->num - 1, 3)] ;
	    if (wlu->scr >= lcritjscr) {
		++nbetter;
		lcritjscr = wlu->scr - vthr;
		JUXT*	wjxt = wlu->jxt;
	 	if (lend.jx > wjxt->jx) lend = *wjxt;
	 	wjxt += wlu->num - 1;
	 	int	l = wjxt->jx + wjxt->jlen;
	 	if (rend.jx < l) {
		    rend.jx = l;
		    rend.jy = wjxt->jy + bbt * wjxt->jlen;
		}
	    } else {
		for (WLUNIT* wwl = wilip->begin(); wwl < wlu; ++wwl) {
		    if (wwl->llmt <= wlu->ulmt && wwl->llmt > wlu->llmt)
			wwl->llmt = wlu->llmt;
		    if (wwl->ulmt >= wlu->llmt && wwl->ulmt < wlu->ulmt)
			wwl->ulmt = wlu->ulmt;
		}
	    }
	}
	if (!nbetter) {delete wilip; return (multi? 2: 0);}

	if (lend.jx && lend.jx < prv.left) {
	    prv.left = lend.jx;
	    if (rvs) {
		if (lend.jx <= min_agap)
		    wrkbp->rb = wrkbp->db;
		else {
		    wrkbp->rb = min(wrkbp->rb + ExtBlockL, wrkbp->zr);
		    for ( ; wrkbp->rb > wrkbp->db; --wrkbp->rb)
			if (bh4->bscr[2][wrkbp->rb]) break;
		}
		wrkbp->db = min(wrkbp->rb + ExtBlock, wrkbp->zr);
	    } else {
		if (lend.jx <= min_agap)
		    wrkbp->lb = wrkbp->ub;
		else {
		    INT	x = wrkbp->lb > ExtBlockL? wrkbp->lb - ExtBlockL: 0;
		    wrkbp->lb = max(x, wrkbp->zl);
		    for ( ; wrkbp->lb < wrkbp->ub; ++wrkbp->lb)
			if (bh4->bscr[0][wrkbp->lb]) break;
		}
		INT	x = wrkbp->lb > ExtBlock? wrkbp->lb - ExtBlock: 0;
		wrkbp->ub = max(x, wrkbp->zl);
	    }
	}
	if (query->right > rend.jx && rend.jx > prv.right) {
	    prv.right = rend.jx;
	    int	dlt = query->right - rend.jx;
	    if (rvs) {
		if (dlt == 1 && wrkbp->lb > wrkbp->ub)
		    --wrkbp->lb;
		else if (dlt < min_agap)
		    wrkbp->lb = wrkbp->ub;
		else {
		    INT	x = wrkbp->lb > ExtBlockL? wrkbp->lb - ExtBlockL: 0;
		    wrkbp->lb = max(x, wrkbp->zl);
		    for ( ; wrkbp->lb < wrkbp->ub; ++wrkbp->lb)
			if (bh4->bscr[3][wrkbp->lb]) break;
		}
		INT	x = wrkbp->lb > ExtBlock? wrkbp->lb - ExtBlock: 0;
		wrkbp->ub = max(x, wrkbp->zl);
	    } else {
		if (dlt == 1 && wrkbp->rb < wrkbp->db)
		    ++wrkbp->rb;
		else if (dlt < min_agap)
		    wrkbp->rb = wrkbp->db;
		else {
		    wrkbp->rb = min(wrkbp->rb + ExtBlockL, wrkbp->zr);
		    for ( ; wrkbp->rb > wrkbp->db; --wrkbp->rb)
			if (bh4->bscr[1][wrkbp->rb]) break;
		}
		wrkbp->db = min(wrkbp->rb + ExtBlock, wrkbp->zr);
	    }
	}
	if (retry_no++ < NoRetry && 
		(wrkbp->lb != orgbp.lb || wrkbp->rb != orgbp.rb)) {
	    orgbp = *wrkbp;
	    goto retry;
	}

	int	lbias = orgbp.lb - wrkbp->lb;
	int	ubias = wrkbp->rb - orgbp.rb;
	if (lbias || ubias) {
	    cursd = setgnmrng(wrkbp);
	    if (!cursd) {delete wilip; return (0);}
	    if (rvs) swap(lbias, ubias);
	    cursd->inex.intr = intr;
	    if (query->isprotein()) cursd->nuc2tron();
	    if (lbias) {
		int	partial_block_len = (rvs && wrkbp->rb == wrkbp->zr)?	 // chromosomal end
		    (wcp.blklen - cursd->len % wcp.blklen): 0;
		lbias = lbias * wcp.blklen - partial_block_len;
	    }
	    wilip->shift_y(lbias, cursd->len);
	}
	cursd->saverange(&grng);
	wilip->sort_on_scr();	// re-sort on score
	wlu = wilip->begin();
	wrkbp->jscr = wlu->scr;
	Seq**	wrkgr;
	for ( ; wlu->num; ++wlu) {
	    if (wlu->scr < lcritjscr) break;
	    JUXT*	wjxt = wlu->jxt;
	    if (wrkbp - bh4->bpair >= (int) OutPrm.MaxOut && wlu->scr < critjscr)
		break;
	    if (wlu > wilip->begin()) {
		if (cursd != *curgr) cursd->aliaseq(*curgr);
		(*curgr)->restrange(&grng);
	    }
	    (*curgr)->jscr = (VTYPE) wlu->scr;
	    (*curgr)->CdsNo = 0;
	    (*curgr)->left = wlu->llmt;
	    if (wlu->ulmt < (*curgr)->right) (*curgr)->right = wlu->ulmt;
	    int	cl = (*curgr)->SiteNo(wjxt->jy);
	    wjxt += wlu->num - 1;
	    int	cr = (*curgr)->SiteNo(wjxt->jy + wjxt->jlen);
	    if (rvs) swap(cl, cr);
	    int	tl = 0;
	    int	tr = sr;
	    for (wrkgr = curgr; --wrkgr >= gener; ) {	// find nearest
		if (((*wrkgr)->did == (*curgr)->did) &&
		    (*wrkgr)->inex.sens == (*curgr)->inex.sens) {
		    wjxt = (*wrkgr)->jxt;
		    int	wl = (*wrkgr)->SiteNo(wjxt->jy);
		    wjxt += (*wrkgr)->CdsNo - 1;
		    int	wr = (*wrkgr)->SiteNo(wjxt->jy + wjxt->jlen);
		    if (rvs) swap(wl, wr);
		    if (cr > wl && cl < wr) break;	// overlap
		    if (cr < wl && wl < tr) tr = wl;
		    if (wr < cl && wr > tl) tl = wr;
		}
	    }
	    if (wrkgr >= gener) continue;		// overlap
	    for (wrkgr = curgr; --wrkgr >= gener; ) {	// sort on score
		if (wrkgr[1]->jscr > (*wrkgr)->jscr)
		    swap(wrkgr[0], wrkgr[1]);
		else break;
	    }
	    if (++wrkgr >= lstgr) break;	// no more locus
	    for ( ; curgr > gener; --curgr) {	// prune low-jscr loci
		if ((*curgr)->jscr >= lcritjscr) break;
		(*curgr)->refresh(1);
	    }
	    if (curgr - gener >= OutPrm.MaxOut - 1) {
		critjscr = gener[OutPrm.MaxOut - 1]->jscr - vthr;
		if (critjscr < 0) critjscr = 0;
	    }
	    if (curgr < lstgr) ++curgr;
	    delete[] (*wrkgr)->jxt;
	    (*wrkgr)->jxt = new JUXT[wlu->num + 1];
	    vcopy((*wrkgr)->jxt, wlu->jxt, wlu->num + 1);
	    (*wrkgr)->CdsNo = wlu->num;
	    if (curgr == lstgr) (*curgr)->refresh(1);
	}
	delete wilip;
	return (1);
}

SHORT SrchBlk::extract_to_work(int d)
{
	SHORT	e = d + 1;
	SHORT	f = (SHORT) (d >> 1);
	if (!bh4->sign[d] && !bh4->sign[e]) return (0);
	SHORT	j = 0;
	BLKTYPE*&	sw = bh4->sigw[f];
	while (j < bh4->sign[d]) {
	    BlkScr&	bs = (*bh4->prqueue_b[d])[j];
	    sw[j++] = bs.key << 1;
	}
	for (SHORT k = 0; k < bh4->sign[e]; ++k) {
	    BlkScr&	bs = (*bh4->prqueue_b[e])[k];
	    sw[j++] = (bs.key << 1) + 1;
	}
	if (j == 1) {
	    sw[0] &= ~1;
	    sw[1] = INT_MAX;
	    return (1);
	}
	qsort((UPTR) sw, j, sizeof(BLKTYPE), (CMPF) cmpf);
	BLKTYPE	p = sw[0] >> 1;
	if (d && j > 1) {
	    for (SHORT i = 1; i < j; ++i) {
		BLKTYPE	q = sw[i] >> 1;
		if (p == q) swap(sw[i - 1], sw[i]);
		else	p = q;
	    }
	}
	p = sw[0] >> 1;
	bool	pr = (sw[0] & 1) ^ f;
	int	cp = findChrNo(p);
	sw[0] &= ~1;		// reset lsb
	SHORT	k = 0;
	SHORT	c = 0;
	for (SHORT i = 1; i < j; ++i) {
	    BLKTYPE	q = sw[i] >> 1;
	    bool	qr = (sw[i] & 1) ^ f;
	    int	cq = findChrNo(q);
	    sw[i] &= ~1;	// reset lsb
	    int	s = q - p;
	    if (cp == cq && (s < 2 ||	// contiguous
		(!pr && qr && s <= int(MaxBlock)) ||
		(pr == qr && s <= int(ExtBlock)))) {
		if (!c++) sw[k++] = sw[i - 1];
	    } else {	// separate
		sw[k++] = sw[i - 1] | (c? 1: 0);
		c = 0;
	    }
	    p = q;
	    pr = qr;
	    cp = cq;
	}
	sw[k++] = sw[j - 1] + (c? 1: 0);
	sw[k] = INT_MAX;
	return (k);
}

int SrchBlk::TestOutput(int force)
{
	int	ReportAln = algmode.nsa != MAP1_FORM && algmode.nsa != MAP2_FORM;
	BPAIR*	wrkbp = bh4->bpair;
	BPAIR*	curbp = bh4->bpair;
	BPAIR*	lstbp = bh4->bpair + Ncand;
	Seq**	wrkgr;
	SHORT	sigm[2];
static	const	char	ofmt[] = 
	"%-7s %c %5d %5d %-7s %5d %5d %6.2f %6.2f %3ld %3ld %2d %2d %2d %2d %d %d\n";
static	const	char	ofmt2[] = 
	"%-7s %c %5d %5d %-7s %5d %5d %6.2f %6.2f %3d %3d %2d %2d %2d %2d\n";

	curgr = gener;
	curbp->bscr = 0;
	for (wrkgr = gener; wrkgr < lstgr; ++wrkgr) (*wrkgr)->len = 0;
	for (SHORT f = 0; f < 2; ++f)		// sort on position
	    sigm[f] = extract_to_work(f << 1);
	for (SHORT f = 0; f < 2; ++f) {
	    SHORT	d = f << 1;
	    SHORT	e = d + 1;
	    BLKTYPE	pu = 0;
	    for (SHORT i = 0; i < sigm[f]; ++i) {
		BLKTYPE	p = bh4->sigw[f][i] >> 1;
		BLKTYPE	q = bh4->sigw[f][i + 1 < sigm[f]? i + 1: i];
		bool	ispair = q & 1;
		q >>= 1;
		if (ispair) ++i;	// block pair
		else	q = p;	// sigleton
		BLKTYPE	qd = bh4->sigw[f][i + 1] >> 1;
		curbp->jscr = 0;
		curbp->rvs = f;
		int	c1 = curbp->chr = findChrNo(q);
		curbp->zl = chrblk(c1);
		curbp->zr = chrblk(c1 + 1) - 1;
		curbp->lb = p;
		curbp->rb = q;
		curbp->bscr = 0;
		for (BLKTYPE r = curbp->lb; r <= curbp->rb; ++r)
		    curbp->bscr += bh4->bscr[d][r] + bh4->bscr[e][r];
		BLKTYPE	exb = ExtBlock;
		BLKTYPE r = curbp->lb;
		BLKTYPE	z = max3((r > exb)? r - exb: 0, curbp->zl, pu);
		while (r && --r >= z && (bh4->bscr[d][r] + bh4->bscr[e][r])) {
		    curbp->lb = r;
		    curbp->bscr += bh4->bscr[d][r] + bh4->bscr[e][r];
		}
		curbp->ub = max((r > exb)? r - exb: 0, curbp->zl);

		r = curbp->rb;
		z = min3(r + exb, curbp->zr, qd);
		while (++r < z && (bh4->bscr[d][r] + bh4->bscr[e][r])) {
		    curbp->rb = r;
		    curbp->bscr += bh4->bscr[d][r] + bh4->bscr[e][r];
		}
		curbp->db = min(r + exb, curbp->zr);
		pu = curbp->rb + 1;
		for (wrkbp = curbp; --wrkbp >= bh4->bpair; ) {
		    if (wrkbp[1].bscr > wrkbp->bscr) {
			swap(wrkbp[0], wrkbp[1]);
		    } else break;	// sort in decending order
		}
		if (curbp < lstbp) ++curbp;
	    }
	}

	if (curbp == bh4->bpair) return (force? ERROR: 0);
	++force;
	lstbp = curbp;
	int	phase1 = 0;	// look for best block pair
	INT	nfail = OutPrm.MaxOut2 + 2;
	for (wrkbp = curbp = bh4->bpair; wrkbp < lstbp && nfail; ++wrkbp) {
	    if (wrkbp->bscr == 0) continue;
	    SHORT	d = wrkbp->rvs << 1;
	    SHORT	e = d + 1;
	    if (force != 2 && (wrkbp->bscr < 
		rdbt->randbs(bh4->mmct[d] + bh4->mmct[e]) + rdbt->Phase1T))
		continue;
	    switch (FindHsp(wrkbp)) {
		case 1: ++phase1; break;
		case 2: --nfail; break;
		default: break;
	    }
	}
	if (phase1) {			// phase1 passed
	    rectiseq(gener, curgr);
	    if (ReportAln) return (curgr - gener);
	} else if (force < 2) return (0);	// next reccurence
	else if (ReportAln && query->isprotein() && algmode.slv) {	// examine all positive blocks
	    curgr = gener;
	    curbp = bh4->bpair + 1;
	    curbp->bscr = 0;
	    for (SHORT d = 0; d < 4; d += 2) {
		bh4->bpair->rvs = d > 1;
		bh4->bpair->bscr = 0;
		bh4->bpair->chr = 0;
		SHORT	s = pbwc->AvrScr;
		int	c2 = 0;
		SHORT	e = d + 1;
		INT	x, y, z;
		for (y = z = x = c2; ++x < nseg; ) {
		    if (bh4->bscr[d][x] || bh4->bscr[e][x]) {
			int	c1 = findChrNo(x);
			if (((!y && z) || c1 != c2) && bh4->bpair->bscr >= s) {
			    FindHsp(bh4->bpair);
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
			    FindHsp(bh4->bpair);
			    if (bh4->bpair->bscr > curbp->bscr) *curbp = *bh4->bpair;
			}
			y = z = bh4->bpair->bscr = 0;
		    }
		}
	    }
	    if (curgr > gener) return (curgr - gener);
	}
	SHORT	d = curbp->rvs << 1;				// map only or phase1 failed
	SHORT	e = d + 1;
	BLKTYPE	p = d? curbp->rb: curbp->lb;
	BLKTYPE	q = d? curbp->lb: curbp->rb;
	char*	cid = dbf->entname(curbp->chr);
	if (algmode.nsa == MAP1_FORM) {
	    printf(ofmt, query->sqname(), dir[d / 2], query->left, query->right, cid, p, q,
		(double) bh4->bscr[d][p] / TFACTOR,
		(double) bh4->bscr[e][q] / TFACTOR,
		bh4->as[d][bh4->maxs[d]] - query->at(0),
		bh4->as[e][bh4->maxs[e]] - query->at(0), 
		bh4->mmct[d], bh4->mmct[e], bh4->nhit[d], bh4->nhit[e], bh4->testword[d], bh4->testword[e]);
	} else if (algmode.nsa == MAP2_FORM) {
	    int	x = p - chrblk(curbp->chr);
	    int	y = q - chrblk(curbp->chr);
	    if (x > 0) x = (x - 1) * wcp.blklen;
	    y *= wcp.blklen;
	    printf("%s\t%s\t%7d\t%7d %c\t%7d\t%7d\n",
		query->sqname(), cid, x + 1, y, dir[d / 2], y - x, query->len);
	} else { 
	    prompt(ofmt2, query->sqname(), dir[d / 2], query->left, query->right,
		cid, p, q,
		(double) bh4->bscr[d][p] / TFACTOR,
		(double) bh4->bscr[e][q] / TFACTOR,
		bh4->as[d][bh4->maxs[d]] - query->at(0),
		bh4->as[e][bh4->maxs[e]] - query->at(0),
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
	as[0] = new const CHAR*[4 * wcp.Nshift];
	sigw[0] = new BLKTYPE[Ncand * 2 + 1];
	sigw[1] = new BLKTYPE[Ncand * 2 + 1];
	for (int d = 1; d < 4; ++d) {
	    bscr[d] = bscr[d-1] + nseg;
	    ascr[d] = ascr[d-1] + nseg;
	    as[d] = as[d-1] + wcp.Nshift;
	}
	bpair = new BPAIR[Ncand + 1];
	INT	na = Nascr + 1;
	INT	nc = Ncand + 1;
	blkscr = new BlkScr[4 * (nc + na)];
	vclear(blkscr, 4 * (nc + na));
	BlkScr*	bs = blkscr;
	for (int d = 0; d < 4; ++d, bs += nc) {
	    prqueue_a[d] = new PrQueue_wh<BlkScr>(bs, Nascr, 0, false, true);
	    prqueue_b[d] = new PrQueue_wh<BlkScr>(bs += na, Ncand, 0, false, true);
	}
}

Bhit4::~Bhit4()
{
	delete[] bscr[0];
	delete[] ascr[0];
	delete[] sigw[0];
	delete[] sigw[1];
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

Qwords::Qwords(int k, int nc, CHAR* ct, ContBlk* pwc, Bitpat** bp, Seq* a) :
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

int Qwords::querywords(const CHAR* ss)
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

int Qwords::querywords(const CHAR* ss, int d, bool rvs)
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

int SrchBlk::findblock(Seq** sqs)
{
const	int	qlen = query->right - query->left;
	if (qlen - (wcp.Nshift + bpp[0]->width) < 1) return (ERROR);
	init4(query);
	int	nohit = 0;
	int	sigpr = 0;
	int	c = qlen / (wcp.Nshift + wcp.Nshift) - 1;
const	bool	is_shorquery = qlen < shortquery;
const	CHAR*	at = is_shorquery? query->at(query->right): 0;
const	CHAR*	ab = is_shorquery? query->at(query->left): 0;
	bool	meet[2] = {false, false};
	Dhash<INT, int>	hh(2 * pbwc->MaxBlk, 0);
	INT	nmmc = 0;
	int	notry = 0;
	int	maxbscr[4];
	vclear(maxbscr, 4);
	Qwords	qwd(kk, DRNA, ConvTab, pbwc, bpp, query);
	INT	totalsign = 0;
	critjscr = 0;
	while (!(meet[0] || meet[1])) {
	  totalsign = 0;
	  for (SHORT d = 0; d < 4; ++d) {
	    if (meet[d / 2]) continue;
	    SHORT	prty = d % 2;		// direct:reverse
	    SHORT	e = prty? d - 1: d + 1;	// tally
	    bool	rvs = d >= 2;
	    int*&	rscr = bh4->bscr[d];
const	    CHAR*	ms = prty? ab: at;
	    SHORT	maxp = 0;
	    for (INT s = 0; s < wcp.Nshift; s++) {
const		CHAR**	ws = bh4->as[d] + s;
		if (!is_shorquery) ms = bh4->as[e][s];
		int	cscr = 0, q = 0, p = 0;
		hh.clear();
		do {				// continueous hits
const		    CHAR*	ss = *ws;
		    if (prty)	*ws -= wcp.Nshift;	// right side
		    else	*ws += wcp.Nshift;	// left side
		    if (prty ^ (ss >= ms)) {		// has scanned
			meet[d / 2] = true;
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
		if  ((c = TestOutput(0))) return (c);		// found
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
	if (c != ERROR) c = TestOutput(1);
	return (c);
}

INT SrchBlk::bestref(Seq* sqs[], KVpair<INT, int>* sh, int n)
{
	Seq*&	a = sqs[1];
	Wilip**	wl = new Wilip*[n];
	GeneRng	gr;
	Mfile	mfd(sizeof(GeneRng));

	swap(sqs[0], sqs[1]);
	for (int i = 0; i < n; ++i, ++sh) {
	    if (!sh->val) continue;
	    setaaseq(a, sh->key);
	    wl[i] = new Wilip((const Seq**) sqs, pwd, -1);
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
	swap(sqs[0], sqs[1]);
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
		setaaseq(sqs[++j], jgr->sid);
		sqs[j]->left = jgr->lend;
		sqs[j]->right = jgr->rend;
		delete[] sqs[j]->jxt;
		sqs[j]->jxt = new JUXT[jgr->num + 1];
		memcpy(sqs[j]->jxt, jgr->jxt, (jgr->num + 1) * sizeof(JUXT));
		sqs[j]->CdsNo = jgr->num;
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
	as[0] = new const CHAR*[2 * wcp.Nshift];
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
	if (sigm < w) w = sigm;
	prqueue_b->hsort();
	return (w);
}

int SrchBlk::findh(Seq** sqs)
{
	Seq*&	b = sqs[0];	// query
	Seq*&	a = sqs[1];	// translated
	Dhash<INT, int>	hh(2 * pbwc->MaxBlk, 0);
	ORF*	orf = b->getorf();
	if (!orf) return (ERROR);
	Dhash<INT, int>	shash(2 * MaxNref, 0);
	Qwords	qwd(kk,  DRNA, ConvTab, pbwc, bpp);
	INT	norf = 0;
	int	c;
	for (ORF* wrf = orf; wrf->len; ++norf, ++wrf) {
	  bool	meet = false;
	  b->translate(a, *wrf);
	  if (a->right - a->left - (wcp.Nshift + bpp[0]->width) < 1) continue;
	  qwd.reset(a);
	  c = (a->right - a->left) / (wcp.Nshift + wcp.Nshift) - 1;
	  init2(a);
const	  CHAR*	rms = a->at(a->right) - wcp.Nshift * ((c < ptpl)? c: ptpl);
	  INT	nmmc = 0;
	  int	sigpr = 0;
	  for ( ; nmmc < maxmmc && !meet; ++nmmc) {
	    for (int d = 0; d < 2; ++d) {
	      int	e = 1 - d;		// tally
	      int	maxp = 0;
	      for (INT s = 0; s < wcp.Nshift; s++) {
const		CHAR**	ws = bh2->as[d] + s;
const		CHAR*	ms = bh2->as[e][s];
		if (ms > rms) ms = rms;
		int	cscr = 0;
		int	p = 0;
		int	q = 0;
		hh.clear();
		do {		// continueous hits
const		    CHAR*	ss = *ws;
		    if (d)	*ws -= wcp.Nshift;	// right side
		    else	*ws += wcp.Nshift;	// left side
		    if (d ^ (ss >= ms)) {		// has scanned
			meet = true;
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
		    setaaseq(sqs[++c], sh->key);
	    } else
		c = bestref(sqs, sh, c);
	} else
	    c = ERROR;
	delete[] orf;
	return (c);
}

// query and database: proteins or DNAs

int SrchBlk::finds(Seq** sqs)
{
	if (query->right - query->left - (wcp.Nshift + bpp[0]->width) < 1) return (ERROR);
	init2(query);
	int	c = (query->right - query->left) / (wcp.Nshift + wcp.Nshift) - 1;
	INT	testword[2] = {0, 0};
	bool	meet = false;
	BLKTYPE	maxb = 0;
	Dhash<INT, int>	hh(int(2 * pbwc->MaxBlk), 0);
	Qwords	qwd(kk,  DRNA, ConvTab, pbwc, bpp, query);
const	CHAR*	rms = query->at(query->right) - wcp.Nshift * ((c < ptpl)? c: ptpl);
#if TESTRAN
	int	sigpr = 0;
#endif
	INT	nmmc = 0;
	for ( ; nmmc < maxmmc && !meet; ++nmmc) {
	    for (SHORT d = 0; d < 2; ++d) {
	      int	e = 1 - d;		// tally
	      int	maxp = 0;
	      for (INT s = 0; s < wcp.Nshift; s++) {
const		CHAR**	ws = bh2->as[d] + s;
const		CHAR*	ms = bh2->as[e][s];
		if (ms > rms) ms =  rms;
		int	cscr = 0;
		int	p = 0;
		int	q = 0;
		hh.clear();
		do {		// continueous hits
const		    CHAR*	ss = *ws;
		    if (d)	*ws -= wcp.Nshift;	// right side
		    else	*ws += wcp.Nshift;	// left side
		    if (d ^ (ss >= ms)) {		// has scanned
			meet = true;
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
	printf("%d\t%d\t%d\t%d\t%d\t%d\n", nmmc, testword[0] + testword[1], query->len,
		chrsize(maxb), bh2->rscr[maxb], rdbt->randbs(nmmc));
#endif
	if (bh2->sigm) {	// heap sort in the order of block score
	    INT	w = bh2->hsort();
	    if (w < bh2->sigm && OutPrm.supself) ++w;
	    Dhash<INT, int>	shash(2 * w, 0);
	    int	j = 1;
	    for (INT i = 0; i < w; ++i) {
		c = findChrNo((*bh2->prqueue_b)[i].key);
		KVpair<INT, int>*	h = shash.incr(c);
		if (h->val < 2) setaaseq(sqs[j], c);
		if (!OutPrm.supself || 
		    strcmp((*query->sname)[0], (*sqs[j]->sname)[0]) < 0) ++j;
	    }
	    c = j - 1;
	} else
	    c = ERROR;
	return (c);
}

