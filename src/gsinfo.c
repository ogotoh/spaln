/*****************************************************************************
*
*	Manage information on gene organization 
*
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
#include "mfile.h"

PFQ	pfqend = {INT_MAX - 1};

static	PFQ*	fusePfqinaGap(PFQ* dst, int igp);
static	int	cmppos(PFQ* a, PFQ* b);
static	void	oneline(Seq* sd, int i, int len, int pp, double sig, double scr);
static	const	int	MaxClm = 80;

VTYPE	SpbFact;

VTYPE	spb_fact() {return (SpbFact = (VTYPE) (alprm.scale * alprm2.spb));}

PfqItr::PfqItr(Seq* sd, int n) :
	pfqnum(sd->sigII? sd->sigII->pfqnum: 0),
	lstnum(sd->sigII? sd->sigII->lstnum: 0),
	step(sd->sigII? sd->sigII->step: 0),
	pfq(sd->sigII? sd->sigII->pfq: 0),
	lst(sd->sigII? sd->sigII->lst: 0),
	tfq(sd->sigII? sd->sigII->pfq + pfqnum: 0) {
#if USE_WEIGHT
	    weight = sd->weight;
#endif
	    reset(n);
}

#if USE_WEIGHT
PfqItr::PfqItr(SigII& sgi, int n, FTYPE* wt) :
	    pfqnum(sgi.pfqnum), lstnum(sgi.lstnum), step(sgi.step),
	    pfq(sgi.pfq), lst(sgi.lst), tfq(pfq + pfqnum), weight(wt) {
		reset(n);
}
#else
PfqItr::PfqItr(SigII& sgi, int n) :
	    pfqnum(sgi.pfqnum), lstnum(sgi.lstnum), step(sgi.step),
	    pfq(sgi.pfq), lst(sgi.lst), tfq(pfq + pfqnum) {
		reset(n);
}
#endif

Iiinfo::Iiinfo(Seq* seqs[], int m, int n, bool save) :
	a(seqs[0]), b(seqs[1]), 
	sgi(0), cpi(0), agap(0), bgap(0), amany(seqs[0]->many)
{
	api = new PfqItr(a, m);
	bpi = new PfqItr(b, n);
	step = a->isprotein()? 3: 1;
	int	igap = m - n;
	if (igap > 0)	bgap =  igap * step;
	else		agap = -igap * step;
	int	npfq = api->pfqnum + bpi->pfqnum;
	int	nlst = (api->lstnum? api->lstnum: api->pfqnum)
			+ (bpi->lstnum? bpi->lstnum: bpi->pfqnum);
	if (npfq && save) {
	    sgi = new SigII(npfq, nlst, step);
	    cpi = new PfqItr(*sgi);
	}
}

SigII::~SigII()
{
	delete[] pfq; delete[] lst;
	if (eijtab) {delete[] *eijtab; delete[] eijtab;}
	delete[] lone;
}

SigII::SigII(Seq* sd) : 
	pfqnum(0), lstnum(0), 
	pfq(0), lst(0), eijtab(0), lone(0)
{
	if (sd) step = sd->isprotein()? 3: 1;
	else	step = 0;
}

SigII::SigII(int p, int l, int s) : 
	pfqnum(p), lstnum(l), step(s),
	eijtab(0), lone(0)
{
	pfq = p? new PFQ[p + 1]: 0;
	lst = l? new int[l]: 0;
	if (p) vclear(pfq, p + 1);
	if (l) vclear(lst, l);
	else	lstnum = pfqnum;
}

SigII::SigII(int* poss, int nn, int s) : 
	pfqnum(nn), lstnum(0), step(s),
	lst(0), eijtab(0), lone(0)
{
	if (nn) {
	    pfq = new PFQ[nn + 1];
	    for (int i = 0; i < nn; ++i) {
		pfq[i].pos = poss[i];
		pfq[i].num = 1;
		pfq[i].gps = 0;
#if USE_WEIGHT
		pfq[i].dns = 1;
#endif
	    }
	    pfq[pfqnum].num = 0;	// mark last
	} else	pfq = 0;
}

SigII::SigII(Seq** slist, GAPS** gsrc, FTYPE* wtlst)
	: pfqnum(0), lstnum(0), step(0), pfq(0), lst(0), eijtab(0), lone(0)
{
	int	i = 0, j = 0, k = 0;
	for (Seq** sq = slist ; *sq; ++sq, ++k) {
	    if ((*sq)->sigII) {
		i += (*sq)->sigII->pfqnum;
		j += (*sq)->sigII->lstnum;
	    }
	}
	if (!i) return;
	PFQ	pfqbuf = {0, 0, 0};
	step = (*slist)->isprotein()? 3: 1;
	pfq = new PFQ[i + 1];
	vclear(pfq, i + 1);
	lst = j? new int[j]: 0;
	if (j) vclear(lst, j);
	PFQ**	pfqs = new PFQ*[k + 1];
	PFQ**	wsqs = new PFQ*[k];
	GAPS**	glst = new GAPS*[k];
	*pfqs = new PFQ[i + k];
	int**	wlst = new int*[k];
	int	len = gaps_span(*gsrc);
	if (step == 3) len *= 3;
	int*	queue = new int[k + 1];
	int*	amany = new int[k + 1];
	int	many = amany[0] = 0;
	j = 0;
	for (Seq** sq = slist; j < k; ++j, ++sq) {
	    wsqs[j] = pfqs[j];
	    glst[j] = gsrc[j] + 1;
	    amany[j + 1] = (many += (*sq)->many);
	    if ((*sq)->sigII) {
		wlst[j] = (*sq)->sigII->lst;
		vcopy(pfqs[j], (*sq)->sigII->pfq, (*sq)->sigII->pfqnum + 1);
		unfoldPfq(pfqs[j], (*sq)->sigII->pfqnum, gsrc[j], step);
		pfqs[j+1] = pfqs[j] + (*sq)->sigII->pfqnum + 1;
	    } else {
		*pfqs[j] = pfqend;
		wlst[j] = 0;
		pfqs[j+1] = pfqs[j] + 1;
	    }
	    queue[j] = pfqs[j]->pos;
	}
	PFQ*	wfq = pfq;
	int*	wst = lst;
	PrQueue_idx<int>	pq(queue, k, k);
	int	pos = queue[j = pq.gettop_idx()];
	while (pos  < len) {
#if USE_WEIGHT
	    FTYPE	wt = (wtlst && many > 1)? wtlst[j]: 1;
#endif
	    for (i = 0; i < wsqs[j]->num; ++i)
		if (wlst[j]) *wst++ = *wlst[j]++;
		else	*wst++ = amany[j];
	    GAPS*&	gj = glst[j];
	    while (gaps_intr(gj) && step * (gj->gps + gj->gln) < pos) ++gj;
	    int play = (gaps_intr(gj) && (step * (gj->gps + gj->gps)) == pos)?
		step * gj->gln: 0;
	    int diffpos = pos - pfqbuf.pos;
	    if (!diffpos || (diffpos <= play && !(diffpos % 3))) {
		pfqbuf.num += wsqs[j]->num;
#if USE_WEIGHT
		pfqbuf.dns += wt * wsqs[j]->dns;
#endif
	    } else {
		if (pfqbuf.num) *wfq++ = pfqbuf;
		pfqbuf.pos = pos;
		pfqbuf.num = wsqs[j]->num;
#if USE_WEIGHT
		pfqbuf.dns = wt * wsqs[j]->dns;
#endif
	    }
	    queue[j] = (++wsqs[j])->pos;	// step forward
	    pos = queue[j = pq.shift_idx(j)];   // reconstruct heap index
	}
	if (pfqbuf.num) *wfq++ = pfqbuf;
	pfqnum = wfq - pfq;
	lstnum = wst - lst;
	*wfq++ = pfqend;
	delete[] queue;
	delete[] *pfqs;
	delete[] pfqs;
	delete[] wsqs;
	delete[] wlst;
	delete[] glst;
	delete[] amany;
}

void SigII::pfqrepos(RANGE* cr)
{
	if (!cr) return;
	int	getrid = cr->left * step;
	int	gpbias = min(pfq[0].gps, pfq[pfqnum].gps);
	PFQ*	nfq = pfq;
	PFQ*	wfq = pfq;
	int*	wst = lst;
	int*	nst = lst;

	do {
	    while (cr->right < wfq->pos / step) {
		if (!neorng(++cr)) goto loopout;
		getrid += step * (cr->left - cr[-1].right);
	    }
	    if (cr->left <= wfq->pos / step) {
		if (nfq != wfq) *nfq = *wfq;
		nfq->gps -= gpbias;
		(nfq++)->pos -= getrid;
		if (lst)
		    for (int i = 0; i < wfq->num; ++i)
			*nst++ = *wst++;
	    } else if (lst) {
		for (int i = 0; i < wfq->num; ++i) ++wst;
	    }
	} while ((wfq++)->num);
loopout:
	pfqnum = nfq - pfq;
	lstnum = nst - lst;
}

void SigII::locate(PFQ*& wfq, int*& wst, int pos)
{
	wfq = pfq;
	wst = lst;
	pos *= step;
	for ( ; wfq->pos < pos; ++wfq)
	    if (lst) wst += wfq->num;
}

void SigII::swaplst(int an, int bn)
{
	int*	tst = lst + lstnum;

	for (int* wst = lst; wst < tst; ++wst) {
	    if (*wst < an) *wst += bn;
	    else	   *wst -= an;
	}
}

// renumber lst according to the permutation table

void SigII::renumlst(int* lst_odr)
{
	if (!lst) return;
	for (int i = 0; i < lstnum; ++i) lst[i] = lst_odr[lst[i]];
}

void SigII::relist(int bias)
{
	if (!pfqnum) return;
	int*	wst = lst;
	int*	tst;

	if (!lst) {
	    lstnum = pfqnum;
	    wst = lst = new int[lstnum];
	    tst = lst + lstnum;
	    while (wst < tst) *wst++ = bias;
	} else if (bias) {
	    tst = lst + lstnum;
	    while (wst < tst) *wst++ += bias;
	}
}

void SigII::mkeijtab(int many)
{
	if (!pfqnum || !lst || eijtab) return;
	int	n = pfqnum;
	eijtab = new int*[many];
	*eijtab = new int[many * n];
	vclear(*eijtab, many * n);
	for (int i = 1; i < many; ++i) eijtab[i] = eijtab[i-1] + n;
	PFQ*	wfq = pfq;
	int*	wst = lst;
	lone = new int[many];
	vclear(lone, many);
	for (int j = 0; j < pfqnum; ++j, ++wfq) {
	    for (int k = 0; k < wfq->num; ++k, ++wst) {
		eijtab[*wst][j] = wfq->pos;
		if (wfq->num == 1) ++lone[*wst];	// lone eij
	    }
	}
}

void SigII::printlones(FILE* fd, int i, int many)
{
	int	intno = 0;
	for (int j = 0; j < pfqnum; ++j) {
	    if (eijtab[i][j] == 0) continue;
	    ++intno;
	    int	mate = 0;
	    for (int k = 0; k < many; ++k)
		if (eijtab[k][j]) ++mate;
	    if (mate == 1) fprintf(fd, " %d", intno);
	}
}

void SigII::printmates(FILE* fd, int i, int many)
{
	for (int j = 0; j < pfqnum; ++j) {
	    if (eijtab[i][j] == 0) fputs(" 0", fd);
	    else {
		int	mate = 0;
		for (int k = 0; k < many; ++k)
		    if (eijtab[k][j]) ++mate;
		fprintf(fd, " %d", min(mate, 9));
	    }
	}
}

float SigII::eij_dist(int i, int j, int* abc)
{
	int	abc_[3];
	if (!abc) abc = abc_;
	vclear(abc, 3);
	int&	a = abc[0];
	int&	b = abc[1];
	int&	c = abc[2];
	for (int k = 0; k < pfqnum; ++k) {
	    if (eijtab[i][k] && eijtab[j][k]) ++c;
	    if (eijtab[i][k]) ++a;
	    if (eijtab[j][k]) ++b;
	}
	return (a + b)? float(a + b - 2 * c) / (a + b): 0.;
}

int SigII::to_gene_end(int m, bool rend)
{
	if (!pfqnum) return (0);
	m *= step;
	PFQ*	tfq = pfq + pfqnum;
	PFQ*	wfq = rend? tfq: pfq;
	int	term = wfq->pos;

	if (rend) {
	    while ((--wfq)->pos >= m) ;
	    if (wfq == tfq - 1) return (term - m);
	    return (abs(tfq->gps - wfq->gps) - m + wfq->pos);
	} else {
	    while (wfq->pos < m) ++wfq;
	    if (wfq == pfq) return (m - term);
	    --wfq;
	    return (abs(wfq->gps) + m - wfq->pos);
	}
}

int SigII::n_common()
{
	int	c = 0;
	PFQ*	tfq = pfq + pfqnum;
	for (PFQ* wfq = pfq; wfq < tfq; ++wfq)
	    if (wfq->num > 1) c += wfq->num * (wfq->num - 1);
	return (c / 2);
}

#if USE_WEIGHT

void SigII::rescale_dns(VTYPE f)
{
	PFQ*	tfq = pfq + pfqnum;
	for (PFQ* wfq = pfq; wfq < tfq; ++wfq)
	    wfq->dns /= f;
}

void SigII::reset_dns(FTYPE wt)
{
	PFQ*	tfq = pfq + pfqnum;
	for (PFQ* wfq = pfq; wfq < tfq; ++wfq)
	    wfq->dns = wt * wfq->num;
}

#endif

FTYPE* eijdmx(Seq* sd)
{
	FTYPE*	dist = new FTYPE[ncomb(sd->many)];
	FTYPE*	d = dist;
	SigII*	sgi = sd->sigII;
	sgi->mkeijtab(sd->many);
	int	abc[3];
	for (int j = 1; j < sd->many; ++j) {
	    for (int i = 0; i < j; ++i)
		*d++ = sgi->eij_dist(i, j, abc);
	}
	return(dist);
}
		
void fouteijdmx(FILE* fd, Seq* sd, bool dmx)
{
	PrintMember	prm(sd->sname, false, dmx? "\n": " ");
	if (dmx) {
	    for (int i = 0; i < sd->many; ++i)
		prm.put_member(fd, i);
	    fputc('\n', fd);
	}
	char	str[MAXL];
	SigII*	sgi = sd->sigII;
	for (int j = 1; j < sd->many; ++j) {
	    int	clm = 0;
	    for (int i = 0; i < j; ++i) {
		int abc[3];
		float	d = 100. * sgi->eij_dist(i, j, abc);
		if (dmx) {
		    sprintf(str, " %7.3f", d);
		    if (clm >= MaxClm) {
			fputc('\n', fd);
			clm = 0;
		    }
		    fputs(str, fd);
		    clm += 8;
		} else {
		    fprintf(fd, "%7.2f %3d %3d %3d\t", d, abc[2], abc[0], abc[1]);
		    prm.put_member(fd, i); prm.put_member(fd, j);
		    fputc('\n', fd);
		}
	    }
	    if (dmx) fputc('\n', fd);
	}
}

static void fouteij_sumary(FILE* fd, Seq* sd, PrintMember& prm)
{
	int	preblank = sd->sname->longest() + 5;
	fputs("Perdeci", fd);
	for (int i = strlen("Perdeci"); i < preblank; ++i) fputc(' ', fd);
	fputc('\t', fd);
	SigII*	sgi = sd->sigII;
	for (int j = 0; j < sgi->pfqnum; ++j) {
	    int	c = 0;
	    for (int i = 0; i < sd->many; ++i)
		if (sgi->eijtab[i][j]) ++c;
	    c = 10 * c / sd->many;
	    if (c < 10) fprintf(fd, "%d", c);
	    else	fputc('*', fd);
	}
	fputs("\nKingdoms", fd);
	for (int i = strlen("Kingdoms"); i < preblank; ++i) fputc(' ', fd);
	fputc('\t', fd);
	for (int j = 0; j < sgi->pfqnum; ++j) {
	    int	kingdom[4];
	    vclear(kingdom, 4);
	    for (int i = 0; i < sd->many; ++i) {
		if (!sgi->eijtab[i][j]) continue;
		switch (*prm[i]) {
		  case 'A': ++kingdom[0]; break;
		  case 'F': ++kingdom[1]; break;
		  case 'P': ++kingdom[2]; break;
		  case 'O': ++kingdom[3]; break;
		}
	    }
	    int	c = 0;
	    for (int i = 0; i < 4; ++i) if (kingdom[i]) ++c;
	    fprintf(fd, "%d", c);
	}
	fputs("\nPhyla", fd);
	for (int i = strlen("Phyla"); i < preblank; ++i) fputc(' ', fd);
	fputc('\t', fd);
	for (int j = 0; j < sgi->pfqnum; ++j) {
	    StrHash<int>	sh(sd->many / 10);
	    for (int i = 0; i < sd->many; ++i)
		if (sgi->eijtab[i][j]) sh.incr(prm[i]);
	    int	c = 0;
	    for (KVpair<INT, int>* kv = sh.begin(); kv < sh.end(); ++kv)
		if (kv->val != sh.undef()) ++c;
	    if (c < 10) c += '0';
	    else if (c < 36) c += 'a' - 10;
	    else if (c < 52) c += 'A' - 36;
	    else c = '*';
	    fputc(c, fd);
	}
	fputc('\n', fd);
}

void fouteij(FILE* fd, Seq* sd)
{
	SigII*	sgi = sd->sigII;
	if (!sgi || !sgi->pfqnum) return;
	if (!fd) fd = qout("");
	if (!fd) return;
	if (!sgi->lst) {
	    fprintf(fd, "%s\t", sd->sqname());
	    PFQ*	pfq = sgi->pfq;
	    for (int i = 0; i < sgi->pfqnum; ++i, ++pfq)
		fprintf(fd, " %d", pfq->pos);
	    fputc('\n', fd);
	    return;
	}
	sgi->mkeijtab(sd->many);
	int	buf[2] = {0, -1};
	Seq*	memsd = 0;
	if (algmode.nsa == 14 || algmode.nsa == 15) {
	    fouteijdmx(fd, sd, algmode.nsa == 14);
	    return;
	} else if (algmode.nsa >= 8 && algmode.nsa <= 11) {
	    memsd = new Seq(1);
	} else if ((algmode.nsa == 0 || algmode.nsa == 2) && OutPrm.deflbl) {
	    PFQ*	pfq = sgi->pfq;
	    fputs("SPB\t", fd);
	    if (algmode.nsa & 4) fputs("No. Ni  Li :", fd);
	    for (int j = 0; j < sgi->pfqnum; ++j, ++pfq)
		fprintf(fd, " %d", pfq->pos / 3);
	    fputc('\n', fd);
	}
	PrintMember	prm(sd->sname, true, "\t");
	for (int i = 0; i < sd->many; ++i) {
	    prm.put_member(fd, i);
	    if (4 <= algmode.nsa && algmode.nsa < 8) {
		int	nint = 0;
		for (int j = 0; j < sgi->pfqnum; ++j)
		    if (sgi->eijtab[i][j]) ++nint;
		fprintf(fd, "%3d %3d %3d : ", i + 1, nint, sgi->lone[i]);
	    }
	    if (algmode.nsa == 4) sgi->printlones(fd, i, sd->many);
	    if (algmode.nsa == 5) sgi->printmates(fd, i, sd->many);
	    PFQ*	pfq = sgi->pfq;
	    Seq*	msd = memsd? memsd->getseq((*sd->sname)[i]): 0;
	    PFQ*	mfq = (msd && msd->sigII)? msd->sigII->pfq: 0;
	    bool	anti = !mfq || (mfq[0].gps >  mfq[1].gps);
	    int		cds = ((algmode.nsa == 8 || algmode.nsa == 9) && !anti)? mfq->pos: 0;
	    int 	prp = 0;
	    for (int j = 0; j < sgi->pfqnum; ++j, ++pfq) {
		switch (algmode.nsa) {
		  case 0:
		    fprintf(fd, "%d", sgi->eijtab[i][j] != 0); break;
		  case 1: 
		    if (sgi->eijtab[i][j]) fprintf(fd, " %d", pfq->pos);
		    break;
		  case 2: case 6:
		    if (sgi->eijtab[i][j]) fprintf(fd, "%d", pfq->pos % 3);
		    else	fprintf(fd, "-");
		    break;
		  case 3: case 7:
		    if (sgi->eijtab[i][j]) fprintf(fd, " %3d %d",
			pfq->pos / 3, pfq->pos % 3);
		    break;
		  case 8: case 9:
		    if (j) fputc('\t', fd);
		    if (mfq && sgi->eijtab[i][j]) {
			int	intlen = (anti? mfq->gps - mfq[1].gps - mfq->pos:
				mfq[1].gps - mfq->gps - mfq[1].pos) + cds;
			if (algmode.nsa == 8)
			    fprintf(fd, "%7d", intlen);
			else
			    fprintf(fd, "%5d:%d", intlen, pfq->pos % 3);
			cds = anti? (mfq++)->pos: (++mfq)->pos;
		    } else	fputs("      -", fd);
		    break;
		  case 10:
		    if (j) fputc('\t', fd);
		    if (mfq && sgi->eijtab[i][j]) {
			fprintf(fd, "%7d", mfq->pos - cds);
			cds = (mfq++)->pos;
		    } else	fputs("      -", fd);
		    break;
		  case 11:
		    if (j) fputc('\t', fd);
		    if (sgi->eijtab[i][j]) {
			int	nogap = sd->countgap(i, (prp + 1) / 3, (pfq->pos + 1) / 3);
			int	exlen = pfq->pos - prp - 3 * nogap;
			if (mfq) {
			    exlen -= (mfq->pos - cds);
			    cds = (mfq++)->pos;
			}
			fprintf(fd, "%7d", exlen);
			prp = pfq->pos;
		    } else	fputs("      -", fd);
		    break;
		  case 12:
		    if (sgi->eijtab[i][j]) {
			*buf = (pfq->pos + 1) / 3;
			sd->pos2num(i, buf);
			fprintf(fd, " %d", 3 * *buf + (pfq->pos + 1) % 3 - 1);
		    }
		    break;
		  case 13:
		    if (sgi->eijtab[i][j]) {
			*buf = (pfq->pos + 1) / 3;
			sd->pos2num(i, buf);
			fprintf(fd, " %3d %d", *buf, pfq->pos % 3);
		    }
		    break;
		  default: break;
		}
	    }
	    if (algmode.nsa == 6) {
		if (mfq && mfq->num) 
		    fprintf(fd, "\t%7d", mfq->pos - cds);
	    }
	    fputc('\n', fd);
	}
	if (OutPrm.taxoncode && (algmode.nsa == 0 || algmode.nsa == 2))
	   fouteij_sumary(fd, sd, prm);
	delete memsd;
}

VTYPE Iiinfo::StoreIIinfo(int m, int n)
{
	int	zero = 0;
	VTYPE	scr = 0;
	bool	anend = !api->end() && *api < m;
	bool	bnend = !bpi->end() && *bpi < n;
	while (anend || bnend) {
	    int	apos = api->wfq?  api->wfq->pos + agap: pfqend.pos;
	    int	bpos = bpi->wfq?  bpi->wfq->pos + bgap: pfqend.pos;
	    if (anend && bnend && apos == bpos) {
#if USE_WEIGHT
		PFQ	pfqbf = {apos, api->wfq->num + bpi->wfq->num, 
			0, api->wfq->dns + bpi->wfq->dns};
		scr += api->wfq->dns * bpi->wfq->dns;
#else
		PFQ	pfqbf = {apos, api->wfq->num + bpi->wfq->num, 0};
		scr += api->wfq->num * bpi->wfq->num;
#endif
		if (cpi) *cpi->wfq++ = pfqbf;
	    } else if (anend && cpi && apos < bpos) {
#if USE_WEIGHT
		PFQ	pfqbf = {apos, api->wfq->num, 0, api->wfq->dns};
#else
		PFQ	pfqbf = {apos, api->wfq->num, 0};
#endif
		*cpi->wfq++ = pfqbf;
	    } else if (bnend && cpi && bpos < apos) {
#if USE_WEIGHT
		PFQ	pfqbf = {bpos, bpi->wfq->num, 0, bpi->wfq->dns};
#else
		PFQ	pfqbf = {bpos, bpi->wfq->num, 0};
#endif
		*cpi->wfq++ = pfqbf;
	    }
	    if (anend && apos <= bpos) {
		if (cpi) {
		    if (api->wst) {
			for (int j = 0; j < api->wfq->num; ++j)
			    *cpi->wst++ = api->wst[j];
		    } else {
			*cpi->wst++ = zero;
		    }
		}
		++(*api);
		anend = !api->end() && *api < m;
	    }
	    if (bnend && bpos <= apos) {
		if (cpi) {
		    if (bpi->wst) {
			for (int j = 0; j < bpi->wfq->num; ++j) {
			    int	lstbf = amany + bpi->wst[j];
			    *cpi->wst++ = lstbf;
			}
		    } else {
			*cpi->wst++ = amany;
		    }
		}
		++(*bpi);
		bnend = !bpi->end() && *bpi < n;
	    }
	}
	return (SpbFact * scr);
}

SigII* Iiinfo::finalize(int len)
{
	if (!cpi) return (0);
	*cpi->wfq = pfqend;
	cpi->wfq->pos = step * len;
	sgi->pfqnum = cpi->pfq? cpi->wfq - cpi->pfq: 0;
	sgi->lstnum = cpi->lst? cpi->wst - cpi->lst: 0;
	SigII*	rv = sgi;
	sgi = 0;
	return (rv);
}

bool Gsinfo::intronless()
{
	if (!eijnc) return false;
	EISCR*	fst = eijnc->begin();
	return ((noeij < 2) && (fstat.unp < IntronPrm.llmt) &&
	    (fst->unp5 + fst->mmc5 <= end_error_thr) &&
	    (fst->unp3 + fst->mmc3 <= end_error_thr));
}

void Gsinfo::SaveGsInfo(Iiinfo* iif, int len)
{
	delete sigII;
	sigII = iif? iif->finalize(len): 0;
}

Gsinfo::Gsinfo(SKL* s) :
	end_error_thr(int(alprm2.jneibr * 0.8)), scr(0), rscr(0),
	skl(s), noeij(0), CDSrng(0), eijnc(0), cigar(0), vlgar(0),
	samfm(0), sigII(0)
{
	vclear(&fstat);
}

Gsinfo::~Gsinfo()
{
	delete[] skl;
	delete[] CDSrng;
	delete eijnc;
	delete sigII;
	delete cigar;
	delete vlgar;
	delete samfm;
}

int Gsinfo::center(int k)
{
	if (!eijnc) return (0);
	EISCR*  eij = eijnc->begin();
	EISCR*  lst = eij + noeij - 1;
	return ((k? eij->left + lst->right: eij->rleft + lst->rright) / 2);
}

RANGE* Gsinfo::eiscr2rng()
{
	if (!eijnc) return (0);
	EISCR*	eij = eijnc->begin();
	CDSrng = new RANGE[noeij + 2];
	RANGE*  wrng = CDSrng;
	(wrng++)->left = noeij;
	while (neoeij(eij)) {
	    wrng->left = eij->left;
	    (wrng++)->right = (eij++)->right;
	}
	wrng->left = eij->left;
	wrng->right = eij->right;
	return (CDSrng);
}

RANGE* Gsinfo::eiscrunfold(GAPS* gp)
{
	int	gap = 0;
	EISCR*	eij = eijnc->begin();
	RANGE*	rng = new RANGE[noeij + 2];
	RANGE*	wrk = rng;

	(wrk++)->left = noeij;
	for ( ; neoeij(eij); ++eij, ++wrk) {
	    while (gp->gps < eij->left) gap += (gp++)->gln;
	    wrk->left = eij->left + gap;
	    while (gp->gps < eij->right) gap += (gp++)->gln;
	    wrk->right = eij->right + gap;
	}
	*wrk = endrng;
	return (rng);
}

// infer exon-exon structure of the reference

RANGE* Gsinfo::querygs(Seq* qry)
{
	int	step = qry->isprotein()? 3: 1;
	RANGE*	rng = new RANGE[noeij + 2];
	RANGE*  wrng = rng;
	EISCR*	eij = eijnc->begin();
	(wrng++)->left = noeij;
	int	gps = 0;
	int	rps = 0;
	for (int cps = rps; neoeij(eij); ++eij) {
	    wrng->left = cps;
	    cps = step * (eij->rright - rps);
	    gps += eij->right - eij->left;
	    if (step == 3) {
		switch (gps % 3) {
		  case 1: ++cps; break;
		  case 2: --cps; break;
		  default: break;
		}
	    }
	    (wrng++)->right = cps;
	}
	*wrng = endrng;
	return (rng);
}

static void oneline(Seq* sd, int i, int len, int pp, double sig, double scr)
{
	char*	decode = sd->inex.molc == TRON? ncodon: nucl;
	int	ll = i - BoundRng;
	int	rr = i + BoundRng;
	int	l = ll;
	int	r = rr;
	if (ll < 0) ll = 0;
	if (rr > sd->len) rr = sd->len;
	CHAR*	ps = sd->at(ll);
	fprintf(out_fd, "%6d  ", sd->SiteNo(ll));
	for ( ; l < ll; ++l) putc(' ', out_fd);
	for ( ; l < i; ++l, ++ps)
	    putc(decode[*ps], out_fd);
	putc(' ', out_fd);
	for ( ; l < rr; ++l, ++ps)
	    putc(decode[*ps], out_fd);
	while (l++ < r) putc(' ', out_fd);
	fprintf(out_fd, "%8d %6d %9.2f ", len, pp, sig);
	if (scr == INT_MIN) putc('\n', out_fd);
	else fprintf(out_fd, "%9.2f\n", scr);
}

void Gsinfo::BoundarySeq(Seq* sd)
{
	RANGE*	wrng = fistrng(CDSrng);
	double	sig;
	int	len = -wrng->left;
	int	clen = 0;

	RANGE*	rng = wrng;
	fputs("//\n", out_fd);
	for ( ; neorng(wrng + 1); ++wrng) {
	    clen += wrng->right - wrng->left;
	    int	intv = wrng[1].left - wrng->right;
	    if (intv > 2) {
		len += wrng->right;
		sig = sd->exin->sig53(wrng->right, 0, IE5);
		oneline(sd, wrng->right, len, clen / 3, sig / sd->exin->fact, INT_MIN);
		len = -wrng[1].left;
	    }
	}
	len += wrng->right;
	putc('\n', out_fd);
	clen = rng->right - rng->left;
	wrng = rng;
	while (neorng(wrng + 1)) {
	    ++wrng;
	    int	intv = wrng->left - wrng[-1].right;
	    if (intv > 2) {
		sig = sd->exin->sig53(wrng[-1].right, wrng->left, IE53);
		oneline(sd, wrng->left, intv, clen % 3, sig / sd->exin->fact, INT_MIN);
	    }
	    clen += wrng->right - wrng->left;
	}
	putc('\n', out_fd);
	sig = sd->exin->sigST(rng->left + 1, true);
	oneline(sd, rng->left, 0, clen % 3, sig, INT_MIN);
	sig = sd->exin->sigST(wrng->right + 1, false);
	oneline(sd, wrng->right, len, clen / 3, sig, INT_MIN);
}

void Gsinfo::BoundaryInf(Seq* sd)
{
	if (!eijnc) {
	    BoundarySeq(sd);
	    return;
	}
	EISCR*	wrng = eijnc->begin();
	EISCR*	trng = wrng + noeij - 1;
	EXIN*	bb = sd->exin->data;
	int	len = -wrng->left;
	int	clen = 0;
	double	tscr = 0.;

	fputs("//\n", out_fd);
	for ( ; wrng != trng; ++wrng) {
	    clen += wrng->right - wrng->left;
	    int	intv = wrng[1].left - wrng->right;
	    if (intv > 2) {
		len += wrng->right;
		double	sig = wrng->sig5;
		tscr += wrng->escr - sig;
		sig /= sd->exin->fact;
		double	scr = (double) wrng->escr / sd->exin->fact;
		oneline(sd, wrng->right, len, clen / 3, sig, scr);
		len = -wrng[1].left;
	    }
	}
	len += wrng->right;
	putc('\n', out_fd);
	wrng = eijnc->begin();
	clen = wrng->right - wrng->left;
	while (wrng != trng) {
	    int	intv = wrng[1].left - wrng->right;
	    if (intv > 2) { 
		double	sig = wrng[1].sig3;
		tscr += wrng->iscr - sig;
		sig /= sd->exin->fact;
		double	scr = (double) wrng->iscr / sd->exin->fact;
		oneline(sd, wrng[1].left, intv, clen % 3, sig, scr);
	    }
	    ++wrng;
	    clen += wrng->right - wrng->left;
	}
	putc('\n', out_fd);
	wrng = eijnc->begin();
	double	sig = bb? bb[wrng->left + 1].sigS / sd->exin->fact: 0.;
	tscr = (tscr + trng->escr) / sd->exin->fact;
	oneline(sd, wrng->left, 0, clen % 3, sig, tscr);
	sig = bb? bb[trng->right + 1].sigT / sd->exin->fact: 0.;
	double	scr = (double) trng->escr / sd->exin->fact;
	oneline(sd, trng->right, len, clen / 3, sig, scr);
}

void cutSigII(Seq* dest, Seq* sorc)
{
	if (dest == sorc) fatal("bad cutSigII!\n");
	SigII*&	src = sorc->sigII;
	SigII*&	dst = dest->sigII; 
	delete dst; dst = 0;
	if (!src || !src->pfq) return;
	PFQ*	pfq = src->pfq;
	int*	lst = src->lst;
	int	bias = sorc->left * src->step;
	int	to = sorc->pfqPos(sorc->left);
	while (pfq->num && pfq->pos < to) {
	    if (lst) lst += pfq->num;
	    ++pfq;
	}
	int	gpbias = pfq->gps;
	to = sorc->pfqPos(sorc->right);
	int	nfq = src->pfqnum - (pfq - src->pfq);
	int	nst = src->lstnum - (lst - src->lst);
	if (!nfq) return;
	dst = new SigII(nfq, nst, src->step);
	PFQ*	wfq = dst->pfq;
	int*	wst = dst->lst;
	while (pfq->num && pfq->pos < to) {
	    PFQ pfqbf = *pfq;
	    pfq->pos -= bias;
	    pfq->gps -= gpbias;
	    *wfq++ = pfqbf;
	    if (lst)
		for (int k = 0; k < pfq->num; ++k)
		    *wst++ = *lst++;
	    ++pfq;
	}
	dst->pfqnum = wfq - dst->pfq;
	if (!dst->pfqnum) {
	    delete dst; dst = 0;
	    return;
	}
	wfq->pos = pfq->pos - bias;
	wfq->num = 0;
	wfq->gps = pfq->gps - gpbias;
	if (lst) dst->lstnum = wst - dst->lst;
}

SigII* copySigII(SigII* src)
{
	if (!src || !src->pfqnum) return (0);
	SigII*	dst = new SigII;
	dst->pfqnum = src->pfqnum;
	dst->lstnum = src->lstnum;
	dst->step = src->step;
	dst->pfq = new PFQ[dst->pfqnum + 1];
	vcopy(dst->pfq, src->pfq, src->pfqnum + 1);
	if (src->lst) {
	    dst->lst = new int[dst->lstnum];
	    vcopy(dst->lst, src->lst, src->lstnum);
	} else
	    dst->lst = 0;
	return (dst);
}

void catSigII(Seq* dest, Seq* sorc, int bias)
{
	SigII*&	head = dest->sigII;
	SigII*&	tail = sorc->sigII;
	bias += dest->right - sorc->left;
	if (!tail) {
	    if (head) {
		bias += sorc->right - sorc->left;
		if (sorc->isprotein()) bias *= 3;
		PFQ*	hfq = head->pfq + head->pfqnum;
		hfq->pos = bias;
	    }
	    return;
	} else if (sorc->isprotein()) bias *= 3;
	if (!head) {
	    head = copySigII(tail);
	    PFQ*	pfq = head->pfq;
	    PFQ*	wfq = pfq + head->pfqnum;
	    for ( ; pfq <= wfq; ++pfq) pfq->pos += bias;
	    return;
	}
	PFQ*	pfq = tail->pfq;
	int*	lst = tail->lst;
	int	to = sorc->pfqPos(sorc->left);
	for ( ; pfq->num && pfq->pos < to; ++pfq)
	    if (lst) lst += pfq->num;
	int*	wst = lst;
	PFQ*	wfq = pfq;
	to = sorc->pfqPos(sorc->right);
	for ( ; wfq->num && wfq->pos < to; ++wfq)
	    if (wst) wst += wfq->num;
	to = head->pfqnum + wfq - pfq;
	if (!to) return;
	int	gpbias = pfq->gps;
	PFQ*	tmp = new PFQ[to + 1];
	vcopy(tmp, head->pfq, head->pfqnum + 1);
	delete[] head->pfq;
	head->pfq = tmp;
	PFQ*	hfq = head->pfq + head->pfqnum;
	head->pfqnum = to;
	for ( ; pfq <= wfq; ++pfq, ++hfq) {
	    *hfq = *pfq;
	    hfq->pos += bias;
	    hfq->gps += gpbias;
	}
	to = head->lstnum;
	head->lstnum += wst - lst;
	if (lst) {
	    int*	tmp = new int[head->lstnum];
	    memcpy(tmp, head->lst, to * sizeof(int));
	    delete[] head->lst;
	    head->lst = tmp;
	    vcopy(head->lst + to, lst, wst - lst);
	}
}

SigII* extSigII(Seq* sorc, int* which, FTYPE nfact, bool renum_lst)
{
	SigII*	src = sorc->sigII;
	if (!src) return (0);
	if (!which) return copySigII(src);

	PFQ*	pfq = src->pfq;
	int*	lst = src->lst;
	int	bias = sorc->left;
	bias *= src->step;
	int	to = sorc->pfqPos(sorc->left);
	while (pfq->num && pfq->pos < to) {
	    if (lst) lst += pfq->num;
	    ++pfq;
	}
	int	nfq = src->pfqnum - (pfq - src->pfq);
	if (!nfq) return (0);

	int*	chosen = new int[sorc->many];
	vset(chosen, -1, sorc->many);
	int	n = 0;
	for (int* w = which; *w >= 0; ++w, ++n)
	    chosen[*w] = renum_lst? n: *w;
	int	nst = lst? src->lstnum - (lst - src->lst): 0;
	SigII*	dst = new SigII(nfq, nst, src->step);
	PFQ*	wfq = dst->pfq;
	int*	wst = dst->lst;

	to = sorc->pfqPos(sorc->right);
	for ( ; pfq->num && pfq->pos < to; ++pfq) {
	    PFQ	pfqbf = {pfq->pos - bias, 0, 0};
#if USE_WEIGHT
	    pfqbf.dns = (lst && sorc->weight)? 0: pfq->num;
#endif
	    if (lst) {
		for (int k = 0; k < pfq->num; ++k, ++lst) {
		    int	m = chosen[*lst];
		    if (m >= 0) {
			*wst++ = m;
			++pfqbf.num;
#if USE_WEIGHT
			if (n == 1) pfqbf.dns = 1;
			else if (sorc->weight) pfqbf.dns += sorc->weight[m];
#endif
		    }
		}
#if USE_WEIGHT
		if (n > 1) pfqbf.dns /= nfact;
#endif
	    } else pfqbf.num = 1;
	    if (pfqbf.num) *wfq++ = pfqbf;
	}
	delete[] chosen;
	dst->pfqnum = wfq - dst->pfq;
	if (!dst->pfqnum) {
	    delete dst;
	    return 0;
	}
	*wfq = pfqend;
	if (lst) dst->lstnum = wst - dst->lst;
	return (dst);
}

/*******************************************
	|		||   |  intron (unfold -> fold)
_______V_____________VV___V__________________
	   ||	|  |	gap (unfold)

*********************************************/

static int cmppos(PFQ* a, PFQ* b)
{
	return (a->pos - b->pos);
}

static PFQ* fusePfqinaGap(PFQ* dst, int igp)
{
	PFQ*	gfq = dst - igp;

	qsort((UPTR) gfq, (INT) igp, sizeof(PFQ), (CMPF) cmppos);
	for (dst = gfq++; --igp; ++gfq) {
	    if (dst->pos == gfq->pos) dst->num += gfq->num;
	    else	(++dst)->num = gfq->num;
	}
	return (dst);
}

void SigII::rmGapPfq(GAPS* gp)
{
	if ((gp++)->gln == 3) return;	/* no gap no action */

	int	gpos = step * gp->gps;
	int	gupp = step * (gp->gps + gp->gln);
	PFQ*	wfq = pfq;
	PFQ*	dst = pfq;
	PFQ*	end = pfq + pfqnum;
	int	glen = 0;
	int	igp = 0;
	int	play = step == 3? 1: 0;
	int	phs = 0;

	while (wfq <= end) {
	    if (!gaps_intr(gp) || gpos + play >= wfq->pos) {
		*dst = *wfq;
		dst->pos -= glen;
		++wfq;
		++dst;
	    } else if (gupp + play >= wfq->pos) {
		if (step == 3) phs = (wfq->pos + 1) % 3 - 1;
		*dst = *wfq;
		dst->pos = gpos + phs - glen;
		++wfq;
		++dst;
		igp++;
	    } else {
		if (igp > 1) dst = fusePfqinaGap(dst, igp);
		glen += step * gp->gln;
		igp = 0;
		++gp;
		gpos = step * gp->gps;
		gupp = gpos + step * gp->gln;
	    }
	}
	pfqnum += dst - wfq;
}

VTYPE spSigII(Seq* sd)
{
	if (!sd->sigII || !sd->sigII->pfqnum || SpbFact <= 0) return (0);
#if USE_WEIGHT
	FTYPE*	pw = sd->pairwt;
#else
	FTYPE*	pw = 0;
#endif
	int	lb = sd->pfqPos(sd->left);
	int	rb = sd->pfqPos(sd->right);

	PFQ*	pfq = sd->sigII->pfq;
	int*	lst = sd->sigII->lst;
	if (!lst || !pw) {
	    for ( ; pfq->pos < lb; ++pfq) ;
	    int	n = 0;
	    for ( ; pfq->pos < rb; ++pfq)
		n += ncomb(pfq->num);
	    return ((VTYPE) (SpbFact * n));
	}
#if USE_WEIGHT
	for ( ; pfq->pos < lb; ++pfq) lst += pfq->num;
	VTYPE	sum = 0;
	for ( ; pfq->pos < rb; ++pfq) {
	    for (int j = 1; j < pfq->num; ++j) {
		for (int i = 0; i < j; ++i)
		    sum += sd->pairwt[elem(lst[i], lst[j])];
	    }
	}
	return ((VTYPE) (SpbFact * sum));
#else
	for ( ; pfq->pos < lb; ++pfq) ;
	int	n = 0;
	for ( ; pfq->pos < rb; ++pfq)
	    n += ncomb(pfq->num);
	return ((VTYPE) (SpbFact * n));
#endif
}

void unfoldPfq(PFQ* pfq, int num, GAPS* gg, int step)
{
	if ((gg++)->gln == 3) return;		// no gap to be inserted
	PFQ*	tfq = pfq + num;
	for (int glen = 0; pfq <= tfq; ++pfq) {
	    while (gaps_intr(gg) && step * gg->gps <= pfq->pos + glen)
		glen += step * (gg++)->gln;
	    pfq->pos += glen;
	}
}

void Vulgar::postproc()
{
	VULGAR*	tm = rec + --num;	// ignore the last dummy record
	for (VULGAR* wk = rec; wk < tm; ++wk) {
	   if (wk->ope == 'M' || wk->ope == 'D') {
		if (wk[1].ope == 'S' && wk[1].blen == 2) {
		    --wk->alen;
		    wk->blen -= 3;
		}
	    } else if (wk->ope == 'F') {
		if (wk->alen == 1) {
		    wk[-1].blen -= 2;
		    wk->blen += 2;
		} else if (wk->alen == 2) {
		    wk[-1].blen -= 1;
		    wk->blen += 1;
		    wk->alen = 1;
		}
	    }
	}
}

Eijnc::Eijnc(bool que) : num(0), rec(0), fstque(0), q(0)
{
	emfd = new Mfile(sizeof(EISCR));
	if (que) {
	    fstque = new FSTAT[alprm2.jneibr];
	    vclear(fstque, alprm2.jneibr);
	}
}

void Eijnc::store(EISCR& rbuf, FSTAT& now, FSTAT& prv, bool nearjnc)
{
	rbuf.mch = (int) (now.mch - prv.mch);
	rbuf.mmc = (int) (now.mmc - prv.mmc);
	rbuf.gap = (int) (now.gap - prv.gap);
	rbuf.unp = (int) (now.unp - prv.unp);
	if (nearjnc) {
	    rbuf.mch5 = rbuf.mch;
	    rbuf.mmc5 = rbuf.mmc;
	    rbuf.gap5 = rbuf.gap;
	    rbuf.unp5 = rbuf.unp;
	}
	rbuf.mch3 = (int) (now.mch - fstque[q].mch);
	rbuf.mmc3 = (int) (now.mmc - fstque[q].mmc);
	rbuf.unp3 = (int) (now.unp - fstque[q].unp);
	rbuf.gap3 = (int) (now.gap - fstque[q].gap);
	vclear(fstque, alprm2.jneibr);
	q = 0;
}

void Eijnc::shift(EISCR& rbuf, FSTAT& now, FSTAT& prv, bool nearjnc)
{
	if (nearjnc) {
	    rbuf.mch5 = (int) (now.mch - prv.mch);
	    rbuf.mmc5 = (int) (now.mmc - prv.mmc);
	    rbuf.unp5 = (int) (now.unp - prv.unp);
	    rbuf.gap5 = (int) (now.gap - prv.gap);
	}
	fstque[q] = now;
	if (++q == alprm2.jneibr) q = 0;
}
