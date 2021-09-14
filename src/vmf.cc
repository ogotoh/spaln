/*****************************************************************************
*
*	A self-expanding memory block
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

#include "cmn.h"
#include "vmf.h"

int	MaxVmfSpace = DefMaxVMF;

void setVmfSpace(int spc)
{
	MaxVmfSpace = spc;
}

VMFBLK* Vmf::newblk()
{
	VMFBLK*	    prv = blk;

	cur = new SKLP[VMFBLKSIZE];
	brk = cur + VMFBLKSIZE;
	blk = new VMFBLK;
	blk->ptr = cur;
	blk->nxt = 0;
	if (prv) prv->nxt = blk;
	return (blk);
}

Vmf::Vmf()
{
	idx = 0;
	blk = 0;
	hed = newblk();
}

void Vmf::freeblk(VMFBLK* blk)
{
	if (blk) {
	    freeblk(blk->nxt);
	    delete[] blk->ptr;
	    delete blk;
	}
}

Vmf::~Vmf()
{
	freeblk(hed);
}

int Vmf::writevmf(SKLP* rec)
{
	if (idx == INT_MAX)
	    fatal("Too many Vmf records !");
	if (cur == brk) newblk();
	movmem(rec, cur, sizeof(SKLP));
	cur++;
	return(idx++);
}

int Vmf::vmfseek(int recno)
{
	if (recno < 0) recno += idx;
	if (0 > recno || recno > idx) {
	    prompt("%ld > %ld: Vmf out of range !\n", recno, idx);
	    return (ERROR);
	}
	int	rn = 0;
	for (blk = hed; blk; blk = blk->nxt)
	    if ((rn += VMFBLKSIZE) > recno) break;
	cur = blk->ptr + (recno - rn + VMFBLKSIZE);
	idx = recno;
	return(recno);
}

int Vmf::readvmf(SKLP* rec, int recno)
{
	if (vmfseek(recno) == ERROR) return(ERROR);
	movmem(cur, rec, sizeof(SKLP));
	cur++;
	return((int) sizeof(SKLP));
}

SKLP& Vmf::operator[](int recno)
{
	if (recno < 0) recno += idx;
	if (0 > recno || recno > idx) {
	    prompt("%ld > %ld: Vmf out of range !\n", recno, idx);
	    return (cur[0]);
	}
	int	rn = 0;
	VMFBLK*	bk = hed;
	for ( ; bk; bk = bk->nxt)
	    if ((rn += VMFBLKSIZE) > recno) break;
	if (!bk) {
	    prompt("%ld > %ld: Vmf out of range !\n", recno, idx);
	    return (cur[0]);
	}
	return (bk->ptr[recno - rn + VMFBLKSIZE]);
}

SKL* Vmf::vmferror(Mfile& mfd) 
{
	SKL*    skl = (SKL*) mfd.flush();
	delete[] skl;
	return (0);
}

SKL* Vmf::traceback(int pp)
{
	SKLP	sv;
	Mfile	mfd(sizeof(SKL));
	mfd.write(&sv);
	if (readvmf(&sv, pp) == ERROR) return (vmferror(mfd));
	mfd.write(&sv);
	while (sv.p) {
	    if (readvmf(&sv, sv.p) == ERROR) return (vmferror(mfd));
	    mfd.write(&sv);
	}
	int	n = mfd.size();
	SKL*	skl = (SKL*) mfd.flush();
	skl->n = --n;
	return (skl);
}
