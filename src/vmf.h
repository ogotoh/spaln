/*****************************************************************************
*
*	Header to vmf.c
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

#ifndef _VMF_H_
#define _VMF_H_

static	const	int	VMFBLKSIZE = 7680;
static	const	long	DefMaxVMF = 16 * 1024 * 1024;

extern	long	MaxVmfSpace;
extern	void	setVmfSpace(long spc);

struct VMFBLK {
	SKLP*	ptr;
	VMFBLK* nxt;
};

class Vmf {
	SKLP*	cur;
	SKLP*	brk;
	VMFBLK* hed;
	VMFBLK*	blk;
	long    idx;
	VMFBLK*	newblk();
	void	freeblk(VMFBLK* blk);
public:
	Vmf();
	~Vmf();
	long	add(int m, int n, long p) {
		SKLP	tmp = {m, n, p};
		return (writevmf(&tmp));
	}
	long	writevmf(SKLP* rec);
	long	vmfseek(long recno);
	int	readvmf(SKLP* rec, long recno);
	SKL*	traceback(long pp);
	SKL*	vmferror(Mfile& mfd);
};
#endif
