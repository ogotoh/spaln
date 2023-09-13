/*****************************************************************************
*
*	Header to vmf.c
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

#ifndef _VMF_H_
#define _VMF_H_

static	const	int	VMFBLKSIZE = 7680;
static	const	int	DefMaxVMF = 32 * 1024 * 1024;

extern	int	MaxVmfSpace;
extern	void	setVmfSpace(int spc);

struct VMFBLK {
	SKLP*	ptr;
	VMFBLK* nxt;
};

class Vmf {
	SKLP*	cur;
	SKLP*	brk;
	VMFBLK* hed;
	VMFBLK*	blk;
	int	idx;
	VMFBLK*	newblk();
	void	freeblk(VMFBLK* blk);
public:
	Vmf();
	~Vmf();
	int	add(int m, int n, int p) {
		SKLP	tmp = {m, n, p};
		return (writevmf(&tmp));
	}
	int	writevmf(SKLP* rec);
	int	vmfseek(int recno);
	int	size() const {return (idx);}
	int	readvmf(SKLP* rec, int recno);
	SKL*	traceback(int pp);
	SKL*	vmferror(Mfile& mfd);
	SKLP&	operator[](int recno);
};
#endif
