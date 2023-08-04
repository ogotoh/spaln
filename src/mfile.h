/*****************************************************************************
*
*	Collection of headers commonly used for sequence comparison
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

#ifndef  _MFILE_
#define  _MFILE_

#define BLOCK_SIZE  1024

class Mfile {
	char	*ptr, *cur, *brk;
	size_t	wwd;
	size_t	recno;
public:
	Mfile(size_t rec_size);
	Mfile(const Mfile& src);
	~Mfile() {delete[] ptr;}
	void	write(const UPTR pi);
	void	undo() {cur -= wwd; --recno;}
	UPTR	flush();
	void	reset(long n = 0);
	size_t	size() {return recno;}
	Mfile& operator=(const Mfile& src);
};

#endif
