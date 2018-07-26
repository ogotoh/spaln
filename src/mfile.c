/*****************************************************************************
*
*	Simulate a block of memory as if a file
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

#include "stdtype.h"
#include "mfile.h"

Mfile::Mfile(size_t rec_size)
{
	wwd = rec_size;
	size_t	offbrk = BLOCK_SIZE * wwd;
	ptr = cur = new char[offbrk];
	brk = ptr + offbrk;
	recno = 0;
}

UPTR Mfile::flush()
{
	char*	tmp = 0;
	if (recno > 0) {
	    tmp = new char[recno * wwd];
	    memcpy(tmp, ptr, recno * wwd);
	}
	delete[] ptr; ptr = 0;
	return ((UPTR) tmp);
}

void Mfile::write(const UPTR pi)
{
	if (cur == brk) {
	    brk += BLOCK_SIZE * wwd;
	    size_t	offset = cur - ptr;
	    size_t	offbrk = brk - ptr;
	    char*	tmp = new char[offbrk];
	    memcpy(tmp, ptr, offset);
	    delete ptr;
	    ptr = tmp;
	    cur = ptr + offset;
	    brk = ptr + offbrk;
	}
	memcpy(cur, pi, wwd);
	cur += wwd;
	++recno;
}
