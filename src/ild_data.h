/*****************************************************************************
*
*	Generate xxx.ild[.dat] from xxx.eij or
*	Convert between text and binary forms
*
*	Osamu Gotoh, ph.D.      (-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.      (2001-2023)
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

#ifndef	_ILD_DATA_H_
#define	_ILD_DATA_H_

#include "stdtype.h"

struct	ILD_HEADER {
	int	ndsize, nsamples;
};

using	var_t = SHORT;
using	Kvp = KVpair<int, var_t>;

struct	sLfp {
	int	len;
	var_t	frq;
};

extern	INT	min_samples;

inline	void	set_min_samples(const INT& n) {min_samples = n;}
inline	int	complen(const sLfp* a, const sLfp* b)
	{return (a->len - b->len);}

class	Ild_data {
	ILD_HEADER	header;
template <typename file_t>
	void read_ild_for_map(file_t fd);
template <typename file_t>
	int	read_text(file_t fd);
template <typename file_t>
	int	read_binary(file_t fd, const char* fname);
template <typename file_t>
	void	write_text(file_t ofd);
template <typename file_t>
	void	write_binary(file_t ofd, const char* oname);
public:
	sLfp*	lenfrq = 0;
	int	size() {return (header.ndsize);}
	int	samples() {return (header.nsamples);}
	int	from_eij(const char* fname);
	int	from_file(const char* fname);
	void	to_file(const char* oname, const int& text, const bool& gzip);
	Ild_data(const char* fname = 0) {
	    header.ndsize = header.nsamples = 0;
	    if (fname && !from_file(fname)) fatal(not_found, fname);
	}
	~Ild_data() {delete[] lenfrq;}
};

template <typename file_t>
void Ild_data::read_ild_for_map(file_t fd)
{
	char	str[MAXL];
	Dhash<int, var_t>	lfhash(28000, 0);

	while (fgets(str, MAXL, fd)) {
	    if (*str == '#' || *str == '\0') continue;
	    Strlist	slt(str, stddelim);
const	    int	nitem = slt.size();
	    int	len = atoi(slt[nitem - 2]);
	    lfhash.incr(len);
	    ++header.nsamples;
	}
	lenfrq = (sLfp*) lfhash.squeeze(&header.ndsize);
	qsort((UPTR) lenfrq, header.ndsize, sizeof(sLfp), (CMPF) complen);
}

template <typename file_t>
int Ild_data::read_text(file_t fd)
{
	char	str[MAXL];
	Mfile	mfd(sizeof(sLfp));

	while (fgets(str, MAXL, fd)) {
	    if (*str == '#' || *str == '\0') continue;
	    Strlist	slt(str, stddelim);
	    if (!isdigit(slt[0][0])) continue;
	    int	len = atoi(slt[0]);
	    var_t	num = atoi(slt[1]);
	    sLfp	lfq = {len, num};
	    mfd.write(&lfq);
	    ++header.ndsize;
	    header.nsamples += num;
	}
	header.ndsize = mfd.size();
	lenfrq = (sLfp*) mfd.flush();
	return (header.ndsize);
}

template <typename file_t>
int Ild_data::read_binary(file_t fd, const char* fname)
{
	if (fread(&header, sizeof(ILD_HEADER), 1, fd) != 1)
	    fatal(read_error, fname);
	lenfrq = new sLfp[header.ndsize];
	if (fread(lenfrq, sizeof(sLfp), header.ndsize, fd) 
		!= (size_t) header.ndsize)
	    fatal(read_error, fname);
	return (header.ndsize);
}

template <typename file_t>
void Ild_data::write_text(file_t ofd)
{
	fputs("ilen\tfreq\n", ofd);
	sLfp*	w = lenfrq;
	sLfp*	t = w + header.ndsize;
	char	str[MAXL];
	for ( ; w < t; ++w) {
	    sprintf(str, "%d\t%d\n", w->len, w->frq);
	    fputs(str, ofd);
	}
}

template <typename file_t>
void Ild_data::write_binary(file_t ofd, const char* oname)
{
	if (fwrite(&header, sizeof(ILD_HEADER), 1, ofd) != 1 ||
	fwrite(lenfrq, sizeof(sLfp), header.ndsize, ofd) != 
		size_t(header.ndsize))
	    fatal(write_error, oname);
}

#endif	// _ILD_DATA_H_
