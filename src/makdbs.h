/*****************************************************************************
*
*	Header file to makdbsc
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

#ifndef	_MAKEDBS_H_
#define	_MAKEDBS_H_

#include "seq.h"

static	const	INT	ddelim = SEQ_DELIM + (SEQ_DELIM << 4);

class Makdbs {
	int	molc;
	int	bias;
	int	ceil;
	bool	cridxf;
	SeqDb*	db;
	SEQ_CODE*	defcode;
	size_t	recnbr;
	FILE*	fgrp;
	FILE*	fseq;
	FILE*	fidx;
	FILE*	fent;
	char	str[MAXL];
	char	prv[MAXL];
	bool	halfway;
#if USE_ZLIB
	gzFile	gzseq;
#endif
	int	encode(int c);
	char*	getDbEntry(DbsRec* rec, int idf);
template <typename file_t, typename ofile_t>
	int	convert(file_t fsrc, ofile_t);
template <typename file_t, typename ofile_t>
	int	convert2(file_t fsrc, ofile_t);
	void	initialize(const char* av);
template <typename file_t>
	char*	get_str(file_t fsrc);
template <typename file_t>
	void	skip_till_nl(file_t fsrc);
template <typename file_t, typename ofile_t>
	void	makdbs(file_t fsrc, ofile_t);

public:
	Makdbs(int ac, const char** av, int mlc);
	Makdbs(const char* dbname);
	~Makdbs();
	void	makdbs(const char* fn);
	void	mkidx();
	void	stamp21();
	void	wrtgrp(const char* ps) {
	    long	fpos = 0L;
	    if (fseq)	fpos = ftell(fseq);
#if USE_ZLIB
	    else	fpos = ftell(gzseq);
#endif
	    fprintf(fgrp, "%8ld %u %s\n", fpos, (INT) recnbr, ps);
	}
};

#endif	// _MAKEDBS_H_
