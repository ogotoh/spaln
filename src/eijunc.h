/*****************************************************************************
*
*	header to read exon/intron sequences indicated in xxx.eij
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

#ifndef _EIJUNC_H_
#define _EIJUNC_H_

#include <vector>
#include "seq.h"

enum	Eijmode {DONOR, ACCPTR, JOIN, INTRON, EXON, BRANCH};

static	const	int	Ref_Column = 7;
static	const	int	def_wing = 25;
static	const	int	def_brch = 120;

class EiJuncSeq {
	Eijmode	eori = INTRON;
	DbsDt*	dbs = 0;
	Seq*	eijseq = 0;
	FILE*	eijfd = 0;
	char	str[MAXL];
	int	wing = 0;
	int	wms1 = 0;
public:
	int	width = 0;
	bool	new_entry = false;
	char	genspc[MAXL];
	char	transcript[MAXL];
	Seq*	nextseq();
	void	typeseq(FILE* fo) {eijseq->typeseq(fo);}
	void	reset() {
	    new_entry = false;
	    if (eijfd) rewind(eijfd);
	    *str = *transcript = '\0';
	}
	char*	goahead() {
	    if (eori == EXON && new_entry) return (str);
	    return (fgets(str, MAXL, eijfd));
	}
	EiJuncSeq(Eijmode mode, const char* eij, 
		const char* genome = 0, int w = 0);
	~EiJuncSeq() {
	    delete eijseq; 
	    if (eijfd) fclose(eijfd);
	    delete dbs;
	}
};

#endif // _EIJUNC_H_
