/*****************************************************************************
*
*	Automatic performance of program "perform" for a given set
*	of protein or nucleotide sequences.
*
*	Subcommand: 'e' all; 		'f' single scan; 	'g' group
*	Subcommand: 'i' internal; 	'j' juxtapose		'k' addon
*	Subcommand: '1' onesequence;
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

#include	"seq.h"
#include	"autocomp.h"

InputSeqTest	zero_InputSeqTest = {0, 0, 0, 0, 0, 0, 0, INT_MAX, 0, 0};

int testInputSeq(CalcServer<Seq>* svr, Seq* sqs[], ThQueue<Seq>* q)
{
	InputSeqTest*	tv = (InputSeqTest*) svr->prm;
	Seq*&	sd = sqs[0];
	if (tv->molc && sd->inex.molc != tv->molc) {
	    prompt("%s is mixed type !\n", sd->spath);
	    ++tv->bad;
	} else if (!sd->many) {
	    prompt("%s is empty !\n", sd->spath);
	    ++tv->bad;
	} else {
	    tv->molc = sd->inex.molc;
	    tv->many += sd->many;
	    if (sd->many > tv->maxmany) tv->maxmany = sd->many;
	    tv->space += sd->many * sd->len;
	    if (sd->inex.dels) tv->dels = 1;
	    if (sd->len > tv->maxlen) tv->maxlen = sd->len;
	    if (sd->len < tv->minlen) tv->minlen = sd->len;
	    if (tv->sname) tv->sname->push(sd->sqname());
	    ++tv->num;
	}
	return (OK);
}

