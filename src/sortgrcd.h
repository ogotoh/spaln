/*******************************************************************************
*
*	sortgrcd.h version 2.2
*
*	Recover the output of spaln -O12
*	Filter the output by several criteria
*	Sort the records on the genomic positions
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

#include "aln.h"

#if M_THREAD
#include <pthread.h>
#include <unistd.h>
#include <sched.h>
#endif

#if MONITOR
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <sys/timeb.h>
#endif

enum InOrder {INPUT_ODR, ALPHABETIC, ABUNDANCE};

static	const	int	max_no_exon = 1024;
static	const	int	MAXRECD = 16 * 1024 * 1024;
static	const	char*	RGB_RED = "255,0,0";
static	const	char*	RGB_BLUE = "0,255,255";

struct GRFn {
	GeneRecord	gr;	// G-record No.
	INT		fn;	// file No.
};

struct GERecN  {
	size_t	grn;		// G-record No.
	size_t	ern;		// E-record No.
	INT	gcn;		// cumulative
	INT	ecn;		// cumulative 
	Strlist* sname;
};

struct ChrInf {
	int	key;
	INT	gcount;
	INT	ecount;
};

typedef int (*CMPCIF)(ChrInf*, ChrInf*);

struct Chash {
	INT	size;
	INT	nelm;
	ChrInf*	hash;
	Chash(INT n);
	~Chash() {delete[] hash;}
	void	clear() {vclear(hash, size);}
	ChrInf*	resize(int wkey);
	ChrInf*	chashv(int wkey, bool incr);
	ChrInf*	scream();
};

struct ExnInf {
	int	left;
	int	right;
	INT	count;
	INT	exonID;
};
 
class Ehash {
	INT	size;
	ExnInf*	hash;
public:
	Ehash(INT n);
	~Ehash() {delete[] hash;}
	void	clear() {memset(hash, '\0', size * sizeof(ExnInf));}
	ExnInf*	ehashv(int left, int right);
	ExnInf*	resize(int left, int right);
};

struct IntronInf {
	int     left;
	int     right;
	INT     count;
	int	Gleft;
	int	Gright;
	int	Rleft;
	int	Rright;
	int	mch;
	int	mmc;
	int	unp;
	int	ilen;
	int	Rid;
	char	intends[8];
};

class Ihash {
	INT	size;
	IntronInf*	hash;
public:
	Ihash(INT n);
	~Ihash() {delete[] hash;}
	void	clear() {vclear(hash, size);}
	IntronInf* ihashv(int left, int right, ExonRecord* ewrk, 
		ExonRecord* fwrk, int mch, int mmc, int unp, 
		char* intends, int rid);
	IntronInf* resize(int left, int right, ExonRecord* ewrk,
		ExonRecord* fwrk, int mch, int mmc, int unp, 
		char* intends, int rid);
	INT	sortihash();
	IntronInf*	begin() {return hash;}
};

struct FiltParam {
	int     bmmc;
	int     bunp;
	int     ncan;
	int     Bmmc;
	int     Bunp;
	int     ng;
	float   Gscore;
	float   Pmatch;
	float   Pcover;
};

class Sortgrcd {
	int	argc;		// total # of .grd files
	INT	ngrcd;		// total # of G-record
	INT	nercd;		// total # of E-record
	INT	nchr;		// total # of chromosomes
	GERecN*	nrcd;		// [argc]
	GRFn*	grcd;		// [ngrcd]
	ExonRecord*	ercd;	// 
	ChrInf*	chrlist;	// [ncrh];
const	char*	grdname;	// input file name
	void	assort_by_chr(Chash* hh);
template <typename file_t>
	int	read_chr_rec(file_t fe, ExonRecord*& ewrk, 
			GRFn* frcd, INT grn, INT fn);
public:
	void	printGrcd();
	void	readGrcd(int ac, const char** av, INT hashsize);
	ExonRecord* ReadRcd(int ac, const char** av);
	ExonRecord* ReadChrRcd(int ac, const char** av, INT nercd, GRFn* frcd, INT grn);
	void	print_cds(GeneRecord* gwrk, RANGE* exon, const char* rname);
	void	print_bed(GeneRecord* gwrk, RANGE* exon, const char* Rname);
	void	Exonform(ExonRecord* ercd, GRFn* grfn = 0, INT gcn = 0);
	void	Intronform(ExonRecord* ercd, GRFn* grfn = 0, INT gcn = 0);
	void	Gff3form(ExonRecord* ercd, GRFn* grfn = 0, INT gcn = 0);
	void	Cdsform(ExonRecord* ercd, GRFn* grfn = 0, INT gcn = 0);
	GRFn*	begin() {return (grcd);}
	INT	nGrcd() {return (ngrcd);}
	INT	nErcd() {return (nercd);}
	ChrInf*	chrbegin() {return (chrlist);}
	ChrInf*	chrend() {return (chrlist + nchr);}
	Sortgrcd(int ac, const char** av);
	~Sortgrcd();
};

#if M_THREAD

struct thread_arg_t {
	int	cpuid;
	GRFn*	grcd;
	int*	occr;
	INT	nchr;
};

#endif	// M_THREAD

