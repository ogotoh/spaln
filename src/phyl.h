/*****************************************************************************
*
*	Collection of headers for calculation of sequence divergence
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

#ifndef  _PHYL_H_
#define  _PHYL_H_

#include "divseq.h"
#include "seq.h"

class mSeq;

#define	EOL	(-1)
#define	LEFT	1
#define	RIGHT	2
#define	CLMWD	80
#define	OTULEN	10
#define	LMARGIN	10
#define	MAXELEM	256
#define	HORLINE	'-'
#define	VERLINE	'|'

enum TreeMet {Def_METHOD, UPG_METHOD, NJ_METHOD, SLINK_METHOD,
		CLINK_METHOD, AB_METHOD, WARD_METHOD};

class	DistTree;

struct ScrRng {FTYPE scr; int left, right;};

static const ScrRng	nullScrRng = {0, INT_MAX, 0};
static const FTYPE	F_UNIF_WEIGHT = 0.0;

struct TOUTMODE {
	TreeMet	method:		4;
	INT	paren:		1;
	INT	tree:		1;
	INT	lroot:		1;
	INT	branch:		1;
	INT	height:		1;
	INT	weight:		1;
	INT	pairwt:		1;
	INT	prnode:		1;
	INT	n_edge:		1;
	INT	balance:	1;
	INT	fixedbase:	1;
	INT	unrooted:	1;
	INT	ddmx:		1;
	INT	path:		1;
	INT	negate:		1;
	INT	dollo:		1;
};

struct LEAF {
	int	tid;
	FTYPE   wheight;
	LEAF*	lnk;
};

struct LIST {
	FTYPE   data;
	int     idx;
};

struct DistMat {
	int	numb;
	DistCal	realign;
	Strlist*	sname;
	double	bias;
	FTYPE*	dist;
	ScrRng*	scrrng;
	DistMat(FILE* fd);
	DistMat(int argc, const char** argv, 
	    const char* catalog = 0, DistCal realin = Composition);
	DistMat(Seq** seqs, int nn);
	DistMat(mSeq** seqs, int nn, DistCal realin = DynAln);
	~DistMat();
};

struct Tnode {
	Tnode*  left;
	Tnode*  right;
	Tnode*  parent;
	int     tid;
	int     ndesc;
	FTYPE   height;
	FTYPE   length;
	char*	tname;
	Tnode() {
	    left = right = parent = 0;
	    height = length = tid = ndesc = 0;
	}
	~Tnode() {}
	bool    isleaf() {return !(left || right);}
	bool	isroot() {return (!parent);}
	int	ndescend() {return (ndesc);}
};

struct  Knode {
	Knode*  left;
	Knode*  right;
	Knode*  parent;
	int     tid;
	int     ndesc;
	FTYPE   height;
	FTYPE   length;
	union	{FTYPE vol; int gain;};
	union 	{FTYPE cur; int loss;};
	FTYPE	res;
	union	{FTYPE ros; int stat;};
	Knode() {
	    left = right = parent = 0;
	    height = length = tid = ndesc = 0;
	    vol = cur = ros = res = 0;
	}
	~Knode() {}
	bool    isleaf() {return !(left || right);}
	bool	isroot() {return (!parent);}
	int	ndescend() {return (ndesc);}
	double	calres();
	Knode*	findroot(double brl);
	Knode*	findcenter();
	int	teachparent();
	void    kirchhof();
	Knode*	balanced();
};

template <class node_t>
class Btree {
public:
	int	members;	// mumber of leaves
	char**	mname;
	char*	mnbuf;
	node_t*	root;
	node_t*	lead;		// array top
	void	GetNHtree(FILE* fd);
	void	GetTXtree(FILE* fd);
	void	fill_tname();
	char*	lname(int i) {return mname[i];}
	char*	lname(node_t* node) {return mname[node->tid];}
	Btree(const char* fname);
	~Btree();
};

template <class node_t>
Btree<node_t>::Btree(const char* fname)
{
        root = lead = 0;
        mname = 0;
        mnbuf = 0;
        members = 0;
        FILE*   fd = fopen(fname, "r");
        if (!fd) return;
        int     c;
        while (isspace(c = getc(fd))) ;
        if (ungetc(c, fd) == '(')
                GetNHtree(fd);
        else    GetTXtree(fd);
        fclose(fd);
}

template <class node_t>
Btree<node_t>::~Btree()
{
	delete[] mnbuf;
	delete[] mname;
	delete[] lead;
}

template <class node_t>
void Btree<node_t>::fill_tname()
{
	node_t*	nd = lead;
	for (int i = 0; i < members; ++i) 
	    (nd++)->tname = mname[i];
}

template <class node_t>
void Btree<node_t>::GetNHtree(FILE* fd)
{
	int	cc, nr = 0, nl = 0, nc = 0;
	int	state = 1;
	int	maxlen = 0;
	int	lenlen = 0;
	int	mnsize = 0;

//	First read

	while ((cc = getc(fd)) != EOF && cc != ';') {
	    if (iscntrl(cc) || isspace(cc)) continue;
	    switch (cc) {
		case '(': ++nl; break;
		case ')': ++nr; if (lenlen > maxlen) maxlen = lenlen; break;
		case ',': ++nc; if (lenlen > maxlen) maxlen = lenlen; state = 1; break;
		case ':': state = 2; lenlen = 0; break;
		default: 
		    if (state == 2)	++lenlen;
		    else	++mnsize;
		    break;
	    }
	}

	if (nl != nr || (nr != nc && nr + 1 != nc) || !nl) {
	    delete[] mnbuf;
	    fprintf(stderr, "Unbalanced: %d( %d, %d)\n", nl, nc, nr);
	    return;
	}
	members = ++nc;
	mname = new char*[members + 1];
	mnbuf = new char[mnsize + members];
	char*	lenbuf = maxlen? new char[maxlen + 1]: 0;
	lead = new node_t[2 * members - 1];
	vclear(lead, 2 * members - 1);
	root = lead;
	Stack<node_t*>   cstk;

//      Parse text and construct tree structure
//      nl: leaf no.; nc: internal node no.;

	rewind(fd);
	FTYPE   hl, hr;
	char*	ps = mnbuf;
	char*	pl = lenbuf;
	state = nl = 0;
	while ((cc = getc(fd)) != EOF && cc != ';') {
	  if (iscntrl(cc) || isspace(cc)) continue;
	  switch (cc) {
	    case '(': break;
	    case ',': 
		if (state == 1) *ps++ = '\0';
		if (pl) {
		    *pl = '\0';
		    root->length = atof(lenbuf);
		} else 
		    root->length = 1.;
		state = 0;
		break;
	    case ':': 
		if (state == 1) *ps++ = '\0';
		state = 2; pl = lenbuf; 
		break;
	    case ')':   // Interna node
		if (state == 1) *ps++ = '\0';
		if (pl) {
		    *pl = '\0';
		    root->length = atof(lenbuf);
		} else 
		    root->length = 1.;
		root = lead + nc;
		root->right = cstk.pop();
		root->left = cstk.pop();
		root->tid = nc++;
		hl = root->left->height + root->left->length;
		hr = root->right->height + root->right->length;
		root->height = max(hl, hr);
		root->left->parent = root->right->parent = root;
		cstk.push(root);
		state = 0;
		break;
	    default:
		if (state == 0) {	// leaf
		    root = lead + nl;
		    mname[nl] = ps;
		    root->tid = nl++;
		    root->height = 0;
		    root->left = root->right = 0;
		    cstk.push(root);
		    state = 1;
		    *ps++ = cc;
		} else if (state == 1) {
		    *ps++ = cc;
		} else if (pl) {
		    *pl++ = cc;
		}
		break;
	  }
	}
	delete[] lenbuf;
	root->length = 0;
}

template <class node_t>
void Btree<node_t>::GetTXtree(FILE* fd)
{
	char    str[MAXL];
	int	nmsize = 0;

	while (fgets(str, MAXL, fd)) {
	    char*       wk = strchr(str, ':');
	    if (wk) {
		++members;
		*wk = '\0';
		char*	pq = str;
		char*	ps = car(pq);
		nmsize += (pq - ps);
	    }
	}
	lead = new node_t[2 * members];
	mname = new char*[members + 1];
	mnbuf = new char[nmsize + members];
	vclear(lead, 2 * members);

	rewind(fd);
	Stack<node_t*>   ctck;
	int	nl = 0;
	int	nc = members;
	char*	ps = mnbuf;
	node_t*	leaf = 0;
	node_t*	node = 0;
	FTYPE	hh = 0.;
	while (fgets(str, MAXL, fd)) {
	    if (isBlankLine(str)) continue;
	    char*	wk = strchr(str, ':');
	    if (wk) {
		*wk = '\0';
		char*	pq = str;
		char*	ss = car(pq);
	        mname[nl] = ps;
		while ((*ps++ = *ss++)) ;
		leaf = lead + nl;
		leaf->tid = nl++;
		leaf->left = leaf->right = 0;
		leaf->height = 0;
	    } else {
		node = lead + nc;
		node->tid = nc++;
		node->height = atof(str);
		if (node->height > hh) {
		    root = node;
		    hh = node->height;
		}
		while (node_t* last = ctck.pop()) {
		    if (last->height > node->height) {
			ctck.push(last);
			break;
		    }
		    last->right = leaf;
		    leaf->parent = last;
		    leaf->length = last->height - leaf->height;
		    leaf = last;
		}
		leaf->parent = node;
		node->left = leaf;
		leaf->length = node->height - leaf->height;
		ctck.push(node);
	    }
	}
	while (node_t* last = ctck.pop()) {
	    last->right = leaf;
	    leaf->parent = last;
	    leaf->length = last->height - leaf->height;
	    leaf = last;
	}
//      Move the root to the last
	if (!root) root = node;
	else if (root != node) {
	    if (node->parent) {
		if (node->parent->left == node) node->parent->left = root;
		if (node->parent->right == node) node->parent->right = root;
	    }
	    if (root->parent) {
		if (root->parent->left == root) root->parent->left = node;
		if (root->parent->right == root) root->parent->right = node;
	    }
	    swap(*node, *root);
	    swap(node->tid, root->tid);
	    node->left->parent = node->right->parent = root;
	    root->left->parent = root->right->parent = node;
	    if (root->left  == root) root->left  = node;
	    if (root->right == root) root->right = node;
	    swap(root, node);
	}
}

template <class node_t>
void fill_ndecent(node_t* nd)
{
	if (nd->isleaf()) nd->ndesc = 1;
	else {
	    fill_ndecent<node_t>(nd->left);
	    fill_ndecent<node_t>(nd->right);
	    nd->ndesc = nd->left->ndesc + nd->right->ndesc;
	}
}

class Ktree {
protected:
	int	members;	// mumber of leaves
	LEAF*	leaves;		// leaves
	int	leaf_no;	// leaf identifier
#if USE_WEIGHT
	FTYPE*	pwt;		// pair weight temporary variable
	FTYPE	bwt;		// basal weight
#endif
	SigII*	sgi;
	LEAF*	repairwt(Knode* node);
	void	plant(DistTree& dt);
	void	dollo_1st(Knode* node, int j);
	void	dollo_2nd(Knode* node, bool peist);
public:
	Knode*	lead;	// array top
	Knode*	root;	// root
	Ktree(int mem = 0, Knode* led = 0);
	Ktree(DistMat* dmat, TreeMet tmethod = Def_METHOD, Knode* led = 0);
	Ktree(Seq* sd, Subset* ss = 0, TreeMet 
		tmethod = Def_METHOD, Knode* led = 0);
	~Ktree() {delete[] lead; delete[] leaves;}
	Knode*	skimtree() {lead = 0; return (root);}
#if USE_WEIGHT
	FTYPE*	calcwt(FTYPE* wt = 0);
	FTYPE* 	calcpw(FTYPE* wt = 0);
	LEAF*	pairwt(Knode* node, double ros);
	FTYPE*	recalcpw(Knode* node = 0, int num = 0);
	FTYPE	sum_of_pairwt(Knode* node = 0);
#endif
	void	mkleaves() {if (!leaves) leaves = new LEAF[members];}
	void	dollo(Seq* sd);
};

class DistTree : public Ktree {
protected:
	FTYPE*	dist;
	int*	row;
	int*	nnbr;
	bool	newdmat;
	char**	mname;
	char*	mnbuf;
	int	ntimes;
	double	lwhi;
	FTYPE	dminidx(int& i, int members);
	int	dminrow(int i, int members);
	void	lowesthi(Knode* node, double hi);
	double	recalhi(Knode* node, double hi);
public:
	FTYPE*	mtrx;
	DistTree(DistMat* dmat, TreeMet tmethod = Def_METHOD, Knode* led = 0);
	DistTree(Seq* sd, Subset* ss = 0, TreeMet tmethod = Def_METHOD,
	    Knode* led = 0, bool nm = false);
	DistTree(const char* user_tree);
	~DistTree() {
	    delete[] nnbr; delete[] row;  delete[] mtrx; delete[] mname;
	    delete[] mnbuf;
	    if (newdmat) {delete[] dist;}
	}
	FTYPE&	distance(int i, int j) {return (dist[elem(i, j)]);}
	Knode*	upg_method();
	Knode*	nj_method();
	void	maketree(TreeMet met = Def_METHOD);
	void	putdmx(FILE* fo = stdout);
	LEAF*	estimatedist(Knode* node);
	char*	lname(int i) {return mname[i];}
	char*	lname(Knode* node) {return mname[node->tid];}
};

class PpPrm : public DistTree {
	char	str[MAXL];
	FILE*	fo;
	int	minelem, lastelem, datp;
	int	otulen, clmpos;
	double	factor;
	double	bias;
	double	maxhi;
	LIST	list[MAXELEM];
	char	sfrmt[10];
	char	ffrmt[10];
	char*	linebuf;
	char*	baseline;
	void	fillstr(int cha, int from, int to);
	void	printcur(Knode* node);
//	Gnm2tab*	g2t;
public:
	PpPrm(FILE* fd, Seq* sd, Subset* ss, TreeMet meth, double bs = 0);
	PpPrm(FILE* fd, DistMat* dmat, TreeMet meth, double bs = 0);
	PpPrm(FILE* fd, const char* user_tree, double bs = 0);
	~PpPrm() {delete[] linebuf;}
	void	prepare();
	void	add(double data);
	void	del(double data);
	double	firstdat();
	double	nextdat();
	void	putwt();
	void	putpw();
	void	putheight(Knode* node);
	void	printnode(Knode* node, double father, int lr);
	void	printpar(Knode* node, bool newpar);
	void	wdadjust();
	void	superimp(int a, int z, double len);
	void	superimp(int a, int z, int gain, int loss);
	void	superimpnn(int a, int z, int nn);
	void	puttree();
	LEAF*	tracepath(Knode* node);
	int	normal(double d, int lroot) {return ((int) ((lroot? 
		(maxhi - d): (d - bias)) * factor));}
};

extern	TOUTMODE	treemode;

extern	void	divseq(FSTAT* stat, Seq* sd, int* group1, int* group2);
extern	FTYPE*	calcdist(mSeq** sqs, int nn, DistCal realn);
extern	FTYPE*	calcdist(Seq* sd, Subset* ss = 0);
extern	FTYPE*	calcdistsum(Seq* sd, Subset* ss = 0, FTYPE* dist = 0);
extern	FTYPE*	calcdist_i(Seq* sd, int k, Subset* ss = 0, FTYPE* dist = 0);
extern	TreeMet	setphyl(int method);
extern	void	setpwfact(double wf);
#if USE_WEIGHT
extern	void	calcweight(Seq* sd);
#endif

#endif	// _PHYL_H_
