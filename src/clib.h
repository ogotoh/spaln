/*****************************************************************************
*
*	Collection of functions for general use
*
*	vclear	vcopy	vset	vreverse	vmax	vmin
*	max3	min3	max4	min4
*	ipower	
*	next_wd	wordcmp	isBlankLine  
*	car	cdr	prefix	chop	replace
*	queue	stack	fbisrch	insort
*	comb2	comb3	decomb
*	Strlist
*
*	Osamu Gotoh, ph.D.	(-2001)
*	Saitama Cancer Center Research Institute
*	818 Komuro, Ina-machi, Saitama 362-0806, Japan
*
*	Osamu Gotoh, Ph.D.	(2001-)
*	National Institute of Advanced Industrial Science and Technology
*	Computational Biology Research Center (CBRC)
*	2-4-7 Aomi, Koutou-ku, Tokyo 135-0064, Japan
*
*	Osamu Gotoh, Ph.D.      (2003-2012)
*	Department of Intelligence Science and Technology
*	Graduate School of Informatics, Kyoto University
*	Yoshida Honmachi, Sakyo-ku, Kyoto 606-8501, Japan
*
*	Copyright(c) Osamu Gotoh <<o.gotoh@aist.go.jp>>
*****************************************************************************/
#ifndef  _CLIB_
#define  _CLIB_

#include <math.h>
#include <algorithm>
#include "adddef.h"

struct ALGMODE {
	INT	nsa :   4;	// no/single/all alignments
	INT	alg :   4;	// rigor in group sequnece
	INT	bnd :   1;	// banded dp space
	INT	mlt :   3;	// multiple output
	INT	aut :	1;	// automatic param set
	INT	lcl :   5;	// local or global
	INT	lsg :   1;	// # of line pieces > 1 ?
	INT	mns :   2;	// minus strand as well
	INT	thr :   1;	// cutoff if < thr
	INT	rng :   1;	// output range
	INT	qck :   2;	// quick calculation
	INT	lvl :	2;	// wilip initial level
	INT	blk :	1;	// rapid genome scan
	INT	any :	2;	// accept non-consensus spj
	INT	crs :	2;	// cross species comparison
	INT	slv :	1;	// salvage all positive blocks
	INT	dim:	1;	// database seq in memory
};

using std::max ;
using std::min ;
using std::swap ;

extern	int	thread_num;
extern	int	cpu_num;
extern	int	max_queue_num;
extern	INT	supprime(INT n);
extern	INT	str2uint(const char* str);
extern	ALGMODE	algmode;

static	const	int	KILO = 1024;
static	const	int	MEGA = 1048576;
static	const	long	GIGA = 1073741824;
static	const	double	PI = 3.14159265358979323846;
static	const	int	DefStackDepth = 128;
static	const	int	defsunit = 128;
static	const	int	defpunit = 16;
static	const	INT	SecondHS = 8;
static	const	float	DefHashFact = 1.2;
static	const	double	epsilon = 1e-6;
static	const	int	HashovLS = 12;
static	const	int	def_un_def = (INT_MIN / 8 * 7);
static	const	float	FACT_QUEUE = 1.5;

template <typename X> X* vcopy(X* dst, const X* src, int n)
{
	if (!dst) dst = new X[n];
	X*	head = dst;
	while (n-- > 0) *dst++ = *src++;
	return (head);
}

template <typename X> X* vset(X* dst, const X& val, size_t n)
{
	if (n == 0) return(dst);
	if (!dst) dst = new X[n];
	X*	w = dst + n;
	while (--w >= dst) *w = val;
	return (dst);
}

template  <typename X> inline void vclear(X* ary, const int n = 1)
{
	memset(ary, '\0', n * sizeof(X));
}

template <typename X> X* vmax(const X* array, int n)
{
	const X* temp = array;
	while (--n > 0) if (*++array > *temp) temp = array;
	return ((X*) temp);
}

template <typename X> X* vmin(const X* array, int n)
{
	const X* temp = array;
	while (--n > 0) if (*++array < *temp) temp = array;
	return ((X*) temp);
}

template <typename X> X* vreverse(X* array, int n)
{
	X*	f = array;
	X*	b = array + n - 1;
	while (f < b) swap(*f++, *b--);
	return (array);
}

template <typename X> X vavsd(X& sd, X* array, int n)
{
	double	av = sd = 0;
	if (n < 1) return ((X) av);
	if (n == 1) return (*array);
	for (int i = 0; i < n; ++i) av += array[i];
	av /= n;
	double	vr = 0;
	for (int i = 0; i < n; ++i) {
	    double	x = array[i] - av;
	    vr += x * x;
	}
	sd = (X) sqrt(vr / (n - 1));
	return ((X) av);
}

template <typename X> inline X max3(X x, X y, X z)
{
	if (y > x) x = y;
	return (max(x, z));
}

template <typename X> inline X min3(X x, X y, X z)
{
	if (y < x) x = y;
	return (min(x, z));
}

template <typename X> inline X max4(X x, X y, X z, X w)
{
	if (y > x) x = y;
	if (w > z) z = w;
	return (max(x, z));
}

template <typename X> inline X min4(X x, X y, X z, X w)
{
	if (y < x) x = y;
	if (w < z) z = w;
	return (min(x, z));
}

// key-value pair used by Dhash

template <class key_t, class val_t>
struct KVpair {key_t key; val_t val;};

// double hash

template <class key_t, class val_t>
class Dhash {
protected:
	INT	size1;
	INT	size2;
	KVpair<key_t, val_t>*	hash;
	KVpair<key_t, val_t>*	hz;
	val_t	un_def;
public:
	Dhash(int n = 0, val_t ud = 0, float hf = DefHashFact);
	~Dhash() {delete[] hash;}
	KVpair<key_t, val_t>*	begin() {return hash;}
	KVpair<key_t, val_t>*	end()	{return hz;}
	KVpair<key_t, val_t>*	map(key_t ky, bool record = true);
	KVpair<key_t, val_t>*	find(key_t ky) {
	    return map(ky, false);
	}
	KVpair<key_t, val_t>*	assign(key_t ky, val_t vl = 1, bool first = false) {
	    KVpair<key_t, val_t>*	sh = map(ky);
	    if (!first || sh->val == un_def) sh->val = vl;
	    return (sh);
	}
	KVpair<key_t, val_t>*	incr(key_t ky, val_t vl = 1){
	    KVpair<key_t, val_t>*	sh = map(ky);
	    if (sh->val == un_def) sh->val = 0;
	    sh->val += vl;
	    return (sh);
	}
	int	count() const;
	KVpair<key_t, val_t>*	press(int* n, val_t* sum = 0);
	void	clear() {
	    KVpair<key_t, val_t> zero = {0, un_def};
	    if (un_def)	vset(hash, zero, size1);
	    else	vclear(hash, size1);
	}
	void	resize(int s = 0);
	void	remove(key_t ky) {
	    KVpair<key_t, val_t>*	sh = map(ky);
	    sh->val = un_def;
	}
	KVpair<key_t, val_t>*	squeeze(int* n) {
	    KVpair<key_t, val_t>* tmp = press(n);
	    hash = 0;
	    return (tmp);
	}
	val_t	undef() {return (un_def);}
};

template <class key_t, class val_t>
Dhash<key_t, val_t>::Dhash(int n, val_t ud, float hf) : un_def(ud)
{
	if (n == 0) {
	    size1 = size2 = 0;
	    hash = hz = 0;
	} else {
	    size1 = supprime(int(hf * n));
	    hash = new KVpair<key_t, val_t>[size1];
	    hz = hash + size1;
	    clear();
	    for (size2 = SecondHS; size2 > size1; size2 /= 2) ;
	}
}

template <class key_t, class val_t>
KVpair<key_t, val_t>* Dhash<key_t, val_t>::map(key_t key, bool record)
{
	INT	v = key % size1;
	INT	u = size2 - key % size2;
	INT	v0 = v;
	KVpair<key_t, val_t>*	sh = hash + v;

	while (sh->val != un_def && sh->key != key) {
	    v = (v + u) % size1;
	    if (v == v0) resize();
	    sh = hash + v;
	}
	if (sh->val == un_def) {
	    if (record) sh->key = key;
	    else	sh = 0;
	}
	return (sh);
}

template <class key_t, class val_t>
int	Dhash<key_t, val_t>::count() const
{
	int	c = 0;
	for (KVpair<key_t, val_t>* sh = hash; sh < hz; ++sh)
	    if (sh->val != un_def) ++c;
	return (c);
}

template <class key_t, class val_t>
KVpair<key_t, val_t>* Dhash<key_t, val_t>::press(int* n, val_t* sum)
{
	if (sum) *sum = 0;
	KVpair<key_t, val_t>*	rh = hash;
	for (KVpair<key_t, val_t>* sh = hash; sh < hz; ++sh) {
	    if (sh->val != un_def) {
		*rh++ = *sh;
		if (sum) *sum += sh->val;
	    }
	}
	hz = rh;
	*n = rh - hash;
	return (hash);
}

template <class key_t, class val_t>
void Dhash<key_t, val_t>::resize(int s)
{
	KVpair<key_t, val_t>*	hold = hash;
	KVpair<key_t, val_t>*	zold = hz;
	size1 = supprime(s? s: 2 * size1);
	if (!size2)
	    for (size2 = SecondHS; size2 > size1; size2 /= 2) ;
	hash = new KVpair<key_t, val_t>[size1];
	hz = hash + size1;
	clear();
	for (KVpair<key_t, val_t>* hw = hold; hw < zold; ++hw)
	    if (hw->val != un_def) assign(hw->key, hw->val);
	delete[] hold;
}

// simple queue

template <typename X>
class Queue {
	int	qsize;
	int	qp;
	X*	queue;
public:
	void	clear(X iv = 0) {
	    qp = 0;
	    vset(queue, iv, qsize);
	}
	Queue(int sz, X iv = 0) : qsize(sz) {
	    queue = new X[sz];
	    clear(iv);
	}
	Queue(Queue& q) {
	    qsize = q.qsize;
	    qp = q.qp;
	    queue = new X[qsize];
	    vcopy(queue, q.queue, qsize);
	}
	~Queue() {delete[] queue;}
	X	shift(X val) {
	    swap(val, queue[qp]);
	    if (++qp == qsize) qp = 0;
	    return (val);
	}
	X	unshift(X val) {
	    swap(val, queue[qp]);
	    if (qp) --qp;
	    else qp = qsize - 1;
	    return (val);
	}
	X	oldest() const {return queue[qp];}
};

// priority queue simple version without heap

template <typename val_t>
class PrQueue {
private:
	val_t*	data;
	int	capacity;
	int	front;
	bool	dec_order;	// default is ascending order
	bool	replace;	// may replace older value
public:				// initial data size, ascending order, never replace data
	PrQueue(val_t* _data, int _size, int fr = 0, bool maxi = false, bool rep = false);
	~PrQueue() {}
	void	downheap(int k, int kk = 0);
	void	upheap(int k);
	val_t	shift();
	val_t	shift(val_t& x);
	void	put(const val_t& x, int p = -1);
	val_t	gettop() {return data[0];}
	val_t&	operator[](int i) {return (data[i]);}
	bool	empty() const{return front == 0;}
	int	size() const {return front;}
	void	reset() {front = 0;}
	void	hsort();
	int	find(const val_t& x);
	bool	lt(const val_t& a, const val_t& b) const {
		     return (dec_order? b < a: a < b);
		}
};

template <typename val_t>
PrQueue<val_t>::PrQueue(val_t* _data, int _size, int fr, bool maxi, bool rep)
	: data(_data), capacity(_size), front(fr), dec_order(maxi), replace(rep)
{
	for (int k = front / 2; --k >= 0; ) downheap(k);
}

template <typename val_t>
void PrQueue<val_t>::downheap(int k, int kmax)
{
	val_t	v = data[k];
	if (kmax == 0) kmax = front;
	while (k < kmax / 2) {
	    int	l = 2 * k + 1;
	    int	r = l + 1;
	    if (r < kmax && lt(data[r], data[l])) ++l;
	    if (!lt(data[l], v)) break;
	    data[k] = data[l];
	    k = l;
	}
	data[k] = v;
}

template <typename val_t>
void PrQueue<val_t>::upheap(int k)
{
	val_t	v = data[k];
	int	h = (k - 1) / 2;
	while (k && lt(v, data[h])) {
	    data[k] = data[h];
	    k = h;
	    h = (h - 1) / 2;
	}
	data[k] = v;
}

template <typename val_t>
val_t PrQueue<val_t>::shift()
{
	val_t v = data[0];
	if (front > 0) {
	    data[0] = data[--front];
	    downheap(0);
	}
	return (v);
}

template <typename val_t>
val_t PrQueue<val_t>::shift(val_t& x)
{
	val_t v = data[0];
	if (front == capacity) {
	    data[0] = data[--front];
	    downheap(0);
	}
	data[front] = x;
	upheap(front++);
	return (v);
}

template <typename val_t>
void PrQueue<val_t>::put(const val_t& x, int p)
{
	if (p < 0) {			// new member
	    if (front < capacity) {
		data[front] = x;	// add
		upheap(front++);
		return;
	    } else	p = 0;
	}
	if (lt(data[p], x)) {		// exam
	    data[p] = x;		// replace
	    downheap(p);
	}

}

template <typename val_t>	// default key finder
int PrQueue<val_t>::find(const val_t& x)
{
        for (int i = 0; i < front; ++i)
            if (x == data[i]) return (i);
        return (-1);
}

template <typename val_t>
void PrQueue<val_t>::hsort()	// heap sort data
{
	for (int k = front; --k > 0; ) {
	    swap(data[0], data[k]);
	    downheap(0, k);
	}
}

// priority queue without hash

template <typename val_t>
class PrQueue_wh {
private:
	val_t*	data;
	int	capacity;
	int	front;
	bool	dec_order;	// default is ascending order
	bool	replace;	// may replace older value
	Dhash<int, int>*	hpos;	// heap position of key
public:				// initial data size, ascending order, never replace data
	PrQueue_wh(val_t* _data, int _size, int fr = 0, bool maxi = false, bool rep = false);
	~PrQueue_wh() {delete hpos;}
	void	downheap(int k, int kk = 0);
	void	upheap(int k);
	val_t	shift();
	val_t	shift(val_t& x);
	void	put(const val_t& x, int p = -1);
	val_t	gettop() {return data[0];}
	val_t&	operator[](int i) {return (data[i]);}
	bool	empty() const{return front == 0;}
	int	size() const {return front;}
	void	reset() {front = 0; if (hpos) hpos->clear();}
	void	hsort();
	int	find(const val_t& x);
	bool	lt(const val_t& a, const val_t& b) const {
		     return (dec_order? b < a: a < b);
		}
};

template <typename val_t>
PrQueue_wh<val_t>::PrQueue_wh(val_t* _data, int _size, int fr, bool maxi, bool rep)
	: data(_data), capacity(_size), front(fr), dec_order(maxi), replace(rep)
{
	hpos = new Dhash<int, int>(2 * capacity, -1);
	for (int k = front / 2; --k >= 0; ) downheap(k);
}

template <typename val_t>
void PrQueue_wh<val_t>::downheap(int k, int kmax)
{
	val_t	v = data[k];
	if (kmax == 0) kmax = front;
	while (k < kmax / 2) {
	    int	l = 2 * k + 1;
	    int	r = l + 1;
	    if (r < kmax && lt(data[r], data[l])) ++l;
	    if (!lt(data[l], v)) break;
	    data[k] = data[l];
	    hpos->assign(data[k].key, k);
	    k = l;
	}
	data[k] = v;
	hpos->assign(data[k].key, k);
}

template <typename val_t>
void PrQueue_wh<val_t>::upheap(int k)
{
	val_t	v = data[k];
	int	h = (k - 1) / 2;
	while (k && lt(v, data[h])) {
	    data[k] = data[h];
	    hpos->assign(data[k].key, k);
	    k = h;
	    h = (h - 1) / 2;
	}
	data[k] = v;
	hpos->assign(data[k].key, k);
}

template <typename val_t>
val_t PrQueue_wh<val_t>::shift()
{
	val_t v = data[0];
	if (front > 0) {
	    hpos->remove(data[0].key);
	    data[0] = data[--front];
	    downheap(0);
	}
	return (v);
}

template <typename val_t>
val_t PrQueue_wh<val_t>::shift(val_t& x)
{
	val_t v = data[0];
	if (front == capacity) {
	    hpos->remove(data[0].key);
	    data[0] = data[--front];
	    downheap(0);
	}
	data[front] = x;
	upheap(front++);
	return (v);
}

template <typename val_t>
void PrQueue_wh<val_t>::put(const val_t& x, int p)
{
	if (p < 0) {			// new member
	    if (front < capacity) {
		data[front] = x;	// add
		upheap(front++);
		return;
	    } else	p = 0;
	}
	if (lt(data[p], x)) {		// exam
	    hpos->remove(data[p].key);
	    data[p] = x;		// replace
	    downheap(p);
	}

}

template <typename val_t>	// default key finder
int PrQueue_wh<val_t>::find(const val_t& x)
{
	KVpair<int, int>* kv = hpos->find(x.key);
	return kv? kv->val: -1;
}

template <typename val_t>
void PrQueue_wh<val_t>::hsort()	// heap sort data
{
	for (int k = front; --k > 0; ) {
	    swap(data[0], data[k]);
	    downheap(0, k);
	}
}

// priority queue with index

template <typename val_t>
class PrQueue_idx {
private:
	int*	aidx;
	val_t*	data;
	int	capacity;
	int	front;
	bool	dec_order;

public:
	PrQueue_idx(val_t* _data, int _size, int fr = 0, bool maxi = false);
	~PrQueue_idx() {delete[] aidx;}
	void	downheap(int k, int kk = 0);
	void	upheap(int k);
	val_t	shift();
	val_t*	shift_ptr();
	val_t	shift(val_t& x);
	int	shift_idx(int x);
	void	put(const val_t& x, int p = -1);
	val_t	gettop() {return data[aidx[0]];}
	int	gettop_idx() {return aidx[0];}
	val_t*	gettop_ptr() {return data + aidx[0];}
	val_t&	operator[](int i) {return (data[aidx[i]]);}
	bool	empty() {return front == 0;}
	int	size()	{return front;}
	bool	lt(const val_t& a, const val_t& b) {
		     if (dec_order) return (b < a);
		     else	return (a < b);
		}
};

template <typename val_t>
PrQueue_idx<val_t>::PrQueue_idx(val_t* _data, int _size, int fr, bool maxi)
	: data(_data), capacity(_size), front(fr), dec_order(maxi)
{
	aidx = new int[capacity];
	for (int i = 0; i < capacity; ++i) aidx[i] = i;
	for (int k = front / 2; --k >= 0; ) downheap(k);
}

template <typename val_t>
void PrQueue_idx<val_t>::downheap(int k, int kmax)
{
	int	v = aidx[k];
	if (kmax == 0) kmax = front;
	while (k < kmax / 2) {
	    int	l = 2 * k + 1;
	    int	r = l + 1;
	    if (r < kmax && lt(data[aidx[r]], data[aidx[l]])) ++l;
	    if (!lt(data[aidx[l]], data[v])) break;
	    aidx[k] = aidx[l];
	    k = l;
	}
	aidx[k] = v;
}

template <typename val_t>
void PrQueue_idx<val_t>::upheap(int k)
{
	int	v = aidx[k];
	int	h = (k - 1) / 2;
	while (k && lt(data[v], data[aidx[h]])) {
	    aidx[k] = aidx[h];
	    k = h;
	    h = (h - 1) / 2;
	}
	aidx[k] = v;
}

template <typename val_t>
val_t PrQueue_idx<val_t>::shift()
{
	val_t v = data[aidx[0]];
	if (front > 0) {
	    aidx[0] = aidx[--front];
	    downheap(0);
	}
	return (v);
}

template <typename val_t>
val_t* PrQueue_idx<val_t>::shift_ptr()
{
	int	v = aidx[0];
	if (front > 0) {
	    aidx[0] = aidx[--front];
	    downheap(0);
	}
	return (data + v);
}

// The following two functions return the top value AFTER shift

template <typename val_t>
val_t PrQueue_idx<val_t>::shift(val_t& x)
{
	if (front == capacity) {
	    aidx[0] = aidx[--front];
	    downheap(0);
	}
	data[aidx[front]] = x;
	upheap(front++);
	return (data[aidx[0]]);
}

template <typename val_t>
int PrQueue_idx<val_t>::shift_idx(int j)
{
	if (front == capacity) {
	    aidx[0] = aidx[--front];
	    downheap(0);
	}
	aidx[front] = j;
	upheap(front++);
	return (aidx[0]);
}

template <typename val_t>
void PrQueue_idx<val_t>::put(const val_t& x, int p)
{
	if (p < 0 && front < capacity) {
	    data[aidx[front]] = x;
	    upheap(front++);
	} else if (lt(data[aidx[p]], x)) {
	    data[aidx[p]] = x;
	    downheap(p);
	}
}

// stack

template <typename X>
class Stack {
	X*	stack;
	X*	sp;
	X*	upplim;
	bool	makehere;
public:
	Stack(int md = DefStackDepth, X* given = 0)
	    : makehere(!given)
	{
	    if (!md) md = DefStackDepth;
	    if (makehere) stack = new X[md];
	    sp = stack;
	    upplim = stack + md;
	}
	~Stack() {if (makehere) delete[] stack;}
	void	push(X val) {
	    if (sp == upplim) {
		int	sz = sp - stack;
		X*	tmp = new X[sz * 2];
		vcopy(tmp, stack, sz);
		delete[] stack;
		stack = tmp;
		sp = stack + sz;
		upplim = stack + (sz *= 2);
	    }
	    *sp++ = val;
	}
	X	pop() {
	    return (sp == stack)? 0: *--sp;
	}
};

inline	char*	prefix(char* dst, const char* src, int n) {
	memcpy(dst, src, --n);
	dst[n] = '\0';
	return dst;
}

// list of strings

class Strlist {
protected:
	char*	strbuf;
	INT*	idxlst;
	INT	sunitsize;
	INT	punitsize;
	INT	totallen;
	INT	lastlen;
	INT	maxlen;
	INT	many;
	bool	filled;
	void	format();
public:
	Strlist(int m = 1, int len = defsunit);
	Strlist(const char* str);
	Strlist(char* str, const char* delim);
	Strlist(FILE* fd, int m = 0);
#if USE_ZLIB
	Strlist(gzFile gzfd, int m = 0);
#endif
	Strlist(Strlist& src);
	~Strlist() {delete[] strbuf; delete[] idxlst;}
	char*	operator[](INT n) const {
	    return ((idxlst && n < many)? strbuf + idxlst[n]: strbuf);
	}
	INT	size() const {return (many);}
	char*	assign(const Strlist& src);
	char*	assign(const char* str);
	char*	push(const char* str);
	void	reset(int m = 0);
	bool	unfilled() const {return !filled;}
	void	setfill() {filled = true;}
	bool	empty() const {return !strbuf || !*strbuf;}
	void	undo() {--many; totallen -= lastlen;}
	char*	squeeze() {char* tmp = strbuf; strbuf = 0; return (tmp);}
	INT	longest() const {return (maxlen? maxlen - 1: 0);}
};

// hash for string key

template <class val_t>
class StrHash {
protected:
	INT	size1;
	INT	size2;
	KVpair<INT, val_t>*	hash;
	KVpair<INT, val_t>*	hz;
	Strlist*	sl;
	val_t	un_def;
public:
	StrHash(int n = 0, val_t udf = def_un_def, float hf = DefHashFact);
	StrHash(const char* fname);
	StrHash(StrHash<val_t>& src) {*this = src;}
	~StrHash() {delete sl; delete[] hash;}
	KVpair<INT, val_t>* begin() {return hash;}
	KVpair<INT, val_t>* end() {return hz;}
	char*	strkey(INT iky) {return ((*sl)[iky]);}
	KVpair<INT, val_t>* map(const char* ky, int record = 3);
	KVpair<INT, val_t>* map(INT iky, int record = 3) {
	    return map(strkey(iky));
	}
	KVpair<INT, val_t>* find(const char* ky) {
	    return map(ky, 0);
	}
	KVpair<INT, val_t>* assign(const char* ky, 
		val_t vl = 1, bool first = false) {
	    KVpair<INT, val_t>*	sh = map(ky);
	    if (!first || sh->val == un_def) sh->val = vl;
	    return (sh);
	}
	KVpair<INT, val_t>* incr(const char* ky, val_t vl = 1) {
	    KVpair<INT, val_t>* sh = map(ky);
	    if (sh->val == un_def) sh->val = vl;
	    else	sh->val += vl;
	    return (sh);
	}
	KVpair<INT, val_t>* incr(INT iky, val_t vl = 1) {
	    return incr(strkey(iky), vl);
	}
	KVpair<INT, val_t>* pile(const char* ky) {
	    int	v = sl->size();
	    KVpair<INT, val_t>* sh = map(ky);
	    if (sh->val == un_def) sh->val = v;
	    return (sh);
	}
	KVpair<INT, val_t>* pile(INT iky) {
	    return push(strkey(iky));
	}
	int	count() const;
	KVpair<INT, val_t>* press(int* n, val_t* sum = 0);
	void	clear() {
	    if (un_def) {
		KVpair<INT, val_t> zero = {0, un_def};
		vset(hash, zero, size1);
	    } else vclear(hash, size1);
	}
	int	size() {return sl->size();}
	void	resize(int s = 0);
	void	erase_sl()	{sl = 0;}
	char*	squeeze() {return (sl->squeeze());}
	val_t	undef() {return (un_def);}
};

template <class val_t>
StrHash<val_t>::StrHash(int n, val_t udf, float hf) : un_def(udf)
{
	if (n == 0) {
	    size1 = size2 = 0;
	    hash = hz = 0;
	    sl = 0;
	} else {
	    size1 = supprime(int(hf * n));
	    sl = new Strlist(size1);
	    hash = new KVpair<INT, val_t>[size1];
	    hz = hash + size1;
	    if (un_def) {
		KVpair<INT, val_t> zero = {0, un_def};
		vset(hash, zero, size1);
	    } else vclear(hash, size1);
	    for (size2 = SecondHS; size2 > size1; size2 /= 2) ;
	}
}

template <class val_t>
KVpair<INT, val_t>* StrHash<val_t>::map(const char* skey, int record)
{
	INT	key = str2uint(skey);
	INT	v = key % size1;
	INT	u = size2 - key % size2;
	INT	v0 = v;
	KVpair<INT, val_t>*	sh = hash + v;

	while (sh->key && strcmp(strkey(sh->key), skey)) {
	    v = (v + u) % size1;
	    if (v == v0) {
		resize();
		v = key % size1;
	    }
	    sh = hash + v;
	}
	if (sh->val == un_def) {
	    if (record & 2) {
		sh->key = sl->size();
		sl->push(skey);
	    } else if (!record)
		sh = 0;
	}
	return (sh);
}

template <class val_t>
int StrHash<val_t>::count() const
{
	int	c = 0;
	for (KVpair<INT, val_t>* sh = hash; sh < hz; ++sh)
	    if (sh->val != un_def) ++c;
	return (c);
}

template <class val_t>
KVpair<INT, val_t>* StrHash<val_t>::press(int* n, val_t* sum)
{
	if (sum) *sum = 0;
	KVpair<INT, val_t>*	rh = hash;
	for (KVpair<INT, val_t>* sh = hash; sh < hz; ++sh) {
	    if (sh->val != un_def) {
		*rh++ = *sh;
		if (sum) *sum += sh->val;
	    }
	}
	hz = rh;
	*n = rh - hash;
	return (hash);
}

template <class val_t>
void StrHash<val_t>::resize(int s)
{
	StrHash<val_t>	org(*this);
	size1 = supprime(s? s: 2 * size1);
	if (!size2)
	    for (size2 = SecondHS; size2 > size1; size2 /= 2) ;
	hash = new KVpair<INT, val_t>[size1];
	hz = hash + size1;
	clear();
	for (KVpair<INT, val_t>* hw = org.hash; hw < org.hz; ++hw) {
	    if (hw->val != un_def) {
		KVpair<INT, val_t>* kv = map(strkey(hw->key), 1);
		*kv = *hw;
	    }
	}
	org.erase_sl();
}

// insert sort; for any 0 <= i < num
// compf(array + i, array + num) < 0 must be satisfied

template <class val_t>
void insort(val_t* array, size_t num, 
	int (*compf)(const val_t* a, const val_t* b))
{
	for (int i = num - 2; i >= 0; --i) {
	    int	j = i;
	    val_t v = array[i];
	    for ( ; compf(array + j + 1, &v) < 0; ++j)
		array[j] = array[j + 1];
	    array[j] = v;
	}
}

class PutIntoBins
{
protected:
	int	ndiv;
	double	llmt;
	double	ulmt;
	double	width;
	double	(*nlitransform)(double x);
	double	(*invtransform)(double x);
	int	n_data;
	int	idx;
	int	sample_size;
public:
	bool	vrtl;
	int*	xval;
	double*	freq;
	void	erase() {
	    if (!vrtl) {
		if (freq) delete[] (--freq);
		if (xval) delete[] (--xval);
	    }
	    xval = 0; freq = 0; 
	}
	PutIntoBins(int n, double l, double u, double (*t)(double x) = 0, 
		double (*it)(double x) = 0, bool savedata = true);
	PutIntoBins(PutIntoBins& src) {
	    *this = src;
	    vrtl = true;
	}
	~PutIntoBins() {erase();}
	void	accumulate(int x, double y);
	void	normalize(double to = 1., double ps = 0.);
	void	reset() {idx = 0;}
	int	size();
	int	step(int n) {return (xval[n] - xval[n-1]);}
	double* squeeze() {double* dt = freq; freq = 0; return dt;}
	double*	begin() {return (freq);}
	double*	end()   {return (freq + n_data);}
	int	samples(int m = 0) {return sample_size;}
};

//	outlier detection by Dixon's method

class Dixon
{
	float	alpha;
	int	elt;	// error limit
public:
	Dixon(float a);
	int	dixon(int* res, float* data, int* odr, int num, float md = 0.);
};

//	generate integer random numbers in proportion to pcmp

class RandNumGen {
	INT*	buf;
	INT	resol;
public:
	RandNumGen(double* pcmp, INT dim, INT resol = 1024);
	~RandNumGen() {delete[] buf;}
	int	get() {return buf[int(resol * drand48())];}
};

class AddExt {
	const	char**	argv_;
	char**	arg_ext;
public:
	AddExt(int argc, const char** argv, const char* ext);
	~AddExt() {if (arg_ext) {delete[] arg_ext[0]; delete[] arg_ext;}}
	const char** add_ext() {return arg_ext? (const char**) arg_ext: argv_;}
};

extern	int	ipower(int x, int n);
extern	long	lpower(long x, int n);
extern	char*	strrealloc(char* dst, const char* src);
extern	bool	isnumber(const char* ps);
extern	bool	isBlankLine(char* str);
extern	char*	replace(char* dest, const char* sorc, const char* to, const char* from);
extern	char*	car(char*& ps);
extern	char*	car(char* token, const char* ps);
extern	char*	car(char* ps, char*& pt, char& cc);
extern	char*	carn(char* token, const char* ps, int n);
extern	char*	cdr(const char* ps);
extern	char*	cdr(const char* ps, const char* delim);
extern	int	chop(char* ps);
extern	char*	next_wd(char** sp);
extern	int	wordcmp(const char* a, const char* b, int n = INT_MAX);
extern	UPTR	fbisrch(UPTR found, const UPTR key, FILE* fd, long left, long right, int width, CMPF cmpf);
extern	int	comb2(int i, int j);
extern	int	comb3(int i, int j, int k);
extern	void	decomb(int* t, int n, int i = 2);
extern	char*	strupr(char* str);
extern	char*	strlwr(char* str);
extern	int	countbit(unsigned long x);
extern	double	rgauss();
extern	int	rpoisson(double mu);
extern	double	ktof(const char* str);
extern	const	char* getarg(int& argc, const char**& argv, bool num = false, int pv = 2);
inline	long	ktol(const char* str) {return (long(ktof(str)));}
inline	int	ktoi(const char* str) {return (int(ktof(str)));}
inline	double	inverse(double x) {return (1. / x);};

/*  End of clib.h  */

#endif
