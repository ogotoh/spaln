/*****************************************************************************
*
*	Collection of functions for general use
*
*	vclear	vcopy	vset	vreverse	vmax	vmin
*	max	min	max3	min3	max4	min4
*	ipower	gswap
*	next_wd	wordcmp	isBlankLine  
*	car	cdr	prefix	chop	replace
*	queue	stack	fbisrch	insort
*	comp2	comp3	decomp
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
*	Copyright(c) Osamu Gotoh <<o.gotoh@i.kyoto-u.ac.jp>>
*****************************************************************************/
#ifndef  _CLIB_
#define  _CLIB_

#include <math.h>

extern	INT	supprime(INT n);

static	const	int	KILO = 1024;
static	const	int	MEGA = 1048576;
static	const	long	GIGA = 1073741824;
static	const	int	DefStackDepth = 128;
static	const	int	defsunit = 128;
static	const	int	defpunit = 16;
static	const	INT	SecondHS = 8;
static	const	float	DefHashFact = 1.2;
static	const	double	epsilon = 1e-6;
static	const	int	HashovLS = 12;

extern	INT	str2uint(const char* str);

template <typename X> X max(X x, X y)
{
	return (x > y)? x: y;
}

template <typename X> X min(X x, X y)
{
	return (x < y)? x: y;
}

template <typename X> X* vcopy(X* dst, const X* src, int n)
{
	if (!dst) dst = new X[n];
	X*	head = dst;
	while (n-- > 0) *dst++ = *src++;
	return (head);
}

template <typename X> X* vset(X* dst, const X& val, int n)
{
	if (!dst) dst = new X[n];
	X*	head = dst;
	while (n-- > 0) *dst++ = val;
	return (head);
}

template  <typename X> void vclear(X* ary, const int n = 1)
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
	for (int i = 0, j = n - 1; i < j; ++i, --j) {
	    X	tmp = array[i];
	    array[i] = array[j];
	    array[j] = tmp;
	}
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

template <typename X> X max3(X x, X y, X z)
{
	if (y > x) x = y;
	return (max(x, z));
}

template <typename X> X min3(X x, X y, X z)
{
	if (y < x) x = y;
	return (min(x, z));
}

template <typename X> X max4(X x, X y, X z, X w)
{
	if (y > x) x = y;
	if (w > z) z = w;
	return (max(x, z));
}

template <typename X> X min4(X x, X y, X z, X w)
{
	if (y < x) x = y;
	if (w < z) z = w;
	return (min(x, z));
}

template <typename X> void gswap(X& x, X& y)
{
	X temp = x;
	x = y;
	y = temp;
}

// key-value pair used by Dhash

template <class key_t, class val_t>
struct KVpair {key_t key; val_t val;};

// double hash

template <class key_t, class val_t>
class Dhash {
protected:
	INT	size;
	INT	size2;
	KVpair<key_t, val_t>*	hash;
	KVpair<key_t, val_t>*	hz;
	val_t	undef;
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
	    if (!first || sh->val == undef) sh->val = vl;
	    return (sh);
	}
	KVpair<key_t, val_t>*	incr(key_t ky, val_t vl = 1){
	    KVpair<key_t, val_t>*	sh = map(ky);
	    if (sh->val == undef) sh->val = 0;
	    sh->val += vl;
	    return (sh);
	}
	int	count();
	KVpair<key_t, val_t>*	press(int* n, val_t* sum = 0);
	void	clear() {
	    KVpair<key_t, val_t> zero = {0, undef};
	    if (undef)	vset(hash, zero, size);
	    else	vclear(hash, size);
	}
	void	resize(int s = 0);
	void	remove(key_t ky) {
	    KVpair<key_t, val_t>*	sh = map(ky);
	    sh->val = undef;
	}
	KVpair<key_t, val_t>*	squeeze(int* n) {
	    KVpair<key_t, val_t>* tmp = press(n);
	    hash = 0;
	    return (tmp);
	}
};

template <class key_t, class val_t>
Dhash<key_t, val_t>::Dhash(int n, val_t ud, float hf) : undef(ud)
{
	if (n == 0) {
	    size = size2 = 0;
	    hash = hz = 0;
	} else {
	    size = supprime(int(hf * n));
	    hash = new KVpair<key_t, val_t>[size];
	    hz = hash + size;
	    clear();
	    for (size2 = SecondHS; size2 > size; size2 /= 2) ;
	}
}

template <class key_t, class val_t>
KVpair<key_t, val_t>* Dhash<key_t, val_t>::map(key_t key, bool record)
{
	INT	v = key % size;
	INT	u = size2 - key % size2;
	INT	v0 = v;
	KVpair<key_t, val_t>*	sh = hash + v;

	while (sh->val != undef && sh->key != key) {
	    v = (v + u) % size;
	    if (v == v0) resize();
	    sh = hash + v;
	}
	if (sh->val == undef) {
	    if (record) sh->key = key;
	    else	sh = 0;
	}
	return (sh);
}

template <class key_t, class val_t>
int	Dhash<key_t, val_t>::count()
{
	int	c = 0;
	for (KVpair<key_t, val_t>* sh = hash; sh < hz; ++sh)
	    if (sh->val != undef) ++c;
	return (c);
}

template <class key_t, class val_t>
KVpair<key_t, val_t>* Dhash<key_t, val_t>::press(int* n, val_t* sum)
{
	if (sum) *sum = 0;
	KVpair<key_t, val_t>*	rh = hash;
	for (KVpair<key_t, val_t>* sh = hash; sh < hz; ++sh) {
	    if (sh->val != undef) {
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
	size = supprime(s? s: 2 * size);
	if (!size2)
	    for (size2 = SecondHS; size2 > size; size2 /= 2) ;
	hash = new KVpair<key_t, val_t>[size];
	hz = hash + size;
	clear();
	for (KVpair<key_t, val_t>* hw = hold; hw < zold; ++hw)
	    if (hw->val != undef) assign(hw->key, hw->val);
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
	    gswap(val, queue[qp]);
	    if (++qp == qsize) qp = 0;
	    return (val);
	}
	X	unshift(X val) {
	    gswap(val, queue[qp]);
	    if (qp) --qp;
	    else qp = qsize - 1;
	    return (val);
	}
	X	oldest() {return queue[qp];}
};

// priority queue

template <typename val_t>
class PrQueue {
private:
	val_t*	data;
	int	capacity;
	int	front;
	bool	dec_order;	// default is ascending order
	bool	replace;	// may replace older value
	Dhash<int, int>*	hpos;	// heap position of key
public:				// initial data size, ascending order, never replace data
	PrQueue(val_t* _data, int _size, int fr = 0, bool maxi = false, bool rep = false);
	~PrQueue() {delete hpos;}
	void	downheap(int k, int kk = 0);
	void	upheap(int k);
	val_t	shift();
	val_t	shift(val_t& x);
	void	put(const val_t& x, int p = -1);
	val_t	gettop() {return data[0];}
	val_t&	operator[](int i) {return (data[i]);}
	bool	empty() {return front == 0;}
	int	size()	{return front;}
	void	reset() {front = 0;}
	void	hsort();
	int	find(const val_t& x);
	bool	lt(const val_t& a, const val_t& b) {
		     return (dec_order? b < a: a < b);
		}
};

template <typename val_t>
PrQueue<val_t>::PrQueue(val_t* _data, int _size, int fr, bool maxi, bool rep)
	: data(_data), capacity(_size), front(fr), dec_order(maxi), replace(rep)
{
	hpos = (rep && capacity > HashovLS)? new Dhash<int, int>(2 * capacity, -1): 0;
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
	    if (hpos) hpos->assign(data[k].key, k);
	    k = l;
	}
	data[k] = v;
	if (hpos) hpos->assign(data[k].key, k);
}

template <typename val_t>
void PrQueue<val_t>::upheap(int k)
{
	val_t	v = data[k];
	int	h = (k - 1) / 2;
	while (k && lt(v, data[h])) {
	    data[k] = data[h];
	    if (hpos) hpos->assign(data[k].key, k);
	    k = h;
	    h = (h - 1) / 2;
	}
	data[k] = v;
	if (hpos) hpos->assign(data[k].key, k);
}

template <typename val_t>
val_t PrQueue<val_t>::shift()
{
	val_t v = data[0];
	if (front > 0) {
	    if (hpos) hpos->remove(data[0].key);
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
	    if (hpos) hpos->remove(data[0].key);
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
	    if (hpos) hpos->remove(data[p].key);
	    data[p] = x;		// replace
	    downheap(p);
	}

}

template <typename val_t>	// default key finder
int PrQueue<val_t>::find(const val_t& x)
{
	if (hpos) {
	    KVpair<int, int>* kv = hpos->find(x.key);
	    return kv? kv->val: -1;
	}
        for (int i = 0; i < front; ++i)
            if (x == data[i]) return (i);
        return (-1);
}

template <typename val_t>
void PrQueue<val_t>::hsort()	// heap sort data
{
	for (int k = front; --k > 0; ) {
	    gswap(data[0], data[k]);
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
	int*	idxlst;
	int	sunitsize;
	int	punitsize;
	int	totallen;
	int	lastlen;
	int	maxlen;
	int	many;
	bool	filled;
public:
	Strlist(int m = 1, int len = defsunit);
	Strlist(const char* str);
	Strlist(char* str, const char* delim);
	Strlist(FILE* fd, int m = 0);
	Strlist(Strlist& src);
	~Strlist() {delete[] strbuf; delete[] idxlst;}
	char*	operator[](int n) {
	    return ((idxlst && n < many)? strbuf + idxlst[n]: strbuf);
	}
	int	terms() {return (many);}
	char*	assign(const Strlist& src);
	char*	assign(const char* str);
	char*	push(const char* str);
	void	reset(int m = 0);
	bool	unfilled() {return !filled;}
	void	setfill() {filled = true;}
	bool	empty() {return !strbuf || !*strbuf;}
	void	undo() {--many; totallen -= lastlen;}
	char*	squeeze() {char* tmp = strbuf; strbuf = 0; return (tmp);}
	int	longest() {return (maxlen? maxlen - 1: 0);}
};

// hash for string key

template <class val_t>
class StrHash : public Strlist {
protected:
	INT	size;
	INT	size2;
	KVpair<INT, val_t>*	hash;
	KVpair<INT, val_t>*	hz;
public:
	StrHash(int n = 0, float hf = DefHashFact);
	StrHash(const char* fname);
	~StrHash() {delete[] hash;}
	KVpair<INT, val_t>* begin() {return hash;}
	KVpair<INT, val_t>* end() {return hz;}
	const char*	memname(INT iky) {return (strbuf + iky);}
	KVpair<INT, val_t>* map(const char* ky, bool record = true);
	KVpair<INT, val_t>* map(INT iky, bool record = true) {
	    return map(strbuf + iky);
	}
	KVpair<INT, val_t>* find(const char* ky) {
	    return map(ky, false);
	}
	KVpair<INT, val_t>* assign(const char* ky, 
		val_t vl = 1, bool first = false) {
	    KVpair<INT, val_t>*	sh = map(ky);
	    if (!first || !sh->val) sh->val = vl;
	    return (sh);
	}
	KVpair<INT, val_t>* incr(const char* ky, val_t vl = 1) {
	    KVpair<INT, val_t>* sh = map(ky);
	    sh->val += vl;
	    return (sh);
	}
	KVpair<INT, val_t>* incr(INT iky, val_t vl = 1) {
	    return incr(strbuf + iky, vl);
	}
	KVpair<INT, val_t>* pile(const char* ky) {
	    int	v = many;
	    KVpair<INT, val_t>* sh = map(ky);
	    if (!sh->val) sh->val = v;
	    return (sh);
	}
	KVpair<INT, val_t>* pile(INT iky) {
	    return push(strbuf + iky);
	}
	int	count();
	KVpair<INT, val_t>* press(int* n, val_t* sum = 0);
	void	clear() {vclear(hash, size);}
	void	resize(int s = 0);
};

template <class val_t>
StrHash<val_t>::StrHash(int n, float hf) : Strlist(n? n: 1)
{
	if (n == 0) {
	    size = size2 = 0;
	    hash = hz = 0;
	} else {
	    size = supprime(int(hf * n));
	    hash = new KVpair<INT, val_t>[size];
	    hz = hash + size;
	    vclear(hash, size);
	    for (size2 = SecondHS; size2 > size; size2 /= 2) ;
	}
	push("");	// 0-th recode
}

template <class val_t>
KVpair<INT, val_t>* StrHash<val_t>::map(const char* strkey, bool record)
{
	INT	key = str2uint(strkey);
	INT	v = key % size;
	INT	u = size2 - key % size2;
	INT	v0 = v;
	KVpair<INT, val_t>*	sh = hash + v;

	while (sh->key && strcmp(strbuf + sh->key, strkey)) {
	    v = (v + u) % size;
	    if (v == v0) resize();
	    sh = hash + v;
	}
	if (!sh->val) {
	    if (record) sh->key = push(strkey) - strbuf;
	    else	sh = 0;
	}
	return (sh);
}

template <class val_t>
int StrHash<val_t>::count()
{
	int	c = 0;
	for (KVpair<INT, val_t>* sh = hash; sh < hz; ++sh)
	    if (sh->val) ++c;
	return (c);
}

template <class val_t>
KVpair<INT, val_t>* StrHash<val_t>::press(int* n, val_t* sum)
{
	if (sum) *sum = 0;
	KVpair<INT, val_t>*	rh = hash;
	for (KVpair<INT, val_t>* sh = hash; sh < hz; ++sh) {
	    if (sh->val) {
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
	KVpair<INT, val_t>*	hold = hash;
	KVpair<INT, val_t>*	zold = hz;
	size = supprime(s? s: 2 * size);
	if (!size2)
	    for (size2 = SecondHS; size2 > size; size2 /= 2) ;
	hash = new KVpair<INT, val_t>[size];
	hz = hash + size;
	clear();
	for (KVpair<INT, val_t>* hw = hold; hw < zold; ++hw)
	    if (hw->val) incr(hw->key, hw->val);
	delete[] hold;
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
	double	(*transform)(double x);
	double	(*invtransform)(double x);
	double*	data;
	int	sample_size;
public:
	bool	vrtl;
	void	erase() {
	    if (!vrtl && data) delete[] (--data); 
	    data = 0; 
	}
	PutIntoBins(int n, double l, double u, double (*t)(double x) = 0, 
		double (*it)(double x) = 0, int m = 1);
	PutIntoBins(PutIntoBins& src) {
	    *this = src;
	    vrtl = true;
	}
	~PutIntoBins() {erase();}
	void	accumulate(int x, double y, int m = 0);
	void	normalize(double to = 1., double ps = 0.);
	int	x2n(double& x) {
	    if (transform) x = transform(x);
	    x = (x - llmt) / width;
	    return int(x + epsilon);
	}
	int	n_in_slice(int n) {
	    if (!invtransform) return (1);
	    double      base = invtransform(llmt + width * n);
	    double      ceil = invtransform(llmt + width * (++n));
	    return (int(ceil) - int(base));
	}
	double	normcoef(int n) {
	    if (!invtransform) return (1.);
	    double	base = invtransform(llmt + width * n);
	    double	ceil = invtransform(llmt + width * (++n));
	    int		nl = int(ceil) - int(base);
	    if (!nl) return (0.);
	    return ((ceil - base) / nl);
	}
	double* squeeze() {double* dt = data; data = 0; return dt;}
	double*	begin() {return (data);}
	double*	end()   {return (data + ndiv);}
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
extern	int	comp2(int i, int j);
extern	int	comp3(int i, int j, int k);
extern	void	decomp(int i, int* t, int n);
extern	char*	strupr(char* str);
extern	char*	strlwr(char* str);
extern	int	countbit(unsigned long x);
extern	double	rgauss();
extern	int	rpoisson(double mu);
extern	double	ktof(const char* str);
extern	const	char* getarg(int& argc, const char**& argv, bool num = false);
inline	long	ktol(const char* str) {return (long(ktof(str)));}
inline	int	ktoi(const char* str) {return (int(ktof(str)));}
inline	double	inverse(double x) {return (1. / x);};

/*  End of clib.h  */

#endif
