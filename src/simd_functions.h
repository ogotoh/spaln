/*****************************************************************************
*
*	Class definition to use Intel SIMD intrinsics
*	wrapper of simd functions
*
*	Assumes capability of store/load in/from unaligned memory
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

#ifndef _SIMD_FUNCTIONS_
#define _SIMD_FUNCTIONS_

#if defined(__SSE4_1__)
#include <x86intrin.h>
#elif defined(__ARM_NEON)
#include <arm_neon.h>
#endif

static	const	CHAR	b32s2c_a[32] = 
	{ 0,  2,  4,  6,  8, 10, 12, 14,  1,  3,  5,  7,  9, 11, 13, 15,
	  0,  2,  4,  6,  8, 10, 12, 14,  1,  3,  5,  7,  9, 11, 13, 15};
static	const	CHAR	b32i2s_a[32] = 
	{ 0,  1,  4,  5,  8,  9, 12, 13,  2,  3,  6,  7, 10, 11, 14, 15,
	  0,  1,  4,  5,  8,  9, 12, 13,  2,  3,  6,  7, 10, 11, 14, 15};
static	const	CHAR	b32i2c_a[32] = 
	{ 0,  4,  8, 12,  1,  5,  9, 13,  2,  6, 10, 14,  3,  7, 11, 15,
	  0,  4,  8, 12,  1,  5,  9, 13,  2,  6, 10, 14,  3,  7, 11, 15};
static	const	CHAR	b64s2c_a[64] = 
	{ 0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
	 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
	  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31,
	 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63};
static	const	CHAR	b64i2c_a[64] = 
	{ 0,  4,  8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60,
	  1,  5,  9, 13, 17, 21, 25, 29, 33, 37, 41, 45, 49, 53, 57, 61,
	  2,  6, 10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62,
	  3,  7, 11, 15, 19, 23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63};
static	const	CHAR	b64i2s_a[64] = 
	{ 0,  1,  4,  5,  8,  9, 12, 13, 16, 17, 20, 21, 24, 25, 28, 29,
	 32, 33, 36, 37, 40, 41, 44, 45, 48, 49, 52, 53, 56, 57, 60, 61,
	  2,  3,  6,  7, 10, 11, 14, 15, 18, 19, 22, 23, 26, 27, 30, 31,
	 34, 35, 38, 39, 42, 43, 46, 47, 50, 51, 54, 55, 58, 59, 62, 63};
static	const	CHAR	b16s2c_m[3][16] =
	{{0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0},
	 {0xff, 0xff, 0, 0, 0xff, 0xff, 0, 0, 0xff, 0xff, 0, 0, 0xff, 0xff, 0, 0},
	 {0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0}
	};
static	const	CHAR	b16i2s_m[3][16] =
	{{0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0, 0xff, 0},
	 {0xff, 0xff, 0, 0, 0xff, 0xff, 0, 0, 0xff, 0xff, 0, 0, 0xff, 0xff, 0, 0},
	 {0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0}
	};
static	const	INT	i32permute[8] = {0, 4, 1, 5, 2, 6, 3, 7};

#define VecClear(d, n) \
	size_t	nn = n / Nelem * Nelem; \
	if (!nn) { \
	    vclear(d, n); \
	    return; \
	} \
	var_v	z_v = clear(); \
	for (var_t* e = d + nn; d < e; d += Nelem) store(d, z_v); \
	nn = n - nn; \
	if (nn) store(d - Nelem + nn, z_v);

#define VecCopy(d, s, n) \
	size_t	nn = n / Nelem * Nelem; \
	if (!nn) return (vcopy(d, s, n)); \
	var_t*	r = d; \
	for (const var_t* e = s + nn; s < e; d += Nelem, s += Nelem) { \
	    var_v	v_v = load(s); \
	    store(d, v_v); \
	} \
	nn = n - nn; \
	if (nn) { \
	    var_v	v_v = load(s - Nelem + nn); \
	    store(d - Nelem + nn, v_v); \
	} \
	return (r);

#define VecSet(d, c, n) \
	size_t	nn = n / Nelem * Nelem; \
	if (!nn) { \
	    vset(d, c, n); \
	    return; \
	} \
	var_v	c_v = splat(c); \
	for (var_t* e = d + nn; d < e; d += Nelem) store(d, c_v); \
	nn = n - nn; \
	if (nn) store(d - Nelem + nn, c_v);

#define VecAdd_c(d, c, n) \
	size_t	nn = n / Nelem * Nelem; \
	var_v	c_v = splat(c); \
	for (var_t* e = d + nn; d < e; d += Nelem) { \
var_v	    v_v = add(load(d), c_v); \
	    store(d, v_v); \
	} \
	nn = n - nn; \
	while (nn--) *d++ += c;

#define VecSub_c(d, c, n) \
	size_t	nn = n / Nelem * Nelem; \
	var_v	c_v = splat(c); \
	for (var_t* e = d + nn; d < e; d += Nelem) { \
var_v	    v_v = sub(load(d), c_v); \
	    store(d, v_v); \
	} \
	nn = n - nn; \
	while (nn--) *d++ -= c;

// O(log_2(Nelem)) algorithm
#define VecMax(d, n) \
	size_t	nn = n / Nelem * Nelem; \
	var_t	buf[Nelem + Nelem / 2]; \
var_v	a_v = splat(d[0]); \
	store(buf, a_v); \
	store(buf + Nelem / 2, a_v); \
	if (n > nn) vcopy(buf, d + nn, n - nn); \
	a_v = load(buf); \
	for (const var_t* e = d + nn; d < e; d += Nelem) { \
var_v	    b_v = load(d); \
	    a_v = max(a_v, b_v); \
	} \
	store(buf, a_v); \
	for (int k = Nelem; (k >>= 1); ) { \
var_v	    b_v = load(buf + k); \
	    a_v = max(a_v, b_v); \
	    store(buf, a_v); \
	} \
	return (buf[0]);

#define VecMin(d, n) \
	size_t	nn = n / Nelem * Nelem; \
	var_t	buf[Nelem + Nelem / 2]; \
var_v	a_v = splat(d[0]); \
	store(buf, a_v);  \
	store(buf + Nelem / 2, a_v);  \
	if (n > nn) vcopy(buf, d + nn, n - nn); \
	a_v = load(buf); \
	for (const var_t* e = d + nn; d < e; d += Nelem) { \
var_v	    b_v = load(d); \
	    a_v = min(a_v, b_v); \
	} \
	store(buf, a_v);  \
	for (int k = Nelem; (k >>= 1); ) { \
var_v	    b_v = load(buf + k); \
	    a_v = min(a_v, b_v); \
	    store(buf, a_v); \
	} \
	return (buf[0]);

#define VecDotP(a, b, n) \
	var_t	dp = 0.; \
	var_t	cbuf[Nelem]; \
	for ( ; n > Nelem; n -= Nelem) { \
	    var_v	a_v = load(a); \
	    var_v	b_v = load(b); \
	    var_v	c_v = mul(a_v, b_v); \
	    store(cbuf, c_v); \
	    for (int i = 0; i < Nelem; ++i) dp += cbuf[i]; \
	    a += Nelem; \
	    b += Nelem; \
	} \
	while (n-- > 0) \
	    dp += *a++ * *b++; \
	return (dp);

/****************************************************************
*	class declaration of elementary simd functions
****************************************************************/


template <typename var_t>
struct Simd_functions {};

// class specialization

#if __AVX512BW__

#define _VecRegSize_ 64

template <>
struct Simd_functions<char> {
using	var_t = char;
using	int_v = __m512i;
using	var_v = int_v;
using	var_m = __mmask64;
const	int	Nelem = _VecRegSize_ / sizeof(char);
	int_v clear() {return _mm512_setzero_si512();}
	int_v splat(const char i) {return _mm512_set1_epi8(i);}
	int_v load(const char* a) {
	    return _mm512_loadu_epi8(a);
	}
	void	store(char* a, int_v v) {_mm512_storeu_epi8(a, v);}
	int_v add(int_v u, int_v v) {
	    return _mm512_adds_epi8(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm512_subs_epi8(u, v);
	}
	int_v max(int_v u, int_v v) {
	    return _mm512_max_epi8(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm512_min_epi8(u, v);
	}
	int_v blend(int_v u, int_v v, var_m m) {
	    return _mm512_mask_blend_epi8(m, v, u);
	}
	var_m cmp_gt(int_v u, int_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 6);
	}
	var_m cmp_ge(int_v u, int_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 5);
	}
	var_m cmp_eq(int_v u, int_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 0);
	}
	int_v bit_and(int_v u, int_v v) {return _mm512_and_si512(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm512_or_si512(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm512_xor_si512(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(int_v v) {return !(_mm512_reduce_or_epi32(v));}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<CHAR> {
using	var_t = CHAR;
using	int_v = __m512i;
using	var_v = int_v;
using	var_m = __mmask64;
const	int	Nelem = _VecRegSize_ / sizeof(CHAR);
	int_v clear() {return _mm512_setzero_si512();}
	int_v splat(const CHAR i) {return _mm512_set1_epi8(i);}
	int_v load(const CHAR* a) {
	    return _mm512_loadu_epi8(a);
	}
	void	store(CHAR* a, int_v v) {_mm512_storeu_epi8(a, v);}
	int_v add(int_v u, int_v v) {
	    return _mm512_adds_epi8(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm512_subs_epi8(u, v);
	}
	int_v max(int_v u, int_v v) {
	    return _mm512_max_epi8(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm512_min_epi8(u, v);
	}
	int_v blend(int_v u, int_v v, var_m m) {
	    return _mm512_mask_blend_epi8(m, v, u);
	}
	var_m cmp_gt(int_v u, int_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 6);
	}
	var_m cmp_ge(int_v u, int_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 5);
	}
	var_m cmp_eq(int_v u, int_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 0);
	}
	int_v bit_and(int_v u, int_v v) {return _mm512_and_si512(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm512_or_si512(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm512_xor_si512(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(int_v v) {return !(_mm512_reduce_or_epi32(v));}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<short> {
using	var_t = short;
using	int_v = __m512i;
using	var_v = int_v;
using	var_m = __mmask32;
const	int	Nelem = _VecRegSize_ / sizeof(short);
	int_v clear() {return _mm512_setzero_si512();}
	int_v splat(const short i) {return _mm512_set1_epi16(i);}
	int_v load(const short* a) {
	    return _mm512_loadu_epi16(a);
	}
	void	store(short* a, int_v v) {_mm512_storeu_epi16(a, v);}
	int_v add(int_v u, int_v v) {
	    return _mm512_adds_epi16(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm512_subs_epi16(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm512_mullo_epi16(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm512_slli_epi16(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm512_srai_epi16(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm512_max_epi16(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm512_min_epi16(u, v);
	}
	int_v blend(int_v u, int_v v, var_m m) {
	    return _mm512_mask_blend_epi16(m, v, u);
	}
	var_m cmp_gt(int_v u, int_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 6);
	}
	var_m cmp_ge(int_v u, int_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 5);
	}
	var_m cmp_eq(int_v u, int_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 0);
	}
	int_v bit_and(int_v u, int_v v) {return _mm512_and_si512(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm512_or_si512(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm512_xor_si512(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(int_v v) {return !(_mm512_reduce_or_epi32(v));}
	int_v cast16to8(int_v v) {
	    int_v	b_v = _mm512_loadu_epi8(b64s2c_a);
	    return _mm512_shuffle_epi8(v, b_v);
	}
	int_v cast32to8(int_v v) {
	    int_v	b_v = _mm512_loadu_epi8(b64i2c_a);
	    return _mm512_shuffle_epi8(v, b_v);
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<SHORT> {
using	var_t = SHORT;
using	int_v = __m512i;
using	var_v = int_v;
using	var_m = __mmask32;
const	int	Nelem = _VecRegSize_ / sizeof(SHORT);
	int_v clear() {return _mm512_setzero_si512();}
	int_v splat(const SHORT i) {return _mm512_set1_epi16(i);}
	int_v load(const SHORT* a) {
	    return _mm512_loadu_epi16(a);
	}
	void	store(SHORT* a, int_v v) {_mm512_storeu_epi16(a, v);}
	int_v add(int_v u, int_v v) {
	    return _mm512_adds_epi16(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm512_subs_epi16(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm512_mullo_epi16(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm512_slli_epi16(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm512_srli_epi16(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm512_max_epi16(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm512_min_epi16(u, v);
	}
	int_v blend(int_v u, int_v v, var_m m) {
	    return _mm512_mask_blend_epi16(m, v, u);
	}
	var_m cmp_gt(int_v u, int_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 6);
	}
	var_m cmp_ge(int_v u, int_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 5);
	}
	var_m cmp_eq(int_v u, int_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 0);
	}
	int_v bit_and(int_v u, int_v v) {return _mm512_and_si512(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm512_or_si512(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm512_xor_si512(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(int_v v) {return !(_mm512_reduce_or_epi32(v));}
	int_v cast16to8(int_v v) {
	    int_v	b_v = _mm512_loadu_epi8(b64s2c_a);
	    return _mm512_shuffle_epi8(v, b_v);
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<int> {
using	var_t = int;
using	int_v = __m512i;
using	var_v = int_v;
using	var_m = __mmask16;
const	int	Nelem = _VecRegSize_ / sizeof(int);
	int_v clear() {return _mm512_setzero_si512();}
	int_v splat(const int i) {return _mm512_set1_epi32(i);}
	int_v load(const int* a) {return _mm512_loadu_epi32(a);}
	void	store(int* a, int_v v) {_mm512_storeu_epi32(a, v);}
	int_v add(int_v u, int_v v) {
	    return _mm512_add_epi32(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm512_sub_epi32(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm512_mullo_epi32(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm512_slli_epi32(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm512_srai_epi32(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm512_max_epi32(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm512_min_epi32(u, v);
	}
	int_v blend(int_v u, int_v v, var_m m) {
	    return _mm512_mask_blend_epi32(m, v, u);
	}
	var_m cmp_gt(int_v u, int_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 6);
	}
	var_m cmp_ge(int_v u, int_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 5);
	}
	var_m cmp_eq(int_v u, int_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 0);
	}
	int_v bit_and(int_v u, int_v v) {return _mm512_and_si512(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm512_or_si512(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm512_xor_si512(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(int_v v) {return !(_mm512_reduce_or_epi32(v));}
	int_v cast32to16(int_v v) {
	    int_v	b_v = _mm512_loadu_epi8(b64i2s_a);
	    return _mm512_shuffle_epi8(v, b_v);
	}
	int_v cast32to8(int_v v) {
	    int_v	b_v = _mm512_loadu_epi8(b64i2c_a);
	    return _mm512_shuffle_epi8(v, b_v);
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<INT> {
using	var_t = INT;
using	int_v = __m512i;
using	var_v = int_v;
using	var_m = __mmask16;
const	int	Nelem = _VecRegSize_ / sizeof(INT);
	int_v clear() {return _mm512_setzero_si512();}
	int_v splat(const INT i) {return _mm512_set1_epi32(i);}
	int_v load(const INT* a) {return _mm512_loadu_epi32(a);}
	void	store(INT* a, int_v v) {_mm512_storeu_epi32(a, v);}
	int_v add(int_v u, int_v v) {
	    return _mm512_add_epi32(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm512_sub_epi32(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm512_mullo_epi32(u, v);
	}
	int_v max(int_v u, int_v v) {
	    return _mm512_max_epi32(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm512_min_epi32(u, v);
	}
	int_v blend(int_v u, int_v v, var_m m) {
	    return _mm512_mask_blend_epi32(m, v, u);
	}
	var_m cmp_gt(int_v u, int_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 6);
	}
	var_m cmp_ge(int_v u, int_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 5);
	}
	var_m cmp_eq(int_v u, int_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 0);
	}
	int_v bit_and(int_v u, int_v v) {return _mm512_and_si512(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm512_or_si512(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm512_xor_si512(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(int_v v) {return !(_mm512_reduce_or_epi32(v));}
	int_v cast32to16(int_v v) {
	    int_v	b_v = _mm512_loadu_epi8(b64i2s_a);
	    return _mm512_shuffle_epi8(v, b_v);
	}
	int_v cast32to8(int_v v) {
	    int_v	b_v = _mm512_loadu_epi8(b64i2c_a);
	    return _mm512_shuffle_epi8(v, b_v);
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<float> {
using	var_t = float;
using	var_v = __m512;
using	int_v = __m512i;
using	var_m = __mmask16;
const	int	Nelem = _VecRegSize_ / sizeof(float);
	var_v clear() {return _mm512_setzero_ps();}
	var_v splat(const float f) {return _mm512_set1_ps(f);}
	var_v load(const float* a) {return _mm512_loadu_ps(a);}
	void	store(float* a, var_v v) {_mm512_storeu_ps(a, v);}
	var_v add(var_v u, var_v v) {
	    return _mm512_add_ps(u, v);
	}
	var_v sub(var_v u, var_v v) {
	    return _mm512_sub_ps(u, v);
	}
	var_v mul(var_v u, var_v v) {
	    return _mm512_mul_ps(u, v);
	}
	var_v max(var_v u, var_v v) {
	    return _mm512_max_ps(u, v);
	}
	var_v min(var_v u, var_v v) {
	    return _mm512_min_ps(u, v);
	}
	var_v blend(var_v u, var_v v, var_m m) {
	    return _mm512_mask_blend_ps(m, v, u);
	}
	var_m cmp_gt(var_v u, var_v v) {
	    return _mm512_cmp_ps_mask(u, v, 0x1e);
	}
	var_m cmp_ge(var_v u, var_v v) {
	    return _mm512_cmp_ps_mask(u, v, 0x1d);
	}
	var_m cmp_eq(var_v u, var_v v) {
	    return _mm512_cmp_ps_mask(u, v, 0x00);
	}
	int_v castf2i(var_v u) {return _mm512_castps_si512(u);}
	var_v casti2f(int_v u) {return _mm512_castsi512_ps(u);}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

#elif __AVX2__

#define _VecRegSize_ 32

template <>
struct Simd_functions<char> {
using	var_t = char;
using	int_v = __m256i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(char);
	int_v clear() {return _mm256_setzero_si256();}
	int_v splat(const char i) {return _mm256_set1_epi8(i);}
	int_v load(const char* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(char* a, int_v v) {
	    _mm256_storeu_si256((__m256i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm256_adds_epi8(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm256_subs_epi8(u, v);
	}
	int_v max(int_v u, int_v v) {
	    return _mm256_max_epi8(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm256_min_epi8(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm256_cmpgt_epi8(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi8(u, v), _mm256_cmpeq_epi8(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm256_cmpeq_epi8(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm256_and_si256(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm256_or_si256(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm256_xor_si256(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(int_v v) {return _mm256_testnzc_si256(v, v);}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<CHAR> {
using	var_t = CHAR;
using	int_v = __m256i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(CHAR);
	int_v clear() {return _mm256_setzero_si256();}
	int_v splat(const CHAR i) {return _mm256_set1_epi8(i);}
	int_v load(const CHAR* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(CHAR* a, int_v v) {
	    _mm256_storeu_si256((__m256i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm256_adds_epi8(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm256_subs_epi8(u, v);
	}
	int_v max(int_v u, int_v v) {
	    return _mm256_max_epi8(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm256_min_epi8(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm256_cmpgt_epi8(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi8(u, v), _mm256_cmpeq_epi8(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm256_cmpeq_epi8(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm256_and_si256(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm256_or_si256(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm256_xor_si256(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(int_v v) {return _mm256_testnzc_si256(v, v);}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<short> {
using	var_t = short;
using	int_v = __m256i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(short);
	int_v clear() {return _mm256_setzero_si256();}
	int_v splat(const short i) {return _mm256_set1_epi16(i);}
	int_v load(const short* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(short* a, int_v v) {
	    _mm256_storeu_si256((__m256i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm256_adds_epi16(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm256_subs_epi16(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm256_mullo_epi16(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm256_slli_epi16(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm256_srai_epi16(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm256_max_epi16(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm256_min_epi16(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm256_cmpgt_epi16(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi16(u, v), _mm256_cmpeq_epi16(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm256_cmpeq_epi16(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm256_and_si256(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm256_or_si256(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm256_xor_si256(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(int_v v) {return _mm256_testnzc_si256(v, v);}
	int_v cast16to8(int_v v) {
	    int_v	b_v = _mm256_loadu_si256((__m256i const*) b32s2c_a);
	    b_v = _mm256_shuffle_epi8(v, b_v);
	    return (_mm256_permute4x64_epi64(b_v, 216));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<SHORT> {
using	var_t = SHORT;
using	int_v = __m256i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(SHORT);
	int_v clear() {return _mm256_setzero_si256();}
	int_v splat(const SHORT i) {return _mm256_set1_epi16(i);}
	int_v load(const SHORT* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(SHORT* a, int_v v) {
	    _mm256_storeu_si256((__m256i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm256_adds_epi16(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm256_subs_epi16(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm256_mullo_epi16(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm256_slli_epi16(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm256_srli_epi16(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm256_max_epi16(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm256_min_epi16(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm256_cmpgt_epi16(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi16(u, v), _mm256_cmpeq_epi16(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm256_cmpeq_epi16(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm256_and_si256(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm256_or_si256(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm256_xor_si256(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(int_v v) {return _mm256_testnzc_si256(v, v);}
	int_v cast16to8(int_v v) {
	    int_v	b_v = _mm256_loadu_si256((__m256i const*) b32s2c_a);
	    b_v = _mm256_shuffle_epi8(v, b_v);
	    return (_mm256_permute4x64_epi64(b_v, 216));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<int> {
using	var_t = int;
using	int_v = __m256i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(int);
	int_v clear() {return _mm256_setzero_si256();}
	int_v splat(const int i) {return _mm256_set1_epi32(i);}
	int_v load(const int* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(int* a, int_v v) {
	    _mm256_storeu_si256((__m256i *) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm256_add_epi32(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm256_sub_epi32(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm256_mullo_epi32(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm256_slli_epi32(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm256_srai_epi32(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm256_max_epi32(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm256_min_epi32(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm256_cmpgt_epi32(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi32(u, v), _mm256_cmpeq_epi32(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm256_cmpeq_epi32(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm256_and_si256(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm256_or_si256(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm256_xor_si256(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(int_v v) {return _mm256_testnzc_si256(v, v);}
	int_v cast32to16(int_v v) {
	    int_v	b_v = _mm256_loadu_si256((__m256i const*) b32i2s_a);
	    b_v = _mm256_shuffle_epi8(v, b_v);
	    return (_mm256_permute4x64_epi64(b_v, 216));	// 0xb8
	}
	int_v cast32to8(int_v v) {
	    int_v	b_v = _mm256_loadu_si256((__m256i const*) b32i2c_a);
	    b_v = _mm256_shuffle_epi8(v, b_v);
	    int_v	i_v = _mm256_loadu_si256((__m256i const*) i32permute);
	    return (_mm256_permutevar8x32_epi32(b_v, i_v));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<INT> {
using	var_t = INT;
using	int_v = __m256i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(INT);
	int_v clear() {return _mm256_setzero_si256();}
	int_v splat(const INT i) {return _mm256_set1_epi32(i);}
	int_v load(const INT* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(INT* a, int_v v) {
	    _mm256_storeu_si256((__m256i *) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm256_add_epi32(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm256_sub_epi32(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm256_mullo_epi32(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm256_slli_epi32(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm256_srli_epi32(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm256_max_epi32(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm256_min_epi32(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm256_cmpgt_epi32(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi32(u, v), _mm256_cmpeq_epi32(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm256_cmpeq_epi32(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm256_and_si256(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm256_or_si256(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm256_xor_si256(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(int_v v) {return _mm256_testnzc_si256(v, v);}
	int_v cast32to16(int_v v) {
	    int_v	b_v = _mm256_loadu_si256((__m256i const*) b32i2s_a);
	    b_v = _mm256_shuffle_epi8(v, b_v);
	    return (_mm256_permute4x64_epi64(b_v, 216));	// 0xb8
	}
	int_v cast32to8(int_v v) {
	    int_v	b_v = _mm256_loadu_si256((__m256i const*) b32i2c_a);
	    b_v = _mm256_shuffle_epi8(v, b_v);
	    int_v	i_v = _mm256_loadu_si256((__m256i const*) i32permute);
	    return (_mm256_permutevar8x32_epi32(b_v, i_v));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<float> {
using	var_t = float;
using	var_v = __m256;
using	var_m = __m256;
using	int_v = __m256i;
const	int	Nelem = _VecRegSize_ / sizeof(float);
	var_v clear() {return _mm256_setzero_ps();}
	var_v splat(const float f) {return _mm256_set1_ps(f);}
	var_v load(const float* a) {
	    return _mm256_loadu_ps(a);
	}
	void	store(float* a, var_v v) {_mm256_storeu_ps(a, v);}
	var_v add(var_v u, var_v v) {
	    return _mm256_add_ps(u, v);
	}
	var_v sub(var_v u, var_v v) {
	    return _mm256_sub_ps(u, v);
	}
	var_v mul(var_v u, var_v v) {
	    return _mm256_mul_ps(u, v);
	}
	var_v max(var_v u, var_v v) {
	    return _mm256_max_ps(u, v);
	}
	var_v min(var_v u, var_v v) {
	    return _mm256_min_ps(u, v);
	}
	var_v blend(var_v u, var_v v, var_m m) {
	    return _mm256_blendv_ps(v, u, m);
	}
	var_m cmp_gt(var_v u, var_v v) {
	    return _mm256_cmp_ps(u, v, 0x1e);
	}
	var_m cmp_ge(var_v u, var_v v) {
	    return _mm256_cmp_ps(u, v, 0x1d);
	}
	var_m cmp_eq(var_v u, var_v v) {
	    return _mm256_cmp_ps(u, v, 0x00);
	}

	int_v castf2i(var_v u) {return _mm256_castps_si256(u);}
	var_v casti2f(int_v u) {return _mm256_castsi256_ps(u);}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
	var_t	vecdotp(const var_t* a, const var_t* b, int n) {
	    var_t	dp = 0.;
	    var_t	cbuf[Nelem];
	    for ( ; n >= Nelem; n -= Nelem) {
		var_v	a_v = _mm256_loadu_ps(a);
		var_v	b_v = _mm256_loadu_ps(b);
		var_v	c_v = _mm256_dp_ps(a_v, b_v, 0xff);
		_mm256_storeu_ps(cbuf, c_v);
		dp += cbuf[0] + cbuf[4];
		a += Nelem;
		b += Nelem;
	    }
	    while (n-- > 0)
		dp += *a++ * *b++;
	    return (dp);
	}
};

#elif __SSE4_1__

#define _VecRegSize_ 16

template <>
struct Simd_functions<char> {
using	var_t = char;
using	int_v = __m128i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(char);
	int_v clear() {return _mm_setzero_si128();}
	int_v splat(const char i) {return _mm_set1_epi8(i);}
	int_v load(const char* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(char* a, int_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm_adds_epi8(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm_subs_epi8(u, v);
	}
	int_v max(int_v u, int_v v) {
	    return _mm_max_epi8(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm_min_epi8(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm_cmpgt_epi8(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi8(u, v), _mm_cmpeq_epi8(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm_cmpeq_epi8(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm_and_si128(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm_or_si128(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm_xor_si128(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(int_v v) {return _mm_test_all_zeros(v, v);}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<CHAR> {
using	var_t = CHAR;
using	int_v = __m128i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(CHAR);
	int_v clear() {return _mm_setzero_si128();}
	int_v splat(const CHAR i) {return _mm_set1_epi8(i);}
	int_v load(const CHAR* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(CHAR* a, int_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm_adds_epi8(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm_subs_epi8(u, v);
	}
	int_v max(int_v u, int_v v) {
	    return _mm_max_epi8(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm_min_epi8(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm_cmpgt_epi8(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi8(u, v), _mm_cmpeq_epi8(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm_cmpeq_epi8(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm_and_si128(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm_or_si128(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm_xor_si128(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(int_v v) {return _mm_test_all_zeros(v, v);}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<short> {
using	var_t = short;
using	int_v = __m128i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(short);
	int_v clear() {return _mm_setzero_si128();}
	int_v splat(const short i) {return _mm_set1_epi16(i);}
	int_v load(const short* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(short* a, int_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm_adds_epi16(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm_subs_epi16(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm_mullo_epi16(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm_slli_epi16(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return  _mm_srai_epi16(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm_max_epi16(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm_min_epi16(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm_cmpgt_epi16(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi16(u, v), _mm_cmpeq_epi16(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm_cmpeq_epi16(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm_and_si128(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm_or_si128(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm_xor_si128(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(int_v v) {return _mm_test_all_zeros(v, v);}
	int_v cast16to8(int_v v) {
	    int_v	b_v = _mm_loadu_si128((__m128i const*) b32s2c_a);
	    return _mm_shuffle_epi8(v, b_v);
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<SHORT> {
using	var_t = SHORT;
using	int_v = __m128i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(SHORT);
	int_v clear() {return _mm_setzero_si128();}
	int_v splat(const SHORT i) {return _mm_set1_epi16(i);}
	int_v load(const SHORT* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(SHORT* a, int_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm_adds_epi16(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm_subs_epi16(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm_mullo_epi16(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm_slli_epi16(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm_srli_epi16(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm_max_epi16(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm_min_epi16(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm_cmpgt_epi16(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi16(u, v), _mm_cmpeq_epi16(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm_cmpeq_epi16(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm_and_si128(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm_or_si128(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm_xor_si128(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(int_v v) {return _mm_test_all_zeros(v, v);}
	int_v cast16to8(int_v v) {
	    int_v	b_v = _mm_loadu_si128((__m128i const*) b32s2c_a);
	    return _mm_shuffle_epi8(v, b_v);
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<int> {
using	var_t = int;
using	int_v = __m128i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(int);
	int_v clear() {return _mm_setzero_si128();}
	int_v splat(const int i) {return _mm_set1_epi32(i);}
	int_v load(const int* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(int* a, int_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm_add_epi32(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm_sub_epi32(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm_mullo_epi32(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm_slli_epi32(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return  _mm_srai_epi32(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm_max_epi32(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm_min_epi32(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm_cmpgt_epi32(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi32(u, v), _mm_cmpeq_epi32(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm_cmpeq_epi32(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm_and_si128(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm_or_si128(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm_xor_si128(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(int_v v) {return _mm_test_all_zeros(v, v);}
	int_v cast32to16(int_v v) {
	    int_v	b_v = _mm_loadu_si128((__m128i const*) b32i2s_a);
	    return _mm_shuffle_epi8(v, b_v);
	}
	int_v cast32to8(int_v v) {
	    int_v	b_v = _mm_loadu_si128((__m128i const*) b32i2c_a);
	    return _mm_shuffle_epi8(v, b_v);
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<INT> {
using	var_t = INT;
using	int_v = __m128i;
using	var_v = int_v;
const	int	Nelem = _VecRegSize_ / sizeof(INT);
	int_v clear() {return _mm_setzero_si128();}
	int_v splat(const INT i) {return _mm_set1_epi32(i);}
	int_v load(const INT* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(INT* a, int_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	int_v add(int_v u, int_v v) {
	    return _mm_add_epi32(u, v);
	}
	int_v sub(int_v u, int_v v) {
	    return _mm_sub_epi32(u, v);
	}
	int_v mul(int_v u, int_v v) {
	    return _mm_mullo_epi32(u, v);
	}
	int_v shiftl(int_v u, const int i) {
	    return _mm_slli_epi32(u, i);
	}
	int_v shiftr(int_v u, const int i) {
	    return _mm_srli_epi32(u, i);
	}
	int_v max(int_v u, int_v v) {
	    return _mm_max_epi32(u, v);
	}
	int_v min(int_v u, int_v v) {
	    return _mm_min_epi32(u, v);
	}
	int_v blend(int_v u, int_v v, int_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	int_v cmp_gt(int_v u, int_v v) {return _mm_cmpgt_epi32(u, v);}
	int_v cmp_ge(int_v u, int_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi32(u, v), _mm_cmpeq_epi32(u, v));
	}
	int_v cmp_eq(int_v u, int_v v) {return _mm_cmpeq_epi32(u, v);}
	int_v bit_and(int_v u, int_v v) {return _mm_and_si128(u, v);}
	int_v bit_or(int_v u, int_v v) {return _mm_or_si128(u, v);}
	int_v bit_xor(int_v u, int_v v) {return _mm_xor_si128(u, v);}
	int_v bit_andnot(int_v u, int_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(int_v v) {return _mm_test_all_zeros(v, v);}
	int_v cast32to16(int_v v) {
	    int_v	b_v = _mm_loadu_si128((__m128i const*) b32i2s_a);
	    return _mm_shuffle_epi8(v, b_v);
	}
	int_v cast32to8(int_v v) {
	    int_v	b_v = _mm_loadu_si128((__m128i const*) b32i2c_a);
	    return _mm_shuffle_epi8(v, b_v);
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<float> {
using	var_t = float;
using	var_v = __m128;
using	var_m = __m128;
using	int_v = __m128i;
const	int	Nelem = _VecRegSize_ / sizeof(float);
	var_v clear() {return _mm_setzero_ps();}
	var_v splat(const float f) {return _mm_set_ps1(f);}
	var_v load(const float* a) {return _mm_loadu_ps(a);}
	void	store(float* a, var_v v) {_mm_storeu_ps(a, v);}
	var_v add(var_v u, var_v v) {
	    return _mm_add_ps(u, v);
	}
	var_v sub(var_v u, var_v v) {
	    return _mm_sub_ps(u, v);
	}
	var_v mul(var_v u, var_v v) {
	    return _mm_mul_ps(u, v);
	}
	var_v max(var_v u, var_v v) {
	    return _mm_max_ps(u, v);
	}
	var_v min(var_v u, var_v v) {
	    return _mm_min_ps(u, v);
	}
	var_v blend(var_v u, var_v v, var_m m) {
	    return _mm_blendv_ps(v, u, m);
	}
	var_m cmp_gt(var_v u, var_v v) {
	    return _mm_cmpgt_ps(u, v);
	}
	var_m cmp_ge(var_v u, var_v v) {
	    return _mm_cmpge_ps(u, v);
	}
	var_m cmp_eq(var_v u, var_v v) {
	    return _mm_cmpeq_ps(u, v);
	}

	int_v castf2i(var_v u) {return _mm_castps_si128(u);}
	var_v casti2f(int_v u) {return _mm_castsi128_ps(u);}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
	var_t	vecdotp(const var_t* a, const var_t* b, int n){
	    var_t	dp = 0.;			// dot product
	    var_t	cbuf[Nelem];
	    for ( ; n > Nelem; n -= Nelem) {
		var_v	a_v = _mm_loadu_ps(a);
		var_v	b_v = _mm_loadu_ps(b);
		var_v	c_v = _mm_dp_ps(a_v, b_v, 0xff);
		_mm_storeu_ps(cbuf, c_v);
		dp += cbuf[0];
		a += Nelem;
		b += Nelem;
	    }
	    while (n-- > 0)
		dp += *a++ * *b++;
	    return (dp);
	}
};

#elif __ARM_NEON

#define _VecRegSize_ 16

template <>
struct Simd_functions<char> {
using	var_t = signed char;
using	int_v = int8x16_t;
using	var_v = int_v;
using	int_m = uint8x16_t;
const	int	Nelem = _VecRegSize_ / sizeof(char);
	int_v clear() {return vdupq_n_s8(0);}
	int_v splat(const int8_t i) {return (vdupq_n_s8(i));}
	int_v load(int8_t const* a) {return vld1q_s8(a);}
	void	store(int8_t* a, int_v v) {vst1q_s8(a, v);}
	int_v add(int_v u, int_v v) {return vqaddq_s8(u, v);}
	int_v sub(int_v u, int_v v) {return vqsubq_s8(u, v);}
	int_v mul(int_v u, int_v v) {return vmulq_s8(u, v);}
	int_v max(int_v u, int_v v) {return vmaxq_s8(u, v);}
	int_v min(int_v u, int_v v) {return vminq_s8(u, v);}
	int_v blend(int_v u, int_v v, int_m m) {return vbslq_s8(m, u, v);}
	int_m cmp_gt(int_v u, int_v v) {return vcgtq_s8(u, v);}
	int_m cmp_ge(int_v u, int_v v) {return vcgeq_s8(u, v);}
	int_m cmp_eq(int_v u, int_v v) {return vceqq_s8(u, v);}
	int_v bit_and(int_v u, int_v v) {return vandq_s8(u, v);}
	int_v bit_or(int_v u, int_v v) {return vorrq_s8(u, v);}
	int_v bit_xor(int_v u, int_v v) {return veorq_s8(u, v);}
	int_v bit_andnot(int_v u, int_v v) {
	    return vandq_s8(u, vmvnq_s8(v));
	}
	int	all_zero(int_v v) {
	    return !(vgetq_lane_s64(vreinterpretq_s64_s8(v), 0) |
		vgetq_lane_s64(vreinterpretq_s64_s8(v), 1));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<CHAR> {
using	var_t = CHAR;
using	int_v = uint8x16_t;
using	var_v = int_v;
using	int_m = uint8x16_t;
const	int	Nelem = _VecRegSize_ / sizeof(CHAR);
	int_v clear() {return vdupq_n_u8(0);}
	int_v splat(const uint8_t i) {return vdupq_n_u8(i);}
	int_v load(uint8_t const* a) {return vld1q_u8(a);}
	void	store(uint8_t* a, int_v v) {vst1q_u8(a, v);}
	int_v add(int_v u, int_v v) {return vqaddq_u8(u, v);}
	int_v sub(int_v u, int_v v) {return vqsubq_u8(u, v);}
	int_v mul(int_v u, int_v v) {return vmulq_u8(u, v);}
	int_v max(int_v u, int_v v) {return vmaxq_u8(u, v);}
	int_v min(int_v u, int_v v) {return vminq_u8(u, v);}
	int_v blend(int_v u, int_v v, int_m m){return vbslq_u8(m, u, v);}
	int	all_zero(int_v v) {
	    return !(vgetq_lane_s64(vreinterpretq_s64_u8(v), 0) |
		vgetq_lane_s64(vreinterpretq_s64_u8(v), 1));
	}
	int_m cmp_gt(int_v u, int_v v) {return vcgtq_u8(u, v);}
	int_m cmp_ge(int_v u, int_v v) {return vcgeq_u8(u, v);}
	int_m cmp_eq(int_v u, int_v v) {return vceqq_u8(u, v);}
	int_v bit_and(int_v u, int_v v) {return vandq_u8(u, v);}
	int_v bit_or(int_v u, int_v v) {return vorrq_u8(u, v);}
	int_v bit_xor(int_v u, int_v v) {return veorq_u8(u, v);}
	int_v bit_andnot(int_v u, int_v v) {
	    return vandq_u8(u, vmvnq_u8(v));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<short> {
using	var_t = short;
using	int_v = int16x8_t;
using	var_v = int_v;
using	int_m = uint16x8_t;
const	int	Nelem = _VecRegSize_ / sizeof(short);
	int_v clear() {return vdupq_n_s16(0);}
	int_v splat(int16_t i) {return vdupq_n_s16(i);}
	int_v load(int16_t const* a) {return vld1q_s16(a);}
	void	store(int16_t* a, int_v v) {vst1q_s16(a, v);}
	int_v add(int_v u, int_v v) {return vqaddq_s16(u, v);}
	int_v sub(int_v u, int_v v) {return vqsubq_s16(u, v);}
	int_v mul(int_v u, int_v v) {return vmulq_s16(u, v);}
	int_v max(int_v u, int_v v) {return vmaxq_s16(u, v);}
	int_v min(int_v u, int_v v) {return vminq_s16(u, v);}
	int_v blend(int_v u, int_v v, int_m m) {
	    return vbslq_s16(m, u, v);
	}
	int_m cmp_gt(int_v u, int_v v) {return vcgtq_s16(u, v);}
	int_m cmp_ge(int_v u, int_v v) {return vcgeq_s16(u, v);}
	int_m cmp_eq(int_v u, int_v v) {return vceqq_s16(u, v);}
	int_v bit_and(int_v u, int_v v) {return vandq_s16(u, v);}
	int_v bit_or(int_v u, int_v v) {return vorrq_s16(u, v);}
	int_v bit_xor(int_v u, int_v v) {return veorq_s16(u, v);}
	int_v bit_andnot(int_v u, int_v v) {
	    return vandq_s16(u, vmvnq_s16(v));
	}
	int	all_zero(int_v v) {
	    return !(vgetq_lane_s64(vreinterpretq_s64_s16(v), 0) |
		vgetq_lane_s64(vreinterpretq_s64_s16(v), 1));
	}
	int_v cast16to8(int_v v) {
	    int8x16_t	w = vreinterpretq_s8_s16(v);
	    uint8x16_t	b_v = vld1q_u8((uint8_t const*) b32s2c_a);
	    return vreinterpretq_s16_s8(vqtbl1q_s8(w, b_v));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<SHORT> {
using	var_t = SHORT;
using	int_v = uint16x8_t;
using	var_v = int_v;
using	int_m = uint16x8_t;
const	int	Nelem = _VecRegSize_ / sizeof(SHORT);
	int_v clear() {return vdupq_n_u16(0);}
	int_v splat(uint16_t i) {return vdupq_n_u16(i);}
	int_v load(uint16_t const* a) {return vld1q_u16(a);}
	void	store(uint16_t* a, int_v v) {vst1q_u16(a, v);}
	int_v add(int_v u, int_v v) {return vqaddq_u16(u, v);}
	int_v sub(int_v u, int_v v) {return vqsubq_u16(u, v);}
	int_v mul(int_v u, int_v v) {return vmulq_u16(u, v);}
	int_v max(int_v u, int_v v) {return vmaxq_u16(u, v);}
	int_v min(int_v u, int_v v) {return vminq_u16(u, v);}
	int_v blend(int_v u, int_v v, int_m m) {
	    return vbslq_u16(m, u, v);
	}
	int_m cmp_gt(int_v u, int_v v) {return vcgtq_u16(u, v);}
	int_m cmp_ge(int_v u, int_v v) {return vcgeq_u16(u, v);}
	int_m cmp_eq(int_v u, int_v v) {return vceqq_u16(u, v);}
	int_v bit_and(int_v u, int_v v) {return vandq_u16(u, v);}
	int_v bit_or(int_v u, int_v v) {return vorrq_u16(u, v);}
	int_v bit_xor(int_v u, int_v v) {return veorq_u16(u, v);}
	int_v bit_andnot(int_v u, int_v v) {
	    return vandq_u16(u, vmvnq_u16(v));
	}
	int	all_zero(int_v v) {
	    return !(vgetq_lane_s64(vreinterpretq_s64_u16(v), 0) |
		vgetq_lane_s64(vreinterpretq_s64_u16(v), 1));
	}
	int_v cast16to8(int_v v) {
	    uint8x16_t	w = vreinterpretq_u8_u16(v);
	    uint8x16_t	b_v = vld1q_u8((uint8_t const*) b32s2c_a);
	    return vreinterpretq_u16_u8(vqtbl1q_u8(w, b_v));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<int> {
using	var_t = int;
using	int_v = int32x4_t;
using	var_v = int_v;
using	int_m = uint32x4_t;
const	int	Nelem = _VecRegSize_ / sizeof(int);
	int_v clear() {return vdupq_n_s32(0);}
	int_v splat(int32_t i) {return vdupq_n_s32(i);}
	int_v load(int32_t const* a) {return vld1q_s32(a);}
	void	store(int32_t* a, int_v v) {vst1q_s32(a, v);}
	int_v add(int_v u, int_v v) {return vqaddq_s32(u, v);}
	int_v sub(int_v u, int_v v) {return vqsubq_s32(u, v);}
	int_v mul(int_v u, int_v v) {return vmulq_s32(u, v);}
	int_v max(int_v u, int_v v) {return vmaxq_s32(u, v);}
	int_v min(int_v u, int_v v) {return vminq_s32(u, v);}
	int_v blend(int_v u, int_v v, int_m m) {
	    return vbslq_s32(m, u, v);
	}
	int_m cmp_gt(int_v u, int_v v) {return vcgtq_s32(u, v);}
	int_m cmp_ge(int_v u, int_v v) {return vcgeq_s32(u, v);}
	int_m cmp_eq(int_v u, int_v v) {return vceqq_s32(u, v);}
	int_v bit_and(int_v u, int_v v) {return vandq_s32(u, v);}
	int_v bit_or(int_v u, int_v v) {return vorrq_s32(u, v);}
	int_v bit_xor(int_v u, int_v v) {return veorq_s32(u, v);}
	int_v bit_andnot(int_v u, int_v v) {
	    return vandq_s32(u, vmvnq_s32(v));
	}
	int	all_zero(int_v v) {
	    return !(vgetq_lane_s64(vreinterpretq_s64_s32(v), 0) |
		vgetq_lane_s64(vreinterpretq_s64_s32(v), 1));
	}
	int_v cast32to16(int_v v) {
	    int8x16_t	w = vreinterpretq_s8_s32(v);
	    uint8x16_t	b_v = vld1q_u8((uint8_t const*) b32i2s_a);
	    return vreinterpretq_s32_s8(vqtbl1q_s8(w, b_v));
	}
	int_v cast32to8(int_v v) {
	    int8x16_t	w = vreinterpretq_s8_s32(v);
	    uint8x16_t	b_v = vld1q_u8((uint8_t const*) b32i2c_a);
	    return vreinterpretq_s32_s8(vqtbl1q_s8(w, b_v));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<INT> {
using	var_t = INT;
using	int_v = uint32x4_t;
using	var_v = int_v;
using	int_m = uint32x4_t;
const	int	Nelem = _VecRegSize_ / sizeof(INT);
	int_v clear() {return vdupq_n_u32(0);}
	int_v splat(uint32_t i) {return vdupq_n_u32(i);}
	int_v load(uint32_t const* a) {return vld1q_u32(a);}
	void	store(uint32_t* a, int_v v) {vst1q_u32(a, v);}
	int_v add(int_v u, int_v v) {return vqaddq_u32(u, v);}
	int_v sub(int_v u, int_v v) {return vqsubq_u32(u, v);}
	int_v mul(int_v u, int_v v) {return vmulq_u32(u, v);}
	int_v max(int_v u, int_v v) {return vmaxq_u32(u, v);}
	int_v min(int_v u, int_v v) {return vminq_u32(u, v);}
	int_v blend(int_v u, int_v v, int_m m) {
	    return vbslq_u32(m, u, v);
	}
	int_m cmp_gt(int_v u, int_v v) {return vcgtq_u32(u, v);}
	int_m cmp_ge(int_v u, int_v v) {return vcgeq_u32(u, v);}
	int_m cmp_eq(int_v u, int_v v) {return vceqq_u32(u, v);}
	int_v bit_and(int_v u, int_v v) {return vandq_u32(u, v);}
	int_v bit_or(int_v u, int_v v) {return vorrq_u32(u, v);}
	int_v bit_xor(int_v u, int_v v) {return veorq_u32(u, v);}
	int_v bit_andnot(int_v u, int_v v) {
	    return vandq_u32(u, vmvnq_u32(v));
	}
	int	all_zero(int_v v) {
	    return !(vgetq_lane_s64(vreinterpretq_s64_u32(v), 0) |
		vgetq_lane_s64(vreinterpretq_s64_u32(v), 1));
	}
	int_v cast32to16(int_v v) {
	    uint8x16_t	w = vreinterpretq_u8_u32(v);
	    uint8x16_t	b_v = vld1q_u8((uint8_t const*) b32i2s_a);
	    return vreinterpretq_u32_u8(vqtbl1q_u8(w, b_v));
	}
	int_v cast32to18(int_v v) {
	    uint8x16_t	w = vreinterpretq_u8_u32(v);
	    uint8x16_t	b_v = vld1q_u8((uint8_t const*) b32i2c_a);
	    return vreinterpretq_u32_u8(vqtbl1q_u8(w, b_v));
	}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n){VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n){VecMin(d, n)}
};

template <>
struct Simd_functions<float> {
using	var_t = float;
using	var_v = float32x4_t;
using	int_v = int32x4_t;
using	var_m = uint32x4_t;
const	int	Nelem = _VecRegSize_ / sizeof(float);
	var_v clear() {return vdupq_n_f32(0.f);}
	var_v splat(float32_t f) {return vdupq_n_f32(f);}
	var_v load(float32_t const* a) {return vld1q_f32(a);}
	void	store(float32_t* a, var_v v) {vst1q_f32(a, v);}
	var_v add(var_v u, var_v v) {
	    return vaddq_f32(u, v);
	}
	var_v sub(var_v u, var_v v) {
	    return vsubq_f32(u, v);
	}
	var_v mul(var_v u, var_v v) {
	    return vmulq_f32(u, v);
	}
	var_v max(var_v u, var_v v) {
	    return vmaxq_f32(u, v);
	}
	var_v min(var_v u, var_v v) {
	    return vminq_f32(u, v);
	}
	var_v blend(var_v u, var_v v, var_m m) {
	    return vbslq_f32(m, u, v);
	}
	var_m cmp_gt(var_v u, var_v v) {
	    return vcgtq_f32(u, v);
	}
	var_m cmp_ge(var_v u, var_v v) {
	    return vcltq_f32(v, u);
	}
	var_m cmp_eq(var_v u, var_v v) {
	    return vceqq_f32(u, v);
	}
	int_v castf2i(var_v u) {return vreinterpretq_s32_f32(u);}
	var_v casti2f(int_v u) {return vreinterpretq_f32_s32(u);}
	void	vecclear(var_t* dst, const size_t n) {VecClear(dst, n)}
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n)
	    {VecCopy(dst, src, n)}
	void	vecset(var_t* dst, const var_t& c, const size_t n)
	    {VecSet(dst, c, n)}
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n)
	    {VecAdd_c(dst, c, n)}
	void	vecsub_c(var_t* dst, const var_t& c, const size_t n)
	    {VecSub_c(dst, c, n)}
	var_t	vecmax(const var_t* d, const size_t n) {VecMax(d, n)}
	var_t	vecmin(const var_t* d, const size_t n) {VecMin(d, n)}
	var_t	vecdotp(const var_t* a, const var_t* b, int n)
	    {VecDotP(a, b, n)}
};

#else	// __ARM_NEON
#define _VecRegSize_ 0;
#endif	// __AVX512BW__

/*************************************************************************
	common public functions
*************************************************************************/

template <typename X>
void vec_clear(X* ary, const size_t n)
{
	Simd_functions<X> sf;
	sf.vecclear(ary, n);
}

template <typename X>
X* vec_copy(X* dst, const X* src, const size_t n)
{
	Simd_functions<X> sf;
	return (sf.veccopy(dst, src, n));
}

template <typename X>
void vec_set(X* ary, const X& c, const size_t n)
{
	Simd_functions<X> sf;
	sf.vecset(ary, c, n);
}

template <typename X>
void vec_add_c(X* ary, const X& c, const size_t n)
{
	Simd_functions<X> sf;
	sf.vecadd_c(ary, c, n);
}

template <typename X>
void vec_sub_c(X* ary, const X& c, const size_t n)
{
	Simd_functions<X> sf;
	sf.vecsub_c(ary, c, n);
}

template <typename X>
X vec_max(const X* ary, const size_t n)
{
	Simd_functions<X> sf;
	return sf.vecmax(ary, n);
}

template <typename X>
X vec_min(const X* ary, const size_t n)
{
	Simd_functions<X> sf;
	return sf.vecmin(ary, n);
}

template <typename X>
X vec_dotp(const X* a_ary, const X* b_ary, const int n)
{
	Simd_functions<X> sf;
	return sf.vecdotp(a_ary, b_ary, n);
}

#endif	// _SIMD_FUNCTIONS_
