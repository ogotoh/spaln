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

#define SIMD_INLINE inline __attribute__((always_inline))

static	const	CHAR	b32_a[32] = 
	{ 0,  2,  4,  6,  8, 10, 12, 14,  1, 3, 5, 7, 9, 11, 13, 15,
	  0,  2,  4,  6,  8, 10, 12, 14,  1, 3, 5, 7, 9, 11, 13, 15};
static	const	CHAR	b64_a[64] = 
	{ 0,  2,  4,  6,  8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30,
	 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62,
	  1,  3,  5,  7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31,
	 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63};
static	const	CHAR	b16_m[3][16] =
	{{0xff, 0, 0, 0xff, 0xff, 0, 0, 0xff, 0xff, 0, 0, 0xff, 0xff, 0, 0, 0xff},
	 {0xff, 0xff, 0, 0, 0, 0, 0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0xff, 0xff},
	 {0xff, 0xff, 0xff, 0xff, 0, 0, 0, 0, 0, 0, 0, 0, 0xff, 0xff, 0xff, 0xff}
	};

// class declaration of elementary simd functions

template <typename var_t, int Nelem, typename regist_v>
struct Simd_functions {
	regist_v clear();
	regist_v splat(const var_t i);
	regist_v load(const var_t* a);
	void	store(var_t* a, regist_v);
	regist_v add(regist_v u, regist_v v);
	regist_v sub(regist_v u, regist_v v);
	regist_v mul(regist_v u, regist_v v);
	regist_v shiftl(regist_v u, const int i);
	regist_v shiftr(regist_v u, const int i);
	regist_v max(regist_v u, regist_v v);
	regist_v min(regist_v u, regist_v v);
	regist_v blend(regist_v u, regist_v v, regist_v m);	// m? u: v
// Warning! Not m? v: u as _mm_blendv_epi8 etc
	regist_v cmp_eq(regist_v u, regist_v v);
	regist_v cmp_ge(regist_v u, regist_v v);
	regist_v cmp_gt(regist_v u, regist_v v);
	regist_v bit_and(regist_v u, regist_v v);
	regist_v bit_or(regist_v u, regist_v v);
	regist_v bit_xor(regist_v u, regist_v v);
	regist_v bit_andnot(regist_v u, regist_v v);	// u & !v
// Warning! Not !u & v as _mm_andnot_si128 etc
	regist_v cast16to8(regist_v v);
	int	ph_v(regist_v v);
};

// class specialization

#if __SSE4_1__

template <typename regist_v>
struct Simd_functions<char, 16, regist_v> {
	regist_v clear() {return _mm_setzero_si128();}
	regist_v splat(const char i) {return _mm_set1_epi8(i);}
	regist_v load(const char* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(char* a, regist_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm_adds_epi8(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm_subs_epi8(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm_max_epi8(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm_min_epi8(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm_cmpgt_epi8(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi8(u, v), _mm_cmpeq_epi8(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm_cmpeq_epi8(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm_and_si128(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm_or_si128(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm_xor_si128(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(regist_v v) {return _mm_test_all_zeros(v, v);}
};

template <typename regist_v>
struct Simd_functions<CHAR, 16, regist_v> {
	regist_v clear() {return _mm_setzero_si128();}
	regist_v splat(const CHAR i) {return _mm_set1_epi8(i);}
	regist_v load(const CHAR* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(CHAR* a, regist_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm_adds_epi8(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm_subs_epi8(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm_max_epi8(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm_min_epi8(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm_cmpgt_epi8(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi8(u, v), _mm_cmpeq_epi8(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm_cmpeq_epi8(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm_and_si128(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm_or_si128(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm_xor_si128(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(regist_v v) {return _mm_test_all_zeros(v, v);}
};

template <typename regist_v>
struct Simd_functions<short, 8, regist_v> {
	regist_v clear() {return _mm_setzero_si128();}
	regist_v splat(const short i) {return _mm_set1_epi16(i);}
	regist_v load(const short* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(short* a, regist_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm_adds_epi16(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm_subs_epi16(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm_mullo_epi16(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm_slli_epi16(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return  _mm_srai_epi16(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm_max_epi16(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm_min_epi16(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm_cmpgt_epi16(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi16(u, v), _mm_cmpeq_epi16(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm_cmpeq_epi16(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm_and_si128(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm_or_si128(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm_xor_si128(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(regist_v v) {return _mm_test_all_zeros(v, v);}
	regist_v cast16to8(regist_v v) {
	    regist_v	b_v = _mm_loadu_si128((__m128i const*) b32_a);
	    return _mm_shuffle_epi8(v, b_v);
	}
};

template <typename regist_v>
struct Simd_functions<SHORT, 8, regist_v> {
	regist_v clear() {return _mm_setzero_si128();}
	regist_v splat(const SHORT i) {return _mm_set1_epi16(i);}
	regist_v load(const SHORT* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(SHORT* a, regist_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm_adds_epi16(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm_subs_epi16(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm_mullo_epi16(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm_slli_epi16(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm_srli_epi16(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm_max_epi16(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm_min_epi16(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm_cmpgt_epi16(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi16(u, v), _mm_cmpeq_epi16(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm_cmpeq_epi16(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm_and_si128(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm_or_si128(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm_xor_si128(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(regist_v v) {return _mm_test_all_zeros(v, v);}
	regist_v cast16to8(regist_v v) {
	    regist_v	b_v = _mm_loadu_si128((__m128i const*) b32_a);
	    return _mm_shuffle_epi8(v, b_v);
	}
};

template <typename regist_v>
struct Simd_functions<int, 4, regist_v> {
	regist_v clear() {return _mm_setzero_si128();}
	regist_v splat(const int i) {return _mm_set1_epi32(i);}
	regist_v load(const int* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(int* a, regist_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm_add_epi32(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm_sub_epi32(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm_mullo_epi32(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm_slli_epi32(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return  _mm_srai_epi32(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm_max_epi32(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm_min_epi32(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm_cmpgt_epi32(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi32(u, v), _mm_cmpeq_epi32(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm_cmpeq_epi32(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm_and_si128(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm_or_si128(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm_xor_si128(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(regist_v v) {return _mm_test_all_zeros(v, v);}
};

template <typename regist_v>
struct Simd_functions<INT, 4, regist_v> {
	regist_v clear() {return _mm_setzero_si128();}
	regist_v splat(const INT i) {return _mm_set1_epi32(i);}
	regist_v load(const INT* a) {
	    return _mm_loadu_si128((__m128i const*) a);
	}
	void	store(INT* a, regist_v v) {
	    _mm_storeu_si128((__m128i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm_add_epi32(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm_sub_epi32(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm_mullo_epi32(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm_slli_epi32(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm_srli_epi32(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm_max_epi32(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm_min_epi32(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm_cmpgt_epi32(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm_or_si128
	    (_mm_cmpgt_epi32(u, v), _mm_cmpeq_epi32(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm_cmpeq_epi32(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm_and_si128(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm_or_si128(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm_xor_si128(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm_andnot_si128(v, u);}
	int	all_zero(regist_v v) {return _mm_test_all_zeros(v, v);}
};

template <typename regist_v>
struct Simd_functions<float, 4, regist_v> {
	regist_v clear() {return _mm_setzero_ps();}
	regist_v splat(const float f) {return _mm_set_ps1(f);}
	regist_v load(const float* a) {return _mm_loadu_ps(a);}
	void	store(float* a, regist_v v) {_mm_storeu_ps(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm_add_ps(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm_sub_ps(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm_mul_ps(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm_max_ps(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm_min_ps(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm_blendv_ps(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm_cmpgt_ps(u, v);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm_cmpge_ps(u, v);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm_cmpeq_ps(u, v);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm_and_ps(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm_or_ps(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm_xor_ps(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm_andnot_ps(v, u);}
};
#endif	// __SSE4_1__

#if __ARM_NEON

template <typename regist_v>
struct Simd_functions<char, 16, regist_v> {
	regist_v clear() {return vdupq_n_s8(0);}
	regist_v splat(const char i) {return vdupq_n_s8(i);}
	regist_v load(const char* a) {return vld1q_s8(a);}
	void	store(char* a, regist_v v) {vst1q_s8(a, v);}
	regist_v add(regist_v u, regist_v v) {return vqaddq_s8(u, v);}
	regist_v sub(regist_v u, regist_v v) {return vqsubq_s8(u, v);}
	regist_v mul(regist_v u, regist_v v) {return vmulq_s8(u, v);}
	regist_v max(regist_v u, regist_v v) {return vmaxq_s8(u, v);}
	regist_v min(regist_v u, regist_v v) {return vminq_s8(u, v);}
	regist_v shiftl(regist_v u, const int i) {
	    return vshlq_n_s8(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return vshrq_n_s8(u, i);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return vbslq_s8(m, v, u);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return vcgtq_s8(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return vcgeq_s8(u, v);}
	regist_v cmp_eq(regist_v u, regist_v v) {return vceqq_s8(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return vandq_s8(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return vorrq_s8(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return veorq_s8(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {
	    return vandq_s8(u, vmvnq_s8(v));
	}
};

template <typename regist_v>
struct Simd_functions<CHAR, 16, regist_v> {
	regist_v clear() {return vdupq_n_u8(0);}
	regist_v splat(const CHAR i) {return vdupq_n_u8(i);}
	regist_v load(const CHAR* a) {return vld1q_u8(a);}
	void	store(CHAR* a, regist_v v) {vst1q_u8(a, v);}
	regist_v add(regist_v u, regist_v v) {return vqaddq_s8(u, v);}
	regist_v sub(regist_v u, regist_v v) {return vqsubq_u8(u, v);}
	regist_v mul(regist_v u, regist_v v) {return vmulq_u8(u, v);}
	regist_v max(regist_v u, regist_v v) {return vmaxq_u8(u, v);}
	regist_v min(regist_v u, regist_v v) {return vminq_u8(u, v);}
	regist_v shiftl(regist_v u, const int i) {
	    return vshlq_n_u8(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return vshrq_n_u8(u, i);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return vbslq_u8(m, v, u);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return vcgtq_u8(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return vcgeq_u8(u, v);}
	regist_v cmp_eq(regist_v u, regist_v v) {return vceqq_u8(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return vandq_u8(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return vorrq_u8(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return veorq_u8(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {
	    return vandq_u8(u, vmvnq_u8(v));
	}
};

template <typename regist_v>
struct Simd_functions<short, 8, regist_v> {
	regist_v clear() {return vdupq_n_s16(0);}
	regist_v splat(const short i) {return vdupq_n_s16(i);}
	regist_v load(const short* a) {return vld1q_s16(a);}
	void	store(short* a, regist_v v) {vst1_s16(a, v);}
	regist_v add(regist_v u, regist_v v) {return vqaddq_s16(u, v);}
	regist_v sub(regist_v u, regist_v v) {return vqsubq_s16(u, v);}
	regist_v mul(regist_v u, regist_v v) {return vmulq_s16(u, v);}
	regist_v max(regist_v u, regist_v v) {return vmaxq_s16(u, v);}
	regist_v min(regist_v u, regist_v v) {return vminq_s16(u, v);}
	regist_v shiftl(regist_v u, const int i) {
	    return vshlq_n_s16(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return vshrq_n_s16(u, i);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return vbslq_s16(m, v, u);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return vcgtq_s16(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return vcgeq_s16(u, v);}
	regist_v cmp_eq(regist_v u, regist_v v) {return vceqq_s16(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return vandq_s16(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return vorrq_s16(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return veorq_s16(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {
	    return vandq_s16(u, vmvnq_s16(v));
	}
	regist_v cast16to8(regist_v v) {
	    uint8x16_t	w = vreinterpretq_u8_s16(v);
	    for (int i = 0, s = 1; i < 3; ++i, s *= 2) {
		uint8x16_t	u = vshlq_n_u8(w, s);
		uint8x16_t	m = vld1q_u8(b16_m[i]);
		w = vbslq_u8(m, w, u);
	    }
	    return vreinterpretq_s16_u8(w);
	}
};

template <typename regist_v>
struct Simd_functions<SHORT, 8, regist_v> {
	regist_v clear() {return vdupq_n_u16(0);}
	regist_v splat(const SHORT i) {return vdupq_n_u16(i);}
	regist_v load(const SHORT* a) {return vld1q_u16(a);}
	void	store(SHORT* a, regist_v v) {vst1(a, v);}
	regist_v add(regist_v u, regist_v v) {return vqaddq_u16(u, v);}
	regist_v sub(regist_v u, regist_v v) {return vqsubq_u16(u, v);}
	regist_v mul(regist_v u, regist_v v) {return vmulq_u16(u, v);}
	regist_v max(regist_v u, regist_v v) {return vmaxq_u16(u, v);}
	regist_v min(regist_v u, regist_v v) {return vminq_u16(u, v);}
	regist_v shiftl(regist_v u, const int i) {
	    return vshlq_n_u16(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return vshrq_n_s16(u, i);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return vbslq_u16(m, v, u);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return vcgtq_u16(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return vcgeq_u16(u, v);}
	regist_v cmp_eq(regist_v u, regist_v v) {return vceqq_u16(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return vandq_u16(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return vorrq_u16(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return veorq_u16(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {
	    return vandq_u16(u, vmvnq_u16(v));
	}
	regist_v cast16to8(regist_v v) {
	    uint8x16_t	w = vreinterpretq_u8_u16(v);
	    for (int i = 0, s = 1; i < 3; ++i, s *= 2) {
		uint8x16_t	u = vshlq_n_u8(w, s);
		uint8x16_t	m = vld1q_u8(b16_m[i]);
		w = vbslq_u8(m, w, u);
	    }
	    return vreinterpretq_u16_u8(w);
	}
};

template <typename regist_v>
struct Simd_functions<int, 4, regist_v> {
	regist_v clear() {return vdupq_n_s32(0);}
	regist_v splat(const int i) {return vdupq_n_s32(i);}
	regist_v load(const int* a) {return vld1q_u16(a);}
	void	store(int* a, regist_v v) {vst1q_s32(a, v);}
	regist_v add(regist_v u, regist_v v) {return vqaddq_s32(u, v);}
	regist_v sub(regist_v u, regist_v v) {return vqsubq_s32(u, v);}
	regist_v mul(regist_v u, regist_v v) {return vmulq_s32(u, v);}
	regist_v max(regist_v u, regist_v v) {return vmaxq_s32(u, v);}
	regist_v min(regist_v u, regist_v v) {return vminq_s32(u, v);}
	regist_v shiftl(regist_v u, const int i) {
	    return vshlq_n_s32(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return vshrq_n_s32(u, i);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return vbslq_s32(m, v, u);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return vcgtq_s32(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return vcgeq_s32(u, v);}
	regist_v cmp_eq(regist_v u, regist_v v) {return vceqq_s32(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return vandq_s32(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return vorrq_s32(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return veorq_s32(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {
	    return vandq_s32(u, vmvnq_s32(v));
	}
};

template <typename regist_v>
struct Simd_functions<INT, 4, regist_v> {
	regist_v clear() {return vdupq_n_u32(0);}
	regist_v splat(const INT i) {return vdupq_n_u32(i);}
	regist_v load(const INT* a) {return vld1q_u32(a);}
	void	store(INT* a, regist_v v) {vst1q_u32(a, v);}
	regist_v add(regist_v u, regist_v v) {return vqaddq_u32(u, v);}
	regist_v sub(regist_v u, regist_v v) {return vqsubq_u32(u, v);}
	regist_v mul(regist_v u, regist_v v) {return vmulq_u32(u, v);}
	regist_v max(regist_v u, regist_v v) {return vmaxq_u32(u, v);}
	regist_v min(regist_v u, regist_v v) {return vminq_u32(u, v);}
	regist_v shiftl(regist_v u, const int i) {
	    return vshlq_n_u32(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return vshrq_n_u32(u, i);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return vbslq_u32(m, v, u);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return vcgtq_u32(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return vcgeq_u32(u, v);}
	regist_v cmp_eq(regist_v u, regist_v v) {return vceqq_u32(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return vandq_u32(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return vorrq_u32(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return veorq_u32(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {
	    return vandq_u32(u, vmvnq_u32(v));
	}
};
#endif	// __ARM_NEON

#if __AVX2__

template <typename regist_v>
struct Simd_functions<char, 32, regist_v> {
	regist_v clear() {return _mm256_setzero_si256();}
	regist_v splat(const char i) {return _mm256_set1_epi8(i);}
	regist_v load(const char* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(char* a, regist_v v) {
	    _mm256_storeu_si256((__m256i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm256_adds_epi8(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm256_subs_epi8(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm256_mullo_epi8(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm256_max_epi8(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm256_min_epi8(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm256_cmpgt_epi8(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi8(u, v), _mm256_cmpeq_epi8(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm256_cmpeq_epi8(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm256_and_si256(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm256_or_si256(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm256_xor_si256(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(regist_v v) {return _mm256_testnzc_si256(v, v);}
};

template <typename regist_v>
struct Simd_functions<CHAR, 32, regist_v> {
	regist_v clear() {return _mm256_setzero_si256();}
	regist_v splat(const CHAR i) {return _mm256_set1_epi8(i);}
	regist_v load(const CHAR* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(CHAR* a, regist_v v) {
	    _mm256_storeu_si256((__m256i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm256_adds_epi8(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm256_subs_epi8(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm256_mullo_epi8(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm256_max_epi8(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm256_min_epi8(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm256_cmpgt_epi8(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi8(u, v), _mm256_cmpeq_epi8(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm256_cmpeq_epi8(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm256_and_si256(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm256_or_si256(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm256_xor_si256(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(regist_v v) {return _mm256_testnzc_si256(v, v);}
};

template <typename regist_v>
struct Simd_functions<short, 16, regist_v> {
	regist_v clear() {return _mm256_setzero_si256();}
	regist_v splat(const short i) {return _mm256_set1_epi16(i);}
	regist_v load(const short* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(short* a, regist_v v) {
	    _mm256_storeu_si256((__m256i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm256_adds_epi16(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm256_subs_epi16(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm256_mullo_epi16(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm256_slli_epi16(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm256_srai_epi16(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm256_max_epi16(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm256_min_epi16(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm256_cmpgt_epi16(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi16(u, v), _mm256_cmpeq_epi16(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm256_cmpeq_epi16(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm256_and_si256(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm256_or_si256(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm256_xor_si256(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(regist_v v) {return _mm256_testnzc_si256(v, v);}
	regist_v cast16to8(regist_v v) {
	    regist_v	b_v = _mm256_loadu_si256((__m256i*) b32_a);
	    b_v = _mm256_shuffle_epi8(v, b_v);
	    return (_mm256_permute4x64_epi64(b_v, 216));
	}
};

template <typename regist_v>
struct Simd_functions<SHORT, 16, regist_v> {
	regist_v clear() {return _mm256_setzero_si256();}
	regist_v splat(const SHORT i) {return _mm256_set1_epi16(i);}
	regist_v load(const SHORT* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(SHORT* a, regist_v v) {
	    _mm256_storeu_si256((__m256i*) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm256_adds_epi16(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm256_subs_epi16(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm256_mullo_epi16(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm256_slli_epi16(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm256_srli_epi16(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm256_max_epi16(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm256_min_epi16(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm256_cmpgt_epi16(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi16(u, v), _mm256_cmpeq_epi16(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm256_cmpeq_epi16(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm256_and_si256(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm256_or_si256(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm256_xor_si256(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(regist_v v) {return _mm256_testnzc_si256(v, v);}
	regist_v cast16to8(regist_v v) {
	    regist_v	b_v = _mm256_loadu_si256((__m256i*) b32_a);
	    b_v = _mm256_shuffle_epi8(v, b_v);
	    return (_mm256_permute4x64_epi64(b_v, 216));
	}
};

template <typename regist_v>
struct Simd_functions<int, 8, regist_v> {
	regist_v clear() {return _mm256_setzero_si256();}
	regist_v splat(const int i) {return _mm256_set1_epi32(i);}
	regist_v load(const int* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(int* a, regist_v v) {
	    _mm256_storeu_si256((__m256i *) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm256_add_epi32(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm256_sub_epi32(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm256_mullo_epi32(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm256_slli_epi32(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm256_srai_epi32(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm256_max_epi32(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm256_min_epi32(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm256_cmpgt_epi32(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi32(u, v), _mm256_cmpeq_epi32(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm256_cmpeq_epi32(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm256_and_si256(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm256_or_si256(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm256_xor_si256(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(regist_v v) {return _mm256_testnzc_si256(v, v);}
};

template <typename regist_v>
struct Simd_functions<INT, 8, regist_v> {
	regist_v clear() {return _mm256_setzero_si256();}
	regist_v splat(const INT i) {return _mm256_set1_epi32(i);}
	regist_v load(const INT* a) {
	    return _mm256_loadu_si256((__m256i const*) a);
	}
	void	store(INT* a, regist_v v) {
	    _mm256_storeu_si256((__m256i *) a, v);
	}
	regist_v add(regist_v u, regist_v v) {
	    return _mm256_add_epi32(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm256_sub_epi32(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm256_mullo_epi32(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm256_slli_epi32(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm256_srli_epi32(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm256_max_epi32(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm256_min_epi32(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm256_blendv_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {return _mm256_cmpgt_epi32(u, v);}
	regist_v cmp_ge(regist_v u, regist_v v) {return _mm256_or_si256
	    (_mm256_cmpgt_epi32(u, v), _mm256_cmpeq_epi32(u, v));
	}
	regist_v cmp_eq(regist_v u, regist_v v) {return _mm256_cmpeq_epi32(u, v);}
	regist_v bit_and(regist_v u, regist_v v) {return _mm256_and_si256(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm256_or_si256(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm256_xor_si256(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm256_andnot_si256(v, u);}
	int	all_zero(regist_v v) {return _mm256_testnzc_si256(v, v);}
};

template <typename regist_v>
struct Simd_functions<float, 8, regist_v> {
	regist_v clear() {return _mm256_setzero_ps();}
	regist_v splat(const float f) {return _mm256_set1_ps(f);}
	regist_v load(const float* a) {
	    return _mm256_loadu_ps(a);
	}
	void	store(float* a, regist_v v) {_mm256_storeu_ps(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm256_add_ps(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm256_sub_ps(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm256_mul_ps(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm256_max_ps(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm256_min_ps(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm256_blendv_ps(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm256_cmp_ps(u, v, 0x1e);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm256_cmp_ps(u, v, 0x1d);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm256_cmp_ps(u, v, 0x00);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm256_and_si256(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm256_or_si256(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm256_xor_si256(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm256_andnot_si256(v, u);}
};

#endif	// __AVX2__

#if __AVX512BW____

template <typename regist_v>
struct Simd_functions<char, 64, regist_v> {
	regist_v clear() {return _mm512_setzero_si512();}
	regist_v splat(const char i) {return _mm512_set1_epi8(i);}
	regist_v load(const char* a) {
	    return _mm512_loadu_epi8(a);
	}
	void	store(char* a, regist_v v) {_mm512_storeu_epi8(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm512_adds_epi8(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm512_subs_epi8(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm512_max_epi8(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm512_min_epi8(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm512_mask_blend_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 6);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 5);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 0);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm512_and_si512(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm512_or_si512(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm512_xor_si512(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(regist_v v) {return !(_mm512_reduce_or_epi32(v));}
};

template <typename regist_v>
struct Simd_functions<CHAR, 64, regist_v> {
	regist_v clear() {return _mm512_setzero_si512();}
	regist_v splat(const CHAR i) {return _mm512_set1_epi8(i);}
	regist_v load(const CHAR* a) {
	    return _mm512_loadu_epi8(a);
	}
	void	store(CHAR* a, regist_v v) {_mm512_storeu_epi8(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm512_adds_epi8(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm512_subs_epi8(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm512_max_epi8(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm512_min_epi8(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm512_mask_blend_epi8(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 6);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 5);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm512_cmp_epi8_mask(u, v, 0);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm512_and_si512(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm512_or_si512(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm512_xor_si512(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(regist_v v) {return !(_mm512_reduce_or_epi32(v));}
};

template <typename regist_v>
struct Simd_functions<short, 32, regist_v> {
	regist_v clear() {return _mm512_setzero_si512();}
	regist_v splat(const short i) {return _mm512_set1_epi16(i);}
	regist_v load(const short* a) {
	    return _mm512_loadu_epi16(a);
	}
	void	store(short* a, regist_v v) {_mm512_storeu_epi16(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm512_adds_epi16(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm512_subs_epi16(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm512_mullo_epi16(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm512_slli_epi16(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm512_srai_epi16(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm512_max_epi16(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm512_min_epi16(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm512_mask_blend_epi16(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 6);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 5);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 0);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm512_and_si512(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm512_or_si512(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm512_xor_si512(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(regist_v v) {return !(_mm512_reduce_or_epi32(v));}
	regist_v cast16to8(regist_v v) {
	    regist_v	b_v = load(b64_a);
	    return _m512i_shuffle_epi8(v, b_v);
	}
};

template <typename regist_v>
struct Simd_functions<SHORT, 32, regist_v> {
	regist_v clear() {return _mm512_setzero_si512();}
	regist_v splat(const SHORT i) {return _mm512_set1_epi16(i);}
	regist_v load(const SHORT* a) {
	    return _mm512_loadu_epi16(a);
	}
	void	store(SHORT* a, regist_v v) {_mm512_storeu_epi16(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm512_adds_epi16(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm512_subs_epi16(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm512_mullo_epi16(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm512_slli_epi16(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm512_srli_epi16(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm512_max_epi16(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm512_min_epi16(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm512_mask_blend_epi16(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 6);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 5);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm512_cmp_epi16_mask(u, v, 0);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm512_and_si512(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm512_or_si512(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm512_xor_si512(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(regist_v v) {return !(_mm512_reduce_or_epi32(v));}
	regist_v cast16to8(regist_v v) {
	    regist_v	b_v = load(b64_a);
	    return _m512i_shuffle_epi8(v, b_v);
	}
};

template <typename regist_v>
struct Simd_functions<int, 16, regist_v> {
	regist_v clear() {return _mm512_setzero_si512();}
	regist_v splat(const int i) {return _mm512_set1_epi32(i);}
	regist_v load(const int* a) {return _mm512_loadu_epi32(a);}
	void	store(int* a, regist_v v) {_mm512_storeu_epi32(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm512_add_epi32(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm512_sub_epi32(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm512_mullo_epi32(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm512_slli_epi32(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm512_srai_epi32(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm512_max_epi32(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm512_min_epi32(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm512_mask_blend_epi32(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 6);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 5);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 0);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm512_and_si512(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm512_or_si512(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm512_xor_si512(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(regist_v v) {return !(_mm512_reduce_or_epi32(v));}
};

template <typename regist_v>
struct Simd_functions<INT, 16, regist_v> {
	regist_v clear() {return _mm512_setzero_si512();}
	regist_v splat(const INT i) {return _mm512_set1_epi32(i);}
	regist_v load(const INT* a) {return _mm512_loadu_epi32(a);}
	void	store(INT* a, regist_v v) {_mm512_storeu_epi32(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm512_add_epi32(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm512_sub_epi32(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm512_mullo_epi32(u, v);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm512_max_epi32(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm512_min_epi32(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm512_mask_blend_epi32(v, u, m);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 6);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 5);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm512_cmp_epi32_mask(u, v, 0);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm512_and_si512(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm512_or_si512(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm512_xor_si512(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(regist_v v) {return !(_mm512_reduce_or_epi32(v));}
};

template <typename regist_v>
struct Simd_functions<float, 16, regist_v> {
	regist_v clear() {return _mm512_setzero_ps();}
	regist_v splat(const float f) {return _mm512_set1_ps(f);}
	regist_v load(const float* a) {return _mm512_loadu_ps(a);}
	void	store(float* a, regist_v v) {_mm512_storeu_ps(a, v);}
	regist_v add(regist_v u, regist_v v) {
	    return _mm512_add_ps(u, v);
	}
	regist_v sub(regist_v u, regist_v v) {
	    return _mm512_sub_ps(u, v);
	}
	regist_v mul(regist_v u, regist_v v) {
	    return _mm512_mul_ps(u, v);
	}
	regist_v shiftl(regist_v u, const int i) {
	    return _mm512_slli_epi32(u, i);
	}
	regist_v shiftr(regist_v u, const int i) {
	    return _mm512_srli_epi32(u, i);
	}
	regist_v max(regist_v u, regist_v v) {
	    return _mm512_max_ps(u, v);
	}
	regist_v min(regist_v u, regist_v v) {
	    return _mm512_min_ps(u, v);
	}
	regist_v blend(regist_v u, regist_v v, regist_v m) {
	    return _mm512_mask_blend_ps(m, v, u);
	}
	regist_v cmp_gt(regist_v u, regist_v v) {
	    return _mm512_cmp_ps_mask(u, v, 0x1e);
	}
	regist_v cmp_ge(regist_v u, regist_v v) {
	    return _mm512_cmp_ps_mask(u, v, 0x1d);
	}
	regist_v cmp_eq(regist_v u, regist_v v) {
	    return _mm512_cmp_ps_mask(u, v, 0x00);
	}
	regist_v bit_and(regist_v u, regist_v v) {return _mm512_and_si512(u, v);}
	regist_v bit_or(regist_v u, regist_v v) {return _mm512_or_si512(u, v);}
	regist_v bit_xor(regist_v u, regist_v v) {return _mm512_xor_si512(u, v);}
	regist_v bit_andnot(regist_v u, regist_v v) {return _mm512_andnot_si512(v, u);}
	int	all_zero(regist_v v) {return !(_mm512_reduce_or_epi32(v));}
};
#endif	//__AVX512BW__

/*************************************************************************
	Array-wise functions using elementary simd functions
*************************************************************************/

template <typename var_t, int Nelem, typename regist_v>
struct  ArrayWiseSimdFunc :
	public Simd_functions<var_t, Nelem, regist_v> {
	void	vecclear(var_t* dst, const size_t n);
	var_t*	veccopy(var_t* dst, const var_t* src, const size_t n);
	void	vecset(var_t* dst, const var_t& c, const size_t n);
	void	vecadd_c(var_t* dst, const var_t& c, const size_t n);
	var_t	vecmax(const var_t* data, const size_t n);
	var_t	vecmin(const var_t* data, const size_t n);
};

template <typename var_t, int Nelem, typename regist_v>
void ArrayWiseSimdFunc<var_t, Nelem, regist_v>::
vecclear(var_t* d, const size_t n)
{
	size_t	nn = n / Nelem * Nelem;
	if (!nn) {
	    vclear(d, n);
	    return;
	}
	regist_v	z_v = this->clear();
	for (var_t* e = d + nn; d < e; d += Nelem) this->store(d, z_v);
	nn = n - nn;
	if (nn) this->store(d - Nelem + nn, z_v);
}

template <typename var_t, int Nelem, typename regist_v>
var_t* ArrayWiseSimdFunc<var_t, Nelem, regist_v>::
veccopy(var_t* d, const var_t* s, const size_t n)
{
	size_t	nn = n / Nelem * Nelem;
	if (!nn) return (vcopy(d, s, n));
	var_t*	r = d;
	for (const var_t* e = s + nn; s < e; d += Nelem, s += Nelem) {
	    regist_v	v_v = this->load(s);
	    this->store(d, v_v);
	}
	nn = n - nn;
	if (nn) {
	    regist_v	v_v = this->load(s - Nelem + nn);
	    this->store(d - Nelem + nn, v_v);
	}
	return (r);
}

template <typename var_t, int Nelem, typename regist_v>
void ArrayWiseSimdFunc<var_t, Nelem, regist_v>::
vecset(var_t* d, const var_t& c, const size_t n)
{
	size_t	nn = n / Nelem * Nelem;
	if (!nn) {
	    vset(d, c, n);
	    return;
	}
	regist_v	c_v = this->splat(c);
	for (var_t* e = d + nn; d < e; d += Nelem) this->store(d, c_v);
	nn = n - nn;
	if (nn) this->store(d - Nelem + nn, c_v);
}

template <typename var_t, int Nelem, typename regist_v>
void ArrayWiseSimdFunc<var_t, Nelem, regist_v>::
vecadd_c(var_t* d, const var_t& c, const size_t n)
{
	size_t	nn = n / Nelem * Nelem;
regist_v    c_v = this->splat(c);
	for (var_t* e = d + nn; d < e; d += Nelem) {
regist_v    v_v = this->add(this->load(d), c_v);
	    this->store(d, v_v);
	}
	nn = n - nn;
	while (nn--) *d++ += c;
}

template <typename var_t, int Nelem, typename regist_v>
var_t ArrayWiseSimdFunc<var_t, Nelem, regist_v>::
vecmax(const var_t* d, const size_t n)
{
	size_t	nn = n / Nelem * Nelem;
	var_t	buf[Nelem + Nelem / 2];
regist_v    a_v = this->splat(d[0]);
	this->store(buf, a_v); 	// d[0]..d[0]
	this->store(buf + Nelem / 2, a_v); 
	if (n > nn) vcopy(buf, d + nn, n - nn);
	a_v = this->load(buf);	// d[nn]..d[n-1]

	for (const var_t* e = d + nn; d < e; d += Nelem) {
regist_v    b_v = this->load(d);
	    a_v = this->max(a_v, b_v);
	}
	this->store(buf, a_v); 

// O(log_2(Nelem)) algorithm
	for (int k = Nelem; (k >>= 1); ) {
regist_v    b_v = this->load(buf + k);
	    a_v = this->max(a_v, b_v);
	    this->store(buf, a_v);
	}
	return (buf[0]);
}

template <typename var_t, int Nelem, typename regist_v>
var_t ArrayWiseSimdFunc<var_t, Nelem, regist_v>::
vecmin(const var_t* d, const size_t n)
{
	size_t	nn = n / Nelem * Nelem;
	var_t	buf[Nelem + Nelem / 2];
regist_v    a_v = this->splat(d[0]);
	this->store(buf, a_v); 
	this->store(buf + Nelem / 2, a_v); 
	if (n > nn) vcopy(buf, d + nn, n - nn);
	a_v = this->load(buf);

	for (const var_t* e = d + nn; d < e; d += Nelem) {
regist_v    b_v = this->load(d);
	    a_v = this->min(a_v, b_v);
	}
	this->store(buf, a_v); 

// O(log_2(Nelem)) algorithm
	for (int k = Nelem; (k >>= 1); ) {
regist_v    b_v = this->load(buf + k);
	    a_v = this->min(a_v, b_v);
	    this->store(buf, a_v);
	}
	return (buf[0]);
}

/*************************************************************************
	common public functions
*************************************************************************/

#if FVAL
template <typename X>
void vec_clear(X* ary, const size_t n)
{
#if !__SSE4_1__
	return (vclear(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	sf.vecclear(ary, n);
#endif	// !__SSE4_1__
}

template <typename X>
X* vec_copy(X* dst, const X* src, const size_t n)
{
	if (!dst) dst = new X[n];
#if !__SSE4_1__
	return (vclear(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	return sf.veccopy(dst, src, n);
#endif	// !__SSE4_1__
}

template <typename X>
void vec_set(X* ary, const X& c, const size_t n)
{
#if !__SSE4_1__
	return (vclear(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	sf.vecset(ary, c, n);
#endif	// !__SSE4_1__
}

template <typename X>
void vec_vecadd_c(X* ary, const X& c, const size_t n)
{
#if !__SSE4_1__
	return (vclear(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	sf.vecadd_c(ary, c, n);
#endif	// !__SSE4_1__
}

template <typename X>
X vec_max(const X* ary, const size_t n)
{
#if !__SSE4_1__
	return (vclear(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	return sf.vecmax(ary, n);
#endif	// !__SSE4_1__
}

template <typename X>
X vec_min(const X* ary, const size_t n)
{
#if !__SSE4_1__
	return (vclear(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	return sf.vecmin(ary, n);
#endif	// !__SSE4_1__
}

#else	// FVAL

template <typename X>
void vec_clear(X* ary, const size_t n)
{
#if !__SSE4_1__
	return (vclear(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512i;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256i;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128i;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	sf.vecclear(ary, n);
#endif	// !__SSE4_1__
}

template <typename X>
X* vec_copy(X* dst, const X* src, const size_t n)
{
	if (!dst) dst = new X[n];
#if !__SSE4_1__
	return (vcopy(dst, src, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512i;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256i;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128i;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	return sf.veccopy(dst, src, n);
#endif	// !__SSE4_1__
}

template <typename X>
void vec_set(X* ary, const X& c, const size_t n)
{
#if !__SSE4_1__
	return (vclear(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512i;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256i;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128i;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	sf.vecset(ary, c, n);
#endif	// !__SSE4_1__
}

template <typename X>
void vec_vecadd_c(X* ary, const X& c, const size_t n)
{
#if !__SSE4_1__
	return (vadd(ary, c, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512i;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256i;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128i;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	sf.vecadd_c(ary, c, n);
#endif	// !__SSE4_1__
}

template <typename X>
X vec_max(const X* ary, const size_t n)
{
#if !__SSE4_1__
	return (*vmax(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512i;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256i;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128i;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	return sf.vecmax(ary, n);
#endif	// !__SSE4_1__
}

template <typename X>
X vec_min(const X* ary, const size_t n)
{
#if !__SSE4_1__
	return (*vmin(ary, n));
#else
#if __AVX512BW__
const	int	nelem = 64 / sizeof(X);
	using	regist_v = __m512i;
#elif __AVX2__
const	int	nelem = 32 / sizeof(X);
	using	regist_v = __m256i;
#else
const	int	nelem = 16 / sizeof(X);
	using	regist_v = __m128i;
#endif
	ArrayWiseSimdFunc<X, nelem, regist_v> sf;
	return sf.vecmin(ary, n);
#endif	// !__SSE4_1__
}
#endif	// FVAL

#endif	// _SIMD_FUNCTIONS_
