#ifndef SIMD_DEFINE
#define SIMD_DEFINE

#ifndef __USE_GNU
#define __USE_GNU
#endif

#ifdef __USE_INTEL
#undef __USE_GNU
#include "immintrin.h"
#include "xmmintrin.h"
#include "emmintrin.h"
#include "pmmintrin.h"
#include "smmintrin.h"
#endif

#ifdef __USE_INTEL
typedef __m256d v4df;
typedef __m256  v8sf;
typedef __m128d v2df;
typedef __m128  v4sf;
typedef __m128i v4si;
// SSE
#define __builtin_prefetch(p,rw,i)                 _mm_prefetch(p,i)
#define __builtin_ia32_shufps(a,b,imm)             _mm_shuffle_ps(a,b,imm)
#define __builtin_ia32_cvtpd2ps(a)                 _mm_cvtpd_ps(a)
#define __builtin_ia32_minpd(a,b)                  _mm_min_pd(a,b)
#define __builtin_ia32_unpckhpd(a,b)               _mm_unpackhi_pd(a,b)
#define __builtin_ia32_unpcklpd(a,b)               _mm_unpacklo_pd(a,b)
#define __builtin_ia32_loaddqu(mem_addr)           _mm_loaddup_pd(mem_addr)
//#define __builtin_ia32_vec_ext_v2df(a,imm)         _mm_store_sd(a,imm)
// AVX
#define __builtin_ia32_rsqrtps256(a)               _mm256_rsqrt_ps(a)
#define __builtin_ia32_vinsertf128_ps256(a,b,imm)  _mm256_insertf128_ps(a,b,imm)
#define __builtin_ia32_cvtpd2ps256(a)              _mm256_cvtpd_ps(a)
#define __builtin_ia32_cvtps2pd256(a)              _mm256_cvtps_pd(a)
#define __builtin_ia32_movntps256(mem_addr,a)      _mm256_stream_ps(mem_addr,a)
#define __builtin_ia32_unpcklps256(a, b)           _mm256_unpacklo_ps(a,b)
#define __builtin_ia32_unpckhps256(a, b)           _mm256_unpackhi_ps(a,b)
#define __builtin_ia32_shufps256(a,b,imm)          _mm256_shuffle_ps(a,b,imm)
#define __builtin_ia32_vextractf128_ps256(a,imm)   _mm256_extractf128_ps(a,imm)
#define __builtin_ia32_haddpd256(a,b)              _mm256_hadd_pd(a,b)
#define __builtin_ia32_vextractf128_pd256(a,imm)   _mm256_extractf128_pd(a,imm)
#define __builtin_ia32_minpd256(a,b)               _mm256_min_pd(a,b)
//#define __builtin_ia32_movntdq(mem_addr,a)         _mm_stream_si128(mem_addr,a)

// SSE
v2df operator + (const v2df& a, const v2df& b) { return _mm_add_pd(a,b); }
v2df operator - (const v2df& a, const v2df& b) { return _mm_sub_pd(a,b); }
v2df operator * (const v2df& a, const v2df& b) { return _mm_mul_pd(a,b); }
v2df operator / (const v2df& a, const v2df& b) { return _mm_div_pd(a,b); }
v4sf operator + (const v4sf& a, const v4sf& b) { return _mm_add_ps(a,b); }
v4sf operator - (const v4sf& a, const v4sf& b) { return _mm_sub_ps(a,b); }
v4sf operator * (const v4sf& a, const v4sf& b) { return _mm_mul_ps(a,b); }
v4sf operator / (const v4sf& a, const v4sf& b) { return _mm_div_ps(a,b); }
v4si operator - (const v4si& a, const v4si& b) { return _mm_sub_epi32(a,b); }

// AVX
v4df operator + (const v4df& a, const v4df& b) { return _mm256_add_pd(a,b); }
v4df operator - (const v4df& a, const v4df& b) { return _mm256_sub_pd(a,b); }
v4df operator * (const v4df& a, const v4df& b) { return _mm256_mul_pd(a,b); }
v4df operator / (const v4df& a, const v4df& b) { return _mm256_div_pd(a,b); }
v4df& operator += (v4df &a, const v4df& b) { a = _mm256_add_pd(a,b); }
v8sf operator + (const v8sf& a, const v8sf& b) { return _mm256_add_ps(a,b); }
v8sf operator - (const v8sf& a, const v8sf& b) { return _mm256_sub_ps(a,b); }
v8sf operator * (const v8sf& a, const v8sf& b) { return _mm256_mul_ps(a,b); }
v8sf operator / (const v8sf& a, const v8sf& b) { return _mm256_div_ps(a,b); }
v8sf& operator += (v8sf &a, const v8sf& b) { a = _mm256_add_ps(a,b); }


#endif

#ifdef __USE_GNU
// SSE
typedef float  v4sf __attribute__((vector_size(16)));
typedef double v2df __attribute__((vector_size(16)));
//typedef int    v4si __attribute__((vector_size(16)));
//typedef long long v2di __attribute__ ((__vector_size__ (16)));
// AVX
typedef float  v8sf __attribute__((vector_size(32)));
typedef double v4df __attribute__((vector_size(32)));
#endif


#define REP4(x) {x,x,x,x}
#define REP8(x) {x,x,x,x,x,x,x,x}

#endif
