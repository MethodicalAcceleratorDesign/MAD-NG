#ifndef MAD_MONO_SSE_H
#ifdef __SSE2__ // minimum requirement
#include "mad_sse_avx.h"

#undef mono_add
#undef mono_order
#undef mono_leq

#define mono_add   mono_add_sse
#define mono_order mono_order_sse
#define mono_leq   mono_leq_sse

static inline void
mono_add_sse(int n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  __m128i ra, rb, rr, rm;
  int i=0, nn=SSE_CRND(n), nm=SSE_CMOD(n);

  for (; i < nn; i+=SSE_CSIZ) {
    ra = _mm_loadu_si128((__m128i*)&a[i]);
    rb = _mm_loadu_si128((__m128i*)&b[i]);
    rr = _mm_adds_epi8(ra,rb);
    _mm_storeu_si128((__m128i*)&r[i],rr);
  }

#if 1
  if (nm) {
    rm = _mm_load_si128 ((__m128i*)mad_sse_msk2[nm]);
    ra = _mm_loadu_si128((__m128i*)&a[i]);
    rb = _mm_loadu_si128((__m128i*)&b[i]);
    rr = _mm_adds_epi8(ra,rb);
    _mm_maskmoveu_si128(rr, rm, (char*)&r[i]);
  }
#else
  for (int j=0; j < nm; ++j)
    r[i+j] = a[i+j] + b[i+j];
#endif
}

static inline int
mono_order_sse(int n, const ord_t a[n])
{
  assert(a);
  __m128i ra, rs, rm, zero = _mm_setzero_si128();
  int i=0, s=0, nn=SSE_CRND(n), nm=SSE_CMOD(n);

  for (; i < nn; i+=SSE_CSIZ) {
    ra = _mm_loadu_si128((__m128i*)&a[i]);
    rs = _mm_sad_epu8(ra, zero);
    s += _mm_cvtsi128_si32(_mm_srli_si128(rs,SSE_CSIZ/2));
    s += _mm_cvtsi128_si32(rs);
  }

#if 1
 if (nm) {
    rm = _mm_load_si128((__m128i*)mad_sse_msk2[nm]);
    ra = _mm_and_si128(rm,_mm_loadu_si128((__m128i*)&a[i]));
    rs = _mm_sad_epu8(ra, zero);
    s += _mm_cvtsi128_si32(_mm_srli_si128(rs,SSE_CSIZ/2));
    s += _mm_cvtsi128_si32(rs);
 }
#else
  for (int j=0; j < nm; j++)
    s += a[i+j];
#endif
  return s;
}

static inline int
mono_leq_sse(int n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  __m128i ra, rb, rr, rm;
  int i=0, nn=SSE_CRND(n), nm=SSE_CMOD(n);

  for (; i < nn; i+=SSE_CSIZ) {
    ra = _mm_loadu_si128((__m128i*)&a[i]);
    rb = _mm_loadu_si128((__m128i*)&b[i]);
    rr = _mm_cmpgt_epi8(ra,rb);
    if (_mm_movemask_epi8(rr)) return 0;
  }

#if 1
  if (nm) {
    rm = _mm_load_si128((__m128i*)mad_sse_msk2[nm]);
    ra = _mm_and_si128(rm,_mm_loadu_si128((__m128i*)&a[i]));
    rb = _mm_and_si128(rm,_mm_loadu_si128((__m128i*)&b[i]));
    rr = _mm_cmpgt_epi8(ra,rb);
    if (_mm_movemask_epi8(rr)) return 0;
  }
#else
  for (int j=0; j < nm; j++)
    if (a[i+j] > b[i+j]) return 0;
#endif
  return 1;
}

#endif
#endif
