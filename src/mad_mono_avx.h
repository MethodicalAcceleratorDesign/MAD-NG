#ifndef MAD_MONO_AVX_H
#ifdef __AVX2__ // minimum requirement
#include "mad_sse_avx.h"

/* TODO

#undef mono_add
#undef mono_order
#undef mono_leq

#define mono_add     mono_add_avx
#define mono_order   mono_order_avx
#define mono_leq     mono_leq_avx
*/

#endif
#endif
