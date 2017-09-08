#ifndef MAD_AVX_H
#define MAD_AVX_H

/*
 o----------------------------------------------------------------------------o
 |
 | AVX & AVX2 interface
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
 */

#include <immintrin.h>

// --- macros ----------------------------------------------------------------o

// char
#define MAD_AVX_CSIZ 32
#define MAD_AVX_CMSK            (MAD_AVX_CSIZ-1)
#define MAD_AVX_CRND(n) ((n) & ~(MAD_AVX_CSIZ-1))
#define MAD_AVX_CMOD(n) ((n) &  (MAD_AVX_CSIZ-1))

// int
#define MAD_AVX_ISIZ 8
#define MAD_AVX_IMSK            (MAD_AVX_ISIZ-1)
#define MAD_AVX_IRND(n) ((n) & ~(MAD_AVX_ISIZ-1))
#define MAD_AVX_IMOD(n) ((n) &  (MAD_AVX_ISIZ-1))

// double
#define MAD_AVX_DSIZ 4
#define MAD_AVX_DMSK            (MAD_AVX_DSIZ-1)
#define MAD_AVX_DRND(n) ((n) & ~(MAD_AVX_DSIZ-1))
#define MAD_AVX_DMOD(n) ((n) &  (MAD_AVX_DSIZ-1))

// --- globals ---------------------------------------------------------------o

extern const unsigned char mad_avx_msk1[32][32];
extern const unsigned char mad_avx_msk2[32][32];

// ---------------------------------------------------------------------------o

#endif // MAD_AVX_H
