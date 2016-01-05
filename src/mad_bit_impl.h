#ifndef MAD_BIT_IMPL_H
#define MAD_BIT_IMPL_H
#else
#error "implementation header, do not include this file directly"
#endif

/*
 o----------------------------------------------------------------------------o
 |
 | Bit module implementation
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

static inline bit_t
mad_bit_set (bit_t b, int n)
{
  return b | (1 << n);
}

static inline bit_t
mad_bit_get (bit_t b, int n)
{
  return b & (1 << n);
}

static inline bit_t
mad_bit_clr (bit_t b, int n)
{
  return b & ~(1 << n);
}

static inline bit_t
mad_bit_add (bit_t a, bit_t b)
{
  return a | b;
}

static inline bit_t
mad_bit_trunc (bit_t b, int n)
{
  return b & ((2 << n) - 1);
}

static inline bit_t
mad_bit_lowest (bit_t b)
{
// find the number of trailing zeros in 32-bit b
// http://stackoverflow.com/questions/757059/position-of-least-significant-bit-that-is-set
  static const bit_t tbl[32] = {
    0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
  };
  return b ? tbl[((b & -b) * 0x077CB531U) >> 27] : 32;
}

static inline bit_t
mad_bit_highest (bit_t b)
{
  bit_t pos = 0;
  while(b >>= 1) ++pos;
  return pos;
}

// ---------------------------------------------------------------------------o
