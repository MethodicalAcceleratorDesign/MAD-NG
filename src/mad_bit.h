#ifndef MAD_BIT_H
#define MAD_BIT_H

/*
 o----------------------------------------------------------------------------o
 |
 | Bit module interface
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
  
  Purpose:
  - provide simple functions to manipulate bits.
 
  Information:
  - all functions are inlined

 o----------------------------------------------------------------------------o
 */

// --- types -----------------------------------------------------------------o

typedef unsigned int bit_t;

// --- interface -------------------------------------------------------------o

static bit_t mad_bit_set     (bit_t b, int n);
static bit_t mad_bit_get     (bit_t b, int n);
static bit_t mad_bit_clr     (bit_t b, int n);
static bit_t mad_bit_add     (bit_t a, bit_t b);
static bit_t mad_bit_trunc   (bit_t b, int n);
static int   mad_bit_lowest  (bit_t b);
static int   mad_bit_highest (bit_t b);

#endif // MAD_BIT_H























































// --- implementation (private) ----------------------------------------------o

#ifndef MAD_BIT_IMPL_H
#define MAD_BIT_IMPL_H

// http://graphics.stanford.edu/~seander/bithacks.html

static inline bit_t
mad_bit_set (bit_t b, int n)
{
  return b | (1u << n);
}

static inline bit_t
mad_bit_get (bit_t b, int n)
{
  return b & (1u << n);
}

static inline bit_t
mad_bit_clr (bit_t b, int n)
{
  return b & ~(1u << n);
}

static inline bit_t
mad_bit_add (bit_t a, bit_t b)
{
  return a | b;
}

static inline bit_t
mad_bit_trunc (bit_t b, int n)
{
  return b & ((2u << n) - 1);
}

static inline int
mad_bit_lowest (bit_t b)
{
  static const unsigned char tbl[32] = {
    0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
    31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
  };
  return b ? tbl[((b & -b) * 0x077CB531u) >> 27] : 32;
}

static inline int
mad_bit_highest (bit_t b)
{
  #define R2(n)    n,     n + 2*64,     n + 1*64,     n + 3*64
  #define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
  #define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )
  static const unsigned char tbl[256] = { R6(0), R6(2), R6(1), R6(3) };
  #undef  R2
  #undef  R4
  #undef  R6
  bit_t r = (tbl[ b        & 0xFF] << 24) | 
            (tbl[(b >>  8) & 0xFF] << 16) | 
            (tbl[(b >> 16) & 0xFF] <<  8) |
            (tbl[(b >> 24) & 0xFF]);
  return 31 - mad_bit_lowest(r);
}

// ---------------------------------------------------------------------------o

#endif // MAD_BIT_IMPL_H
