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

// ---------------------------------------------------------------------------o

#endif // MAD_BIT_H

/*\

 |
 |
 |
 |
 |

 P
 R
 I
 V
 A
 T
 E

 I
 M
 P
 L
 E
 M
 E
 N
 T
 A
 T
 I
 O
 N

 |
 |
 |
 |
 |

\*/

#ifndef MAD_BIT_IMPL_H
#define MAD_BIT_IMPL_H

// --- implementation (private) ----------------------------------------------o

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
  extern const unsigned char mad_bit_lowest_tbl_[32];
  return b ? mad_bit_lowest_tbl_[((b & -b) * 0x077CB531u) >> 27] : 32;
}

static inline int
mad_bit_highest (bit_t b)
{
  extern const unsigned char mad_bit_highest_tbl_[256];
  bit_t r = (mad_bit_highest_tbl_[ b        & 0xFF] << 24) | 
            (mad_bit_highest_tbl_[(b >>  8) & 0xFF] << 16) | 
            (mad_bit_highest_tbl_[(b >> 16) & 0xFF] <<  8) |
            (mad_bit_highest_tbl_[(b >> 24) & 0xFF]);
  return 31 - mad_bit_lowest(r);
}

// ---------------------------------------------------------------------------o

#endif // MAD_BIT_IMPL_H
