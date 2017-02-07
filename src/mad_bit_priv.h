#ifndef MAD_BIT_PRIV_H
#define MAD_BIT_PRIV_H

/*
 o----------------------------------------------------------------------------o
 |
 | Bit module private implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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

// http://graphics.stanford.edu/~seander/bithacks.html

static inline bit_t  __attribute__((const))
mad_bit_set (bit_t b, int n)
{
  return b | (1u << n);
}

static inline bit_t  __attribute__((const))
mad_bit_get (bit_t b, int n)
{
  return b & (1u << n);
}

static inline bit_t __attribute__((const))
mad_bit_clr (bit_t b, int n)
{
  return b & ~(1u << n);
}

static inline bit_t __attribute__((const))
mad_bit_add (bit_t a, bit_t b)
{
  return a | b;
}

static inline bit_t __attribute__((const))
mad_bit_trunc (bit_t b, int n)
{
  return b & ((2u << n) - 1);
}

int mad_bit_lowest  (bit_t b) __attribute__((const));
int mad_bit_highest (bit_t b) __attribute__((const));

// ---------------------------------------------------------------------------o

#endif // MAD_BIT_PRIV_H
