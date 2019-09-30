#ifndef MAD_BIT_H
#define MAD_BIT_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Bit module interface
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - provide simple functions to manipulate bits.

  Information:
  - all functions are inlined

 o-----------------------------------------------------------------------------o
*/

// --- types ------------------------------------------------------------------o

typedef unsigned int bit_t;

// --- interface --------------------------------------------------------------o

static bit_t mad_bit_set     (bit_t b, int n);
static bit_t mad_bit_get     (bit_t b, int n);
static bit_t mad_bit_clr     (bit_t b, int n);
static bit_t mad_bit_lcut    (bit_t b, int n);
static bit_t mad_bit_hcut    (bit_t b, int n);

#ifndef __SSE2__
extern int   mad_bit_lowest  (bit_t b);
extern int   mad_bit_highest (bit_t b);
#else
static int   mad_bit_lowest  (bit_t b);
static int   mad_bit_highest (bit_t b);
#endif

// ----------------------------------------------------------------------------o
// --- implementation (private) -----------------------------------------------o
// ----------------------------------------------------------------------------o

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
mad_bit_lcut (bit_t b, int n) // clear bits < n
{
  return b & ~((1u << n) - 1);
}

static inline bit_t __attribute__((const))
mad_bit_hcut (bit_t b, int n) // clear bits > n
{
  return b & ((2u << n) - 1);
}

#ifndef __SSE2__
// #warning "Default (<SSE2) selected"

int mad_bit_lowest  (bit_t b) __attribute__((const));
int mad_bit_highest (bit_t b) __attribute__((const));

// --- optimized versions -----------------------------------------------------o

#elif defined(__AVX2__)
// #warning "AVX2 selected"
#include "sse/mad_bit_avx2.tc"
#elif defined(__SSE2__)
// #warning "SSE selected"
#include "sse/mad_bit_sse2.tc"
#else
#error "unsupported architecture"
#endif // __SSE2__ || __AVX2__

// --- end --------------------------------------------------------------------o

#endif // MAD_BIT_H
