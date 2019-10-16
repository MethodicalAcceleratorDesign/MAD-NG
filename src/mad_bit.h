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

#include <stdint.h>

// --- types ------------------------------------------------------------------o

typedef uint64_t bit_t;

// --- interface --------------------------------------------------------------o

// bit mask
static _Bool mad_bit_mtst    (bit_t b, bit_t m);
static bit_t mad_bit_mget    (bit_t b, bit_t m);
static bit_t mad_bit_mset    (bit_t b, bit_t m);
static bit_t mad_bit_mflp    (bit_t b, bit_t m);
static bit_t mad_bit_mclr    (bit_t b, bit_t m);

// bit nth
static _Bool mad_bit_tst     (bit_t b, int n);
static bit_t mad_bit_get     (bit_t b, int n);
static bit_t mad_bit_set     (bit_t b, int n);
static bit_t mad_bit_flp     (bit_t b, int n);
static bit_t mad_bit_clr     (bit_t b, int n);

// other
static bit_t mad_bit_lcut    (bit_t b, int n);
static bit_t mad_bit_hcut    (bit_t b, int n);
static int   mad_bit_lowest  (bit_t b);
static int   mad_bit_highest (bit_t b);

// ----------------------------------------------------------------------------o
// --- implementation (private) -----------------------------------------------o
// ----------------------------------------------------------------------------o

// #undef __AVX2__
// #undef __SSE2__

#ifndef __SSE2__
extern int mad_bit_lowest32  (uint32_t b) __attribute__((const));
extern int mad_bit_lowest64  (uint64_t b) __attribute__((const));
extern int mad_bit_highest32 (uint32_t b) __attribute__((const));
extern int mad_bit_highest64 (uint64_t b) __attribute__((const));
#else
static int mad_bit_lowest32  (uint32_t b) __attribute__((const));
static int mad_bit_lowest64  (uint64_t b) __attribute__((const));
static int mad_bit_highest32 (uint32_t b) __attribute__((const));
static int mad_bit_highest64 (uint64_t b) __attribute__((const));
#endif

static inline _Bool __attribute__((const))
mad_bit_mtst (bit_t b, bit_t m)
{
  return (b & m) != 0;
}

static inline bit_t  __attribute__((const))
mad_bit_mget (bit_t b, bit_t m)
{
  return b & m;
}

static inline bit_t  __attribute__((const))
mad_bit_mset (bit_t b, bit_t m)
{
  return b | m;
}

static inline bit_t __attribute__((const))
mad_bit_mflp (bit_t b, bit_t m)
{
  return b ^ m;
}

static inline bit_t __attribute__((const))
mad_bit_mclr (bit_t b, bit_t m)
{
  return b & ~m;
}

static inline _Bool __attribute__((const))
mad_bit_tst (bit_t b, int n)
{
  return mad_bit_mtst(b, 1ull << n);
}

static inline bit_t  __attribute__((const))
mad_bit_get (bit_t b, int n)
{
  return mad_bit_mget(b, 1ull << n);
}

static inline bit_t  __attribute__((const))
mad_bit_set (bit_t b, int n)
{
  return mad_bit_mset(b, 1ull << n);
}

static inline bit_t __attribute__((const))
mad_bit_flp (bit_t b, int n)
{
  return mad_bit_mflp(b, 1ull << n);
}

static inline bit_t __attribute__((const))
mad_bit_clr (bit_t b, int n)
{
  return mad_bit_mclr(b, 1ull << n);
}

static inline bit_t __attribute__((const))
mad_bit_lcut (bit_t b, int n) // clear bits < n
{
  return mad_bit_mclr(b, (1ull << n)-1);
}

static inline bit_t __attribute__((const))
mad_bit_hcut (bit_t b, int n) // clear bits > n
{
  return mad_bit_mget(b, (2ull << n)-1);
}

static inline int __attribute__((const))
mad_bit_lowest (bit_t b) // 0..64 (0x0 -> 64)
{
  return sizeof b == sizeof(uint32_t)
         ? mad_bit_lowest32(b) : mad_bit_lowest64(b);
}

static inline int __attribute__((const))
mad_bit_highest (bit_t b) // -1..63 (0x0 -> -1)
{
  return sizeof b == sizeof(uint32_t)
         ? mad_bit_highest32(b) : mad_bit_highest64(b);
}

// --- optimized versions -----------------------------------------------------o

#if defined(__AVX2__)
// #warning("AVX2 selected")
#include "sse/mad_bit_avx2.tc"
#elif defined(__SSE2__)
// #warning("SSE2 selected")
#include "sse/mad_bit_sse2.tc"
#endif // __SSE2__ || __AVX2__

// --- debug/print lcut, lowest, hcut, highest --------------------------------o

void mad_bit_check (void);

// --- end --------------------------------------------------------------------o

#endif // MAD_BIT_H
