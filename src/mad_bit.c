/*
 o-----------------------------------------------------------------------------o
 |
 | Bit module implementation
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
*/

#include "mad_bit.h"

// --- default versions -------------------------------------------------------o

#ifndef __SSE2__

// http://graphics.stanford.edu/~seander/bithacks.html

static const unsigned char mad_bit_lowest_tbl_[32] = {
  0, 1, 28, 2, 29, 14, 24, 3, 30, 22, 20, 15, 25, 17, 4, 8,
  31, 27, 13, 23, 21, 19, 16, 7, 26, 12, 18, 6, 11, 5, 10, 9
};

#define R2(n)    n,     n + 2*64,     n + 1*64,     n + 3*64
#define R4(n) R2(n), R2(n + 2*16), R2(n + 1*16), R2(n + 3*16)
#define R6(n) R4(n), R4(n + 2*4 ), R4(n + 1*4 ), R4(n + 3*4 )

static const unsigned char mad_bit_highest_tbl_[256] = {
  R6(0), R6(2), R6(1), R6(3)
};

int
mad_bit_lowest (bit_t b)
{
  return b ? mad_bit_lowest_tbl_[((b & -b) * 0x077CB531u) >> 27] : 32;
}

int
mad_bit_highest (bit_t b)
{
  bit_t r = (mad_bit_highest_tbl_[ b        & 0xFF] << 24) |
            (mad_bit_highest_tbl_[(b >>  8) & 0xFF] << 16) |
            (mad_bit_highest_tbl_[(b >> 16) & 0xFF] <<  8) |
            (mad_bit_highest_tbl_[(b >> 24) & 0xFF]);
  return 31 - mad_bit_lowest(r);
}

#endif // __SSE2__

// --- end --------------------------------------------------------------------o

