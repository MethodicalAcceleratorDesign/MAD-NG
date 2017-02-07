#ifndef MAD_BIT_H
#define MAD_BIT_H

/*
 o----------------------------------------------------------------------------o
 |
 | Bit module interface
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
extern int   mad_bit_lowest  (bit_t b);
extern int   mad_bit_highest (bit_t b);

// --- implementation (private) ----------------------------------------------o

#include "mad_bit_priv.h"

// ---------------------------------------------------------------------------o

#endif // MAD_BIT_H
