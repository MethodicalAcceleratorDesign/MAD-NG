#ifndef MAD_TPSA_PRIV_H
#define MAD_TPSA_PRIV_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module implementation (private)
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
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
#include "mad_tpsa.h"

// --- types ------------------------------------------------------------------o

struct tpsa { // warning: must be kept identical to LuaJIT definition (cmad.lua)
  desc_t *d;
  ord_t   lo, hi, mo; // lowest/highest used ord, max ord (allocated)
  bit_t   nz;
  num_t   coef[];
};

// --- helpers ----------------------------------------------------------------o

#ifndef MAD_TPSA_NOHELPER

#define T           tpsa_t
#define NUM         num_t
#define FUN(name)   MKNAME(mad_tpsa_,name)
#define PFX(name)   name
#define VAL(num)    num
#define FMT         " %+6.4lE"
#define SELECT(R,C) R

#endif

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_PRIV_H
