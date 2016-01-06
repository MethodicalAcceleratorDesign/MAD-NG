#ifndef MAD_CTPSA_IMPL_H
#define MAD_CTPSA_IMPL_H

/*
 o----------------------------------------------------------------------------o
 |
 | Complex Truncated Power Series Algebra module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
*/

#include <complex.h>

#include "mad_bit.h"
#include "mad_ctpsa.h"

// --- types -----------------------------------------------------------------o

struct ctpsa { // warning: must be kept identical to LuaJIT definition (cmad.lua)
  struct tpsa_desc *desc;
  ord_t             lo, hi, mo; // lowest/highest used ord, trunc ord
  bit_t             nz;
  int               tmp;
  cnum_t            coef[];
};

#define T struct ctpsa
#define NUM cnum_t
#define FUN(name) MKNAME(mad_ctpsa_,name)
#define PFX(name) MKNAME(c,name)
#define FMT "%g%+gi"
#define VAL(num) real(num), imag(num) 

// ---------------------------------------------------------------------------o

#endif // MAD_CTPSA_IMPL_H
