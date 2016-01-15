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

#include "mad_bit.h"
#include "mad_ctpsa.h"

// --- types -----------------------------------------------------------------o

struct ctpsa { // warning: must be kept identical to LuaJIT definition (cmad.lua)
  desc_t *desc;
  ord_t   lo, hi, mo; // lowest/highest used ord, trunc ord
  bit_t   nz;
  cnum_t  coef[];
};

// --- helpers ---------------------------------------------------------------o

#define T           ctpsa_t
#define NUM         cnum_t
#define FUN(name)   MKNAME(mad_ctpsa_,name)
#define PFX(name)   MKNAME(c,name)
#define VAL(num)    creal(num), cimag(num)
#define FMT         "%g%+gi"
#define SELECT(R,C) C

#define CNUM(a) (* (cnum_t*) & (num_t[2]) { MKNAME(a,_re), MKNAME(a,_im) })

#include <tgmath.h>
#include <complex.h>

// ---------------------------------------------------------------------------o

#endif // MAD_CTPSA_IMPL_H
