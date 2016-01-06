#ifndef MAD_TPSA_IMPL_H
#define MAD_TPSA_IMPL_H

/*
 o----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module implementation
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

#include "mad.h"
#include "mad_bit.h"
#include "mad_mono.h"

struct tpsa { // warning: must be kept identical to LuaJIT definition (cmad.lua)
  struct tpsa_desc *desc;
  ord_t             lo, hi, mo; // lowest/highest used ord, trunc ord
  bit_t             nz;
  int               tmp;
  num_t             coef[];
};

// ---------------------------------------------------------------------------o

#endif // MAD_TPSA_IMPL_H
