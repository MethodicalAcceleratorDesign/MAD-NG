#ifndef MAD_POLY_H
#define MAD_POLY_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Polygon module interface
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
  - wrappers around functions of polygons for LuaJIT

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"

// --- interface --------------------------------------------------------------o

// polygon contains a point? (winding number algo)
log_t mad_pol_inside (num_t px, num_t py, const num_t *vx, const num_t *vy, ssz_t n);

// ----------------------------------------------------------------------------o

#endif // MAD_POLY_H
