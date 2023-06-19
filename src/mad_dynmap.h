#ifndef MAD_DYNMAP_H
#define MAD_DYNMAP_H

/*
 o-----------------------------------------------------------------------------o
 |
 | C interface to dynamic maps
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
  - Provide interface to dynamic maps in C.

 o-----------------------------------------------------------------------------o
 */

#include "mad_tpsa.h"

// --- interface --------------------------------------------------------------o

void mad_trk_strex_drift (void);
void mad_trk_strex_kick  (void);

// --- end --------------------------------------------------------------------o

#endif // MAD_DYNMAP_H
