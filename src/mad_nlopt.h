#ifndef MAD_NLOPT_H
#define MAD_NLOPT_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Non Linear Optimization module interface
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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
  - interface to NLOpt library for LuaJIT

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"

void mad_nlopt (void);

#endif // MAD_NLOPT_H
