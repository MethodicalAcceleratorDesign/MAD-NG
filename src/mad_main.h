#ifndef MAD_MAIN_H
#define MAD_MAIN_H

/*
 o----------------------------------------------------------------------------o
 |
 | Main module interface
 |
 | Methodical Accelerator Design - Copyright CERN 2015
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
  - provide some global and portable definitions commonly used in MAD

  Information:
  - formal parameters ending with an underscope can be null (optional).

 o----------------------------------------------------------------------------o
 */

#include "mad.h"

// --- interface -------------------------------------------------------------o

void mad_lua_setloc  (int level);
void mad_fatal_extrn (void) __attribute__((noreturn));

// ---------------------------------------------------------------------------o

#endif // MAD_MAIN_H
