#ifndef MAD_DESC_H
#define MAD_DESC_H

/*
 o----------------------------------------------------------------------------o
 |
 | Descriptor (TPSA) module interface
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
  
  Purpose:
  - provide a full feathered Generalized TPSA package
 
  Information:
  - parameters ending with an underscope can be null.

  Errors:
  - TODO

 o----------------------------------------------------------------------------o
 */

#include "mad.h"
#include "mad_mono.h"

// --- types -----------------------------------------------------------------o

typedef struct desc desc_t;

// --- globals ---------------------------------------------------------------o

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;
extern       int   mad_tpsa_strict;

// --- interface -------------------------------------------------------------o

#define D desc_t

// ctors, dtor
D*    mad_desc_new (int nv, const ord_t var_ords[nv], const ord_t map_ords_[nv], str_t var_nam_[nv]);
D*    mad_desc_newk(int nv, const ord_t var_ords[nv], const ord_t map_ords_[nv], str_t var_nam_[nv],
                    int nk, const ord_t knb_ords[nk], ord_t dk); // knobs
void  mad_desc_del (D *d);

// introspection
int   mad_desc_maxsize (const D *d);
ord_t mad_desc_maxord  (const D *d);
ord_t mad_desc_gtrunc  (      D *d, ord_t to);

#undef D

// ---------------------------------------------------------------------------o

#endif // MAD_DESC_H
