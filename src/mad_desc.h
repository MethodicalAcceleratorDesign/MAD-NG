#ifndef MAD_DESC_H
#define MAD_DESC_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Descriptor (TPSA) module interface
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

  Purpose:
  - provide a full feathered Generalized TPSA package

  Information:
  - parameters ending with an underscope can be null.

  Errors:
  - TODO

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"
#include "mad_mono.h"

// --- types -----------------------------------------------------------------o

typedef struct desc desc_t;

// --- globals ---------------------------------------------------------------o

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;
extern       int   mad_tpsa_strict;

// --- interface -------------------------------------------------------------o

// ctors, dtor
desc_t* mad_desc_new (int nmv, const ord_t mvar_ords[nmv], str_t mvar_names_[nmv]);
desc_t* mad_desc_newv(int nmv, const ord_t mvar_ords[nmv], str_t mvar_names_[nmv],
                      int nv , const ord_t  var_ords[nv ], ord_t dk);
void    mad_desc_del (desc_t *d);

// introspection
int     mad_desc_maxsize (const desc_t *d);
ord_t   mad_desc_maxord  (const desc_t *d);
ord_t   mad_desc_gtrunc  (      desc_t *d, ord_t to);

// ---------------------------------------------------------------------------o

#endif // MAD_DESC_H
