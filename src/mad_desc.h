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

#include "mad_mono.h"

// --- types ------------------------------------------------------------------o

typedef struct desc desc_t;

// --- globals ----------------------------------------------------------------o

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;

// --- interface --------------------------------------------------------------o

// ctors, dtor
desc_t* mad_desc_newn (int nmv, ord_t mvo); // nmv mvars of order mvo
desc_t* mad_desc_newm (int nmv, const ord_t mvar_ords[nmv]);
desc_t* mad_desc_newv (int nmv, const ord_t mvar_ords[nmv],
                       int nv , const ord_t  var_ords[nv ], ord_t dk); // knobs X-order
void    mad_desc_del  (desc_t *d);

// introspection
ord_t   mad_desc_maxord (const desc_t *d);
ssz_t   mad_desc_maxlen (const desc_t *d); // ordlen(maxord) == maxlen
ssz_t   mad_desc_ordlen (const desc_t *d, ord_t mo);
ord_t   mad_desc_gtrunc (      desc_t *d, ord_t to);

// --- end --------------------------------------------------------------------o

#endif // MAD_DESC_H
