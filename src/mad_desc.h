#ifndef MAD_DESC_H
#define MAD_DESC_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Descriptor (TPSA) module interface
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
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

extern const  ord_t  mad_tpsa_default;
extern const  ord_t  mad_tpsa_same;
extern const desc_t *mad_desc_curr;

// --- interface --------------------------------------------------------------o

// -- ctor: new general interface
// mo = max(1,mo_)
const desc_t* mad_desc_newn(int nv, ord_t mo_);

// if nk == 0, same as mad_desc_newn, otherwise
// mo = max(1, mo_)
// ko = ko_ ? min(mo,ko_) : mo
const desc_t* mad_desc_newk(int nv, ord_t mo_, int nk, ord_t ko_);

// mo = max(vars[0:nv-1])
// ko = nk>0 ? min(mo, max(ko_, max( vars[nv-nk:nv-1] ))) : mo
const desc_t* mad_desc_newv(int nv, const ord_t vars[nv], int nk, ord_t ko_);

// -- dtor
void  mad_desc_del    (const desc_t *d);

// -- introspection
int   mad_desc_nvmok  (const desc_t *d, ord_t *mo_, int *nk_, ord_t *ko_);
ord_t mad_desc_getvar (const desc_t *d, int nv, ord_t vars_[nv]); // return mo
ssz_t mad_desc_maxlen (const desc_t *d); // ordlen(mo) == maxlen
ssz_t mad_desc_ordlen (const desc_t *d, ord_t mo);
ord_t mad_desc_gtrunc (const desc_t *d, ord_t to);

// global cleanup (warning: no GTSPA must still be in use!)
void  mad_desc_cleanup(void);

// --- end --------------------------------------------------------------------o

#endif // MAD_DESC_H
