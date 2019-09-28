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
// all arguments can be zero or null for default except nv>0 and nk>=0.

// mo = max(1,mo_)
const desc_t* mad_desc_newn(int nv, ord_t mo_);

// if nk == 0, same as newn, otherwise
// mo = max(1, mo_)
// vo = vo_ ? min(mo,vo_) : mo
// ko = ko_ ? min(mo,ko_) : mo
const desc_t* mad_desc_newk(int nv, ord_t mo_, int nk, ord_t vo_, ord_t ko_);

// if vars_ == null, same as newk, otherwise
// to = ord(vars[0:nv-1])
// mo = min(to, max(mo_, max( vars[0:nv-1]     )))
// vo = min(mo, max(vo_, max( vars[0:nv-nk-1]  )))
// ko = min(mo, max(ko_, max( vars[nv-nk:nv-1] )))
const desc_t* mad_desc_newv(int nv, ord_t mo_, const ord_t vars_[nv],
                            int nk, ord_t vo_, ord_t ko_);

// -- dtor
void  mad_desc_del    (const desc_t *d);

// -- introspection
int   mad_desc_nv     (const desc_t *d);
int   mad_desc_nk     (const desc_t *d);
int   mad_desc_nmv    (const desc_t *d); // nv-nk
ord_t mad_desc_maxord (const desc_t *d);
ssz_t mad_desc_maxlen (const desc_t *d); // ordlen(maxord) == maxlen
ssz_t mad_desc_ordlen (const desc_t *d, ord_t mo);
ord_t mad_desc_gtrunc (const desc_t *d, ord_t to);

// global cleanup (warning: no GTSPA must still be in use!)
void  mad_desc_cleanup (void);

// --- end --------------------------------------------------------------------o

#endif // MAD_DESC_H
