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
  - provide a full featured parametric Generalized TPSA package

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
// mo = max(1, mo)
const desc_t* mad_desc_newv(int nv, ord_t mo);

// if np == 0, same as mad_desc_newv, otherwise
// mo = max(1, mo)
// po = max(1, po_)
const desc_t* mad_desc_newvp(int nv, ord_t mo, int np_, ord_t po_);

// mo = max(mo , no[0 :nn-1]), nn = nv+np
// po = max(po_, no[nv:nn-1]), po <= mo
const desc_t* mad_desc_newvpo(int nv, ord_t mo, int np_, ord_t po_, const ord_t no_[nv+np_]);

// -- dtor
void  mad_desc_del    (const desc_t *d);

// -- introspection
int   mad_desc_getnv  (const desc_t *d, ord_t *mo_, int *np_, ord_t *po_); // return nv
ord_t mad_desc_maxord (const desc_t *d, int nn, ord_t no_[nn]); // return mo
ssz_t mad_desc_maxlen (const desc_t *d, ord_t mo);
ord_t mad_desc_gtrunc (const desc_t *d, ord_t to);

// -- indexes / monomials
log_t mad_desc_isvalids  (const desc_t *d,          ssz_t n,       str_t s    ); // string
log_t mad_desc_isvalidm  (const desc_t *d,          ssz_t n, const ord_t m [n]); // mono
log_t mad_desc_isvalidsm (const desc_t *d,          ssz_t n, const idx_t m [n]); // sparse mono
idx_t mad_desc_idxs      (const desc_t *d,          ssz_t n,       str_t s    ); // string
idx_t mad_desc_idxm      (const desc_t *d,          ssz_t n, const ord_t m [n]); // mono
idx_t mad_desc_idxsm     (const desc_t *d,          ssz_t n, const idx_t m [n]); // sparse mono
idx_t mad_desc_nxtbyvar  (const desc_t *d,          ssz_t n,       ord_t m [n]);
idx_t mad_desc_nxtbyord  (const desc_t *d,          ssz_t n,       ord_t m [n]);
ord_t mad_desc_mono      (const desc_t *d, idx_t i, ssz_t n,       ord_t m_[n]);

// for debugging
void  mad_desc_info      (const desc_t *d, FILE *fp_);

// global cleanup (warning: no GTSPA must still be in use!)
void  mad_desc_cleanup(void);

// --- end --------------------------------------------------------------------o

#endif // MAD_DESC_H
