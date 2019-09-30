#ifndef MAD_DESC_PRIV_H
#define MAD_DESC_PRIV_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Descriptor (TPSA) module implementation (private)
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
*/

#include "mad_log.h"
#include "mad_bit.h"
#include "mad_desc.h"
#include "mad_tpsa.h"
#include "mad_ctpsa.h"

// --- types ------------------------------------------------------------------o

enum { DESC_MAX_TMP = 6 };

struct desc { // warning: must be identical to LuaJIT def (see mad_gtpsa.mad)
  int   nv, nk, nmv; // #vars (max 100000), #knobs, #mvars = nv-nk (main variables)
  ord_t mo, ko, to;  // max order, max order of vars[nmv:nv-1] and global order of truncation
              // end of compatibility with LuaJIT FFI

  int   id, nth;     // index in list of registered descriptors, max #threads or 1
  int   uvars;       // user provided vars
  ssz_t nc;          // number of coefs (max length of TPSA)

  const ord_t *vars; // orders of vars[nv] (max order for each monomial variable)

  ord_t *monos,      // 'matrix' storing the monomials (sorted by var)
        *ords,       // order of each mono of To
       **To,         // Table by orders -- pointers to monos, sorted by order
       **Tv,         // Table by vars   -- pointers to monos, sorted by vars
       **ocs;        // ocs[t,i] -> o; in mul, compute o on thread t; 3 <= o <= mo; terminated with 0

  idx_t *ord2idx,    // order to polynomial start index in To (i.e. in TPSA coef[])
        *tv2to,      // lookup tv->to
        *to2tv,      // lookup to->tv
        *H,          // indexing matrix in Tv
       **L,          // multiplication indexes: L[oa,ob]->L_ord; L_ord[ia,ib]->ic
      ***L_idx;      // L_idx[oa,ob]->[start] [split] [end] idxs in L

  size_t size;       // bytes used by desc

  // permanent temporaries per thread for internal use
   tpsa_t ** t;      // tmp for  tpsa
  ctpsa_t **ct;      // tmp for ctpsa
  idx_t *ti, *cti;   // index of tmp used
};

// --- interface --------------------------------------------------------------o

#define D desc_t

ord_t    mad_desc_get_mono        (const D *d, ssz_t n,       ord_t m_[n], idx_t i);
idx_t    mad_desc_get_idx_s       (const D *d, ssz_t n,       str_t s    );
idx_t    mad_desc_get_idx_m       (const D *d, ssz_t n, const ord_t m [n]);
idx_t    mad_desc_get_idx_sm      (const D *d, ssz_t n, const idx_t m [n]);
int      mad_desc_mono_isvalid_s  (const D *d, ssz_t n,       str_t s    );
int      mad_desc_mono_isvalid_m  (const D *d, ssz_t n, const ord_t m [n]);
int      mad_desc_mono_isvalid_sm (const D *d, ssz_t n, const idx_t m [n]);
int      mad_desc_mono_nxtbyvar   (const D *d, ssz_t n,       ord_t m [n]);

tpsa_t*  mad_tpsa_newd  (const D *d, ord_t mo);
void     mad_tpsa_del   (const tpsa_t *t);

ctpsa_t* mad_ctpsa_newd (const D *d, ord_t mo);
void     mad_ctpsa_del  (const ctpsa_t *t);

// --- helpers ----------------------------------------------------------------o

static inline idx_t
hpoly_idx (idx_t ib, idx_t ia, idx_t ia_size)
{
  return ib*ia_size + ia;
}

#ifdef DEBUG
#define CHECK_VALIDITY(t) if (!FUN(is_valid)(t)) FUN(debug)(t,__func__,0)
#else
#define CHECK_VALIDITY(t)
#endif

// --- macros for temporaries -------------------------------------------------o

#define TRC_TMPX(a) (void)func // a

#define GET_TMPX(t)       FUN(gettmp)(t, __func__)
#define REL_TMPX(t)       FUN(reltmp)(t, __func__)
#define GET_TMPC(t) mad_ctpsa_gettmp (t, __func__)
#define REL_TMPC(t) mad_ctpsa_reltmp (t, __func__)
#define GET_TMPR(t) mad_ctpsa_gettmpr(t, __func__)

// --- end --------------------------------------------------------------------o

#endif // MAD_DESC_PRIV_H
