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

// --- constants --------------------------------------------------------------o

#include <limits.h>

enum { DESC_WARN_MONO  = 1000000, // warn if tpsa can have 1e6 coefs or more
       DESC_MAX_ORD    = CHAR_BIT*sizeof(bit_t)-1, // max ord of a tpsa (=63)
       DESC_MAX_VAR    = 100000,  // max number of variables in a tpsa
       DESC_MAX_ARR    = 100,     // max number of simultaneous descriptors
       DESC_MAX_TMP    = 8,       // max number of temp. per thread in each desc
       TPSA_STRICT_NZ  = 1 // enforce strict nz (chk can be quadratic at level 3)
};

// --- types ------------------------------------------------------------------o

struct desc { // warning: must be identical to LuaJIT def (see mad_gtpsa.mad)
  int   id;          // index in list of registered descriptors
  int   nv, nk, nmv; // #vars (max 100000), #knobs, #mvars = nv-nk (main variables)
  ord_t mo, ko, to;  // max order, max order of vo[nmv:nv-1] and global order of truncation
  const ord_t *vo;   // orders of vars, vo[nv] (max order for each monomial variable)
              // end of compatibility with LuaJIT FFI

  int   uvo, nth;    // user provided vo, max #threads or 1
  ssz_t nc;          // number of coefs (max length of TPSA)

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

tpsa_t*  mad_tpsa_newd  (const desc_t *d, ord_t mo);
void     mad_tpsa_del   (const tpsa_t *t);

ctpsa_t* mad_ctpsa_newd (const desc_t *d, ord_t mo);
void     mad_ctpsa_del  (const ctpsa_t *t);

// --- TPSA sanity checks -----------------------------------------------------o

#define DEBUG 0

#if DEBUG > 0
#  define DBGTPSA(t) FUN(debug)(t,#t,__func__,__LINE__,0)
#else
#  define DBGTPSA(t)
#endif

// --- trace functions --------------------------------------------------------o

#if DEBUG > 1
#  define DBGFUN(a) printf(#a " %s:%d:\n",__func__,__LINE__)
#else
#  define DBGFUN(a)
#endif

// --- helpers ----------------------------------------------------------------o

static inline idx_t
hpoly_idx (idx_t ib, idx_t ia, ssz_t ia_size)
{
  return ib*ia_size + ia;
}

// --- macros for temporaries -------------------------------------------------o

#define TRC_TMPX(a) (void)func // a

#define GET_TMPX(t)       FUN(gettmp)(t, __func__)
#define REL_TMPX(t)       FUN(reltmp)(t, __func__)
#define GET_TMPC(t) mad_ctpsa_gettmpt(t, __func__)
#define REL_TMPC(t) mad_ctpsa_reltmp (t, __func__)
#define GET_TMPR(t)  mad_tpsa_gettmpt(t, __func__)
#define REL_TMPR(t)  mad_tpsa_reltmp (t, __func__)

// --- end --------------------------------------------------------------------o

#endif // MAD_DESC_PRIV_H
