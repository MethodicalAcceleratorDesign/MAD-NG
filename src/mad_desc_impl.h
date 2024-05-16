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
       DESC_MAX_ORD    = 250,     // max ord of a tpsa
       DESC_MAX_VAR    = 100000,  // max number of variables in a tpsa
       DESC_MAX_ARR    = 250,     // max number of simultaneous descriptors
       DESC_MAX_TMP    = 8,       // max number of temp. per thread in each desc
};

#define TPSA_STRICT  1 // see calls to update
#define TPSA_DEBUG   0 // 0-1: print fname in/out, call mad_tpsa_debug, more I/O
#define DESC_DEBUG   0 // 0-3: print debug info during descriptor construction
#define DESC_USE_TMP 0 // 0: use new, 1: use TMP

// --- types ------------------------------------------------------------------o

struct desc_ { // warning: must be identical to LuaJIT def (see mad_gtpsa.mad)
  int   id;          // index in list of registered descriptors
  int   nn, nv, np;  // #variables, #parameters, nn=nv+np <= 100000
  ord_t mo, po, sh;  // max order of vars & params, shared with id
  const ord_t *no;   // orders of each vars & params, no[nn]
              // end of compatibility with LuaJIT FFI

  int   uno, nth;    // user provided no, max #threads or 1
  ssz_t nc;          // number of coefs (max length of TPSA)

  int   *shared;     // counter of shared desc (all tables below except prms)
  ord_t *monos,      // 'matrix' storing the monomials (sorted by var)
        *ords,       // order of each mono of To
        *prms,       // order of parameters in each mono of To (zero = no prms)
       **To,         // Table by orders -- pointers to monos, sorted by order
       **Tv,         // Table by vars   -- pointers to monos, sorted by vars
       **ocs;        // ocs[t,i] -> o; in mul, compute o on thread t; 3 <= o <= mo; terminated with 0

  idx_t *ord2idx,    // order to polynomial start index in To (i.e. in TPSA coef[])
        *tv2to,      // lookup tv->to
        *to2tv,      // lookup to->tv
        *H,          // indexing matrix in Tv
       **L,          // multiplication indexes: L[oa,ob]->L_ord; L_ord[ia,ib]->ic
      ***L_idx;      // L_idx[oa,ob]->[start] [split] [end] idxs in L

  size_t size;       // bytes used by tables

  // permanent temporaries per thread for internal use (not shared)
#if DESC_USE_TMP
   tpsa_t ** t;      // tmp for  tpsa
  ctpsa_t **ct;      // tmp for ctpsa
  idx_t *ti, *cti;   // index of tmp used
#endif
};

// --- interface --------------------------------------------------------------o

#define D desc_t

// --- TPSA sanity checks -----------------------------------------------------o

#if TPSA_DEBUG
#  define DBGTPSA(t) ((void)(mad_tpsa_dbga && FUN(debug)(t,#t,__func__,__LINE__,0)))
#else
#  define DBGTPSA(t)
#endif

// --- trace functions --------------------------------------------------------o

#if TPSA_DEBUG
#  define DBGFUN(a) ((void)(mad_tpsa_dbgf && printf(#a " %s:%d:\n",__func__,__LINE__)))
#else
#  define DBGFUN(a)
#endif

// --- helpers ----------------------------------------------------------------o

static inline idx_t
hpoly_idx (idx_t ib, idx_t ia, ssz_t ia_size)
{
  return ib*ia_size + ia;
}

#define IS_COMPAT(...) MKNAME(IS_COMPAT_,NARG(__VA_ARGS__))(__VA_ARGS__)
#define IS_COMPAT_2(t1,t2)          ((t1)->d->monos == (t2)->d->monos)
#define IS_COMPAT_3(t1,t2,t3)       (IS_COMPAT_2(t1,t3)    && IS_COMPAT_2(t2,t3))
#define IS_COMPAT_4(t1,t2,t3,t4)    (IS_COMPAT_3(t1,t2,t4) && IS_COMPAT_2(t3,t4))
#define IS_COMPAT_5(t1,t2,t3,t4,t5) (IS_COMPAT_3(t1,t2,t5) && IS_COMPAT_3(t3,t4,t5))

// --- macros for temporaries -------------------------------------------------o

#if DESC_USE_TMP
#define TRC_TMPX(a) (void)func // a

#define GET_TMPX(t)       FUN(gettmp)(t, __func__)
#define REL_TMPX(t)       FUN(reltmp)(t, __func__)
#define GET_TMPC(t) mad_ctpsa_gettmpt(t, __func__)
#define GET_TMPR(t)  mad_tpsa_gettmpt(t, __func__)
#define REL_TMPC(t) mad_ctpsa_reltmp (t, __func__)
#define REL_TMPR(t)  mad_tpsa_reltmp (t, __func__)
#else
#define GET_TMPX(t)          FUN(new)(          t, mad_tpsa_same)
#define REL_TMPX(t)          FUN(del)(          t)
#define GET_TMPC(t)     mad_ctpsa_new((ctpsa_t*)t, mad_tpsa_same)
#define GET_TMPR(t)      mad_tpsa_new( (tpsa_t*)t, mad_tpsa_same)
#define REL_TMPC(t)     mad_ctpsa_del(          t)
#define REL_TMPR(t)      mad_tpsa_del(          t)
#endif // DESC_USE_TMP

// --- end --------------------------------------------------------------------o

#endif // MAD_DESC_PRIV_H
