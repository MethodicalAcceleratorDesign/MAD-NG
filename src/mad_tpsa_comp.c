/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA map composition module implementation
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

#include <string.h>

#include "mad_mem.h"
#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

#define DEBUG_COMPOSE 0
#define TC (const T**)

// --- helpers ----------------------------------------------------------------o

static inline void
check_same_desc (ssz_t sa, const T *ma[sa])
{
  FOR(i,1,sa)
    ensure(ma[i]->d == ma[i-1]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
check_compose (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa],
               log_t chk_sa)
{
  ensure(sa > 0 && sb > 0, "invalid map sizes (zero or negative sizes)");
  if (chk_sa) {
    ensure(sa <= sb          , "incompatibles damap #A > #B");
    ensure(sa <= ma[0]->d->nv, "incompatibles damap #A > NV(A)");
  }
  ensure(  sb <= ma[0]->d->nn, "incompatibles damap #B > NV(A)+NP(A)");
  check_same_desc(sa, ma);
  check_same_desc(sb, mb);
  check_same_desc(sa, TC mc);
  ensure(IS_COMPAT(*ma,*mb,*mc), "incompatibles GTPSA (descriptors differ)");
}

static inline log_t
is_aliased (const T *ma, ssz_t sa, const T *mb[sa])
{
  FOR(ia,sa) if (ma == mb[ia]) return TRUE;
  return FALSE;
}

static inline void
print_damap (ssz_t sa, const T *ma[sa], FILE *fp_)
{
  char nam[NAMSZ];

  if (!fp_) fp_ = stdout;

  FOR(ia,sa) {
    strcpy(nam, ma[ia]->nam);
    if (!nam[0]) snprintf(nam, sizeof(nam), "#%d", ia+1);
    FUN(print)(ma[ia], nam, -1, 0, fp_);
  }

  (void)print_damap;
}

// --- implementation ---------------------------------------------------------o

typedef struct {
  ssz_t sa, sb;
  ord_t hi_ord;
  ord_t mo_ord;
  log_t *required;
  const T **ma, **mb;
        T **mc, **ords;
  const D *d;
} cmpctx_t;

static inline void
compose_ord1 (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa])
{
  FOR(ia,sa) {
    FUN(seti)(mc[ia], 0, 0, ma[ia]->coef[0]);
    FOR(ib,sb) {
      NUM coef = FUN(geti)(ma[ia],ib+1);
      if (coef) FUN(acc)(mb[ib], coef, mc[ia]);
    }
  }
}

static inline void
compose_mono (int ib, idx_t idx, ord_t o, ord_t mono[], cmpctx_t *ctx)
{
  // ib  : current variable index (in mb)
  // idx : current monomial index
  // o   : current monomial order
  // mono: current monomial
  if (idx < 0 || !ctx->required[idx]) return;

  const D *d = ctx->d;
#if DEBUG_COMPOSE
  printf("compose: ib=%d, o=%d, req[% 3d]->", ib, o, idx);
  mad_mono_print(d->nn, mono, 0,0);
  printf("\n");
#endif

  // compute monomial tpsa from mb
  if (o > 0) FUN(mul)(ctx->ords[o-1], ctx->mb[ib], ctx->ords[o]);

  // compose monomial tpsa with ma
  FOR(ia,ctx->sa) {
    NUM coef = FUN(geti)(ctx->ma[ia],idx);
    if (coef) FUN(acc)(ctx->ords[o], coef, ctx->mc[ia]);
  }

  if (o < ctx->hi_ord)
    // continue for each var & prm in mb
    for(; ib < ctx->sb; ++ib) {
      mono[ib]++;
      idx = mad_desc_idxm(d, d->nn, mono);
      compose_mono(ib, idx, o+1, mono, ctx); // recursive call
      mono[ib]--;
    }
}

static inline void
init_required(ssz_t sa, const T *ma[sa], log_t required[], ord_t hi_ord)
{
  assert(ma && required);
  const D *d = ma[0]->d;

  // root is always required
  required[0] = 1;

  // primary nodes (non-zero coefs)
  FOR(ia,sa) {
    TPSA_SCAN(ma[ia]) if (ma[ia]->coef[i]) required[i] = 1;
  }

  // fathers of primary nodes
  ord_t mono[d->nn];
  idx_t j, father;
  for (ord_t o=hi_ord; o > 1; --o) {
    TPSA_SCAN(ma[0],o,o) {
      if (required[i]) {
        mad_mono_copy(d->nn, d->To[i], mono);
        for (j = d->nn-1; j >= 0 && !mono[j]; --j) ;
        mono[j]--;
        father = mad_desc_idxm(d, d->nn, mono);

#if DEBUG_COMPOSE
  printf("compini: ");
  mad_mono_print(d->nn, d->To[i], 0,0);
  printf("->");
  mad_mono_print(d->nn, mono, 0,0);
  printf(" req[% 3d]%c\n", father, required[father] ? '*' : ' ');
#endif
        required[father] = 1;
      }
    }
  }
}

static inline void
compose (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa],
         ord_t hi_ord, ord_t mo_ord)
{
  const D *d = ma[0]->d;

  ssz_t nc = mad_desc_maxlen(d, hi_ord);
  mad_alloc_tmp(log_t, required, nc);
  init_required(sa, ma, memset(required, 0, nc*sizeof *required), hi_ord);

  // initialization
  T *ords[hi_ord+1]; // one for each order [0,hi_ord]
  FOR(o,hi_ord+1) ords[o] = FUN(newd)(d, mo_ord);
  FUN(seti)(ords[0],0,0,1);
  FOR(ia,sa) FUN(clear)(mc[ia]);

  cmpctx_t ctx = { .d=d, .sa=sa, .sb=sb, .ma=ma, .mb=mb, .mc=mc,
                   .ords=ords, .hi_ord=hi_ord, .mo_ord=mo_ord,
                   .required=required };

  // compose starting at root of tree: ib 0, idx 0, ord 0, mono 0
  ord_t mono[d->nn]; mad_mono_fill(d->nn, mono, 0);
  compose_mono(0, 0, 0, mono, &ctx);

  // cleanup
  FOR(o,hi_ord+1) FUN(del)(ords[o]);
  mad_free_tmp(required);
}

// --- public -----------------------------------------------------------------o

void //             sa <= nv                   sb <= nn
FUN(compose) (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa])
{
  assert(ma && mb && mc); DBGFUN(->);
  log_t chk_sa = TRUE;
  if (sa < 0) chk_sa = FALSE, sa = -sa; // special case sa > nv (not for damap)
  check_compose(sa, ma, sb, mb, mc, chk_sa);

  const D *d = ma[0]->d; (void)d;

  // handle aliasing
  log_t amc[sa];
  mad_alloc_tmp(T*, mc_, sa);
  FOR(ia,sa) {
    amc[ia] = mc[ia] == ma[ia] || is_aliased(mc[ia], sa, mb);
    mc_[ia] = amc[ia] ? FUN(new)(mc[ia], mad_tpsa_same) : FUN(reset0)(mc[ia]);
  }

  ord_t hi_ord = FUN(mord)(sa, TC ma, TRUE );
  ord_t mo_ord = FUN(mord)(sa, TC mc, FALSE);

#if DEBUG_COMPOSE
  printf("hi: %d\n", hi_ord);
  printf("mo: %d\n", mo_ord);
  printf("ma:\n"); print_damap(sa, ma, 0);
  printf("mb:\n"); print_damap(sb, mb, 0);
#endif

  if (hi_ord == 1) compose_ord1(sa,ma, sb,mb, mc);

#ifdef _OPENMP
  else if (d->pcomp && hi_ord >= 6 && d->ord2idx[hi_ord+1] >= d->pcomp) {
    #pragma omp parallel for
    FOR(ia,sa) {
#if DEBUG_COMPOSE
    printf("compose: thread no %d\n", omp_get_thread_num());
#endif
      compose(1,&ma[ia], sb,mb, &mc_[ia], ma[ia]->hi, mc_[ia]->mo);
    }
  }
  #endif // _OPENMP

  else compose(sa,ma, sb,mb, mc_, hi_ord, mo_ord);

  // copy back
  FOR(ia,sa) if (amc[ia]) {
    FUN(copy)(mc_[ia], mc[ia]);
    FUN(del )(mc_[ia]);
  }
  mad_free_tmp(mc_);
  DBGFUN(<-);
}

void
FUN(translate) (ssz_t sa, const T *ma[sa], ssz_t sb, const NUM tb[sb], T *mc[sa])
{
  assert(ma && tb && mc); DBGFUN(->);
  ensure(sa > 0 && sb > 0, "invalid map/vector sizes (zero or negative sizes)");
  ensure(sa <= sb        , "incompatibles map/vector #A > #B");

  // transform translation vector into damap of order 1
  mad_alloc_tmp(const T*, mb, sb);
  FOR(ib,sb) {
    T *t = FUN(newd)(ma[0]->d,1);
    FUN(setvar)(t, tb[ib], ib+1, 0);
    mb[ib] = t;
  }

  FUN(compose)(sa, ma, sb, mb, mc);

  // cleanup
  FOR(ib,sb) FUN(del)(mb[ib]);
  mad_free_tmp(mb);
  DBGFUN(<-);
}

void
FUN(eval) (ssz_t sa, const T *ma[sa], ssz_t sb, const NUM tb[sb], NUM tc[sa])
{
  assert(ma && tb && tc); DBGFUN(->);
  ensure(sa > 0 && sb > 0, "invalid map/vector sizes (zero or negative sizes)");
  ensure(sa <= sb        , "incompatibles map/vector #A > #B");

  // transform vectors into damap of order 0
  mad_alloc_tmp(const T*, mb, sb);
  mad_alloc_tmp(      T*, mc, sa);
  FOR(ib,sb) {
    T *t = FUN(newd)(ma[0]->d,0);
    FUN(setval)(t, tb[ib]);
    mb[ib] = t;
  }
  FOR(ic,sa) {
    T *t = FUN(newd)(ma[0]->d,0);
    FUN(setval)(t, tc[ic]);
    mc[ic] = t;
  }

  FUN(compose)(sa, ma, sb, mb, mc);

  // cleanup, save result
  FOR(ib,sb) FUN(del)(mb[ib]);
  FOR(ic,sa) {
    tc[ic] = FUN(geti)(mc[ic],0);
    FUN(del)(mc[ic]);
  }
  mad_free_tmp(mb);
  mad_free_tmp(mc);
  DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
