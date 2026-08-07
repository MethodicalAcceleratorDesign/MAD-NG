/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA division module implementation
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
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
#include <float.h>
#include <complex.h>

#include "mad_mem.h"
#include "mad_num.h"
#include "mad_vec.h"
#include "mad_mat.h"
#include "mad_tpsa_impl.h"
#include "mad_ctpsa_impl.h"

// --- locals -----------------------------------------------------------------o

#define MAD_TPSA_DIV_RCOND    1e-12
#define MAD_TPSA_DIVC_SVDMAX (512u*512u)
#define MAD_TPSA_DIVC_RATIO   1e6   // trigger cancellation if ||[b]_k|| > |b0|*RATIO
#define MAD_TPSA_DIVC_AUTO    0     // >0 enable auto divc in div if np>0
                                    // >1 enable auto divc in div, >2 enable trace

// Counters for div and div->divc calls (diagnostic, not thread-safe).
static ssz_t div_cnt   = 0;
static ssz_t divc_cnt  = 0;
static ssz_t divc_fail = 0;

void
FUN(divc_clrcnt) (void)
{
  div_cnt = divc_cnt = divc_fail = 0;
}

void
FUN(divc_getcnt) (ssz_t *cnt_, ssz_t *fail_)
{
  if (cnt_ ) *cnt_  = divc_cnt;
  if (fail_) *fail_ = divc_fail;
#if MAD_TPSA_DIVC_AUTO > 0
  printf("divc auto: div=%d, divc=%d, ok=%d, fail=%d\n",
            div_cnt, divc_cnt, divc_cnt-divc_fail, divc_fail);
#endif
}

// --- helpers ----------------------------------------------------------------o
// 1/b = (1/b0) * sum_{k=0}^{to} (-(b-b0)/b0)^k.
// ord_coef[k] = (-1/b0)^k for k>=1; ord_coef[0] is repurposed as b0 for divn.
// ----------------------------------------------------------------------------o

static inline void
div_taylor (const T *a, const T *b, T *c, ord_t n, const NUM ord_coef[n+1])
{
  assert(a && b && c && ord_coef && n >= 1);

  T *acp = GET_TMPX(c);
  T *bcp = GET_TMPX(c);  // b - b0
  T *tmp = GET_TMPX(c);

  FUN(copy)(a, acp);
  FUN(copy)(b, bcp);
  FUN(seti)(bcp, 0, 0, 0);            // bcp = b - b0

  FUN(divn)(acp, ord_coef[0], c);     // c = a/b0
  FUN(mul) (acp, bcp, tmp);
  FUN(acc) (tmp, ord_coef[1], c);

  if (n >= 2) {
    T *pow = GET_TMPX(c), *t;
    FUN(mul)(bcp, bcp, pow);
    FUN(mul)(acp, pow, tmp);
    FUN(acc)(tmp, ord_coef[2], c);

    for (ord_t i = 3; i <= n; ++i) {
      FUN(mul)(bcp, pow, tmp);
      FUN(mul)(acp, tmp, pow);
      FUN(acc)(pow, ord_coef[i], c);
      SWAP(pow, tmp, t);
    }
    if (n & 1) SWAP(pow, tmp, t);
    REL_TMPX(pow);
  }

  REL_TMPX(tmp); REL_TMPX(bcp); REL_TMPX(acp);
}

static inline void
div_ (const T *a, const T *b, T *c, NUM b0)
{
  ord_t to = c->mo;
  if (!to || FUN(isval)(b)) { FUN(divn)(a,b0,c);        return; }
  if (FUN(isval)(a))        { FUN(inv)(b,a->coef[0],c); return; }

  NUM ord_coef[to+1];
  ord_coef[0] = NUMF(inv)(b0);
  for (ord_t o = 1; o <= to; ++o)
    ord_coef[o] = -NUMF(div)(ord_coef[o-1],b0);

  ord_coef[0] = b0;
  div_taylor(a,b,c,to,ord_coef);
}

// --- cancel -----------------------------------------------------------------o
// Polynomial cancellation division: c = a/b, b0 may be zero or negligible.
//
// Assumes b = [b]_rb + [b]_{rb+1} + ... with rb = first order where
// ||[b]_rb||_inf > tol; a = b*c + r with ||r|| ~ 0.
// When rb > 0 (b0 negligible vs [b]_rb), applies L'Hopital in the TPSA
// monomial sense: uses [b]_rb as the effective pivot via SVD at each order.
//
// Homogeneous recurrence:
//   [b]_rb * [c]_n = [a - b * c_prefix]_{rb+n},   n = 0,..,to-rb
//
// Each order is an (nr x nc) linear system solved by SVD (minimum-norm).
// If [b]_rb is rank-deficient, the greedy per-order choice may not be the
// globally compatible quotient; the final residual check decides success.
// SVD cost is bounded by MAD_TPSA_DIVC_SVDMAX to prevent runaway as
// singularities amplification generally occurs at low parametric orders.
// System is scaled by 1/||[b]_rb||_inf to improve the SVD condition number.
//
// Returns TRUE iff ||a - b*c||_1 <= tol * max(||a||_1, 1); c unchanged on fail.
// ----------------------------------------------------------------------------o

// Infinity norm of a TPSA at a single homogeneous order.
static inline num_t
divc_nrm (const T *a, ord_t o)
{
  assert(a);
  if (!o) return fabs(a->coef[0]);
  if (!a->hi || o < a->lo || o > a->hi) return 0;

  const idx_t *o2i = a->d->ord2idx;
  ssz_t nr = o2i[o+1] - o2i[o];
  mad_alloc_tmp(NUM, v, nr);
  FUN(getv)(a, o2i[o], nr, v);

  num_t nrm = 0;
  FOR(i, nr) { num_t vi = fabs(v[i]); if (vi > nrm) nrm = vi; }
  mad_free_tmp(v);
  return nrm;
}

// First order of a with ||[a]_o||_inf > tol (absolute) or > |a0|*RATIO
// (amplification risk). Starts scanning at lo; lo=0 includes the scalar.
// Returns 0 if b0 is the safe pivot, rb>0 if a higher order dominates,
// -1 if a is negligible everywhere in [lo,hi].
static inline ssz_t
divc_fst (const T *a, ord_t lo, ord_t hi, num_t tol)
{
  assert(a);
  hi = MIN(hi, a->mo);

  num_t a0 = fabs(a->coef[0]);
  num_t th = MAX(tol, a0 * MAD_TPSA_DIVC_RATIO);

  if (!lo && a0 > tol) return 0;

  const idx_t *o2i = a->d->ord2idx;
  const ord_t *i2o = a->d->ords;
  lo = MAX(lo, 1, a->lo);
  hi = MIN(hi, a->hi);
  for (idx_t i = o2i[lo]; i < o2i[hi+1]; ++i)
    if (fabs(a->coef[i]) > th) return i2o[i];

  return a0 > tol ? 0 : -1;
}

// Cancellation solver: [b]_rb * [c]_n = [a - b*c_prefix]_{rb+n}.
// Called with rb>0 and rbn = ||[b]_rb||_inf.
static log_t
divc_ (const T *a, const T *b, T *c, ord_t rb, num_t rbn, num_t tol)
{
  assert(a && b && c && rb > 0 && rbn > 0);

  ord_t to  = c->mo;
  ord_t qto = to >= rb ? to - rb : 0;

  T *prf = GET_TMPX(c);
  T *tmp = GET_TMPX(c);
  T *res = GET_TMPX(c);

  const D     *d   = b->d;
  const idx_t *o2i = d->ord2idx;
  ssz_t        nn  = d->nn;

  for (ord_t n = 0; n <= qto; ++n) {
    // rhs = [a - b*prf]_{rb+n}; for n=0 prf=0 so rhs=[a]_rb directly
    if (!n) FUN(getord)(a, res, rb);
    else {
      FUN(mul)(b, prf, tmp);
      FUN(sub)(a, tmp, res);
      FUN(getord)(res, res, rb+n);
    }

    ssz_t nr = o2i[rb+n+1] - o2i[rb+n];
    ssz_t nc = o2i[n+1]    - o2i[n];
    if (!nr || !nc) continue;
    if ((size_t)nr * nc > MAD_TPSA_DIVC_SVDMAX) goto fail;

    mad_alloc_tmp(NUM, mat, nr*nc);
    mad_alloc_tmp(NUM, rhs, nr);
    mad_alloc_tmp(NUM, sol, nc);

    FUN(getv)(res, o2i[rb+n], nr, rhs);
    SELECT(mad_vec_divn, mad_cvec_divn)(rhs, rbn, rhs, nr);

    memset(mat, 0, sizeof(NUM)*nr*nc);
    ord_t mk[nn];
    FOR(j,nc) {
      const ord_t *mj = d->To[o2i[n]+j];
      FOR(i,nr) {
        const ord_t *mi = d->To[o2i[rb+n]+i];
        if (!mad_mono_le(nn, mj, mi)) continue;
        mad_mono_sub(nn, mi, mj, mk);
        idx_t ik = mad_desc_idxm(d, nn, mk);
        if (ik >= 0) mat[i*nc+j] = FUN(geti)(b, ik) / rbn;
      }
    }

    int rnk = SELECT(mad_mat_ssolve, mad_cmat_ssolve)
                    (mat, rhs, sol, nr, nc, 1, MAD_TPSA_DIV_RCOND, 0, NULL);
    if (rnk >= 0) FUN(setv)(prf, o2i[n], nc, sol);

    mad_free_tmp(sol); mad_free_tmp(rhs); mad_free_tmp(mat);
    if (rnk < 0) goto fail;
  }

  // Global residual check
  FUN(mul)(b, prf, tmp);
  FUN(sub)(a, tmp, res);
  if (FUN(nrm)(res) > tol * MAX(FUN(nrm)(a), 1)) goto fail;

  if (qto < to) FUN(cutord)(prf, prf, qto+1);
  FUN(copy)(prf, c);

  REL_TMPX(res); REL_TMPX(tmp); REL_TMPX(prf);
  return TRUE;

fail:
  REL_TMPX(res); REL_TMPX(tmp); REL_TMPX(prf);
  return FALSE;
}

// --- public interface -------------------------------------------------------o

void
FUN(div) (const T *a, const T *b, T *c)  // c = a/b, b0 != 0
{
  assert(a && b && c); DBGFUN(->);
  ensure(IS_COMPAT(a,b,c), "incompatibles GTPSA (descriptors differ)");

  NUM b0 = b->coef[0];
  ensure(b0 != 0, "invalid domain div(.,"FMT")", VAL(b0));

  // Detect amplification risk from higher-order terms (orders 1+).
  // If found, attempt cancellation; fall back to Taylor on failure.
#if MAD_TPSA_DIVC_AUTO > 0
#if MAD_TPSA_DIVC_AUTO < 2
  if (a->d->np > 0)
#endif
  {
    ++div_cnt;
    ssz_t rb = divc_fst(b,1,c->mo,mad_tpsa_eps);
    if (rb > 0) {
      ++divc_cnt;
      num_t rbn = divc_nrm(b,rb);
      if (divc_(a,b,c,rb,rbn,mad_tpsa_eps)) { DBGFUN(<-); return; }
      ++divc_fail;
      trace(2, "div: divc auto failed (b0="FMT"), falling back to Taylor", VAL(b0));
    }
  }
#endif

  div_(a,b,c,b0);
  DBGFUN(<-);
}

log_t
FUN(divc) (const T *a, const T *b, T *c, num_t tol_)  // c = a/b, b0 may be 0
{
  assert(a && b && c); DBGFUN(->);
  ensure(IS_COMPAT(a,b,c), "incompatibles GTPSA (descriptors differ)");

  num_t tol = tol_ > 0 ? tol_ : mad_tpsa_eps;

  // rb: first homogeneous order of b suitable as pivot (not amplifying).
  // rb=0: b0 dominates, use Taylor. rb>0: cancellation via [b]_rb.
  ssz_t rb = divc_fst(b,0,c->mo,tol);
  ensure(rb >= 0, "invalid domain divc: b is identically zero");

  if (!rb) {
    div_(a,b,c,b->coef[0]);    // b0 is the safe pivot
    DBGFUN(<-); return TRUE;
  }

  num_t rbn = divc_nrm(b,rb); // guaranteed > tol by divc_fst
  log_t ok = divc_(a,b,c,rb,rbn,tol);
  DBGFUN(<-); return ok;
}
