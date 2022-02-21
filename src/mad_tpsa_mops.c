/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA map operation module implementation
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

#include <math.h>
#include <float.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

// --- local ------------------------------------------------------------------o

static inline void
check_same_desc (ssz_t sa, const T *ma[sa])
{
  assert(ma);
  for (idx_t i=1; i < sa; ++i)
    ensure(ma[i]->d == ma[i-1]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
check_compat (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa])
{
  assert(ma && mb && mc);
  ensure(sa>0 && sb>0, "invalid map sizes (zero or negative sizes)");
  ensure(sb == ma[0]->d->nv, "incompatibles GTPSA (number of map variables differ)");
  check_same_desc(sa,ma);
  check_same_desc(sb,mb);
  check_same_desc(sa,(const T**)mc);
  ensure(ma[0]->d == mb[0]->d, "incompatibles GTPSA (descriptors differ)");
  ensure(ma[0]->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
print_damap (ssz_t sa, const T *ma[sa], FILE *fp)
{
  char nam[12];
  for (ssz_t i = 0; i < sa; ++i) {
    snprintf(nam, sizeof(nam), "#%d", i+1);
    FUN(print)(ma[i], nam, 1e-15, 0, fp);
  }
  (void)print_damap;
}

#include <stdio.h>
extern int mad_trace_fortid;

static inline void
exppb1 (ssz_t sa, const T *ma[sa], const T *b, T *c, T *t[4], log_t inv, int n)
{
  const num_t nrm_min1 = 1e-10, nrm_min2 = 4*DBL_EPSILON*sa;
  const int imax = 100;
  num_t nrm0 = INFINITY;
  log_t conv = FALSE;

  FILE *fp = NULL, *fp2 = NULL;
  char nam[100];
  if (mad_trace_fortid > 0) {
    snprintf(nam, 100, "fort/fort_n.%d.dat", mad_trace_fortid+100);
    fp = fopen(nam, "a");
    assert(fp);

    snprintf(nam, 100, "fort/fort_n.%d.dat", mad_trace_fortid+200);
    fp2 = fopen(nam, "a");
    assert(fp2);
    fprintf(fp2, "\nvar(exp)=%d\n", n);
  }

  FUN(copy)(b, t[0]);                                    // b1=x
  FUN(copy)(b, c);                                       // b4=x

  idx_t i;
  for (i = 1; i <= imax; ++i) {

    FUN(scl)(t[0], 1.0/i, t[1]);                         // b2=coe*b1

    if (fp) {
      fprintf(fp2, "itr(*)=%d\nd(0)=\n", i);
      FUN(print) (t[1], "s2", 1e-16, FALSE, fp2);
    }

    FUN(clear)(t[0]);
    for (idx_t j = 0; j < sa; ++j) {                     // b1=h*b2
      FUN(deriv)(t[1], t[2], j+1);
      FUN(mul)(ma[j], t[2], t[3]);

      if (fp) {
        fprintf(fp2,"d(i)=%d\n", j+1);
        FUN(print)(ma[j], "s1.v(i)", 1e-16, FALSE, fp2);
        FUN(print)(t[0] , "s22 before add", 1e-16, FALSE, fp2);
        FUN(print)(t[2] , "s2.d.i", 1e-16, FALSE, fp2);
        FUN(print)(t[3] , "s1.v(i)*s2.d.i", 1e-16, FALSE, fp2);
      }

      (inv ? FUN(sub) : FUN(add))(t[0], t[3], t[0]);

      if (fp && !FUN(isnul)(t[0])) {
        snprintf(nam, sizeof(nam), "t[0].%d.%d.%d", n, i, j+1);
        FUN(print)(t[0], nam, 1e-16, FALSE, fp);
        FUN(print)(t[0], nam, 1e-16, FALSE, fp2);
      }
    }
    FUN(add)(t[0], c, c);                                // b3=b1+b4

    // check for convergence
    const num_t nrm = FUN(nrm)(t[0]);                    // r=full_abs(b1)

    // avoid numerical oscillations around very small values
    if (nrm <= nrm_min2 || (conv && nrm >= nrm0))
      break;

    // assume convergence is ok, just refine
    if (nrm <= nrm_min1) conv = TRUE;

    nrm0 = nrm;
  }

  if (i > imax)
    warn("exppb did not converged after %d iterations for variable %d", imax,n);

  if (fp) { fclose(fp); fclose(fp2); }
}

static inline void // TODO!!!
logpb1 (ssz_t sa, const T *ma[sa], const T *b, T *c, T *t[4], int n)
{
  const num_t nrm_min1 = 1e-10, nrm_min2 = 4*DBL_EPSILON*sa;
  const int imax = 100;
  num_t nrm0 = INFINITY;
  log_t conv = FALSE;

  FILE *fp = NULL, *fp2 = NULL;
  char nam[100];
  if (mad_trace_fortid > 0) {
    snprintf(nam, 100, "fort/fort_n.%d.dat", mad_trace_fortid+100);
    fp = fopen(nam, "a");
    assert(fp);

    snprintf(nam, 100, "fort/fort_n.%d.dat", mad_trace_fortid+200);
    fp2 = fopen(nam, "a");
    assert(fp2);
    fprintf(fp2, "\nvar(exp)=%d\n", n);
  }

  FUN(copy)(b, t[0]);                                    // b1=x
  FUN(copy)(b, c);                                       // b4=x

  idx_t i;
  for (i = 1; i <= imax; ++i) {

    FUN(scl)(t[0], 1.0/i, t[1]);                         // b2=coe*b1

    if (fp) {
      fprintf(fp2, "itr(*)=%d\nd(0)=\n", i);
      FUN(print) (t[1], "s2", 1e-16, FALSE, fp2);
    }

    FUN(clear)(t[0]);
    for (idx_t j = 0; j < sa; ++j) {                     // b1=h*b2
      FUN(deriv)(t[1], t[2], j+1);
      FUN(mul)(ma[j], t[2], t[3]);

      if (fp) {
        fprintf(fp2,"d(i)=%d\n", j+1);
        FUN(print)(ma[j], "s1.v(i)", 1e-16, FALSE, fp2);
        FUN(print)(t[0] , "s22 before add", 1e-16, FALSE, fp2);
        FUN(print)(t[2] , "s2.d.i", 1e-16, FALSE, fp2);
        FUN(print)(t[3] , "s1.v(i)*s2.d.i", 1e-16, FALSE, fp2);
      }

      FUN(add)(t[0], t[3], t[0]);

      if (fp && !FUN(isnul)(t[0])) {
        snprintf(nam, sizeof(nam), "t[0].%d.%d.%d", n, i, j+1);
        FUN(print)(t[0], nam, 1e-16, FALSE, fp);
        FUN(print)(t[0], nam, 1e-16, FALSE, fp2);
      }
    }
    FUN(add)(t[0], c, c);                                // b3=b1+b4

    // check for convergence
    const num_t nrm = FUN(nrm)(t[0]);                    // r=full_abs(b1)

    // avoid numerical oscillations around very small values
    if (nrm <= nrm_min2 || (conv && nrm >= nrm0))
      break;

    // assume convergence is ok, just refine
    if (nrm <= nrm_min1) conv = TRUE;

    nrm0 = nrm;
  }

  if (i > imax)
    warn("logpb did not converged after %d iterations for variable %d", imax,n);

  if (fp) { fclose(fp); fclose(fp2); }
}

// --- public -----------------------------------------------------------------o

void // compute M x = exp(:f(x;0):) x (eq. 32, 33 & 38 and inverse)
FUN(exppb) (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa], int inv)
{
  DBGFUN(->);
  ensure(inv == 1 || inv == -1, "invalid inv value, -1 or 1 expected, got %d", inv);
  check_compat(sa, ma, sb, mb, mc);

  // handle aliasing
  mad_alloc_tmp(T*, mc_, sa);
  for (idx_t ic = 0; ic < sa; ++ic) {
    DBGTPSA(ma[ic]); DBGTPSA(mc[ic]);
    mc_[ic] = FUN(new)(mc[ic], mad_tpsa_same);
  }

  for (idx_t ib = 0; ib < sb; ++ib) DBGTPSA(mb[ib]);

  // temporaries
  T *t[4];
  for (int i = 0; i < 4; ++i) t[i] = FUN(new)(mc[0], mad_tpsa_same);

  for (idx_t i = 0; i < sa; ++i)
    exppb1(sa, ma, mb[i], mc_[i], t, inv == -1, i);

  // temporaries
  for (int i = 0; i < 4; i++) FUN(del)(t[i]);

  // copy back
  for (idx_t ic = 0; ic < sa; ++ic) {
    FUN(copy)(mc_[ic], mc[ic]);
    FUN(del )(mc_[ic]);
    DBGTPSA(mc[ic]);
  }
  mad_free_tmp(mc_);
  DBGFUN(<-);
}

void // compute log(M) x = log(exp(:f(x;0):)) x => f
FUN(logpb) (ssz_t sa, const T *ma[sa], ssz_t sb, const T *mb[sb], T *mc[sa]) // TODO!!!
{
  DBGFUN(->);
  check_compat(sa, ma, sb, mb, mc);

  // handle aliasing
  mad_alloc_tmp(T*, mc_, sa);
  for (idx_t ic = 0; ic < sa; ++ic) {
    DBGTPSA(ma[ic]); DBGTPSA(mc[ic]);
    mc_[ic] = FUN(new)(mc[ic], mad_tpsa_same);
  }

  for (idx_t ib = 0; ib < sb; ++ib) DBGTPSA(mb[ib]);

  // temporaries
  T *t[4];
  for (int i = 0; i < 4; ++i) t[i] = FUN(new)(mc[0], mad_tpsa_same);

  for (idx_t i = 0; i < sa; ++i) {
    logpb1(sa, ma, mb[i], mc_[i], t, i);
  }

  // temporaries
  for (int i = 0; i < 4; i++) FUN(del)(t[i]);

  // copy back
  for (idx_t ic = 0; ic < sa; ++ic) {
    FUN(copy)(mc_[ic], mc[ic]);
    FUN(del )(mc_[ic]);
    DBGTPSA(mc[ic]);
  }
  mad_free_tmp(mc_);
  DBGFUN(<-);
}

void // compute G(x;0) = -J grad.f(x;0) (eq. 34),
FUN(vec2fld) (ssz_t sc, const T *a, T *mc[sc]) // pbbra
{
  DBGFUN(->);
  assert(a && mc);
  check_same_desc(sc,(const T**)mc);
  ensure(a->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");

  T *t = FUN(new)(a, mad_tpsa_same);

  for (idx_t ic = 0; ic < sc; ++ic) {
    FUN(setvar)(t, 0, ic+1, 0);
    FUN(poisson)(a, t, mc[ic], 0);
  }

  FUN(del)(t);
  DBGFUN(<-);
}

void // compute f(x;0) = \int_0^x J G(x';0) dx' = x^t J phi G(x;0) (eq. 34, 36 & 37)
FUN(fld2vec) (ssz_t sa, const T *ma[sa], T *c) // getpb
{
  DBGFUN(->);
  assert(ma && c);
  check_same_desc(sa, ma);
  ensure(ma[0]->d == c->d, "incompatibles GTPSA (descriptors differ)");

  FUN(reset0)(c);

  T *t1 = FUN(new)(c, mad_tpsa_same);
  T *t2 = FUN(new)(c, mad_tpsa_same);

  for (idx_t ia = 0; ia < sa; ++ia) {
    idx_t iv = ia & 1 ? ia : ia+2;

    FUN(setvar)(t2, 0, iv, 0); // q_i -> p_i monomial, p_i -> q_i monomial
    FUN(mul)(ma[ia], t2, t1);  // integrate by monomial of "paired" canon. var.
    FUN(sclord)(t1, t1, TRUE); // scale coefs by orders, i.e. integrate order

    (ia & 1 ? FUN(add) : FUN(sub))(c, t1, c); // \sum p_i - q_i to c
  }

  FUN(del)(t2);
  FUN(del)(t1);
  DBGFUN(<-);
}

void // convert maps to another maps using tpsa conversion.
FUN(mconv) (ssz_t sa, const T *ma[sa], ssz_t sc, T *mc[sc], ssz_t n, idx_t t2r_[n], int pb)
{
  DBGFUN(->);
  assert(ma && mc);

  if (!t2r_) {
    ssz_t nn = MIN(sa,sc);
    for (idx_t i=0; i < nn; ++i) FUN(convert)(ma[i], mc[i], 0,0,0);
    DBGFUN(<-); return;
  }

  for (idx_t i=0; i < n; ++i) {
    if (t2r_[i] < 0) continue; // discard vars
    idx_t ii = t2r_[i];
    ensure(0 <= ii && ii < sc, "translation index out of range 0 <= %d < %d", ii, sc);
    FUN(convert)(ma[i], mc[ii], n, t2r_, pb);
    int ss = (ii-i)%2 * pb;
    if (ss == -1) FUN(scl)(mc[ii], -1, mc[ii]);
#if DEBUG > 2
    printf("cvt: % 2d -> % d x % 2d [pb=% d]\n", i, ss, ii, pb);
#endif
  }

  DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
