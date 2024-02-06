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
#include <string.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

#define DEBUG_CONVERT 0

// --- local ------------------------------------------------------------------o

static inline void
check_same_desc (ssz_t sa, const T *ma[sa])
{
  assert(ma);
  FOR(i,1,sa)
    ensure(ma[i]->d == ma[i-1]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
check_compat (ssz_t sa, const T *ma[sa], const T *mb[sa], T *mc[sa])
{
  assert(ma && mc);
  ensure(sa>0, "invalid map sizes (zero or negative sizes)");
  check_same_desc(sa,ma);
  check_same_desc(sa,(const T**)mc);
  ensure(ma[0]->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");

  if (mb) {
    check_same_desc(sa,mb);
    ensure(ma[0]->d == mb[0]->d, "incompatibles GTPSA (descriptors differ)");
  }
}

static inline void
print_damap (ssz_t sa, const T *ma[sa], FILE *fp_)
{
  char nam[NAMSZ];

  if (!fp_) fp_ = stdout;

  FOR(i,sa) {
    strcpy(nam, ma[i]->nam);
    if (!nam[0]) snprintf(nam, sizeof(nam), "#%d", i+1);
    FUN(print)(ma[i], nam, 1e-15, 0, fp_);
  }

  (void)print_damap;
}

static inline void
fgrad (ssz_t sa, const T *ma[sa], const T *b, T *c, T *t[2])
{
  FUN(clear)(c);
  FOR(i,sa) {
    FUN(deriv)(b    , t[0], i+1 );
    FUN(mul)  (ma[i], t[0], t[1]);
    FUN(add)  (c    , t[1], c   );
  }
}

static inline void
liebra (ssz_t sa, const T *ma[sa], const T *mb[sa], T *mc[sa], T *t[3])
{
  FOR(i,sa) {
    fgrad(sa, mb, ma[i], mc[i], t);
    fgrad(sa, ma, mb[i], t[2] , t);
    FUN(sub)(t[2], mc[i], mc[i]);
  }
}

static inline num_t
mnrm (ssz_t sa, const T *ma[sa])
{
  num_t nrm = 0;
  FOR(i,sa) nrm += FUN(nrm)(ma[i]);
  return nrm;
}

static inline void
exppb (ssz_t sa, const T *ma[sa], const T *mb[sa], T *mc[sa], T *t[4])
{
  const int nmax = 100;
  const num_t nrm_min1 = 1e-9 , nrm_min2 = 100*DBL_EPSILON*sa;
//const num_t nrm_min1 = 1e-10, nrm_min2 =   4*DBL_EPSILON*sa;

  FOR(i,sa) {
    num_t nrm_ = INFINITY;
    log_t conv = FALSE;

    FUN(copy)(mb[i], t[0]);
    FUN(copy)(mb[i], mc[i]);

    idx_t n;
    for (n=1; n <= nmax; ++n) {
      if (n==nmax/4) warn("exppb: n=%d (slow convergence)", n);
      FUN(scl)(t[0], 1.0/n, t[1]);
      fgrad(sa, ma, t[1], t[0], &t[2]);
      FUN(add)(mc[i], t[0], mc[i]);

      // check for convergence (avoid oscillations around very small values)
      const num_t nrm = FUN(nrm)(t[0]);
      if (nrm <= nrm_min2 || (conv && nrm >= nrm_))
        break;

      // convergence looks ok, just refine
      if (nrm <= nrm_min1) conv = TRUE;

      nrm_ = nrm;
    }

    if (n > nmax)
      warn("exppb did not converged after %d iterations for variable %d",nmax,i);
  }
}

#define TC (const T**)

static inline void // see 2nd Etienne's book, ch11
logpb (ssz_t sa, const T *ma[sa], T *mc[sa], T *t[4+5*sa], num_t eps)
{
  const int nmax = 100;
  const num_t nrm_min1 = 1e-9 , nrm_min2 = 100*DBL_EPSILON*sa;
//const num_t nrm_min1 = 1e-10, nrm_min2 =   4*DBL_EPSILON*sa;
  num_t nrm = 1e4, nrm_ = INFINITY, nrm0 = mnrm(sa, ma);
  num_t epsone = eps ? eps : nrm0/1000;
  log_t conv = FALSE;

  // temporary damaps
  T **t0 = &t[4+0*sa];
  T **t1 = &t[4+1*sa];
  T **t2 = &t[4+2*sa];
  T **t3 = &t[4+3*sa];
  T **t4 = &t[4+4*sa];

  idx_t n;
  for (n=1; n <= nmax; ++n) {
    if (n==nmax/4) warn("logpb: n=%d (slow convergence)", n);

    FOR(i,sa) FUN(scl) (mc[i], -1, t1[i]);     // t1 = -mc
    exppb(sa, TC t1, ma, t0, t);               // t0 = exp(:-mc:) ma
    FOR(i,sa) FUN(seti)(t0[i], i+1, 1, -1);    // t0 = t0-Id

    if (nrm < epsone) {
      FOR(i,sa) {  // t2 = -0.5*fgrad(t0, t0_i)
        fgrad(sa, TC t0, t0[i], t2[i], &t[2]);
        FUN(scl)(t2[i], -0.5, t2[i]);
      }
      FOR(i,sa) {  // t3 = -0.5*fgrad(t2, t0_i) - 1/6*fgrad(t0, t2_i)
        fgrad(sa, TC t2, t0[i], t3[i], &t[2]);
        fgrad(sa, TC t0, t2[i], t4[i], &t[2]);
        FUN(axpbypc)(-0.5, t3[i], -1.0/6, t4[i], 0, t3[i]);
      }
      FOR(i,sa) {  // t0 = t0+t2+t3
        FUN(add)(t0[i], t2[i], t0[i]);
        FUN(add)(t0[i], t3[i], t0[i]);
      }

      liebra(sa, TC mc, TC t0, t1, t); // t1 = <mc, t0>
      liebra(sa, TC mc, TC t1, t2, t); // t2 = <mc, <mc, t0>>
      liebra(sa, TC t0, TC t1, t3, t); // t3 = <t0, <mc, t0>>
      liebra(sa, TC t0, TC t2, t4, t); // t4 = <t0, <mc,  <mc, t0>>>

      FOR(i,sa) { // t0 = t0 + 0.5 t1 + 1/12 (t2-t3) - 1/24 t4
        FUN(axpbypc)(1, t0[i],  1.0/ 2, t1[i], 0, t0[i]);
        FUN(axpbypc)(1, t0[i],  1.0/12, t2[i], 0, t0[i]);
        FUN(axpbypc)(1, t0[i], -1.0/12, t3[i], 0, t0[i]);
        FUN(axpbypc)(1, t0[i], -1.0/24, t4[i], 0, t0[i]);
      }
    }
    FOR(i,sa) FUN(add)(mc[i], t0[i], mc[i]); // mc = mc + t0

    // check for convergence
    nrm = mnrm(sa, TC t0)/nrm0;

    // avoid numerical oscillations around very small values
    if (nrm <= nrm_min2 || (conv && nrm >= nrm_))
      break;

    // assume convergence is ok, just refine
    if (nrm <= nrm_min1) conv = TRUE;

    nrm_ = nrm;
  }

  if (n > nmax)
    warn("logpb did not converged after %d iterations", nmax);
}

#undef TC

// --- public -----------------------------------------------------------------o

void // compute M x = exp(:f(x;0):) x (eq. 32, 33 & 38 and inverse)
FUN(exppb) (ssz_t sa, const T *ma[sa], const T *mb[sa], T *mc[sa])
{
  DBGFUN(->); assert(mb);
  check_compat(sa, ma, mb, mc);
  const D *d = ma[0]->d;

  // handle aliasing
  mad_alloc_tmp(T*, mc_, sa);
  FOR(i,sa) {
    DBGTPSA(ma[i]); DBGTPSA(mb[i]); DBGTPSA(mc[i]);
    mc_[i] = FUN(newd)(d, d->to);
  }

  // temporaries
  T *t[4];
  FOR(i,4) t[i] = FUN(newd)(d, d->to);

  // main call
  exppb(sa, ma, mb, mc_, t);

  // temporaries
  FOR(i,4) FUN(del)(t[i]);

  // copy back
  FOR(i,sa) {
    FUN(copy)(mc_[i], mc[i]);
    FUN(del )(mc_[i]);
    DBGTPSA(mc[i]);
  }
  mad_free_tmp(mc_);
  DBGFUN(<-);
}

void // compute log(M) x = exp(log(:f(x;0):)) x => f
FUN(logpb) (ssz_t sa, const T *ma[sa], const T *mb[sa], T *mc[sa])
{
  DBGFUN(->);
  check_compat(sa, ma, mb, mc);
  const D *d = ma[0]->d;

  // handle aliasing
  mad_alloc_tmp(T*, mc_, sa);
  FOR(i,sa) {
    DBGTPSA(ma[i]); DBGTPSA(mc[i]);
    mc_[i] = FUN(newd)(d, d->to);
  }

  // initial guess provided
  if (mb) FOR(i,sa) {
    DBGTPSA(mb[i]);
    FUN(copy)(mb[i], mc_[i]);
  }

  // temporaries: 4 tpsa + 5 damap
  const int nt = 4+5*sa;
  T *t[nt];
  FOR(i,nt) t[i] = FUN(newd)(d, d->to);

  // main call
  logpb(sa, ma, mc_, t, 0);

  // temporaries
  FOR(i,nt) FUN(del)(t[i]);

  // copy back
  FOR(i,sa) {
    FUN(copy)(mc_[i], mc[i]);
    FUN(del )(mc_[i]);
    DBGTPSA(mc[i]);
  }
  mad_free_tmp(mc_);
  DBGFUN(<-);
}

void
FUN(liebra) (ssz_t sa, const T *ma[sa], const T *mb[sa], T *mc[sa])
{
  DBGFUN(->); assert(mb);
  check_compat(sa, ma, mb, mc);
  const D *d = ma[0]->d;

  // handle aliasing
  mad_alloc_tmp(T*, mc_, sa);
  FOR(i,sa) {
    DBGTPSA(ma[i]); DBGTPSA(mb[i]); DBGTPSA(mc[i]);
    mc_[i] = FUN(newd)(d, d->to);
  }

  // temporaries: 3 tpsa
  T *t[3];
  FOR(i,3) t[i] = FUN(newd)(d, d->to);

  // main call
  liebra(sa, ma, mb, mc_, t);

  // temporaries
  FOR(i,3) FUN(del)(t[i]);

  // copy back
  FOR(i,sa) {
    FUN(copy)(mc_[i], mc[i]);
    FUN(del )(mc_[i]);
    DBGTPSA(mc[i]);
  }
  mad_free_tmp(mc_);
  DBGFUN(<-);
}

void
FUN(fgrad) (ssz_t sa, const T *ma[sa], const T *b, T *c)
{
  DBGFUN(->);
  assert(ma && b && c);
  check_same_desc(sa,(const T**)ma);
  ensure(ma[0]->d == b->d, "incompatibles GTPSA (descriptors differ)");
  ensure(ma[0]->d == c->d, "incompatibles GTPSA (descriptors differ)");
  const D *d = ma[0]->d;

  // temporaries: 3 tpsa
  T *t[2];
  FOR(i,2) t[i] = FUN(newd)(d, d->to);

  // main call
  fgrad(sa, ma, b, c, t);

  // temporaries
  FOR(i,2) FUN(del)(t[i]);

  DBGFUN(<-);
}

void // compute G(x;0) = -J grad.f(x;0) (eq. 34),
FUN(vec2fld) (ssz_t sc, const T *a, T *mc[sc]) // pbbra
{
  DBGFUN(->);
  assert(a && mc);
  check_same_desc(sc,(const T**)mc);
  ensure(a->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");
  const D *d = a->d;

  T *t = FUN(newd)(d, d->to);

  FOR(i,sc) {
    FUN(setvar)(t, 0, i+1, 0);
    FUN(poisbra)(a, t, mc[i], 0);
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
  const D *d = ma[0]->d;

  FUN(reset0)(c);

  T *t1 = FUN(newd)(d, d->to);
  T *t2 = FUN(newd)(d, d->to);

  FOR(i,sa) {
    idx_t iv = i & 1 ? i : i+2;

    FUN(setvar)(t2, 0, iv, 0); // q_i -> p_i monomial, p_i -> q_i monomial
    FUN(mul)(ma[i], t2, t1);   // integrate by monomial of "paired" canon. var.
    FUN(sclord)(t1, t1, TRUE, FALSE); // scale coefs by orders, i.e. integ. var.

    (i & 1 ? FUN(add) : FUN(sub))(c, t1, c); // \sum p_i - q_i to c
  }

  FUN(del)(t2);
  FUN(del)(t1);
  DBGFUN(<-);
}

num_t
FUN(mnrm) (ssz_t sa, const T *ma[sa])
{
  DBGFUN(->);
  assert(ma);

  num_t nrm = 0;
  FOR(i,sa) {
    DBGTPSA(ma[i]);
    nrm += FUN(nrm)(ma[i]);
  }

  DBGFUN(<-);
  return nrm;
}

void // convert maps to another maps using tpsa conversion.
FUN(mconv) (ssz_t sa, const T *ma[sa], ssz_t sc, T *mc[sc], ssz_t n, idx_t t2r_[n], int pb)
{
  DBGFUN(->);
  assert(ma && mc);

  if (!t2r_) {
    ssz_t nn = MIN(sa,sc);
    FOR(i,nn) FUN(convert)(ma[i], mc[i], 0,0,0);
    DBGFUN(<-); return;
  }

  FOR(i,MIN(n,sa)) {
    idx_t ii = t2r_[i];
    if (ii < 0) continue; // discard var
    ensure(0 <= ii && ii < sc, "translation index out of range 0 <= %d < %d", ii, sc);
    FUN(convert)(ma[i], mc[ii], n, t2r_, pb);
    int ss = (ii-i)%2 * pb;
    if (ss == -1) FUN(scl)(mc[ii], -1, mc[ii]);
#if DEBUG_CONVERT
    printf("cvt: % 2d -> % d x % 2d [pb=% d]\n", i, ss, ii, pb);
#endif
  }

  DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
