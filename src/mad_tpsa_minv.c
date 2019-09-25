/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA map inversion module implementation
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

#include <assert.h>

#include "mad_mem.h"
#include "mad_vec.h"
#include "mad_mat.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

// --- local ------------------------------------------------------------------o

static inline void
check_same_desc(ssz_t sa, const T *ma[sa])
{
  assert(ma);
  for (idx_t i=1; i < sa; ++i)
    ensure(ma[i]->d == ma[i-1]->d, "incompatibles GTPSA (descriptors differ)");
}

static inline void
check_minv(ssz_t sa, const T *ma[sa], T *mc[sa])
{
  ensure(sa == ma[0]->d->nmv, "non-square system"); // 'square' matrix, ignoring knobs
  check_same_desc(sa,ma);
  check_same_desc(sa,(const T**)mc);
  ensure(ma[0]->d == mc[0]->d, "incompatibles GTPSA (descriptors differ)");
}

/* GSL BASED REAL VERSION (NO COMPLEX)

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_blas.h>

static void
split_and_inv(const D *d, const T *ma[], T *lininv[], T *nonlin[])
{
  int nv = d->nv, cv = d->nmv, nk = nv - cv;                    // #vars, #canonical vars, #knobs
  gsl_matrix *mat_var  =      gsl_matrix_calloc(cv,cv),        // canonical vars
             *mat_vari =      gsl_matrix_alloc (cv,cv),        // inverse of vars
             *mat_knb  = nk ? gsl_matrix_calloc(cv,nk) : NULL, // knobs
             *mat_knbi = nk ? gsl_matrix_alloc (cv,nk) : NULL; // 'inverse' of knobs

  // split linear, (-1 * nonlinear)
  for (int i = 0; i < cv; ++i) {
    int v = 0;
    for (; v < cv; ++v) gsl_matrix_set(mat_var, i,v   , ma[i]->coef[v+1]);
    for (; v < nv; ++v) gsl_matrix_set(mat_knb, i,v-cv, ma[i]->coef[v+1]);

    FUN(copy)(ma[i], nonlin[i]);

    // clear constant and linear part
    for (int c = 0; c < d->ord2idx[2]; ++c)
      nonlin[i]->coef[c] = 0;
    nonlin[i]->nz = bclr(nonlin[i]->nz,0);
    nonlin[i]->nz = bclr(nonlin[i]->nz,1);
    FUN(scl)(nonlin[i],-1,nonlin[i]);
  }

  // invert linear
  gsl_permutation *p = gsl_permutation_alloc(cv);
  int signum;
  gsl_linalg_LU_decomp(mat_var, p, &signum);
  gsl_linalg_LU_invert(mat_var, p, mat_vari);
  gsl_permutation_free(p);

  if (nk != 0) {
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans,
                   -1, mat_vari, mat_knb,
                    0, mat_knbi);
  }

  // copy result into TPSA
  for (int i = 0; i < cv; ++i) {
    for (int v = 0; v < cv; ++v)
      FUN(seti)(lininv[i], v    +1, 0, gsl_matrix_get(mat_vari, i,v));
    for (int k = 0; k < nk; ++k)
      FUN(seti)(lininv[i], k+cv +1, 0, gsl_matrix_get(mat_knbi, i,k));
  }

  gsl_matrix_free(mat_var);
  gsl_matrix_free(mat_knb);
  gsl_matrix_free(mat_vari);
  gsl_matrix_free(mat_knbi);
}
*/

static void
split_and_inv(const D *d, const T *ma[], T *lininv[], T *nonlin[])
{
  ssz_t nv = d->nv, cv = d->nmv, nk = nv - cv; // #vars, #canonical vars, #knobs
  mad_alloc_tmp(NUM, mat_var , cv*cv); // canonical vars
  mad_alloc_tmp(NUM, mat_vari, cv*cv); // inverse of vars
  mad_alloc_tmp(NUM, mat_knb , cv*nk); // knobs
  mad_alloc_tmp(NUM, mat_knbi, cv*nk); // 'inverse' of knobs

  // split linear, (-1 * nonlinear)
  for (idx_t i = 0; i < cv; ++i) {
    idx_t v = 0;
    for (; v < cv; ++v) mat_var[i*cv +  v    ] = ma[i]->coef[v+1];
    for (; v < nv; ++v) mat_knb[i*nk + (v-cv)] = ma[i]->coef[v+1];

    FUN(copy)(ma[i], nonlin[i]);

    // clear constant and linear part
    for (idx_t c = 0; c < d->ord2idx[2]; ++c)
      nonlin[i]->coef[c] = 0;

    nonlin[i]->nz = mad_bit_clr(nonlin[i]->nz,0);
    nonlin[i]->nz = mad_bit_clr(nonlin[i]->nz,1);
    FUN(scl)(nonlin[i],-1,nonlin[i]);
  }

  // invert linear part: mat_vari = mat_var^-1
# ifndef MAD_CTPSA_IMPL
  mad_mat_invn(mat_var, 1, mat_vari, cv, cv, -1);
  if (nk != 0) {
    // mat_knbi = - mat_vari * mat_knb
    mad_mat_mul(mat_vari, mat_knb, mat_knbi, cv, nk, cv);
    mad_vec_muln(mat_knbi, -1, mat_knbi, cv*nk);
  }
# else
  mad_cmat_invn(mat_var, 1, mat_vari, cv, cv, -1);
  if (nk != 0) {
    // mat_knbi = - mat_vari * mat_knb
    mad_cmat_mul(mat_vari, mat_knb, mat_knbi, cv, nk, cv);
    mad_cvec_muln(mat_knbi, -1, mat_knbi, cv*nk);
  }
# endif

  // copy result into TPSA
  for (idx_t i = 0; i < cv; ++i) {
    for (idx_t v = 0; v < cv; ++v)
      FUN(seti)(lininv[i], v    +1, 0, mat_vari[i*cv + v]);
    for (idx_t k = 0; k < nk; ++k)
      FUN(seti)(lininv[i], k+cv +1, 0, mat_knbi[i*nk + k]);
  }

  mad_free_tmp(mat_var );
  mad_free_tmp(mat_vari);
  mad_free_tmp(mat_knb );
  mad_free_tmp(mat_knbi);
}

// --- public -----------------------------------------------------------------o

void
FUN(minv) (ssz_t sa, const T *ma[sa], T *mc[sa])
{
  assert(ma && mc);
  check_minv(sa,ma,mc);
  for (idx_t i = 0; i < sa; ++i)
    ensure(mad_bit_get(ma[i]->nz,1), "invalid domain");

  const D *d = ma[0]->d;
  T *lininv[sa], *nonlin[sa], *tmp[sa];
  for (idx_t i = 0; i < sa; ++i) {
    lininv[i] = FUN(newd)(d,1);
    nonlin[i] = FUN(new)(ma[i], mad_tpsa_same);
    tmp[i]    = FUN(new)(ma[i], mad_tpsa_same);
  }

  split_and_inv(d, ma, lininv, nonlin);

  // iteratively compute higher orders of the inverse
  // MC (OF ORDER I) = AL^-1 o [ I - ANL (NONLINEAR) o MC (OF ORDER I-1) ]

  for (idx_t i = 0; i < sa; ++i)
    FUN(copy)(lininv[i], mc[i]);

  ord_t o_prev = mad_desc_gtrunc(d, 2);
  for (ord_t o = 2; o <= d->mo; ++o) {
    mad_desc_gtrunc(d, o);
    FUN(compose)(sa, (const T**)nonlin, sa, (const T**)mc, tmp);

    for (idx_t v = 0; v < sa; ++v)
      FUN(seti)(tmp[v], v+1, 1,1);    // add I

    FUN(compose)(sa, (const T**)lininv, sa, (const T**)tmp, mc);
  }
  mad_desc_gtrunc(d, o_prev);

  // cleanup
  for (idx_t i = 0; i < sa; ++i) {
    FUN(del)(lininv[i]);
    FUN(del)(nonlin[i]);
    FUN(del)(tmp[i]);
  }
}

void
FUN(pminv) (ssz_t sa, const T *ma[sa], T *mc[sa], idx_t select[sa])
{
  assert(ma && mc && select);
  check_minv(sa,ma,mc);
  for (idx_t i = 0; i < sa; ++i)
    if (select[i])
      ensure(mad_bit_get(ma[i]->nz,1), "invalid domain");

  const D *d = ma[0]->d;
  // split input map into rows that are inverted and rows that are not
  T *mUsed[sa], *mUnused[sa], *mInv[sa];
  for (idx_t i = 0; i < sa; ++i) {
    if (select[i]) {
      mUsed  [i] = FUN(new) (ma[i], mad_tpsa_same);
      mInv   [i] = FUN(new) (ma[i], mad_tpsa_same);
      mUnused[i] = FUN(newd)(d,1);
      FUN(copy)(ma[i],mUsed[i]);
      FUN(seti)(mUnused[i], i+1,  0,1);
    }
    else {
      mUsed  [i] = FUN(newd)(d,1);
      mInv   [i] = FUN(newd)(d,1);
      mUnused[i] = FUN(new) (ma[i], mad_tpsa_same);
      FUN(copy)(ma[i],mUnused[i]);
      FUN(seti)(mUsed[i], i+1, 0,1);
    }
    FUN(set0)(mUsed  [i], 0,0);
    FUN(set0)(mUnused[i], 0,0);
  }

  FUN(minv)   (sa,(const T**)mUsed  ,              mInv);
  FUN(compose)(sa,(const T**)mUnused,sa,(const T**)mInv,mc);

  for (idx_t i = 0; i < sa; ++i) {
    FUN(del)(mUsed[i]);
    FUN(del)(mUnused[i]);
    FUN(del)(mInv[i]);
  }
}

// --- end --------------------------------------------------------------------o
