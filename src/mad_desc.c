/*
 o-----------------------------------------------------------------------------o
 |
 | Descriptor (TPSA) module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

// --- globals ----------------------------------------------------------------o

const  ord_t  mad_tpsa_default = -1;
const  ord_t  mad_tpsa_same    = -2;
const desc_t *mad_desc_curr    =  0;

// --- constants --------------------------------------------------------------o

enum { DESC_MAX_ORD = CHAR_BIT * sizeof(bit_t) };

// --- helpers ----------------------------------------------------------------o

static inline ssz_t
max_nc(ssz_t nv, ssz_t no)
{
  // #coeff(nv,no) = (nv+no)! / (nv! no!)
  ssz_t max = MAX(nv,no);
  u64_t num = 1, den = 1;
  for (ssz_t i=max+1; i <= nv+no; ++i) {
    num *= i;
    den *= i - max;
  }
  return num / den;
}

// --- monomials --------------------------------------------------------------o

static inline void
nxt_mono_by_unk(ssz_t n, const ord_t a[n], const idx_t sort[n],
                idx_t r, int o, ord_t m[n])
{
  assert(a && sort && m);
  mad_mono_fill(n, m, 0);
  for (idx_t k=r; k < n; ++k) {
    idx_t v = sort[k];
    m[v] = a[v];
    o -= a[v];
    if (o <  0) { m[v] += o; break; }
    if (o == 0) {            break; }
  }
}

// --- tables -----------------------------------------------------------------o

static inline void
tbl_print(ssz_t w, ssz_t h, ord_t **t)
{
  for (idx_t i=0; i < MIN(h,50); ++i) {
    printf("(%2d) ", i);
    mad_mono_print(w, t[i]);
    printf(" o=%d\n", mad_mono_ord(w, t[i]));
  }
  if (h > 50) printf("... [ %d more rows ] ...\n", h - 50);
}

static inline void
tbl_realloc_monos(D *d, ssz_t new_nc_)
{
  // reallocates array of monomials to fit at least new_nc_ items
  // if new_nc_ is 0, it increases the capacity of exiting array
  assert(d);
  if (!new_nc_) new_nc_ = d->nc * 2;
  d->nc = new_nc_;
  d->monos = mad_realloc(d->monos, d->nc * d->nv * sizeof *d->monos);
  d->ords  = mad_realloc(d->ords,  d->nc         * sizeof *d->monos);
}

static inline idx_t  // idx to end of d->monos
make_higher_ord_monos(D *d, idx_t curr_mono_idx, int need_realloc, idx_t var_at_idx[])
{
  // ords 2..mo
  int nv = d->nv;
  ord_t m[nv], *curr_mono = d->monos + curr_mono_idx*nv;
  idx_t *pi = d->ord2idx;

  for (ord_t o=2; o <= d->mo; o++) {          // to build ord o:
    for   (idx_t i=pi[ 1 ]; i < pi[2]; ++i) { // i goes through ord  1
      for (idx_t j=pi[o-1]; j < pi[o]; ++j) { // j goes through ord (o-1)
        mad_mono_add(nv, d->monos + i*nv, d->monos + j*nv, m);
        if (mad_desc_mono_isvalid_m(d, nv, m)) {
          // ensure there is space for m
          if (need_realloc && curr_mono_idx >= d->nc) {
            tbl_realloc_monos(d, 0);
            curr_mono = d->monos + nv*curr_mono_idx;
          }
          mad_mono_copy(nv, m, curr_mono);
          d->ords[curr_mono_idx] = o;
          curr_mono_idx++;
          curr_mono += nv;
        }
        // TODO: find better way to stop
        if (m[var_at_idx[i]] > d->var_ords[var_at_idx[i]] || m[var_at_idx[i]] >= o)
          break;
      }
      pi[o+1] = curr_mono_idx;
    }
  }

#if DEBUG > 1
  printf("ps={ ");
  for (int o=0; o <= d->mo+1; ++o)
    printf("[%d]=%d ", o, pi[o]);
  printf("}\n");
#endif

  return curr_mono_idx;
}

static inline void
make_monos(D *d)
{
  // builds the monomials matrix in To order
  assert(d && d->var_ords);
  const int max_init_alloc = 20000;  // to fit (6,12)
  d->nc = max_nc(d->nv, d->mo);
  int nv = d->nv, need_realloc = 0;
  if (d->nc > max_init_alloc || d->nc < 0) {  // nc < 0 when overflow in max_nc
    need_realloc = 1;
    d->nc = max_init_alloc;
  }

  d->monos   = mad_malloc(d->nc*nv  * sizeof *d->monos);
  d->ords    = mad_malloc(d->nc     * sizeof *d->ords);
  d->ord2idx = mad_malloc((d->mo+2) * sizeof *d->ord2idx);

  // ord 0
  mad_mono_fill(nv, d->monos, 0);
  d->ord2idx[0] = d->ords[0] = 0;
  d->ord2idx[1] = 1;
  idx_t curr_mono_idx = 1;

  // ord 1
  idx_t var_at_idx[nv+1];
  if (d->mo >= 1) {
    mad_mono_fill(nv*nv, d->monos+nv, 0);  // TODO: check d->monos has nv*nv slots
    for (idx_t i=0; i < nv; ++i)
      if (d->var_ords[i]) {
        d->monos  [curr_mono_idx*nv + i] = 1;
        d->ords   [curr_mono_idx       ] = 1;
        var_at_idx[curr_mono_idx       ] = i;
        curr_mono_idx++;
      }
    d->ord2idx[2] = curr_mono_idx;
  }

  int real_nc = curr_mono_idx;
  if (d->mo >= 2)
    real_nc = make_higher_ord_monos(d, curr_mono_idx, need_realloc, var_at_idx);

  tbl_realloc_monos(d, real_nc);
  d->size += real_nc*d->nv * sizeof *d->monos;
  d->size += real_nc       * sizeof *d->ords;
  d->size += (d->mo+2)     * sizeof *d->ord2idx;
}

/* kept for debugging binary search below.
static inline int
find_index_lin(ssz_t n, ord_t **T_, const ord_t m[n], idx_t start, idx_t stop)
{
  assert(T_ && m);
  const ord_t **T = (const ord_t **)T_;
  for (idx_t i=start; i < stop; ++i)
    if (mad_mono_eq(n, T[i], m)) return i;

  // error
  mad_mono_print(n, m);
  error("this monomial was not found in table (unexpected)");
  return -1; // never reached
}
*/

enum tbl_ordering {BY_ORD, BY_VAR};

static inline int
find_index_bin(ssz_t n, ord_t **T_, const ord_t m[n], idx_t from, idx_t to,
               enum tbl_ordering tbl_ord)
{
  const ord_t **T = (const ord_t **)T_;
  idx_t start = from, count = to-from, i = 0, step = 0;
  int ord_m = mad_mono_ord(n, m);
  while (count > 0) {
    step = count / 2;
    i = start + step;
    if ((tbl_ord == BY_ORD && mad_mono_ord (n, T[i]) < ord_m)
                           || mad_mono_rcmp(n, T[i], m) < 0) {
      start = ++i;
      count -= step + 1;
    }
    else
      count = step;
  }
  if (start < to && mad_mono_eq(n, T[start], m))
    return start;

  // error
  tbl_print(n, to, T_);
  mad_mono_print(n, m);
  error("this monomial was not found in table_by_%s from %d to %d (unexpected)",
         tbl_ord == BY_ORD ? "ord" : "var", from, to);
  return -1; // never reached
}

static inline void
tbl_by_ord(D *d)
{
  assert(d && d->var_ords && d->ord2idx && d->monos);

  d->To    = mad_malloc(d->nc * sizeof *(d->To));
  d->size += d->nc * sizeof *(d->To);

  ord_t *monos = d->monos;
  for (int i=0; i < d->nc; ++i, monos += d->nv)
    d->To[i] = monos;

#if DEBUG > 1
  tbl_print(d->nv, d->nc, d->To);
#endif
}

/**
 Tv depends on To, as it needs to search its monomials there.
 But H matrix wants Tv to be built according to sorted d->var_ords.
 So we pretend we build it with sorted d->var_ords, whereas actually we just
 generate nxt_mono_by_var differently, such that the variable with the highest
 order varies the fastest. And H will access it based on the sorted order
 */
static inline void
tbl_by_var(D *d)
{
  assert(d && d->var_ords && d->sort_var && d->monos && d->ord2idx);

  d->Tv    = mad_malloc(d->nc * sizeof *d->Tv);
  d->tv2to = mad_malloc(d->nc * sizeof *d->tv2to);
  d->to2tv = mad_malloc(d->nc * sizeof *d->to2tv);
  d->size +=   d->nc * sizeof *d->Tv;
  d->size += 2*d->nc * sizeof *d->tv2to;  // tv2to + to2tv

  int mi = 0, nv = d->nv;
  ord_t m[nv];
  mad_mono_fill(nv, m, 0);
  do {
    int o = mad_mono_ord(nv, m);
    int idx = find_index_bin(nv, d->To, m, d->ord2idx[o], d->ord2idx[o+1], BY_ORD);
    d->tv2to[mi]  = idx;
    d->to2tv[idx] = mi;
    d->Tv[mi]     = d->To[idx];
    ++mi;
  } while (mad_desc_mono_nxtbyvar(d, nv, m));
  assert(mi == d->nc);

#if DEBUG > 1
  printf("sort:  ");
  for (int i=0; i < d->nv; ++i)
    printf("%d ", d->sort_var[i]);
  printf("\n");
  tbl_print(nv, d->nc, d->Tv);
  printf("tv2to=[ ");
  for (int i=0; i < MIN(d->nc,50); ++i) printf("%d ", d->tv2to[i]);
  printf("%s]\n", d->nc > 50 ? " ... " : "");
  printf("to2tv=[ ");
  for (int i=0; i < MIN(d->nc,50); ++i) printf("%d ", d->to2tv[i]);
  printf("%s]\n", d->nc > 50 ? " ... " : "");
#endif
}

// --- H indexing matrix ------------------------------------------------------

static inline void
tbl_print_H(const D *d)
{
  assert(d && d->H);
  int cols = d->mo+2;
  for (int i=0; i < d->nv; ++i) {
    printf("%2d | ", d->sort_var[i]);
    for (int j=0; j < cols; ++j)
      printf("%2d ", d->H[i*cols + j]);
    printf("\n");
  }
}

static inline int
tbl_index_H(const D *d, int n, const ord_t m[n])
{
  assert(d && n <= d->nv);
  int s = 0, I = 0, cols = d->mo+2;
  const idx_t *H = d->H, *sort = d->sort_var;
  idx_t v;

  for (int r = n-1; r >= 0; --r) {
    v = sort[r];
    I += H[r*cols + s+m[v]] - H[r*cols + s];
    s += m[v];
  }
  assert(I > -1);
  return I;
}

static inline int
tbl_index_H_sm(const D *d, int n, const idx_t m[n])
{
  assert(d && n/2 <= d->nv && !(n & 1));
  int s = 0, I = 0, cols = d->mo+2, idx, o;
  const int *H = d->H;
  // indexes in m should be in ascending order
  for (int i = n-2; i >= 0; i -= 2) {
    idx = m[i] - 1;
    o = m[i+1];
    I += H[idx*cols + s+o] - H[idx*cols + s];
    s += o;
  }
  assert(I > -1);
  return I;
}

static inline void
tbl_clear_H(D *d)
{
  int cols = d->mo+2, accum = 0;
  idx_t *sort = d->sort_var, si;

  for (int r=d->nv-1; r >= 0; --r) {  // rows (sorted vars)
    si = sort[r];
    accum += d->var_ords[si];
    for (int o=1 + MIN(accum, d->mo); o < cols; ++o)  // ords
      d->H[r*cols + o] = -1;
  }
}

static inline void
tbl_solve_H(D *d)
{
  idx_t *sort = d->sort_var, v;
  int nv = d->nv, cols = d->mo+2, accum = d->var_ords[sort[nv-1]];
  ord_t mono[nv];
  const ord_t *vo = d->var_ords;

  // solve system of equations
  // r goes through H rows (vars) from lowest to highest ord
  for (int r=nv-2; r >= 1; --r) {
    v = sort[r];  // the var at this row
    accum += vo[v];
    for (int o=vo[v]+2; o <= MIN(accum, d->mo); o++) { // orders
      nxt_mono_by_unk(nv, vo, sort, r, o, mono);  // NOTE: it's r, not v
      if (mad_desc_mono_isvalid_m(d, nv, mono)) {
        idx_t idx0 = tbl_index_H(d, nv, mono);
        idx_t idx1 = find_index_bin(nv, d->Tv, mono, idx0, d->nc, BY_VAR);
        d->H[r*cols + o] = idx1 - idx0;
      }
      else
        d->H[r*cols + o] = 0;
    }
  }
}

static inline void
tbl_build_H(D *d)
{
  idx_t *H = d->H;
  int rows = d->nv, cols = d->mo + 2, nc = d->nc;
  ord_t **Tv = d->Tv;
  idx_t *sort = d->sort_var;
  const ord_t *vo = d->var_ords;

  // minimal constants for 1st row
  for (int c = 0; c < cols; ++c)
    H[0*cols + c] = c;

  // remaining rows == vars in descending order
  for (int r = 1; r < rows; ++r) {
    H[r*cols + 0] = 0;
    int curr_col = 1;
    idx_t var_idx = sort[r];

    // initial congruence from Tv
    for (int m = 1; m < nc; m++) { // monomials
      if (Tv[m][var_idx] != Tv[m-1][var_idx]) {
        H[r*cols + curr_col] = m;
        curr_col++;
        if (Tv[m][var_idx] == 0) break;
      }
    }

    // complete row with zeros
    while(curr_col < cols) H[r*cols + curr_col] = 0, curr_col++;
  }

  // close congruence of the last var
  idx_t var_idx = sort[rows-1];
  H[(rows-1)*cols + vo[var_idx] + 1] = nc;

#if DEBUG > 1
  tbl_print_H(d);
#endif
}

static inline void
tbl_set_H(D *d)
{
  assert(d && d->var_ords && d->Tv && d->To);
  assert(d->nv != 0);

  d->H = mad_malloc(d->nv * (d->mo+2) * sizeof *(d->H));
  d->size += d->nv * (d->mo+2) * sizeof *(d->H);

  tbl_build_H(d);
  tbl_solve_H(d);
  tbl_clear_H(d);

#if DEBUG > 1
  printf("H = {\n");
  tbl_print_H(d);
#endif
}

// --- L indexing matrix ------------------------------------------------------

static inline void
tbl_print_LC(const idx_t *lc, int oa, int ob, int *pi)
{
  int iao = pi[oa], ibo = pi[ob], cols = pi[oa+1] - pi[oa];
  for (int ib = pi[ob]; ib < pi[ob+1]; ++ib) {
    printf("\n  ");
    for (int ia = pi[oa]; ia < pi[oa+1]; ++ia) {
      int ic = lc[hpoly_idx(ib-ibo,ia-iao,cols)];
      printf("%3d ", ic);
    }
  }
  printf("\n");
}

static inline void
tbl_print_L(const D *d)
{
  int ho = d->mo / 2;
  for (int oc = 2; oc <= MIN(d->mo,5); ++oc)
    for (int j = 1; j <= oc/2; ++j) {
      int oa = oc - j, ob = j;
      printf("L[%d][%d] = {", ob, oa);
      tbl_print_LC(d->L[oa*ho + ob], oa, ob, d->ord2idx);
    }
  if (d->mo > 5)
    printf("Orders 5 to %d omitted...\n", d->mo);
}

static inline idx_t*
tbl_build_LC(int oa, int ob, D *d)
{
#if DEBUG > 1
  printf("tbl_set_LC oa=%d ob=%d\n", oa, ob);
#endif
  assert(d && d->To && d->ord2idx && d->tv2to);
  assert(oa < d->mo && ob < d->mo);

  ord_t **To = d->To;
  const int *pi   = d->ord2idx,      *tv2to = d->tv2to,          // shorter names
             iao  = pi[oa],            ibo  = pi[ob],            // offsets
             cols = pi[oa+1] - pi[oa], rows = pi[ob+1] - pi[ob]; // sizes

  int mat_size = rows * cols;
  idx_t *lc = mad_malloc(mat_size * sizeof *lc);
  d->size += mat_size * sizeof *lc;
  for (int i = 0; i < mat_size; ++i)
    lc[i] = -1;

  ord_t m[d->nv];
  int ic, idx_lc;
  for (int ib = pi[ob]; ib < pi[ob+1]; ++ib) {
    int lim_a = oa == ob ? ib+1 : pi[oa+1];   // triangular is lower left
    for (int ia = pi[oa]; ia < lim_a; ++ia) {
      mad_mono_add(d->nv, To[ia], To[ib], m);
      if (mad_desc_mono_isvalid_m(d,d->nv,m)) {
        ic = tv2to[tbl_index_H(d,d->nv,m)];
        idx_lc = hpoly_idx(ib-ibo, ia-iao, cols);
        lc[idx_lc] = ic;
#if DEBUG > 2
        printf(" ib=%d ", ib); mad_mono_print(d->nv, To[ib]);
        printf(" ia=%d ", ia); mad_mono_print(d->nv, To[ia]);
        printf(" ic=%d ", ic); mad_mono_print(d->nv, m);
        printf(" ilc=%d\n", idx_lc);
#endif
      }
    }
  }

  return lc;
}

static inline int**
get_LC_idxs(int oa, int ob, D *d)
{
  int oc = oa + ob;
  const idx_t *pi = d->ord2idx, *lc = d->L[d->mo/2*oa + ob],
                T = (pi[oc+1] + pi[oc] - 1) / 2;  // splitting threshold of oc
  const int  cols =  pi[oa+1] - pi[oa],
             rows =  pi[ob+1] - pi[ob];

  int **LC_idx = mad_malloc(3 * sizeof *LC_idx);
  d->size += 3 * sizeof *LC_idx;

  int *limits = mad_malloc(3 * rows * sizeof *limits); // rows: [start split end]
  d->size += 3 * rows * sizeof *limits;

  const int START = 0, SPLIT = 1, END = 2;

  LC_idx[START] = limits;
  LC_idx[SPLIT] = limits +   rows;
  LC_idx[END  ] = limits + 2*rows;

  for (int ib = 0; ib < rows; ++ib) {
    int ia;
    for (ia = 0; lc[hpoly_idx(ib,ia,cols)] == -1; ++ia)
      ;  // shift ia to first valid entry
    LC_idx[START][ib] = ia;

    for (ia = oa == ob ? ib : cols-1; lc[hpoly_idx(ib,ia,cols)] == -1; --ia)
      ;  // shift ia to last valid entry
    LC_idx[END  ][ib] = ia + 1;

    LC_idx[SPLIT][ib] = oa == ob ? ib+1 : cols;
    for (ia = LC_idx[START][ib]; ia < LC_idx[END][ib]; ++ia)
      if (lc[hpoly_idx(ib,ia,cols)] >= T) {
        LC_idx[SPLIT][ib] = ia;
        break;
      }
  }

#if DEBUG > 1
  if (oc <= 5) {
    printf("LC_idx[%d][%d] = { [T=%d]\n", ob, oa, T);
    printf("  -->\t  //\t<--\n");
    for (int r = 0; r < rows; ++r)
      printf("  [%3d\t%3d\t%3d]\n", LC_idx[START][r],LC_idx[SPLIT][r],LC_idx[END][r]);
  }
#endif

  return LC_idx;
}

static inline void
tbl_set_L(D *d)
{
  ord_t o = d->mo, ho = d->mo / 2;
  int size_L = (o*ho + 1) * sizeof *(d->L);
  d->L = mad_malloc(size_L);
  d->size += size_L;

  int size_lci = (o*ho + 1) * sizeof *(d->L_idx);
  d->L_idx = mad_malloc(size_lci);
  d->size += size_lci;

  memset(d->L,     0, size_L);
  memset(d->L_idx, 0, size_lci);
  // #ifdef _OPENMP
  // #pragma omp parallel for schedule(guided,1)
  // #endif
  for (int oc = 2; oc <= d->mo; ++oc)
    for (int j = 1; j <= oc / 2; ++j) {
      int oa = oc - j, ob = j;

      d->L    [oa*ho + ob] = tbl_build_LC(oa, ob, d);
      d->L_idx[oa*ho + ob] = get_LC_idxs (oa, ob, d);
    }

#if DEBUG > 1
  tbl_print_L(d);
#endif
}

// --- descriptor internal checks ---------------------------------------------o

static int
tbl_check_L(D *d)
{
  assert(d && d->ord2idx && d->L && d->var_ords && d->To && d->H);
  int ho = d->mo / 2, *pi = d->ord2idx;
  ord_t m[d->nv];
  for (int oc = 2; oc <= d->mo; ++oc)
    for (int j = 1; j <= oc / 2; ++j) {
      int oa = oc - j, ob = j;
      idx_t *lc = d->L[oa*ho + ob];
      if (!lc)                                     return  1e7 + oa*1e3 + ob;

      int sa = pi[oa+1]-pi[oa], sb = pi[ob+1]-pi[ob];

      for (int ibl = 0; ibl < sb; ++ibl) {
        int lim_a = oa == ob ? ibl+1 : sa;
        for (int ial = 0; ial < lim_a; ++ial) {
          int ib = ibl + pi[ob], ia = ial + pi[oa];
          int il = hpoly_idx(ibl,ial,sa);
          if (il < 0)                              return -2e7 - ia*1e5 - ib;
          if (il >= sa * sb)                       return  2e7 + ia*1e5 + ib;

          int ic = lc[il];
          if (ic >= pi[oc+1])                      return  3e7 + ic*1e5 + 11;
          if (ic >= 0 && ic < d->ord2idx[oc])      return  3e7 + ic*1e5 + 12;

          mad_mono_add(d->nv, d->To[ia], d->To[ib], m);
          if (ic < 0 && mad_desc_mono_isvalid_m(d,d->nv,m))
                                                   return -3e7          - 13;
        }
      }
    }
  return 0;
}

static int  // error code
tbl_check(D *d)
{
  const idx_t
      *tv2to = d->tv2to,
      *to2tv = d->to2tv,
          *H = d->H,
       *sort = d->sort_var,
          nv = d->nv,
        cols = d->mo+2;
  ord_t **Tv = d->Tv,
        **To = d->To,
      *monos = d->monos;
  const ord_t *vo = d->var_ords;

  // check H
  for (int i = 0; i < nv; ++i)
    if (H[i*cols + 0] != 0)                   return 1e6 + i;

  for (int r = 1; r < nv; ++r) {
    if (vo[sort[r]] == 0) {
      for (int o = 1; o <= d->mo + 1; ++o)
        if (H[r*cols + o] != -1)              return 2e6 + r;
    } else // TODO: restore test (LD:?)
    if (H[r*cols + 1] != (H[(r-1)*cols + vo[sort[r-1]]  + 1] != -1 ?
                          H[(r-1)*cols + vo[sort[r-1]]  + 1] :
                          H[(r-1)*cols + vo[sort[r-1]]] + 1) )
                                              return 3e6 + r;
  }

  for (int i=0; i < d->nc; ++i) {
    if (to2tv[tv2to[i]] != i)                 return 4e6 + i;
    if (tv2to[tbl_index_H(d,nv,To[i])] != i)  return 5e6 + i;
    if (!mad_mono_eq(nv,To[tv2to[i]],Tv[i]))  return 6e6 + i;
    if (!mad_mono_eq(nv,To[i],monos + nv*i))  return 7e6 + i;
  }

  return tbl_check_L(d);
}

// --- thread dispatch --------------------------------------------------------o

static inline void
get_ops(D *d, long long int ops[])
{
  int *pi = d->ord2idx;
  ops[0] = ops[1] = ops[2] = 0;
  for (int o = 3; o <= d->mo; ++o) {
    ops[o] = 0;
    for (int j = 1; j <= (o-1)/2; ++j) {
      int oa = o-j, ob = j;            // oa > ob >= 1
      int na = pi[oa+1] - pi[oa], nb = pi[ob+1] - pi[ob];
      ops[o] += 2 * na * nb;
    }
    if (!(o & 1)) {
      int ho = o/2;
      ops[o] += (pi[ho+1]-pi[ho]) * (pi[ho+1]-pi[ho]);
    }
  }
  ops[d->mo+1] = ops[d->mo]/2;
  ops[d->mo]  -= ops[d->mo+1];
}

static inline int
get_min_dispatched_idx(int nb_threads, long long int dops[])
{
  long long int min_disp = dops[nb_threads-1];
  int min_disp_idx = nb_threads - 1;
  for (int t = nb_threads-1; t >= 0; --t)
    if (dops[t] <= min_disp) {
      min_disp = dops[t];
      min_disp_idx = t;
    }
  return min_disp_idx;
}

static inline void
build_dispatch (D *d)
{
  // [0] serial(all), [1..nth] parallel(split)
  int nth = d->nth + (d->nth > 1);

  d->ocs = mad_malloc(nth * sizeof *(d->ocs));
  d->size += nth * sizeof *(d->ocs);

  int sizes[nth];
  for (int t = 0; t < nth; ++t) {
    d->ocs[t] = mad_calloc(d->mo, sizeof *d->ocs[0]);
    d->size += (d->mo) * sizeof *d->ocs[0];
    sizes[t] = 0;
  }

  long long int ops[d->mo+2], dops[nth];
  memset(dops, 0, nth * sizeof *dops);
  get_ops(d, ops);

  // serial
  for (int o = d->mo; o > 2; --o) {
    d->ocs[0][sizes[0]++] = o;
    dops[0] += ops[o];
  }
  dops[0] += ops[d->mo+1];

  // parallel
  if (nth > 1) {
    for (int o = d->mo+1; o > 2; --o) {
      int idx = get_min_dispatched_idx(d->nth,dops+1) + 1;
      assert(idx > 0 && idx <= d->nth);
      d->ocs[idx][sizes[idx]++] = o;
      dops[idx] += ops[o];
    }
  }

#ifdef DEBUG
  printf("\nTHREAD DISPATCH:\n");
  for (int t = 0; t < nth; ++t) {
    printf("[%d]: ", t);
    for (int i = 0; d->ocs[t][i]; ++i)
      printf("%d ", d->ocs[t][i]);
    printf("[ops:%lld] \n", dops[t]);
  }
  printf("\n");
#endif
}

static inline void
set_temps (D *d)
{
  d->  t = mad_malloc(DESC_MAX_TMP * d->nth * sizeof *d-> t );
  d-> ct = mad_malloc(DESC_MAX_TMP * d->nth * sizeof *d->ct );
  d-> ti = mad_malloc(               d->nth * sizeof *d-> ti);
  d->cti = mad_malloc(               d->nth * sizeof *d->cti);

  for(int j = 0; j < d->nth; ++j) {
  for(int i = 0; i < DESC_MAX_TMP; ++i) {
    d-> t[j*DESC_MAX_TMP+i] = mad_tpsa_newd (d,d->mo);
    d->ct[j*DESC_MAX_TMP+i] = mad_ctpsa_newd(d,d->mo); }
    d->ti[j] = d->cti[j] = 0;
  }

#if DEBUG > 1
  printf("\nTEMPS #TPSA = 2 (R&C) x %d (#TMPS) x %d (Threads) = %d\n"
         "TEMPS TMEM  = %d (TPSA) x %d (nc) = %llu bytes\n",
          DESC_MAX_TMP, d->nth, 2*DESC_MAX_TMP*d->nth, 2*DESC_MAX_TMP*d->nth,
          d->nc, 2*DESC_MAX_TMP*d->nth* 3ull*d->nc*sizeof(num_t)/2);
#endif
}

static inline void
del_temps (D *d)
{
  for(int j = 0; j < d->nth; ++j) {
  for(int i = 0; i < DESC_MAX_TMP; ++i) {
    mad_tpsa_del (d-> t[j*DESC_MAX_TMP+i]);
    mad_ctpsa_del(d->ct[j*DESC_MAX_TMP+i]);
  }}
  mad_free(d->  t);
  mad_free(d-> ct);
  mad_free(d-> ti);
  mad_free(d->cti);
}

// --- descriptor management --------------------------------------------------o

enum { MAX_TPSA_DESC = 50 };    // max number of simultaneous descriptors

static D  *Ds[MAX_TPSA_DESC];

static inline void
set_var_ords (D *d, const ord_t var_ords[])
{
  assert(d && var_ords);

  d->sort_var = mad_malloc(d->nv * sizeof *d->sort_var);
  mad_mono_sort(d->nv, var_ords, d->sort_var);
  d->size += d->nv * sizeof *d->sort_var;

  ord_t *vo = mad_malloc(d->nv * sizeof *d->var_ords);
  mad_mono_copy(d->nv, var_ords, vo);
  d->var_ords = vo;
  d->size += d->nv * sizeof *d->var_ords;

#if DEBUG > 1
  printf("var_ords sorting: [");
  for(int i = 0; i < d->nv; ++i)
    printf("%d ", d->sort_var[i]);
  printf("]\n");
#endif
}

static inline D*
desc_init (int nmv, const ord_t mvar_ords[nmv], int nv, ord_t ko)
{
  assert(mvar_ords);

  D *d = mad_malloc(sizeof *d);
  memset(d, 0, sizeof *d);
  d->size = sizeof *d;
  d->nmv = nmv;
  d->nv = nv;
  d->ko = ko;
  d->mo = d->to = mad_mono_max(nmv, mvar_ords);

  ord_t *mo = mad_malloc(nmv * sizeof *d->mvar_ords);
  mad_mono_copy(nmv, mvar_ords, mo);
  d->mvar_ords = mo;
  d->size += nmv * sizeof *d->mvar_ords;
  d->nth = omp_get_max_threads();

  return d;
}

static D*
desc_build (int nmv, const ord_t mvar_ords[nmv],
            int nv , const ord_t  var_ords[nv ], ord_t ko)
{
  assert(mvar_ords && var_ords);

  // input validation
  ensure(nmv > 0 && nmv < 1000, "invalid number of map variables [%d]", nmv);
  ensure(nv >= nmv, "#map variables [%d] exceeds #variables [%d]", nmv, nv);

  // variables max orders validation
  ord_t mo = mad_mono_max(nmv, mvar_ords);
  ensure(mo <= DESC_MAX_ORD,
         "some map variables order exceeds maximum order [%u]", DESC_MAX_ORD);
  ensure(mad_mono_le(nmv, var_ords, mvar_ords),
         "some variables orders exceed map variables orders");
  ensure(mad_mono_max(nv, var_ords) <= mo,
         "some variables orders exceed map variables orders");
  ensure(mad_mono_min(nv, var_ords) > 0,
         "some variables have invalid zero order");

  D *d = desc_init(nmv, mvar_ords, nv, ko);

  set_var_ords(d, var_ords);
  make_monos(d);
  tbl_by_ord(d);
  tbl_by_var(d);  // requires To
  tbl_set_H(d);
  tbl_set_L(d);
  build_dispatch(d);
  set_temps(d);

#if DEBUG > 1
  printf("nc = %d ---- Total desc size: %ld bytes\n", d->nc, d->size);
#endif

  int err = tbl_check(d);
  if (err != 0) {
    printf("\nA= ");   mad_mono_print(d->nv, d->var_ords);
    printf("\nH=\n");  tbl_print_H(d);
    printf("\nTv=\n"); tbl_print(d->nv, d->nc, d->Tv);
    printf("\nTo=\n"); tbl_print(d->nv, d->nc, d->To);
    printf("\nChecking table consistency ... %d\n", err);
    assert(NULL);
  }

  return d;
}

static inline int
desc_equiv (const D *d, int nmv, const ord_t mvar_ords[nmv],
                        int nv , const ord_t  var_ords[nv ], ord_t ko)
{
  return d->nmv == nmv && d->nv == nv && d->ko == ko
         && mad_mono_eq(nmv, d->mvar_ords, mvar_ords)
         && mad_mono_eq(nv , d-> var_ords,  var_ords);
}

static inline D*
get_desc (int nmv, const ord_t mvar_ords[nmv],
          int nv , const ord_t  var_ords[nv ], ord_t ko)
{
  assert(mvar_ords && var_ords);

  for (int i=0; i < MAX_TPSA_DESC; ++i)
    if (Ds[i] && desc_equiv(Ds[i], nmv, mvar_ords, nv, var_ords, ko))
      return mad_desc_curr=Ds[i], Ds[i];

  for (int i=0; i < MAX_TPSA_DESC; ++i)
    if (!Ds[i]) {
      Ds[i] = desc_build(nmv, mvar_ords, nv, var_ords, ko);
      Ds[i]->id = i;
      return mad_desc_curr=Ds[i], Ds[i];
    }

  error("Too many descriptors in concurrent use");
}

// --- public -----------------------------------------------------------------o

int
mad_desc_mono_isvalid_m (const D *d, ssz_t n, const ord_t m[n])
{
  assert(d && m);
  return n <= d->nv
         && mad_mono_ord(n, m) <= d->mo
         && mad_mono_le (n, m, d->var_ords)
         && mad_mono_ord(n-d->nmv, m+d->nmv) <= d->ko;
}

int
mad_desc_mono_isvalid_s (const D *d, ssz_t n, str_t s)
{
  assert(s && n <= 1000000);
  if (n <= 0) n = strlen(s);
  ord_t m[n];
  n = mad_mono_str(n, m, s);
  return mad_desc_mono_isvalid_m(d, n, m);
}

int
mad_desc_mono_isvalid_sm (const D *d, ssz_t n, const idx_t m[n])
{
  assert(d && m);
  if (n & 1) return 0;

  int mo = 0, ko = 0;
  for (idx_t i = 0; i < n; i += 2) {
    idx_t mono_idx = m[i] - 1; // translate from var idx to mono idx
    if (mono_idx >= d->nv)         return 0;

    ord_t o = m[i+1];
    if (o > d->var_ords[mono_idx]) return 0;
    mo += o;
    if (mono_idx > d->nmv) ko += o;
  }
  return mo <= d->mo && ko <= d->ko;
}

int
mad_desc_mono_nxtbyvar (const D *d, ssz_t n, ord_t m[n])
{
  assert(d && m);
  const idx_t *sort = d->sort_var;
  for (idx_t i = 0; i < n; ++i) {
    ++m[sort[i]];
    if (mad_desc_mono_isvalid_m(d, n, m)) return 1;
    m[sort[i]] = 0;
  }
  return 0;
}

ord_t
mad_desc_get_mono (const D *d, ssz_t n, ord_t m_[n], idx_t i)
{
  assert(d);
  ensure(0 <= n && n < d->nv, "invalid monomial length");
  ensure(0 <= i && i < d->nc, "index out of bounds");
  if (m_ && n) mad_mono_copy(n, d->To[i], m_);
  return d->ords[i];
}

idx_t
mad_desc_get_idx_m (const D *d, ssz_t n, const ord_t m[n])
{
  assert(d && m);
  ensure(mad_desc_mono_isvalid_m(d, n, m), "invalid monomial");
  return d->tv2to[tbl_index_H(d, n, m)];
}

idx_t
mad_desc_get_idx_s (const D *d, ssz_t n, str_t s)
{
  assert(s);
  if (n <= 0) n = strlen(s);
  ord_t m[n];
  n = mad_mono_str(n, m, s);
  return mad_desc_get_idx_m(d, n, m);
}

idx_t
mad_desc_get_idx_sm (const D *d, ssz_t n, const idx_t m[n])
{
  assert(d && m);
  ensure(mad_desc_mono_isvalid_sm(d, n, m), "invalid monomial");
  return d->tv2to[tbl_index_H_sm(d, n, m)];
}

ord_t
mad_desc_maxord (const D *d)
{
  assert(d);
  return d->mo;
}

ssz_t
mad_desc_maxlen (const D *d)
{
  assert(d);
  return d->nc;
}

ssz_t
mad_desc_ordlen (const D *d, ord_t mo)
{
  assert(d);
  ensure(mo <= d->mo, "invalid order (exceeds maximum order)");
  return d->ord2idx[mo+1];
}

ord_t
mad_desc_gtrunc (const D *d_, ord_t to)
{
  D* d = (void*)d_;
  assert(d);
  if (to == mad_tpsa_same)
    return d->to;

  ord_t orig = d->to;

  if (to == mad_tpsa_default)
    return d->to = d->mo, orig;

  ensure(to <= d->mo, "invalid order (exceeds maximum order)");
  return d->to = to, orig;
}

// --- ctors, dtor ------------------------------------------------------------o

const D*
mad_desc_newn (int nmv, ord_t mvo)
{
  ensure(nmv > 0     , "invalid map variables orders specification");
  ensure(nmv < 100000, "too many variables");

  ord_t mvar_ords[nmv];
  mad_mono_fill(nmv, mvar_ords, mvo);

  return get_desc(nmv, mvar_ords, nmv, mvar_ords, 0);
}

const D*
mad_desc_newm (int nmv, const ord_t mvar_ords[nmv])
{
  assert(mvar_ords);
  ensure(nmv > 0     , "invalid map variables orders specification");
  ensure(nmv < 100000, "too many variables");

  return get_desc(nmv, mvar_ords, nmv, mvar_ords, 0);
}

const D*
mad_desc_newk (int nmv, ord_t mvo, int nk, ord_t ko, ord_t dk)
{
  int nv = nmv+nk;
  ensure(nmv > 0     , "invalid map variables orders specification");
  ensure(nk  >= 0    , "invalid knob variables orders specification");
  ensure(nv  < 100000, "too many variables");

  ord_t var_ords[nv];
  mad_mono_fill(nmv, var_ords    , mvo);
  mad_mono_fill(nk , var_ords+nmv,  ko);

  if (nk > 0) {
    if (!dk || dk > ko) dk = ko;
  } else dk = 0;

  return get_desc(nmv, var_ords, nv, var_ords, dk);
}

const D*
mad_desc_newv (int nmv, const ord_t mvar_ords[nmv],
               int nv_, const ord_t _var_ords[nv_], ord_t dk)
{
  int nv = MAX(nmv, nv_);
  assert(mvar_ords);
  ensure(nmv > 0     , "invalid map variables orders specification");
  ensure(nv  < 100000, "too many variables");

  ord_t var_ords[nv];
  mad_mono_copy(nmv, mvar_ords, var_ords);

  if (_var_ords) {
    ensure(nv_ >= 0, "invalid variable orders specification");
    mad_mono_copy(nv_, _var_ords, var_ords); // override overlapping part
  }

  int nk = nv-nmv;
  if (nk > 0) {
    ord_t ko = mad_mono_max(nk, var_ords+nmv);
    if (!dk || dk > ko) dk = ko;
  } else dk = 0;

  return get_desc(nmv, mvar_ords, nv, var_ords, dk);
}

const D*
mad_desc_newkv (int nmv, const ord_t mvar_ords[nmv],
                int nk , const ord_t kvar_ords[nk ],
                int nv_, const ord_t _var_ords[nv_], ord_t dk)
{
  int nv = MAX(nmv+nk, nv_);
  assert(mvar_ords && kvar_ords);
  ensure(nmv > 0     , "invalid map variables orders specification");
  ensure(nk  >= 0    , "invalid knob variables orders specification");
  ensure(nv  < 100000, "too many variables");

  ord_t var_ords[nv];
  mad_mono_copy(nmv, mvar_ords, var_ords    );
  mad_mono_copy(nk , kvar_ords, var_ords+nmv);

  if (_var_ords) {
    ensure(nv_ >= 0, "invalid variable orders specification");
    mad_mono_copy(nv_, _var_ords, var_ords); // override overlapping part
  }

  if (nk > 0) {
    ord_t ko = mad_mono_max(nk, var_ords+nmv);
    if (!dk || dk > ko) dk = ko;
  } else dk = 0;

  return get_desc(nmv, mvar_ords, nv, var_ords, dk);
}

void
mad_desc_del (const D *d_)
{
  D *d = (void*)d_;
  assert(d);
  mad_free((void*)d->mvar_ords);
  mad_free((void*)d->var_ords);
  mad_free(d->monos);
  mad_free(d->ords);
  mad_free(d->To);
  mad_free(d->Tv);
  mad_free(d->ord2idx);
  mad_free(d->tv2to);
  mad_free(d->to2tv);
  mad_free(d->H);

  if (d->L) {  // if L exists, then L_idx exists too
    for (idx_t i=0; i < 1 + d->mo * (d->mo/2); ++i) {
      mad_free(d->L[i]);
      if (d->L_idx[i]) {
        mad_free(*d->L_idx[i]);  // allocated as single block
        mad_free( d->L_idx[i]);
      }
    }
    mad_free(d->L);
    mad_free(d->L_idx);
  }

  if (d->ocs) {
    int nth = d->nth + (d->nth > 1);
    for (int t=0; t < nth; ++t)
      mad_free(d->ocs[t]);
    mad_free(d->ocs);
  }

  // destroy temporaries
  del_temps(d);

  // remove descriptor from global array
  Ds[d->id] = NULL;
  mad_free(d);
}

void
mad_desc_cleanup (void)
{
  for (idx_t i = 0; i < MAX_TPSA_DESC; ++i)
    if (Ds[i]) mad_desc_del(Ds[i]);
}

// --- end --------------------------------------------------------------------o
