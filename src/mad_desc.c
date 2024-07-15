/*
 o-----------------------------------------------------------------------------o
 |
 | Descriptor (TPSA) module implementation
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: C. Tomoiaga
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - provide a full featured Generalized TPSA package

  Information:
  - parameters ending with an underscope can be null.

  Errors:
  - TODO

 o-----------------------------------------------------------------------------o
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

// --- globals ----------------------------------------------------------------o

// must be global variables for access from LuaJIT FFI.
const  ord_t  mad_tpsa_dflt  = -1;
const  ord_t  mad_tpsa_same  = -2;
       num_t  mad_tpsa_eps   =  1e-40;
       ord_t  mad_tpsa_dbgo  = -1; // effective only with TPSA_DEBUG > 0
       int    mad_tpsa_dbgf  =  0; // effective only with TPSA_DEBUG > 0
       int    mad_tpsa_dbga  =  0; // effective only with TPSA_DEBUG > 0

// last descriptor created or searched (used to create GTPSA when d is NULL)
const desc_t *mad_desc_curr = NULL;

// --- sizes ------------------------------------------------------------------o

static inline ssz_t
max_nc(ssz_t nn, ssz_t mo)
{
  // #coeff(nn,mo) = (nn+mo)! / (nn! mo!)
  ssz_t sum = nn+mo, max = MAX(nn,mo);
  u64_t num = 1, den = 1, lim = ULLONG_MAX/sum;
  for (ssz_t i=max+1; i <= sum && num < lim; ++i) {
    num *= i;
    den *= i - max;
  }
  return num < lim ? num / den : 0;
}

// --- monomials --------------------------------------------------------------o

static inline log_t
mono_isvalid (const D *d, ssz_t n, const ord_t m[n])
{
  assert(d && m);
  ensure(0 <= n && n <= d->nn, "invalid monomial length, 0<= %d <=%d", n,d->nn);
  return mad_mono_le (n, m, d->no)
      && mad_mono_ord(n, m)             <= d->mo
      && mad_mono_ord(n-d->nv, m+d->nv) <= d->po;
}

static inline log_t
mono_isvalidsm (const D *d, ssz_t n, const idx_t m[n])
{
  assert(d && m);
  ensure(0 <= n && n/2 <= d->nn && !(n&1),
                             "invalid monomial length, 0<= %d <=%d", n/2,d->nn);
  idx_t prev = -1;
  ord_t mo = 0, po = 0;

  FOR(i,0,n,2) {
    idx_t idx = m[i]-1; // translate from var idx to mono idx
    ord_t ord = m[i+1];

    // check index
    if (idx >= d->nn || idx < prev) return FALSE;

    // check order
    mo += ord;
    if (mo > d->mo || ord > d->no[idx]) return FALSE;

    // check param x-order
    if (idx > d->nv) {
      po += ord;
      if (po > d-> po) return FALSE;
    }

    // backup index
    prev = idx;
  }
  return TRUE;
}

static inline log_t
mono_nxtbyvar (const D *d, ssz_t n, ord_t m[n])
{
  assert(d && m && n <= d->nn);
  const ord_t *no = d->no;
  FOR(i,n) {
    if (++m[i] <= no[i] && mono_isvalid(d, n, m)) return TRUE;
    m[i] = 0;
  }
  return FALSE; // no more valid monomial
}

static inline idx_t
mono_findidx (ord_t **T_, ssz_t n, const ord_t m[n], idx_t from, idx_t to)
{
  const ord_t *const *T = (const ord_t *const *)T_;
  idx_t start = from, count = to-from, step, i;

  while (count > 0) {
    step = count / 2;
    i = start + step;
    if (mad_mono_rcmp(n, T[i], m) < 0) {
      start = ++i;
      count -= step + 1;
    }
    else
      count = step;
  }
  if (start < to && mad_mono_eq(n, T[start], m))
    return start;

  // error
  printf("** ERR: "); mad_mono_print(n, m, 0,0); printf(" <- not found\n");
  error("monomial not found in [%d,%d) (unexpected)", from, to);
  return -1; // never reached
}

static inline void
mono_realloc (D *d, ssz_t nc)
{
  assert(d);
  if (nc > 10*DESC_WARN_MONO)
    error("gtpsa are too large (%d > %d monomials)", nc, 10*DESC_WARN_MONO);

#if DESC_DEBUG > 1
  printf("desc nc: %d\n", nc);
#endif

  d->nc = nc;
  d->monos = mad_realloc(d->monos, nc*(d->nn * sizeof *d->monos));
}

static inline void
set_monos (D *d) // builds the monomials matrix in Tv order
{
  DBGFUN(->);
  assert(d && d->no);

  int n = d->nn; // TODO: avoid overflow for large np and/or po...
  u64_t nc = max_nc(n, d->mo); // upper bound
  // u64_t nc = max_nc(nv, mo) + max_nc(nv, mo-po) * max_nc(np, po);

  // check for overflow, then start with (6,8)=3003
  d->nc = nc <= 0 || nc > DESC_WARN_MONO ? (ssz_t)max_nc(6,8) : (ssz_t)nc;
  mono_realloc(d, d->nc);

  ord_t m[n];
  mad_mono_fill(n, m, 0);
  mad_mono_copy(n, m, d->monos);

  // fill the matrix
  idx_t i = 1;
  for (; mono_nxtbyvar(d, n, m); ++i) {
    if (i >= d->nc) mono_realloc(d, 2*d->nc);
    mad_mono_copy(n, m, d->monos + n*i);
  }

  // resize the matrix (shrink)
  mono_realloc(d, i);
  d->size += d->nc*d->nn * sizeof *d->monos;

  if (d->nc > DESC_WARN_MONO)
    warn("gtpsa are very large (%d monomials),\n\t"
         "building descriptor may take some time,\n\t"
         "consider that multiplication is quadratic in size.", d->nc);

  DBGFUN(<-);
}

// --- tables -----------------------------------------------------------------o

static inline void
hsm_print (ssz_t n, const idx_t m[n], int l0_)
{
  assert(m);
  ssz_t len = printf("[ ") + l0_;
  FOR(i,0,n,2) {
    len += printf("%d^%2hhu ", m[i]-1, (ord_t)m[i+1]);
    if (len >= 80) { printf("\n"); len = 0; }
  }
  printf("]\n");
}

static inline void
ord_print (ssz_t n, ord_t v[n], int l0_)
{
  assert(v);
  ssz_t len=printf("[ ") + l0_;
  FOR(i,n) {
    len += printf("%d ", v[i]);
    if (len >= 80) { printf("\n"); len = 0; }
  }
  printf("]\n");
}

static inline void
idx_print (ssz_t n, idx_t v[n], int l0_)
{
  assert(v);
  ssz_t len = printf("[ ") + l0_;
  FOR(i,n) {
    len += printf("%d ", v[i]);
    if (len >= 80) { printf("\n"); len = 0; }
  }
  printf("]\n");
}

static inline void
tbl_print (ssz_t n, ssz_t h, ord_t **t) // t[h][n]
{
  assert(t);
  idx_t i=0;
#if DESC_DEBUG > 2
  for (; i < h; ++i) {
    printf("(%2d) ",i); mad_mono_print(n,t[i],0,0); printf(" o=%d\n", mad_mono_ord(n,t[i]));
  }
#else
  for (; i < MIN(h,50); ++i) {
    printf("(%2d) ",i); mad_mono_print(n,t[i],0,0); printf(" o=%d\n", mad_mono_ord(n,t[i]));
  }
  if (h > 100) {
    printf("... [ %d more rows ] ...\n", h - 100);
  }
  for (i=MAX(i,h-50); i < h; ++i) {
    printf("(%2d) ",i); mad_mono_print(n,t[i],0,0); printf(" o=%d\n", mad_mono_ord(n,t[i]));
  }
#endif
}

static inline void
tbl_by_var(D *d)
{
  DBGFUN(->);
  assert(d && d->monos);

  d->Tv = mad_malloc(d->nc * sizeof *d->Tv);
  d->size += d->nc * sizeof *d->Tv;

  FOR(i,d->nc) d->Tv[i] = d->monos + i*d->nn;

#if DESC_DEBUG > 1
  printf("Tv =\n");
  tbl_print(d->nn, d->nc, d->Tv);
#endif
  DBGFUN(<-);
}

static const D *cmp_d = NULL; // not re-entrant...

static int
cmp_mono (const void *a, const void *b)
{
  const D *d = cmp_d;
  idx_t i1 = *(const idx_t*)a;
  idx_t i2 = *(const idx_t*)b;

  int o1 = mad_mono_ord(d->nn, d->Tv[i1]);
  int o2 = mad_mono_ord(d->nn, d->Tv[i2]);
  if (o1 != o2) return o1 - o2;

  return mad_mono_rcmp(d->nn, d->Tv[i1], d->Tv[i2]);
}

static inline void
tbl_by_ord(D *d)
{
  DBGFUN(->);
  assert(d && d->monos);

  d->To      = mad_malloc( d->nc    * sizeof *d->To     );
  d->tv2to   = mad_malloc( d->nc    * sizeof *d->tv2to  );
  d->to2tv   = mad_malloc( d->nc    * sizeof *d->to2tv  );
  d->ords    = mad_malloc( d->nc    * sizeof *d->ords   );
  d->prms    = mad_malloc( d->nc    * sizeof *d->prms   );
  d->ord2idx = mad_malloc((d->mo+2) * sizeof *d->ord2idx);

  d->size +=  d->nc    * sizeof *d->To;
  d->size +=  d->nc*2  * sizeof *d->tv2to; // tv2to + to2tv
  d->size +=  d->nc*2  * sizeof *d->ords;  // ords  + prms
  d->size += (d->mo+2) * sizeof *d->ord2idx;

  assert(!cmp_d);
  cmp_d = d;
  FOR(i,d->nc) d->to2tv[i] = i;
  qsort(d->to2tv, d->nc, sizeof *d->to2tv, cmp_mono);
  cmp_d = NULL;

  d->To[0] = d->monos;
  d->tv2to[0] = 0, d->ord2idx[0] = 0, d->ords[0] = d->prms[0] = 0;

  idx_t j=0;
  FOR(i,1,d->nc) {
    d->tv2to[d->to2tv[i]] = i;
    d->To[i] = d->monos + d->to2tv[i]*d->nn;
    d->ords[i] = mad_mono_ord(d->nn, d->To[i]);
    d->prms[i] = mad_mono_ord(d->np, d->To[i]+d->nv);
    if (d->ords[i] > d->ords[i-1]) d->ord2idx[++j] = i;
//  printf("i=%d, o=%d, p=%d\n", i, d->ords[i], d->prms[i]);
  }
  d->ord2idx[d->mo+1] = d->nc;

#if DESC_DEBUG > 1
  printf("To =\n"); tbl_print(d->nn, d->nc, d->To);
  ord_print(d->nc  , d->ords,    printf("ords = "));
  idx_print(d->mo+2, d->ord2idx, printf("ord2idx = "));
#endif
  DBGFUN(<-);
}

// --- H indexing matrix ------------------------------------------------------o

static inline void
tbl_print_H(const D *d)
{
  assert(d && d->H);
  const idx_t *H = d->H;
  ssz_t nj = d->nn, ni = d->mo+2;

  FOR(j,nj) {
    printf("%2d | ", j);
    FOR(i,ni) printf("%2d ", H[j*ni+i]);
    printf("\n");
  }
}

static inline int
tbl_index_H(const D *d, ssz_t n, const ord_t m[n])
{
//  DBGFUN(->);
  assert(d && n <= d->nn);
  const idx_t *H = d->H;
  ssz_t ni = d->mo+2;
  idx_t I = 0, s = 0;

  // eq. 10 from IPAC'15 paper
  RFOR(j,n) {
    idx_t i0 = j*ni + s;
    idx_t i1 = i0 + m[j];
#if DESC_DEBUG > 0
    assert(m[j] <= d->no[j] && i1 < (d->mo+2)*d->nn);
#endif
    I += H[i1] - H[i0];
    s += m[j];
  }
  if (I < 0) {
    printf("%s: I=%d for monomial ", __func__, I);
    mad_mono_print(n, m, 0,0); printf("\n");
    assert(I > -1);
  }
//  DBGFUN(<-);
  return I;
}

static inline int
tbl_index_Hsm(const D *d, ssz_t n, const idx_t m[n])
{
//  DBGFUN(->);
  assert(d && n/2 <= d->nn && !(n & 1));
  const idx_t *H = d->H;
  ssz_t ni = d->mo+2;
  idx_t I = 0;

  // indexes in m expected to be in ascending order
  for (idx_t i=n-1, j=0, s=0; i > 0; i -= 2) {
    ensure(j <= m[i-1], "sparse monomial must be in ascending indexes");
    j = m[i-1]-1;
    idx_t i0 = j*ni + s;
    idx_t i1 = i0 + m[i];
#if DESC_DEBUG > 0
    assert(m[i] <= d->no[j] && i0 <= i1 && i1 < (d->mo+2)*d->nn);
#endif
    I += H[i1] - H[i0];
    s += m[i];
  }
  if (I < 0) {
    hsm_print(n, m, printf("%s: I=%d for monomial ", __func__, I));
    assert(I > -1);
  }
//  DBGFUN(<-);
  return I;
}

static inline void
tbl_solve_H(D *d)
{
  DBGFUN(->);
  assert(d);
  const idx_t *o2i = d->ord2idx, *o2v = d->to2tv;
  idx_t *H = d->H;
  ord_t **To = d->To;
  ssz_t nj = d->nn, ni = d->mo+2, n = d->nn;
  ord_t m[n];

  // solve unknowns for j=1..nn-2 variables in reverse order (skip 1st and last)
  for (idx_t i, j=nj-2; j > 0; --j) {
    if (H[(j+1)*ni-2] > 0) continue;      // no unknown for this variable
    for (i=ni-2; i > 0 && H[j*ni+i-1] <= 0; --i) ; // find lowest unknown  // 1.
    for (; i < ni && !H[j*ni+i]; ++i) {   // for each unknown..
      mad_mono_copy(n, To[o2i[i]-1], m);  // monomial without the unknown
      m[j]++;                             // add the unknown               // 2.
      if (!mono_isvalid(d, n, m)) {       // monomial blocked by po
        for (; i < ni && !H[j*ni+i]; ++i) H[j*ni+i] = -1;
        break;
      }
      idx_t i0 = tbl_index_H(d, n, m);                                     // 3.
      idx_t i1 = o2v[mono_findidx(To, n, m, o2i[i], o2i[i+1])];            // 4.
      H[j*ni+i] = i1 - i0;                                                 // 5.
    }
  }

#if DESC_DEBUG > 1
  printf("H =\n");
  tbl_print_H(d);
#endif
  DBGFUN(<-);
}

static inline void
tbl_bound_H(D *d)
{
  DBGFUN(->);
  assert(d);
  ssz_t nj = d->nn, ni = d->mo+2;
  idx_t *H = d->H, s = 0;

  RFOR(j,nj) {
    s += d->no[j];           // fill unreacheable orders in H with -1
    FOR(i,1+MIN(s,d->mo),ni) H[j*ni+i] = -1;
  }

#if DESC_DEBUG > 1
  printf("H =\n");
  tbl_print_H(d);
#endif
  DBGFUN(<-);
}

static inline void
tbl_build_H(D *d)
{
  DBGFUN(->);
  assert(d);
  idx_t *H = d->H;
  ord_t **Tv = d->Tv;
  ssz_t nj = d->nn, ni = d->mo+2, nc = d->nc;

  // unit congruence for j=0 (1st) variable
  FOR(i,ni-1) H[i] = i;
  H[ni-1] = 0;                         // complete row with zero

  // congruence for j=1..nn-1 variables (skip 1st)
  idx_t m = 2;
  FOR(j,1,nj) {
    H[j*ni] = 0;                       // first column
    idx_t i = 1;
    for (; m < nc; ++m) {              // scan monomials
      if (Tv[m][j] != Tv[m-1][j]) {    // transition for variable j
        H[j*ni + i++] = m;             // save index in Tv
        if (Tv[m][j] == 0) break;      // congruence for variable j closed
      }
    }
    for (; i < ni; ++i) H[j*ni+i] = 0; // complete row with zeros
  }

#if DESC_DEBUG > 1
  printf("H =\n");
  tbl_print_H(d);
#endif
  DBGFUN(<-);
}

static inline void
tbl_set_H(D *d)
{
  DBGFUN(->);
  assert(d && d->no && d->Tv && d->To);

  d->H = mad_malloc((d->mo+2)*d->nn * sizeof *d->H);
  d->size += d->nn *(d->mo+2) * sizeof *d->H;

  tbl_build_H(d);
  tbl_bound_H(d);
  tbl_solve_H(d);
  DBGFUN(<-);
}

// --- L indexing matrix ------------------------------------------------------o

static inline void
tbl_print_LC(const idx_t *lc, ord_t oa, ord_t ob, const idx_t *o2i)
{
  ssz_t cols = o2i[oa+1] - o2i[oa];
  FOR(ib, o2i[ob], o2i[ob+1]) { printf("\n  ");
  FOR(ia, o2i[oa], o2i[oa+1]) {
    idx_t ic = lc[hpoly_idx(ib-o2i[ob],ia-o2i[oa],cols)];
    printf("%3d ", ic);
  }}
  printf("\n");
}

static inline void
tbl_print_L(const D *d)
{
  ssz_t ho = d->mo/2;
  for (ord_t oc=2; oc <= MIN(d->mo,5); ++oc)
    for (ord_t j=1; j <= oc/2; ++j) {
      ord_t oa = oc-j, ob = j;
      printf("L[%d][%d] = {", ob, oa);
      tbl_print_LC(d->L[oa*ho + ob], oa, ob, d->ord2idx);
    }
  if (d->mo > 5) printf("Orders 5 to %d omitted...\n", d->mo);
}

static inline idx_t*
tbl_build_LC (ord_t oa, ord_t ob, D *d)
{
  DBGFUN(->);
#if DESC_DEBUG > 2
  printf("tbl_set_LC oa=%d ob=%d\n", oa, ob);
#endif
  assert(d && d->To && d->ord2idx && d->tv2to);
  assert(oa < d->mo && ob < d->mo);

  ssz_t nn   = d->nn;
  ord_t **To = d->To, m[nn];
  const idx_t *o2i = d->ord2idx,
            *tv2to = d->tv2to;
  const ssz_t cols = o2i[oa+1] - o2i[oa], // sizes of orders
              rows = o2i[ob+1] - o2i[ob];

  // allocation lc[rows,cols]: lc[ib,ia] = lc[(ib-o2i[ob])*cols + ia-o2i[oa]]
  size_t mat_size = (size_t)rows*cols;

#if DESC_DEBUG > 2
  printf("LC[%d,%d]=%zu index slots\n", rows, cols, mat_size);
#endif

  if (mat_size > INT32_MAX)
    error("gtpsa are too large (%d monomials),\n\t"
          "indexing matrix #slots for orders %d x %d = %zu > 2^31",
          d->nc, oa, ob, mat_size);

  idx_t *lc = mad_malloc(mat_size * sizeof *lc);
  d->size += mat_size * sizeof *lc;

  // initialisation
  for (size_t i=0; i < mat_size; ++i) lc[i] = -1;

  // loop over indexes of order ob
  FOR(ib, o2i[ob], o2i[ob+1]) {
    int lim_a = oa == ob ? ib+1 : o2i[oa+1];   // triangular is lower left

    // loop over indexes of order oa
    FOR(ia, o2i[oa], lim_a) {
      // get the resulting monomial
      mad_mono_add(nn, To[ia], To[ib], m);
      // check for validity
      if (mono_isvalid(d,nn,m)) {
        // get its index in To
        idx_t ic = tv2to[tbl_index_H(d,nn,m)];
        // get its index in lc
        idx_t ilc = hpoly_idx(ib-o2i[ob], ia-o2i[oa], cols);
        // fill lc
        lc[ilc] = ic;
#if DESC_DEBUG > 2
        printf(" ib=%d ", ib); mad_mono_print(nn, To[ib], 0,0);
        printf(" ia=%d ", ia); mad_mono_print(nn, To[ia], 0,0);
        printf(" ic=%d ", ic); mad_mono_print(nn, m     , 0,0);
        printf(" ilc=%d\n", ilc);
#endif
      }
    }
  }

#if DESC_DEBUG > 2
  tbl_print_LC(lc, oa, ob, o2i);
#endif
  DBGFUN(<-);
  return lc;
}

static inline idx_t**
get_LC_idxs (ord_t oa, ord_t ob, D *d)
{
  DBGFUN(->);
  ord_t oc = oa + ob;
  ssz_t ho = d->mo/2;

  const idx_t *o2i = d->ord2idx,
               *lc = d->L[oa*ho + ob];
  const idx_t    T = (o2i[oc+1]+o2i[oc]-1) / 2;  // splitting threshold of oc  (???)
  const ssz_t cols = o2i[oa+1] - o2i[oa],
              rows = o2i[ob+1] - o2i[ob];

  idx_t **LC_idx = mad_malloc(3 * sizeof *LC_idx);
  d->size += 3 * sizeof *LC_idx;

  idx_t *limits = mad_malloc(3*rows * sizeof *limits); // rows: [start split end]
  d->size += 3*rows * sizeof *limits;

  const idx_t START = 0, SPLIT = 1, END = 2;

  LC_idx[START] = limits;
  LC_idx[SPLIT] = limits +   rows;
  LC_idx[END  ] = limits + 2*rows;

  FOR(ib,rows) {
    // shift ia to first valid entry
    idx_t ia = 0;
    for (; ia < cols && lc[hpoly_idx(ib,ia,cols)] == -1; ++ia) ;
    LC_idx[START][ib] = ia;

    // shift ia to last valid entry
    ia = oa == ob ? ib : cols-1;
    for (; ia >= 0 && lc[hpoly_idx(ib,ia,cols)] == -1; --ia) ;
    LC_idx[END  ][ib] = ia + 1;

    LC_idx[SPLIT][ib] = oa == ob ? ib+1 : cols;
    for (ia = LC_idx[START][ib]; ia < LC_idx[END][ib]; ++ia)
      if (lc[hpoly_idx(ib,ia,cols)] >= T) {
        LC_idx[SPLIT][ib] = ia;
        break;
      }
  }

#if DESC_DEBUG > 2
  if (oc <= 5) {
    printf("LC_idx[%d][%d] = { [T=%d]\n", ob, oa, T);
    printf("  -->\t  //\t<--\n");
    FOR(i,rows)
      printf("  [%3d\t%3d\t%3d]\n", LC_idx[START][i],LC_idx[SPLIT][i],LC_idx[END][i]);
  }
#endif

  DBGFUN(<-);
  return LC_idx;
}

static inline void
tbl_set_L (D *d)
{
  DBGFUN(->);
  ssz_t ho = d->mo/2;

  size_t L_sz = (ho*d->mo+1) * sizeof *d->L;
  d->L = mad_malloc(L_sz); memset(d->L, 0, L_sz);
  d->size += L_sz;

  size_t Li_sz = (ho*d->mo+1) * sizeof *d->L_idx;
  d->L_idx = mad_malloc(Li_sz); memset(d->L_idx, 0, Li_sz);
  d->size += Li_sz;

  #ifdef _OPENMP
  if (d->mo > 6) {
    #pragma omp parallel for schedule(guided)
    for (ord_t oc=2; oc <= d->mo; ++oc) {
      for (ord_t j=1; j <= oc/2; ++j) {
        ord_t oa = oc-j, ob = j;
        d->L    [oa*ho + ob] = tbl_build_LC(oa, ob, d);
        d->L_idx[oa*ho + ob] = get_LC_idxs (oa, ob, d);
      }
    }
  } else
  #endif
    for (ord_t oc=2; oc <= d->mo; ++oc) {
      for (ord_t j=1; j <= oc/2; ++j) {
        ord_t oa = oc-j, ob = j;
        d->L    [oa*ho + ob] = tbl_build_LC(oa, ob, d);
        d->L_idx[oa*ho + ob] = get_LC_idxs (oa, ob, d);
      }
    }

#if DESC_DEBUG > 1
  tbl_print_L(d);
#endif
  DBGFUN(<-);
}

// --- descriptor internal checks ---------------------------------------------o

static int
tbl_check_L (D *d)
{
  DBGFUN(->);
  assert(d && d->ord2idx && d->L && d->no && d->To && d->H);
  const idx_t *o2i = d->ord2idx;
  int ho = d->mo / 2;
  ord_t m[d->nn];
  for (int oc = 2; oc <= d->mo; ++oc)
    for (int j = 1; j <= oc / 2; ++j) {
      int oa = oc - j, ob = j;
      idx_t *lc = d->L[oa*ho + ob];
      if (!lc)                                     return  1e7 + oa*1e3 + ob;

      int sa = o2i[oa+1]-o2i[oa], sb = o2i[ob+1]-o2i[ob];

      FOR(ibl,sb) {
        int lim_a = oa == ob ? ibl+1 : sa;
        FOR(ial,lim_a) {
          int ib = ibl + o2i[ob], ia = ial + o2i[oa];
          int il = hpoly_idx(ibl,ial,sa);
          if (il < 0)                              return -2e7 - ia*1e5 - ib;
          if (il >= sa * sb)                       return  2e7 + ia*1e5 + ib;

          int ic = lc[il];
          if (ic >= o2i[oc+1])                     return  3e7 + ic*1e5 + 11;
          if (ic >= 0 && ic < d->ord2idx[oc])      return  3e7 + ic*1e5 + 12;

          mad_mono_add(d->nn, d->To[ia], d->To[ib], m);
          if (ic < 0 && mono_isvalid(d,d->nn,m))   return -3e7          - 13;
        }
      }
    }
  DBGFUN(<-);
  return 0;
}

static int  // error code
tbl_check_H (D *d)
{
  DBGFUN(->);
  const idx_t *tv2to = d->tv2to,
                  *H = d->H,
                  nn = d->nn;
        ord_t   **To = d->To;

  ssz_t nj = d->nn, ni = d->mo+2, nc = d->nc;

  FOR(j,nj) // check for zeros at order 0
    if (H[j*ni+0])                           return 4e6 + j*ni;

  FOR(j,nj) // check for -1 at order mo+1
    if (H[(j+1)*ni-1] != -1)                 return 5e6 + (j+1)*ni-1;

  FOR(i,ni-1) // check for 1..n for first variable
    if (H[i] != i)                           return 6e6 + i;

  FOR(j,1,nj) FOR(i,1,ni) // check for no more zeros
    if (!H[j*ni+i])                          return 7e6 + j*ni+i;

  FOR(i,nc)
    if (tv2to[tbl_index_H(d,nn,To[i])] != i) return 8e6 + i;

  DBGFUN(<-);
  return 0;
}

static int  // error code
tbl_check_T (D *d)
{
  DBGFUN(->);
  const idx_t *tv2to = d->tv2to,
              *to2tv = d->to2tv,
                  nn = d->nn;
        ord_t   **Tv = d->Tv,
                **To = d->To,
              *monos = d->monos;

  FOR(i,d->nc) {
    if (!mad_mono_eq(nn,Tv[i],monos + nn*i)) return 1e6 + i;
    if (!mad_mono_eq(nn,To[tv2to[i]],Tv[i])) return 2e6 + i;
    if (to2tv[tv2to[i]] != i)                return 3e6 + i;
  }
  DBGFUN(<-);
  return 0;
}

// --- thread dispatch --------------------------------------------------------o

static inline void
get_ops(D *d, long long int ops[])
{
  DBGFUN(->);
  const idx_t *o2i = d->ord2idx;
  ops[0] = ops[1] = ops[2] = 0;
  for (ord_t o = 3; o <= d->mo; ++o) {
    ops[o] = 0;
    for (ord_t j = 1; j <= (o-1)/2; ++j) {
      ord_t oa = o-j, ob = j;            // oa > ob >= 1
      idx_t na = o2i[oa+1] - o2i[oa], nb = o2i[ob+1] - o2i[ob];
      ops[o] += 2 * na * nb;
    }
    if (!(o & 1)) {
      ord_t ho = o/2;
      ops[o] += (o2i[ho+1]-o2i[ho]) * (o2i[ho+1]-o2i[ho]);
    }
  }
  ops[d->mo+1] = ops[d->mo]/2;
  ops[d->mo]  -= ops[d->mo+1];
  DBGFUN(<-);
}

static inline int
get_min_dispatched_idx(int nb_threads, long long int dops[])
{
  DBGFUN(->);
  long long int min_disp = dops[nb_threads-1];
  int min_disp_idx = nb_threads - 1;
  RFOR(t,nb_threads)
    if (dops[t] <= min_disp) {
      min_disp = dops[t];
      min_disp_idx = t;
    }
  DBGFUN(<-);
  return min_disp_idx;
}

static inline void
set_thread (D *d)
{
  DBGFUN(->);
  // [0] serial(all), [1..nth] parallel(split)
  int nth = d->nth + (d->nth > 1);

  d->ocs = mad_malloc(nth * sizeof *d->ocs);
  d->size += nth * sizeof *(d->ocs);

  int sizes[nth];
  FOR(t,nth) {
    d->ocs[t] = mad_calloc(d->mo, sizeof *d->ocs[0]);
    d->size += (d->mo) * sizeof *d->ocs[0];
    sizes[t] = 0;
  }

  long long int ops[d->mo+2], dops[nth];
  memset(dops, 0, nth * sizeof *dops);
  get_ops(d, ops);

  // serial
  for (ord_t o = d->mo; o > 2; --o) {
    d->ocs[0][sizes[0]++] = o;
    dops[0] += ops[o];
  }
  dops[0] += ops[d->mo+1];

  // parallel
  if (nth > 1) {
    for (ord_t o = d->mo+1; o > 2; --o) {
      int idx = get_min_dispatched_idx(d->nth,dops+1) + 1;
      assert(idx > 0 && idx <= d->nth);
      d->ocs[idx][sizes[idx]++] = o;
      dops[idx] += ops[o];
    }
  }

#if DESC_DEBUG > 1
  printf("\nTHREAD DISPATCH:\n");
  FOR(t,nth) {
    printf("[%d]: ", t);
    for (idx_t i=0; d->ocs[t][i]; ++i)
      printf("%d ", d->ocs[t][i]);
    printf("[ops:%lld] \n", dops[t]);
  }
  printf("\n");
#endif
  DBGFUN(<-);
}

#if DESC_USE_TMP

static inline void
set_temp (D *d)
{
  DBGFUN(->);
  d->  t = mad_malloc(DESC_MAX_TMP * d->nth * sizeof *d-> t );
  d-> ct = mad_malloc(DESC_MAX_TMP * d->nth * sizeof *d->ct );
  d-> ti = mad_malloc(               d->nth * sizeof *d-> ti);
  d->cti = mad_malloc(               d->nth * sizeof *d->cti);

  FOR(j,d->nth) {
  FOR(i,DESC_MAX_TMP) {
    d-> t[j*DESC_MAX_TMP+i] = mad_tpsa_newd (d,d->mo);
    d->ct[j*DESC_MAX_TMP+i] = mad_ctpsa_newd(d,d->mo); }
    d->ti[j] = d->cti[j] = 0;
  }

#if DESC_DEBUG > 1
  printf("\nTEMPS #TPSA = 2 (R&C) x %d (#TMPS) x %d (Threads) = %d\n"
         "TEMPS TMEM  = %d (TPSA) x %d (nc) = %llu bytes\n",
          DESC_MAX_TMP, d->nth, 2*DESC_MAX_TMP*d->nth, 2*DESC_MAX_TMP*d->nth,
          d->nc, 2*DESC_MAX_TMP*d->nth* 3ull*d->nc*sizeof(num_t)/2);
#endif
  DBGFUN(<-);
}

static inline void
del_temps (D *d)
{
  DBGFUN(->);
  if (d->t) {
    FOR(j,d->nth) {
    FOR(i,DESC_MAX_TMP) {
      mad_tpsa_del (d-> t[j*DESC_MAX_TMP+i]);
      mad_ctpsa_del(d->ct[j*DESC_MAX_TMP+i]);
    }}
    mad_free(d-> t );
    mad_free(d->ct );
    mad_free(d-> ti);
    mad_free(d->cti);
  }
  DBGFUN(<-);
}

#endif // DESC_USE_TMP

// --- descriptor management --------------------------------------------------o

static int desc_max = 0;
static D *Ds[DESC_MAX_ARR];

static inline log_t
desc_equiv (const D *d, int nn, ord_t mo, int np, ord_t po, const ord_t no_[nn])
{
  int same = d->nn == nn && d->mo == mo && d->np == np && (!np || d->po == po);

  if (same) return no_ ? mad_mono_eq (nn   , d->no      , no_) : !d->uno ||
                       ( mad_mono_eqn(nn-np, d->no      , mo ) &&
                         mad_mono_eqn(np   , d->no+nn-np, po ) );
  return 0;
}

static inline log_t
desc_compat (const D *dc, const D *d)
{
  return dc->nn == d->nn && dc->mo == d->mo && dc->np <= d->np &&
        (!d->np || (dc->po == d->po && dc->mo >= d->po)) &&
         mad_mono_eq(d->nn, dc->no, d->no);
}

static inline D*
desc_init (int nn, ord_t mo, int np, ord_t po, const ord_t no_[nn])
{
  DBGFUN(->);
  ensure(mo <= DESC_MAX_ORD, // variables max orders validation
         "invalid gtpsa order exceeds maximum order, %u < %u", mo, DESC_MAX_ORD);

#if DESC_DEBUG > 1
  printf("desc in: nn=%d, mo=%d, np=%d, po=%d\n", nn, mo, np, po);
#endif

  D *d = mad_malloc(sizeof *d); memset(d, 0, sizeof *d);

  d->nn = nn;
  d->nv = nn-np;
  d->np = np;
  d->mo = mo;
  d->po = po;

  ord_t *no = mad_malloc(nn * sizeof *d->no);
  if (no_) {
    mad_mono_copy(nn, no_, no);
    d->uno = !mad_mono_eqn(nn-np, no      , mo) ||
             !mad_mono_eqn(np   , no+nn-np, po);
  } else {
    mad_mono_fill(nn-np, no, mo);
    mad_mono_fill(np, no+nn-np, po);
    d->uno = 0;
  }
  d->no = no;

#if DESC_DEBUG > 1
  printf("desc no: "); mad_mono_print(nn,d->no, 0,0); printf("\n");
#endif

  d->nth  = omp_get_max_threads();
  d->size = 0;

  DBGFUN(<-);
  return d;
}

static D*
desc_build (int nn, ord_t mo, int np, ord_t po, const ord_t no_[nn], log_t share)
{
  DBGFUN(->);
  D *d = desc_init(nn, mo, np, po, no_);
  int err = 0, eid=0;

  D *dc = NULL; // search for compatible descriptor
  if (share)    // never true, disable sharing which is suboptimal vs convert
    FOR(i,desc_max)
      if (Ds[i] && desc_compat(Ds[i], d)) { dc = Ds[i]; break; }

  if (!dc) {
    d->shared = mad_malloc(sizeof *d->shared); *d->shared=0, d->sh=mad_tpsa_dflt;
    set_monos (d);
    tbl_by_var(d);
    tbl_by_ord(d); if (DESC_DEBUG && (err = tbl_check_T(d))) { eid=1; goto error; }
    tbl_set_H (d); if (DESC_DEBUG && (err = tbl_check_H(d))) { eid=2; goto error; }
    tbl_set_L (d); if (DESC_DEBUG && (err = tbl_check_L(d))) { eid=3; goto error; }
    set_thread(d);
  } else {
    d->shared  = dc->shared; ++*d->shared, d->sh = dc->id;
    d->nc      = dc->nc;      // set in set_monos
    d->monos   = dc->monos;
    d->ords    = dc->ords;
    d->To      = dc->To;
    d->Tv      = dc->Tv;
    d->ocs     = dc->ocs;
    d->ord2idx = dc->ord2idx;
    d->tv2to   = dc->tv2to;
    d->to2tv   = dc->to2tv;
    d->H       = dc->H;
    d->L       = dc->L;
    d->L_idx   = dc->L_idx;
    d->prms    = mad_malloc(d->nc * sizeof *d->prms);
    d->size   += d->nc * sizeof *d->prms;

    if (d->np) FOR(i,d->nc) d->prms[i] = mad_mono_ord(d->np, d->To[i]+d->nv);
    else       memset(d->prms, 0, d->nc * sizeof *d->prms);
  }

#if DESC_USE_TMP
  set_temp  (d);
#endif

#if DESC_DEBUG > 1
  printf("desc nc: %d ---- Total desc size: %ld bytes\n", d->nc, d->size);
#endif

  DBGFUN(<-); return d;

error:
  printf(eid==1 ? "** >>>>> T BUG <<<<<\n" : "");
  printf("no= ");  mad_mono_print(d->nn, d->no, 0,0); printf("\n");
  printf("Tv=\n"); tbl_print(d->nn, d->nc, d->Tv);    printf("\n");
  printf("To=\n"); tbl_print(d->nn, d->nc, d->To);    printf("\n");
  idx_print(d->nc  , d->tv2to,   printf("tv2to= "  ));
  idx_print(d->nc  , d->to2tv,   printf("to2tv= "  ));
  ord_print(d->nc  , d->ords ,   printf("ords = "  ));
  idx_print(d->mo+2, d->ord2idx, printf("ord2idx= "));
  assert(--eid);

  printf(eid==1 ? "** >>>>> H BUG <<<<<\n" : "");
  printf("\nH=\n"); tbl_print_H(d);
  printf("\nH consistency reported error %d\n", err);
  assert(--eid);

  printf(eid==1 ? "** >>>>> L BUG <<<<<\n" : "");
  printf("\nL=\n"); tbl_print_L(d);
  assert(0);
}

static inline const D*
get_desc (int nn, ord_t mo, int np, ord_t po, const ord_t no_[nn])
{
  DBGFUN(->);
  FOR(i,desc_max)
    if (Ds[i] && desc_equiv(Ds[i], nn, mo, np, po, no_)) {
      DBGFUN(<-); return mad_desc_curr=Ds[i];
    }

  FOR(i,DESC_MAX_ARR)
    if (!Ds[i]) {
      Ds[i] = desc_build(nn, mo, np, po, no_, FALSE);
      Ds[i]->id = i;
      if (i == desc_max) ++desc_max;
      DBGFUN(<-); return mad_desc_curr=Ds[i];
    }

  error("Too many descriptors in concurrent use (max %d)", DESC_MAX_ARR);
}

// --- public -----------------------------------------------------------------o

log_t
mad_desc_isvalidm (const D *d, ssz_t n, const ord_t m[n])
{
  DBGFUN(->);
  assert(d && m);
  const log_t ret = 0 <= n && n <= d->nn && mono_isvalid(d, n, m);
  DBGFUN(<-); return ret;
}

log_t
mad_desc_isvalids (const D *d, ssz_t n, str_t s)
{
  DBGFUN(->);
  assert(d && s);
  if (n <= 0) n = strlen(s);

  ord_t m[n];
  n = mad_mono_str(n, m, s); // n can be shrinked by '\0'
  const log_t ret = 0 <= n && n <= d->nn && mono_isvalid(d, n, m);
  DBGFUN(<-); return ret;
}

log_t
mad_desc_isvalidsm (const D *d, ssz_t n, const idx_t m[n])
{
  DBGFUN(->);
  assert(d && m);
  const log_t ret = 0 <= n && n <= d->nn && mono_isvalidsm(d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
mad_desc_nxtbyvar (const D *d, ssz_t n, ord_t m[n])
{
  DBGFUN(->);
  assert(d && m);

  if (!mono_isvalid(d,n,m)) { DBGFUN(<-); return -1; }

  idx_t idx = tbl_index_H(d,n,m)+1;

  if (d->sh != mad_tpsa_dflt)
    for (; idx < d->nc; idx++) if (mono_isvalid(d,n,d->Tv[idx])) break;

  if (idx == d->nc) { DBGFUN(<-); return -1; }

  mad_mono_copy(n, d->Tv[idx], m);
  DBGFUN(<-); return idx;
}

idx_t
mad_desc_nxtbyord (const D *d, ssz_t n, ord_t m[n])
{
  DBGFUN(->);
  assert(d && m);

  if (!mono_isvalid(d,n,m)) { DBGFUN(<-); return -1; }

  idx_t idx = d->tv2to[tbl_index_H(d,n,m)]+1;

  if (d->sh != mad_tpsa_dflt)
    for (; idx < d->nc; idx++) if (mono_isvalid(d,n,d->To[idx])) break;

  if (idx == d->nc) { DBGFUN(<-); return -1; }

  mad_mono_copy(n, d->To[idx], m);
  DBGFUN(<-); return idx;
}

ord_t
mad_desc_mono (const D *d, idx_t i, ssz_t n, ord_t m_[n], ord_t *p_)
{
  DBGFUN(->);
  assert(d);
  ensure(0 <= i && i < d->nc, "index out of bounds");
  if (m_ && n > 0) mad_mono_copy(MIN(n,d->nn), d->To[i], m_);
  if (p_) *p_ = d->prms[i];
  const ord_t ret = d->ords[i];
  DBGFUN(<-); return ret;
}

idx_t
mad_desc_idxm (const D *d, ssz_t n, const ord_t m[n])
{
  DBGFUN(->);
  assert(d && m);
  const idx_t ret = mono_isvalid(d,n,m) ? d->tv2to[tbl_index_H(d,n,m)] : -1;
  DBGFUN(<-); return ret;
}

idx_t
mad_desc_idxs (const D *d, ssz_t n, str_t s)
{
  assert(d && s); DBGFUN(->);
  if (n <= 0) n = strlen(s);

  ord_t m[n];
  n = mad_mono_str(n, m, s); // n can be shrinked by '\0'
  const idx_t ret = mad_desc_idxm(d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
mad_desc_idxsm (const D *d, ssz_t n, const idx_t m[n])
{
  assert(d && m); DBGFUN(->);
  const idx_t ret = mono_isvalidsm(d,n,m) ? d->tv2to[tbl_index_Hsm(d,n,m)] : -1;
  DBGFUN(<-); return ret;
}

int
mad_desc_getnv (const D *d, ord_t *mo_, int *np_, ord_t *po_)
{
  assert(d);
  if (mo_) *mo_ = d->mo;
  if (np_) *np_ = d->np;
  if (po_) *po_ = d->po;
  return d->nv;
}

ord_t
mad_desc_maxord (const D *d, int n, ord_t no_[n])
{
  assert(d); DBGFUN(->);
  if (no_) {
    ensure(0 <= n && n <= d->nn, "invalid monomial length, 0<= %d <=%d", n,d->nn);
    mad_mono_copy(n, d->no, no_);
  }
  DBGFUN(<-); return d->mo;
}

ssz_t
mad_desc_maxlen (const D *d, ord_t mo)
{
  assert(d); DBGFUN(->);
  if (mo == mad_tpsa_dflt) mo = d->mo; else
  ensure(mo <= d->mo, "invalid order %d (exceeds maximum order %d)", mo,d->mo);
  DBGFUN(<-); return d->ord2idx[mo+1];
}

void
mad_desc_paropsth (const D *d, ssz_t *mult_, ssz_t *comp_)
{
  assert(d); DBGFUN(->);
  D *dd = (D*)d;
  ssz_t t;
  if (mult_) SWAP(dd->pmul , *mult_, t);
  if (comp_) SWAP(dd->pcomp, *comp_, t);
  DBGFUN(<-);
}

void
mad_desc_info (const D *d, FILE *fp_)
{
  assert(d); DBGFUN(->);
  char s[d->nn+1];
  fprintf(fp_ ? fp_ : stdout,
          "id=%d, nn=%d, nv=%d, np=%d, mo=%d, po=%d, uno=%d, no=[%s]\n",
           d->id, d->nn, d->nv, d->np, d->mo, d->po, d->uno,
           mad_mono_prt(d->nn, d->no, s));
  DBGFUN(<-);
}

// --- ctors, dtor ------------------------------------------------------------o

const D*
mad_desc_newv (int nv, ord_t mo)
{
  DBGFUN(->);
  ensure(0 < nv && nv <= DESC_MAX_VAR,
         "invalid #variables, 0< %d <=%d", nv, DESC_MAX_VAR);
  ensure(0 < mo && mo <= DESC_MAX_ORD,
         "invalid maximum order, 0< %d <=%d", mo, DESC_MAX_ORD);

#if DESC_DEBUG > 1
  printf(">> nv=%d, mo=%d\n", nv, mo);
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnonnull"
  const desc_t* ret = get_desc(nv, mo, 0, 0, NULL);
#pragma GCC diagnostic pop

  DBGFUN(<-); return ret;
}

const desc_t*
mad_desc_newvp(int nv, ord_t mo, int np_, ord_t po_)
{
  if (np_ <= 0) return mad_desc_newv(nv, mo);

  DBGFUN(->);
  int np = MAX(np_,0);
  int nn = nv+np;
  ensure(0 < nn && nn <= DESC_MAX_VAR,
         "invalid #variables+#parameters, 0< %d <=%d", nn, DESC_MAX_VAR);
  ensure(0 < mo && mo <= DESC_MAX_ORD,
         "invalid maximum order, 0< %d <=%d", mo, DESC_MAX_ORD);

  ord_t po = MAX(po_,1);
  ensure(0 < po && po <= mo, "invalid parameter order, 0< %d <=%d", po, mo);

#if DESC_DEBUG > 1
  printf(">> nn=%d, mo=%d, np=%d, po=%d[%d]\n", nn, mo, np, po, po_);
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wnonnull"
  const desc_t* ret = get_desc(nn, mo, np, po, NULL);
#pragma GCC diagnostic pop

  DBGFUN(<-); return ret;
}

const desc_t*
mad_desc_newvpo(int nv, ord_t mo, int np_, ord_t po_, const ord_t no_[nv+np_])
{
  if (!no_) return mad_desc_newvp(nv, mo, np_, po_);

  DBGFUN(->);
  int np = MAX(np_,0);
  int nn = nv+np;
  ensure(0 < nn && nn <= DESC_MAX_VAR,
         "invalid #variables+#parameters, 0< %d <=%d", nn, DESC_MAX_VAR);
  ensure(mad_mono_min(nn, no_) > 0,
         "some variables (or parameters) have invalid zero order");

  ord_t mo_ = mad_mono_max(nn, no_); mo = MAX(mo, mo_);
  ensure(0 < mo && mo <= DESC_MAX_ORD,
         "invalid maximum order, 0< %d <=%d", mo, DESC_MAX_ORD);

  ord_t po = MAX(po_,1);
  if (np) {
    ord_t po_ = mad_mono_max(np, no_+nv); po = MAX(po, po_);
    ensure(0 < po && po <= mo, "invalid parameter order, 0< %d <=%d", po, mo);
  }

#if DESC_DEBUG > 1
  printf(">> nn=%d, mo=%d, np=%d, po=%d[%d]\n", nn, mo, np, po, po_);
#endif

  const desc_t* ret = get_desc(nn, mo, np, po, no_);
  DBGFUN(<-); return ret;
}

static void
mad_desc_cleanup (void)
{
  FOR(i,desc_max) if (Ds[i]) mad_desc_del(Ds[i]);
}

void
mad_desc_del (const D *d_)
{
  if (d_ == NULL) { mad_desc_cleanup(); return; }

  DBGFUN(->);
  D *d = (void*)d_;

  mad_free((void*)d->no);
  mad_free(d->prms);

  if (*d->shared > 0) --*d->shared;
  else {
    mad_free(d->shared);
    mad_free(d->monos);
    mad_free(d->ords);
    mad_free(d->To);
    mad_free(d->Tv);
    mad_free(d->ord2idx);
    mad_free(d->tv2to);
    mad_free(d->to2tv);
    mad_free(d->H);

    if (d->L) {  // if L exists, then L_idx exists too
      FOR(i, 1+d->mo*(d->mo/2)) {
        mad_free(d_->L[i]);
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
      FOR(t,nth) mad_free(d->ocs[t]);
      mad_free(d->ocs);
    }
  }

#if DESC_USE_TMP
  del_temps(d); // destroy temporaries
#endif

  // remove descriptor from global array
  if (d == mad_desc_curr) mad_desc_curr = NULL;
  if (d->id+1 == desc_max) {
    int i = d->id;
    while (i > 0 && !Ds[i-1]) --i;
    desc_max = i;
  }
  Ds[d->id] = NULL;
  mad_free(d);
  DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
