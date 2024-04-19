/*
 o-----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module implementation
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

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <limits.h>

#include "mad_cst.h"
#include "mad_mem.h"
#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

/* --- definition --------------------------------------------------------------
  GTPSA are *defined* in [0] U [lo,hi] with 0 < lo <= hi <= mo <= d->mo <= 253
    - scalar GTPSA have lo=1, hi=0, the only case where hi=0 and i.e. lo > hi.
    - new/clear/reset GTPSA are scalar GTPSA with coef[0]=0 (see reset0).
    - always true: coef[0] is defined and lo >= 1, even if coef[0] != 0.
------------------------------------------------------------------------------*/

// --- debugging --------------------------------------------------------------o

// histogram by deciles 0..10
static long long ratio_nz[11] = {0};
static long long ratio_nn[11] = {0};

static inline num_t
ratio (const T *t, num_t eps_)
{
  if (!t->hi) return -1;
  num_t eps = eps_ ? eps_ : mad_tpsa_eps;

  long  nz = 0;
  TPSA_SCAN(t) if (fabs(t->coef[i]) >= eps) ++nz;
  long  nn = o2i[t->hi+1] - o2i[t->lo];
  num_t rr = (num_t)nz / nn;

  int i = rr*10;
  assert(0 <= i && i <= 10);
  ratio_nz[i] += nz;
  ratio_nn[i] += nn;
  return rr;
}

static inline log_t
check (const T *t, ord_t *o_, num_t *r_)
{
  ord_t _o = 0;

  if (!t->d || t->mo > t->d->mo || t->hi > t->mo ||
      (t->lo > t->hi && t->lo != 1)) goto ret;

  if (t->hi) {
    if (FUN(nzero0)(t,t->lo,t->lo) < 0) {_o = t->lo; goto ret;}
    if (FUN(nzero0)(t,t->hi,t->hi) < 0) {_o = t->hi; goto ret;}
  }

  if (r_) *r_ = ratio(t,0);
  return TRUE;

ret:
  if (o_) *o_ = _o;
  return FALSE;
}

/* Enabled only if TPSA_DEBUG > 0 (compile time)
   mad_tpsa_dbga  effects (runtime)
         0        none
         1        check ;
         2        check ; density
         3        check ; density ; print header +tpsa(ok) or +raw(!ok)
*/

int
FUN(debug) (const T *t, str_t name_, str_t fname_, int line_, FILE *stream_)
{
  assert(t);
  ord_t o; num_t r = 0;
  log_t ok = check(t,&o, mad_tpsa_dbga >= 2 ? &r : 0);  // dbga = 2..3 -> ratio

  if (ok && mad_tpsa_dbga <= 2) return ok;              // dbga = 1..2 -> no prn

  const D* d = t->d;
  if (!stream_) stream_ = stdout;

  fprintf(stream_, "%s:%d: '%s' { lo=%d hi=%d mo=%d ao=%d uid=%d did=%d",
          fname_ ? fname_ : "??", line_, name_ ? name_ : "?",
          t->lo, t->hi, t->mo, t->ao, t->uid, d ? d->id : -1);

  if (ok) {
    fprintf(stream_," r=%.2f } 0x%p\n", r, (void*)t); fflush(stream_);

    char name[48];
    strncpy(name, name_ ? name_ : t->nam, 48); name[47] = '\0';
    FUN(print)(t, name, 1e-40, 0, stream_);
    return ok;
  }

  fprintf(stream_," ** bug @ o=%d } 0x%p\n", o, (void*)t); fflush(stream_);

  if (d) {
    const idx_t *o2i = d->ord2idx;
    idx_t ni = o2i[MIN(t->mo,t->ao)+1]; // corrupted TPSA cannot use print
    FOR(i,ni) fprintf(stream_," [%d:%d]=" FMT "\n", i,d->ords[i],VAL(t->coef[i]));
    fprintf(stream_,"\n"); fflush(stream_);
  }

  ensure(ok, "corrupted TPSA detected"); return ok;
}

#ifndef MAD_CTPSA_IMPL
void
mad_tpsa_clrdensity (void)
{
  FOR(i,11) ratio_nz[i] = ratio_nn[i] = 0;
}

void
mad_tpsa_prtdensity (FILE *stream_)
{
  if (!stream_) stream_ = stdout;
  long long sum_nn = 0, sum_nz = 0;
  FOR(i,11) sum_nz += ratio_nz[i], sum_nn += ratio_nn[i];
  if (!sum_nn) { fprintf(stream_,"no tpsa density available.\n"); return; }

  fprintf(stream_,"tpsa average density:\n");
  FOR(i,11) {
    ensure(ratio_nz[i] <= ratio_nn[i], "unexpect ratio > 1");
    fprintf(stream_,"i=%2d: nz=%15lld, nn=%15lld, r=%6.2f, p=%6.2f%%\n",
      i, ratio_nz[i], ratio_nn[i],
      ratio_nn[i] ? 10*(num_t)ratio_nz[i]/ratio_nn[i] : 0,
      100*(num_t)ratio_nn[i]/sum_nn);
  }
  fprintf(stream_,"total nz=%15lld, nn=%15lld, r=%6.2f\n",
          sum_nz, sum_nn, 10*(num_t)sum_nz/sum_nn);
  fflush(stream_);
}
#endif

// --- introspection ----------------------------------------------------------o

const D*
FUN(desc) (const T *t)
{
  assert(t);
  return t->d;
}

ord_t
FUN(mo) (T *t, ord_t mo_)
{
  assert(t);
  if (mo_ == mad_tpsa_same) return t->mo;

  ord_t ret = t->mo;
  if (mo_ < t->mo)
    t->lo = MIN(t->lo, mo_), t->hi = MIN(t->hi, mo_), t->mo = mo_;
  else
    t->mo = MIN(mo_, t->ao);
  return ret;
}

int32_t
FUN(uid) (T *t, int32_t uid_) // set uid if uid != 0
{
  assert(t);
  int32_t ret = t->uid;
  if (uid_) t->uid = uid_;
  return ret;
}

str_t
FUN(nam) (T *t, str_t nam_)
{
  assert(t);
  if (nam_) strncpy(t->nam, nam_, NAMSZ), t->nam[NAMSZ-1] = 0;
  return t->nam;
}

ssz_t
FUN(len) (const T *t)
{
  assert(t);
  return t->d->ord2idx[t->mo+1];
}

ord_t
FUN(ord) (const T *t, log_t hi_)
{
  assert(t);
  return hi_ ? t-> hi : t->mo;
}

log_t
FUN(isnul) (const T *t)
{
  assert(t);
  return !(t->hi || t->coef[0]);
}

log_t
FUN(isvalid) (const T *t)
{
  assert(t);
  return check(t,0,0);
}

num_t
FUN(density) (const T *t, num_t eps)
{
  assert(t);
  return ratio(t,eps);
}

// --- init (unsafe) ----------------------------------------------------------o

T*
FUN(init) (T *t, const D *d, ord_t mo)
{
  assert(t && d && mo <= d->mo); DBGFUN(->);
  t->d = d, t->ao = t->mo = mo, t->uid = 0, t->nam[0] = 0, FUN(reset0)(t);
  DBGTPSA(t); DBGFUN(<-); return t;
}

// --- ctors, dtor ------------------------------------------------------------o

T*
FUN(newd) (const D *d, ord_t mo)
{
  assert(d); DBGFUN(->);
  mo = MIN(mo, d->mo);
  T *r = mad_malloc(sizeof(T) + d->ord2idx[mo+1] * sizeof(NUM)); assert(r);
  r->d = d, r->ao = r->mo = mo, r->uid = 0, r->nam[0] = 0, FUN(reset0)(r);
  DBGTPSA(r); DBGFUN(<-); return r;
}

T*
FUN(new) (const T *t, ord_t mo)
{
  assert(t);
  return FUN(newd)(t->d, mo == mad_tpsa_same ? t->mo : mo);
}

void
FUN(del) (const T *t)
{
  DBGFUN(->);
  mad_free((void*)t);
  DBGFUN(<-);
}

// --- clear, setvar, setprm, setval, update ----------------------------------o

void
FUN(clear) (T *t)
{
  assert(t); DBGFUN(->);
  FUN(reset0)(t);
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(update) (T *t)
{
  assert(t); DBGFUN(->);
  if (t->hi) {
    const ord_t *ords = t->d->ords;
    idx_t j;
    if ((j=FUN(nzero0 )(t,t->lo,        t->hi)) > 0 &&
        (j=FUN(nzero0r)(t,t->lo=ords[j],t->hi)) > 0) t->hi=ords[j];
    else t->lo=1, t->hi=0;
  }
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setvar) (T *t, NUM v, idx_t iv, NUM scl)
{
  assert(t); DBGFUN(->);
  ensure(t->mo >= 1, "variables must have an order >= 1, got %d", t->mo);
  ensure(0 < iv && iv <= t->d->nv,
         "index 1<= %d <=%d is not a GTPSA variable", iv, t->d->nv);

  t->lo = t->hi = 1, t->coef[0] = v;
  FUN(clear0)(t,1,1), t->coef[iv] = scl ? scl : 1;
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setprm) (T *t, NUM v, idx_t ip)
{
  assert(t); DBGFUN(->);
  ensure(t->mo == 1 , "parameters must be a GPTSA of order 1, got %d", t->mo);
  ensure(0 < ip && ip <= t->d->np,
         "index 1<= %d <=%d is not a GPTSA parameter", ip, t->d->np);

  t->lo = t->hi = 1, t->coef[0] = v;
  FUN(clear0)(t,1,1), t->coef[ip+t->d->nv] = 1;
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setval) (T *t, NUM v)
{
  assert(t); DBGFUN(->);
  t->lo = 1, t->hi = 0, t->coef[0] = v;
  DBGTPSA(t); DBGFUN(<-);
}

// --- copy, convert, swap ----------------------------------------------------o

void
FUN(copy) (const T *t, T *r)
{
  assert(t && r); DBGFUN(->);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", t->d, r->d);

  if (t != r) {
    FUN(copy0)(t, r);
    r->coef[0] = t->coef[0];
    TPSA_SCAN(r) r->coef[i] = t->coef[i];
  }
  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(sclord) (const T *t, T *r, log_t inv, log_t prm)
{
  assert(t && r); DBGFUN(->);

  FUN(copy)(t,r);

  const int    np = !prm; // do not consider parameters orders
  const ord_t *po = r->d->prms, lo = MAX(r->lo,2);
  if (inv) { // scale coefs by 1/order
    TPSA_SCAN_Z(r,lo,r->hi) TPSA_SCAN_O(r) r->coef[i] /= o - po[i] * np;
  } else {   // scale coefs by order
    TPSA_SCAN_Z(r,lo,r->hi) TPSA_SCAN_O(r) r->coef[i] *= o - po[i] * np;
  }
  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(clrord) (T *t, ord_t o)
{
  assert(t); DBGFUN(->);
  idx_t j = 0;
  if (!o) t->coef[0] = 0;                                                         else
  if (o  > t->lo && o < t->hi) FUN(clear0)(t, o, o);                              else
  if (o == t->lo && (j=FUN(nzero0 )(t,t->lo+1,t->hi)) > 0) t->lo = t->d->ords[j]; else
  if (o == t->hi && (j=FUN(nzero0r)(t,t->lo,t->hi-1)) > 0) t->hi = t->d->ords[j]; else
  if (j < 0) t->lo = 1, t->hi = 0;
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(getord) (const T *t, T *r, ord_t o)
{
  assert(t && r); DBGFUN(->);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", t->d, r->d);

  if (o < t->lo || o > MIN(t->hi, r->mo)) {
    FUN(setval)(r, o ? 0 : t->coef[0]);
    DBGFUN(<-); return;
  }

  r->lo = r->hi = o, r->coef[0] = 0;
  if (t != r) { TPSA_SCAN_O(r, o) r->coef[i] = t->coef[i]; }

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(cutord) (const T *t, T *r, int o)
{
  assert(t && r); DBGFUN(->);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", t->d, r->d);

  if (o <= 0) {    // cut 0..|o|, see copy0 with t->lo = |o|+1
    r->lo = 1-o;             // min 1 -> keep 1..
    r->hi = MIN(t->hi, r->mo);
    r->coef[0] = 0;
  } else {         // cut |o|..mo, see copy0 with t->hi = |o|-1
    r->lo = t->lo;
    r->hi = MIN(o-1, r->mo); // min 0 -> cut 1..
    r->coef[0] = t->coef[0];
  }

  if (r->lo > r->hi) FUN(setval)(r, r->coef[0]); else
  if (r != t) { TPSA_SCAN(r) r->coef[i] = t->coef[i]; }

  DBGTPSA(r); DBGFUN(<-);
}

idx_t
FUN(maxord) (const T *t, ssz_t n, idx_t idx_[n])
{
  assert(t); DBGFUN(->);

  if (idx_) { FOR(i,n) idx_[i] = -1; idx_[0] = 0; }

  idx_t mi = 0;                       // idx of mv
  num_t mv = fabs(t->coef[0]);        // max of all values
  ord_t hi = MIN(n-1, t->hi);
  TPSA_SCAN_Z(t,t->lo,hi) {
    num_t mo = 0;                     // max of this order
    TPSA_SCAN_O(t)
      if (mo < fabs(t->coef[i])) {
        mo = fabs(t->coef[i]);        // save max for this order
        if (idx_) idx_[o] = i;        // save idx for this order
        if (mv < mo) mv = mo, mi = i; // save max and idx for all orders
      }
  }
  DBGFUN(<-); return mi;
}

idx_t
FUN(cycle) (const T *t, idx_t i, ssz_t n, ord_t m_[n], NUM *v_)
{
  assert(t); DBGFUN(->);
  if (++i <= 0 && t->coef[0]) { i = 0; goto ret; }

  const idx_t *o2i = t->d->ord2idx;
  for (i=MAX(i,o2i[t->lo]); i < o2i[t->hi+1]; i++) if (t->coef[i]) goto ret;
  DBGFUN(<-); return -1;

ret:
  if (v_) *v_ = t->coef[i];
  if (m_) mad_mono_copy(MIN(n,t->d->nn), t->d->To[i], m_);
  DBGFUN(<-); return i;
}

void
FUN(convert) (const T *t, T *r_, ssz_t n, idx_t t2r_[n], int pb)
{
  assert(t && r_); DBGFUN(->);
  ensure(pb >= -1 && pb <= 1,
         "invalid Poisson bracket direction %d, {-1, 0, 1} expected", pb);
  const D *td = t->d, *rd = r_->d;

  if (td == rd && !t2r_) { FUN(copy)(t,r_); DBGFUN(<-); return; }

  T *r = t == r_ ? GET_TMPX(r_) : FUN(reset0)(r_);
  ssz_t rn = rd->nv, tn = td->nv;
  ord_t rm[rn], tm[tn];
  idx_t t2r[tn]; // if t2r[i]>=0 then rm[t2r[i]] = tm[i] for i=0..tn-1
  int   pbs[tn];

  // prepare t2r and pbs
  { idx_t i = 0;
    if (!t2r_)
      for (; i < MIN(tn,rn); ++i) t2r[i] = i, pbs[i] = 0; // identity
    else
      for (; i < MIN(tn, n); ++i) {
        t2r[i] = t2r_[i] >= 0 && t2r_[i] < rn ? t2r_[i] : -1; // -1 discard var
        pbs[i] = pb*(t2r[i]-i)%2 < 0; // pb sign, ignored for discarded vars
      } // fromptc: pt => 1*(6-5)<0=0, t => 1*(5-6)<0=1, x,px,y,py => 1*(i-i)<0=0
    for (; i < tn; i++) t2r[i] = -1;  // discard remaining vars
  }

  // convert t -> r
  const ord_t *ords = rd->ords;
  const ord_t t_hi = MIN(t->hi, r->mo);
  r->coef[0] = t->coef[0];
  TPSA_SCAN(t,t_hi) {
    if (!t->coef[i]) continue;
    mad_desc_mono(td, i, tn, tm, NULL);             // get tm mono at index i
    mad_mono_fill(rn, rm, 0);
    int sgn = 0;
    FOR(j,tn) {                                     // set rm mono
      if (t2r[j] < 0 && tm[j]) goto skip;           // discard coef
      rm[t2r[j]] = tm[j];                           // translate tm to rm
      sgn = sgn - pbs[j] * (tm[j] & 1);             // poisson bracket
    }
    idx_t ri = mad_desc_idxm(rd, rn, rm);           // get index ri of mono rm
#if TPSA_DEBUG > 2
    printf("cvt %d -> %d %c : ", i+1, ri+1, i==ri?' ' : SIGN1(sgn%2)<0?'-':'+');
    mad_mono_print(tn, tm, 0); printf(" -> "); mad_mono_print(rn, rm, 0);
#ifndef MAD_CTPSA_IMPL
    printf(" : %-.16e\n", t->coef[i]);
#else
    printf(" : %-.16e%+.16ei\n", creal(t->coef[i]), cimag(t->coef[i]));
#endif
#endif
    if (ri >= 0) {
      r->coef[ri] = SIGN1(sgn%2)*t->coef[i];
      if (ri) {
        ord_t o = ords[ri];
        if (o < r->lo) r->lo = o;
        if (o > r->hi) r->hi = o;
      }
    }
  skip: ;
  }

  if (r != r_) { FUN(copy)(r,r_); REL_TMPX(r); } else
  DBGTPSA(r_); DBGFUN(<-);
}

// --- indexing / monomials ---------------------------------------------------o

ord_t
FUN(mono) (const T *t, idx_t i, ssz_t n, ord_t m_[n], ord_t *p_)
{
  assert(t);
  ord_t ret = mad_desc_mono(t->d, i, n, m_, p_);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxs) (const T *t, ssz_t n, str_t s)
{
  assert(t && s);
  idx_t ret = mad_desc_idxs(t->d, n, s);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m);
  idx_t ret = mad_desc_idxm(t->d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxsm) (const T *t, ssz_t n, const idx_t m[n])
{
  assert(t && m); DBGFUN(->);
  idx_t ret = mad_desc_idxsm(t->d, n, m);
  DBGFUN(<-); return ret;
}

// --- getters ----------------------------------------------------------------o

static inline NUM
geti (const T *t, idx_t i)
{
  const ord_t o = t->d->ords[i];
  return !o || (t->lo <= o && o <= t->hi) ? t->coef[i] : 0;
}

NUM
FUN(geti) (const T *t, idx_t i)
{
  assert(t);
  if (!i) return t->coef[0];

  DBGFUN(->);
  ensure(0 < i && i < t->d->nc, "index %d out of bounds", i);
  NUM ret = geti(t,i);
  DBGFUN(<-); return ret;
}

NUM
FUN(gets) (const T *t, ssz_t n, str_t s)
{ // --- mono is a string; represented as "[0-9A-Z]*"
  assert(t && s); DBGFUN(->);
  idx_t i = mad_desc_idxs(t->d,n,s);
  ensure(i >= 0, "invalid monomial");
  NUM ret = geti(t,i);
  DBGFUN(<-); return ret;
}

NUM
FUN(getm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGFUN(->);
  idx_t i = mad_desc_idxm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  NUM ret = geti(t,i);
  DBGFUN(<-); return ret;
}

NUM
FUN(getsm) (const T *t, ssz_t n, const idx_t m[n])
{ // --- mono is sparse; represented as [(i,o)]
  assert(t && m); DBGFUN(->);
  idx_t i = mad_desc_idxsm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  NUM ret = geti(t,i);
  DBGFUN(<-); return ret;
}

/* getv cases
   0   1     lo=2      hi=3        mo=4
  [.|.....|........|..........|............]
    |i000n|        |          |               lo=2, hi=1
    | i0n |        |          |               lo=2, hi=1
    | i00n|        |          |               lo=2, hi=1
    | i000|....n   |          |               lo=2, hi=2
    | i000|........|..........|0000000n       lo=2, hi=3
    |     |  i.....|......n   |               lo=2, hi=3
    |     |        |   i......|0000000n       lo=2, hi=3
    |     |        |          |i000000n       lo=2, hi=3
    |     |        |          |  i0000n       lo=2, hi=3
    |     |        |          |i0000000000n   lo=2, hi=3
*/

void
FUN(getv) (const T *t, idx_t i, ssz_t n, NUM v[n])
{
  if (n <= 0) return;
  assert(t && v); DBGFUN(->);
  ssz_t nn = i+n;
  ensure(0 <= i && nn <= t->d->nc, "indexes %d:%d out of bounds", i, nn);

  const ord_t *ord = t->d->ords;
  const idx_t *o2i = t->d->ord2idx;
  ord_t lo = t->lo;
  ord_t hi = MIN(ord[nn-1], t->hi);

  ssz_t n0 = MIN(o2i[lo  ], nn);
  ssz_t nj = MAX(o2i[lo  ], i );
  ssz_t ni = MIN(o2i[hi+1], nn);

//printf("getv: i=%d, n=%d, lo=%d, hi=%d, n0=%d, ni=%d, nj=%d, nn=%d %c\n",
//              i   , n   , lo   , hi   , n0   , ni   , nj   , nn,
//              ni == i+n ? ' ' : '*');

  FOR(j, i,n0) v[j-i] = 0;
  FOR(j,nj,ni) v[j-i] = t->coef[j];
  FOR(j,ni,nn) v[j-i] = 0;

  if (!i) v[0] = t->coef[0];

  DBGFUN(<-);
}

// --- setters ----------------------------------------------------------------o

void
FUN(seti) (T *t, idx_t i, NUM a, NUM b)
{
  assert(t);
  if (!i) { t->coef[0] = a*t->coef[0]+b; return; }

  DBGFUN(->);
  const D *d = t->d;
  ensure(0 < i && i < d->nc, "index %d out of bounds", i);

  // discard value
  const ord_t o = d->ords[i];
  if (o > t->mo) { DBGFUN(<-); return; }

  NUM v = t->lo <= o && o <= t->hi ? a*t->coef[i]+b : b;

  if (v) {
    if (   !t->hi) { FUN(clear0)(t,o,o); t->lo = t->hi = o; } else
    if (o < t->lo) { FUN(clear0)(t,o,t->lo-1);   t->lo = o; } else
    if (o > t->hi) { FUN(clear0)(t,t->hi+1,o);   t->hi = o; }
    t->coef[i] = v;
  } else {
    t->coef[i] = 0;
    idx_t j = 0;
    if (o == t->lo && (j=FUN(nzero0 )(t,t->lo,t->hi)) > 0) t->lo = d->ords[j]; else
    if (o == t->hi && (j=FUN(nzero0r)(t,t->lo,t->hi)) > 0) t->hi = d->ords[j]; else
    if (j <= 0) t->lo = 1, t->hi = 0;
  }
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(sets) (T *t, ssz_t n, str_t s, NUM a, NUM b)
{
  assert(t && s); DBGFUN(->);
  idx_t i = mad_desc_idxs(t->d,n,s);
  ensure(i >= 0, "invalid monomial");
  FUN(seti)(t,i,a,b);
  DBGFUN(<-);
}

void
FUN(setm) (T *t, ssz_t n, const ord_t m[n], NUM a, NUM b)
{
  assert(t && m); DBGFUN(->);
  idx_t i = mad_desc_idxm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  FUN(seti)(t,i,a,b);
  DBGFUN(<-);
}

void
FUN(setsm) (T *t, ssz_t n, const idx_t m[n], NUM a, NUM b)
{
  assert(t && m); DBGFUN(->);
  const idx_t i = mad_desc_idxsm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  FUN(seti)(t,i,a,b);
  DBGFUN(<-);
}

/* setv cases
   0   1     lo=2      hi=3        mo=4
  [.|.....|........|..........|............]
    |i...n|        |          |               lo=1, hi=1
    |0i.n0|        |          |               lo=1, hi=1
    |0i..n|        |          |               lo=1, hi=1
    |0i...|....n   |          |               lo=1, hi=2
    |0i...|........|..........|.......n0000   lo=1, hi=4
    |     |  i.....|......n   |               lo=2, hi=3
    |     |        |   i......|.......n0000   lo=3, hi=4
    |     |        |          |i......n0000   lo=4, hi=4
    |     |        |          |00i....n0000   lo=4, hi=4
    |     |        |          |i..........n   lo=4, hi=4
*/

void
FUN(setv) (T *t, idx_t i, ssz_t n, const NUM v[n])
{
  if (n <= 0) return;
  assert(t && v); DBGFUN(->);
  ssz_t nn = i+n;
  ensure(0 <= i && nn <= t->d->nc, "indexes %d:%d out of bounds", i, nn);

  const ord_t *ord = t->d->ords;
  const idx_t *o2i = t->d->ord2idx;
  ord_t lo = MAX(ord[i],1);
  ord_t hi = MIN(ord[nn-1], t->mo);

  ssz_t n0 = t->lo > lo ? o2i[lo  ] : i;
  ssz_t ni = MIN(o2i[hi+1], nn);
  ssz_t nj = t->hi < hi ? o2i[hi+1] : MAX(o2i[lo], nn);

//printf("setv: i=%d, n=%d, lo=%d, hi=%d, n0=%d, ni=%d, nj=%d, nn=%d %c\n",
//              i   , n   , lo   , hi   , n0   , ni   , nj   , nn,
//              ni == i+n ? ' ' : '*');

  FOR(j,n0, i) t->coef[j] = 0;
  FOR(j, i,ni) t->coef[j] = v[j-i];
  FOR(j,ni,nj) t->coef[j] = 0;

  if (!i) t->coef[0] = v[0];

  if (t->lo > lo) t->lo = lo;
  if (t->hi < hi) t->hi = hi;

  FUN(update)(t); // v may contain zeros...
  DBGFUN(<-); return;
}

// --- copiers ----------------------------------------------------------------o

void
FUN(cpyi) (const T *t, T *r, idx_t i)
{
  assert(t && r);
  if (!i) { FUN(setval)(r, t->coef[0]); return; }

  DBGFUN(->);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p",t->d,r->d);
  ensure(0 <= i && i < t->d->nc, "index %d out of bounds", i);
  NUM v = geti(t,i); FUN(reset0)(r); if (v) FUN(seti)(r,i,0,v);
  DBGFUN(<-);
}

void
FUN(cpys) (const T *t, T *r, ssz_t n, str_t s)
{
  assert(t && r && s); DBGFUN(->);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p",t->d,r->d);
  idx_t i = mad_desc_idxs(t->d,n,s);
  ensure(i >= 0, "invalid monomial");
  NUM v = geti(t,i); FUN(reset0)(r); if (v) FUN(seti)(r,i,0,v);
  DBGFUN(<-);
}

void
FUN(cpym) (const T *t, T *r, ssz_t n, const ord_t m[n])
{
  assert(t && r && m); DBGFUN(->);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p",t->d,r->d);
  idx_t i = mad_desc_idxm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  NUM v = geti(t,i); FUN(reset0)(r); if (v) FUN(seti)(r,i,0,v);
  DBGFUN(<-);
}

void
FUN(cpysm) (const T *t, T *r, ssz_t n, const idx_t m[n])
{
  assert(t && m); DBGFUN(->);
  ensure(t->d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p",t->d,r->d);
  idx_t i = mad_desc_idxsm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  NUM v = geti(t,i); FUN(reset0)(r); if (v) FUN(seti)(r,i,0,v);
  DBGFUN(<-);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(geti_r) (const T *t, idx_t i, NUM *r)
{ assert(r); *r = FUN(geti)(t, i); }

void FUN(gets_r) (const T *t, ssz_t n, str_t s, NUM *r)
{ assert(r); *r = FUN(gets)(t, n, s); }

void FUN(getm_r) (const T *t, ssz_t n, const ord_t m[n], NUM *r)
{ assert(r); *r = FUN(getm)(t, n, m); }

void FUN(getsm_r) (const T *t, ssz_t n, const idx_t m[n], NUM *r)
{ assert(r); *r = FUN(getsm)(t, n, m); }

void FUN(seti_r) (T *t, idx_t i, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(seti)(t, i, CPX(a), CPX(b)); }

void FUN(sets_r) (T *t, ssz_t n, str_t s, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(sets)(t, n, s, CPX(a), CPX(b)); }

void FUN(setm_r) (T *t, ssz_t n, const ord_t m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(setm)(t, n, m, CPX(a), CPX(b)); }

void FUN(setsm_r) (T *t, ssz_t n, const idx_t m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(setsm)(t, n, m, CPX(a), CPX(b)); }

void FUN(setvar_r) (T *t, num_t v_re, num_t v_im, idx_t iv, num_t scl_re, num_t scl_im)
{ FUN(setvar)(t, CPX(v), iv, CPX(scl)); }

void FUN(setval_r) (T *t, num_t v_re, num_t v_im)
{ FUN(setval)(t, CPX(v)); }

#endif // MAD_CTPSA_IMPL

// --- end --------------------------------------------------------------------o
