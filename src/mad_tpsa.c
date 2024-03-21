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

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

/* --- definition --------------------------------------------------------------
  GTPSA are *defined* in [0] U [lo,hi] with 0 <= hi <= mo <= d->mo <= 63
    - new/clear/reset GTPSA have lo=1, hi=0, nz=0, coef[0]=0 (see reset0).
    - coef[0] is always defined and lo >= 1 even if coef[0] != 0.
    - nz[o] == 0 for o in [0,lo) U (hi,mo], i.e. nz[0] == 0 because lo >= 1.
    - in [lo,hi]: nz[o] == 0 <=> all coef[[o]] == 0 by definition (not by value)
                  nz[o] == 1 <=> one coef[[o]] != 0 by value      (at least one)
------------------------------------------------------------------------------*/

// --- debugging --------------------------------------------------------------o

static inline log_t
FUN(check) (const T *t, ord_t *o_, idx_t *i_)
{
  ord_t _o = 0;
  idx_t _i = -1;

  if (!t->d || t->mo > t->d->mo || t->hi > t->mo ||
      (t->lo > t->hi && t->lo != 1)) goto ret;

  FOR(o,t->lo)
    if (mad_bit_tst(t->nz,o)) { _o = o; goto ret; }

  FOR(o,t->hi+1,t->mo+1)
    if (mad_bit_tst(t->nz,o)) { _o = o; goto ret; }

  TPSA_SCAN(t) {
    if (t->coef[i]) continue; // ensure one non-zero
    if (i+1 == o2i[o+1]) { _o = o, _i = i; goto ret; }
  }
  return TRUE;

ret:
  if (o_) *o_ = _o;
  if (i_) *i_ = _i;
  return FALSE;
}

log_t
FUN(isvalid) (const T *t)
{
  assert(t); DBGFUN(->);
  log_t ret = FUN(check)(t,0,0);
  DBGFUN(<-);
  return ret;
}

void
FUN(debug) (const T *t, str_t name_, str_t fname_, int line_, FILE *stream_)
{
  assert(t);

  ord_t o; idx_t i;
  if (FUN(check)(t,&o,&i)) return;

  const D* d = t->d;
  if (!stream_) stream_ = stdout;

  fprintf(stream_, "%s:%d: '%s' { lo=%d hi=%d mo=%d uid=%d, did=%d",
          fname_ ? fname_ : "??", line_, name_ ? name_ : "?",
          t->lo, t->hi, t->mo, t->uid, d ? d->id : -1);
  fflush(stream_);

  if (!d) { fprintf(stream_," }\n"); fflush(stream_); assert(d); }

  char bnz[DESC_MAX_ORD+2] = {0};
  for (ord_t b = 0; b <= t->mo; ++b)
    bnz[b] = '0' + mad_bit_tst(t->nz,b);
  fprintf(stream_," nz=%s ** bug @ o=%d i=%d }\n", bnz, o, i); fflush(stream_);

  const idx_t *o2i = d->ord2idx;
  idx_t ni = o2i[MIN(t->mo,d->to)+1];
  FOR(i,ni) fprintf(stream_," [%d:%d]=" FMT "\n", i,d->ords[i],VAL(t->coef[i]));
  fprintf(stream_,"\n"); fflush(stream_);
}

// --- introspection ----------------------------------------------------------o

const D*
FUN(desc) (const T *t)
{
  assert(t);
  return t->d;
}

int32_t
FUN(uid) (T *t, int32_t uid_) // set uid if uid != 0
{
  assert(t);
  const int32_t ret = t->uid;
  if (uid_) t->uid = uid_;
  return ret;
}

void
FUN(setnam) (T *t, str_t nam)
{
  assert(t); DBGFUN(->);
  strncpy(t->nam, nam, NAMSZ-1), t->nam[NAMSZ-1] = 0;
  DBGFUN(<-);
}

str_t
FUN(nam) (const T *t)
{
  assert(t);
  return t->nam;
}

ssz_t
FUN(len) (const T *t)
{
  assert(t); DBGFUN(->);
  const ssz_t ret = t->d->ord2idx[t->mo+1];
  DBGFUN(<-); return ret;
}

ord_t
FUN(ord) (const T *t)
{
  assert(t);
  return t->mo;
}

ord_t
(FUN(ordv)) (const T *t, ...)
{
  assert(t); DBGFUN(->);
  ord_t mo = t->mo;
  va_list args;
  va_start(args,t);
  while ((t = va_arg(args,T*)) != NULL) if (t->mo > mo) mo = t->mo;
  va_end(args);
  DBGFUN(<-); return mo;
}

ord_t
FUN(ordn) (ssz_t n, const T *t[n], log_t hi)
{
  assert(t); DBGFUN(->);
  ord_t mx = 0;
  if (hi) { FOR(i,n) if (t[i]->hi > mx) mx = t[i]->hi; }
  else    { FOR(i,n) if (t[i]->mo > mx) mx = t[i]->mo; }
  DBGFUN(<-); return mx;
}

log_t
FUN(isnul) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const log_t ret = FUN(isnul0)(t);
  DBGFUN(<-); return ret;
}

// --- init (unsafe) ----------------------------------------------------------o

T*
FUN(init) (T *t, const D *d, ord_t mo)
{
  assert(t && d && mo <= d->mo); DBGFUN(->);

  t->d = d, t->uid = 0, t->mo = mo, t->nam[0] = 0;
  FUN(reset0)(t);

  DBGFUN(<-); return t;
}

// --- ctors, dtor ------------------------------------------------------------o

T*
FUN(newd) (const D *d, ord_t mo)
{
  assert(d); DBGFUN(->);

  mo = MIN(mo, d->mo);
  ssz_t nc = d->ord2idx[mo+1]; // i.e. mad_desc_maxlen(d, mo);
  T *ret = mad_malloc(sizeof(T) + nc * sizeof(NUM));
  FUN(init)(ret, d, mo);

  DBGFUN(<-); return ret;
}

T*
FUN(new) (const T *t, ord_t mo)
{
  assert(t); DBGFUN(->);
  const D *d = t->d;

  mo = mo == mad_tpsa_same ? t->mo : MIN(mo, d->mo);
  ssz_t nc = d->ord2idx[mo+1]; // i.e. mad_desc_maxlen(d, mo);
  T *ret = mad_malloc(sizeof(T) + nc * sizeof(NUM));
  FUN(init)(ret, d, mo);

  DBGFUN(<-); return ret;
}

void
FUN(del) (const T *t)
{
  DBGFUN(->);
  if (t) mad_free((void*)t);
  DBGFUN(<-);
}

// --- clear, setvar, setprm, setval, update ----------------------------------o

void
FUN(clear) (T *t)
{
  assert(t); DBGFUN(->);
  FUN(reset0)(t);
  DBGFUN(<-);
}

void
FUN(setvar) (T *t, NUM v, idx_t iv, NUM scl)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;

  ensure(t->mo >= 1, "variables must have an order >= 1, got %d", t->mo);
  ensure(0 < iv && iv <= d->nv,
         "index 1<= %d <=%d is not a GTPSA variable", iv, d->nv);

  // set zero order
  t->coef[0] = v;

  // clear first order
  if (mad_bit_tst(t->nz,1)) FUN(clear0)(t,1);

  // set lo, hi, nz, coef[iv]
  t->lo = t->hi = 1, t->nz = 2, t->coef[iv] = scl ? scl : 1;

  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setprm) (T *t, NUM v, idx_t ip)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;

  ensure(t->mo == 1 , "parameters must be a GPTSA of order 1, got %d", t->mo);
  ensure(0 < ip && ip <= d->np,
         "index 1<= %d <=%d is not a GPTSA parameter", ip, d->np);

  // set zero order
  t->coef[0] = v;

  // clear first order
  if (mad_bit_tst(t->nz,1)) FUN(clear0)(t,1);

  // set lo, hi, nz, coef[ip]
  t->lo = t->hi = 1, t->nz = 2, t->coef[ip+d->nv] = 1;

  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setval) (T *t, NUM v)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  FUN(reset0)(t); t->coef[0] = v;
  DBGTPSA(t); DBGFUN(<-);
}

log_t
FUN(update) (T *t, num_t eps)
{
  assert(t); DBGFUN(->);
  const bit_t nz = t->nz;
  if (eps <= 0) { TPSA_SCAN_Z(t,t->lo,t->hi) FUN(stabilize0)(t,o,eps); }
  else          { TPSA_SCAN_Z(t,t->lo,t->hi) FUN(update0   )(t,o    ); }
  const log_t up = t->nz != nz;
  if (up) FUN(adjust0)(t);
  DBGTPSA(t); DBGFUN(<-); return up;
}

// --- copy, update, convert, swap --------------------------------------------o

void
FUN(copy) (const T *t, T *r)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t);

  if (t != r) {
    const D *d = t->d;
    ensure(d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", d, r->d);

    // copy lo, hi, nz
    FUN(copy0)(t, r);

    // copy coef[0]
    r->coef[0] = t->coef[0];

    // copy coefs
    TPSA_SCAN(r) r->coef[i] = t->coef[i];
  } else
    FUN(trunc0)(r);

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(sclord) (const T *t, T *r, log_t inv, log_t prm)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t);

  FUN(copy)(t,r);

  const int    np = !prm;
  const ord_t *po = r->d->prms, lo = MAX(r->lo,2);

  if (inv) { // scale coefs by 1/order
    TPSA_SCAN(r,lo,r->hi) r->coef[i] /= o - po[i] * np;
  } else {   // scale coefs by order
    TPSA_SCAN(r,lo,r->hi) r->coef[i] *= o - po[i] * np;
  }

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(getord) (const T *t, T *r, ord_t ord)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t); DBGTPSA(r);
  const D *d = t->d;
  ensure(d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", d, r->d);

  if (!mad_bit_tst(t->nz,ord) || ord > MIN(r->mo,d->to)) {
    FUN(setval)(r, ord ? 0 : t->coef[0]);
    DBGFUN(<-); return;
  }

  // set coef[0], lo, hi, nz
  r->coef[0] = 0, r->lo = r->hi = ord, r->nz = mad_bit_set(0, ord);

  // copy data
  if (t != r) { TPSA_SCAN_O(r,ord) r->coef[i] = t->coef[i]; }

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(cutord) (const T *t, T *r, int ord)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t); DBGTPSA(r);
  const D *d = t->d;
  ensure(d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", d, r->d);

  if (ord <= 0) { // cut 0..|ord|, see copy0 with t->lo = |ord|+1
    r->hi = MIN(t->hi, r->mo, d->to);
    r->lo = MIN(1-ord, r->hi);
    r->nz = mad_bit_hcut(mad_bit_lcut(t->nz, r->lo), r->hi);
    r->coef[0] = 0;
  } else {        // cut |ord|..mo, see copy0 with t->hi = |ord|-1
    r->lo = t->lo;
    r->hi = MIN(ord-1, r->mo, d->to);
    r->nz = mad_bit_hcut(t->nz, r->hi);
    r->coef[0] = t->coef[0];
  }

  if (!r->nz) FUN(setval)(r, t->coef[0]); else
  if (r != t) { TPSA_SCAN(r) r->coef[i] = t->coef[i]; }

  DBGTPSA(r); DBGFUN(<-);
}

idx_t
FUN(maxord) (const T *t, ssz_t n, idx_t idx_[n])
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  if (idx_) FOR(i,n) idx_[i] = -1;

  num_t mv =  0; // max of all values
  idx_t mi = -1; // idx of max for all
  ord_t hi = MIN(n-1,t->hi,t->d->to);
  TPSA_SCAN_Z(t,t->lo,hi) {
    num_t mo = 0; // max of this order
    TPSA_SCAN_O(t,o)
      if (mo < fabs(t->coef[i])) {
        mo = fabs(t->coef[i]);        // save max for this order
        if (idx_) idx_[o] = i;        // save idx for this order
        if (mv < mo) mv = mo, mi = i; // save max and idx for all orders
      }
  }
  DBGFUN(<-); return mi;
}

void
FUN(convert) (const T *t, T *r_, ssz_t n, idx_t t2r_[n], int pb)
{
  assert(t && r_); DBGFUN(->); DBGTPSA(t); DBGTPSA(r_);
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
  ord_t t_hi = MIN(t->hi, r->mo, rd->to);         // see copy0
  r->coef[0] = t->coef[0], r->nz = 0;
  TPSA_SCAN(t,t_hi) {
    if (!t->coef[i]) continue;
    mad_desc_mono(td, i, tn, tm, NULL);           // get tm mono at index i
    mad_mono_fill(rn, rm, 0);
    int sgn = 0;
    FOR(j,tn) {                                   // set rm mono
      if (t2r[j] < 0 && tm[j]) goto skip;         // discard coef
      rm[t2r[j]] = tm[j];                         // translate tm to rm
      sgn = sgn - pbs[j] * (tm[j] & 1);           // poisson bracket
    }
    idx_t ri = mad_desc_idxm(rd, rn, rm);         // get index ri of mono rm
#if DEBUG > 2
    printf("cvt %d -> %d %c : ", i+1, ri+1, i==ri?' ' : SIGN1(sgn%2)<0?'-':'+');
    mad_mono_print(tn, tm, 0); printf(" -> "); mad_mono_print(rn, rm, 0);
    printf(" : %-.16e\n", t->coef[i]); // works only for real, warn for complex
#endif
    if (ri >= 0) {
      r->nz = mad_bit_set(r->nz, ords[ri]);
      r->coef[ri] = SIGN1(sgn%2)*t->coef[i];
    }
  skip: ;
  }
  FUN(adjust0)(r);

  if (r != r_) { FUN(copy)(r,r_); REL_TMPX(r); }

  DBGTPSA(r_); DBGFUN(<-);
}

// --- indexing / monomials ---------------------------------------------------o

ord_t
FUN(mono) (const T *t, idx_t i, ssz_t n, ord_t m_[n], ord_t *p_)
{
  assert(t); DBGTPSA(t);
  const ord_t ret = mad_desc_mono(t->d, i, n, m_, p_);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxs) (const T *t, ssz_t n, str_t s)
{
  assert(t && s); DBGTPSA(t);
  const idx_t ret = mad_desc_idxs(t->d, n, s);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGTPSA(t);
  const idx_t ret = mad_desc_idxm(t->d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxsm) (const T *t, ssz_t n, const idx_t m[n])
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  const idx_t ret = mad_desc_idxsm(t->d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
FUN(cycle) (const T *t, idx_t i, ssz_t n, ord_t m_[n], NUM *v_)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  const ord_t mo = MIN(t->mo, d->to);
  if (++i >= d->ord2idx[mo+1] || i < 0) { DBGFUN(<-); return -1; }

  if (i || !t->coef[0]) {
    ord_t lo = MAX(t->lo, d->ords[i]);
    idx_t hi = MIN(t->hi, d->to);
    TPSA_SCAN_Z(t,lo,hi) {
      for (i = MAX(i,o2i[o]); i < o2i[o+1]; i++)
        if (t->coef[i]) goto ret;
    }
    DBGFUN(<-); return -1;
  }

ret:
  if (v_) *v_ = t->coef[i];
  if (m_) mad_mono_copy(MIN(n,d->nn), d->To[i], m_);
  DBGFUN(<-); return i;
}

// --- getters ----------------------------------------------------------------o

static inline NUM
geti (const T *t, idx_t i)
{
  const D *d = t->d;
  const ord_t o = d->ords[i] > d->to ? 0 : d->ords[i];
  return !mad_bit_tst(t->nz, o) && i ? 0 : t->coef[i];
}

NUM
FUN(get0) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const NUM ret = t->coef[0];
  DBGFUN(<-); return ret;
}

NUM
FUN(geti) (const T *t, idx_t i)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  ensure(0 <= i && i < t->d->nc, "index %d out of bounds", i);
  const NUM ret = geti(t,i);
  DBGFUN(<-); return ret;
}

NUM
FUN(gets) (const T *t, ssz_t n, str_t s)
{
  // --- mono is a string; represented as "[0-9A-Z]*"
  assert(t && s); DBGFUN(->); DBGTPSA(t);
  const idx_t i = mad_desc_idxs(t->d,n,s);
  ensure(i >= 0, "invalid monomial");
  const NUM ret = geti(t,i);
  DBGFUN(<-); return ret;
}

NUM
FUN(getm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  const idx_t i = mad_desc_idxm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  const NUM ret = geti(t,i);
  DBGFUN(<-); return ret;
}

NUM
FUN(getsm) (const T *t, ssz_t n, const idx_t m[n])
{
  // --- mono is sparse; represented as [(i,o)]
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  const idx_t i = mad_desc_idxsm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  const NUM ret = geti(t,i);
  DBGFUN(<-); return ret;
}

void
FUN(getv) (const T *t, idx_t i, ssz_t n, NUM v[n])
{
  assert(t && v); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  ensure(0 <= i && i+n <= d->nc, "indexes %d:%d out of bounds", i, i+n);

  const ord_t *ord = d->ords+i;
  const NUM  *coef = t->coef+i;
  const ord_t   hi = MIN(t->hi, d->to);
  const idx_t   nj = MIN(d->ord2idx[hi+1], i+n)-i;
  FOR(j,nj)   v[j] = mad_bit_tst(t->nz, ord[j]) ? coef[j] : 0;
  FOR(j,nj,n) v[j] = 0; if (!i && n) v[0] = coef[0];
  DBGFUN(<-);
}

// --- setters ----------------------------------------------------------------o

void
FUN(set0) (T *t, NUM a, NUM b)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  t->coef[0] = a*t->coef[0] + b;
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(seti) (T *t, idx_t i, NUM a, NUM b)
{
  assert(t); DBGFUN(->);
  if (!i) { FUN(set0)(t,a,b); DBGFUN(<-); return; }

  DBGTPSA(t);
  const D *d = t->d;
  ensure(0 < i && i < d->nc, "index order exceeds GPTSA maximum order");

  // discard value
  const ord_t o = d->ords[i];
  if (o > MIN(t->mo, d->to)) { DBGFUN(<-); return; }

  NUM v = mad_bit_tst(t->nz,o) ? a*t->coef[i]+b : b;

  if (!v && mad_bit_tst(t->nz,o)) {
    t->coef[i] = 0, FUN(update0)(t,o);
    if (!mad_bit_tst(t->nz,o)) FUN(adjust0)(t);
  } else
  if (v && !mad_bit_tst(t->nz,o)) {
    FUN(clear0)(t,o)->coef[i] = v;
    t->nz = mad_bit_set(t->nz,o), FUN(adjust0)(t);
  } else
    t->coef[i] = v;

  DBGTPSA(t); DBGFUN(<-);
}

/*
       0   1      2           3                4                   mo=5
      [.|.....|.....[..|............|...)...............|......................]
             vlo=2  i          v        i+n            vhi=4
lo=hi=1 [.....]00000)                   [000000000000000]
lo=hi=2       [........]                [000000000000000]
lo=hi=3       [00000)  [............]   [000000000000000]
lo=hi=4       [00000)               [...................]
lo=hi=5       [00000)                   [000000000000000[......................]
lo=2,hi=3     [........|............]   [000000000000000]
lo=3,hi=4     [00000)  [............|...................]
*/

void
FUN(setv) (T *t, idx_t i, ssz_t n, const NUM v[n])
{
  assert(t && v); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  ensure(0 <= i && i+n <= d->nc, "indexes %d:%d out of bounds", i, i+n);

  // compute boundaries
  const idx_t *o2i = d->ord2idx;
  const ord_t *ord = d->ords+i;
  ord_t vlo = MIN(t->hi+1, ord[ 0 ]);
  ord_t vhi = MAX(t->lo-1, ord[n-1]);
        vhi = MIN(vhi, t->mo, d->to);
  idx_t j = o2i[vlo], nj = MIN(i+n, o2i[vhi+1]);

  for (; j < i         ; j++) t->coef[j] = 0;      // right gap (if any)
  for (; j < nj        ; j++) t->coef[j] = v[j-i]; // copy data
  for (; j < o2i[vhi+1]; j++) t->coef[j] = 0;      // left gap  (if any)

  // enlarge [lo,hi] if needed
  if (t->lo > vlo) t->lo = MAX(1,vlo);
  if (t->hi < vhi) t->hi = vhi;

  // activate nz bits that might be affected by vector content
  t->nz |= mad_bit_mask(~0ull, t->lo, t->hi);

  FUN(update)(t,0);
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(sets) (T *t, ssz_t n, str_t s, NUM a, NUM b)
{
  assert(t && s); DBGFUN(->); DBGTPSA(t);
  idx_t i = mad_desc_idxs(t->d,n,s);
  ensure(i >= 0, "invalid monomial");
  FUN(seti)(t,i,a,b);
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setm) (T *t, ssz_t n, const ord_t m[n], NUM a, NUM b)
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  idx_t i = mad_desc_idxm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  FUN(seti)(t,i,a,b);
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setsm) (T *t, ssz_t n, const idx_t m[n], NUM a, NUM b)
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  idx_t i = mad_desc_idxsm(t->d,n,m);
  ensure(i >= 0, "invalid monomial");
  FUN(seti)(t,i,a,b);
  DBGTPSA(t); DBGFUN(<-);
}

// --- without complex-by-value version ---------------------------------------o

#ifdef MAD_CTPSA_IMPL

void FUN(get0_r) (const T *t, NUM *r)
{ assert(r); *r = FUN(get0)(t); }

void FUN(geti_r) (const T *t, idx_t i, NUM *r)
{ assert(r); *r = FUN(geti)(t, i); }

void FUN(gets_r) (const T *t, ssz_t n, str_t s, NUM *r)
{ assert(r); *r = FUN(gets)(t, n, s); }

void FUN(getm_r) (const T *t, ssz_t n, const ord_t m[n], NUM *r)
{ assert(r); *r = FUN(getm)(t, n, m); }

void FUN(getsm_r) (const T *t, ssz_t n, const idx_t m[n], NUM *r)
{ assert(r); *r = FUN(getsm)(t, n, m); }

void FUN(set0_r) (T *t, num_t a_re, num_t a_im, num_t b_re, num_t b_im)
{ FUN(set0)(t, CPX(a), CPX(b)); }

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
