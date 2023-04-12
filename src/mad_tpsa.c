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
#include <assert.h>
#include <limits.h>
#include <math.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

// --- debugging --------------------------------------------------------------o

/*
GTPSA are defined in [0] U [lo,hi], and nz[o] == 0 in [0,lo) U (hi,mo]
GTPSA just initialized have lo=0, hi=0, nz=0, coef[0]=0 (see reset0)
in [lo,hi]: nz[o] == 0 <=> all coef[[o]] == 0
            nz[o] == 1 <=> any coef[[o]] != 0 (or none: conservative)
            nz[0] == 0 <=>     coef[ 0 ] == 0 && lo >= 1
            nz[0] == 1 <=>     coef[ 0 ] != 0 && lo == 0

Improvement: remove constraint t->lo > 0 => t->coef[0] == 0
  t->coef[0] could be handled separately and let t->lo > 0 to start higher.
*/

static inline log_t
FUN(check) (const T *t, ord_t *o_, idx_t *i_)
{
  ord_t o  = 0;
  idx_t i  = -1;
  log_t ok = t->d != NULL;

  if (!ok) goto ret; // i == -1

  ok &= t->mo <= t->d->mo;
  ok &= t->hi <= t->mo;
  ok &= t->lo <= t->hi;
  ok &= t->lo ? t->coef[0] == 0 : TRUE;

  if (!ok) goto ret; // i == -1

  for (; o < t->lo; ++o) {
    ok &= !mad_bit_tst(t->nz,o);
    if (!ok) goto ret;
  }

  const idx_t *o2i = t->d->ord2idx;

  for (; o <= t->hi; ++o) {
    if (!mad_bit_tst(t->nz,o))
      for (i = o2i[o]; i < o2i[o+1]; ++i) {
        ok &= t->coef[i] == 0; // ensure all zero
        if (!ok) goto ret;
      }
    else if (TPSA_STRICT_NZ > 2) { // max level required for validation...
      for (i = o2i[o]; i < o2i[o+1] && t->coef[i] == 0; ++i) ;
      ok &= i < o2i[o+1]; // ensure any non-zero
      if (!ok) goto ret;
    }
  }

  for (; o <= t->mo; ++o) {
    ok &= !mad_bit_tst(t->nz,o);
    if (!ok) goto ret;
  }

ret:
  if (o_) *o_ = o;
  if (i_) *i_ = i;
  return ok;
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

  fprintf(stream_, "%s:%d: '%s' { lo=%d hi=%d mo=%d id=%d",
          fname_ ? fname_ : "??", line_, name_ ? name_ : "?",
          t->lo, t->hi, t->mo, d ? d->id : -1);
  fflush(stream_);

  if (!d) { fprintf(stream_," }\n"); fflush(stream_); assert(d); }

  char bnz[DESC_MAX_ORD+2] = {0};
  for (ord_t b = 0; b <= t->mo; ++b)
    bnz[b] = '0' + mad_bit_tst(t->nz,b);
  fprintf(stream_," nz=%s ** o=%d i=%d }", bnz, o, i); fflush(stream_);

  const idx_t *o2i = d->ord2idx;
  idx_t ni = DEBUG > 1 ? o2i[t->mo+1] : MIN(25,o2i[t->mo+1]);
  for (idx_t i = 0; i < ni; ++i)
    fprintf(stream_," [%d]=" FMT, i, VAL(t->coef[i]));
  fprintf(stream_,"\n"); fflush(stream_);
}

// --- introspection ----------------------------------------------------------o

const D*
FUN(desc) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D* ret = t->d;
  DBGFUN(<-); return ret;
}

int32_t
FUN(uid) (T *t, int32_t uid_) // set uid if != 0
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  int32_t ret = t->uid;
  if (uid_) t->uid = uid_;
  DBGFUN(<-); return ret;
}

str_t
FUN(nam) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  str_t ret = t->nam;
  DBGFUN(<-); return ret;
}

ssz_t
FUN(len) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  ssz_t ret = t->d->ord2idx[t->mo+1];
  DBGFUN(<-); return ret;
}

ord_t
FUN(ord) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  ord_t ret = t->mo;
  DBGFUN(<-); return ret;
}

ord_t
(FUN(ordv)) (const T *t, ...)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  ord_t mo = t->mo;
  va_list args;
  va_start(args,t);
  while ((t = va_arg(args,T*)) != NULL) {
    DBGTPSA(t);
    if (t->mo > mo) mo = t->mo;
  }
  va_end(args);
  DBGFUN(<-); return mo;
}

ord_t
FUN(ordn) (ssz_t n, const T *t[n])
{
  assert(t); DBGFUN(->);
  ord_t mo = 0;
  for (idx_t i = 0; i < n; i++){
    DBGTPSA(t[i]);
    if (t[i]->mo > mo) mo = t[i]->mo;
  }
  DBGFUN(<-); return mo;
}

// --- init (unsafe) ----------------------------------------------------------o

T*
FUN(init) (T *t, const D *d, ord_t mo)
{
  assert(t); DBGFUN(->);

  if (!d) d = mad_desc_curr;
  ensure(d, "GTPSA descriptor not found (no current one?)");

  if (mo == mad_tpsa_default) mo = d->mo;
  else ensure(mo <= d->mo, "GTPSA order exceeds descriptor maximum order");

  t->d = d, t->uid = 0, t->mo = mo, t->nam[0] = 0;
  FUN(reset0)(t);

  DBGFUN(<-); return t;
}

// --- ctors, dtor ------------------------------------------------------------o

T*
FUN(newd) (const D *d, ord_t mo)
{
  DBGFUN(->);
  if (!d) d = mad_desc_curr;
  ensure(d, "GTPSA descriptor not found (no current one?)");

  if (mo == mad_tpsa_default) mo = d->mo;
  else ensure(mo <= d->mo, "GTPSA order exceeds descriptor maximum order");

  ssz_t nc = mad_desc_maxlen(d, mo);
  T *t = mad_malloc(sizeof(T) + nc * sizeof(NUM));
  t->d = d, t->uid = 0, t->mo = mo, t->nam[0] = 0;
  FUN(reset0)(t);

  DBGFUN(<-); return t;
}

T*
FUN(new) (const T *t, ord_t mo)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  T *ret = FUN(newd)(t->d, mo == mad_tpsa_same ? t->mo : mo);
  DBGFUN(<-); return ret;
}

void
FUN(del) (const T *t)
{
  DBGFUN(->);
  if (t) mad_free((void*)t);
  DBGFUN(<-);
}

// --- clear, scalar ----------------------------------------------------------o

void
FUN(clear) (T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  FUN(reset0)(t);
  DBGFUN(<-);
}

log_t
FUN(isnul) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  const idx_t *o2i = t->d->ord2idx;
  for (ord_t o = t->lo; o <= t->hi; o++) {
    if (mad_bit_tst(t->nz, o))
      for (idx_t i = o2i[o]; i < o2i[o+1]; ++i)
        if (t->coef[i]) { DBGFUN(<-); return FALSE; }
  }

  DBGFUN(<-);
  return TRUE;
}

void
FUN(setnam) (T *t, str_t nam)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  strncpy(t->nam, nam, NAMSZ-1), t->nam[NAMSZ-1] = 0;
  DBGFUN(<-);
}

void
FUN(setvar) (T *t, NUM v, idx_t iv, NUM scl)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  const int nfo = !(iv && t->mo && d->to);

  // v=0 and no first order: reset
  if (nfo && !v) {
    FUN(reset0)(t); DBGFUN(<-); return;
  }

  t->coef[0] = v;

  // no first order: set lo, hi, nz
  if (nfo) {
    t->lo = t->hi = 0, t->nz = !!v;
    DBGTPSA(t); DBGFUN(<-); return;
  }

  ensure(0 <= iv && iv <= d->nv,
         "index %d exceeds GPTSA number of variables %d", iv, d->nv);

  // clear first order
  const idx_t *o2i = d->ord2idx;
  for (idx_t i = o2i[1]; i < o2i[2]; ++i) t->coef[i] = 0;

  // set lo, hi, nz, coef[iv]
  t->hi = 1, t->lo = !v, t->nz = v ? 3 : 2;
  t->coef[iv] = scl ? scl : 1;

  DBGTPSA(t); DBGFUN(<-);
}

// --- copy, convert, swap ----------------------------------------------------o

void
FUN(copy) (const T *t, T *r)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t);
  if (t != r) {
    const D *d = t->d;
    ensure(d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", d, r->d);

    // copy lo, hi(mo,to), nz(hi)
    FUN(copy0)(t, r);

    // copy name if unnamed
    if (!r->nam[0]) strcpy(r->nam, t->nam);

    // copy coefs
    const idx_t *o2i = d->ord2idx;
    for (idx_t i = o2i[r->lo]; i < o2i[r->hi+1]; ++i) r->coef[i] = t->coef[i];
  }
  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(sclord) (const T *t, T *r, log_t inv)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t);

  FUN(copy)(t,r);

  // scale coefs
  const idx_t *o2i = r->d->ord2idx;
  if (inv) {
    for (ord_t o = MAX(r->lo,2); o <= r->hi; ++o)
      if (mad_bit_tst(r->nz,o))
        for (idx_t i = o2i[o]; i < o2i[o+1]; ++i)
          r->coef[i] /= o; // scale coefs by 1/order
  } else {
    for (ord_t o = MAX(r->lo,2); o <= r->hi; ++o)
      if (mad_bit_tst(r->nz,o))
        for (idx_t i = o2i[o]; i < o2i[o+1]; ++i)
          r->coef[i] *= o; // scale coefs by order
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
    FUN(reset0)(r); DBGFUN(<-); return;
  }

  // set lo, hi, nz, coef[0]
  r->lo = r->hi = ord;
  r->nz = mad_bit_set(0, ord);
  if (r->lo) r->coef[0] = 0;

  // copy data
  if (t != r) {
    const idx_t *o2i = d->ord2idx;
    for (idx_t i = o2i[r->lo]; i < o2i[r->hi+1]; ++i) r->coef[i] = t->coef[i];
  }

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(cutord) (const T *t, T *r, int ord)
{
  assert(t && r); DBGFUN(->); DBGTPSA(t); DBGTPSA(r);
  const D *d = t->d;
  ensure(d == r->d, "incompatible GTPSAs descriptors 0x%p vs 0x%p", d, r->d);

  if (ord < 0) { // cut 0..|ord|, see copy0 with t->lo = |ord|+1
    r->hi = MIN(t->hi, r->mo, d->to);
    r->nz = mad_bit_hcut(mad_bit_lcut(t->nz, -ord+1), r->hi);
    if (!r->nz) { FUN(reset0)(r); DBGFUN(<-); return; }
    r->lo = -ord+1; r->coef[0] = 0;
  } else {      // cut |ord|..mo, see copy0 with t->hi = |ord|-1
    r->hi = MIN(ord-1, r->mo, d->to);
    r->nz = mad_bit_hcut(t->nz, r->hi);
    if (!r->nz) { FUN(reset0)(r); DBGFUN(<-); return; }
    if ((r->lo=t->lo)) r->coef[0] = 0;
  }

  // copy data
  if (r != t) {
    const idx_t *o2i = d->ord2idx;
    for (idx_t i = o2i[r->lo]; i < o2i[r->hi+1]; ++i) r->coef[i] = t->coef[i];
  }

  DBGTPSA(r); DBGFUN(<-);
}

idx_t
FUN(maxord) (const T *t, ssz_t n, idx_t idx_[n])
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  if (idx_) for (ord_t o=0; o < n; ++o) idx_[o] = -1;

  const idx_t *o2i = t->d->ord2idx;
  num_t mv =  0; // max of all values
  idx_t mi = -1; // idx of max for all
  for (ord_t o = t->lo; o < MIN(n,t->hi+1); ++o)
    if (mad_bit_tst(t->nz,o)) {
      num_t mo = 0; // max of this order
      for (idx_t i = o2i[o]; i < o2i[o+1]; ++i)
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

  // fast branch for (almost) compatible cases avoiding monomials translation
  if (!t2r_) {
    const D *td = t->d, *rd = r_->d;

    if (td == rd) { FUN(copy)(t,r_); DBGFUN(<-); return; }

    if (td->nn == rd->nn && td->np == rd->np && td->po == rd->po && !td->uno && !rd->uno) {
      FUN(copy0)(t, r_); // copy lo, hi(mo,to), nz(hi)
      const idx_t *o2i = rd->ord2idx, *to2i = td->ord2idx;
      if (o2i[r_->lo] == to2i[r_->lo] && o2i[r_->hi] == to2i[r_->hi]) {
        for (idx_t i = o2i[r_->lo]; i < o2i[r_->hi+1]; ++i) // copy coefs
          r_->coef[i] = t->coef[i];
        DBGTPSA(r_); DBGFUN(<-); return;
      }
    }
  }

  // slow branch for non-compatible cases with monomials translation
  T *r = t == r_ ? GET_TMPX(r_) : FUN(reset0)(r_);
  ssz_t rn = r->d->nv, tn = t->d->nv;
  ord_t rm[rn], tm[tn];
  idx_t t2r[tn]; // if t2r[i]>=0 then rm[t2r[i]] = tm[i] for i=0..tn-1
  int   pbs[tn];

  idx_t i = 0;
  if (!t2r_)
    for (; i < MIN(tn,rn); ++i) t2r[i] = i, pbs[i] = 0; // identity
  else
    for (; i < MIN(tn, n); ++i) {
      t2r[i] = t2r_[i] >= 0 && t2r_[i] < rn ? t2r_[i] : -1; // -1 discard var
      pbs[i] = pb*(t2r[i]-i)%2 < 0; // pb sign, ignored for discarded vars
    } // fromptc: pt => 1*(6-5)<0=0, t => 1*(5-6)<0=1, x,px,y,py => 1*(i-i)<0=0
  for (; i < tn; i++) t2r[i] = -1;  // discard remaining vars

  const idx_t *o2i = t->d->ord2idx;
  ord_t t_hi = MIN(t->hi, r->mo, t->d->to);
  for (idx_t ti = o2i[t->lo]; ti < o2i[t_hi+1]; ++ti) {
    if (t->coef[ti] == 0) goto skip;
    mad_desc_mono(t->d, ti, tn, tm);              // get tm mono at index ti
    mad_mono_fill(rn, rm, 0);
    int sgn = 0;
    for (idx_t i = 0; i < tn; ++i) {              // set rm mono
      if (t2r[i] < 0 && tm[i]) goto skip;         // discard coef
      rm[t2r[i]] = tm[i];                         // translate tm to rm
      sgn = sgn - pbs[i] * (tm[i] & 1);           // poisson bracket
    }
    idx_t ri = mad_desc_idxm(r->d, rn, rm);       // get index ri of mono rm
#if DEBUG > 2
    printf("cvt %d -> %d %c : ", ti+1, ri+1, ti==ri?' ' : SIGN1(sgn%2)<0?'-':'+');
    mad_mono_print(tn, tm, 0); printf(" -> "); mad_mono_print(rn, rm, 0);
    printf(" : %-.16e\n", t->coef[ti]); // works only for real, warn for complex
#endif
    if (ri >= 0) FUN(seti)(r, ri, 0, SIGN1(sgn%2)*t->coef[ti]);
  skip: ;
  }

  if (r != r_) { FUN(copy)(r,r_); REL_TMPX(r); }

  DBGTPSA(r_); DBGFUN(<-);
}

// --- indexing / monomials ---------------------------------------------------o

ord_t
FUN(mono) (const T *t, idx_t i, ssz_t n, ord_t m_[n])
{
  assert(t); DBGTPSA(t);
  ord_t ret = mad_desc_mono(t->d, i, n, m_);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxs) (const T *t, ssz_t n, str_t s)
{
  assert(t && s); DBGTPSA(t);
  idx_t ret = mad_desc_idxs(t->d, n, s);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGTPSA(t);
  idx_t ret = mad_desc_idxm(t->d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxsm) (const T *t, ssz_t n, const idx_t m[n])
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  idx_t ret = mad_desc_idxsm(t->d, n, m);
  DBGFUN(<-); return ret;
}

idx_t
FUN(cycle) (const T *t, idx_t i, ssz_t n, ord_t m_[n], NUM *v_)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  i += 1;
  ensure(0 <= i && i <= d->nc, "index %d out of bounds", i);

  const idx_t *o2i = d->ord2idx;
  idx_t ni = o2i[MIN(t->hi,d->to)+1];
  for (i = MAX(i,o2i[t->lo]); i < ni && t->coef[i] == 0; ++i) ;

  if (i >= ni) { DBGFUN(<-); return -1; }

  if (m_) {
    ensure(0 <= n && n <= d->nn, "invalid monomial length %d", n);
    mad_mono_copy(n, d->To[i], m_);
  }
  if (v_) *v_ = t->coef[i];

  DBGFUN(<-); return i;
}

// --- getters ----------------------------------------------------------------o

NUM
FUN(get0) (const T *t)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  NUM ret = t->coef[0];
  DBGFUN(<-); return ret;
}

NUM
FUN(geti) (const T *t, idx_t i)
{
  assert(t); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  ensure(0 <= i && i < d->nc, "index %d out of bounds", i);
  ord_t o = d->ords[i];
  NUM ret = t->lo <= o && o <= MIN(t->hi,d->to) ? t->coef[i] : 0;
  DBGFUN(<-); return ret;
}

NUM
FUN(gets) (const T *t, ssz_t n, str_t s)
{
  // --- mono is a string; represented as "[0-9A-Z]*"
  assert(t && s); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  idx_t i = mad_desc_idxs(d,n,s);
  ensure(i >= 0, "invalid monomial");
  ord_t o = d->ords[i];
  NUM ret = t->lo <= o && o <= MIN(t->hi,d->to) ? t->coef[i] : 0;
  DBGFUN(<-); return ret;
}

NUM
FUN(getm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  idx_t i = mad_desc_idxm(d,n,m);
  ensure(i >= 0, "invalid monomial");
  ord_t o = d->ords[i];
  NUM ret = t->lo <= o && o <= MIN(t->hi,d->to) ? t->coef[i] : 0;
  DBGFUN(<-); return ret;
}

NUM
FUN(getsm) (const T *t, ssz_t n, const idx_t m[n])
{
  // --- mono is sparse; represented as [(i,o)]
  assert(t && m); DBGFUN(->); DBGTPSA(t);
  const D *d = t->d;
  idx_t i = mad_desc_idxsm(d,n,m);
  ensure(i >= 0, "invalid monomial");
  ord_t o = d->ords[i];
  NUM ret = t->lo <= o && o <= MIN(t->hi,d->to) ? t->coef[i] : 0;
  DBGFUN(<-); return ret;
}

void
FUN(getv) (const T *t, idx_t i, ssz_t n, NUM v[n])
{
  assert(t && v); DBGFUN(->); DBGTPSA(t);

  const D *d = t->d;
  ensure(0 <= i && i+n <= d->nc, "indexes %d:%d out of bounds", i, i+n);

  ord_t hi = MIN(t->hi, d->to);
  const ord_t *ord = d->ords+i;
  const NUM  *coef = t->coef+i;
  for (idx_t j = 0; j < n; ++j)
    v[j] = t->lo <= ord[j] && ord[j] <= hi ? coef[j] : 0;

  DBGFUN(<-);
}

// --- setters ----------------------------------------------------------------o

void
FUN(set0) (T *t, NUM a, NUM b)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  NUM v = a*t->coef[0]+b;

  if (!v) { // clear coef[0], set lo and nz
    t->nz = mad_bit_clr(t->nz, 0);
    if (!t->nz) { FUN(reset0)(t); DBGFUN(<-); return; }
    t->lo = mad_bit_lowest(t->nz);
    t->coef[0] = 0;
  } else {  // set coef[0], clear (0,lo)
    const idx_t *o2i = t->d->ord2idx;
    for (idx_t c = o2i[1]; c < o2i[t->lo]; ++c) t->coef[c] = 0;
    t->nz = mad_bit_set(t->nz, 0);
    t->lo = 0;
    t->coef[0] = v;
  }

  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(seti) (T *t, idx_t i, NUM a, NUM b)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  if (!i) { FUN(set0)(t,a,b); DBGTPSA(t); DBGFUN(<-); return; }

  const D *d = t->d;
  ensure(0 <= i && i < d->nc, "index order exceeds GPTSA maximum order");

  ord_t o = d->ords[i];

  // discard value
  if (o > MIN(t->mo, d->to)) { DBGTPSA(t); DBGFUN(<-); return; }

  NUM v = mad_bit_tst(t->nz, o) ? a*t->coef[i]+b : b;

  if (!v) { // clear coef[i], might clear order
    t->coef[i] = 0;
    // update0 below is stricter but can be quadratic with consecutive seti(),
    // like in the pattern of clearing an order minus few ending monomials, so
    // I choosed to be lazy and fast with such consecutive seti() as there is
    // little chance that seti() clears the last non-zero coef of the hpoly.
    if (TPSA_STRICT_NZ > 2 && mad_bit_tst(t->nz, o)) FUN(update0)(t,o,o);
    DBGTPSA(t); DBGFUN(<-); return;
  }

  const idx_t *o2i = d->ord2idx;
  if (o < t->lo) { // extend left
    for (idx_t c = o2i[o]; c < o2i[t->lo]; ++c) t->coef[c] = 0;
    t->lo = o;
  } else
  if (o > t->hi) { // extend right
    for (idx_t c = o2i[t->hi+1]; c < o2i[o+1]; ++c) t->coef[c] = 0;
    t->hi = o;
  }
  t->nz = mad_bit_set(t->nz, o);
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
  ensure(0 <= i && i+n <= d->nc, "index order exceeds GPTSA maximum order");

  // compute boundaries
  const idx_t *o2i = d->ord2idx;
  const ord_t *ord = d->ords+i;
  ord_t vlo =         MIN(t->hi+1, ord[ 0 ]);
  ord_t vhi = t->lo ? MAX(t->lo-1, ord[n-1]) : ord[n-1];
        vhi = MIN(vhi, t->mo, d->to);
  idx_t j = o2i[vlo], nj = MIN(i+n, o2i[vhi+1]);

  for (; j < i         ; j++) t->coef[j] = 0;      // right gap (if any)
  for (; j < nj        ; j++) t->coef[j] = v[j-i]; // copy data
  for (; j < o2i[vhi+1]; j++) t->coef[j] = 0;      // left gap  (if any)

  // enlarge [lo,hi] if needed
  if (t->lo > vlo) t->lo = vlo;
  if (t->hi < vhi) t->hi = vhi;

  // activate nz bits that might be affected by vector content
  t->nz |= mad_bit_hcut(mad_bit_lcut(~0ull, vlo), vhi);

  // update nz for zeros hpoly (see seti for discussion)
  if (TPSA_STRICT_NZ > 1) FUN(update0)(t, vlo, vhi);

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

#endif // MAD_CTPSA_IMPL

// --- end --------------------------------------------------------------------o
