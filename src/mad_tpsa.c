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
  GTPSA are *defined* in [0,nc] with:
    - 0 < nc <= mc=idx(mo+1) ; mo <= ao <= d->mo <= 250
    - scalar GTPSA have nc=1.
    - new/clear/reset GTPSA have nc=1 with val[0]=0 (see reset0).
    - always true: val[0] is defined and nc >= 1, even if val[0] == 0.
------------------------------------------------------------------------------*/

// --- debugging --------------------------------------------------------------o

static inline log_t
check (const T *t)
{
  return !(!t->d || t->nc > t->mc || t->mo > t->ao || t->ao > t->d->mo ||
            t->mc > t->d->ord2idx[t->mo+1] || t->mo > mad_tpsa_dbgo ||
            t->idx[0] != 0);
}

/* Enabled only if TPSA_DEBUG > 0 (compile time)
   mad_tpsa_dbga  effects (runtime)
         0        none
         1        check ;
         2        check ; print header +tpsa(ok) or +raw(!ok)
         3        check ; print header +tpsa +raw
*/

int
FUN(debug) (const T *t, str_t name_, str_t fname_, int line_, FILE *stream_)
{
  assert(t);
  ord_t o; num_t r = 0;
  log_t ok = check(t);

  if (ok && mad_tpsa_dbga <= 1) return ok;

  const D* d = t->d;
  if (!stream_) stream_ = stdout;

  fprintf(stream_, "%s:%d: '%s' { nc=%d mc=%d mo=%d(%d) ao=%d uid=%d did=%d",
          fname_ ? fname_ : "??", line_, name_ ? name_ : "?",
          t->nc, t->mc, t->mo, mad_tpsa_dbgo, t->ao, t->uid, d ? d->id : -1);

  if (ok) {
    fprintf(stream_," } 0x%p\n", (void*)t); fflush(stream_);

    char name[48];
    strncpy(name, name_ ? name_ : t->nam, 48); name[47] = '\0';
    FUN(print)(t, name, 1e-40, 0, stream_);
    if (mad_tpsa_dbga <= 2) return ok;
  } else {
    fprintf(stream_," ** bug } 0x%p\n", (void*)t); fflush(stream_);
  }

  if (d) {
    const idx_t *o2i = d->ord2idx;
    idx_t ni = o2i[t->ao+1]; // corrupted TPSA cannot use print, display full memory
    FOR(i,ni) fprintf(stream_," [%d:%d]=" FMT "\n", i,d->ords[i],VAL(t->val[i]));
    fprintf(stream_,"\n"); fflush(stream_);
  }

  if (mad_tpsa_dbga <= 2) ensure(ok, "corrupted TPSA detected");
  return ok;
}

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

  const idx_t* o2i = t->d->ord2idx;
  ord_t ret = t->mo;
  if (mo_ < t->mo)
    t->mo = mo_, t->mc = o2i[t->mo+1], t->nc = MIN(t->nc, t->mc),
  else
    t->mo = MIN(mo_, t->ao), t->mc = o2i[t->mo+1];
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
FUN(len) (const T *t, log_t nc_)
{
  assert(t);
  return nc_ ? t->nc : t->mc;
}

ord_t
FUN(ord) (const T *t, log_t nc_)
{
  assert(t);
  return nc_ ? t->d->ords[t->nc] : t->mo;
}

log_t
FUN(isnul) (const T *t)
{
  assert(t);
  return t->nc == 1 && !t->val[0];
}

log_t
FUN(isval) (const T *t)
{
  assert(t);
  return t->nc == 1;
}

log_t
FUN(isvalid) (const T *t)
{
  assert(t);
  return check(t);
}

// --- init (unsafe) ----------------------------------------------------------o

T*
FUN(init) (T *t, const D *d, ord_t mo)
{
  assert(t && d && mo <= d->mo); DBGFUN(->);
  idx_t mc = d->ord2idx[mo+1];
  t->d = d, t->ao = t->mo = mo, t->mc = mc, t->uid=0, *t->nam=0;
  t->idx = (idx_t*)(t->val + t->mc), t->idx[0] = 0, FUN(reset0)(t);
#if TPSA_DEBUG
  if (mad_tpsa_dbga >= 3) FOR(i,1,mc) t->val[i] = M_PI;
#endif
  DBGTPSA(t); DBGFUN(<-); return t;
}

// --- ctors, dtor ------------------------------------------------------------o

T*
FUN(newd) (const D *d, ord_t mo)
{
  assert(d); DBGFUN(->);
  mo = MIN(mo, d->mo);
  idx_t mc = d->ord2idx[mo+1];
  T *r = mad_malloc(sizeof(T) + mc * (sizeof(NUM)+sizeof(idx_t))); assert(r);
  r->d = d, r->ao = r->mo = mo, t->mc = mc, r->uid=0, *r->nam=0;
  r->idx = (idx_t*)(r->val + t->mc), r->idx[0] = 0, FUN(reset0)(r);
#if TPSA_DEBUG
  if (mad_tpsa_dbga >= 3) FOR(i,1,mc) r->val[i] = M_PI;
#endif
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
FUN(setvar) (T *t, NUM v, idx_t iv, NUM scl)
{
  assert(t); DBGFUN(->);
  ensure(t->mo >= 1, "variables must have an order >= 1, got %d", t->mo);
  ensure(0 < iv && iv <= t->d->nv,
         "index 1<= %d <=%d is not a GTPSA variable", iv, t->d->nv);

  t->nc = 2, t->idx[1] = iv, t->val[0] = v, t->val[1] = scl ? scl : 1;
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setprm) (T *t, NUM v, idx_t ip)
{
  assert(t); DBGFUN(->);
  ensure(t->mo == 1 , "parameters must be a GPTSA of order 1, got %d", t->mo);
  ensure(0 < ip && ip <= t->d->np,
         "index 1<= %d <=%d is not a GPTSA parameter", ip, t->d->np);

  t->nc = 2, t->idx[1] = ip+t->d->nv, t->val[0] = v, t->val[1] = 1;
  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(setval) (T *t, NUM v)
{
  assert(t); DBGFUN(->);
  t->nc = 1, t->val[0] = v;
  DBGTPSA(t); DBGFUN(<-);
}

// --- copy, convert, swap ----------------------------------------------------o

void
FUN(copy) (const T *t, T *r)
{
  assert(t && r); DBGFUN(->);
  ensure(IS_COMPAT(t,r), "incompatibles GTPSA (descriptors differ)");

  if (t != r) {
    r->nc = MIN(t->nc, r->mc);
    FOR(i,r->nc) r->val[i] = t->val[i], r->idx[i] = t->idx[i];
  }
  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(sclord) (const T *t, T *r, log_t inv, log_t prm)
{
  assert(t && r); DBGFUN(->);

  FUN(copy)(t,r);

  const ord_t *co = r->d->ords, *po = r->d->prms;
  r->nc = MIN(t->nc, r->mc);
  if (inv) { // scale coefs by 1/order
    FOR(i,1,r->nc) r->val[i] /= co[r->idx[i]] - po[r->idx[i]] * !prm;
  } else {   // scale coefs by order
    FOR(i,1,r->nc) r->val[i] *= co[r->idx[i]] - po[r->idx[i]] * !prm;
  }
  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(clrord) (T *t, ord_t o)
{
  assert(t); DBGFUN(->);

  idx_t idx0 = t->d->ord2idx[o+0];
  idx_t idx1 = t->d->ord2idx[o+1];
  idx_t i0 = geti(t,idx0,0);
  idx_t i1 = geti(t,idx1,i0);
  idx_t ni = i0-i1;

  // left shift
  FOR(i,i1,t->nc) t->val[i+ni] = t->val[i], t->idx[i+ni] = t->idx[i];
  t->nc += ni;

  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(getord) (const T *t, T *r, ord_t o)
{
  assert(t && r); DBGFUN(->);
  ensure(IS_COMPAT(t,r), "incompatibles GTPSA (descriptors differ)");

  idx_t idx0 = t->d->ord2idx[o+0];
  idx_t idx1 = t->d->ord2idx[o+1];
  idx_t i0 = geti(t,idx0,0);
  idx_t i1 = geti(t,idx1,i0);

  FOR(i,i0,i1) r->val[i-i0] = t->val[i], r->idx[i-i0] = t->idx[i];
  r->nc = i1-i0;

  DBGTPSA(r); DBGFUN(<-);
}

void
FUN(cutord) (const T *t, T *r, int o)
{
  assert(t && r); DBGFUN(->);
  ensure(IS_COMPAT(t,r), "incompatibles GTPSA (descriptors differ)");

  if (o <= 0) {
    idx_t i0 = geti(t,t->d->ord2idx[1-o],0);
    FOR(i,i0,t->nc) r->val[i-i0] = t->val[i], r->idx[i-i0] = t->idx[i];
    r->nc = t->nc-i0;
  } else {
    idx_t i0 = geti(t,t->d->ord2idx[0+o],0);
    if (r != t) FOR(i,i0) r->val[i] = t->val[i], r->idx[i] = t->idx[i];
    r->nc = i0;
  }

  DBGTPSA(r); DBGFUN(<-);
}

idx_t
FUN(maxord) (const T *t, ssz_t n, idx_t idx_[n])
{
  assert(t); DBGFUN(->);
  ensure(FALSE, "NYI"); (void)n, (void)idx_;
/*
  if (idx_) { FOR(i,n) idx_[i] = -1; idx_[0] = 0; }

  idx_t mi = 0;                       // idx of mv
  num_t mv = fabs(t->coef[0]);        // max of all values
  ord_t hi = MIN(n-1, t->hi);
  TPSA_SCAN_O(t,t->lo,hi) {
    num_t mo = 0;                     // max of this order
    TPSA_SCAN(t,o,o) {
      if (mo < fabs(t->coef[i])) {
        mo = fabs(t->coef[i]);        // save max for this order
        if (idx_) idx_[o] = i;        // save idx for this order
        if (mv < mo) mv = mo, mi = i; // save max and idx for all orders
      }
    }
  }
*/
  DBGFUN(<-); return 0; // mi;
}

idx_t
FUN(cycle) (const T *t, idx_t i, ssz_t n, ord_t m_[n], idx_t *i_, NUM *v_)
{
  assert(t); DBGFUN(->);
  ensure(i >= -1 && i < t->d->nc, "index %d out of bounds", i);
  if (++i < t->nc) {
   if (i < 0) i = 0;
   if (!i && !t->val[0]) ++i;
   if (i_) *i_ = t->idx[i];
   if (v_) *v_ = t->val[t->idx[i]];
   if (m_) mad_mono_copy(MIN(n,t->d->nn), t->d->To[t->idx[i]], m_);
  } else i = -1;
  DBGFUN(<-); return i;
}

void
FUN(convert) (const T *t, T *r_, ssz_t n, idx_t t2r_[n], int pb)
{
  assert(t && r_); DBGFUN(->);
  ensure(pb >= -1 && pb <= 1,
         "invalid Poisson bracket direction %d, {-1, 0, 1} expected", pb);

  if (IS_COMPAT(t,r_) && !t2r_) { FUN(copy)(t,r_); DBGFUN(<-); return; }

  ensure(FALSE, "NYI"); (void)n;
/*
  T *r = t == r_ ? GET_TMPX(r_) : FUN(reset0)(r_);

  const D *td = t->d, *rd = r_->d;
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
    mad_mono_print(tn, tm, 0,0); printf(" -> "); mad_mono_print(rn, rm, 0,0);
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

  if (r != r_) { FUN(copy)(r,r_); REL_TMPX(r); }
  else         { DBGTPSA(r_); }
*/
  DBGFUN(<-);
}

// --- indexing / monomials ---------------------------------------------------o

ord_t
FUN(mono) (const T *t, idx_t idx, ssz_t n, ord_t m_[n], ord_t *p_)
{
  assert(t); DBGFUN(->);
  ensure(0 <= idx && idx < t->d->nc, "index %d out of bounds", idx);
  ord_t ret = mad_desc_mono(t->d, idx, n, m_, p_);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxs) (const T *t, ssz_t n, str_t s)
{
  assert(t && s); DBGFUN(->);
  idx_t ret = mad_desc_idxs(t->d, n, s);
  DBGFUN(<-); return ret;
}

idx_t
FUN(idxm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGFUN(->);
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

static inline idx_t
geti (const T *t, idx_t idx, idx_t i0) // TODO: use dichotomy
{
  idx_t i = i0;
  while (i < t->nc && t->idx[i] < idx) i++;
  return i;
}

static inline NUM
getv (const T *t, idx_t idx, idx_t i0)
{
  idx_t i = geti(t,idx,i0);
  return i < t->nc && t->idx[i] == idx ? t->val[i] : 0;
}

idx_t
FUN(getidx) (const T *t, idx_t i)
{
  assert(t); DBGFUN(->);
  ensure(0 <= i && i < t->d->nc, "index %d out of bounds", i);
  idx_t ret = i < t->nc ? t->idx[i] : -1;
  DBGFUN(<-); return ret;
}

NUM
FUN(geti) (const T *t, idx_t idx)
{
  assert(t); DBGFUN(->);
  NUM ret = !idx ? t->val[0] : getv(t,idx,1);
  DBGFUN(<-); return ret;
}

NUM
FUN(gets) (const T *t, ssz_t n, str_t s)
{ // --- mono is a string; represented as "[0-9A-Z]*"
  assert(t && s); DBGFUN(->);
  idx_t idx = mad_desc_idxs(t->d,n,s);
  ensure(idx >= 0, "invalid monomial");
  NUM ret = !idx ? t->val[0] : getv(t,idx,1);
  DBGFUN(<-); return ret;
}

NUM
FUN(getm) (const T *t, ssz_t n, const ord_t m[n])
{
  assert(t && m); DBGFUN(->);
  idx_t idx = mad_desc_idxm(t->d,n,m);
  ensure(idx >= 0, "invalid monomial");
  NUM ret = !idx ? t->val[0] : getv(t,idx,1);
  DBGFUN(<-); return ret;
}

NUM
FUN(getsm) (const T *t, ssz_t n, const idx_t m[n])
{ // --- mono is sparse; represented as [(i,o)]
  assert(t && m); DBGFUN(->);
  idx_t idx = mad_desc_idxsm(t->d,n,m);
  ensure(idx >= 0, "invalid monomial");
  NUM ret = !idx ? t->val[0] : getv(t,idx,1);
  DBGFUN(<-); return ret;
}

void
FUN(getv) (const T *t, idx_t idx, ssz_t n, NUM v[n])
{
  if (n <= 0) return;
  assert(t && v); DBGFUN(->);
  ssz_t nn = idx+n;
  ensure(0 <= idx && nn <= t->d->nc, "indexes %d:%d out of bounds", idx, nn);

  // clear vector
  FOR(i,n) v[i] = 0;

  // fill vector
  for (idx_t i = geti(t,idx,0); i < t->nc && t->idx[i] < nn; i++)
    v[t->idx[i]-idx] = t->val[i];

  DBGFUN(<-);
}

// --- setters ----------------------------------------------------------------o

void
FUN(seti) (T *t, idx_t idx, NUM a, NUM b)
{
  assert(t); DBGFUN(->);
  if (!idx) { t->val[0] = a*t->val[0]+b; DBGFUN(<-); return; }

  ensure(0 < idx && idx < t->d->nc, "index %d out of bounds", idx);

  // discard value
  if (idx >= t->mc) { DBGFUN(<-); return; }

  idx_t i = geti(t,idx,1);
  if (i < t->nc && t->idx[i] == idx) { // index exists
    t->val[i] = a*t->val[i]+b;
    if (!t->val[i]) {
      FOR(j,i,t->nc) t->val[j] = t->val[j+1], t->idx[j] = t->idx[j+1]; // shift left
      --t->nc;
    }
  } else if (b) {                      // insert index/value
    RFOR(j,t->nc,i) t->val[j+1] = t->val[j], t->idx[j+1] = t->idx[j]; // shift right
    t->val[i] = b, t->idx[i] = idx, t->nc += i == t->nc;
  }

  DBGTPSA(t); DBGFUN(<-);
}

void
FUN(sets) (T *t, ssz_t n, str_t s, NUM a, NUM b)
{
  assert(t && s); DBGFUN(->);
  idx_t idx = mad_desc_idxs(t->d,n,s);
  ensure(idx >= 0, "invalid monomial");
  FUN(seti)(t,idx,a,b);
  DBGFUN(<-);
}

void
FUN(setm) (T *t, ssz_t n, const ord_t m[n], NUM a, NUM b)
{
  assert(t && m); DBGFUN(->);
  idx_t idx = mad_desc_idxm(t->d,n,m);
  ensure(idx >= 0, "invalid monomial");
  FUN(seti)(t,idx,a,b);
  DBGFUN(<-);
}

void
FUN(setsm) (T *t, ssz_t n, const idx_t m[n], NUM a, NUM b)
{
  assert(t && m); DBGFUN(->);
  const idx_t idx = mad_desc_idxsm(t->d,n,m);
  ensure(idx >= 0, "invalid monomial");
  FUN(seti)(t,idx,a,b);
  DBGFUN(<-);
}

void
FUN(setv) (T *t, idx_t idx, ssz_t n, const NUM v[n])
{
  if (n <= 0) return;
  assert(t && v); DBGFUN(->);
  ssz_t nn = idx+n;
  ensure(0 <= idx && nn <= t->d->nc, "indexes %d:%d out of bounds", idx, nn);

  // truncate v to t
  if (nn > t->mc) {
    n -= nn-t->mc, nn = idx+n;
    if (n <= 0) { DBGFUN(<-); return; }
  }

  // find indexes of v != 0 => v' (linear in v)
  idx_t vi[n], nc = 0;
  FOR(i,n) if (v[i]) vi[nc++] = i;
  if (!nc) { DBGFUN(<-); return; }

  // find indexes of v in t and shift data (linear in t)
  idx_t i0 = geti(t,idx,0), i1 = geti(t,nn,i0), ni = i0-i1+nc;
  if (ni > 0) RFOR(i,t->nc,i1) t->val[i+ni]=t->val[i], t->idx[i+ni]=t->idx[i]; else
  if (ni < 0)  FOR(i,i1,t->nc) t->val[i+ni]=t->val[i], t->idx[i+ni]=t->idx[i];

  // copy v' in t (linear in v')
  FOR(i,nc) t->val[i0+i] = v[vi[i]], t->idx[i0+i] = idx+vi[i];
  t->nc += ni;

  DBGFUN(<-); return;
}

// --- copiers ----------------------------------------------------------------o

void
FUN(cpyi) (const T *t, T *r, idx_t idx)
{
  assert(t && r); DBGFUN(->);
  ensure(IS_COMPAT(t,r), "incompatibles GTPSA (descriptors differ)");
  if (!idx) { FUN(setval)(r, t->val[0]); DBGFUN(<-); return; }
  ensure(0 < idx && idx < t->d->nc, "index %d out of bounds", idx);
  NUM v = getv(t,idx,1); FUN(reset0)(r);
  if (v) t->val[1] = v, t->idx[1] = idx, t->nc = 2;
  DBGFUN(<-);
}

void
FUN(cpys) (const T *t, T *r, ssz_t n, str_t s)
{
  assert(t && r && s); DBGFUN(->);
  ensure(IS_COMPAT(t,r), "incompatibles GTPSA (descriptors differ)");
  idx_t idx = mad_desc_idxs(t->d,n,s);
  if (!idx) { FUN(setval)(r, t->val[0]); DBGFUN(<-); return; }
  ensure(idx > 0, "invalid monomial");
  NUM v = getv(t,idx,1); FUN(reset0)(r);
  if (v) t->val[1] = v, t->idx[1] = idx, t->nc = 2;
  DBGFUN(<-);
}

void
FUN(cpym) (const T *t, T *r, ssz_t n, const ord_t m[n])
{
  assert(t && r && m); DBGFUN(->);
  ensure(IS_COMPAT(t,r), "incompatibles GTPSA (descriptors differ)");
  idx_t idx = mad_desc_idxm(t->d,n,m);
  if (!idx) { FUN(setval)(r, t->val[0]); DBGFUN(<-); return; }
  ensure(idx > 0, "invalid monomial");
  NUM v = geti(t,idx,1); FUN(reset0)(r);
  if (v) t->val[1] = v, t->idx[1] = idx, t->nc = 2;
  DBGFUN(<-);
}

void
FUN(cpysm) (const T *t, T *r, ssz_t n, const idx_t m[n])
{
  assert(t && m); DBGFUN(->);
  ensure(IS_COMPAT(t,r), "incompatibles GTPSA (descriptors differ)");
  idx_t idx = mad_desc_idxsm(t->d,n,m);
  if (!idx) { FUN(setval)(r, t->val[0]); DBGFUN(<-); return; }
  ensure(idx > 0, "invalid monomial");
  NUM v = geti(t,idx,1); FUN(reset0)(r);
  if (v) t->val[1] = v, t->idx[1] = idx, t->nc = 2;
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
