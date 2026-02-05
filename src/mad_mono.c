/*
 o-----------------------------------------------------------------------------o
 |
 | Monomial module implementation
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

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_num.h"
#include "mad_mono.h"

// --- implementation ---------------------------------------------------------o

str_t
mad_mono_prt (ssz_t n, const ord_t a[n], char s[n+1])
{
  assert(a && s);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) s[i] = a[i] + (a[i] < 10 ? '0' : (a[i] < 36 ? 'A'-10 : 'a'-36));
  s[n] = '\0';
  return s;
}

void
mad_mono_str (ssz_t n, ord_t a[n], str_t s)
{
  assert(a && s);
  ensure(n >= 0, "invalid monomial length");
  idx_t i = 0;
  for (; i < n && s[i]; ++i)
    a[i] = s[i] - (s[i] < 'A' ? '0' : (s[i] < 'a' ? 'A'-10 : 'a'-36));
  for (; i < n        ; ++i)
    a[i] = 0;
}

void
mad_mono_fill (ssz_t n, ord_t a[n], ord_t v)
{
  assert(a);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) a[i] = v;
}

void
mad_mono_copy (ssz_t n, const ord_t a[n], ord_t r[n])
{
  assert(a && r);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) r[i] = a[i];
}

ord_t
mad_mono_min (ssz_t n, const ord_t a[n])
{
  assert(a);
  ensure(n >= 0, "invalid monomial length");
  ord_t mo = ~0;
  FOR(i,n) if (a[i] < mo) mo = a[i];
  return mo;
}

ord_t
mad_mono_max (ssz_t n, const ord_t a[n])
{
  assert(a);
  ensure(n >= 0, "invalid monomial length");
  ord_t mo = 0;
  FOR(i,n) if (a[i] > mo) mo = a[i];
  return mo;
}

u32_t
mad_mono_ord (ssz_t n, const ord_t a[n])
{
  assert(a);
  ensure(n >= 0, "invalid monomial length");
  u32_t s = 0;
  FOR(i,n) s += a[i];
  return s;
}

num_t
mad_mono_ordp (ssz_t n, const ord_t a[n], idx_t stp)
{
  assert(a);
  ensure(n >= 0, "invalid monomial length");
  ensure(stp >= 1, "invalid step %d (>= 1)", stp);
  num_t p = 1;
  FOR(i,0,n,stp) p *= a[i];
  return p;
}

num_t
mad_mono_ordpf (ssz_t n, const ord_t a[n], idx_t stp)
{
  assert(a);
  ensure(n >= 0, "invalid monomial length");
  ensure(stp >= 1, "invalid step %d (>= 1)", stp);
  num_t p = 1;
  FOR(i,0,n,stp) p *= mad_num_fact(a[i]);
  return p;
}

log_t
mad_mono_eqn (ssz_t n, const ord_t a[n], ord_t b)
{
  assert(a);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) if (a[i] != b) return FALSE;
  return TRUE;
}

log_t
mad_mono_eq (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) if (a[i] != b[i]) return FALSE;
  return TRUE;
}

log_t
mad_mono_lt (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) if (a[i] >= b[i]) return FALSE;
  return TRUE;
}

log_t
mad_mono_le (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) if (a[i] > b[i]) return FALSE;
  return TRUE;
}

log_t
mad_mono_ok (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);       // b is the ref !
  ensure(n >= 0, "invalid monomial length");
  u32_t sa = 0; ord_t mb = 0; // sum(a), max(b)
  FOR(i,n) {
    if (a[i] > b[i]) return FALSE;
    sa += a[i], mb = MAX(mb, b[i]);
  }
  return sa <= mb;
}

log_t
mad_mono_ok_ (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);       // b is the ref !
  ensure(n >= 0, "invalid monomial length");
  u32_t sa = 0; ord_t mb = 0; // sum(a), max(b) when a[i] != 0 (new version)
  FOR(i,n) if (a[i]) {
    if (a[i] > b[i]) return FALSE;
    sa += a[i], mb = MAX(mb, b[i]);
  }
  return sa <= mb;
}

int
mad_mono_cmp (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) if (a[i] != b[i]) return (int)a[i] - b[i];
  return 0;
}

int
mad_mono_rcmp (ssz_t n, const ord_t a[n], const ord_t b[n])
{
  assert(a && b);
  ensure(n >= 0, "invalid monomial length");
  RFOR(i,n) if (a[i] != b[i]) return (int)a[i] - b[i];
  return 0;
}

void
mad_mono_add (ssz_t n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) r[i] = a[i] + b[i];
}

void
mad_mono_sub (ssz_t n, const ord_t a[n], const ord_t b[n], ord_t r[n])
{
  assert(a && b && r);
  ensure(n >= 0, "invalid monomial length");
  FOR(i,n) r[i] = a[i] > b[i] ? a[i] - b[i] : 0;
}

void
mad_mono_cat (ssz_t n, const ord_t a[n],
              ssz_t m, const ord_t b[m], ord_t r[n+m])
{
  mad_mono_copy(n, a, r  );
  mad_mono_copy(m, b, r+n);
}

void
mad_mono_rev (ssz_t n, const ord_t a[n], ord_t r[n])
{
  assert(a && r);
  ensure(n >= 0, "invalid monomial length");
  ord_t t;
  if (a != r) FOR(i,n)           r[i] = a[n-1-i];
  else        FOR(i,n/2) t=r[i], r[i] = r[n-1-i], r[n-1-i]=t; // swap
}


// -- printing

void
mad_mono_print (ssz_t n, const ord_t a[n], str_t sep_, FILE *fp_)
{
  assert(a);
  ensure(n >= 0, "invalid monomial length");
  char sep[6] = {' ', 0, '[', 0, ']', 0};

  if (!fp_) fp_ = stdout;

  if (sep_) FOR(i,3) if (*sep_) sep[i*2] = *sep_++ ;

  if (!n)  { fprintf(fp_, "%s%s", sep+2, sep+4); return; }

  n -= 1;    fprintf(fp_,               sep+2);
  FOR(i,n) { fprintf(fp_, "%d%s", a[i], sep+0); }
             fprintf(fp_, "%d%s", a[n], sep+4);
}

// --- end --------------------------------------------------------------------o
