/*
 o-----------------------------------------------------------------------------o
 |
 | Vector module implementation
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

#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_mem.h"
#include "mad_num.h"
#include "mad_vec.h"

// --- implementation ---------------------------------------------------------o

#define CHKR   assert( r )
#define CHKX   assert( x )
#define CHKXY  assert( x && y )
#define CHKXR  assert( x && r )
#define CHKYR  assert( y && r )
#define CHKXYR assert( x && y && r )

// --- vec

log_t mad_vec_isnul (const num_t x[], ssz_t n)
{ CHKX; return mad_vec_isval(x, 0, n); }

log_t mad_vec_isval (const num_t x[], num_t v, ssz_t n)
{ CHKX; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) { if (x[i] != v) return FALSE; } return n>0; }

void mad_vec_fill (num_t x, num_t r[], ssz_t n)
{ CHKR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x; }

void mad_vec_copy (const num_t x[], num_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  if (x == r) return;
  uintptr_t x_ = (uintptr_t)x, r_ = (uintptr_t)r;
  if (x_ > r_) FOR(i,n) r[i] = x[i];
  else        RFOR(i,n) r[i] = x[i]; }

void mad_vec_copyv (const num_t x[], cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i]; } // x and r cannot overlap

void mad_vec_cplx (const num_t re[], const num_t im[], cpx_t r[], ssz_t n)
{ assert(r && (re || im)); ensure(n>=0, "invalid vector/matrix size %d", n);
  if (re && im) FOR(i,n) r[i] = CPX2(re[i],im[i]);
  else if  (re) FOR(i,n) r[i] =      re[i]       ;
  else          FOR(i,n) r[i] = CPX2(0    ,im[i]);
}

void mad_vec_abs (const num_t x[], num_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = fabs(x[i]); }

num_t mad_vec_sum (const num_t x[], ssz_t n)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t r=x[0]; FOR(i,1,n) r += x[i]; return r; }

num_t mad_vec_mean (const num_t x[], ssz_t n)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  return mad_vec_sum(x,n)/n; }

num_t mad_vec_var (const num_t x[], ssz_t n)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  if (n == 1) return 0;
  num_t m = mad_vec_mean(x,n);
  num_t s=0, s2=0; FOR(i,n) s += x[i]-m, s2 += SQR(x[i]-m);
  return (s2 - SQR(s)/n)/(n-1); // Bessel's correction on centered values.
}

num_t mad_vec_sdev (const num_t x[], ssz_t n)
{ CHKX; return sqrt(mad_vec_var(x,n)); }

num_t mad_vec_dot (const num_t x[], const num_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t r=0; FOR(i,n) r += x[i] * y[i]; return r; }

num_t mad_vec_nrm (const num_t x[], ssz_t n)
{ CHKX; return sqrt(mad_vec_dot(x,x,n)); }

num_t mad_vec_rms (const num_t x[], ssz_t n)
{ CHKX; return sqrt(mad_vec_dot(x,x,n)/n); }

num_t mad_vec_dst (const num_t x[], const num_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t r=0; FOR(i,n) r += SQR(x[i] - y[i]);
  return sqrt(r);
}

num_t mad_vec_dstv (const num_t x[], const cpx_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t r=0; FOR(i,n) r += creal(conj(x[i]-y[i])*(x[i]-y[i]));
  return sqrt(r);
}

void mad_vec_add (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y[i]; }

void mad_vec_addn (const num_t x[], num_t y, num_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y; }

void mad_vec_addc (const num_t x[], cpx_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y; }

void mad_vec_addc_r (const num_t x[], num_t y_re, num_t y_im, cpx_t r[], ssz_t n)
{ CHKXR; mad_vec_addc(x, CPX(y), r, n); }

void mad_vec_sub (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] - y[i]; }

void mad_vec_subv (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] - y[i]; }

void mad_vec_subn (const num_t y[], num_t x, num_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x - y[i]; }

void mad_vec_subc (const num_t y[], cpx_t x, cpx_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x - y[i]; }

void mad_vec_subc_r (const num_t y[], num_t x_re, num_t x_im, cpx_t r[], ssz_t n)
{ CHKYR; mad_vec_subc(y, CPX(x), r, n); }

void mad_vec_mul (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y[i]; }

void mad_vec_muln (const num_t x[], num_t y, num_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y; }

void mad_vec_mulc (const num_t x[], cpx_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y; }

void mad_vec_mulc_r (const num_t x[], num_t y_re, num_t y_im, cpx_t r[], ssz_t n)
{ CHKXR; mad_vec_mulc(x, CPX(y), r, n); }

void mad_vec_div (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] / y[i]; }

void mad_vec_divv (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] / y[i]; }

void mad_vec_divn (const num_t x[], num_t y, num_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  ensure(y!=0, "invalid domain (zero divide)");
  FOR(i,n) r[i] = x[i] / y; }

void mad_vec_divc (const num_t x[], cpx_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  ensure(y!=0, "invalid domain (zero divide)");
  FOR(i,n) r[i] = x[i] / y; }

void mad_vec_divc_r (const num_t x[], num_t y_re, num_t y_im, cpx_t r[], ssz_t n)
{ CHKXR; mad_vec_divc(x, CPX(y), r, n); }

void mad_vec_invn (const num_t y[], num_t x, num_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x / y[i]; }

void mad_vec_invc (const num_t y[], cpx_t x, cpx_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x / y[i]; }

void mad_vec_invc_r (const num_t y[], num_t x_re, num_t x_im, cpx_t r[], ssz_t n)
{ CHKYR; mad_vec_invc(y, CPX(x), r, n); }

void mad_vec_dif (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = (x[i] - y[i]) / MAX(fabs(x[i]),1); }

void mad_vec_difv (const num_t x[], const cpx_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = (x[i] - y[i]) / MAX(fabs(x[i]),1); }

num_t mad_vec_eval (const num_t x[], num_t x0, ssz_t n) // Horner scheme
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t v=x[n-1]; FOR(i,1,n) v = v*x0 + x[n-1-i]; return v;
}

void mad_vec_minmax(const num_t x[], log_t absf, idx_t r[2], ssz_t n)
{ CHKXR; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t v[2];
  r[0] = r[1] = 0;
  if (absf) {
    v[0] = v[1] = fabs(x[0]);
    FOR(i,1,n) {
      num_t a = fabs(x[i]);
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  } else {
    v[0] = v[1] = x[0];
    FOR(i,1,n) {
      num_t a = x[i];
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  }
}

// Neumaier variants of Kahan sum

__attribute__((optimize("no-fast-math")))
num_t mad_vec_ksum (const num_t x[], ssz_t n)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t s = x[0], c = 0, t;
  FOR(i,1,n) {
    num_t xi = x[i];
    t = s + xi;
    if (fabs(s) >= fabs(xi))
      c += (s-t) + xi;
    else
      c += (xi-t) + s;
    s = t;
  }
  return s + c;
}

__attribute__((optimize("no-fast-math")))
cpx_t mad_cvec_ksum (const cpx_t x[], ssz_t n)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  cpx_t s = x[0], c = 0, t;
  FOR(i,1,n) {
    cpx_t xi = x[i];
    t = s + xi;
    if (cabs(s) >= cabs(xi))
      c += (s-t) + xi;
    else
      c += (xi-t) + s;
    s = t;
  }
  return s + c;
}

__attribute__((optimize("no-fast-math")))
num_t mad_vec_kdot (const num_t x[], const num_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t s = x[0]*y[0], c = 0, t, v;
  FOR(i,1,n) {
    v = x[i]*y[i];
    t = s + v;
    if (fabs(s) >= fabs(v))
      c += (s-t) + v;
    else
      c += (v-t) + s;
    s = t;
  }
  return s + c;
}

__attribute__((optimize("no-fast-math")))
cpx_t mad_cvec_kdot (const cpx_t x[], const cpx_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  cpx_t s = conj(x[0])*y[0], c = 0, t, v;
  FOR(i,1,n) {
    v = conj(x[i])*y[i];
    t = s + v;
    if (cabs(s) >= cabs(v))
      c += (s-t) + v;
    else
      c += (v-t) + s;
    s = t;
  }
  return s + c;
}

__attribute__((optimize("no-fast-math")))
cpx_t mad_cvec_kdotv (const cpx_t x[], const num_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  cpx_t s = conj(x[0])*y[0], c = 0, t, v;
  FOR(i,1,n) {
    v = conj(x[i])*y[i];
    t = s + v;
    if (cabs(s) >= cabs(v))
      c += (s-t) + v;
    else
      c += (v-t) + s;
    s = t;
  }
  return s + c;
}

// End of Neumaier variants of Kahan sum

void mad_vec_roll (num_t x[], ssz_t n, int nroll)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  nroll %= n;
  ssz_t nsz = iabs(nroll);
  mad_alloc_tmp(num_t, a, nsz);
  if (nroll > 0) {
    mad_vec_copy(x+n-nsz, a    ,   nsz); // end of x to a
    mad_vec_copy(x      , x+nsz, n-nsz); // shift x down (or right)
    mad_vec_copy(a      , x    ,   nsz); // a to beginning of x
  } else
  if (nroll < 0) {
    mad_vec_copy(x    , a      ,   nsz); // beginning of x to a
    mad_vec_copy(x+nsz, x      , n-nsz); // shift x up (or left)
    mad_vec_copy(a    , x+n-nsz,   nsz); // a to end of x
  }
  mad_free_tmp(a);
}

void mad_vec_kadd (int k, const num_t a[], const num_t *x[], num_t r[], ssz_t n)
{ assert(a && x && r);
  ensure(n>=0, "invalid vector/matrix size %d", n);
  ensure(k> 0, "invalid number of vectors");
  int j = k%8;

  switch(j) {
  case 0:
    assert(x[0] && x[1] && x[2] && x[3] && x[4] && x[5] && x[6] && x[7]);
    j = 8;
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i] + a[6]*x[6][i] + a[7]*x[7][i];
    break;
  case 1:
    assert(x[0]);
    FOR(i,n)
      r[i] = a[0]*x[0][i];
    break;
  case 2:
    assert(x[0] && x[1]);
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i];
    break;
  case 3:
    assert(x[0] && x[1] && x[2]);
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i];
    break;
  case 4:
    assert(x[0] && x[1] && x[2] && x[3]);
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i];
    break;
  case 5:
    assert(x[0] && x[1] && x[2] && x[3] && x[4]);
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i];
    break;
  case 6:
    assert(x[0] && x[1] && x[2] && x[3] && x[4] && x[5]);
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i];
    break;
  case 7:
    assert(x[0] && x[1] && x[2] && x[3] && x[4] && x[5] && x[6]);
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i] + a[6]*x[6][i];
    break;
  }

  for(; j < k; j+=8) {
    assert(x[j] && x[j+1] && x[j+2] && x[j+3] && x[j+4] && x[j+5] && x[j+6] && x[j+7]);
    FOR(i,n)
      r[i] += a[j+0]*x[j+0][i] + a[j+1]*x[j+1][i] + a[j+2]*x[j+2][i] + a[j+3]*x[j+3][i]
            + a[j+4]*x[j+4][i] + a[j+5]*x[j+5][i] + a[j+6]*x[j+6][i] + a[j+7]*x[j+7][i];
  }
}

void mad_vec_sort (num_t x[], idx_t r[], ssz_t n)
{
  CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t t;
  if (!r) {
    FOR(i,1,n) RFOR(j,i+1) if (x[j-1] > x[j]) SWAP(x[j-1], x[j], t);
    return;
  }

  // Set permutation indexes.
  FOR(i,n) r[i] = i;

  // Sort values by ascending order.
  idx_t tr;
  FOR(i,1,n) RFOR(j,i+1) if (x[j-1] > x[j]) {
    SWAP(x[j-1], x[j], t );
    SWAP(r[j-1], r[j], tr);
  }
}

// --- cvec

log_t mad_cvec_isnul (const cpx_t x[], ssz_t n)
{ CHKX; return mad_cvec_isval(x, 0, n); }

log_t mad_cvec_isval (const cpx_t x[], cpx_t v, ssz_t n)
{ CHKX; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) { if (x[i] != v) return FALSE; } return n>0; }

log_t mad_cvec_isval_r (const cpx_t x[], num_t v_re, num_t v_im, ssz_t n)
{ CHKX; return mad_cvec_isval (x, CPX(v), n); }

void mad_cvec_fill (cpx_t x, cpx_t r[], ssz_t n)
{ CHKR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x; }

void mad_cvec_fill_r (num_t x_re, num_t x_im, cpx_t r[], ssz_t n)
{ CHKR; mad_cvec_fill(CPX(x), r, n); }

void mad_cvec_copy (const cpx_t x[], cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  if (x == r) return;
  uintptr_t x_ = (uintptr_t)x, r_ = (uintptr_t)r;
  if (x_ > r_) FOR(i,n) r[i] = x[i];
  else        RFOR(i,n) r[i] = x[i]; }

void mad_cvec_reim (const cpx_t x[], num_t re[], num_t ri[], ssz_t n)
{ assert(x && (re || ri)); ensure(n>=0, "invalid vector/matrix size %d", n);
  if (re && ri) FOR(i,n) re[i]=creal(x[i]),
                         ri[i]=cimag(x[i]);
  else if (re)  FOR(i,n) re[i]=creal(x[i]);
  else          FOR(i,n) ri[i]=cimag(x[i]);
}

void mad_cvec_abs (const cpx_t x[], num_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = cabs(x[i]); }

void mad_cvec_conj (const cpx_t x[], cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = conj(x[i]); }

cpx_t mad_cvec_sum (const cpx_t x[], ssz_t n)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  cpx_t r=x[0]; FOR(i,1,n) r += x[i]; return r; }

void mad_cvec_sum_r (const cpx_t x[], cpx_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_sum(x,n); }

void mad_cvec_ksum_r (const cpx_t x[], cpx_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_ksum(x,n); }

cpx_t mad_cvec_dot (const cpx_t x[], const cpx_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  cpx_t r=0; FOR(i,n) r += conj(x[i])*y[i]; return r; }

void mad_cvec_dot_r (const cpx_t x[], const cpx_t y[], cpx_t *r, ssz_t n)
{ CHKR; *r = mad_cvec_dot(x,y,n); }

cpx_t mad_cvec_dotv (const cpx_t x[], const num_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  cpx_t r=0; FOR(i,n) r += conj(x[i])*y[i]; return r; }

void mad_cvec_dotv_r (const cpx_t x[], const num_t y[], cpx_t *r, ssz_t n)
{ CHKR; *r = mad_cvec_dotv(x,y,n); }

void mad_cvec_kdot_r (const cpx_t x[], const cpx_t y[], cpx_t *r, ssz_t n)
{ CHKR; *r = mad_cvec_kdot(x,y,n); }

void mad_cvec_kdotv_r (const cpx_t x[], const num_t y[], cpx_t *r, ssz_t n)
{ CHKR; *r = mad_cvec_kdotv(x,y,n); }

cpx_t mad_cvec_mean (const cpx_t x[], ssz_t n)
{ CHKX; return mad_cvec_sum(x,n)/n; }

void mad_cvec_mean_r (const cpx_t x[], cpx_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_mean(x,n); }

cpx_t mad_cvec_var (const cpx_t x[], ssz_t n)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  if (n == 1) return 0;
  cpx_t m = mad_cvec_mean(x,n);
  cpx_t s=0, s2=0; FOR(i,n) s += x[i]-m, s2 += SQR(x[i]-m);
  return (s2 - SQR(s)/n)/(n-1); // Bessel's correction on centered values.
}

void mad_cvec_var_r (const cpx_t x[], cpx_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_var(x,n); }

cpx_t mad_cvec_sdev (const cpx_t x[], ssz_t n)
{ CHKX; return sqrt(mad_cvec_var(x,n)); }

void mad_cvec_sdev_r (const cpx_t x[], cpx_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_sdev(x,n); }

cpx_t mad_cvec_rms (const cpx_t x[], ssz_t n)
{ CHKX; return sqrt(mad_cvec_dot(x,x,n)/n); }

void mad_cvec_rms_r (const cpx_t x[], cpx_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_rms(x,n); }

num_t mad_cvec_nrm (const cpx_t x[], ssz_t n)
{ CHKX; return sqrt(creal(mad_cvec_dot(x,x,n))); }

num_t mad_cvec_dst (const cpx_t x[], const cpx_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t r=0; FOR(i,n) r += creal(conj(x[i]-y[i])*(x[i]-y[i]));
  return sqrt(r);
}

num_t mad_cvec_dstv (const cpx_t x[], const num_t y[], ssz_t n)
{ CHKXY; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t r=0; FOR(i,n) r += creal(conj(x[i]-y[i])*(x[i]-y[i]));
  return sqrt(r);
}

void mad_cvec_add (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y[i]; }

void mad_cvec_addv (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y[i]; }

void mad_cvec_addn (const cpx_t x[], num_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y; }

void mad_cvec_addc (const cpx_t x[], cpx_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y; }

void mad_cvec_addc_r (const cpx_t x[], num_t y_re, num_t y_im, cpx_t r[], ssz_t n)
{ CHKXR; mad_cvec_addc(x, CPX(y), r, n); }

void mad_cvec_sub (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] - y[i]; }

void mad_cvec_subv (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] - y[i]; }

void mad_cvec_subn (const cpx_t y[], num_t x, cpx_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x - y[i]; }

void mad_cvec_subc (const cpx_t y[], cpx_t x, cpx_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x - y[i]; }

void mad_cvec_subc_r (const cpx_t y[], num_t x_re, num_t x_im, cpx_t r[], ssz_t n)
{ CHKYR; mad_cvec_subc(y, CPX(x), r, n); }

void mad_cvec_mul (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y[i]; }

void mad_cvec_mulv (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y[i]; }

void mad_cvec_muln (const cpx_t x[], num_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y; }

void mad_cvec_mulc (const cpx_t x[], cpx_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y; }

void mad_cvec_mulc_r (const cpx_t x[], num_t y_re, num_t y_im, cpx_t r[], ssz_t n)
{ CHKXR; mad_cvec_mulc(x, CPX(y), r, n); }

void mad_cvec_div (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] / y[i]; }

void mad_cvec_divv (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] / y[i]; }

void mad_cvec_divn (const cpx_t x[], num_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  ensure(y!=0, "invalid domain (zero divide)");
  FOR(i,n) r[i] = x[i] / y; }

void mad_cvec_divc (const cpx_t x[], cpx_t y, cpx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  ensure(y!=0, "invalid domain (zero divide)");
  FOR(i,n) r[i] = x[i] / y; }

void mad_cvec_divc_r (const cpx_t x[], num_t y_re, num_t y_im, cpx_t r[], ssz_t n)
{ CHKXR; mad_cvec_divc(x, CPX(y), r, n); }

void mad_cvec_invn (const cpx_t y[], num_t x, cpx_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x / y[i]; }

void mad_cvec_invc (const cpx_t y[], cpx_t x, cpx_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x / y[i]; }

void mad_cvec_invc_r (const cpx_t y[], num_t x_re, num_t x_im, cpx_t r[], ssz_t n)
{ CHKYR; mad_cvec_invc(y, CPX(x), r, n); }

void mad_cvec_dif  (const cpx_t x[], const cpx_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = (x[i] - y[i]) / MAX(cabs(x[i]),1); }

void mad_cvec_difv (const cpx_t x[], const num_t y[], cpx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = (x[i] - y[i]) / MAX(cabs(x[i]),1); }

cpx_t mad_cvec_eval (const cpx_t x[], cpx_t x0, ssz_t n) // Horner scheme
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  cpx_t v=x[n-1]; FOR(i,1,n) v = v*x0 + x[n-1-i]; return v; }

void mad_cvec_eval_r (const cpx_t x[], num_t x0_re, num_t x0_im, cpx_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_eval(x, CPX(x0), n); }

void mad_cvec_minmax(const cpx_t x[], idx_t r[2], ssz_t n)
{ CHKXR; ensure(n>0, "invalid vector/matrix size %d", n);
  num_t v[2];
  r[0] = r[1] = 0;
  v[0] = v[1] = cabs(x[0]);
  FOR(i,1,n) {
    num_t a = cabs(x[i]);
         if (a < v[0]) v[0]=a, r[0]=i;
    else if (a > v[1]) v[1]=a, r[1]=i;
  }
}

void mad_cvec_roll (cpx_t x[], ssz_t n, int nroll)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  nroll %= n;
  ssz_t nsz = iabs(nroll);
  mad_alloc_tmp(cpx_t, a, nsz);
  if (nroll > 0) {
    mad_cvec_copy(x+n-nsz, a    ,   nsz); // end of x to a
    mad_cvec_copy(x      , x+nsz, n-nsz); // shift x down (or right)
    mad_cvec_copy(a      , x    ,   nsz); // a to beginning of x
  } else
  if (nroll < 0) {
    mad_cvec_copy(x    , a      ,   nsz); // beginning of x to a
    mad_cvec_copy(x+nsz, x      , n-nsz); // shift x up (or left)
    mad_cvec_copy(a    , x+n-nsz,   nsz); // a to end of x
  }
  mad_free_tmp(a);
}

void mad_cvec_kadd (int k, const cpx_t a[], const cpx_t *x[], cpx_t r[], ssz_t n)
{ assert(a && x && r);
  ensure(n>=0, "invalid vector/matrix size %d", n);
  ensure(k> 0, "invalid number of vectors");
  int j = k%4;

  switch(j) {
  case 0:
    assert(x[0] && x[1] && x[2] && x[3]);
    j = 4;
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i];
    break;
  case 1:
    assert(x[0]);
    FOR(i,n)
      r[i] = a[0]*x[0][i];
    break;
  case 2:
    assert(x[0] && x[1]);
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i];
    break;
  case 3:
    assert(x[0] && x[1] && x[2]);
    FOR(i,n)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i];
    break;
  }

  for(; j < k; j+=4) {
    assert(x[j] && x[j+1] && x[j+2] && x[j+3]);
    FOR(i,n)
      r[i] += a[j+0]*x[j+0][i] + a[j+1]*x[j+1][i] + a[j+2]*x[j+2][i] + a[j+3]*x[j+3][i];
  }
}

// --- ivec

log_t mad_ivec_isnul (const idx_t x[], ssz_t n)
{ CHKX; return mad_ivec_isval(x, 0, n); }

log_t mad_ivec_isval (const idx_t x[], idx_t v, ssz_t n)
{ CHKX; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) { if (x[i] != v) return FALSE; } return n>0; }

void mad_ivec_fill (idx_t x, idx_t r[], ssz_t n)
{ CHKR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x; }

void mad_ivec_copy (const idx_t x[], idx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  if (x == r) return;
  uintptr_t x_ = (uintptr_t)x, r_ = (uintptr_t)r;
  if (x_ >  r_) FOR(i,n) r[i] = x[i];
  else         RFOR(i,n) r[i] = x[i]; }

void mad_ivec_add (const idx_t x[], const idx_t y[], idx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y[i]; }

void mad_ivec_addn (const idx_t x[], idx_t y, idx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] + y; }

void mad_ivec_sub (const idx_t x[], const idx_t y[], idx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] - y[i]; }

void mad_ivec_subn (const idx_t y[], idx_t x, idx_t r[], ssz_t n)
{ CHKYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x - y[i]; }

void mad_ivec_mul (const idx_t x[], const idx_t y[], idx_t r[], ssz_t n)
{ CHKXYR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y[i]; }

void mad_ivec_muln (const idx_t x[], idx_t y, idx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] * y; }

void mad_ivec_divn (const idx_t x[], idx_t y, idx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] / y; }

void mad_ivec_modn (const idx_t x[], idx_t y, idx_t r[], ssz_t n)
{ CHKXR; ensure(n>=0, "invalid vector/matrix size %d", n);
  FOR(i,n) r[i] = x[i] % y; }

void mad_ivec_minmax(const idx_t x[], log_t absf, idx_t r[2], ssz_t n)
{ CHKXR; ensure(n>0, "invalid vector/matrix size %d", n);
  idx_t v[2];
  r[0] = r[1] = 0;
  if (absf) {
    v[0] = v[1] = iabs(x[0]);
    FOR(i,1,n) {
      idx_t a = iabs(x[i]);
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  } else {
    v[0] = v[1] = x[0];
    FOR(i,1,n) {
      idx_t a = x[i];
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  }
}

void mad_ivec_roll (idx_t x[], ssz_t n, int nroll)
{ CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  nroll %= n;
  ssz_t nsz = iabs(nroll);
  mad_alloc_tmp(idx_t, a, nsz);
  if (nroll > 0) {
    mad_ivec_copy(x+n-nsz, a    ,   nsz); // end of x to a
    mad_ivec_copy(x      , x+nsz, n-nsz); // shift x down (or right)
    mad_ivec_copy(a      , x    ,   nsz); // a to beginning of x
  } else
  if (nroll < 0) {
    mad_ivec_copy(x    , a      ,   nsz); // beginning of x to a
    mad_ivec_copy(x+nsz, x      , n-nsz); // shift x up (or left)
    mad_ivec_copy(a    , x+n-nsz,   nsz); // a to end of x
  }
  mad_free_tmp(a);
}

ssz_t mad_ivec_sort (idx_t x[], ssz_t n, log_t rmdup)
{
  CHKX; ensure(n>0, "invalid vector/matrix size %d", n);
  idx_t ti;

  // Sort indexes by ascending order.
  FOR(i,1,n) RFOR(j,i+1) if (x[j-1] > x[j]) SWAP(x[j-1], x[j], ti);

  // Remove duplicates.
  if (rmdup) {
    idx_t k = 1;
    FOR(i,1,n) if (x[k-1] < x[i]) x[k++] = x[i];
    return k;
  }

  return n;
}
