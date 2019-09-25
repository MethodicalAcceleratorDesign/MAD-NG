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

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_mem.h"
#include "mad_vec.h"

// --- implementation ---------------------------------------------------------o

#define CHKR    assert( r )
#define CHKX    assert( x )
#define CHKXY   assert( x && y )
#define CHKXR   assert( x && r )
#define CHKYR   assert( y && r )
#define CHKXYR  assert( x && y && r )

#define SQR(a)      ((a)*(a))
#define CNUM(re,im) (* (cnum_t*) & (num_t[2]) { re, im })

// --- vec

void mad_vec_zero (num_t r[], ssz_t n)
{ CHKR; for (idx_t i=0; i < n; i++) r[i] = 0; }

void mad_vec_fill (num_t x, num_t r[], ssz_t n)
{ CHKR; for (idx_t i=0; i < n; i++) r[i] = x; }

void mad_vec_copy (const num_t x[], num_t r[], ssz_t n)
{ CHKXR; if (x > r) for (idx_t i=0; i <  n; i++) r[  i] = x[  i];
  else   if (x < r) for (idx_t i=1; i <= n; i++) r[n-i] = x[n-i]; }

void mad_vec_copyv (const num_t x[], cnum_t r[], ssz_t n)
{ CHKXR; if (x > (num_t*)r) for (idx_t i=0; i <  n; i++) r[  i] = x[  i];
  else   if (x < (num_t*)r) for (idx_t i=1; i <= n; i++) r[n-i] = x[n-i]; }

void mad_vec_cvec (const num_t x[], const num_t y[], cnum_t r[], ssz_t n)
{ assert( r && (x || y) );
  if (x && y) for (idx_t i=0; i < n; i++) r[i] = CNUM(x[i],y[i]);
  else if (x) for (idx_t i=0; i < n; i++) r[i] =      x[i]      ;
  else        for (idx_t i=0; i < n; i++) r[i] = CNUM(0   ,y[i]);
}

num_t mad_vec_dot (const num_t x[], const num_t y[], ssz_t n)
{ CHKXY; num_t r=0; for (idx_t i=0; i < n; i++) r += x[i] * y[i]; return r; }

cnum_t mad_vec_dotv (const num_t x[], const cnum_t y[], ssz_t n)
{ CHKXY; cnum_t r=0; for (idx_t i=0; i < n; i++) r += x[i] * y[i]; return r; }

void mad_vec_dotv_r (const num_t x[], const cnum_t y[], cnum_t *r, ssz_t n)
{ CHKR; *r = mad_vec_dotv(x, y, n); }

num_t mad_vec_norm (const num_t x[], ssz_t n)
{ return sqrt(mad_vec_dot(x, x, n)); }

num_t mad_vec_dist (const num_t x[], const num_t y[], ssz_t n)
{ CHKXY; num_t r=0;
  for (idx_t i=0; i < n; i++) r += SQR(x[i] - y[i]);
  return sqrt(r);
}

num_t mad_vec_distv (const num_t x[], const cnum_t y[], ssz_t n)
{ CHKXY; num_t r=0;
  for (idx_t i=0; i < n; i++) r += SQR(x[i] - creal(y[i])) + SQR(cimag(y[i]));
  return sqrt(r);
}

void mad_vec_add (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] + y[i]; }

void mad_vec_addn (const num_t x[], num_t y, num_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = x[i] + y; }

void mad_vec_addc (const num_t x[], cnum_t y, cnum_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = x[i] + y; }

void mad_vec_addc_r (const num_t x[], num_t y_re, num_t y_im, cnum_t r[], ssz_t n)
{ mad_vec_addc(x, CNUM(y_re,y_im), r, n); }

void mad_vec_sub (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] - y[i]; }

void mad_vec_subv (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] - y[i]; }

void mad_vec_subn (const num_t y[], num_t x, num_t r[], ssz_t n)
{ CHKYR; for (idx_t i=0; i < n; i++) r[i] = x - y[i]; }

void mad_vec_subc (const num_t y[], cnum_t x, cnum_t r[], ssz_t n)
{ CHKYR; for (idx_t i=0; i < n; i++) r[i] = x - y[i]; }

void mad_vec_subc_r (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t n)
{ mad_vec_subc(y, CNUM(x_re,x_im), r, n); }

void mad_vec_mul (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] * y[i]; }

void mad_vec_muln (const num_t x[], num_t y, num_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = x[i] * y; }

void mad_vec_mulc (const num_t x[], cnum_t y, cnum_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = x[i] * y; }

void mad_vec_mulc_r (const num_t x[], num_t y_re, num_t y_im, cnum_t r[], ssz_t n)
{ mad_vec_mulc(x, CNUM(y_re,y_im), r, n); }

void mad_vec_div (const num_t x[], const num_t y[], num_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] / y[i]; }

void mad_vec_divv (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] / y[i]; }

void mad_vec_divn (const num_t y[], num_t x, num_t r[], ssz_t n)
{ CHKYR; for (idx_t i=0; i < n; i++) r[i] = x / y[i]; }

void mad_vec_divc (const num_t y[], cnum_t x, cnum_t r[], ssz_t n)
{ CHKYR; for (idx_t i=0; i < n; i++) r[i] = x / y[i]; }

void mad_vec_divc_r (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t n)
{  mad_vec_divc(y, CNUM(x_re,x_im), r, n); }

void mad_vec_center (const num_t x[], num_t r[], ssz_t n)
{ CHKXR;
  num_t mu = 0;
  for (idx_t i=0; i < n; i++) mu += x[i];
  mu /= n;
  for (idx_t i=0; i < n; i++) r[i] = x[i] - mu;
}

num_t mad_vec_eval (const num_t x[], num_t x0, ssz_t n) // Horner scheme
{ CHKX;
  num_t v = x[n-1];
  for (idx_t i=n-2; i >= 0; i++) v = v*x0 + x[i];
  return v;
}

void mad_vec_minmax(const num_t x[], log_t abs, idx_t r[2], ssz_t n)
{ CHKXR; num_t v[2];
  r[0] = r[1] = 0;
  if (abs) {
    v[0] = v[1] = fabs(x[0]);
    for (idx_t i=1; i < n; i++) {
      num_t a = fabs(x[i]);
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  } else {
    v[0] = v[1] = x[0];
    for (idx_t i=1; i < n; i++) {
      num_t a = x[i];
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  }
}

#pragma GCC push_options
#pragma GCC optimize ("O3")
num_t mad_vec_sum (const num_t x[], ssz_t n) // Neumaier variant
{ CHKX;
  num_t s = x[0], c = 0, t;
  for (idx_t i=1; i < n; i++) {
    t = s + x[i];
    if (fabs(s) >= fabs(t))
      c = c + ((s-t) + x[i]);
    else
      c = c + ((x[i]-t) + s);
    s = t;
  }
  return s + c;
}
#pragma GCC pop_options

void mad_vec_shift (num_t x[], ssz_t n, int nshft)
{ CHKX;
  if (nshft > 0) mad_vec_copy(x, x+nshft, n-nshft); // shift x down (or right)
  else
  if (nshft < 0) mad_vec_copy(x-nshft, x, n+nshft); // shift x up (or left)
}

void mad_vec_roll (num_t x[], ssz_t n, int nroll)
{ CHKX; nroll %= n;
  ssz_t nsz = abs(nroll);
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
  if (k == 0) return;
  int j = k%8;

  switch(j) {
  case 0: j = 8;
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i] + a[6]*x[6][i] + a[7]*x[7][i];
    break;
  case 1:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i];
    break;
  case 2:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i];
    break;
  case 3:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i];
    break;
  case 4:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i];
    break;
  case 5:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i];
    break;
  case 6:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i];
    break;
  case 7:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i] + a[6]*x[6][i];
    break;
  }

  for(; j < k; j+=8) {
    for (idx_t i=0; i < n; i++)
      r[i] += a[j+0]*x[j+0][i] + a[j+1]*x[j+1][i] + a[j+2]*x[j+2][i] + a[j+3]*x[j+3][i]
            + a[j+4]*x[j+4][i] + a[j+5]*x[j+5][i] + a[j+6]*x[j+6][i] + a[j+7]*x[j+7][i];
  }
}

// --- cvec

void mad_cvec_zero (cnum_t r[], ssz_t n)
{ CHKR; for (idx_t i=0; i < n; i++) r[i] = 0; }

void mad_cvec_fill (cnum_t x, cnum_t r[], ssz_t n)
{ CHKR; for (idx_t i=0; i < n; i++) r[i] = x; }

void mad_cvec_fill_r (num_t x_re, num_t x_im, cnum_t r[], ssz_t n)
{ mad_cvec_fill(CNUM(x_re,x_im), r, n); }

void mad_cvec_copy (const cnum_t x[], cnum_t r[], ssz_t n)
{ mad_vec_copy((const num_t*)x, (num_t*)r, 2*n); }

void mad_cvec_shift (cnum_t x[], ssz_t n, int nshft)
{ mad_vec_shift((num_t*)x, 2*n, 2*nshft); }

void mad_cvec_roll (cnum_t x[], ssz_t n, int nroll)
{ mad_vec_roll((num_t*)x, 2*n, 2*nroll); }

void mad_cvec_vec (const cnum_t x[], num_t re[], num_t ri[], ssz_t n)
{ assert( x && (re || ri) );
  if (re && ri) for (idx_t i=0; i < n; i++) re[i]=creal(x[i]),
                                            ri[i]=cimag(x[i]);
  else if (re)  for (idx_t i=0; i < n; i++) re[i]=creal(x[i]);
  else          for (idx_t i=0; i < n; i++) ri[i]=cimag(x[i]);
}

void mad_cvec_conj (const cnum_t x[], cnum_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = conj(x[i]); }

cnum_t mad_cvec_dot (const cnum_t x[], const cnum_t y[], ssz_t n)
{ CHKXY; cnum_t r=0; for (idx_t i=0; i < n; i++) r += conj(x[i])*y[i]; return r; }

void mad_cvec_dot_r (const cnum_t x[], const cnum_t y[], cnum_t *r, ssz_t n)
{ CHKR; *r = mad_cvec_dot(x, y, n); }

cnum_t mad_cvec_dotv (const cnum_t x[], const num_t y[], ssz_t n)
{ CHKXY; cnum_t r=0; for (idx_t i=0; i < n; i++) r += conj(x[i])*y[i]; return r; }

void mad_cvec_dotv_r (const cnum_t x[], const num_t y[], cnum_t *r, ssz_t n)
{ CHKR; *r = mad_cvec_dotv(x, y, n); }

num_t mad_cvec_norm (const cnum_t x[], ssz_t n)
{ return mad_vec_norm((const num_t*)x, 2*n); }

num_t mad_cvec_dist (const cnum_t x[], const cnum_t y[], ssz_t n)
{ return mad_vec_dist((const num_t*)x, (const num_t*)y, 2*n); }

num_t mad_cvec_distv (const cnum_t x[], const num_t y[], ssz_t n)
{ CHKXY; num_t r=0;
  for (idx_t i=0; i < n; i++) r += SQR(creal(x[i]) - y[i]) + SQR(cimag(x[i]));
  return sqrt(r);
}

void mad_cvec_add (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] + y[i]; }

void mad_cvec_addv (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] + y[i]; }

void mad_cvec_addn (const cnum_t x[], num_t y, cnum_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = x[i] + y; }

void mad_cvec_addc (const cnum_t x[], cnum_t y, cnum_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = x[i] + y; }

void mad_cvec_addc_r (const cnum_t x[], num_t y_re, num_t y_im, cnum_t r[], ssz_t n)
{ mad_cvec_addc(x, CNUM(y_re,y_im), r, n); }

void mad_cvec_sub (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] - y[i]; }

void mad_cvec_subv (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] - y[i]; }

void mad_cvec_subn (const cnum_t y[], num_t x, cnum_t r[], ssz_t n)
{ CHKYR; for (idx_t i=0; i < n; i++) r[i] = x - y[i]; }

void mad_cvec_subc (const cnum_t y[], cnum_t x, cnum_t r[], ssz_t n)
{ CHKYR; for (idx_t i=0; i < n; i++) r[i] = x - y[i];; }

void mad_cvec_subc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t n)
{ mad_cvec_subc(y, CNUM(x_re,x_im), r, n); }

void mad_cvec_mul (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] * y[i]; }

void mad_cvec_mulv (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] * y[i]; }

void mad_cvec_muln (const cnum_t x[], num_t y, cnum_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = x[i] * y; }

void mad_cvec_mulc (const cnum_t x[], cnum_t y, cnum_t r[], ssz_t n)
{ CHKXR; for (idx_t i=0; i < n; i++) r[i] = x[i] * y; }

void mad_cvec_mulc_r (const cnum_t x[], num_t y_re, num_t y_im, cnum_t r[], ssz_t n)
{ mad_cvec_mulc(x, CNUM(y_re,y_im), r, n); }

void mad_cvec_div (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] / y[i]; }

void mad_cvec_divv (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t n)
{ CHKXYR; for (idx_t i=0; i < n; i++) r[i] = x[i] / y[i]; }

void mad_cvec_divn (const cnum_t y[], num_t x, cnum_t r[], ssz_t n)
{ CHKYR; for (idx_t i=0; i < n; i++) r[i] = x / y[i]; }

void mad_cvec_divc (const cnum_t y[], cnum_t x, cnum_t r[], ssz_t n)
{ CHKYR; for (idx_t i=0; i < n; i++) r[i] = x / y[i]; }

void mad_cvec_divc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t n)
{ mad_cvec_divc(y, CNUM(x_re,x_im), r, n); }

void mad_cvec_center (const cnum_t x[], cnum_t r[], ssz_t n)
{ CHKXR;
  cnum_t mu = 0;
  for (idx_t i=0; i < n; i++) mu += x[i];
  mu /= n;
  for (idx_t i=0; i < n; i++) r[i] = x[i] - mu;
}

cnum_t mad_cvec_eval (const cnum_t x[], cnum_t x0, ssz_t n) // Horner scheme
{ CHKX;
  cnum_t v = x[n-1];
  for (idx_t i=n-2; i >= 0; i++) v = v*x0 + x[i];
  return v;
}

void mad_cvec_eval_r (const cnum_t x[], num_t x0_re, num_t x0_im, cnum_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_eval(x, CNUM(x0_re, x0_im), n); }

void mad_cvec_minmax(const cnum_t x[], idx_t r[], ssz_t n)
{ CHKXR; num_t v[2];
  r[0] = r[1] = 0;
  v[0] = v[1] = cabs(x[0]);
  for (idx_t i=1; i < n; i++) {
    num_t a = cabs(x[i]);
         if (a < v[0]) v[0]=a, r[0]=i;
    else if (a > v[1]) v[1]=a, r[1]=i;
  }
}

#pragma GCC push_options
#pragma GCC optimize ("O3")
cnum_t mad_cvec_sum (const cnum_t x[], ssz_t n) // Neumaier variant
{ CHKX;
  cnum_t s = x[0], c = 0, t;
  for (idx_t i=1; i < n; i++) {
    t = s + x[i];
    if (fabs(creal(s))+fabs(cimag(s)) >= fabs(creal(t))+fabs(cimag(t)))
      c = c + ((s-t) + x[i]);
    else
      c = c + ((x[i]-t) + s);
    s = t;
  }
  return s + c;
}
#pragma GCC pop_options

void mad_cvec_sum_r (const cnum_t x[], cnum_t *r, ssz_t n)
{ CHKXR; *r = mad_cvec_sum(x, n); }

void mad_cvec_kadd (int k, const cnum_t a[], const cnum_t *x[], cnum_t r[], ssz_t n)
{ assert(a && x && r);
  if (k == 0) return;
  int j = k%4;
  switch(j) {
  case 1:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i];
    break;
  case 2:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i];
    break;
  case 3:
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i];
    break;
  case 0: j = 4;
    for (idx_t i=0; i < n; i++)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i];
  }
  for(; j < k; j+=4) {
    for (idx_t i=0; i < n; i++)
      r[i] += a[j+0]*x[j+0][i] + a[j+1]*x[j+1][i] + a[j+2]*x[j+2][i] + a[j+3]*x[j+3][i];
  }
}

