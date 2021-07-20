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

#define CHKD    (d=MAX(1,d), n*=d)

#define CNUM(re,im) (* (cnum_t*) & (num_t[2]) { re, im })

// --- vector, cvector, ivector

struct  matrix { ssz_t nr, nc;  num_t data[]; };
struct cmatrix { ssz_t nr, nc; cnum_t data[]; };
struct imatrix { ssz_t nr, nc;  idx_t data[]; };

void mad_vec_append (struct matrix *x, num_t v)
{ CHKX; if (x->nc == 1) x->data[x->nr++] = v; else x->data[x->nc++] = v; }

void mad_ivec_append (struct imatrix *x, idx_t v)
{ CHKX; if (x->nc == 1) x->data[x->nr++] = v; else x->data[x->nc++] = v; }

void mad_cvec_append (struct cmatrix *x, cnum_t v)
{ CHKX; if (x->nc == 1) x->data[x->nr++] = v; else x->data[x->nc++] = v; }

void mad_cvec_append_r (struct cmatrix *x, num_t v_re, num_t v_im)
{ CHKX; mad_cvec_append(x, CNUM(v_re, v_im)); }

// --- vec

void mad_vec_zero (num_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = 0; }

void mad_vec_seq (num_t x, num_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = i+x; }

void mad_vec_fill (num_t x, num_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x; }

void mad_vec_copy (const num_t x[], num_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; if (x > r) for (idx_t i=0; i <  n; i+=d) r[  i] = x[  i];
          else if (x < r) for (idx_t i=d; i <= n; i+=d) r[n-i] = x[n-i]; }

void mad_vec_copyv (const num_t x[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; if (x > (num_t*)r) for (idx_t i=0; i <  n; i+=d) r[  i] = x[  i];
          else if (x < (num_t*)r) for (idx_t i=d; i <= n; i+=d) r[n-i] = x[n-i]; }

void mad_vec_cvec (const num_t x[], const num_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ assert( r && (x || y) ); CHKD;
  if (x && y) for (idx_t i=0; i < n; i+=d) r[i] = CNUM(x[i],y[i]);
  else if (x) for (idx_t i=0; i < n; i+=d) r[i] =      x[i]      ;
  else        for (idx_t i=0; i < n; i+=d) r[i] = CNUM(0   ,y[i]);
}

num_t mad_vec_abs (const num_t x[], num_t r_[], ssz_t n, ssz_t d)
{ CHKX; CHKD; num_t s=0;
  if (r_) for (idx_t i=0; i < n; i+=d) r_[i] = fabs(x[i]), s += r_[i];
  else    for (idx_t i=0; i < n; i+=d) s    += fabs(x[i]);
  return s;
}

num_t mad_vec_sum (const num_t x[], ssz_t n, ssz_t d)
{ CHKX; CHKD; num_t r=0; for (idx_t i=0; i < n; i+=d) r += x[i]; return r; }

num_t mad_vec_mean (const num_t x[], ssz_t n, ssz_t d)
{ CHKX; return mad_vec_sum(x,n,d)/n; }

num_t mad_vec_var (const num_t x[], ssz_t n, ssz_t d)
{ num_t xb = mad_vec_mean(x,n,d); CHKD;
  num_t r=0; for (idx_t i=0; i < n; i+=d) r += SQR(x[i]-xb);
  return sqrt(r/n);
}

num_t mad_vec_dot (const num_t x[], const num_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD; num_t r=0; for (idx_t i=0; i < n; i+=d) r += x[i] * y[i]; return r; }

cnum_t mad_vec_dotv (const num_t x[], const cnum_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD; cnum_t r=0; for (idx_t i=0; i < n; i+=d) r += x[i] * y[i]; return r; }

void mad_vec_dotv_r (const num_t x[], const cnum_t y[], cnum_t *r, ssz_t n, ssz_t d)
{ CHKR; *r = mad_vec_dotv(x,y,n,d); }

num_t mad_vec_norm (const num_t x[], ssz_t n, ssz_t d)
{ return sqrt(mad_vec_dot(x,x,n,d)); }

num_t mad_vec_dist (const num_t x[], const num_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD; num_t r=0; for (idx_t i=0; i < n; i+=d) r += SQR(x[i] - y[i]);
  return sqrt(r);
}

num_t mad_vec_distv (const num_t x[], const cnum_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD; num_t r=0;
  for (idx_t i=0; i < n; i+=d) r += conj(x[i]-y[i])*(x[i]-y[i]);
  return sqrt(r);
}

void mad_vec_add (const num_t x[], const num_t y[], num_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y[i]; }

void mad_vec_addn (const num_t x[], num_t y, num_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y; }

void mad_vec_addc (const num_t x[], cnum_t y, cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y; }

void mad_vec_addc_r (const num_t x[], num_t y_re, num_t y_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_vec_addc(x, CNUM(y_re,y_im), r, n, d); }

void mad_vec_sub (const num_t x[], const num_t y[], num_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] - y[i]; }

void mad_vec_subv (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] - y[i]; }

void mad_vec_subn (const num_t y[], num_t x, num_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x - y[i]; }

void mad_vec_subc (const num_t y[], cnum_t x, cnum_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x - y[i]; }

void mad_vec_subc_r (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_vec_subc(y, CNUM(x_re,x_im), r, n, d); }

void mad_vec_mul (const num_t x[], const num_t y[], num_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] * y[i]; }

void mad_vec_muln (const num_t x[], num_t y, num_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] * y; }

void mad_vec_mulc (const num_t x[], cnum_t y, cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] * y; }

void mad_vec_mulc_r (const num_t x[], num_t y_re, num_t y_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_vec_mulc(x, CNUM(y_re,y_im), r, n, d); }

void mad_vec_div (const num_t x[], const num_t y[], num_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] / y[i]; }

void mad_vec_divv (const num_t x[], const cnum_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] / y[i]; }

void mad_vec_divn (const num_t y[], num_t x, num_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x / y[i]; }

void mad_vec_divc (const num_t y[], cnum_t x, cnum_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x / y[i]; }

void mad_vec_divc_r (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t n, ssz_t d)
{  mad_vec_divc(y, CNUM(x_re,x_im), r, n, d); }

void mad_vec_center (const num_t x[], num_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD;
  num_t mu  = 0; for (idx_t i=0; i < n; i+=d) mu  += x[i];
        mu /= n; for (idx_t i=0; i < n; i+=d) r[i] = x[i] - mu;
}

num_t mad_vec_eval (const num_t x[], num_t x0, ssz_t n, ssz_t d) // Horner scheme
{ CHKX; CHKD; num_t v=x[n-d]; for (idx_t i=n-2*d; i >= 0; i-=d) v = v*x0 + x[i];
  return v;
}

void mad_vec_minmax(const num_t x[], log_t abs, idx_t r[2], ssz_t n, ssz_t d)
{ CHKXR; CHKD; num_t v[2];
  r[0] = r[1] = 0;
  if (abs) {
    v[0] = v[1] = fabs(x[0]);
    for (idx_t i=d; i < n; i+=d) {
      num_t a = fabs(x[i]);
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  } else {
    v[0] = v[1] = x[0];
    for (idx_t i=d; i < n; i+=d) {
      num_t a = x[i];
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  }
}

#pragma GCC push_options
#pragma GCC optimize ("O3")
// Neumaier variants of Kahan sum

num_t mad_vec_ksum (const num_t x[], ssz_t n, ssz_t d)
{ CHKX; CHKD;
  num_t s = x[0], c = 0, t;
  for (idx_t i=d; i < n; i+=d) {
    t = s + x[i];
    if (fabs(s) >= fabs(t))
      c = c + ((s-t) + x[i]);
    else
      c = c + ((x[i]-t) + s);
    s = t;
  }
  return s + c;
}

num_t mad_vec_knorm (const num_t x[], ssz_t n, ssz_t d)
{ CHKX; CHKD;
  num_t s = x[0]*x[0], c = 0, t, v;
  for (idx_t i=d; i < n; i+=d) {
    v = x[i]*x[i];
    t = s + v;
    if (s >= t)
      c = c + ((s-t) + v);
    else
      c = c + ((v-t) + s);
    s = t;
  }
  return sqrt(s + c);
}

num_t mad_vec_kdot (const num_t x[], const num_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD;
  num_t s = x[0]*y[0], c = 0, t, v;
  for (idx_t i=d; i < n; i+=d) {
    v = x[i]*y[i];
    t = s + v;
    if (fabs(s) >= fabs(t))
      c = c + ((s-t) + v);
    else
      c = c + ((v-t) + s);
    s = t;
  }
  return s + c;
}
#pragma GCC pop_options

void mad_vec_shift (num_t x[], ssz_t n, ssz_t d, int nshft)
{ CHKX; ssz_t dd=MAX(1,d);
  if (nshft > 0) mad_vec_copy(x, x+nshft*dd, n-nshft, d); // shift x down (or right)
  else
  if (nshft < 0) mad_vec_copy(x-nshft*dd, x, n+nshft, d); // shift x up (or left)
}

void mad_vec_roll (num_t x[], ssz_t n, ssz_t d, int nroll)
{ CHKX; ssz_t dd=MAX(1,d); nroll %= n;
  ssz_t nsz = abs(nroll);
  mad_alloc_tmp(num_t, a, nsz*dd);
  if (nroll > 0) {
    mad_vec_copy(x+(n-nsz)*dd, a       ,   nsz, d); // end of x to a
    mad_vec_copy(x           , x+nsz*dd, n-nsz, d); // shift x down (or right)
    mad_vec_copy(a           , x       ,   nsz, d); // a to beginning of x
  } else
  if (nroll < 0) {
    mad_vec_copy(x       , a           ,   nsz, d); // beginning of x to a
    mad_vec_copy(x+nsz*dd, x           , n-nsz, d); // shift x up (or left)
    mad_vec_copy(a       , x+(n-nsz)*dd,   nsz, d); // a to end of x
  }
  mad_free_tmp(a);
}

void mad_vec_kadd (int k, const num_t a[], const num_t *x[], num_t r[], ssz_t n, ssz_t d)
{ assert(a && x && r);
  if (k == 0) return;
  CHKD; int j = k%8;

  switch(j) {
  case 0: j = 8;
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i] + a[6]*x[6][i] + a[7]*x[7][i];
    break;
  case 1:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i];
    break;
  case 2:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i];
    break;
  case 3:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i];
    break;
  case 4:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i];
    break;
  case 5:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i];
    break;
  case 6:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i];
    break;
  case 7:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i]
           + a[4]*x[4][i] + a[5]*x[5][i] + a[6]*x[6][i];
    break;
  }

  for(; j < k; j+=8) {
    for (idx_t i=0; i < n; i+=d)
      r[i] += a[j+0]*x[j+0][i] + a[j+1]*x[j+1][i] + a[j+2]*x[j+2][i] + a[j+3]*x[j+3][i]
            + a[j+4]*x[j+4][i] + a[j+5]*x[j+5][i] + a[j+6]*x[j+6][i] + a[j+7]*x[j+7][i];
  }
}

// --- cvec

void mad_cvec_zero (cnum_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = 0; }

void mad_cvec_seq (cnum_t x, cnum_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = i+x; }

void mad_cvec_seq_r (num_t x_re, num_t x_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_cvec_seq(CNUM(x_re,x_im), r, n, d); }

void mad_cvec_fill (cnum_t x, cnum_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x; }

void mad_cvec_fill_r (num_t x_re, num_t x_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_cvec_fill(CNUM(x_re,x_im), r, n, d); }

void mad_cvec_copy (const cnum_t x[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; if (x > r) for (idx_t i=0; i <  n; i+=d) r[  i] = x[  i];
          else if (x < r) for (idx_t i=d; i <= n; i+=d) r[n-i] = x[n-i]; }

void mad_cvec_vec (const cnum_t x[], num_t re[], num_t ri[], ssz_t n, ssz_t d)
{ assert( x && (re || ri) ); CHKD;
  if (re && ri) for (idx_t i=0; i < n; i+=d) re[i]=creal(x[i]),
                                             ri[i]=cimag(x[i]);
  else if (re)  for (idx_t i=0; i < n; i+=d) re[i]=creal(x[i]);
  else          for (idx_t i=0; i < n; i+=d) ri[i]=cimag(x[i]);
}

num_t mad_cvec_abs (const cnum_t x[], num_t r_[], ssz_t n, ssz_t d)
{ CHKX; CHKD; num_t s=0;
  if (r_) for (idx_t i=0; i < n; i+=d) r_[i] = cabs(x[i]), s += r_[i];
  else    for (idx_t i=0; i < n; i+=d) s    += cabs(x[i]);
  return s;
}

void mad_cvec_conj (const cnum_t x[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = conj(x[i]); }

cnum_t mad_cvec_sum (const cnum_t x[], ssz_t n, ssz_t d)
{ CHKX; CHKD; cnum_t r=0; for (idx_t i=0; i < n; i+=d) r += x[i]; return r; }

void mad_cvec_sum_r (const cnum_t x[], cnum_t *r, ssz_t n, ssz_t d)
{ CHKXR; *r = mad_cvec_sum(x,n,d); }

cnum_t mad_cvec_mean (const cnum_t x[], ssz_t n, ssz_t d)
{ CHKX; return mad_cvec_sum(x,n,d)/n; }

void mad_cvec_mean_r (const cnum_t x[], cnum_t *r, ssz_t n, ssz_t d)
{ CHKXR; *r = mad_cvec_mean(x,n,d); }

cnum_t mad_cvec_var (const cnum_t x[], ssz_t n, ssz_t d)
{ cnum_t xb = mad_cvec_mean(x,n,d); CHKD;
  cnum_t r=0; for (idx_t i=0; i < n; i+=d) r += SQR(x[i]-xb);
  return csqrt(r/n);
}

void mad_cvec_var_r (const cnum_t x[], cnum_t *r, ssz_t n, ssz_t d)
{ CHKXR; *r = mad_cvec_var(x,n,d); }

cnum_t mad_cvec_dot (const cnum_t x[], const cnum_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD; cnum_t r=0; for (idx_t i=0; i < n; i+=d) r += conj(x[i])*y[i]; return r; }

void mad_cvec_dot_r (const cnum_t x[], const cnum_t y[], cnum_t *r, ssz_t n, ssz_t d)
{ CHKR; *r = mad_cvec_dot(x,y,n,d); }

cnum_t mad_cvec_dotv (const cnum_t x[], const num_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD; cnum_t r=0; for (idx_t i=0; i < n; i+=d) r += conj(x[i])*y[i]; return r; }

void mad_cvec_dotv_r (const cnum_t x[], const num_t y[], cnum_t *r, ssz_t n, ssz_t d)
{ CHKR; *r = mad_cvec_dotv(x,y,n,d); }

num_t mad_cvec_norm (const cnum_t x[], ssz_t n, ssz_t d)
{ return sqrt(mad_cvec_dot(x,x,n,d)); }

num_t mad_cvec_dist (const cnum_t x[], const cnum_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD; num_t r=0; for (idx_t i=0; i < n; i+=d) r += conj(x[i]-y[i])*(x[i]-y[i]);
  return sqrt(r);
}

num_t mad_cvec_distv (const cnum_t x[], const num_t y[], ssz_t n, ssz_t d)
{ CHKXY; CHKD; num_t r=0; for (idx_t i=0; i < n; i+=d) r += conj(x[i]-y[i])*(x[i]-y[i]);
  return sqrt(r);
}

void mad_cvec_add (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y[i]; }

void mad_cvec_addv (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y[i]; }

void mad_cvec_addn (const cnum_t x[], num_t y, cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y; }

void mad_cvec_addc (const cnum_t x[], cnum_t y, cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y; }

void mad_cvec_addc_r (const cnum_t x[], num_t y_re, num_t y_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_cvec_addc(x, CNUM(y_re,y_im), r, n, d); }

void mad_cvec_sub (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] - y[i]; }

void mad_cvec_subv (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] - y[i]; }

void mad_cvec_subn (const cnum_t y[], num_t x, cnum_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x - y[i]; }

void mad_cvec_subc (const cnum_t y[], cnum_t x, cnum_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x - y[i]; }

void mad_cvec_subc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_cvec_subc(y, CNUM(x_re,x_im), r, n, d); }

void mad_cvec_mul (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] * y[i]; }

void mad_cvec_mulv (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] * y[i]; }

void mad_cvec_muln (const cnum_t x[], num_t y, cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] * y; }

void mad_cvec_mulc (const cnum_t x[], cnum_t y, cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] * y; }

void mad_cvec_mulc_r (const cnum_t x[], num_t y_re, num_t y_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_cvec_mulc(x, CNUM(y_re,y_im), r, n, d); }

void mad_cvec_div (const cnum_t x[], const cnum_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] / y[i]; }

void mad_cvec_divv (const cnum_t x[], const num_t y[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] / y[i]; }

void mad_cvec_divn (const cnum_t y[], num_t x, cnum_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x / y[i]; }

void mad_cvec_divc (const cnum_t y[], cnum_t x, cnum_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x / y[i]; }

void mad_cvec_divc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], ssz_t n, ssz_t d)
{ mad_cvec_divc(y, CNUM(x_re,x_im), r, n, d); }

void mad_cvec_center (const cnum_t x[], cnum_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD;
  cnum_t mu  = 0; for (idx_t i=0; i < n; i+=d) mu  += x[i];
         mu /= n; for (idx_t i=0; i < n; i+=d) r[i] = x[i] - mu;
}

cnum_t mad_cvec_eval (const cnum_t x[], cnum_t x0, ssz_t n, ssz_t d) // Horner scheme
{ CHKX; CHKD; cnum_t v=x[n-d]; for (idx_t i=n-2*d; i >= 0; i-=d) v = v*x0 + x[i];
  return v;
}

void mad_cvec_eval_r (const cnum_t x[], num_t x0_re, num_t x0_im, cnum_t *r, ssz_t n, ssz_t d)
{ CHKXR; *r = mad_cvec_eval(x, CNUM(x0_re, x0_im), n, d); }

void mad_cvec_minmax(const cnum_t x[], idx_t r[2], ssz_t n, ssz_t d)
{ CHKXR; num_t v[2];
  r[0] = r[1] = 0;
  v[0] = v[1] = cabs(x[0]);
  for (idx_t i=d; i < n; i+=d) {
    num_t a = cabs(x[i]);
         if (a < v[0]) v[0]=a, r[0]=i;
    else if (a > v[1]) v[1]=a, r[1]=i;
  }
}

void mad_cvec_shift (cnum_t x[], ssz_t n, ssz_t d, int nshft)
{ CHKX; ssz_t dd=MAX(1,d);
  if (nshft > 0) mad_cvec_copy(x, x+nshft*dd, n-nshft, d); // shift x down (or right)
  else
  if (nshft < 0) mad_cvec_copy(x-nshft*dd, x, n+nshft, d); // shift x up (or left)
}

void mad_cvec_roll (cnum_t x[], ssz_t n, ssz_t d, int nroll)
{ CHKX; ssz_t dd=MAX(1,d); nroll %= n;
  ssz_t nsz = abs(nroll);
  mad_alloc_tmp(cnum_t, a, nsz*dd);
  if (nroll > 0) {
    mad_cvec_copy(x+(n-nsz)*dd, a       ,   nsz, d); // end of x to a
    mad_cvec_copy(x           , x+nsz*dd, n-nsz, d); // shift x down (or right)
    mad_cvec_copy(a           , x       ,   nsz, d); // a to beginning of x
  } else
  if (nroll < 0) {
    mad_cvec_copy(x       , a           ,   nsz, d); // beginning of x to a
    mad_cvec_copy(x+nsz*dd, x           , n-nsz, d); // shift x up (or left)
    mad_cvec_copy(a       , x+(n-nsz)*dd,   nsz, d); // a to end of x
  }
  mad_free_tmp(a);
}

void mad_cvec_kadd (int k, const cnum_t a[], const cnum_t *x[], cnum_t r[], ssz_t n, ssz_t d)
{ assert(a && x && r);
  if (k == 0) return;
  CHKD; int j = k%4;
  switch(j) {
  case 1:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i];
    break;
  case 2:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i];
    break;
  case 3:
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i];
    break;
  case 0: j = 4;
    for (idx_t i=0; i < n; i+=d)
      r[i] = a[0]*x[0][i] + a[1]*x[1][i] + a[2]*x[2][i] + a[3]*x[3][i];
  }
  for(; j < k; j+=4) {
    for (idx_t i=0; i < n; i+=d)
      r[i] += a[j+0]*x[j+0][i] + a[j+1]*x[j+1][i] + a[j+2]*x[j+2][i] + a[j+3]*x[j+3][i];
  }
}

// --- ivec

void mad_ivec_zero (idx_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = 0; }

void mad_ivec_seq (idx_t x, idx_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = i+x; }

void mad_ivec_fill (idx_t x, idx_t r[], ssz_t n, ssz_t d)
{ CHKR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x; }

void mad_ivec_copy (const idx_t x[], idx_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; if (x > r) for (idx_t i=0; i <  n; i+=d) r[  i] = x[  i];
          else if (x < r) for (idx_t i=d; i <= n; i+=d) r[n-i] = x[n-i]; }

void mad_ivec_copyv (const idx_t x[], num_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; if (x > (idx_t*)r) for (idx_t i=0; i <  n; i+=d) r[  i] = x[  i];
          else if (x < (idx_t*)r) for (idx_t i=d; i <= n; i+=d) r[n-i] = x[n-i]; }

void mad_ivec_add (const idx_t x[], const idx_t y[], idx_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y[i]; }

void mad_ivec_addn (const idx_t x[], idx_t y, idx_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y; }

void mad_ivec_sub (const idx_t x[], const idx_t y[], idx_t r[], ssz_t n, ssz_t d)
{ CHKXYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] - y[i]; }

void mad_ivec_subn (const idx_t y[], idx_t x, idx_t r[], ssz_t n, ssz_t d)
{ CHKYR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x - y[i]; }

void mad_ivec_muln (const idx_t x[], idx_t y, idx_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] * y; }

void mad_ivec_divn (const idx_t x[], idx_t y, idx_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] + y; }

void mad_ivec_modn (const idx_t x[], idx_t y, idx_t r[], ssz_t n, ssz_t d)
{ CHKXR; CHKD; for (idx_t i=0; i < n; i+=d) r[i] = x[i] % y; }

void mad_ivec_minmax(const idx_t x[], log_t absf, idx_t r[2], ssz_t n, ssz_t d)
{ CHKXR; CHKD; idx_t v[2];
  r[0] = r[1] = 0;
  if (absf) {
    v[0] = v[1] = abs(x[0]);
    for (idx_t i=d; i < n; i+=d) {
      idx_t a = abs(x[i]);
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  } else {
    v[0] = v[1] = x[0];
    for (idx_t i=d; i < n; i+=d) {
      idx_t a = x[i];
           if (a < v[0]) v[0]=a, r[0]=i;
      else if (a > v[1]) v[1]=a, r[1]=i;
    }
  }
}

void mad_ivec_shift (idx_t x[], ssz_t n, ssz_t d, int nshft)
{ CHKX; ssz_t dd=MAX(1,d);
  if (nshft > 0) mad_ivec_copy(x, x+nshft*dd, n-nshft, d); // shift x down (or right)
  else
  if (nshft < 0) mad_ivec_copy(x-nshft*dd, x, n+nshft, d); // shift x up (or left)
}

void mad_ivec_roll (idx_t x[], ssz_t n, ssz_t d, int nroll)
{ CHKX; ssz_t dd=MAX(1,d); nroll %= n;
  ssz_t nsz = abs(nroll);
  mad_alloc_tmp(idx_t, a, nsz*dd);
  if (nroll > 0) {
    mad_ivec_copy(x+(n-nsz)*dd, a       ,   nsz, d); // end of x to a
    mad_ivec_copy(x           , x+nsz*dd, n-nsz, d); // shift x down (or right)
    mad_ivec_copy(a           , x       ,   nsz, d); // a to beginning of x
  } else
  if (nroll < 0) {
    mad_ivec_copy(x       , a          ,   nsz, d); // beginning of x to a
    mad_ivec_copy(x+nsz*dd, x          , n-nsz, d); // shift x up (or left)
    mad_ivec_copy(a       , x+(n-nsz)*d,   nsz, d); // a to end of x
  }
  mad_free_tmp(a);
}
