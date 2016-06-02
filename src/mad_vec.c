/*
 o----------------------------------------------------------------------------o
 |
 | Vector module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
*/

#include <math.h>
#include <complex.h>
#include <assert.h>

#include "mad_log.h"
#include "mad_mem.h"
#include "mad_vec.h"

// --- implementation --------------------------------------------------------o

#define CHKR    assert( r )
#define CHKX    assert( x )
#define CHKXY   assert( x && y )
#define CHKXR   assert( x && r )
#define CHKYR   assert( y && r )
#define CHKXYR  assert( x && y && r )

#define NO(a) 
#define ID(a) a

#define CNUM2(re,im) (* (cnum_t*) & (num_t[2]) { re, im })
#define CNUMR(re)    CNUM2(re,0)
#define CNUMI(im)    CNUM2(0,im)
#define CNUM(a)      cnum_t a = CNUM2(MKNAME(a,_re), MKNAME(a,_im))

#define SET(OP)      for (size_t i=0; i < n; i++)  r[i] OP##= x[0]
#define CPY(OP)      for (size_t i=0; i < n; i++)  r[i] OP##= x[i]
#define FOLD(OP)     for (size_t i=0; i < n; i++)  r[0] OP##= x[i]
#define MAP(FN)      for (size_t i=0; i < n; i++)  r[i]  = FN(x[i])
#define MAP2(FN)     for (size_t i=0; i < n; i++)  r[i]  = FN(x[i], y[i])
#define VEC(OP)      for (size_t i=0; i < n; i++)  r[i]  = x[i] OP y[i]
#define VECS(OP)     for (size_t i=0; i < n; i++)  r[i]  = x[i] OP y
#define SVEC(OP)     for (size_t i=0; i < n; i++)  r[i]  = x    OP y[i]
#define DOT(C)       for (size_t i=0; i < n; i++)  r[0] += C(x[i]) * y[i]

// --- vec

void mad_vec_fill (num_t xx, num_t r[], size_t n)
{ CHKR; num_t *x=&xx; SET(); }

void mad_vec_copy (const num_t x[], num_t r[], size_t n)
{ CHKXR; CPY(); }

void mad_vec_copyv (const num_t x[], cnum_t r[], size_t n)
{ CHKXR; CPY(); }

void mad_vec_cvec (const num_t x[], const num_t y[], cnum_t r[], size_t n)
{ assert( r && (x || y) );
       if (x && y)       { MAP2(CNUM2); }
  else if (x)            { MAP (CNUMR); }
  else { const num_t *x=y; MAP (CNUMI); }
}

num_t mad_vec_dot (const num_t x[], const num_t y[], size_t n)
{ CHKXY; num_t r_=0, *r=&r_; DOT(ID); return *r; }

cnum_t mad_vec_dotv (const  num_t x[], const cnum_t y[], size_t n)
{ CHKXY; cnum_t r_=0, *r=&r_; DOT(ID); return *r;}

void mad_vec_dotv_r (const  num_t x[], const cnum_t y[], cnum_t *r, size_t n)
{ CHKXYR; *r=0; DOT(ID); }

void mad_vec_add (const num_t x[], const num_t y[], num_t r[], size_t n)
{ CHKXYR; VEC(+); }

void mad_vec_addn (const num_t x[], num_t y, num_t r[], size_t n)
{ CHKXR; VECS(+); }

void mad_vec_addc (const num_t x[], cnum_t y, cnum_t r[], size_t n)
{ CHKXR; VECS(+); }

void mad_vec_addc_r (const num_t x[], num_t y_re, num_t y_im, cnum_t r[], size_t n)
{ CHKXR; CNUM(y); VECS(+); }

void mad_vec_sub (const num_t x[], const num_t y[], num_t r[], size_t n)
{ CHKXYR; VEC(-); }

void mad_vec_subv (const num_t x[], const cnum_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(-); }

void mad_vec_subn (const num_t y[], num_t x, num_t r[], size_t n)
{ CHKYR; SVEC(-); }

void mad_vec_subc (const num_t y[], cnum_t x, cnum_t r[], size_t n)
{ CHKYR; SVEC(-); }

void mad_vec_subc_r (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], size_t n)
{ CHKYR; CNUM(x); SVEC(-); }

void mad_vec_mul (const num_t x[], const num_t y[], num_t r[], size_t n)
{ CHKXYR; VEC(*); }

void mad_vec_muln (const num_t x[], num_t y, num_t r[], size_t n)
{ CHKXR; VECS(*); }

void mad_vec_mulc (const num_t x[], cnum_t y, cnum_t r[], size_t n)
{ CHKXR; VECS(*); }

void mad_vec_mulc_r (const num_t x[], num_t y_re, num_t y_im, cnum_t r[], size_t n)
{ CHKXR; CNUM(y); VECS(*); }

void mad_vec_div (const num_t x[], const num_t y[], num_t r[], size_t n)
{ CHKXYR; VEC(/); }

void mad_vec_divv (const num_t x[], const cnum_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(/); }

void mad_vec_divn (const num_t y[], num_t x, num_t r[], size_t n)
{ CHKYR; SVEC(/); }

void mad_vec_divc (const num_t y[], cnum_t x, cnum_t r[], size_t n)
{ CHKYR; SVEC(/); }

void mad_vec_divc_r (const num_t y[], num_t x_re, num_t x_im, cnum_t r[], size_t n)
{ CHKYR; CNUM(x); SVEC(/); }

void mad_vec_center (const num_t x[], num_t r[], size_t n)
{ CHKXR; num_t m=0; { num_t *r=&m; FOLD(+); }
                    { num_t *x=&m; SET (-); } }

// --- cvec

void mad_cvec_fill (cnum_t xx, cnum_t r[], size_t n)
{ CHKR; cnum_t *x=&xx; SET(); }

void mad_cvec_fill_r (num_t xx_re, num_t xx_im, cnum_t r[], size_t n)
{ CHKR; CNUM(xx); { cnum_t *x=&xx; SET(); } }

void mad_cvec_copy (const cnum_t x[], cnum_t r[], size_t n)
{ CHKXR; CPY(); }

void mad_cvec_vec (const cnum_t x[], num_t re[], num_t ri[], size_t n)
{ CHKX; num_t *r;
  if (re) { r=re; MAP(creal); }
  if (ri) { r=ri; MAP(cimag); }
}

void mad_cvec_conj (const cnum_t x[], cnum_t r[], size_t n)
{ CHKXR; MAP(conj); }

cnum_t mad_cvec_dot (const cnum_t x[], const cnum_t y[], size_t n)
{ CHKXY; cnum_t r_=0, *r=&r_; DOT(conj); return *r; }

void mad_cvec_dot_r (const cnum_t x[], const cnum_t y[], cnum_t *r, size_t n)
{ CHKXYR; *r=0; DOT(conj); }

cnum_t mad_cvec_dotv (const cnum_t x[], const num_t y[], size_t n)
{ CHKXY; cnum_t r_=0, *r=&r_; DOT(conj); return *r; }

void mad_cvec_dotv_r (const cnum_t x[], const num_t y[], cnum_t *r, size_t n)
{ CHKXYR; *r=0; DOT(conj); }

void mad_cvec_add (const cnum_t x[], const cnum_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(+); }

void mad_cvec_addv (const cnum_t x[], const num_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(+); }

void mad_cvec_addn (const cnum_t x[], num_t y, cnum_t r[], size_t n)
{ CHKXR; VECS(+); }

void mad_cvec_addc (const cnum_t x[], cnum_t y, cnum_t r[], size_t n)
{ CHKXR; VECS(+); }

void mad_cvec_addc_r (const cnum_t x[], num_t y_re, num_t y_im, cnum_t r[], size_t n)
{ CHKXR; CNUM(y); VECS(+); }

void mad_cvec_sub (const cnum_t x[], const cnum_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(-); }

void mad_cvec_subv (const cnum_t x[], const num_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(-); }

void mad_cvec_subn (const cnum_t y[], num_t x, cnum_t r[], size_t n)
{ CHKYR; SVEC(-); }

void mad_cvec_subc (const cnum_t y[], cnum_t x, cnum_t r[], size_t n)
{ CHKYR; SVEC(-); }

void mad_cvec_subc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], size_t n)
{ CHKYR; CNUM(x); SVEC(-); }

void mad_cvec_mul (const cnum_t x[], const cnum_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(*); }

void mad_cvec_mulv (const cnum_t x[], const num_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(*); }

void mad_cvec_muln (const cnum_t x[], num_t y, cnum_t r[], size_t n)
{ CHKXR; VECS(*); }

void mad_cvec_mulc (const cnum_t x[], cnum_t y, cnum_t r[], size_t n)
{ CHKXR; VECS(*); }

void mad_cvec_mulc_r (const cnum_t x[], num_t y_re, num_t y_im, cnum_t r[], size_t n)
{ CHKXR; CNUM(y); VECS(*); }

void mad_cvec_div (const cnum_t x[], const cnum_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(/); }

void mad_cvec_divv (const cnum_t x[], const num_t y[], cnum_t r[], size_t n)
{ CHKXYR; VEC(/); }

void mad_cvec_divn (const cnum_t y[], num_t x, cnum_t r[], size_t n)
{ CHKYR; SVEC(/); }

void mad_cvec_divc (const cnum_t y[], cnum_t x, cnum_t r[], size_t n)
{ CHKYR; SVEC(/); }

void mad_cvec_divc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t r[], size_t n)
{ CHKYR; CNUM(x); SVEC(/); }

void mad_cvec_center (const cnum_t x[], cnum_t r[], size_t n)
{ CHKXR; cnum_t m=0; { cnum_t *r=&m; FOLD(+); }
                     { cnum_t *x=&m; SET (-); } }

// -- FFT ---------------------------------------------------------------------o

#include <fftw3.h>

void // x [n] -> r [n]
mad_vec_fft (const num_t x[], cnum_t r[], size_t n)
{
  CHKXR;
  mad_alloc_tmp(cnum_t, cx, n);
  mad_vec_copyv(x, cx, n);
  mad_cvec_fft(cx, r, n);
  mad_free_tmp(cx);
}

void // x [n] -> r [n/2+1]
mad_vec_rfft (const num_t x[], cnum_t r[], size_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_r2c_1d(n, (num_t*)x, r, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void // x [n/2+1] -> r [n]
mad_vec_irfft (const cnum_t x[], num_t r[], size_t n)
{
  CHKXR;
  size_t nn = n/2+1;
  mad_alloc_tmp(cnum_t, cx, nn);
  mad_cvec_copy(x, cx, nn);
  fftw_plan p = fftw_plan_dft_c2r_1d(n, cx, r, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  mad_free_tmp(cx);
  mad_vec_muln(r, 1.0/n, r, n);
}

void
mad_cvec_fft (const cnum_t x[], cnum_t r[], size_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_1d(n, (cnum_t*)x, r, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
}

void
mad_cvec_ifft(const cnum_t x[], cnum_t r[], size_t n)
{
  CHKXR;
  fftw_plan p = fftw_plan_dft_1d(n, (cnum_t*)x, r, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  mad_cvec_muln(r, 1.0/n, r, n);
}

// -- NFFT --------------------------------------------------------------------o

#if 0
#include <nfft3.h>

void
mad_cvec_nfft (const cnum_t x[], const num_t x_pos[], cnum_t r[], size_t n, size_t n_pos)
{
  assert( x && x_pos && r );
  nfft_plan p;
  nfft_init_1d(&p, n, n_pos);
  memcpy(p.x, x_pos, p.M_total * sizeof *x_pos); // TODO:  resample from n_pos to p.M_total?
  if(p.nfft_flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);
  memcpy(p.f_hat, x, p.N_total * sizeof *x); // TODO: resample from n to p.N_total?
  nfft_trafo(&p);
  memcpy(r, p.f, p.N_total * sizeof *r); // TODO: resample from p.N_total to n? 
  nfft_finalize(&p);
}

void
mad_cvec_infft (const cnum_t x[], const num_t x_pos[], cnum_t r[], size_t n, size_t n_pos)
{
  assert( x && x_pos && r );
  nfft_plan p;
  nfft_init_1d(&p, n, n_pos);
  memcpy(p.x, x_pos, p.M_total * sizeof *x_pos); // TODO:  resample from n_pos to p.M_total
  if(p.nfft_flags & PRE_ONE_PSI) nfft_precompute_one_psi(&p);
  memcpy(p.f, x, p.N_total * sizeof *x); // TODO: resample from n to p.N_total
  nfft_adjoint(&p);
  memcpy(r, p.f_hat, p.N_total * sizeof *r); // TODO: resample from p.N_total to n? 
  nfft_finalize(&p);
}

void
mad_cvec_nnfft (const cnum_t x[], const num_t x_pos[], const num_t f_pos[], cnum_t r[], size_t n, size_t n_pos)
{
  assert( x && x_pos && f_pos && r );
}

void
mad_cvec_innfft (const cnum_t x[], const num_t x_pos[], const num_t f_pos[], cnum_t r[], size_t n, size_t n_pos)
{
  assert( x && x_pos && f_pos && r );
}
#endif
