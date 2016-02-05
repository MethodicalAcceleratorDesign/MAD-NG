#ifndef MAD_MAT_H
#define MAD_MAT_H

/*
 o----------------------------------------------------------------------------o
 |
 | Matrix module interface
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
  
  Purpose:
  - wrappers around functions of real and complex matrices for LuaJIT

 o----------------------------------------------------------------------------o
 */

#include "mad.h"

// --- interface -------------------------------------------------------------o

void   mad_mat_ident   (                                           num_t  r[], size_t m, size_t n,             size_t ldr); // ident-> mat
void   mad_mat_set     (                         num_t x  ,        num_t  r[], size_t m, size_t n,             size_t ldr); //  num -> mat
void   mad_mat_copy    (const  num_t x[],                          num_t  r[], size_t m, size_t n, size_t ldx, size_t ldr); //  mat -> mat
void   mad_mat_copym   (const  num_t x[],                         cnum_t  r[], size_t m, size_t n, size_t ldx, size_t ldr); //  mat ->cmat
void   mad_mat_trans   (const  num_t x[],                          num_t  r[], size_t m, size_t n);                         //  mat.t()
num_t  mad_mat_dot     (const  num_t x[], const  num_t y[],                    size_t m, size_t n, size_t p);               // <mat ,  mat>
cnum_t mad_mat_dotm    (const  num_t x[], const cnum_t y[],                    size_t m, size_t n, size_t p);               // <mat , cmat>
void   mad_mat_dotm_r  (const  num_t x[], const cnum_t y[],       cnum_t *r  , size_t m, size_t n, size_t p);               // <mat , cmat>
void   mad_mat_mul     (const  num_t x[], const  num_t y[],        num_t  r[], size_t m, size_t n, size_t p);               //  mat *  mat
void   mad_mat_mulm    (const  num_t x[], const cnum_t y[],       cnum_t  r[], size_t m, size_t n, size_t p);               //  mat * cmat
int    mad_mat_invn    (const  num_t y[],        num_t x  ,        num_t  r[], size_t m, size_t n,           num_t rcond);  //  num /  mat
int    mad_mat_invc    (const  num_t y[],       cnum_t x  ,       cnum_t  r[], size_t m, size_t n,           num_t rcond);  // cnum /  mat
int    mad_mat_invc_r  (const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t m, size_t n,           num_t rcond);  // cnum /  mat
int    mad_mat_div     (const  num_t x[], const  num_t y[],        num_t  r[], size_t m, size_t n, size_t p, num_t rcond);  //  mat /  mat
int    mad_mat_divm    (const  num_t x[], const cnum_t y[],       cnum_t  r[], size_t m, size_t n, size_t p, num_t rcond);  //  mat / cmat
int    mad_mat_svd     (const  num_t x[], num_t u[], num_t s[],    num_t  v[], size_t m, size_t n);                         // u * s * v.t
int    mad_mat_eigen   (const  num_t x[], cnum_t w[], num_t vl[],  num_t vr[],           size_t n);                         //  w, vl, vr
void   mad_mat_fft     (const  num_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  mat ->cmat
void   mad_mat_rfft    (const  num_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  mat ->cmat
void   mad_mat_irfft   (const cnum_t x[],                          num_t  r[], size_t m, size_t n);                         // cmat -> mat
void   mad_mat_center  (const  num_t x[],                          num_t  r[], size_t m, size_t n);                         //  mat -> mat-<mat>_r

void   mad_cmat_ident  (                                          cnum_t  r[], size_t m, size_t n,             size_t ldr); //  ident->cmat
void   mad_cmat_set    (                        cnum_t x  ,       cnum_t  r[], size_t m, size_t n,             size_t ldr); //  cnum ->cmat
void   mad_cmat_set_r  (                  num_t x_re, num_t x_im, cnum_t  r[], size_t m, size_t n,             size_t ldr); //  cnum ->cmat
void   mad_cmat_copy   (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n, size_t ldx, size_t ldr); //  cmat ->cmat
void   mad_cmat_trans  (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  cmat.t()
void   mad_cmat_ctrans (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  cmat.ct()
cnum_t mad_cmat_dot    (const cnum_t x[], const cnum_t y[],                    size_t m, size_t n, size_t p);               // <cmat , cmat>
cnum_t mad_cmat_dotm   (const cnum_t x[], const  num_t y[],                    size_t m, size_t n, size_t p);               // <cmat ,  mat>
void   mad_cmat_dot_r  (const cnum_t x[], const cnum_t y[],       cnum_t *r  , size_t m, size_t n, size_t p);               // <cmat , cmat>
void   mad_cmat_dotm_r (const cnum_t x[], const  num_t y[],       cnum_t *r  , size_t m, size_t n, size_t p);               // <cmat ,  mat>
void   mad_cmat_mul    (const cnum_t x[], const cnum_t y[],       cnum_t  r[], size_t m, size_t n, size_t p);               //  cmat * cmat
void   mad_cmat_mulm   (const cnum_t x[], const  num_t y[],       cnum_t  r[], size_t m, size_t n, size_t p);               //  cmat *  mat
int    mad_cmat_invn   (const cnum_t y[],        num_t x  ,       cnum_t  r[], size_t m, size_t n,           num_t rcond);  //   num / cmat
int    mad_cmat_invc   (const cnum_t y[],       cnum_t x  ,       cnum_t  r[], size_t m, size_t n,           num_t rcond);  //  cnum / cmat
int    mad_cmat_invc_r (const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t m, size_t n,           num_t rcond);  //  cnum / cmat
int    mad_cmat_div    (const cnum_t x[], const cnum_t y[],       cnum_t  r[], size_t m, size_t n, size_t p, num_t rcond);  //  cmat / cmat
int    mad_cmat_divm   (const cnum_t x[], const  num_t y[],       cnum_t  r[], size_t m, size_t n, size_t p, num_t rcond);  //  cmat /  mat
int    mad_cmat_svd    (const cnum_t x[], cnum_t u[], num_t s[],  cnum_t  v[], size_t m, size_t n);                         // u * s * v.t
int    mad_cmat_eigen  (const cnum_t x[], cnum_t w[], cnum_t vl[],cnum_t vr[],           size_t n);                         // w, vl, vr
void   mad_cmat_fft    (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  cmat ->cmat
void   mad_cmat_ifft   (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  cmat ->cmat
void   mad_cmat_center (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  cmat ->cmat-<cmat>_r

// ---------------------------------------------------------------------------o

#endif // MAD_MAT_H