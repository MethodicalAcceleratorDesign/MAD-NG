#ifndef MAD_MAT_H
#define MAD_MAT_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Matrix module interface
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

  Purpose:
  - wrappers around functions of real and complex matrices for LuaJIT

 o-----------------------------------------------------------------------------o
 */

#include "mad_def.h"

// --- interface --------------------------------------------------------------o

void  mad_mat_rev      (      num_t x[],                                   ssz_t m, ssz_t n, int d);                //  mat -> rev(mat)
void  mad_mat_center   (      num_t x[],                                   ssz_t m, ssz_t n, int d);                //  mat -> mat-<mat>
void  mad_mat_roll     (      num_t x[],                                   ssz_t m, ssz_t n, int mroll, int nroll); //  mat -> mat
void  mad_mat_eye      (      num_t x[],                        num_t v  , ssz_t m, ssz_t n, ssz_t ldx);            //  mat -> mat
void  mad_mat_copy     (const num_t x[],                        num_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr); //  mat -> mat
void  mad_mat_copym    (const num_t x[],                        cpx_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr); //  mat ->cmat
void  mad_mat_trans    (const num_t x[],                        num_t r[], ssz_t m, ssz_t n);                       //  mat.t()
void  mad_mat_mul      (const num_t x[], const num_t y[],       num_t r[], ssz_t m, ssz_t n, ssz_t p);              //  mat *  mat
void  mad_mat_mulm     (const num_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  mat * cmat
void  mad_mat_tmul     (const num_t x[], const num_t y[],       num_t r[], ssz_t m, ssz_t n, ssz_t p);              //  mat'*  mat
void  mad_mat_tmulm    (const num_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  mat'* cmat
void  mad_mat_mult     (const num_t x[], const num_t y[],       num_t r[], ssz_t m, ssz_t n, ssz_t p);              //  mat *  mat'
void  mad_mat_multm    (const num_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  mat * cmat'
void  mad_mat_dmul     (const num_t x[], const num_t y[],       num_t r[], ssz_t m, ssz_t n, ssz_t p);              //  diag(mat) *  mat
void  mad_mat_dmulm    (const num_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  diag(mat) * cmat
void  mad_mat_muld     (const num_t x[], const num_t y[],       num_t r[], ssz_t m, ssz_t n, ssz_t p);              //  mat * diag( mat)
void  mad_mat_muldm    (const num_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  mat * diag(cmat)
int   mad_mat_det      (const num_t x[],                        num_t *r ,          ssz_t n);                       //  det(mat)
int   mad_mat_invn     (const num_t y[],       num_t x,         num_t r[], ssz_t m, ssz_t n,          num_t rcond); //  num /  mat
int   mad_mat_invc     (const num_t y[],       cpx_t x,         cpx_t r[], ssz_t m, ssz_t n,          num_t rcond); //  cpx /  mat
int   mad_mat_invc_r   (const num_t y[], num_t x_re,num_t x_im, cpx_t r[], ssz_t m, ssz_t n,          num_t rcond); //  cpx /  mat
int   mad_mat_div      (const num_t x[], const num_t y[],       num_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond); //  mat /  mat
int   mad_mat_divm     (const num_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond); //  mat / cmat
int   mad_mat_solve    (const num_t a[], const num_t b[],       num_t x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond); //  min|b-ax| (QR)
int   mad_mat_nsolve   (const num_t a[], const num_t b[],       num_t x[], ssz_t m, ssz_t n, ssz_t N, num_t rcond, num_t r_[]); // min|b-ax| (MICADO)
int   mad_mat_ssolve   (const num_t a[], const num_t b[],       num_t x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond, num_t s_[]); // min|b-ax| (SVD)
int   mad_mat_gsolve   (const num_t a[], const num_t b[], const num_t c[], const num_t d[],
                                                                num_t x[], ssz_t m, ssz_t n, ssz_t p, num_t *nrm_);
int   mad_mat_gmsolve  (const num_t a[], const num_t b[], const num_t d[],
                                               num_t x[],       num_t y[], ssz_t m, ssz_t n, ssz_t p);
int   mad_mat_pcacnd   (const num_t a[],       idx_t c[],                  ssz_t m, ssz_t n, ssz_t N, num_t cut, num_t s_[]);
int   mad_mat_svdcnd   (const num_t a[],       idx_t c[],                  ssz_t m, ssz_t n, ssz_t N, num_t cut, num_t s_[], num_t tol);
int   mad_mat_svd      (const num_t x[], num_t u[], num_t s[],  num_t v[], ssz_t m, ssz_t n);                       //  u * s * v.t
int   mad_mat_eigen    (const num_t x[], cpx_t w[], num_t vl[], num_t vr[],         ssz_t n);                       //  w, vl, vr
void  mad_mat_fft      (const num_t x[],                        cpx_t r[], ssz_t m, ssz_t n);                       //  mat ->cmat
void  mad_mat_rfft     (const num_t x[],                        cpx_t r[], ssz_t m, ssz_t n);                       //  mat ->cmat
void  mad_mat_nfft     (const num_t x[], const num_t x_node[],  cpx_t r[], ssz_t m, ssz_t n, ssz_t nr);
void  mad_mat_sympconj (const num_t x[],                        num_t r[],          ssz_t n);                       //  -J M' J
num_t mad_mat_symperr  (const num_t x[],                        num_t r[],          ssz_t n, num_t *tol_);          //  M' J M - J

void  mad_cmat_rev     (      cpx_t x[],                                   ssz_t m, ssz_t n, int d);                //  cmat -> rev(cmat)
void  mad_cmat_center  (      cpx_t x[],                                   ssz_t m, ssz_t n, int d);                //  cmat ->cmat-<cmat>
void  mad_cmat_roll    (      cpx_t x[],                                   ssz_t m, ssz_t n, int mroll, int nroll); //  cmat ->cmat
void  mad_cmat_eye     (      cpx_t x[],                        cpx_t v  , ssz_t m, ssz_t n, ssz_t ldx);            //  cmat ->cmat
void  mad_cmat_eye_r   (      cpx_t x[], num_t v_re,num_t v_im,            ssz_t m, ssz_t n, ssz_t ldx);            //  cmat ->cmat
void  mad_cmat_copy    (const cpx_t x[],                        cpx_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr); //  cmat ->cmat
void  mad_cmat_trans   (const cpx_t x[],                        cpx_t r[], ssz_t m, ssz_t n);                       //  cmat.t()
void  mad_cmat_ctrans  (const cpx_t x[],                        cpx_t r[], ssz_t m, ssz_t n);                       //  cmat.ct()
void  mad_cmat_mul     (const cpx_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  cmat * cmat
void  mad_cmat_mulm    (const cpx_t x[], const num_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  cmat *  mat
void  mad_cmat_tmul    (const cpx_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  cmat'* cmat
void  mad_cmat_tmulm   (const cpx_t x[], const num_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  cmat'*  mat
void  mad_cmat_mult    (const cpx_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  cmat * cmat'
void  mad_cmat_multm   (const cpx_t x[], const num_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  cmat *  mat'
void  mad_cmat_dmul    (const cpx_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  diag(cmat) * cmat
void  mad_cmat_dmulm   (const cpx_t x[], const num_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  diag(cmat) *  mat
void  mad_cmat_muld    (const cpx_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  cmat * diag(cmat)
void  mad_cmat_muldm   (const cpx_t x[], const num_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p);              //  cmat * diag( mat)
int   mad_cmat_det     (const cpx_t x[],                        cpx_t *r ,          ssz_t n);                       //  det(cmat)
int   mad_cmat_invn    (const cpx_t y[],       num_t x  ,       cpx_t r[], ssz_t m, ssz_t n,          num_t rcond); //  num  / cmat
int   mad_cmat_invc    (const cpx_t y[],       cpx_t x  ,       cpx_t r[], ssz_t m, ssz_t n,          num_t rcond); //  cpx  / cmat
int   mad_cmat_invc_r  (const cpx_t y[], num_t x_re,num_t x_im, cpx_t r[], ssz_t m, ssz_t n,          num_t rcond); //  cpx  / cmat
int   mad_cmat_div     (const cpx_t x[], const cpx_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond); //  cmat / cmat
int   mad_cmat_divm    (const cpx_t x[], const num_t y[],       cpx_t r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond); //  cmat /  mat
int   mad_cmat_solve   (const cpx_t a[], const cpx_t b[],       cpx_t x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond); //  min|b-ax| (QR)
int   mad_cmat_ssolve  (const cpx_t a[], const cpx_t b[],       cpx_t x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond, num_t s_[]); // min|b-ax| (SVD)
int   mad_cmat_gsolve  (const cpx_t a[], const cpx_t b[], const cpx_t c[], const cpx_t d[],
                                                                cpx_t x[], ssz_t m, ssz_t n, ssz_t p, num_t *nrm_);
int   mad_cmat_gmsolve (const cpx_t a[], const cpx_t b[], const cpx_t d[],
                                               cpx_t x[],       cpx_t y[], ssz_t m, ssz_t n, ssz_t p);
int   mad_cmat_pcacnd  (const cpx_t a[],       idx_t c[],                  ssz_t m, ssz_t n, ssz_t N, num_t cut, num_t s_[]);
int   mad_cmat_svdcnd  (const cpx_t a[],       idx_t c[],                  ssz_t m, ssz_t n, ssz_t N, num_t cut, num_t s_[], num_t tol);
int   mad_cmat_svd     (const cpx_t x[], cpx_t u[], num_t s[],  cpx_t v[], ssz_t m, ssz_t n);                       //  u * s * v.t
int   mad_cmat_eigen   (const cpx_t x[], cpx_t w[], cpx_t vl[], cpx_t vr[],         ssz_t n);                       //  w, vl, vr
void  mad_cmat_fft     (const cpx_t x[],                        cpx_t r[], ssz_t m, ssz_t n);                       //  cmat ->cmat
void  mad_cmat_nfft    (const cpx_t x[], const num_t x_node[],  cpx_t r[], ssz_t m, ssz_t n, ssz_t nr);
void  mad_cmat_ifft    (const cpx_t x[],                        cpx_t r[], ssz_t m, ssz_t n);                       //  cmat ->cmat
void  mad_cmat_irfft   (const cpx_t x[],                        num_t r[], ssz_t m, ssz_t n);                       //  cmat -> mat
void  mad_cmat_infft   (const cpx_t x[], const num_t r_node[],  cpx_t r[], ssz_t m, ssz_t n, ssz_t nx);
void  mad_cmat_sympconj(const cpx_t x[],                        cpx_t r[],          ssz_t n);                       //  -J M' J
num_t mad_cmat_symperr (const cpx_t x[],                        cpx_t r[],          ssz_t n, num_t *tol_);          //  M' J M - J

void  mad_imat_rev     (      idx_t x[],                                   ssz_t m, ssz_t n, int d);                //  imat -> rev(imat)
void  mad_imat_roll    (      idx_t x[],                                   ssz_t m, ssz_t n, int mroll, int nroll); //  imat ->imat
void  mad_imat_eye     (      idx_t x[],                        idx_t v,   ssz_t m, ssz_t n, ssz_t ldx);            //  imat ->imat
void  mad_imat_copy    (const idx_t x[],                        idx_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr); //  imat ->imat
void  mad_imat_copym   (const idx_t x[],                        num_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr); //  imat -> mat
void  mad_imat_trans   (const idx_t x[],                        idx_t r[], ssz_t m, ssz_t n);                       //  imat ->imat

// global fft cleanup
void  mad_fft_cleanup  (void);

// 2D/3D geometry -------------------------------------------------------------o

// rotations
void  mad_mat_rot      (      num_t x[2*2], num_t a);  // R
void  mad_mat_rotx     (      num_t x[3*3], num_t ax); // Rx
void  mad_mat_roty     (      num_t x[3*3], num_t ay); // Ry
void  mad_mat_rotz     (      num_t x[3*3], num_t az); // Rz
void  mad_mat_rotxy    (      num_t x[3*3], num_t ax, num_t ay, log_t inv); // Rx.Ry
void  mad_mat_rotxz    (      num_t x[3*3], num_t ax, num_t az, log_t inv); // Rx.Rz
void  mad_mat_rotyz    (      num_t x[3*3], num_t ay, num_t az, log_t inv); // Ry.Rz
void  mad_mat_rotxyz   (      num_t x[3*3], num_t ax, num_t ay, num_t az, log_t inv); // Rx.Ry.Rz
void  mad_mat_rotxzy   (      num_t x[3*3], num_t ax, num_t ay, num_t az, log_t inv); // Rx.Rz.Ry
void  mad_mat_rotyxz   (      num_t x[3*3], num_t ax, num_t ay, num_t az, log_t inv); // Ry.Rx.Rz
void  mad_mat_torotxyz (const num_t x[3*3], num_t r[3]                  , log_t inv); // ax, ay, az from rotxyz
void  mad_mat_torotxzy (const num_t x[3*3], num_t r[3]                  , log_t inv); // ax, ay, az from rotxzy
void  mad_mat_torotyxz (const num_t x[3*3], num_t r[3]                  , log_t inv); // ax, ay, az from rotyxz

// vector
void  mad_mat_rotv     (      num_t x[3*3], const num_t v [3], num_t a, log_t inv); // Rv
num_t mad_mat_torotv   (const num_t x[3*3],       num_t v_[3]         , log_t inv); // av from rotv

// quaternion
void  mad_mat_rotq     (      num_t x[3*3], const num_t q[4], log_t inv);
void  mad_mat_torotq   (const num_t x[3*3],       num_t q[4], log_t inv);

// misalignments
void  mad_mat_rtbar    (      num_t Rb[3*3],       num_t Tb[3], num_t el, num_t ang, num_t tlt,
                        const num_t R_[3*3], const num_t T [3]);

// special flags
extern int mad_use_madx_micado;
extern int mad_use_madx_svdcnd;

// matrix functions -----------------------------------------------------------o

struct  matrix;
struct cmatrix;
struct imatrix;

// unsafe functions, assume that matrices reshaped sizes are valid (no check!).
void mad_mat_reshape  (struct  matrix *x, ssz_t m, ssz_t n);
void mad_cmat_reshape (struct cmatrix *x, ssz_t m, ssz_t n);
void mad_imat_reshape (struct imatrix *x, ssz_t m, ssz_t n);

// ----------------------------------------------------------------------------o

#endif // MAD_MAT_H
