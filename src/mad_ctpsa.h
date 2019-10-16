#ifndef MAD_CTPSA_H
#define MAD_CTPSA_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Complex Truncated Power Series Algebra module interface
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
  - provide a full feathered Generalized Complex TPSA package

  Information:
  - parameters ending with an underscope can be null.

  Errors:
  - TODO

 o-----------------------------------------------------------------------------o
 */

#include <stdio.h>

#include "mad_mono.h"
#include "mad_desc.h"
#include "mad_tpsa.h"

// --- types -----------------------------------------------------------------o

typedef struct ctpsa ctpsa_t;

// --- globals ---------------------------------------------------------------o

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;

// --- interface -------------------------------------------------------------o

// ctors, dtor
ctpsa_t* mad_ctpsa_newd    (const  desc_t *d, ord_t mo); // if mo > d_mo, mo = d_mo
ctpsa_t* mad_ctpsa_new     (const ctpsa_t *t, ord_t mo); // ok with t=(ctpsa_t*)tpsa
void     mad_ctpsa_del     (const ctpsa_t *t);

// introspection
const
desc_t*  mad_ctpsa_desc    (const ctpsa_t *t);
int32_t  mad_ctpsa_uid     (      ctpsa_t *t, int32_t uid_); // set uid if != 0
ssz_t    mad_ctpsa_len     (const ctpsa_t *t);
ord_t    mad_ctpsa_ord     (const ctpsa_t *t);
ord_t    mad_ctpsa_ordv    (const ctpsa_t *t, ...);        // max order of all
ord_t    mad_ctpsa_ordn    (ssz_t n, const ctpsa_t *t[n]); // max order of all

// initialization
void     mad_ctpsa_copy    (const ctpsa_t *t, ctpsa_t *r);
void     mad_ctpsa_getord  (const ctpsa_t *t, ctpsa_t *r, ord_t ord);
void     mad_ctpsa_cutord  (const ctpsa_t *t, ctpsa_t *r, int   ord); // ord..mo = 0 or 0..-ord=0
void     mad_ctpsa_convert (const ctpsa_t *t, ctpsa_t *r, ssz_t n, idx_t t2r_[n]);
void     mad_ctpsa_setvar  (      ctpsa_t *t, cnum_t v, idx_t iv_, cnum_t scl_);
void     mad_ctpsa_setvar_r(      ctpsa_t *t, num_t v_re, num_t v_im, idx_t iv_, num_t scl_re_, num_t scl_im_);
void     mad_ctpsa_clear   (      ctpsa_t *t);

// conversion
void     mad_ctpsa_real    (const ctpsa_t *t, tpsa_t *r);
void     mad_ctpsa_imag    (const ctpsa_t *t, tpsa_t *r);
void     mad_ctpsa_complex (const  tpsa_t *re_, const tpsa_t *im_, ctpsa_t *r);

// indexing / monomials (return idx_t = -1 if invalid)
ord_t    mad_ctpsa_mono    (const ctpsa_t *t, ssz_t n,       ord_t m_[n], idx_t i);
idx_t    mad_ctpsa_idxs    (const ctpsa_t *t, ssz_t n,       str_t s    ); // string mono "[0-9]*"
idx_t    mad_ctpsa_idxm    (const ctpsa_t *t, ssz_t n, const ord_t m [n]);
idx_t    mad_ctpsa_idxsm   (const ctpsa_t *t, ssz_t n, const int   m [n]); // sparse mono [(i,o)]
idx_t    mad_ctpsa_cycle   (const ctpsa_t *t, ssz_t n,       ord_t m_[n], idx_t i, num_t *v_);

// accessors
cnum_t   mad_ctpsa_get0    (const ctpsa_t *t);
cnum_t   mad_ctpsa_geti    (const ctpsa_t *t, idx_t i);
cnum_t   mad_ctpsa_gets    (const ctpsa_t *t, ssz_t n,       str_t s   ); // string w orders in '0'-'9'
cnum_t   mad_ctpsa_getm    (const ctpsa_t *t, ssz_t n, const ord_t m[n]);
cnum_t   mad_ctpsa_getsm   (const ctpsa_t *t, ssz_t n, const int   m[n]); // sparse mono [(i,o)]
void     mad_ctpsa_set0    (      ctpsa_t *t, /* i = 0 */                cnum_t a, cnum_t b); // a*x[0]+b
void     mad_ctpsa_seti    (      ctpsa_t *t, idx_t i,                   cnum_t a, cnum_t b); // a*x[i]+b
void     mad_ctpsa_sets    (      ctpsa_t *t, ssz_t n,       str_t s   , cnum_t a, cnum_t b); // a*x[m]+b
void     mad_ctpsa_setm    (      ctpsa_t *t, ssz_t n, const ord_t m[n], cnum_t a, cnum_t b); // a*x[m]+b
void     mad_ctpsa_setsm   (      ctpsa_t *t, ssz_t n, const int   m[n], cnum_t a, cnum_t b); // a*x[m]+b

// accessors without complex-by-value
void     mad_ctpsa_get0_r  (const ctpsa_t *t, cnum_t *r);
void     mad_ctpsa_geti_r  (const ctpsa_t *t, idx_t i, cnum_t *r);
void     mad_ctpsa_gets_r  (const ctpsa_t *t, ssz_t n,       str_t s   , cnum_t *r);
void     mad_ctpsa_getm_r  (const ctpsa_t *t, ssz_t n, const ord_t m[n], cnum_t *r);
void     mad_ctpsa_getsm_r (const ctpsa_t *t, ssz_t n, const int   m[n], cnum_t *r);
void     mad_ctpsa_set0_r  (      ctpsa_t *t, /* i = 0 */                num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_seti_r  (      ctpsa_t *t, idx_t i,                   num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_sets_r  (      ctpsa_t *t, ssz_t n,       str_t s   , num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_setm_r  (      ctpsa_t *t, ssz_t n, const ord_t m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_setsm_r (      ctpsa_t *t, ssz_t n, const int   m[n], num_t a_re, num_t a_im, num_t b_re, num_t b_im);

// accessors vector based
void     mad_ctpsa_getv    (const ctpsa_t *t, idx_t i, ssz_t n,          cnum_t v[n]);
void     mad_ctpsa_setv    (      ctpsa_t *t, idx_t i, ssz_t n,    const cnum_t v[n]);

// operators
log_t    mad_ctpsa_equ     (const ctpsa_t *a, const ctpsa_t *b, num_t eps_);
void     mad_ctpsa_add     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_sub     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_mul     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_div     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_pow     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_powi    (const ctpsa_t *a, int            n, ctpsa_t *c);
void     mad_ctpsa_pown    (const ctpsa_t *a, cnum_t         v, ctpsa_t *c);

// operators without complex-by-value arguments
void     mad_ctpsa_pown_r    (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);

// operators with internal real-to-complex conversion
log_t    mad_ctpsa_equt    (const ctpsa_t *a, const  tpsa_t *b, num_t tol);
void     mad_ctpsa_addt    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_subt    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_tsub    (const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_mult    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_divt    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_tdiv    (const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_powt    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_tpow    (const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_poisst  (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c, int nv);
void     mad_ctpsa_tpoiss  (const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int nv);

// functions
void     mad_ctpsa_abs     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_arg     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_conj    (const ctpsa_t *a, ctpsa_t *c);
cnum_t   mad_ctpsa_nrm1    (const ctpsa_t *a, const ctpsa_t *b_);
cnum_t   mad_ctpsa_nrm2    (const ctpsa_t *a, const ctpsa_t *b_);
void     mad_ctpsa_deriv   (const ctpsa_t *a, ctpsa_t *c, int iv); // TODO: check functions that rely on it
void     mad_ctpsa_derivm  (const ctpsa_t *a, ctpsa_t *c, ssz_t n, const ord_t m[n]);
void     mad_ctpsa_poisson (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int nv);
void     mad_ctpsa_taylor  (const ctpsa_t *a, ssz_t n, const cnum_t coef[n], ctpsa_t *c);

void     mad_ctpsa_acc     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c += v*a, aliasing OK
void     mad_ctpsa_scl     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v*a
void     mad_ctpsa_inv     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v/a
void     mad_ctpsa_invsqrt (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v/sqrt(a)

void     mad_ctpsa_sqrt    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_exp     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_log     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sincos  (const ctpsa_t *a, ctpsa_t *s, ctpsa_t *c);
void     mad_ctpsa_sin     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cos     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_tan     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cot     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sinc    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sincosh (const ctpsa_t *a, ctpsa_t *s, ctpsa_t *c);
void     mad_ctpsa_sinh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cosh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_tanh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_coth    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sinhc   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_asin    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acos    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_atan    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acot    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_asinh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acosh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_atanh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acoth   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_erf     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_erfc    (const ctpsa_t *a, ctpsa_t *c);

// functions without complex-by-value arguments
void     mad_ctpsa_nrm1_r    (const ctpsa_t *a, const ctpsa_t *b_, cnum_t *r);
void     mad_ctpsa_nrm2_r    (const ctpsa_t *a, const ctpsa_t *b_, cnum_t *r);
void     mad_ctpsa_acc_r     (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_scl_r     (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_inv_r     (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_invsqrt_r (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);

// high level functions (aliasing OK)
void     mad_ctpsa_axpb        (cnum_t a, const ctpsa_t *x,
                                cnum_t b, ctpsa_t *r);
void     mad_ctpsa_axpbypc     (cnum_t a, const ctpsa_t *x,
                                cnum_t b, const ctpsa_t *y,
                                cnum_t c, ctpsa_t *r);
void     mad_ctpsa_axypb       (cnum_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cnum_t b, ctpsa_t *r);
void     mad_ctpsa_axypbzpc    (cnum_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cnum_t b, const ctpsa_t *z,
                                cnum_t c, ctpsa_t *r);
void     mad_ctpsa_axypbvwpc   (cnum_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cnum_t b, const ctpsa_t *v, const ctpsa_t *w,
                                cnum_t c, ctpsa_t *r);
void     mad_ctpsa_ax2pby2pcz2 (cnum_t a, const ctpsa_t *x,
                                cnum_t b, const ctpsa_t *y,
                                cnum_t c, const ctpsa_t *z, ctpsa_t *r);

void     mad_ctpsa_axpsqrtbpcx2    (const ctpsa_t *x, cnum_t a, cnum_t b, cnum_t c, ctpsa_t *r);
void     mad_ctpsa_logaxpsqrtbpcx2 (const ctpsa_t *x, cnum_t a, cnum_t b, cnum_t c, ctpsa_t *r);
void     mad_ctpsa_logxdy          (const ctpsa_t *x, const ctpsa_t *y, ctpsa_t *r);

// high level functions without complex-by-value
void     mad_ctpsa_axpb_r        (num_t a_re, num_t a_im, const ctpsa_t *x,
                                  num_t b_re, num_t b_im, ctpsa_t *r);
void     mad_ctpsa_axpbypc_r     (num_t a_re, num_t a_im, const ctpsa_t *x,
                                  num_t b_re, num_t b_im, const ctpsa_t *y,
                                  num_t c_re, num_t c_im, ctpsa_t *r);
void     mad_ctpsa_axypb_r       (num_t a_re, num_t a_im, const ctpsa_t *x, const ctpsa_t *y,
                                  num_t b_re, num_t b_im, ctpsa_t *r);
void     mad_ctpsa_axypbzpc_r    (num_t a_re, num_t a_im, const ctpsa_t *x, const ctpsa_t *y,
                                  num_t b_re, num_t b_im, const ctpsa_t *z,
                                  num_t c_re, num_t c_im, ctpsa_t *r);
void     mad_ctpsa_axypbvwpc_r   (num_t a_re, num_t a_im, const ctpsa_t *x, const ctpsa_t *y,
                                  num_t b_re, num_t b_im, const ctpsa_t *v, const ctpsa_t *w,
                                  num_t c_re, num_t c_im, ctpsa_t *r);
void     mad_ctpsa_ax2pby2pcz2_r (num_t a_re, num_t a_im, const ctpsa_t *x,
                                  num_t b_re, num_t b_im, const ctpsa_t *y,
                                  num_t c_re, num_t c_im, const ctpsa_t *z, ctpsa_t *r);

void     mad_ctpsa_axpsqrtbpcx2_r    (const ctpsa_t *x, num_t a_re, num_t a_im,
                                                        num_t b_re, num_t b_im,
                                                        num_t c_re, num_t c_im, ctpsa_t *r);
void     mad_ctpsa_logaxpsqrtbpcx2_r (const ctpsa_t *x, num_t a_re, num_t a_im,
                                                        num_t b_re, num_t b_im,
                                                        num_t c_re, num_t c_im, ctpsa_t *r);

// to check for non-homogeneous maps & knobs
void     mad_ctpsa_minv     (ssz_t na, const ctpsa_t *ma[na],                                  ctpsa_t *mc[na]);
void     mad_ctpsa_pminv    (ssz_t na, const ctpsa_t *ma[na],                                  ctpsa_t *mc[na], idx_t select[na]);
void     mad_ctpsa_compose  (ssz_t na, const ctpsa_t *ma[na], ssz_t nb, const ctpsa_t *mb[nb], ctpsa_t *mc[na]);
void     mad_ctpsa_translate(ssz_t na, const ctpsa_t *ma[na], ssz_t nb, const cnum_t   tb[nb], ctpsa_t *mc[na]);
void     mad_ctpsa_eval     (ssz_t na, const ctpsa_t *ma[na], ssz_t nb, const cnum_t   tb[nb], cnum_t   tc[nb]);

// I/O
void     mad_ctpsa_print    (const ctpsa_t *t, str_t name_, num_t eps_, int nohdr_, FILE *stream_);
ctpsa_t* mad_ctpsa_scan     (                                                       FILE *stream_);
const
desc_t*  mad_ctpsa_scan_hdr (                  int  *kind_,                         FILE *stream_);
void     mad_ctpsa_scan_coef(      ctpsa_t *t,                                      FILE *stream_);
void     mad_ctpsa_debug    (const ctpsa_t *t, str_t name_, str_t fnam_, int line_, FILE *stream_);
log_t    mad_ctpsa_isvalid  (const ctpsa_t *t);

// unsafe operation (mo vs allocated!!)
ctpsa_t* mad_ctpsa_init (ctpsa_t *t, const desc_t *d, ord_t mo);

// macro wrapper for safe use
#define  mad_ctpsa_ordv(...) mad_ctpsa_ordv(__VA_ARGS__,NULL)

// ---------------------------------------------------------------------------o

#endif // MAD_CTPSA_H
