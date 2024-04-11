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
  - provide a full featured parametric Generalized Complex TPSA package

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

typedef struct ctpsa_ ctpsa_t;

// --- interface -------------------------------------------------------------o

// ctors, dtor
ctpsa_t* mad_ctpsa_newd    (const  desc_t *d, ord_t mo); // if mo > d_mo, mo = d_mo
ctpsa_t* mad_ctpsa_new     (const ctpsa_t *t, ord_t mo); // ok with t=(ctpsa_t*)tpsa
void     mad_ctpsa_del     (const ctpsa_t *t);

// introspection
const
desc_t*  mad_ctpsa_desc    (const ctpsa_t *t);
ssz_t    mad_ctpsa_len     (const ctpsa_t *t);
int32_t  mad_ctpsa_uid     (      ctpsa_t *t, int32_t uid_); // set uid if != 0
str_t    mad_ctpsa_nam     (      ctpsa_t *t, str_t   nam_); // set nam if != null
ord_t    mad_ctpsa_ord     (const ctpsa_t *t, log_t   hi_ ); // mo or hi
log_t    mad_ctpsa_isnul   (const ctpsa_t *t);
log_t    mad_ctpsa_isvalid (const ctpsa_t *t);

// initialization / manipulation
void     mad_ctpsa_copy    (const ctpsa_t *t, ctpsa_t *r);
void     mad_ctpsa_convert (const ctpsa_t *t, ctpsa_t *r, ssz_t n, idx_t t2r_[], int pb);
idx_t    mad_ctpsa_maxord  (const ctpsa_t *t,             ssz_t n, idx_t idx_[]);
void     mad_ctpsa_sclord  (const ctpsa_t *t, ctpsa_t *r, log_t inv, log_t prm); // t[i]*o[i]
void     mad_ctpsa_getord  (const ctpsa_t *t, ctpsa_t *r, ord_t ord);
void     mad_ctpsa_cutord  (const ctpsa_t *t, ctpsa_t *r, int   ord); // ord..mo = 0 or 0..-ord=0
void     mad_ctpsa_clrord  (      ctpsa_t *t, ord_t ord);
void     mad_ctpsa_setvar  (      ctpsa_t *t, cpx_t v, idx_t iv, cpx_t scl_);
void     mad_ctpsa_setprm  (      ctpsa_t *t, cpx_t v, idx_t ip);
void     mad_ctpsa_setval  (      ctpsa_t *t, cpx_t v);
void     mad_ctpsa_update  (      ctpsa_t *t);
void     mad_ctpsa_clear   (      ctpsa_t *t);

// initialization without complex-by-value
void     mad_ctpsa_setvar_r(      ctpsa_t *t, num_t v_re, num_t v_im, idx_t iv, num_t scl_re_, num_t scl_im_);
void     mad_ctpsa_setprm_r(      ctpsa_t *t, num_t v_re, num_t v_im, idx_t ip);
void     mad_ctpsa_setval_r(      ctpsa_t *t, num_t v_re, num_t v_im);

// real, imaginary, norm, phase, conversion
void     mad_ctpsa_cplx    (const  tpsa_t *re_, const tpsa_t *im_, ctpsa_t *r);
void     mad_ctpsa_real    (const ctpsa_t *t,  tpsa_t *r);
void     mad_ctpsa_imag    (const ctpsa_t *t,  tpsa_t *r);
void     mad_ctpsa_cabs    (const ctpsa_t *t,  tpsa_t *r);
void     mad_ctpsa_carg    (const ctpsa_t *t,  tpsa_t *r);
void     mad_ctpsa_rect    (const ctpsa_t *t, ctpsa_t *r);
void     mad_ctpsa_polar   (const ctpsa_t *t, ctpsa_t *r);

// indexing / monomials (return idx_t = -1 if invalid)
ord_t    mad_ctpsa_mono    (const ctpsa_t *t, idx_t i, ssz_t n,       ord_t m_[], ord_t *p_);
idx_t    mad_ctpsa_idxs    (const ctpsa_t *t,          ssz_t n,       str_t s   ); // string mono "[0-9]*"
idx_t    mad_ctpsa_idxm    (const ctpsa_t *t,          ssz_t n, const ord_t m []);
idx_t    mad_ctpsa_idxsm   (const ctpsa_t *t,          ssz_t n, const idx_t m []); // sparse mono [(i,o)]
idx_t    mad_ctpsa_cycle   (const ctpsa_t *t, idx_t i, ssz_t n,       ord_t m_[], cpx_t *v_);

// accessors
cpx_t    mad_ctpsa_get0    (const ctpsa_t *t);
cpx_t    mad_ctpsa_geti    (const ctpsa_t *t, idx_t i);
cpx_t    mad_ctpsa_gets    (const ctpsa_t *t, ssz_t n,       str_t s  ); // string w orders in '0'-'9'
cpx_t    mad_ctpsa_getm    (const ctpsa_t *t, ssz_t n, const ord_t m[]);
cpx_t    mad_ctpsa_getsm   (const ctpsa_t *t, ssz_t n, const idx_t m[]); // sparse mono [(i,o)]
void     mad_ctpsa_set0    (      ctpsa_t *t, /* i = 0 */               cpx_t a, cpx_t b); // a*x[0]+b
void     mad_ctpsa_seti    (      ctpsa_t *t, idx_t i,                  cpx_t a, cpx_t b); // a*x[i]+b
void     mad_ctpsa_sets    (      ctpsa_t *t, ssz_t n,       str_t s  , cpx_t a, cpx_t b); // a*x[m]+b
void     mad_ctpsa_setm    (      ctpsa_t *t, ssz_t n, const ord_t m[], cpx_t a, cpx_t b); // a*x[m]+b
void     mad_ctpsa_setsm   (      ctpsa_t *t, ssz_t n, const idx_t m[], cpx_t a, cpx_t b); // a*x[m]+b
void     mad_ctpsa_cpy0    (const ctpsa_t *t, ctpsa_t *r);
void     mad_ctpsa_cpyi    (const ctpsa_t *t, ctpsa_t *r,          idx_t i);
void     mad_ctpsa_cpys    (const ctpsa_t *t, ctpsa_t *r, ssz_t n, str_t s); // string mono "[0-9]*"
void     mad_ctpsa_cpym    (const ctpsa_t *t, ctpsa_t *r, ssz_t n, const ord_t m[]);
void     mad_ctpsa_cpysm   (const ctpsa_t *t, ctpsa_t *r, ssz_t n, const idx_t m[]); // sparse mono [(i,o)]

// accessors without complex-by-value
void     mad_ctpsa_get0_r  (const ctpsa_t *t, cpx_t *r);
void     mad_ctpsa_geti_r  (const ctpsa_t *t, idx_t i, cpx_t *r);
void     mad_ctpsa_gets_r  (const ctpsa_t *t, ssz_t n,       str_t s  , cpx_t *r);
void     mad_ctpsa_getm_r  (const ctpsa_t *t, ssz_t n, const ord_t m[], cpx_t *r);
void     mad_ctpsa_getsm_r (const ctpsa_t *t, ssz_t n, const idx_t m[], cpx_t *r);
void     mad_ctpsa_set0_r  (      ctpsa_t *t, /* i = 0 */               num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_seti_r  (      ctpsa_t *t, idx_t i,                  num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_sets_r  (      ctpsa_t *t, ssz_t n,       str_t s  , num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_setm_r  (      ctpsa_t *t, ssz_t n, const ord_t m[], num_t a_re, num_t a_im, num_t b_re, num_t b_im);
void     mad_ctpsa_setsm_r (      ctpsa_t *t, ssz_t n, const idx_t m[], num_t a_re, num_t a_im, num_t b_re, num_t b_im);

// accessors vector based
void     mad_ctpsa_getv    (const ctpsa_t *t, idx_t i, ssz_t n,       cpx_t v[]); // return copied length
void     mad_ctpsa_setv    (      ctpsa_t *t, idx_t i, ssz_t n, const cpx_t v[]); // return copied length

// operators
log_t    mad_ctpsa_equ     (const ctpsa_t *a, const ctpsa_t *b, num_t tol_);
void     mad_ctpsa_dif     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c); // (a_i-b_i)/max(|a_i|,1)
void     mad_ctpsa_add     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_sub     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_mul     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_div     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_pow     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_powi    (const ctpsa_t *a,       int      n, ctpsa_t *c);
void     mad_ctpsa_pown    (const ctpsa_t *a,       cpx_t    v, ctpsa_t *c);

// operators without complex-by-value
void     mad_ctpsa_pown_r  (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);

// operators with internal real-to-complex conversion
log_t    mad_ctpsa_equt    (const ctpsa_t *a, const  tpsa_t *b, num_t tol);
void     mad_ctpsa_dift    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_tdif    (const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_addt    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_subt    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_tsub    (const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_mult    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_divt    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_tdiv    (const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_powt    (const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_tpow    (const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c);

// functions
num_t    mad_ctpsa_nrm     (const ctpsa_t *a);
void     mad_ctpsa_unit    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_conj    (const ctpsa_t *a, ctpsa_t *c);
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
void     mad_ctpsa_asinc   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_asinh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acosh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_atanh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acoth   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_asinhc  (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_erf     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_erfc    (const ctpsa_t *a, ctpsa_t *c);

void     mad_ctpsa_acc     (const ctpsa_t *a, cpx_t v, ctpsa_t *c); // c += v*a, aliasing OK
void     mad_ctpsa_scl     (const ctpsa_t *a, cpx_t v, ctpsa_t *c); // c  = v*a
void     mad_ctpsa_inv     (const ctpsa_t *a, cpx_t v, ctpsa_t *c); // c  = v/a
void     mad_ctpsa_invsqrt (const ctpsa_t *a, cpx_t v, ctpsa_t *c); // c  = v/sqrt(a)

void     mad_ctpsa_hypot   (const ctpsa_t *x, const ctpsa_t *y, ctpsa_t *r);
void     mad_ctpsa_hypot3  (const ctpsa_t *x, const ctpsa_t *y, const ctpsa_t *z, ctpsa_t *r);

// functions without complex-by-value arguments
void     mad_ctpsa_acc_r    (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_scl_r    (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_inv_r    (const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);
void     mad_ctpsa_invsqrt_r(const ctpsa_t *a, num_t v_re, num_t v_im, ctpsa_t *c);

// functions for differential algebra
void     mad_ctpsa_integ   (const ctpsa_t *a, ctpsa_t *c, idx_t iv);
void     mad_ctpsa_deriv   (const ctpsa_t *a, ctpsa_t *c, idx_t iv);
void     mad_ctpsa_derivm  (const ctpsa_t *a, ctpsa_t *c, ssz_t n, const ord_t m[]);
void     mad_ctpsa_poisbra (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int nv);
void     mad_ctpsa_taylor  (const ctpsa_t *a, ssz_t n, const cpx_t coef[], ctpsa_t *c);

// functions for differential algebra with internal real-to-complex conversion
void     mad_ctpsa_poisbrat(const ctpsa_t *a, const  tpsa_t *b, ctpsa_t *c, int nv);
void     mad_ctpsa_tpoisbra(const  tpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int nv);

// high level functions
void     mad_ctpsa_axpb        (cpx_t a, const ctpsa_t *x,
                                cpx_t b, ctpsa_t *r);
void     mad_ctpsa_axpbypc     (cpx_t a, const ctpsa_t *x,
                                cpx_t b, const ctpsa_t *y,
                                cpx_t c, ctpsa_t *r);
void     mad_ctpsa_axypb       (cpx_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cpx_t b, ctpsa_t *r);
void     mad_ctpsa_axypbzpc    (cpx_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cpx_t b, const ctpsa_t *z,
                                cpx_t c, ctpsa_t *r);
void     mad_ctpsa_axypbvwpc   (cpx_t a, const ctpsa_t *x, const ctpsa_t *y,
                                cpx_t b, const ctpsa_t *v, const ctpsa_t *w,
                                cpx_t c, ctpsa_t *r);
void     mad_ctpsa_ax2pby2pcz2 (cpx_t a, const ctpsa_t *x,
                                cpx_t b, const ctpsa_t *y,
                                cpx_t c, const ctpsa_t *z, ctpsa_t *r);

void     mad_ctpsa_axpsqrtbpcx2    (const ctpsa_t *x, cpx_t a, cpx_t b, cpx_t c, ctpsa_t *r);
void     mad_ctpsa_logaxpsqrtbpcx2 (const ctpsa_t *x, cpx_t a, cpx_t b, cpx_t c, ctpsa_t *r);
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

// map functions (to check for non-homogeneous maps & parameters)
void     mad_ctpsa_vec2fld  (ssz_t na, const ctpsa_t *a   ,                      ctpsa_t *mc[]); // F . grad
void     mad_ctpsa_fld2vec  (ssz_t na, const ctpsa_t *ma[],                      ctpsa_t *c   );
void     mad_ctpsa_fgrad    (ssz_t na, const ctpsa_t *ma[], const ctpsa_t * b  , ctpsa_t *c   );
void     mad_ctpsa_liebra   (ssz_t na, const ctpsa_t *ma[], const ctpsa_t *mb[], ctpsa_t *mc[]);
void     mad_ctpsa_exppb    (ssz_t na, const ctpsa_t *ma[], const ctpsa_t *mb[], ctpsa_t *mc[]); // exp(:F:) K
void     mad_ctpsa_logpb    (ssz_t na, const ctpsa_t *ma[], const ctpsa_t *mb[], ctpsa_t *mc[]); // exp(log(:F:)) K

ord_t    mad_ctpsa_mord     (ssz_t na, const ctpsa_t *ma[], log_t hi); // max mo (or max hi)
num_t    mad_ctpsa_mnrm     (ssz_t na, const ctpsa_t *ma[]);
void     mad_ctpsa_minv     (ssz_t na, const ctpsa_t *ma[], ssz_t nb,                      ctpsa_t *mc[]);
void     mad_ctpsa_pminv    (ssz_t na, const ctpsa_t *ma[], ssz_t nb,                      ctpsa_t *mc[], idx_t select[]);
void     mad_ctpsa_compose  (ssz_t na, const ctpsa_t *ma[], ssz_t nb, const ctpsa_t *mb[], ctpsa_t *mc[]);
void     mad_ctpsa_translate(ssz_t na, const ctpsa_t *ma[], ssz_t nb, const cpx_t    tb[], ctpsa_t *mc[]);
void     mad_ctpsa_eval     (ssz_t na, const ctpsa_t *ma[], ssz_t nb, const cpx_t    tb[], cpx_t    tc[]);
void     mad_ctpsa_mconv    (ssz_t na, const ctpsa_t *ma[], ssz_t nc,                      ctpsa_t *mc[], ssz_t n, idx_t t2r_[], int pb);

// I/O
void     mad_ctpsa_print    (const ctpsa_t *t, str_t name_, num_t eps_, int nohdr_, FILE *stream_);
ctpsa_t* mad_ctpsa_scan     (                                                       FILE *stream_);
const
desc_t*  mad_ctpsa_scan_hdr (      int *kind_, char  name_[NAMSZ],                  FILE *stream_);
void     mad_ctpsa_scan_coef(      ctpsa_t *t,                                      FILE *stream_);

// unsafe operation (mo vs allocated!!)
ctpsa_t* mad_ctpsa_init     (      ctpsa_t *t, const desc_t *d, ord_t mo);

// debug
void     mad_ctpsa_debug    (const ctpsa_t *t, str_t name_, str_t fnam_, int line_, FILE *stream_);

// ---------------------------------------------------------------------------o

#endif // MAD_CTPSA_H
