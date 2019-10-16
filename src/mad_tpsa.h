#ifndef MAD_TPSA_H
#define MAD_TPSA_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module interface
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

  Purpose:
  - provide a full feathered Generalized TPSA package

  Information:
  - parameters ending with an underscope are optional (i.e. can be null).

  Errors:
  - TODO

 o-----------------------------------------------------------------------------o
 */

#include <stdio.h>

#include "mad_mono.h"
#include "mad_desc.h"

// --- types ------------------------------------------------------------------o

typedef struct tpsa tpsa_t;

// --- globals ----------------------------------------------------------------o

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;

// --- interface --------------------------------------------------------------o

// ctors, dtor, shape
tpsa_t* mad_tpsa_newd    (const desc_t *d, ord_t mo); // if mo > d_mo, mo = d_mo
tpsa_t* mad_tpsa_new     (const tpsa_t *t, ord_t mo); // ok with t=(tpsa_t*)ctpsa
void    mad_tpsa_del     (const tpsa_t *t);

// introspection
const
desc_t* mad_tpsa_desc    (const tpsa_t *t);
int32_t mad_tpsa_uid     (      tpsa_t *t, int32_t uid_); // set uid if != 0
ssz_t   mad_tpsa_len     (const tpsa_t *t);
ord_t   mad_tpsa_ord     (const tpsa_t *t);
ord_t   mad_tpsa_ordv    (const tpsa_t *t, ...);        // max order of all
ord_t   mad_tpsa_ordn    (ssz_t n, const tpsa_t *t[n]); // max order of all

// initialization
void    mad_tpsa_copy    (const tpsa_t *t, tpsa_t *r);
void    mad_tpsa_getord  (const tpsa_t *t, tpsa_t *r, ord_t ord);
void    mad_tpsa_cutord  (const tpsa_t *t, tpsa_t *r, int   ord); // ord..mo = 0 or 0..-ord=0
void    mad_tpsa_convert (const tpsa_t *t, tpsa_t *r, ssz_t n, idx_t t2r_[n]);
void    mad_tpsa_setvar  (      tpsa_t *t, num_t v, idx_t iv_, num_t scl_);
void    mad_tpsa_clear   (      tpsa_t *t);

// indexing / monomials (return idx_t = -1 if invalid)
ord_t   mad_tpsa_mono    (const tpsa_t *t, ssz_t n,       ord_t m_[n], idx_t i);
idx_t   mad_tpsa_idxs    (const tpsa_t *t, ssz_t n,       str_t s    ); // string mono "[0-9]*"
idx_t   mad_tpsa_idxm    (const tpsa_t *t, ssz_t n, const ord_t m [n]);
idx_t   mad_tpsa_idxsm   (const tpsa_t *t, ssz_t n, const int   m [n]); // sparse mono [(i,o)]
idx_t   mad_tpsa_cycle   (const tpsa_t *t, ssz_t n,       ord_t m_[n], idx_t i, num_t *v_);

// accessors
num_t   mad_tpsa_get0    (const tpsa_t *t);
num_t   mad_tpsa_geti    (const tpsa_t *t, idx_t i);
num_t   mad_tpsa_gets    (const tpsa_t *t, ssz_t n,       str_t s   ); // string mono "[0-9]*"
num_t   mad_tpsa_getm    (const tpsa_t *t, ssz_t n, const ord_t m[n]);
num_t   mad_tpsa_getsm   (const tpsa_t *t, ssz_t n, const int   m[n]); // sparse mono [(i,o)]
void    mad_tpsa_set0    (      tpsa_t *t, /* i = 0 */                   num_t a, num_t b);
void    mad_tpsa_seti    (      tpsa_t *t, idx_t i,                      num_t a, num_t b);
void    mad_tpsa_sets    (      tpsa_t *t, ssz_t n,       str_t s   ,    num_t a, num_t b);
void    mad_tpsa_setm    (      tpsa_t *t, ssz_t n, const ord_t m[n],    num_t a, num_t b);
void    mad_tpsa_setsm   (      tpsa_t *t, ssz_t n, const int   m[n],    num_t a, num_t b);

// accessors vector based
void    mad_tpsa_getv    (const tpsa_t *t, idx_t i, ssz_t n,             num_t v[n]);
void    mad_tpsa_setv    (      tpsa_t *t, idx_t i, ssz_t n,    const    num_t v[n]);

// operators
log_t   mad_tpsa_equ     (const tpsa_t *a, const tpsa_t *b, num_t eps_);
void    mad_tpsa_add     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_sub     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_mul     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_div     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_pow     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_powi    (const tpsa_t *a, int           n, tpsa_t *c);
void    mad_tpsa_pown    (const tpsa_t *a, num_t         v, tpsa_t *c);

// functions
void    mad_tpsa_abs     (const tpsa_t *a, tpsa_t *c);
num_t   mad_tpsa_nrm1    (const tpsa_t *a, const tpsa_t *b_);
num_t   mad_tpsa_nrm2    (const tpsa_t *a, const tpsa_t *b_);
void    mad_tpsa_deriv   (const tpsa_t *a, tpsa_t *c, int iv);
void    mad_tpsa_derivm  (const tpsa_t *a, tpsa_t *c, ssz_t n, const ord_t m[n]);
void    mad_tpsa_poisson (const tpsa_t *a, const tpsa_t *b, tpsa_t *c, int nv);
void    mad_tpsa_taylor  (const tpsa_t *a, ssz_t n, const num_t coef[n], tpsa_t *c);

void    mad_tpsa_acc     (const tpsa_t *a, num_t v, tpsa_t *c);  // c += v*a, aliasing OK
void    mad_tpsa_scl     (const tpsa_t *a, num_t v, tpsa_t *c);  // c  = v*a
void    mad_tpsa_inv     (const tpsa_t *a, num_t v, tpsa_t *c);  // c  = v/a
void    mad_tpsa_invsqrt (const tpsa_t *a, num_t v, tpsa_t *c);  // c  = v/sqrt(a)

void    mad_tpsa_sqrt    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_exp     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_log     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_sincos  (const tpsa_t *a, tpsa_t *s, tpsa_t *c);
void    mad_tpsa_sin     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_cos     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_tan     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_cot     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_sinc    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_sincosh (const tpsa_t *a, tpsa_t *s, tpsa_t *c);
void    mad_tpsa_sinh    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_cosh    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_tanh    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_coth    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_sinhc   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_asin    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acos    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_atan    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acot    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_asinh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acosh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_atanh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acoth   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_erf     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_erfc    (const tpsa_t *a, tpsa_t *c);

// high level functions (aliasing OK)
void    mad_tpsa_axpb       (num_t a, const tpsa_t *x,
                             num_t b, tpsa_t *r);
void    mad_tpsa_axpbypc    (num_t a, const tpsa_t *x,
                             num_t b, const tpsa_t *y,
                             num_t c, tpsa_t *r);
void    mad_tpsa_axypb      (num_t a, const tpsa_t *x, const tpsa_t *y,
                             num_t b, tpsa_t *r);
void    mad_tpsa_axypbzpc   (num_t a, const tpsa_t *x, const tpsa_t *y,
                             num_t b, const tpsa_t *z,
                             num_t c, tpsa_t *r);
void    mad_tpsa_axypbvwpc  (num_t a, const tpsa_t *x, const tpsa_t *y,
                             num_t b, const tpsa_t *v, const tpsa_t *w,
                             num_t c, tpsa_t *r);
void    mad_tpsa_ax2pby2pcz2(num_t a, const tpsa_t *x,
                             num_t b, const tpsa_t *y,
                             num_t c, const tpsa_t *z, tpsa_t *r);

void    mad_tpsa_axpsqrtbpcx2    (const tpsa_t *x, num_t a, num_t b, num_t c, tpsa_t *r);
void    mad_tpsa_logaxpsqrtbpcx2 (const tpsa_t *x, num_t a, num_t b, num_t c, tpsa_t *r);
void    mad_tpsa_logxdy          (const tpsa_t *x, const tpsa_t *y, tpsa_t *r);

// to check for non-homogeneous maps & knobs
void    mad_tpsa_minv     (ssz_t n , const tpsa_t *ma[n ],                                 tpsa_t *mc[n ]);
void    mad_tpsa_pminv    (ssz_t n , const tpsa_t *ma[n ],                                 tpsa_t *mc[n ], idx_t select[n]);
void    mad_tpsa_compose  (ssz_t na, const tpsa_t *ma[na], ssz_t nb, const tpsa_t *mb[nb], tpsa_t *mc[na]);
void    mad_tpsa_translate(ssz_t na, const tpsa_t *ma[na], ssz_t nb, const num_t   tb[nb], tpsa_t *mc[na]);
void    mad_tpsa_eval     (ssz_t na, const tpsa_t *ma[na], ssz_t nb, const num_t   tb[nb], num_t   tc[nb]);

// I/O
void    mad_tpsa_print    (const tpsa_t *t, str_t name_, num_t eps_, int nohdr_, FILE *stream_);
tpsa_t* mad_tpsa_scan     (                                                      FILE *stream_);
const
desc_t* mad_tpsa_scan_hdr (                 int  *kind_,                         FILE *stream_);
void    mad_tpsa_scan_coef(      tpsa_t *t,                                      FILE *stream_);
void    mad_tpsa_debug    (const tpsa_t *t, str_t name_, str_t fnam_, int line_, FILE *stream_);
log_t   mad_tpsa_isvalid  (const tpsa_t *t);

// unsafe operation (mo vs allocated!!)
tpsa_t* mad_tpsa_init (tpsa_t *t, const desc_t *d, ord_t mo);

// macro wrapper for safe use
#define mad_tpsa_ordv(...) mad_tpsa_ordv(__VA_ARGS__,NULL)

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_H
