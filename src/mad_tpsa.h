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
  - provide a full featured parametric Generalized TPSA package

  Information:
  - parameters ending with an underscope are optional (i.e. can be null).

  Errors:
  - TODO

 o-----------------------------------------------------------------------------o
 */

#include <stdio.h>

#include "mad_mono.h"
#include "mad_desc.h"

// --- constants --------------------------------------------------------------o

enum { NAMSZ=16 };

// --- globals ----------------------------------------------------------------o

extern const ord_t mad_tpsa_dflt;
extern const ord_t mad_tpsa_same;

// --- types ------------------------------------------------------------------o

typedef struct tpsa_ tpsa_t;

// --- interface --------------------------------------------------------------o

// ctors, dtor, shape
tpsa_t* mad_tpsa_newd    (const desc_t *d, ord_t mo); // if mo > d_mo, mo = d_mo
tpsa_t* mad_tpsa_new     (const tpsa_t *t, ord_t mo); // ok with t=(tpsa_t*)ctpsa
void    mad_tpsa_del     (const tpsa_t *t);

// introspection
const
desc_t* mad_tpsa_desc    (const tpsa_t *t);
ord_t   mad_tpsa_mo      (      tpsa_t *t, ord_t   mo  ); // set mo
int32_t mad_tpsa_uid     (      tpsa_t *t, int32_t uid_); // set uid if != 0
str_t   mad_tpsa_nam     (      tpsa_t *t, str_t   nam_); // set nam if != null
ssz_t   mad_tpsa_len     (const tpsa_t *t, log_t   hi_ ); // get mo or hi
ord_t   mad_tpsa_ord     (const tpsa_t *t, log_t   hi_ ); // get mo or hi
log_t   mad_tpsa_isnul   (const tpsa_t *t);
log_t   mad_tpsa_isval   (const tpsa_t *t);
log_t   mad_tpsa_isvalid (const tpsa_t *t);
num_t   mad_tpsa_density (const tpsa_t *t); // ratio nz/nc in [lo,hi]

// initialization / manipulation
void    mad_tpsa_copy    (const tpsa_t *t, tpsa_t *r);
void    mad_tpsa_convert (const tpsa_t *t, tpsa_t *r, ssz_t n, idx_t t2r_[], int pb);
idx_t   mad_tpsa_maxord  (const tpsa_t *t,            ssz_t n, idx_t idx_[]);
void    mad_tpsa_sclord  (const tpsa_t *t, tpsa_t *r, log_t inv, log_t prm); // t[i]*o[i]
void    mad_tpsa_getord  (const tpsa_t *t, tpsa_t *r, ord_t ord);
void    mad_tpsa_cutord  (const tpsa_t *t, tpsa_t *r, int   ord); // ord..mo = 0 or 0..-ord=0
void    mad_tpsa_clrord  (      tpsa_t *t, ord_t ord);
void    mad_tpsa_setvar  (      tpsa_t *t, num_t v, idx_t iv, num_t scl_);
void    mad_tpsa_setprm  (      tpsa_t *t, num_t v, idx_t ip);
void    mad_tpsa_setval  (      tpsa_t *t, num_t v);
void    mad_tpsa_update  (      tpsa_t *t);
void    mad_tpsa_clear   (      tpsa_t *t);

// indexing / monomials (return idx_t = -1 if invalid)
ord_t   mad_tpsa_mono    (const tpsa_t *t, idx_t i, ssz_t n,       ord_t m_[], ord_t *p_);
idx_t   mad_tpsa_idxs    (const tpsa_t *t,          ssz_t n,       str_t s   ); // string mono "[0-9]*"
idx_t   mad_tpsa_idxm    (const tpsa_t *t,          ssz_t n, const ord_t m []);
idx_t   mad_tpsa_idxsm   (const tpsa_t *t,          ssz_t n, const idx_t m []); // sparse mono [(i,o)]
idx_t   mad_tpsa_cycle   (const tpsa_t *t, idx_t i, ssz_t n,       ord_t m_[], num_t *v_);

// accessors
num_t   mad_tpsa_geti    (const tpsa_t *t, idx_t i);
num_t   mad_tpsa_gets    (const tpsa_t *t, ssz_t n,       str_t s  ); // string mono "[0-9]*"
num_t   mad_tpsa_getm    (const tpsa_t *t, ssz_t n, const ord_t m[]);
num_t   mad_tpsa_getsm   (const tpsa_t *t, ssz_t n, const idx_t m[]); // sparse mono [(i,o)]
void    mad_tpsa_seti    (      tpsa_t *t, idx_t i,                  num_t a, num_t b);
void    mad_tpsa_sets    (      tpsa_t *t, ssz_t n,       str_t s  , num_t a, num_t b);
void    mad_tpsa_setm    (      tpsa_t *t, ssz_t n, const ord_t m[], num_t a, num_t b);
void    mad_tpsa_setsm   (      tpsa_t *t, ssz_t n, const idx_t m[], num_t a, num_t b);
void    mad_tpsa_cpyi    (const tpsa_t *t, tpsa_t *r,          idx_t i);
void    mad_tpsa_cpys    (const tpsa_t *t, tpsa_t *r, ssz_t n, str_t s); // string mono "[0-9]*"
void    mad_tpsa_cpym    (const tpsa_t *t, tpsa_t *r, ssz_t n, const ord_t m[]);
void    mad_tpsa_cpysm   (const tpsa_t *t, tpsa_t *r, ssz_t n, const idx_t m[]); // sparse mono [(i,o)]

// accessors vector based
void    mad_tpsa_getv    (const tpsa_t *t, idx_t i, ssz_t n,       num_t v[]); // return copied length
void    mad_tpsa_setv    (      tpsa_t *t, idx_t i, ssz_t n, const num_t v[]); // return copied length

// operators
log_t   mad_tpsa_equ     (const tpsa_t *a, const tpsa_t *b, num_t tol_);
void    mad_tpsa_dif     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c); // (a_i-b_i)/max(|a_i|,1)
void    mad_tpsa_add     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_sub     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_mul     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_div     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_pow     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_powi    (const tpsa_t *a, int           n, tpsa_t *c);
void    mad_tpsa_pown    (const tpsa_t *a, num_t         v, tpsa_t *c);

// functions
num_t   mad_tpsa_nrm     (const tpsa_t *a);
void    mad_tpsa_unit    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_abs     (const tpsa_t *a, tpsa_t *c);
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
void    mad_tpsa_asinc   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_asinh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acosh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_atanh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acoth   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_asinhc  (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_erf     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_erfc    (const tpsa_t *a, tpsa_t *c);

void    mad_tpsa_acc     (const tpsa_t *a, num_t v, tpsa_t *c); // c += v*a, aliasing OK
void    mad_tpsa_scl     (const tpsa_t *a, num_t v, tpsa_t *c); // c  = v*a
void    mad_tpsa_inv     (const tpsa_t *a, num_t v, tpsa_t *c); // c  = v/a
void    mad_tpsa_invsqrt (const tpsa_t *a, num_t v, tpsa_t *c); // c  = v/sqrt(a)

void    mad_tpsa_atan2   (const tpsa_t *y, const tpsa_t *x, tpsa_t *r);
void    mad_tpsa_hypot   (const tpsa_t *x, const tpsa_t *y, tpsa_t *r);
void    mad_tpsa_hypot3  (const tpsa_t *x, const tpsa_t *y, const tpsa_t *z, tpsa_t *r);

// functions for differential algebra
void    mad_tpsa_integ   (const tpsa_t *a, tpsa_t *c, idx_t iv);
void    mad_tpsa_deriv   (const tpsa_t *a, tpsa_t *c, idx_t iv);
void    mad_tpsa_derivm  (const tpsa_t *a, tpsa_t *c, ssz_t n, const ord_t m[]);
void    mad_tpsa_poisbra (const tpsa_t *a, const tpsa_t *b, tpsa_t *c, int nv);
void    mad_tpsa_taylor  (const tpsa_t *a, ssz_t n, const num_t coef[], tpsa_t *c);

// high level functions
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

// to check for non-homogeneous maps & parameters
void    mad_tpsa_vec2fld  (ssz_t na, const tpsa_t *a   ,                     tpsa_t *mc[]);
void    mad_tpsa_fld2vec  (ssz_t na, const tpsa_t *ma[],                     tpsa_t *c   );
void    mad_tpsa_fgrad    (ssz_t na, const tpsa_t *ma[], const tpsa_t * b  , tpsa_t *c   );
void    mad_tpsa_liebra   (ssz_t na, const tpsa_t *ma[], const tpsa_t *mb[], tpsa_t *mc[]);
void    mad_tpsa_exppb    (ssz_t na, const tpsa_t *ma[], const tpsa_t *mb[], tpsa_t *mc[]); // exp(:F:) K
void    mad_tpsa_logpb    (ssz_t na, const tpsa_t *ma[], const tpsa_t *mb[], tpsa_t *mc[]); // exp(log(:F:))K

ord_t   mad_tpsa_mord     (ssz_t na, const tpsa_t *ma[], log_t hi); // max mo (or max hi)
num_t   mad_tpsa_mnrm     (ssz_t na, const tpsa_t *ma[]);
void    mad_tpsa_minv     (ssz_t na, const tpsa_t *ma[], ssz_t nb,                     tpsa_t *mc[]);
void    mad_tpsa_pminv    (ssz_t na, const tpsa_t *ma[], ssz_t nb,                     tpsa_t *mc[], idx_t select[]);
void    mad_tpsa_compose  (ssz_t na, const tpsa_t *ma[], ssz_t nb, const tpsa_t *mb[], tpsa_t *mc[]);
void    mad_tpsa_translate(ssz_t na, const tpsa_t *ma[], ssz_t nb, const num_t   tb[], tpsa_t *mc[]);
void    mad_tpsa_eval     (ssz_t na, const tpsa_t *ma[], ssz_t nb, const num_t   tb[], num_t   tc[]);
void    mad_tpsa_mconv    (ssz_t na, const tpsa_t *ma[], ssz_t nc,                     tpsa_t *mc[], ssz_t n, idx_t t2r_[], int pb);

// I/O
void    mad_tpsa_print    (const tpsa_t *t, str_t name_, num_t eps_, int nohdr_, FILE *stream_);
tpsa_t* mad_tpsa_scan     (                                                      FILE *stream_);
const
desc_t* mad_tpsa_scan_hdr (     int *kind_, char  name_[NAMSZ],                  FILE *stream_);
void    mad_tpsa_scan_coef(      tpsa_t *t,                                      FILE *stream_);

// unsafe operation (mo vs allocated!!)
tpsa_t* mad_tpsa_init     (      tpsa_t *t, const desc_t *d, ord_t mo);

// debug
int     mad_tpsa_debug    (const tpsa_t *t, str_t name_, str_t fnam_, int line_, FILE *stream_);
void    mad_tpsa_prtdensity(FILE *stream_);
void    mad_tpsa_clrdensity(void);

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_H
