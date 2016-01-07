#ifndef MAD_CTPSA_H
#define MAD_CTPSA_H

/*
 o----------------------------------------------------------------------------o
 |
 | Complex Truncated Power Series Algebra module interface
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
  
  Purpose:
  - provide a full feathered Generalized Complex TPSA package
 
  Information:
  - parameters ending with an underscope can be null.

  Errors:
  - TODO

 o----------------------------------------------------------------------------o
 */

#include <stdio.h>

#include "mad.h"
#include "mad_mono.h"

// --- types -----------------------------------------------------------------o

struct ctpsa;
struct tpsa_desc;

// --- globals ---------------------------------------------------------------o

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;
extern       int   mad_tpsa_strict;

// --- interface -------------------------------------------------------------o

#define T struct ctpsa

// ctors, dtor
T*     mad_ctpsa_newd    (struct tpsa_desc
                                  *d, ord_t mo); // if mo > d_mo, mo = d_mo
T*     mad_ctpsa_new     (const T *t, ord_t mo);
void   mad_ctpsa_del     (      T *t);
void   mad_ctpsa_delv    (      T *t1, T *t2, ...);

// introspection
struct tpsa_desc*
       mad_ctpsa_desc    (const T *t);
ord_t  mad_ctpsa_ord     (const T *t);
ord_t  mad_ctpsa_ordv    (const T *t1, const T *t2, ...);  // max order of all

// initialization
T*     mad_ctpsa_copy    (const T *t, T *dst);
void   mad_ctpsa_clear   (      T *t);
void   mad_ctpsa_scalar  (      T *t, cnum_t v);

// indexing / monomials
const ord_t*
       mad_ctpsa_mono    (const T *t, int i, int *n, ord_t *total_ord_);
int    mad_ctpsa_midx    (const T *t, int n, const ord_t m[n]);
int    mad_ctpsa_midx_sp (const T *t, int n, const int   m[n]); // sparse mono [(i,o)]

// accessors
cnum_t mad_ctpsa_get0    (const T *t);
cnum_t mad_ctpsa_geti    (const T *t, int i);
cnum_t mad_ctpsa_getm    (const T *t, int n, const ord_t m[n]);
cnum_t mad_ctpsa_getm_sp (const T *t, int n, const int   m[n]); // sparse mono [(i,o)]
void   mad_ctpsa_set0    (      T *t, /* i = 0 */              cnum_t a, cnum_t b);
void   mad_ctpsa_seti    (      T *t, int i,                   cnum_t a, cnum_t b);
void   mad_ctpsa_setm    (      T *t, int n, const ord_t m[n], cnum_t a, cnum_t b);
void   mad_ctpsa_setm_sp (      T *t, int n, const int   m[n], cnum_t a, cnum_t b);

// tranformations
T*     mad_ctpsa_map     (const T *a,             T *c, cnum_t (*f)(cnum_t v, int i_));
T*     mad_ctpsa_map2    (const T *a, const T *b, T *c, cnum_t (*f)(cnum_t va, cnum_t vb, int i_));

// operations
void   mad_ctpsa_abs     (const T *a, T *c);
cnum_t mad_ctpsa_nrm1    (const T *t, const T *t2_);
cnum_t mad_ctpsa_nrm2    (const T *t, const T *t2_);
void   mad_ctpsa_der     (const T *a, T *c, int var);  // TODO: check functions that rely on it
void   mad_ctpsa_mder    (const T *a, T *c, int n, const ord_t m[n]);

void   mad_ctpsa_add     (const T *a, const T *b, T *c);
void   mad_ctpsa_sub     (const T *a, const T *b, T *c);
void   mad_ctpsa_mul     (const T *a, const T *b, T *c);
void   mad_ctpsa_div     (const T *a, const T *b, T *c);

void   mad_ctpsa_acc     (const T *a, cnum_t v, T *c);  // c += v*a, aliasing OK
void   mad_ctpsa_scl     (const T *a, cnum_t v, T *c);  // c  = v*a
void   mad_ctpsa_inv     (const T *a, cnum_t v, T *c);  // c  = v/a
void   mad_ctpsa_invsqrt (const T *a, cnum_t v, T *c);  // c  = v/sqrt(a)
void   mad_ctpsa_sqrt    (const T *a, T *c);
void   mad_ctpsa_exp     (const T *a, T *c);
void   mad_ctpsa_log     (const T *a, T *c);
void   mad_ctpsa_sin     (const T *a, T *c);
void   mad_ctpsa_cos     (const T *a, T *c);
void   mad_ctpsa_sinh    (const T *a, T *c);
void   mad_ctpsa_cosh    (const T *a, T *c);
void   mad_ctpsa_sincos  (const T *a, T *s, T *c);
void   mad_ctpsa_sincosh (const T *a, T *s, T *c);
void   mad_ctpsa_sinc    (const T *a, T *c);
void   mad_ctpsa_sirx    (const T *a, T *c);
void   mad_ctpsa_corx    (const T *a, T *c);

void   mad_ctpsa_tan     (const T *a, T *c);
void   mad_ctpsa_cot     (const T *a, T *c);
void   mad_ctpsa_asin    (const T *a, T *c);
void   mad_ctpsa_acos    (const T *a, T *c);
void   mad_ctpsa_atan    (const T *a, T *c);
void   mad_ctpsa_acot    (const T *a, T *c);
void   mad_ctpsa_tanh    (const T *a, T *c);
void   mad_ctpsa_coth    (const T *a, T *c);
void   mad_ctpsa_asinh   (const T *a, T *c);
void   mad_ctpsa_acosh   (const T *a, T *c);
void   mad_ctpsa_atanh   (const T *a, T *c);
void   mad_ctpsa_acoth   (const T *a, T *c);

void   mad_ctpsa_erf     (const T *a, T *c);

void   mad_ctpsa_ipow    (const T *a, T *c, int n);

void   mad_ctpsa_axpb       (cnum_t a, const T *x, cnum_t b,                       T *r);  // aliasing OK
void   mad_ctpsa_axpbypc    (cnum_t a, const T *x, cnum_t b, const T *y, cnum_t c, T *r);  // aliasing OK
void   mad_ctpsa_axypb      (cnum_t a, const T *x,           const T *y, cnum_t b, T *r);  // aliasing OK
void   mad_ctpsa_axypbzpc   (cnum_t a, const T *x,           const T *y, cnum_t b,
                                                             const T *z, cnum_t c, T *r);  // aliasing OK
void   mad_ctpsa_axypbvwpc  (cnum_t a, const T *x,           const T *y,
                             cnum_t b, const T *v,           const T *w, cnum_t c, T *r);  // aliasing OK
void   mad_ctpsa_ax2pby2pcz2(cnum_t a, const T *x, cnum_t b, const T *y, cnum_t c, const T *z, T *r); // aliasing OK

// to check for non-homogeneous maps & knobs
void   mad_ctpsa_poisson (const T *a, const T *b, T *c, int n);  // TO CHECK n
void   mad_ctpsa_compose (int sa, const T *ma[], int sb, const T *mb[], int sc, T *mc[]);
void   mad_ctpsa_minv    (int sa, const T *ma[],                        int sc, T *mc[]);
void   mad_ctpsa_pminv   (int sa, const T *ma[],                        int sc, T *mc[], int row_select[sa]);

// I/O
void   mad_ctpsa_scan_coef(      T *t, FILE *stream_); // TODO
T*     mad_ctpsa_scan     (            FILE *stream_); // TODO
void   mad_ctpsa_print    (const T *t, str_t name_, FILE *stream_);
void   mad_ctpsa_debug    (const T *t);

#define mad_ctpsa_ordv(...) mad_ctpsa_ordv(__VA_ARGS__,NULL)
#define mad_ctpsa_delv(...) mad_ctpsa_delv(__VA_ARGS__,NULL)

#undef T

// ---------------------------------------------------------------------------o

#endif // MAD_CTPSA_H
