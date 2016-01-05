#ifndef MAD_TPSA_H
#define MAD_TPSA_H

/*
 o----------------------------------------------------------------------------o
 |
 | Truncated Power Series Algebra module interface
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
  - provide a full feathered Generalized TPSA package
 
  Information:
  - parameters ending with an underscope can be null.

  Errors:
  - TODO

 o----------------------------------------------------------------------------o
 */

#include <stdio.h>
#include "mad_mono.h"

// --- types -----------------------------------------------------------------o

struct tpsa;
struct tpsa_desc;

// --- interface -------------------------------------------------------------o

#define T struct tpsa
#define D struct tpsa_desc

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;
extern       int   mad_tpsa_strict;

// descriptors (tpsa factories, bounded to maps)
D*    mad_tpsa_desc_new (int nv, const ord_t var_ords[nv], const ord_t map_ords_[nv], str_t var_nam_[nv]);
D*    mad_tpsa_desc_newk(int nv, const ord_t var_ords[nv], const ord_t map_ords_[nv], str_t var_nam_[nv],
                         int nk, const ord_t knb_ords[nk], ord_t dk); // knobs
void  mad_tpsa_desc_del (D *d);

// Introspection
ord_t mad_tpsa_gtrunc  (      D *d, ord_t to);
int   mad_tpsa_maxsize (const D *d);
ord_t mad_tpsa_maxord  (const D *d);
D*    mad_tpsa_desc    (const T *t);
ord_t mad_tpsa_ord     (const T *t);
ord_t mad_tpsa_ordv    (const T *t1, const T *t2, ...);  // max order of all

// ctors, dtor
T*    mad_tpsa_newd    (      D *d, ord_t mo); // if mo > d_mo, mo = d_mo
T*    mad_tpsa_new     (const T *t, ord_t mo);
T*    mad_tpsa_copy    (const T *t, T *dst);
void  mad_tpsa_clear   (      T *t);
void  mad_tpsa_scalar  (      T *t, num_t v);
void  mad_tpsa_del     (      T *t);
void  mad_tpsa_delv    (      T *t1, T *t2, ...);

// indexing / monomials
const ord_t*
      mad_tpsa_mono    (const T *t, int i, int *n, ord_t *total_ord_);
int   mad_tpsa_midx    (const T *t, int n, const ord_t m[n]);
int   mad_tpsa_midx_sp (const T *t, int n, const int   m[n]); // sparse mono [(i,o)]

// accessors
num_t mad_tpsa_get0    (const T *t);
num_t mad_tpsa_geti    (const T *t, int i);
num_t mad_tpsa_getm    (const T *t, int n, const ord_t m[n]);
num_t mad_tpsa_getm_sp (const T *t, int n, const int   m[n]); // sparse mono [(i,o)]
void  mad_tpsa_set0    (      T *t, /* i = 0 */              num_t a, num_t b);
void  mad_tpsa_seti    (      T *t, int i,                   num_t a, num_t b);
void  mad_tpsa_setm    (      T *t, int n, const ord_t m[n], num_t a, num_t b);
void  mad_tpsa_setm_sp (      T *t, int n, const int   m[n], num_t a, num_t b);

// tranformations
T*    mad_tpsa_map     (const T *a,             T *c, num_t (*f)(num_t v, int i_));
T*    mad_tpsa_map2    (const T *a, const T *b, T *c, num_t (*f)(num_t va, num_t vb, int i_));

// operations
void  mad_tpsa_abs     (const T *a, T *c);
num_t mad_tpsa_nrm1    (const T *t, const T *t2_);
num_t mad_tpsa_nrm2    (const T *t, const T *t2_);
void  mad_tpsa_der     (const T *a, T *c, int var);  // TODO: check functions that rely on it
void  mad_tpsa_mder    (const T *a, T *c, int n, const ord_t m[n]);

void  mad_tpsa_add     (const T *a, const T *b, T *c);
void  mad_tpsa_sub     (const T *a, const T *b, T *c);
void  mad_tpsa_mul     (const T *a, const T *b, T *c);
void  mad_tpsa_div     (const T *a, const T *b, T *c);

void  mad_tpsa_scl     (const T *a, num_t v, T *c);  // c = v*a
void  mad_tpsa_inv     (const T *a, num_t v, T *c);  // c = v/a
void  mad_tpsa_invsqrt (const T *a, num_t v, T *c);  // c = v/sqrt(a)
void  mad_tpsa_sqrt    (const T *a, T *c);
void  mad_tpsa_exp     (const T *a, T *c);
void  mad_tpsa_log     (const T *a, T *c);
void  mad_tpsa_sin     (const T *a, T *c);
void  mad_tpsa_cos     (const T *a, T *c);
void  mad_tpsa_sinh    (const T *a, T *c);
void  mad_tpsa_cosh    (const T *a, T *c);
void  mad_tpsa_sincos  (const T *a, T *s, T *c);
void  mad_tpsa_sincosh (const T *a, T *s, T *c);
void  mad_tpsa_sinc    (const T *a, T *c);
void  mad_tpsa_sirx    (const T *a, T *c);
void  mad_tpsa_corx    (const T *a, T *c);

void  mad_tpsa_tan     (const T *a, T *c);
void  mad_tpsa_cot     (const T *a, T *c);
void  mad_tpsa_asin    (const T *a, T *c);
void  mad_tpsa_acos    (const T *a, T *c);
void  mad_tpsa_atan    (const T *a, T *c);
void  mad_tpsa_acot    (const T *a, T *c);
void  mad_tpsa_tanh    (const T *a, T *c);
void  mad_tpsa_coth    (const T *a, T *c);
void  mad_tpsa_asinh   (const T *a, T *c);
void  mad_tpsa_acosh   (const T *a, T *c);
void  mad_tpsa_atanh   (const T *a, T *c);
void  mad_tpsa_acoth   (const T *a, T *c);

void  mad_tpsa_erf     (const T *a, T *c);

void  mad_tpsa_ipow    (const T *a, T *c, int n);

void  mad_tpsa_axpb       (num_t a, const T *x, num_t b,                      T *r);  // aliasing OK
void  mad_tpsa_axpbypc    (num_t a, const T *x, num_t b, const T *y, num_t c, T *r);  // aliasing OK
void  mad_tpsa_axypb      (num_t a, const T *x,          const T *y, num_t b, T *r);  // aliasing OK
void  mad_tpsa_axypbzpc   (num_t a, const T *x,          const T *y, num_t b,
                                                         const T *z, num_t c, T *r);  // aliasing OK
void  mad_tpsa_axypbvwpc  (num_t a, const T *x,          const T *y,
                           num_t b, const T *v,          const T *w, num_t c, T *r);  // aliasing OK
void  mad_tpsa_ax2pby2pcz2(num_t a, const T *x, num_t b, const T *y, num_t c, const T *z, T *r); // aliasing OK

// to check for non-homogeneous maps & knobs
void  mad_tpsa_poisson (const T *a, const T *b, T *c, int n);  // TO CHECK n
void  mad_tpsa_compose (int sa, const T *ma[], int sb, const T *mb[], int sc, T *mc[]);
void  mad_tpsa_minv    (int sa, const T *ma[],                        int sc, T *mc[]);
void  mad_tpsa_pminv   (int sa, const T *ma[],                        int sc, T *mc[], int row_select[sa]);

// I/O
void  mad_tpsa_scan_coef(      T *t, FILE *stream_);
T*    mad_tpsa_scan     (            FILE *stream_); // TODO
void  mad_tpsa_print    (const T *t, FILE *stream_);
D*    mad_tpsa_scan_desc(            FILE *stream_);
void  mad_tpsa_debug    (const T *t);

#define mad_tpsa_ordv(...) mad_tpsa_ordv(__VA_ARGS__,NULL)
#define mad_tpsa_delv(...) mad_tpsa_delv(__VA_ARGS__,NULL)

#undef T
#undef D

// ---------------------------------------------------------------------------o

#endif // MAD_TPSA_H
