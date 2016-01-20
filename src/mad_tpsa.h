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
#include "mad_desc.h"

// --- types -----------------------------------------------------------------o

typedef struct tpsa tpsa_t;

// --- globals ---------------------------------------------------------------o

extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;
extern       int   mad_tpsa_strict;

// --- interface -------------------------------------------------------------o

// ctors, dtor
tpsa_t* mad_tpsa_newd    (desc_t *d, ord_t mo); // if mo > d_mo, mo = d_mo
tpsa_t* mad_tpsa_new     (const tpsa_t *t, ord_t mo);
void    mad_tpsa_del     (      tpsa_t *t);
void    mad_tpsa_delv    (      tpsa_t *t1, tpsa_t *t2, ...);

// introspection
desc_t* mad_tpsa_desc    (const tpsa_t *t);
ord_t   mad_tpsa_ord     (const tpsa_t *t);
ord_t   mad_tpsa_ordv    (const tpsa_t *t1, const tpsa_t *t2, ...);  // max order of all

// initialization
tpsa_t* mad_tpsa_copy    (const tpsa_t *t, tpsa_t *dst);
void    mad_tpsa_clear   (      tpsa_t *t);
void    mad_tpsa_scalar  (      tpsa_t *t, num_t v);

// indexing / monomials
int     mad_tpsa_mono    (const tpsa_t *t, int n,       ord_t m_[], idx_t i);
int     mad_tpsa_midx    (const tpsa_t *t, int n, const ord_t m []);
int     mad_tpsa_midx_sp (const tpsa_t *t, int n, const int   m []); // sparse mono [(i,o)]

// accessors
num_t   mad_tpsa_get0    (const tpsa_t *t);
num_t   mad_tpsa_geti    (const tpsa_t *t, idx_t i);
num_t   mad_tpsa_getm    (const tpsa_t *t, int n, const ord_t m[]);
num_t   mad_tpsa_getm_sp (const tpsa_t *t, int n, const int   m[]); // sparse mono [(i,o)]
void    mad_tpsa_set0    (      tpsa_t *t, /* i = 0 */             num_t a, num_t b);
void    mad_tpsa_seti    (      tpsa_t *t, idx_t i,                num_t a, num_t b);
void    mad_tpsa_setm    (      tpsa_t *t, int n, const ord_t m[], num_t a, num_t b);
void    mad_tpsa_setm_sp (      tpsa_t *t, int n, const int   m[], num_t a, num_t b);

// operations
void    mad_tpsa_abs     (const tpsa_t *a, tpsa_t *c);
num_t   mad_tpsa_nrm1    (const tpsa_t *t, const tpsa_t *t2_);
num_t   mad_tpsa_nrm2    (const tpsa_t *t, const tpsa_t *t2_);
void    mad_tpsa_der     (const tpsa_t *a, tpsa_t *c, int var);  // TODO: check functions that rely on it
void    mad_tpsa_mder    (const tpsa_t *a, tpsa_t *c, int n, const ord_t m[]);

void    mad_tpsa_add     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_sub     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_mul     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);
void    mad_tpsa_div     (const tpsa_t *a, const tpsa_t *b, tpsa_t *c);

void    mad_tpsa_acc     (const tpsa_t *a, num_t v, tpsa_t *c);  // c += v*a, aliasing OK
void    mad_tpsa_scl     (const tpsa_t *a, num_t v, tpsa_t *c);  // c  = v*a
void    mad_tpsa_inv     (const tpsa_t *a, num_t v, tpsa_t *c);  // c  = v/a
void    mad_tpsa_invsqrt (const tpsa_t *a, num_t v, tpsa_t *c);  // c  = v/sqrt(a)

void    mad_tpsa_sqrt    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_exp     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_log     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_sin     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_cos     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_sinh    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_cosh    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_sincos  (const tpsa_t *a, tpsa_t *s, tpsa_t *c);
void    mad_tpsa_sincosh (const tpsa_t *a, tpsa_t *s, tpsa_t *c);
void    mad_tpsa_sinc    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_sirx    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_corx    (const tpsa_t *a, tpsa_t *c);

void    mad_tpsa_tan     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_cot     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_asin    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acos    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_atan    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acot    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_tanh    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_coth    (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_asinh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acosh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_atanh   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_acoth   (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_erf     (const tpsa_t *a, tpsa_t *c);
void    mad_tpsa_ipow    (const tpsa_t *a, tpsa_t *c, int n);

// high level functions
void    mad_tpsa_axpb       (num_t a, const tpsa_t *x,
                             num_t b, tpsa_t *r);  // aliasing OK
void    mad_tpsa_axpbypc    (num_t a, const tpsa_t *x,
                             num_t b, const tpsa_t *y,
                             num_t c, tpsa_t *r);  // aliasing OK
void    mad_tpsa_axypb      (num_t a, const tpsa_t *x, const tpsa_t *y,
                             num_t b, tpsa_t *r);  // aliasing OK
void    mad_tpsa_axypbzpc   (num_t a, const tpsa_t *x, const tpsa_t *y,
                             num_t b, const tpsa_t *z,
                             num_t c, tpsa_t *r);  // aliasing OK
void    mad_tpsa_axypbvwpc  (num_t a, const tpsa_t *x, const tpsa_t *y,
                             num_t b, const tpsa_t *v, const tpsa_t *w,
                             num_t c, tpsa_t *r);  // aliasing OK
void    mad_tpsa_ax2pby2pcz2(num_t a, const tpsa_t *x,
                             num_t b, const tpsa_t *y,
                             num_t c, const tpsa_t *z, tpsa_t *r); // aliasing OK

// to check for non-homogeneous maps & knobs
void    mad_tpsa_poisson (const tpsa_t *a, const tpsa_t *b, tpsa_t *c, int n);  // TO CHECK n
void    mad_tpsa_compose (int sa, const tpsa_t *ma[], int sb, const tpsa_t *mb[], int sc, tpsa_t *mc[]);
void    mad_tpsa_minv    (int sa, const tpsa_t *ma[],                             int sc, tpsa_t *mc[]);
void    mad_tpsa_pminv   (int sa, const tpsa_t *ma[],                             int sc, tpsa_t *mc[], int row_select[]);

// I/O
void    mad_tpsa_print    (const tpsa_t *t, str_t name_, FILE *stream_);
tpsa_t* mad_tpsa_scan     (                              FILE *stream_); // TODO
desc_t* mad_tpsa_scan_hdr (                              FILE *stream_);
void    mad_tpsa_scan_coef(      tpsa_t *t,              FILE *stream_); // TODO
void    mad_tpsa_debug    (const tpsa_t *t);

#define mad_tpsa_ordv(...) mad_tpsa_ordv(__VA_ARGS__,NULL)
#define mad_tpsa_delv(...) mad_tpsa_delv(__VA_ARGS__,NULL)

// ---------------------------------------------------------------------------o

#endif // MAD_TPSA_H
