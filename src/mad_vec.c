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

#include "mad_vec.h"

// --- implementation --------------------------------------------------------o

#define CHKR    assert( r )
#define CHKX    assert( x )
#define CHKXY   assert( x && y )
#define CHKXR   assert( x && r )
#define CHKYR   assert( y && r )
#define CHKXYR  assert( x && y && r )

#define CNUM2(re,im) (* (cnum_t*) & (num_t[2]) { re, im })
#define CNUMR(re)    CNUM2(re,0)
#define CNUMI(im)    CNUM2(0,im)
#define CNUM(a)      cnum_t a = CNUM2(MKNAME(a,_re), MKNAME(a,_im))

#define SET(OP)      for (size_t i=0; i < n; i++)  r[i] OP##= x
#define CPY(OP)      for (size_t i=0; i < n; i++)  r[i] OP##= x[i]
#define MAP(FN)      for (size_t i=0; i < n; i++)  r[i]  = FN(x[i])
#define MAP2(FN)     for (size_t i=0; i < n; i++)  r[i]  = FN(x[i], y[i])
#define VEC(OP)      for (size_t i=0; i < n; i++)  r[i]     = x[i] OP y[i]
#define VECS(OP)     for (size_t i=0; i < n; i++)  r[i]     = x[i] OP y
#define SVEC(OP)     for (size_t i=0; i < n; i++)  r[i]     = x    OP y[i]
#define DOT(C) *r=0; for (size_t i=0; i < n; i++) *r       += C(x[i]) * y[i]

// --- vec

void mad_vec_fill (num_t x, num_t r[], size_t n)
{ CHKR; SET(); }

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
{ CHKXY; num_t r_, *r=&r_; DOT(); return *r; }

cnum_t mad_vec_dotv (const  num_t x[], const cnum_t y[], size_t n)
{ CHKXY; cnum_t r_, *r=&r_; DOT(); return *r;}

void mad_vec_dotv_r (const  num_t x[], const cnum_t y[], cnum_t *r, size_t n)
{ CHKXYR; DOT(); }

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

// --- cvec
 
void mad_cvec_fill (cnum_t x, cnum_t r[], size_t n)
{ CHKR; SET(); }

void mad_cvec_fill_r (num_t x_re, num_t x_im, cnum_t r[], size_t n)
{ CHKR; CNUM(x); SET(); }

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
{ CHKXY; cnum_t r_, *r=&r_; DOT(conj); return *r; }

void mad_cvec_dot_r (const cnum_t x[], const cnum_t y[], cnum_t *r, size_t n)
{ CHKXYR; DOT(conj); }

cnum_t mad_cvec_dotv (const cnum_t x[], const num_t y[], size_t n)
{ CHKXY; cnum_t r_, *r=&r_; DOT(conj); return *r; }

void mad_cvec_dotv_r (const cnum_t x[], const num_t y[], cnum_t *r, size_t n)
{ CHKXYR; DOT(conj); }

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
