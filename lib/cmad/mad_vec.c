#include "mad_vec.h"
#include <math.h>
#include <complex.h>
#include <assert.h>

typedef double           num_t;
typedef double _Complex cnum_t;

#define CHKXY   assert( x && y )
#define CHKXR   assert( x && r )
#define CHKYR   assert( y && r )
#define CHKXYR  assert( x && y && r )

#define CNUM(a)  cnum_t a = (* (cnum_t*) & (num_t[2]) { (a##_re), (a##_im) })

#define VEC(OP)        for (size_t i=0; i < n; i++) r[i] = x[i] OP y[i]
#define VECS(OP)       for (size_t i=0; i < n; i++) r[i] = x[i] OP y
#define SVEC(OP)       for (size_t i=0; i < n; i++) r[i] = x    OP y[i]
#define DOT(C)   *r=0; for (size_t i=0; i < n; i++) *r += C(x[i]) * y[i]

num_t mad_vec_dot (const num_t *x, const num_t *y, size_t n)
{ CHKXY; num_t r_, *r=&r_; DOT(); return *r; }

void mad_vec_dotv (const  num_t *x, const cnum_t *y, cnum_t *r, size_t n)
{ CHKXYR; DOT(); }

void mad_vec_add (const num_t *x, const num_t *y, num_t *r, size_t n)
{ CHKXYR; VEC(+); }

void mad_vec_addn (const num_t *x, num_t y, num_t *r, size_t n)
{ CHKXR; VECS(+); }

void mad_vec_addc (const num_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n)
{ CHKXR; CNUM(y); VECS(+); }

void mad_vec_sub (const num_t *x, const num_t *y, num_t *r, size_t n)
{ CHKXYR; VEC(-); }

void mad_vec_subv (const num_t *x, const cnum_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(-); }

void mad_vec_subn (const num_t *y, num_t x, num_t *r, size_t n)
{ CHKYR; SVEC(-); }

void mad_vec_subc (const num_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n)
{ CHKYR; CNUM(x); SVEC(-); }

void mad_vec_mul (const num_t *x, const num_t *y, num_t *r, size_t n)
{ CHKXYR; VEC(*); }

void mad_vec_muln (const num_t *x, num_t y, num_t *r, size_t n)
{ CHKXR; VECS(*); }

void mad_vec_mulc (const num_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n)
{ CHKXR; CNUM(y); VECS(*); }

void mad_vec_div (const num_t *x, const num_t *y, num_t *r, size_t n)
{ CHKXYR; VEC(/); }

void mad_vec_divv (const num_t *x, const cnum_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(/); }

void mad_vec_divn (const num_t *y, num_t x, num_t *r, size_t n)
{ CHKYR; SVEC(/); }

void mad_vec_divc (const num_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n)
{ CHKYR; CNUM(x); SVEC(/); }


void mad_cvec_dot (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t n)
{ CHKXYR; DOT(conj); }

void mad_cvec_dotv (const cnum_t *x, const num_t *y, cnum_t *r, size_t n)
{ CHKXYR; DOT(conj); }

void mad_cvec_add (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(+); }

void mad_cvec_addv (const cnum_t *x, const num_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(+); }

void mad_cvec_addn (const cnum_t *x, num_t y, cnum_t *r, size_t n)
{ CHKXR; VECS(+); }

void mad_cvec_addc (const cnum_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n)
{ CHKXR; CNUM(y); VECS(+); }

void mad_cvec_sub (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(-); }

void mad_cvec_subv (const cnum_t *x, const num_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(-); }

void mad_cvec_subn (const cnum_t *y, num_t x, cnum_t *r, size_t n)
{ CHKYR; SVEC(-); }

void mad_cvec_subc (const cnum_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n)
{ CHKYR; CNUM(x); SVEC(-); }

void mad_cvec_mul (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(*); }

void mad_cvec_mulv (const cnum_t *x, const num_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(*); }

void mad_cvec_muln (const cnum_t *x, num_t y, cnum_t *r, size_t n)
{ CHKXR; VECS(*); }

void mad_cvec_mulc (const cnum_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n)
{ CHKXR; CNUM(y); VECS(*); }

void mad_cvec_div (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(/); }

void mad_cvec_divv (const cnum_t *x, const num_t *y, cnum_t *r, size_t n)
{ CHKXYR; VEC(/); }

void mad_cvec_divn (const cnum_t *y, num_t x, cnum_t *r, size_t n)
{ CHKYR; SVEC(/); }

void mad_cvec_divc (const cnum_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n)
{ CHKYR; CNUM(x); SVEC(/); }
