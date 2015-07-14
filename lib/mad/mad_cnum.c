// gcc -std=c11 -W -Wall -Wextra -pedantic -O3 -ffast-math -ftree-vectorize -shared -fPIC -static-libgcc *.c -o libmad-OSX.so -lm

#include <complex.h>
#include "mad_cnum.h"

typedef double           num_t;
typedef double _Complex cnum_t;

#define CNUM(re,im) (* (cnum_t*) & (num_t[2]) { (re), (im) })

num_t mad_cnum_abs  (num_t x_re, num_t x_im)  { return cabs( CNUM(x_re, x_im) ); }
num_t mad_cnum_arg  (num_t x_re, num_t x_im)  { return carg( CNUM(x_re, x_im) ); }

void mad_cnum_exp   (num_t x_re, num_t x_im, cnum_t *r)  { *r = cexp   ( CNUM(x_re, x_im) ); }
void mad_cnum_log   (num_t x_re, num_t x_im, cnum_t *r)  { *r = clog   ( CNUM(x_re, x_im) ); }
void mad_cnum_sqrt  (num_t x_re, num_t x_im, cnum_t *r)  { *r = csqrt  ( CNUM(x_re, x_im) ); }
void mad_cnum_proj  (num_t x_re, num_t x_im, cnum_t *r)  { *r = cproj  ( CNUM(x_re, x_im) ); }
void mad_cnum_sin   (num_t x_re, num_t x_im, cnum_t *r)  { *r = csin   ( CNUM(x_re, x_im) ); }
void mad_cnum_cos   (num_t x_re, num_t x_im, cnum_t *r)  { *r = ccos   ( CNUM(x_re, x_im) ); }
void mad_cnum_tan   (num_t x_re, num_t x_im, cnum_t *r)  { *r = ctan   ( CNUM(x_re, x_im) ); }
void mad_cnum_sinh  (num_t x_re, num_t x_im, cnum_t *r)  { *r = csinh  ( CNUM(x_re, x_im) ); }
void mad_cnum_cosh  (num_t x_re, num_t x_im, cnum_t *r)  { *r = ccosh  ( CNUM(x_re, x_im) ); }
void mad_cnum_tanh  (num_t x_re, num_t x_im, cnum_t *r)  { *r = ctanh  ( CNUM(x_re, x_im) ); }
void mad_cnum_asin  (num_t x_re, num_t x_im, cnum_t *r)  { *r = casin  ( CNUM(x_re, x_im) ); }
void mad_cnum_acos  (num_t x_re, num_t x_im, cnum_t *r)  { *r = cacos  ( CNUM(x_re, x_im) ); }
void mad_cnum_atan  (num_t x_re, num_t x_im, cnum_t *r)  { *r = catan  ( CNUM(x_re, x_im) ); }
void mad_cnum_asinh (num_t x_re, num_t x_im, cnum_t *r)  { *r = casinh ( CNUM(x_re, x_im) ); }
void mad_cnum_acosh (num_t x_re, num_t x_im, cnum_t *r)  { *r = cacosh ( CNUM(x_re, x_im) ); }
void mad_cnum_atanh (num_t x_re, num_t x_im, cnum_t *r)  { *r = catanh ( CNUM(x_re, x_im) ); }

void mad_cnum_pow   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r) {
  *r = cpow( CNUM(x_re, x_im), CNUM(y_re, y_im) );
}
