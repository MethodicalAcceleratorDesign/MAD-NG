local oss = jit.os
local ffi = require 'ffi'
local mad = ffi.load("lib/mad/libmad-" .. oss .. ".so")

ffi.cdef [[

// Wrappers around functions of complex numbers for LuaJIT

typedef double num_t;
typedef double _Complex cnum_t;

num_t mad_cnum_abs   (num_t x_re, num_t x_im);
num_t mad_cnum_arg   (num_t x_re, num_t x_im);

void  mad_cnum_exp   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_log   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_sqrt  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_proj  (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_sin   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_cos   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_tan   (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_sinh  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_cosh  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_tanh  (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_asin  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_acos  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_atan  (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_asinh (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_acosh (num_t x_re, num_t x_im, cnum_t *r);
void  mad_cnum_atanh (num_t x_re, num_t x_im, cnum_t *r);

void  mad_cnum_pow   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);

]]

return mad