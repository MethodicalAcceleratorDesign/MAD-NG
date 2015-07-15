local oss  = jit.os
local ffi  = require 'ffi'
local cmad = ffi.load("lib/cmad/libmad-" .. oss .. ".so")

ffi.cdef[[
typedef double num_t;
typedef double _Complex cnum_t;
]]

-- functions for complex numbers for LuaJIT

ffi.cdef [[
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

void  mad_cnum_div   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);
void  mad_cnum_pow   (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cnum_t *r);
]]

-- functions for vector-vector, vector-scalar operatrions for LuaJIT
-- note: matrices can be treated as vector for element wise operations

ffi.cdef[[
num_t mad_vec_dot   (const  num_t *x, const  num_t *y       ,            size_t n); // <vec, vec>
void  mad_vec_dotv  (const  num_t *x, const cnum_t *y       , cnum_t *r, size_t n); // <vec,cvec>
void  mad_vec_add   (const  num_t *x, const  num_t *y       ,  num_t *r, size_t n); // vec +  vec
void  mad_vec_addn  (const  num_t *x,        num_t  y       ,  num_t *r, size_t n); // vec +  num
void  mad_vec_addc  (const  num_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n); // vec +  cpx
void  mad_vec_sub   (const  num_t *x, const  num_t *y       ,  num_t *r, size_t n); // vec -  vec
void  mad_vec_subv  (const  num_t *x, const cnum_t *y       , cnum_t *r, size_t n); // vec - cvec
void  mad_vec_subn  (const  num_t *y,        num_t  x       ,  num_t *r, size_t n); // num -  vec
void  mad_vec_subc  (const  num_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n); // cpx -  vec
void  mad_vec_mul   (const  num_t *x, const  num_t *y       ,  num_t *r, size_t n); // vec *  vec
void  mad_vec_muln  (const  num_t *x,        num_t  y       ,  num_t *r, size_t n); // vec *  num
void  mad_vec_mulc  (const  num_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n); // vec *  cpx
void  mad_vec_div   (const  num_t *x, const  num_t *y       ,  num_t *r, size_t n); // vec /  vec
void  mad_vec_divv  (const  num_t *x, const cnum_t *y       , cnum_t *r, size_t n); // vec / cvec
void  mad_vec_divn  (const  num_t *y,        num_t  x       ,  num_t *r, size_t n); // num /  vec
void  mad_vec_divc  (const  num_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n); // cpx /  vec 

void  mad_cvec_dot  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); // <cvec,cvec>
void  mad_cvec_dotv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); // <cvec, vec>
void  mad_cvec_add  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); // cvec + cvec
void  mad_cvec_addv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); // cvec +  vec
void  mad_cvec_addn (const cnum_t *x,        num_t  y       , cnum_t *r, size_t n); // cvec +  num
void  mad_cvec_addc (const cnum_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n); // cvec +  cpx
void  mad_cvec_sub  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); // cvec - cvec
void  mad_cvec_subv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); // cvec -  vec
void  mad_cvec_subn (const cnum_t *y,        num_t  x       , cnum_t *r, size_t n); // num  - cvec
void  mad_cvec_subc (const cnum_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n); // cpx  - cvec
void  mad_cvec_mul  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); // cvec * cvec
void  mad_cvec_mulv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); // cvec *  vec
void  mad_cvec_muln (const cnum_t *x,        num_t  y       , cnum_t *r, size_t n); // cvec *  num
void  mad_cvec_mulc (const cnum_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n); // cvec *  cpx
void  mad_cvec_div  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); // cvec / cvec
void  mad_cvec_divv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); // cvec /  vec
void  mad_cvec_divn (const cnum_t *y,        num_t  x       , cnum_t *r, size_t n); // num  / cvec
void  mad_cvec_divc (const cnum_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n); // cpx  / cvec
]]

-- functions for matrix-matrix, vector-matrix and matrix-vector multiplications for LuaJIT

ffi.cdef[[
void mad_mat_trans (const  num_t *x,                   num_t *r, size_t m, size_t n);           //  mat.t()
void mad_mat_mul   (const  num_t *x, const  num_t *y,  num_t *r, size_t m, size_t n, size_t p); //  mat *  mat
void mad_mat_mulm  (const  num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p); //  mat * cmat
void mad_mat_muln  (const  num_t *x, const  num_t *y,  num_t *r, size_t m, size_t p);           //  mat *  vec
void mad_mat_mulc  (const  num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t p);           //  mat * cvec
void mad_mat_nmul  (const  num_t *x, const  num_t *y,  num_t *r, size_t n, size_t p);           //  vec *  mat
void mad_mat_cmul  (const cnum_t *x, const  num_t *y, cnum_t *r, size_t n, size_t p);           // cvec *  mat

void mad_cmat_trans(const cnum_t *x,                  cnum_t *r, size_t m, size_t n);           // cmat.t()
void mad_cmat_mul  (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p); // cmat * cmat
void mad_cmat_mulm (const cnum_t *x, const  num_t *y, cnum_t *r, size_t m, size_t n, size_t p); // cmat *  mat
void mad_cmat_muln (const cnum_t *x, const  num_t *y, cnum_t *r, size_t m, size_t p);           // cmat *  vec
void mad_cmat_mulc (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t p);           // cmat * cvec
void mad_cmat_nmul (const  num_t *x, const cnum_t *y, cnum_t *r, size_t n, size_t p);           //  vec * cmat
void mad_cmat_cmul (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t n, size_t p);           // cvec * cmat
]]

return cmad