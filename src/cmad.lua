local ffi = require 'ffi'

-- types

ffi.cdef[[
typedef const char*      str_t;
typedef int              idx_t;
typedef double           num_t;
typedef double _Complex cnum_t;
]]

-- functions for logging (mad_log.h)

ffi.cdef [[
void mad_fatal (str_t msg);
void mad_error (str_t msg);
void mad_warn  (str_t msg);
void mad_info  (int lvl, str_t msg);
void mad_debug (int lvl, str_t msg);

void mad_log_setloc1 (str_t file, int line);
]]

-- functions for memory management (mad_mem.h)

ffi.cdef [[
void*  mad_malloc  (size_t size_);
void*  mad_realloc (void  *ptr_ , size_t size_);
void   mad_free    (void  *ptr_);

size_t mad_mem_size    (void *ptr_);
size_t mad_mem_cached  (void);
size_t mad_mem_collect (void);

// alternate for memcheck
void*  malloc  (size_t size_);
void*  realloc (void  *ptr_ , size_t size_);
void   free    (void  *ptr_);
]]

-- functions for real and complex numbers for LuaJIT (mad_num.h)

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

-- functions for vector-vector and vector-scalar operations for LuaJIT (mad_vec.h)
-- note: matrices can be treated as vector for element wise operations

ffi.cdef[[
void  mad_vec_set   (                        num_t  x       ,  num_t *r, size_t n); //  num -> vec
void  mad_vec_cpy   (const  num_t *x,                          num_t *r, size_t n); //  vec -> vec
void  mad_vec_cpyv  (const  num_t *x,                         cnum_t *r, size_t n); //  vec ->cvec
num_t mad_vec_dot   (const  num_t *x, const  num_t *y       ,            size_t n); // <vec ,  vec>
void  mad_vec_dotv  (const  num_t *x, const cnum_t *y       , cnum_t *r, size_t n); // <vec , cvec>
void  mad_vec_add   (const  num_t *x, const  num_t *y       ,  num_t *r, size_t n); //  vec +  vec
void  mad_vec_addn  (const  num_t *x,        num_t  y       ,  num_t *r, size_t n); //  vec +  num
void  mad_vec_addc  (const  num_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n); //  vec +  cpx
void  mad_vec_sub   (const  num_t *x, const  num_t *y       ,  num_t *r, size_t n); //  vec -  vec
void  mad_vec_subv  (const  num_t *x, const cnum_t *y       , cnum_t *r, size_t n); //  vec - cvec
void  mad_vec_subn  (const  num_t *y,        num_t  x       ,  num_t *r, size_t n); //  num -  vec
void  mad_vec_subc  (const  num_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n); //  cpx -  vec
void  mad_vec_mul   (const  num_t *x, const  num_t *y       ,  num_t *r, size_t n); //  vec *  vec
void  mad_vec_muln  (const  num_t *x,        num_t  y       ,  num_t *r, size_t n); //  vec *  num
void  mad_vec_mulc  (const  num_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n); //  vec *  cpx
void  mad_vec_div   (const  num_t *x, const  num_t *y       ,  num_t *r, size_t n); //  vec /  vec
void  mad_vec_divv  (const  num_t *x, const cnum_t *y       , cnum_t *r, size_t n); //  vec / cvec
void  mad_vec_divn  (const  num_t *y,        num_t  x       ,  num_t *r, size_t n); //  num /  vec
void  mad_vec_divc  (const  num_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n); //  cpx /  vec 

void  mad_cvec_set  (                 num_t x_re, num_t x_im, cnum_t *r, size_t n); //  cnum ->cvec
void  mad_cvec_cpy  (const cnum_t *x,                         cnum_t *r, size_t n); //  cvec ->cvec
void  mad_cvec_cpyv (const cnum_t *x,                          num_t *r, size_t n); //  cvec -> vec
void  mad_cvec_dot  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); // <cvec , cvec>
void  mad_cvec_dotv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); // <cvec ,  vec>
void  mad_cvec_add  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); //  cvec + cvec
void  mad_cvec_addv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); //  cvec +  vec
void  mad_cvec_addn (const cnum_t *x,        num_t  y       , cnum_t *r, size_t n); //  cvec +  num
void  mad_cvec_addc (const cnum_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n); //  cvec +  cpx
void  mad_cvec_sub  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); //  cvec - cvec
void  mad_cvec_subv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); //  cvec -  vec
void  mad_cvec_subn (const cnum_t *y,        num_t  x       , cnum_t *r, size_t n); //  num  - cvec
void  mad_cvec_subc (const cnum_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n); //  cpx  - cvec
void  mad_cvec_mul  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); //  cvec * cvec
void  mad_cvec_mulv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); //  cvec *  vec
void  mad_cvec_muln (const cnum_t *x,        num_t  y       , cnum_t *r, size_t n); //  cvec *  num
void  mad_cvec_mulc (const cnum_t *x, num_t y_re, num_t y_im, cnum_t *r, size_t n); //  cvec *  cpx
void  mad_cvec_div  (const cnum_t *x, const cnum_t *y       , cnum_t *r, size_t n); //  cvec / cvec
void  mad_cvec_divv (const cnum_t *x, const  num_t *y       , cnum_t *r, size_t n); //  cvec /  vec
void  mad_cvec_divn (const cnum_t *y,        num_t  x       , cnum_t *r, size_t n); //  num  / cvec
void  mad_cvec_divc (const cnum_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t n); //  cpx  / cvec
]]

-- functions for matrix-matrix, vector-matrix and matrix-vector operations for LuaJIT (mad_mat.h)

ffi.cdef[[
void  mad_mat_ident   (                                          num_t *r, size_t m, size_t n,             size_t ldr); // ident-> mat
void  mad_mat_set     (                        num_t  x,         num_t *r, size_t m, size_t n,             size_t ldr); //  num -> mat
void  mad_mat_cpy     (const  num_t *x,                          num_t *r, size_t m, size_t n, size_t ldx, size_t ldr); //  mat -> mat
void  mad_mat_cpym    (const  num_t *x,                         cnum_t *r, size_t m, size_t n, size_t ldx, size_t ldr); //  mat ->cmat
void  mad_mat_trans   (const  num_t *x,                          num_t *r, size_t m, size_t n);                         //  mat.t()
num_t mad_mat_dot     (const  num_t *x, const  num_t *y,                   size_t m, size_t n, size_t p);               // <mat ,  mat>
void  mad_mat_dotm    (const  num_t *x, const cnum_t *y,        cnum_t *r, size_t m, size_t n, size_t p);               // <mat , cmat>
void  mad_mat_mul     (const  num_t *x, const  num_t *y,         num_t *r, size_t m, size_t n, size_t p);               //  mat *  mat
void  mad_mat_mulm    (const  num_t *x, const cnum_t *y,        cnum_t *r, size_t m, size_t n, size_t p);               //  mat * cmat
int   mad_mat_invn    (const  num_t *y,        num_t  x,         num_t *r, size_t m, size_t n,           num_t rcond);  //  num /  mat
int   mad_mat_invc    (const  num_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t m, size_t n,           num_t rcond);  // cnum /  mat
int   mad_mat_div     (const  num_t *x, const  num_t *y,         num_t *r, size_t m, size_t n, size_t p, num_t rcond);  //  mat /  mat
int   mad_mat_divm    (const  num_t *x, const cnum_t *y,        cnum_t *r, size_t m, size_t n, size_t p, num_t rcond);  //  mat / cmat

void  mad_cmat_ident  (                                         cnum_t *r, size_t m, size_t n,             size_t ldr); //  ident->cmat
void  mad_cmat_set    (                num_t x_re, num_t x_im,  cnum_t *r, size_t m, size_t n,             size_t ldr); //  cnum ->cmat
void  mad_cmat_cpy    (const cnum_t *x,                         cnum_t *r, size_t m, size_t n, size_t ldx, size_t ldr); //  cmat ->cmat
void  mad_cmat_trans  (const cnum_t *x,                         cnum_t *r, size_t m, size_t n);                         //  cmat.t()
void  mad_cmat_ctrans (const cnum_t *x,                         cnum_t *r, size_t m, size_t n);                         //  cmat.ct()
void  mad_cmat_dot    (const cnum_t *x, const cnum_t *y,        cnum_t *r, size_t m, size_t n, size_t p);               // <cmat , cmat>
void  mad_cmat_dotm   (const cnum_t *x, const  num_t *y,        cnum_t *r, size_t m, size_t n, size_t p);               // <cmat ,  mat>
void  mad_cmat_mul    (const cnum_t *x, const cnum_t *y,        cnum_t *r, size_t m, size_t n, size_t p);               //  cmat * cmat
void  mad_cmat_mulm   (const cnum_t *x, const  num_t *y,        cnum_t *r, size_t m, size_t n, size_t p);               //  cmat *  mat
int   mad_cmat_invn   (const cnum_t *y,        num_t  x,        cnum_t *r, size_t m, size_t n,           num_t rcond);  //   num / cmat
int   mad_cmat_invc   (const cnum_t *y, num_t x_re, num_t x_im, cnum_t *r, size_t m, size_t n,           num_t rcond);  //  cnum / cmat
int   mad_cmat_div    (const cnum_t *x, const cnum_t *y,        cnum_t *r, size_t m, size_t n, size_t p, num_t rcond);  //  cmat / cmat
int   mad_cmat_divm   (const cnum_t *x, const  num_t *y,        cnum_t *r, size_t m, size_t n, size_t p, num_t rcond);  //  cmat /  mat
]]

return ffi.C