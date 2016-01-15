local ffi = require 'ffi'

-- types

ffi.cdef[[
typedef  const char*       str_t; // mad.h
typedef  int               idx_t; // mad.h
typedef  double            num_t; // mad.h
typedef  double _Complex  cnum_t; // mad.h

typedef  unsigned int      bit_t; // mad_bit.h
typedef  unsigned char     ord_t; // mad_mono.h

typedef  struct _IO_FILE    FILE; // stdio.h
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
void  mad_vec_set   (                         num_t x        ,  num_t  r[], size_t n); //  num -> vec
void  mad_vec_cpy   (const  num_t x[],                          num_t  r[], size_t n); //  vec -> vec
void  mad_vec_cpyv  (const  num_t x[],                         cnum_t  r[], size_t n); //  vec ->cvec
num_t mad_vec_dot   (const  num_t x[], const  num_t y[]      ,              size_t n); // <vec ,  vec>
void  mad_vec_dotv  (const  num_t x[], const cnum_t y[]      , cnum_t *r  , size_t n); // <vec , cvec>
void  mad_vec_add   (const  num_t x[], const  num_t y[]      ,  num_t  r[], size_t n); //  vec +  vec
void  mad_vec_addn  (const  num_t x[],        num_t y        ,  num_t  r[], size_t n); //  vec +  num
void  mad_vec_addc  (const  num_t x[], num_t y_re, num_t y_im, cnum_t  r[], size_t n); //  vec +  cpx
void  mad_vec_sub   (const  num_t x[], const  num_t y[]      ,  num_t  r[], size_t n); //  vec -  vec
void  mad_vec_subv  (const  num_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  vec - cvec
void  mad_vec_subn  (const  num_t y[],        num_t x        ,  num_t  r[], size_t n); //  num -  vec
void  mad_vec_subc  (const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cpx -  vec
void  mad_vec_mul   (const  num_t x[], const  num_t y[]      ,  num_t  r[], size_t n); //  vec *  vec
void  mad_vec_muln  (const  num_t x[],        num_t y        ,  num_t  r[], size_t n); //  vec *  num
void  mad_vec_mulc  (const  num_t x[], num_t y_re, num_t y_im, cnum_t  r[], size_t n); //  vec *  cpx
void  mad_vec_div   (const  num_t x[], const  num_t y[]      ,  num_t  r[], size_t n); //  vec /  vec
void  mad_vec_divv  (const  num_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  vec / cvec
void  mad_vec_divn  (const  num_t y[],        num_t x        ,  num_t  r[], size_t n); //  num /  vec
void  mad_vec_divc  (const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cpx /  vec 

void  mad_cvec_set  (                  num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cnum ->cvec
void  mad_cvec_cpy  (const cnum_t x[],                         cnum_t  r[], size_t n); //  cvec ->cvec
void  mad_cvec_cpyv (const cnum_t x[],                          num_t  r[], size_t n); //  cvec -> vec
void  mad_cvec_dot  (const cnum_t x[], const cnum_t y[]      , cnum_t *r  , size_t n); // <cvec , cvec>
void  mad_cvec_dotv (const cnum_t x[], const  num_t y[]      , cnum_t *r  , size_t n); // <cvec ,  vec>
void  mad_cvec_add  (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  cvec + cvec
void  mad_cvec_addv (const cnum_t x[], const  num_t y[]      , cnum_t  r[], size_t n); //  cvec +  vec
void  mad_cvec_addn (const cnum_t x[],        num_t y        , cnum_t  r[], size_t n); //  cvec +  num
void  mad_cvec_addc (const cnum_t x[], num_t y_re, num_t y_im, cnum_t  r[], size_t n); //  cvec +  cpx
void  mad_cvec_sub  (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  cvec - cvec
void  mad_cvec_subv (const cnum_t x[], const  num_t y[]      , cnum_t  r[], size_t n); //  cvec -  vec
void  mad_cvec_subn (const cnum_t y[],        num_t x        , cnum_t  r[], size_t n); //  num  - cvec
void  mad_cvec_subc (const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cpx  - cvec
void  mad_cvec_mul  (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  cvec * cvec
void  mad_cvec_mulv (const cnum_t x[], const  num_t y[]      , cnum_t  r[], size_t n); //  cvec *  vec
void  mad_cvec_muln (const cnum_t x[],        num_t y        , cnum_t  r[], size_t n); //  cvec *  num
void  mad_cvec_mulc (const cnum_t x[], num_t y_re, num_t y_im, cnum_t  r[], size_t n); //  cvec *  cpx
void  mad_cvec_div  (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  cvec / cvec
void  mad_cvec_divv (const cnum_t x[], const  num_t y[]      , cnum_t  r[], size_t n); //  cvec /  vec
void  mad_cvec_divn (const cnum_t y[],        num_t x        , cnum_t  r[], size_t n); //  num  / cvec
void  mad_cvec_divc (const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cpx  / cvec
]]

-- functions for matrix-matrix, vector-matrix and matrix-vector operations for LuaJIT (mad_mat.h)

ffi.cdef[[
void  mad_mat_ident   (                                           num_t  r[], size_t m, size_t n,             size_t ldr); // ident-> mat
void  mad_mat_set     (                         num_t x,          num_t  r[], size_t m, size_t n,             size_t ldr); //  num -> mat
void  mad_mat_cpy     (const  num_t x[],                          num_t  r[], size_t m, size_t n, size_t ldx, size_t ldr); //  mat -> mat
void  mad_mat_cpym    (const  num_t x[],                         cnum_t  r[], size_t m, size_t n, size_t ldx, size_t ldr); //  mat ->cmat
void  mad_mat_trans   (const  num_t x[],                          num_t  r[], size_t m, size_t n);                         //  mat.t()
num_t mad_mat_dot     (const  num_t x[], const  num_t y[],                    size_t m, size_t n, size_t p);               // <mat ,  mat>
void  mad_mat_dotm    (const  num_t x[], const cnum_t y[],       cnum_t *r  , size_t m, size_t n, size_t p);               // <mat , cmat>
void  mad_mat_mul     (const  num_t x[], const  num_t y[],        num_t  r[], size_t m, size_t n, size_t p);               //  mat *  mat
void  mad_mat_mulm    (const  num_t x[], const cnum_t y[],       cnum_t  r[], size_t m, size_t n, size_t p);               //  mat * cmat
int   mad_mat_invn    (const  num_t y[],        num_t x,          num_t  r[], size_t m, size_t n,           num_t rcond);  //  num /  mat
int   mad_mat_invc    (const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t m, size_t n,           num_t rcond);  // cnum /  mat
int   mad_mat_div     (const  num_t x[], const  num_t y[],        num_t  r[], size_t m, size_t n, size_t p, num_t rcond);  //  mat /  mat
int   mad_mat_divm    (const  num_t x[], const cnum_t y[],       cnum_t  r[], size_t m, size_t n, size_t p, num_t rcond);  //  mat / cmat

void  mad_cmat_ident  (                                          cnum_t  r[], size_t m, size_t n,             size_t ldr); //  ident->cmat
void  mad_cmat_set    (                  num_t x_re, num_t x_im, cnum_t  r[], size_t m, size_t n,             size_t ldr); //  cnum ->cmat
void  mad_cmat_cpy    (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n, size_t ldx, size_t ldr); //  cmat ->cmat
void  mad_cmat_trans  (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  cmat.t()
void  mad_cmat_ctrans (const cnum_t x[],                         cnum_t  r[], size_t m, size_t n);                         //  cmat.ct()
void  mad_cmat_dot    (const cnum_t x[], const cnum_t y[],       cnum_t *r  , size_t m, size_t n, size_t p);               // <cmat , cmat>
void  mad_cmat_dotm   (const cnum_t x[], const  num_t y[],       cnum_t *r  , size_t m, size_t n, size_t p);               // <cmat ,  mat>
void  mad_cmat_mul    (const cnum_t x[], const cnum_t y[],       cnum_t  r[], size_t m, size_t n, size_t p);               //  cmat * cmat
void  mad_cmat_mulm   (const cnum_t x[], const  num_t y[],       cnum_t  r[], size_t m, size_t n, size_t p);               //  cmat *  mat
int   mad_cmat_invn   (const cnum_t y[],        num_t x,         cnum_t  r[], size_t m, size_t n,           num_t rcond);  //   num / cmat
int   mad_cmat_invc   (const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t m, size_t n,           num_t rcond);  //  cnum / cmat
int   mad_cmat_div    (const cnum_t x[], const cnum_t y[],       cnum_t  r[], size_t m, size_t n, size_t p, num_t rcond);  //  cmat / cmat
int   mad_cmat_divm   (const cnum_t x[], const  num_t y[],       cnum_t  r[], size_t m, size_t n, size_t p, num_t rcond);  //  cmat /  mat
]]

-- functions for GTPSA descriptor (mad_desc.h)

ffi.cdef[[
typedef struct desc desc_t;

// globals
extern const ord_t mad_tpsa_default;
extern const ord_t mad_tpsa_same;
extern       int   mad_tpsa_strict;

// ctors, dtor
desc_t* mad_desc_new  (int nv, const ord_t var_ords[], const ord_t map_ords_[], str_t var_nam_[]);
desc_t* mad_desc_newk (int nv, const ord_t var_ords[], const ord_t map_ords_[], str_t var_nam_[],
                       int nk, const ord_t knb_ords[], ord_t dk); // knobs
void    mad_desc_del  (desc_t *d);

// introspection
int     mad_desc_maxsize (const desc_t *d);
ord_t   mad_desc_maxord  (const desc_t *d);
ord_t   mad_desc_gtrunc  (      desc_t *d, ord_t to);
]]

-- functions for real GTPSA (mad_tpsa.h)

ffi.cdef[[
typedef struct tpsa tpsa_t;

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
const ord_t*
        mad_tpsa_mono    (const tpsa_t *t, int i, int *n, ord_t *total_ord_);
int     mad_tpsa_midx    (const tpsa_t *t, int n, const ord_t m[]);
int     mad_tpsa_midx_sp (const tpsa_t *t, int n, const int   m[]); // sparse mono [(i,o)]

// accessors
num_t   mad_tpsa_get0    (const tpsa_t *t);
num_t   mad_tpsa_geti    (const tpsa_t *t, int i);
num_t   mad_tpsa_getm    (const tpsa_t *t, int n, const ord_t m[]);
num_t   mad_tpsa_getm_sp (const tpsa_t *t, int n, const int   m[]); // sparse mono [(i,o)]
void    mad_tpsa_set0    (      tpsa_t *t, /* i = 0 */             num_t a, num_t b);
void    mad_tpsa_seti    (      tpsa_t *t, int i,                  num_t a, num_t b);
void    mad_tpsa_setm    (      tpsa_t *t, int n, const ord_t m[], num_t a, num_t b);
void    mad_tpsa_setm_sp (      tpsa_t *t, int n, const int   m[], num_t a, num_t b);

// tranformations
tpsa_t* mad_tpsa_map     (const tpsa_t *a,                  tpsa_t *c, num_t (*f)(num_t v, int i_));
tpsa_t* mad_tpsa_map2    (const tpsa_t *a, const tpsa_t *b, tpsa_t *c, num_t (*f)(num_t va, num_t vb, int i_));

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

void    mad_tpsa_axpb       (num_t a, const tpsa_t *x, num_t b,                           tpsa_t *r);  // aliasing OK
void    mad_tpsa_axpbypc    (num_t a, const tpsa_t *x, num_t b, const tpsa_t *y, num_t c, tpsa_t *r);  // aliasing OK
void    mad_tpsa_axypb      (num_t a, const tpsa_t *x,          const tpsa_t *y, num_t b, tpsa_t *r);  // aliasing OK
void    mad_tpsa_axypbzpc   (num_t a, const tpsa_t *x,          const tpsa_t *y, num_t b,
                                                                const tpsa_t *z, num_t c, tpsa_t *r);  // aliasing OK
void    mad_tpsa_axypbvwpc  (num_t a, const tpsa_t *x,          const tpsa_t *y,
                             num_t b, const tpsa_t *v,          const tpsa_t *w, num_t c, tpsa_t *r);  // aliasing OK
void    mad_tpsa_ax2pby2pcz2(num_t a, const tpsa_t *x, num_t b, const tpsa_t *y, num_t c, const tpsa_t *z, tpsa_t *r); // aliasing OK

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
]]

-- functions for complex GTPSA (mad_ctpsa.h)

ffi.cdef[[
typedef struct ctpsa ctpsa_t;

// ctors, dtor
ctpsa_t* mad_ctpsa_newd    (desc_t *d, ord_t mo); // if mo > d_mo, mo = d_mo
ctpsa_t* mad_ctpsa_new     (const ctpsa_t *t, ord_t mo);
void     mad_ctpsa_del     (      ctpsa_t *t);
void     mad_ctpsa_delv    (      ctpsa_t *t1, ctpsa_t *t2, ...);

// introspection
desc_t*  mad_ctpsa_desc    (const ctpsa_t *t);
ord_t    mad_ctpsa_ord     (const ctpsa_t *t);
ord_t    mad_ctpsa_ordv    (const ctpsa_t *t1, const ctpsa_t *t2, ...);  // max order of all

// initialization
ctpsa_t* mad_ctpsa_copy    (const ctpsa_t *t, ctpsa_t *dst);
void     mad_ctpsa_clear   (      ctpsa_t *t);
void     mad_ctpsa_scalar  (      ctpsa_t *t, cnum_t v);

// indexing / monomials
const ord_t*
         mad_ctpsa_mono    (const ctpsa_t *t, int i, int *n, ord_t *total_ord_);
int      mad_ctpsa_midx    (const ctpsa_t *t, int n, const ord_t m[]);
int      mad_ctpsa_midx_sp (const ctpsa_t *t, int n, const int   m[]); // sparse mono [(i,o)]

// accessors
cnum_t   mad_ctpsa_get0    (const ctpsa_t *t);
cnum_t   mad_ctpsa_geti    (const ctpsa_t *t, int i);
cnum_t   mad_ctpsa_getm    (const ctpsa_t *t, int n, const ord_t m[]);
cnum_t   mad_ctpsa_getm_sp (const ctpsa_t *t, int n, const int   m[]); // sparse mono [(i,o)]
void     mad_ctpsa_set0    (      ctpsa_t *t, /* i = 0 */             cnum_t a, cnum_t b);
void     mad_ctpsa_seti    (      ctpsa_t *t, int i,                  cnum_t a, cnum_t b);
void     mad_ctpsa_setm    (      ctpsa_t *t, int n, const ord_t m[], cnum_t a, cnum_t b);
void     mad_ctpsa_setm_sp (      ctpsa_t *t, int n, const int   m[], cnum_t a, cnum_t b);

// tranformations
ctpsa_t* mad_ctpsa_map     (const ctpsa_t *a,                   ctpsa_t *c, cnum_t (*f)(cnum_t v, int i_));
ctpsa_t* mad_ctpsa_map2    (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c, cnum_t (*f)(cnum_t va, cnum_t vb, int i_));

// operations
void     mad_ctpsa_abs     (const ctpsa_t *a, ctpsa_t *c);
cnum_t   mad_ctpsa_nrm1    (const ctpsa_t *t, const ctpsa_t *t2_);
cnum_t   mad_ctpsa_nrm2    (const ctpsa_t *t, const ctpsa_t *t2_);
void     mad_ctpsa_der     (const ctpsa_t *a, ctpsa_t *c, int var);  // TODO: check functions that rely on it
void     mad_ctpsa_mder    (const ctpsa_t *a, ctpsa_t *c, int n, const ord_t m[]);

void     mad_ctpsa_add     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_sub     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_mul     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);
void     mad_ctpsa_div     (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c);

void     mad_ctpsa_acc     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c += v*a, aliasing OK
void     mad_ctpsa_scl     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v*a
void     mad_ctpsa_inv     (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v/a
void     mad_ctpsa_invsqrt (const ctpsa_t *a, cnum_t v, ctpsa_t *c);  // c  = v/sqrt(a)
void     mad_ctpsa_sqrt    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_exp     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_log     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sin     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cos     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sinh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cosh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sincos  (const ctpsa_t *a, ctpsa_t *s, ctpsa_t *c);
void     mad_ctpsa_sincosh (const ctpsa_t *a, ctpsa_t *s, ctpsa_t *c);
void     mad_ctpsa_sinc    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_sirx    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_corx    (const ctpsa_t *a, ctpsa_t *c);

void     mad_ctpsa_tan     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_cot     (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_asin    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acos    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_atan    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acot    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_tanh    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_coth    (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_asinh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acosh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_atanh   (const ctpsa_t *a, ctpsa_t *c);
void     mad_ctpsa_acoth   (const ctpsa_t *a, ctpsa_t *c);

void     mad_ctpsa_erf     (const ctpsa_t *a, ctpsa_t *c);

void     mad_ctpsa_ipow    (const ctpsa_t *a, ctpsa_t *c, int n);

void     mad_ctpsa_axpb       (cnum_t a, const ctpsa_t *x, cnum_t b,                             ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axpbypc    (cnum_t a, const ctpsa_t *x, cnum_t b, const ctpsa_t *y, cnum_t c, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypb      (cnum_t a, const ctpsa_t *x,           const ctpsa_t *y, cnum_t b, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypbzpc   (cnum_t a, const ctpsa_t *x,           const ctpsa_t *y, cnum_t b,
                                                                     const ctpsa_t *z, cnum_t c, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_axypbvwpc  (cnum_t a, const ctpsa_t *x,           const ctpsa_t *y,
                               cnum_t b, const ctpsa_t *v,           const ctpsa_t *w, cnum_t c, ctpsa_t *r);  // aliasing OK
void     mad_ctpsa_ax2pby2pcz2(cnum_t a, const ctpsa_t *x, cnum_t b, const ctpsa_t *y, cnum_t c, const ctpsa_t *z, ctpsa_t *r); // aliasing OK

// to check for non-homogeneous maps & knobs
void     mad_ctpsa_poisson (const ctpsa_t *a, const ctpsa_t *b, ctpsa_t *c, int n);  // TO CHECK n
void     mad_ctpsa_compose (int sa, const ctpsa_t *ma[], int sb, const ctpsa_t *mb[], int sc, ctpsa_t *mc[]);
void     mad_ctpsa_minv    (int sa, const ctpsa_t *ma[],                              int sc, ctpsa_t *mc[]);
void     mad_ctpsa_pminv   (int sa, const ctpsa_t *ma[],                              int sc, ctpsa_t *mc[], int row_select[]);

// I/O
void     mad_ctpsa_print    (const ctpsa_t *t, str_t name_, FILE *stream_);
ctpsa_t* mad_ctpsa_scan     (                               FILE *stream_); // TODO
desc_t*  mad_ctpsa_scan_hdr (                               FILE *stream_);
void     mad_ctpsa_scan_coef(      ctpsa_t *t,              FILE *stream_); // TODO
void     mad_ctpsa_debug    (const ctpsa_t *t);
]]

return ffi.C