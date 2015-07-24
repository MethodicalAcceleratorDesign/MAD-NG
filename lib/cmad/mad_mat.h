#ifndef MAD_MAT_H
#define MAD_MAT_H

#include <stddef.h>

#define  num_t double
#define cnum_t double _Complex

void  mad_mat_trans   (const  num_t *x,                   num_t *r, size_t m, size_t n);           //  mat.t()
num_t mad_mat_dot     (const  num_t *x, const  num_t *y,            size_t m, size_t n, size_t p); // <mat ,  mat>
void  mad_mat_dotm    (const  num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p); // <mat , cmat>
void  mad_mat_mul     (const  num_t *x, const  num_t *y,  num_t *r, size_t m, size_t n, size_t p); //  mat *  mat
void  mad_mat_mulm    (const  num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p); //  mat * cmat

void  mad_cmat_trans  (const cnum_t *x,                  cnum_t *r, size_t m, size_t n);           //  cmat.t()
void  mad_cmat_ctrans (const cnum_t *x,                  cnum_t *r, size_t m, size_t n);           //  cmat.ct()
void  mad_cmat_dot    (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p); // <cmat , cmat>
void  mad_cmat_dotm   (const cnum_t *x, const  num_t *y, cnum_t *r, size_t m, size_t n, size_t p); // <cmat ,  mat>
void  mad_cmat_mul    (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p); //  cmat * cmat
void  mad_cmat_mulm   (const cnum_t *x, const  num_t *y, cnum_t *r, size_t m, size_t n, size_t p); //  cmat *  mat

#undef  num_t
#undef cnum_t

#endif // MAD_MAT_H