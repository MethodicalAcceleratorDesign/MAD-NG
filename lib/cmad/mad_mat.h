#ifndef MAD_MAT_H
#define MAD_MAT_H

#include <stddef.h>

#define  num_t double
#define cnum_t double _Complex

void mad_mat_trans  (const  num_t *x,                   num_t *r, size_t m, size_t n);           //  mat.t()
void mad_mat_mul    (const  num_t *x, const  num_t *y,  num_t *r, size_t m, size_t n, size_t p); //  mat *  mat
void mad_mat_mulm   (const  num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p); //  mat * cmat
void mad_mat_mulv   (const  num_t *x, const  num_t *y,  num_t *r, size_t m, size_t p);           //  mat *  vec
void mad_mat_mulc   (const  num_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t p);           //  mat * cvec
void mad_mat_vmul   (const  num_t *x, const  num_t *y,  num_t *r, size_t n, size_t p);           //  vec *  mat
void mad_mat_cmul   (const cnum_t *x, const  num_t *y, cnum_t *r, size_t n, size_t p);           // cvec *  mat

void mad_cmat_trans (const cnum_t *x,                  cnum_t *r, size_t m, size_t n);           // cmat.t()
void mad_cmat_mul   (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t n, size_t p); // cmat * cmat
void mad_cmat_mulm  (const cnum_t *x, const  num_t *y, cnum_t *r, size_t m, size_t n, size_t p); // cmat *  mat
void mad_cmat_mulv  (const cnum_t *x, const  num_t *y, cnum_t *r, size_t m, size_t p);           // cmat *  vec
void mad_cmat_mulc  (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t m, size_t p);           // cmat * cvec
void mad_cmat_vmul  (const  num_t *x, const cnum_t *y, cnum_t *r, size_t n, size_t p);           //  vec * cmat
void mad_cmat_cmul  (const cnum_t *x, const cnum_t *y, cnum_t *r, size_t n, size_t p);           // cvec * cmat

#undef  num_t
#undef cnum_t

#endif // MAD_MAT_H