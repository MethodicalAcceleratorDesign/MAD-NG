#ifndef MAD_VEC_H
#define MAD_VEC_H

/*
 o----------------------------------------------------------------------------o
 |
 | Vector module interface
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
  
  Purpose:
  - wrappers around functions of real and complex vectors for LuaJIT

 o----------------------------------------------------------------------------o
 */

#include "mad.h"

// --- interface -------------------------------------------------------------o

void   mad_vec_set   (                         num_t x        ,  num_t  r[], size_t n); //  num -> vec
void   mad_vec_cpy   (const  num_t x[],                          num_t  r[], size_t n); //  vec -> vec
void   mad_vec_cpyv  (const  num_t x[],                         cnum_t  r[], size_t n); //  vec ->cvec
num_t  mad_vec_dot   (const  num_t x[], const  num_t y[]                   , size_t n); // <vec ,  vec>
cnum_t mad_vec_dotv  (const  num_t x[], const cnum_t y[]                   , size_t n); // <vec , cvec>
void   mad_vec_dotv_r(const  num_t x[], const cnum_t y[]      , cnum_t *r  , size_t n); // <vec , cvec>
void   mad_vec_add   (const  num_t x[], const  num_t y[]      ,  num_t  r[], size_t n); //  vec +  vec
void   mad_vec_addn  (const  num_t x[],        num_t y        ,  num_t  r[], size_t n); //  vec +  num
void   mad_vec_addc  (const  num_t x[],       cnum_t y        , cnum_t  r[], size_t n); //  vec +  cpx
void   mad_vec_addc_r(const  num_t x[], num_t y_re, num_t y_im, cnum_t  r[], size_t n); //  vec +  cpx
void   mad_vec_sub   (const  num_t x[], const  num_t y[]      ,  num_t  r[], size_t n); //  vec -  vec
void   mad_vec_subv  (const  num_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  vec - cvec
void   mad_vec_subn  (const  num_t y[],        num_t x        ,  num_t  r[], size_t n); //  num -  vec
void   mad_vec_subc  (const  num_t y[],       cnum_t x        , cnum_t  r[], size_t n); //  cpx -  vec
void   mad_vec_subc_r(const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cpx -  vec
void   mad_vec_mul   (const  num_t x[], const  num_t y[]      ,  num_t  r[], size_t n); //  vec *  vec
void   mad_vec_muln  (const  num_t x[],        num_t y        ,  num_t  r[], size_t n); //  vec *  num
void   mad_vec_mulc  (const  num_t x[],       cnum_t y        , cnum_t  r[], size_t n); //  vec *  cpx
void   mad_vec_mulc_r(const  num_t x[], num_t y_re, num_t y_im, cnum_t  r[], size_t n); //  vec *  cpx
void   mad_vec_div   (const  num_t x[], const  num_t y[]      ,  num_t  r[], size_t n); //  vec /  vec
void   mad_vec_divv  (const  num_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  vec / cvec
void   mad_vec_divn  (const  num_t y[],        num_t x        ,  num_t  r[], size_t n); //  num /  vec
void   mad_vec_divc  (const  num_t y[],       cnum_t x        , cnum_t  r[], size_t n); //  cpx /  vec 
void   mad_vec_divc_r(const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cpx /  vec 

void   mad_cvec_set   (                        cnum_t x        , cnum_t  r[], size_t n); //  cnum ->cvec
void   mad_cvec_set_r (                  num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cnum ->cvec
void   mad_cvec_cpy   (const cnum_t x[],                         cnum_t  r[], size_t n); //  cvec ->cvec
void   mad_cvec_cpyv  (const cnum_t x[],                          num_t  r[], size_t n); //  cvec -> vec
cnum_t mad_cvec_dot   (const cnum_t x[], const cnum_t y[]                   , size_t n); // <cvec , cvec>
cnum_t mad_cvec_dotv  (const cnum_t x[], const  num_t y[]                   , size_t n); // <cvec ,  vec>
void   mad_cvec_dot_r (const cnum_t x[], const cnum_t y[]      , cnum_t *r  , size_t n); // <cvec , cvec>
void   mad_cvec_dotv_r(const cnum_t x[], const  num_t y[]      , cnum_t *r  , size_t n); // <cvec ,  vec>
void   mad_cvec_add   (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  cvec + cvec
void   mad_cvec_addv  (const cnum_t x[], const  num_t y[]      , cnum_t  r[], size_t n); //  cvec +  vec
void   mad_cvec_addn  (const cnum_t x[],        num_t y        , cnum_t  r[], size_t n); //  cvec +  num
void   mad_cvec_addc  (const cnum_t x[],       cnum_t y        , cnum_t  r[], size_t n); //  cvec +  cpx
void   mad_cvec_addc_r(const cnum_t x[], num_t y_re, num_t y_im, cnum_t  r[], size_t n); //  cvec +  cpx
void   mad_cvec_sub   (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  cvec - cvec
void   mad_cvec_subv  (const cnum_t x[], const  num_t y[]      , cnum_t  r[], size_t n); //  cvec -  vec
void   mad_cvec_subn  (const cnum_t y[],        num_t x        , cnum_t  r[], size_t n); //  num  - cvec
void   mad_cvec_subc  (const cnum_t y[],       cnum_t x        , cnum_t  r[], size_t n); //  cpx  - cvec
void   mad_cvec_subc_r(const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cpx  - cvec
void   mad_cvec_mul   (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  cvec * cvec
void   mad_cvec_mulv  (const cnum_t x[], const  num_t y[]      , cnum_t  r[], size_t n); //  cvec *  vec
void   mad_cvec_muln  (const cnum_t x[],        num_t y        , cnum_t  r[], size_t n); //  cvec *  num
void   mad_cvec_mulc  (const cnum_t x[],       cnum_t y        , cnum_t  r[], size_t n); //  cvec *  cpx
void   mad_cvec_mulc_r(const cnum_t x[], num_t y_re, num_t y_im, cnum_t  r[], size_t n); //  cvec *  cpx
void   mad_cvec_div   (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], size_t n); //  cvec / cvec
void   mad_cvec_divv  (const cnum_t x[], const  num_t y[]      , cnum_t  r[], size_t n); //  cvec /  vec
void   mad_cvec_divn  (const cnum_t y[],        num_t x        , cnum_t  r[], size_t n); //  num  / cvec
void   mad_cvec_divc  (const cnum_t y[],       cnum_t x        , cnum_t  r[], size_t n); //  cpx  / cvec
void   mad_cvec_divc_r(const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], size_t n); //  cpx  / cvec

// ---------------------------------------------------------------------------o

#endif // MAD_VEC_H