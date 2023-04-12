#ifndef MAD_MONO_H
#define MAD_MONO_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Monomial module interface
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: C. Tomoiaga
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - provide some functions to handle monomials as array of ord_t

  Information:
  - parameters ending with an underscope can be null.

  Errors:
  - TODO

 o-----------------------------------------------------------------------------o
 */

#include <stdio.h>

#include "mad_def.h"

// --- types ------------------------------------------------------------------o

typedef unsigned char ord_t;

// --- interface --------------------------------------------------------------o

ssz_t mad_mono_str   (ssz_t n,       ord_t a[n], str_t s);     // from string
str_t mad_mono_prt   (ssz_t n, const ord_t a[n], char s[n+1]); // to   string

void  mad_mono_fill  (ssz_t n,       ord_t a[n], ord_t v);
void  mad_mono_copy  (ssz_t n, const ord_t a[n], ord_t r[n]);

ord_t mad_mono_min   (ssz_t n, const ord_t a[n]);
ord_t mad_mono_max   (ssz_t n, const ord_t a[n]);
int   mad_mono_ord   (ssz_t n, const ord_t a[n]); // sum of orders
num_t mad_mono_ordp  (ssz_t n, const ord_t a[n], idx_t stp); // product of orders
num_t mad_mono_ordpf (ssz_t n, const ord_t a[n], idx_t stp); // product of factorials of orders

log_t mad_mono_eq    (ssz_t n, const ord_t a[n], const ord_t b[n]);
log_t mad_mono_lt    (ssz_t n, const ord_t a[n], const ord_t b[n]);
log_t mad_mono_le    (ssz_t n, const ord_t a[n], const ord_t b[n]);

int   mad_mono_cmp   (ssz_t n, const ord_t a[n], const ord_t b[n]);
int   mad_mono_rcmp  (ssz_t n, const ord_t a[n], const ord_t b[n]);

void  mad_mono_add   (ssz_t n, const ord_t a[n], const ord_t b[n], ord_t r[n]);
void  mad_mono_sub   (ssz_t n, const ord_t a[n], const ord_t b[n], ord_t r[n]);
void  mad_mono_cat   (ssz_t n, const ord_t a[n], ssz_t m, const ord_t b[m], ord_t r[n+m]);
void  mad_mono_rev   (ssz_t n, const ord_t a[n], ord_t r[n]);

void  mad_mono_print (ssz_t n, const ord_t a[n], FILE *fp_);

// --- end --------------------------------------------------------------------o

#endif //  MAD_MONO_H
