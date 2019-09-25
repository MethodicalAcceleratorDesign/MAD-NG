/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA I/O module implementation
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 |          C. Tomoiaga
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o
*/

#include <math.h>
#include <string.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#define  SPC "                      "
#else
#include "mad_tpsa_impl.h"
#define  SPC
#endif

// --- local ------------------------------------------------------------------o

static inline void
print_ords(int n, const ord_t ords[n], FILE *stream)
{
  assert(ords && stream);
  for (int i=0; i < n-1; i += 2)
    fprintf(stream, "  %hhu %hhu", ords[i], ords[i+1]);
  if (n % 2)
    fprintf(stream, "  %hhu"     , ords[n-1]);
}

static inline void
read_ords(int n, ord_t ords[n], FILE *stream)
{
  assert(ords && stream);
  for (int i=0; i < n; ++i) {
    int cnt = fscanf(stream, "%hhu", &ords[i]);
    ensure(cnt == 1, "invalid input (missing order?)");
  }
}

// --- public -----------------------------------------------------------------o

#ifdef MAD_CTPSA_IMPL

extern const D* mad_tpsa_scan_hdr(FILE*);

const D*
FUN(scan_hdr) (FILE *stream_)
{
  return mad_tpsa_scan_hdr(stream_);
}

#else

const D*
FUN(scan_hdr) (FILE *stream_)
{
  enum { BUF_SIZE=256 };
  char buf[BUF_SIZE];

  int nmv=0, nk=0, cnt=0;
  ord_t mo, ko;

  if (!stream_) stream_ = stdin;

  // discard leading white space and the name (which is 10 chars and comma)
  ensure(!fscanf(stream_, " %*11c"), "unexpected fscanf returned value");
  ensure(!feof(stream_) && !ferror(stream_), "invalid input (file error?)");

  // 1st line
  cnt = fscanf(stream_, " NO =%5hhu, NV =%5d, KO =%5hhu, NK =%5d",
                              &mo,       &nmv,    &ko,       &nk);

  if (cnt == 2) {
    // TPSA  -- ignore rest of lines; default values
    ensure(fgets(buf, BUF_SIZE, stream_), "invalid input (file error?)"); // finish  1st line
    ensure(fgets(buf, BUF_SIZE, stream_), "invalid input (file error?)"); // discard 2nd line
    ensure(fgets(buf, BUF_SIZE, stream_), "invalid input (file error?)"); // discard 3rd line
    ensure(fgets(buf, BUF_SIZE, stream_), "invalid input (file error?)"); // discard coeff header

    ord_t mvar_ords[nmv];
    mad_mono_fill(nmv, mvar_ords, mo);
    return mad_desc_newm(nmv, mvar_ords);
  }

  if (cnt == 4) {
    // GTPSA -- process rest of lines
    int nv = nmv+nk;
    ord_t mvar_ords [nmv], var_ords[nv];

    // read mvars orders
    ensure(!fscanf(stream_, " MAP ORDS: "), "unexpected fscanf returned value");
    ensure(!feof(stream_) && !ferror(stream_), "invalid input (file error?)");
    read_ords(nmv, mvar_ords, stream_);

    // read var orders
    ensure(!fscanf(stream_, " ||| VAR ORDS: "), "unexpected fscanf returned value");
    ensure(!feof(stream_) && !ferror(stream_), "invalid input (file error?)");
    read_ords(nv, var_ords, stream_);

    return mad_desc_newv(nmv, mvar_ords, nv, var_ords, ko);
  }

  if (cnt <  2) error("could not read (NO,NV) from header");
  if (cnt == 3) error("could not read NK from header");
  return NULL; // never reached
}

#endif // !MAD_CTPSA_IMPL

void
FUN(scan_coef) (T *t, FILE *stream_)
{
  assert(t);
  if (!stream_) stream_ = stdin;

  NUM c;
  int nv = t->d->nv, cnt = -1;
  ord_t o, ords[nv];
  FUN(reset0)(t);

#ifndef MAD_CTPSA_IMPL
  while ((cnt = fscanf(stream_, "%*d %lG %hhu", &c, &o)) == 2) {
#else
  while ((cnt = fscanf(stream_, "%*d %lG%lGi %hhu", (num_t*)&c, (num_t*)&c+1, &o)) == 2) {
#endif

    #ifdef DEBUG
      printf("c=" FMT ", o=%d\n", VAL(c), o);
    #endif
    read_ords(nv,ords,stream_); // sanity check
    ensure(mad_mono_ord(nv,ords) == o, "invalid input (bad order?)");
    if (o <= t->mo)             // discard too high mononial
     FUN(setm)(t,nv,ords, 0.0,c);
  }
}

T*
FUN(scan) (FILE *stream_)
{
  const D *d = FUN(scan_hdr)(stream_);
  T *t = FUN(newd)(d, mad_tpsa_default);
  FUN(scan_coef)(t, stream_);
  return t;
}

void
FUN(print) (const T *t, str_t name_, num_t eps_, FILE *stream_)
{
  assert(t);

  if (!name_  ) name_ = "-UNNAMED--";
  if (eps_ < 0) eps_ = 1e-16;
  if (!stream_) stream_ = stdout;

  const D *d = t->d;

  // print header
  fprintf(stream_, "\n %10s, NO =%5hhu, NV =%5d, KO =%5hhu, NK =%5d\n MAP ORDS:",
                       name_,    d->mo,  d->nmv,     d->ko, d->nv - d->nmv);
  print_ords(d->nmv, d->mvar_ords, stream_);
  fprintf(stream_, " ||| VAR ORDS: ");
  print_ords(d->nv, d->var_ords, stream_);
  fprintf(stream_, "\n *******************************************************");

  // print coefficients
  fprintf(stream_, "\n     I   COEFFICIENT         " SPC "  ORDER   EXPONENTS");
  idx_t *pi = d->ord2idx, idx = 0;
  for (ord_t o = t->lo; o <= t->hi ; ++o) {
    if (!mad_bit_get(t->nz,o)) continue;
    for (idx_t i = pi[o]; i < pi[o+1]; ++i) {
      if (fabs(t->coef[i]) < eps_) continue;
#ifndef MAD_CTPSA_IMPL
      fprintf(stream_, "\n%6d  %21.14lE%5hhu   "          , ++idx, VAL(t->coef[i]), d->ords[i]);
#else
      fprintf(stream_, "\n%6d  %21.14lE%+21.14lEi%5hhu   ", ++idx, VAL(t->coef[i]), d->ords[i]);
#endif
      print_ords(d->nv, d->To[i], stream_);
    }
  }

  if (!idx)
    fprintf(stream_, "\n          ALL COMPONENTS ZERO \n");
  else
    fprintf(stream_, "\n\n");
}

// --- end --------------------------------------------------------------------o
