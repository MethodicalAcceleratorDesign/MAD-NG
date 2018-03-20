/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA I/O module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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

extern D* mad_tpsa_scan_hdr(FILE*);

D*
FUN(scan_hdr) (FILE *stream_)
{
  return mad_tpsa_scan_hdr(stream_);
}

#else

D*
FUN(scan_hdr) (FILE *stream_)
{
  enum { BUF_SIZE=256 };
  char buf[BUF_SIZE];

  int nmv=0, nk=0, cnt=0;
  ord_t mo, ko;

  if (!stream_) stream_ = stdin;

  // discard leading white space and the name (which is 10 chars and comma)
  fscanf(stream_, " %*11c");
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
    fscanf(stream_, " MAP ORDS: ");
    ensure(!feof(stream_) && !ferror(stream_), "invalid input (file error?)");
    read_ords(nmv, mvar_ords, stream_);

    // read var orders
    fscanf(stream_, " ||| VAR ORDS: ");
    ensure(!feof(stream_) && !ferror(stream_), "invalid input (file error?)");
    read_ords(nv, var_ords, stream_);

    desc_t *d = mad_desc_newv(nmv, mvar_ords, nv, var_ords, ko);
    return d;
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
  FUN(clear)(t);

#ifndef MAD_CTPSA_IMPL
  while ((cnt = fscanf(stream_, "%*d %lG %hhu", &c, &o)) == 2) {
#else
  while ((cnt = fscanf(stream_, "%*d %lG%lGi %hhu", (num_t*)&c, (num_t*)&c+1, &o)) == 2) {
#endif

    #ifdef DEBUG
      printf("c=%.2f o=%d\n", c, o);
    #endif
    read_ords(nv,ords,stream_); // sanity check
    ensure(mad_mono_ord(nv,ords) == o, "invalid input (bad order?)");
    if (o <= t->mo) FUN(setm)(t,nv,ords, 0.0,c); // discard too high mononial
  }
}

T*
FUN(scan) (FILE *stream_)
{
  desc_t *d = FUN(scan_hdr)(stream_);
  T *t = FUN(newd)(d, mad_tpsa_default);
  FUN(scan_coef)(t, stream_);
  return t;
}

void
FUN(print) (const T *t, str_t name_, FILE *stream_)
{
  assert(t);
  // TODO: print map vars and name

  if (!stream_) stream_ = stdout;
  D *d = t->d;

  // print header
  if (!name_) name_ = "-UNNAMED--";
  fprintf(stream_, "\n %10s, NO =%5hhu, NV =%5d, KO =%5hhu, NK =%5d\n MAP ORDS:",
                       name_,    d->mo,  d->nmv,     d->ko, d->nv - d->nmv);
  print_ords(d->nmv, d->mvar_ords, stream_);
  fprintf(stream_, " ||| VAR ORDS: ");
  print_ords(d->nv, d->var_ords, stream_);
  fprintf(stream_, "\n *******************************************************");

  if (!t->nz) {
    fprintf(stream_, "\n   ALL COMPONENTS ZERO \n");
    return;
  }

  fprintf(stream_, "\n     I   COEFFICIENT         " SPC "  ORDER   EXPONENTS");
  int idx = 1;
  ssz_t nc = mad_desc_tpsa_len(d, t->mo);
  for (int c = 0; c < nc; ++c) {
    if (mad_bit_get(t->nz,d->ords[c]) && fabs(t->coef[c]) > 1e-10) {
#ifndef MAD_CTPSA_IMPL
      fprintf(stream_, "\n%6d  %21.14lE%5hhu   "          , idx, VAL(t->coef[c]), d->ords[c]);
#else
      fprintf(stream_, "\n%6d  %21.14lE%+21.14lEi%5hhu   ", idx, VAL(t->coef[c]), d->ords[c]);
#endif
      print_ords(d->nv, d->To[c], stream_);
      idx++;
    }
  }
  fprintf(stream_, "\n\n");
}

// --- end --------------------------------------------------------------------o
