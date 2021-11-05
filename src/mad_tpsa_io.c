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
#include <ctype.h>
#include <string.h>
#include <assert.h>

#include "mad_mem.h"
#include "mad_desc_impl.h"

#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#define  SPC "                       "
#else
#include "mad_tpsa_impl.h"
#define  SPC
#endif

//#undef  DEBUG
//#define DEBUG 3

// --- local ------------------------------------------------------------------o

static inline int
skip_line(FILE *stream)
{
  int c;
  while ((c = fgetc(stream)) != '\n' && c != EOF) ;
  return c;
}

static inline void
print_ords_sm(int n, const ord_t ords[n], FILE *stream)
{
  assert(ords && stream);
  for (int i=0; i < n; i++)
    if (ords[i]) fprintf(stream, "  %d^%hhu", i+1, ords[i]);
}

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
read_ords(int ci, str_t name, int n, ord_t ords[n], FILE *stream)
{
  assert(ords && stream);
  idx_t idx;
  ord_t ord;
  char  chr;

  if (!name[0]) name = "-UNNAMED-";

  mad_mono_fill(n, ords, 0);
  for (int i=0; i < n; i++) {
    idx = chr = ord = 0;
    int cnt = fscanf(stream, " %d%c%hhu", &idx, &chr, &ord);

#if DEBUG > 2
    printf("mono: ci=%d, i=%d[%d], cnt=%d, idx=%d, chr='%c', ord=%d\n",
           ci, i, n, cnt, idx, chr, ord);
#endif

    if (cnt == 3 && chr == '^') {
      ensure(0 < idx && idx <= n,
             "invalid index (expecting 0 < %d <= %d) at index %d of '%s'",
              idx, n, ci, name);
      ords[idx-1] = ord, i = idx-1;
    } else
    if (cnt == 3 && chr == ' ') {
      ords[i] = idx, ords[++i] = ord;
    } else
    if (cnt == 2 && chr == '\n') {
      ords[i] = idx; ungetc(chr, stream);
    } else
      error("invalid monomial input at index %d of '%s'", ci, name);
  }
}

// --- public -----------------------------------------------------------------o

#ifdef MAD_CTPSA_IMPL

extern const D* mad_tpsa_scan_hdr(int*, char[NAMSZ], FILE*);

const D*
FUN(scan_hdr) (int *kind_, char name_[NAMSZ], FILE *stream_)
{
  DBGFUN(->); // complex and real header are the same...
  const D* ret = mad_tpsa_scan_hdr(kind_, name_, stream_);
  DBGFUN(<-);
  return ret;
}

#else

const D*
FUN(scan_hdr) (int *kind_, char name_[NAMSZ], FILE *stream_)
{
  DBGFUN(->);
  if (!stream_) stream_ = stdin;

  // backup stream position
  fpos_t fpos;
  fgetpos(stream_, &fpos);

  // eat white space
  int c;
  while ((c=getc(stream_)) != EOF && isspace(c)) ;
  ungetc(c, stream_);

  // check the name (which is 10 chars) and the type
  char name[NAMSZ]="", knd='?';
  int nc=0, n = fscanf(stream_, "%16[^ \t:]%*[ ]: %c%n", name, &knd, &nc);
  name[NAMSZ-1]=0;

#if DEBUG > 2
    printf("header: n=%d, name='%s', knd='%c', nc=%d\n", n,name,knd,nc);
#endif

  if (n != 2
      || nc < 4 || !strchr(" RC", knd)
      || (kind_ && *kind_ != -1 && (*kind_ != (knd == 'C'))) ) {

         if (nc < 4)              strncpy(name, "-INVALIDNAME-", NAMSZ);
    else if (!strchr(" RC", knd)) strncpy(name, "-INVALIDTYPE-", NAMSZ);
    else if (kind_ && *kind_ != -1 && (*kind_ != (knd == 'C')))
                                  strncpy(name, "-UNXPCTDTYPE-", NAMSZ);

    if (name_) strncpy(name_, name, NAMSZ-1), name_[NAMSZ-1] = '\0';

    warn("unable to parse GTPSA header: %s", name);

    fsetpos(stream_, &fpos); // may fail for non-seekable stream (e.g. pipes)...
    return NULL;
  }

  ensure(!feof(stream_) && !ferror(stream_), "invalid input (file error?)");

  if (kind_) *kind_ = knd == 'C';
  if (name_) strncpy(name_, name, NAMSZ-1), name_[NAMSZ-1] = '\0';

  ord_t mo=0, po=0;
  int nv=0, np=0, cnt=0;

  // 1st line (cnt includes knd)
  cnt = 1+fscanf(stream_, ", NV = %d, MO = %hhu, NP = %d, PO = %hhu%n",
                                  &nv,     &mo,       &np,     &po,    &nc);

#if DEBUG > 2
    printf("header: cnt=%d, nv=%d, mo=%d, np=%d, po=%d\n", cnt,nv,mo,np,po);
#endif

  // sanity checks
  ensure(0 <  nv && nv <= DESC_MAX_VAR, "invalid NV=%d", nv);
  ensure(           mo <= DESC_MAX_ORD, "invalid MO=%d", mo);

  if (cnt == 3) {
    // TPSA -- ignore rest of lines
    ensure(skip_line(stream_) != EOF, "invalid input (file error?)"); // finish NV,MO line
    ensure(fscanf(stream_, "%*[*]\n") != 1, "unexpected input (invalid header?)");
    ensure(skip_line(stream_) != EOF, "invalid input (file error?)"); // discard coeff header

    const D* ret = mad_desc_newv(nv, mo);
    DBGFUN(<-);
    return ret;
  }

  if (cnt == 5) {
    int nn = nv+np;

    // sanity checks
    ensure(0 <= np && nn <= DESC_MAX_VAR, "invalid NP=%d", np);
    ensure(           po <= DESC_MAX_ORD, "invalid PO=%d", po);

    // GTPSA -- process rest of lines
    ord_t no[nn];

    // read variables orders if present
    if ((cnt += fscanf(stream_, ", N%*[O] = ")) == 6) {
      read_ords(-1,name,nn,no,stream_);
      ensure(skip_line(stream_) != EOF, "invalid input (file error?)"); // finish VO line
    }

    ensure(fscanf(stream_, "%*[*]\n") != 1, "unexpected input (invalid header?)");
    ensure(skip_line(stream_) != EOF, "invalid input (file error?)"); // discard coeff header
    ensure(!feof(stream_) && !ferror(stream_), "invalid input (file error?)");

    const D* ret = cnt == 5 ? mad_desc_newvp (nv, np, mo, po)
                            : mad_desc_newvpo(nv, np, no, po);
    DBGFUN(<-);
    return ret;
  }

       if (cnt < 3) warn("could not read (NV,NO) from header");
  else if (cnt < 5) warn("could not read (NP,PO) from header");
  else              warn("unable to parse GTPSA header for '%s'",
                         name[0] ? name : "-UNNAMED-");

  fsetpos(stream_, &fpos); // may fail for non-seekable stream (e.g. sockets)...
  return NULL;
}

#endif // !MAD_CTPSA_IMPL

void
FUN(scan_coef) (T *t, FILE *stream_)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  if (!stream_) stream_ = stdin;

  NUM c;
  int nn = t->d->nn, cnt = -1, i = -1;
  ord_t o, ords[nn];
  FUN(reset0)(t);

#ifndef MAD_CTPSA_IMPL
  while ((cnt = fscanf(stream_, "%d %lG %hhu", &i, &c, &o)) == 3) {
#else
  while ((cnt = fscanf(stream_, "%d %lG%lGi %hhu", &i, (num_t*)&c, (num_t*)&c+1, &o)) == 4) {
#endif

    #if DEBUG > 2
      printf("coef: i=%d, c=" FMT ", o=%d\n", i, VAL(c), o);
    #endif
    read_ords(i,t->nam,nn,ords,stream_); // sanity check
    ensure(mad_mono_ord(nn,ords) == o,
           "invalid monomial order at index %d of '%s'", i, t->nam);
    // discard too high mononial
    if (o <= t->mo) FUN(setm)(t,nn,ords,0,c);
  }
  #if DEBUG > 2
    printf("coef: i=%d, cnt=%d, c=" FMT ", o=%d\n", i, cnt, VAL(c), o);
  #endif
  FUN(update0)(t, t->lo, t->hi);
  DBGTPSA(t); DBGFUN(<-);
}

T*
FUN(scan) (FILE *stream_)
{
  DBGFUN(->);
#ifndef MAD_CTPSA_IMPL
  int knd = 0;
#else
  int knd = 1;
#endif
  T *t = NULL;
  char name[NAMSZ];
  const D *d = FUN(scan_hdr)(&knd, name, stream_);
  if (d) {
    t = FUN(newd)(d, mad_tpsa_default);
    FUN(scan_coef)(t, stream_);
    FUN(setnam)   (t, name   );
  }
  DBGFUN(<-);
  return t;
}

void
FUN(print) (const T *t, str_t name_, num_t eps_, int nohdr_, FILE *stream_)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  if (!name_  ) name_   = t->nam[0] ? t->nam : "-UNNAMED-";
  if (eps_ < 0) eps_    = 1e-16;
  if (!stream_) stream_ = stdout;

#ifndef MAD_CTPSA_IMPL
  const char typ = 'R';
#else
  const char typ = 'C';
#endif

  const D *d = t->d;

  if (nohdr_) goto coeffonly;

  // print header
  fprintf(stream_, d->np || d->uno
                 ? "\n %-8s:  %c, NV = %3d, MO = %2hhu, NP = %3d, PO = %2hhu"
                 : "\n %-8s:  %c, NV = %3d, MO = %2hhu",
                      name_, typ,    d->nn,      d->mo,    d->np,      d->po);

  if (d->uno) {
    fprintf(stream_, ", NO = ");
    print_ords(d->nn, d->no, stream_);
  }
  fprintf(stream_, "\n********************************************************");
#ifdef MAD_CTPSA_IMPL
  fprintf(stream_, "***********************");
#endif

coeffonly:

  // print coefficients
  fprintf(stream_, "\n     I   COEFFICIENT         " SPC "  ORDER   EXPONENTS");
  const idx_t *o2i = d->ord2idx;
  idx_t idx = 0;
  for (ord_t o = t->lo; o <= t->hi ; ++o) {
    if (!mad_bit_tst(t->nz,o)) continue;
    for (idx_t i = o2i[o]; i < o2i[o+1]; ++i) {
#ifndef MAD_CTPSA_IMPL
      if (fabs(t->coef[i]) < eps_) continue;
      fprintf(stream_, "\n%6d  %21.14lE   %2hhu   "           , ++idx, VALEPS(t->coef[i],eps_), d->ords[i]);
#else
      if (fabs(creal(t->coef[i])) < eps_ && fabs(cimag(t->coef[i])) < eps_) continue;
      fprintf(stream_, "\n%6d  %21.14lE %+21.14lEi   %2hhu   ", ++idx, VALEPS(t->coef[i],eps_), d->ords[i]);
#endif
      (d->nn > 20 ? print_ords_sm : print_ords)(d->nn, d->To[i], stream_);
    }
  }

  if (!idx) fprintf(stream_, "\n          ALL COMPONENTS ZERO");

  fprintf(stream_, "\n\n");

  DBGTPSA(t); DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
