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

#include <ctype.h>
#include <string.h>

#include "mad_mem.h"
#include "mad_str.h"
#ifdef    MAD_CTPSA_IMPL
#include "mad_ctpsa_impl.h"
#else
#include "mad_tpsa_impl.h"
#endif

//#undef  DEBUG
//#define DEBUG 3

// --- local ------------------------------------------------------------------o

static inline int
skip_line_dbg(FILE *stream)
{
  int c;
  printf("skip line\n");
  while ((c = fgetc(stream)) != '\n' && c != EOF) putchar(c);
  putchar(c);
  return c;
}

static inline int
skip_line(FILE *stream)
{
#if DEBUG > 2
  return skip_line_dbg(stream);
#else
  int c;
  while ((c = fgetc(stream)) != '\n' && c != EOF) ;
  return c;
#endif
}

static inline int
skip_spaces(FILE *stream)
{
  int c;
  while ((c = fgetc(stream)) == ' ' || c == '\t') ;
  ungetc(c, stream);
  return c;
}

static inline int
skip_wspaces(FILE *stream)
{
  int c;
  while ((c=getc(stream)) != EOF && isspace(c)) ;
  ungetc(c, stream);
  return c;
}

// static inline void
// print_ords_sm(int n, const ord_t ords[n], FILE *stream)
// {
//   assert(ords && stream);
//   for (int i=0; i < n; i++)
//     if (ords[i]) fprintf(stream, "  %d^%hhu", i+1, ords[i]);
// }

static inline void
print_ords(int nv, int np, ord_t po, const ord_t ords[nv+np], FILE *stream)
{
  assert(ords && stream);

  // print variables
  FOR(i,0,nv-1,2) fprintf(stream, "  %hhu %hhu", ords[i], ords[i+1]);
  if (nv % 2)     fprintf(stream, "  %hhu"     , ords[nv-1]);

  // print parameters
  FOR(i,nv,nv+np) if (ords[i] != po) fprintf(stream, "  %d^%hhu", i+1, ords[i]);
}

static inline void
read_ords(int nv, int np, ord_t po, ord_t ords[nv+np], FILE *stream, int ci, str_t name)
{
  assert(ords && stream);

  if (!name[0]) name = "-UNNAMED-";

  mad_mono_fill(nv, ords   ,  0);
  mad_mono_fill(np, ords+nv, po);

  // read variables
  FOR(i,0,nv-1,2)
    if (fscanf(stream, "%*[ ]%hhu%*[ ]%hhu", &ords[i], &ords[i+1]) != 2)
      error("invalid monomial input at index %d of '%s'", ci, name);
  if (nv % 2)
    if (fscanf(stream, "%*[ ]%hhu"         , &ords[nv-1]         ) != 1)
      error("invalid monomial input at index %d of '%s'", ci, name);

  // read parameters
  idx_t idx;
  ord_t ord;
  FOR(i,nv,nv+np) {
    idx = 0, ord = -1;
    int cnt = fscanf(stream, "%*[ ]%d^%hhu", &idx, &ord);

#if DEBUG > 2
    int chr = getc(stream); ungetc(chr, stream);
    printf("mono: ci=%d, oi=%d[%d], cnt=%d, idx=%d, ord=%d, nxtchr='%c'\n",
           ci, i, nv+np, cnt, idx, ord, chr);
#endif

         if (cnt == 0) break;
    else if (cnt == 2) ;
    else error("invalid monomial input at index %d of '%s'", ci, name);

    ensure(nv < idx && idx <= nv+np,
           "invalid parameter index (expecting %d < %d <= %d) at index %d of '%s'",
            nv, idx, nv+np, ci, name);
    ensure(0 < ord && ord <= DESC_MAX_ORD,
           "invalid order (expecting 0 < %d <= %d) at index %d of '%s'",
            ord, DESC_MAX_ORD, ci, name);

    ords[idx-1] = ord;
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
  skip_wspaces(stream_);

  // read the name (which is 15+'\0' chars)
  char name[NAMSZ]= {0}, sep='?';
  int cnt = fscanf(stream_, "%15[^:,\t\n]%c", name, &sep);

#if DEBUG > 2
    printf("header: cnt=%d, name='%s', sep='%c'\n", cnt, name, sep);
#endif

  if (!(cnt == 2 && (sep == ':' || sep == ','))) {
    warn("unable to parse TPSA header: '%s'", name);
    fsetpos(stream_, &fpos); // may fail for non-seekable stream (e.g. pipes)...
    return NULL;
  }

  ensure(!feof(stream_) && !ferror(stream_), "invalid input (file error?)");

  char knd=0;
  int nv=0, np=0, nc;
  ord_t mo=0, po=0;
  log_t ptc=sep == ',';

  if (ptc) // n = 2, swap input for PTC/BERZ
    cnt = fscanf(stream_, "%*[ ]NO%*[ ]=%hhu,%*[ ]NV%*[ ]=%d",
                                        &mo,              &nv);
  else     // n = 3 or 5
    cnt = fscanf(stream_,"%*[ ]%c,%*[ ]NV%*[ ]=%d,%*[ ]MO%*[ ]=%hhu"
                                ",%*[ ]NP%*[ ]=%d,%*[ ]PO%*[ ]=%hhu",
                               &knd,           &nv,            &mo,
                                               &np,            &po);
#if DEBUG > 2
    printf("header: cnt=%d, knd='%c', nv=%d, mo=%d, np=%d, po=%d\n",
                       cnt, knd?knd:'?', nv,    mo,    np,    po);
#endif

  // sanity checks
  ensure(0 <  nv && nv <= DESC_MAX_VAR, "invalid NV=%d", nv);
  ensure(           mo <= DESC_MAX_ORD, "invalid MO=%d", mo);
  ensure(strchr("RC ", knd), "invalid kind='%c' (expecting R or C)", knd);

  if (kind_) {
    ensure(*kind_ >= -1 && *kind_ <= 1, "invalid kind (expecting -1, 0, 1)");
    if (*kind_ == -1) *kind_ = knd == 'C';
    else if (knd && knd != "RC"[*kind_])
      warn("kind specification '%c' differs from input '%c'", "RC"[*kind_], knd);

#if DEBUG > 2
    printf("header: kind='%c'\n", "RC"[*kind_]);
#endif
  }

  if (name_) {
    int bnd[2] = {0, strlen(name)};
    mad_str_trim(name, bnd);
    memcpy(name_, name+bnd[0], bnd[1]);
    name_[bnd[1]] = '\0';

#if DEBUG > 2
    printf("header: name='%s'\n", name_);
#endif
  }

  if (cnt == 2 || cnt == 3) {
    // TPSA -- ignore rest of lines
    ensure(skip_line(stream_) != EOF, "invalid input (file error?)"); // finish NV,MO line
    ensure(skip_line(stream_) != EOF, "invalid input (file error?)"); // discard *****

    const D* ret = mad_desc_newv(nv, mo);
    DBGFUN(<-);
    return ret;
  }

  if (cnt == 5) {
    int nn = nv+np;

    // sanity checks
    ensure(0 <= np && nn <= DESC_MAX_VAR, "invalid NP=%d", np);
    ensure(           po <= DESC_MAX_ORD, "invalid PO=%d", po);

    // read variables orders if present (GTPSA)
    ord_t no[nn];
    (void)fscanf(stream_, ",%*[ ]NO%*[ ]=%n", &nc);
    if (nc >= 6) {
      read_ords(nv,np,po,no,stream_,-1,name); ++cnt;
    }
    ensure(skip_line(stream_) != EOF, "invalid input (file error?)"); // finish NO line
    ensure(skip_line(stream_) != EOF, "invalid input (file error?)"); // discard *****

    const D* ret = cnt == 5 ? mad_desc_newvp (nv, mo, np, po)
                            : mad_desc_newvpo(nv, mo, np, po, no);
    DBGFUN(<-);
    return ret;
  }

       if (cnt < 2) warn("could not read (NV,%s) from header", ptc?"NO":"MO");
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

  NUM   v = 0;
  int   nn = t->d->nn, nv = t->d->nv, np = t->d->np;
  int   cnt = -1, i = -1, nc;
  ord_t o = 0, ords[nn];
  FUN(reset0)(t);

  // process coef header or summary
  int c = skip_wspaces(stream_);

  #if DEBUG > 2
    printf("coef: c='%c'\n", c);
  #endif

  if (c == 'I') {
    (void)fscanf(stream_, "I%*[ ]COEFFICIENT%*[ ]ORDER%*[ ]EXPONENTS%n", &nc);
    if (nc < 29) warn("unable to parse GTPSA coefficients for '%s'",
                       t->nam[0] ? t->nam : "-UNNAMED-");
    #if DEBUG > 2
      printf("coef: header 'I COEF...' parsed\n");
    #endif
    ensure((c=skip_wspaces(stream_)) != EOF, "invalid input (file error?)");
  }

  if (c == 'A') {
    // works for 0_dp, ZERO and EPS
    (void)fscanf(stream_, "ALL%*[ ]COMPONENTS%n", &nc);
    if (nc < 14) warn("unable to parse GTPSA coefficients for '%s'",
                      t->nam[0] ? t->nam : "-UNNAMED-");
    #if DEBUG > 2
      printf("coef: 'ALL COMP...' parsed\n");
    #endif
    ensure(skip_line(stream_) != EOF, "invalid input (file error?)");
    DBGTPSA(t); DBGFUN(<-); return;
  }

  for (;;) {
    skip_spaces(stream_);

    // read index (avoid %d in case next TPSA name is a number!)
    char idx[16] = {0};
    cnt = fscanf(stream_, "%15[0-9]", idx);
    if (cnt != 1) break;
    i = strtol(idx, 0, 0);

    // read coef and order
#ifndef MAD_CTPSA_IMPL
    cnt = fscanf(stream_, "%lG%*[ ]%hhu", &v, &o);
    if (cnt != 2) break;
#else
    char chr;
    cnt = fscanf(stream_, "%lG%lG%c%*[ ]%hhu", (num_t*)&v, (num_t*)&v+1, &chr, &o);
    if (cnt != 4) break;
    ensure(chr == ' ' || chr == 'i',
           "invalid complex number format (' ' or 'i' expected ending)"
           " at index %d of '%s'", i, t->nam);
#endif

    #if DEBUG > 2
      printf("coef: i=%d, v=" FMT ", o=%d\n", i, VAL(v), o);
    #endif

    read_ords(nv,np,0,ords,stream_,i,t->nam); // sanity check
    ensure(mad_mono_ord(nn,ords) == o,
           "invalid monomial order at index %d of '%s'", i, t->nam);

    // discard too high mononial and zeros
    if (o <= t->mo && v) FUN(setm)(t,nn,ords,0,v);

    // finish line (handle no '\n' before EOF)
    skip_line(stream_);
  }

  #if DEBUG > 2
    printf("endc: i=%d, cnt=%d, v=" FMT ", o=%d\n\n", i, cnt, VAL(v), o);
  #endif

  if (i == -1) // no coef read
    warn("unable to parse GTPSA coefficients for '%s'",
         t->nam[0] ? t->nam : "-UNNAMED-");

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
    t = FUN(newd)(d, mad_tpsa_dflt);
    FUN(scan_coef)(t, stream_);
    FUN(nam)      (t, name   );
  }
  DBGFUN(<-);
  return t;
}

void
FUN(print) (const T *t, str_t name_, num_t eps_, int nohdr_, FILE *stream_)
{
  assert(t); DBGFUN(->); DBGTPSA(t);

  if (!name_  ) name_   = t->nam[0] ? t->nam : "-UNNAMED-";
  if (eps_ < 0) eps_    = 0;
  if (!stream_) stream_ = stdout;

#ifndef MAD_CTPSA_IMPL
  const char typ = 'R';
#else
  const char typ = 'C';
#endif
  const D *d = t->d;
  if (nohdr_) goto onlycoef;

  // print header
  fprintf(stream_, d->np || d->uno
                 ? "\n %-8s:  %c, NV = %3d, MO = %2hhu, NP = %3d, PO = %2hhu"
                 : "\n %-8s:  %c, NV = %3d, MO = %2hhu",
                      name_, typ,    d->nv,      d->mo,    d->np,      d->po);

  if (d->uno) {
    fprintf(stream_, ", NO = ");
    print_ords(d->nv, d->np, d->po, d->no, stream_);
  }
  fprintf(stream_, "\n *******************************************************");
#ifdef MAD_CTPSA_IMPL
  fprintf(stream_, "***********************");
#endif

onlycoef: ;
  idx_t idx = 0;
  bit_t tnz = mad_bit_set(t->nz,0);
  TPSA_SCAN(t,0,t->hi,tnz) {
#ifndef MAD_CTPSA_IMPL
    if (fabs(t->coef[i]) < eps_) continue;
    if (!idx) fprintf(stream_, "\n     I   COEFFICIENT             ORDER   EXPONENTS");
    fprintf(stream_, "\n%6d  %23.16lE   %2hhu   "           , ++idx, VALEPS(t->coef[i],eps_), d->ords[i]);
#else
    if (fabs(creal(t->coef[i])) < eps_ && fabs(cimag(t->coef[i])) < eps_) continue;
    if (!idx) fprintf(stream_, "\n     I   COEFFICIENT                                      ORDER   EXPONENTS");
    fprintf(stream_, "\n%6d  %23.16lE %+23.16lEi   %2hhu   ", ++idx, VALEPS(t->coef[i],eps_), d->ords[i]);
#endif
    print_ords(d->nv, d->np, 0, d->To[i], stream_);
  }
  if (!idx) fprintf(stream_, "\n\n         ALL COMPONENTS ZERO (EPS=%.1lE)", eps_);
  fprintf(stream_, "\n");
  DBGFUN(<-);
}

// --- end --------------------------------------------------------------------o
