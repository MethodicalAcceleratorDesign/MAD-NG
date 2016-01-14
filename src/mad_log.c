/*
 o----------------------------------------------------------------------------o
 |
 | Logging module implementation
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
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>

#include "mad_main.h"
#include "mad_log.h"

#undef fatal
#undef error
#undef warn
#undef info
#undef debug
#undef ensure

/*
Note about thread-safety:
  - Globals  are shared ints, lockless atomic access is safe
  - Counters are shared ints, lockless atomic access is safe
  - Locations must be thread specific for consistent output
*/

// --- globals ---------------------------------------------------------------o

int mad_info_level     = 0;
int mad_debug_level    = 0;
int mad_trace_location = 0;

// --- locals ----------------------------------------------------------------o

static long error_count = 0; // only incremented
static long warn_count  = 0; // only incremented
static long info_count  = 0; // only incremented
static long debug_count = 0; // only incremented

static struct location {
  str_t file;
  int   line;
  int   extrn;
} location[1] = {{ "", 0, 0 }};

#ifdef _OPENMP
#pragma omp threadprivate(location)
#endif

// --- implementation --------------------------------------------------------o

void
(mad_fatal) (str_t msg) {
  mad_fatalf(msg);
}

void
(mad_error) (str_t msg) {
  mad_errorf(msg);
}

void
(mad_warn) (str_t msg) {
  mad_warnf(msg);
}

void
(mad_info) (int lvl, str_t msg) {
  mad_infof(lvl, msg);
}

void
(mad_debug) (int lvl, str_t msg) {
  mad_debugf(lvl, msg);
}

void
mad_log_setloc (str_t file, int line)
{
  struct location *loc = location;
  loc->file = file;
  loc->line = line;
  loc->extrn = 0;
}

void
mad_log_setloc1 (str_t file, int line)
{
  struct location *loc = location;
  loc->file = file;
  loc->line = line;
  loc->extrn = 1;
}

void
mad_log_getstat (long *error_, long *warn_, long *info_, long *debug_)
{
  if (error_) *error_  = error_count;
  if (warn_ ) *warn_   =  warn_count;
  if (info_ ) *info_   =  info_count;
  if (debug_) *debug_  = debug_count;
}

void
mad_fatalf (str_t fmt, ...)
{
  assert(fmt);
  fflush(stdout);

  error_count += 1;

  struct location *loc = location;
  fprintf(stderr, "ERROR:%s:%d: ", loc->file, loc->line);

  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);

  putc('\n', stderr);
  fflush(stderr);

  // if (location->extrn) mad_fatal_extrn();

  exit(EXIT_FAILURE);
}

void
mad_errorf (str_t fmt, ...)
{
  assert(fmt);
  fflush(stdout);

  error_count += 1;

  if (mad_trace_location >= 2) {
    struct location *loc = location;
    fprintf(stderr, "ERROR:%s:%d: ", loc->file, loc->line);
  } else
    fputs("ERROR: ", stderr);

  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);

  putc('\n', stderr);
  fflush(stderr);
}

void
mad_warnf (str_t fmt, ...)
{
  assert(fmt);
  fflush(stdout);

  warn_count += 1;

  if (mad_trace_location >= 2) {
    struct location *loc = location;
    fprintf(stderr, "WARNING:%s:%d: ", loc->file, loc->line);
  } else
    fputs("WARNING: ", stderr);

  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);

  putc('\n', stderr);
  fflush(stderr);
}

void
mad_infof (int level, str_t fmt, ...)
{
  if (mad_info_level < level) return;

  assert(fmt);
  fflush(stderr);

  info_count += 1;

  if (mad_trace_location >= 3) {
    struct location *loc = location;
    fprintf(stdout, "INFO%d:%s:%d: ", mad_info_level, loc->file, loc->line);
  }

  va_list ap;
  va_start(ap, fmt);
  vfprintf(stdout, fmt, ap);
  va_end(ap);

  putc('\n', stdout);
  fflush(stdout);
}

void
mad_debugf (int level, str_t fmt, ...)
{
  if (mad_debug_level < level) return;

  assert(fmt);
  fflush(stdout);

  debug_count += 1;

  if (mad_trace_location >= 1) {
    struct location *loc = location;
    fprintf(stderr, "DEBUG%d:%s:%d: ", mad_debug_level, loc->file, loc->line);
  } else
    fprintf(stderr, "DEBUG%d: ", mad_debug_level);

  va_list ap;
  va_start(ap, fmt);
  vfprintf(stderr, fmt, ap);
  va_end(ap);

  putc('\n', stderr);
  fflush(stderr);
}
