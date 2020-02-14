/*
 o-----------------------------------------------------------------------------o
 |
 | Logging module implementation
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  WARNING: this file is not part of MAD, it must be used only for standalone
           library build.

*/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "mad_log.h"

int mad_warn_count     = 0;
int mad_trace_level    = 0;
int mad_trace_location = 0;

LUALIB_API void (mad_error) (str_t fn, str_t fmt, ...)
{
  va_list va;
  va_start(va, fmt);
   fprintf(stderr, fn ? "error: %s: " : "error: ", fn);
  vfprintf(stderr, fmt, va);
  va_end(va);
  fputc('\n', stderr);
  exit(EXIT_FAILURE); /* never reached */
}

LUALIB_API void (mad_warn) (str_t fn, str_t fmt, ...)
{
  ++mad_warn_count;
  va_list va;
  va_start(va, fmt);
   fprintf(stderr, fn ? "warning: %s: " : "warning: ", fn);
  vfprintf(stderr, fmt, va);
  va_end(va);
  fputc('\n', stderr);
}

LUALIB_API void (mad_trace) (int lvl, str_t fn, str_t fmt, ...)
{
  if (mad_trace_level < lvl) return;
  va_list va;
  va_start(va, fmt);
  if (fn) fprintf(stderr, "%s", fn);
  vfprintf(stderr, fmt, va);
  va_end(va);
  fputc('\n', stderr);
}
