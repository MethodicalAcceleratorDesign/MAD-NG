#ifndef MAD_MAIN_H
#define MAD_MAIN_H

/*
 o-----------------------------------------------------------------------------o
 |
 | MAD frontend
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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

  Purpose:
  - Frontend (main) of the MAD application.

  Comment:
  - MAD is embedding the LuaJIT library and frontend developped by Mike Pall
    modified for the purpose of MAD. See COPYRIGHT.luajit in the lib/patches
    directory.
 o-----------------------------------------------------------------------------o
 */

// --- interface --------------------------------------------------------------o

#ifndef LUALIB_API // Only valid in a Lua environment
#define LUALIB_API
#endif

LUALIB_API void mad_error (     const char *fn, const char *fmt, ...);
LUALIB_API void mad_warn  (     const char *fn, const char *fmt, ...);
LUALIB_API void mad_trace (int, const char *fn, const char *fmt, ...);

// --- globals ----------------------------------------------------------------o

extern int mad_trace_level;
extern int mad_trace_location;

// ----------------------------------------------------------------------------o
// --- implementation (private) -----------------------------------------------o
// ----------------------------------------------------------------------------o

LUALIB_API void mad_error(const char*,const char*,...) __attribute((noreturn));

// ----------------------------------------------------------------------------o

#endif // MAD_MAIN_H