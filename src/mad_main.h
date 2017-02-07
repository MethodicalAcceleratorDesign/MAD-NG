#ifndef MAD_MAIN_H
#define MAD_MAIN_H

/*
 o----------------------------------------------------------------------------o
 |
 | MAD frontend
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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
  - Frontend (main) of the MAD application.

  Comment:
  - MAD is embedding the LuaJIT library and frontend developped by Mike Pall
    modified for the purpose of MAD. See COPYRIGHT.luajit in the lib/patches
    directory.
 o----------------------------------------------------------------------------o
 */

// --- interface -------------------------------------------------------------o

void mad_error (const char *fname, const char *fmt, ...);
void mad_warn  (const char *fname, const char *fmt, ...);
void mad_trace (const char *fname, const char *fmt, ...);

#ifdef LUALIB_API // Only valid in a Lua environment
LUALIB_API int luaL_error (lua_State *L, const char *fmt, ...);
LUALIB_API int luaL_warn  (lua_State *L, const char *fmt, ...);
LUALIB_API int luaL_trace (lua_State *L, const char *fmt, ...);
#endif

// --- globals ---------------------------------------------------------------o

extern int mad_trace_level;
extern int mad_trace_location;

// ---------------------------------------------------------------------------o

#endif // MAD_MAIN_H