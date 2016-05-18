#ifndef MAD_MAIN_H
#define MAD_MAIN_H

/*
 o----------------------------------------------------------------------------o
 |
 | MAD frontend
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
  - Frontend (main) of the MAD application.

  Comment:
    MAD is embedding the LuaJIT library and frontend developped by Mike Pall.
    See COPYRIGHT.luajit in this directory.
 o----------------------------------------------------------------------------o
 */

// --- interface -------------------------------------------------------------o

void mad_luaerr (const char *msg);
void mad_luawrn (const char *msg);
void mad_luatrc (const char *msg);

#ifdef LUALIB_API // Only valid in a Lua environment
LUALIB_API int luaL_error (lua_State *L, const char *fmt, ...);
LUALIB_API int luaL_warn  (lua_State *L, const char *fmt, ...);
LUALIB_API int luaL_trace (lua_State *L, const char *fmt, ...);
#endif

// ---------------------------------------------------------------------------o

#endif // MAD_MAIN_H