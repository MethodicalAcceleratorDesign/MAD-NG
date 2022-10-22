/*
 o-----------------------------------------------------------------------------o
 |
 | MAD frontend
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

  Purpose:
  - Frontend (main) of the MAD application.

  Comment:
  - MAD is embedding the LuaJIT library and frontend developped by Mike Pall
    modified for the purpose of MAD. See COPYRIGHT.luajit in the lib/patches
    directory.
 o-----------------------------------------------------------------------------o
 */

/*
** LuaJIT frontend. Runs commands, scripts, read-eval-print (REPL) etc.
** Copyright (C) 2005-2017 Mike Pall. See Copyright Notice in luajit.h
**
** Major portions taken verbatim or adapted from the Lua interpreter.
** Copyright (C) 1994-2008 Lua.org, PUC-Rio. See Copyright Notice in lua.h
*/

#define _GNU_SOURCE 1
#define _POSIX_C_SOURCE 1
#define _XOPEN_UNIX 1
#define _XOPEN_VERSION 700
#define _DARWIN_BETTER_REALPATH 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// to shut up gcc 8 invalid warnings!
int snprintf(char *restrict str, size_t size, const char* restrict format, ...);

#define luajit_c

#include "lua.h"
#include "lauxlib.h"
#include "lualib.h"
#include "luajit.h"

#include "lj_arch.h"

#if LJ_TARGET_POSIX
#include <unistd.h>
#define lua_stdin_is_tty()  isatty(0)
#elif LJ_TARGET_WINDOWS
#include <io.h>
#ifdef __BORLANDC__
#define lua_stdin_is_tty()  isatty(_fileno(stdin))
#else
#define lua_stdin_is_tty()  _isatty(_fileno(stdin))
#endif
#else
#define lua_stdin_is_tty()  1
#endif

#if !LJ_TARGET_CONSOLE
#include <signal.h>
#endif

static lua_State *globalL = NULL;
static const char *progname = "mad";

/* --- MAD (start) -----------------------------------------------------------*/

/* Assume Posix: MacOSX, Linux, Mingw32/64 or Cygwin */
#include <unistd.h>
#include <limits.h>
#include <assert.h>
#include <time.h>
#include "lj_def.h"
#include "mad_ver.h"
#include "mad_log.h"

/* globals */
int mad_warn_count     = 0;
int mad_stdin_is_tty   = 0;
int mad_is_interactive = 0;

int mad_trace_level    = 0;
int mad_trace_location = 0;
int mad_trace_fortid   = 0;

/* Forward declarations. */
static int dofile   (lua_State *L, str_t name);
static int dostring (lua_State *L, str_t s, str_t name);
static int report   (lua_State *L, int status);

#if defined(__MINGW32__) || defined(__MINGW64__)
#include <sys/stat.h>

/* Missing declaration in Mingw */
int setenv(str_t name, str_t value, int overwrite);

/* Better tty detection under Mingw */
#undef     lua_stdin_is_tty
static int lua_stdin_is_tty (void)
{
  struct stat stats;
  fstat(0, &stats);
  return S_ISFIFO(stats.st_mode) || isatty(0);
}
#endif /* __MINGW32/64__ */

static void print_mad_version (void)
{
  str_t rel = MKSTR(MAD_VERSION)
  #ifdef _OPENMP
  "_P"
  #endif
   " (" LJ_OS_NAME " " MKSTR(LJ_ARCH_BITS) ")";
  str_t msg =
  "    ____  __   ______    ______     |   Methodical Accelerator Design\n"
  "     /  \\/  \\   /  _  \\   /  _  \\   |   release: %s\n"
  "    /  __   /  /  /_/ /  /  /_/ /   |   support: http://cern.ch/mad\n"
  "   /__/  /_/  /__/ /_/  /_____ /    |   licence: GPL3 (C) CERN 2016+\n"
  "                                    |   started: %s\n"
  "\n";

  time_t t = time(NULL);
  struct tm *tm = localtime(&t);
  char buf[80];

  strftime(buf, sizeof buf, "%Y-%m-%d %H:%M:%S", tm);
  printf(msg, rel, buf);
}

LUALIB_API void (mad_error) (str_t fn, str_t fmt, ...)
{
  va_list va;
  va_start(va, fmt);
  fflush(stdout);
  fprintf(stderr, fn ? "error: %s: " : "error: ", fn);
  vfprintf(stderr, fmt, va);
  va_end(va);
  fputc('\n', stderr);
  lua_pushstring(globalL, "");  // NOT SAFE! (but let's try...)
  lua_error(globalL);
  exit(EXIT_FAILURE); /* never reached */
}

LUALIB_API void (mad_warn) (str_t fn, str_t fmt, ...)
{
  ++mad_warn_count;
  va_list va;
  va_start(va, fmt);
  fflush(stdout);
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
  fflush(stdout);
  if (fn) fprintf(stderr, "%s", fn);
  vfprintf(stderr, fmt, va);
  va_end(va);
  fputc('\n', stderr);
}

static int mad_luawarn (lua_State *L)
{
  (mad_warn)(NULL, "%s", luaL_checkstring(L,1));
  return 0;
}

static int mad_luatrace (lua_State *L)
{
  int level = luaL_checknumber(L,1);
  str_t where = NULL;

  if (mad_trace_location && mad_trace_level >= level) {
    luaL_where(L, 2);
    where = lua_tostring(L, -1);
  }
  (mad_trace)(level, where, "%s", luaL_checkstring(L,2));
  return 0;
}

/* Handle paths */

/* Windows: not declared by any mean but provided by libgettextlib */
extern char *realpath (str_t restrict fname, char *restrict rname);

#if LUAJIT_OS == LUAJIT_OS_WINDOWS
static char* winpath (char *buf)
{
  char *p = buf;
  while ((p = strchr(p, '/'))) *p++ = '\\';
  return buf;
}

static char* winexe (char *buf)
{
  int len = strlen(buf);
  if (len > 0 && (len < 4 || strcmp(buf+len-4, ".exe")))
    strcpy(buf+len, ".exe");
  return buf;
}
#endif

static void mad_setenv (lua_State *L, int no_env)
{
  char prog_name[PATH_MAX+1] = "";
  char prog_path[PATH_MAX+1] = "";
  char curr_path[PATH_MAX+1] = "";
  char home_path[PATH_MAX+1] = "";
  char nam[PATH_MAX+4+1], buf[PATH_MAX+4+1];
  char *path, *p, nul = 0, psep = ':', dsep = '/';

  /* get curr_path [=getenv(unix:"PWD" or win:"CD")+cannonize] */
  if (!(path = getcwd(curr_path, sizeof curr_path))) *curr_path = nul;
  if (!path && !(path = getenv("PATH"))) path = &nul;

  /* retrieve separators */
  p = strchr(path, '/'); if (!p) p = strchr(path, '\\');
  if (p) dsep = *p, psep = *p == '\\' ? ';' : ':';

  /* get home_path */
#if LUAJIT_OS == LUAJIT_OS_WINDOWS
  if ((path = getenv("HOMEDRIVE"))) strcpy(home_path, path);
  if ((path = getenv("HOMEPATH" ))) strcat(home_path, path);
  winpath(home_path);
#else
  if ((path = getenv("HOME"     ))) strcpy(home_path, path);
#endif

  /* get prog_name */
  strcpy(nam, progname);
#if LUAJIT_OS == LUAJIT_OS_WINDOWS
  winexe(winpath(prog_path));  /* canonize */
#endif
  p = strrchr(nam, dsep);
  strcpy(prog_name, p ? p+1 : nam);

  /* get prog_path */
  if (realpath(nam, prog_path)) goto found; /* absolute path */
  if (snprintf(buf, sizeof buf, "%s%c%s", curr_path, dsep, nam) > 0 &&
      realpath(buf, prog_path)) goto found; /* relative path */
  if ((path = getenv("PATH")) && *path) {   /* $PATH path */
    for(;;) {
      p = strchr(path, psep);
      if (p) *p = nul;
      if (snprintf(buf, sizeof buf, "%s%c%s", path, dsep, prog_name) > 0 &&
          realpath(buf, prog_path)) { if (p) *p = psep; goto found; }
      if (p) *p = psep, path = p+1; else break;
    }
  }
  *prog_path = nul;

found:
#if LUAJIT_OS == LUAJIT_OS_WINDOWS
  winexe(winpath(prog_path));  /* canonize */
#endif
  if ((p=strrchr(prog_path, dsep))) {
    *p = nul;
    if (p[1]) strcpy(prog_name, p+1);
  }

  /* add trailing dir separator */
  if (*curr_path && (p=curr_path+strlen(curr_path))[-1] != dsep)
    p[0] = dsep, p[1] = nul;
  if (*home_path && (p=home_path+strlen(home_path))[-1] != dsep)
    p[0] = dsep, p[1] = nul;
  if (*prog_path && (p=prog_path+strlen(prog_path))[-1] != dsep)
    p[0] = dsep, p[1] = nul;

  buf[0] = dsep, buf[1] = nul;
  str_t dir_sep = buf;

  /* set global table '_M' */
  enum { nitem = 9 };
  str_t list[2][nitem] = {
    { "currpath", "homepath", "progpath", "progname", "dirsep",
      "version", "version_date", "os", "arch" },
    { curr_path , home_path , prog_path , prog_name , dir_sep ,
      MKSTR(MAD_VERSION), MKSTR(MAD_VERSION_DATE),
      LJ_OS_NAME, MKSTR(LJ_ARCH_BITS) }
  };
  lua_createtable(L, 0, nitem);
  for (int i = 0; i < nitem; i++) {
    assert(list[0][i]); lua_pushstring(L, list[0][i]);
    assert(list[1][i]); lua_pushstring(L, list[1][i]);
    lua_rawset(L, -3);
  }
  lua_setglobal(L, "_M");

  if (no_env) return;

  str_t mpath = getenv("MAD_PATH"); str_t cpath  = getenv("MAD_CPATH");
  str_t lpath = getenv("LUA_PATH"); str_t lcpath = getenv("LUA_CPATH");
  int len, mlen, clen, marg, carg;

  if (mpath) marg = 0, mlen = strlen(mpath);
  else
    mpath = "./?.mad;./?.lua;"                        // CURR_PATH
            /* MAD modules (user's home) */
            "%s.mad/?.mad;"                           // HOME_PATH  (1)
            "%s.mad/?.lua;"                           // HOME_PATH  (1)
            /* MAD modules (relative) */
            "%s?.mad;"                                // MAD_PATH   (1)
            "%s?.lua;"                                // MAD_PATH   (1)
#if LUAJIT_OS == LUAJIT_OS_WINDOWS
            /* MAD modules (relative) */
            "%s../mad/?.mad;"                         // MAD_PATH   (1)
            "%s../mad/?.lua;"                         // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (relative) */
            "%s../Lua/5.1/?.lua;"                     // MAD_PATH   (1)
            "%s../Lua/5.1/?/init.lua;"                // MAD_PATH   (1)
            "%s../Lua/5.1/lua/?.lua;"                 // MAD_PATH   (1)
            "%s../Lua/5.1/lua/?/init.lua;"            // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (absolute) */
            "C:/Program Files/Lua/5.1/?.lua;"
            "C:/Program Files/Lua/5.1/?/init.lua;"
            "C:/Program Files/Lua/5.1/lua/?.lua;"
            "C:/Program Files/Lua/5.1/lua/?/init.lua;",
    marg = 4+6,
#else
            /* MAD modules (relative) */
            "%s../share/mad/?.mad;"                   // MAD_PATH   (1)
            "%s../share/mad/?.lua;"                   // MAD_PATH   (1)
            "%s../lib/mad/?.mad;"                     // MAD_PATH   (1)
            "%s../lib/mad/?.lua;"                     // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (relative) */
            "%s../share/lua/5.1/?.lua;"               // MAD_PATH   (1)
            "%s../share/lua/5.1/?/init.lua;"          // MAD_PATH   (1)
            "%s../lib/lua/5.1/?.lua;"                 // MAD_PATH   (1)
            "%s../lib/lua/5.1/?/init.lua;"            // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (/opt/local absolute) */
            "/opt/local/share/lua/5.1/?.lua;"
            "/opt/local/share/lua/5.1/?/init.lua;"
            "/opt/local/lib/lua/5.1/?.lua;"
            "/opt/local/lib/lua/5.1/?/init.lua;"
            /* From Lua unofficial uFAQ (/usr/local absolute) */
            "/usr/local/share/lua/5.1/?.lua;"
            "/usr/local/share/lua/5.1/?/init.lua;"
            "/usr/local/lib/lua/5.1/?.lua;"
            "/usr/local/lib/lua/5.1/?/init.lua;",
    marg = 4+8,
#endif
    mlen = strlen(mpath)-(2*marg)
           + 2*strlen(home_path) + (marg-2)*strlen(prog_path);

  if (cpath) carg = 0, clen = strlen(cpath);
  else
    cpath =
#if LUAJIT_OS == LUAJIT_OS_WINDOWS
            "./?.dll;./?.so;"                         // CURR_PATH
            /* MAD modules (user's home) */
            "%s.mad/?.dll;"                           // HOME_PATH  (1)
            "%s.mad/?.so;"                            // HOME_PATH  (1)
            /* MAD modules (relative) */
            "%s?.dll;"                                // MAD_PATH   (1)
            "%s?.so;"                                 // MAD_PATH   (1)
            "%s../lib/mad/?.dll;"                     // MAD_PATH   (1)
            "%s../lib/mad/?.so;"                      // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (relative) */
            "%s../lib/lua/5.1/?.dll;"                 // MAD_PATH   (1)
            "%s../lib/lua/5.1/?.so;"                  // MAD_PATH   (1)
            "C:/Program Files/Lua/5.1/?.dll;"
            "C:/Program Files/Lua/5.1/?.so;"
            "C:/Program Files/Lua/5.1/clibs/?.dll;"
            "C:/Program Files/Lua/5.1/clibs/?.so;",
    carg = 8,
#elif LUAJIT_OS == LUAJIT_OS_OSX
            "./?.dylib;./?.so;"                       // CURR_PATH
            /* MAD modules (user's home) */
            "%s.mad/?.dylib;"                         // HOME_PATH  (1)
            "%s.mad/?.so;"                            // HOME_PATH  (1)
            /* MAD modules (relative) */
            "%s?.dylib;"                              // MAD_PATH   (1)
            "%s?.so;"                                 // MAD_PATH   (1)
            "%s../lib/mad/?.dylib;"                   // MAD_PATH   (1)
            "%s../lib/mad/?.so;"                      // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (relative) */
            "%s../lib/lua/5.1/?.dylib;"               // MAD_PATH   (1)
            "%s../lib/lua/5.1/?.so;"                  // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (/opt/local absolute) */
            "/opt/local/lib/lua/5.1/?.dylib;"
            "/opt/local/lib/lua/5.1/?.so;"
            /* From Lua unofficial uFAQ (/usr/local absolute) */
            "/usr/local/lib/lua/5.1/?.dylib;"
            "/usr/local/lib/lua/5.1/?.so;",
    carg = 8,
#else
            "./?.so;"                                 // CURR_PATH
            /* MAD modules (user's home) */
            "%s.mad/?.so;"                            // HOME_PATH  (1)
            /* MAD modules (relative) */
            "%s?.so;"                                 // MAD_PATH   (1)
            "%s../lib/mad/?.so;"                      // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (relative) */
            "%s../lib/lua/5.1/?.so;"                  // MAD_PATH   (1)
            /* From Lua unofficial uFAQ (/opt/local absolute) */
            "/opt/local/lib/lua/5.1/?.so;"
            /* From Lua unofficial uFAQ (/usr/local absolute) */
            "/usr/local/lib/lua/5.1/?.so;",
    carg = 4,
#endif
    clen = strlen(cpath)-(2*carg)
           + (carg/4)*strlen(home_path) + (carg-carg/4)*strlen(prog_path);

  /* LUA_PATH = $MAD_PATH;$LUA_PATH */
  char env[(mlen>clen?mlen:clen) + strlen(lpath ? lpath : "") + 1];
  ensure(marg == 10 || marg == 12, "invalid number of path argument");
  len = snprintf(env, sizeof env, mpath, /* next args discarded without %s */
                  home_path, home_path, prog_path, prog_path, prog_path,
                  prog_path, prog_path, prog_path, prog_path, prog_path,
                  prog_path, prog_path);
  if (lpath) strcat(env, lpath); else if (len>0) env[len-1] = nul;
#if LUAJIT_OS == LUAJIT_OS_WINDOWS
  winpath(env); /* canonize */
#endif
  setenv("LUA_PATH", env, 1);

  /* LUA_CPATH = $MAD_CPATH;$LUA_CPATH */
  ensure(carg == 4 || carg == 8, "invalid number of path argument");
  if (carg == 4)
    len = snprintf(env, sizeof env, cpath, /* next args discarded without %s */
                   home_path, prog_path, prog_path, prog_path);
  else
    len = snprintf(env, sizeof env, cpath, /* next args discarded without %s */
                   home_path, home_path, prog_path, prog_path,
                   prog_path, prog_path, prog_path, prog_path);
  if (lcpath) strcat(env, lcpath); else if (len>0) env[len-1] = nul;
#if LUAJIT_OS == LUAJIT_OS_WINDOWS
  winpath(env); /* canonize */
#endif
  setenv("LUA_CPATH", env, 1);
}

/* MAD initialization. */
static int handle_madinit (lua_State *L)
{
#if LJ_TARGET_CONSOLE
  str_t init = NULL;
#else
  str_t init = getenv("MAD_INIT");
#endif
  if (init == NULL)
    return 0;  /* status OK */
  else if (init[0] == '@')
    return dofile(L, init+1);
  else
    return dostring(L, init, "=MAD_INIT");
}

/* Extra integrated libs to load. */
LUALIB_API int luaopen_lpeg (lua_State *L);
LUALIB_API int luaopen_lfs (lua_State *L);

static void mad_openlibs (lua_State *L)
{
  static struct { str_t name; int(*func)(lua_State*); } libs[] = {
    { "lpeg", luaopen_lpeg },
    { "lfs", luaopen_lfs },
    { NULL  , NULL },
  };

  for (int i=0; libs[i].name; i++) {
    lua_pushcfunction(L, libs[i].func);
    lua_pushstring(L, libs[i].name);
    lua_call(L, 1, 0);
  }
}

/* Handle signals */
static const int   sig_i[] = { SIGSEGV , SIGFPE , SIGPIPE , SIGBUS , SIGILL , SIGABRT };
static const str_t sig_s[] = {"SIGSEGV","SIGFPE","SIGPIPE","SIGBUS","SIGILL","SIGABRT"};
static const str_t sig_m[] = {"segmentation fault!", "floating-point exception!",
                              "broken pipe!", "bus error!", "illegal instruction!",
                              "abort!"};
enum { sig_ni = sizeof sig_i / sizeof *sig_i,
       sig_ns = sizeof sig_s / sizeof *sig_s,
       sig_nm = sizeof sig_m / sizeof *sig_m,
       static_assert__sig_s_len = 1/(sig_ni == sig_ns),
       static_assert__sig_m_len = 1/(sig_ni == sig_nm) };

static void mad_signal(int sig)
{

  int i = 0;
  while (i < sig_ni && sig_i[i] != sig) ++i;

  // NOT SAFE! (but let's try...)
  fflush(stdout);
  fprintf(stderr, "%s: %s\n", sig_s[i], sig_m[i]);
  lua_pushstring(globalL, "");
  lua_error(globalL);
  exit(EXIT_FAILURE); /* never reached */
}

// SIGINT forward decl
static void laction(int i);

static void mad_setsignal (void)
{
  ensure(signal(SIGINT, laction) != SIG_ERR,
         "unable to set signal hanlder %s", "SIGINT");

  for (int i=0; i < sig_ni; i++)
    ensure(signal(sig_i[i], mad_signal) != SIG_ERR,
           "unable to set signal hanlder %s", sig_s[i]);
}

static void mad_regfunvar (void)
{
  mad_stdin_is_tty = lua_stdin_is_tty();
  lua_register(globalL, "warn" , mad_luawarn );
  lua_register(globalL, "trace", mad_luatrace);
}

/* --- MAD (end) -------------------------------------------------------------*/

#if !LJ_TARGET_CONSOLE
static int sigint_count = 0;

static void lstop(lua_State *L, lua_Debug *ar)
{
  (void)ar;  /* unused arg. */
  sigint_count = 0;
  lua_sethook(L, NULL, 0, 0);
  /* Avoid luaL_error -- a C hook doesn't add an extra frame. */
  luaL_where(L, 0);
  lua_pushfstring(L, "%sinterrupted!", lua_tostring(L, -1));
  lua_error(L);
}

static void laction(int i)
{
  /* protect against multiple Ctrl-C in interactive mode */
  signal(i, laction);
  if (sigint_count++)     // 2nd SIGINT, sethook was called but not efficient
    lstop(globalL, NULL); // so we are certainly in compiled code, call lstop

  lua_sethook(globalL, lstop, LUA_MASKCALL | LUA_MASKRET | LUA_MASKCOUNT, 1);
}
#endif

static void print_usage(void)
{
  fputs("usage: ", stderr);
  fputs(progname, stderr);
  fputs(" [options]... [script [args]...].\n"
  "Available options are:\n"
  "  -e chunk  Execute string " LUA_QL("chunk") ".\n"
  "  -l name   Require library " LUA_QL("name") ".\n"
  "  -b ...    Save or list bytecode.\n"
  "  -j cmd    Perform JIT control command.\n"
  "  -O[opt]   Control JIT optimizations.\n"
  "  -i        Enter interactive mode after executing " LUA_QL("script") ".\n"
  "  -q        Do not show version information.\n"
  "  -M        Do not load MAD environment.\n"
  "  -Mt[=num] Set initial MAD trace level to " LUA_QL("num") ".\n"
  "  -MT[=num] Set initial MAD trace level to " LUA_QL("num") " and location.\n"
  "  -E        Ignore environment variables.\n"
  "  --        Stop handling options.\n"
  "  -         Execute stdin and stop handling options.\n", stderr);
  fflush(stderr);
}

static void l_message(const char *pname, const char *msg)
{
  if (pname) { fputs(pname, stderr); fputc(':', stderr); fputc(' ', stderr); }
  fputs(msg, stderr); fputc('\n', stderr);
  fflush(stderr);
}

static int report(lua_State *L, int status)
{
  if (status && !lua_isnil(L, -1)) {
    const char *msg = lua_tostring(L, -1);
    if (msg == NULL) msg = "(error object is not a string)";
    l_message(progname, msg);
    lua_pop(L, 1);
  }
  return status;
}

static int traceback(lua_State *L)
{
  if (!lua_isstring(L, 1)) { /* Non-string error object? Try metamethod. */
    if (lua_isnoneornil(L, 1) ||
        !luaL_callmeta(L, 1, "__tostring") ||
        !lua_isstring(L, -1))
      return 1;  /* Return non-string error object. */
    lua_remove(L, 1);  /* Replace object by result of __tostring metamethod. */
  }
  luaL_traceback(L, L, lua_tostring(L, 1), 1);
  return 1;
}

static int docall(lua_State *L, int narg, int clear)
{
  int status;
  int base = lua_gettop(L) - narg;  /* function index */
  lua_pushcfunction(L, traceback);  /* push traceback function */
  lua_insert(L, base);  /* put it under chunk and args */
#if !LJ_TARGET_CONSOLE
//  signal(SIGINT, laction);
#endif
  status = lua_pcall(L, narg, (clear ? 0 : LUA_MULTRET), base);
#if !LJ_TARGET_CONSOLE
//  signal(SIGINT, SIG_DFL);
#endif
  lua_remove(L, base);  /* remove traceback function */
  /* force a complete garbage collection in case of errors */
  if (status != LUA_OK) lua_gc(L, LUA_GCCOLLECT, 0);
  return status;
}

static void print_version(void)
{
  print_mad_version();
}

static void print_jit_status(lua_State *L)
{
  int n;
  const char *s;
  lua_getfield(L, LUA_REGISTRYINDEX, "_LOADED");
  lua_getfield(L, -1, "jit");  /* Get jit.* module table. */
  lua_remove(L, -2);
  lua_getfield(L, -1, "status");
  lua_remove(L, -2);
  n = lua_gettop(L);
  lua_call(L, 0, LUA_MULTRET);
  fputs(lua_toboolean(L, n) ? "JIT: ON" : "JIT: OFF", stdout);
  for (n++; (s = lua_tostring(L, n)); n++) {
    putc(' ', stdout);
    fputs(s, stdout);
  }
  putc('\n', stdout);
}

static void createargtable(lua_State *L, char **argv, int argc, int argf)
{
  int i;
  lua_createtable(L, argc - argf, argf);
  for (i = 0; i < argc; i++) {
    lua_pushstring(L, argv[i]);
    lua_rawseti(L, -2, i - argf);
  }
  lua_setglobal(L, "arg");
}

static int dofile(lua_State *L, const char *name)
{
  int status = luaL_loadfile(L, name) || docall(L, 0, 1);
  return report(L, status);
}

static int dostring(lua_State *L, const char *s, const char *name)
{
  int status = luaL_loadbuffer(L, s, strlen(s), name) || docall(L, 0, 1);
  return report(L, status);
}

static int dolibrary(lua_State *L, const char *name)
{
  lua_getglobal(L, "require");
  lua_pushstring(L, name);
  return report(L, docall(L, 1, 1));
}

static void write_prompt(lua_State *L, int firstline)
{
  const char *p;
  lua_getfield(L, LUA_GLOBALSINDEX, firstline ? "_PROMPT" : "_PROMPT2");
  p = lua_tostring(L, -1);
  if (p == NULL) p = firstline ? LUA_PROMPT : LUA_PROMPT2;
  fputs(p, stdout);
  fflush(stdout);
  lua_pop(L, 1);  /* remove global */
}

static int incomplete(lua_State *L, int status)
{
  if (status == LUA_ERRSYNTAX) {
    size_t lmsg;
    const char *msg = lua_tolstring(L, -1, &lmsg);
    const char *tp = msg + lmsg - (sizeof(LUA_QL("<eof>")) - 1);
    if (strstr(msg, LUA_QL("<eof>")) == tp) {
      lua_pop(L, 1);
      return 1;
    }
  }
  return 0;  /* else... */
}

static int pushline(lua_State *L, int firstline)
{
  char buf[LUA_MAXINPUT];
  write_prompt(L, firstline);
  if (fgets(buf, LUA_MAXINPUT, stdin)) {
    size_t len = strlen(buf);
    if (len > 0 && buf[len-1] == '\n')
      buf[len-1] = '\0';
    if (firstline && buf[0] == '=')
      lua_pushfstring(L, "return %s", buf+1);
    else
      lua_pushstring(L, buf);
    return 1;
  }
  return 0;
}

static int loadline(lua_State *L)
{
  int status;
  lua_settop(L, 0);
  if (!pushline(L, 1))
    return -1;  /* no input */
  for (;;) {  /* repeat until gets a complete line */
    status = luaL_loadbuffer(L, lua_tostring(L, 1), lua_strlen(L, 1), "=stdin");
    if (!incomplete(L, status)) break;  /* cannot try to add lines? */
    if (!pushline(L, 0)) /* no more input? */
      return -1;
    lua_pushliteral(L, "\n");  /* add a new line... */
    lua_insert(L, -2);  /* ...between the two lines */
    lua_concat(L, 3);  /* join them */
  }
  lua_remove(L, 1);  /* remove line */
  return status;
}

static void dotty(lua_State *L)
{
  int status;
  const char *oldprogname = progname;
  progname = NULL;
  mad_is_interactive = 1;                                              /* MAD */
  while ((status = loadline(L)) != -1) {
    if (status == LUA_OK) status = docall(L, 0, 0);
    report(L, status);
    if (status == LUA_OK && lua_gettop(L) > 0) {  /* any result to print? */
      lua_getglobal(L, "print");
      lua_insert(L, 1);
      if (lua_pcall(L, lua_gettop(L)-1, 0, 0) != 0)
        l_message(progname,
                  lua_pushfstring(L, "error calling " LUA_QL("print") " (%s)",
                  lua_tostring(L, -1)));
    }
  }
  lua_settop(L, 0);  /* clear stack */
  fputs("\n", stdout);
  fflush(stdout);
  progname = oldprogname;
}

static int handle_script(lua_State *L, char **argx)
{
  int status;
  const char *fname = argx[0];
  if (strcmp(fname, "-") == 0 && strcmp(argx[-1], "--") != 0)
    fname = NULL;  /* stdin */
  status = luaL_loadfile(L, fname);
  if (status == LUA_OK) {
    /* Fetch args from arg table. LUA_INIT or -e might have changed them. */
    int narg = 0;
    lua_getglobal(L, "arg");
    if (lua_istable(L, -1)) {
      do {
        narg++;
        lua_rawgeti(L, -narg, narg);
      } while (!lua_isnil(L, -1));
      lua_pop(L, 1);
      lua_remove(L, -narg);
      narg--;
    } else {
      lua_pop(L, 1);
    }
    status = docall(L, narg, 0);
  }
  return report(L, status);
}

/* Load add-on module. */
static int loadjitmodule(lua_State *L)
{
  lua_getglobal(L, "require");
  lua_pushliteral(L, "ljit_");
  lua_pushvalue(L, -3);
  lua_concat(L, 2);
  if (lua_pcall(L, 1, 1, 0)) {
    const char *msg = lua_tostring(L, -1);
    if (msg && !strncmp(msg, "module ", 7))
      goto nomodule;
    return report(L, 1);
  }
  lua_getfield(L, -1, "start");
  if (lua_isnil(L, -1)) {
  nomodule:
    l_message(progname,
      "unknown MAD command or ljit.* modules not installed");
    return 1;
  }
  lua_remove(L, -2);  /* Drop module table. */
  return 0;
}

/* Run command with options. */
static int runcmdopt(lua_State *L, const char *opt)
{
  int narg = 0;
  if (opt && *opt) {
    for (;;) {  /* Split arguments. */
      const char *p = strchr(opt, ',');
      narg++;
      if (!p) break;
      if (p == opt)
        lua_pushnil(L);
      else
        lua_pushlstring(L, opt, (size_t)(p - opt));
      opt = p + 1;
    }
    if (*opt)
      lua_pushstring(L, opt);
    else
      lua_pushnil(L);
  }
  return report(L, lua_pcall(L, narg, 0, 0));
}

/* JIT engine control command: try jit library first or load add-on module. */
static int dojitcmd(lua_State *L, const char *cmd)
{
  const char *opt = strchr(cmd, '=');
  lua_pushlstring(L, cmd, opt ? (size_t)(opt - cmd) : strlen(cmd));
  lua_getfield(L, LUA_REGISTRYINDEX, "_LOADED");
  lua_getfield(L, -1, "jit");  /* Get jit.* module table. */
  lua_remove(L, -2);
  lua_pushvalue(L, -2);
  lua_gettable(L, -2);  /* Lookup library function. */
  if (!lua_isfunction(L, -1)) {
    lua_pop(L, 2);  /* Drop non-function and jit.* table, keep module name. */
    if (loadjitmodule(L))
      return 1;
  } else {
    lua_remove(L, -2);  /* Drop jit.* table. */
  }
  lua_remove(L, -2);  /* Drop module name. */
  return runcmdopt(L, opt ? opt+1 : opt);
}

/* Optimization flags. */
static int dojitopt(lua_State *L, const char *opt)
{
  lua_getfield(L, LUA_REGISTRYINDEX, "_LOADED");
  lua_getfield(L, -1, "jit.opt");  /* Get jit.opt.* module table. */
  lua_remove(L, -2);
  lua_getfield(L, -1, "start");
  lua_remove(L, -2);
  return runcmdopt(L, opt);
}

/* Save or list bytecode. */
static int dobytecode(lua_State *L, char **argv)
{
  int narg = 0;
  lua_pushliteral(L, "bcsave");
  if (loadjitmodule(L))
    return 1;
  if (argv[0][2]) {
    narg++;
    argv[0][1] = '-';
    lua_pushstring(L, argv[0]+1);
  }
  for (argv++; *argv != NULL; narg++, argv++)
    lua_pushstring(L, *argv);
  report(L, lua_pcall(L, narg, 0, 0));
  return -1;
}

/* check that argument has no extra characters at the end */
#define notail(x,i) {if ((x)[i] != '\0') return -1;}

#define FLAGS_INTERACTIVE 1
#define FLAGS_VERSION     2
#define FLAGS_EXEC        4
#define FLAGS_OPTION      8
#define FLAGS_NOENV       16
#define FLAGS_MADENV      128

static int collectargs(char **argv, int *flags)
{
  int i;
  *flags |= FLAGS_VERSION|FLAGS_MADENV;
  for (i = 1; argv[i] != NULL; i++) {
    if (argv[i][0] != '-')  /* Not an option? */
      return i;
    switch (argv[i][1]) {  /* Check option. */
    case '-':
      notail(argv[i],2);
      return i+1;
    case '\0':
      return i;
    case 'i':
      notail(argv[i],2);
      *flags |= FLAGS_INTERACTIVE;
      break;
    case 'q':
      notail(argv[i],2);
      *flags &= ~FLAGS_VERSION;
      break;
    case 'e':
      *flags |= FLAGS_EXEC; /* FALLTHRU */
    case 'j':  /* LuaJIT extension */
    case 'l':
      *flags |= FLAGS_OPTION;
      if (argv[i][2] == '\0') {
        i++;
        if (argv[i] == NULL) return -1;
      }
      break;
    case 'O': break;  /* LuaJIT extension */
    case 'b':  /* LuaJIT extension */
      *flags &= ~(FLAGS_VERSION|FLAGS_MADENV);
      if (*flags) return -1;
      *flags |= FLAGS_EXEC;
      return i+1;
    case 'E':
      notail(argv[i],2);
      *flags |= FLAGS_NOENV;
      break;
    case 'M': /* MAD options */
      if (argv[i][2] == 't' || argv[i][2] == 'T') { /* MAD trace level */
        mad_trace_location = argv[i][2] == 'T';
        mad_trace_level    = argv[i][3] == '=' ? strtol(argv[i]+4, NULL, 0) : 1;
      } else { /* don't load MAD environment */
        notail(argv[i],2);
        *flags &= ~FLAGS_MADENV;
      }
      break;
    default: return -1;  /* invalid option */
    }
  }
  return i;
}

static int runargs(lua_State *L, char **argv, int argn)
{
  int i;
  for (i = 1; i < argn; i++) {
    if (argv[i] == NULL) continue;
    lua_assert(argv[i][0] == '-');
    switch (argv[i][1]) {
    case 'e': {
      const char *chunk = argv[i] + 2;
      if (*chunk == '\0') chunk = argv[++i];
      lua_assert(chunk != NULL);
      if (dostring(L, chunk, "=(command line)") != 0)
        return 1;
      break;
    }
    case 'l': {
      const char *filename = argv[i] + 2;
      if (*filename == '\0') filename = argv[++i];
      lua_assert(filename != NULL);
      if (dolibrary(L, filename))
        return 1;
      break;
    }
    case 'j': {  /* LuaJIT extension. */
      const char *cmd = argv[i] + 2;
      if (*cmd == '\0') cmd = argv[++i];
      lua_assert(cmd != NULL);
      if (dojitcmd(L, cmd))
        return 1;
      break;
    }
    case 'O':  /* LuaJIT extension. */
      if (dojitopt(L, argv[i] + 2))
        return 1;
      break;
    case 'b':  /* LuaJIT extension. */
      return dobytecode(L, argv+i);
    default: break;
    }
  }
  return LUA_OK;
}

static int handle_luainit(lua_State *L)
{
#if LJ_TARGET_CONSOLE
  const char *init = NULL;
#else
  const char *init = getenv(LUA_INIT);
#endif
  if (init == NULL)
    return LUA_OK;
  else if (init[0] == '@')
    return dofile(L, init+1);
  else
    return dostring(L, init, "=" LUA_INIT);
}

static struct Smain {
  char **argv;
  int argc;
  int status;
} smain;

static int pmain(lua_State *L)
{
  struct Smain *s = &smain;
  char **argv = s->argv;
  int argn;
  int flags = 0;
  globalL = L;
  if (argv[0] && argv[0][0]) progname = argv[0];

  LUAJIT_VERSION_SYM();  /* Linker-enforced version check. */

  argn = collectargs(argv, &flags);
  if (argn < 0) {  /* Invalid args? */
    print_usage();
    s->status = 1;
    return 0;
  }

  if ((flags & FLAGS_NOENV)) {
    lua_pushboolean(L, 1);
    lua_setfield(L, LUA_REGISTRYINDEX, "LUA_NOENV");
  }

  /* Set MAD env _before_ libraries are open. */
  mad_setenv(L, flags & FLAGS_NOENV);

  /* Stop collector during library initialization. */
  lua_gc(L, LUA_GCSTOP, 0);
  luaL_openlibs(L);
  mad_openlibs(L);
  lua_gc(L, LUA_GCRESTART, -1);

  createargtable(L, argv, s->argc, argn);

  if (!(flags & FLAGS_NOENV)) {
    s->status = handle_luainit(L);
    if (s->status != LUA_OK) return 0;
  }

  /* MAD section. */
  mad_setsignal();
  mad_regfunvar();

  if ((flags & FLAGS_MADENV))
    dolibrary(L, "madl_main");
  if (!(flags & FLAGS_NOENV)) {
    s->status = handle_madinit(L);
    if (s->status != LUA_OK) return 0;
  }

  if ((flags & FLAGS_VERSION)) print_version();

  s->status = runargs(L, argv, argn);
  if (s->status != LUA_OK) return 0;

  if (s->argc > argn) {
    s->status = handle_script(L, argv + argn);
    if (s->status != LUA_OK) return 0;
  }

  if ((flags & FLAGS_INTERACTIVE)) {
    (void)print_jit_status;
    if ((flags & FLAGS_MADENV))
      dostring(L, "MAD.strict(false)", "=(interactive setup)");
    dotty(L);
  } else if (s->argc == argn && !(flags & FLAGS_EXEC)) {
    if (lua_stdin_is_tty()) {
      (void)print_version;
      (void)print_jit_status;
      if ((flags & FLAGS_MADENV))
        dostring(L, "MAD.strict(false)", "=(interactive setup)");
      dotty(L);
    } else {
      dofile(L, NULL);  /* Executes stdin as a file. */
    }
  }
  return 0;
}

int main(int argc, char **argv)
{
//  fprintf(stderr, "main=0x%p\n", (void*)main);
  int status;

  lua_State *L = lua_open();
  if (L == NULL) {
    l_message(argv[0], "cannot create state: not enough memory");
    return EXIT_FAILURE;
  }
  smain.argc = argc;
  smain.argv = argv;
  status = lua_cpcall(L, pmain, NULL);
  report(L, status);
  lua_close(L);
  return (status || smain.status > 0) ? EXIT_FAILURE : EXIT_SUCCESS;
}

