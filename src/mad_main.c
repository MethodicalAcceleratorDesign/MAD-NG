/*
 o----------------------------------------------------------------------------o
 |
 | Main module implementation
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

#include <lua.h>
#include <lualib.h>
#include <lauxlib.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#include "mad_main.h"
#include "mad_log.h"

#ifdef POSIX_VERSION
#include <unistd.h>
#define stdin_isatty()   isatty( fileno(stdin))
#else // WIN_VERSION
#include <io.h>
#define stdin_isatty()  _isatty(_fileno(stdin))
#endif

// --- constants -------------------------------------------------------------o

#define VERSION "0.0.0"
#define OSNAME  "Linux64"
// "Linux64"
// "Linux32"
// "MacOSX"
// "Win64"
// "Win32"

#define MAD_MAXINPUTLEN 512

// --- helpers ---------------------------------------------------------------o

static lua_State *mad_luastate;
static str_t      mad_progname;

void
mad_lua_setloc (int level)
{
  static char buf[FILENAME_MAX];
  lua_State *L = mad_luastate;
  ensure(L, "invalid script location");
  lua_Debug ar;
  lua_getstack(L, level, &ar);
  lua_getinfo (L, "Sl", &ar);

  strncpy(buf, ar.source+1, sizeof buf);
  int line = ar.currentline;
  debug(4, "lua_setloc: file='%s', line='%d'", buf, line);
  mad_log_setloc1(buf, line);
}

void
mad_fatal_extrn (void)
{
  debug(4, "lua_fatal: entering mad_fatal_extrn");
  lua_State *L = mad_luastate;
  ensure(L, "invalid script error");
  lua_pushstring(L, "");
  debug(4, "lua_fatal: calling lua_error");
  lua_error(L);  // no return, but crash!!!
  exit(EXIT_FAILURE); // make gcc happy (never reached)
}

// --- locals ----------------------------------------------------------------o

static struct option {
  int iscr, nopt, narg, argc;
  int quiet, interact;
  str_t input;
} option[1] = {{0}};

static void
banner(void)
{
  static str_t msg = 
  "    ____  __   ______    ______     |   Methodical Accelerator Design\n"
  "     /  \\/  \\   /  _  \\   /  _  \\   |   release: " VERSION " (" OSNAME ")\n"
  "    /  __   /  /  /_/ /  /  /_/ /   |   support: http://cern.ch/mad\n"
  "   /__/  /_/  /__/ /_/  /_____ /    |   licence: GPLv3 (C) CERN 2015\n"
  "                                    |   started: %s\n"
  "\n";

  time_t t = time(NULL);
  struct tm *tm = localtime(&t);
  char buf[80];

  strftime(buf, sizeof buf, "%F %T", tm);
  printf(msg, buf);
}

static void
usage(void)
{
  static str_t msg =
  "usage: %s [options]... [script [args]...]\n"
  "  -e chunck   execute string 'chunck'\n"
  "  -l name     load library 'name'\n"
  "  -i          run interactive mode\n"
  "  -q          set quiet mode (no banner)\n"
  "  -v[level]   set info  'level', default is 0\n"
  "  -d[level]   set debug 'level', default is 0\n"
  "  -t[level]   set trace location 'level', default is 0\n"
  "  -           set input to <stdin>\n"
  "\n";

  fprintf(stderr, msg, mad_progname);
  exit(EXIT_FAILURE);
}

static void
set_level(str_t level, int *target)
{
  char *end;
  int   lvl = *level ? strtol(level, &end, 10) : 1;
  if ((*level && end && *end) || lvl < 0 || lvl > 100 ) usage();
  *target = lvl;
}

static void
check_args(int argc, char *argv[])
{
  int i, j=0;

  for (i=1; argv[i] && !j; i++) {
    if (argv[i][0] != '-') { j=i; continue; } // not an option

    switch(argv[i][1]) {
    case '\0': j=i; break;// stdin
    case 'v' : set_level(&argv[i][2], &mad_info_level    ); break;
    case 'd' : set_level(&argv[i][2], &mad_debug_level   ); break;
    case 't' : set_level(&argv[i][2], &mad_trace_location); break;
    default  :
      if (argv[i][2]) usage();      
      switch(argv[i][1]) {
      case 'q': option->quiet    = 1; break;
      case 'i': option->interact = 1; break;
      case 'e':
      case 'l': if (argv[++i]) break;
      default : usage();
      }
    }
  }

  option->argc  = argc;
  option->iscr  = j ? j : argc;
  option->nopt  = j ? j-1 : argc-1;
  option->narg  = j ? argc - (j+1) : 0;
  option->input = j ? argv[j] : "-";

  if (!option->quiet) banner();

  if (mad_debug_level >= 3) {
    debug(4, "prog: %s", mad_progname);
    debug(4, "#opt: %d", option->nopt);
    for (int i=1; i <= option->nopt; i++)
      debug(4, "opt%d: %s", i, argv[i]);
  
    debug(4, "scr%d: %s", option->iscr, option->input);
    debug(4, "#arg: %d", option->narg);
    for (int i=1; i <= option->narg; i++)
      debug(4, "arg%d: %s", i, argv[option->iscr+i]);

    debug(4, "verbose: level %d", mad_info_level);
    debug(4, "debug  : level %d", mad_debug_level);
    debug(4, "trace  : level %d", mad_trace_location);
  }
}

static void
make_arg(lua_State *L, char *argv[])
{
  int argc = option->argc;
  int iscr = option->iscr;

  // ensure enough space on the stack
  luaL_checkstack(L, option->narg+3, "too many arguments to script");

  // push args of script for use as '...'
  for (int i=iscr+1; i < argc; i++)
    lua_pushstring(L, argv[i]);

  lua_createtable(L, option->narg+1, option->nopt+1);

  // push all args for use as 'arg[]'
  for (int i=0; i < argc; i++) {
    lua_pushstring(L, argv[i]);
    lua_rawseti(L, -2, i - iscr);
  }
  lua_setglobal(L, "arg");
}

static void
load_file(lua_State *L, str_t file)
{
  debug(4, "loading script '%s'", file);

  if (!strcmp(file, "-")) file = NULL; // stdin

  if (luaL_loadfile(L, file)) {
    fprintf(stderr, "%s\n", lua_tostring(L, -1));
    lua_pop(L, 1);
    exit(EXIT_FAILURE);
  }
}

static void
load_string(lua_State *L, str_t str)
{
  debug(4, "loading string '%s'", str);

  if (luaL_loadstring(L, str)) {
    fprintf(stderr, "%s\n", lua_tostring(L, -1));
    lua_pop(L, 1);
    exit(EXIT_FAILURE);
  }
}

static void
run_chunk(lua_State *L, str_t chunk, int narg)
{
  debug(4, "running script '%s' with %d arg", chunk, narg);

  if (lua_pcall(L, narg, 0, 0)) {
    fprintf(stderr, "%s\n", lua_tostring(L, -1));
    lua_pop(L, 1);
    exit(EXIT_FAILURE);
  } 
}

#if 0
static void
interactive(lua_State *L)
{
  printf("NYI\n");
  (void)L;
}

static void
exec_args(lua_State *L, char *argv[])
{
  (void)n;
  int do_interactive = FALSE;
  if (do_interactive) banner();

  for (int i=1; argv[i]; i++) {

    if (!strcmp(argv[i], "-v")) {
      continue;
    }

    if (!strcmp(argv[i], "-i")) {
      interactive(L);
      continue;
    }

    if (!strcmp(argv[i], "-e")) {
      load_string(L, argv[++i]);
      continue;
    }

    if (!strcmp(argv[i], "-l")) {
      char buf[MAXINPUTLEN+25];
      snprintf(buf, sizeof(buf), "require '%s'", argv[++i]);
      load_string(L, buf);
      continue;
    }

    if (!strcmp(argv[i], "--"))
      ++i;

    load_file(L, argv[i]);
  }
}
#endif

static void
run_init(lua_State *L)
{
  str_t init = getenv("MAD_INIT");
  if (init) {
    if (*init == '@')
      load_file  (L, init+1);
    else
      load_string(L, init);
  }
}

static void
at_exit(void)
{
  debug(4, "closing Lua");
  lua_close(mad_luastate);
}

int
main(int argc, char *argv[])
{
  lua_State *L = luaL_newstate();
  mad_luastate = L;
  mad_progname = argv[0];

  if (!L) 
    fatal("%s - unable to create Lua state (out of memory)\n", argv[0]);
  if (atexit(at_exit))
    fatal("%s - unable to save exit handler\n", argv[0]);

  luaL_openlibs(L);

  check_args(argc, argv);
  run_init  (L);
  load_file (L, option->input);
  make_arg  (L, argv);
  run_chunk (L, option->input, option->narg);
//  exec_args(L   , argv);

  return EXIT_SUCCESS;
}