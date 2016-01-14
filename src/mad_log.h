#ifndef MAD_LOG_H
#define MAD_LOG_H

/*
 o----------------------------------------------------------------------------o
 |
 | Logging module interface
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
  - MAD log handlers: mad_(fatal, error, warn, info, debug, ensure)
    shortcuts: fatal, error, warn, info, debug, ensure
    all are macros that set automatically the location of the caller
  - non-variadic log handlers for foreign languages:
    mad_(fatal, error, warn, info, debug)  (same name but not macros)
    use mad_log_setloc  to set caller location from the application
    use mad_log_setloc1 to set caller location from a scripting language

  Information:
  - fatal, error, warn, debug print to stderr; info prints to stdout
  - fatal, error, warn and debug use tags 'Fatal', Error', 'Warning' and 'Debug'
  - info uses tag 'Info' only if the location is traced
  - only MAD handlers save implicitly the location of the caller
  - mad_info_level =
     0 fatal, error and warn are active (default)
    >0 fatal, error, warn and info with level <= mad_info_level are active
  - mad_debug_level =
  	 0 fatal, error and warn are active (default)
  	>0 fatal, error, warn and debug with level <= mad_debug_level are active
  - mad_trace_location =
  	0  location is traced by fatal (default)
  	1  location is traced by debug and fatal
    2  location is traced by debug, warn, error and fatal
    >2 location is traced by debug, info, warn, error and fatal

 o----------------------------------------------------------------------------o
 */

#include "mad.h"

// --- interface -------------------------------------------------------------o

#define fatal(...)         mad_fatal (     __VA_ARGS__)
#define error(...)         mad_error (     __VA_ARGS__)
#define warn(...)          mad_warn  (     __VA_ARGS__)
#define info(lvl,...)      mad_info  (lvl, __VA_ARGS__)
#define debug(lvl, ...)    mad_debug (lvl, __VA_ARGS__)
#define ensure(cond, ...)  mad_ensure(cond,__VA_ARGS__)

// interface for manual call, don't set location
void mad_fatalf  (str_t fmt, ...);
void mad_errorf  (str_t fmt, ...);
void mad_warnf   (str_t fmt, ...);
void mad_infof   (int lvl, str_t fmt, ...);
void mad_debugf  (int lvl, str_t fmt, ...);

// interface only for ffi without variadic support (!= macros)
void mad_fatal   (str_t msg);
void mad_error   (str_t msg);
void mad_warn    (str_t msg);
void mad_info    (int lvl, str_t msg);
void mad_debug   (int lvl, str_t msg);

// utils
void mad_log_setloc  (str_t file, int line);
void mad_log_setloc1 (str_t file, int line); // foreign location
void mad_log_getstat (long *error_, long *warn_, long *info_, long *debug_);

// automatic location
#define mad_log_saveloc() mad_log_setloc(__FILE__, __LINE__)

// --- globals ---------------------------------------------------------------o

extern int mad_info_level;
extern int mad_debug_level;
extern int mad_inside_script;
extern int mad_trace_location;

// ---------------------------------------------------------------------------o

#endif // MAD_LOG_H

/*\

 |
 |
 |
 |
 |

 P
 R
 I
 V
 A
 T
 E

 I
 M
 P
 L
 E
 M
 E
 N
 T
 A
 T
 I
 O
 N

 |
 |
 |
 |
 |

\*/

#ifndef MAD_LOG_IMPL_H
#define MAD_LOG_IMPL_H

// --- implementation (private) ----------------------------------------------o

#define mad_fatal(...)         (mad_savtrcloc(0),mad_fatalf(  __VA_ARGS__))
#define mad_error(...)         (mad_savtrcloc(2),mad_errorf(  __VA_ARGS__))
#define mad_warn(...)          (mad_savtrcloc(2),mad_warnf (  __VA_ARGS__))
#define mad_info(l,...) ((void)(mad_chkmsglvl(l) && \
                               (mad_savtrcloc(3),mad_infof (l,__VA_ARGS__),0)))
#define mad_debug(l,...)((void)(mad_chkerrlvl(l) && \
                               (mad_savtrcloc(1),mad_debugf(l,__VA_ARGS__),0)))

#define mad_ensure(c,...)((void)(!(c) && (mad_fatal(__VA_ARGS__),0)))

void (mad_fatal)   (str_t)      __attribute__((noreturn));
void mad_fatalf    (str_t, ...) __attribute__((format(printf,1,2),noreturn));
void mad_errorf    (str_t, ...) __attribute__((format(printf,1,2)));
void mad_warnf     (str_t, ...) __attribute__((format(printf,1,2)));
void mad_infof (int,str_t, ...) __attribute__((format(printf,2,3)));
void mad_debugf(int,str_t, ...) __attribute__((format(printf,2,3)));

void mad_log_setloc (str_t, int) __attribute__((hot));

#define mad_savtrcloc(l) \
  ((void)(mad_trace_location >= (l) && (mad_log_saveloc(),0)))

#define mad_chkmsglvl(l) \
  (mad_info_level >= (l))

#define mad_chkerrlvl(l) \
  (mad_debug_level >= (l))

// ---------------------------------------------------------------------------o

#endif // MAD_LOG_IMPL_H
