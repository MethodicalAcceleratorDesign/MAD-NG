#ifndef MAD_LOG_H
#define MAD_LOG_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Logging module interface
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
  - MAD log handlers: mad_(error, warn, trace, ensure)
    shortcuts: error, warn, trace, ensure
    all are macros that set automatically the location of the caller

  Information:
  - error, warn, trace print to stderr; ensure = cond + error
  - warn and trace use tags 'warning:' and 'trace:'
  - mad_trace_level:
  	trace with level <= mad_trace_level are active
  - mad_trace_location:
    >0 location is traced by error, warn and trace
     0 location is traced by error (only)

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"
#include "mad_main.h"

// --- interface -------------------------------------------------------------o

#define error(...)        mad_error (     __VA_ARGS__)
#define warn(...)         mad_warn  (     __VA_ARGS__)
#define trace(lvl,...)    mad_trace (lvl, __VA_ARGS__)
#define ensure(cond,...)  mad_ensure(cond,__VA_ARGS__)

// --- implementation (private) ----------------------------------------------o

#define mad_error(...)                  (mad_error)(  mad_logloc_,__VA_ARGS__)
#define mad_warn(...)                   (mad_warn )(  mad_logloc_,__VA_ARGS__)
#define mad_trace(l,...)  mad_loglvl_(l,(mad_trace)(l,mad_logloc_,__VA_ARGS__))
#define mad_ensure(c,...) mad_logcnd_(c,(mad_error)(  mad_logloc_,__VA_ARGS__))

#define mad_loglvl_(l,f) ((void)(mad_trace_level >= (l) && (f,0)))
#define mad_logcnd_(c,f) ((void)(                  !(c) && (f,0)))
#define mad_logloc_      __FILE__ ":" MKSTR(__LINE__) ": "

// ---------------------------------------------------------------------------o

#endif // MAD_LOG_H
