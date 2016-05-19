#ifndef MAD_LOG_PRIV_H
#define MAD_LOG_PRIV_H

/*
 o----------------------------------------------------------------------------o
 |
 | Logging module private implementation
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

#define mad_error(...)     mad_savtrcfun(  mad_error(__VA_ARGS__))
#define mad_warn(...)      mad_savtrcfun(  mad_warn (__VA_ARGS__))
#define mad_trace(l,...)   mad_savchkfun(l,mad_trace(__VA_ARGS__))
#define mad_ensure(c,...)  mad_savtstfun(c,mad_error(__VA_ARGS__))

#define mad_savtrcfun(f) \
  (mad_trace_function = __func__, f)

#define mad_savchkfun(l,f) \
  ((void)(mad_trace_level >= (l) && (mad_savtrcfun(f),0)))

#define mad_savtstfun(c,f) \
  ((void)(!(c) && (mad_savtrcfun(f),0)))

void (mad_error) (str_t, ...) __attribute__((format(printf,1,2),noreturn));
void (mad_warn)  (str_t, ...) __attribute__((format(printf,1,2)));
void (mad_trace) (str_t, ...) __attribute__((format(printf,1,2)));

// ---------------------------------------------------------------------------o

#endif // MAD_LOG_PRIV_H
