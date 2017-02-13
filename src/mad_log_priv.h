#ifndef MAD_LOG_PRIV_H
#define MAD_LOG_PRIV_H

/*
 o----------------------------------------------------------------------------o
 |
 | Logging module private implementation
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
*/

#define mad_error(...)                     (mad_error)(  __func__, __VA_ARGS__)
#define mad_warn(...)                      (mad_warn )(  __func__, __VA_ARGS__)
#define mad_trace(l,...)   mad_logchkfun(l,(mad_trace)(l,__func__, __VA_ARGS__))
#define mad_ensure(c,...)  mad_logtstfun(c,(mad_error)(  __func__, __VA_ARGS__))

#define mad_logchkfun(l,f) ((void)(mad_trace_level >= (l) && (f,0)))
#define mad_logtstfun(c,f) ((void)(                  !(c) && (f,0)))

void (mad_error) (    str_t,str_t,...) __attribute__((format(printf,2,3),noreturn));
void (mad_warn)  (    str_t,str_t,...) __attribute__((format(printf,2,3)));
void (mad_trace) (int,str_t,str_t,...) __attribute__((format(printf,3,4)));

// ---------------------------------------------------------------------------o

#endif // MAD_LOG_PRIV_H
