#ifndef MAD_LOG_IMPL_H
#define MAD_LOG_IMPL_H
#else
#error "implementation header, do not include this file directly"
#endif

/*
 o----------------------------------------------------------------------------o
 |
 | Logging module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2015
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o----------------------------------------------------------------------------o
 | You can redistr_tibute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distr_tibuted in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o----------------------------------------------------------------------------o
*/

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
