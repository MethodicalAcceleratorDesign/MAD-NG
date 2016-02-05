#ifndef MAD_MEM_PRIV_H
#define MAD_MEM_PRIV_H

/*
 o----------------------------------------------------------------------------o
 |
 | Memory module private implementation
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

#include "mad_log.h"

#ifndef MAD_MEM_STD

#define mad_malloc(s)    (mad_savtrcloc(2), mad_malloc (s)  )
#define mad_calloc(c,s)  (mad_savtrcloc(2), mad_calloc (c,s))
#define mad_realloc(p,s) (mad_savtrcloc(2), mad_realloc(p,s))
#define mad_free(p)      (mad_savtrcloc(2), mad_free   (p)  )

#else

#include <stdlib.h>
#define mad_malloc(s)    (mad_savtrcloc(2), mad_mem_check(malloc (s))  )
#define mad_calloc(c,s)  (mad_savtrcloc(2), mad_mem_check(calloc (c,s)))
#define mad_realloc(p,s) (mad_savtrcloc(2), mad_mem_check(realloc(p,s)))
#define mad_free(p)      (                                free   (p)   )

#endif // MAD_MEM_STD

void* (mad_malloc)   (size_t)         __attribute__((hot,malloc));
void* (mad_calloc)   (size_t, size_t) __attribute__((hot,malloc));
void* (mad_realloc)  (void* , size_t) __attribute__((hot));
void  (mad_free)     (void*)          __attribute__((hot));
size_t mad_mem_size  (void*)          __attribute__((hot,const));
void*  mad_mem_check (void*)          __attribute__((hot,const));

#undef  mad_alloc_tmp
#define mad_alloc_tmp(T,NAME,L) \
  T NAME##_local_tmp__[8192/sizeof(T)]; \
  T *NAME = ((size_t)(L) > (8192/sizeof(T)) ? \
            mad_malloc((size_t)(L) * sizeof(T)) : NAME##_local_tmp__)

#undef  mad_free_tmp
#define mad_free_tmp(NAME) \
  (NAME != NAME##_local_tmp__ ? mad_free(NAME) : (void)0)

// ---------------------------------------------------------------------------o

#endif // MAD_MEM_PRIV_H