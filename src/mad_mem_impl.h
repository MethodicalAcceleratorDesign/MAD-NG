#ifndef MAD_MEM_IMPL_H
#define MAD_MEM_IMPL_H
#else
#error "implementation header, do not include this file directly"
#endif

/*
 o----------------------------------------------------------------------------o
 |
 | Memory module implementation
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

#define mad_malloc(s)    (mad_savtrcloc(2), mad_malloc(s))
#define mad_realloc(p,s) (mad_savtrcloc(2), mad_realloc(p,s))
#define mad_free(p)      (mad_savtrcloc(2), mad_free(p))

void* (mad_malloc)  (size_t)        __attribute__((hot,malloc));
void* (mad_realloc) (void*, size_t) __attribute__((hot));
void  (mad_free)    (void*)         __attribute__((hot));
size_t mad_mem_size (void*)         __attribute__((hot,const));

// note: 2048 == mblk_max from mad_mem.c

#define mad_alloc_tmp(T,NAME,L) \
  T NAME##_local_tmp__[2048/sizeof(T)]; \
  T *NAME = ((size_t)(L) > (2048/sizeof(T)) ? \
            mad_malloc((size_t)(L) * sizeof(T)) : NAME##_local_tmp__)

#define mad_free_tmp(NAME) \
  (NAME != NAME##_local_tmp__ ? mad_free(NAME) : (void)0)


// ---------------------------------------------------------------------------o
