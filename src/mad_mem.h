#ifndef MAD_MEM_H
#define MAD_MEM_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Memory module interface
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
  - fast memory allocator (per-thread pool).
  - MAD memory handlers: mad_malloc, mad_realloc, mad_free
  - allocated memory can be used-by/moved-to any thread (global allocator).
  - temporary allocation are only for local use (scoped), and the length
    corresponds to the number of elements of type 'type' in the buffer.

  Information:
  - parameters ending with an underscope can be null (i.e. optional).
  - mad_calloc calls mad_malloc and set to zeros the allocated memory.
  - mad_mcached returns the amount of memory cached (slow).
  - mad_mcollect frees the cached memory and returns its amount (slow).
  - temporay buffers can be either on the stack or allocated with mad_malloc
    depending on their size, and must _always_ be locally freed.
  - defining MAD_MEM_STD replaces mad allocator by C allocator

  Errors:
  - mad_realloc or mad_free on pointers not returned by mad_malloc or
    mad_realloc is undefined.

 o-----------------------------------------------------------------------------o
 */

#include "mad_log.h"

// --- interface --------------------------------------------------------------o

// local buffer (macros)
#define mad_alloc_tmp(type,name,length)
#define mad_free_tmp(      name)

// allocators
void*  mad_malloc  (size_t size);
void*  mad_calloc  (size_t count, size_t size );
void*  mad_realloc (void*  ptr_ , size_t size_);
void   mad_free    (void*  ptr_);

// utils
size_t mad_msize    (void*  ptr_);
size_t mad_mcached  (void);
size_t mad_mcollect (void);

// ----------------------------------------------------------------------------o
// --- implementation (private) -----------------------------------------------o
// ----------------------------------------------------------------------------o

#define mad_malloc(s)    mad_mcheck(__func__, mad_malloc (s)  )
#define mad_calloc(c,s)  mad_mcheck(__func__, mad_calloc (c,s))
#define mad_realloc(p,s) mad_mcheck(__func__, mad_realloc(p,s))

void*  (mad_malloc)  (size_t)         __attribute__((hot,malloc));
void*  (mad_calloc)  (size_t, size_t) __attribute__((hot,malloc));
void*  (mad_realloc) (void* , size_t) __attribute__((hot));
void    mad_free     (void*)          __attribute__((hot));

static inline void*
mad_mcheck (str_t fname, void *ptr_)
{
  if (!ptr_)
    (mad_error)(fname, "invalid null pointer (out of memory?)");

  return ptr_;
}

#undef  mad_alloc_tmp
#define mad_alloc_tmp(T,NAME,L) \
  T NAME##_local_tmp__[sizeof(T)*(L) < 8192 ? (L) : 1]; \
  T *NAME = (sizeof(T)*(L) < 8192 ? \
             NAME##_local_tmp__ : mad_malloc(sizeof(T)*(L)) )

#undef  mad_free_tmp
#define mad_free_tmp(NAME) \
  (NAME != NAME##_local_tmp__ ? mad_free(NAME) : (void)0)

// --- end --------------------------------------------------------------------o

#endif // MAD_MEM_H
