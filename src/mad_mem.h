#ifndef MAD_MEM_H
#define MAD_MEM_H

/*
 o----------------------------------------------------------------------------o
 |
 | Memory module interface
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
  - fast memory allocator (per-thread pool).
  - MAD memory handlers: mad_malloc, mad_realloc, mad_free
  - allocated memory can be used-by/moved-to any thread (global allocator).
  - temporary allocation are only for local use (scoped), and the length
    corresponds to the number of elements of type 'type' in the buffer.
 
  Information:
  - parameters ending with an underscope can be null.
  - mad_malloc and mad_realloc call mad_fatal with caller location if
    (re)allocation fails instead of returning a NULL pointer.
  - mad_calloc calls mad_malloc and set to zeros the allocated memory.
  - mad_mem_cached returns the amount of memory cached (slow).
  - mad_mem_collect frees the cached memory and returns its amount (slow).
  - temporay buffers can be either on the stack or allocated with mad_malloc
    depending on their size, and must _always_ be locally freed.
  - defining MAD_MEM_STD replaces mad allocator by C allocator

  Errors:
  - mad_realloc or mad_free on pointers not returned by mad_malloc or
    mad_realloc is undefined.

 o----------------------------------------------------------------------------o
 */

#include "mad.h"

// --- interface -------------------------------------------------------------o

// local buffer
#define mad_alloc_tmp(type,name,length)
#define mad_free_tmp(name)

// allocator (note: no calloc!)
void*  mad_malloc   (size_t size);
void*  mad_calloc   (size_t count, size_t size );
void*  mad_realloc  (void  *ptr_ , size_t size_);
void   mad_free     (void  *ptr_);

// utils
size_t mad_msize    (void *ptr_);
size_t mad_mcached  (void);
size_t mad_mcollect (void);

// --- implementation (private) ----------------------------------------------o

#include "mad_mem_priv.h"

// ---------------------------------------------------------------------------o

#endif // MAD_MEM_H
