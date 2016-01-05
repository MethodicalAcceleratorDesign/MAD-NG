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
    shortcuts: malloc, realloc, free
  - allocated memory can be used-by/moved-to any thread (global allocator).
  - temporary buffer are only for local use (scoped), and the length
    corresponds to the number of elements of type 'type' in the buffer.
 
  Information:
  - parameters ending with an underscope can be null.
  - mad_malloc and mad_realloc call mad_fatal with caller location if
    (re)allocation fails instead of returning a NULL pointer.
  - mad_mem_cached returns the amount of memory cached (slow).
  - mad_mem_collect frees the cached memory and returns its amount (slow).
  - temporay buffers can be either on the stack or allocated with mad_malloc
    depending on their size, and must always be locally freed.

  Errors:
  - mad_realloc or mad_free on pointers not returned by mad_malloc or
    mad_realloc is undefined.

 o----------------------------------------------------------------------------o
 */

#include "mad.h"

// --- interface -------------------------------------------------------------o

#define malloc(size_)         mad_malloc  (size_)
#define realloc(ptr_, size_)  mad_realloc (ptr_, size_)
#define free(ptr_)            mad_free    (ptr_)

// local temporary buffer
#define alloc_tmp(type,name,length) mad_alloc_tmp (type,name,length)
#define free_tmp(name)              mad_free_tmp  (name)

// allocator
void*  mad_malloc  (size_t size_);
void*  mad_realloc (void  *ptr_ , size_t size_);
void   mad_free    (void  *ptr_);

// utils
size_t mad_mem_size    (void *ptr_);
size_t mad_mem_cached  (void);
size_t mad_mem_collect (void);

// --- implementation (private) ----------------------------------------------o

#include "mad_mem_impl.h"

// ---------------------------------------------------------------------------o

#endif // MAD_MEM_H
