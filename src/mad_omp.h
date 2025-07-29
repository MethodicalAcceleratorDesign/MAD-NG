#ifndef MAD_OMP_H
#define MAD_OMP_H

/*
 o-----------------------------------------------------------------------------o
 |
 | MAD definitions
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
  - separate omp.h inclusion from other includes for C++ due to a bug in omp.h

 o-----------------------------------------------------------------------------o
 */

// --- Open Multi-Processing -------------------------------------------------o

#ifdef _OPENMP
#include <omp.h>
#else
#define __thread
#define omp_get_num_procs()   1
#define omp_get_num_threads() 1
#define omp_get_max_threads() 1
#define omp_get_thread_num()  0
#endif

// ---------------------------------------------------------------------------o

#endif // MAD_OMP_H
