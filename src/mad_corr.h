#ifndef MAD_CORR_H
#define MAD_CORR_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Orbit Correction module interface
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
  - todo

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"

void mad_corr_mic(num_t A[], num_t B[], num_t X[], ssz_t m, ssz_t n,
                  num_t tol, int ncorr);

void mad_corr_lsq(void);
void mad_corr_svd(void);


// ----------------------------------------------------------------------------o

#endif // MAD_CORR_H
