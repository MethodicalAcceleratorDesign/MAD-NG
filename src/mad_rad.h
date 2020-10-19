#ifndef MAD_RAD_H
#define MAD_RAD_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Radiation module interface
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
  - Routines for synchrotron radiation from Placet (courtesy to A. Latina)

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"

// --- interface --------------------------------------------------------------o

num_t mad_rad_nrjloss_quantum (num_t particle_nrj, num_t kick, num_t length);
num_t mad_rad_nrjloss_average (num_t particle_nrj, num_t kick, num_t length);
num_t mad_rad_freepath        (num_t particle_nrj, num_t kick, num_t length);
num_t mad_rad_freepath_mean   (num_t particle_nrj, num_t kick, num_t length);
num_t mad_rad_synrad_prob     (num_t particle_nrj, num_t kick);

// ----------------------------------------------------------------------------o

#endif // MAD_RAD_H
