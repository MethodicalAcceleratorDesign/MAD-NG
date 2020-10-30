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
  - Routines for synchrotron radiation.

  Reference:
  - H. Burkhardt, CERN CLIC-Note-709 of 08-06-2007 for routine InvSynFracInt.

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"

// --- interface --------------------------------------------------------------o

num_t mad_rad_InvSynFracInt (num_t x); // HBU 2007

// ----------------------------------------------------------------------------o

#if 0
// Obsolete code not used, adapted from Placet (courtesy A. Latina)

num_t mad_rad_nrjloss_average (num_t gamma , num_t kick, num_t length);
num_t mad_rad_nrjloss_quantum (num_t gamma , num_t kick, num_t length);
num_t mad_rad_freepath        (num_t betgam, num_t kick, num_t length);
num_t mad_rad_synrad_prob     (num_t betgam, num_t kick);

num_t mad_rad_randexp_seed    (num_t val);
num_t mad_rad_randexp         (num_t mu);
#endif

// ----------------------------------------------------------------------------o

#endif // MAD_RAD_H
