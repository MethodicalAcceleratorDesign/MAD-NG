/*
 o-----------------------------------------------------------------------------o
 |
 | Dynamic maps implementation in mixed C/C++
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
*/

extern "C" {
#include "mad_dynmap.h"
}

#include "mad_tpsa.hpp"

// --- implementation ---------------------------------------------------------o

void
mad_trk_strex_drift (void)
{
  mad::tpsa a(mad_tpsa_newd(mad_desc_curr, mad_tpsa_default));
}

void
mad_trk_strex_kick (void)
{

}

// ----------------------------------------------------------------------------o
