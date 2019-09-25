/*
 o-----------------------------------------------------------------------------o
 |
 | Constants module implementation
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

#include <math.h>
#include <float.h>

#include "mad_cst.h"

// --- implementation ---------------------------------------------------------o

// constants for maths

const num_t mad_cst_EPS      = DBL_EPSILON;
const num_t mad_cst_TINY     = DBL_MIN;
const num_t mad_cst_HUGE     = DBL_MAX;
const num_t mad_cst_INFINITY = INFINITY;

const num_t mad_cst_E        = M_E;
const num_t mad_cst_LOG2E    = M_LOG2E;
const num_t mad_cst_LOG10E   = M_LOG10E;
const num_t mad_cst_LN2      = M_LN2;
const num_t mad_cst_LN10     = M_LN10;
const num_t mad_cst_PI       = M_PI;
const num_t mad_cst_2PI      = M_2PI;
const num_t mad_cst_PI_2     = M_PI_2;
const num_t mad_cst_PI_4     = M_PI_4;
const num_t mad_cst_1_PI     = M_1_PI;
const num_t mad_cst_2_PI     = M_2_PI;
const num_t mad_cst_2_SQRTPI = M_2_SQRTPI;
const num_t mad_cst_SQRT2    = M_SQRT2;
const num_t mad_cst_SQRT1_2  = M_SQRT1_2;
const num_t mad_cst_DEGRAD   = M_DEGRAD;
const num_t mad_cst_RADDEG   = M_RADDEG;

// constants for physics

const num_t mad_cst_MINLEN   = 1e-10;
const num_t mad_cst_MINANG   = 1e-10;
const num_t mad_cst_MINSTR   = 1e-12;

const num_t mad_cst_CLIGHT   = P_CLIGHT;
const num_t mad_cst_MU0      = P_MU0;
const num_t mad_cst_EPSILON0 = P_EPSILON0;
const num_t mad_cst_QELECT   = P_QELECT;
const num_t mad_cst_HBAR     = P_HBAR;
const num_t mad_cst_AMASS    = P_AMASS;
const num_t mad_cst_EMASS    = P_EMASS;
const num_t mad_cst_NMASS    = P_NMASS;
const num_t mad_cst_PMASS    = P_PMASS;
const num_t mad_cst_MUMASS   = P_MUMASS;
const num_t mad_cst_DEUMASS  = P_DEUMASS;
const num_t mad_cst_ERADIUS  = P_ERADIUS;
