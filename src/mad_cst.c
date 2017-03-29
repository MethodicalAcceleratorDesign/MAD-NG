/*
 o-----------------------------------------------------------------------------o
 |
 | Constants module implementation
 |
 | Methodical Accelerator Design - Copyright CERN 2016+
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

#include "mad_cst.h"

// --- implementation ---------------------------------------------------------o

const num_t mad_cst_minlen   = 1e-12     ; // m   minimum length tolerance
const num_t mad_cst_minang   = 1e-12     ; // rad minimum angle  tolerance

const num_t mad_cst_E        = M_E       ; // e
const num_t mad_cst_LOG2E    = M_LOG2E   ; // log_2 e
const num_t mad_cst_LOG10E   = M_LOG10E  ; // log_10 e
const num_t mad_cst_LN2      = M_LN2     ; // log_e 2
const num_t mad_cst_LN10     = M_LN10    ; // log_e 10
const num_t mad_cst_PI       = M_PI      ; // pi
const num_t mad_cst_PI_2     = M_PI_2    ; // pi/2
const num_t mad_cst_PI_4     = M_PI_4    ; // pi/4
const num_t mad_cst_1_PI     = M_1_PI    ; // 1/pi
const num_t mad_cst_2_PI     = M_2_PI    ; // 2/pi
const num_t mad_cst_2_SQRTPI = M_2_SQRTPI; // 2/sqrt(pi)
const num_t mad_cst_SQRT2    = M_SQRT2   ; // sqrt(2)
const num_t mad_cst_1_SQRT2  = M_1_SQRT2 ; // 1/sqrt(2)

const num_t mad_cst_CLIGHT   = P_CLIGHT  ; // m/s
const num_t mad_cst_MU0      = P_MU0     ; // T.m/A (or N/A^2)
const num_t mad_cst_EPSILON0 = P_EPSILON0; // F/m
const num_t mad_cst_QELECT   = P_QELECT  ; // C
const num_t mad_cst_HBAR     = P_HBAR    ; // GeV.s
const num_t mad_cst_EMASS    = P_EMASS   ; // GeV
const num_t mad_cst_PMASS    = P_PMASS   ; // GeV
const num_t mad_cst_NMASS    = P_NMASS   ; // GeV
const num_t mad_cst_MUMASS   = P_MUMASS  ; // GeV
const num_t mad_cst_DEUMASS  = P_DEUMASS ; // GeV
const num_t mad_cst_ERADIUS  = P_ERADIUS ; // m (elec. mag. radius)
