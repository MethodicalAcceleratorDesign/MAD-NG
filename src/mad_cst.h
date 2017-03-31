#ifndef MAD_CST_H
#define MAD_CST_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Constants module interface
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

  Purpose:
  - Provide a unique place to define constants.

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"

// --- Math constants ---------------------------------------------------------o

#ifndef M_E
#define M_E         2.7182818284590452353602874713526625  // e
#define M_LOG2E     1.4426950408889634073599246810018921  // log_2 e
#define M_LOG10E    0.4342944819032518276511289189166051  // log_10 e
#define M_LN2       0.6931471805599453094172321214581766  // log_e 2
#define M_LN10      2.3025850929940456840179914546843642  // log_e 10
#define M_PI        3.1415926535897932384626433832795029  // pi
#define M_PI_2      1.5707963267948966192313216916397514  // pi/2
#define M_PI_4      0.7853981633974483096156608458198757  // pi/4
#define M_1_PI      0.3183098861837906715377675267450287  // 1/pi
#define M_2_PI      0.6366197723675813430755350534900574  // 2/pi
#define M_2_SQRTPI  1.1283791670955125738961589031215452  // 2/sqrt(pi)
#define M_SQRT2     1.4142135623730950488016887242096981  // sqrt(2)
#define M_1_SQRT2   0.7071067811865475244008443621048490  // 1/sqrt(2)
#endif

// --- Physics constants ------------------------------------------------------o

#ifndef P_CLIGHT                           // Source: CODATA 2014
#define P_CLIGHT    299792458.0            // [m/s]   Speed of light in vacuum
#define P_CLIGHT2   (P_CLIGHT*P_CLIGHT)    //         c^2
#define P_MU0       (4e-7*M_PI)            // [T.m/A] Permeability of vacuum
#define P_EPSILON0  (1/(P_MU0*P_CLIGHT2))  // [F/m]   Permittivity of vacuum
#define P_QELECT    1.602176620898e-19     // [C]     Elementary electric charge
#define P_HBAR      6.58211951440e-25      // [GeV.s] Reduced Plack's constant
#define P_EMASS     5.10998946131e-4       // [GeV]   Electron mass
#define P_PMASS     0.938272081358         // [GeV]   Proton mass
#define P_NMASS     0.939565413358         // [GeV]   Neutron mass
#define P_MUMASS    0.105658374524         // [GeV]   Muon mass
#define P_DEUMASS   1.87561292812          // [GeV]   Deuteron mass
#define P_ERADIUS   2.817940322719e-15     // [m]     Classical electron radius
#endif

// --- interface --------------------------------------------------------------o

extern const num_t mad_cst_minlen  ;  // [m]   Minimum length tolerance
extern const num_t mad_cst_minang  ;  // [rad] Minimum angle  tolerance

extern const num_t mad_cst_E       ;  // e
extern const num_t mad_cst_LOG2E   ;  // log_2 e
extern const num_t mad_cst_LOG10E  ;  // log_10 e
extern const num_t mad_cst_LN2     ;  // log_e 2
extern const num_t mad_cst_LN10    ;  // log_e 10
extern const num_t mad_cst_PI      ;  // pi
extern const num_t mad_cst_PI_2    ;  // pi/2
extern const num_t mad_cst_PI_4    ;  // pi/4
extern const num_t mad_cst_1_PI    ;  // 1/pi
extern const num_t mad_cst_2_PI    ;  // 2/pi
extern const num_t mad_cst_2_SQRTPI;  // 2/sqrt(pi)
extern const num_t mad_cst_SQRT2   ;  // sqrt(2)
extern const num_t mad_cst_1_SQRT2 ;  // 1/sqrt(2)

extern const num_t mad_cst_CLIGHT  ;  // [m/s]
extern const num_t mad_cst_MU0     ;  // [T.m/A] or [N/A^2] or [V.s/(A.m)]
extern const num_t mad_cst_EPSILON0;  // [F/m]
extern const num_t mad_cst_QELECT  ;  // [C]
extern const num_t mad_cst_HBAR    ;  // [GeV.s]
extern const num_t mad_cst_EMASS   ;  // [GeV]
extern const num_t mad_cst_PMASS   ;  // [GeV]
extern const num_t mad_cst_NMASS   ;  // [GeV]
extern const num_t mad_cst_MUMASS  ;  // [GeV]
extern const num_t mad_cst_DEUMASS ;  // [GeV]
extern const num_t mad_cst_ERADIUS ;  // [m]

// ----------------------------------------------------------------------------o

#endif // MAD_CST_H