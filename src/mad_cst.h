#ifndef MAD_CST_H
#define MAD_CST_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Constants module interface
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
  - Provide a unique place to define constants.

 o-----------------------------------------------------------------------------o
 */

#include "mad_defs.h"

// --- interface --------------------------------------------------------------o

// constants for maths

extern const num_t mad_cst_EPS     ;  // minimum value such that 1+eps > 1
extern const num_t mad_cst_TINY    ;  // minimum value above 0
extern const num_t mad_cst_HUGE    ;  // maximum value below inf
extern const num_t mad_cst_INFINITY;  // maximum representable value inf

extern const num_t mad_cst_E       ;  // e
extern const num_t mad_cst_LOG2E   ;  // log_2 e
extern const num_t mad_cst_LOG10E  ;  // log_10 e
extern const num_t mad_cst_LN2     ;  // log_e 2
extern const num_t mad_cst_LN10    ;  // log_e 10
extern const num_t mad_cst_PI      ;  // pi
extern const num_t mad_cst_2PI     ;  // 2*pi
extern const num_t mad_cst_PI_2    ;  // pi/2
extern const num_t mad_cst_PI_4    ;  // pi/4
extern const num_t mad_cst_SQRTPI  ;  // sqrt(pi)
extern const num_t mad_cst_1_PI    ;  // 1/pi
extern const num_t mad_cst_1_SQRTPI;  // 1/sqrt(pi)
extern const num_t mad_cst_2_PI    ;  // 2/pi
extern const num_t mad_cst_2_SQRTPI;  // 2/sqrt(pi)
extern const num_t mad_cst_SQRT2   ;  // sqrt(2)
extern const num_t mad_cst_SQRT1_2 ;  // 1/sqrt(2), sqrt(1/2)
extern const num_t mad_cst_SQRT3   ;  // sqrt(3)
extern const num_t mad_cst_SQRT1_3 ;  // 1/sqrt(3), sqrt(1/3)

extern const num_t mad_cst_DEGRAD  ;  // degree to radian factor
extern const num_t mad_cst_RADDEG  ;  // radian to degree factor

// constants for physics

extern const num_t mad_cst_MINLEN  ;  // [m]   Minimum tolerance on lengths
extern const num_t mad_cst_MINANG  ;  // [rad] Minimum tolerance on angles
extern const num_t mad_cst_MINSTR  ;  // [1/m] Minimum tolerance on strengths

extern const num_t mad_cst_CLIGHT  ;  // [m/s]
extern const num_t mad_cst_MU0     ;  // [T.m/A] or [N/A^2] or [V.s/(A.m)]
extern const num_t mad_cst_EPSILON0;  // [F/m]
extern const num_t mad_cst_QELECT  ;  // [C]
extern const num_t mad_cst_HBAR    ;  // [GeV.s]
extern const num_t mad_cst_AMASS   ;  // [GeV]
extern const num_t mad_cst_EMASS   ;  // [GeV]
extern const num_t mad_cst_NMASS   ;  // [GeV]
extern const num_t mad_cst_PMASS   ;  // [GeV]
extern const num_t mad_cst_MUMASS  ;  // [GeV]
extern const num_t mad_cst_DEUMASS ;  // [GeV]
extern const num_t mad_cst_ERADIUS ;  // [m]

// --- math constants ---------------------------------------------------------o

#ifndef M_E         // standard constants
#define M_E         2.71828182845904523536028747135266250   // e
#define M_LOG2E     1.44269504088896340735992468100189214   // log_2 e
#define M_LOG10E    0.434294481903251827651128918916605082  // log_10 e
#define M_LN2       0.693147180559945309417232121458176568  // log_e 2
#define M_LN10      2.30258509299404568401799145468436421   // log_e 10
#define M_PI        3.14159265358979323846264338327950288   // pi
#define M_PI_2      1.57079632679489661923132169163975144   // pi/2
#define M_PI_4      0.785398163397448309615660845819875721  // pi/4
#define M_1_PI      0.318309886183790671537767526745028724  // 1/pi
#define M_2_PI      0.636619772367581343075535053490057448  // 2/pi
#define M_2_SQRTPI  1.12837916709551257389615890312154517   // 2/sqrt(pi)
#define M_SQRT2     1.41421356237309504880168872420969808   // sqrt(2)
#define M_SQRT1_2   0.707106781186547524400844362104849039  // 1/sqrt(2)
#endif
                    // extra constants
#define M_2PI       6.28318530717958647692528676655900577   // 2*pi
#define M_SQRTPI    1.77245385090551602729816748334114518   // sqrt(pi)
#define M_1_SQRTPI  0.564189583547756286948079451560772586  // 1/sqrt(pi)
#define M_SQRT3     1.73205080756887729352744634150587237   // sqrt(3)
#define M_SQRT1_3   0.577350269189625764509148780501957456  // 1/sqrt(3)

                    // constants for conversion
#define M_RADDEG    57.2957795130823208767981548141051703   // 180/pi
#define M_DEGRAD    0.0174532925199432957692369076848861271 // pi/180

// --- physics constants ------------------------------------------------------o

// https://physics.nist.gov/cuu/pdf/wall_2018.pdf
#ifndef P_CLIGHT                           // Source: CODATA 2018
#define P_CLIGHT     299792458.0           // [m/s]   Speed of light in vacuum
#define P_CLIGHT2   (P_CLIGHT*P_CLIGHT)    //         c^2
#define P_MU0       (4e-7*M_PI)            // [T.m/A] Permeability of vacuum
#define P_EPSILON0  (1/(P_MU0*P_CLIGHT2))  // [F/m]   Permittivity of vacuum
#define P_QELECT     1.602176634e-19       // [C]     Elementary electric charge
#define P_HBAR      (6.582119569e-16*1e-9) // [GeV.s] Reduced Plack's constant
#define P_AMASS     (931.49410242   *1e-3) // [GeV]   Unified atomic mass
#define P_EMASS     (0.51099895000  *1e-3) // [GeV]   Electron energy-mass
#define P_PMASS     (938.27208816   *1e-3) // [GeV]   Proton energy-mass
#define P_NMASS     (939.56542052   *1e-3) // [GeV]   Neutron energy-mass
#define P_MUMASS    (105.6583755    *1e-3) // [GeV]   Muon energy-mass
#define P_DEUMASS   (1875.61294257  *1e-3) // [GeV]   Deuteron energy-mass
#define P_ERADIUS    2.8179403262e-15      // [m]     Classical electron radius
#define P_ALPHAEM    7.2973525693e-3       //         Fine-structure constant
#endif

// https://physics.nist.gov/cuu/pdf/wall_2014.pdf
#ifndef P_CLIGHT                           // Source: CODATA 2014
#define P_CLIGHT    299792458.0            // [m/s]   Speed of light in vacuum
#define P_CLIGHT2   (P_CLIGHT*P_CLIGHT)    //         c^2
#define P_MU0       (4e-7*M_PI)            // [T.m/A] Permeability of vacuum
#define P_EPSILON0  (1/(P_MU0*P_CLIGHT2))  // [F/m]   Permittivity of vacuum
#define P_QELECT    1.6021766208e-19       // [C]     Elementary electric charge
#define P_HBAR      6.582119514e-25        // [GeV.s] Reduced Plack's constant
#define P_AMASS     0.9314940954           // [GeV]   Unified atomic mass
#define P_EMASS     5.109989461e-4         // [GeV]   Electron mass
#define P_PMASS     0.9382720813           // [GeV]   Proton mass
#define P_NMASS     0.9395654133           // [GeV]   Neutron mass
#define P_MUMASS    0.1056583745           // [GeV]   Muon mass
#define P_DEUMASS   1.875612928            // [GeV]   Deuteron mass
#define P_ERADIUS   2.8179403227e-15       // [m]     Classical electron radius
#define P_ALPHAEM   7.2973525693e-3        //         Fine-structure constant
#endif

// ----------------------------------------------------------------------------o

#endif // MAD_CST_H
