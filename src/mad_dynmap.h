#ifndef MAD_DYNMAP_H
#define MAD_DYNMAP_H

/*
 o-----------------------------------------------------------------------------o
 |
 | C interface to dynamic maps
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
  - Provide interface to dynamic maps in C.

 o-----------------------------------------------------------------------------o
 */

#include "mad_def.h"

typedef union mflw_ mflw_t;
typedef void (trkfun) (mflw_t*, num_t, int);

// --- interface --------------------------------------------------------------o

// -- track one slice
void mad_trk_slice_dkd (mflw_t *m, num_t lw, trkfun *dft, trkfun *kck, int ord);// yoshida
void mad_trk_slice_tkt (mflw_t *m, num_t lw, trkfun *dft, trkfun *kck, int ord);// yoshida
void mad_trk_slice_kmk (mflw_t *m, num_t lw, trkfun *dft, trkfun *kck, int ord);// boole
void mad_trk_slice_tpt (mflw_t *m, num_t lw, trkfun *dft, trkfun *kck, int knd);// teapot
void mad_trk_slice_one (mflw_t *m, num_t lw, trkfun *dft_or_kck);               // single

// -- tilt & misalignment
void mad_trk_tilt_r         (mflw_t *m, num_t lw);
void mad_trk_tilt_t         (mflw_t *m, num_t lw);
void mad_trk_tilt_p         (mflw_t *m, num_t lw);

void mad_trk_misalign_r     (mflw_t *m, num_t lw);
void mad_trk_misalign_t     (mflw_t *m, num_t lw);
void mad_trk_misalign_p     (mflw_t *m, num_t lw);

// -- fringe maps
void mad_trk_strex_fringe_r (mflw_t *m, num_t lw);
void mad_trk_curex_fringe_r (mflw_t *m, num_t lw);
void mad_trk_rfcav_fringe_r (mflw_t *m, num_t lw);

void mad_trk_strex_fringe_t (mflw_t *m, num_t lw);
void mad_trk_curex_fringe_t (mflw_t *m, num_t lw);
void mad_trk_rfcav_fringe_t (mflw_t *m, num_t lw);

void mad_trk_strex_fringe_p (mflw_t *m, num_t lw);
void mad_trk_curex_fringe_p (mflw_t *m, num_t lw);
void mad_trk_rfcav_fringe_p (mflw_t *m, num_t lw);

// -- patches
void mad_trk_xrotation_r    (mflw_t *m, num_t lw, int _);
void mad_trk_yrotation_r    (mflw_t *m, num_t lw, int _);
void mad_trk_srotation_r    (mflw_t *m, num_t lw, int _);
void mad_trk_translate_r    (mflw_t *m, num_t lw, int _);
void mad_trk_changeref_r    (mflw_t *m, num_t lw, int _);

void mad_trk_xrotation_t    (mflw_t *m, num_t lw, int _);
void mad_trk_yrotation_t    (mflw_t *m, num_t lw, int _);
void mad_trk_srotation_t    (mflw_t *m, num_t lw, int _);
void mad_trk_translate_t    (mflw_t *m, num_t lw, int _);
void mad_trk_changeref_t    (mflw_t *m, num_t lw, int _);

void mad_trk_xrotation_p    (mflw_t *m, num_t lw, int _);
void mad_trk_yrotation_p    (mflw_t *m, num_t lw, int _);
void mad_trk_srotation_p    (mflw_t *m, num_t lw, int _);
void mad_trk_translate_p    (mflw_t *m, num_t lw, int _);
void mad_trk_changeref_p    (mflw_t *m, num_t lw, int _);

// -- DKD maps
void mad_trk_strex_drift_r  (mflw_t *m, num_t lw, int _);
void mad_trk_strex_kick_r   (mflw_t *m, num_t lw, int _);
void mad_trk_strex_kickhs_r (mflw_t *m, num_t lw, int _);
void mad_trk_curex_drift_r  (mflw_t *m, num_t lw, int _);
void mad_trk_curex_kick_r   (mflw_t *m, num_t lw, int _);

void mad_trk_strex_drift_t  (mflw_t *m, num_t lw, int _);
void mad_trk_strex_kick_t   (mflw_t *m, num_t lw, int _);
void mad_trk_strex_kickhs_t (mflw_t *m, num_t lw, int _);
void mad_trk_curex_drift_t  (mflw_t *m, num_t lw, int _);
void mad_trk_curex_kick_t   (mflw_t *m, num_t lw, int _);

void mad_trk_strex_drift_p  (mflw_t *m, num_t lw, int _);
void mad_trk_strex_kick_p   (mflw_t *m, num_t lw, int _);
void mad_trk_strex_kickhs_p (mflw_t *m, num_t lw, int _);
void mad_trk_curex_drift_p  (mflw_t *m, num_t lw, int _);
void mad_trk_curex_kick_p   (mflw_t *m, num_t lw, int _);

// -- TKT maps
void mad_trk_sbend_thick_r  (mflw_t *m, num_t lw, int _);
void mad_trk_sbend_kick_r   (mflw_t *m, num_t lw, int _);
void mad_trk_rbend_thick_r  (mflw_t *m, num_t lw, int _);
void mad_trk_rbend_kick_r   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thick_r   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thicks_r  (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thickh_r  (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kick_r    (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kicks_r   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kickh_r   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kick__r   (mflw_t *m, num_t lw, int is);
void mad_trk_quad_kicks__r  (mflw_t *m, num_t lw, int is);
void mad_trk_quad_kickh__r  (mflw_t *m, num_t lw, int is);

void mad_trk_sbend_thick_t  (mflw_t *m, num_t lw, int _);
void mad_trk_sbend_kick_t   (mflw_t *m, num_t lw, int _);
void mad_trk_rbend_thick_t  (mflw_t *m, num_t lw, int _);
void mad_trk_rbend_kick_t   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thick_t   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thicks_t  (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thickh_t  (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kick_t    (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kicks_t   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kickh_t   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kick__t   (mflw_t *m, num_t lw, int is);
void mad_trk_quad_kicks__t  (mflw_t *m, num_t lw, int is);
void mad_trk_quad_kickh__t  (mflw_t *m, num_t lw, int is);

void mad_trk_sbend_thick_p  (mflw_t *m, num_t lw, int _);
void mad_trk_sbend_kick_p   (mflw_t *m, num_t lw, int _);
void mad_trk_rbend_thick_p  (mflw_t *m, num_t lw, int _);
void mad_trk_rbend_kick_p   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thick_p   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thicks_p  (mflw_t *m, num_t lw, int _);
void mad_trk_quad_thickh_p  (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kick_p    (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kicks_p   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kickh_p   (mflw_t *m, num_t lw, int _);
void mad_trk_quad_kick__p   (mflw_t *m, num_t lw, int is);
void mad_trk_quad_kicks__p  (mflw_t *m, num_t lw, int is);
void mad_trk_quad_kickh__p  (mflw_t *m, num_t lw, int is);

// -- other maps
void mad_trk_solen_thick_r  (mflw_t *m, num_t lw, int _);
void mad_trk_esept_thick_r  (mflw_t *m, num_t lw, int _);
void mad_trk_rfcav_kick_r   (mflw_t *m, num_t lw, int _);
void mad_trk_rfcav_kickn_r  (mflw_t *m, num_t lw, int _);

void mad_trk_solen_thick_t  (mflw_t *m, num_t lw, int _);
void mad_trk_esept_thick_t  (mflw_t *m, num_t lw, int _);
void mad_trk_rfcav_kick_t   (mflw_t *m, num_t lw, int _);
void mad_trk_rfcav_kickn_t  (mflw_t *m, num_t lw, int _);

void mad_trk_solen_thick_p  (mflw_t *m, num_t lw, int _);
void mad_trk_esept_thick_p  (mflw_t *m, num_t lw, int _);
void mad_trk_rfcav_kick_p   (mflw_t *m, num_t lw, int _);
void mad_trk_rfcav_kickn_p  (mflw_t *m, num_t lw, int _);

// -- do nothing
void mad_trk_fnil           (mflw_t *m, num_t lw, int _);

// -- benchmark
void mad_trk_spdtest (int n, int k);
void mad_trk_cpptest (void);

// -- cross-check with LJFFI --

extern const size_t mad_cflw_rsize;
extern const size_t mad_cflw_tsize;
extern const size_t mad_cflw_psize;
extern const size_t mad_cflw_size;

// --- end --------------------------------------------------------------------o

#endif // MAD_DYNMAP_H
