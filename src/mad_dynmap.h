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

typedef struct elem elem_t;
typedef struct mflw mflw_t;

typedef struct part  par_t;
typedef struct damap map_t;

// --- interface --------------------------------------------------------------o

typedef void (trkfun) (elem_t*, mflw_t*, num_t, int);

void mad_trk_strex_drift_r (elem_t *e, mflw_t *m, num_t lw, int i);
void mad_trk_strex_drift_t (elem_t *e, mflw_t *m, num_t lw, int i);
void mad_trk_strex_kick_r  (elem_t *e, mflw_t *m, num_t lw, int i);
void mad_trk_strex_kick_t  (elem_t *e, mflw_t *m, num_t lw, int i);

void mad_trk_curex_drift_r (elem_t *e, mflw_t *m, num_t lw, int i);
void mad_trk_curex_drift_t (elem_t *e, mflw_t *m, num_t lw, int i);
void mad_trk_curex_kick_r  (elem_t *e, mflw_t *m, num_t lw, int i);
void mad_trk_curex_kick_t  (elem_t *e, mflw_t *m, num_t lw, int i);

void mad_trk_slice_r       (elem_t *e, mflw_t *m, num_t lw, trkfun *dft, trkfun *kck);
void mad_trk_slice_t       (elem_t *e, mflw_t *m, num_t lw, trkfun *dft, trkfun *kck);

// --- end --------------------------------------------------------------------o

#endif // MAD_DYNMAP_H
