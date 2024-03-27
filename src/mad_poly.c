/*
 o-----------------------------------------------------------------------------o
 |
 | Polygon module implementation
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

#include <assert.h>

#include "mad_poly.h"

// Check if a point p is inside a closed polygon v (vector of points)
//   based on the winding number test for a point in a polygon
// Input:   px,py = a point,
//          vx,vy = vector of n points of a polygon v with v[n-1]=v[0]
// Return:  wn = the winding number (=0 only when p is outside)
// Ref: http://geomalgorithms.com/a03-_inclusion.html

static inline num_t __attribute__((pure))
is_left (num_t px, num_t py, const num_t vx[], const num_t vy[], ssz_t i)
{
  return (vx[i+1]-vx[i]) * (py-vy[i]) - (px-vx[i]) * (vy[i+1]-vy[i]);
}

log_t
mad_pol_inside(num_t px, num_t py, ssz_t n, const num_t vx[n], const num_t vy[n])
{
  assert(vx && vy);
  int wn = 0; // the winding number counter

  for (ssz_t i=0; i < n-1; i++) {
    if (vy[i  ] <= py && vy[i+1] > py && is_left(px, py, vx, vy, i) > 0) ++wn;
    else
    if (vy[i+1] <= py &&                 is_left(px, py, vx, vy, i) < 0) --wn;
  }

  return !wn;
}
