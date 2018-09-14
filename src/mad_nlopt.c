/*
 o-----------------------------------------------------------------------------o
 |
 | Non Linear Optimization module implementation
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

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include "mad_log.h"
#include "mad_nlopt.h"

void mad_nlopt (nlopt_args_t *a)
{
  assert(a);
  ensure(          a->n > 0 && a->f    && a->x       , "invalid objective function");
  ensure(!a->p || (a->p > 0 && a->c_eq && a->ctol_eq), "invalid equality constraints");
  ensure(!a->q || (a->q > 0 && a->c_le && a->ctol_le), "invalid inequality constraints");

  // create optimizer, set algorithm and problem dimension
  nlopt_opt opt = nlopt_create(a->algorithm, a->n);

  // set objective function to minimize
  nlopt_set_min_objective(opt, a->f, a->fdat);

  // set contraint functions to satisfy withing tolerances
  if (a->p)
    nlopt_add_equality_mconstraint  (opt, a->p, a->c_eq, a->cdat_eq, a->ctol_eq);
  if (a->q)
    nlopt_add_inequality_mconstraint(opt, a->q, a->c_le, a->cdat_le, a->ctol_le);

  // set variables boundary constraints
  if (a->xmin) nlopt_set_lower_bounds(opt, a->xmin);
  if (a->xmax) nlopt_set_upper_bounds(opt, a->xmax);

  // set variables intial steps
  if (a->xstep) nlopt_set_initial_step(opt, a->xstep);

  // set variables tolerances
  if (a->xtol_rel) nlopt_set_xtol_rel(opt, *a->xtol_rel);
  if (a->xtol_abs) nlopt_set_xtol_abs(opt,  a->xtol_abs);

  // set objective function tolerances
  if (a->ftol_rel) nlopt_set_ftol_rel(opt, *a->ftol_rel);
  if (a->ftol_abs) nlopt_set_ftol_abs(opt, *a->ftol_abs);

  // set abnormal stopping criteria
  if (a->maxcall) nlopt_set_maxeval(opt, *a->maxcall);
  if (a->maxtime) nlopt_set_maxtime(opt, *a->maxtime);

  a->status = nlopt_optimize(opt, a->x, &a->fmin);

  nlopt_destroy(opt);
}
