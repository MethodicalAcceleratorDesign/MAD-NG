/*
 o-----------------------------------------------------------------------------o
 |
 | Non Linear Optimization module implementation
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
#include <stdio.h>
#include <assert.h>
#include "mad_log.h"
#include "mad_nlopt.h"

void mad_nlopt (nlopt_args_t *a)
{
  assert(a);
  ensure(          a->n > 0 && a->fn && a->x    , "invalid objective function");
  ensure(!a->p || (a->p > 0 && a->eq && a->etol), "invalid equality constraints");
  ensure(!a->q || (a->q > 0 && a->le && a->ltol), "invalid inequality constraints");

  // create optimizer, set algorithm and problem dimension
  nlopt_opt opt = nlopt_create(a->algo, a->n);

  // set objective function to minimize/maximize
  if (a->fdir > 0) nlopt_set_max_objective(opt, a->fn, NULL);
  else             nlopt_set_min_objective(opt, a->fn, NULL);

  // set objective function tolerances
  nlopt_set_ftol_abs    (opt, a->ftol);
  // set objective function stop value
  nlopt_set_stopval     (opt, a->fstop);
  // set variables intial steps
  nlopt_set_initial_step(opt, a->xstp);
  // set variables tolerances
  nlopt_set_xtol_abs    (opt, a->xtol);
  // set variables relative tolerance
  nlopt_set_xtol_rel    (opt, a->rtol);
  // set variables boundary constraints
  nlopt_set_lower_bounds(opt, a->xmin);
  nlopt_set_upper_bounds(opt, a->xmax);

  // set constraint functions to satisfy withing tolerances
  if (a->p) nlopt_add_equality_mconstraint  (opt, a->p, a->eq, NULL, a->etol);
  if (a->q) nlopt_add_inequality_mconstraint(opt, a->q, a->le, NULL, a->ltol);

  // set stopping criteria
  nlopt_set_maxeval(opt, a->maxcall);
  nlopt_set_maxtime(opt, a->maxtime);

  // seach for optimum
  a->status = nlopt_optimize(opt, a->x, &a->fval);

  // destroy optimizer
  nlopt_destroy(opt);
}

void mad_nlopt_srand (u64_t seed)
{
  nlopt_srand(seed);
  nlopt_srand_time();
}
