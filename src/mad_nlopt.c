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

#define DBG(s) if (a->debug) puts("nlopt: " #s " set")
#define CHK(a) ensure(a == 1, "unable to configure nlopt")

static inline
num_t min(ssz_t n, const num_t x[n])
{
  num_t mx = INFINITY;
  for (ssz_t i=0; i<n; i++) if (x[i] < mx) mx = x[i];
  return mx;
}

static inline
num_t max(int n, const num_t x[n])
{
  num_t mx = -INFINITY;
  for (ssz_t i=0; i<n; i++) if (x[i] > mx) mx = x[i];
  return mx;
}

void mad_nlopt (nlopt_args_t *a)
{
  assert(a);

  // retrieve algorithm full name (and sanity check)
  a->algonam = nlopt_algorithm_name(a->algo);
  ensure(a->algonam, "invalid nlopt algorithm id=%d", a->algo);

  // create optimizer, set algorithm and problem dimension
  a->opt = nlopt_create(a->algo, a->n);

  // set objective function to minimize
  CHK(nlopt_set_min_objective(a->opt, a->fun, NULL)); DBG(fun);

  // set objective function stop value
  if (a->fmin > -INFINITY) { CHK(nlopt_set_stopval(a->opt, a->fmin)); DBG(fmin); }

  // set objective function tolerance (value change)
  if (a->ftol > 0) { CHK(nlopt_set_ftol_abs(a->opt, a->ftol)); DBG(ftol); }

  // set objective function relative tolerance (relative value change)
  if (a->frtol > 0) { CHK(nlopt_set_ftol_rel(a->opt, a->frtol)); DBG(frtol); }

  // set variables tolerances
  if (max(a->n,a->xtol) > 0) { CHK(nlopt_set_xtol_abs(a->opt, a->xtol)); DBG(xtol); }

  // set variables relative tolerance (same for all)
  if (a->xrtol > 0) { CHK(nlopt_set_xtol_rel(a->opt, a->xrtol)); DBG(xrtol); }

  // set variables initial step size
  if (min(a->n,a->dx) > 0) { CHK(nlopt_set_initial_step(a->opt, a->dx)); DBG(xstp); }

  // set variables boundary constraints
  if (max(a->n,a->xmin) > -INFINITY) { CHK(nlopt_set_lower_bounds(a->opt, a->xmin)); DBG(xmin); }
  if (min(a->n,a->xmax) <  INFINITY) { CHK(nlopt_set_upper_bounds(a->opt, a->xmax)); DBG(xmax); }

  // check constraints tolerances
  if (! (max(a->p,a->etol) > 0) ) a->etol = NULL; else DBG(etol);
  if (! (max(a->q,a->ltol) > 0) ) a->ltol = NULL; else DBG(ltol);

  // set constraints to satisfy (within tolerances)
  if (a->efun) { CHK(nlopt_add_equality_mconstraint  (a->opt, a->p, a->efun, NULL, a->etol)); DBG(efun); }
  if (a->lfun) { CHK(nlopt_add_inequality_mconstraint(a->opt, a->q, a->lfun, NULL, a->ltol)); DBG(lfun); }

  // set extra stopping criteria
  if (a->maxcall > 0) { CHK(nlopt_set_maxeval(a->opt, a->maxcall)); DBG(maxcall); }
  if (a->maxtime > 0) { CHK(nlopt_set_maxtime(a->opt, a->maxtime)); DBG(maxtime); }

  // seach for minimum
  a->status = nlopt_optimize(a->opt, a->x, &a->fval);

  // destroy optimizer
  nlopt_destroy(a->opt); a->opt = NULL;
}

void mad_nlopt_srand (u64_t seed)
{
  nlopt_srand(seed);
  nlopt_srand_time();
}
