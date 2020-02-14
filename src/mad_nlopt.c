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

#define CHK_(  c) ensure((c) == 1, "unable to configure nlopt");
#define CHK( s,c) do { CHK_(c) if (a->debug) printf("nlopt: %s set\n",#s); } while(0)
#define CHKI(s,c) do { CHK_(c) if (a->debug) printf("nlopt: %s set to %d\n",#s,a->s); } while(0)
#define CHKN(s,c) do { CHK_(c) if (a->debug) printf("nlopt: %s set to % -.6e\n",#s,a->s); } while(0)

#define INF INFINITY

static inline
num_t min (ssz_t n, const num_t x[n])
{
  num_t mx = INF;
  if (x) for (ssz_t i=0; i<n; i++) if (x[i] < mx) mx = x[i];
  return mx;
}

static inline
num_t max (ssz_t n, const num_t x[n])
{
  num_t mx = -INF;
  if (x) for (ssz_t i=0; i<n; i++) if (x[i] > mx) mx = x[i];
  return mx;
}

void mad_nlopt (nlopt_args_t *a)
{
  assert(a);

  // retrieve algorithm full name (and sanity check)
  a->algonam = nlopt_algorithm_name(a->algo);
  ensure(a->algonam, "invalid nlopt algorithm id=%d", a->algo);

  // create optimizer, set algorithm and problem dimension
  nlopt_opt opt = nlopt_create(a->algo, a->n); a->opt = opt;

  // set objective function to minimize
  CHK(fun, nlopt_set_min_objective(opt, a->fun, a->fdat));

  // set objective function stop value
  if (a->fmin > -INF) CHKN(fmin, nlopt_set_stopval(opt, a->fmin));

  // set objective function tolerance (value change)
  if (a->ftol > 0   ) CHKN(ftol, nlopt_set_ftol_abs(opt, a->ftol));

  // set objective function relative tolerance (relative value change)
  if (a->frtol > 0  ) CHKN(frtol, nlopt_set_ftol_rel(opt, a->frtol));

  // set variables relative tolerance (same for all)
  if (a->xrtol > 0  ) CHKN(xrtol, nlopt_set_xtol_rel(opt, a->xrtol));

  // set variables tolerances
  if (max(a->n,a->xtol) > 0) CHK(xtol, nlopt_set_xtol_abs(opt, a->xtol));

  // set variables initial step size
  if (max(a->n,a->xstp) > 0) CHK(xstp, nlopt_set_initial_step(opt, a->xstp));

  // set variables boundary constraints
  if (max(a->n,a->xmin) > -INF) CHK(xmin, nlopt_set_lower_bounds(opt, a->xmin));
  if (min(a->n,a->xmax) <  INF) CHK(xmax, nlopt_set_upper_bounds(opt, a->xmax));

  // check constraints tolerances
  if (! (max(a->p,a->etol) > 0) ) a->etol = NULL; else CHK(etol,1);
  if (! (max(a->q,a->ltol) > 0) ) a->ltol = NULL; else CHK(ltol,1);

  // set constraints to satisfy (within tolerances)
  if (a->efun) CHK(efun, nlopt_add_equality_mconstraint  (opt, a->p, a->efun, a->edat, a->etol));
  if (a->lfun) CHK(lfun, nlopt_add_inequality_mconstraint(opt, a->q, a->lfun, a->ldat, a->ltol));

  // set extra stop criteria
  if (a->maxcall > 0) CHKI(maxcall, nlopt_set_maxeval(opt, a->maxcall));
  if (a->maxtime > 0) CHKN(maxtime, nlopt_set_maxtime(opt, a->maxtime));

  // seach for minimum
  a->status = nlopt_optimize(opt, a->x, &a->fval);

  // destroy optimizer
  nlopt_destroy(opt); a->opt = NULL;
}

void mad_nlopt_srand (u64_t seed) { nlopt_srand(seed);  }
void mad_nlopt_srand_time (void)  { nlopt_srand_time(); }
