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

#define CHK_(  c) ensure((c) == 1, "NLOPT config error '%s': %s.",\
                                   nlopt_result_to_string(c), nlopt_get_errmsg(opt))
#define CHK( s,c) do { CHK_(c); if (a->debug) printf("nlopt: %s set.\n",#s); } while(0)
#define CHKS(s,c) do { CHK_(c); if (a->debug) printf("nlopt: %s set to '%s'.\n"  ,#s,a->s); } while(0)
#define CHKI(s,c) do { CHK_(c); if (a->debug) printf("nlopt: %s set to %d.\n"    ,#s,a->s); } while(0)
#define CHKN(s,c) do { CHK_(c); if (a->debug) printf("nlopt: %s set to % -.6e.\n",#s,a->s); } while(0)

#define INF INFINITY

static inline
num_t min (ssz_t n, const num_t x[n])
{
  assert(x);
  num_t mx = INF;
  for (ssz_t i=0; i<n; i++) if (x[i] < mx) mx = x[i];
  return mx;
}

static inline
num_t max (ssz_t n, const num_t x[n])
{
  assert(x);
  num_t mx = -INF;
  for (ssz_t i=0; i<n; i++) if (x[i] > mx) mx = x[i];
  return mx;
}

static void
dbg_args (nlopt_args_t *a)
{
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wpedantic"
  printf("** algorithm **\n");
  printf("algo        = %d\n"   , a->algo);
  printf("subalgo     = %d\n"   , a->algo);

  printf("** objective function **\n");
  printf("fun         = %p\n"   , (void*)a->fun);
  printf("fval        = %.16e\n", a->fval );
  printf("fmin        = %.16e\n", a->fmin );
  printf("ftol        = %.16e\n", a->ftol );
  printf("frtol       = %.16e\n", a->frtol);
  printf("fdat        = %p\n"   , a->fdat );

  printf("** state variables **\n");
  printf("n           = %d\n"   , a->n    );
  if (a->x)    FOR(i,a->n) printf("x   [%d]    = %.16e\n", i, a->x   [i]);
  if (a->xstp) FOR(i,a->n) printf("xstp[%d]    = %.16e\n", i, a->xstp[i]);
  if (a->xmin) FOR(i,a->n) printf("xmin[%d]    = %.16e\n", i, a->xmin[i]);
  if (a->xmax) FOR(i,a->n) printf("xmax[%d]    = %.16e\n", i, a->xmax[i]);
  if (a->xtol) FOR(i,a->n) printf("xtol[%d]    = %.16e\n", i, a->xtol[i]);
  printf("xrtol       = %.16e\n", a->xrtol);

  printf("** equality constraints **\n");
  printf("p           = %d\n"   , a->p    );
  printf("efun        = %p\n" , (void*)a->efun);
  if (a->etol) FOR(i,a->p) printf("etol[%d]    = %.16e\n", i, a->etol[i]);
  printf("edat        = %p\n" , a->edat   );

  printf("** inequality constraints **\n");
  printf("q           = %d\n"   , a->q    );
  printf("lfun        = %p\n" , (void*)a->lfun);
  if (a->ltol) FOR(i,a->q) printf("ltol[%d]    = %.16e\n", i, a->ltol[i]);
  printf("ldat        = %p\n" , a->ldat   );

  printf("** stop criteria **\n");
  printf("maxcall     = %d\n"   , a->maxcall);
  printf("maxtime     = %.16e\n", a->maxtime);
#pragma GCC diagnostic pop
}

void mad_nlopt (nlopt_args_t *a)
{
  assert(a);

  if (a->debug >= 4) dbg_args(a);

  // create optimizer, set algorithm and problem dimension
  nlopt_opt opt = nlopt_create(a->algo, a->n); a->opt = opt;

  // retrieve algorithm full name (and sanity check)
  a->algonam = nlopt_algorithm_name(a->algo);
  ensure(a->algonam, "invalid nlopt algorithm id=%d.", a->algo);
  CHKS(algonam, 1);

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
  if (a->xtol && max(a->n,a->xtol) > 0) CHK(xtol, nlopt_set_xtol_abs(opt, a->xtol));

  // set variables initial step size
  if (a->xstp && max(a->n,a->xstp) > 0) CHK(xstp, nlopt_set_initial_step(opt, a->xstp));

  // set variables boundary constraints
  if (a->xmin && max(a->n,a->xmin) > -INF) CHK(xmin, nlopt_set_lower_bounds(opt, a->xmin));
  if (a->xmax && min(a->n,a->xmax) <  INF) CHK(xmax, nlopt_set_upper_bounds(opt, a->xmax));

  // check constraints tolerances
  if (a->etol && !(max(a->p,a->etol) > 0) ) a->etol = NULL; else CHK(etol,1);
  if (a->ltol && !(max(a->q,a->ltol) > 0) ) a->ltol = NULL; else CHK(ltol,1);

  // set constraints to satisfy (within tolerances)
  if (a->efun) CHK(efun, nlopt_add_equality_mconstraint  (opt, a->p, a->efun, a->edat, a->etol));
  if (a->lfun) CHK(lfun, nlopt_add_inequality_mconstraint(opt, a->q, a->lfun, a->ldat, a->ltol));

  // set extra stop criteria
  if (a->maxcall > 0) CHKI(maxcall, nlopt_set_maxeval(opt, a->maxcall));
  if (a->maxtime > 0) CHKN(maxtime, nlopt_set_maxtime(opt, a->maxtime));

  // local optimiser
  if (a->algo >= NLOPT_AUGLAG && a->algo <= NLOPT_G_MLSL_LDS) {
    // retrieve sub algorithm full name (and sanity check)
    a->subalgonam = nlopt_algorithm_name(a->subalgo);
    ensure(a->subalgonam, "invalid nlopt local algorithm id=%d", a->subalgo);

    // create local optimizer, set algorithm and problem dimension
    nlopt_opt sopt = nlopt_create(a->subalgo, a->n);

    // set objective function stop value
    if (a->fmin > -INF) CHKN(fmin, nlopt_set_stopval(sopt, a->fmin));

    // set objective function tolerance (value change)
    if (a->ftol > 0   ) CHKN(ftol, nlopt_set_ftol_abs(sopt, a->ftol));

    // set objective function relative tolerance (relative value change)
    if (a->frtol > 0  ) CHKN(frtol, nlopt_set_ftol_rel(sopt, a->frtol));

    // set variables relative tolerance (same for all)
    if (a->xrtol > 0  ) CHKN(xrtol, nlopt_set_xtol_rel(sopt, a->xrtol));

    // set variables tolerances
    if (a->xtol && max(a->n,a->xtol) > 0) CHK(xtol, nlopt_set_xtol_abs(sopt, a->xtol));

    // set variables initial step size
    if (a->xstp && max(a->n,a->xstp) > 0) CHK(xstp, nlopt_set_initial_step(sopt, a->xstp));

    // set extra stop criteria
    if (a->maxcall > 0) CHKI(maxcall, nlopt_set_maxeval(sopt, a->maxcall));
    if (a->maxtime > 0) CHKN(maxtime, nlopt_set_maxtime(sopt, a->maxtime));

    // set local optimizer
    CHK(subalgo, nlopt_set_local_optimizer(opt, sopt));

    // destroy local optimizer
    nlopt_destroy(sopt);
  }

  // seach for minimum
  a->status = nlopt_optimize(opt, a->x, &a->fval);

  // destroy optimizer
  nlopt_destroy(opt); a->opt = NULL;
}

void mad_nlopt_srand (u64_t seed) { nlopt_srand(seed);  }
void mad_nlopt_srand_time (void)  { nlopt_srand_time(); }
