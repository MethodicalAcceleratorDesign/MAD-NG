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
#include <nlopt.h>
#include "mad_nlopt.h"

// FIRST EXAMPLE

typedef struct {
    double a, b;
} my_constraint_data;

static double myfunc(unsigned n, const double *x, double *grad, void *data)
{
  assert(n == 2);
  assert(x && data);
  ++*(int*)data;

  if (grad) {
      grad[0] = 0.0;
      grad[1] = 0.5 / sqrt(x[1]);
  }

  return sqrt(x[1]);
}

static double myconstraint(unsigned n, const double *x, double *grad, void *data)
{
  assert(n == 2);
  assert(x && data);
  my_constraint_data *d = (my_constraint_data *) data;

  double a = d->a, b = d->b;

  if (grad) {
      grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
      grad[1] = -1.0;
  }

  return ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
 }

void mad_nlopt (void)
{
  nlopt_opt opt = nlopt_create(NLOPT_LD_MMA, 2); /* algorithm and dimensionality */

  double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
  nlopt_set_lower_bounds(opt, lb);

  int count = 0;
  nlopt_set_min_objective(opt, myfunc, &count);

  my_constraint_data data[2] = { {2,0}, {-1,1} };
  nlopt_add_inequality_constraint(opt, myconstraint, &data[0], 1e-12);
  nlopt_add_inequality_constraint(opt, myconstraint, &data[1], 1e-12);

  nlopt_set_xtol_rel(opt, 1e-12);

  double x[2] = { 1.234, 5.678 };  /* `*`some` `initial` `guess`*` */
  double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */

  if (nlopt_optimize(opt, x, &minf) < 0) {
      printf("nlopt failed!\n");
  }
  else {
      printf("found minimum at f(%g,%g) = %0.10g in %d iterations\n", x[0], x[1], minf, count);
  }

  nlopt_destroy(opt);
}
