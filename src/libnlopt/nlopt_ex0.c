/*
compile & run:
gcc -std=c99 -Wall -W -pedantic -O3 nlopt_ex0.c \
    -o nlopt_ex0 -I../../lib/nlopt-2.6.1/src/api/ ../../bin/macosx/libnlopt.a
./nlopt_ex0

output:
NLOPT_LD_MMA:
found minimum after 11 evaluations
found minimum at f(0.333333,0.296296) = 0.5443310474
NLOPT_LN_COBYLA:
found minimum after 31 evaluations
found minimum at f(0.333329,0.2962) = 0.5442423017
*/

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <nlopt.h>

typedef struct {
  double a, b;
} my_constraint_data;

static int count = 0;

double myfunc(unsigned n, const double *x, double *grad, void *data)
{
  assert(n == 2); assert(data == NULL);
  ++count;
  if (grad) {
    grad[0] = 0.0;
    grad[1] = 0.5 / sqrt(x[1]);
  }
  return sqrt(x[1]);
}

void myconstraints(unsigned m, double *r, unsigned n, const double *x, double *grad, void *data)
{
  assert(m == 2); assert(n == 2); assert(data);
  my_constraint_data *d = (my_constraint_data*) data;

  { double a = d[0].a, b = d[0].b;
    if (grad) {
      grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
      grad[1] = -1.0;
    }
    r[0] = ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
  }
  { double a = d[1].a, b = d[1].b;
    if (grad) {
      grad[2] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
      grad[3] = -1.0;
    }
    r[1] = ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
  }
}

int main(void)
{
  double lb[2] = { -HUGE_VAL, 0 }; /* lower bounds */
  double tol[2] = { 1e-8, 1e-8 }; /* lower bounds */
  my_constraint_data data[2] = { {2,0}, {-1,1} };
  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, 2); // NLOPT_LN_COBYLA, NLOPT_LD_MMA

  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_min_objective(opt, myfunc, NULL);
  nlopt_add_inequality_mconstraint(opt, 2, myconstraints, data, tol);

  nlopt_set_xtol_rel(opt, 1e-4);

  double x[2] = { 1.234, 5.678 };  /* `*`some` `initial` `guess`*` */
  double minf; /* `*`the` `minimum` `objective` `value,` `upon` `return`*` */
  if (nlopt_optimize(opt, x, &minf) < 0) {
      printf("nlopt failed!\n");
  }
  else {
      printf("found minimum after %d evaluations\n", count);
      printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], minf);
  }

  nlopt_destroy(opt);
}
