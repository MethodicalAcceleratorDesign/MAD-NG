/*
compile & run:
gcc -std=c99 -Wall -W -pedantic -O3 nlopt_ex0.c \
    -o nlopt_ex0 -I../../lib/nlopt-git/src/api/ ../../bin/macosx/libnlopt.a
./nlopt_ex0

output:
NLOPT_LD_MMA:
found minimum after 11 evaluations, reason: 4
found minimum at f(0.333333,0.296296) = 0.5443310477
relative errors: x0=0.000000, x1=-0.000000, f(x0,x1)=-0.000000
NLOPT_LN_COBYLA:
found minimum after 31 evaluations, reason: 4
found minimum at f(0.333329,0.296200) = 0.5442423017
relative errors: x0=-0.000012, x1=-0.000326, f(x0,x1)=-0.000163
*/

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <nlopt.h>

typedef struct {
  double a, b;
} my_constraint_data;

static int count = 0;

static inline double sqr(double x) { return x*x; }

static double
myfunc(unsigned n, const double *x, double *grad, void *data)
{
  assert(n == 2); assert(!data);
  ++count;
  if (grad) {
    grad[0] = 0.0;
    grad[1] = 0.5 / sqrt(x[1]);
  }
  return sqrt(x[1]);
}

static void
myconstraints(unsigned m, double *r, unsigned n, const double *x, double *grad, void *data)
{
  assert(m == 2 && n == 2); assert(x && data);
  my_constraint_data *d = (my_constraint_data*) data;

  { double a = d[0].a, b = d[0].b;
    if (grad) {
      grad[0] = 3*a*sqr(a*x[0]+b);
      grad[1] = -1.0;
    }
    r[0] = sqr(a*x[0]+b)*(a*x[0]+b) - x[1];
  }
  { double a = d[1].a, b = d[1].b;
    if (grad) {
      grad[2] = 3*a*sqr(a*x[0]+b);
      grad[3] = -1.0;
    }
    r[1] = sqr(a*x[0]+b)*(a*x[0]+b) - x[1];
  }
}

int main(void)
{
  double x  [2] = { 1.234, 5.678 };
  double lb [2] = { -HUGE_VAL, 0 };
  double tol[2] = { 1e-8, 1e-8 };
  my_constraint_data data[2] = { {2,0}, {-1,1} };

  nlopt_opt opt = nlopt_create(NLOPT_LN_COBYLA, 2); // NLOPT_LN_COBYLA, NLOPT_LD_MMA

  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_min_objective(opt, myfunc, NULL);
  nlopt_add_inequality_mconstraint(opt, 2, myconstraints, data, tol);
  nlopt_set_xtol_rel(opt, 1e-4);

  double minf;
  int status = nlopt_optimize(opt, x, &minf);
  if (status < 0) {
    printf("nlopt failed! reason: %d, count: %d\n", status, count);
  }
  else {
      printf("found minimum after %d evaluations, reason: %d\n", count, status);
      printf("found minimum at f(%.6f,%.6f) = %.10f\n", x[0], x[1], minf);
      printf("relative errors: x0=%.6f, x1=%.6f, f(x0,x1)=%.6f\n",
             (x[0]-1.0/3)/x[0], (x[1]-8.0/27)/x[1], (minf-sqrt(8.0/27))/minf);
  }

  nlopt_destroy(opt);
}
