/* NLOPT direct
compile & run:
gcc -std=c99 -Wall -W -pedantic -O3 nlopt_ex1.c \
    -o nlopt_ex1 -I../../lib/nlopt/src/api/ ../../bin/macosx/libnlopt.a
./nlopt_ex1 [algorithm]
*/

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <nlopt.h>

static int debug = 1;
static int count = 0;

typedef struct {
  double a, b;
} my_constraint_data;

static double
myfunc(unsigned n, const double *x, double *grad, void *data)
{
  assert(n == 2); assert(!data);
  double fval = sqrt(x[1]);
  ++count;

  if (debug) {
    printf("objective (%s)\n", grad != NULL ? "D" : "N");
    for (unsigned i=0; i<n; i++) printf("x[%d]=%.16e\n", i+1, x[i]);
    printf("fval=%.16e\n", fval);
  }

  if (grad) {
    grad[0] = 0.0;
    grad[1] = 0.5 / fval;
  }
  return fval;
}

static void
myconstraints(unsigned m, double *r, unsigned n, const double *x, double *grad, void *data)
{
  assert(m == 2 && n == 2); assert(x && data);
  my_constraint_data *d = (my_constraint_data*) data;

  if (debug) {
    printf("constraint (%s)\n", grad != NULL ? "D" : "N");
    for (unsigned i=0; i<n; i++) printf("x[%d]=%.16e\n", i+1, x[i]);
  }

  for (int i=0; i<=1; i++) {
    double a = d[i].a, b = d[i].b;
    if (grad) {
      grad[2*i+0] = 3*a*(a*x[0]+b)*(a*x[0]+b);
      grad[2*i+1] = -1.0;
    }
    r[i] = (a*x[0]+b)*(a*x[0]+b)*(a*x[0]+b) - x[1];
  }
}

int main(int argc, const char *argv[])
{
  double fmin   = sqrt(8.0/27);
  double r  [2] = { 1.0/3, 8.0/27 };
  double x  [2] = { 1.234, 5.678 };
  double lb [2] = { -HUGE_VAL, 0 };
  double tol[2] = { 1e-8, 1e-8 };
  my_constraint_data data[2] = { {2,0}, {-1,1} };

  int algo = nlopt_algorithm_from_string(argc>1 ? argv[1] : "LN_COBYLA");
  if (algo == -1) {
    printf("invalid algorithm: %s (set to LN_COBYLA)\n", argv[1]);
    algo = NLOPT_LN_COBYLA;
  }

  nlopt_opt opt = nlopt_create(algo, 2);

  nlopt_set_lower_bounds(opt, lb);
  nlopt_set_min_objective(opt, myfunc, NULL);
  nlopt_add_inequality_mconstraint(opt, 2, myconstraints, data, tol);
  nlopt_set_xtol_rel(opt, 1e-4);

  double fval;
  int status = nlopt_optimize(opt, x, &fval);

  printf("method: %s\n", nlopt_algorithm_name(algo));
  if (status < 0) {
    printf("nlopt failed! reason: %s, count: %d\n",
           nlopt_result_to_string(status), count);
  }
  else {
    printf("found minimum after %d evaluations, reason: %s\n",
           count, nlopt_result_to_string(status));
    printf("found minimum at f(%.6e,%.6e) = %.8e [%.8e]\n",
           x[0], x[1], fval, fmin);
    printf("relative errors: x0=%.6e, x1=%.6e, f(x0,x1)=%.6e\n",
           (x[0]-r[0])/x[0], (x[1]-r[1])/x[1], (fval-fmin)/fval);
  }

  nlopt_destroy(opt);
}
