/*
compile & run:
gcc -std=c99 -Wall -W -pedantic -O3 nlopt_ex2.c mad_log.c ../mad_nlopt.c \
    -o nlopt_ex2 -I.. -I../../lib/nlopt/src/api/ ../../bin/macosx/libnlopt.a
./nlopt_ex2 [algorithm]
*/

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <nlopt.h>

static int debug = 1;
static int count = 0;

static inline double sqr(double x) { return x*x; }

static double
myfunc(unsigned n, const double *x, double *grad, void *data)
{
  assert(n == 5); assert(x && !grad && !data);
  ++count;

  double prod = 1;
  for (unsigned i=0; i<5; i++) prod *= x[i];

  double fval = exp(prod) - 0.5*sqr(x[0]*x[0]*x[0] + x[1]*x[1]*x[1] + 1);

  if (debug) {
    printf("objective (%s)\n", grad != NULL ? "D" : "N");
    for (unsigned i=0; i<n; i++) printf("x[%d]=%.16e\n", i+1, x[i]);
    printf("fval=%.16e\n", fval);
  }

  return fval;
}

static void
myconstraints(unsigned m, double *r, unsigned n, const double *x, double *grad, void *data)
{
  assert(m == 3 && n == 5); assert(r && x && !grad && !data);

  if (debug) {
    printf("constraint (%s)\n", grad != NULL ? "D" : "N");
    for (unsigned i=0; i<5; i++) printf("x[%d]=%.16e\n", i+1, x[i]);
  }

  double sumsq = 0;
  for (unsigned i=0; i<5; i++) sumsq += sqr(x[i]);

  r[0] = sumsq - 10;
  r[1] = x[1]*x[2] - 5*x[3]*x[4];
  r[2] = x[0]*x[0]*x[0] + x[1]*x[1]*x[1] + 1;
}

int main(int argc, const char *argv[])
{
  double fmin = 0.053950;
  double r[5] = { -1.71714, 1.59571, 1.82725, -0.763643, -0.763643 };
  double x[5] = {-1.8, 1.7, 1.9, -0.8, -0.8};

  int algo = nlopt_algorithm_from_string(argc>1 ? argv[1] : "LN_COBYLA");
  if (algo == -1) {
    printf("invalid algorithm: %s (set to LN_COBYLA)\n", argv[1]);
    algo = NLOPT_LN_COBYLA;
  }

  nlopt_opt opt = nlopt_create(algo, 5);

  nlopt_set_min_objective(opt, myfunc, NULL);
  nlopt_add_equality_mconstraint(opt, 3, myconstraints, NULL, NULL);

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
    printf("found minimum at f(%.6e,%.6e,%.6e,%.6e,%.6e) = %.8e [%.8e]\n",
           x[0], x[1], x[2], x[3], x[4], fval, fmin);
    printf("relative errors: x0=%.6e, x1=%.6e, x2=%.6e, x3=%.6e,\n"
           "                 x4=%.6e, f(x0,x1)=%.6e\n",
       (x[0]-r[0])/x[0], (x[1]-r[1])/x[1], (x[2]-r[2])/x[2], (x[3]-r[3])/x[3],
       (x[4]-r[4])/x[4], (fval-fmin)/fval);
  }
}
