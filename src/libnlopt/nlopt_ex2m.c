/* NLOPT from MAD C interface
compile & run:
gcc -std=c99 -Wall -W -pedantic -O3 nlopt_ex2m.c mad_log.c ../mad_nlopt.c \
    -o nlopt_ex2m -I.. -I../../lib/nlopt/src/api/ ../../bin/macosx/libnlopt.a
./nlopt_ex2m [algorithm]

output:
NLOPT_LN_COBYLA:
found minimum after 13 evaluations
found minimum at f(-1.85006,1.61757,1.52261,-0.801924,-0.811983) = -0.5534111746
*/

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <mad_nlopt.h>

static int debug = 1;
static int count = 0;

static inline num_t sqr(num_t x) { return x*x; }

static num_t
myfunc(u32_t n, const num_t *x, num_t *grad, void *data)
{
  assert(n == 5); assert(x && !grad && !data);
  ++count;

  if (debug) {
    printf("objective (%s)\n", grad != NULL ? "D" : "N");
    for (unsigned i=0; i<n; i++) printf("x[%d]=%.16e\n", i+1, x[i]);
  }

  num_t prod = 1;
  for (u32_t i=0; i<n; i++) prod *= x[i];

  return exp(prod) - 0.5*sqr(x[0]*sqr(x[0]) + x[1]*sqr(x[1]) + 1);
}

static void
myconstraints(u32_t m, num_t *r, u32_t n, const num_t *x, num_t *grad, void *data)
{
  assert(m == 3 && n == 5); assert(r && x && !grad && !data);

  if (debug) {
    printf("constraint (%s)\n", grad != NULL ? "D" : "N");
    for (unsigned i=0; i<n; i++) printf("x[%d]=%.16e\n", i+1, x[i]);
  }

  num_t sumsq = 0;
  for (u32_t i=0; i<n; i++) sumsq += sqr(x[i]);

  r[0] = sumsq - 10;
  r[1] = x[1]*x[2] - 5*x[3]*x[4];
  r[2] = x[0]*sqr(x[0]) + x[1]*sqr(x[1]) + 1;
}

int main(int argc, const char *argv[])
{
  num_t fmin = 0.053950;
  num_t r[5] = { -1.71714, 1.59571, 1.82725, -0.763643, -0.763643 };
  num_t x[5] = {-1.8, 1.7, 1.9, -0.8, -0.8};

  int algo = nlopt_algorithm_from_string(argc>1 ? argv[1] : "LN_COBYLA");
  if (algo == -1) {
    printf("invalid algorithm: %s (set to LN_COBYLA)\n", argv[1]);
    algo = NLOPT_LN_COBYLA;
  }

  nlopt_args_t arg = { 0 };

  // algorithm
  arg.algo  = algo;
  // objective
  arg.fun   = myfunc;
  arg.fmin  = -INFINITY;
  arg.debug = debug;
  // variables
  arg.n     = 5;
  arg.x     = x;
  // equalities
  arg.p     = 3;
  arg.efun  = myconstraints;

  mad_nlopt(&arg);

  num_t fval = arg.fval;
  printf("method: %s\n", nlopt_algorithm_name(algo));
  if (arg.status < 0) {
    printf("nlopt failed! reason: %d, count: %d\n", arg.status, count);
  }
  else {
    printf("found minimum after %d evaluations, reason: %d\n", count, arg.status);
    printf("found minimum at f(%g,%g,%g,%g,%g) = %0.10g\n",
           x[0], x[1], x[2], x[3], x[4], fval);
    printf("relative errors: x0=%.6e, x1=%.6e, x2=%.6e, x3=%.6e,\n"
           "                 x4=%.6e, f(x0,x1)=%.6e\n",
       (x[0]-r[0])/x[0], (x[1]-r[1])/x[1], (x[2]-r[2])/x[2], (x[3]-r[3])/x[3],
       (x[4]-r[4])/x[4], (fval-fmin)/fval);
  }
}
