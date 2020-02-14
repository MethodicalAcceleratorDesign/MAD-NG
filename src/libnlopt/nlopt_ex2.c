/*
compile & run:
gcc -std=c99 -Wall -W -pedantic -O3 nlopt_ex2.c mad_log.c ../mad_nlopt.c \
    -o nlopt_ex2 -I.. -I../../lib/nlopt-2.6.1/src/api/ ../../bin/macosx/libnlopt.a
./nlopt_ex2

output:
NLOPT_LN_COBYLA:
found minimum after 13 evaluations
found minimum at f(-1.85006,1.61757,1.52261,-0.801924,-0.811983) = -0.5534111746
*/

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <mad_nlopt.h>

static int count = 0;

static inline num_t sqr(num_t x) { return x*x; }

static num_t
myfunc(u32_t n, const num_t *x, num_t *grad, void *data)
{
  assert(n == 5); assert(x && !grad && !data);
  ++count;

  num_t prod = 1;
  for (u32_t i=0; i<n; i++) prod *= x[i];

  return exp(prod) - 0.5*sqr(x[0]*sqr(x[0]) + x[1]*sqr(x[1]) + 1);
}

static void
myconstraints(u32_t m, num_t *r, u32_t n, const num_t *x, num_t *grad, void *data)
{
  assert(m == 3 && n == 5); assert(r && x && !grad && !data);

  num_t sumsq = 0;
  for (u32_t i=0; i<n; i++) sumsq += sqr(x[i]);

  r[0] = sumsq - 10;
  r[1] = x[1]*x[2] - 5*x[3]*x[4];
  r[2] = x[0]*sqr(x[0]) + x[1]*sqr(x[1]) + 1;
}

int main(void)
{
  num_t x[5] = {-1.8, 1.7, 1.9, -0.8, -0.8};

  nlopt_args_t arg = { 0 };

  // algorithm
  arg.algo = NLOPT_LN_COBYLA;
  // objective
  arg.fun  = myfunc;
  // variables
  arg.n    = 5;
  arg.x    = x;
  // equalities
  arg.p    = 3;
  arg.efun = myconstraints;

  mad_nlopt(&arg);

  if (arg.status < 0) {
    printf("nlopt failed! reason: %d, count: %d\n", arg.status, count);
  }
  else {
    printf("found minimum after %d evaluations\n", count);
    printf("found minimum at f(%g,%g,%g,%g,%g) = %0.10g\n",
           x[0], x[1], x[2], x[3], x[4], arg.fval);
  }
}
