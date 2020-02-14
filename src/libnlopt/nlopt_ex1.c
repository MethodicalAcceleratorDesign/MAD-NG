/*
compile & run:
gcc -std=c99 -Wall -W -pedantic -O3 nlopt_ex1.c mad_log.c ../mad_nlopt.c \
    -o nlopt_ex1 -I.. -I../../lib/nlopt-2.6.1/src/api/ ../../bin/macosx/libnlopt.a
./nlopt_ex1

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
#include <mad_nlopt.h>

typedef struct {
  num_t a, b;
} my_constraint_data;

static int count = 0;

static num_t
myfunc(u32_t n, const num_t *x, num_t *grad, void *data)
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
myconstraints(u32_t m, num_t *r, u32_t n, const num_t *x, num_t *grad, void *data)
{
  assert(m == 2 && n == 2); assert(x && data);
  my_constraint_data *d = (my_constraint_data*) data;

  { num_t a = d[0].a, b = d[0].b;
    if (grad) {
      grad[0] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
      grad[1] = -1.0;
    }
    r[0] = ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
  }
  { num_t a = d[1].a, b = d[1].b;
    if (grad) {
      grad[2] = 3 * a * (a*x[0] + b) * (a*x[0] + b);
      grad[3] = -1.0;
    }
    r[1] = ((a*x[0] + b) * (a*x[0] + b) * (a*x[0] + b) - x[1]);
  }
}

int main(void)
{
  num_t x   [2] = { 1.234, 5.678 };
  num_t lb  [2] = { -HUGE_VAL, 0 };
  num_t ltol[2] = { 1e-8, 1e-8 };
  my_constraint_data data[2] = { {2,0}, {-1,1} };

  nlopt_args_t arg = { 0 };

  // algorithm
  arg.algo  = NLOPT_LN_COBYLA; // NLOPT_LN_COBYLA, NLOPT_LD_MMA
  // objective
  arg.fun   = myfunc;
  // variables
  arg.n     = 2;
  arg.x     = x;
  arg.xmin  = lb;
  arg.xrtol = 1e-4;
  // inequalities
  arg.q     = 2;
  arg.lfun  = myconstraints;
  arg.ltol  = ltol;
  arg.ldat  = data;

  mad_nlopt(&arg);

  if (arg.status < 0) {
    printf("nlopt failed! reason: %d, count: %d\n", arg.status, count);
  }
  else {
    printf("found minimum after %d evaluations\n", count);
    printf("found minimum at f(%g,%g) = %0.10g\n", x[0], x[1], arg.fval);
  }
}
