/* NLOPT from MAD C interface
compile & run:
gcc -std=c99 -Wall -W -pedantic -O3 nlopt_ex1m.c mad_log.c ../mad_nlopt.c \
    -o nlopt_ex1m -I.. -I../../lib/nlopt/src/api/ ../../bin/macosx/libnlopt.a
./nlopt_ex1m [algorithm]
*/

#include <math.h>
#include <stdio.h>
#include <assert.h>
#include <mad_nlopt.h>

static int debug = 1;
static int count = 0;

static inline num_t sqr(num_t x) { return x*x; }

typedef struct {
  num_t a, b;
} my_constraint_data;

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
      grad[0] = 3*a*sqr(a*x[0]+b);
      grad[1] = -1.0;
    }
    r[0] = sqr(a*x[0]+b)*(a*x[0]+b) - x[1];
  }
  { num_t a = d[1].a, b = d[1].b;
    if (grad) {
      grad[2] = 3*a*sqr(a*x[0]+b);
      grad[3] = -1.0;
    }
    r[1] = sqr(a*x[0]+b)*(a*x[0]+b) - x[1];
  }
}

int main(int argc, const char *argv[])
{
  num_t fmin    = sqrt(8.0/27);
  num_t r   [2] = { 1.0/3, 8.0/27 };
  num_t x   [2] = { 1.234, 5.678 };
  num_t lb  [2] = { -HUGE_VAL, 0 };
  num_t ltol[2] = { 1e-8, 1e-8 };
  my_constraint_data data[2] = { {2,0}, {-1,1} };

  int algo = nlopt_algorithm_from_string(argc>1 ? argv[1] : "LN_COBYLA");
  if (algo == -1) {
    printf("invalid algorithm: %s (set to LN_COBYLA)\n", argv[1]);
    algo = NLOPT_LN_COBYLA;
  }

  nlopt_args_t arg = { 0 };

  // algorithm
  arg.algo   = algo
  // objective
  arg.fun    = myfunc;
  arg.fmin   = -INFINITY;
  arg.debug  = debug;
  // variables
  arg.n      = 2;
  arg.x      = x;
  arg.xmin   = lb;
  arg.xrtol  = 1e-4;
  // inequalities
  arg.q      = 2;
  arg.lfun   = myconstraints;
  arg.ltol   = ltol;
  arg.ldat   = data;

  mad_nlopt(&arg);

  num_t fval = arg.fval;
  printf("method: %s\n", nlopt_algorithm_name(algo));
  if (arg.status < 0) {
    printf("nlopt failed! reason: %d, count: %d\n", arg.status, count);
  }
  else {
      printf("found minimum after %d evaluations, reason: %d\n", count, arg.status);
      printf("found minimum at f(%.6e,%.6e) = %.10e (%.10e)\n", x[0], x[1], fval, fmin);
      printf("relative errors: x0=%.6e, x1=%.6e, f(x0,x1)=%.6e\n",
             (x[0]-r[0])/x[0], (x[1]-r[1])/x[1], (fval-fmin)/fval);

  }
}
