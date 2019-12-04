#ifndef MAD_NLOPT_H
#define MAD_NLOPT_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Non Linear Optimization module interface
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

  Purpose:
  - interface to NLOpt library for LuaJIT

 o-----------------------------------------------------------------------------o
 */

#include <nlopt.h>
#include "mad_defs.h"

// -- interface ---------------------------------------------------------------o

typedef struct nlopt_args nlopt_args_t;

void mad_nlopt (nlopt_args_t *args);
void mad_nlopt_srand (u64_t seed);

// -- types -------------------------------------------------------------------o

// Define explicitly nlopt_func and nlopt_mfunc here to ensure FFI x-check
typedef num_t (nlopt_obj_t)(u32_t n, const num_t *x, num_t *grad, void *data);
typedef void  (nlopt_cts_t)(u32_t m, num_t *res,
                            u32_t n, const num_t* x, num_t* grad, void* data);

struct nlopt_args {
  // algorithm
        int    algo;

  // objective function (required)
  nlopt_obj_t *fun;
        int    fdir; // 1: maximize, -1: minimize (default)
        num_t  fval;
        num_t  ftol;
        num_t  fstop;

  // state variables [n] (required)
        ssz_t  n;
        num_t *x;

  // initial value, step, tolerances and bounds [n] (optional)
  const num_t *x0;
  const num_t *xstp;
  const num_t *xtol;
        num_t  rtol;
  const num_t *xmin;
  const num_t *xmax;

  // equality constraints [p] (optional)
        ssz_t  p;
  nlopt_cts_t *efun;
        num_t *edat;
  const num_t *etol;

  // inequality constraints [q] (optional)
        ssz_t  q;
  nlopt_cts_t *lfun;
        num_t *ldat;
  const num_t *ltol;

  // stop criteria (if >0)
        int    maxcall;
        num_t  maxtime;

  // returned values
        int    status;
        int    ncall;
};

// -- end ---------------------------------------------------------------------o

#endif // MAD_NLOPT_H
