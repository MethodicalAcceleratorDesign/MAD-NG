#ifndef MAD_NLOPT_H
#define MAD_NLOPT_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Non Linear Optimization module interface
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

  Purpose:
  - interface to NLOpt library for LuaJIT

 o-----------------------------------------------------------------------------o
 */

#include <nlopt.h>
#include "mad_defs.h"

// -- interface ---------------------------------------------------------------o

typedef struct nlopt_args nlopt_args_t;

void mad_nlopt (nlopt_args_t *args);

// -- types -------------------------------------------------------------------o

// Define explicitly nlopt_func and nlopt_mfunc here to ensure FFI x-check
typedef num_t (objective_t )(u32_t n, const num_t *x, num_t *grad, void *data);
typedef void  (constraint_t)(u32_t m, num_t *res,
                             u32_t n, const num_t* x, num_t* grad, void* data);

struct nlopt_args {
  // algorithm
  int             algorithm;

  // variables [n] (required)
  num_t          *x;
  const num_t    *xtol_rel;
  const num_t    *xtol_abs;
  const num_t    *xstep;

  // objective function [n] (required)
  ssz_t           n;
  objective_t    *f;
  void           *fdat;
  const num_t    *ftol_rel;
  const num_t    *ftol_abs;

  // equality constraints [p] (optional)
  ssz_t           p;
  constraint_t   *c_eq;
  void           *cdat_eq;
  const num_t    *ctol_eq;

  // inequality constraints [q] (optional)
  ssz_t           q;
  constraint_t   *c_le;
  void           *cdat_le;
  const num_t    *ctol_le;

  // bounds constraint [n] (optional)
  const num_t    *xmin;
  const num_t    *xmax;

  // stop criteria
  const int      *maxcall;
  const num_t    *maxtime;

  // returned values
  num_t           fmin;
  int             status;
};

// -- end ---------------------------------------------------------------------o

#endif // MAD_NLOPT_H
