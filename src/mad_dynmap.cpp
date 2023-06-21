/*
 o-----------------------------------------------------------------------------o
 |
 | Dynamic maps implementation in mixed C/C++
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
*/

extern "C" {
#include "mad_cst.h"
#include "mad_dynmap.h"
}

#include "mad_tpsa.hpp"

// --- types ------------------------------------------------------------------o

extern "C" {

struct mflw {
  num_t el, eld;
  int npar, T;
  par_t *par;
  map_t *map;
};

struct part {
  num_t x, px, y, py, t, pt;
  num_t beta;
};

struct damap {
  tpsa_t *x, *px, *y, *py, *t, *pt;
  num_t beta;
};

}

// --- implementation ---------------------------------------------------------o

using namespace mad;

// --- helpers ---

const num_t minlen = mad_cst_MINLEN;
const num_t minang = mad_cst_MINANG;
const num_t minstr = mad_cst_MINSTR;

inline num_t sqr (       num_t  a) { return a*a; }
inline tpsa  sqr (const tpsa_t &a) { return a*a; }
inline tpsa  sqr (const tpsa_t *a) { return sqr(*a); }
inline tpsa  sqr (const tpsa   &a) { return sqr(*a); }

inline num_t dp2 (par_t &p) {
  return 1 + (2/p.beta)*p.pt + sqr(p.pt);
}

inline tpsa dp2 (map_t &p) {
  return 1 + (2/p.beta)**p.pt + sqr(p.pt);
}

inline num_t pz2 (par_t &p) {
  return 1 + (2/p.beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py);
}

inline tpsa pz2 (map_t &p) {
  return 1 + (2/p.beta)**p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py);
}

inline num_t dp (par_t &p) { return std::sqrt(dp2(p)); }
inline num_t pz (par_t &p) { return std::sqrt(pz2(p)); }
inline tpsa  dp (map_t &p) { return mad::sqrt(dp2(p)); }
inline tpsa  pz (map_t &p) { return mad::sqrt(pz2(p)); }

// --- DKD maps ---

void
mad_trk_strex_drift_r (elem_t *e, mflw_t *m, num_t lw, int istp)
{
  (void)e; (void)istp;
  if (std::abs(m->el*lw) < minlen) return;

  num_t l = m->el*lw, ld = m->eld*lw;
  int T = m->T;

  FOR(i,m->npar) {
    par_t &p = m->par[i];
    num_t l_pz = l/pz(p);

    p.x += p.px*l_pz;
    p.y += p.py*l_pz;
    p.t -= l_pz*(1/p.beta+p.pt) + (1-T)*(ld/p.beta);
  }
}

void
mad_trk_strex_drift_t (elem_t *e, mflw_t *m, num_t lw, int istp)
{
  (void)e; (void)istp;
  if (std::abs(m->el*lw) < minlen) return;

  num_t l = m->el*lw, ld = m->eld*lw;
  int T = m->T;

  FOR(i,m->npar) {
    map_t &p = m->map[i];
    const tpsa l_pz = l/pz(p);

    ref(p.x) += p.px*l_pz;
    ref(p.y) += p.py*l_pz;
    ref(p.t) -= l_pz*(1/p.beta+*p.pt) + (1-T)*(ld/p.beta);
  }
}

void
mad_trk_strex_kick_t (elem_t *e, mflw_t *m, num_t lw, int istp)
{
  (void)e; (void)m; (void)lw; (void)istp;
}

void
mad_trk_curex_drift_t (elem_t *e, mflw_t *m, num_t lw, int istp)
{
  (void)e; (void)m; (void)lw; (void)istp;
}

void
mad_trk_curex_kick_r (elem_t *e, mflw_t *m, num_t lw, int istp)
{
  (void)e; (void)m; (void)lw; (void)istp;
}

void
mad_trk_curex_kick_t (elem_t *e, mflw_t *m, num_t lw, int istp)
{
  (void)e; (void)m; (void)lw; (void)istp;
}

void mad_trk_slice_r (elem_t *e, mflw_t *m, num_t lw, trkfun *dft, trkfun *kck)
{
  (void)e; (void)m; (void)lw; (void)dft; (void)kck;
}

void mad_trk_slice_t (elem_t *e, mflw_t *m, num_t lw, trkfun *dft, trkfun *kck)
{
  (void)e; (void)m; (void)lw; (void)dft; (void)kck;
}

// ----------------------------------------------------------------------------o
