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

// comment to remove mdump code (+2.5% code size), (basic) speed and C++ tests
#define TPSA_DBGMDUMP 1
#define TPSA_SPDTESTS 1
#define TPSA_CPPTESTS 1

// --- includes ---------------------------------------------------------------o

#include <type_traits>
#include "mad_tpsa.hpp"

extern "C" {
#include "mad_log.h"
#include "mad_cst.h"
#include "mad_mem.h"
#include "mad_mat.h"
#include "mad_dynmap.h"
}

// --- types ------------------------------------------------------------------o

extern "C" {
enum { nmul_max=22, snm_max=(nmul_max+1)*(nmul_max+2)/2 };
}

template <typename M, typename MT=M::MT, typename MP=M::MP>
struct cflw { // must be identical to def of 3 variants in madl_etck.mad!!
  str_t name;
  int dbg;

  // beam
  num_t pc, beta, betgam;
  int charge;

  // directions, path
  int sdir, edir, pdir, T, Tbak;

// start of polymorphic section

  // element data
  MP el, eld, elc, lrad;
  MP eh, ang, mang;

  // quad, solenoid, multipole, esptum, rfcav, sine & cosine
  MP k1, ks, volt, freq, lag, sa, ca;
  int nbsl;

  // fringes
  int frng, fmax;
  MP e, h, a, fint, hgap, f1, f2;

  // patches, misalignments & tilt
  bool rot, trn;
  MP dx,   dy,   ds;
  MP dthe, dphi, dpsi, tlt;

  // multipoles
  int nmul;
  MP  knl[nmul_max];
  MP  ksl[nmul_max];

  // curved multipoles
  int snm;
  MP  bfx[snm_max];
  MP  bfy[snm_max];

  // particles/damaps/parametric_damaps (must be last!!)
  int npar;
  MT **par;
};

struct par_t { // particles
  // traits
  using T  = num_t;         // type of variables  in maps (num_t  or tpsa)
  using P  = num_t;         // type of parameters in maps (num_t  or tpsa)
  using R  = num_t&;        // type of prms refs  in maps (num_t  or tpsa_ref)
  using A  = num_t*;        // type of prms array in maps (num_t& or tpsa_refs)
  using MT = num_t;         // type of variables  in mflw (num_t  or tpsa_t*)
  using MP = num_t;         // type of parameters in mflw (num_t  or tpsa_t*)
  // ctor
  par_t(struct cflw<par_t> &m, int i)
    : x(m.par[i][0]), px(m.par[i][1]),
      y(m.par[i][2]), py(m.par[i][3]),
      t(m.par[i][4]), pt(m.par[i][5]) {}
  // members
  num_t &x, &px, &y, &py, &t, &pt;
};

struct map_t { // damaps
  // traits
  using T  = mad::tpsa;     // type of variables  in maps (num_t  or tpsa)
  using P  = num_t;         // type of parameters in maps (num_t  or tpsa)
  using R  = num_t&;        // type of prms refs  in maps (num_t  or tpsa_ref)
  using A  = num_t*;        // type of prms array in maps (num_t& or tpsa_refs)
  using MT = tpsa_t*;       // type of variables  in mflw (num_t  or tpsa_t*)
  using MP = num_t;         // type of parameters in mflw (num_t  or tpsa_t*)
  // ctor
  map_t(struct cflw<map_t> &m, int i)
    : x(m.par[i][0]), px(m.par[i][1]),
      y(m.par[i][2]), py(m.par[i][3]),
      t(m.par[i][4]), pt(m.par[i][5]) {}
  // members
  mad::tpsa_ref x, px, y, py, t, pt;
};

struct prm_t { // parametric damaps
  // traits
  using T  = mad::tpsa;     // type of variables  in maps (num_t  or tpsa)
  using P  = mad::tpsa;     // type of parameters in maps (num_t  or tpsa)
  using R  = mad::tpsa_ref; // type of prms refs  in maps (num_t  or tpsa_ref)
  using A  = mad::tpsa_refs;// type of prms array in maps (num_t& or tpsa_refs)
  using MT = tpsa_t*;       // type of variables  in mflw (num_t  or tpsa_t*)
  using MP = tpsa_t*;       // type of parameters in mflw (num_t  or tpsa_t*)
  // ctor
  prm_t(struct cflw<prm_t> &m, int i)
    : x(m.par[i][0]), px(m.par[i][1]),
      y(m.par[i][2]), py(m.par[i][3]),
      t(m.par[i][4]), pt(m.par[i][5]) {}
  // members
  mad::tpsa_ref x, px, y, py, t, pt;
};

extern "C" {
union mflw_ {
  struct cflw<par_t> rflw;
  struct cflw<map_t> tflw;
  struct cflw<prm_t> pflw;
};

const size_t mad_cflw_rsize = sizeof(struct cflw<par_t>);
const size_t mad_cflw_tsize = sizeof(struct cflw<map_t>);
const size_t mad_cflw_psize = sizeof(struct cflw<prm_t>);
const size_t mad_cflw_size  = sizeof(union  mflw_      );
} // extern "C"

// --- implementation ---------------------------------------------------------o

using namespace mad;

// --- debug ------------------------------------------------------------------o

#if TPSA_DBGMDUMP // set to 0 to remove debug code, ~2.5% of code size
#define mdump(n) if (m.dbg) mdump<M>(m, __func__, n)
#else
#define mdump(n)
#endif

template <typename M, typename T=M::T>
inline void (mdump) (cflw<M> &m, str_t s, int n)
{
  if (!m.dbg) return;

  char fun[30];
  snprintf(fun, 30, "%s:%d", s, n);
  printf("@@ %-15s %-15s ", m.name, fun);

  if (!m.npar) { printf("no particle found\n"); return; }

  M p(m,0);
  if constexpr (std::is_floating_point<T>::value)
    printf("% -.16e  % -.16e  % -.16e  % -.16e  % -.16e  % -.16e\n",
                p.x,    p.px,     p.y,    p.py,     p.t,    p.pt);
  else
    printf("% -.16e  % -.16e  % -.16e  % -.16e  % -.16e  % -.16e\n",
             p.x[0], p.px[0],  p.y[0], p.py[0],  p.t[0], p.pt[0]);
}

// --- constants --------------------------------------------------------------o

const num_t one    = 1;
const num_t zero   = 0;
const num_t minlen = mad_cst_MINLEN;
const num_t minang = mad_cst_MINANG;
const num_t minstr = mad_cst_MINSTR;

const num_t    pi_clight = mad_cst_PI /mad_cst_CLIGHT;
const num_t twopi_clight = mad_cst_2PI/mad_cst_CLIGHT;

// --- multipoles -------------------------------------------------------------o

template <typename M, typename T=M::T, typename V, typename R=M::R>
inline void bxby (const cflw<M> &m, const V &x, const V &y, T &bx, T &by)
{
  bx = R(m.ksl[m.nmul-1]);
  by = R(m.knl[m.nmul-1]);

  if (m.nmul > 1) {
    T byt;
    RFOR(i,m.nmul-1) {
      byt = by*x - bx*y + R(m.knl[i]);
      bx  = by*y + bx*x + R(m.ksl[i]);
      by  = byt;
    }
  }
}

template <typename M, typename T=M::T, typename V, typename R=M::R>
inline void bxbyh (const cflw<M> &m, const V &x, const V &y, T &bx, T &by)
{
  bx = 0., by = 0.;

  int k = -1;
  T btx, bty;

  RFOR(i,m.snm) {
    btx = 0., bty = 0.;

    RFOR(j,m.snm-i-1) { ++k; // must skip the first iteration
      btx = (btx + R(m.bfx[k])) * y;
      bty = (bty + R(m.bfy[k])) * y;
    }

    ++k;
    bx = (bx + btx + R(m.bfx[k])) * x;
    by = (by + bty + R(m.bfy[k])) * x;
  }

  btx = 0., bty = 0.;
  RFOR(i,m.snm) { ++k;
    btx = (btx + R(m.bfx[k])) * y;
    bty = (bty + R(m.bfy[k])) * y;
  }

  bx += btx + R(m.bfx[k+1]); // better to enforce associativity in Lua.
  by += bty + R(m.bfy[k+1]);
}

// --- patches ----------------------------------------------------------------o

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void xrotation (cflw<M> &m, num_t lw, const V &dphi_)
{
  if (fabs(dphi_)+fabs(m.ang) < minang) return;
  mdump(0);
  lw *= m.edir;
  P a;
  if (fval(dphi_)) a = lw*dphi_; else a = lw*R(m.ang);
  P sa=sin(a), ca=cos(a), ta=tan(a);

  FOR(i,m.npar) {
    M p(m,i);
    T   pz = sqrt(1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = 1/pz;
    T _ptt = p.y/(1 - ta*p.py*_pz);
    T _pzt = ta*_pz*_ptt;

    // eq. 127 in Forest06
    p.y   = _ptt/ca;
    p.py  = ca*p.py + sa*pz;
    p.x  += _pzt*p.px;
    p.t  -= _pzt*(1/m.beta+p.pt);
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void yrotation (cflw<M> &m, num_t lw, const V &dthe_)
{
  if (fabs(dthe_)+fabs(m.ang) < minang) return;
  mdump(0);
  lw *= -m.edir;
  P a;
  if (fval(dthe_)) a = lw*dthe_; else a = lw*R(m.ang);
  P sa=sin(a), ca=cos(a), ta=tan(a);

  FOR(i,m.npar) {
    M p(m,i);
    T   pz = sqrt(1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = inv(pz);
    T _ptt = p.x/(1 - ta*p.px*_pz);
    T _pzt = ta*_pz*_ptt;

    // eq. 127 in Forest06
    p.x   = _ptt/ca;
    p.px  = ca*p.px + sa*pz;
    p.y  += _pzt*p.py;
    p.t  -= _pzt*(1/m.beta+p.pt);
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void srotation (cflw<M> &m, num_t lw, const V &dpsi_)
{
  if (fabs(dpsi_)+fabs(m.ang) < minang) return;
  mdump(0);
  lw *= m.edir;
  P a;
  if (fval(dpsi_)) a = lw*dpsi_; else a = lw*R(m.ang);
  P sa=sin(a), ca=cos(a);

  FOR(i,m.npar) {
    M p(m,i);
    T nx  = ca*p.x  + sa*p.y;
    T npx = ca*p.px + sa*p.py;

    p.y  = ca*p.y  - sa*p.x;
    p.py = ca*p.py - sa*p.px;
    p.x  = nx;
    p.px = npx;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void translate (cflw<M> &m, num_t lw, const V &dx_, const V &dy_, const V &ds_)
{
  if (fabs(dx_)+fabs(m.dx) < minlen &&
      fabs(dy_)+fabs(m.dy) < minlen &&
      fabs(ds_)+fabs(m.ds) < minlen) return;
  mdump(0);
  P dx, dy, ds;
  if (fval(dx_)) dx = lw*m.edir*dx_; else dx = lw*m.edir*R(m.dx);
  if (fval(dy_)) dy = lw*m.edir*dy_; else dy = lw*m.edir*R(m.dy);
  if (fval(ds_)) ds = lw       *ds_; else ds = lw*       R(m.ds);

  if (fabs(ds) < minlen)
    FOR(i,m.npar) {
      M p(m,i);
      p.x -= dx;
      p.y -= dy;
    }
  else
    FOR(i,m.npar) {
      M p(m,i);
      T l_pz = ds*invsqrt(1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));

      p.x += l_pz*p.px - dx;
      p.y += l_pz*p.py - dy;
      p.t -= l_pz*(1/m.beta+p.pt);
    }
  mdump(1);
}

template <typename M, typename R=M::R>
inline void changeref (cflw<M> &m, num_t lw)
{
  bool trn = fabs(m.dx  )+fabs(m.dy  )+fabs(m.ds  ) >= minlen;
  bool rot = fabs(m.dthe)+fabs(m.dphi)+fabs(m.dpsi) >= minang;

  if (!trn && !rot) return;
  mdump(0);

  if (rot && lw > 0) {
    yrotation<M>(m,  1, R(m.dthe));
    xrotation<M>(m, -1, R(m.dphi));
    srotation<M>(m,  1, R(m.dpsi));
  }

  if (trn) translate<M>(m, m.sdir, R(m.dx), R(m.dy), R(m.ds));

  if (rot && lw < 0) {
    srotation<M>(m, -1, R(m.dpsi));
    xrotation<M>(m,  1, R(m.dphi));
    yrotation<M>(m, -1, R(m.dthe));
  }
  mdump(1);
}

// --- misalignments ----------------------------------------------------------o

template <typename M, typename R=M::R>
inline void misalignent (cflw<M> &m)
{
  mdump(0);
  num_t tb[3], t[3]={m.edir*fval(m.dx), m.edir*fval(m.dy), fval(m.ds)}; // t needs to be weighted before rotation

  if (m.rot && m.trn){
    num_t r[3*3];
    mad_mat_rotyxz(r, m.edir*fval(m.dphi),-m.edir*fval(m.dthe),-m.edir*fval(m.dpsi), false);
    mad_mat_mul(r, t, tb, 3, 1, 3);
  }

  if (m.rot && m.sdir > 0) {
    yrotation<M>(m,  1, R(m.dthe));
    xrotation<M>(m, -1, R(m.dphi));
    srotation<M>(m,  1, R(m.dpsi));
  }

  if (m.rot && m.trn) {
    translate<M>(m, m.sdir,
      R(m.dx)-m.edir*(t[0]-tb[0]), R(m.dy)-m.edir*(t[1]-tb[1]), R(m.ds)-(t[2]-tb[2]));
  }
  else if (m.trn) {
    translate<M>(m, m.sdir, R(m.dx), R(m.dy), R(m.ds));
  }

  if (m.rot && m.sdir < 0) {
    srotation<M>(m, -1, R(m.dpsi));
    xrotation<M>(m,  1, R(m.dphi));
    yrotation<M>(m, -1, R(m.dthe));
  }
  mdump(1);
}

template <typename M, typename R=M::R>
inline void misalignexi (cflw<M> &m)
{
  mdump(0);
  num_t rb[3*3], r[3*3],
        tb[3]  , t[3]={m.edir*fval(m.dx)  , m.edir*fval(m.dy)  , fval(m.ds)  },
                 a[3]={m.edir*fval(m.dphi), m.edir*fval(m.dthe), m.edir*fval(m.dpsi)};

  if (m.rot)
    mad_mat_rotyxz(r, a[0], -a[1], -a[2], true);

  // compute Rbar, Tbar
  mad_mat_rtbar(rb, tb, fabs(m.el), m.edir*fval(m.mang), fval(m.tlt), m.rot ? r:0, t);

  if (m.rot && m.sdir > 0) {
    num_t v[3];
    mad_mat_torotyxz(rb, v, true);
    srotation<M>(m,  1, R(m.dpsi)-m.edir*(a[2]-v[2]));
    xrotation<M>(m,  1,-R(m.dphi)+m.edir*(a[0]+v[0]));
    yrotation<M>(m,  1, R(m.dthe)-m.edir*(a[1]-v[1]));
  }

  translate<M>(m, -m.sdir,
     R(m.dx)-m.edir*(t[0]-tb[0]), R(m.dy)-m.edir*(t[1]-tb[1]), R(m.ds)-(t[2]-tb[2]));

  if (m.rot && m.sdir < 0) {
    num_t v[3];
    mad_mat_torotyxz(rb, v, true);
    yrotation<M>(m,  1,-R(m.dthe)+m.edir*(a[1]-v[1]));
    xrotation<M>(m,  1, R(m.dphi)-m.edir*(a[0]+v[0]));
    srotation<M>(m,  1,-R(m.dpsi)+m.edir*(a[2]-v[2]));
  }
  mdump(1);
}

template <typename M>
inline void misalign (cflw<M> &m, num_t lw)
{
  if (lw >= 0) misalignent<M>(m);
  else         misalignexi<M>(m);
}

// --- special maps -----------------------------------------------------------o

template <typename M, typename T=M::T, typename P=M::P>
inline void drift_adj (cflw<M> &m, const P &l)
{
  mdump(0);
  FOR(i,m.npar) {
    M p(m,i);
    T l_pz = l*invsqrt(1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));

    p.x += p.px*(l_pz-l);
    p.y += p.py*(l_pz-l);
    p.t  = p.t - (1/m.beta+p.pt)*l_pz + (1-m.T)/m.beta*l;
  }
  mdump(1);
}

// --- DKD maps ---------------------------------------------------------------o

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void strex_drift (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  if (fabs(m.el) < minlen) return;

  mdump(0);
  P l  = R(m.el)*lw;
  P ld = (fval(m.eld) ? R(m.eld) : R(m.el))*lw;

  FOR(i,m.npar) {
    M p(m,i);
    T l_pz = l*invsqrt(1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));

    p.x += p.px*l_pz;
    p.y += p.py*l_pz;
    p.t  = p.t - (1/m.beta+p.pt)*l_pz + (1-m.T)/m.beta*ld;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void strex_kick (cflw<M> &m, num_t lw, int is, bool no_k0l=false)
{                                           (void)is;
  if (!m.nmul || !m.charge) return;

  mdump(0);
  num_t wchg = lw*m.edir*m.charge;
  P dby;
  if (no_k0l) dby = R(m.knl[0]); else dby = 0.;
  T bx, by;

  FOR (i,m.npar) {
    M p(m,i);
    bxby(m, p.x, p.y, bx, by);

    p.px -= wchg*(by-dby);
    p.py += wchg* bx;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void strex_kicks (cflw<M> &m, num_t lw, M &p, T &pz)
{
  if (fabs(m.ks) < minstr || fabs(m.lrad) < minlen) return;

  num_t wchg = lw*m.edir*m.charge;
  P hss  = lw*R(m.lrad)*sqr(R(m.ks));
  T _dpp = inv(pz);
  T  ang = (0.5*wchg)*R(m.ks)*R(m.lrad)*_dpp;
  T  ca  = cos(ang), sa = sin(ang);

  T nx  = ca*p. x + sa*p. y;
  T npx = ca*p.px + sa*p.py;
  T ny  = ca*p. y - sa*p. x;
  T npy = ca*p.py - sa*p.px;
  T nt  = p.t - ang*(1/m.beta+p.pt)*(p.y*p.px - p.x*p.py)*sqr(_dpp);

  p.x  = nx;
  p.px = npx - 0.25 *hss*nx*_dpp;
  p.y  = ny;
  p.py = npy - 0.25 *hss*ny*_dpp;
  p.t  = nt  - 0.125*hss*(1/m.beta+p.pt)*(sqr(nx)+sqr(ny))*pow(_dpp,3);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void strex_kickhs (cflw<M> &m, num_t lw, int is)
{                                             (void)is;
  if ((!m.nmul && fabs(m.ks) < minstr) || !m.charge) return;

  mdump(0);
  num_t wchg = lw*m.edir*m.charge;
  T bx, by;

  FOR(i,m.npar) {
    M p(m,i);
    T pz = sqrt(1 + 2/m.beta*p.pt + sqr(p.pt));

    if (m.sdir == -1) strex_kicks(m, lw, p, pz);

    if (m.nmul > 0) {
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*by;
      p.py += wchg*bx;

      if (fabs(m.knl[0])+fabs(m.ksl[0]) >= minstr) {
        p.px += wchg* R(m.knl[0])*pz;
        p.py -= wchg* R(m.ksl[0])*pz;
        p.t  -= wchg*(R(m.knl[0])*p.x - R(m.ksl[0])*p.y)*(1/m.beta+p.pt)/pz;

        if (fabs(m.lrad) >= minlen) {
          p.px -= lw*sqr(R(m.knl[0]))/R(m.lrad)*p.x;
          p.py -= lw*sqr(R(m.ksl[0]))/R(m.lrad)*p.y;
        }
      }
    }

    if (m.sdir == 1) strex_kicks(m, lw, p, pz);
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void curex_drift (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  if (fabs(R(m.el)*R(m.eh)) < minang) return strex_drift<M>(m, lw, is);

  mdump(0);
  P ld  = (fval(m.eld) ? R(m.eld) : R(m.el))*lw;
  P ang = R(m.eh)*R(m.el)*lw*m.edir, rho = 1/R(m.eh)*m.edir;
  P ca  = cos(ang), sa = sin(ang), sa2 = sin(ang/2);

  FOR(i,m.npar) {
    M p(m,i);
    T   pz = sqrt(1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = inv(pz);
    T  pxt = p.px*_pz;
    T _ptt = inv(ca - sa*pxt);
    T  pst = (p.x+rho)*sa*_pz*_ptt;

    p.x  = (p.x + rho*(2*sqr(sa2) + sa*pxt))*_ptt;
    p.px = ca*p.px + sa*pz;
    p.y += pst*p.py;
    p.t  = p.t - pst*(1/m.beta+p.pt) + (1-m.T)/m.beta*ld;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void curex_kick (cflw<M> &m, num_t lw, int is, bool no_k0l=false)
{                                           (void)is;
  if (!m.nmul || !m.charge) return;
  if (fabs(R(m.el)*R(m.eh)) < minang) return strex_kick<M>(m, lw, is);

  mdump(0);
  num_t wchg = lw*m.edir*m.charge;
  T bx, by; bx = 0., by = R(m.knl[0]);

  FOR(i,m.npar) {
    M p(m,i);
    T r = 1+R(m.eh)*p.x*m.edir;

    if (m.snm) bxbyh(m, p.x, p.y, bx, by);

    p.px -= wchg*by*r;
    p.py += wchg*bx*r;

    if (no_k0l) p.px += (wchg*R(m.knl[0]))*r;
  }
  mdump(1);
}

// --- TKT maps ---------------------------------------------------------------o

// --- sbend ---

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void sbend_thick (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  if (fabs(m.knl[0]) < minstr || !m.charge) return curex_drift<M>(m, lw, is);

  mdump(0);
  P ld  = (fval(m.eld) ? R(m.eld) : R(m.el))*lw;
  P ang = R(m.eh)*R(m.el)*lw*m.edir, rho=1/R(m.eh)*m.edir;
  P k0q = R(m.knl[0])/R(m.el)*(m.edir*m.charge);
  P ca  = cos(ang), sa = sin(ang);

  FOR(i,m.npar) {
    M p(m,i);
    T  pw2 = 1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.py);
    T  pzx = sqrt(pw2 - sqr(p.px)) - k0q*(rho+p.x); // can be numerically unstable
    T  npx = sa*pzx + ca*p.px;
    T  dpx = ca*pzx - sa*p.px;
    T _ptt = invsqrt(pw2);
    T  dxs = (ang + asin(p.px*_ptt) - asin(npx*_ptt))/k0q;

    // eq. 126 in Forest06
    p.x  = (sqrt(pw2 - sqr(npx)) - dpx)/k0q - rho;  // can be numerically unstable
    p.px = npx;
    p.y += dxs*p.py;
    p.t  = p.t - dxs*(1/m.beta+p.pt) + (1-m.T)/m.beta*ld;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void sbend_thick_new (cflw<M> &m, num_t lw, int is)
{                                                (void)is;
  if (fabs(m.knl[0]) < minstr || !m.charge) return curex_drift<M>(m, lw, is);

  mdump(0);
  P ld  = (fval(m.eld) ? R(m.eld) : R(m.el))*lw;
  P ang = R(m.eh)*R(m.el)*lw*m.edir, rho=1/R(m.eh)*m.edir;
  P k0q = R(m.knl[0])/R(m.el)*(m.edir*m.charge);
  P ca  = cos(ang), sa = sin(ang), s2a = sin(2*ang);

  FOR(i,m.npar) {
    M p(m,i);
    T  pw2 = 1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.py);
    T   pz = sqrt(pw2 - sqr(p.px));
    T   xr = p.x+rho;
    T  pzx = pz - k0q*xr;
    T  npx = sa*pzx + ca*p.px;
    T  dpx = ca*pzx - sa*p.px;
    T _ptt = invsqrt(pw2);

    T  xt1 = -k0q*sqr(p.x) + 2*(pz*xr - (k0q*rho)*p.x) - k0q*sqr(rho);
    T   xi = p.px*_ptt;
    T zeta =  npx*_ptt;
    T  sxi = sqrt(1-sqr(xi));
    T    w = (ca*xi + sa*sxi) * sqrt(1-sqr(zeta));
    T    v = (sa*xi - ca*sxi) * zeta;
    T  xt2 = (s2a*p.px + sqr(sa)*(2*pz - k0q*xr)) * xr*sqr(_ptt) / (w - v);
    T  dxs = asinc(xt2*k0q)*xt2;

    // eq. 126 in Forest06 with modif. from Sagan
    p.x  = xt1/(dpx + sqrt(pw2 - sqr(npx))) - rho;
    p.px = npx;
    p.y += dxs*p.py;
    p.t  = p.t - dxs*(1/m.beta+p.pt) + (1-m.T)/m.beta*ld;
  }
  mdump(1);
}

// --- rbend ---

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void rbend_thick (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  if (fabs(m.knl[0]) < minstr || !m.charge) return strex_drift<M>(m, lw, is);

  mdump(0);
  P ld   = (fval(m.eld) ? R(m.eld) : R(m.el))*lw;
  P k0q  = R(m.knl[0])/R(m.el)*(m.edir*m.charge);
  P k0lq = R(m.knl[0])*(lw    * m.edir*m.charge);

  FOR(i,m.npar) {
    M p(m,i);
    T  npx = p.px - k0lq;
    T  pw2 = 1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.py);
    T _ptt = invsqrt(pw2);
    T   pz = sqrt(pw2 - sqr(p.px));
    T  pzs = sqrt(pw2 - sqr(npx));
    T  dxs = (asin(p.px*_ptt) - asin(npx*_ptt))/k0q;

    // eq. 126 in Forest06
    p.x += (pzs-pz)/k0q;
    p.px = npx;
    p.y += dxs*p.py;
    p.t  = p.t - dxs*(1/m.beta+p.pt) + (1-m.T)/m.beta*ld;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void rbend_thick_new (cflw<M> &m, num_t lw, int is)
{                                                (void)is;
  if (fabs(m.knl[0]) < minstr || !m.charge) return strex_drift<M>(m, lw, is);

  mdump(0);
  P l    = R(m.el)*lw;
  P ld   = (fval(m.eld) ? R(m.eld) : R(m.el))*lw;
  P k0q  = R(m.knl[0])/R(m.el)*(m.edir*m.charge);
  P k0lq = R(m.knl[0])*(lw    * m.edir*m.charge);

  FOR(i,m.npar) {
    M p(m,i);
    T  npx = p.px - k0lq;
    T  pw2 = 1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.py);
    T _ptt = invsqrt(pw2);
    T   xi = p.px*_ptt;
    T zeta =  npx*_ptt;
    T  xtd = xi*sqrt(1-sqr(zeta)) + zeta*sqrt(1-sqr(xi));
    T   xt = l*(2*p.px - k0lq)*sqr(_ptt) / xtd;
    T  dxs = asinc(xt*k0q)*xt;

    // eq. 126 in Forest06 with modif. from Sagan
    p.x += l*(2*p.px - k0lq) / (sqrt(pw2 - sqr(p.px)) + sqrt(pw2 - sqr(npx)));
    p.px = npx;
    p.y += dxs*p.py;
    p.t  = p.t - dxs*(1/m.beta+p.pt) + (1-m.T)/m.beta*ld;
  }
  mdump(1);
}

// --- quadrupole ---

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void quad_thick (cflw<M> &m, num_t lw, int is)
{                                           (void)is;
  if (fabs(m.k1) < minstr || !m.charge) return strex_drift<M>(m, lw, is);

  mdump(0);
  P    l = R(m.el)*lw;
  int ws = fval(m.k1)*m.edir < 0 ? -1 : 1;

  P cx, sx, mx1, mx2;
  P cy, sy, my1, my2;

  if (fabs(m.k1) >= minstr) {
    P w = sqrt(abs(R(m.k1)))*(ws*m.sdir*m.edir);
    cx = cos (w*l), sx  = sin (w*l);
    cy = cosh(w*l), sy  = sinh(w*l);
    mx1 = sx/w    , mx2 = -sx*w;
    my1 = sy/w    , my2 =  sy*w;
  } else {
    cx = 1., sx = 0., mx1 = l, mx2 = 0.;
    cy = 1., sy = 0., my1 = l, my2 = 0.;
  }

  if (ws != m.charge) // swap x <-> y
    swap(cx,cy), swap(sx,sy), swap(mx1,my1), swap(mx2,my2);

  FOR(i,m.npar) {
    M p(m,i);
    T nx  = p.x*cx  + p.px*mx1;
    T npx = p.x*mx2 + p.px*cx;
    T ny  = p.y*cy  + p.py*my1;
    T npy = p.y*my2 + p.py*cy;

    p.x  = nx;
    p.px = npx;
    p.y  = ny;
    p.py = npy;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void quad_kick (cflw<M> &m, num_t lw, int is)
{                                          (void)is;
  if (fabs(m.k1) < minstr) return strex_kick<M>(m, lw, is);

  P l = R(m.el)*lw;
  num_t dw = is == 0 ? 1./2 : 1.; // drift weight
  if (is >= 0) drift_adj<M>(m, dw*l);

  mdump(0);
  if (m.nmul > 0) {
    num_t wchg = lw*m.edir*m.charge;
    T bx, by;

    FOR (i,m.npar) {
      M p(m,i);
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - R(m.knl[1])*p.x);
      p.py += wchg*(bx - R(m.knl[1])*p.y);
    }
  }
  mdump(1);

  if (is <= 0) drift_adj<M>(m, dw*l);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void quad_thicks (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  if (fabs(m.k1) < minstr || !m.charge) return strex_drift<M>(m, lw, is);

  mdump(0);
  P l   = R(m.el)*lw;
  int ws = fval(m.k1)*m.edir < 0 ? -1 : 1;
  P w   = sqrt(abs(R(m.k1)))*m.sdir*m.edir*ws;
  P cx  = cos (w*l), sx  = sin (w*l);
  P cy  = cosh(w*l), sy  = sinh(w*l);
  P mx1 = sx/w     , mx2 = -sx*w;
  P my1 = sy/w     , my2 =  sy*w;
  R ca  = m.ca, sa = m.sa;

  if (ws != m.charge) // swap x <-> y
    swap(cx,cy), swap(sx,sy), swap(mx1,my1), swap(mx2,my2);

  FOR(i,m.npar) {
    M p(m,i);
    // srotation
    T rx  = ca*p.x  + sa*p.y;
    T rpx = ca*p.px + sa*p.py;
    T ry  = ca*p.y  - sa*p.x;
    T rpy = ca*p.py - sa*p.px;

    T nx  = rx*cx  + rpx*mx1;
    T npx = rx*mx2 + rpx*cx;
    T ny  = ry*cy  + rpy*my1;
    T npy = ry*my2 + rpy*cy;

    // srotation^-1
    p.x   = ca*nx  - sa*ny;
    p.px  = ca*npx - sa*npy;
    p.y   = ca*ny  + sa*nx;
    p.py  = ca*npy + sa*npx;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void quad_kicks (cflw<M> &m, num_t lw, int is)
{                                           (void)is;
  if (fabs(m.k1) < minstr) return strex_kick<M>(m, lw, is);

  P l = R(m.el)*lw;
  num_t dw = is == 0 ? 1./2 : 1.; // drift weight
  if (is >= 0) drift_adj<M>(m, dw*l);

  mdump(0);
  if (m.nmul > 0) {
    num_t wchg = lw*m.edir*m.charge;
    T bx, by;

    FOR (i,m.npar) {
      M p(m,i);
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - R(m.knl[1])*p.x + R(m.ksl[1])*p.y);
      p.py += wchg*(bx - R(m.knl[1])*p.y - R(m.ksl[1])*p.x);
    }
  }
  mdump(1);

  if (is <= 0) drift_adj<M>(m, dw*l);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void quad_thickh (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  mdump(0);
  P l     = R(m.el)*lw;
  P kx    = (R(m.knl[1]) + R(m.eh)*R(m.knl[0]))/R(m.el);
  P ky    = -R(m.knl[1])/R(m.el);
  int wxs = fval(kx)*m.edir < 0 ? -1 : 1;
  int wys = fval(ky)*m.edir < 0 ? -1 : 1;

  P wx, cx, sx, wy, cy, sy;
  P mx11, mx12, mx13, mx21, mx22, mx23, mx31, mx32, mx33;
  P my11, my12,       my21, my22;

  if (fabs(kx) >= minstr) {
    wx = sqrt(abs(kx))*(wxs*m.charge); wxs *= -m.charge;
    if (fval(wx) > 0) cx = cos (wx*l), sx = sin (wx*l);
    else             cx = cosh(wx*l), sx = sinh(wx*l);
    mx11 = cx       , mx12 = sx/wx, mx13 =      R(m.eh) *  (cx-1)*wxs/sqr(wx);
    mx21 = wxs*wx*sx, mx22 = cx   , mx23 =      R(m.eh) *   mx12;
    mx31 = mx23     , mx32 = mx13 , mx33 = -sqr(R(m.eh))*(l-mx12)*wxs/sqr(wx);
  } else {
    wx = 0.;
    mx11 = 1.  , mx12 = l   , mx13 = R(m.eh)*sqr(l)/2;
    mx21 = 0.  , mx22 = 1.  , mx23 = R(m.eh)*l;
    mx31 = mx23, mx32 = mx13, mx33 = mx13*mx23/3;
  }

  if (fabs(ky) >= minstr) {
    wy = sqrt(abs(ky))*(wys*m.charge); wys *= -m.charge;
    if (fval(wy) > 0) cy = cos (wy*l), sy = sin (wy*l);
    else             cy = cosh(wy*l), sy = sinh(wy*l);
    my11 = cy       , my12 = sy/wy;
    my21 = wys*wy*sy, my22 = cy;
  } else {
    wy = 0.;
    my11 = 1., my12 = l ;
    my21 = 0., my22 = 1.;
  }

  FOR (i,m.npar) {
    M p(m,i);
    T nx  = p.x*mx11 + p.px*mx12 + p.pt*(mx13/m.beta);
    T npx = p.x*mx21 + p.px*mx22 + p.pt*(mx23/m.beta);
    T ny  = p.y*my11 + p.py*my12;
    T npy = p.y*my21 + p.py*my22;
    T dt  = p.x*(mx31/m.beta) + p.px*(mx32/m.beta) + p.pt*mx33;

    p.x   = nx;
    p.y   = ny;
    p.px  = npx;
    p.py  = npy;
    p.t  -= dt;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void quad_kickh (cflw<M> &m, num_t lw, int is)
{                                           (void)is;
  P l  = R(m.el)*lw;
  P lh = R(m.eh)*l;
  num_t dw = is == 0 ? 1./2 : 1.; // drift weight

  if (is >= 0) drift_adj<M>(m, dw*l);

  mdump(0);
  if (m.nmul > 0) {
    num_t wchg = lw*m.edir*m.charge;
    T bx, by;

    FOR (i,m.npar) {
      M p(m,i);
      T pz = sqrt(1 + 2/m.beta*p.pt + sqr(p.pt));
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - R(m.knl[1])*p.x) - lh*(pz-(1/m.beta*p.pt));
      p.py += wchg*(bx - R(m.knl[1])*p.y);
      p.t  -= lh*((1/m.beta+p.pt)/pz - 1/m.beta)*p.x;
    }
  }
  mdump(1);

  if (is <= 0) drift_adj<M>(m, dw*l);
}

// --- solenoid ---

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void solen_thick (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  if (fabs(m.ks) < minstr || !m.charge) return strex_drift<M>(m, lw, is);

  mdump(0);
  P l    = R(m.el)*lw;
  P bsol = R(m.ks)*0.5*m.edir*m.charge;

  FOR (i,m.npar) {
    M p(m,i);
    T   xp  = p.px + bsol*p.y;
    T   yp  = p.py - bsol*p.x;
    T  l_pz = l*invsqrt(1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(xp) - sqr(yp));
    T  ang  = bsol*l_pz;

    T ca = cos(ang), sa = sin(ang), sc = sinc(ang);

    T lsc = sc*l_pz;
    T xt  = ca*p.x  + lsc*p.px;
    T pxt = ca*p.px - lsc*p.x *sqr(bsol);
    T yt  = ca*p.y  + lsc*p.py;
    T pyt = ca*p.py - lsc*p.y *sqr(bsol);

    p.x  = ca*xt  + sa*yt;
    p.px = ca*pxt + sa*pyt;
    p.y  = ca*yt  - sa*xt;
    p.py = ca*pyt - sa*pxt;
    p.t  = p.t - (1/m.beta+p.pt)*l_pz + (1-m.T)/m.beta*l;
  }
  mdump(1);
}

// --- eseptum ---

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void esept_thick (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  if (!fval(m.volt) || !m.charge) return strex_drift<M>(m, lw, is);

  mdump(0);
  P l  = R(m.el)*lw;
  P k1 = R(m.volt)*(m.edir*m.charge/m.pc);
  R ca = m.ca, sa = m.sa;

  FOR (i,m.npar) {
    M p(m,i);
    // srotation
    T  nx  = ca*p.x  + sa*p.y;
    T  npx = ca*p.px + sa*p.py;
    T  ny  = ca*p.y  - sa*p.x;
    T  npy = ca*p.py - sa*p.px;

    T   e1 = 1/m.beta+p.pt;
    T   dp = e1 + k1*ny;
    T l_pz = l*invsqrt(sqr(dp) - 1/sqr(m.betgam) - sqr(npx) - sqr(npy));
    T  arg = k1*l_pz;
    T  shx = sinhc(arg)*l_pz;
    T   ch = cosh(arg), sh = sinh(arg);
    T  chm = sqr(sinh(0.5*arg))*(2/k1);
    T   dt = chm*npy + sh *ny  + e1*shx;
    T   yt = ch *ny  + shx*npy + e1*chm;
    T  pyt = ch *npy + sh*dp;

    nx += npx*l_pz;
    ny  = yt, npy = pyt;

    // srotation^-1
    p.x  = ca*nx  - sa*ny;
    p.px = ca*npx - sa*npy;
    p.y  = ca*ny  + sa*nx;
    p.py = ca*npy + sa*npx;
    p.t  = p.t - dt + (1-m.T)/m.beta*l;
  }
  mdump(1);
}

// --- rfcavity ---

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void rfcav_kick (cflw<M> &m, num_t lw, int is)
{                                           (void)is;
  if (!fabs(m.volt) || !m.charge) return;

  mdump(0);
  P w  = R(m.freq)*twopi_clight;
  P vl = R(m.volt)*(lw*m.edir*m.charge/m.pc);

  FOR (i,m.npar) {
    M p(m,i);
    p.pt += vl*sin(R(m.lag) - w*p.t);
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void rfcav_kickn (cflw<M> &m, num_t lw, int is)
{                                            (void)is;
  if ((!m.nmul && !fabs(m.volt)) || !m.charge) return;

  mdump(0);
  num_t wchg = lw*m.edir*m.charge;
  P     w    = R(m.freq)*twopi_clight;
  P     vl   = R(m.volt)*(wchg/m.pc);

  T bx, by, byt;

  FOR (i,m.npar) {
    M p(m,i);
    T ph = R(m.lag) - w*p.t;
    T sa = sin(ph), ca = cos(ph);
    T f; f = 1.;

    if (fval(m.volt)) {
      if (m.nbsl > 0) {
        T df, r2; df = 0., r2 = 1.;

        FOR(i,1,m.nbsl+1) {
          r2  = -r2*(sqr(w)/(4*sqr(i+1)));
          df +=  r2*(i*2);
          r2  =  r2*(sqr(p.x)+sqr(p.y));
          f  +=  r2;
        }

        T c1 = vl/w*df*ca;
        p.px += p.x*c1;
        p.py += p.y*c1;
      }

      p.pt += f*sa*vl;
    }

    if (m.nmul > 0) {
      bxby(m, p.x, p.y, bx, by);

      p.px += wchg/m.pc*by*ca;
      p.py -= wchg/m.pc*bx*ca;

      by = -R(m.knl[m.nmul-1])/m.nmul;
      bx = -R(m.ksl[m.nmul-1])/m.nmul;
      RFOR(i,m.nmul-1) {
        byt = p.x*by - p.y*bx - R(m.knl[i])/(i+1);
        bx  = p.y*by + p.x*bx - R(m.ksl[i])/(i+1);
        by  = byt;
      }
      byt = p.x*by - p.y*bx;
      bx  = p.y*by + p.x*bx;
      by  = byt;

      p.pt -= wchg/m.pc*w*by*sa;
    }
  }
  mdump(1);
}

// --- fringe maps ------------------------------------------------------------o

// must be identical to M.fringe in madl_dynamp.mad
enum {
 fringe_none  = 0,
 fringe_bend  = 1, fringe_mult  = 2 , fringe_qsad = 2+4,
 fringe_rfcav = 8, fringe_solen = 16, fringe_comb = 1+2, fringe_combqs = 1+2+4
};

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void adjust_time (cflw<M> &m, num_t lw)
{                                    (void)lw;
  if (fabs(m.el) < minlen) return;

  P Tl = (m.T-m.Tbak)/m.beta*R(m.el)*m.sdir;

  FOR (i,m.npar) {
    M p(m,i);
    p.t += Tl;
  }
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void cav_fringe (cflw<M> &m, num_t lw)
{
  if (fabs(m.el) < minlen || !fval(m.volt)) return;

  mdump(0);
  P w  = R(m.freq)*twopi_clight;
  P vl = R(m.volt)/R(m.el)*(0.5*lw*m.edir*m.charge/m.pc);

  FOR (i,m.npar) {
    M p(m,i);
    T s1 = sin(R(m.lag) - w*p.t);
    T c1 = cos(R(m.lag) - w*p.t);

    p.px -= vl*s1*p.x;
    p.py -= vl*s1*p.y;
    p.pt += 0.5*vl*w*c1*(sqr(p.x) + sqr(p.y));
  }
  mdump(1);
}


#if 0  // version from L.D. & PTC (requires correct inversion for backtracking)
template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void bend_face (cflw<M> &m, num_t lw, const V &h)
{                                  (void)lw;
  if (!fval(h) || fabs(m.el) < minlen || fabs(m.knl[0]) < minstr) return;

  mdump(0);
  P k0hq = R(m.knl[0])/R(m.el)*(0.5*h*m.sdir*m.edir*m.charge);

  FOR (i,m.npar) {
    M p(m,i);
    if (m.sdir == 1) p.px += k0hq*sqr(p.x);

    T dpp      =  1 + 2/m.beta*p.pt + sqr(p.pt);
    T _pt2     =  1/(dpp - sqr(p.px));
    T xi       =  2*k0hq * sqrt(dpp)*_pt2;
    T dxi_px   =  2*p.px*xi         *_pt2;
    T dxi_ddel = -2     *xi*(1+p.pt)*_pt2;
    T y2       = sqr(p.y);

    p.x  /= 1-dxi_px*y2;
    p.px -= xi*y2;
    p.py -= 2*xi*p.x*p.y;
    p.t  += dxi_ddel*p.x*y2;

    if (m.sdir == -1) p.px += k0hq*sqr(p.x);
  }
  mdump(1);
}
#else  // version from J.Gray, slightly improves backtracking (review needed)
template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void bend_face (cflw<M> &m, num_t lw, const V &h)
{                                  (void)lw;
  if (!fval(h) || fabs(m.el) < minlen || fabs(m.knl[0]) < minstr) return;

  mdump(0);
  P k0hq = R(m.knl[0])/R(m.el)*(0.5*h*m.sdir*m.edir*m.charge);

  FOR (i,m.npar) {
    M p(m,i);
    if (m.sdir == 1) p.px += k0hq*sqr(p.x);

    T dpp    =  1 + 2/m.beta*p.pt + sqr(p.pt);
    T y2     = sqr(p.y);
    T _pt2   =  1/(dpp - sqr(p.px));
    T xi     =  2*k0hq * sqrt(dpp)*_pt2;
    T dxi_px =  2*p.px*xi         *_pt2;
    
    if (m.sdir == -1) {
      T npx  = p.px - xi*y2;
      _pt2   = 1/(dpp - sqr(npx));
      xi     = 2*k0hq * sqrt(dpp)*_pt2;
      dxi_px = 2*npx*xi          *_pt2;
    }

    T dxi_ddel = -2*xi*(1+p.pt)*_pt2;

    if (m.sdir == 1) {
      p.x  /= 1-dxi_px*y2;
      p.px -= xi*y2;
      p.py -= 2*xi*p.x*p.y;
      p.t  += dxi_ddel*p.x*y2;
    } else {
      p.t  += dxi_ddel*p.x*y2;
      p.py -= 2*xi*p.x*p.y;
      p.x  /= 1-dxi_px*y2;
      p.px -= xi*y2 - k0hq*sqr(p.x);
    }
  }
  mdump(1);
}
#endif

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void bend_ptch (cflw<M> &m, num_t lw, const V &a)
{                                  (void)lw;
  if (!fval(a) || !fval(m.elc)) return;

  P dx = R(m.elc)*sin(m.sdir*m.edir*a/2);

  FOR (i,m.npar) {
    M p(m,i);
    p.x += dx;
  }
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void bend_wedge (cflw<M> &m, num_t lw, const V &e)
{                                   (void)lw;
  if (!fval(e)) return;
  if (fabs(m.knl[0]) < minstr) return yrotation<M>(m,m.sdir,-e);

  mdump(0);
  P b1 = R(m.knl[0])/R(m.el)*(m.edir*m.charge);
  P we = e*m.sdir*m.edir;
  P sa = sin(we), ca = cos(we), s2a = sin(2*we);

  FOR (i,m.npar) {
    M p(m,i);
    T pzy = 1 + 2/m.beta*p.pt + sqr(p.pt) - sqr(p.py);
    T _pt = invsqrt(pzy);
    T  pz = sqrt(pzy - sqr(p.px));
    T pzx = pz - b1*p.x;
    T npx = p.px*ca + pzx*sa;
    T pzs = sqrt(pzy - sqr(npx));
    T dxs = (we + asin(p.px*_pt) - asin(npx*_pt))/b1;

    p.x  *= ca + (p.px*s2a + sqr(sa)*(pz+pzx))/(pzs + pz*ca - p.px*sa);
    p.px  = npx;
    p.y  += dxs*p.py;
    p.t  -= dxs*(1/m.beta+p.pt);
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R, typename V>
inline void mad8_wedge (cflw<M> &m, num_t lw, const V &e)
{                                   (void)lw;
  if (!fval(e) || fabs(m.knl[1]) < minstr) return;

  mdump(0);
  num_t wc = m.frng <= 1 ? 0 : 0.25;
  P    k1e = R(m.knl[1])/R(m.el)*e;
  P     c1 = (1+wc)*m.charge*k1e;
  P     c2 = (1-wc)*m.charge*k1e;

  FOR (i,m.npar) {
    M p(m,i);
    p.px +=   c1*sqr(p.x) - c2*sqr(p.y);
    p.py -= 2*c2*p.x*p.y;
  }
  mdump(1);
}

#if 0  // version from L.D. & PTC (requires correct inversion for backtracking)
template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void bend_fringe (cflw<M> &m, num_t lw)
{
  if (fabs(m.knl[0]) < minang) return;

  mdump(0);
  P   fh = 2*R(m.fint)*R(m.hgap);
  P   b0 = R(m.knl[0])/abs(R(m.el))*(lw*m.sdir*m.edir*m.charge);
  P   c2 = fh*b0;
  P fsad;

  if (fval(fh)) fsad=1/(36*fh); else fsad = 0.;

  FOR (i,m.npar) {
    M p(m,i);
    T   dpp = 1 + 2/m.beta*p.pt + sqr(p.pt);
    T    pz = sqrt(dpp - sqr(p.px) - sqr(p.py));
    T   _pz = 1/pz;
    T  _pz2 = sqr(_pz);
    T  relp = invsqrt(dpp);
    T  tfac = -1/m.beta - p.pt;
    T    c3 = fsad*sqr(b0)*relp;

    T xp  = p.px/pz,  yp  = p.py/pz;
    T xyp = xp*yp  ,  yp2 = 1+sqr(yp);
    T xp2 = sqr(xp), _yp2 = 1/yp2;

    T fi0 = atan((xp*_yp2)) - c2*(1 + xp2*(1+yp2))*pz;
    T co2 = b0/sqr(cos(fi0));
    T co1 = co2/(1 + sqr(xp*_yp2))*_yp2;
    T co3 = co2*c2;

    T fi1 =    co1          - co3*2*xp*(1+yp2)*pz;
    T fi2 = -2*co1*xyp*_yp2 - co3*2*xp*xyp    *pz;
    T fi3 =                 - co3*(1 + xp2*(1+yp2));

    T kx = fi1*(1+xp2)*_pz  + fi2*xyp*_pz      - fi3*xp;
    T ky = fi1*xyp*_pz      + fi2*yp2*_pz      - fi3*yp;
    T kz = fi1*tfac*xp*_pz2 + fi2*tfac*yp*_pz2 - fi3*tfac*_pz;

    if (m.sdir == 1) {
      T y  = 2*p.y / (1 + sqrt(1-2*ky*p.y));
      T y2 = sqr(y);

      p.x  += 0.5*kx*y2;
      p.py -= (4*c3*y2 + b0*tan(fi0))*y;
      p.t  += (0.5*kz  + c3*y2*sqr(relp)*tfac)*y2;
      p.y   = y;
    } else { // need to reverse y-dependence
      T y2 = sqr(p.y);

      p.x  -= 0.5*kx*y2;
      p.py += (4*c3*y2 + b0*tan(fi0))*p.y;
      p.t  -= (0.5*kz  + c3*y2*sqr(relp)*tfac)*y2;
      p.y  *= 0.5*(1 + sqrt(1-2*ky*p.y));
    }
  }
  mdump(1);
}
#else  // version from J.Gray, slightly improves backtracking (review needed)
template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void bend_fringe (cflw<M> &m, num_t lw)
{
  if (fabs(m.knl[0]) < minang) return;

  mdump(0);
  P   fh = 2*R(m.fint)*R(m.hgap);
  P   b0 = R(m.knl[0])/abs(R(m.el))*(lw*m.sdir*m.edir*m.charge);
  P   c2 = fh*b0;
  P fsad;

  if (fval(fh)) fsad=1/(36*fh); else fsad = 0.;

  FOR (i,m.npar) {
    M p(m,i);
    T   dpp = 1 + 2/m.beta*p.pt + sqr(p.pt);
    T    pz = sqrt(dpp - sqr(p.px) - sqr(p.py));
    T   _pz = 1/pz;
    T  _pz2 = sqr(_pz);
    T  relp = invsqrt(dpp);
    T  tfac = -1/m.beta - p.pt;
    T    c3 = fsad*sqr(b0)*relp;

    T xp  = p.px/pz,  yp  = p.py/pz;
    T xyp = xp*yp  ,  yp2 = 1+sqr(yp);
    T xp2 = sqr(xp), _yp2 = 1/yp2;

    T fi0 = atan((xp*_yp2)) - c2*(1 + xp2*(1+yp2))*pz;
    T co2 = b0/sqr(cos(fi0));
    T co1 = co2/(1 + sqr(xp*_yp2))*_yp2;
    T co3 = co2*c2;

    T fi1 =    co1          - co3*2*xp*(1+yp2)*pz;
    T fi2 = -2*co1*xyp*_yp2 - co3*2*xp*xyp    *pz;
    T fi3 =                 - co3*(1 + xp2*(1+yp2));

    T ky = fi1*xyp*_pz      + fi2*yp2*_pz      - fi3*yp;
    T y  = 2*p.y / (1 + sqrt(1-2*ky*p.y));

    if (m.sdir == 1) {
      // Only calculate kx and kz once
      T kx = fi1*(1+xp2)*_pz  + fi2*xyp*_pz      - fi3*xp;
      T kz = fi1*tfac*xp*_pz2 + fi2*tfac*yp*_pz2 - fi3*tfac*_pz;
      T y2 = sqr(y);

      p.x  += 0.5*kx*y2;
      p.py -= (4*c3*y2 + b0*tan(fi0))*y;
      p.t  += (0.5*kz  + c3*y2*sqr(relp)*tfac)*y2;
      p.y   = y;
    } else { // need to reverse y-dependence
      T y2 = sqr(p.y);
      // recalculate kx, ky and kz with a new py (Work out what it probably was)
      T   npy = p.py + (4*c3*y2 + b0*tan(fi0))*y;
      
        pz = sqrt(dpp - sqr(p.px) - sqr(npy));
       _pz = 1/pz;
      _pz2 = sqr(_pz);
      
      xp  = p.px/pz,  yp  = npy/pz;
      xyp = xp*yp  ,  yp2 = 1+sqr(yp);
      xp2 = sqr(xp), _yp2 = 1/yp2;

      fi0 = atan((xp*_yp2)) - c2*(1 + xp2*(1+yp2))*pz;
      co2 = b0/sqr(cos(fi0));
      co1 = co2/(1 + sqr(xp*_yp2))*_yp2;
      co3 = co2*c2;

      fi1 =    co1          - co3*2*xp*(1+yp2)*pz;
      fi2 = -2*co1*xyp*_yp2 - co3*2*xp*xyp    *pz;
      fi3 =                 - co3*(1 + xp2*(1+yp2));

      T kx = fi1*(1+xp2)*_pz  + fi2*xyp*_pz      - fi3*xp;
        ky = fi1*xyp*_pz      + fi2*yp2*_pz      - fi3*yp;
      T kz = fi1*tfac*xp*_pz2 + fi2*tfac*yp*_pz2 - fi3*tfac*_pz;

      p.x  -= 0.5*kx*y2;
      p.py += (4*c3*y2 + b0*tan(fi0))*p.y;
      p.t  -= (0.5*kz  + c3*y2*sqr(relp)*tfac)*y2;
      p.y  -= y2*ky/2;
    }
  }
  mdump(1);
}
#endif

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void qsad_fringe (cflw<M> &m, num_t lw)
{
  if (fabs(m.knl[1])+fabs(m.ksl[1]) < minstr) return;
  if (fabs(m.f1)    +fabs(m.f2)     < minstr) return;

  mdump(0);
  P     a    = -0.5*atan2(R(m.ksl[1]), R(m.knl[1]));
  P     b2   = hypot(R(m.knl[1]), R(m.ksl[1]))/R(m.el)*m.edir;
  P     ca   = cos(a), sa = sin(a);
  P     bf1  = (abs(R(m.f1))*R(m.f1)/-24)*b2;
  P     bf2  =               R(m.f2)     *b2*m.sdir;

  // Lee-Whiting formula, E. Forest ch 13.2.3, eq 13.33
  FOR (i,m.npar) {
    M p(m,i);
    T _pz = invsqrt(1 + 2/m.beta*p.pt + sqr(p.pt));
    T  dt = (1/m.beta+p.pt)*_pz;

    T  f1 = lw*m.charge*bf1*_pz;
    T  f2 =    m.charge*bf2*_pz;

    T nx  = ca*p.x  + sa*p.y;
    T npx = ca*p.px + sa*p.py;
    T ny  = ca*p.y  - sa*p.x;
    T npy = ca*p.py - sa*p.px;

    p.t  -= dt*((f1*nx + (1+f1/2)*exp(-f1)*f2*npx*_pz)*npx -
                (f1*ny + (1-f1/2)*exp( f1)*f2*npy*_pz)*npy)*_pz;

    nx  =  nx*exp( f1) + f2*npx*_pz;
    ny  =  ny*exp(-f1) - f2*npy*_pz;
    npx = npx*exp(-f1);
    npy = npy*exp( f1);

    p.x  = ca*nx  - sa*ny;
    p.px = ca*npx - sa*npy;
    p.y  = ca*ny  + sa*nx;
    p.py = ca*npy + sa*npx;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void mult_fringe (cflw<M> &m, num_t lw)
{
  int n = MIN(m.nmul,m.fmax);
  if (!n) return;

  mdump(0);
  num_t  wchg = lw*m.charge;
  int   no_k1 = m.frng & fringe_bend;
  P        _l;

  if (fval(m.el)) _l = m.edir/R(m.el); else _l = m.edir;

  FOR (i,m.npar) {
    M p(m,i);
    T rx, ix;       rx = 1., ix =0.;
    T fx, fxx, fxy; fx = 0., fxx=0., fxy=0.;
    T fy, fyy, fyx; fy = 0., fyy=0., fyx=0.;

    T _pz = invsqrt(1 + 2/m.beta*p.pt + sqr(p.pt));

    FOR (j,1,n+1) {
      T drx, dix; drx = rx, dix = ix;
      rx = drx*p.x - dix*p.y;
      ix = drx*p.y + dix*p.x;

      num_t nj = -wchg/(4*(j+1)), nf = (j+2.)/j;
      P kj = R(m.knl[j-1])*_l, ksj = R(m.ksl[j-1])*_l;

      T u, v, du, dv;
      if (j == 1 && no_k1) {
        u  = nj*(       - ksj*ix );
        v  = nj*(       + ksj*rx );
        du = nj*(       - ksj*dix);
        dv = nj*(       + ksj*drx);
      } else {
        u  = nj*(kj*rx  - ksj*ix );
        v  = nj*(kj*ix  + ksj*rx );
        du = nj*(kj*drx - ksj*dix);
        dv = nj*(kj*dix + ksj*drx);
      }

      T dux =  j*du, dvx = j*dv;
      T duy = -j*dv, dvy = j*du;

      fx  +=   u*p.x + nf*   v*p.y;
      fy  +=   u*p.y - nf*   v*p.x;
      fxx += dux*p.x + nf* dvx*p.y + u;
      fyy += duy*p.y - nf* dvy*p.x + u;
      fxy += duy*p.x + nf*(dvy*p.y + v);
      fyx += dux*p.y - nf*(dvx*p.x + v);
    }

    T    a = 1 - fxx*_pz;
    T    b =   - fyx*_pz;
    T    c =   - fxy*_pz;
    T    d = 1 - fyy*_pz;
    T _det = 1/(a*d - b*c);
    T npx = (d*p.px - b*p.py)*_det; // Create variables so they do not affect each other
    T npy = (a*p.py - c*p.px)*_det;

    p.x  -= fx*_pz;
    p.y  -= fy*_pz;
    p.px  = npx;
    p.py  = npy;
    p.t  += (1/m.beta+p.pt)*(p.px*fx + p.py*fy)*sqr(_pz)*_pz;
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void curex_fringe (cflw<M> &m, num_t lw)
{
  mdump(0);
  if (m.sdir*lw == 1) { // 'forward entry' or 'backward exit'
    yrotation<M> (m, m.sdir,-R(m.e));
    bend_face<M> (m, 1     , R(m.h));
    if (m.frng) {
      if (m.frng & fringe_bend) bend_fringe<M>(m, 1);
      if (m.frng & fringe_mult) mult_fringe<M>(m, 1);
      if (m.frng & fringe_qsad) qsad_fringe<M>(m, 1);
    }
    mad8_wedge<M>(m, 1, R(m.e));
    bend_wedge<M>(m, 1,-R(m.e));
  } else {
    bend_wedge<M>(m,-1,-R(m.e));
    mad8_wedge<M>(m,-1, R(m.e));
    if (m.frng) {
      if (m.frng & fringe_qsad) qsad_fringe<M>(m,-1);
      if (m.frng & fringe_mult) mult_fringe<M>(m,-1);
      if (m.frng & fringe_bend) bend_fringe<M>(m,-1);
    }
    bend_face<M> (m,-1     , R(m.h));
    yrotation<M> (m,-m.sdir, R(m.e));
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void strex_fringe (cflw<M> &m, num_t lw)
{
  mdump(0);
  bool p =  m.pdir == lw;
  P    a = !m.pdir*(0.5*R(m.eh)*(fval(m.eld) ? R(m.eld) : R(m.el)) - R(m.e));

  if (m.sdir*lw == 1) { // 'forward entry' or 'backward exit'
    yrotation<M> (m, m.sdir,-R(m.e)  );
    bend_ptch<M> (m, 1     , R(m.a)*p);
    bend_face<M> (m, 1     , R(m.h)  );
    if (m.frng) {
      if (m.frng & fringe_bend) bend_fringe<M>(m, 1);
      if (m.frng & fringe_mult) mult_fringe<M>(m, 1);
      if (m.frng & fringe_qsad) qsad_fringe<M>(m, 1);
    }
    bend_wedge<M>(m, 1, a);
  } else {
    bend_wedge<M>(m,-1, a);
    if (m.frng) {
      if (m.frng & fringe_qsad) qsad_fringe<M>(m,-1);
      if (m.frng & fringe_mult) mult_fringe<M>(m,-1);
      if (m.frng & fringe_bend) bend_fringe<M>(m,-1);
    }
    bend_face<M> (m,-1     , R(m.h)  );
    bend_ptch<M> (m,-1     , R(m.a)*p);
    yrotation<M> (m,-m.sdir, R(m.e)  );
  }
  mdump(1);
}

template <typename M, typename T=M::T, typename P=M::P, typename R=M::R>
inline void rfcav_fringe (cflw<M> &m, num_t lw)
{
  mdump(0);
  if (m.sdir*lw == 1) { // 'forward entry' or 'backward exit'
    ensure(m.Tbak < 0, "inconsistent totalpath when entering rfcavity");
    if (m.Tbak == -1)
      m.Tbak = m.T, m.T = 1, R(m.lag) -= R(m.freq)*abs(R(m.el))*(pi_clight/m.beta);
    else
      m.Tbak = m.T, m.T = 0;
    if (lw == -1 && m.T != m.Tbak) adjust_time<M>(m, 1);
    if (m.frng & fringe_rfcav)      cav_fringe<M>(m, 1);
  } else {
    if (m.frng & fringe_rfcav)      cav_fringe<M>(m,-1);
    if (lw == -1 && m.T != m.Tbak) adjust_time<M>(m,-1);
    m.T = m.Tbak, m.Tbak = -1;
  }
  mdump(1);
}

// --- specializations --------------------------------------------------------o

// --- tilt & misalignment ---
void mad_trk_tilt_r (mflw_t *m, num_t lw) {
  srotation<par_t>(m->rflw, lw*m->rflw.sdir, m->rflw.tlt);
}
void mad_trk_tilt_t (mflw_t *m, num_t lw) {
  srotation<map_t>(m->tflw, lw*m->tflw.sdir, m->tflw.tlt);
}
void mad_trk_tilt_p (mflw_t *m, num_t lw) {
  srotation<prm_t>(m->pflw, lw*m->pflw.sdir, tpsa_ref(m->pflw.tlt));
}

void mad_trk_misalign_r (mflw_t *m, num_t lw) {
  misalign<par_t>(m->rflw, lw);
}
void mad_trk_misalign_t (mflw_t *m, num_t lw) {
  misalign<map_t>(m->tflw, lw);
}
void mad_trk_misalign_p (mflw_t *m, num_t lw) {
  misalign<prm_t>(m->pflw, lw);
}

// -- fringe maps ---

void mad_trk_strex_fringe_r (mflw_t *m, num_t lw) {
  strex_fringe<par_t>(m->rflw, lw);
}
void mad_trk_curex_fringe_r (mflw_t *m, num_t lw) {
  curex_fringe<par_t>(m->rflw, lw);
}
void mad_trk_rfcav_fringe_r (mflw_t *m, num_t lw) {
  rfcav_fringe<par_t>(m->rflw, lw);
}

void mad_trk_strex_fringe_t (mflw_t *m, num_t lw) {
  strex_fringe<map_t>(m->tflw, lw);
}
void mad_trk_curex_fringe_t (mflw_t *m, num_t lw) {
  curex_fringe<map_t>(m->tflw, lw);
}
void mad_trk_rfcav_fringe_t (mflw_t *m, num_t lw) {
  rfcav_fringe<map_t>(m->tflw, lw);
}

void mad_trk_strex_fringe_p (mflw_t *m, num_t lw) {
  strex_fringe<prm_t>(m->pflw, lw);
}
void mad_trk_curex_fringe_p (mflw_t *m, num_t lw) {
  curex_fringe<prm_t>(m->pflw, lw);
}
void mad_trk_rfcav_fringe_p (mflw_t *m, num_t lw) {
  rfcav_fringe<prm_t>(m->pflw, lw);
}

// --- patches ---

void mad_trk_xrotation_r (mflw_t *m, num_t lw, int is) {
  xrotation<par_t>(m->rflw, lw, zero); (void)is;
}
void mad_trk_yrotation_r (mflw_t *m, num_t lw, int is) {
  yrotation<par_t>(m->rflw, lw, zero); (void)is;
}
void mad_trk_srotation_r (mflw_t *m, num_t lw, int is) {
  srotation<par_t>(m->rflw, lw, zero); (void)is;
}
void mad_trk_translate_r (mflw_t *m, num_t lw, int is) {
  translate<par_t>(m->rflw, lw, zero, zero, zero); (void)is;
}
void mad_trk_changeref_r (mflw_t *m, num_t lw, int is) {
  changeref<par_t>(m->rflw, lw); (void)is;
}

void mad_trk_xrotation_t (mflw_t *m, num_t lw, int is) {
  xrotation<map_t>(m->tflw, lw, zero); (void)is;
}
void mad_trk_yrotation_t (mflw_t *m, num_t lw, int is) {
  yrotation<map_t>(m->tflw, lw, zero); (void)is;
}
void mad_trk_srotation_t (mflw_t *m, num_t lw, int is) {
  srotation<map_t>(m->tflw, lw, zero); (void)is;
}
void mad_trk_translate_t (mflw_t *m, num_t lw, int is) {
  translate<map_t>(m->tflw, lw, zero, zero, zero); (void)is;
}
void mad_trk_changeref_t (mflw_t *m, num_t lw, int is) {
  changeref<map_t>(m->tflw, lw); (void)is;
}

void mad_trk_xrotation_p (mflw_t *m, num_t lw, int is) {
  xrotation<prm_t>(m->pflw, lw, zero); (void)is;
}
void mad_trk_yrotation_p (mflw_t *m, num_t lw, int is) {
  yrotation<prm_t>(m->pflw, lw, zero); (void)is;
}
void mad_trk_srotation_p (mflw_t *m, num_t lw, int is) {
  srotation<prm_t>(m->pflw, lw, zero); (void)is;
}
void mad_trk_translate_p (mflw_t *m, num_t lw, int is) {
  translate<prm_t>(m->pflw, lw, zero, zero, zero); (void)is;
}
void mad_trk_changeref_p (mflw_t *m, num_t lw, int is) {
  changeref<prm_t>(m->pflw, lw); (void)is;
}

// --- DKD straight ---

void mad_trk_strex_drift_r (mflw_t *m, num_t lw, int is) {
  strex_drift<par_t>(m->rflw,lw,is);
}
void mad_trk_strex_kick_r (mflw_t *m, num_t lw, int is) {
  strex_kick<par_t>(m->rflw,lw,is);
}
void mad_trk_strex_kickhs_r (mflw_t *m, num_t lw, int is) {
  strex_kickhs<par_t>(m->rflw,lw,is);
}

void mad_trk_strex_drift_t (mflw_t *m, num_t lw, int is) {
  strex_drift<map_t>(m->tflw,lw,is);
}
void mad_trk_strex_kick_t (mflw_t *m, num_t lw, int is) {
  strex_kick<map_t>(m->tflw,lw,is);
}
void mad_trk_strex_kickhs_t (mflw_t *m, num_t lw, int is) {
  strex_kickhs<map_t>(m->tflw,lw,is);
}

void mad_trk_strex_drift_p (mflw_t *m, num_t lw, int is) {
  strex_drift<prm_t>(m->pflw,lw,is);
}
void mad_trk_strex_kick_p (mflw_t *m, num_t lw, int is) {
  strex_kick<prm_t>(m->pflw,lw,is);
}
void mad_trk_strex_kickhs_p (mflw_t *m, num_t lw, int is) {
  strex_kickhs<prm_t>(m->pflw,lw,is);
}

// --- DKD curved ---

void mad_trk_curex_drift_r (mflw_t *m, num_t lw, int is) {
  curex_drift<par_t>(m->rflw,lw,is);
}
void mad_trk_curex_kick_r (mflw_t *m, num_t lw, int is) {
  curex_kick<par_t>(m->rflw,lw,is);
}

void mad_trk_curex_drift_t (mflw_t *m, num_t lw, int is) {
  curex_drift<map_t>(m->tflw,lw,is);
}
void mad_trk_curex_kick_t (mflw_t *m, num_t lw, int is) {
  curex_kick<map_t>(m->tflw,lw,is);
}

void mad_trk_curex_drift_p (mflw_t *m, num_t lw, int is) {
  curex_drift<prm_t>(m->pflw,lw,is);
}
void mad_trk_curex_kick_p (mflw_t *m, num_t lw, int is) {
  curex_kick<prm_t>(m->pflw,lw,is);
}

// --- sbend ---

void mad_trk_sbend_thick_r (mflw_t *m, num_t lw, int is) {
  sbend_thick<par_t>(m->rflw,lw,is);
}
void mad_trk_sbend_kick_r (mflw_t *m, num_t lw, int is) {
  curex_kick<par_t>(m->rflw,lw,is,true);
}

void mad_trk_sbend_thick_t (mflw_t *m, num_t lw, int is) {
  sbend_thick<map_t>(m->tflw,lw,is);
}
void mad_trk_sbend_kick_t (mflw_t *m, num_t lw, int is) {
  curex_kick<map_t>(m->tflw,lw,is,true);
}

void mad_trk_sbend_thick_p (mflw_t *m, num_t lw, int is) {
  sbend_thick<prm_t>(m->pflw,lw,is);
}
void mad_trk_sbend_kick_p (mflw_t *m, num_t lw, int is) {
  curex_kick<prm_t>(m->pflw,lw,is,true);
}

// --- rbend ---

void mad_trk_rbend_thick_r (mflw_t *m, num_t lw, int is) {
  rbend_thick<par_t>(m->rflw,lw,is);
}
void mad_trk_rbend_kick_r (mflw_t *m, num_t lw, int is) {
  strex_kick<par_t>(m->rflw,lw,is,true);
}

void mad_trk_rbend_thick_t (mflw_t *m, num_t lw, int is) {
  rbend_thick<map_t>(m->tflw,lw,is);
}
void mad_trk_rbend_kick_t (mflw_t *m, num_t lw, int is) {
  strex_kick<map_t>(m->tflw,lw,is,true);
}

void mad_trk_rbend_thick_p (mflw_t *m, num_t lw, int is) {
  rbend_thick<prm_t>(m->pflw,lw,is);
}
void mad_trk_rbend_kick_p (mflw_t *m, num_t lw, int is) {
  strex_kick<prm_t>(m->pflw,lw,is,true);
}

// --- quadrupole ---

void mad_trk_quad_thick_r (mflw_t *m, num_t lw, int is) {
  quad_thick<par_t>(m->rflw,lw,is);
}
void mad_trk_quad_thicks_r (mflw_t *m, num_t lw, int is) {
  quad_thicks<par_t>(m->rflw,lw,is);
}
void mad_trk_quad_thickh_r (mflw_t *m, num_t lw, int is) {
  quad_thickh<par_t>(m->rflw,lw,is);
}
void mad_trk_quad_kick_r (mflw_t *m, num_t lw, int is) {
  quad_kick<par_t>(m->rflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kicks_r (mflw_t *m, num_t lw, int is) {
  quad_kicks<par_t>(m->rflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kickh_r (mflw_t *m, num_t lw, int is) {
  quad_kickh<par_t>(m->rflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kick__r (mflw_t *m, num_t lw, int is) {
  quad_kick<par_t>(m->rflw,lw,is);
}
void mad_trk_quad_kicks__r (mflw_t *m, num_t lw, int is) {
  quad_kicks<par_t>(m->rflw,lw,is);
}
void mad_trk_quad_kickh__r (mflw_t *m, num_t lw, int is) {
  quad_kickh<par_t>(m->rflw,lw,is);
}

void mad_trk_quad_thick_t (mflw_t *m, num_t lw, int is) {
  quad_thick<map_t>(m->tflw,lw,is);
}
void mad_trk_quad_thicks_t (mflw_t *m, num_t lw, int is) {
  quad_thicks<map_t>(m->tflw,lw,is);
}
void mad_trk_quad_thickh_t (mflw_t *m, num_t lw, int is) {
  quad_thickh<map_t>(m->tflw,lw,is);
}
void mad_trk_quad_kick_t (mflw_t *m, num_t lw, int is) {
  quad_kick<map_t>(m->tflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kicks_t (mflw_t *m, num_t lw, int is) {
  quad_kicks<map_t>(m->tflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kickh_t (mflw_t *m, num_t lw, int is) {
  quad_kickh<map_t>(m->tflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kick__t (mflw_t *m, num_t lw, int is) {
  quad_kick<map_t>(m->tflw,lw,is);
}
void mad_trk_quad_kicks__t (mflw_t *m, num_t lw, int is) {
  quad_kicks<map_t>(m->tflw,lw,is);
}
void mad_trk_quad_kickh__t (mflw_t *m, num_t lw, int is) {
  quad_kickh<map_t>(m->tflw,lw,is);
}

void mad_trk_quad_thick_p (mflw_t *m, num_t lw, int is) {
  quad_thick<prm_t>(m->pflw,lw,is);
}
void mad_trk_quad_thicks_p (mflw_t *m, num_t lw, int is) {
  quad_thicks<prm_t>(m->pflw,lw,is);
}
void mad_trk_quad_thickh_p (mflw_t *m, num_t lw, int is) {
  quad_thickh<prm_t>(m->pflw,lw,is);
}
void mad_trk_quad_kick_p (mflw_t *m, num_t lw, int is) {
  quad_kick<prm_t>(m->pflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kicks_p (mflw_t *m, num_t lw, int is) {
  quad_kicks<prm_t>(m->pflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kickh_p (mflw_t *m, num_t lw, int is) {
  quad_kickh<prm_t>(m->pflw,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kick__p (mflw_t *m, num_t lw, int is) {
  quad_kick<prm_t>(m->pflw,lw,is);
}
void mad_trk_quad_kicks__p (mflw_t *m, num_t lw, int is) {
  quad_kicks<prm_t>(m->pflw,lw,is);
}
void mad_trk_quad_kickh__p (mflw_t *m, num_t lw, int is) {
  quad_kickh<prm_t>(m->pflw,lw,is);
}

// --- solenoid ---

void mad_trk_solen_thick_r (mflw_t *m, num_t lw, int is) {
  solen_thick<par_t>(m->rflw,lw,is);
}
void mad_trk_solen_thick_t (mflw_t *m, num_t lw, int is) {
  solen_thick<map_t>(m->tflw,lw,is);
}
void mad_trk_solen_thick_p (mflw_t *m, num_t lw, int is) {
  solen_thick<prm_t>(m->pflw,lw,is);
}

// --- eseptum ---

void mad_trk_esept_thick_r (mflw_t *m, num_t lw, int is) {
  esept_thick<par_t>(m->rflw,lw,is);
}
void mad_trk_esept_thick_t (mflw_t *m, num_t lw, int is) {
  esept_thick<map_t>(m->tflw,lw,is);
}
void mad_trk_esept_thick_p (mflw_t *m, num_t lw, int is) {
  esept_thick<prm_t>(m->pflw,lw,is);
}

// --- rfcavity ---

void mad_trk_rfcav_kick_r (mflw_t *m, num_t lw, int is) {
  rfcav_kick<par_t>(m->rflw,lw,is);
}
void mad_trk_rfcav_kickn_r (mflw_t *m, num_t lw, int is) {
  rfcav_kickn<par_t>(m->rflw,lw,is);
}

void mad_trk_rfcav_kick_t (mflw_t *m, num_t lw, int is) {
  rfcav_kick<map_t>(m->tflw,lw,is);
}
void mad_trk_rfcav_kickn_t (mflw_t *m, num_t lw, int is) {
  rfcav_kickn<map_t>(m->tflw,lw,is);
}

void mad_trk_rfcav_kick_p (mflw_t *m, num_t lw, int is) {
  rfcav_kick<prm_t>(m->pflw,lw,is);
}
void mad_trk_rfcav_kickn_p (mflw_t *m, num_t lw, int is) {
  rfcav_kickn<prm_t>(m->pflw,lw,is);
}

// --- do nothing ---

void mad_trk_fnil (mflw_t *m, num_t lw, int is) {
  (void)m, (void)lw, (void)is;
}

// --- track one thick or thin slice ------------------------------------------o

void mad_trk_slice_one (mflw_t *m, num_t lw, trkfun *fun)
{
  fun(m, lw, zero);
}

// --- track one Yoshida slice ------------------------------------------------o

const ssz_t yosh2_n   = 1;
const num_t yosh2_d[] = {0.5};
const num_t yosh2_k[] = {1};
const ssz_t yosh4_n   = 2;
const num_t yosh4_d[] = { 0x1.59e8b6eb96339p-1,-0x1.67a2dbae58ce4p-3 };
const num_t yosh4_k[] = { 0x1.59e8b6eb96339p+0,-0x1.b3d16dd72c672p+0 };
const ssz_t yosh6_n   = 4;
const num_t yosh6_d[] = { 0x1.91abc4988937bp-2, 0x1.052468fb75c74p-1,
                         -0x1.e25bd194051b9p-2, 0x1.199cec1241558p-4 };
const num_t yosh6_k[] = { 0x1.91abc4988937bp-1, 0x1.e2743579895b4p-3,
                         -0x1.2d7c6f7933b93p+0, 0x1.50b00cfb7be3ep+0 };
const ssz_t yosh8_n   = 8;
const num_t yosh8_d[] = { 0x1.d466770cfb237p-2, 0x1.2b25476e416dap-1,
                         -0x1.30efca291a66ep-1,-0x1.9a644b62ac4e7p-1,
                          0x1.c7a76da161edap-1,-0x1.702a9ae94c280p-7,
                         -0x1.db997617a90dfp-1, 0x1.cfae4578f406ep-1 };
const num_t yosh8_k[] = { 0x1.d466770cfb237p-1, 0x1.03c82f9f0f6fbp-2,
                         -0x1.71e1d610de42dp+0,-0x1.4413aa8e705cep-3,
                          0x1.f029e2f32ff94p+0,-0x1.f5ea8d5ed529ep+0,
                          0x1.a5117472c1becp-4, 0x1.b55d2e31c7eafp+0 };

const struct {
  const ssz_t  n;
  const num_t *d;
  const num_t *k;
} yosh[] = {
  {yosh2_n, yosh2_d, yosh2_k},   // yosh2 -> 0..0, ord=2, j=0, n=1, k=-1
  {yosh4_n, yosh4_d, yosh4_k},   // yosh4 -> 0..1, ord=4, j=1, n=2, k=-3
  {yosh6_n, yosh6_d, yosh6_k},   // yosh6 -> 0..3, ord=6, j=2, n=4, k=-7
  {yosh8_n, yosh8_d, yosh8_k},   // yosh8 -> 0..7, ord=8, j=3, n=8, k=-15
};

void mad_trk_slice_dkd (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick, int ord)
{
  ensure(ord >= 2 && ord <= 8, "invalid dkd/tkt order 2..8");
  int j = ord/2-1;
  int n = 1<<j;
  int k = -2*n;
  int d = m->rflw.sdir; // this is above polymorphic data, hence same for all

  FOR(i,n) {
    thick(m, lw*yosh[j].d[i  ], ++k*d);
     kick(m, lw*yosh[j].k[i  ], ++k*d);
  } thick(m, lw*yosh[j].d[--n], ++k*d);
  RFOR(i,n) {
     kick(m, lw*yosh[j].k[i  ], ++k*d);
    thick(m, lw*yosh[j].d[i  ], ++k*d);
  }
}

void mad_trk_slice_tkt (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick, int n)
{
  mad_trk_slice_dkd(m, lw, thick, kick, n);
}

// --- track one Boole slice --------------------------------------------------o

const ssz_t boole2_n    = 1;
const num_t boole2_d    = 1.;
const num_t boole2_k[]  = {1./2};
const ssz_t boole4_n    = 2;
const num_t boole4_d    = 1./2;
const num_t boole4_k[]  = {1./6, 4./6};
const ssz_t boole6_n    = 3;
const num_t boole6_d    = 1./4;
const num_t boole6_k[]  = {7./90, 32./90, 12./90};
const ssz_t boole8_n    = 4;
const num_t boole8_d    = 1./6;
const num_t boole8_k[]  = {41./840, 216./840, 27./840, 272./840};
const ssz_t boole10_n   = 5;
const num_t boole10_d   = 1./8;
const num_t boole10_k[] = {  989./28350, 5888./28350, -928./28350, 10496./28350,
                           -4540./28350};
const ssz_t boole12_n   = 6;
const num_t boole12_d   = 1./10;
const num_t boole12_k[] = { 16067./598752,  106300./598752, -48525./598752,
                           272400./598752, -260550./598752, 427368./598752};

const struct {
  const ssz_t  n;
  const num_t  d;
  const num_t *k;
} boole[] = {
  {boole2_n , boole2_d , boole2_k }, // bool2  -> 0..0, ord=2,  j=0, n=0, k=-1
  {boole4_n , boole4_d , boole4_k }, // bool4  -> 0..1, ord=4,  j=1, n=1, k=-2
  {boole6_n , boole6_d , boole6_k }, // bool6  -> 0..2, ord=6,  j=2, n=2, k=-4
  {boole8_n , boole8_d , boole8_k }, // bool8  -> 0..3, ord=8,  j=3, n=3, k=-6
  {boole10_n, boole10_d, boole10_k}, // bool10 -> 0..4, ord=10, j=4, n=4, k=-8
  {boole12_n, boole12_d, boole12_k}, // bool12 -> 0..5, ord=12, j=5, n=5, k=-10
} ;

void mad_trk_slice_kmk (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick, int ord)
{
  ensure(ord >= 2 && ord <= 12, "invalid kmk order 2..12");
  int j = ord/2-1;
  int n = j;
  int k = -2*n;                      if (!k) --k;
  int d = m->rflw.sdir; // this is above polymorphic data, hence same for all

  FOR(i,n) {
     kick(m, lw*boole[j].k[i], k++*d);
    thick(m, lw*boole[j].d   , k++*d);
  }  kick(m, lw*boole[j].k[n], k++*d); if (!n) ++n;
  RFOR(i,n) {
    thick(m, lw*boole[j].d   , k++*d);
     kick(m, lw*boole[j].k[i], k++*d);
  }
}

// --- track one Teapot slice -------------------------------------------------o

const struct {
  const num_t d;
  const num_t D;
  const num_t k;
} teapot[] = {
  {1./6 , 4./6 , 1./2 }, // teapot2 -> 0..0, knd=2, j=0, n=1, k=-2
  {1./8 , 3./8 , 1./3 }, // teapot3 -> 0..1, knd=3, j=1, n=2, k=-3
  {3./30, 8./30, 1./4 }, // teapot4 -> 0..2, knd=4, j=2, n=3, k=-4
} ;

void mad_trk_slice_tpt (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick, int knd)
{
  ensure(knd >= 2 && knd <= 4, "invalid teapot kind 2..4");
  int j = knd-2;
  int n = knd-1;
  int k = knd-1;
  int d = m->rflw.sdir; // this is above polymorphic data, hence same for all

    thick(m, lw*teapot[j].d, k++*d);
  FOR(i,n) {
     kick(m, lw*teapot[j].k, k++*d);
    thick(m, lw*teapot[j].D, k++*d);
  }  kick(m, lw*teapot[j].k, k++*d);
    thick(m, lw*teapot[j].d, k++*d);
}

// --- speed tests ------------------------------------------------------------o

#if TPSA_SPDTEST
void mad_trk_spdtest (int n, int k)
{
  mad_desc_newv(6, 1);

  tpsa x ( "X"); x .set( 0   , 1);
  tpsa px("PX"); px.set( 1e-7, 2);
  tpsa y ( "Y"); y .set( 0   , 3);
  tpsa py("PY"); py.set(-1e-7, 4);
  tpsa t ( "T"); t .set( 0   , 5);
  tpsa pt("PT"); pt.set( 0   , 6);

  tpsa_t*   par[] = { x.ptr(), px.ptr(), y.ptr(), py.ptr(), t.ptr(), pt.ptr() };
  tpsa_t* *pars[] = {par};

  union mflw_ m = { .tflw = {
    .name="spdtest", .dbg=0,

    .el=1, .eld=1, .elc=0, .lrad=0,
    .eh=0, .ang=0, .mang=0,

    .pc=1, .beta=1, .betgam=0, .charge=1,

    .sdir=1, .edir=1, .pdir=0, .T=0, .Tbak=-1,

    .k1=0, .ks=0, .volt=0, .freq=0, .lag=0, .sa=0, .ca=1, .nbsl=0,

    .frng=0, .fmax=2, .e=0, .h=0, .a=0, .fint=0, .hgap=0, .f1=0, .f2=0,

    .rot=false, .trn=false,
    .dx=0, .dy=0, .ds=0, .dthe=0, .dphi=0, .dpsi=0, .tlt=0,

    .nmul=1, .knl={1e-7}, .ksl={0},
    .snm=0,  .bfx={0}   , .bfy={0},

    .npar=1, .par=pars  ,// .map=maps,
  }};

  switch(k) {
//  case 0: {
//    FOR(i,n) mad_trk_strex_drift_r (&m, 1, 1);
//    par_t p(&m,0);
//    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
//            p.x, p.px, p.y, p.py, p.t, p.pt);
//  } break;

  case 1: {
    FOR(i,n) mad_trk_strex_drift_t (&m, 1, 1);
    map_t p(m.tflw,0);
    stdout << p.x << p.px << p.y << p.py << p.t << p.pt;
  } break;

//  case 2: {
//    FOR(i,n) {
//      mad_trk_strex_drift_r (&m, 0.5, 1);
//      mad_trk_strex_kick_r  (&m,   1, 1);
//      mad_trk_strex_drift_r (&m, 0.5, 1);
//    }
//    par_t p(&m,0);
//    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
//            p.x, p.px, p.y, p.py, p.t, p.pt);
//  } break;

  case 3: {
    FOR(i,n) {
      mad_trk_strex_drift_t (&m, 0.5, 1);
      mad_trk_strex_kick_t  (&m,   1, 1);
      mad_trk_strex_drift_t (&m, 0.5, 1);
    }
    map_t p(m.tflw,0);
    stdout << p.x << p.px << p.y << p.py << p.t << p.pt;
  } break;

//  case 4: {
//    FOR(i,n) mad_trk_slice_tkt(&m, 1, mad_trk_strex_drift_r, mad_trk_strex_kick_r, 2);
//    par_t p(&m,0);
//    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
//            p.x, p.px, p.y, p.py, p.t, p.pt);
//  } break;

  case 5: {
    FOR(i,n) mad_trk_slice_tkt(&m, 1, mad_trk_strex_drift_t, mad_trk_strex_kick_t, 2);
    map_t p(m.tflw,0);
    stdout << p.x << p.px << p.y << p.py << p.t << p.pt;
  } break;
/*
  case 6: {
    FOR(i,n) mad_trk_slice_tkt(&m, 1, mad_trk_curex_drift_r, mad_trk_curex_kick_r, 2);
    par_t p(&m,0);
    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
            p.x, p.px, p.y, p.py, p.t, p.pt);
  } break;

  case 7: {
    FOR(i,n) mad_trk_slice_tkt(&m, 1, mad_trk_curex_drift_t, mad_trk_curex_kick_t, 2);
    map_t p(&m,0);
    stdout << p.x << p.px << p.y << p.py << p.t << p.pt;
  } break;
*/
  default:
    printf("unknown use case %d\n", k);
  }
}
#endif // TPSA_SPDTEST

/*
time: 0.001262 sec
local t=os.clock() MAD._C.mad_trk_spdtest(1e6,0) print(os.clock()-t, "sec")

time: 0.397203 sec
MAD._C.mad_mcollect()
local t=os.clock() MAD._C.mad_trk_spdtest(1e6,1) print(os.clock()-t, "sec")
MAD._C.mad_mdump(nil)

time: 0.005795 sec
do
local m = {el=1, eld=1, beam={beta=1}, T=0, atdebug=\->(), npar=1,
{x=0,px=1e-7,y=0,py=-1e-7,t=0,pt=0}}
local f =\n=> for i=1,n do MAD.dynmap.strex_drift(nil,m,1,1) end end
m[1].px, m[1].py = 1e-7, -1e-7
local t=os.clock() f(1e6) print(os.clock()-t, "sec")
MAD.utility.printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
m[1].x, m[1].px, m[1].y, m[1].py, m[1].t, m[1].pt)
end

time: 2.679041 sec
do
local m = {el=1, eld=1, beam={beta=1}, T=0, atdebug=\->(), npar=1, MAD.damap()}
m[1].px, m[1].py = 1e-7, -1e-7
local f =\n=> for i=1,n do MAD.dynmap.strex_drift(nil,m,1,1) end end
local t=os.clock() f(1e6) print(os.clock()-t, "sec")
m[1]:print()
end
*/

// --- unit tests -------------------------------------------------------------o

#if TPSA_CPPTEST
#include "mad_ctpsa.hpp"

void mad_trk_cpptest (void)
{
  mad_desc_newv(6, 1);

#if TPSA_USE_TRC
#define TRC(...) printf(#__VA_ARGS__ "\n"); __VA_ARGS__
#else
#define TRC(...) __VA_ARGS__
#endif

// Real version
{ TRC(tpsa a("A");                             )
  TRC(tpsa_ref ar(a.ptr());                    )

  TRC( a   = 1.;                               )
  TRC( a  += 1.;                               )
  TRC( ar  = 1.;                               )
  TRC( ar += 1.;                               )

  TRC( a  = a;                                 )
  TRC( a += a;                                 )
  TRC( a  = ar;                                )
  TRC( a += ar;                                )
  TRC( a  = a+a;                               )
  TRC( a += a+a;                               )
  TRC( a  = 2*a;                               )
  TRC( a += 2*a;                               )
  TRC( a  = 1.;                                )
  TRC( a += 1.;                                )
  TRC( a  = tpsa();                            ) // tpsa(a);       // error
  TRC( a += tpsa();                            ) // tpsa(a);       // error
  TRC( a  = tpsa(a+a);                         )
  TRC( a += tpsa(a+a);                         )
  TRC( a  = tpsa_ref(a.ptr());                 ) // tpsa_ref(a);   // error
  TRC( a += tpsa_ref(a.ptr());                 ) // tpsa_ref(a+b); // error

  TRC( ar  = a;                                )
  TRC( ar += a;                                )
  TRC( ar  = ar;                               )
  TRC( ar += ar;                               )
  TRC( ar  = a+a;                              )
  TRC( ar += a+a;                              )
  TRC( ar  = 2*a;                              )
  TRC( ar += 2*a;                              )
  TRC( ar  = 1.;                               )
  TRC( ar += 1.;                               )
  TRC( ar  = tpsa();                           )  // tpsa(a);       // error
  TRC( ar += tpsa();                           )  // tpsa(a);       // error
  TRC( ar  = tpsa(a+a);                        )
  TRC( ar += tpsa(a+a);                        )
  TRC( ar  = tpsa_ref(a.ptr());                ) // tpsa_ref(a);   // error
  TRC( ar += tpsa_ref(a.ptr());                ) // tpsa_ref(a+b); // error

  TRC( tpsa()  = a;                            )
  TRC( tpsa() += a;                            )
  TRC( tpsa()  = ar;                           )
  TRC( tpsa() += ar;                           )
  TRC( tpsa()  = a+a;                          )
  TRC( tpsa() += a+a;                          )
  TRC( tpsa()  = 2*a;                          )
  TRC( tpsa() += 2*a;                          )
  TRC( tpsa()  = 1.;                           )
  TRC( tpsa() += 1.;                           )
  TRC( tpsa()  = tpsa();                       )  // tpsa(a); // error
  TRC( tpsa() += tpsa();                       )  // tpsa(a); // error
  TRC( tpsa()  = tpsa(a+a);                    )
  TRC( tpsa() += tpsa(a+a);                    )
  TRC( tpsa()  = tpsa_ref(a.ptr());            ) // tpsa_ref(a);   // error
  TRC( tpsa() += tpsa_ref(a.ptr());            ) // tpsa_ref(a+b); // error

  TRC( tpsa_ref(a.ptr())  = a;                 )
  TRC( tpsa_ref(a.ptr()) += a;                 )
  TRC( tpsa_ref(a.ptr())  = ar;                )
  TRC( tpsa_ref(a.ptr()) += ar;                )
  TRC( tpsa_ref(a.ptr())  = a+a;               )
  TRC( tpsa_ref(a.ptr()) += a+a;               )
  TRC( tpsa_ref(a.ptr())  = 2*a;               )
  TRC( tpsa_ref(a.ptr()) += 2*a;               )
  TRC( tpsa_ref(a.ptr())  = 1.;                )
  TRC( tpsa_ref(a.ptr()) += 1.;                )
  TRC( tpsa_ref(a.ptr())  = tpsa();            ) // tpsa(a);       // error
  TRC( tpsa_ref(a.ptr()) += tpsa();            ) // tpsa(a);       // error
  TRC( tpsa_ref(a.ptr())  = tpsa(a+a);         )
  TRC( tpsa_ref(a.ptr()) += tpsa(a+a);         )
  TRC( tpsa_ref(a.ptr())  = tpsa_ref(a.ptr()); ) // tpsa_ref(a);   // error
  TRC( tpsa_ref(a.ptr()) += tpsa_ref(a.ptr()); ) // tpsa_ref(a+b); // error

  TRC( const tpsa b {  a+1+a+2+a+3 };          )
  TRC( const tpsa c { (a+1)*sqr(a+2)+a*2 };    )

  TRC( const tpsa d =  a+1+a+2+a+2;            )
  TRC( const tpsa e = (a+1)*sqr(a+2)+a*2;      )

  TRC(       tpsa f =  a+1+a+2+a+2;            )
  TRC(       tpsa g = (a+1)*sqr(a+2)+a*2;      )
}

// Complex version
{ TRC(ctpsa a("A");                            )
  TRC(ctpsa_ref ar(a.ptr());                   )

  TRC( a   = 1.;                               )
  TRC( a  += 1.;                               )
  TRC( ar  = 1.;                               )
  TRC( ar += 1.;                               )

  TRC( a  = a;                                 )
  TRC( a += a;                                 )
  TRC( a  = ar;                                )
  TRC( a += ar;                                )
  TRC( a  = a+a;                               )
  TRC( a += a+a;                               )
  TRC( a  = 2*a;                               )
  TRC( a += 2*a;                               )
  TRC( a  = 1.;                                )
  TRC( a += 1.;                                )
  TRC( a  = ctpsa();                           ) // ctpsa(a);       // error
  TRC( a += ctpsa();                           ) // ctpsa(a);       // error
  TRC( a  = ctpsa(a+a);                        )
  TRC( a += ctpsa(a+a);                        )
  TRC( a  = ctpsa_ref(a.ptr());                ) // ctpsa_ref(a);   // error
  TRC( a += ctpsa_ref(a.ptr());                ) // ctpsa_ref(a+b); // error

  TRC( ar  = a;                                )
  TRC( ar += a;                                )
  TRC( ar  = ar;                               )
  TRC( ar += ar;                               )
  TRC( ar  = a+a;                              )
  TRC( ar += a+a;                              )
  TRC( ar  = 2*a;                              )
  TRC( ar += 2*a;                              )
  TRC( ar  = 1.;                               )
  TRC( ar += 1.;                               )
  TRC( ar  = ctpsa();                          )  // ctpsa(a);       // error
  TRC( ar += ctpsa();                          )  // ctpsa(a);       // error
  TRC( ar  = ctpsa(a+a);                       )
  TRC( ar += ctpsa(a+a);                       )
  TRC( ar  = ctpsa_ref(a.ptr());               ) // ctpsa_ref(a);   // error
  TRC( ar += ctpsa_ref(a.ptr());               ) // ctpsa_ref(a+b); // error

  TRC( ctpsa()  = a;                           )
  TRC( ctpsa() += a;                           )
  TRC( ctpsa()  = ar;                          )
  TRC( ctpsa() += ar;                          )
  TRC( ctpsa()  = a+a;                         )
  TRC( ctpsa() += a+a;                         )
  TRC( ctpsa()  = 2*a;                         )
  TRC( ctpsa() += 2*a;                         )
  TRC( ctpsa()  = 1.;                          )
  TRC( ctpsa() += 1.;                          )
  TRC( ctpsa()  = ctpsa();                     )  // ctpsa(a); // error
  TRC( ctpsa() += ctpsa();                     )  // ctpsa(a); // error
  TRC( ctpsa()  = ctpsa(a+a);                  )
  TRC( ctpsa() += ctpsa(a+a);                  )
  TRC( ctpsa()  = ctpsa_ref(a.ptr());          ) // ctpsa_ref(a);   // error
  TRC( ctpsa() += ctpsa_ref(a.ptr());          ) // ctpsa_ref(a+b); // error

  TRC( ctpsa_ref(a.ptr())  = a;                )
  TRC( ctpsa_ref(a.ptr()) += a;                )
  TRC( ctpsa_ref(a.ptr())  = ar;               )
  TRC( ctpsa_ref(a.ptr()) += ar;               )
  TRC( ctpsa_ref(a.ptr())  = a+a;              )
  TRC( ctpsa_ref(a.ptr()) += a+a;              )
  TRC( ctpsa_ref(a.ptr())  = 2*a;              )
  TRC( ctpsa_ref(a.ptr()) += 2*a;              )
  TRC( ctpsa_ref(a.ptr())  = 1.;               )
  TRC( ctpsa_ref(a.ptr()) += 1.;               )
  TRC( ctpsa_ref(a.ptr())  = ctpsa();          ) // ctpsa(a);       // error
  TRC( ctpsa_ref(a.ptr()) += ctpsa();          ) // ctpsa(a);       // error
  TRC( ctpsa_ref(a.ptr())  = ctpsa(a+a);       )
  TRC( ctpsa_ref(a.ptr()) += ctpsa(a+a);       )
  TRC( ctpsa_ref(a.ptr())  = ctpsa_ref(a.ptr());) // ctpsa_ref(a);   // error
  TRC( ctpsa_ref(a.ptr()) += ctpsa_ref(a.ptr());) // ctpsa_ref(a+b); // error

  TRC( const ctpsa b {  a+1+a+2+a+3 };         )
  TRC( const ctpsa c { (a+1)*sqr(a+2)+a*2 };   )

  TRC( const ctpsa d =  a+1+a+2+a+2;           )
  TRC( const ctpsa e = (a+1)*sqr(a+2)+a*2;     )

  TRC(       ctpsa f =  a+1+a+2+a+2;           )
  TRC(       ctpsa g = (a+1)*sqr(a+2)+a*2;     )
}

// Real <-> Complex conversion
{ TRC( tpsa a("A");                            )
  TRC( tpsa_ref ar(a.ptr());                   )
  TRC( tpsa b("B");                            )
  TRC( tpsa_ref br(b.ptr());                   )

  TRC(ctpsa c(a);                              )
  TRC(ctpsa_ref cr(c.ptr());                   )
  TRC(ctpsa d(a,b);                            )
  TRC(ctpsa_ref dr(d.ptr());                   )

  c = a;
  c = ar;
  a = real(d);
  b = imag(d);
  a = real(dr);
  b = imag(dr);
}
}
#endif // TPSA_CPPTEST
