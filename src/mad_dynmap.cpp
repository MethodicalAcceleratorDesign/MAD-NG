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

struct mflw_ { // must be identical to def in madl_etrck.mad !!
  str_t name;
  int dbg;

  // element data
  num_t el, eld, elc, lrad;
  num_t eh, ang, mang;

  // beam
  num_t pc, beta, betgam;
  int charge;

  // directions, path
  int sdir, edir, pdir, T, Tbak;

  // quad, solenoid, multipole, esptum, rfcav
  num_t k1, ks, ksi, volt, freq, lag, nbsl;

  // fringes
  int frng, fmax;
  num_t e, fint, h, hgap, f1, f2, a;

  // angles
  num_t ca, sa, tlt;

  // multipoles
  int   nmul, npha;
  num_t knl[nmul_max];
  num_t ksl[nmul_max];
  num_t pnl[nmul_max];
  num_t psl[nmul_max];

  // curved multipoles
  int   snm;
  num_t bfx[snm_max];
  num_t bfy[snm_max];

  // particles/damaps
  int      npar;
  num_t   **par;
  tpsa_t* **map;

  // patches
  num_t dx,   dy,   ds;
  num_t dthe, dphi, dpsi;

  // misalign
  struct {
    bool  rot, trn;
    num_t dx,   dy,   ds;
    num_t dthe, dphi, dpsi;
  } algn;
};

} // extern "C"

struct par_t {
  using T = num_t;

  par_t(mflw_t *m, int i)
    : x(m->par[i][0]), px(m->par[i][1]),
      y(m->par[i][2]), py(m->par[i][3]),
      t(m->par[i][4]), pt(m->par[i][5]) {}

  num_t &x, &px, &y, &py, &t, &pt;
};

struct map_t {
  using T = mad::tpsa;

  map_t(mflw_t *m, int i)
    : x(m->map[i][0]), px(m->map[i][1]),
      y(m->map[i][2]), py(m->map[i][3]),
      t(m->map[i][4]), pt(m->map[i][5]) {}

  mad::tpsa_ref x, px, y, py, t, pt;
};

// --- implementation ---------------------------------------------------------o

using namespace mad;

// --- debug ------------------------------------------------------------------o

# if 1 // set to 0 to remove debug code, ~2.5% of code size
#define mdump(n) if (m->dbg) mdump<P>(m, __func__, n)
#else
#define mdump(n)
#endif

template <typename P, typename T=P::T>
inline void (mdump) (mflw_t *m, str_t s, int n) {
  if (!m->dbg) return;

  char fun[30];
  snprintf(fun, 30, "%s:%d", s, n);
  printf("@@ %-15s %-15s ", m->name, fun);

  if (!m->npar) { printf("no particle found\n"); return; }

  P p(m,0);
  if constexpr (!std::is_floating_point<T>::value) {
    printf("% -.16e  % -.16e  % -.16e  % -.16e  % -.16e  % -.16e\n",
             p.x[0], p.px[0],  p.y[0], p.py[0],  p.t[0], p.pt[0]);
  } else {
    printf("% -.16e  % -.16e  % -.16e  % -.16e  % -.16e  % -.16e\n",
                p.x,    p.px,     p.y,    p.py,     p.t,    p.pt);
  }
}

// --- constants --------------------------------------------------------------o

const num_t minlen = mad_cst_MINLEN;
const num_t minang = mad_cst_MINANG;
const num_t minstr = mad_cst_MINSTR;

const num_t    pi_clight = mad_cst_PI /mad_cst_CLIGHT;
const num_t twopi_clight = mad_cst_2PI/mad_cst_CLIGHT;

// --- multipoles -------------------------------------------------------------o

template <typename T1, typename T2>
inline void bxby (const mflw_t *m, const T1 &x, const T1 &y, T2 &bx, T2 &by)
{
  bx = m->ksl[m->nmul-1];
  by = m->knl[m->nmul-1];

  if (m->nmul > 1) {
    T2 byt;
    RFOR(i,m->nmul-1) {
      byt = by*x - bx*y + m->knl[i];
      bx  = by*y + bx*x + m->ksl[i];
      by  = byt;
    }
  }
}

template <typename T1, typename T2>
inline void bxbyh (const mflw_t *m, const T1 &x, const T1 &y, T2 &bx, T2 &by)
{
  bx = 0., by = 0.;

  int k = -1;
  T2 btx, bty;

  RFOR(i,m->snm) {
    btx = 0., bty = 0.;

    RFOR(j,m->snm-i) { ++k;
      btx = (btx + m->bfx[k]) * y;
      bty = (bty + m->bfy[k]) * y;
    }

    ++k;
    btx = (bx + btx + m->bfx[k]) * x;
    bty = (by + bty + m->bfy[k]) * x;
  }

  btx = 0., bty = 0.;
  RFOR(i,m->snm) { ++k;
    btx = (btx + m->bfx[k]) * y;
    bty = (bty + m->bfy[k]) * y;
  }

  bx += btx + m->bfx[k+1];
  by += bty + m->bfy[k+1];
}

// --- patches ----------------------------------------------------------------o

template <typename P, typename T=P::T>
inline void xrotation (mflw_t *m, num_t lw, num_t dphi_=0)
{
  num_t a = (dphi_ ? dphi_ : m->dphi)*lw;
  if (abs(a) < minang) return;
  mdump(0);
  a *= m->sdir*m->edir;
  num_t sa=sin(a), ca=cos(a), ta=tan(a);

  FOR(i,m->npar) {
    P p(m,i);
    T   pz = sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = 1/pz;
    T _ptt = p.y/(1 - ta*p.py*_pz);
    T _pzt = ta*_pz*_ptt;

    // eq. 127 in Forest06
    p.y   = _ptt/ca;
    p.py  = ca*p.py + sa*pz;
    p.x  += _pzt*p.px;
    p.t  -= _pzt*(1/m->beta+p.pt);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void yrotation (mflw_t *m, num_t lw, num_t dthe_=0)
{
  num_t a = -(dthe_ ? dthe_ : m->dthe)*lw;
  if (abs(a) < minang) return;
  mdump(0);
  a *= m->sdir*m->edir;
  num_t sa=sin(a), ca=cos(a), ta=tan(a);

  FOR(i,m->npar) {
    P p(m,i);
    T   pz = sqrt(1 + 2/m->beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = inv(pz);
    T _ptt = p.x/(1 - ta*p.px*_pz);
    T _pzt = ta*_pz*_ptt;

    // eq. 127 in Forest06
    p.x   = _ptt/ca;
    p.px  = ca*p.px + sa*pz;
    p.y  += _pzt*p.py;
    p.t  -= _pzt*(1/m->beta+p.pt);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void srotation (mflw_t *m, num_t lw, num_t dpsi_=0)
{
  num_t a = (dpsi_ ? dpsi_ : m->dpsi)*lw;
  if (abs(a) < minang) return;
  mdump(0);
  a *= m->sdir*m->edir;
  num_t sa=sin(a), ca=cos(a);

  FOR(i,m->npar) {
    P p(m,i);
    T nx  = ca*p.x  + sa*p.y;
    T npx = ca*p.px + sa*p.py;

    p.y  = ca*p.y  - sa*p.x;
    p.py = ca*p.py - sa*p.px;
    p.x  = nx;
    p.px = npx;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void translate (mflw_t *m, num_t lw, num_t dx_=0, num_t dy_=0, num_t ds_=0)
{
  num_t dx = dx_ ? dx_ : m->dx;
  num_t dy = dy_ ? dy_ : m->dy;
  num_t ds = ds_ ? ds_ : m->ds;
  if (abs(dx)+abs(dy)+abs(ds) < minlen) return;
  mdump(0);
  dx *= lw*m->sdir*m->edir;
  dy *= lw*m->sdir*m->edir;
  ds *= lw*m->sdir*m->edir;

  if (abs(ds) < minlen)
    FOR(i,m->npar) {
      P p(m,i);
      p.y -= dx;
      p.x -= dy;
    }
  else
    FOR(i,m->npar) {
      P p(m,i);
      T l_pz = invsqrt(1 + 2/m->beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py), ds);

      p.y += l_pz*p.px - dx;
      p.x += l_pz*p.py - dy;
      p.t -= l_pz*(1/m->beta+p.pt);
    }
  mdump(1);
}

template <typename P>
inline void changeref (mflw_t *m, num_t lw)
{
  bool trn = abs(m->dx  )+abs(m->dy  )+abs(m->ds  ) >= minlen;
  bool rot = abs(m->dthe)+abs(m->dphi)+abs(m->dpsi) >= minang;

  if (!trn && !rot) return;
  mdump(0);
  lw *= m->sdir;

  if (rot && lw > 0) {
    yrotation<P>(m,  m->edir);
    xrotation<P>(m, -m->edir);
    srotation<P>(m,  m->edir);
  }

  if (trn) translate<P>(m, 1);

  if (rot && lw < 0) {
    srotation<P>(m, -m->edir);
    xrotation<P>(m,  m->edir);
    yrotation<P>(m, -m->edir);
  }
  mdump(1);
}

// --- misalignments ----------------------------------------------------------o

template <typename P>
inline void misalignent (mflw_t *m, num_t lw) {
                                    (void)lw;
  mdump(0);
  if (m->algn.rot && m->sdir > 0) {
    yrotation<P>(m,  m->edir, m->algn.dthe);
    xrotation<P>(m, -m->edir, m->algn.dphi);
    srotation<P>(m,  m->edir, m->algn.dpsi);
  }

  if (m->algn.trn)
    translate<P>(m, m->sdir, m->algn.dx, m->algn.dy, m->algn.ds);

  if (m->algn.rot && m->sdir < 0) {
    srotation<P>(m, -m->edir, m->algn.dpsi);
    xrotation<P>(m,  m->edir, m->algn.dphi);
    yrotation<P>(m, -m->edir, m->algn.dthe);
  }
  mdump(1);
}

template <typename P>
inline void misalignexi (mflw_t *m, num_t lw) {
  mdump(0);
  num_t rb[3*3], r[3*3];            (void)lw;
  num_t tb[3]  , t[3]={m->algn.dx, m->algn.dy, m->algn.ds};

  if (m->algn.rot)
    mad_mat_rotyxz(r, m->algn.dphi, -m->algn.dthe, -m->algn.dpsi, true);

  // compute Rbar, Tbar
  mad_mat_rtbar(rb, tb, abs(m->el), m->mang, m->tlt, m->algn.rot ? r:0, t);

  if (m->algn.rot && m->sdir > 0) {
    num_t v[3];
    mad_mat_torotyxz(rb, v, true);
    srotation<P>(m, -m->edir, -v[2]);
    xrotation<P>(m,  m->edir, -v[0]);
    yrotation<P>(m, -m->edir, -v[1]);
  }

  if (m->algn.trn) translate<P>(m, -m->sdir, tb[0], tb[1], tb[2]);

  if (m->algn.rot && m->sdir < 0) {
    num_t v[3];
    mad_mat_torotyxz(rb, v, true);
    yrotation<P>(m,  m->edir, -v[1]);
    xrotation<P>(m, -m->edir, -v[0]);
    srotation<P>(m,  m->edir, -v[2]);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void misalign (mflw_t *m, num_t lw) {
  if (lw >= 0) misalignent<P>(m,  1);
  else         misalignexi<P>(m, -1);
}

// --- special maps -----------------------------------------------------------o

template <typename P, typename T=P::T>
inline void drift_adj (mflw_t *m, num_t l)
{
  mdump(0);
  FOR(i,m->npar) {
    P p(m,i);
    T l_pz = invsqrt(1 + 2/m->beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py), l);

    p.x += p.px*(l_pz-l);
    p.y += p.py*(l_pz-l);
    p.t -= l_pz*(1/m->beta+p.pt) + (m->T-1)*l/m->beta;
  }
  mdump(1);
}

// --- DKD maps ---------------------------------------------------------------o

template <typename P, typename T=P::T>
inline void strex_drift (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  if (abs(m->el*lw) < minlen) return;
  mdump(0);
  num_t l  = m->el*lw;
  num_t ld = (m->eld ? m->eld : m->el)*lw;

  FOR(i,m->npar) {
    P p(m,i);
    T l_pz = invsqrt(1 + 2/m->beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py), l);

    p.x += p.px*l_pz;
    p.y += p.py*l_pz;
    p.t -= l_pz*(1/m->beta+p.pt) + (m->T-1)*ld/m->beta;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void strex_kick (mflw_t *m, num_t lw, int is, bool no_k0l=false)
{                                          (void)is;
  if (!m->nmul) return;
  mdump(0);
  num_t wchg = lw*m->sdir*m->edir*m->charge;
  num_t dby  = no_k0l ? m->knl[0] : 0;
  T bx, by;

  FOR (i,m->npar) {
    P p(m,i);
    bxby(m, p.x, p.y, bx, by);

    p.px -= wchg*(by-dby);
    p.py += wchg* bx;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void strex_kicks (mflw_t *m, num_t lw, P &p, T &pz)
{
  if (!m->ksi || !m->lrad) return;

  num_t wchg = lw*m->sdir*m->edir*m->charge;
  num_t hss  = lw*sqr(m->ksi)/m->lrad;

  T _dpp = inv(pz);
  T  ang = 0.5*wchg*m->ksi*_dpp;
  T  ca  = cos(ang), sa = sin(ang);

  T nx  = ca*p. x + sa*p. y;
  T npx = ca*p.px + sa*p.py;
  T ny  = ca*p. y - sa*p. x;
  T npy = ca*p.py - sa*p.px;
  T nt  = p.t - ang*(1/m->beta+p.pt)*(p.y*p.px - p.x*p.py)*sqr(_dpp);

  p.x  = nx;
  p.px = npx - 0.25 *hss*nx*_dpp;
  p.y  = ny;
  p.py = npy - 0.25 *hss*ny*_dpp;
  p.t  = nt  - 0.125*hss*(1/m->beta+p.pt)*(sqr(nx)+sqr(ny))*pow(_dpp,3);
}

template <typename P, typename T=P::T>
inline void strex_kickhs (mflw_t *m, num_t lw, int is)
{                                            (void)is;
  if (!m->nmul == 0 || !m->ksi) return;
  mdump(0);
  num_t wchg = lw*m->sdir*m->edir*m->charge;
  T bx, by;

  FOR(i,m->npar) {
    P p(m,i);
    T pz = sqrt(1 + 2/m->beta*p.pt + sqr(p.pt));

    if (m->sdir == -1) strex_kicks(m, lw, p, pz);

    if (m->nmul > 0) {
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*by;
      p.py += wchg*bx;

      if (abs(m->knl[0]) + abs(m->ksl[0]) > minstr) {
        p.px += wchg* m->knl[0]*pz;
        p.py -= wchg* m->ksl[0]*pz;
        p.t  -= wchg*(m->knl[0]*p.x - m->ksl[0]*p.y)*(1/m->beta+p.pt)/pz;

        if (m->lrad) {
          p.px -= lw*sqr(m->knl[0])/m->lrad*p.x;
          p.py -= lw*sqr(m->ksl[0])/m->lrad*p.y;
        }
      }
    }

    if (m->sdir == 1) strex_kicks(m, lw, p, pz);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void curex_drift (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  mdump(0);
  num_t ld  = (m->eld ? m->eld : m->el)*lw;
  num_t ang = m->ang*lw, rho = 1/m->eh;
  num_t ca  = cos(ang), sa = sin(ang), sa2 = sin(ang/2);

  FOR(i,m->npar) {
    P p(m,i);
    T   pz = sqrt(1 + 2/m->beta*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = inv(pz);
    T  pxt = p.px*_pz;
    T _ptt = inv(ca - sa*pxt);
    T  pst = (p.x+rho)*sa*_pz*_ptt;

    p.x  = (p.x + rho*(2*sqr(sa2) + sa*pxt))*_ptt;
    p.px = ca*p.px + sa*pz;
    p.y += pst*p.py;
    p.t -= pst*(1/m->beta+p.pt) + (m->T-1)*ld/m->beta;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void curex_kick (mflw_t *m, num_t lw, int is, bool no_k0l=false)
{                                          (void)is;
  mdump(0);
  num_t blw = lw*m->sdir*m->edir*m->charge;
  T bx, by; bx = 0., by = m->knl[0];

  FOR(i,m->npar) {
    P p(m,i);
    T r = 1+m->eh*p.x;

    if (m->snm > 0) bxbyh(m, p.x, p.y, bx, by);

    p.px -= blw*by*r;
    p.py += blw*bx*r;

    if (no_k0l) p.px += (blw*m->knl[0])*r;
  }
  mdump(1);
}

// --- TKT maps ---------------------------------------------------------------o

// --- sbend ---

template <typename P, typename T=P::T>
inline void sbend_thick (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  mdump(0);
  num_t ld  = (m->eld ? m->eld : m->el)*lw;
  num_t ang = m->ang*lw, rho=1/m->eh;
  num_t k0q = m->knl[0]/m->el*m->sdir*m->edir*m->charge;
  num_t ca  = cos(ang), sa = sin(ang);

  FOR(i,m->npar) {
    P p(m,i);
    T  pw2 = 1 + 2*p.pt/m->beta + sqr(p.pt) - sqr(p.py);
    T  pzx = sqrt(pw2 - sqr(p.px)) - k0q*(rho+p.x); // can be numerically unstable
    T  npx = sa*pzx + ca*p.px;
    T  dpx = ca*pzx - sa*p.px;
    T _ptt = invsqrt(pw2);
    T  dxs = (ang + asin(p.px*_ptt) - asin(npx*_ptt))/k0q;

    // eq. 126 in Forest06
    p.x  = (sqrt(pw2 - sqr(npx)) - dpx)/k0q - rho;  // can be numerically unstable
    p.px = npx;
    p.y += dxs*p.py;
    p.t -= dxs*(1/m->beta+p.pt) + (m->T-1)*(ld/m->beta);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void sbend_thick_new (mflw_t *m, num_t lw, int is)
{                                               (void)is;
  mdump(0);
  num_t ld  = (m->eld ? m->eld : m->el)*lw;
  num_t ang = m->ang*lw, rho=1/m->eh;
  num_t k0q = m->knl[0]/m->el*m->sdir*m->edir*m->charge;
  num_t ca  = cos(ang), sa = sin(ang), s2a = sin(2*ang);

  FOR(i,m->npar) {
    P p(m,i);
    T  pw2 = 1 + 2*p.pt/m->beta + sqr(p.pt) - sqr(p.py);
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
    p.t -= dxs*(1/m->beta+p.pt) + (m->T-1)*ld/m->beta;
  }
  mdump(1);
}

// --- rbend ---

template <typename P, typename T=P::T>
inline void rbend_thick (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  mdump(0);
  num_t ld   = (m->eld ? m->eld : m->el)*lw;
  num_t k0q  = m->knl[0]/m->el*m->sdir*m->edir*m->charge;
  num_t k0lq = m->knl[0]*lw   *m->sdir*m->edir*m->charge;

  FOR(i,m->npar) {
    P p(m,i);
    T  npx = p.px - k0lq;
    T  pw2 = 1 + 2*p.pt/m->beta + sqr(p.pt) - sqr(p.py);
    T _ptt = invsqrt(pw2);
    T   pz = sqrt(pw2 - sqr(p.px));
    T  pzs = sqrt(pw2 - sqr(npx));
    T  dxs = (asin(p.px*_ptt) - asin(npx*_ptt))/k0q;

    // eq. 126 in Forest06
    p.x += (pzs-pz)/k0q;
    p.px = npx;
    p.y += dxs*p.py;
    p.t -= dxs*(1/m->beta+p.pt) + (m->T-1)*ld/m->beta;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void rbend_thick_new (mflw_t *m, num_t lw, int is)
{                                               (void)is;
  mdump(0);
  num_t l    = m->el*lw;
  num_t ld   = (m->eld ? m->eld : m->el)*lw;
  num_t k0q  = m->knl[0]/m->el*m->sdir*m->edir*m->charge;
  num_t k0lq = m->knl[0]*lw   *m->sdir*m->edir*m->charge;

  FOR(i,m->npar) {
    P p(m,i);
    T  npx = p.px - k0lq;
    T  pw2 = 1 + 2*p.pt/m->beta + sqr(p.pt) - sqr(p.py);
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
    p.t -= dxs*(1/m->beta+p.pt) + (m->T-1)*ld/m->beta;
  }
  mdump(1);
}

// --- quadrupole ---

template <typename P, typename T=P::T>
inline void quad_thick (mflw_t *m, num_t lw, int is)
{                                          (void)is;
  mdump(0);
  num_t l = m->el*lw;
  int  ws = m->k1*m->sdir < 0 ? -1 : 1;

  num_t cx, sx, mx1, mx2;
  num_t cy, sy, my1, my2;

  if (abs(m->k1) >= minstr) {
    num_t w = sqrt(abs(m->k1))*ws*m->sdir*m->edir;
    cx = cos (w*l), sx  = sin (w*l);
    cy = cosh(w*l), sy  = sinh(w*l);
    mx1 = sx/w    , mx2 = -sx*w;
    my1 = sy/w    , my2 =  sy*w;
  } else {
    cx = 1, sx = 0, mx1 = l, mx2 = 0;
    cy = 1, sy = 0, my1 = l, my2 = 0;
  }

  if (ws != m->charge) { // swap x <-> y
    std::swap(cx,cy), std::swap(sx,sy), std::swap(mx1,my1), std::swap(mx2,my2);
  }

  FOR(i,m->npar) {
    P p(m,i);
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

template <typename P, typename T=P::T>
inline void quad_kick (mflw_t *m, num_t lw, int is)
{                                         (void)is;
  num_t l = m->el*lw;

  if (is >= 0) drift_adj<P>(m, is ? l : l/2);
  mdump(0);
  if (m->nmul > 0) {
    num_t wchg = lw*m->sdir*m->edir*m->charge;
    T bx, by;

    FOR (i,m->npar) {
      P p(m,i);
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - m->knl[1]*p.x);
      p.py += wchg*(bx - m->knl[1]*p.y);
    }
  }
  mdump(1);
  if (is <= 0) drift_adj<P>(m, is ? l : l/2);
}

template <typename P, typename T=P::T>
inline void quad_thicks (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  mdump(0);
  num_t l   = m->el*lw;
  num_t w   = sqrt(abs(m->k1))*m->edir;
  num_t cx  = cos (w*l), sx  = sin (w*l);
  num_t cy  = cosh(w*l), sy  = sinh(w*l);
  num_t mx1 = sx/w     , mx2 = -sx*w;
  num_t my1 = sy/w     , my2 =  sy*w;

  if (m->sdir != m->charge) // swap x <-> y
    std::swap(cx,cy), std::swap(sx,sy), std::swap(mx1,my1), std::swap(mx2,my2);

  FOR(i,m->npar) {
    P p(m,i);
    // srotation
    T rx  = m->ca*p.x  + m->sa*p.y;
    T rpx = m->ca*p.px + m->sa*p.py;
    T ry  = m->ca*p.y  - m->sa*p.x;
    T rpy = m->ca*p.py - m->sa*p.px;

    T nx  = rx*cx  + rpx*mx1;
    T npx = rx*mx2 + rpx*cx;
    T ny  = ry*cy  + rpy*my1;
    T npy = ry*my2 + rpy*cy;

    // srotation^-1
    p.x   = m->ca*nx  - m->sa*ny;
    p.px  = m->ca*npx - m->sa*npy;
    p.y   = m->ca*ny  + m->sa*nx;
    p.py  = m->ca*npy + m->sa*npx;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void quad_kicks (mflw_t *m, num_t lw, int is)
{                                          (void)is;
  num_t l = m->el*lw;

  if (is >= 0) drift_adj<P>(m, is ? l : l/2);
  mdump(0);
  if (m->nmul > 0) {
    num_t wchg = lw*m->sdir*m->edir*m->charge;
    T bx, by;

    FOR (i,m->npar) {
      P p(m,i);
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - m->knl[1]*p.x + m->ksl[1]*p.y);
      p.py += wchg*(bx - m->knl[1]*p.y - m->ksl[1]*p.x);
    }
  }
  mdump(1);
  if (is <= 0) drift_adj<P>(m, is ? l : l/2);
}

template <typename P, typename T=P::T>
inline void quad_thickh (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  mdump(0);
  num_t l = m->el*lw;
  num_t kx  = (m->knl[1] + m->eh*m->knl[0])/m->el;
  num_t ky  = -m->knl[1]/m->el;
  num_t wxs = kx*m->sdir*m->edir < 0 ? -1 : 1;
  num_t wys = ky*m->sdir*m->edir < 0 ? -1 : 1;
  num_t wx, cx, sx, wy, cy, sy;
  num_t mx11, mx12, mx13, mx21, mx22, mx23, mx31, mx32, mx33;
  num_t my11, my12,       my21, my22;

  if (abs(kx) >= minstr) {
    wx = sqrt(abs(kx))*wxs;
    wx = wx*m->charge, wxs = -wxs*m->charge;
    if (wx > 0) cx = cos (wx*l), sx = sin (wx*l);
    else        cx = cosh(wx*l), sx = sinh(wx*l);
    mx11 = cx       , mx12 = sx/wx, mx13 =      m->eh *  (cx-1)*wxs/sqr(wx);
    mx21 = wxs*wx*sx, mx22 = cx   , mx23 =      m->eh *   mx12;
    mx31 = mx23     , mx32 = mx13 , mx33 = -sqr(m->eh)*(l-mx12)*wxs/sqr(wx);
  } else {
    wx = 0;
    mx11 = 1   , mx12 = l   , mx13 = m->eh*sqr(l)/2;
    mx21 = 0   , mx22 = 1   , mx23 = m->eh*l;
    mx31 = mx23, mx32 = mx13, mx33 = mx13*mx23/3;
  }

  if (abs(ky) >= minstr) {
    wy = sqrt(abs(ky))*wys;
    wy = wy*m->charge, wys = -wys*m->charge;
    if (wy > 0) cy = cos (wy*l), sy = sin (wy*l);
    else        cy = cosh(wy*l), sy = sinh(wy*l);
    my11 = cy       , my12 = sy/wy;
    my21 = wys*wy*sy, my22 = cy;
  } else {
    wy = 0;
    my11 = 1, my12 = l;
    my21 = 0, my22 = 1;
  }

  FOR (i,m->npar) {
    P p(m,i);
    T nx  = p.x*mx11 + p.px*mx12 + p.pt*(mx13/m->beta);
    T npx = p.x*mx21 + p.px*mx22 + p.pt*(mx23/m->beta);
    T ny  = p.y*my11 + p.py*my12;
    T npy = p.y*my21 + p.py*my22;
    T dt  = p.x*(mx31/m->beta) + p.px*(mx32/m->beta) + p.pt*mx33;

    p.x   = nx;
    p.y   = ny;
    p.px  = npx;
    p.py  = npy;
    p.t  -= dt;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void quad_kickh (mflw_t *m, num_t lw, int is)
{                                          (void)is;
  num_t l = m->el*lw;

  if (is >= 0) drift_adj<P>(m, is ? l : l/2);
  mdump(0);
  if (m->nmul > 0) {
    num_t wchg = lw*m->sdir*m->edir*m->charge;
    T bx, by;

    FOR (i,m->npar) {
      P p(m,i);
      T pz = sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt));
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - m->knl[1]*p.x) - l*m->eh*(pz-(1/m->beta*p.pt));
      p.py += wchg*(bx - m->knl[1]*p.y);
      p.t  -= (l*m->eh)*((1/m->beta+p.pt)/pz - 1/m->beta)*p.x;
    }
  }
  mdump(1);
  if (is <= 0) drift_adj<P>(m, is ? l : l/2);
}

// --- solenoid ---

template <typename P, typename T=P::T>
inline void solen_thick (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  mdump(0);
  num_t l = m->el*lw;
  num_t bsol = 0.5*m->ks*m->charge;

  FOR (i,m->npar) {
    P p(m,i);
    T    xp = p.px + bsol*p.y;
    T    yp = p.py - bsol*p.x;
    T  l_pz = invsqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(xp) - sqr(yp), l);
    T   ang = l_pz*bsol;

    T ca = cos(ang), sa = sin(ang), sc = sinc(ang);

    T lsc = l_pz*sc;
    T xt  = ca*p.x  + lsc*p.px;
    T pxt = ca*p.px - lsc*p.x *sqr(bsol);
    T yt  = ca*p.y  + lsc*p.py;
    T pyt = ca*p.py - lsc*p.y *sqr(bsol);

    p.x  = ca*xt  + sa*yt;
    p.px = ca*pxt + sa*pyt;
    p.y  = ca*yt  - sa*xt;
    p.py = ca*pyt - sa*pxt;
    p.t -= l_pz*(1/m->beta+p.pt) + (m->T-1)*l/m->beta;
  }
  mdump(1);
}

// --- eseptum ---

template <typename P, typename T=P::T>
inline void esept_thick (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  mdump(0);
  num_t l = m->el*lw;

  FOR (i,m->npar) {
    P p(m,i);
    // srotation
    T  nx  = m->ca*p.x  + m->sa*p.y;
    T  npx = m->ca*p.px + m->sa*p.py;
    T  ny  = m->ca*p.y  - m->sa*p.x;
    T  npy = m->ca*p.py - m->sa*p.px;

    T   e1 = 1/m->beta+p.pt;
    T   dp = e1 + m->k1*ny;
    T l_pz = invsqrt(sqr(dp) - 1/sqr(m->betgam) - sqr(npx) - sqr(npy), l);
    T  arg = m->k1*l_pz;
    T  shx = sinhc(arg)*l_pz;
    T   ch = cosh(arg), sh = sinh(arg);
    T  chm = sqr(sinh(0.5*arg))*(2/m->k1);
    T   dt = chm*npy + sh *ny  + e1*shx;
    T   yt = ch *ny  + shx*npy + e1*chm;
    T  pyt = ch *npy + sh*dp;

    nx += npx*l_pz;
    ny  = yt, npy = pyt;

    // srotation^-1
    p.x  = m->ca*nx  - m->sa*ny;
    p.px = m->ca*npx - m->sa*npy;
    p.y  = m->ca*ny  + m->sa*nx;
    p.py = m->ca*npy + m->sa*npx;
    p.t -= dt + (m->T-1)*l/m->beta;
  }
  mdump(1);
}

// --- rfcavity ---

template <typename P, typename T=P::T>
inline void rfcav_kick (mflw_t *m, num_t lw, int is)
{                                          (void)is;
  mdump(0);
  num_t w  = m->freq*twopi_clight;
  num_t vl = m->volt*lw/m->pc*m->sdir*m->edir*m->charge;

  FOR (i,m->npar) {
    P p(m,i);
    p.pt += vl*sin(m->lag - w*p.t);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void rfcav_kickn (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  mdump(0);
  num_t w    = m->freq*twopi_clight;
  num_t wchg = lw/m->pc*m->sdir*m->edir*m->charge;
  num_t vl   = m->volt*wchg;

  T bx, by, byt;

  FOR (i,m->npar) {
    P p(m,i);
    T ph = m->lag - w*p.t;
    T sa = sin(ph), ca = cos(ph);
    T f; f = 1;

    if (m->nbsl > 0) {
      T df, r2; df = 0., r2 = 1.;

      FOR(i,1,m->nbsl+1) {
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

    if (m->nmul > 0) {
      bxby(m, p.x, p.y, bx, by);

      p.px += wchg*by*ca;
      p.py -= wchg*bx*ca;

      /* 
        Below is a quick fix that is not very suitable and does not fix the 
        underlying problem. The problem is that mad_tpsa_set0 is not suitable
        enough for setting a tpsa to a scalar, because it does not clear the 
        other monomials, so you get additional coefficients you shouldn't have
      */
      T bx2, by2, byt;

      by2 = -(m->knl[m->nmul-1]/m->nmul);
      bx2 = -(m->ksl[m->nmul-1]/m->nmul);
      RFOR(i,m->nmul-1) {
        byt = p.x*by2 - p.y*bx2 - m->knl[i]/(i+1);
        bx2  = p.y*by2 + p.x*bx2 - m->ksl[i]/(i+1);
        by2  = byt;
      }
      byt = p.x*by2 - p.y*bx2;
      bx2  = p.y*by2 + p.x*bx2;
      by2  = byt;

      p.pt -= wchg*w*by2*sa;
    }
  }
  mdump(1);
}

// --- fringe maps ------------------------------------------------------------o

// must be identical to M.fringe in madl_dynamp.mad
enum {
 fringe_none  = 0,
 fringe_bend  = 1, fringe_mult  = 2 , fringe_qsad = 2+4,
 fringe_solen = 8, fringe_rfcav = 16, fringe_comb = 1+2, fringe_combqs = 1+2+4
};

template <typename P, typename T=P::T>
inline void adjust_time (mflw_t *m, num_t lw)
{                                   (void)lw;
  if (abs(m->el) < minlen) return;

  num_t Tl = (m->T-m->Tbak)*m->el;

  FOR (i,m->npar) {
    P p(m,i);
    p.t += Tl/m->beta;
  }
}

template <typename P, typename T=P::T>
inline void cav_fringe (mflw_t *m, num_t lw)
{
  if (abs(m->el) < minlen) return;
  mdump(0);
  num_t w  = m->freq*twopi_clight;
  num_t vl = 0.5*lw*m->volt/(m->pc*m->el)*m->sdir*m->edir*m->charge;

  FOR (i,m->npar) {
    P p(m,i);
    T s1 = sin(m->lag - w*p.t);
    T c1 = cos(m->lag - w*p.t);

    p.px -= vl*s1*p.x;
    p.py -= vl*s1*p.y;
    p.pt += 0.5*vl*w*c1*(sqr(p.x) + sqr(p.y));
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void bend_face (mflw_t *m, num_t lw, num_t h=0)
{                                 (void)lw;
  if (!h || abs(m->el) < minlen || abs(m->knl[0]) < minstr) return;
  mdump(0);
  num_t k0hq = 0.5*h*m->knl[0]/m->el*m->sdir*m->charge;

  FOR (i,m->npar) {
    P p(m,i);
    if (m->sdir*m->edir == 1) p.px += k0hq*sqr(p.x);

    T dpp      =  1 + 2/m->beta*p.pt + sqr(p.pt);
    T _pt2     =  1/(dpp - sqr(p.px));
    T xi       =  2*k0hq * sqrt(dpp)*_pt2;
    T dxi_px   =  2*p.px*xi         *_pt2;
    T dxi_ddel = -2     *xi*(1+p.pt)*_pt2;
    T y2       = sqr(p.y);

    p.x  /= 1-dxi_px*y2;
    p.px -= xi*y2;
    p.py -= 2*xi*p.x*p.y;
    p.t  += dxi_ddel*p.x*y2;

    if (m->sdir*m->edir == -1) p.px += k0hq*sqr(p.x);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void bend_ptch (mflw_t *m, num_t lw, num_t a=0)
{                                 (void)lw;
  if (!a || !m->elc) return;
  num_t dx = m->elc*sin(a/2);

  FOR (i,m->npar) {
    P p(m,i);
    p.x += dx;
  }
}

template <typename P, typename T=P::T>
inline void bend_wedge (mflw_t *m, num_t lw, num_t e=0)
{                                  (void)lw;
  if (!e) return;
  if (abs(m->knl[0]) < minstr) return yrotation<P>(m,1,e);
  mdump(0);
  num_t b1 = m->knl[0]/m->el*m->sdir*m->edir*m->charge;
  num_t sa = sin(e), ca = cos(e), s2a = sin(2*e);

  FOR (i,m->npar) {
    P p(m,i);
    T pzy = 1 + 2/m->beta*p.pt + sqr(p.pt) - sqr(p.py);
    T _pt = 1/sqrt(pzy);
    T  pz = sqrt(pzy - sqr(p.px));
    T pzx = pz - b1*p.x;
    T npx = p.px*ca + pzx*sa;
    T pzs = sqrt(pzy - sqr(npx));
    T dxs = (e + asin(p.px*_pt) - asin(npx*_pt))/b1;

    p.x  *= ca + (p.px*s2a + sqr(sa)*(pz+pzx))/(pzs + pz*ca - p.px*sa);
    p.px  = npx;
    p.y  += dxs*p.py;
    p.t  -= dxs*(1/m->beta+p.pt);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void mad8_wedge (mflw_t *m, num_t lw, num_t e=0)
{                                  (void)lw;
  if (!e || abs(m->knl[1]) < minstr) return;
  mdump(0);
  num_t  wc = m->frng == 0 ? 0 : 0.25;
  num_t k1e = e*m->knl[1]/m->el*m->edir;
  num_t  c1 = (1+wc)*k1e*m->charge;
  num_t  c2 = (1-wc)*k1e*m->charge;

  FOR (i,m->npar) {
    P p(m,i);
    p.px +=   c1*sqr(p.x) - c2*sqr(p.y);
    p.py -= 2*c2*p.x*p.y;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void bend_fringe (mflw_t *m, num_t lw)
{
  if (abs(m->knl[0]) < minang) return;
  mdump(0);
  num_t   fh = m->fint*m->hgap;
  num_t fsad = fh ? 1/(72*fh) : 0;
  num_t   b0 = lw*m->knl[0]/abs(m->el)*m->sdir*m->edir*m->charge;
  num_t   c2 = b0*fh*2;

  FOR (i,m->npar) {
    P p(m,i);
    T   dpp = 1 + 2/m->beta*p.pt + sqr(p.pt);
    T    pz = sqrt(dpp - sqr(p.px) - sqr(p.py));
    T   _pz = 1/pz;
    T  _pz2 = sqr(_pz);
    T  relp = 1/sqrt(dpp);
    T  tfac = -1/m->beta - p.pt;
    T    c3 = sqr(b0)*fsad*relp;

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

    if (m->sdir == 1) {
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

template <typename P, typename T=P::T>
inline void qsad_fringe (mflw_t *m, num_t lw)
{
  if (abs(m->knl[1])+abs(m->ksl[1]) < minstr) return;
  if (abs(m->f1)    +abs(m->f2)     < minstr) return;
  mdump(0);
  num_t wchg = lw*m->charge;
  num_t  a   = -0.5*atan2(m->ksl[1], m->knl[1]);
  num_t b2   = hypot(m->knl[1], m->ksl[1])/m->el*m->edir;
  num_t ca   = cos(a), sa = sin(a);
  num_t bf1  = -abs(m->f1)*m->f1*b2/24;
  num_t bf2  =             m->f2*b2;

  // Lee-Whiting formula, E. Forest ch 13.2.3, eq 13.33
  FOR (i,m->npar) {
    P p(m,i);
    T _pz = 1/sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt));
    T  dt = (1/m->beta+p.pt)*_pz;

    T  f1 = wchg*bf1*_pz;
    T  f2 =      bf2*_pz;

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

template <typename P, typename T=P::T>
inline void mult_fringe (mflw_t *m, num_t lw)
{
  int n = MIN(m->nmul,m->fmax);
  if (!n) return;
  mdump(0);
  num_t    _l = m->el ? m->sdir*m->edir/m->el : m->edir;
  num_t  wchg = lw*m->charge;
  num_t no_k1 = m->frng & fringe_bend;

  FOR (i,m->npar) {
    P p(m,i);
    T rx, ix;       rx = 1., ix=0.;
    T fx, fxx, fxy; fx = 0., fxx=0., fxy=0.;
    T fy, fyy, fyx; fy = 0., fyy=0., fyx=0.;

    T _pz = 1/sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt));

    FOR (j,1,n+1) {
      T drx, dix; drx = rx, dix = ix;
      rx = drx*p.x - dix*p.y;
      ix = drx*p.y + dix*p.x;

      num_t nj = -wchg/(4*(j+1)), nf = (j+2)/j;
      num_t kj = m->knl[j-1]*_l, ksj = m->ksl[j-1]*_l;

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

    p.x  -= fx*_pz;
    p.y  -= fy*_pz;
    p.px  = (d*p.px - b*p.py)*_det;
    p.py  = (a*p.py - c*p.px)*_det;
    p.t  += (1/m->beta+p.pt)*(p.px*fx + p.py*fy)*sqr(_pz)*_pz;
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void curex_fringe (mflw_t *m, num_t lw)
{
  mdump(0);
  if (m->sdir*lw == 1) { // 'forward entry' or 'backward exit'
    yrotation<P> (m, 1,-m->e);
    bend_face<P> (m, 1, m->h);
    if (m->frng) {
      if (m->frng & fringe_bend) bend_fringe<P>(m, 1);
      if (m->frng & fringe_mult) mult_fringe<P>(m, 1);
      if (m->frng & fringe_qsad) qsad_fringe<P>(m, 1);
    }
    mad8_wedge<P>(m, 1, m->e);
    bend_wedge<P>(m, 1,-m->e);
  } else {
    bend_wedge<P>(m,-1,-m->e);
    mad8_wedge<P>(m,-1, m->e);
    if (m->frng) {
      if (m->frng & fringe_qsad) qsad_fringe<P>(m,-1);
      if (m->frng & fringe_mult) mult_fringe<P>(m,-1);
      if (m->frng & fringe_bend) bend_fringe<P>(m,-1);
    }
    bend_face<P> (m,-1, m->h);
    yrotation<P> (m,-1, m->e);
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void strex_fringe (mflw_t *m, num_t lw)
{
  mdump(0);
  num_t a = !m->pdir*(0.5*m->eh*(m->eld ? m->eld : m->el) - m->e);
  bool  p =  m->pdir == lw;

  if (m->sdir*lw == 1) { // 'forward entry' or 'backward exit'
    yrotation<P> (m, 1,-m->e  );
    bend_ptch<P> (m, 1, m->a*p);
    bend_face<P> (m, 1, m->h  );
    if (m->frng) {
      if (m->frng & fringe_bend) bend_fringe<P>(m, 1);
      if (m->frng & fringe_mult) mult_fringe<P>(m, 1);
      if (m->frng & fringe_qsad) qsad_fringe<P>(m, 1);
    }
    bend_wedge<P>(m, 1, a);
  } else {
    bend_wedge<P>(m,-1, a);
    if (m->frng) {
      if (m->frng & fringe_qsad) qsad_fringe<P>(m,-1);
      if (m->frng & fringe_mult) mult_fringe<P>(m,-1);
      if (m->frng & fringe_bend) bend_fringe<P>(m,-1);
    }
    bend_face<P> (m,-1, m->h  );
    bend_ptch<P> (m,-1, m->a*p);
    yrotation<P> (m,-1, m->e  );
  }
  mdump(1);
}

template <typename P, typename T=P::T>
inline void rfcav_fringe (mflw_t *m, num_t lw)
{
  mdump(0);
  if (m->sdir*lw == 1) { // 'forward entry' or 'backward exit'
    ensure(m->Tbak < 0, "inconsistent totalpath when entering rfcavity");
    if (m->Tbak == -1)
      m->Tbak = m->T, m->T = 1, m->lag -= m->freq*abs(m->el)*pi_clight/m->beta;
    else
      m->Tbak = m->T, m->T = 0;
    if (lw == -1 && m->T != m->Tbak) adjust_time<P>(m,  1);
    if (m->frng & fringe_rfcav)      cav_fringe<P> (m,  1);
  } else {
    if (m->frng & fringe_rfcav)      cav_fringe<P> (m, -1);
    if (lw == -1 && m->T != m->Tbak) adjust_time<P>(m, -1);
    m->T = m->Tbak, m->Tbak = -1;
  }
  mdump(1);
}

// --- specializations --------------------------------------------------------o

// --- patches ---

void mad_trk_xrotation_r (mflw_t *m, num_t lw, num_t phi) {
  xrotation<par_t>(m, lw, phi);
}
void mad_trk_xrotation_t (mflw_t *m, num_t lw, num_t phi) {
  xrotation<map_t>(m, lw, phi);
}

void mad_trk_yrotation_r (mflw_t *m, num_t lw, num_t the) {
  yrotation<par_t>(m, lw, the);
}
void mad_trk_yrotation_t (mflw_t *m, num_t lw, num_t the) {
  yrotation<map_t>(m, lw, the);
}

void mad_trk_srotation_r (mflw_t *m, num_t lw, num_t psi) {
  srotation<par_t>(m, lw, psi);
}
void mad_trk_srotation_t (mflw_t *m, num_t lw, num_t psi) {
  srotation<map_t>(m, lw, psi);
}

void mad_trk_translate_r (mflw_t *m, num_t lw) {
  srotation<par_t>(m, lw);
}
void mad_trk_translate_t (mflw_t *m, num_t lw) {
  srotation<map_t>(m, lw);
}

void mad_trk_changeref_r (mflw_t *m, num_t lw) {
  changeref<par_t>(m, lw);
}
void mad_trk_changeref_t (mflw_t *m, num_t lw) {
  changeref<map_t>(m, lw);
}

// --- misalignment ---

void mad_trk_misalign_r (mflw_t *m, num_t lw) {
  misalign<par_t>(m, lw);
}
void mad_trk_misalign_t (mflw_t *m, num_t lw) {
  misalign<map_t>(m, lw);
}

// -- fringe maps ---

void mad_trk_strex_fringe_r (mflw_t *m, num_t lw) {
  strex_fringe<par_t>(m, lw);
}
void mad_trk_strex_fringe_t (mflw_t *m, num_t lw) {
  strex_fringe<map_t>(m, lw);
}

void mad_trk_curex_fringe_r (mflw_t *m, num_t lw) {
  curex_fringe<par_t>(m, lw);
}
void mad_trk_curex_fringe_t (mflw_t *m, num_t lw) {
  curex_fringe<map_t>(m, lw);
}

void mad_trk_rfcav_fringe_r (mflw_t *m, num_t lw) {
  rfcav_fringe<par_t>(m, lw);
}
void mad_trk_rfcav_fringe_t (mflw_t *m, num_t lw) {
  rfcav_fringe<map_t>(m, lw);
}

// --- DKD straight ---

void mad_trk_strex_drift_r (mflw_t *m, num_t lw, int is) {
  strex_drift<par_t>(m,lw,is);
}
void mad_trk_strex_drift_t (mflw_t *m, num_t lw, int is) {
  strex_drift<map_t>(m,lw,is);
}

void mad_trk_strex_kick_r (mflw_t *m, num_t lw, int is) {
  strex_kick<par_t>(m,lw,is);
}
void mad_trk_strex_kick_t (mflw_t *m, num_t lw, int is) {
  strex_kick<map_t>(m,lw,is);
}

void mad_trk_strex_kickhs_r (mflw_t *m, num_t lw, int is) {
  strex_kickhs<par_t>(m,lw,is);
}
void mad_trk_strex_kickhs_t (mflw_t *m, num_t lw, int is) {
  strex_kickhs<map_t>(m,lw,is);
}

// --- DKD curved ---

void mad_trk_curex_drift_r (mflw_t *m, num_t lw, int is) {
  curex_drift<par_t>(m,lw,is);
}
void mad_trk_curex_drift_t (mflw_t *m, num_t lw, int is) {
  curex_drift<map_t>(m,lw,is);
}

void mad_trk_curex_kick_r (mflw_t *m, num_t lw, int is) {
  curex_kick<par_t>(m,lw,is);
}
void mad_trk_curex_kick_t (mflw_t *m, num_t lw, int is) {
  curex_kick<map_t>(m,lw,is);
}

// --- sbend ---

void mad_trk_sbend_thick_r (mflw_t *m, num_t lw, int is) {
  sbend_thick<par_t>(m,lw,is);
}
void mad_trk_sbend_thick_t (mflw_t *m, num_t lw, int is) {
  sbend_thick<map_t>(m,lw,is);
}

void mad_trk_sbend_kick_r (mflw_t *m, num_t lw, int is) {
  curex_kick<par_t>(m,lw,is,true);
}
void mad_trk_sbend_kick_t (mflw_t *m, num_t lw, int is) {
  curex_kick<map_t>(m,lw,is,true);
}

// --- rbend ---

void mad_trk_rbend_thick_r (mflw_t *m, num_t lw, int is) {
  rbend_thick<par_t>(m,lw,is);
}
void mad_trk_rbend_thick_t (mflw_t *m, num_t lw, int is) {
  rbend_thick<map_t>(m,lw,is);
}

void mad_trk_rbend_kick_r (mflw_t *m, num_t lw, int is) {
  strex_kick<par_t>(m,lw,is,true);
}
void mad_trk_rbend_kick_t (mflw_t *m, num_t lw, int is) {
  strex_kick<map_t>(m,lw,is,true);
}

// --- quadrupole ---

void mad_trk_quad_thick_r (mflw_t *m, num_t lw, int is) {
  quad_thick<par_t>(m,lw,is);
}
void mad_trk_quad_thick_t (mflw_t *m, num_t lw, int is) {
  quad_thick<map_t>(m,lw,is);
}

void mad_trk_quad_kick_r (mflw_t *m, num_t lw, int is) {
  quad_kick<par_t>(m,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kick_t (mflw_t *m, num_t lw, int is) {
  quad_kick<map_t>(m,lw,0); (void)is; // always yoshida
}

void mad_trk_quad_thicks_r (mflw_t *m, num_t lw, int is) {
  quad_thicks<par_t>(m,lw,is);
}
void mad_trk_quad_thicks_t (mflw_t *m, num_t lw, int is) {
  quad_thicks<map_t>(m,lw,is);
}

void mad_trk_quad_kicks_r (mflw_t *m, num_t lw, int is) {
  quad_kicks<par_t>(m,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kicks_t (mflw_t *m, num_t lw, int is) {
  quad_kicks<map_t>(m,lw,0); (void)is; // always yoshida
}

void mad_trk_quad_thickh_r (mflw_t *m, num_t lw, int is) {
  quad_thickh<par_t>(m,lw,is);
}
void mad_trk_quad_thickh_t (mflw_t *m, num_t lw, int is) {
  quad_thickh<map_t>(m,lw,is);
}

void mad_trk_quad_kickh_r (mflw_t *m, num_t lw, int is) {
  quad_kickh<par_t>(m,lw,0); (void)is; // always yoshida
}
void mad_trk_quad_kickh_t (mflw_t *m, num_t lw, int is) {
  quad_kickh<map_t>(m,lw,0); (void)is; // always yoshida
}

// --- solenoid ---

void mad_trk_solen_thick_r (mflw_t *m, num_t lw, int is) {
  solen_thick<par_t>(m,lw,is);
}
void mad_trk_solen_thick_t (mflw_t *m, num_t lw, int is) {
  solen_thick<map_t>(m,lw,is);
}

// --- eseptum ---

void mad_trk_esept_thick_r (mflw_t *m, num_t lw, int is) {
  esept_thick<par_t>(m,lw,is);
}
void mad_trk_esept_thick_t (mflw_t *m, num_t lw, int is) {
  esept_thick<map_t>(m,lw,is);
}

// --- rfcavity ---

void mad_trk_rfcav_kick_r (mflw_t *m, num_t lw, int is) {
  rfcav_kick<par_t>(m,lw,is);
}
void mad_trk_rfcav_kick_t (mflw_t *m, num_t lw, int is) {
  rfcav_kick<map_t>(m,lw,is);
}

void mad_trk_rfcav_kickn_r (mflw_t *m, num_t lw, int is) {
  rfcav_kickn<par_t>(m,lw,is);
}
void mad_trk_rfcav_kickn_t (mflw_t *m, num_t lw, int is) {
  rfcav_kickn<map_t>(m,lw,is);
}

// --- do nothing ---

void mad_trk_fnil (mflw_t *m, num_t lw, int is) {
  (void)m, (void)lw, (void)is;
}

// --- track one thick or thin ------------------------------------------------o

void mad_trk_slice_thk (mflw_t *m, num_t lw, trkfun *thick)
{
  thick(m, lw, 0);
}

void mad_trk_slice_thn (mflw_t *m, num_t lw, trkfun *kick)
{
  kick(m, lw, 0);
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
  FOR(i,n) {
    thick(m, lw*yosh[j].d[i  ], ++k);
     kick(m, lw*yosh[j].k[i  ], ++k);
  } thick(m, lw*yosh[j].d[--n], ++k);
  RFOR(i,n) {
     kick(m, lw*yosh[j].k[i  ], ++k);
    thick(m, lw*yosh[j].d[i  ], ++k);
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
  FOR(i,n) {
     kick(m, lw*boole[j].k[i], k++);
    thick(m, lw*boole[j].d   , k++);
  }  kick(m, lw*boole[j].k[n], k++); if (!n) ++n;
  RFOR(i,n) {
    thick(m, lw*boole[j].d   , k++);
     kick(m, lw*boole[j].k[i], k++);
  }
}

// --- track one Teapot slice -------------------------------------------------o

const struct {
  const num_t d;
  const num_t D;
  const num_t k;
} teapot[] = {
  {1./6 , 1./6 , 1./2 }, // teapot2 -> 0..0, knd=2, j=0, n=1, k=-2
  {1./8 , 3./8 , 1./3 }, // teapot3 -> 0..1, knd=3, j=1, n=2, k=-3
  {3./30, 8./30, 1./4 }, // teapot4 -> 0..2, knd=4, j=2, n=3, k=-4
} ;

void mad_trk_slice_tpt (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick, int knd)
{
  ensure(knd >= 2 && knd <= 4, "invalid teapot kind 2..4");
  int j = knd-2;
  int n = knd-1;
  int k = knd-1;
    thick(m, lw*teapot[j].d, k++);
  FOR(i,n) {
     kick(m, lw*teapot[j].k, k++);
    thick(m, lw*teapot[j].D, k++);
  }  kick(m, lw*teapot[j].k, k++);
    thick(m, lw*teapot[j].d, k++);
}

// --- speed tests ------------------------------------------------------------o

#if 1

#if TPSA_USE_TRC
#define TRC(...) printf(#__VA_ARGS__ "\n"); __VA_ARGS__
#else
#define TRC(...) __VA_ARGS__
#endif

void mad_trk_spdtest (int n, int k)
{
  mad_desc_newv(6, 1);

  tpsa x ( "X"); x .set( 0   , 1);
  tpsa px("PX"); px.set( 1e-7, 2);
  tpsa y ( "Y"); y .set( 0   , 3);
  tpsa py("PY"); py.set(-1e-7, 4);
  tpsa t ( "T"); t .set( 0   , 5);
  tpsa pt("PT"); pt.set( 0   , 6);

  num_t     par[] = { x[0], px[0], y[0], py[0], t[0], pt[0] };
  tpsa_t*   map[] = { x.ptr(), px.ptr(), y.ptr(), py.ptr(), t.ptr(), pt.ptr() };

  num_t   *pars[] = {par};
  tpsa_t* *maps[] = {map};

  struct mflw_ m = {
    .name="spdtest", .dbg=0,

    .el=1, .eld=1, .elc=0, .lrad=0,
    .eh=0, .ang=0, .mang=0,

    .pc=1, .beta=1, .betgam=0, .charge=1,

    .sdir=1, .edir=1, .pdir=0, .T=0, .Tbak=-1,

    .k1=0, .ks=0, .ksi=0, .volt=0, .freq=0, .lag=0, .nbsl=0,

    .frng=0, .fmax=2, .e=0, .fint=0, .h=0, .hgap=0, .f1=0, .f2=0, .a=0,

    .ca=0, .sa=0, .tlt=0,

    .nmul=1, .npha=0, .knl={1e-7}, .ksl={0}, .pnl={0}, .psl={0},
    .snm=0,  .bfx={0}   , .bfy={0},
    .npar=1, .par=pars  , .map=maps,

    .dx=0, .dy=0, .ds=0, .dthe=0, .dphi=0, .dpsi=0,

    .algn = {.rot=false, .trn=false,
    .dx=0, .dy=0, .ds=0, .dthe=0, .dphi=0, .dpsi=0},
  };

  switch(k) {
  case 0: {
    FOR(i,n) mad_trk_strex_drift_r (&m, 1, 1);
    par_t p(&m,0);
    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
            p.x, p.px, p.y, p.py, p.t, p.pt);
  } break;

  case 1: {
    FOR(i,n) mad_trk_strex_drift_t (&m, 1, 1);
    map_t p(&m,0);
    stdout << p.x << p.px << p.y << p.py << p.t << p.pt;
  } break;

  case 2: {
    FOR(i,n) {
      mad_trk_strex_drift_r (&m, 0.5, 1);
      mad_trk_strex_kick_r  (&m,   1, 1);
      mad_trk_strex_drift_r (&m, 0.5, 1);
    }
    par_t p(&m,0);
    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
            p.x, p.px, p.y, p.py, p.t, p.pt);
  } break;

  case 3: {
    FOR(i,n) {
      mad_trk_strex_drift_t (&m, 0.5, 1);
      mad_trk_strex_kick_t  (&m,   1, 1);
      mad_trk_strex_drift_t (&m, 0.5, 1);
    }
    map_t p(&m,0);
    stdout << p.x << p.px << p.y << p.py << p.t << p.pt;
  } break;

  case 4: {
    FOR(i,n) mad_trk_slice_tkt(&m, 1, mad_trk_strex_drift_r, mad_trk_strex_kick_r, 2);
    par_t p(&m,0);
    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
            p.x, p.px, p.y, p.py, p.t, p.pt);
  } break;

  case 5: {
    FOR(i,n) mad_trk_slice_tkt(&m, 1, mad_trk_strex_drift_t, mad_trk_strex_kick_t, 2);
    map_t p(&m,0);
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
#endif

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

#if 1
void mad_trk_cpptest (void)
{
  mad_desc_newv(6, 1);

  TRC(tpsa a("A");                             )
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
#endif