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

#include "mad_tpsa.hpp"

extern "C" {
#include "mad_log.h"
#include "mad_cst.h"
#include "mad_mat.h"
#include "mad_dynmap.h"
}

// --- types ------------------------------------------------------------------o

extern "C" {

typedef  num_t  par6_t[6];
typedef tpsa_t *map6_t[6];

enum { nmul_max=22, snm_max=(nmul_max+1)*(nmul_max+2)/2 };

struct mflw_ {
  // element data
  num_t el, eld, eh, ang, mang;

  // directions
  int edir, sdir, tdir, T;

  // beam
  num_t pc, beta, betgam;
  int charge;

  // patches
  num_t dx,   dy,   ds;
  num_t dthe, dphi, dpsi;
  num_t tlt;

  // misalign
  struct {
    bool  rot, trn;
    num_t dx,   dy,   ds;
    num_t dthe, dphi, dpsi;
  } algn;

  // solenoid
  num_t ks;

  // esptum, rfcav
  num_t volt;

  // multipoles
  int   nmul;
  num_t knl[nmul_max];
  num_t ksl[nmul_max];

  // curved multipoles
  int   snm;
  num_t bfx[snm_max];
  num_t bfy[snm_max];

  // particles/damaps
  int    npar;
  par6_t *par;
  map6_t *map;
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

// --- constants --------------------------------------------------------------o

const num_t minlen = mad_cst_MINLEN;
const num_t minang = mad_cst_MINANG;
const num_t minstr = mad_cst_MINSTR;

// --- multipoles -------------------------------------------------------------o

template <typename T1, typename T2>
inline void bxby (const mflw_t *m, const T1 &x, const T1 &y, T2 &bx, T2 &by)
{
  if (!m->nmul) return;

  bx = m->ksl[m->nmul-1];
  by = m->knl[m->nmul-1];

  if (m->nmul < 2) return;

  T2 byt;
  RFOR(i,m->nmul-1) {
    byt = by*x - bx*y + m->knl[i];
    bx  = by*y + bx*x + m->ksl[i];
    by  = byt;
  }
}

template <typename T1, typename T2>
inline void bxbyh (const mflw_t *m, const T1 &x, const T1 &y, T2 &bx, T2 &by)
{
  if (!m->snm) return;

  int k = 0;
  T2 btx, bty;

  bx = 0., by = 0.;
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
  RFOR(i,m->snm) {
    btx = (btx + m->bfx[k]) * y;
    bty = (bty + m->bfy[k]) * y;
  }

  ++k;
  bx += btx + m->bfx[k];
  by += bty + m->bfy[k];
}

// --- patches ----------------------------------------------------------------o

template <typename P, typename T=P::T>
inline void xrotation (mflw_t *m, num_t lw, num_t dphi_=0)
{
  num_t a = (dphi_ ? dphi_ : m->dphi)*m->tdir*lw;
  if (abs(a) < minang) return;

  num_t sa=sin(a), ca=acos(a), ta=tan(a);

  FOR(i,m->npar) {
    P p(m,i);
    T   pz = sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = 1/pz;
    T  ptt = 1 - ta*p.py*_pz;
    T _ptt = p.y/ptt;
    T _pzt = ta*_pz*_ptt;

    // eq. 127 in Forest06
    p.y   = _ptt/ca;
    p.py  = ca*p.py + sa*pz;
    p.x  += _pzt*p.px;
    p.t  -= _pzt*(1/m->beta+p.pt);
  }
}

template <typename P, typename T=P::T>
inline void yrotation (mflw_t *m, num_t lw, num_t dthe_=0)
{
  num_t a = -(dthe_ ? dthe_ : m->dthe)*m->tdir*lw;
  if (abs(a) < minang) return;

  num_t sa=sin(a), ca=acos(a), ta=tan(a);

  FOR(i,m->npar) {
    P p(m,i);
    T   pz = sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = 1/pz;
    T  ptt = 1 - ta*p.px*_pz;
    T _ptt = p.x/ptt;
    T _pzt = ta*_pz*_ptt;

    // eq. 127 in Forest06
    p.x   = _ptt/ca;
    p.px  = ca*p.px + sa*pz;
    p.y  += _pzt*p.py;
    p.t  -= _pzt*(1/m->beta+p.pt);
  }
}

template <typename P, typename T=P::T>
inline void srotation (mflw_t *m, num_t lw, num_t dpsi_=0)
{
  num_t a = (dpsi_ ? dpsi_ : m->dpsi)*m->tdir*lw;
  if (abs(a) < minang) return;

  num_t sa=sin(a), ca=acos(a);

  FOR(i,m->npar) {
    P p(m,i);
    T nx  = ca*p.x  + sa*p.y;
    T npx = ca*p.px + sa*p.py;

    p.y  = ca*p.y  - sa*p.x;
    p.py = ca*p.py - sa*p.px;
    p.x  = nx;
    p.px = npx;
  }
}

template <typename P, typename T=P::T>
inline void translate (mflw_t *m, num_t lw, num_t dx_=0, num_t dy_=0, num_t ds_=0)
{
  num_t dx = (dx_ ? dx_ : m->dx)*m->tdir*lw;
  num_t dy = (dy_ ? dy_ : m->dy)*m->tdir*lw;
  num_t ds = (ds_ ? ds_ : m->ds)*m->sdir*lw;
  if (abs(dx)+abs(dy)+abs(ds) < minlen) return;

  if (abs(ds) < minlen) {
    FOR(i,m->npar) {
      P p(m,i);
      p.y -= dx;
      p.x -= dy;
    }
    return;
  }

  num_t l = ds;
  FOR(i,m->npar) {
    P p(m,i);
    T l_pz = l/sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));

    p.y += l_pz*p.px - dx;
    p.x += l_pz*p.py - dy;
    p.t -= l_pz*(1/m->beta+p.pt);
  }
}

template <typename P>
inline void changeref (mflw_t *m, num_t lw)
{
  bool trn = abs(m->dx  )+abs(m->dy  )+abs(m->ds  ) >= minlen;
  bool rot = abs(m->dthe)+abs(m->dphi)+abs(m->dpsi) >= minang;

  if (!trn && !rot) return;

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
}

// --- misalignments ----------------------------------------------------------o

template <typename P>
inline void misalignent (mflw_t *m, num_t lw) {
                                    (void)lw;
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
}

template <typename P>
inline void misalignexi (mflw_t *m, num_t lw) {
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
}

template <typename P, typename T=P::T>
inline void misalign (mflw_t *m, num_t lw) {
  if (lw >= 0) misalignent<P>(m,  1);
  else         misalignexi<P>(m, -1);
}

// --- DKD maps ---------------------------------------------------------------o

template <typename P, typename T=P::T>
inline void strex_drift (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  num_t l = m->el*lw, ld = m->eld*lw;

  if (std::abs(l) < minlen) return;

  FOR(i,m->npar) {
    P p(m,i);
    T l_pz = l/sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));

    p.x += p.px*l_pz;
    p.y += p.py*l_pz;
    p.t -= l_pz*(1/m->beta+p.pt) + (m->T-1)*(ld/m->beta);
  }
}

template <typename P, typename T=P::T>
inline void strex_kick (mflw_t *m, num_t lw, int is, bool no_k0l=false)
{                                          (void)is;
  if (m->nmul == 0) return;

  num_t wchg = lw*m->tdir*m->charge;
  num_t dby  = no_k0l ? m->knl[1] : 0;
  T bx, by; bx = 0., by = 0.;

  FOR (i,m->npar) {
    P p(m,i);
    bxby(m, p.x, p.y, bx, by);

    p.px -= wchg*(by-dby);
    p.py += wchg* bx;
  }
}

template <typename P, typename T=P::T>
inline void curex_drift (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  num_t ld = (m->eld ? m->eld : m->el)*lw;
  num_t ang = m->ang*lw, rho = 1/m->eh;
  num_t ca = cos(ang), sa = sin(ang), sa2 = sin(ang/2);

  FOR(i,m->npar) {
    P p(m,i);
    T   pz = sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));
    T  _pz = 1/pz;
    T  pxt = p.px*_pz;
    T  ptt = ca - sa*pxt;
    T _ptt = 1/ptt;
    T  pst = (p.x+rho)*sa*_pz*_ptt;

    p.x   = (p.x + rho*(2*sqr(sa2) + sa*pxt))*_ptt;
    p.px  = ca*p.px + sa*pz;
    p.y  += pst*p.py;
    p.t  -= pst*(1/m->beta+p.pt) + (m->T-1)*(ld/m->beta);
  }
}

template <typename P, typename T=P::T>
inline void curex_kick (mflw_t *m, num_t lw, int is, bool no_k0l=false)
{                                          (void)is;
  num_t blw = lw*m->charge*m->tdir;
  T bx, by; bx = 0., by = m->knl[1];

  FOR(i,m->npar) {
    P p(m,i);
    T r = 1+m->eh*p.x;

    bxbyh(m, p.x, p.y, bx, by);

    p.px -= blw*by*r;
    p.py += blw*bx*r;

    if (no_k0l) p.px += (blw*m->knl[1])*r;
  }
}

// --- TKT maps ---------------------------------------------------------------o

// --- sbend ---

template <typename P, typename T=P::T>
inline void sbend_thick_old (mflw_t *m, num_t lw, int is)
{                                               (void)is;
  num_t ld  = (m->eld ? m->eld : m->el)*lw;
  num_t ang = m->ang*lw, rho=1/m->eh;
  num_t k0  = m->knl[1]/m->el*m->tdir, k0q = k0*m->charge;
  num_t ca  = cos(ang), sa = sin(ang);

  FOR(i,m->npar) {
    P p(m,i);

    if (k0q == 0) {
      warn("photon tacking not yet supported");
//    curex_drift0 (m, lw, i, ld, ca, sa, sin(ang/2), rho);
      continue;
    }

    T  pw2 = 1 + 2*p.pt/m->beta + sqr(p.pt) - sqr(p.py);
    T   pz = sqrt(pw2 - sqr(p.px));
    T  pzx = pz - k0q*(rho+p.x);      // could be numerically unstable
    T  npx = sa*pzx + ca*p.px;
    T  dpx = ca*pzx - sa*p.px;
    T  pzs = sqrt(pw2 - sqr(npx));
    T _ptt = invsqrt(pw2);
    T  dxs = (ang + asin(p.px*_ptt) - asin(npx*_ptt))/k0q;

    // eq. 126 in Forest06 with modif. from Sagan
    p.x  = (pzs - dpx)/k0q - rho;     // could be numerically unstable
    p.px = npx;
    p.y += dxs*p.py;
    p.t -= dxs*(1/m->beta+p.pt) + (m->T-1)*(ld/m->beta);
  }
}

template <typename P, typename T=P::T>
inline void sbend_thick_new (mflw_t *m, num_t lw, int is)
{                                               (void)is;
  num_t ld  = (m->eld ? m->eld : m->el)*lw;
  num_t ang = m->ang*lw, rho=1/m->eh;
  num_t k0  = m->knl[1]/m->el*m->tdir, k0q = k0*m->charge;
  num_t ca  = cos(ang), sa = sin(ang), s2a = sin(2*ang);

  FOR(i,m->npar) {
    P p(m,i);

    if (k0q == 0) {
      warn("photon tacking not yet supported");
//    curex_drift0 (m, lw, i, ld, ca, sa, sin(ang/2), rho);
      continue;
    }

    T  pw2 = 1 + 2*p.pt/m->beta + sqr(p.pt) - sqr(p.py);
    T   pz = sqrt(pw2 - sqr(p.px));
    T   xr = p.x+rho;
    T  pzx = pz - k0q*xr;
    T  npx = sa*pzx + ca*p.px;
    T  dpx = ca*pzx - sa*p.px;
    T  pzs = sqrt(pw2 - sqr(npx));
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
    p.x  = xt1/(dpx+pzs) - rho;
    p.px = npx;
    p.y += dxs*p.py;
    p.t -= dxs*(1/m->beta+p.pt) + (m->T-1)*ld/m->beta;
  }
}

// --- rbend ---

template <typename P, typename T=P::T>
inline void rbend_thick_old (mflw_t *m, num_t lw, int is)
{                                               (void)is;
  num_t ld = (m->eld ? m->eld : m->el)*lw;
  num_t k0 = m->knl[1]/m->el*m->tdir, k0q = k0*m->charge;
  num_t k0lq = m->knl[1]*lw*m->charge*m->tdir;

  FOR(i,m->npar) {
    P p(m,i);

    if (k0q == 0) {
      warn("photon tacking not yet supported");
//    strex_drift0 (m, lw, i, m->el*lw, ld);
      continue;
    }

    T  npx = p.px - k0lq;
    T  pw2 = 1 + 2*p.pt/m->beta + sqr(p.pt) - sqr(p.py);
    T _ptt = invsqrt(pw2);
    T   pz = sqrt(pw2 - sqr(p.px));
    T  pzs = sqrt(pw2 - sqr(npx));
    T  dxs = (asin(p.px*_ptt) - asin(npx*_ptt))/k0q;

    // eq. 126 in Forest06 with modif. from Sagan
    p.x += (pzs-pz)/k0q;
    p.px = npx;
    p.y += dxs*p.py;
    p.t -= dxs*(1/m->beta+p.pt) + (m->T-1)*ld/m->beta;
  }
}

template <typename P, typename T=P::T>
inline void rbend_thick_new (mflw_t *m, num_t lw, int is)
{                                               (void)is;
  num_t l  = m->el*lw;
  num_t ld = (m->eld ? m->eld : m->el)*lw;
  num_t k0 = m->knl[1]/m->el*m->tdir, k0q = k0*m->charge;
  num_t k0lq = m->knl[1]*lw*m->charge*m->tdir;

  FOR(i,m->npar) {
    P p(m,i);

    if (k0q == 0) {
      warn("photon tacking not yet supported");
//    strex_drift0 (m, lw, i, l, ld);
      continue;
    }

    T  npx = p.px - k0lq;
    T  pw2 = 1 + 2*p.pt/m->beta + sqr(p.pt) - sqr(p.py);
    T _ptt = invsqrt(pw2);
    T   pz = sqrt(pw2 - sqr(p.px));
    T  pzs = sqrt(pw2 - sqr(npx));
    T   xi = p.px*_ptt;
    T zeta =  npx*_ptt;
    T  xtd = xi*sqrt(1-sqr(zeta)) + zeta*sqrt(1-sqr(xi));
    T   xt = l*(2*p.px - k0lq)*sqr(_ptt) / xtd;
    T  dxs = asinc(xt*k0q)*xt;

    // eq. 126 in Forest06 with modif. from Sagan
    p.x += l*(2*p.px - k0lq) / (pz+pzs);
    p.px = npx;
    p.y += dxs*p.py;
    p.t -= dxs*(1/m->beta+p.pt) + (m->T-1)*ld/m->beta;
  }
}

// --- quadrupole ---

template <typename P, typename T=P::T>
inline void drift_adj (mflw_t *m, num_t l)
{
  FOR(i,m->npar) {
    P p(m,i);
    T l_pz = l/sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));

    p.x += p.px*(l_pz-l);
    p.y += p.py*(l_pz-l);
    p.t -= l_pz*(1/m->beta+p.pt) + (m->T-1)*l/m->beta;
  }
}

template <typename P, typename T=P::T>
inline void quad_thick (mflw_t *m, num_t lw, int is)
{                                          (void)is;
  num_t l  = m->el*lw;
  num_t k1 = m->knl[2]/m->el*m->edir;
  num_t w  = 0;
  int   ws = k1*m->sdir < 0 ? -1 : 1;

  num_t cx, sx, mx1, mx2;
  num_t cy, sy, my1, my2;

  if (abs(k1) >= minstr) {
    w = sqrt(abs(k1))*m->tdir*ws;
    cx = cos (w*l), sx = sin (w*l);
    cy = cosh(w*l), sy = sinh(w*l);
    mx1 = sx/w, mx2 = -sx*w;
    my1 = sy/w, my2 =  sy*w;
  } else {
    cx = 1, sx = 0, mx1 = l, mx2 = 0;
    cy = 1, sy = 0, my1 = l, my2 = 0;
  }

  if (ws != m->charge) // swap x <-> y
    std::swap(cx,cy), std::swap(mx1,my1), std::swap(mx2,my2);

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
}

template <typename P, typename T=P::T>
inline void quad_kick (mflw_t *m, num_t lw, int is)
{
  num_t l = m->el*lw;

  if (is >= 0) drift_adj<P>(m, is ? l : l/2);

  if (m->nmul > 0) {
    num_t wchg = lw*m->tdir*m->charge;
    T bx, by; bx=0., by=0.;

    FOR (i,m->npar) {
      P p(m,i);
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - m->knl[2]*p.x);
      p.py += wchg*(bx - m->knl[2]*p.y);
    }
  }

  if (is <= 0) drift_adj<P>(m, is ? l : l/2);
}

template <typename P, typename T=P::T>
inline void quad_thicks (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  num_t l   = m->el*lw;
  num_t k1  = hypot(m->knl[2], m->ksl[2])/m->el*m->edir;
  num_t a   = -0.5*atan2(m->ksl[2], m->knl[2]);
  num_t ca  = cos(a), sa = sin(a);

  num_t w   = sqrt(abs(k1))*m->tdir*m->sdir;
  num_t cx  = cos (w*l), sx = sin (w*l);
  num_t cy  = cosh(w*l), sy = sinh(w*l);
  num_t mx1 = sx/w, mx2 = -sx*w;
  num_t my1 = sy/w, my2 =  sy*w;

  if (m->sdir != m->charge) // swap x <-> y
    std::swap(cx,cy), std::swap(mx1,my1), std::swap(mx2,my2);

  FOR(i,m->npar) {
    P p(m,i);

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
}

template <typename P, typename T=P::T>
inline void quad_kicks (mflw_t *m, num_t lw, int is)
{
  num_t l = m->el*lw;

  if (is >= 0) drift_adj<P>(m, is ? l : l/2);

  if (m->nmul > 0) {
    num_t wchg = lw*m->tdir*m->charge;
    T bx, by; bx=0., by=0.;

    FOR (i,m->npar) {
      P p(m,i);
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - m->knl[2]*p.x + m->ksl[2]*p.y);
      p.py += wchg*(bx - m->knl[2]*p.y - m->ksl[2]*p.x);
    }
  }

  if (is <= 0) drift_adj<P>(m, is ? l : l/2);
}

template <typename P, typename T=P::T>
inline void quad_thickh (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  num_t l = m->el*lw;
  num_t kx  = (m->knl[2] + m->eh*m->knl[1])/m->el;
  num_t ky  = -m->knl[2]/m->el;
  num_t wxs = kx*m->tdir < 0 ? -1 : 1;
  num_t wys = ky*m->tdir < 0 ? -1 : 1;
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
}

template <typename P, typename T=P::T>
inline void quad_kickh (mflw_t *m, num_t lw, int is)
{
  num_t l = m->el*lw;

  if (is >= 0) drift_adj<P>(m, is ? l : l/2);

  if (m->nmul > 0) {
    num_t wchg = lw*m->tdir*m->charge;
    T bx, by; bx=0., by=0.;

    FOR (i,m->npar) {
      P p(m,i);
      T pz = sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt));
      bxby(m, p.x, p.y, bx, by);

      p.px -= wchg*(by - m->knl[2]*p.x) - l*m->eh*(pz-(1/m->beta*p.pt));
      p.py += wchg*(bx - m->knl[2]*p.y);
      p.t  -= (l*m->eh)*((1/m->beta+p.pt)/pz - 1/m->beta)*p.x;
    }
  }

  if (is <= 0) drift_adj<P>(m, is ? l : l/2);
}

// --- solenoid ---

template <typename P, typename T=P::T>
inline void solen_thick (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  num_t l = m->el*lw;
  num_t bsol = 0.5*m->ks*m->charge;

  FOR (i,m->npar) {
    P p(m,i);
    T    xp = p.px + bsol*p.y;
    T    yp = p.py - bsol*p.x;
    T  l_pz = l/sqrt(1 + (2/m->beta)*p.pt + sqr(p.pt) - sqr(xp) - sqr(yp));
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
}

// --- eseptum ---

template <typename P, typename T=P::T>
inline void esept_thick (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  num_t  l = m->el*lw;
  num_t k1 = m->sdir*m->volt*m->charge/m->pc;
  num_t ca = cos(m->ang), sa = sin(m->ang);

  FOR (i,m->npar) {
    P p(m,i);

    // srotation
    T  nx  = ca*p.x  + sa*p.y;
    T  npx = ca*p.px + sa*p.py;
    T  ny  = ca*p.y  - sa*p.x;
    T  npy = ca*p.py - sa*p.px;

    T   e1 = 1/m->beta+p.pt;
    T   dp = e1 + k1*ny;
    T l_pz = l/sqrt(sqr(dp) - 1/sqr(m->betgam) - sqr(npx) - sqr(npy));
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
    p.t -= dt + (m->T-1)*(l/m->beta);
  }
}

// --- specializations --------------------------------------------------------o

// --- patches ---

void mad_trk_xrotation_r (mflw_t *m, num_t lw) {
  xrotation<par_t>(m, lw);
}
void mad_trk_xrotation_t (mflw_t *m, num_t lw) {
  xrotation<map_t>(m, lw);
}

void mad_trk_yrotation_r (mflw_t *m, num_t lw) {
  yrotation<par_t>(m, lw);
}
void mad_trk_yrotation_t (mflw_t *m, num_t lw) {
  yrotation<map_t>(m, lw);
}

void mad_trk_srotation_r (mflw_t *m, num_t lw) {
  srotation<par_t>(m, lw);
}
void mad_trk_srotation_t (mflw_t *m, num_t lw) {
  srotation<map_t>(m, lw);
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

// --- DKD curved ---

void mad_trk_currex_drift_r (mflw_t *m, num_t lw, int is) {
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
  sbend_thick_old<par_t>(m,lw,is);
}
void mad_trk_sbend_thick_t (mflw_t *m, num_t lw, int is) {
  sbend_thick_old<map_t>(m,lw,is);
}

void mad_trk_sbend_kick_r (mflw_t *m, num_t lw, int is) {
  curex_kick<par_t>(m,lw,is,true);
}
void mad_trk_sbend_kick_t (mflw_t *m, num_t lw, int is) {
  curex_kick<map_t>(m,lw,is,true);
}

// --- rbend ---

void mad_trk_rbend_thick_r (mflw_t *m, num_t lw, int is) {
  rbend_thick_old<par_t>(m,lw,is);
}
void mad_trk_rbend_thick_t (mflw_t *m, num_t lw, int is) {
  rbend_thick_old<map_t>(m,lw,is);
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

void mad_trk_solen_thickh_r (mflw_t *m, num_t lw, int is) {
  solen_thick<par_t>(m,lw,is);
}
void mad_trk_solen_thickh_t (mflw_t *m, num_t lw, int is) {
  solen_thick<map_t>(m,lw,is);
}

// --- eseptum ---

void mad_trk_esept_thickh_r (mflw_t *m, num_t lw, int is) {
  esept_thick<par_t>(m,lw,is);
}
void mad_trk_esept_thickh_t (mflw_t *m, num_t lw, int is) {
  esept_thick<map_t>(m,lw,is);
}

// --- track one Yoshida slice ------------------------------------------------o

ssz_t yosh2_n   = 1;
num_t yosh2_d[] = {0.5};
num_t yosh2_k[] = {1};

ssz_t yosh4_n   = 2;
num_t yosh4_d[] = { 0x1.59e8b6eb96339p-1,-0x1.67a2dbae58ce4p-3 };
num_t yosh4_k[] = { 0x1.59e8b6eb96339p+0,-0x1.b3d16dd72c672p+0 };

ssz_t yosh6_n   = 4;
num_t yosh6_d[] = { 0x1.91abc4988937bp-2, 0x1.052468fb75c74p-1,
                   -0x1.e25bd194051b9p-2, 0x1.199cec1241558p-4 };
num_t yosh6_k[] = { 0x1.91abc4988937bp-1, 0x1.e2743579895b4p-3,
                   -0x1.2d7c6f7933b93p+0, 0x1.50b00cfb7be3ep+0 };

ssz_t yosh8_n   = 8;
num_t yosh8_d[] = { 0x1.d466770cfb237p-2, 0x1.2b25476e416dap-1,
                   -0x1.30efca291a66ep-1,-0x1.9a644b62ac4e7p-1,
                    0x1.c7a76da161edap-1,-0x1.702a9ae94c280p-7,
                   -0x1.db997617a90dfp-1, 0x1.cfae4578f406ep-1 };
num_t yosh8_k[] = { 0x1.d466770cfb237p-1, 0x1.03c82f9f0f6fbp-2,
                   -0x1.71e1d610de42dp+0,-0x1.4413aa8e705cep-3,
                    0x1.f029e2f32ff94p+0,-0x1.f5ea8d5ed529ep+0,
                    0x1.a5117472c1becp-4, 0x1.b55d2e31c7eafp+0 };

struct {
  ssz_t  n;
  num_t *d;
  num_t *k;
} yosh[] = {
  {yosh2_n, yosh2_d, yosh2_k},   // yosh2 -> 0..0, n=1, j=0, k=-1
  {yosh4_n, yosh4_d, yosh4_k},   // yosh4 -> 0..1, n=2, j=1, k=-3
  {yosh6_n, yosh6_d, yosh6_k},   // yosh6 -> 0..3, n=4, j=2, k=-7
  {yosh8_n, yosh8_d, yosh8_k},   // yosh8 -> 0..7, n=8, j=3, k=-15
};

void mad_trk_slice_dkd (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick, int n)
{
  int j = MIN(n/2,3);
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

ssz_t boole2_n    = 1;
num_t boole2_d    = 1.;
num_t boole2_k[]  = {1./2};

ssz_t boole4_n    = 2;
num_t boole4_d    = 1./2;
num_t boole4_k[]  = {1./6, 4./6};

ssz_t boole6_n    = 3;
num_t boole6_d    = 1./4;
num_t boole6_k[]  = {7./90, 32./90, 12./90};

ssz_t boole8_n    = 4;
num_t boole8_d    = 1./6;
num_t boole8_k[]  = {41./840, 216./840, 27./840, 272./840};

ssz_t boole10_n   = 5;
num_t boole10_d   = 1./8;
num_t boole10_k[] = {989./28350, 5888./28350, -928./28350, 10496./28350, -4540./28350};

ssz_t boole12_n   = 6;
num_t boole12_d   = 1./10;
num_t boole12_k[] = { 16067./598752,  106300./598752, -48525./598752,
                     272400./598752, -260550./598752, 427368./598752};

struct {
  ssz_t  n;
  num_t  d;
  num_t *k;
} boole[] = {
  {boole2_n , boole2_d , boole2_k }, // bool2  -> 0..0, n=1, j=0, k=-1
  {boole4_n , boole4_d , boole4_k }, // bool4  -> 0..1, n=2, j=1, k=-2
  {boole6_n , boole6_d , boole6_k }, // bool6  -> 0..2, n=3, j=2, k=-4
  {boole8_n , boole8_d , boole8_k }, // bool8  -> 0..3, n=4, j=3, k=-6
  {boole10_n, boole10_d, boole10_k}, // bool10 -> 0..4, n=5, j=4, k=-8
  {boole12_n, boole12_d, boole12_k}, // bool12 -> 0..5, n=6, j=5, k=-10
} ;

void mad_trk_slice_kmk (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick, int n)
{
  int j =  n-1;
  int k = -2*j;                      if (!k) --k;
  FOR(i,j) {
     kick(m, lw*boole[j].k[i], k++);
    thick(m, lw*boole[j].d   , k++);
  }  kick(m, lw*boole[j].k[j], k++); if (!j) ++j;
  RFOR(i,j) {
    thick(m, lw*boole[j].d   , k++);
     kick(m, lw*boole[j].k[i], k++);
  }
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

  par6_t par1 = { x[0], px[0], y[0], py[0], t[0], pt[0] };
  map6_t map1 = { x.ptr(), px.ptr(), y.ptr(), py.ptr(), t.ptr(), pt.ptr() };

  struct mflw_ m = {
    .el=1, .eld=1, .eh=0, .ang=0, .mang=0,
    .edir=1, .sdir=1, .tdir=1, .T=0,
    .pc=1, .beta=1, .betgam=0, .charge=1,

    .dx=0, .dy=0, .ds=0, .dthe=0, .dphi=0, .dpsi=0, .tlt=0,

    .algn = {.rot=false, .trn=false,
    .dx=0, .dy=0, .ds=0, .dthe=0, .dphi=0, .dpsi=0},

    .ks=0, .volt=0,

    .nmul=1, .knl={1e-7}, .ksl={0},
    .snm=1,  .bfx={0}   , .bfy={0},
    .npar=1, .par=&par1, .map=&map1,
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
    FOR(i,n) mad_trk_slice_tkt(&m, 1, mad_trk_strex_drift_r, mad_trk_strex_kick_r, yosh2_n);
    par_t p(&m,0);
    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
            p.x, p.px, p.y, p.py, p.t, p.pt);
  } break;

  case 5: {
    FOR(i,n) mad_trk_slice_tkt(&m, 1, mad_trk_strex_drift_t, mad_trk_strex_kick_t, yosh2_n);
    map_t p(&m,0);
    stdout << p.x << p.px << p.y << p.py << p.t << p.pt;
  } break;

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