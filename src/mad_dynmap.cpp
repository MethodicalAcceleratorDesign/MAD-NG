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
  num_t beta;
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
//par_t(num_t *a)
//  : x(a[0]), px(a[1]), y(a[2]), py(a[3]), t(a[4]), pt(a[5]) {}
  par_t(mflw_t *m, int i)
    : x(m->par[i][0]), px(m->par[i][1]),
      y(m->par[i][2]), py(m->par[i][3]),
      t(m->par[i][4]), pt(m->par[i][5]) {}
  num_t &x, &px, &y, &py, &t, &pt;
};

struct map_t {
//map_t(tpsa_t **a)
//  : x(a[0]), px(a[1]), y(a[2]), py(a[3]), t(a[4]), pt(a[5]) {}
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

template <typename P, typename T>
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

template <typename P, typename T>
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

template <typename P, typename T>
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

template <typename P, typename T>
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

template <typename P, typename T>
inline void changeref (mflw_t *m, num_t lw)
{
  bool trn = abs(m->dx  )+abs(m->dy  )+abs(m->ds  ) >= minlen;
  bool rot = abs(m->dthe)+abs(m->dphi)+abs(m->dpsi) >= minang;

  if (!trn && !rot) return;

  if (rot && lw > 0) {
    yrotation<P,T>(m,  1);
    xrotation<P,T>(m, -1);
    srotation<P,T>(m,  1);
  }

  if (trn) translate<P,T>(m, lw);

  if (rot && lw < 0) {
    srotation<P,T>(m, -1);
    xrotation<P,T>(m,  1);
    yrotation<P,T>(m, -1);
  }
}

// --- misalignments ----------------------------------------------------------o

template <typename P, typename T>
inline void misalignent (mflw_t *m, num_t lw) {
                                    (void)lw;
  if (m->algn.rot && m->sdir > 0) {
    yrotation<P,T>(m,  m->edir, m->algn.dthe);
    xrotation<P,T>(m, -m->edir, m->algn.dphi);
    srotation<P,T>(m,  m->edir, m->algn.dpsi);
  }

  if (m->algn.trn)
    translate<P,T>(m, m->sdir, m->algn.dx, m->algn.dy, m->algn.ds);

  if (m->algn.rot && m->sdir < 0) {
    srotation<P,T>(m, -m->edir, m->algn.dpsi);
    xrotation<P,T>(m,  m->edir, m->algn.dphi);
    yrotation<P,T>(m, -m->edir, m->algn.dthe);
  }
}

template <typename P, typename T>
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
    srotation<P,T>(m, -m->edir, -v[2]);
    xrotation<P,T>(m,  m->edir, -v[0]);
    yrotation<P,T>(m, -m->edir, -v[1]);
  }

  if (m->algn.trn) translate<P,T>(m, -m->sdir, tb[0], tb[1], tb[2]);

  if (m->algn.rot && m->sdir < 0) {
    num_t v[3];
    mad_mat_torotyxz(rb, v, true);
    yrotation<P,T>(m,  m->edir, -v[1]);
    xrotation<P,T>(m, -m->edir, -v[0]);
    srotation<P,T>(m,  m->edir, -v[2]);
  }
}

template <typename P, typename T>
inline void misalign (mflw_t *m, num_t lw) {
  if (lw >= 0) misalignent<P,T>(m,  1);
  else         misalignexi<P,T>(m, -1);
}

// --- DKD maps ---------------------------------------------------------------o

template <typename P, typename T>
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

template <typename P, typename T>
inline void strex_kick (mflw_t *m, num_t lw, int is, bool no_k0l=false)
{                                          (void)is;
  if (m->nmul == 0) return;

  num_t wchg = lw*m->tdir*m->charge;
  num_t dby  = no_k0l ? m->knl[1] : 0;
  T bx=0, by=0;

  FOR (i,m->npar) {
    P p(m,i);
    bxby(m, p.x, p.y, bx, by);

    p.px -= wchg*(by-dby);
    p.py += wchg* bx;
  }
}

template <typename P, typename T>
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

template <typename P, typename T>
inline void curex_kick (mflw_t *m, num_t lw, int is, bool no_k0l=false)
{                                          (void)is;
  num_t blw = lw*m->charge*m->tdir;
  T bx = 0, by = m->knl[1];

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

template <typename P, typename T>
inline void sbend_thick (mflw_t *m, num_t lw, int is)
{                                           (void)is;
  num_t ld = (m->eld ? m->eld : m->el)*lw;
  num_t ang = m->ang*lw, rho=1/m->eh;
  num_t k0 = m->knl[1]/m->el*m->tdir, k0q = k0*m->charge;
  num_t ca = cos(ang), sa = sin(ang), s2a = sin(2*ang);

  FOR(i,m->npar) {
    P p(m,i);

    if (k0q == 0) {
      warn("photon tacking not yet supported");
//    curex_drift0 (elm, m, lw, i, ld, ca, sa, sin(ang/2), rho)
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
    p.x   = xt1/(dpx+pzs) - rho;
    p.px  = npx;
    p.y  += dxs*p.py;
    p.t  -= dxs*(1/m->beta+p.pt) + (m->T-1)*(ld/m->beta);
  }
}

// --- specializations --------------------------------------------------------o

// --- patches ---

void mad_trk_xrotation_r (mflw_t *m, num_t lw) {
  xrotation<par_t, num_t>(m, lw);
}
void mad_trk_xrotation_t (mflw_t *m, num_t lw) {
  xrotation<map_t, tpsa>(m, lw);
}

void mad_trk_yrotation_r (mflw_t *m, num_t lw) {
  yrotation<par_t, num_t>(m, lw);
}
void mad_trk_yrotation_t (mflw_t *m, num_t lw) {
  yrotation<map_t, tpsa>(m, lw);
}

void mad_trk_srotation_r (mflw_t *m, num_t lw) {
  srotation<par_t, num_t>(m, lw);
}
void mad_trk_srotation_t (mflw_t *m, num_t lw) {
  srotation<map_t, tpsa>(m, lw);
}

void mad_trk_translate_r (mflw_t *m, num_t lw) {
  srotation<par_t, num_t>(m, lw);
}
void mad_trk_translate_t (mflw_t *m, num_t lw) {
  srotation<map_t, tpsa>(m, lw);
}

void mad_trk_changeref_r (mflw_t *m, num_t lw) {
  changeref<par_t, num_t>(m, lw);
}
void mad_trk_changeref_t (mflw_t *m, num_t lw) {
  changeref<map_t, tpsa>(m, lw);
}

// --- misalignment ---

void mad_trk_misalign_r (mflw_t *m, num_t lw) {
  misalign<par_t, num_t>(m, lw);
}
void mad_trk_misalign_t (mflw_t *m, num_t lw) {
  misalign<par_t, num_t>(m, lw);
}

// --- DKD straight ---

void mad_trk_strex_drift_r (mflw_t *m, num_t lw, int is) {
  strex_drift<par_t, num_t>(m,lw,is);
}
void mad_trk_strex_drift_t (mflw_t *m, num_t lw, int is) {
  strex_drift<map_t, tpsa>(m,lw,is);
}

void mad_trk_strex_kick_r (mflw_t *m, num_t lw, int is) {
  strex_kick<par_t, num_t>(m,lw,is);
}
void mad_trk_strex_kick_t (mflw_t *m, num_t lw, int is) {
  strex_kick<map_t, tpsa>(m,lw,is);
}

// --- DKD curved ---

void mad_trk_currex_drift_r (mflw_t *m, num_t lw, int is) {
  curex_drift<par_t, num_t>(m,lw,is);
}
void mad_trk_curex_drift_t (mflw_t *m, num_t lw, int is) {
  curex_drift<map_t, tpsa>(m,lw,is);
}

void mad_trk_curex_kick_r (mflw_t *m, num_t lw, int is) {
  curex_kick<par_t, num_t>(m,lw,is);
}
void mad_trk_curex_kick_t (mflw_t *m, num_t lw, int is) {
  curex_kick<map_t, tpsa>(m,lw,is);
}

// --- sbend ---

void mad_trk_sbend_thick_r (mflw_t *m, num_t lw, int is) {
  sbend_thick<par_t, num_t>(m,lw,is);
}
void mad_trk_sbend_thick_t (mflw_t *m, num_t lw, int is) {
  sbend_thick<map_t, tpsa>(m,lw,is);
}

// --- track Slice ------------------------------------------------------------o

void mad_trk_slice_dkd (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick,
                        int n, num_t *yosh_d, num_t *yosh_k)
{
  // yosh2 -> 0..0, n=1, k=-1
  // yosh4 -> 0..1, n=2, k=-3
  // yosh6 -> 0..3, n=4, k=-7
  // yosh8 -> 0..7, n=8, k=-15

  int k=-2*n;
  FOR(i,n) {
    thick(m, lw*yosh_d[i  ], ++k);
     kick(m, lw*yosh_k[i  ], ++k);
  } thick(m, lw*yosh_d[n-1], ++k);
  RFOR(i,n-1) {
     kick(m, lw*yosh_k[i  ], ++k);
    thick(m, lw*yosh_d[i  ], ++k);
  }
}

void mad_trk_slice_tkt (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick,
                        int n, num_t *yosh_d, num_t *yosh_k) {
  mad_trk_slice_dkd(m, lw, thick, kick, n, yosh_d, yosh_k);
}

void mad_trk_slice_kmk (mflw_t *m, num_t lw, trkfun *thick, trkfun *kick,
                        int n, num_t bool_d, num_t *bool_k)
{
  // bool2  -> 0..0, n=1, k=-1
  // bool4  -> 0..1, n=2, k=-2
  // bool6  -> 0..2, n=3, k=-4
  // bool8  -> 0..3, n=4, k=-6
  // bool10 -> 0..4, n=5, k=-8
  // bool12 -> 0..5, n=6, k=-10

  int k=-2*(n-1);                  if (n==1) --k;
  FOR(i,n-1) {
     kick(m, lw*bool_k[i  ], k++);
    thick(m, lw*bool_d     , k++);
  }  kick(m, lw*bool_k[n-1], k++); if (n==1) ++n;
  RFOR(i,n-1) {
    thick(m, lw*bool_d     , k++);
     kick(m, lw*bool_k[i  ], k++);
  }
}

// --- speed tests ------------------------------------------------------------o

#if 1

#if TPSA_USE_TRC
#define TRC(...) printf(#__VA_ARGS__ "\n"); __VA_ARGS__
#else
#define TRC(...) __VA_ARGS__
#endif

ssz_t yosh2_n   = 1;
num_t yosh2_d[] = {0.5};
num_t yosh2_k[] = {1};

ssz_t yosh4_n   = 2;
num_t yosh4_d[] = { 6.7560359597982889e-01, -1.7560359597982889e-01 };
num_t yosh4_k[] = { 1.3512071919596578e+00, -1.7024143839193155e+00 };

ssz_t yosh6_n   = 4;
num_t yosh6_d[] = { 3.9225680523877998e-01,  5.1004341191845848e-01,
                   -4.7105338540975655e-01,  6.8753168252518093e-02 };
num_t yosh6_k[] = { 7.8451361047755996e-01,  2.3557321335935699e-01,
                   -1.1776799841788701e+00,  1.3151863206839063e+00 };

ssz_t yosh8_n   = 8;
num_t yosh8_d[] = { 4.5742212311487002e-01,  5.8426879139798449e-01,
                   -5.9557945014712543e-01, -8.0154643611436149e-01,
                    8.8994925112725842e-01, -1.1235547676365032e-02,
                   -9.2890519179175246e-01,  9.0562646008949144e-01 };
num_t yosh8_k[] = { 9.1484424622974003e-01,  2.5369333656622900e-01,
                   -1.4448522368604799e+00, -1.5824063536824301e-01,
                    1.9381391376227599e+00, -1.9606102329754900e+00,
                    1.0279984939198500e-01,  1.7084530707869978e+00 };

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
    .beta=1, .charge=1,

    .dx=0, .dy=0, .ds=0, .dthe=0, .dphi=0, .dpsi=0, .tlt=0,

    .algn = {.rot=false, .trn=false,
    .dx=0, .dy=0, .ds=0, .dthe=0, .dphi=0, .dpsi=0},

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
    FOR(i,n)
      mad_trk_slice_tkt(&m, 1, mad_trk_strex_drift_r, mad_trk_strex_kick_r, yosh2_n, yosh2_d, yosh2_k);
    par_t p(&m,0);
    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
            p.x, p.px, p.y, p.py, p.t, p.pt);
  } break;

  case 5: {
    FOR(i,n)
      mad_trk_slice_tkt(&m, 1, mad_trk_strex_drift_t, mad_trk_strex_kick_t, yosh2_n, yosh2_d, yosh2_k);
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

#if 0
void mad_trk_cpptest (void)
{
  mad_desc_newv(6, 1);

  TRC(tpsa a("A");                             )
  TRC(tpsa_ref ar(a.ptr());                    )

  TRC( a  = a;                                 )
  TRC( a += a;                                 )
  TRC( a  = ar;                                )
  TRC( a += ar;                                )
  TRC( a  = a+a;                               )
  TRC( a += a+a;                               )
  TRC( a  = 2*a;                               )
  TRC( a += 2*a;                               )
  TRC( a  = 1;                                 )
  TRC( a += 1;                                 )
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
  TRC( ar  = 1;                                )
  TRC( ar += 1;                                )
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
  TRC( tpsa()  = 1;                            )
  TRC( tpsa() += 1;                            )
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
  TRC( tpsa_ref(a.ptr())  = 1;                 )
  TRC( tpsa_ref(a.ptr()) += 1;                 )
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