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
#include "mad_cst.h"
#include "mad_dynmap.h"
}

// --- types ------------------------------------------------------------------o

extern "C" {

typedef  num_t  par6_t[6];
typedef tpsa_t *map6_t[6];

struct mflw {
  num_t el, eld, beta, T;
  int    npar;
  par6_t *par;
  map6_t *map;
};

}

struct par_t {
  par_t(num_t *a)
    : x(a[0]), px(a[1]), y(a[2]), py(a[3]), t(a[4]), pt(a[5]) {}
  num_t &x, &px, &y, &py, &t, &pt;
};

struct map_t {
  map_t(tpsa_t **a)
    : x(a[0]), px(a[1]), y(a[2]), py(a[3]), t(a[4]), pt(a[5]) {}
  mad::tpsa_ref x, px, y, py, t, pt;
};

// --- implementation ---------------------------------------------------------o

using namespace mad;

// --- constants ---

const num_t minlen = mad_cst_MINLEN;
const num_t minang = mad_cst_MINANG;
const num_t minstr = mad_cst_MINSTR;

// --- DKD maps ---

void
mad_trk_strex_drift_r (elem_t *e, mflw_t *m, num_t lw, int istp)
{
  (void)e; (void)istp;
  if (std::abs(m->el*lw) < minlen) return;

  num_t l = m->el*lw, ld = m->eld*lw, beta = m->beta;
  int T = m->T;

  FOR(i,m->npar) {
    par_t p { m->par[i] };
    num_t l_pz = l/sqrt(1 + (2/beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py));

    p.x += p.px*l_pz;
    p.y += p.py*l_pz;
    p.t -= l_pz*(1/beta+p.pt) - (1-T)*(ld/beta);
  }
}

void
mad_trk_strex_drift_t (elem_t *e, mflw_t *m, num_t lw, int istp)
{
  (void)e; (void)istp;
  if (std::abs(m->el*lw) < minlen) return;

  num_t l = m->el*lw, ld = m->eld*lw, beta = m->beta;
  int T = m->T;

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
    yrotation<P,T>(m, 1);
    xrotation<P,T>(m,-1);
    srotation<P,T>(m, 1);
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
  mad_mat_rtbar(rb, tb, abs(m->el), m->algn.ang, m->algn.tlt, m->algn.rot ? r:0, t);

  if (m->algn.rot && m->sdir > 0) {
    num_t v[3];
    mad_mat_torotxyz(rb, v, false);
    srotation<P,T>(m, -m->edir, v[2]);
    xrotation<P,T>(m,  m->edir, v[0]);
    yrotation<P,T>(m, -m->edir, v[1]);
  }

  if (m->algn.trn) translate<P,T>(m, -m->sdir, t[0], t[1], t[2]);

  if (m->algn.rot && m->sdir < 0) {
    num_t v[3];
    mad_mat_torotxyz(rb, v, false);
    yrotation<P,T>(m,  m->edir, v[1]);
    xrotation<P,T>(m, -m->edir, v[0]);
    srotation<P,T>(m,  m->edir, v[2]);
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
    p.t -= l_pz*(1/beta+p.pt) - (1-T)*(ld/beta);
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

// --- speed tests ------------------------------------------------------------o

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

  struct mflw m = {
    .el=1, .eld=1, .beta=1, .T=0,
    .npar=1, .par=&par1, .map=&map1,
  };

  switch(k) {
  case 0: {
    FOR(i,n) mad_trk_strex_drift_r (nullptr, &m, 1, 1);
    par_t p { m.par[0] };
    printf("x =% -.16e\npx=% -.16e\ny =% -.16e\npy=% -.16e\nt =% -.16e\npt=% -.16e\n",
            p.x, p.px, p.y, p.py, p.t, p.pt);
  } break;

  case 1: {
    FOR(i,n) mad_trk_strex_drift_t (nullptr, &m, 1, 1);
    map_t p { m.map[0] };
    stdout << p.x << p.px << p.y << p.py << p.t << p.pt;
  } break;

  default:
    printf("unknown use case %d\n", k);
  }
}

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
