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
//#include "mad_tpsa_notmp.hpp"

// --- types ------------------------------------------------------------------o

extern "C" {

struct mflw {
  num_t el, eld, beta, T;
  int    npar;
  num_t   *par[6];
  tpsa_t **map[6];
};

}

struct par_t {
  par_t(num_t *a, num_t b)
    : x(a[0]), px(a[1]), y(a[2]), py(a[3]), t(a[4]), pt(a[5]), beta(b) {}
  num_t x, px, y, py, t, pt, beta;
};

struct map_t {
  map_t(tpsa_t **a, num_t b)
    : x(a[0]), px(a[1]), y(a[2]), py(a[3]), t(a[4]), pt(a[5]), beta(b) {}
  mad::tpsa_ref x, px, y, py, t, pt;
  num_t beta;
};

// --- implementation ---------------------------------------------------------o

using namespace mad;

// --- helpers ---

const num_t minlen = mad_cst_MINLEN;
const num_t minang = mad_cst_MINANG;
const num_t minstr = mad_cst_MINSTR;

inline num_t dp2 (par_t &p) {
  return 1 + (2/p.beta)*p.pt + sqr(p.pt);
}

inline tpsa dp2 (map_t &p) {
  return 1 + (2/p.beta)*p.pt + sqr(p.pt);
}

inline num_t pz2 (par_t &p) {
  return 1 + (2/p.beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py);
}

inline tpsa pz2 (map_t &p) {
  return 1 + (2/p.beta)*p.pt + sqr(p.pt) - sqr(p.px) - sqr(p.py);
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
    par_t p { m->par[i], m->beta };
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
    map_t p { m->map[i], m->beta };
    const tpsa l_pz = l/pz(p);

    p.x += p.px*l_pz;
    p.y += p.py*l_pz;
    p.t -= l_pz*(1/p.beta+p.pt) + (1-T)*(ld/p.beta);
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

void mad_trk_test (int n)
{
  mad_desc_newv(6, 1);

  tpsa a(newt());    mad_tpsa_setvar(&*a, 10, 1, 0);
  tpsa b(newt());    mad_tpsa_setvar(&*b, 20, 1, 0);

//mad_tpsa_print(&*a, "A", 0,0,0);
//mad_tpsa_print(&*b, "B", 0,0,0);

  FOR(i,n) {
//    tpsa c = a+1+b+2+a+2;
    tpsa c = (a+1)*sqr(b+2)+a*2;
  }

/*
-----------------------------------------------------------
tpsa c = (a+1)*sqr(b+2)+a*2;  WITH TMP
-> a
mad_tpsa.hpp:116:        newt: void                   => t1
mad_tpsa.hpp: 60:   tpsa_tmp_: tpsa_t* 0x6000029c8008
-> b
mad_tpsa.hpp:116:        newt: void                   => t2
mad_tpsa.hpp: 60:   tpsa_tmp_: tpsa_t* 0x6000029c8078
-> a*2
mad_tpsa.hpp:325:   operator*: tpa,num
mad_tpsa.hpp: 68:    tpsa_ref: tpsa_t& 0x6000029c8008
mad_tpsa.hpp:288:   operator*: ref,num
mad_tpsa.hpp:126:        newt: ref                    => t3
mad_tpsa.hpp: 60:   tpsa_tmp_: tpsa_t* 0x6000029c80e8
-> b+2
mad_tpsa.hpp:206:   operator+: tpa,num
mad_tpsa.hpp: 68:    tpsa_ref: tpsa_t& 0x6000029c8078
mad_tpsa.hpp:189:   operator+: ref,num
mad_tpsa.hpp:126:        newt: ref                    => t4
mad_tpsa.hpp: 60:   tpsa_tmp_: tpsa_t* 0x6000029c8158
-> sqr(b+2)
mad_tpsa.hpp:501:         sqr: tpa
mad_tpsa.hpp:322:   operator*: tpa,tpa
mad_tpsa.hpp: 68:    tpsa_ref: tpsa_t& 0x6000029c8158
mad_tpsa.hpp: 68:    tpsa_ref: tpsa_t& 0x6000029c8158
mad_tpsa.hpp:281:   operator*: ref,ref
mad_tpsa.hpp:126:        newt: ref                    => t5
mad_tpsa.hpp: 60:   tpsa_tmp_: tpsa_t* 0x6000029c81c8
-> a+1
mad_tpsa.hpp:206:   operator+: tpa,num
mad_tpsa.hpp: 68:    tpsa_ref: tpsa_t& 0x6000029c8008
mad_tpsa.hpp:189:   operator+: ref,num
mad_tpsa.hpp:126:        newt: ref                    => t6
mad_tpsa.hpp: 60:   tpsa_tmp_: tpsa_t* 0x6000029c8238
-> (a+1)*sqr(b+2)
mad_tpsa.hpp:295:   operator*: tmp,tmp
mad_tpsa.hpp: 60:   tpsa_tmp_: tpsa_t* 0x6000029c8238 -> t6
-> (a+1)*sqr(b+2)+a*2
mad_tpsa.hpp:168:   operator+: tmp,tmp
mad_tpsa.hpp: 60:   tpsa_tmp_: tpsa_t* 0x6000029c8238 -> t6
-> dtor
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000029c81c8
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000029c8158
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000029c80e8
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000029c8238
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000029c8078
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000029c8008

-----------------------------------------------------------
tpsa c = (a+1)*sqr(b+2)+a*2;  WITHOUT TMP
-> a
mad_tpsa_notmp.hpp:109:        newt: void             => t1
-> b
mad_tpsa_notmp.hpp:109:        newt: void             => t2
-> a*2
mad_tpsa.hpp:211:   operator*: tpa,num
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e8008
mad_tpsa.hpp:202:   operator*: ref,num
mad_tpsa.hpp:119:        newt: ref                    => t3
-> b+2
mad_tpsa.hpp:157:   operator+: tpa,num
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e8078
mad_tpsa.hpp:147:   operator+: ref,num
mad_tpsa.hpp:119:        newt: ref                    => t4
-> sqr(b+2)
mad_tpsa.hpp:302:         sqr: tpa
mad_tpsa.hpp:208:   operator*: tpa,tpa
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e8158
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e8158
mad_tpsa.hpp:195:   operator*: ref,ref
mad_tpsa.hpp:119:        newt: ref                    => t5
-> a+1
mad_tpsa.hpp:157:   operator+: tpa,num
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e8008
mad_tpsa.hpp:147:   operator+: ref,num
mad_tpsa.hpp:119:        newt: ref                    => t6
-> (a+1)*sqr(b+2)
mad_tpsa.hpp:208:   operator*: tpa,tpa
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e81c8
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e8238
mad_tpsa.hpp:195:   operator*: ref,ref
mad_tpsa.hpp:119:        newt: ref                    => t7
-> (a+1)*sqr(b+2)+a*2
mad_tpsa.hpp:154:   operator+: tpa,tpa
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e80e8
mad_tpsa.hpp: 61:    tpsa_ref: tpsa_t& 0x6000020e82a8
mad_tpsa.hpp:140:   operator+: ref,ref
mad_tpsa.hpp:119:        newt: ref                    => t8
-> dtor
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000020e82a8
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000020e8238
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000020e81c8
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000020e8158
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000020e80e8
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000020e8318
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000020e8078
mad_tpsa.hpp: 51:  operator(): tpsa_t* 0x6000020e8008
*/
}

// ----------------------------------------------------------------------------o
