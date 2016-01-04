#include <stdio.h>
#include <assert.h>
#include "mad_tpsa.h"

#define T struct tpsa
#define num_t double

// #define TRACE

void
mad_track_drift(T * restrict m[], num_t L, num_t B, num_t E)
{
#ifdef TRACE
  printf("track_drift\n");
#endif

  assert(m);

  T *x = m[0], *px = m[1];
  T *y = m[2], *py = m[3];
  T *s = m[4], *ps = m[5];

  assert(x ); assert(y ); assert(s ); 
  assert(px); assert(py); assert(ps);

  T *t1 = mad_tpsa_new(x, mad_tpsa_default);
  T *t2 = mad_tpsa_new(s, mad_tpsa_same   );

  mad_tpsa_ax2pby2pcz2(1,ps, -1,px, -1,py, t1); // ps^2 - px^2 - py^2
  mad_tpsa_axpbypc(2/B,ps, 1,t1, 1, t1);        // 1 + 2/B*m.ps + ps^2 - px^2 - py^2
  mad_tpsa_invsqrt(t1,L, t1);                   // L/sqrt(1 + 2/B*m.ps + ps^2 - px^2 - py^2)

  T *l_pz = t1;

  mad_tpsa_axypbzpc(1,px,l_pz, 1,x, 0, x);      // x + px*l_pz -> x
  mad_tpsa_axypbzpc(1,py,l_pz, 1,y, 0, y);      // y + py*l_pz -> y

  mad_tpsa_copy(ps, t2);
  mad_tpsa_set0(t2,1,1/B);                      // ps + 1/B
  mad_tpsa_axypbzpc(1,t2,l_pz, 1,s, -E/B, s);   // s + (ps + 1/B)*l_pz -> s

  mad_tpsa_del(t2);
  mad_tpsa_del(t1);
}

void
mad_track_kick(T * restrict m[], num_t L, num_t B, int n, num_t Bn[n], num_t An[n])
{
#ifdef TRACE
  printf("track_drift\n");
#endif

  assert(m);

  T *x = m[0], *px = m[1];
  T *y = m[2], *py = m[3];
  T *s = m[4], *ps = m[5];

  assert(x ); assert(y ); assert(s ); 
  assert(px); assert(py); assert(ps);

  int dir = 1; // TODO: (m.dir or 1) * (m.charge or 1)

  T* bbxtw = mad_tpsa_new(px, mad_tpsa_same);
  T* bbytw = mad_tpsa_new(py, mad_tpsa_same);

  mad_tpsa_scalar(bbxtw, Bn[n-1]);
  mad_tpsa_scalar(bbytw, An[n-1]);

  if (n > 2) {
    T* bbytwt = mad_tpsa_new(py, mad_tpsa_same);

    for (int j = n-2; j >= 0; j--) {
      mad_tpsa_axypbvwpc(1,x,bbytw, -1,y,bbxtw, Bn[j], bbytwt);
      mad_tpsa_axypbvwpc(1,y,bbytw,  1,x,bbxtw, An[j], bbxtw);
      T* tmp = bbytw; bbytw = bbytwt; bbytwt = tmp;
    }

    mad_tpsa_del(bbytwt);
  }

  mad_tpsa_axpbypc(1,px, -L*dir,bbytw, 0, px); 
  mad_tpsa_axpbypc(1,py,  L*dir,bbxtw, 0, py);
  
  mad_tpsa_axypbzpc(1,ps,ps, 2/B,ps, 1, ps);
  mad_tpsa_sqrt(ps, ps);
  mad_tpsa_set0(ps,1,-1);

  mad_tpsa_del(bbytw);
  mad_tpsa_del(bbxtw);
}

#undef T
