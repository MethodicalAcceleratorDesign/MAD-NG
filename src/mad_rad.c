/*
 o-----------------------------------------------------------------------------o
 |
 | Radiation module implementation
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

#include <math.h>

#include "mad_log.h"
#include "mad_rad.h"

// --- InvSynFracInt ----------------------------------------------------------o

static inline num_t
Chebyshev (num_t a, num_t b, const num_t c[], int m, num_t x)
{
  num_t y = (2*x-a-b)/(b-a);    // Change of variable.
  num_t y2 = 2*y, d = 0, dd = 0, sv;

  for (int j=m-1; j>=1; j--) {  // Clenshaw's recurrence.
    sv = d;
    d = y2*d-dd+c[j];
    dd = sv;
  }

  return y*d - dd + 0.5*c[0];
}

num_t
mad_rad_InvSynFracInt (num_t x)
{
  ensure(x >= 0 && x < 1, "invalid argument #1 (0 <= x < 1 expected)");

  // from 0 to 0.7
  static const num_t aa1=0, aa2=0.7;
  static const num_t cheb1[] = {
    1.22371665676046468821,0.108956475422163837267,0.0383328524358594396134,0.00759138369340257753721,
    0.00205712048644963340914,0.000497810783280019308661,0.000130743691810302187818,0.0000338168760220395409734,
    8.97049680900520817728e-6,2.38685472794452241466e-6,6.41923109149104165049e-7,1.73549898982749277843e-7,
    4.72145949240790029153e-8,1.29039866111999149636e-8,3.5422080787089834182e-9,9.7594757336403784905e-10,
    2.6979510184976065731e-10,7.480422622550977077e-11,2.079598176402699913e-11,5.79533622220841193e-12,
    1.61856011449276096e-12,4.529450993473807e-13,1.2698603951096606e-13,3.566117394511206e-14,1.00301587494091e-14,
    2.82515346447219e-15,7.9680747949792e-16};

  // from 0.7 to 0.9132260271183847
  static const num_t aa3=0.9132260271183847;
  static const num_t cheb2[] = {
    1.1139496701107756,0.3523967429328067,0.0713849171926623,0.01475818043595387,0.003381255637322462,
    0.0008228057599452224,0.00020785506681254216,0.00005390169253706556,0.000014250571923902464,3.823880733161044e-6,
    1.0381966089136036e-6,2.8457557457837253e-7,7.86223332179956e-8,2.1866609342508474e-8,6.116186259857143e-9,
    1.7191233618437565e-9,4.852755117740807e-10,1.3749966961763457e-10,3.908961987062447e-11,1.1146253766895824e-11,
    3.1868887323415814e-12,9.134319791300977e-13,2.6211077371181566e-13,7.588643377757906e-14,2.1528376972619e-14,
    6.030906040404772e-15,1.9549163926819867e-15};

  // Chebyshev with exp/log  scale
  // a = -Log[1 - SynFracInt[1]]; b = -Log[1 - SynFracInt[7]];
  static const num_t aa4=2.4444485538746025480, aa5=9.3830728608909477079;
  static const num_t cheb3[] = {
    1.2292683840435586977,0.160353449247864455879,-0.0353559911947559448721,0.00776901561223573936985,
    -0.00165886451971685133259,0.000335719118906954279467,-0.0000617184951079161143187,9.23534039743246708256e-6,
    -6.06747198795168022842e-7,-3.07934045961999778094e-7,1.98818772614682367781e-7,-8.13909971567720135413e-8,
    2.84298174969641838618e-8,-9.12829766621316063548e-9,2.77713868004820551077e-9,-8.13032767247834023165e-10,
    2.31128525568385247392e-10,-6.41796873254200220876e-11,1.74815310473323361543e-11,-4.68653536933392363045e-12,
    1.24016595805520752748e-12,-3.24839432979935522159e-13,8.44601465226513952994e-14,-2.18647276044246803998e-14,
    5.65407548745690689978e-15,-1.46553625917463067508e-15,3.82059606377570462276e-16,-1.00457896653436912508e-16};

  static const num_t aa6=33.122936966163038145;
  static const num_t cheb4[] = {
    1.69342658227676741765,0.0742766400841232319225,-0.019337880608635717358,0.00516065527473364110491,
    -0.00139342012990307729473,0.000378549864052022522193,-0.000103167085583785340215,0.0000281543441271412178337,
    -7.68409742018258198651e-6,2.09543221890204537392e-6,-5.70493140367526282946e-7,1.54961164548564906446e-7,
    -4.19665599629607704794e-8,1.13239680054166507038e-8,-3.04223563379021441863e-9,8.13073745977562957997e-10,
    -2.15969415476814981374e-10,5.69472105972525594811e-11,-1.48844799572430829499e-11,3.84901514438304484973e-12,
    -9.82222575944247161834e-13,2.46468329208292208183e-13,-6.04953826265982691612e-14,1.44055805710671611984e-14,
    -3.28200813577388740722e-15,6.96566359173765367675e-16,-1.294122794852896275e-16};

  enum { // avoid array size error
    ncheb1 = sizeof cheb1/sizeof cheb1[0],
    ncheb2 = sizeof cheb2/sizeof cheb2[0],
    ncheb3 = sizeof cheb3/sizeof cheb3[0],
    ncheb4 = sizeof cheb4/sizeof cheb4[0] };

  if(x<aa2)        return x*x*x*Chebyshev(aa1,aa2,cheb1,ncheb1,x);
  if(x<aa3)        return       Chebyshev(aa2,aa3,cheb2,ncheb2,x);

  num_t y = -log1p(-x);
  if(x<1-0.0000841363) return y*Chebyshev(aa4,aa5,cheb3,ncheb3,y);
  else                 return y*Chebyshev(aa5,aa6,cheb4,ncheb4,y);
}

#if 0
// Obsolete code not used, adapted from Placet (courtesy A. Latina)

#include <float.h>
#include <time.h>
#include <assert.h>

#include "mad_cst.h"
#include "mad_num.h"

// --- RNG for radiation ------------------------------------------------------o

static __thread rng_state_t rng[1];
static __thread num_t seed = 0;

static void
rngseed (num_t val)
{
  mad_num_randseed(rng, val);
  mad_num_randjump(rng);
  seed = val;
}

static num_t
rnguni (void)
{
  if (!seed) rngseed(12345678.9*time(NULL));
  return mad_num_rand(rng); // [0.,1.)
}

static num_t
rngexp (num_t mu)
{
  return -mu * log1p(-rnguni());
}

// --- SynRad -----------------------------------------------------------------o

//  x :    energy normalized to the critical energy
//  returns function value SynRad photon spectrum dn/dx
//  (integral of modified 1/3 order Bessel function)
//  Reference:
//  - Chebyshev series see H.H.Umstaetter CERN/PS/SM/81-13 10-3-1981
//  - H. Burkhardt LEP Note 632 of 12-1990
//  - H. Burkhardt CLIC-Note-709 of 08-06-2007

static num_t
SynRad(num_t x)
{
  num_t synrad=0;

  if (x > 0 && x < 800) {  // otherwise result synrad remains 0
    if (x < 6) {
      num_t a,b,z;
      z=x*x/16-2;
      b=        0.00000000000000000012;
      a=z*b  +  0.00000000000000000460;
      b=z*a-b+  0.00000000000000031738;
      a=z*b-a+  0.00000000000002004426;
      b=z*a-b+  0.00000000000111455474;
      a=z*b-a+  0.00000000005407460944;
      b=z*a-b+  0.00000000226722011790;
      a=z*b-a+  0.00000008125130371644;
      b=z*a-b+  0.00000245751373955212;
      a=z*b-a+  0.00006181256113829740;
      b=z*a-b+  0.00127066381953661690;
      a=z*b-a+  0.02091216799114667278;
      b=z*a-b+  0.26880346058164526514;
      a=z*b-a+  2.61902183794862213818;
      b=z*a-b+ 18.65250896865416256398;
      a=z*b-a+ 92.95232665922707542088;
      b=z*a-b+308.15919413131586030542;
      a=z*b-a+644.86979658236221700714;

      num_t p;
      p=0.5*z*a-b+414.56543648832546975110;
      a=       0.00000000000000000004;
      b=z*a+   0.00000000000000000289;
      a=z*b-a+ 0.00000000000000019786;
      b=z*a-b+ 0.00000000000001196168;
      a=z*b-a+ 0.00000000000063427729;
      b=z*a-b+ 0.00000000002923635681;
      a=z*b-a+ 0.00000000115951672806;
      b=z*a-b+ 0.00000003910314748244;
      a=z*b-a+ 0.00000110599584794379;
      b=z*a-b+ 0.00002581451439721298;
      a=z*b-a+ 0.00048768692916240683;
      b=z*a-b+ 0.00728456195503504923;
      a=z*b-a+ 0.08357935463720537773;
      b=z*a-b+ 0.71031361199218887514;
      a=z*b-a+ 4.26780261265492264837;
      b=z*a-b+17.05540785795221885751;
      a=z*b-a+41.83903486779678800040;

      num_t q;
      q=0.5*z*a-b+28.41787374362784178164;

      num_t y;
      y=pow(x,2.0/3);
      synrad=(p/y-q*y-1)*1.81379936423421784215530788143;
    } else { // 6 < x < 174
      num_t a,b,z;
      z=20/x-2;
      a=      0.00000000000000000001;
      b=z*a  -0.00000000000000000002;
      a=z*b-a+0.00000000000000000006;
      b=z*a-b-0.00000000000000000020;
      a=z*b-a+0.00000000000000000066;
      b=z*a-b-0.00000000000000000216;
      a=z*b-a+0.00000000000000000721;
      b=z*a-b-0.00000000000000002443;
      a=z*b-a+0.00000000000000008441;
      b=z*a-b-0.00000000000000029752;
      a=z*b-a+0.00000000000000107116;
      b=z*a-b-0.00000000000000394564;
      a=z*b-a+0.00000000000001489474;
      b=z*a-b-0.00000000000005773537;
      a=z*b-a+0.00000000000023030657;
      b=z*a-b-0.00000000000094784973;
      a=z*b-a+0.00000000000403683207;
      b=z*a-b-0.00000000001785432348;
      a=z*b-a+0.00000000008235329314;
      b=z*a-b-0.00000000039817923621;
      a=z*b-a+0.00000000203088939238;
      b=z*a-b-0.00000001101482369622;
      a=z*b-a+0.00000006418902302372;
      b=z*a-b-0.00000040756144386809;
      a=z*b-a+0.00000287536465397527;
      b=z*a-b-0.00002321251614543524;
      a=z*b-a+0.00022505317277986004;
      b=z*a-b-0.00287636803664026799;
      a=z*b-a+0.06239591359332750793;

      num_t p;
      p=0.5*z*a-b+1.06552390798340693166;
      synrad=p*sqrt(M_PI_2/x)/exp(x);
    }
  }
  return synrad;
}

struct syngen {
  num_t last_xmin, xlow, ratio, a1, a2, c1;
};

static void
syngen_init (struct syngen *sg, num_t xmin)
{
  sg->last_xmin = xmin;
  sg->xlow = xmin > 1 ? xmin : 1;

  // initialize constants used in the approximate expressions
  // for SYNRAD (integral over the modified Bessel function K5/3)
  sg->a1 = SynRad(1e-38) / pow(1e-38, -2.0/3); // = 2**2/3 GAMMA(2/3)
  sg->a2 = SynRad(sg->xlow) / exp(-sg->xlow);
  sg->c1 = pow(xmin, 1.0/3);

  // calculate the integrals of the approximate expressions
  if (xmin < 1) {               // low and high approx needed
    num_t sum1ap = 3*sg->a1*(1-pow(xmin,1.0/3));       // integral xmin --> 1
    num_t sum2ap = sg->a2*exp(-1);                     // integral 1 --> infin
    sg->ratio = sum1ap / (sum1ap+sum2ap);
  } else                        // only high approx needed
    sg->ratio = 0;              // generate only high energies using approx. 2
}

static num_t
syngen (num_t xmin)
{
  ensure(xmin >= 0, "invalid negative input value %.5g (>= expected)", xmin);

  static __thread struct syngen sg = {.last_xmin = -1};
  num_t approx, exact, result;

  if (sg.last_xmin != xmin) syngen_init(&sg, xmin);

  do {
    if(rnguni() < sg.ratio) { // use low energy approximation
      result = sg.c1 + (1-sg.c1)*rnguni();
      exact  = SynRad(CUB(result));
      approx = sg.a1/SQR(result);
    } else {                 // use high energy approximation
      result = sg.xlow - log(rnguni());
      exact  = SynRad(result);
      approx = sg.a2*exp(-result);
    }

  } while (exact < approx*rnguni());  // reject in proportion of approx

  return result;  // result now exact spectrum with unity weight
}

// --- public routines --------------------------------------------------------o

// e*e/6/pi/epsilon0/m/((Gev/GeV)**4) GeV
#define ENERGY_LOSS 9.5997636523126797e-19

num_t
mad_rad_randexp_seed (num_t val)
{
  ensure(val > 0, "invalid argument #1 (positive number expected)");
  num_t ret = seed;
  rngseed(val);
  return ret;
}

num_t
mad_rad_randexp (num_t mu)
{
  ensure(mu > 0, "invalid argument #1 (positive number expected)");
  return rngexp(mu);
}

num_t
mad_rad_synrad_prob (num_t betgam, num_t kick)
{
  ensure(betgam >= 0, "invalid betgam %.5g (>=0 expected)", betgam);
  ensure(kick   >= 0, "invalid kick %.5g (>=0 expected)"  , kick  );
  const num_t c1 = 2.5/M_SQRT3*P_ALPHAEM;
  return c1*kick*betgam;
}

num_t
mad_rad_freepath (num_t betgam, num_t kick, num_t length)
{
  ensure(betgam >= 0, "invalid betgam %.5g (>=0 expected)", betgam);
  ensure(kick   >= 0, "invalid kick %.5g (>=0 expected)"  , kick  );
  ensure(length >  0, "invalid length %.5g (>0 expected)" , length);
  if (kick < DBL_EPSILON) return INFINITY;
  num_t mu = length / mad_rad_synrad_prob(betgam, kick);
  return rngexp(mu);
}

num_t
mad_rad_nrjloss_quantum (num_t gamma, num_t kick, num_t length)
{
  ensure(gamma  >= 1, "invalid gamma %.5g (>=1 expected)", gamma );
  ensure(kick   >= 0, "invalid kick %.5g (>=0 expected)" , kick  );
  ensure(length >  0, "invalid length %.5g (>0 expected)", length);
  if (kick < DBL_EPSILON) return 0;
  const num_t c1 = 1.5*P_HBAR*P_CLIGHT;
  num_t nrj_crit = c1*CUB(gamma)*kick / length;
  num_t nrj_loss = nrj_crit*syngen(0);
  return nrj_loss;
}

num_t
mad_rad_nrjloss_average (num_t gamma, num_t kick, num_t length)
{
  ensure(gamma  >= 1, "invalid gamma %.5g (>=1 expected)", gamma );
  ensure(kick   >= 0, "invalid kick %.5g (>=0 expected)" , kick  );
  ensure(length >  0, "invalid length %.5g (>0 expected)", length);
  if (kick < DBL_EPSILON) return 0;
  num_t EEangle = SQR(gamma) * kick;
  num_t nrj_loss = ENERGY_LOSS * SQR(EEangle) / length;
  return nrj_loss;
}

#endif
