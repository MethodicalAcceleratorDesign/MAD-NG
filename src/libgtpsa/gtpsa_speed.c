/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA Example Speed
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

  WARNING: this file is not part of MAD, it must be used as an example with the
           standalone library, see mad/src/REAME.GTPSA for more info.
*/

#include <time.h>
#include <complex.h>
#include "mad_ctpsa.h"

int main(void)
{ 
   // descriptor for TPSA with 4 variables of order 4 (max)
   const desc_t *d = mad_desc_newv(4, 4);

   // create three TPSAs of order 0 (as scalars) and set values
   ctpsa_t *t1 = mad_ctpsa_newd(d,4); mad_ctpsa_setvar(t1,csqrt(-1.0),1,0);
   ctpsa_t *t2 = mad_ctpsa_newd(d,4); mad_ctpsa_setvar(t2,csqrt(-1.5),2,0);
   ctpsa_t *t3 = mad_ctpsa_newd(d,4); mad_ctpsa_setvar(t3,csqrt(-2.0),3,0);
   ctpsa_t *t4 = mad_ctpsa_newd(d,4); mad_ctpsa_setvar(t4,csqrt(-2.5),4,0);
   ctpsa_t *t5 = mad_ctpsa_newd(d,4);

   mad_ctpsa_sqrt(t1,t1);                     // mad_ctpsa_print(t1,"T1",0,0,stdout);
   mad_ctpsa_sqrt(t2,t2);                     // mad_ctpsa_print(t2,"T2",0,0,stdout);
   mad_ctpsa_sqrt(t3,t3);                     // mad_ctpsa_print(t3,"T3",0,0,stdout);
   mad_ctpsa_sqrt(t4,t4);                     // mad_ctpsa_print(t4,"T4",0,0,stdout);

   mad_ctpsa_axypbvwpc(1,t1,t2,1,t3,t4,0,t5); // mad_ctpsa_print(t5,"T5",0,0,stdout);
   mad_ctpsa_sqrt(t5,t4);                     // mad_ctpsa_print(t4,"T4",0,0,stdout);

   // multiply the scalars
   clock_t ft0 = clock();
   FOR(i,5000000) mad_ctpsa_sqrt(t4,t5);
   clock_t ft1 = clock();
   printf("CPU time: %.6f s\n", (num_t)(ft1-ft0)/CLOCKS_PER_SEC);

   // destroy the three TPSAs
   mad_ctpsa_del(t1);
   mad_ctpsa_del(t2);
   mad_ctpsa_del(t3);
   mad_ctpsa_del(t4);
   mad_ctpsa_del(t5);

   // destroy all created descriptors (optional cleanup)
   mad_desc_del(0);
   return 0;
}
