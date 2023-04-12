/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA Example 8
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

#include "mad_tpsa.h"

int main(void)
{ 
   // descriptor for TPSA with 6 variables of order 2 (max)
   const desc_t *d  = mad_desc_newv(6, 2);

   // const desc_t *d0 = mad_desc_newv(6, 0); (void)d0;
   // -> error: mad_desc.c:1230: : invalid maximum order: 0 (0<?<=63)

   // create three TPSAs of order 0 (as scalars) and set values
   tpsa_t *t1 = mad_tpsa_newd(d,0); mad_tpsa_set0(t1,0,1.5);
   tpsa_t *t2 = mad_tpsa_newd(d,0); mad_tpsa_set0(t2,0,2.0);
   tpsa_t *t3 = mad_tpsa_newd(d,0);

   mad_tpsa_print(t1,"T1",0,0,stdout);
   mad_tpsa_print(t2,"T2",0,0,stdout);

   // add the scalars
   mad_tpsa_add(t1,t2,t3);
   mad_tpsa_print(t3,"T1+T2",0,0,stdout);

   // multiply the scalars
   mad_tpsa_mul(t1,t2,t3);
   mad_tpsa_print(t3,"T1*T2",0,0,stdout);

   // log of the scalar
   mad_tpsa_log(t3,t3);
   mad_tpsa_print(t3,"log(T1*T2)",0,0,stdout);

   // destroy the three TPSAs
   mad_tpsa_del(t1);
   mad_tpsa_del(t2);
   mad_tpsa_del(t3);

   // destroy all created descriptors (optional cleanup)
   mad_desc_cleanup();
   return 0;
} 
