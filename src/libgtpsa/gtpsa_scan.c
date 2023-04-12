/*
 o-----------------------------------------------------------------------------o
 |
 | TPSA regression tests on descriptors
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

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "mad_log.h"
#include "mad_tpsa.h"

// gtpsa_scan [..|--] nv mo np po 'm:' o1 o2 o3
int main(int argc, const char *argv[])
{
  // regression test that scan all possible descriptors for nv,mo.
  int  argi = 1;
  int   nvl = argc > argi ? strcmp(argv[argi],"..") ? 0 : 1 : 0;   argi+=nvl;
  int   nol = argc > argi ? strcmp(argv[argi],"--") ? 1 : 0 : 1;   argi+=!nol;
  int    nv = argc > argi ? strtol(argv[argi],0,10) : 6;           argi+=1;
  ord_t  mo = argc > argi ? strtol(argv[argi],0,10) : 4*(2-nvl);   argi+=1;
  ord_t  np = argc > argi ? strtol(argv[argi],0,10) : 0;           argi+=1;
  ord_t  po = argc > argi ? strtol(argv[argi],0,10) : mo;          argi+=1;
  int   mol = argc > argi ? strcmp(argv[argi],"m:") ? 0 : 1 : 0;   argi+=mol;

  ssz_t n0 = mol ? argc-argi : 0;
  int nn = nv+np;
  ord_t m0[nn];

  if (!nvl)
    printf("**** nv=%d, mo=%d, np=%d, po=%d ****\n", nv, mo, np, po);

  if (!n0) {
    idx_t i=0;
    for (; i < nv; i++) m0[i] = mo;
    for (; i < nn; i++) m0[i] = MIN(mo,po);
  } else {
    ensure(n0 == nn, "inconsistency between monomial len %d and nn %d", n0, nn);
    for (idx_t i=0; i < nn; i++)
      m0[i] = strtol(argv[argi++],0,10);
    printf("**** m0: "); mad_mono_print(nn,m0); printf("\n");
  }

  for (idx_t nv=nvl?np+1:nn; nv <= nn; nv++) {
    printf("**** nv=%d, mo=%d, np=%d, po=%d ****\n", nv, mo, np, po);

    ord_t m[nv];
    if (nol) mad_mono_fill(nv, m,  1);
    else     mad_mono_copy(nv, m0, m);

    ord_t o = mad_mono_ord(nv, m);
    printf("** "); mad_mono_print(nv,m); printf(", o=%2d | ", o);
    mad_desc_del(mad_desc_newvpo(nv, 0, np, po, m)); printf("\n");

    for(idx_t k=1; nol; k++) {
      if (k == 1000) { fprintf(stderr, "."); k=0; }

      idx_t i;
      for (i=0; i < nv; ++i) {
        if (++m[i] <= m0[i]) break;
        m[i] = 1;
      }
      if (i == nv) break; // no more case

      ord_t o = mad_mono_ord(nv, m);
      printf("** "); mad_mono_print(nv,m); printf(", o=%2d | ", o);
      mad_desc_del(mad_desc_newvpo(nv, 0, np, po, m)); printf("\n");
    }
  }

  fprintf(stderr, "\n");
  return 0;
}
