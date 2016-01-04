#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "mad_tpsa.h"
#include "mad_track.h"

int
main(void)
{
	typedef unsigned char ord_t;

	typedef struct tpsa tpsa_t;
	typedef struct tpsa_desc desc_t;
	typedef struct { tpsa_t *x,*px,*y,*py,*s,*ps; } map_t;

	desc_t *d = mad_tpsa_desc_new(6, (ord_t[]){1,1,1,1,1,1}, 0, 0);
	map_t m = {
		.x  = mad_tpsa_newd(d, mad_tpsa_default),
		.px = mad_tpsa_newd(d, mad_tpsa_default),
		.y  = mad_tpsa_newd(d, mad_tpsa_default),
		.py = mad_tpsa_newd(d, mad_tpsa_default),
		.s  = mad_tpsa_newd(d, mad_tpsa_default),
		.ps = mad_tpsa_newd(d, mad_tpsa_default),
	};

	mad_tpsa_set0(m.x,  0,0);
	mad_tpsa_set0(m.y,  0,0);
	mad_tpsa_set0(m.s,  0,0);
	mad_tpsa_set0(m.px, 0,0.001);
	mad_tpsa_set0(m.py, 0,0.001);
	mad_tpsa_set0(m.ps, 0,1e-6);
	mad_tpsa_setm(m.ps, 3,(ord_t[]){1,0,0}, 0.0,1.0);
	mad_tpsa_setm(m.ps, 3,(ord_t[]){0,0,1}, 0.0,1.0);

	for (int i=0; i < 10000000; i++)
	 	mad_track_drift((tpsa_t**)&m, 1, 1, 0);

  mad_tpsa_print(m.x, stdout);
  mad_tpsa_print(m.y, stdout);
  mad_tpsa_print(m.s, stdout);
  mad_tpsa_print(m.px, stdout);
  mad_tpsa_print(m.py, stdout);
  mad_tpsa_print(m.ps, stdout);

  mad_tpsa_del(m.x );
  mad_tpsa_del(m.y );
  mad_tpsa_del(m.s );
  mad_tpsa_del(m.px);
  mad_tpsa_del(m.py);
  mad_tpsa_del(m.ps);
  mad_tpsa_desc_del(d);

	return 0;
}
