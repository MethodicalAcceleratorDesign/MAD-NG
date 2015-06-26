// compilation
// gcc -std=c11 -W -Wall -Wextra -pedantic -static-libgcc -O3 -ffast-math complex.c -o complex
// usage: ./complex 100000001

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

int
main(int argc, char *argv[])
{
	double complex a = (1+I)/(sqrt(2)+1e-8), b = 1;
	int n = argc > 1 ? strtol(argv[1],0,0) : 100000000;

	for(int i=1; i <= n; i++)
		b *= a ;

	printf("%d: a=%.14g%+.14gi, b=%.14g%+.14gi\n", n, creal(a), cimag(a), creal(b), cimag(b));
}
