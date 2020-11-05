#ifndef KERNINTTABLE_H
#define KERNINTTABLE_H
#include "include_all.h"
using namespace std;

float KernIntTable[NINTERP + 1][2];

int InitKernIntTable()
{
	int i;
	float xw, dx, sum, kint, kern;
	float Kernel();

	sum = kint = 0;
	dx = NSRCHRAD / NINTERP;
	for (i = NINTERP - 1; i >= 0; i--)
	{
		xw = i * dx;
		if (xw <= 0.5)
		{
			kern = 1 - 6 * xw * xw + 6 * xw * xw * xw;
		}
		else
		{
			kern = 2 * (1 - xw) * (1 - xw) * (1 - xw);
		}
		kern *= 8. / Pi;
		sum += kern * 4  * xw * xw * dx *Pi;
		kint += kern * dx;
		KernIntTable[i][0] = sum;
		KernIntTable[i][1] = kint;
		if (i % 1000 == 0)
			fprintf(stdout, "KERNTABLE BUILD: %d %g %g %g\n", i, xw, sum, kint);
	}
	return 0;
}


#endif