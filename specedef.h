#ifndef SPECEDEF_H
#define SPECEDEF_H
#include "math.h"

#define Pi 3.14159265358979323846;
#define km 1.e5;
#define Mpc 3.086E24;
#define m_p 1.6726231E-24; /* proton mass */
#define k_B 1.380622E-16;  /* Boltzman constant */

#define NDIM 3
#define BOXSIZE 1.0
#define HALFBOX 0.5
#define HBEXTRA 0.55
#define NSRCHRAD 2.0

#define NIONS 29
#define SOLAR_METALS 0.0189

/* redshift and velocity resolution */
//#define ZRES            6.0e-07
//#define VRES            3.0e-06
#define ZRES 3.0e-06
#define VRES 1.5e-05
//#define ZRES            3.0e-07

/* Physical constants */
#define XH 0.76
#define KBOLTZ 1.381e-16
#define MHYDR 1.673e-24
#define CLIGHT 2.99792458e10
#define PI 3.14159265

#define NMETALS 4

//Parameter spece_para;
//Setting spece_set;

#define NINTERP 10000

#define NHILIM 1.58489e17
//#define NHILIM 1.e18	// agrees better w/Faucher-Giguere,Keres 2010 Fig 3
#define P0BLITZ 3.5e+04	 // Blitz & Rosolowski 2006
#define ALPHA0BLITZ 0.92 // Blitz & Rosolowski 2006

//model set
#define VELOCITY_UNIT_CORRECTION

int binarysearch(double key, double *array, int nbins)
{
	int low = 0;
	int high = nbins - 1;
	int middle;

	if (key >= array[high])
		return high;
	if (key < array[low])
		return low;

	while (low < high)
	{
		middle = floor(low + high) / 2;
		if (key >= array[middle] && key < array[middle + 1])
			return middle;
		else if (key < array[middle])
			high = middle; //search low end of array
		else
			low = middle + 1; //search high end of array
	}
	return -1; //search key not found
}

int harxsearch(double key, double *array, int nbins)
{
	int low = 0;
	int high = nbins - 1;
	int middle;
	double dbin = array[1] - array[0];

	if (key >= array[high])
		return high;
	if (key < array[low])
		return low;

	middle = (int)((key - array[0]) / dbin);
	if (middle >= low && middle <= high)
		return middle;
	else
		return -1;
}

#endif