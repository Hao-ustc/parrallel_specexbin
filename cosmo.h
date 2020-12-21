#ifndef COSMO_H
#define COSMO_H
#include "include_all.h"
static int startflag = 1;
static double etaold;


class Cosmo
{
public:
	double hubble;
	double aex;
	double redshift;
};
Cosmo spece_cosmo;

int cosmopar(float &t,Setting &spece_set) //here t may have a special unit please pay attention to it
{
	float t1;
	double tol = 1.e-6;
	double a0, astar, eta, etalast;
	double f, fprime;
	int it;

	double t0 = spece_set.spece_para.t0;
	double H0 = spece_set.spece_para.H0;
	double aex3 = spece_set.spece_para.aex3;
	double aexhub = spece_set.spece_para.aexhub;
	double hubble = spece_cosmo.hubble;
	double aex = spece_cosmo.aex;
	double redshift = spece_cosmo.redshift;
	double totMass = spece_set.spece_para.totMass;
	double Lambda = spece_set.spece_para.Lambda;
	if (startflag)
	{
		
		startflag = 0;
		etaold = 1.;
	}
    
	if (fabs(spece_set.spece_para.totMass - 1.0) < tol)
	{
		t1 = t / t0;
		aex3 = t1 * t1;
		aex = pow((aex3), 1. / 3.);
		hubble = 2.0 / 3.0 / t; //what is the unit of hubble is it meet us?
		aexhub = aex * hubble;
		redshift = 1. / aex - 1.;
	}
	else if (totMass < 1.0)
	{
		if (Lambda > 0.0)
		{ /* Flat, low-density universe */
			eta = sqrt(1. - totMass) * 1.5 * H0 * t;
			aex = pow(sqrt(totMass / (1. - totMass)) * sinh(eta), 2. / 3);
			aex3 = aex * aex * aex;
			hubble = H0 * sqrt(totMass / aex3 + Lambda);
			aexhub = aex * hubble;
			redshift = 1. / aex - 1.;
		}
		else
		{ /* Open universe */
			a0 = 1. / H0 / sqrt(1. - totMass);
			astar = .5 * H0 * H0 * totMass * a0 * a0 * a0;
			it = 0;
			eta = etaold;
			do
			{
				f = astar * (sinh(eta) - eta) - t;
				fprime = astar * (cosh(eta) - 1.);
				etalast = eta;
				eta = eta - f / fprime;
				if ((it++) > 20)
				{
					fprintf(stderr, "Overiterated in cosmopar %d %g %g\n", it, eta, etalast);
					break;
				}
			} while (fabs(eta - etalast) / etalast > tol);

			aex = astar * (cosh(eta) - 1.) / a0;
			aex3 = aex * aex * aex;
			etaold = eta;
			redshift = 1. / aex - 1.;
			hubble = H0 * (1. + redshift) * sqrt(1. + totMass * redshift);
			aexhub = aex * hubble;
		}
	}
	else if (totMass > 1.0)
	{
		a0 = 1. / H0 / sqrt(totMass - 1.);
		astar = .5 * H0 * H0 * totMass * a0 * a0 * a0;
		it = 0;
		eta = etaold;
		do
		{
			f = astar * (eta - sin(eta)) - t;
			fprime = astar * (1. - cos(eta));
			etalast = eta;
			eta = eta - f / fprime;
			if ((it++) > 20)
			{
				fprintf(stderr, "Overiterated in cosmopar %d %g %g\n", it, eta, etalast);
				break;
			}
		} while (fabs(eta - etalast) / etalast > tol);

		aex = astar * (1. - cos(eta)) / a0;
		aex3 = aex * aex * aex;
		etaold = eta;
		redshift = 1. / aex - 1.;
		hubble = H0 * (1. + redshift) * sqrt(1. + totMass * redshift);
		aexhub = aex * hubble;
	}
	spece_set.spece_para.t0 = t0;
	spece_set.spece_para.H0 = H0;
	spece_set.spece_para.aex3 = aex3;
	spece_set.spece_para.aexhub = aexhub;
	spece_cosmo.hubble = hubble;
	spece_cosmo.aex = aex;
	spece_cosmo.redshift = redshift;
	spece_set.spece_para.totMass = totMass;
	spece_set.spece_para.Lambda = Lambda;
	return 0;
}

float CosmicTime(float &z,Setting &spece_set) /* returns system time at redshift z */ //unit=40Gyr;
	
{
	double tol = 1.e-6;
	double a0, astar, eta, etalast, t = 1.;
	double f, fprime, aextemp;
	int it;

	double t0 = spece_set.spece_para.t0;
	double H0 = spece_set.spece_para.H0;
	double aex3 = spece_set.spece_para.aex3;
	double aexhub = spece_set.spece_para.aexhub;
	double hubble = spece_cosmo.hubble;
	double aex = spece_cosmo.aex;
	double redshift = spece_cosmo.redshift;
	double totMass = spece_set.spece_para.totMass;
	double Lambda = spece_set.spece_para.Lambda;
	if (startflag)
	{
		startflag = 0;
		etaold = 1.;
	}

	if (fabs(totMass - 1.0) < tol)
		t = t0 * pow(1. + z, - 1.5);
	else if (totMass < 1.0)
	{
		if (Lambda > 0.0)
		{
			it = 0;
			aextemp = 1. / (1. + z);
			t = etaold;
			do
			{
				f = sqrt(Lambda / totMass) * pow(aextemp, 1.5) +
					sqrt(aextemp * aextemp * aextemp * Lambda / totMass + 1) -
					exp(1.5 * sqrt(Lambda) * H0 * t);
				fprime = -1.5 * sqrt(Lambda) * H0 * exp(1.5 * sqrt(Lambda) * H0 * t);
				etalast = t;
				t = t - f / fprime;
				if ((it++) > 20)
					break;
			} while (fabs(t - etalast) / etalast > tol);
			etaold = t;
		}
		else
		{
			a0 = 1. / H0 / sqrt(1. - totMass);
			astar = .5 * H0 * H0 * totMass * a0 * a0 * a0;
			aextemp = 1. / (1. + z);
			it = 0;
			eta = etaold;
			do
			{
				f = astar * (cosh(eta) - 1.) / a0 - aextemp;
				fprime = sinh(eta) * astar / a0;
				etalast = eta;
				eta = eta - f / fprime;
				if ((it++) > 20)
				{
					fprintf(stderr, "Overiterated in CosmicTime %d %g %g\n", it, eta, etalast);
					break;
				}
			} while (fabs(eta - etalast) / etalast > tol);

			t = astar * (sinh(eta) - eta);
			etaold = eta;
		}
	}
	else if (totMass > 1.0)
	{
		a0 = 1. / H0 / sqrt(totMass - 1.);
		astar = .5 * H0 * H0 * totMass * a0 * a0 * a0;
		aextemp = 1. / (1. + z);
		it = 0;
		eta = etaold;
		do
		{
			f = astar * (1. - cos(eta)) / a0 - aextemp;
			fprime = sin(eta) * astar / a0;
			etalast = eta;
			eta = eta - f / fprime;
			if ((it++) > 20)
			{
				fprintf(stderr, "Overiterated in CosmicTime %d %g %g\n", it, eta, etalast);
				break;
			}
		} while (fabs(eta - etalast) / etalast > tol);

		t = astar * (eta - sin(eta));
		etaold = eta;
	}

	spece_set.spece_para.t0 = t0;
	spece_set.spece_para.H0 = H0;
	spece_set.spece_para.aex3 = aex3;
	spece_set.spece_para.aexhub = aexhub;
	spece_cosmo.hubble = hubble;
	spece_cosmo.aex = aex;
	spece_cosmo.redshift = redshift;
	spece_set.spece_para.totMass = totMass;
	spece_set.spece_para.Lambda = Lambda;
	
	return t;
}

#endif