#ifndef PARAMETER_H
#define PARAMETER_H
#include "include_all.h"
using namespace std;

class Parameter
{
public:
    double totMass;
    double Lambda;
    double Omega_b;
    double H_0;
    double h;

    double unit_Time, unit_Density, unit_Length, unit_Mass, unit_Velocity, unit_Temperature, unit_DUDT;
    double H0, t0, aex3, aexhub;

    void load();
};

void Parameter::load()
{
    float L;

    spece_set.boxsize /= h;
    L = spece_set.boxsize;
    if (Lambda > 0.01 && totMass < 1.0)
    {
        if (Lambda != 1. - totMass)
        {
            /*			fprintf(stderr,"Setting Lambda = %g\n",1.-totMass);*/
            Lambda = 1. - totMass;
        }
    }
    else
        Lambda = 0.0;

    H0 = (8.0 / 3.0) * Pi;
    H0 = sqrt(H0);
    t0 = 2. / (3 * H0);

    unit_Time = H0 * Mpc;
    unit_Time /= 100.0 * h * km;
    unit_Density = 6.8451e-29 * h * h;
    unit_Length = L * Mpc;//??is the unit is calculate by special melogy????
    unit_Mass = unit_Density * unit_Length * unit_Length * unit_Length;
    unit_Velocity = unit_Length / unit_Time;
    unit_Temperature = unit_Velocity * unit_Velocity * m_p;
    unit_Temperature /=k_B;
    unit_DUDT = unit_Density * unit_Velocity * unit_Velocity / unit_Time;
}

#endif