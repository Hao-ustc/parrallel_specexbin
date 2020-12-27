#ifndef _DATA_RESOLVE_H
#define _DATA_RESOLVE_H
#include "include_all.h"

using namespace std;

class my_EW
{
public:
    vector<double> redshift_track;
    double  *ew;
    int Nions;
    void data_resolve( Ion_all &spece_ionall, Setting spece_set);
};

void my_EW::data_resolve( Ion_all &spece_ionall, Setting spece_set)
{
    Nions = spece_ionall.nions;
    ew = new double[Nions];
    for (int i = 0; i < Nions; i++)
    {
        ew[i] = 0.0;
    }
    for (int j = 0; j < Nions; j++)
    {
        for (int i = 1; i < redshift_track.size(); i++)
        {
            double dlambta = fabs(spece_ionall.ions[j].lambda * (redshift_track[i] - redshift_track[i - 1]));
            double flux = pow(2.71828, (-1.0 * spece_ionall.ions[j].vbins[i]));
            flux += pow(2.71828, (-1.0 * spece_ionall.ions[j].vbins[i - 1]));
            flux /= 2.0;
            double dwidth = dlambta * (1 - flux);
            ew[j] += dwidth;
        } //now width is the EW of the line
    }
}

#endif