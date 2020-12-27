#ifndef _OUT_H
#define _OUT_H
#include "include_all.h"

int OutTau(LOS &los,Ion_all &spece_ionall,Setting spece_set,vector<double> &redshift_track_EW)
{
    int pointuseful=0;



    char filenametau[200];
    sprintf(filenametau, "./result/%s/galaxy%05d/tau_%f_%f_%f_%f.txt", spece_set.TIME,spece_set.galaxyindex_spec,spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
    ofstream out(filenametau);
    sprintf(filenametau, "./result/%s/galaxy%05d/flux_%f_%f_%f_%f.txt", spece_set.TIME,spece_set.galaxyindex_spec ,spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
    ofstream out2(filenametau);
    out << "redshift "
        << "all_rho "
        << "all_temp "
        << "all_metal ";
    for (int i = 0; i < spece_ionall.nions; i++)
    {
        out << spece_ionall.ions[i].name << "_rho "
            << spece_ionall.ions[i].name << "_temp "
            << spece_ionall.ions[i].name << "_metal "
            << spece_ionall.ions[i].name << "_tau ";
    }
    out << endl;
    out2 << "redshift "
         << "all_rho "
         << "all_temp "
         << "all_metal ";
    for (int i = 0; i < spece_ionall.nions; i++)
    {
        out2 << spece_ionall.ions[i].name << "_rho "
             << spece_ionall.ions[i].name << "_temp "
             << spece_ionall.ions[i].name << "_metal "
             << spece_ionall.ions[i].name << "_flux ";
    }
    out2 << endl;
    cerr << los.nzbins << " " << los.nvbins << endl;
    double redshift_track = spece_ionall.ion_total.redshift[0];
    cerr<<"sdjkfaaskjdf "<<redshift_track<<endl;
    double redshift_begin = spece_set.redshift_begin;
    double redshift_end = spece_set.redshift_end;
    int i = 0;
    while (redshift_track > redshift_begin - 0.01)
    {
        redshift_track -= VRES;
        if(redshift_track==0) pointuseful=i;
        if (redshift_track > redshift_begin && redshift_track < redshift_end)
        {
            redshift_track_EW.push_back(redshift_track);
            
            double h = spece_set.spece_para.h;
            double aex = spece_cosmo.aex;
            double rhomean = XH * 1.88e-29 * spece_set.spece_para.Omega_b * h * h / (aex * aex * aex);
            if (rhomean <= 0)
                cerr << rhomean << " " << spece_set.spece_para.Omega_b << " " << h << " " << aex << endl;
            for (int k = -1; k < spece_ionall.nions; k++)
            {
                if (k == -1)
                {
                    spece_ionall.ion_total.rhobins[i] = log10(spece_ionall.ion_total.rhobins[i] / rhomean);
                    spece_ionall.ion_total.tbins[i] = log10(spece_ionall.ion_total.tbins[i]);
                    if (isnan(spece_ionall.ion_total.rhobins[i]) || isinf(spece_ionall.ion_total.rhobins[i]))
                        spece_ionall.ion_total.rhobins[i] = 0.0;
                    if (isnan(spece_ionall.ion_total.tbins[i]) || isinf(spece_ionall.ion_total.tbins[i]))
                        spece_ionall.ion_total.tbins[i] = 0.0;
                }
                else
                {
                    spece_ionall.ions[k].rhobins[i] = log10(spece_ionall.ions[k].rhobins[i] / rhomean);
                    spece_ionall.ions[k].tbins[i] = log10(spece_ionall.ions[k].tbins[i]);
                    if (isnan(spece_ionall.ions[k].rhobins[i]) || isinf(spece_ionall.ions[k].rhobins[i]))
                        spece_ionall.ions[k].rhobins[i] = 0.0;
                    if (isnan(spece_ionall.ions[k].tbins[i]) || isinf(spece_ionall.ions[k].tbins[i]))
                        spece_ionall.ions[k].tbins[i] = 0.0;
                    if (isnan(spece_ionall.ions[k].Zbins[i]) || isinf(spece_ionall.ions[k].Zbins[i]))
                        spece_ionall.ions[k].Zbins[i] = 0.0;
                    if (isnan(spece_ionall.ions[k].vbins[i]) || isinf(spece_ionall.ions[k].vbins[i]))
                        spece_ionall.ions[k].vbins[i] = 0.0;
                }
                if (k == -1)
                {
                    out
                        << redshift_track << " "
                        << spece_ionall.ion_total.rhobins[i] << " "
                        << spece_ionall.ion_total.tbins[i] << " "
                        << spece_ionall.ion_total.Zbins[i] << " ";
                    out2
                        << redshift_track << " "
                        << spece_ionall.ion_total.rhobins[i] << " "
                        << spece_ionall.ion_total.tbins[i] << " "
                        << spece_ionall.ion_total.Zbins[i] << " ";
                }
                else
                {
                    out
                        << spece_ionall.ions[k].rhobins[i] << " "
                        << spece_ionall.ions[k].tbins[i] << " "
                        << spece_ionall.ions[k].Zbins[i] << " "
                        << spece_ionall.ions[k].vbins[i] << " ";
                    out2
                        << spece_ionall.ions[k].rhobins[i] << " "
                        << spece_ionall.ions[k].tbins[i] << " "
                        << spece_ionall.ions[k].Zbins[i] << " "
                        << pow(2.71828, (-1.0 * spece_ionall.ions[k].vbins[i])) << " ";
                }
            }
            out
                << spece_ionall.xbins[i] << " "
                << spece_ionall.ybins[i] << " "
                << spece_ionall.zbins[i]

                << endl;
            out2
                << spece_ionall.xbins[i] << " "
                << spece_ionall.ybins[i] << " "
                << spece_ionall.zbins[i]

                << endl;
        }
        i++;
    }
    time_t now_time = time(NULL);
    tm *t_tm = localtime(&now_time);
    out << "local time    :" << asctime(t_tm) << "        " << endl;
    out2 << "local time    :" << asctime(t_tm) << "        " << endl;
    out.close();
    out2.close();
    cerr << "over" << endl;
    return pointuseful;
}

#endif