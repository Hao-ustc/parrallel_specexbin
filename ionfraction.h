#ifndef IONFRACTION_H
#define IONFRACTION_H
//#include "/data6/Hao_L/code/head/coma/include.h"
#include "include_list.h"
#include "DEAL.h"
#include "include_all.h"
#define N 25
#define XH 0.76
#define MH_YDR 1.673e-24

#define NHPTS 240
#define TPTS 140
#define NHLOW -9.0
#define TLOW 2.5
#define DELTANH 0.05
#define DELTAT 0.05

using namespace std;

class ionfraction
{
private:
    double iontable[240 * 140];
    float nhl[NHPTS], tl[TPTS];
    float y1, y2, y3, y4, t, u;
    float n_h;
    float n_hmin, n_hmax, tmin, tmax; //all after log10();
    double file_h;
    double unit_density;
    double aex;

    float realrho;

public:
    double fraction;
    //   float rho;
    //  float T;

    void loadtable(int ionid, Setting spece_set);
    void loadparam(datahead &b, Setting spece_set);

    int pointtrace(int i, int j);
    double interplate(float &rho, float &T);
    void test();
    ionfraction()
    {
        n_hmin = -9.0;
        n_hmax = 3.0;
        tmin = 2.5;
        tmax = 9.5;
        //  rho = 0.0;
        // T = 0.0;
    }
};

void ionfraction::loadtable(int ionid, Setting spece_set)
{
    double ionredshift = spece_set.redshift_center;
    FILE *table1, *table2;
    char prefix[200];
    sprintf(prefix, "./ionfiles/");
    int z10 = (int)(ionredshift * 10.);
    int zlo = z10;
    char filename[200];
    sprintf(filename, "%slt%.2dHM12_i9", prefix, z10);
    while ((table1 = fopen(filename, "r")) == NULL && zlo > 0)
    {
        zlo--;
        sprintf(filename, "%slt%.2dHM12_i9", prefix, zlo);
    }
    fclose(table1);
    ifstream ionin1(filename);
    int zhi = z10 + 1;
    sprintf(filename, "%slt%.2dHM12_i9", prefix, zhi);
    while ((table2 = fopen(filename, "r")) == NULL && zhi < 80)
    {
        zhi++;
        sprintf(filename, "%slt%.2dHM12_i9", prefix, zhi);
    }
    fclose(table2);
    ifstream ionin2(filename);
    if (zlo < 0 || zhi >= 80)
        exit(-1);
    cerr << "LOAD ION" << ionid << endl;
    float redzlo = 0.1 * zlo;
    float redzhi = 0.1 * zhi;
    int count = 0;
    int i = 0;
    double fraction1;
    double fraction2;
    while (!ionin1.eof())
    {
        char buffer1[100];
        ionin1 >> buffer1;
        char buffer2[100];
        ionin2 >> buffer2;
        count++;
        if ((count - ionid) % 9 == 0)
        {
            R_assii(buffer1, fraction1);
            R_assii(buffer1, fraction2);

            iontable[i] = ((redzhi - ionredshift) * fraction1 + (ionredshift - redzlo) * fraction2) / (redzhi - redzlo);
            i++;
        }
    }
    ionin1.close();
    ionin2.close();
}
void ionfraction::loadparam(datahead &b, Setting spece_set)
{
    file_h = b.HubbleParam;
    unit_density = spece_set.spece_para.unit_Density;
    aex = 1. / (b.Redshift + 1.0);

    for (int inh = 0; inh < NHPTS; inh++)
        nhl[inh] = NHLOW + inh * DELTANH;
    for (int itemp = 0; itemp < TPTS; itemp++)
        tl[itemp] = TLOW + itemp * DELTAT;
}

int ionfraction::pointtrace(int i, int j)
{
    int temppoint = i * TPTS + j;
    if (temppoint >= NHPTS * TPTS)
        exit(-1);
    else
        return temppoint;
}
double ionfraction::interplate(float &rho, float &T)
{

    realrho = rho * XH * unit_density / (aex * aex * aex);
    n_h = realrho / MH_YDR;
    int inh = (int)((log10(n_h) - NHLOW) / DELTANH);
    int itemp = (int)((log10(T) - TLOW) / DELTAT);

    if (inh > NHPTS - 2)
        inh = NHPTS - 2;
    if (itemp > TPTS - 2)
        itemp = TPTS - 2;
    if (inh < 0)
        inh = 0;
    if (itemp < 0)
        itemp = 0;
    t = (log10(n_h) - nhl[inh]) / DELTANH;
    u = (log10(T) - tl[itemp]) / DELTAT;

    y1 = iontable[pointtrace(inh, itemp)];
    y2 = iontable[pointtrace((inh + 1), itemp)];
    y3 = iontable[pointtrace((inh + 1), (itemp + 1))];
    y4 = iontable[pointtrace(inh, (itemp + 1))];
    //cerr <<y1<<" "<<y2<<" "<<y3<<" "<<y4<<" "<<t<<" " <<itemp<< endl;
    fraction = pow(10., (1. - t) * (1. - u) * y1 + t * (1. - u) * y2 + t * u * y3 + (1. - t) * u * y4);
    if (fraction < 1e-10)
        cerr << "T : " << T << " rho : " << rho << " realrho : " << realrho << endl;
    if (fraction > 0.01)
        // cerr << "good fraction ";
        return realrho;
}
void ionfraction::test()
{
    if (isnan(fraction) && isinf(t))
    {

        cerr << " unitdensity:  " << unit_density;
        cerr << " aex:   " << aex;
        cerr << " rho : " << realrho << endl;
    }
}
#endif
