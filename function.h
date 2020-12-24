#ifndef _FUNCTION_H
#define _FUNCTION_H
#include "include_all.h"
#include <algorithm>

#define Nmax 64
using namespace std;

class pixel
{
public:
    vector<double> value;
};
class galaxy_map
{
public:
    double dx[Nmax];
    double dy[Nmax];
    char filedestination[200];
    pixel colum_density[9 * Nmax * Nmax];

    galaxy_map()
    {
        for (int i = 0; i < Nmax; i++)
        {
            for (int k = 0; k < Nmax; k++)
            {
                dx[i] = 0.001 * (i - 315);
                dy[k] = 0.001 * (k - 315); //temppoint=ionid*Nmax*Nmax+i*Nmax+j;
            }
        }
    }
    void load(int &ionid, int &count,double &value);
    void pixel_sort();
    void output();
};
void galaxy_map::load(int &ionid,int &count,double &value)
{
    int temppoint=ionid*Nmax*Nmax+count%(Nmax*Nmax);
    colum_density[temppoint].value.push_back(value);

}
#ifdef para

void galaxy_map::pixel_sort()
{
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel
    {
        int id = 0;
        //id = omp_get_thread_num();
        for (int i = 0; i < 9 * Nmax * Nmax; i += NUM_THREADS)
        {
            sort(colum_density[i].value.begin(), colum_density[i].value.end());
        }
    }
}

#endif
void galaxy_map::output()
{
    char file[300];
    for (int ionid = 0; ionid < 9; ionid++)
    {
        sprintf(file,"%s_%d,txt",filedestination,ionid);
        ofstream tf(file);

        for (int i = 0; i < Nmax; i++)
        {
            for (int k = 0; k < Nmax; k++)
            {

                int temppoint = ionid * Nmax * Nmax + i * Nmax+ k;
                int medianpoint = colum_density[temppoint].value.size() / 2;
                tf << dx[i] << " "
                   << dy[k] << " "
                   << colum_density[temppoint].value[medianpoint] << endl;
            }
        }
        tf.close();
    }
}

void cleanworkplace(LOS &los,Ion_all spece_ionall)
{
    //spece_set.clear();
    cerr << "spece_set has been cleared" << endl;
    los.clear();
    cerr << "los has been cleared" << endl;
    spece_ionall.Freeions();
    cerr << "spece_ionall has been cleared" << endl;
}

void function_test(vector<Gas_1> &gp)
{

    double max = 0.0;
    double min = 0.0;
    char ffilenametest[200];
    for (int j = 0; j < 1; j++)
    {

        FILE *tf;
        sprintf(ffilenametest, "/data6/Hao_L/specexbin_huang/test%d.txt", j);
        if ((tf = fopen(ffilenametest, "wb")) == NULL)
        {
            cerr << "wrong" << endl;
            //return -1;
        }
        int n = gp.size();
        cerr << n << endl;
        W_binary(tf, n, 1);
        for (int i = 0; i < gp.size(); i++)
        {
            if (i % 100000 == 0)
                cerr << i << endl;
            W_binary(tf, gp[i].id, 1);

            W_binary(tf, gp[i].pos[0], 1);
            W_binary(tf, gp[i].pos[1], 1);
            W_binary(tf, gp[i].pos[2], 1);
            W_binary(tf, gp[i].pos[3], 1);
            //cerr<<gp[i].pos[3]<<endl;
            W_binary(tf, gp[i].vel[0], 1);
            W_binary(tf, gp[i].vel[2], 1);
            W_binary(tf, gp[i].vel[2], 1);

            W_binary(tf, gp[i].mass, 1);
            W_binary(tf, gp[i].density, 1);
            W_binary(tf, gp[i].T, 1);

            W_binary(tf, gp[i].Z[0], 1);
            W_binary(tf, gp[i].Z[1], 1);
            W_binary(tf, gp[i].Z[2], 1);
            W_binary(tf, gp[i].Z[3], 1);

            W_binary(tf, gp[i].Hsml, 1);
            W_binary(tf, gp[i].sfr, 1);
        }

        fclose(tf);
    }
    exit(-1);
}

#endif