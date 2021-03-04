#include <iostream>
#include "include_list.h"


#define RRR
#define Nions 9

using namespace std;


class columnDatas
{
public:
    int t_lines;
    int galaxy_index;

    double r;
    double r2rvir;
    double r2rspsh;
    double ys[Nions];
    columnDatas()
    {
        r = -1.0;
        r2rvir = 0.0;
        r2rspsh = 0.0;
        for (int i = 0; i < Nions; i++)
        {
            ys[i] = 0.0;
        }
    }
};

int main(int argc, char *argv[])
{
    double h = 0.72;

    char result_path[100];
    sprintf(result_path, "/data6/Hao_L/my_specexbin/result");
    float redshift;
    R_assii(argv[2],redshift);
    double deltarin = 10; //kpc
    double deltarout = 100;
#ifdef RRR
    double rin = 200;
    int Nin = (int)rin / deltarin;
    int Nout = (int)(1000.0 - rin) / deltarout;
    int N = Nin + Nout + 1;
#endif
#ifdef R2Rv

    double deltalogin = 0.02;
    double deltalogout = 0.3;
    int N = 2.0 / deltalogin + 3.0 / deltalogout + 1;
#endif
#ifdef R2Rs
    double deltalogin = 0.02;
    double deltalogout = 0.3;
    int N = 2.0 / deltalogin + 3.0 / deltalogout + 1;
#endif
    int *averagenumber = new int[N + 1];
    columnDatas *result = new columnDatas[N + 1];
    for (int i = 0; i < N + 1; i++)
    {
        averagenumber[i] = 0;
    }
    vector<columnDatas> pp;

    char path[1000];
    sprintf(path, "%s/0_a_width/%s/lines_rp_rper_ews.txt", result_path, argv[1]);
    ifstream in(path);

    char path2[100];
    sprintf(path2, "%s/0_a_width/%s/average_EW.txt", result_path, argv[1]);
    ofstream out(path2);
    cerr << 1 << endl;
    while (!in.eof())
    {
        columnDatas p;
        char buffer[80];
        double values;
        int value;

        in >> buffer;
        R_assii(buffer, value);
        p.t_lines = value; //1
        R_assii(buffer, value);
        p.galaxy_index = value; //2
        R_assii(buffer, values);
        p.r = value; //3
        p.r = 1000.0*p.r / 0.72 / (1 + redshift);
        R_assii(buffer, value);
        p.r2rvir = value; //4

        for (int i = 0; i < Nions; i++)
        {
            in >> buffer; //5-13
            R_assii(buffer, values);
            p.ys[i] = values;
        }

        pp.push_back(p);
    }
    pp.pop_back();
    cerr << pp.size() << endl;
    in.close();

    for (int i = 0; i < pp.size(); i++)
    {
#ifdef RRR
        if (pp[i].r > 0 && pp[i].r < rin)
        {
            int n = (int)(pp[i].r / (deltarin));
            averagenumber[n]++;
            for(int ion=0;ion<Nions;ion++)
            {
                result[n].ys[ion] += pp[i].ys[ion];
            }
            
        }
        else if (pp[i].r > rin && pp[i].r < 1000.0)
        {
            int n = Nin + (int)((pp[i].r - rin) / (deltarout));
            averagenumber[n]++;
            for(int ion=0;ion<Nions;ion++)
            {
                result[n].ys[ion] += pp[i].ys[ion];
            }
            
        }
#endif
#ifdef R2Rv
        if (pp[i].r2rvir > 0 && pp[i].r2rvir < 2.0)
        {
            int n = (int)(pp[i].r2rvir / deltalogin);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }
        else if (pp[i].r2rvir > 2.0 && pp[i].r2rvir < 5.0)
        {
            int n = 2.0 / deltalogin + (int)((pp[i].r2rvir - 2.0) / deltalogout);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }
#endif
#ifdef R2Rs
        if (pp[i].r2rspsh > 0 && pp[i].r2rspsh < 2.0)
        {
            int n = (int)(pp[i].r2rspsh / deltalogin);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }
        else if (pp[i].r2rspsh > 2.0 && pp[i].r2rspsh < 5.0)
        {
            int n = 2.0 / deltalogin + (int)((pp[i].r2rspsh - 2.0) / deltalogout);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }

#endif
        else
            cerr << "wrong!!!!!!!!" << pp[i].r << endl;
    }

    for (int i = 0; i < N + 1; i++)
    {
#ifdef RRR
        if (i < Nin)
            result[i].r = deltarin / 2.0 + (double)(i)*deltarin;
        else if (i >= Nin)
            result[i].r = rin + deltarout / 2.0 + (double)(i - Nin) * deltarout;
#else
        if (i < 2.0 / deltalogin)
            result[i].r = deltalogin / 2.0 + (double)(i)*deltalogin;
        else if (i >= 2.0 / deltalogin)
            result[i].r = 2.0 + deltalogout / 2.0 + (double)(i - 2.0 / deltalogin) * deltalogout;
#endif
        if (averagenumber[1] != 0)
        {
            out << result[i].r << " ";
            for(int ion=0;ion<Nions;ion++)
            {
                result[i].ys[ion] /= averagenumber[i];
                out<<result[i].ys[ion]<<" ";
            }
             out<< endl;
        }
        else if (averagenumber[i] == 0)
        {
            cerr << " 00000:   " << result[i].ys[1] << " " << result[i].ys[2] << " " << result[i].ys[3] << endl;
        }
    }
    out.close();
    delete[] averagenumber;
    delete[] result;
}
