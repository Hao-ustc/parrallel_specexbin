#ifndef DEAL_H
#define DEAL_H
#include "include_list.h"
//mademychange
#define N 25
using namespace std;

class Data
{
public:
    double galaxyid_column;
    int particalid;
    double mass;
    double pos[4];
    double vel[4];
    double galaxymass;
    Data()
    {
        galaxyid_column = -1.0;
        particalid = -1;
        mass = 0.0;
        pos[3] = -1.0;
        vel[3] = -100.0;
        galaxymass = 0.0;
    };
    //all get from class Gas
};

void Readfile(ifstream &in, vector<Data> &d)
{
    Data dpoint;
    int count = 0;
    while (!in.eof())
    {
        char buffer[100];
        int tempdata1 = 0;
        double tempdata2 = 0.0;
        in >> buffer;
        count++;
        if (count > 9)
        {
            if ((count - 10) % 9 == 0)
            {
                R_assii(buffer, dpoint.galaxyid_column);
                //cerr<<count<<endl;
            }
            if ((count - 11) % 9 == 0)
            {
                R_assii(buffer, dpoint.particalid);
                //cerr<<count<<endl;
            }

            else if ((count - 12) % 9 == 0)
            {
                R_assii(buffer, dpoint.mass);
            }
            else if ((count - 13) % 9 == 0)
            {
                R_assii(buffer, dpoint.pos[0]);
            }
            else if ((count - 14) % 9 == 0)
            {
                R_assii(buffer, dpoint.pos[1]);
            }
            else if ((count - 15) % 9 == 0)
            {
                R_assii(buffer, dpoint.pos[2]);
            }
            else if ((count - 16) % 9 == 0)
            {
                R_assii(buffer, dpoint.vel[0]);
            }
            else if ((count - 17) % 9 == 0)
            {
                R_assii(buffer, dpoint.vel[1]);
            }
            else if ((count - 18) % 9 == 0)
            {
                R_assii(buffer, dpoint.vel[2]);
                d.push_back(dpoint);
            }
        }
    }
}
void Readfile_B(FILE *in, vector<Data> &b)
{
    Data bpoint;
    while (!feof(in))
    {
        //if(b.size()<20) cerr<<bpoint.galaxyid_column<<" "<<bpoint.particalid<<" "<<bpoint.pos[3]<<endl;
        R_binary(in, bpoint.galaxyid_column, 1);
        R_binary(in, bpoint.particalid, 1);
        R_binary(in, bpoint.pos[3], 1);
        b.push_back(bpoint);
    }
    b.pop_back();
}

void Calculate(ofstream &out, vector<Data> &c)
{
    out << "rp "
        << "density "
        << endl;

    double density[2 * N];
    double hx = 0.25 / N;
    for (int i = 0; i < 2 * N; i++)
    {
        density[i] = 0.0;
    }
    for (int i = 0; i < c.size(); i++)
    {
        int x = (int)(c[i].pos[3] / hx);
        if (x > 24)
            cerr << x << endl;
        density[x] += c[i].mass;
    }
    for (int i = 0; i < N; i++)
    {
        double s = (double)i;
        double r1 = s * 10.0;
        double r2 = (s + 1.0) * 10.0; //kpc
        s = 3.14 * (DD(r2) - DD(r1));
        double rp = 5 + 10 * (double)(i);
        density[i] /= s;

        out << rp << " " << log10(rp) << " " << density[i] << " " << log10(density[i]) << endl;
    }
}
void Calculate(ofstream &out, vector<Data> &c, double range, double vhx)
{
    out << "rp "
        << "density "
        << endl;

    double density[2 * N];
    double hx = 0.25 / N;
    for (int i = 0; i < 2 * N; i++)
    {
        density[i] = 0.0;
    }
    for (int i = 0; i < c.size(); i++)
    {
        double judge = c[i].vel[2] - range;
        if (judge > 0 && judge < vhx)
        {
            int x = (int)(c[i].pos[3] / hx);
            if (x > 24)
                cerr << x << endl;
            density[x] += c[i].mass;
        }
    }
    for (int i = 0; i < N; i++)
    {
        double s = (double)i;
        double r1 = s * 10.0;
        double r2 = (s + 1.0) * 10.0; //kpc
        s = 3.14 * (DD(r2) - DD(r1));
        double rp = 5 + 10 * (double)(i);
        density[i] /= s;

        out << rp << " " << log10(rp) << " " << density[i] << " " << log10(density[i]) << endl;
    }
}
void Calculate2(ofstream &out, vector<Data> &c, double range, double rphx)
{
    out << "vz "
        << "density "
        << endl;

    double density[2 * N];
    double hx = 2000 / N;
    for (int i = 0; i < 2 * N; i++)
    {
        density[i] = 0.0;
    }
    for (int i = 0; i < c.size(); i++)
    {
        double judge = c[i].pos[3] - range;
        if (judge > 0 && judge < rphx && c[i].vel[2] < 1000.0 && c[i].vel[2] > -1000.0)
        {
            int x = (int)((c[i].vel[2] + 1000.0) / hx);
            if (x > 24)
                cerr << x << endl;
            density[x] += c[i].mass;
        }
    }
    for (int i = 0; i < N; i++)
    {
        double s = (double)i;
        double r1 = range * 1000.0;
        double r2 = (range + rphx) * 1000.0; //kpc
        s = 3.14 * (DD(r2) - DD(r1));
        double vz = hx / 2.0 + hx * (double)(i)-1000.0;
        density[i] /= s;

        out << vz << " " << log10(vz) << " " << density[i] << " " << log10(density[i]) << endl;
    }
}
void Deal(vector<Data> &c, vector<Data> &s, int filenum)
{
    char filename[100];
    sprintf(filename, "/data6/Hao_L/risk/column_density/data/columndensty_rp/1c_%03d_local.txt", filenum);
    ofstream out1(filename);
    sprintf(filename, "/data6/Hao_L/risk/column_density/data/columndensty_rp/1s_%03d_local.txt", filenum);
    ofstream out2(filename);
    Calculate(out1, c);
    Calculate(out2, s);
    out1.close();
    out2.close();
}
void Deal(vector<Data> &c, vector<Data> &s, int filenum, float mass_bin_start)
{
    char filename[100];
    sprintf(filename, "/data6/Hao_L/risk/column_density/data/columndensty_rp/halomass/1c_%03d_%.2f_local.txt", filenum, mass_bin_start);
    ofstream out1(filename);
    sprintf(filename, "/data6/Hao_L/risk/column_density/data/columndensty_rp/halomass/1s_%03d_%.2f_local.txt", filenum, mass_bin_start);
    ofstream out2(filename);
    Calculate(out1, c);
    Calculate(out2, s);
    out1.close();
    out2.close();
}
void Deal(vector<Data> &c, vector<Data> &s, int filenum, double range, double vhx)
{
    char filename[100];
    int label = (int)range;
    sprintf(filename, "/data6/Hao_L/risk/column_density/data/columndensty_rp/velocity/%03d/%d1c_%03d_local.txt", filenum, label, filenum);
    ofstream out1(filename);
    sprintf(filename, "/data6/Hao_L/risk/column_density/data/columndensty_rp/velocity/%03d/%d1s_%03d_local.txt", filenum, label, filenum);
    ofstream out2(filename);
    Calculate(out1, c, range, vhx);
    Calculate(out2, s, range, vhx);
    out1.close();
    out2.close();
}
void Deal2(vector<Data> &c, vector<Data> &s, int filenum, double range, double rphx)
{
    char filename[100];
    int label = (int)range;
    sprintf(filename, "/data6/Hao_L/risk/column_density/data/columndensty_rp/density_to_velocity/%03d/%f1c_local.txt", filenum, range);
    ofstream out1(filename);
    sprintf(filename, "/data6/Hao_L/risk/column_density/data/columndensty_rp/density_to_velocity/%03d/%f1s_local.txt", filenum, range);
    ofstream out2(filename);
    Calculate2(out1, c, range, rphx);
    Calculate2(out2, s, range, rphx);
    out1.close();
    out2.close();
}
#endif
