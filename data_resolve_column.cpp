#include <iostream>
#include "include_all.h"
#define N2 30
//#define deltah 0.01;
//#define PALACE
#define Nions 9
//#define LINES 5000
#define DELTR 0.005

using namespace std;

class Ion_resolve
{
public:
    char name[10];
    float lambda, fraction, Xsec, atomwt, bsys, alpha;
    int Zcolumn;
};

class Ion_resolves
{
public:
    char path_data[200];
    vector<Ion_resolve> ions;
    Ion_resolve ion;
    int load();

    Ion_resolves()
    {
        sprintf(path_data, "./prepare/ionfiles/specions_i9.dat");
    }
};
int Ion_resolves::load()
{
    int nions;
    FILE *specfile;
    char line[80];
    if ((specfile = fopen(path_data, "r")) == NULL)
    {
        fprintf(stderr, "cannot find specion file anywhere\n");
        exit(-1);
    }
    int i = 0;
    while (fgets(line, 80, specfile) != NULL)
    {
        if (strstr(line, "#") != NULL)
            continue;
        if (i >= MAXIONS)
            break;
        ions.push_back(ion);
        sscanf(line, "%10s %g %g %g %g %d %g", ions[i].name, &ions[i].lambda, &ions[i].Xsec, &ions[i].atomwt, &ions[i].fraction, &ions[i].Zcolumn, &ions[i].alpha);
        i++;
    }
    nions = i;
    fclose(specfile);

    fprintf(stderr, "Processing %d ions from specions.dat:\n", nions);
    for (i = 0; i < nions; i++)
    {

        ions[i].bsys = sqrt(2. * KBOLTZ / (MHYDR * ions[i].atomwt)) / 1.e5;
        ions[i].Xsec *= 2.648e-2 * ions[i].lambda * 1.e-13;
        fprintf(stderr, "%5d %10s %12.6g %12.6g %10.5g %10.5g %10.5g % 3.1f % 2d\n", i, ions[i].name, ions[i].lambda, ions[i].Xsec, ions[i].atomwt, ions[i].fraction, ions[i].bsys * sqrt(1.e4), ions[i].alpha, ions[i].Zcolumn);
    }
}

class Palace
{
public:
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> width;

    double deltah;
    void load(int &galaxy,int &LINES);
    void clear();
};
void Palace::clear()
{
    x.clear();
    y.clear();
    z.clear();
    width.clear();
}
void Palace::load(int &galaxy,int &LINES)
{

    char path[200];
    sprintf(path, "./result/%s/palace%03d.txt", spece_set.TIME, galaxy);
    ifstream in(path);
    int row = 2;
    int countpalace = 0;
    while (!in.eof())
    {
        char buffer[100];
        in >> buffer;
        countpalace++;
        int judge = countpalace % 2;
        if (judge == 1)
        {
            double dx = 0;
            R_assii(buffer, dx);
            x.push_back(dx);
        }
        else if (judge == 0)
        {
            double dy = 0;
            R_assii(buffer, dy);
            y.push_back(dy);
        }
    }
    in.close();
    while (x.size() - y.size())
    {
        x.pop_back();
        cerr << "get palace coordinate:     "
             << x.size() - y.size() << endl;
    }
    if(x.size()!=LINES)
    {
        cerr<<"size wrrong!!!!!!"<<x.size()<<endl;
        exit(-1);
    }
    deltah = 3.0 / N2;

    for (int i = 10; i < N2 + 1; i++)
    {

        z.push_back(deltah * i);
        width.push_back(0.0);
    }
}

class Data_resolve
{

public:
    double Redshift;
    double tau;
    double lambta;
    double flux;
};

void deal(int &ion, Ion_resolves &ion_info, Palace &tab, int &galaxymin, int &galaxymax,int &LINES)
{

    vector<Data_resolve> data;
    Data_resolve p;

    int row = 4 * (Nions + 1) + 3;

    char result_path[100];
    sprintf(result_path, "./result/%s", spece_set.TIME);

    sprintf(spece_set.file_redshift, "/data6/Hao_L/result/pure_redshift.txt");
    sprintf(spece_set.spec_ion_filename, "%sspecions_i9.dat", spece_set.prefix);
    char filename_column[200];
    sprintf(filename_column, "%s/width/column.data", result_path);
    ofstream out1(filename_column);
    sprintf(filename_column, "%s/width/rate_3.data", result_path);
    ofstream out2(filename_column);

    double galaxies = (double)(galaxymax - galaxymin) + 1.0;
    for (int galaxy = galaxymin; galaxy <= galaxymax; galaxy++)
    {
        cerr << galaxy << endl;
        tab.load(galaxy,LINES);
        int read_flag = 0;

        char los[200];
        sprintf(los, "%s/los%03d.txt", result_path, galaxy);
        if ((spece_set.LOSfile = fopen(los, "r")) == NULL)
        {
            cerr << "Could not open file " << los << endl;
            exit(-100);
        }
        while (!feof(spece_set.LOSfile) && read_flag < LINES) //&& read_flag < 1
        {

            spece_set.boxsize = 60; //in Mpc_h
            spece_set.flux_fac = -1.0;
            spece_set.spece_para.totMass = 0.258; //0.28; //0.238; //0.28;
            spece_set.spece_para.Lambda = 0.742;  //0.72; //0.762; //0.72;
            spece_set.spece_para.Omega_b = 0.045; //0.046; //0.0418; //0.046;
            spece_set.spece_para.H_0 = 73;        //73;
            spece_set.spece_para.h = 0.01 * spece_set.spece_para.H_0;
            spece_set.load();
            //cerr << spece_set.short_filenum << endl;
            double r = 1000 * tab.x[read_flag];
            out1 << r << " " << r / 0.72 / 2.0 << " ";
            out2 << r << " ";
            read_flag++;

            double redshift = spece_set.redshift_center;
            //spece_set.redshift_center = 0.0;
            char filename_efc[200];
            sprintf(filename_efc, "%s/galaxy%03d/efele_%f_%f_%f_%f.txt", result_path, galaxy, spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
            ifstream in1(filename_efc);
            
            /*
            sprintf(filename_efc, "%s/galaxy%03d/efc_%f_%f_%f_%f.txt", result_path, galaxy, spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
            ifstream in(filename_efc);
            
            sprintf(filename_efc, "%s/galaxy%03d/efc2_%f_%f_%f_%f.txt", result_path, galaxy, spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
            ifstream in2(filename_efc);
            */
            int incount = 0;
            double colh = 0.0;
            double colMg = 0.0;
            double colc = 0.0;
            double coH = 0.0;
            double coHe = 0.0;
            double coC3 = 0.0;
            double coC4 = 0.0;
            double coO4 = 0.0;
            double coO6 = 0.0;
            double coMg2 = 0.0;
            double coNe8 = 0.0;
            double coSi4 = 0.0;
            int oldpoint=0;
            while (!in1.eof())
            {
                char buffer[80];

                double mass = -1.0;
                ///*
                in1 >> buffer; //1
                R_assii(buffer, mass);
                coH += mass;
                in1 >> buffer; //2
                R_assii(buffer, mass);
                coHe += mass;
                in1 >> buffer; //3
                R_assii(buffer, mass);
                coC3 += mass;
                
                in1 >> buffer; //4
                R_assii(buffer, mass);
                coC4 += mass;
                in1 >> buffer; //5
                R_assii(buffer, mass);
                coO4 += mass;
                in1 >> buffer; //6
                R_assii(buffer, mass);
                coO6 += mass;
                in1 >> buffer; //7
                R_assii(buffer, mass);
                coNe8 += mass;
                in1 >> buffer; //8
                R_assii(buffer, mass);
                coMg2 += mass;
                in1 >> buffer; //9
                R_assii(buffer, mass);
                coSi4 += mass;
                
                /*
                in2 >> buffer;
                incount++;
                R_assii(buffer, mass);

                //int judge = incount % 2;

                if (  mass >= 0)//judge == 1&&
                {
                    //mass /= ion_info.ions[0].atomwt;
                    //mass *= 6.02e23;
                    col += mass;
                }
                else if (mass == -999)
                {
                    out1<<log10(col) << " ";
                    out2<<log10((col*0.72/60)/3.086e24/6e23) << " ";
                    col = 0.0;
                }
                */
            }
            //cerr<<"right"<<endl;
            //in.close();
            in1.close();
            //in2.close();
            out1 << colh << " " << coHe << " " << colc << " " << coC4 << " " << coO4 << " " << coO6 << " " << coNe8 << " " << colMg << " " << coSi4;
            //out1 << coH << " " << coHe << " " << coC3 << " " << coC4 << " " << coO4 << " " << coO6 << " " << //coNe8 << " " << coMg2 << " " << coSi4;
            //out1 << coH << " " << coHi << " " << coMg << " " << coMg2 << " " << coc << " " << coc3;
            //out2 << 1 / (coH / coHi) << " " << 1 / (coMg / coMg2) << " " << 1 / (coc / coc3);
            out1 << endl;
            out2 << endl;
            //out3.close();
            if (read_flag % 100 == 0)
                cerr << read_flag << endl;
        }

        fclose(spece_set.LOSfile);
        tab.clear();
    }

    out1.close();
    out2.close();
}
//./data_resolve_column 0629 0 55 64
int main(int argc, char *argv[])
{
    sprintf(spece_set.TIME, "%s", argv[1]);
    //sprintf(spece_set.galaxy, "%s", argv[2]);
    //sprintf(los, "./result/%s/los%s.txt",argv[1],argv[2]);
    Ion_resolves ion_info;
    ion_info.load();
    cerr << "ion load completed" << endl;

    int ion = 0;
    R_assii(argv[2], ion);
    cerr << ion << endl;
    int galaxymin = 0;
    int galaxymax = 0;
    R_assii(argv[3], galaxymin);
    R_assii(argv[4], galaxymax);
    int lines;
    R_assii(argv[5],lines);
    if (ion == -1)
    {
        cerr << " 1" << endl;
        for (ion = 0; ion < 9; ion++)
        {
            Palace tab;
            deal(ion, ion_info, tab, galaxymin, galaxymax,lines);
        }
    }
    else if (ion >= 0)
    {
        cerr << " 11" << endl;
        Palace tab;
        cerr << " 111" << endl;
        deal(ion, ion_info, tab, galaxymin, galaxymax,lines);
    }
}
