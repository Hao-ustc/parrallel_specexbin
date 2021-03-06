#include <iostream>
#include "include_all.h"
#define N2 30
//#define deltah 0.01;
//#define PALACE
#define Nions 9
#define LINES 2495
#define DELTR 0.005
#define NUM_THREADS 100

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
    void load(int galaxy);
    void clear();
};
void Palace::clear()
{
    x.clear();
    y.clear();
    z.clear();
    width.clear();
}
void Palace::load(int galaxy)
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

    deltah = 3.0 / N2;

    for (int i = 10; i < N2 + 1; i++)
    {

        z.push_back(deltah * i);
        width.push_back(0.0);
    }
}

class Galaxy_info
{
public:
    vector<double> R_Vir;
    vector<double> rspash2RVir;
    void load(char file[]);
};

void Galaxy_info::load(char file[])
{
    ifstream in(file);
    while (!in.eof())
    {
        char buffer[80];
        int point = 0;
        double r = 0.0;
        in >> buffer;
        R_assii(buffer, r);
        R_Vir.push_back(r);
        in >> buffer;
        R_assii(buffer, r);
        rspash2RVir.push_back(r);
    }
    in.close();
}

class Data_resolve
{

public:
    double Redshift;
    double tau;
    double lambta;
    double flux;
};

class Out
{

public:
    double r;
    double r2rvir;
    double r2rspash;
    double width;
};

double deal(int &ion, Ion_resolves &ion_info, Palace &tab, int &galaxymin, int &galaxymax)
{
    int countrun=0;
    int countrunall=0;

    cerr << "1" << endl;
    vector<Out> resultout;
    resultout.clear();
    Out ropoint;

    vector<Data_resolve> data;
    Data_resolve p;

    int row = 4 * (Nions + 1) + 3;

    char result_path[100];
    sprintf(result_path, "./result/%s", spece_set.TIME);

    sprintf(spece_set.file_redshift, "/data6/Hao_L/result/pure_redshift.txt");
    sprintf(spece_set.spec_ion_filename, "%sspecions_i9.dat", spece_set.prefix);

    sprintf(filename, "%s/width/test.data", result_path);
    ofstream out1(filename);

    double EW_x[LINES];
    double EW_y[LINES];
    for (int i = 0; i < LINES; i++)
    {
        EW_x[i] = 0.0;
        EW_y[i] = 0.0;
    }
    double galaxies = (double)(galaxymax - galaxymin) + 1.0;

    //Galaxy_info galaxyR;
    char filegalaxyR[100];
    //sprintf(filegalaxyR, "%s/galaxypara1.txt", result_path);
    cerr << "here" << endl;
    //galaxyR.load(filegalaxyR);
    for (int galaxy = galaxymin; galaxy <= galaxymax; galaxy++)
    {

        cerr << galaxy << endl;
        tab.load(galaxy);
        int read_flag = 0;
        int count = 0;
        int count_all = 0;
        char los[200];
        sprintf(los, "%s/los%03d.txt", result_path, galaxy);
        if ((spece_set.LOSfile = fopen(los, "r")) == NULL)
        {
            cerr << "Could not open file " << los << endl;
            exit(-100);
        }
        while (!feof(spece_set.LOSfile)) //&& read_flag < 1
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

            double redshift = spece_set.redshift_center;
            double r = 1000 * tab.x[read_flag];
            ropoint.r = r / 0.74 / (1 + redshift);
            //ropoint.r2rvir = r/0.72 / (1000 * galaxyR.R_Vir[galaxy]); //output the radius
            //ropoint.r2rspash = ropoint.r2rvir / galaxyR.rspash2RVir[galaxy];
            //spece_set.redshift_center = 0.0;

            sprintf(filename, "%s/galaxy%03d/tau_%f_%f_%f_%f.txt", result_path, galaxy, spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
            ifstream in(filename);
            count_all++;
            countrunall++;

            if (!(in.is_open()))
            {
                count++;
                countrun++;
                //out3 << galaxy << " " << read_flag << endl
                ;
            }
            else
                in.close();

            read_flag++;
        }
        out1 << galaxy << " " << count << " " << count_all << endl;
        fclose(spece_set.LOSfile);
    }
    out1.close();
    double completeness=(double)(countrun)/(double)(countrunall);
    return completeness;
}
//./data_resolve 0910 -1 10 30 will be run 10 11 ......29.
int main(int argc, char *argv[])
{
    char filetest[200];
    sprintf(filetest,"./%s_%s.test",argv[1],argv[2]);
    ofstream outtest(filetest);


    char rundatec[5];
    int start;
    int end;
    R_assii(argv[1], start);
    R_assii(argv[2], end);
    for (int rundate = start; rundate <= end; rundate++)
    {
        sprintf(rundatec, "%04d", rundate);
        sprintf(spece_set.TIME, "%s", rundatec);
        //sprintf(spece_set.galaxy, "%s", argv[2]);
        //sprintf(los, "./result/%s/los%s.txt",argv[1],argv[2]);
        Ion_resolves ion_info;
        ion_info.load();
        cerr << "ion load completed" << endl;

        int ion = 0;
      
        cerr << ion << endl;
        int galaxymin = 0;
        int galaxymax = 0;
        R_assii(argv[3], galaxymin);
        R_assii(argv[4], galaxymax);
        if (ion == -1)
        {
            cerr << " 1" << endl;
            for (ion = 0; ion < 9; ion++)
            {
                Palace tab;
                deal(ion, ion_info, tab, galaxymin, galaxymax);
            }
        }
        else if (ion >= 0)
        {
            cerr << " 11" << endl;
            Palace tab;
            cerr << " 111" << endl;
            outtest<<rundate<<" "<< 1.0-deal(ion, ion_info, tab, galaxymin, galaxymax)<<endl;;
        }
    }
    outtest.close();
}
