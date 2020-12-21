#include <iostream>
#include "include_all.h"
#define N2 30
//#define deltah 0.01;
//#define PALACE
#define Nions 9
#define LINES 4995
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
    vector<double> halomass;
    vector<int> flag;
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

        in >> buffer; //1
        in >> buffer; //2 SFR
        in >> buffer; //3 x
        in >> buffer; //4 y
        in >> buffer; //5 z
        in >> buffer; //6 halomass
        R_assii(buffer, r);
        halomass.push_back(log10(r));
        in >> buffer; //7 stellarmass
        in >> buffer; //8 flag
        R_assii(buffer, point);
        flag.push_back(point);
        in >> buffer; //9 redshift
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
    double width[Nions] ;
};

void deal(int &ion, Ion_resolves &ion_info, Palace &tab, int &galaxymin, int &galaxymax)
{

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

    double galaxies = (double)(galaxymax - galaxymin) + 1.0;

    Galaxy_info massga;
    char filenamemassga[100];
    sprintf(filenamemassga, "%s/galaxypara.txt", result_path);
    massga.load(filenamemassga);

    char fileionall[100];
    sprintf(fileionall, "%s/width/ew_r_all.data", result_path);
    ofstream out1(fileionall);

    for (int galaxy = galaxymin; galaxy <= galaxymax; galaxy++)
    {
        /*
        if (massga.flag[galaxy] == 1)
        {
            continue;
        }
        */
        resultout.clear();
        char fileions[100];
        sprintf(fileions, "%s/width/ew_r_all_%03d.data", result_path, galaxy);
        ifstream in1(fileions);

        while (!in1.eof())
        {
            char buffer[100];
            double value;
            in1 >> buffer; //1
            R_assii(buffer, value);
            ropoint.r = value;

            for (int i = 0; i < Nions; i++)
            {
                in1 >> buffer; //2-10
                R_assii(buffer, value);
                ropoint.width[i] = value;
            }

            resultout.push_back(ropoint);
        }
        if (resultout.size() != 4995)
        {
            cerr << "yes " << resultout.size() << endl;
            resultout.pop_back();
        }
        cerr << galaxy << endl;

        for (int i = 0; i < resultout.size(); i++)
        {
            out1 << resultout[i].r << " ";
            for (int j = 0; j < Nions; j++)
            {
                out1 << resultout[i].width[j] << " ";
            }
            out1 << endl;
        }
    }
    out1.close();
}
//./assemble 0924 7 1 60
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
    
        cerr << " 11" << endl;
        Palace tab;
        cerr << " 111" << endl;
        deal(ion, ion_info, tab, galaxymin, galaxymax);
    
}
