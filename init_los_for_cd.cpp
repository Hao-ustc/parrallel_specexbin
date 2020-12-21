#include "include_list.h"
#include "gdata.h"
using namespace std;
#define N2 30
#define Pi 3.1415926
//#define deltah 0.01;
#define PALACE
#define RAN

int main(int argc, char *argv[])
{

    double min = 0;
    double max = 0;
    int filenum = 200;
    int galaxy_flag = 0;
    R_assii(argv[2], min);
    R_assii(argv[3], max);
    R_assii(argv[4], filenum);
    R_assii(argv[5], galaxy_flag); //1:central
    double deltah = 3.0 / N2;

    char result_path[100];

    int direction;
    double redshift_center;
    double xspec, yspec, zspec;
    double boxsize = 60;

    //try to be a function
    vector<Galaxy> ga;

    Readdata(filenum, ga);

    double mass_max = 0;
    vector<int> ga_choose;
    //ga_choose.push_back(-1);
    for (int i = 0; i < ga.size(); i++)
    {
        if (log10(ga[i].halomass) > min && log10(ga[i].halomass) < max)
        {
            if ((ga[i].center_coordinates[2] - posx) < 5 && (ga[i].center_coordinates[2] - posx) > -5)
            {
                if (galaxy_flag != 0)
                {
                    if (log10(ga[i].Mass) > 10.2 && log10(ga[i].Mass) < 10.6 && ga[i].flag == galaxy_flag)
                    {
                        ga_choose.push_back(i);
                    }
                }
                else if (galaxy_flag == 0)
                {
                    if (log10(ga[i].Mass) > 10.2 && log10(ga[i].Mass) < 10.6)
                    {
                        ga_choose.push_back(i);
                    }
                }
                
            }
        }
    }
    cerr << ga_choose.size() << endl;
    //cerr << "galaxy " << ga[ga_choose[0]].id << " is choosen. The mass of it is " << ga[ga_choose[0]].halomass << endl;

    sprintf(result_path, "/data6/Hao_L/risk/column_density/data/map/%s", argv[1]);
    char los[100];
    char fileplace[200];

    sprintf(fileplace, "%s/palace%03d.txt", result_path, 0);
    ofstream tf(fileplace);
    if (!(tf.is_open()))
    {
        cerr << mkdir(result_path, S_IRWXU) << endl;
        for (int i = 0; i < ga_choose.size(); i++)
        {
            sprintf(result_path, "/data6/Hao_L/risk/column_density/data/map/%s/galaxy%03d", argv[1], i);
            cerr << mkdir(result_path, S_IRWXU) << endl;
        }
    }
    tf.close();
    sprintf(result_path, "/data6/Hao_L/risk/column_density/data/map/%s", argv[1]);

#ifndef PALACE
    for (int i = 0; i < ga_choose.size(); i++)
    {
        FILE *LOSfile;
        sprintf(los, "%s/los%03d.txt", result_path, i);
        cerr << los << endl;

        if ((LOSfile = fopen(los, "w")) == NULL)
        {
            cerr << "Could not open file " << los << endl;
            return 0;
        }
        sprintf(fileplace, "%s/palace%03d.txt", result_path, i);
        ofstream tf(fileplace);

        for (int j = 0; j < N2 + 1; j++)
        {
            double dx = (j)*deltah;
            double dy = 0;
            for (int l = -1; l < 2; l += 2)
            {
                tf << l * dx << " " << dy << endl;
                xspec = ga[ga_choose[i]].center_coordinates[0] + l * (pow(10, fabs(dx)) / 1000.0);
                yspec = ga[ga_choose[i]].center_coordinates[1];
                zspec = ga[ga_choose[i]].center_coordinates[2];
                direction = 2;
                redshift_center = 0.0;
                xspec -= 250;
                yspec -= 250;
                zspec -= 250;

                xspec /= boxsize;
                yspec /= boxsize;
                zspec /= boxsize;
                fprintf(LOSfile, "%lf %lf %lf %lf %d\n", redshift_center, xspec, yspec, zspec, direction);
            }
        }
        for (int k = 0; k < N2 + 1; k++)
        {
            double dx = 0;
            double dy = (k)*deltah;

            for (int l = -1; l < 2; l += 2)
            {
                tf << dx << " " << l * dy << endl;
                xspec = ga[ga_choose[i]].center_coordinates[0];
                yspec = ga[ga_choose[i]].center_coordinates[1] + l * (pow(10, fabs(dy)) / 1000.0);
                zspec = ga[ga_choose[i]].center_coordinates[2];
                direction = 2;
                redshift_center = 0.0;
                xspec -= 250;
                yspec -= 250;
                zspec -= 250;

                xspec /= boxsize;
                yspec /= boxsize;
                zspec /= boxsize;
                fprintf(LOSfile, "%lf %lf %lf %lf %d\n", redshift_center, xspec, yspec, zspec, direction);
            }
        }
        tf.close();
        fclose(LOSfile);
    }
#else
#ifdef RAN
    srand((unsigned)time(NULL));
    for (int i = 0; i < ga_choose.size(); i++)
    {
        FILE *LOSfile;
        sprintf(los, "%s/los%03d.txt", result_path, i);
        cerr << los << endl;

        if ((LOSfile = fopen(los, "w")) == NULL)
        {
            cerr << "Could not open file " << los << endl;
            return 0;
        }
        sprintf(fileplace, "%s/palace%03d.txt", result_path, i);
        ofstream tf(fileplace);
        double rangep = 0.25; //Mpc
        int N_1 = 200;
        double deltar = (double)(rangep / N_1);
        int N_2 = 10;
        double deltatheta = (double)360 / N_2;

        for (int k = 1; k < N_1; k++)
        {
            double r = ((double)rand()) / (double)RAND_MAX; //Mpc
            r = r * 0.25;
            for (int j = 0; j < N_2; j++)
            {

                double theta = (double)j * deltatheta / 180 * Pi;
                double dx = r * cos(theta);
                double dy = r * sin(theta);
                tf << r << " " << theta << endl;

                xspec = ga[ga_choose[i]].center_coordinates[0] + dx;
                yspec = ga[ga_choose[i]].center_coordinates[1] + dy;
                zspec = ga[ga_choose[i]].center_coordinates[2];
                direction = 2;
                redshift_center = 0.0;
                xspec -= 250;
                yspec -= 250;
                zspec -= 250;

                xspec /= boxsize;
                yspec /= boxsize;
                zspec /= boxsize;
                fprintf(LOSfile, "%lf %lf %lf %lf %d\n", redshift_center, xspec, yspec, zspec, direction);
            }
        }

        tf.close();
        fclose(LOSfile);
    }
#else
    for (int i = 0; i < ga_choose.size(); i++)
    {
        FILE *LOSfile;
        sprintf(los, "%s/los%03d.txt", result_path, i);
        cerr << los << endl;

        if ((LOSfile = fopen(los, "w")) == NULL)
        {
            cerr << "Could not open file " << los << endl;
            return 0;
        }
        sprintf(fileplace, "%s/palace%03d.txt", result_path, i);
        ofstream tf(fileplace);
        double rangep = 1; //Mpc
         int N_1 = 200;
        R_assii(argv[6], N_1);
        double deltar = (double)(rangep / N_1);
        int N_2 = 5;
        R_assii(argv[7], N_2);
        double deltatheta = (double)360.0 / N_2;

        for (int k = 1; k < N_1; k++)
        {
            for (int j = 0; j < N_2; j++)
            {
                double r = (double)k * deltar; //Mpc
                double theta = (double)j * deltatheta / 180.0 * Pi;
                double dx = r * cos(theta);
                double dy = r * sin(theta);
                tf << r << " " << theta << endl;

                xspec = ga[ga_choose[i]].center_coordinates[0] + dx;
                yspec = ga[ga_choose[i]].center_coordinates[1] + dy;
                zspec = ga[ga_choose[i]].center_coordinates[2];
                direction = 2;
                redshift_center = 0.0;
                xspec -= 250;
                yspec -= 250;
                zspec -= 250;

                xspec /= boxsize;
                yspec /= boxsize;
                zspec /= boxsize;
                fprintf(LOSfile, "%lf %lf %lf %lf %d\n", redshift_center, xspec, yspec, zspec, direction);
            }
        }

        tf.close();
        fclose(LOSfile);
    }
#endif
#endif

}