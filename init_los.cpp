#include "include_list.h"
#include <set>

using namespace std;

#define Pi 3.1415926
//#define deltah 0.01;
//#define PALACE
#define RAN //ROUND
class RUNNING_init_los
{
};
//./init_los 0629 12.0 12.5 12.0 12.5 100 0 0.25 99 1 5 
int main(int argc, char *argv[])
{
    double min = 0;
    double max = 0;
    double min1 = 0.0;
    double max1 = 0.0;
    int N_lines = 0;
    int filenum = 200;
    int galaxy_flag = 0;

    double rangepmin = 0.0; //Mpc
    double rangepmax = 1.0; //Mpc
    int effective_number = 5;
    //double effctive_controler=0.0;

    R_assii(argv[2], min); //halomass_min
    R_assii(argv[3], max);
    R_assii(argv[4], min1); //stellarmass_min
    R_assii(argv[5], max1);
    R_assii(argv[6], N_lines); //n_rbins
    R_assii(argv[7], rangepmin);
    R_assii(argv[8], rangepmax); //range of rp
    R_assii(argv[9], filenum);
    R_assii(argv[10], galaxy_flag); //1:central
    R_assii(argv[11], effective_number);
    double r_range = rangepmax - rangepmin;
    if (r_range <= 0)
    {
        cerr << "what had happened on r???????" << endl;
        exit(-1);
    }

    effective_number = (double)pow(10, effective_number);

    char result_path[100];
    sprintf(result_path, "./result/%s", argv[1]);

    int direction;
    double redshift_center;
    Z_File t;
    sprintf(t.file_redshift, "/data6/Hao_L/result/simulation/voidA/redshiftmovie.txt");
    redshift_center = t.Check_File_Z(filenum);
    double rhocz = (2.77e11) * (1 - 0.258 + 0.258 * pow((1 + redshift_center), 3.0)) / (pow(1.0 + redshift_center, 3.0));
    double xspec, yspec, zspec;
    double boxsize = 80;

    //try to be a function
    GALAXYS galaxys;
    galaxys.load(filenum);
    vector<GALAXY> ga;
    ga.clear();
    ga.swap(galaxys.galaxys);



    char fileplace[200];

    sprintf(fileplace, "%s/palace%03d.txt", result_path, 0);
    ofstream tf(fileplace);
    if (!(tf.is_open()))
    {
        char folder[100];
        cerr << mkdir(result_path, S_IRWXU) << endl;
        sprintf(folder, "./result/0_a_width/%s", argv[1]);
        cerr << mkdir(folder, S_IRWXU) << endl;
        sprintf(folder, "%s/a_los", result_path);
        cerr << mkdir(folder, S_IRWXU) << endl;
       
    }
    tf.close();//prepare mkdir

    char filetest[100];
    sprintf(filetest, "%s/paratest.txt", result_path);
    ofstream outtest(filetest);

    vector<int> ga_choose;
    //ga_choose.push_back(-1);
    for (int i = 0; i < ga.size(); i++)
    {
        //output
        double rc3 = 3.0 / (4.0 * Pi) * ga[i].halo.submass / (200.0 * rhocz);
        double rc = pow(rc3, 1.0 / 3.0);
        outtest << ga[i].id << " "
                << ga[i].halo.submass << " "
                << ga[i].Rvir << " "
                << rc << endl;
        if (ga[i].SFR < 1.0 || ga[i].SFR > 20.0)
            continue;
        if (log10(ga[i].halo.submass) > min && log10(ga[i].halo.submass) < max)
        {
            //if ((ga[i].center_coordinates[2] - posx) < 5 && (ga[i].center_coordinates[2] - posx) > -5)
            //{
            if (galaxy_flag != 0)
            {
                if (log10(ga[i].stellarmass) > min1 && log10(ga[i].stellarmass) < max1 && ga[i].flag == galaxy_flag)
                {
                    ga_choose.push_back(i);
                }
            }
            else if (galaxy_flag == 0)
            {
                if (log10(ga[i].stellarmass) > min1 && log10(ga[i].stellarmass) < max1)
                {
                    ga_choose.push_back(i);
                }
            }
            // }
        }
    }
    outtest.close();
    cerr << ga_choose.size() << endl;
    //cerr << "galaxy " << ga[ga_choose[0]].id << " is choosen. The mass of it is " << ga[ga_choose[0]].halo.submass << endl;

   
    
     for (int i = 0; i < ga_choose.size(); i++)
        {
            char folder[100];
            int j = ga_choose[i];
            sprintf(folder, "./result/%s/galaxy%05d", argv[1], i);
            cerr << mkdir(folder, S_IRWXU) << endl;
        }

#ifdef PA
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
#ifdef PALACE
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
        for (int l = -1; l < 2; l += 2)
        {
            for (int j = 1; j < N2 + 1; j++)
            {

                for (int ll = -1; ll < 2; ll += 2)
                {
                    for (int k = 1; k < N2 + 1; k++)
                    {
                        double dx = l * (j)*deltah;
                        double dy = ll * (k)*deltah;

                        tf << dx << " " << dy << endl;

                        xspec = ga[ga_choose[i]].center_coordinates[0] + l * (pow(10, fabs(dx)) / 1000.0);
                        yspec = ga[ga_choose[i]].center_coordinates[1] + ll * (pow(10, fabs(dy)) / 1000.0);
                        zspec = ga[ga_choose[i]].center_coordinates[2];
                        direction = 2;

                        xspec -= posx;
                        yspec -= posy;
                        zspec -= posz;

                        xspec /= boxsize;
                        yspec /= boxsize;
                        zspec /= boxsize;
                        fprintf(LOSfile, "%lf %lf %lf %lf %d\n", redshift_center, xspec, yspec, zspec, direction);
                    }
                }
            }
        }
        tf.close();
        fclose(LOSfile);
    }
#endif
#ifdef ROUND
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
        double deltatheta = (double)360 / N_2;

        for (int k = 1; k < N_1; k++)
        {
            for (int j = 0; j < N_2; j++)
            {
                double r = (double)k * deltar; //Mpc
                double theta = (double)j * deltatheta / 180 * Pi;
                double dx = r * cos(theta);
                double dy = r * sin(theta);
                tf << r << " " << theta << endl;

                xspec = ga[ga_choose[i]].center_coordinates[0] + dx;
                yspec = ga[ga_choose[i]].center_coordinates[1] + dy;
                zspec = ga[ga_choose[i]].center_coordinates[2];
                direction = 2;

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
#ifdef RAN
    set<int> xset;
    set<int> yset;
    
    // ser<>
    xset.clear();
    yset.clear();
    N_lines *= ga_choose.size();
    srand((unsigned)time(NULL));

    int setcountx = xset.size();
    int setcounty = yset.size();

    if (setcountx != 0 || setcounty != 0)
    {
        exit(-1);
    }
    FILE *LOSfile;
     char los[100];
    sprintf(los, "%s/los.txt", result_path);
    cerr << los << endl;
    if ((LOSfile = fopen(los, "w")) == NULL)
    {
        cerr << "Could not open file " << los << endl;
        return 0;
    }
    sprintf(fileplace, "%s/palace.txt", result_path);
    ofstream tfpa(fileplace);

    for (int i = 0; i < N_lines; i++)
    {
        double r = 0.0;
        r = ((double)rand()) / (double)RAND_MAX; //Mpc

        r = pow(r, 0.5);
        r = r * rangepmax;
        if (r < rangepmin)
        {
            i--;
            continue;
        }

        double theta = ((double)rand()) / (double)RAND_MAX;
        theta *= (2.0 * Pi);

        int gacindex = rand() % ga_choose.size();

        double dx = r * cos(theta);
        double dy = r * sin(theta);

        xspec = ga[ga_choose[gacindex]].center_coordinates[0] + dx;
        yspec = ga[ga_choose[gacindex]].center_coordinates[1] + dy;
        zspec = ga[ga_choose[gacindex]].center_coordinates[2];
        direction = 2;
        xspec -= 250;
        yspec -= 250;
        zspec -= 250;

        xspec /= boxsize;
        yspec /= boxsize;
        zspec /= boxsize;

        int tempx = xspec * effective_number;
        int tempy = yspec * effective_number;
        xset.insert(tempx);
        yset.insert(tempy);
        if (xset.size() == setcountx && yset.size() == setcounty)
        {
            i--;
            continue;
        }
        else if (xset.size() == (setcountx + 1) || yset.size() == (setcounty + 1))
        {
            setcountx = xset.size();
            setcounty = yset.size();
        }
        else
        {
            cerr << "we must exit now!!!!" << endl;
            cerr << xset.size() << " " << yset.size() << " "
                 << setcountx << " " << setcounty << endl;
            exit(-1);
        }

        fprintf(LOSfile, "%lf %lf %lf %lf %d\n", redshift_center, xspec, yspec, zspec, direction);
        tfpa << r << " "
           << r / ga[ga_choose[gacindex]].halo.Rvir << " " << gacindex << endl;
    }
    tfpa.close();
    fclose(LOSfile);
#endif

    char filename_para[200];
    sprintf(filename_para, "%s/para.txt", result_path);
    ofstream outpara(filename_para);
    outpara 
            << "galaxy number " << ga_choose.size() << endl
            << "workpalace " << argv[1] << endl
            << "halomass bin  " << min << "~" << max << endl
            << "stellarmass bin " << min1 << "~" << max1 << endl
            << "rangep " << rangepmin << "~" << rangepmax << endl
            << "lines " << N_lines << endl
            << "filenum " << filenum << " " << redshift_center << endl
            << "galaxyflag " << galaxy_flag << endl;

    outpara << "./build_pbs " << argv[1] << " "  
            <<N_lines << endl;
    outpara.close();
    sprintf(filename_para, "%s/galaxypara.txt", result_path);
    ofstream outgpara(filename_para);

    for (int i = 0; i < ga_choose.size(); i++)
    {
        outgpara << i << " "
                 << ga[ga_choose[i]].SFR << " "
                 << ga[ga_choose[i]].center_coordinates[0] << " "
                 << ga[ga_choose[i]].center_coordinates[1] << " "
                 << ga[ga_choose[i]].center_coordinates[2] << " "
                 << ga[ga_choose[i]].halo.submass << " "
                 << ga[ga_choose[i]].stellarmass << " "
                 << ga[ga_choose[i]].flag << " "
                 << redshift_center << endl;
    }

    outgpara.close();
}