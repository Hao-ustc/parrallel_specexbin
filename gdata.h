#ifndef G_DATA_H
#define G_DATA_H
#include <iostream>
#include "include_list.h"

using namespace std;

int galaxypos[GALAXYIDMAX];

void Readdata(int filenum, vector<Galaxy> &galaxy)
{
    for (int i = 0; i < GALAXYIDMAX; i++)
    {
        galaxypos[i] = -1;
    }
    char filename1[100];
    sprintf(filename1, "/data11/huiwang/ComaGG/so_z%03d.sovcirc", filenum);
    ifstream in1(filename1);
    cerr << "halo" << endl;

    vector<Halo> HALO;
    Halo hal;
    HALO.clear();
    int counthalo = 0;
    double htemp;

    while (!in1.eof())
    {
        char halo[20];
        in1 >> halo;
        counthalo++;

        if (counthalo > 9)
        {
            sscanf(halo, "%lf", &htemp);
            if ((counthalo - 10) % 8 == 0)
                hal.id = htemp;
            if ((counthalo - 11) % 8 == 0)
                hal.Mass = htemp/0.72;
            if ((counthalo - 12) % 8 == 0)
                hal.Rvir = htemp/0.72;

            if ((counthalo - 16) % 8 == 0)
            {
                hal.subMass = htemp/0.72;
                HALO.push_back(hal);
            }
        }
    }

    char filename2[100];
    sprintf(filename2, "/data11/huiwang/ComaGG/so_z%03d.par", filenum);
    ifstream in2(filename2);

    int countpair = 0;
    int pairtemp;

    while (!in2.eof())
    {
        int haloid;
        char pair[20];
        in2 >> pair;

        countpair++;
        sscanf(pair, "%d", &pairtemp);

        if ((countpair - 1) % 2 == 0)
            haloid = pairtemp;
        else if ((countpair - 2) % 2 == 0)
            HALO[haloid - 1].hostid = pairtemp;
    }

    //---galaxy----------------------------------------------------------------------
    char filename3[100];
    sprintf(filename3, "/data11/huiwang/ComaGG/gal_z%03d.stat", filenum);
    ifstream in3(filename3);

    int countgalaxy = 0;
    double gtemp;
    Galaxy gal;
    galaxy.clear();

    while (!in3.eof())
    {
        int i = 0;
        in3 >> gtemp;
        countgalaxy++;

        if ((countgalaxy - 1) % 22 == 0)
            gal.id = gtemp;
        else if ((countgalaxy - 2) % 22 == 0)
            gal.Npartical = gtemp;
        else if ((countgalaxy - 5) % 22 == 0)
            gal.Mass = gtemp * 4.8179e10;
        else if ((countgalaxy - 16) % 22 == 0)
            gal.vel[0] = gtemp * 1.7274e1;
        else if ((countgalaxy - 17) % 22 == 0)
            gal.vel[1] = gtemp * 1.7274e1;
        else if ((countgalaxy - 18) % 22 == 0)
            gal.vel[2] = gtemp * 1.7274e1;

        else if ((countgalaxy - 19) % 22 == 0)
        {
            gal.center_coordinates[0] = gtemp * 500.0 + 250.0;
        }

        else if ((countgalaxy - 20) % 22 == 0)
            gal.center_coordinates[1] = gtemp * 500.0 + 250.0;
        else if ((countgalaxy - 21) % 22 == 0)
            gal.center_coordinates[2] = gtemp * 500.0 + 250.0;

        else if ((countgalaxy - 22) % 22 == 0)
        {
            gal.SFR = gtemp;
            
            //if (HALO[gal.id - 1].subMass != 0)
            // {
            double gal_x = gal.center_coordinates[0] - posx;
            double gal_y = gal.center_coordinates[1] - posy;
            double gal_z = gal.center_coordinates[2] - posz;
            double r = gal_x * gal_x + gal_y * gal_y + gal_z * gal_z;
            r = pow(r, 0.5);
            gal.hostid = HALO[gal.id - 1].hostid;
            gal.halomass = HALO[gal.hostid-1].Mass;
            gal.Rvir = HALO[gal.hostid-1].Rvir;
            if (gal.hostid == gal.id)
                gal.flag = 1;
            else if (gal.hostid != gal.id)
                gal.flag = -1;
            if (r > 0 && r < 30)
            {
                galaxypos[(int)gal.id] = galaxy.size();
                galaxy.push_back(gal);
                //cerr << r << endl;
            }

            //else cerr << "error" << endl;
            // }
        }
    }

    cerr << "1" << endl;
    in1.close();
    in2.close();
    in3.close();
}

#endif