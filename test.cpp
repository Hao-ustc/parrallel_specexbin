#include "include_all.h"

using namespace std;

int main()
{
    vector<Gas_1> gp;
    //Readdata(200, gp);
    double max = 0.0;
    double min = 0.0;
    char filenametest[200];

    char file_redshift[200];
    sprintf(file_redshift, "/data6/Hao_L/result/pure_redshift.txt");
    int File = 99;
    double Z = -1.0;
    ifstream in(file_redshift);
    double redshift_movie = 0.0;

    for (int i = 200; i >= File; i--)
    {

        char buffer[80];
        in >> buffer;
        R_assii(buffer, redshift_movie);
    }
    cerr<<redshift_movie<<endl;
    in.close();
    /*
    for (int j = 0; j < -1; j++)
    {

        FILE *tf;
        sprintf(filenametest, "/data6/Hao_L/specexbin_huang/test%d.txt", j);
        if ((tf = fopen(filenametest, "wb")) == NULL)
        {
            cerr << "wrong" << endl;
            return -1;
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
            cerr << gp[i].pos[3] << endl;
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
    */
}