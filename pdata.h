#ifndef P_DATA_H
#define P_DATA_H
//#include "/data6/Hao_L/code/head/coma/include.h"
#include "include_all.h"
#define TRY
#define UNIT_M 1.989e43 /* 10^10 Msolar */
using namespace std;

class Gas_1
{
public:
    int id;
    int galaxyid_column;
    int stellarmassbin;
    float mass,
        pos[4],
        vel[3],
        Hsml,
        sfr,
        //  delaytime,
        T;
    //pos[3] is the impact factor
    float U;
    float density;
    float um;
    float Z[4];
    double galaxymass_column;
    Gas_1()
    {
        id = -1;
        galaxyid_column = -1;
        stellarmassbin = -1;

        mass = 0.0;
        //  delaytime = 0.0;
        Hsml = 0.0;
        sfr = 0.0;
        T = 0.0;
        U = 0.0;
        um = 0.0;
        galaxymass_column = 0.0;
    }

    ;
};

bifstream fi;
datahead b;
int blksize;

void skipblock(bifstream &fi)
{
    int blksize;
    fi >> blksize;
    //  CERR(blksize);
    fi.skip(blksize);
    fi >> blksize;
    return;
}

void skipblock(bifstream &fi, int num)
{
    for (int i = 0; i < num; i++)
        skipblock(fi);
    return;
}

void readhead(bifstream &fi, datahead &b)
{
    int a;
    fi >> a;
    //CERR(a);
    for (int i = 0; i < 6; i++)
        fi >> b.Npart[i];
    for (int i = 0; i < 6; i++)
        fi >> b.Massarr[i];
    fi >> b.Time;
    fi >> b.Redshift;
    fi >> b.FlagSfr;
    fi >> b.FlagFeedback;
    for (int i = 0; i < 6; i++)
        fi >> b.Nall[i];
    fi >> b.FlagCooling;
    fi >> b.NumFiles;
    fi >> b.BoxSize;
    fi >> b.Omega0;
    fi >> b.OmegaLambda;
    fi >> b.HubbleParam;
    fi >> b.FlagAge;
    fi >> b.FlagMetals;
    for (int i = 0; i < 6; i++)
        fi >> b.NallHW[i];
    fi >> b.flag_entr_ics;
    fi.skip(b.fill);
    fi >> a;
    //CERR(a);
    return;
}
int gasid[IDMAX];

void Prepare_partical(vector<Gas_1> &gp,Setting &spece_set)
{
    double h = spece_set.spece_para.h;
    double vel_corr;
    vel_corr = sqrt(b.Time * b.Time * b.Time) / h;
    double unit_vel = spece_set.spece_para.unit_Velocity;
    double unit_mass = spece_set.spece_para.unit_Mass;
    double pos_corr = spece_set.boxsize;
    for (int i = 0; i < gp.size(); i++)
    {
#ifndef TRY
        //gp[i].density = gp[i].density * h * h;
        gp[i].vel[0] = gp[i].vel[0] / sqrt(b.Time) * 1e+05 / unit_vel;
        gp[i].vel[1] = gp[i].vel[1] / sqrt(b.Time) * 1e+05 / unit_vel;
        gp[i].vel[2] = gp[i].vel[2] / sqrt(b.Time) * 1e+05 / unit_vel;
#ifdef VELOCITY_UNIT_CORRECTION
        gp[i].vel[0] = gp[i].vel[0] / vel_corr;
        gp[i].vel[1] = gp[i].vel[1] / vel_corr;
        gp[i].vel[2] = gp[i].vel[2] / vel_corr;
#endif

        gp[i].pos[0] = gp[i].pos[0] / pos_corr;
        gp[i].pos[1] = gp[i].pos[1] / pos_corr;
        gp[i].pos[2] = gp[i].pos[2] / pos_corr;
        gp[i].pos[3] = gp[i].pos[3] / pos_corr;

        gp[i].mass = gp[i].mass * 1e-10;
        gp[i].mass = gp[i].mass * 6.77e-31 / (1.88e-29 * pow(pos_corr, 3));

        gp[i].Hsml = gp[i].Hsml / pos_corr;

        gp[i].Hsml = gp[i].Hsml * 0.5 / 1.2275; //?????????????????????????????????????????

        gp[i].density = gp[i].density * (6.77e-31 / 1.88e-29);
#else
       // gp[i].density = gp[i].density * h * h;
        gp[i].vel[0] = gp[i].vel[0] / sqrt(b.Time) * 1e+05 / unit_vel;
        gp[i].vel[1] = gp[i].vel[1] / sqrt(b.Time) * 1e+05 / unit_vel;
        gp[i].vel[2] = gp[i].vel[2] / sqrt(b.Time) * 1e+05 / unit_vel;
#ifdef VELOCITY_UNIT_CORRECTION
        gp[i].vel[0] = gp[i].vel[0] / vel_corr;
        gp[i].vel[1] = gp[i].vel[1] / vel_corr;
        gp[i].vel[2] = gp[i].vel[2] / vel_corr;
#endif

        gp[i].pos[0] = gp[i].pos[0] / pos_corr;
        gp[i].pos[1] = gp[i].pos[1] / pos_corr;
        gp[i].pos[2] = gp[i].pos[2] / pos_corr;
        gp[i].pos[3] = gp[i].pos[3] / pos_corr;

        gp[i].mass = gp[i].mass * 1e-10;
        gp[i].mass = gp[i].mass * 6.77e-31 / (1.88e-29 * pow(pos_corr, 3));

        gp[i].Hsml = gp[i].Hsml / pos_corr;

        gp[i].Hsml = gp[i].Hsml * 0.5 / 1.2275; //?????????????????????????????????????????

        gp[i].density = gp[i].density * (6.77e-31 / 1.88e-29); //right
#endif
    }

    //Rotate90Box();
    int i, ri;
    double hold;

    ri = 0;

    for (i = 0; i < gp.size(); i++)
    {
        if (spece_set.direction == 1)
        {
            hold = gp[i].pos[2];
            gp[i].pos[2] = gp[i].pos[1];
            gp[i].pos[1] = gp[i].pos[0];
            gp[i].pos[0] = hold;

            hold = gp[i].vel[2];
            gp[i].vel[2] = gp[i].vel[1];
            gp[i].vel[1] = gp[i].vel[0];
            gp[i].vel[0] = hold;
            ri++;
        }

        if (spece_set.direction == 0)
        {
            hold = gp[i].pos[0];
            gp[i].pos[0] = gp[i].pos[1];
            gp[i].pos[1] = gp[i].pos[2];
            gp[i].pos[2] = hold;

            hold = gp[i].vel[0];
            gp[i].vel[0] = gp[i].vel[1];
            gp[i].vel[1] = gp[i].vel[2];
            gp[i].vel[2] = hold;
            ri++;
        }
    }
    fprintf(stderr, "Physically rotated %d particles, because direction= %d\n", ri, spece_set.direction);
}

void Readdata(int filenum, vector<Gas_1> &gas,Setting &spece_set)
{

    for (int i = 0; i < IDMAX; i++)
    {
        gasid[i] = -1;
    }
    int tempgasn = 0;

    for (int i = 0; i < 32; i++)
    {

        int id;
        float mass;
        float density;
        float age;
        float pos[3], vel[3];
        float delaytime;
        float Hsml;
        float sfr;
        float z[4];
        Gas_1 gaspoint;

        sprintf(filename, "/data11/huiwang/Coma/snapdir_%03d/snap_gw_%03d.%d", filenum, filenum, i);
        fi.open(filename, 0);
        fi.clear();
        fi.rewind();
        readhead(fi, b);

        //read id
        skipblock(fi, 1);
        skipblock(fi, 1);
        fi >> blksize;

        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            fi >> id;
            gaspoint.id = id;
            gas.push_back(gaspoint);
            gasid[id] = gas.size() - 1; //just the position in list
        }

        fi.skip(b.Npart[1] * 4);
        fi.skip(b.Npart[2] * 4);
        fi.skip(b.Npart[3] * 4);
        fi.skip(b.Npart[4] * 4);
        fi.skip(blksize - (b.Npart[0] + b.Npart[1] + b.Npart[2] + b.Npart[3] + b.Npart[4]) * 4);
        fi.skip(4);
        fi.rewind();
        readhead(fi, b);

        //flag
        char flag_yes = 1;
        char flag_no = 0;
        int flag_initvector = 0;

        //read pos

        fi >> blksize;

        for (int j = 0; j < b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> pos[0] >> pos[1] >> pos[2];

            gas[p].pos[0] = pos[0] - 250.0;
            gas[p].pos[1] = pos[1] - 250.0;
            gas[p].pos[2] = pos[2] - 250.0;
            double px = pos[0] - posx;
            double py = pos[1] - posy;
            double pz = pos[2] - posz;
            //impact factor
            float temppos = DD(px, py, pz);
            temppos = pow(temppos, 0.5);
            gas[p].pos[3] = temppos;
        }

        fi.skip(blksize - 3 * (b.Npart[0]) * sizeof(float));
        fi.skip(4);

        //skip vel

        fi >> blksize;

        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> vel[0] >> vel[1] >> vel[2];
            gas[p].vel[0] = vel[0];
            gas[p].vel[1] = vel[1];
            gas[p].vel[2] = vel[2];
        }
        fi.skip(blksize - 3 * (b.Npart[0]) * sizeof(float));
        fi.skip(4);

        cerr << i << endl;
        //id
        skipblock(fi, 1);

        //read mass

        fi >> blksize;

        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> mass;
            mass *= 1e10;
            gas[p].mass = mass;
        }

        fi.skip(blksize - (b.Npart[0]) * sizeof(float));
        fi.skip(4);
        //read U
        fi >> blksize;
        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> gas[p].U;
        }
        fi.skip(blksize - b.Npart[0] * sizeof(float));
        fi.skip(4);
        //read density
        fi >> blksize;
        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> gas[p].density;
        }
        fi.skip(blksize - b.Npart[0] * sizeof(float));
        fi.skip(4);
        //read um  (ne , not um)
        fi >> blksize;
        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> gas[p].um;
        }
        fi.skip(blksize - b.Npart[0] * sizeof(float));
        fi.skip(4);
        //skip NH
        skipblock(fi, 1);
        //read Hsml
        fi >> blksize;
        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> Hsml;
            gas[p].Hsml = Hsml;
        }
        fi.skip(blksize - b.Npart[0] * sizeof(float));
        fi.skip(4);
        //read SFR
        fi >> blksize;
        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> sfr;
            gas[p].sfr = sfr;
        }
        fi.skip(blksize - b.Npart[0] * sizeof(float));
        fi.skip(4);
        //read delaytime
        /*fi >> blksize;
        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> delaytime;
            gas[p].delaytime = delaytime;
        }
        fi.skip(blksize - b.Npart[0] * sizeof(float));
        fi.skip(4);
        */
        skipblock(fi, 1);
        //read age(2)
        //skipblock(fi, 1);
        if (b.Npart[4] != 0)
            skipblock(fi, 1);

        //read Z
        fi >> blksize;
        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            fi >> z[0] >> z[1] >> z[2] >> z[3];
            gas[p].Z[0] = z[0];
            gas[p].Z[1] = z[1];
            gas[p].Z[2] = z[2];
            gas[p].Z[3] = z[3];
        }
        fi.skip(blksize - b.Npart[0] * 4 * sizeof(float));
        fi.skip(4);
        //get T
        compute_T h2;
        for (int j = 0; j < 0 + b.Npart[0]; j++)
        {
            int p = j + tempgasn;
            gas[p].T = h2.get_T(gas[p].Z, gas[p].U, gas[p].um, b.Time);
            //	testf << gas[p].t << endl;
        }

        //get Tmax
        skipblock(fi, 1);

        fi.close();
        tempgasn = gas.size();
    }

    Prepare_partical(gas,spece_set);
}

#endif