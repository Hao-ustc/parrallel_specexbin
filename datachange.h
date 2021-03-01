#ifndef _DATAFIT_H
#define _DATAFIT_H
#include "voidA_include_all.h"

void datafit(vector<GAS> &gp, Setting &spece_set, datahead b)
{
    double h = spece_set.spece_para.h;
    double vel_corr;
    vel_corr = sqrt(b.Time * b.Time * b.Time) / h;
    double unit_vel = spece_set.spece_para.unit_Velocity;
    double unit_mass = spece_set.spece_para.unit_Mass;
    double pos_corr = spece_set.boxsize;

    for (int i = 0; i < gp.size(); i++)
    {

        for (int j = 0; j < 3; j++)
        {
            gp[i].vel[j] = gp[i].vel[j] / sqrt(b.Time) * 1e+05 / unit_vel;
#ifdef VELOCITY_UNIT_CORRECTION
            gp[i].vel[j] = gp[i].vel[j] / vel_corr;

#endif
            gp[i].pos[j] /= _UNIT_L_;
            gp[i].pos[j] -= 250.0;
        }

        double px = gp[i].pos[0] ;
        double py = gp[i].pos[1] ;
        double pz = gp[i].pos[2] ;
        //impact factor
        float temppos = DD(px, py, pz);
        temppos = pow(temppos, 0.5);
        gp[i].pos[3] = temppos;
        
        gp[i].pos[0] = gp[i].pos[0] / pos_corr;
        gp[i].pos[1] = gp[i].pos[1] / pos_corr;
        gp[i].pos[2] = gp[i].pos[2] / pos_corr;
        gp[i].pos[3] = gp[i].pos[3] / pos_corr;
        gp[i].mass = gp[i].mass * 1e-10;
        gp[i].mass = gp[i].mass * 6.77e-31 / (1.88e-29 * pow(pos_corr, 3));

        gp[i].Hsml = gp[i].Hsml / pos_corr;

        gp[i].Hsml = gp[i].Hsml * 0.5 / 1.2275; //?????????????????????????????????????????

        gp[i].density = gp[i].density * (6.77e-31 / 1.88e-29); //right

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
}

#endif