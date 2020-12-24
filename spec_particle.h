#ifndef SPEC_PARTICLE_H
#define SPEC_PARTICLE_H
#include "include_all.h"
#define TEST
using namespace std;
//class
class Spec_particle
{
public:
    int ID;
    float mass;
    float rho;
    float temp;
    float metals[NMETALS];
    float hsmooth;
    float pos[NDIM];
    float vel[NDIM];
    float sfr;
    float vlaunch;
    float ageaway;
    int nrec;
    float mgal_launch;
    float mgal;
    float dtravel;
    float dgal;
    float dpeculiar;
    float vrel;
};

class Spec_particles
{
public:
    vector<Spec_particle> spec_particles;
    Spec_particle sp;

    int Check_Partical_Los(vector<Gas_1> &gp,Setting &spece_set);
    void SmoothSpec(LOS &los,Ion_all &spece_ionall,Setting &spece_set);
};

int Spec_particles::Check_Partical_Los(vector<Gas_1> &gp,Setting &spece_set)
{
    spec_particles.clear();
    double xspec = spece_set.xspec;
    double yspec = spece_set.yspec;
    double zspec = spece_set.zspec;
    for (int i = 0; i < gp.size(); i++)
    {
        double dx, dy, dr2;
        double irep[NDIM];
        int int_hold;

        if (gp[i].pos[3] > 0.5)
            continue;

        //#ifdef SHORTSPEC
        dx = fabs(gp[i].pos[0] - spece_set.xspec);
        dy = fabs(gp[i].pos[1] - spece_set.yspec);
        /*#else
  RotateCoords(gp[i].pos[0], gp[i].pos[1], gp[i].pos[2], theta, phi);
  dx = fabs(xrot - xrotline);
  dy = fabs(yrot - yrotline);
#endif
*/
        for (irep[0] = -1; irep[0] <= 1; irep[0]++)
        {
            for (irep[1] = -1; irep[1] <= 1; irep[1]++)
            {
                for (irep[2] = -1; irep[2] <= 1; irep[2]++)
                {
                    if (fabs(gp[i].pos[0] + irep[0] * BOXSIZE + 2 * gp[i].Hsml) < HALFBOX || fabs(gp[i].pos[0] + irep[0] * BOXSIZE - 2 * gp[i].Hsml) < HALFBOX)
                    {
                        if (fabs(gp[i].pos[1] + irep[1] * BOXSIZE + 2 * gp[i].Hsml) < HALFBOX || fabs(gp[i].pos[1] + irep[1] * BOXSIZE - 2 * gp[i].Hsml) < HALFBOX)
                        {
                            //#ifdef SHORTSPEC
                            if (fabs(gp[i].pos[0] + irep[0] * BOXSIZE - xspec) < dx)
                            {
                                dx = fabs(gp[i].pos[0] + irep[0] * BOXSIZE - xspec);
                                // gp[i].pos[0] += irep[0]*BOXSIZE;
                            }
                            if (fabs(gp[i].pos[1] + irep[1] * BOXSIZE - yspec) < dy)
                            {
                                dy = fabs(gp[i].pos[1] + irep[1] * BOXSIZE - yspec);
                                // gp[i].pos[1] += irep[1]*BOXSIZE;
                            }
                            /*#else
            if (fabs(gp[i].pos[2] + irep[2] * BOXSIZE + 2 * gp[i].hsmooth) < HALFBOX || fabs(gp[i].pos[2] + irep[2] * BOXSIZE - 2 * gp[i].hsmooth) < HALFBOX)
            {
              RotateCoords(gp[i].pos[0] + irep[0] * BOXSIZE, gp[i].pos[1] + irep[1] * BOXSIZE, gp[i].pos[2] + irep[2] * BOXSIZE, theta, phi);
              if (fabs(xrot - xrotline) < dx)
                dx = fabs(xrot - xrotline);
              if (fabs(yrot - yrotline) < dy)
                dy = fabs(yrot - yrotline);
            }
#endif
*/
                        }
                    }
                }
            }
        }

        dr2 = dx * dx + dy * dy;
        if (dr2 < NSRCHRAD * gp[i].Hsml * NSRCHRAD * gp[i].Hsml)
        {
            sp.ID = gp[i].id;
            sp.mass = gp[i].mass;
            sp.temp = gp[i].T;
            sp.rho = gp[i].density;
            sp.hsmooth = gp[i].Hsml;
            sp.sfr = gp[i].sfr;

            sp.pos[0] = gp[i].pos[0];
            sp.pos[1] = gp[i].pos[1];
            sp.pos[2] = gp[i].pos[2];
            sp.vel[0] = gp[i].vel[0];
            sp.vel[1] = gp[i].vel[1];
            sp.vel[2] = gp[i].vel[2];

            sp.metals[0] = gp[i].Z[0];
            sp.metals[1] = gp[i].Z[1];
            sp.metals[2] = gp[i].Z[2];
            sp.metals[3] = gp[i].Z[3];

            //#ifdef SHORTSPEC
            sp.pos[0] -= xspec;
            sp.pos[1] -= yspec;
            sp.pos[2] -= zspec;
            //#endif
            spec_particles.push_back(sp);
        }
    }
    
    return (0);
}
void Spec_particles::SmoothSpec(LOS &los,Ion_all &spece_ionall,Setting &spece_set)

{
    char fileoutsmooth[200];
    sprintf(fileoutsmooth, "./process/smooth_120.txt");
    ofstream spf(fileoutsmooth);

    double irep[NDIM], part_pos[NDIM];
    float bound_min[NDIM], bound_max[NDIM];

    spece_set.redshift_center = los.bin[(int)((los.nzbins) / 2)].redshift;
    float redshift_center = spece_set.redshift_center;
    float t = CosmicTime(redshift_center,spece_set);
    cosmopar(t,spece_set);

    
    int n = los.nzbins;
    double *bin_coord;
    bin_coord = (double *)malloc((n + 1) * sizeof(double));
    for (int l = 0; l < n; l++)
    {
        bin_coord[l] = los.bin[l].relative_coord;
    }
    //int tempn = spece_ionall.nions;
    vector<ionfraction> ionfrac;
    ionfrac.clear();
   
    for (int i = 0; i < spece_ionall.nions; i++)
    {
         
        ionfraction tempionfrac;
        ionfrac.push_back(tempionfrac);
        ionfrac[i].loadtable(i + 1,spece_set);
        ionfrac[i].loadparam(b,spece_set);
        

    } //load fraction_lookback_table
    //maybe we want to output some files here
   
    bound_min[0] = bound_min[1] = bound_min[2] = -HALFBOX;
    bound_max[0] = bound_max[1] = bound_max[2] = HALFBOX;
    InitKernIntTable ktable;
    //omp_set_num_threads(NUM_THREADS);
    int countscp = 0;
    int countid = 0;
    double kernelcount = 0.0;
    int countscpall = 0;
    
    //#pragma omp parallel
    {
        int id = 0;
        //id = omp_get_thread_num();
        int counttest = 0;
        for (int i = id; i < spec_particles.size(); i += 1)
        {
            if ((i % 1000) == 0)
            {
                //#pragma omp critical
                {
                    time_t now_time = time(NULL);
                    tm *t_tm = localtime(&now_time);
                    spf << "local time is    :" << asctime(t_tm) << "        ";

                    spf << (float)(i) / (float)(spec_particles.size()) << endl;
                }
            }

            int j, k;
            Spec_particle scp;

            double vz;

            float distnorm;
            float unit_vel;
            int comove = 1;
            float kernel, vkernel;
            float zlower, zupper;
            float abs_zlower, abs_zupper;
            float radius2, radius;
            float zo2, zo, zi2, zi;
            float d, d2;
            int bin, bin_min, bin_max;
            double zcoord;

            FILE *binfile, *partfile;
            char binname[80], partname[80];

            double ion_weight;

            int ilo, ihi;
            float frh;
            float coldphasemassfrac, coldphasetemp;
            /*
        i = floor(spece_set.theta * 180.0 / PI + 0.001);
        k = floor(spece_set.redshift_begin * 100 + 0.0001);
        //k = floor((redshift_low-delta_redshift+0.01)*100+0.0001);
        j = floor(spece_set.redshift_end * 100 + 0.0001);
*/

            scp = spec_particles[i];

            for (irep[0] = -1; irep[0] <= 1; irep[0]++)
            {
                for (irep[1] = -1; irep[1] <= 1; irep[1]++)
                {
                    for (irep[2] = -1; irep[2] <= 1; irep[2]++)
                    {
                        if (irep[0] * BOXSIZE + scp.pos[0] > bound_max[0] && irep[0] * BOXSIZE + scp.pos[0] - 2 * scp.hsmooth > bound_max[0])
                            continue;
                        if (irep[0] * BOXSIZE + scp.pos[0] < bound_min[0] && irep[0] * BOXSIZE + scp.pos[0] + 2 * scp.hsmooth < bound_min[0])
                            continue;
                        if (irep[1] * BOXSIZE + scp.pos[1] > bound_max[1] && irep[1] * BOXSIZE + scp.pos[1] - 2 * scp.hsmooth > bound_max[1])
                            continue;
                        if (irep[1] * BOXSIZE + scp.pos[1] < bound_min[1] && irep[1] * BOXSIZE + scp.pos[1] + 2 * scp.hsmooth < bound_min[1])
                            continue;
                        if (irep[2] * BOXSIZE + scp.pos[2] > bound_max[2] && irep[2] * BOXSIZE + scp.pos[2] - 2 * scp.hsmooth > bound_max[2])
                            continue;
                        if (irep[2] * BOXSIZE + scp.pos[2] < bound_min[2] && irep[2] * BOXSIZE + scp.pos[2] + 2 * scp.hsmooth < bound_min[2])
                            continue;
                        for (k = 0; k < NDIM; k++)
                            part_pos[k] = scp.pos[k] + irep[k] * BOXSIZE;
                        los.RotateCoords(part_pos[0], part_pos[1], part_pos[2], spece_set.theta, spece_set.phi);
                        distnorm = 1. / (scp.hsmooth * scp.hsmooth);
                        radius2 = (pow(los.xrot - los.xrotline, 2) + pow(los.yrot - los.yrotline, 2)) * distnorm;
                        if (radius2 >= 4.0)
                            continue;
                        radius = sqrt(radius2);
                        zo2 = 4. - radius2;
                        zo = sqrt(zo2);
                        if (radius2 < 1.)
                        {
                            zi2 = 1. - radius2;
                            zi = sqrt(zi2);
                        }
                        else
                        {
                            zi2 = zi = 0.0;
                        }
                        zcoord = los.zrot; // hold zrot because we need to do velocity

                        los.RotateCoords(scp.vel[0], scp.vel[1], scp.vel[2], spece_set.theta, spece_set.phi);
                        vz = los.zrot;

                        counttest++;

                        coldphasemassfrac = 1;
                        coldphasetemp = scp.temp;
                        if (scp.sfr > 0)
                        {
                            //coldphasemassfrac = (1e+08-cp.temp)/1e+08; // WRONG TO ASSUME TEMPERATURE 1E+08 FOR HOT PHASE, IS A COMPLICATED FUNCTION.
                            coldphasemassfrac = 0.9; // MOTIVATED BY HONG ET AL.  TO GET AROUND COMPLICATED FUNCTION OF TEMPERATURE.
                            coldphasetemp = 1e+03;   /* Modified 7-28-11 to split SH03 two phase. 10^3 and 10^8 */
                        }

                        bin_min = binarysearch((zcoord - zo * scp.hsmooth - los.zbeginline), bin_coord, n);
                        bin_max = binarysearch((zcoord + zo * scp.hsmooth - los.zbeginline), bin_coord, n);

                        for (k = 0; k < spece_ionall.nions; k++)
                        {
                            ionfrac[k].interplate(scp.rho, scp.temp);

                            if (k == 0 || k == 7)
                            { // For HI and MgII and assuming i9

                                double rhocgs = ionfrac[k].interplate(scp.rho, coldphasetemp);
                                ilo = 0;
                                ihi = NINTERP - 1;
                                frh = ionfrac[k].fraction * rhocgs / MHYDR * scp.hsmooth * spece_set.spece_para.unit_Length * spece_cosmo.aex; // ???
                                while (ihi - ilo > 1)
                                {
                                    if (ktable.KernIntTable[(ilo + ihi) / 2][1] * frh < NHILIM)
                                        ihi = (ilo + ihi) / 2;
                                    else
                                        ilo = (ilo + ihi) / 2;
                                }
                                if (ilo > 0)
                                    ionfrac[k].fraction = coldphasemassfrac * (ionfrac[k].fraction * (ktable.KernIntTable[(ilo + ihi) / 2][0]) + 1 / (1 + pow((rhocgs / MHYDR * coldphasetemp) / P0BLITZ, ALPHA0BLITZ)) * (1 - ktable.KernIntTable[(ilo + ihi) / 2][0]));
                                //if( ktable.KernIntTable[(ilo+ihi)/2][0] < 0.999999 )
                                //fprintf(stdout,"INTKERNELNHLIMIT: % 5.3f %5.3f % 5.3f %6.2f %5d %5d %5.3e %7.5f % 7.5f % 5.3e % 5.3e % 5.3e %5.3e % 5.3e % 5.3e % 5.3e %5.3e\n",log10(rhocgs/MHYDR),log10(cp.temp),log10(ionfrac[k]),log10(frh),ilo,ihi,ktable.KernIntTable[(ilo+ihi)/2][1]*frh,ktable.KernIntTable[(ilo+ihi)/2][0],ktable.KernIntTable[(ilo+ihi)/2][1],log10(1-ktable.KernIntTable[(ilo+ihi)/2][1]),cp.hsmooth*unit_Length,cp.hsmooth,cp.sfr,coldphasemassfrac,1/(1+pow((rhocgs/MHYDR*coldphasetemp)/P0BLITZ,ALPHA0BLITZ)),(rhocgs/MHYDR*cp.temp),IonFrac(coldphasetemp,rhocgs,k));
                            }
                        }
                        if ((bin_min <= 0 && bin_max == 0) || (bin_min >= n - 1 && bin_max >= n - 1))
                            continue;
                        for (bin = bin_min; bin <= bin_max; bin++)
                        {
                            countscpall += bin_max - bin_min + 1;
                            if (bin_min >= bin_max)
                                fprintf(stderr, "(%d)ALERT!!! bin_min (%d) >= bin_max (%d)\n", countscp, bin_min, bin_max);

                            /*#ifndef SHORTSPEC
                        if (bin >= nzbins)
                            cosmopar(CosmicTime(bin_redshift[nzbins - 1]));
                        else if (bin < bin_min)
                            cosmopar(CosmicTime(bin_redshift[0]));
                        else
                            cosmopar(CosmicTime(bin_redshift[bin]));
#endif
*/
                            zlower = (bin_coord[bin] + los.zbeginline - zcoord) / scp.hsmooth;
                            zupper = (bin_coord[bin + 1] + los.zbeginline - zcoord) / scp.hsmooth;

                            //if(bin==0)fprintf(stdout,"BINMIN: irep= % 5.3f zcoord= %7.5f zhsmooth= %7.5f zbeginline= %7.5f bin_coord[0]= %5.3e bin_min= %5d bin_max= %5d nzbins= %5d %5d bin= %5d\n",irep[2],zcoord, zo*cp.hsmooth, zbeginline,bin_coord[0],bin_min,bin_max,bin_max-bin_min,nzbins,bin);
                            //if(bin==nzbins-1)fprintf(stdout,"BINMAX: irep= % 5.3f zcoord= %7.5f zhsmooth= %7.5f zbeginline= %7.5f bin_coord[0]= %5.3e bin_min= %5d bin_max= %5d nzbins= %5d %5d bin= %5d\n",irep[2],zcoord, zo*cp.hsmooth, zbeginline,bin_coord[0],bin_min,bin_max,bin_max-bin_min,nzbins,bin);

                            kernel = 0.0;
                            vkernel = 0.0;
                            abs_zlower = fabs(zlower);
                            if (abs_zlower >= zo)
                            {
                                kernel += 1.3125 * radius2 * zo - 1.5 * zo +
                                          0.5 * zo2 * zo - (1.5 * radius2 + 0.09375 * radius2 * radius2) * log(zo + 2.0);
                                /*kernel at zo in region 2 */
                                if (comove == 1)
                                {
                                    vkernel -= zo2 + 0.75 * radius2 * zo2 +
                                               0.375 * zo2 * zo2 - 9.6;
                                }
                            }
                            else if (abs_zlower < zo && (abs_zlower > zi ||
                                                         zlower == zi))
                            {
                                d2 = radius2 + zlower * zlower;
                                d = sqrt(d2);
                                kernel -= copysign(1.0, zlower) * ((2.0 +
                                                                    1.5 * radius2) *
                                                                       abs_zlower +
                                                                   0.5 * abs_zlower * abs_zlower * abs_zlower -
                                                                   0.0625 * abs_zlower * d * d2 - (1.5 + 0.09375 * radius2) * abs_zlower * d - (1.5 * radius2 + 0.09375 * radius2 * radius2) * log(abs_zlower + d));
                                /* sign(zlower)*kernel at abs_zlower in region 2 */
                                if (comove == 1)
                                {
                                    vkernel -= zlower * zlower * (1. + 0.75 * radius2) + 0.375 * zlower * zlower * zlower * zlower - d * d2 -
                                               0.05 * d2 * d2 * d;
                                }
                            }
                            else
                            {
                                d2 = radius2 + zlower * zlower;
                                d = sqrt(d2);
                                kernel -= copysign(1.0, zlower) * ((1.0 -
                                                                    1.5 * radius2) *
                                                                       abs_zlower -
                                                                   0.5 *
                                                                       abs_zlower * abs_zlower * abs_zlower +
                                                                   0.1875 * abs_zlower * d2 * d + 0.28125 * radius2 * abs_zlower * d + 0.28125 * radius2 * radius2 * log(abs_zlower + d));
                                /* sign(zlower)* kernel at abs_zlower in region 1 */
                                if (comove == 1)
                                {
                                    vkernel -= (0.5 - 0.75 * radius2) * zlower *
                                                   zlower -
                                               0.375 * zlower * zlower *
                                                   zlower * zlower +
                                               0.15 * d2 * d2 * d;
                                }
                            }
                            abs_zupper = fabs(zupper);
                            if (zlower * zupper <= 0.0)
                            { /* bin stradles zero */
                                if (radius2 < 1.)
                                {
                                    kernel -= 0.5625 * radius2 * radius2 *
                                              log(radius);
                                    /* 2.0 * kernel at zero in region 1 */
                                }
                                else
                                {
                                    kernel -= -(3.0 * radius2 + 0.1875 * radius2 * radius2) *
                                              log(radius);
                                    /* 2.0 *kernel at zero in region 2 */
                                }
                            }
                            if (radius2 < 1.)
                            {
                                if (zlower * zupper < 0.0 && abs_zupper > zi &&
                                    abs_zlower > zi)
                                {
                                    kernel += 2.375 * zi - 2.4375 * radius2 * zi -
                                              zi2 * zi + 0.5625 * radius2 * radius2 * log(zi + 1);
                                    /* 2.0 * kernel at zi in region 1 */
                                    kernel -= 0.875 * zi + 2.8125 * radius2 * zi +
                                              zi2 * zi - (3.0 * radius2 + 0.1875 * radius2 * radius2) * log(zi + 1);
                                    /* 2.0 * kernel at zi in region 2 */
                                }
                                else
                                {
                                    if ((zlower < -zi && zupper > -zi) ||
                                        (zlower < zi && zupper > zi))
                                    {
                                        kernel += 1.1875 * zi - 1.21875 * radius2 * zi - 0.5 * zi2 * zi + 0.28125 * radius2 * radius2 * log(zi + 1);
                                        /* kernel at zi in region 1 */
                                        if (comove == 1)
                                        {
                                            vkernel += copysign(1.0, zlower) * ((0.5 - 0.75 * radius2) *
                                                                                    zi2 -
                                                                                0.375 * zi2 * zi2 + 0.15);
                                        }
                                        kernel -= 0.4375 * zi + 1.40625 * radius2 * zi + 0.5 * zi2 * zi - (1.5 * radius2 + 0.09375 * radius2 * radius2) * log(zi + 1);
                                        /* kernel at zi in region 2 */
                                        if (comove == 1)
                                        {
                                            vkernel -= copysign(1.0, zlower) * ((1.0 + 0.75 * radius2) *
                                                                                    zi2 +
                                                                                0.375 * zi2 * zi2 - 1.05);
                                        }
                                    }
                                }
                            }

                            if (abs_zupper >= zo)
                            {
                                kernel += 1.3125 * radius2 * zo - 1.5 * zo +
                                          0.5 * zo2 * zo - (1.5 * radius2 + 0.09375 * radius2 * radius2) * log(zo + 2.0);
                                /* kernel at zo in region 2 */
                                if (comove == 1)
                                {
                                    vkernel += zo2 + 0.75 * radius2 * zo2 +
                                               0.375 * zo2 * zo2 - 9.6;
                                }
                            }
                            else if (abs_zupper < zo && abs_zupper >= zi)
                            {
                                d2 = radius2 + zupper * zupper;
                                d = sqrt(d2);
                                kernel += copysign(1.0, zupper) * ((2.0 + 1.5 * radius2) *
                                                                       abs_zupper +
                                                                   0.5 * abs_zupper *
                                                                       abs_zupper * abs_zupper -
                                                                   0.0625 *
                                                                       abs_zupper * d * d2 -
                                                                   (1.5 + 0.09375 *
                                                                              radius2) *
                                                                       abs_zupper * d -
                                                                   (1.5 * radius2 +
                                                                    0.09375 * radius2 * radius2) *
                                                                       log(abs_zupper + d));
                                /* sign(zupper)*kernel at abs_zupper in region 2 */
                                if (comove == 1)
                                {
                                    vkernel += zupper * zupper * (1. + 0.75 * radius2) + 0.375 * zupper * zupper * zupper * zupper - d * d2 -
                                               0.05 * d2 * d2 * d;
                                }
                            }
                            else
                            {
                                d2 = radius2 + zupper * zupper;
                                d = sqrt(d2);
                                kernel += copysign(1.0, zupper) * ((1.0 -
                                                                    1.5 * radius2) *
                                                                       abs_zupper -
                                                                   0.5 *
                                                                       abs_zupper * abs_zupper * abs_zupper +
                                                                   0.1875 * abs_zupper * d2 * d + 0.28125 * radius2 * abs_zupper * d + 0.28125 * radius2 * radius2 * log(abs_zupper + d));
                                /* sign(zupper)* kernel at abs_zupper in region 1 */
                                if (comove == 1)
                                {
                                    vkernel += (0.5 - 0.75 * radius2) * zupper *
                                                   zupper -
                                               0.375 * zupper * zupper *
                                                   zupper * zupper +
                                               0.15 * d2 * d2 * d;
                                }
                            }
                            kernel *= distnorm / PI;
                            if (comove == 1)
                            {
                                vkernel /= PI * scp.hsmooth;
                            }

                            //unit_vel = 1;
                            //#ifndef OWLSFORMAT /* Not sure about this... */
                            unit_vel = spece_set.spece_para.unit_Velocity * spece_cosmo.aex / 1.e5;
                            //#endif
                            //if(kernel > -1.0){ /* It is an unexplained occurence why some kernel values are anomolously very negative... but they are rare so we will ignore them */
                            if (kernel > 0.0)
                            { /* Really, negative kernel values are bad.  They should be ignored at all costs.  7-9-11 */
                                if (countid != scp.ID)
                                {
                                    countid = scp.ID;
                                    countscp++;
                                    // cerr << countid << "kernelcount  " << kernelcount << " " << endl;
                                    //kernelcount = 0.0;
                                }
                                //kernelcount += kernel * scp.mass;//test
                                spece_ionall.ion_total.mass[bin + 0] += kernel * scp.mass;
                                //if(bin==nzbins-2 || bin==nzbins-1){
                                //fprintf(stdout,"LASTBIN: nzbins= %5d kernel= % 5.3e scp.mass= %5.3e mass= %5.3e tot= %5.3e zlower= % 5.3e zupper= % 5.3e radius2= % 5.3e radius= % 5.3e\n",bin,kernel,scp.mass,kernel*scp.mass,spece_ionall.ion_total.mass[bin+0],zlower,zupper,radius2,radius);
                                //}
                                //if(bin%3000==0)fprintf(stderr,"(%d)= %g %g %g\n",bin+0,spece_ionall.ion_total.vel[bin+0],kernel,scp.mass);
#ifdef ZEROVEL
                                spece_ionall.ion_total.vel[bin + 0] = 0;
#else
                                spece_ionall.ion_total.vel[bin + 0] += kernel * (scp.mass) * unit_vel * vz;
#endif
                                if (comove == 1)
                                {
#ifndef ZEROVEL
                                    spece_ionall.ion_total.vel[bin + 0] += vkernel * (scp.mass) * spece_cosmo.hubble * unit_vel;
#endif
                                }
                                spece_ionall.ion_total.temp[bin + 0] += kernel * (scp.mass) * (scp.temp);
                                spece_ionall.ion_total.rho[bin + 0] += kernel * (scp.mass) * (scp.rho * XH * spece_set.spece_para.unit_Density / pow(spece_cosmo.aex, 3.0));

                                for (int m = 0; m < NMETALS; m++)
                                {
                                    spece_ionall.ion_total.metals[m][bin + 0] += kernel * (scp.mass) * (scp.metals[m]);
                                    if (scp.metals[m] > 10)
                                    {
                                        printf("BAD scp METAL!!\n");
                                        exit(-1);
                                    }
                                }
#ifdef PHYSSPEC
                                spece_ionall.ion_total.sfr[bin + 0] += kernel * (scp.mass) * (scp.sfr);
                                if (scp.mgal > 0 && scp.dgal > 0 && scp.ageaway > 0)
                                {
                                    spece_ionall.ion_total.wtmass[bin + 0] += kernel * scp.mass;
                                    spece_ionall.ion_total.age[bin + 0] += kernel * scp.mass * (scp.ageaway);
                                    spece_ionall.ion_total.dgal[bin + 0] += kernel * scp.mass * (log10(scp.dgal));
                                    spece_ionall.ion_total.mgal[bin + 0] += kernel * scp.mass * (scp.mgal);
                                    spece_ionall.ion_total.vlaunch[bin + 0] += kernel * scp.mass * (scp.vlaunch);
                                    spece_ionall.ion_total.nrec[bin + 0] += kernel * scp.mass * (scp.nrec);
                                }

#endif

                                for (k = 0; k < spece_ionall.nions; k++)
                                {
                                    
                                    if (spece_ionall.ions[k].Zcolumn == -1)
                                    {
                                        ion_weight = ionfrac[k].fraction;
                                    }
                                    else
                                    {
                                        if (spece_ionall.ions[k].Zcolumn < -1)
                                        {
                                            //ion_weight = ionfrac[k]*scp.metals[3]/0.001267*spece_ionall.ion[k].fraction*pow(10,spece_ionall.ion[k].alpha);
                                            ion_weight = ionfrac[k].fraction * scp.metals[1] / 0.009618 * spece_ionall.ions[k].fraction * pow(10, spece_ionall.ions[k].alpha); /* 2-11-10 */
                                        }
                                        else
                                        {
                                            ion_weight = ionfrac[k].fraction * scp.metals[spece_ionall.ions[k].Zcolumn];
                                        }
                                    }

                                    spece_ionall.ions[k].mass[bin + 0] += kernel * ion_weight * scp.mass;
                                    if (counttest <= 10)
                                    {
                                        //cerr<<k<<" mass: "<< scp.mass<<" T: "<<scp.temp<<" v: "<<ion_weight<<endl;;
                                    }
#ifdef ZEROVEL
                                    spece_ionall.ions[k].vel[bin + 0] = 0;
#else
                                    spece_ionall.ions[k].vel[bin + 0] += kernel * ion_weight * (scp.mass) * vz * unit_vel;
#endif
                                    if (comove == 1)
                                    {
                                        spece_ionall.ions[k].vel[bin + 0] += vkernel * ion_weight * (scp.mass) * spece_cosmo.hubble * unit_vel;
                                    }

                                    spece_ionall.ions[k].temp[bin + 0] += kernel * ion_weight * scp.mass * scp.temp;
                                    spece_ionall.ions[k].rho[bin + 0] += kernel * ion_weight * scp.mass * (scp.rho * XH * spece_set.spece_para.unit_Density / pow(spece_cosmo.aex, 3.0));
                                    for (int m = 0; m < NMETALS; m++)
                                        spece_ionall.ions[k].metals[m][bin + 0] += kernel * ion_weight * (scp.mass) * (scp.metals[m]);
#ifdef PHYSSPEC
                                    spece_ionall.ions[k].sfr[bin + 0] += kernel * ion_weight * (scp.mass) * (scp.sfr);
                                    if (scp.mgal > 0 && scp.dgal > 0 && scp.ageaway > 0)
                                    {
                                        spece_ionall.ions[k].wtmass[bin + 0] += kernel * ion_weight * scp.mass;
                                        spece_ionall.ions[k].age[bin + 0] += kernel * ion_weight * scp.mass * scp.ageaway;
                                        spece_ionall.ions[k].dgal[bin + 0] += kernel * ion_weight * scp.mass * (log10(scp.dgal));
                                        spece_ionall.ions[k].mgal[bin + 0] += kernel * ion_weight * scp.mass * scp.mgal;
                                        spece_ionall.ions[k].vlaunch[bin + 0] += kernel * ion_weight * scp.mass * (scp.vlaunch);
                                        spece_ionall.ions[k].nrec[bin + 0] += kernel * ion_weight * scp.mass * (scp.nrec);
                                    }
#endif
                                }
                            }
                            else
                            {

                                {

                                    //cerr << "    kernl wrrong" << kernel << endl;
                                }

                                //exit(-1);
                            }

                            //if(bin==0) printf("MIN r = %9.7e k = % 9.7e m = % 9.7e r = % 9.7e b_max = %5d d = %5d scp.rho = %5.3e z = %7.5f aex = %7.5e\n",radius,kernel,spece_ionall.ion_total.mass[bin+0],spece_ionall.ion_total.rho[bin+0],bin_max,bin_max-bin_min,rhocgs,bin_redshift[bin],XH*unit_Density/(aex*aex*aex));
                            //if(bin==nzbins-1) printf("MAX r = %9.7e k = % 9.7e m = % 9.7e r = %9.7e b_min = %5d d = %5d scp.rho = %5.3e z = %7.5f aex = %7.5e\n",radius,kernel,spece_ionall.ion_total.mass[bin+0],spece_ionall.ion_total.rho[bin+0],bin_min,bin_max-bin_min,rhocgs,bin_redshift[bin],XH*unit_Density/(aex*aex*aex));
                        }
                    }
                }
            }
        }
    }
    cerr << endl
         << "special:  " << countscp << " " << kernelcount << "  " << endl;
    free(bin_coord);
    int counteffective = 0;
#ifdef TEST
    char filenameeffective[100];
    sprintf(filenameeffective, "./result/%s/galaxy%s/efc_%f_%f_%f_%f.txt", spece_set.TIME, spece_set.galaxy, spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
    ofstream efc(filenameeffective);
    sprintf(filenameeffective, "./result/%s/galaxy%s/efele_%f_%f_%f_%f.txt", spece_set.TIME, spece_set.galaxy, spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
    ofstream efc1(filenameeffective);
#endif
    for (int i = 0; i < los.nzbins; i++)
    {

        if (spece_ionall.ion_total.mass[i] != 0.0)
        {
            counteffective++;
            spece_ionall.ion_total.vel[i] /= spece_ionall.ion_total.mass[i];
            spece_ionall.ion_total.temp[i] /= spece_ionall.ion_total.mass[i];
            spece_ionall.ion_total.rho[i] /= spece_ionall.ion_total.mass[i];
            for (int m = 0; m < NMETALS; m++)
            {
                spece_ionall.ion_total.metals[m][i] /= spece_ionall.ion_total.mass[i];
                if (spece_ionall.ion_total.metals[m][i] < 0)
                    spece_ionall.ion_total.metals[m][i] = 0.0e+00;
            }
        }
        spece_ionall.ion_total.redshift[i] = los.bin[i].redshift;

        spece_ionall.ion_total.binsize[i] = los.bin[i].size;
        spece_ionall.ion_total.bincoord[i] = los.bin[i].relative_coord; // + hold_coord;
        spece_ionall.x[i] = los.bin[i].x;
        spece_ionall.y[i] = los.bin[i].y;
        spece_ionall.z[i] = los.bin[i].z;
        //#ifdef PARTIONFRAC
        for (int k = 0; k < spece_ionall.nions; k++)
        {
            if (spece_ionall.ions[k].mass[i] > 0.0)
            {

                spece_ionall.ions[k].vel[i] /= spece_ionall.ions[k].mass[i];
                spece_ionall.ions[k].temp[i] /= spece_ionall.ions[k].mass[i];
                spece_ionall.ions[k].rho[i] /= spece_ionall.ions[k].mass[i];
                for (int m = 0; m < NMETALS; m++)
                {
                    spece_ionall.ions[k].metals[m][i] /= spece_ionall.ions[k].mass[i];
                    if (spece_ionall.ions[k].metals[m][i] < 0)
                        spece_ionall.ions[k].metals[m][i] = 0.0e+00;
                }
            }
        }
#ifdef TEST
        double unit_mass = spece_set.spece_para.unit_Mass;
        double unit_binsize = spece_cosmo.aex * spece_set.spece_para.unit_Length; //(spece_set.spece_para.unit_Length / (double)los.nzbins);
        double unittemp = unit_mass / DD(unit_binsize) / MHYDR;
        efc
            << spece_ionall.ion_total.mass[i] * 0.76 * unittemp / spece_ionall.ions[0].atomwt << " "
            << spece_ionall.ion_total.mass[i] * spece_ionall.ion_total.metals[1][i] * 1.0 / 0.009618 * spece_ionall.ions[7].fraction * pow(10, spece_ionall.ions[7].alpha) * unittemp / spece_ionall.ions[7].atomwt << " "
            << spece_ionall.ion_total.mass[i] * spece_ionall.ion_total.metals[0][i] * unittemp / spece_ionall.ions[2].atomwt << endl;

        efc1
            << spece_ionall.ions[0].fraction * spece_ionall.ions[0].mass[i] * unittemp / spece_ionall.ions[0].atomwt << " "
            << spece_ionall.ions[1].fraction * spece_ionall.ions[1].mass[i] * unittemp / spece_ionall.ions[1].atomwt << " "
            << spece_ionall.ions[2].mass[i] * unittemp / spece_ionall.ions[2].atomwt << " "
            << spece_ionall.ions[3].mass[i] * unittemp / spece_ionall.ions[3].atomwt << " "
            << spece_ionall.ions[4].mass[i] * unittemp / spece_ionall.ions[4].atomwt << " "
            << spece_ionall.ions[5].mass[i] * unittemp / spece_ionall.ions[5].atomwt << " "
            << spece_ionall.ions[6].mass[i] * unittemp / spece_ionall.ions[6].atomwt << " "
            << spece_ionall.ions[7].mass[i] * unittemp / spece_ionall.ions[7].atomwt << " "
            << spece_ionall.ions[8].mass[i] * unittemp / spece_ionall.ions[8].atomwt << endl;
#endif
    }
#ifdef TEST
    efc.close();
    efc1.close();
#endif
    spf.close();
}

#endif