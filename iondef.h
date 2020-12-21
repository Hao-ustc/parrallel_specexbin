#ifndef IONDEF_H
#define IONDEF_H
#include "include_all.h"
using namespace std;
#define NMETALS 4
/* Ion lookup table definitions */
#define MAXIONS 35
#define NHPTS 240 /* Number of n_h gridpoints in lookup table */
#define TPTS 140  /* Number of T gridpoints in lookup table */
#define GALPTS 11
/* Table limits */
#define NHLOW -9.0
#define DELTANH 0.05
#define TLOW 2.5
#define DELTAT 0.05
#define GALLOW 6.0
#define DELTAGAL 0.5

#define NSUBBINMIN 4
#define NSUBBINMAX 1000000
#define NBSMOOTH 5.0
#define SIGN(x) ((x) < 0.0 ? -1.0 : 1.0)
//#define clight 2.99792458e10TEST_ION
#define ZRES_ION_TABLE 0.05

//#define TEST_ION

class Ion_special
{

public:
  char name[10];
  double *mass;
  double *vel;
  double *temp;
  double *rho;
  double *metals[NMETALS];
#ifdef PHYSSPEC
  double *sfr;
  double *wtmass;
  double *mgal;
  double *dgal;
  double *age;
  double *nrec;
  double *vlaunch;
#endif
  double *redshift;
  double *binsize;
  double *bincoord;
  double *vbins;
  double *tbins;
  double *rhobins;
  double *Zbins;
#ifdef PHYSSPEC
  double *sfrbins;
  double *mgalbins;
  double *dgalbins;
  double *agebins;
  double *nrecbins;
  double *vlaunchbins;
#endif
  float lambda, fraction, Xsec, atomwt, bsys, alpha;
  int Zcolumn;
};

class Ion_all
{
private:
  Ion_special ion;
  void Load_ion(LOS &los, Setting &spece_set);
  void Load_ion_total(LOS &los);

public:
  vector<Ion_special> ions;
  Ion_special ion_total;
  int nions;
  double *redshift;
  double *x;
  double *y;
  double *z;
  double *xbins;
  double *ybins;
  double *zbins;
  double *gal_field;

  void Load(LOS &los, Setting &spece_set); //
  int Tau(LOS &los, Setting &spece_set);
  void Freeions();
};

void Ion_all::Load_ion(LOS &los, Setting &spece_set)
{
  FILE *specfile;
  char line[80];
  int nzbins = los.nzbins;
  int nvbins = los.nvbins;

  if ((specfile = fopen(spece_set.spec_ion_filename, "r")) == NULL)
  {
    fprintf(stderr, "cannot find specion file anywhere\n");
    exit(-1);
  }

  for (int i = 0; i < MAXIONS; i++)
  {

    ions.push_back(ion);
    ions[i].mass = (double *)malloc(nzbins * sizeof(double));
    ions[i].vel = (double *)malloc(nzbins * sizeof(double));
    ions[i].temp = (double *)malloc(nzbins * sizeof(double));
    ions[i].rho = (double *)malloc(nzbins * sizeof(double));
    for (int j = 0; j < NMETALS; j++)
      ions[i].metals[j] = (double *)malloc(nzbins * sizeof(double));

    ions[i].vbins = (double *)malloc(nvbins * sizeof(double));
    ions[i].tbins = (double *)malloc(nvbins * sizeof(double));
    ions[i].rhobins = (double *)malloc(nvbins * sizeof(double));
    ions[i].Zbins = (double *)malloc(nvbins * sizeof(double));
  }
  int i = 0;
  while (fgets(line, 80, specfile) != NULL)
  {
    if (strstr(line, "#") != NULL)
      continue;
    if (i >= MAXIONS)
      break;
    sscanf(line, "%10s %g %g %g %g %d %g", ions[i].name, &ions[i].lambda, &ions[i].Xsec, &ions[i].atomwt, &ions[i].fraction, &ions[i].Zcolumn, &ions[i].alpha);
    i++;
  }
  nions = i;
  fclose(specfile);

  fprintf(stderr, "Processing %d ions from specions.dat:\n", nions);
  for (i = 0; i < nions; i++)
  {
    for (int j = 0; j < nzbins; j++)
    {
      ions[i].mass[j] = ions[i].vel[j] = ions[i].temp[j] = ions[i].rho[j] = 0.0;
      for (int m = 0; m < NMETALS; m++)
        ions[i].metals[m][j] = 0;
    }
    ions[i].bsys = sqrt(2. * KBOLTZ / (MHYDR * ions[i].atomwt)) / 1.e5;
    ions[i].Xsec *= 2.648e-2 * ions[i].lambda * 1.e-13;
    fprintf(stderr, "%5d %10s %12.6g %12.6g %10.5g %10.5g %10.5g % 3.1f % 2d\n", i, ions[i].name, ions[i].lambda, ions[i].Xsec, ions[i].atomwt, ions[i].fraction, ions[i].bsys * sqrt(1.e4), ions[i].alpha, ions[i].Zcolumn);
  }
}
void Ion_all::Load_ion_total(LOS &los)
{
  int nzbins = los.nzbins;
  int nvbins = los.nvbins;

  ion_total.mass = (double *)malloc(nzbins * sizeof(double));
  ion_total.vel = (double *)malloc(nzbins * sizeof(double));
  ion_total.temp = (double *)malloc(nzbins * sizeof(double));
  ion_total.rho = (double *)malloc(nzbins * sizeof(double));
  for (int j = 0; j < NMETALS; j++)
    ion_total.metals[j] = (double *)malloc(nzbins * sizeof(double));
  ion_total.redshift = (double *)malloc(nzbins * sizeof(double));
  ion_total.binsize = (double *)malloc(nzbins * sizeof(double));
  ion_total.bincoord = (double *)malloc(nzbins * sizeof(double));

  ion_total.vbins = (double *)malloc(nvbins * sizeof(double));
  ion_total.tbins = (double *)malloc(nvbins * sizeof(double));
  ion_total.rhobins = (double *)malloc(nvbins * sizeof(double));
  ion_total.Zbins = (double *)malloc(nvbins * sizeof(double));

  redshift = (double *)malloc(nvbins * sizeof(double));
  x = (double *)malloc(nzbins * sizeof(double));
  y = (double *)malloc(nzbins * sizeof(double));
  z = (double *)malloc(nzbins * sizeof(double));
  xbins = (double *)malloc(nvbins * sizeof(double));
  ybins = (double *)malloc(nvbins * sizeof(double));
  zbins = (double *)malloc(nvbins * sizeof(double));
  gal_field = (double *)malloc(nvbins * sizeof(float));
  for (int j = 0; j < nzbins; j++)
  {
    ion_total.mass[j] = ion_total.vel[j] = ion_total.temp[j] = ion_total.rho[j] = 0;
    for (int m = 0; m < NMETALS; m++)
      ion_total.metals[m][j] = 0;
  }
}
void Ion_all::Load(LOS &los, Setting &spece_set)
{
  Load_ion(los, spece_set);
  Load_ion_total(los);
}

int Ion_all::Tau(LOS &los, Setting &spece_set)
{
  int counttest = 0;
  /* Calculate optical depth for given ion along line of sight */

  int ionid;
  int i, j, k, l;
  int Zcol, imet;
  //float z;
  double b;
  double mass_interp, v_interp;
  double irepz;
  double t_interp, rho_interp, Z_interp;
#ifdef PHYSSPEC
  float mgal_interp, dgal_interp, age_interp, nrec_interp, vlaunch_interp, sfr_interp;
#endif

  int bin, bin_min, bin_max, bin_cen;
  double vlower, vupper;
  float abs_vlower, abs_vupper;
  float dvcol;
  double unit_col;

  i = 0;
  k = 0;
  j = 0;
  l = 0;
  Zcol = 0.0;
  imet = 0.0;
  b = 0.0;
  mass_interp = 0.0;
  v_interp = 0.0;
  irepz = 0;
  t_interp = rho_interp = Z_interp = 0.0;

  double voffset, vmin;

  float *norm;
  float *norm1; //test hao
  int *norm_field;

  double *vbin_size, *vbin_coord, *vbin_zsize, *vbin_zcoord;
  double zstep, vcoord, vstep, zcoord;

  float hubble_expansion;
  double hubble = spece_cosmo.hubble;
  double aex = spece_cosmo.aex;
  double unit_Velocity = spece_set.spece_para.unit_Velocity;
  double unit_Mass = spece_set.spece_para.unit_Mass;
  double unit_Lenth = spece_set.spece_para.unit_Length;

  Ion_special I;
  double x_interp, y_interp, z_interp, weight_sub;
  int nsubbinvar;

  /* we are out of the loop so reassign nzbins and nvbins 
	nzbins = nzloopbins;
	nvbins = nvloopbins;
  */

  float floatredshift = spece_set.redshift_center;
  float t = CosmicTime(floatredshift, spece_set);
  cosmopar(t, spece_set);

  //redshift_track = IonTotal.redshift[0];
  int nzbins = los.nzbins;
  int nvbins = los.nvbins;

  vbin_size = (double *)malloc(nzbins * sizeof(double));
  vbin_coord = (double *)malloc(nzbins * sizeof(double));
  vbin_zsize = (double *)malloc(nzbins * sizeof(double));
  vbin_zcoord = (double *)malloc(nzbins * sizeof(double));
  norm = (float *)malloc(nvbins * sizeof(float));
  norm1 = (float *)malloc(nvbins * sizeof(float)); //test hao
  norm_field = (int *)malloc(nvbins * sizeof(int));

  i = 0;
  vcoord = 0;
  zcoord = 0;
  //	fprintf(stderr, "TAUBEGIN: redshift_track= %g  IonTotal.redshift[nzbins-1]= %g nzbins= %d\n", redshift_track, IonTotal.redshift[nzbins - 1], nzbins);

  zstep = VRES / (BOXSIZE * spece_cosmo.hubble * spece_set.spece_para.unit_Velocity / CLIGHT); //const

  vstep = zstep * spece_cosmo.aex * spece_cosmo.hubble * spece_set.spece_para.unit_Velocity; //const
                                                                                             //if (vstep == 0)
  cerr << "vstep  " << vstep << endl;
  double redshift_track = ion_total.redshift[0];
  while (redshift_track >= ion_total.redshift[nzbins - 1])
  {
    redshift_track -= VRES;
    zcoord += zstep;
    vcoord += vstep;

    i++;

    vbin_size[i - 1] = vstep / 1.e5;   //const
    vbin_coord[i - 1] = vcoord / 1.e5; //i
    vbin_zsize[i - 1] = zstep;         //const
    vbin_zcoord[i - 1] = zcoord;       //i
  }
  cerr << "v_coord.size=" << i << endl;
  //column setting
  voffset = 0;
  vmin = 0.0;
//0706  start
#ifdef TEST_ION
  char filenameeffective[200];
  sprintf(filenameeffective, "./result/%s/galaxy%s/efc2_%f_%f_%f_%f.txt", spece_set.TIME, spece_set.galaxy, spece_set.redshift_center, spece_set.xspec, spece_set.yspec, spece_set.zspec);
  ofstream efc2(filenameeffective);
#endif
  //0706 end
  for (ionid = -1; ionid < nions; ionid++)
  {
    counttest = 0;
    if (ionid == -1)
      fprintf(stderr, "Computing tau for ion ");
    fprintf(stderr, "%d ", ionid);
    if (ionid == nions - 1)
      //    fprintf(stderr, "\n");
      fflush(stderr);
    if (ionid == -1)
    { /* Outputing in physical space */
      I = ion_total;
      I.atomwt = ions[0].atomwt;
      I.fraction = 0.0122;
      I.Zcolumn = -1;
      I.alpha = 0;
      I.bsys = 1e-10;
      I.Xsec = ions[0].Xsec;
    }
    else
    {
      I = ions[ionid];
    }
    for (i = 0; i < nvbins; i++)
    {
      I.vbins[i] = I.tbins[i] = I.rhobins[i] = I.Zbins[i] = norm[i] = norm1[i] = 0.0;
    }
    Zcol = I.Zcolumn;
    hubble_expansion = 0;

    double redshift = ion_total.redshift[0];

    for (i = 0; i < nzbins; i++)
    {

      hubble_expansion += hubble * aex * unit_Velocity / 1.e5 * ion_total.binsize[i]; //binsize with hubble

      nsubbinvar = NSUBBINMIN;
      if (i > 0 && i < nzbins - 1)
      {
        nsubbinvar = abs((int)((I.vel[i + 1] - I.vel[i - 1]) / ((ion_total.redshift[i - 1] - ion_total.redshift[i + 1]) * CLIGHT / 1e+05 / (1 + ion_total.redshift[i]))));
      }
      //if(i%100==0 && ionid==1)fprintf(stderr,"nsubbinvar= %3d %5.3e %5.3e",nsubbinvar,(I.vel[i+1]-I.vel[i-1]),(ion_total.redshift[i-1]-ion_total.redshift[i+1])*clight/1e+05/(1+ion_total.redshift[i]));
      if (nsubbinvar < NSUBBINMIN)
        nsubbinvar = NSUBBINMIN;
      if (nsubbinvar > NSUBBINMAX)
        nsubbinvar = NSUBBINMAX;
      //if(i%100==0 && ionid==1)fprintf(stderr," %3d\n",nsubbinvar);

      //if(i!=0)hubble_expansion += hubble*aex*unit_Velocity/1.e5*vbin_zsize[vi]*ZRES/VRES;
      if (ion_total.mass[i] == 0.0)
        continue;
      for (j = 0; j < nsubbinvar; j++)
      {
        counttest++;
        //z = (los.bin[i].relative_coord + ((double)(j)) / ((double)(nsubbinvar)) * los.bin[i].size);
        if (2 * j < nsubbinvar)
        {
          k = i - 1;
        }
        else
        {
          k = i;
        }
        if (k < 0)
          k = 0;
        if (k > nzbins - 2)
          k = nzbins - 2;

        //if(ionid==-1){
        //v_interp = 0;
        //}else{

        weight_sub = (i + ((double)(j)) / ((double)(nsubbinvar)) - (k + 0.5));
        if (I.mass[k + 1] < I.mass[k])
          weight_sub *= 2 * I.mass[k + 1] / (I.mass[k + 1] + I.mass[k]);
        else
          weight_sub = 1. - 2 * weight_sub * I.mass[k] / (I.mass[k + 1] + I.mass[k]);
        if (ionid == -1)
        {
          x_interp = (x[k + 1] - x[k]) * weight_sub + x[k];
          y_interp = (y[k + 1] - y[k]) * weight_sub + y[k];
          z_interp = (z[k + 1] - z[k]) * weight_sub + z[k];
          v_interp = hubble_expansion + hubble * aex * unit_Velocity / 1.e5 * ion_total.binsize[k] * ((double)(j)) / ((double)(nsubbinvar)) + voffset;
          bin_cen = binarysearch(v_interp - vmin, vbin_coord, nvbins);
          bin_cen = binarysearch(v_interp - vmin + vbin_size[bin_cen] / 2, vbin_coord, nvbins);
          if (bin_cen < 0 || bin_cen >= nvbins)
            continue;
          xbins[bin_cen] = x_interp;
          ybins[bin_cen] = y_interp;
          zbins[bin_cen] = z_interp;
        }
        mass_interp = (I.mass[k + 1] - I.mass[k]) * weight_sub + I.mass[k];
        //if(ionid>=0 && mass_interp==0) continue;
        v_interp = (I.vel[k + 1] - I.vel[k]) * weight_sub + I.vel[k]; /* No ionization weighting here: NOT TRUE NOW THERE IS */
                                                                      //}
        t_interp = (I.temp[k + 1] - I.temp[k]) * weight_sub + I.temp[k];
        rho_interp = (I.rho[k + 1] - I.rho[k]) * weight_sub + I.rho[k];

        //if( ionid==2 && (k==5062||k==5061)) fprintf(stderr,"%d %g %g %g %g %g %g\n",k,log10(I.temp[k]),log10(I.temp[k+1]),I.mass[k],I.mass[k+1],weight_sub,log10(t_interp));

        if (Zcol == -1)
        {
          Z_interp = 0;
          for (imet = 0; imet < NMETALS; imet++)
            Z_interp += (I.metals[imet][k + 1] - I.metals[imet][k]) * (i + ((double)(j)) / ((double)(nsubbinvar)) - (k + 0.5)) + I.metals[imet][k]; // Sum to total metals
          Z_interp *= 1.28;
        }
        else
        {
          if (Zcol < -1)
          {
            Z_interp = (I.metals[1][k + 1] - I.metals[1][k]) * (i + ((double)(j)) / ((double)(nsubbinvar)) - (k + 0.5)) + I.metals[1][k];
            //if(Z_interp>0)fprintf(stderr,"IONCALC %d: %g %g ",ionid,I.metals[3][k],Z_interp);
            //Z_interp *= I.fraction/0.001267*pow(10,I.alpha);
            Z_interp *= I.fraction / 0.009618 * pow(10, I.alpha); /* 2-11-10 */
                                                                  //if(Z_interp>0)fprintf(stderr," %g %g \n",I.fraction/0.001267*pow(10,I.alpha),Z_interp);
          }
          else
          {
            Z_interp = (I.metals[Zcol][k + 1] - I.metals[Zcol][k]) * (i + ((double)(j)) / ((double)(nsubbinvar)) - (k + 0.5)) + I.metals[Zcol][k]; // Fraction by mass of species
          }
        }

        b = I.bsys * sqrt(t_interp);

        v_interp += hubble_expansion + hubble * aex * unit_Velocity / 1.e5 * ion_total.binsize[k] * ((double)(j)) / ((double)(nsubbinvar)) + voffset;

        l = 0;

        irepz = 0;

        for (irepz = -1; irepz <= 1; irepz++)
        {

          bin_min = binarysearch((v_interp + irepz * vstep / 1.e5 * nvbins - NBSMOOTH * b - vmin), vbin_coord, nvbins);
          bin_max = binarysearch((v_interp + irepz * vstep / 1.e5 * nvbins + NBSMOOTH * b - vmin), vbin_coord, nvbins);

          if (counttest <= 10)
          {
            if (ionid != -1)
            {
              cerr << "0v_interp: " << I.vel[k] << " " << hubble_expansion << " ";
              cerr << ionid << " bin   " << bin_min << " " << bin_max << endl;
            }
          }

          if (bin_min < 0)
            bin_min = 0;
          if (bin_max >= nvbins)
            bin_max = nvbins - 1;

          //printf("ionid = %d i = %d vi = %5.3e bin_min = %d bin_max = %d nvbins = %d nzbins = %d\n",ionid,i,z,bin_min,bin_max,nvbins,nzbins);

          if ((bin_min == 0 && bin_max == 0) || (bin_min >= nvbins - 1 && bin_max >= nvbins - 1))
          {

            continue;
          }

          //if((bin_min == bin_max) && (irepz<-0.5 || irepz>0.5)) continue;

          for (bin = bin_min; bin <= bin_max; bin++)
          {
            /*if (ionid != -1)
            {
              cerr << ionid << " bin   " << bin_min << " " << bin_max << endl;
            }
            */
            //if(bin == bin_min){
            //vlower = -NBSMOOTH*b ; /* this appears to be screwing things up! 7-9-11 */
            //}
            //else{
            vlower = vbin_coord[bin] + vmin - v_interp - irepz * vstep / 1.e5 * nvbins;
            //}
            //if(bin == bin_max){
            //vupper = NBSMOOTH*b ;  /* along with this! 7-9-11 */
            //}
            //else{
            vupper = vbin_coord[bin + 1] + vmin - v_interp - irepz * vstep / 1.e5 * nvbins;
            //}

            vlower /= b;
            vupper /= b;
            abs_vlower = fabs(vlower);
            abs_vupper = fabs(vupper);

            if (vupper * vlower < 0)
            {
              dvcol = 1. / (double)nsubbinvar * 0.5 * (erf(abs_vlower) + erf(vupper));
            }
            else
            {
              if (abs_vlower < abs_vupper)
              {
                dvcol = 1. / (double)nsubbinvar * 0.5 * (erf(abs_vupper) - erf(abs_vlower));
              }
              else
              {
                dvcol = 1. / (double)nsubbinvar * 0.5 * (erf(abs_vlower) - erf(abs_vupper));
              }
            }
            if (mass_interp > 0)
            {

              /*
              if (ionid != -1)
              {
                cerr << ionid << " mass_interp   " << endl;
              }
              */
              if (I.Zcolumn == -1)
              {
                I.vbins[bin] += dvcol * mass_interp * I.fraction;
                norm1[bin] += mass_interp * I.fraction; //test hao
                if (I.vbins[bin] == 0)
                  cerr << "I.vbins[i] " << endl;
                if (dvcol == 0)
                  cerr << "col;" << endl;
                if (mass_interp == 0)
                  cerr << "mass_interp;" << endl;
              }
              else
              {
                I.vbins[bin] += dvcol * mass_interp;
                norm1[bin] += mass_interp; //test hao
                if (I.vbins[bin] == 0)
                  cerr << "I.vbins[i] " << endl;
                if (dvcol == 0)
                  cerr << "col;" << endl;
                if (mass_interp == 0)
                  cerr << "mass_interp;" << endl;
              }

              I.rhobins[bin] += dvcol * rho_interp * mass_interp;
              /*if (I.rhobins[bin] == 0 && rho_interp == 0)
              {
                cerr << "weight_sub : " << weight_sub << endl;
              }
              */
              I.tbins[bin] += dvcol * t_interp * mass_interp;
              I.Zbins[bin] += dvcol * Z_interp * mass_interp;

              norm[bin] += dvcol * mass_interp;

              /*if (isnan(norm[bin]) && isnan(dvcol))
              {
                counttest++;
                if (counttest <= 10)
                {
                  cerr << "v_interp : " << v_interp << " t_interp : " << t_interp << " b : " << b << endl;
                  cerr << "up : " << vupper << " low : " << vlower << endl;
                }
              }
              */
              //if(ionid==6 || ionid==1)fprintf(stdout,"WRAPAROUNDVEL: ionid= %2d irepz= % 5.3f bin_min= %5d bin_max= %5d bin= %5d v_interp= % 7.2f irepz*vstep/1.e5*nvbins= % 7.2f b= %7.2f vmin= %4.2f vbin_coord[bin]= %7.2f dvcol= % 5.3e vlower= % 7.2f vupper= % 7.2f b= %7.3f t_interp= %5.3e I.temp[k]= %5.3e I.temp[k+1]= %5.3e weight_sub= %5.3e (I.temp[k+1]-I.temp[k])*weight_sub= %5.3e\n",ionid,irepz,bin_min,bin_max,bin,v_interp,irepz*vstep/1.e5*nvbins,b,vmin,vbin_coord[bin],dvcol,vlower,vupper,I.bsys,t_interp,I.temp[k],I.temp[k+1],weight_sub,(I.temp[k+1] - I.temp[k])*weight_sub);
              //if(weight_sub>1.0)fprintf(stdout,"WRAPAROUNDVELWRONG: i= %5d j= %5d k= %5d nsubbinvar= %5d mass= %5.3e %5.3e ionid= %2d irepz= % 5.3f bin_min= %5d bin_max= %5d bin= %5d v_interp= % 7.2f irepz*vstep/1.e5*nvbins= % 7.2f b= %7.2f vmin= %4.2f vbin_coord[bin]= %7.2f dvcol= % 5.3e vlower= % 7.2f vupper= % 7.2f b= %7.3f t_interp= %5.3e I.temp[k]= %5.3e I.temp[k+1]= %5.3e weight_sub= %5.3e (I.temp[k+1]-I.temp[k])*weight_sub= %5.3e\n",i,j,k,nsubbinvar,I.mass[k+1],I.mass[k],ionid,irepz,bin_min,bin_max,bin,v_interp,irepz*vstep/1.e5*nvbins,b,vmin,vbin_coord[bin],dvcol,vlower,vupper,I.bsys,t_interp,I.temp[k],I.temp[k+1],weight_sub,(I.temp[k+1] - I.temp[k])*weight_sub);

              //if(ionid==6  && (irepz<-0.5 || irepz>0.5))fprintf(stdout,"WRAPAROUNDVEL: ionid= %2d irepz= % 5.3f bin_min= %5d bin_max= %5d bin= %5d v_interp= % 7.2f irepz*vstep/1.e5*nvbins= % 7.2f b= %7.2f vmin= %4.2f vbin_coord[bin]= %7.2f dvcol= % 5.3e vlower= % 7.2f vupper= % 7.2f b= %7.3f t_interp= %5.3e\n",ionid,irepz,bin_min,bin_max,bin,v_interp,irepz*vstep/1.e5*nvbins,b,vmin,vbin_coord[bin],dvcol,vlower,vupper,I.bsys,sqrt(t_interp));
            }
          }
        }
      }
    }
    int wrong = 0;
    for (i = 0; i < nvbins; i++)
    {
      if (I.vbins[i] == 0)
        wrong++;
      unit_col = I.Xsec / (vbin_size[i]) / (aex * aex * spece_set.spece_para.unit_Length * spece_set.spece_para.unit_Length) / MHYDR;

      I.vbins[i] *= unit_col;

      I.vbins[i] *= unit_Mass / I.atomwt;

      norm1[i] *= unit_col / I.Xsec;
      norm1[i] *= unit_Mass / I.atomwt;
#ifdef TEST_ION
      efc2 << norm1[i] << endl; //test hao
#endif
      if (isinf(norm1[i]))
        cerr << unit_col << endl;
      else if (norm1[i] < 0)
        cerr << "-0.0000000" << endl;

      if (norm[i] > 0)
      {

        I.tbins[i] /= norm[i];
        I.rhobins[i] /= norm[i];
        //if (I.rhobins[i] == 0)
        //cerr << I.rhobins[i] << endl;
        I.Zbins[i] /= norm[i];
      }
      //
      //cerr << " norm: " << norm[i] << endl;
    }
#ifdef TEST_ION
    efc2 << -999 << endl; //test hao
#endif
    cerr << "vbins0 " << wrong << endl;
    if (ionid == -1)
    {
      ion_total = I;
    }
    else
    {
      ions[ionid] = I;
    }
  }
#ifdef TEST_ION
  efc2.close();
#endif
  free(vbin_size);
  free(vbin_coord);
  free(vbin_zsize);
  free(vbin_zcoord);
  free(norm);
  free(norm1);
  free(norm_field);

  return 0;
}

void Ion_all::Freeions()
{

  for (int i = 0; i < MAXIONS; i++)
  {
    free(ions[i].mass);
    free(ions[i].vel);
    free(ions[i].temp);
    free(ions[i].rho);
    for (int j = 0; j < NMETALS; j++)
      free(ions[i].metals[j]);
    free(ions[i].vbins);
    free(ions[i].tbins);
    free(ions[i].rhobins);
    free(ions[i].Zbins);
  }

  ions.clear();
  free(ion_total.mass);
  free(ion_total.vel);
  free(ion_total.temp);
  free(ion_total.rho);
  for (int j = 0; j < NMETALS; j++)
    free(ion_total.metals[j]);
  free(ion_total.redshift);
  free(ion_total.binsize);
  free(ion_total.bincoord);
  free(ion_total.vbins);
  free(ion_total.tbins);
  free(ion_total.rhobins);
  free(ion_total.Zbins);
  free(redshift);
  free(x);
  free(y);
  free(z);
  free(gal_field);
}

#endif