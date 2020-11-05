#ifndef SETTING_H
#define SETTING_H
#include "include_all.h"

using namespace std;

class Parameter
{
public:
  double totMass;
  double Lambda;
  double Omega_b;
  double H_0;
  double h;

  double unit_Time, unit_Density, unit_Length, unit_Mass, unit_Velocity, unit_Temperature, unit_DUDT;
  double H0, t0, aex3, aexhub;

  void load(double set_boxsize);
  void clear();
};

class Setting
{
public:
  char TIME[5];
  char galaxy[5];
  //int runtime;
  FILE *LOSfile;
  //char *los;            //while use we need to give it a adreess
  double xspec;         //load range -0.5 to 0.5
  double yspec;         //load
  double zspec;         //load
  int direction;        //load
  char namesuffix[200]; //load
  double theta;         //to be 0
  double phi;           //to be 0

  double redshift_center; //we will set it while use it
  double redshift_begin;
  double redshift_end;
  double redshift;
  int short_filenum;
  int last_filenum;

  double boxsize;
  double flux_fac;
  double taufact0;

  char prefix[50];
  char spec_ion_filename[100];
  char file_redshift[100];

  Parameter spece_para;

  int load(); //temply for short
  //void setpara(Parameter &para);
  int Check_Z_File(double redshift);
  void clear();
  Setting()
  {
    theta = 0;
    phi = 0;
    last_filenum = 0;
    taufact0 = 1.0;
    sprintf(prefix, "./prepare/ionfiles/");
  }
};

int Setting::Check_Z_File(double redshift)
{
  int File = -1;
  ifstream in(file_redshift);
  for (int i = 200; i >= 0; i--)
  {
    double redshift_movie;
    char buffer[80];
    in >> buffer;
    R_assii(buffer, redshift_movie);
    if (fabs(redshift_movie - redshift) <= 0.0001)
    {
      File = i;
      in.close();
      return File;
    }
  }
  in.close();
  return File;
}

int Setting::load()
{

  fscanf(LOSfile, "%lf %lf %lf %lf %d", &redshift_center, &xspec, &yspec, &zspec, &direction);
  //cerr << "read file  " << xspec<< " " << yspec <<" "<<zspec<<" "<<redshift_center<<endl; 
  redshift_begin = redshift_center;
  redshift_end = redshift_center;
  redshift = redshift_center;

  short_filenum = Check_Z_File(redshift_center);

  if (xspec > HALFBOX) //-------------why???????????????________________________________
    xspec -= BOXSIZE;
  if (xspec < -HALFBOX)
    xspec += BOXSIZE;
  if (yspec > HALFBOX)
    yspec -= BOXSIZE;
  if (yspec < -HALFBOX)
    yspec += BOXSIZE;
  if (zspec > HALFBOX)
    zspec -= BOXSIZE;
  if (zspec < -HALFBOX)
    zspec += BOXSIZE; //___________________________end____________________________________

  if (direction == 1)
  {
    double hold = zspec;
    zspec = yspec;
    yspec = xspec;
    xspec = hold;
  }
  if (direction == 0)
  {
    double hold = xspec;
    xspec = yspec;
    yspec = zspec;
    zspec = hold;
  }
  spece_para.load(boxsize);
  if (flux_fac > 0)
    cerr << "the setting and parasetting hasfinished. This should happen only once" << endl;
}

//Setting spece_set;
void Parameter::load(double set_boxsize)
{
  float L;
  double boxsize=set_boxsize;

  set_boxsize /= h;
  L = set_boxsize;
  //set_boxsize=boxsize;
  if (Lambda > 0.01 && totMass < 1.0)
  {
    if (Lambda != 1. - totMass)
    {
      /*			fprintf(stderr,"Setting Lambda = %g\n",1.-totMass);*/
      Lambda = 1. - totMass;
    }
  }
  else
    Lambda = 0.0;

  H0 = (8.0 / 3.0) * Pi;
  H0 = sqrt(H0);
  t0 = 2. / (3 * H0);

  unit_Time = H0 * Mpc;
  unit_Time /= 100.0 * h * km;

  unit_Density = 1.879E-29 * h * h;
  unit_Length = L * Mpc; //??is the unit is calculate by special melogy????
  unit_Mass = unit_Density * unit_Length * unit_Length * unit_Length;
  unit_Velocity = unit_Length / unit_Time; //cm/s
  unit_Temperature = unit_Velocity * unit_Velocity * m_p;
  unit_Temperature /= k_B;
  unit_DUDT = unit_Density * unit_Velocity * unit_Velocity / unit_Time;
}

void Setting::clear()
{
  xspec == -200.0;  //load range -0.5 to 0.5
  yspec == -200.0;  //load
  zspec == -200.0;  //load
  direction == -20; //load

  theta == -200.0; //to be 0
  phi == -200.0;   //to be 0

  redshift_center == -200.0; //we will set it while use it
  redshift_begin == -200.0;
  redshift_end == -200.0;
  redshift == -200.0;
  short_filenum == -1;
  last_filenum == -1;

  //boxsize == -200.0;
  //flux_fac == -200.0;
  //taufact0 == -200.0;

  spece_para.clear();
}
void Parameter::clear()
{

  unit_Time = -20;
  unit_Density = -20;
  unit_Length = -20;
  unit_Mass = -20;
  unit_Velocity = -20;
  unit_Temperature = -20;
  unit_DUDT = -20;
  H0 = -20;
  t0 = -20;
  aex3 = -20;
  aexhub = -20;
}
#endif