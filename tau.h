#ifndef IONDEF_H
#define IONDEF_H
#include "include_list.h"
using namespace std;
#define NMETALS 4
#define ZRES 3.0e-06
#define VRES 1.5e-05
#define BOXSIZE 1.0
#define HALFBOX 0.5

class ion_special
{
private:
  int i;

public:
  char name[10];
  vector<double> mass;
  vector<double> vel;
  vector<double> temp;
  vector<double> rho;
  //vector<double> metals[NMETALS];
#ifdef PHYSSPEC
  vector<double> sfr;
  vector<double> wtmass;
  vector<double> mgal;
  vector<double> dgal;
  vector<double> age;
  vector<double> nrec;
  vector<double> vlaunch;
#endif
  
  vector<double> vbins;
  vector<double> tbins;
  vector<double> rhobins;
  vector<double> Zbins;
#ifdef PHYSSPEC
  vector<double> sfrbins;
  vector<double> mgalbins;
  vector<double> dgalbins;
  vector<double> agebins;
  vector<double> nrecbins;
  vector<double> vlaunchbins;
#endif
  float lambda, fraction, Xsec, atomwt, bsys, alpha;
  int Zcolumn;
};
class ion_all
{
private:
  ion_special ion;
  
public:
  vector<ion_special> ions;
  vector<double> redshift;
  vector<double> binsize;
  vector<double> bincoord;
  vector<double> metals[NMETALS];
  void Load_ion(int a);
  void Load_ions();
};
void ion_all::Load_ion(int a)
{
  //read ionfile
}
void ion_all::Load_ions()
{
  //read ionfile
}
#endif