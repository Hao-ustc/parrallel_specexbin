#ifndef LINE_OF_SIGHT
#define LINE_OF_SIGHT
#include "include_all.h"
using namespace std;

class BIN
{
public:
    double x;
    double y;
    double z;
    double size;
    double relative_coord;
    double redshift;
};
class LOS
{
public:
    int nzbins;
    int nvbins;
    double dz;
    double zstep;
    double zbeginline;
    double xrot;
    double yrot;
    double zrot;
    double xorig;
    double yorig;
    double zorig;
    double xrotline;
    double yrotline;
    BIN *bin;

    void shortLOS(Setting &spece_set);
    void clear();
    int RotateCoords(double x, double y, double z, double theta, double phi);
    int InverseRotateCoords(double x, double y, double z, double theta, double phi);
    LOS()
    {
        xrotline=0;
        yrotline=0;
    }
};

void LOS::shortLOS(Setting &spece_set)
{
    float redshift_center;
    redshift_center = spece_set.redshift_center;
    float t = CosmicTime(redshift_center,spece_set);
    cosmopar(t,spece_set);
    float hubble = spece_cosmo.hubble;
    float unit_Velocity = spece_set.spece_para.unit_Velocity;
    spece_set.redshift = redshift_center;

    dz = BOXSIZE * hubble * unit_Velocity / CLIGHT;
    zstep = ZRES / (BOXSIZE * hubble * unit_Velocity / CLIGHT);
    nzbins = (BOXSIZE * hubble * unit_Velocity / CLIGHT) / ZRES;
    nvbins = floor(nzbins / (VRES / ZRES));
    zbeginline=0;

    bin = (BIN *)malloc((nzbins + 1) * sizeof(BIN));
    for (int i = 0; i < nzbins; i++)
    {
        bin[i].size = zstep;
        bin[i].relative_coord = ((float)i - (float)nzbins / 2) / ((float)nzbins); // (float)i/((float)nbins);
        bin[i].redshift = (nzbins / 2 - i) * ZRES + redshift_center;

        bin[i].x = spece_set.xspec;
        bin[i].y = spece_set.yspec;
        bin[i].z = spece_set.zspec + ((float)i - (float)nzbins / 2) / ((float)nzbins);
        if (bin[i].z > HALFBOX)
            bin[i].z -= BOXSIZE;
        if (bin[i].z < -HALFBOX)
            bin[i].z += BOXSIZE;
    }

    spece_set.redshift_end = bin[0].redshift;
    spece_set.redshift_begin = bin[nzbins - 1].redshift;
    double tx=spece_set.xspec;
    double ty=spece_set.yspec;
    double tz=spece_set.zspec;
    RotateCoords(tx,ty,tz,0,0);
    /*xrotline=xrot;
    yrotline=yrot;
    zbeginline = zrot+zstep;
    */
}
void LOS::clear()
{
    
    free(bin);
}
int LOS::RotateCoords(double x, double y, double z, double theta, double phi)
{
    xrot = x * cos(phi) + y * sin(phi);
    yrot = -x * sin(phi) * cos(theta) + y * cos(phi) * cos(theta) + z * sin(theta);
    zrot = x * sin(phi) * sin(theta) - y * sin(theta) * cos(phi) + z * cos(theta);
    return 0;
}
int LOS::InverseRotateCoords(double x, double y, double z, double theta, double phi)
{

  xorig =  x*cos(phi) - y*sin(phi)*cos(theta) + z*sin(phi)*sin(theta);
  yorig =  x*sin(phi) + y*cos(phi)*cos(theta) - z*cos(phi)*sin(theta);
  zorig =  y*sin(theta) + z*cos(theta);
  return 0;	
}
#endif