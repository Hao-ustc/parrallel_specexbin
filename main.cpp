#include <iostream>
//#define para
#include "function.h"
using namespace std;

int main(int argc,char *argv[])
{
    
    //sprintf(spece_set.TIME,"0419");
    
    //spece_set.runtime=0;
    int test = -1;
    char los[100];

    Setting spece_set;
    sprintf(spece_set.TIME,"%s",argv[1]);
    sprintf(spece_set.galaxy,"%s",argv[2]);


    int galaxyposion=0;
    int parallelcore=0;
    R_assii(argv[2],galaxyposion);
    R_assii(argv[3],parallelcore);
    sprintf(los, "./result/%s/parallel/los%03d_%01d.txt",argv[1],galaxyposion,parallelcore);
    sprintf(spece_set.file_redshift, "/data6/Hao_L/result/pure_redshift.txt");
    sprintf(spece_set.spec_ion_filename, "%sspecions_i9.dat", spece_set.prefix);

    if ((spece_set.LOSfile = fopen(los, "r")) == NULL)
    {
        cerr << "Could not open file " << los << endl;
        return 0;
    }

    vector<Gas_1> gp;

    galaxy_map density_map;
    sprintf(density_map.filedestination,"./analyze/density");
    int count=0;
    while (!feof(spece_set.LOSfile) )
    {
        //if(test<1510) continue;
        spece_set.boxsize = 60; //in Mpc_h
        spece_set.flux_fac = 1.0;
        spece_set.spece_para.totMass = 0.258; //0.28; //0.238; //0.28;
        spece_set.spece_para.Lambda = 0.742;  //0.72; //0.762; //0.72;
        spece_set.spece_para.Omega_b = 0.045; //0.046; //0.0418; //0.046;
        spece_set.spece_para.H_0 = 73;        //73;
        spece_set.spece_para.h = 0.01 * spece_set.spece_para.H_0;
        
        spece_set.load();
        cerr << spece_set.short_filenum << endl;
        LOS shortlos;
        shortlos.shortLOS(spece_set);

        spece_ionall.Load(shortlos);

        if ((spece_set.short_filenum - spece_set.last_filenum) != 0)
        {
            gp.clear();
            Readdata(spece_set.short_filenum, gp);
            //function_test(gp);
            spece_set.last_filenum = spece_set.short_filenum;
        }

        Spec_particles spece_particles;
        spece_particles.Check_Partical_Los(gp);
        spece_particles.SmoothSpec(shortlos);

        spece_ionall.Tau(shortlos);
        int pospoint=OutTau(shortlos);

        count++;
       /* for(int i=0;i<9;i++)
        {
            
            double density_value=spece_ionall.ions[i].rhobins[pospoint];
            density_map.load(i,count,density_value);
        }
        */
        cleanworkplace(shortlos);
    }
    //density_map.pixel_sort();
    //density_map.output();
    fclose(spece_set.LOSfile);
}