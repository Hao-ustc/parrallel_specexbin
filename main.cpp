#include <iostream>
//#define para
#include "function.h"
#include "datachange.h"

using namespace std;

int main(int argc, char *argv[])
{

    //sprintf(spece_set.TIME,"0419");

    //spece_set.runtime=0;
    int NUM_THREADS = 64;
    int test = -1;

    Setting spece_set;
    sprintf(spece_set.TIME, "%s", argv[1]);
    sprintf(spece_set.galaxy, "%s", argv[2]);
    sprintf(spece_set.file_redshift, "/data6/Hao_L/result/pure_redshift.txt");
    sprintf(spece_set.spec_ion_filename, "%sspecions_i9.dat", spece_set.prefix);
    spece_set.boxsize = 80; //in Mpc_h
    spece_set.flux_fac = 1.0;
    spece_set.spece_para.totMass = 0.258; //0.28; //0.238; //0.28;
    spece_set.spece_para.Lambda = 0.742;  //0.72; //0.762; //0.72;
    spece_set.spece_para.Omega_b = 0.045; //0.046; //0.0418; //0.046;
    spece_set.spece_para.H_0 = 72;        //73;
    spece_set.spece_para.h = 0.01 * spece_set.spece_para.H_0;

    spece_set.phi = -234.0;

    int galaxyposion = 0;
    int parallelcore = 0;
    R_assii(argv[2], galaxyposion);
    R_assii(argv[3], parallelcore);

    Los_poss los_poss;
    R_assii(argv[4], los_poss.los_numbers);
    char los[100];
    sprintf(los, "./result/%s/parallel/los%03d_%01d.txt", argv[1], galaxyposion, parallelcore);
    if ((spece_set.LOSfile = fopen(los, "r")) == NULL)
    {
        cerr << "Could not open file " << los << endl;
        return 0;
    }
    los_poss.in_los = spece_set.LOSfile;
    fclose(spece_set.LOSfile);
    sprintf(los, "./result/%s/parallel/palace%03d_%01d.txt", argv[1], galaxyposion, parallelcore);
    if ((spece_set.LOSfile = fopen(los, "r")) == NULL)
    {
        cerr << "Could not open file " << los << endl;
        return 0;
    }
    los_poss.in_palace = spece_set.LOSfile;
    fclose(spece_set.LOSfile);
    los_poss.load();
    fclose(los_poss.in_los);
    fclose(los_poss.in_palace); //finished

    spece_set.load(los_poss.los_pos_vector[0]); //to get filenum
    GASs gass;
    gass.load(spece_set.short_filenum);
    datafit(gass.void_gas,spece_set,gass.b);
    vector<GAS> gp;
    gp.clear();
    gp.swap(gass.void_gas);

    sprintf(los, "./result/0_a_width/%s/lines_rp_rper_ews.txt", argv[1]); //output
    ofstream out_ew(los);

    R_assii(argv[5], NUM_THREADS);
    omp_set_num_threads(NUM_THREADS);
#pragma omp parallel firstprivate(spece_set)
    {
        int thread_id;
        thread_id = omp_get_thread_num();

        for (int losnum = thread_id; losnum < 100; losnum += NUM_THREADS)
        {
            spece_set.load(los_poss.los_pos_vector[losnum]);

            LOS shortlos;

            shortlos.shortLOS(spece_set);
            //cerr << "&&& 1" << endl;
            Ion_all spece_ionall;
            spece_ionall.Load(shortlos, spece_set);
            //cerr << "&&& 2" << endl;
            Spec_particles spece_particles;
            spece_particles.Check_Partical_Los(gp, spece_set);
            //cerr << "&&& 3" << endl;
            spece_particles.SmoothSpec(shortlos, spece_ionall, spece_set,gass.b);
            //cerr << "&&& 4" << endl;
            spece_ionall.Tau(shortlos, spece_set);
            //cerr << "&&& 5" << endl;
            my_EW my_ew;
            int pospoint = OutTau(shortlos, spece_ionall, spece_set, my_ew.redshift_track);

            my_ew.data_resolve(spece_ionall, spece_set);
            int t_main_lines = spece_set.lines_spec;
#pragma omp critical
            {
                out_ew << t_main_lines << " "
                       << los_poss.los_pos_vector[t_main_lines].rp << " "
                       << los_poss.los_pos_vector[t_main_lines].rp_rvir << " ";
                for (int ion = 0; ion < my_ew.Nions; ion++)
                {
                    out_ew << my_ew.ew[ion] << " ";
                }
                out_ew << endl;
            }

            cleanworkplace(shortlos, spece_ionall);
        }
    }
    out_ew.close();
    /*
    galaxy_map density_map;
    sprintf(density_map.filedestination, "./analyze/density");
    int count = 0;
    while (!feof(spece_set.LOSfile))
    {
        //if(test<1510) continue;

        
        cerr << spece_set.short_filenum << endl;
        

        ///*
        if ((spece_set.short_filenum - spece_set.last_filenum) != 0)
        {
            gp.clear();
            Readdata(spece_set.short_filenum, gp);
            //function_test(gp);
            spece_set.last_filenum = spece_set.short_filenum;
        }
        //

        

        count++;
        
        
    }
    //density_map.pixel_sort();
    //density_map.output();
    */
}