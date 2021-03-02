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
    int argnumber = 1;
    Los_poss los_poss;

    Setting spece_set;
    sprintf(spece_set.TIME, "%s", argv[argnumber]);
    argnumber++;
    R_assii(argv[argnumber], NUM_THREADS);
    argnumber++;
    
    R_assii(argv[argnumber], los_poss.los_numbers); //line 的总数
    cerr<<1111<<endl;
    sprintf(spece_set.file_redshift, "/data6/Hao_L/result/simulation/voidA/redshiftmovie.txt");
    sprintf(spece_set.spec_ion_filename, "%sspecions_i9.dat", spece_set.prefix);
    spece_set.boxsize = 80; //in Mpc_h
    spece_set.flux_fac = 1.0;
    spece_set.spece_para.totMass = 0.258; //0.28; //0.238; //0.28;
    spece_set.spece_para.Lambda = 0.742;  //0.72; //0.762; //0.72;
    spece_set.spece_para.Omega_b = 0.045; //0.046; //0.0418; //0.046;
    spece_set.spece_para.H_0 = 72;        //73;
    spece_set.spece_para.h = 0.01 * spece_set.spece_para.H_0;

    spece_set.phi = -234.0;

    
    char los[100];
    sprintf(los, "./result/%s/los.txt", argv[1]);
    if ((los_poss.in_los = fopen(los, "r")) == NULL)
    {
        cerr << "Could not open file " << los << endl;
        return 0;
    }
    
    sprintf(los, "./result/%s/palace.txt", argv[1]);
    if ((los_poss.in_palace = fopen(los, "r")) == NULL)
    {
        cerr << "Could not open file " << los << endl;
        return 0;
    }
    
    cerr<<111<<endl;
    los_poss.load();
    fclose(los_poss.in_los);
    fclose(los_poss.in_palace); //finished
    cerr<<1<<endl;
    //line has already been read in

    spece_set.load(los_poss.los_pos_vector[0]); //to get filenum
    GASs gass;
    gass.load(spece_set.short_filenum);
    datafit(gass.void_gas, spece_set, gass.b);
    vector<GAS> gp;
    gp.clear();
    gp.swap(gass.void_gas);
    cerr<<11<<endl;
    sprintf(los, "./result/0_a_width/%s/lines_rp_rper_ews.txt", argv[1]); //output
    ofstream out_ew(los);

    omp_set_num_threads(NUM_THREADS);
    
#pragma omp parallel firstprivate(spece_set)
    {
        int thread_id;
        thread_id = omp_get_thread_num();

        for (int losnum = thread_id; losnum < los_poss.los_numbers; losnum += NUM_THREADS)
        {
            spece_set.load(los_poss.los_pos_vector[losnum]);

            LOS shortlos;

            shortlos.shortLOS(spece_set);

            Ion_all spece_ionall;
            spece_ionall.Load(shortlos, spece_set);

            Spec_particles spece_particles;
            spece_particles.Check_Partical_Los(gp, spece_set);

            spece_particles.SmoothSpec(shortlos, spece_ionall, spece_set, gass.b);

            spece_ionall.Tau(shortlos, spece_set);

            my_EW my_ew;
            int pospoint = OutTau(shortlos, spece_ionall, spece_set, my_ew.redshift_track);

            my_ew.data_resolve(spece_ionall, spece_set);
            int t_main_lines = spece_set.lines_spec;
#pragma omp critical
            {
                out_ew << t_main_lines << " "
                       << spece_set.galaxyindex_spec << " "
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