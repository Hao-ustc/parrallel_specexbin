#include "include_list.h"
#include "gdata.h"

using namespace std;
//./build_pbs_resolve 0629 0 65 squeue 20
int main(int argc, char *argv[])
{
    char rundate[5];
    int galaxymin;
    int galaxymax;
    int memory;
    int lines;
    //int core;
    sprintf(rundate, "%s", argv[1]);
    sscanf(argv[2], "%d", &galaxymin);
    sscanf(argv[3], "%d", &galaxymax);
    sscanf(argv[5], "%d", &memory);
    sscanf(argv[6], "%d", &lines);
    //sscanf(argv[6], "%d", &core);

    char filepalace[100];
    char result_path[100];
    sprintf(result_path, "./pbs/%s", argv[1]);
    sprintf(filepalace, "%s/palace%03d.test", result_path, 0);
    ofstream tf(filepalace);

    if (!(tf.is_open()))
    {
        cerr << mkdir(result_path, S_IRWXU) << endl;
    }
    tf.close();

    char path[100];
    sprintf(path, "%s/data_resolve.sh", result_path);
    ofstream tfsh(path);
    for (int i = galaxymin; i <= galaxymax; i++)
    {

        char galaxy[5];
        sprintf(galaxy, "%03d", i);
        /*
            char path[100];
            sprintf(path, "%s/data_resolve%03d.pbs", result_path, i);
            ofstream tf(path);
            
            tf << "#PBS -N main" << endl
               << "#PBS -l nodes=1:ppn=1" << endl
               << "#PBS -q " << argv[4] << endl
               << "#PBS -j oe" << endl
               << "#PBS -l walltime=10000:00:00" << endl
               << "#PBS -l mem=" << memory << "G" << endl
               << "cd /data6/Hao_L/my_specexbin" << endl
               << "./data_resolve_tau"
               << " " << rundate << " 7 " << galaxy << " " << galaxy << endl;
               tf.close();
            */
        tfsh << "./data_resolve_tau"
             << " " << rundate << " 7 " << galaxy << " " << galaxy <<" "<<lines<< " &" << endl;
    }
}