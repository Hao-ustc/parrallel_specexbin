#include "include_list.h"
#include "gdata.h"

using namespace std;
//./build_pbs 0629 0 65 fat2 20 20
int main(int argc, char *argv[])
{
    char rundate[5];
    int galaxymin;
    int galaxymax;
    int memory;
    int core;
    int wait;
    sprintf(rundate, "%s", argv[1]);
    sscanf(argv[2], "%d", &galaxymin);
    sscanf(argv[3], "%d", &galaxymax);
    sscanf(argv[5], "%d", &memory);
    sscanf(argv[6], "%d", &core);
    wait=90/core;

    char filepalace[100];
    char result_path[100];
    sprintf(result_path, "./bash/%s", argv[1]);
    sprintf(filepalace, "%s/palace%03d.test", result_path, 0);
    ofstream tf(filepalace);

    if (!(tf.is_open()))
    {
        cerr << mkdir(result_path, S_IRWXU) << endl;
    }
    tf.close();
    char path3[100];
        sprintf(path3, "%s/run.sh", result_path);
        ofstream tfsh(path3);
    for (int i = galaxymin; i <= galaxymax; i++)
    {

        /*
        char path1[100];
        sprintf(path1, "%s/galaxymain%03d.sh", result_path, i);
        ofstream tfg(path1);
        char path2[100];
        sprintf(path2,"/data6/Hao_L/my_specexbin/bash/0910/main%03d",i);
        tfg << "for file in " <<  path2 << "*.sh;do bash $file;done" << endl;
        tfg<<"wait"<<endl;
        tfg.close();
        sprintf(path2,"/data6/Hao_L/my_specexbin/bash/0910/galaxymain%03d.sh",i);
        tfsh << "for file in " <<  path2 << ";do bash $file;done" << endl;
        */

        for (int k = 0; k < core; k++)
        {
            char galaxy[5];
            sprintf(galaxy, "%03d", i);
            /*
            char path[100];
            sprintf(path, "%s/main%03d_%01d.sh", result_path, i, k);
            ofstream tf(path);
            */
            /*
            tf << "#PBS -N main" << endl
               << "#PBS -l nodes=1:ppn=1" << endl
               << "#PBS -q " << argv[4] << endl
               << "#PBS -j oe" << endl
               << "#PBS -l walltime=10000:00:00" << endl
               << "#PBS -l mem=" << memory << "G" << endl
               << "cd /data6/Hao_L/my_specexbin" << endl
               << "./main"
               << " " << rundate << " " << galaxy << " " << k << endl;
            */
            tfsh << "./main"
               << " " << rundate << " " << galaxy << " " << k
               << " &" << endl;
        }
        if((i+1)%wait!=0)continue;
        tfsh << "wait"<<endl;
    }
    tfsh.close();
}