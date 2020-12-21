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
    int coresnumber;
    int day=16;
    int hour=21;
    int minute=10;
    int month=10;
    int year=20;
    sprintf(rundate, "%s", argv[1]);
    sscanf(argv[2], "%d", &galaxymin);
    sscanf(argv[3], "%d", &galaxymax);
    sscanf(argv[5], "%d", &memory);
    sscanf(argv[6], "%d", &core);
    sscanf(argv[7], "%d", &coresnumber);
    sscanf(argv[8], "%d", &hour);
    sscanf(argv[9], "%d", &minute);
    sscanf(argv[10], "%d", &month);
    sscanf(argv[11], "%d", &year);

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

    int countspecial=0;
    for (int i = galaxymin; i <= galaxymax; i++)
    {
        
        if(minute>=60)
        {
            minute-=60;
            hour++;
        }
        
        if(hour==24)
        {
            hour-=24;
            day++;
        }
        char date_time[100];
        sprintf(date_time, "%02d%02d%02d",day, hour,minute);
        hour++;
        
        for (int k = 0; k < core; k++)
        {
            /*
            char node[2];
            int nodenumber=(((i-galaxymin)*core+k)%8+1);
            if(nodenumber==3) 
            {
                countspecial++;
                nodenumber=(countspecial%5+4);
            }
            */
            //sprintf(node,"%02d",nodenumber);
            char galaxy[5];
            sprintf(galaxy, "%03d", i);
            char path[100];
            sprintf(path, "%s/main%03d_%01d.pbs", result_path, i, k);
            ofstream tf(path);
            tf << "#PBS -N main" << endl
               
               << "#PBS -l nodes=cu09"<<":ppn=" <<coresnumber<< endl
               << "#PBS -q " << argv[4] << endl
               << "#PBS -j oe" << endl
               << "#PBS -l walltime=10000:00:00" << endl
               << "#PBS -l mem=" << memory << "G" << endl
               << "cd /data6/Hao_L/my_specexbin" << endl
               << "./main"
               << " " << rundate << " " << galaxy << " " << k << endl;

            //<< "#PBS -a " << date_time << endl
            /*
            tf << "nohup ./main"
               << " " << rundate << " " << galaxy << " " << k
               << " >main" << i << "_" << k << ".out"
               << "2>&1 &" << endl;
               */
        }
    }
}
