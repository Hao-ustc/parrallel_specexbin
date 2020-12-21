#include "include_list.h"
#include "gdata.h"

using namespace std;
//./bashdeal 1111 1134 0 29 1000
int main(int argc, char *argv[])
{
    char rundate1[5];
    char rundate2[5];
    int start;
    int end;
    int galaxymin;
    int galaxymax;
    int lines;


    sprintf(rundate1, "%s", argv[1]);
    sprintf(rundate2, "%s", argv[2]);
    R_assii(rundate1,start);
    R_assii(rundate2,end);
    R_assii(argv[3],galaxymin);
    R_assii(argv[4],galaxymax);
    R_assii(argv[5],lines);

    

    char filepalace[100];
    char result_path[100];
    sprintf(result_path, "./bash", argv[1]);
    
    char path3[100];
        sprintf(path3, "%s/run_multi.sh", result_path);
        ofstream tfsh(path3);
    for (int i = start; i <= end; i++)
    {
        char rundate[5];
        sprintf(rundate,"%04d",i);
        char command[200];
         
        for(int galaxy=galaxymin;galaxy<=galaxymax;galaxy++)
        {
            tfsh << "./data_resolve_tau"
             << " " << rundate << " 7 " << galaxy << " " << galaxy <<" "<<lines<< " &" << endl;
        }
        
        tfsh<<"wait"<<endl;
        sprintf(command,"./assemble %s 7 %03d %03d &",rundate,galaxymin,galaxymax);
        tfsh<<command<<endl;
        tfsh<<"wait"<<endl;
        sprintf(command,"./average_EW %s &",rundate);
        tfsh<<command<<endl;
        tfsh<<"wait"<<endl;

       
    }
    tfsh.close();
    
}