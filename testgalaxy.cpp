#include "include_list.h"
#include "gdata.h"
using namespace std;
#define N2 30
#define Pi 3.1415926
//#define deltah 0.01;
//#define PALACE
#define ROUND

int main(int argc, char *argv[])
{

    double min = 0;
    double max = 0;
    int filenum = 200;
    int galaxy_flag = 0;
    R_assii(argv[2], min);
    R_assii(argv[3], max);
    R_assii(argv[4], filenum);
    R_assii(argv[5], galaxy_flag); //1:central
    double deltah = 3.0 / N2;

    char result_path[100];

    int direction;
    double redshift_center;
    double xspec, yspec, zspec;
    double boxsize = 60;

    //try to be a function
    vector<Galaxy> ga;

    Readdata(filenum, ga);

    double mass_max = 0;
    int count = 0;
    int countcentral1 = 0;
    int countcentral2 = 0;
    int countcentral3 = 0;
    int countsatelite1 = 0;
    int countsatelite2 = 0;
    int countsatelite3 = 0;
    vector<int> ga_choose;
    //ga_choose.push_back(-1);
    char fileout[200];
    sprintf(fileout, "./testgalaxy.txt");
    ofstream tf(fileout);
    for (int i = 0; i < ga.size(); i++)
    {
        tf << log10(ga[i].halomass * 0.72) << " " << log10(ga[i].center_coordinates[1] * 0.72) << " " << log10(ga[i].center_coordinates[2] * 0.72) << " " << log10(ga[i].Mass * 0.72 * 0.72) << " " << ga[i].flag << endl;
        if (log10(ga[i].halomass) > min && log10(ga[i].halomass) < max)
        {
            count++;
            if (log10(ga[i].Mass) > 10.6 && log10(ga[i].Mass) < 10.7 && ga[i].flag == 1)
            {
                countcentral1++;
            }
            if (log10(ga[i].Mass) > 10.7 && log10(ga[i].Mass) < 10.8 && ga[i].flag == 1)
            {
                countcentral1++;
            }
            if (log10(ga[i].Mass) > 10.8 && ga[i].flag == 1)
            {
                countcentral3++;
            }
            if (log10(ga[i].Mass) < 10.2 && ga[i].flag == -1)
            {
                countsatelite1++;
            }
            if (log10(ga[i].Mass) > 10.2 && log10(ga[i].Mass) < 10.6 && ga[i].flag == -1)
            {
                countsatelite2++;
            }
            if (log10(ga[i].Mass) > 10.6 && ga[i].flag == -1)
            {
                countcentral3++;
            }
        }
    }
    tf.close();
    cerr << count << endl;
    cerr << "central1 " << countcentral1 << endl;
    cerr << "central2 " << countcentral2 << endl;
    cerr << "central3 " << countcentral3 << endl;
    cerr << "satelite1 " << countsatelite1 << endl;
    cerr << "satelite2 " << countsatelite2 << endl;
    cerr << "satelite3 " << countsatelite3 << endl;
    cerr << ga_choose.size() << endl;
}
