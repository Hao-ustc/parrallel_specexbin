#include <iostream>
#include "include_list.h"
using namespace std;
#define RRR
class columnData
{
public:
    double r;
    double r2rvir;
    double r2rspsh;
    double y1;
    double y2;
    double y3;
    double y4;
    double y5;
    double y6;
    double y7;
    double y8;
    double y9;
    columnData()
    {
        r = -1.0;
        r2rvir = 0.0;
        r2rspsh = 0.0;
        y1 = 0.0;
        y2 = 0.0;
        y3 = 0.0;
        y4 = 0.0;
        y5 = 0.0;
        y6 = 0.0;
        y7 = 0.0;
        y8 = 0.0;
        y9 = 0.0;
    }
};

int main(int argc, char *argv[])
{
    double h = 0.72;

    double deltarin = 10; //kpc
    double deltarout = 100;
#ifdef RRR
    double rin=200;
    int Nin=(int)rin/deltarin;
    int Nout=(int)(1000.0-rin)/deltarout;
    int N = Nin+Nout + 1;
#endif
#ifdef R2Rv

    double deltalogin = 0.02;
    double deltalogout = 0.3;
    int N = 2.0 / deltalogin + 3.0 / deltalogout + 1;
#endif
#ifdef R2Rs
    double deltalogin = 0.02;
    double deltalogout = 0.3;
    int N = 2.0 / deltalogin + 3.0 / deltalogout + 1;
#endif
    int *averagenumber = new int[N + 1];
    columnData *result = new columnData[N + 1];
    for (int i = 0; i < N + 1; i++)
    {
        averagenumber[i] = 0;
    }
    vector<columnData> pp;

    char path[1000];
    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/        HI_ew_r.data", argv[1]);
    ifstream in1(path);

    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/      HeII_ew_r.data", argv[1]);
    ifstream in2(path);
    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/      CIII_ew_r.data", argv[1]);
    ifstream in3(path);
    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/       CIV_ew_r.data", argv[1]);
    ifstream in4(path);
    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/       OIV_ew_r.data", argv[1]);
    ifstream in5(path);
    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/       OVI_ew_r.data", argv[1]);
    ifstream in6(path);
    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/    NeVIII_ew_r.data", argv[1]);
    ifstream in7(path);
    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/      MgII_ew_r_all.data", argv[1]);
    ifstream in8(path);
    sprintf(path, "/data6/Hao_L/my_specexbin/result/%s/width/      SiIV_ew_r.data", argv[1]);
    ifstream in9(path);

    char path2[100];
    sprintf(path2, "/data6/Hao_L/my_specexbin/result/%s/width/average_EW.txt", argv[1]);
    ofstream out(path2);
    cerr << 1 << endl;
    while (!in8.eof())
    {
        columnData p;
        char buffer[80];
        double value;

        in8 >> buffer;
        R_assii(buffer, value);
        p.r = value; //1
        in8 >> buffer;
        R_assii(buffer, value);
        p.r2rvir = value; //2
        in8 >> buffer;
        R_assii(buffer, value);
        p.r2rspsh = value; //3
        in8 >> buffer;
        R_assii(buffer, value);
        p.y1 = value; //4

        /*
        in2 >> buffer;
        in2 >> buffer;
        in2 >> buffer;
        R_assii(buffer, value);
        p.y2 = value;
        in2 >> buffer;

        in3 >> buffer;
        in3 >> buffer;
        in3 >> buffer;
        R_assii(buffer, value);
        p.y3 = value;
        in3 >> buffer;

        in4 >> buffer;
        in4 >> buffer;
        in4 >> buffer;
        R_assii(buffer, value);
        p.y4 = value;
        in4 >> buffer;

        in5 >> buffer;
        in5 >> buffer;
        in5 >> buffer;
        R_assii(buffer, value);
        p.y5 = value;
        in5 >> buffer;

        in6 >> buffer;
        in6 >> buffer;
        in6 >> buffer;
        R_assii(buffer, value);
        p.y6 = value;
        in6 >> buffer;

        in7 >> buffer;
        in7 >> buffer;
        in7 >> buffer;
        R_assii(buffer, value);
        p.y7 = value;
        in7 >> buffer;

        in8 >> buffer;
        in8 >> buffer;
        in8 >> buffer;
        R_assii(buffer, value);
        p.y8 = value;
        in8 >> buffer;

        in9 >> buffer;
        in9 >> buffer;
        in9 >> buffer;
        R_assii(buffer, value);
        p.y9 = value;
        in9 >> buffer;
*/
        pp.push_back(p);
    }
    pp.pop_back();
    cerr << pp.size() << endl;
    in1.close();

    in2.close();
    in3.close();
    in4.close();
    in5.close();
    in6.close();
    in7.close();
    in8.close();
    in9.close();

    for (int i = 0; i < pp.size(); i++)
    {
#ifdef RRR
        if (pp[i].r > 0 && pp[i].r < rin)
        {
            int n = (int)(pp[i].r / (deltarin ));
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }
        else if (pp[i].r > rin  && pp[i].r < 1000.0 )
        {
            int n = Nin + (int)((pp[i].r - rin ) / (deltarout ));
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }
#endif
#ifdef R2Rv
        if (pp[i].r2rvir > 0 && pp[i].r2rvir < 2.0)
        {
            int n = (int)(pp[i].r2rvir / deltalogin);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }
        else if (pp[i].r2rvir > 2.0 && pp[i].r2rvir < 5.0)
        {
            int n = 2.0 / deltalogin + (int)((pp[i].r2rvir - 2.0) / deltalogout);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }
#endif
#ifdef R2Rs
        if (pp[i].r2rspsh > 0 && pp[i].r2rspsh < 2.0)
        {
            int n = (int)(pp[i].r2rspsh / deltalogin);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }
        else if (pp[i].r2rspsh > 2.0 && pp[i].r2rspsh < 5.0)
        {
            int n = 2.0 / deltalogin + (int)((pp[i].r2rspsh - 2.0) / deltalogout);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
            result[n].y5 += pp[i].y5;
            result[n].y6 += pp[i].y6;
            result[n].y7 += pp[i].y7;
            result[n].y8 += pp[i].y8;
            result[n].y9 += pp[i].y9;
        }

#endif
        else
            cerr << "wrong!!!!!!!!" << pp[i].r << endl;
    }

    for (int i = 0; i < N + 1; i++)
    {
#ifdef RRR
        if (i < Nin)
            result[i].r = deltarin  / 2.0 + (double)(i)*deltarin ;
        else if (i >= Nin)
            result[i].r = rin + deltarout  / 2.0 + (double)(i - Nin) * deltarout ;
#else
        if (i < 2.0 / deltalogin)
            result[i].r = deltalogin / 2.0 + (double)(i)*deltalogin;
        else if (i >= 2.0 / deltalogin)
            result[i].r = 2.0 + deltalogout / 2.0 + (double)(i-2.0 / deltalogin)*deltalogout;
#endif
        if (averagenumber[1] != 0)
        {
            result[i].y1 /= averagenumber[i];
            result[i].y2 /= averagenumber[i];
            result[i].y3 /= averagenumber[i];
            result[i].y4 /= averagenumber[i];
            result[i].y5 /= averagenumber[i];
            result[i].y6 /= averagenumber[i];
            result[i].y7 /= averagenumber[i];
            result[i].y8 /= averagenumber[i];
            result[i].y9 /= averagenumber[i];
            out << result[i].r << " "
                << result[i].y1 << " "
                << result[i].y2 << " "
                << result[i].y3 << " "
                << result[i].y4 << " "
                << result[i].y5 << " "
                << result[i].y6 << " "
                << result[i].y7 << " "
                << result[i].y8 << " "
                << result[i].y9 << endl;
        }
        else if (averagenumber[i] == 0)
        {
            cerr << " 00000:   " << result[i].y1 << " " << result[i].y2 << " " << result[i].y3 << endl;
        }
    }
    out.close();
    delete[] averagenumber;
    delete[] result;
}
