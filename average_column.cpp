#include <iostream>
#include "include_all.h"
using namespace std;

class columnData
{
public:
    double r;
    double y1;
    double y2;
    double y3;
    double y4;

    columnData()
    {
        r = -1.0;
        y1 = 0.0;
        y2 = 0.0;
        y3 = 0.0;
        y4 = 0.0;
    }
};

int main(int argc, char *argv[])
{
    double deltar = 5; //kpc
    int averagenumber[52];
    columnData result[52];
    for (int i = 0; i < 52; i++)
    {
        averagenumber[i] = 0;
    }
    vector<columnData> pp;

    char path[100];
    sprintf(path, "./result/%s/width/column.data", argv[1]);
    ifstream in(path);

    while (!in.eof())
    {
        columnData p;
        char buffer[80];
        double value;
        in >> buffer;//1
        

        in >> buffer;//2
        R_assii(buffer, value);
        p.r = value;
        in >> buffer;//3
        R_assii(buffer, value);
        p.y1 = value;
        in >> buffer;//4
        in >> buffer;//5
        in >> buffer;//6
        in >> buffer;//7
        in >> buffer;//8
        R_assii(buffer, value);
        p.y2 = value;
        in >> buffer;//9
        R_assii(buffer, value);
        p.y3 = value;
        in >> buffer;//10
        in >> buffer;//11
        R_assii(buffer, value);
        p.y4 = value;
        pp.push_back(p);
    }
    pp.pop_back();
    in.close();
    for (int i = 0; i < pp.size(); i++)
    {
        if (pp[i].r > 0 && pp[i].r < 250)
        {
            int n = (int)(pp[i].r / deltar);
            averagenumber[n]++;
            result[n].y1 += pp[i].y1;
            result[n].y2 += pp[i].y2;
            result[n].y3 += pp[i].y3;
            result[n].y4 += pp[i].y4;
        }
        else if(pp[i].r<0)
            cerr << "wrong!!!!!!!!" << endl;
    }
    char path2[100];
    sprintf(path2, "./result/%s/width/average.txt", argv[1]);
    ofstream out(path2);

    for (int i = 0; i < 52; i++)
    {
        result[i].r = deltar / 2.0 + (double)(i)*deltar;
        if (averagenumber[1] != 0)
        {
            result[i].y1 /= averagenumber[i];
            result[i].y2 /= averagenumber[i];
            result[i].y3 /= averagenumber[i];
            result[i].y4 /= averagenumber[i];
            out << result[i].r << " "
                << log10(result[i].y1) << " "
                << log10(result[i].y2) << " "
                << log10(result[i].y3) << " "
                << log10(result[i].y4) << endl;
        }
        else if (averagenumber[i] == 0)
        {
            cerr << " 00000:   " << result[i].y1 << " " << result[i].y2 << " " << result[i].y3 << " " << result[i].y4 << endl;
        }
    }
    out.close();
}