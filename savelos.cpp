#include <iostream>
#include "include_list.h"
#include <set>
#include "gdata.h"
using namespace std;

int main(int argc, char *argv[])
{
    char result_path[100];
    sprintf(result_path, "./result/%s", argv[1]);
    char fileall[100];
    for (int i = 0; i <= 5; i++)
    {
        sprintf(fileall, "%s/los%03d.txt", result_path, i);
        ofstream out(fileall);
        //out<<1<<endl;
        for (int j = 0; j < 10; j++)
        {
            char files[100];
            sprintf(files, "%s/parallel/los%03d_%01d.txt", result_path, i,j);
            ifstream in(files);
            int count=0;
            while (!(in.eof())&&count<100)
            {
                //cerr<<1<<endl;
                char buffer[100];
                double number1;
                int number2;
                in>>buffer;//1
                R_assii(buffer,number1);
                out<<number1<<" ";
                in>>buffer;//2
                R_assii(buffer,number1);
                out<<number1<<" ";
                in>>buffer;//3
                R_assii(buffer,number1);
                out<<number1<<" ";
                in>>buffer;//4
                R_assii(buffer,number1);
                out<<number1<<" ";
                in>>buffer;//5
                R_assii(buffer,number2);
                out<<number2<<endl;
                count++;
            }
            //cerr<<files<<endl;
            in.close();

        }
        //cerr<<fileall<<endl;
        out.close();
    }
}