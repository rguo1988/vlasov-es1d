#include<iostream>
#include<iomanip>
#include<fstream>
#include"diagnose.h"
using namespace std;
void OutputMatrix(string filename, MatrixXd f)
{
    ofstream ofile;
    ofile.open(filename.c_str());
    ofile << f;
    ofile.clear();
    ofile.close();
}
void OutputVector(string filename, vector<double> v)
{
    ofstream ofile;
    ofile.open(filename.c_str());
    for(auto out_v : v)
    {
        ofile << setprecision(13) << out_v << endl;
    }
    ofile.clear();
    ofile.close();
}
