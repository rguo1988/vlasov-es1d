#ifndef _diagnose_h
#define _diagnose_h
#include<eigen3/Eigen/Core>
#include<string>
#include<vector>
using namespace std;
using namespace Eigen;
void OutputMatrix(string filename, MatrixXd f);
void OutputVector(string filename, vector<double> v);
#endif
