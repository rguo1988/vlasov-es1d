//To Use this input by changing the filename to input.h
#ifndef _input_h
#define _input_h
#include<cmath>
#include<string>
using namespace std;
class Input
{
  protected:
    //general parameters
    const double k = 0.6;
    const double L = 2 * M_PI / k; //simulaiton length
    const double T = 1; //temperature
    const double m = 1.0;
    const double vmax = 5;//10.0 * sqrt(temperature / m);
    const double e = -1.0;

    //definition of simulation constant
    static const int nx = 201;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 201;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.1;
    const int max_steps = 2000;
    //
    //data saveing
    const string data_path = "./data/";
    const int data_steps = 100;
    //special parameter
    const double delta_e = 0.1;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double f = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v, 2) / (2 * T)) * (1.0 + delta_e * cos(k * x));
        return f;
    }
    double GetIonInitDistrib(double x, double v)
    {
        return 1.0;
    }
};
#endif
