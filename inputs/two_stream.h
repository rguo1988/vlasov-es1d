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
    const double L = 10.0; //simulaiton length
    const double T = 0.1; //temperature
    const double m = 1.0;
    const double vmax = 5;//10.0 * sqrt(temperature / m);
    const double e = -1.0;

    //definition of simulation constant
    static const int nx = 301;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 301;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 10000;
    //
    //data saveing
    const string data_path = "./data/";
    const int data_steps = 10;

    //special parameters
    const double u1 = 1.6; //drift speed of stream 1
    const double u2 = -1.4; //drift speed of stream 2

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double r1 = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v - u1, 2) / (2 * T));
        double r2 = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v - u2, 2) / (2 * T)) * (1 + 0.02 * cos(2 * M_PI * x / L));
        return 0.5 * (r1 + r2);
    }
    double GetIonInitDistrib(double x, double v)
    {
        return 1.0;
    }
};
#endif
