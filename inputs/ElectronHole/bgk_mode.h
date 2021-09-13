//To Use this input by changing the filename to input.h
/***********************************
 * BGK mode
 ***********************************/
#ifndef _input_h
#define _input_h
#include<cmath>
#include<string>
using namespace std;
class Input
{
  protected:
    //general parameters
    const double L = 15.0; //simulaiton length
    const double T = 1.0; //temperature
    const double m = 1.0;
    const double vmax = 5;//10.0 * sqrt(temperature / m);
    const double e = -1.0;

    //definition of simulation constant
    static const int nx = 101;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 2001;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.02;
    const int max_steps = 300000;

    //special parameters
    const double uae = 0.05;
    //const double uai = 0.6;
    //data recording
    const string data_path = "./data/";
    const int data_steps = 30000;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double rv = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v, 2) / (2 * T));
        double rx = 1.0 + uae * cos(2 * M_PI * x / L) ;
        return rv * rx;
    }
    double GetIonInitDistrib(double x, double v)
    {
        //double r =   1.0 + uai * cos(2 * M_PI * x / L) ;
        double r =   1.0;
        return r;
    }
};
#endif
