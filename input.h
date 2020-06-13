//To Use this input by changing the filename to input.h
/***********************************
 * kappa stationary distribution with non-uniform temperature
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
    const double L = 20.0; //simulaiton length
    const double T = 1.0; //temperature
    const double m = 1.0;
    const double vmax = 20;
    const double e = -1.0;

    //definition of simulation constant
    static const int nx = 101;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 1001;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 10000;

    //special parameters
    const double uae = 0.1;
    const double uai = 0.1;
    const double kappa = 5.0;
    //data recording
    const string data_path = "./data/";
    const int data_steps = 1000;
    const int data_num = max_steps / data_steps;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double rx = 1.0 + uae * cos(2 * M_PI * x / L) ;
        //double Tx = T * rx;
        //double rv = sqrt(1.0 / (2 * M_PI * Tx * kappa)) * tgamma(kappa + 1.5) / tgamma(kappa + 1.0) * pow(1 + pow(v, 2) / (2 * Tx * kappa), -kappa-1.5);
        double rv = sqrt(1.0 / (2 * M_PI * T * kappa)) * tgamma(kappa + 1.5) / tgamma(kappa + 1.0) * pow(1 + pow(v, 2) / (2 * T * kappa), -kappa-1.5);
        //double rv = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v, 2) / (2 * T));
        return rx*rv;
    }
    double GetIonInitDistrib(double x, double v)
    {
        double r =   1.0 + uai * cos(2 * M_PI * x / L) ;
        //double r = 1.0;
        return r;
    }
};
#endif
