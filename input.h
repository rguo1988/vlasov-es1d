//To Use this input by changing the filename to input.h
/***********************************
 * temperature wave in uniform kappa plasma
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
    const double k = 0.05;
    const double L = 1 * 2 * M_PI / k; //simulaiton length
    //const double L = 80; //simulaiton length
    //const double k = 5 * 2 * M_PI / L;
    const double T = 1.0; //temperature
    const double m = 1.0;
    const double vmax = 50;
    const double e = -1.0;
    const double lambda = sqrt(T / e / e);

    //definition of simulation constant
    static const int nx = 1001;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 1001;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 15000;

    //special parameters
    const double d = 0.001;
    const double kappa = 5;
    //data recording
    const string data_path = "./data/";
    const int data_steps = 100;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double fx = 1 + d * cos(k * x);
        double fv = 0.0;
        if (kappa == 0.0)
            fv = sqrt(1.0 / (2 * M_PI * T)) * exp(-v * v / (2 * T));
        else
            //fv = sqrt(1.0 / (2 * M_PI * T * kappa)) * tgamma(kappa + 1.5) / tgamma(kappa + 1.0) * pow(1 + pow(v, 2) / (2 * T * kappa), -kappa - 1.5);
            fv = sqrt(1.0 / (2 * M_PI * T * (kappa-1.5))) * tgamma(kappa) / tgamma(kappa - 0.5) * pow(1 + pow(v, 2) / (2 * T * (kappa - 1.5)), -kappa);
        return fx * fv;
    }
    double GetIonInitDistrib(double x, double v, double t)
    {
        double fx = 1.0;
        return fx;
    }
};
#endif
