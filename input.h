//To Use this input by changing the filename to input.h
/***********************************
 * formation of kappa stationary distribution in an inhomogeneous plasma
 ***********************************/
#ifndef _input_h
#define _input_h
#include<cmath>
#include<iostream>
#include<iomanip>
#include<string>
using namespace std;
class Input
{
  protected:
    //general parameters
    const double L = 80.0; //simulaiton length
    const double k = 1 * 2.0 * M_PI / L;
    const double T = 1.0; //temperature
    const double m = 1.0;
    const double vmax = 20;
    const double e = -1.0;

    //definition of simulation constant
    static const int nx = 201;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 2001;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.02;
    const int max_steps = 100000;

    //special parameters
    const double uae = 0.52;
    const double uai = 0.5;
    const double kappa = 0.0;
    const double u = 0.0;
    //data recording
    const string data_path = "./data/";
    const int data_steps = 10000;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double rx = 1.0 + uae * cos(k * x) ;
        //double rx = 1.0;
        //double rv = sqrt(1.0 / (2 * M_PI * T * kappa)) * tgamma(kappa + 1.5) / tgamma(kappa + 1.0) * pow(1 + pow(v, 2) / (2 * T * kappa), -kappa-1.5);
        double rv = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v - u, 2) / (2 * T));
        return rx * rv;
    }
    double GetIonInitDistrib(double x, double v, double t)
    {
        double r = 1.0;
        r = 1.0 + uai * cos(k * x);
        return r;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << "Special Parameters: " << endl;
        cout << " kappa = " << setw(8) << kappa
             << "   uae = " << setw(8) << uae
             << "   uai = " << setw(8) << uai << endl;
        cout << "************************************" << endl;
    }
};
#endif
