//To Use this input by changing the filename to input.h
/***********************************
 * simulation of $Pn^{-3}$ in Maxwellian plasma
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
    const double k = 0.1;
    const double L = 2.0 * M_PI / k; //simulaiton length
    const double T = 1.0; //temperature
    const double m = 1.0;
    const double vmax = 5;
    const double e = -1.0;

    //definition of simulation constant
    static const int nx = 301;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 301;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.02;
    const int max_steps = 5000;

    //special parameters
    const double delta = 0.003;
    const double kappa = 0.0;
    //data recording
    const string data_path = "./data/";
    const int data_steps = 10;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        double rx = 1.0 + delta * cos(k * x) ;

        //non-adiabatic initial condition
        double rv = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v, 2) / (2 * T));
        return rx * rv;
        
        //adiabatic initial condition
        //double rv = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v/rx, 2) / (2 * T));
        //return 1.0 * rv;
    }
    double GetIonInitDistrib(double x, double v, double t)
    {
        double r = 1.0;
        return r;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << "Special Parameters: " << endl;
        cout << " kappa = " << setw(8) << kappa
             << "   delta = " << setw(8) << delta << endl;
        cout << "************************************" << endl;
    }
};
#endif
