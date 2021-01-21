//To Use this input by changing the filename to input.h
/***********************************
 * formation of 4 hole BGK wave
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
    const double k = 0.452;
    const double L = 4 * 2.0 * M_PI / k; //simulaiton length
    const double vth = 1.0; //temperature
    const double T = (5 - 2 * xi) / (3 - 2 * xi) * vth * vth;
    const double m = 1.0;
    const double vmax = 4;
    const double e = -1.0;

    //definition of simulation constant
    static const int nx = 1001;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 501;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 10000;

    //special parameters
    const double xi = 0.95;
    const double delta = 0.1;
    //data recording
    const string data_path = "./data/";
    const int data_steps = 10;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double rx = 1.0 + delta * cos(2 * k * x) ;
        double W = 0.5 * pow(v / vth, 2);
        double rv = 1.0 / (sqrt(2 * M_PI) * vth) * (2 - 2 * xi) / (3 - 2 * xi) * (1 + W / (1 - xi) ) * exp(-W);
        return rx * rv;
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
        cout << "    xi = " << setw(8) << xi
             << " delta = " << setw(8) << delta << endl;
        cout << "************************************" << endl;
    }
};
#endif
