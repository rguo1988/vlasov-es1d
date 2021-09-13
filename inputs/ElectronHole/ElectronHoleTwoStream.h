//To Use this input by changing the filename to input.h
#ifndef _input_h
#define _input_h
#include<cmath>
#include<string>
#include<iostream>
#include<iomanip>
using namespace std;
class Input
{
  protected:
    //title
    const string title = "formation of electron hole due to two stream instabilities";

    //general parameters
    const double k = 0.2;
    const double L = 2 * M_PI / k; //simulaiton length
    const double T = 0.1; //temperature
    const double m = 1.0;
    const double vmax = 5;//10.0 * sqrt(temperature / m);
    const double e = -1.0;
    const double n = 1.0;
    const double w_p = sqrt(n*e*e / m);
    const double l_D = sqrt(T / n*e*e);
    const bool if_E_External = false;

    //special parameters
    const double u1 = 1.0; //drift speed of stream 1
    const double u2 = -1.0; //drift speed of stream 2
    const double d = 0.01;

    //definition of simulation constant
    static const int nx = 500;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 500;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 70000;
    const double dt_max = min(dx / vmax, dv * m * k / abs(e * d));
    //
    //data saveing
    const string data_path = "./data/";
    const int data_steps = max_steps;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double r1 = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v - u1, 2) / (2 * T));
        double r2 = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v - u2, 2) / (2 * T));
        return 0.5 * (r1 + r2) * (1 + d * cos(k * x));
        //double r = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v, 2) / (2 * T));
        //return r;
    }
    double GetIonInitDistrib(double x, double v, double t)
    {
        return 1.0;
    }
    double E_External(double x, double t)
    {
        double r = 0.0;
        return r;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "       u1 = " << setw(6) << u1
             << "       u2 = " << setw(6) << u2 << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(6) << dt_max << endl;
        cout << "************************************" << endl;
    }
};
#endif
