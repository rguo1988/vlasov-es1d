//To Use this input by changing the filename to input.h
/***********************************
 * formation of Schamel electron hole
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
    //title
    const string title = "formation of electron hole";

    //general parameters
    const double k = 0.1;
    const double L = 2.0 * M_PI / k; //simulaiton length
    const double T = 1.0; //temperature
    const double m = 1.0;
    const double vmax = 5;
    const double e = -1.0;
    const double n = 1.0;
    const double w_p = sqrt(n*e*e / m);
    const double l_D = sqrt(T / n*e*e);
    const bool if_E_External = false;

    //special parameters
    const double d = 0.2;//initial disturbance
    const double del = 5.0;

    //simulation constant
    static const int nx = 500;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 500;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 30000;
    const double dt_max = min(dx / vmax, dv * m * k / abs(e * d));


    //data recording
    const string data_path = "./data/";
    const int data_steps = max_steps;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        double rv = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v, 2) / (2 * T));
        double xp = (x - L / 2) / del;
        //double rx = 1 + d * sin(k * x);
        double rx = 1.0 - d * exp(-xp * xp);
        return rx * rv;
    }

    double GetIonInitDistrib(double x, double v, double t)
    {
        double r = 1.0;
        return r;
    }

    double E_External(double x, double t)
    {
        double r = 0.0;
        if (if_E_External)
        {
            //double w = 1.0;
            //double tau1 = 3 * dtau;
            //double tau2 = 8 * dtau;
            //double g1 = 1.0 / ( 1 + exp(-5 * (t - tau1) / dtau));
            //double g2 = 1.0 / ( 1 + exp( 5 * (t - tau2) / dtau));
        }
        return r;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "        d = " << setw(6) << d << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(6) << dt_max << endl;
        cout << "************************************" << endl;
    }
};
#endif
