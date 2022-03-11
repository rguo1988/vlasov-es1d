//To Use this input by changing the filename to input.h
/***********************************
 * formation of electron hole due to delta disturbance
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
    const string title = "formation of electron hole due to disturbance";

    //general parameters
    const double L = 80.0;
    const double k = 2.0 * M_PI / L; //simulaiton length
    const double vmax = 5.0;
    const double e = -1.0;
    const double n = 1.0;
#define _ions_motion false

    const double Te = 1.0; //temperature
    const double me = 1.0;
    const double vt_e = sqrt(Te / me);
    const double w_pe = sqrt(n*e*e / me);
    const double l_e = sqrt(Te / n*e*e); //Debye length

    const double Ti = 1.0; //temperature
    const double mi = 100.0;
    const double vt_i = sqrt(Ti / mi);
    const double w_pi = sqrt(n*e*e / mi);
    const double l_i = sqrt(Ti / n*e*e);

    //special parameters
    const double d = 0.05;//initial disturbance
    const double u = 0.0;
    const double delta_x = 1.0;
    const double delta_v = 0.1;

    //simulation constant
    static const int nx = 500;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 500;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.02;
    const int max_steps = 10000;
    const double dt_max = min(dx / vmax, dv * me * k / abs(e * d));


    //data recording
    const string data_path = "./data/";
    const int data_steps = 1000;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        double rv = exp(-0.5 * pow((v + u) / vt_e, 2)) / sqrt(2 * M_PI) / vt_e;
        double xp = (x - L / 2) / delta_x;
        double vp = v / delta_v;
        double rv1_x = (20 * pow(cosh(xp), -6) - 16 * pow(cosh(xp), -4)) / delta_x / delta_x;
        //double rv1_x = (6 * pow(cosh(xp), -4) - 4 * pow(cosh(xp), -2)) / delta_x / delta_x;
        double rv1 = -d / cosh(vp) * rv1_x;
        return rv + rv1;
    }

    double GetIonInitDistrib(double x, double v, double t)
    {
        double r = 1.0;
        return r;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "        d = " << setw(6) << d
             << "  delta_x = " << setw(6) << delta_x
             << "  delta_v = " << setw(6) << delta_v << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(6) << dt_max << endl;
        cout << "************************************" << endl;
    }
};
#endif
