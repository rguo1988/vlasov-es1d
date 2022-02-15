//To Use this input by changing the filename to input.h
/***********************************
 ***********************************/
#ifndef _input_h
#define _input_h
#include<cmath>
#include<iostream>
#include<iomanip>
#include<string>
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/QR>
#include<eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class Input
{
  protected:
    //title
    const string title = "linear ion acoustic waves";

    //general parameters
    const double k = 0.5; //simulaiton length
    const double L = 2.0 * M_PI / k;
    const double vmax = 5;
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
    const double d = 0.01;//initial disturbance

    //simulation constant
    static const int nx = 200;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 500;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 20000;
    const double dt_max = min(dx / vmax, dv * me * k / abs(e * d));

    //data recording
    const string data_path = "./data/";
    const int data_steps = max_steps;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        double rv = exp(-0.5 * pow(v / vt_e, 2)) / sqrt(2.0 * M_PI) / vt_e;
        double rx = 1.0 + d * cos(k * x);
        return rx*rv;
    }

    double GetIonInitDistrib(double x, double v)
    {
        double rv = exp(-0.5 * pow(v / vt_i, 2)) / sqrt(2.0 * M_PI) / vt_i;
        //double rx = 1.0 + d * cos(k * x);
        return rv;
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
