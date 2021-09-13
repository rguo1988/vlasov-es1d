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
    //title
    const string title = "formation of kappa stationary distribution";

    //general parameters
    const double L = 80.0; //simulaiton length
    const double k = 1 * 2.0 * M_PI / L;
    const double T = 1.0; //temperature
    const double m = 1.0;
    const double vmax = 10;
    const double e = -1.0;
    const double n = 1.0;
    const double w_p = sqrt(n*e*e / m);
    const double l_D = sqrt(T / n*e*e);
    const bool if_E_External = false;

    //special parameters
    const double uae = 0.8;
    const double uai = 0.8;
    const double kappa = 0.0;
    const double u = 0.5;
    const double dtau = 100;
    const double tau = 3 * dtau;

    //simulation constant
    static const int nx = 201;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 1001;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.02;
    const int max_steps = 100000;
    const double dt_max = dv * m * k / abs(e * (uae - uai));


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
        //double g = 1.0 / (1 + exp(-5 * (t - tau) / dtau));
        double r = 1.0 + uai * cos(k * x);
        return r;
    }

    double E_External(double x, double t)
    {
        double r = 0.0;
        if (if_E_External)
        {
            double w = 1.009;
            double tau1 = 3 * dtau;
            double tau2 = 8 * dtau;
            double g1 = 1.0 / ( 1 + exp(-5 * (t - tau1) / dtau));
            double g2 = 1.0 / ( 1 + exp( 5 * (t - tau2) / dtau));
            r = 0.01 * (g1 + g2 - 1.0) * cos(k * x + w * t);
        }
        return r;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "    kappa = " << setw(6) << kappa
             << "      uae = " << setw(6) << uae
             << "      uai = " << setw(6) << uai << endl;
        cout << "        u = " << setw(6) << u << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(6) << dt_max << endl;
        cout << "************************************" << endl;
    }
};
#endif
