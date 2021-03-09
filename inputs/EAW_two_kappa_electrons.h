//To Use this input by changing the filename to input.h
/***********************************
    electron acoustic wave: two-temperature electrons
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
    const string title = "electron acoustic wave: two-kappa electrons";

    //general parameters
    const double k = 0.6;
    const double L = 2 * M_PI / k; //simulaiton length
    const double T = 5.0; //average temperature for all electrons
    const double m = 1.0;
    const double vmax = 20;
    const double e = -1.0;
    const double n = 1.0;
    const double w_p = sqrt(n*e*e / m);
    const double l_D = sqrt(T / n*e*e);
    const bool if_E_External = false;

    //special parameters
    const double ns = 0.79;
    const double nf = n - ns;
    const double kappa_s = 1.51;
    const double kappa_f = 100.0;
    const double l_D_c = sqrt(T / ns);
    const double l_D_h = sqrt(T / nf);
    const double d = 1e-3;

    //simulation constant
    static const int nx = 401;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 1001;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 10000;
    const double dt_max = dv * m * k / abs(e * d);


    //data recording
    const string data_path = "./data/";
    const int data_steps = 2000;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double rx = 1.0 + d * cos(k * x) ;
        double As = ns / sqrt(2 * M_PI * T * (kappa_s - 1.5)) * tgamma(kappa_s) / tgamma(kappa_s - 0.5);
        double Af = nf / sqrt(2 * M_PI * T * (kappa_f - 1.5)) * tgamma(kappa_f) / tgamma(kappa_f - 0.5);
        double rvs = pow(1 + v * v / (2 * T * (kappa_s - 1.5)), -kappa_s);
        double rvf = pow(1 + v * v / (2 * T * (kappa_f - 1.5)), -kappa_f);

        return rx * (As * rvs + Af * rvf);
    }

    double GetIonInitDistrib(double x, double v, double t)
    {
        double r = 1.0;
        return r;
    }

    double E_External(double x, double t)
    {
        double r = 0.0;
        //if (if_E_External)
        //{
        //double w = 1.009;
        //double tau1 = 3 * dtau;
        //double tau2 = 8 * dtau;
        //double g1 = 1.0 / ( 1 + exp(-5 * (t - tau1) / dtau));
        //double g2 = 1.0 / ( 1 + exp( 5 * (t - tau2) / dtau));
        //r = 0.01 * (g1 + g2 - 1.0) * cos(k * x + w * t);
        //}
        return r;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "       nc = " << setw(8) << ns
             << "       nh = " << setw(8) << nf << endl;
        cout << "  kappa_c = " << setw(8) << kappa_s
             << "  kappa_h = " << setw(8) << kappa_f << endl;
        cout << "    l_D_c = " << setw(8) << l_D_c
             << "    l_D_h = " << setw(8) << l_D_h << endl;
        cout << "  k*l_D_c = " << setw(8) << k*l_D_c
             << "  k*l_D_h = " << setw(8) << k*l_D_h << endl;
        cout << "        d = " << setw(8) << d << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(8) << dt_max << endl;
        cout << "************************************" << endl;
    }
};
#endif
