//To Use this input by changing the filename to input.h
/***********************************
    electron acoustic instability: two-temperature electrons
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
    const string title = "electron acoustic instability: two-kappa electrons";

    //general parameters
    const double k = 0.6;
    const double L = 2 * M_PI / k; //simulaiton length
    const double vmax = 150;
    const double e = -1.0;
    const double n = 1.0;
#define _ions_motion false

    const double me = 1.0;
    const double Te = 1.0; //average temperature for all electrons
    const double w_pe = sqrt(n*e*e / me);
    const double l_e = sqrt(Te / n*e*e);

    //special parameters
    const double ns = 0.5;
    const double nf = n - ns;
    const double v_s = 0.1;
    const double v_f = 1.0;
    //const double kappa_s = 2.0;
    //const double kappa_f = 2.0;
    //const double Ts = 0.5 * me * kappa_s / (kappa_s - 1.5) * v_s * v_s;
    //const double Tf = 0.5 * me * kappa_f / (kappa_f - 1.5) * v_f * v_f;
    const double Ts = 0.5 * me * v_s * v_s;
    const double Tf = 0.5 * me * v_f * v_f;
    const double l_s = sqrt(Ts / ns);
    const double l_f = sqrt(Tf / nf);
    const double d = 1e-6;
    const double a = 10;
    const double vt_e = v_f;
    const double u_f = a * v_s;
    const double u_s = 0.0;

    //simulation constant
    static const int nx = 501;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 7501;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.1;
    const int max_steps = 2500;
    const double dt_max = dv * me * k / abs(e * d);


    //data recording
    const string data_path = "./data/";
    const int data_steps = max_steps;
    const int data_num = max_steps / data_steps + 1;

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double rx = 1.0 + d * cos(k * x) ;

        //double As = ns / sqrt(M_PI * kappa_s) / v_s * tgamma(kappa_s) / tgamma(kappa_s - 0.5);
        //double rvs = pow(1 + (v - u_s) * (v - u_s) / kappa_s / v_s / v_s, -kappa_s);

        double As = ns / sqrt(M_PI) / v_s;
        double rvs = exp(-(v - u_s) * (v - u_s) / v_s / v_s);

        double Af = nf / sqrt(M_PI) / v_f;
        double rvf = exp(-(v - u_f) * (v - u_f) / v_f / v_f);

        //double Af = nf / sqrt(M_PI * kappa_f) / v_f * tgamma(kappa_f) / tgamma(kappa_f - 0.5);
        //double rvf = pow(1 + (v - u_f) * (v - u_f) / kappa_f / v_f / v_f, -kappa_f);

        return rx * (As * rvs + Af * rvf);
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
        cout << "       ns = " << setw(8) << ns
             << "       nf = " << setw(8) << nf << endl;
        //cout << "  kappa_s = " << setw(8) << kappa_s
        cout << "  kappa_s = " << setw(8) << "inf"
             << "  kappa_f = " << setw(8) << "inf" << endl;
             //<< "  kappa_f = " << setw(8) << kappa_f << endl;
        cout << "      v_s = " << setw(8) << v_s
             << "      v_f = " << setw(8) << v_f << endl;
        cout << "      u_s = " << setw(8) << u_s
             << "      u_f = " << setw(8) << u_f
             << "        a = " << setw(8) << a << endl;
        cout << "    k*l_s = " << setw(8) << k*l_s
             << "    k*l_f = " << setw(8) << k*l_f << endl;
        cout << "        d = " << setw(8) << d << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(8) << dt_max << endl;
        cout << "************************************" << endl;
    }
};
#endif
