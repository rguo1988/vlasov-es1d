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
    const string title = "formation of electron hole in two stream plasmas";

    //general parameters
    //const double k = 0.5;
    //const double L = 2 * M_PI / k; //simulaiton length
    const double L = 100; //simulaiton length
    const double k = 2 * M_PI / L;
    const double Te = 1.0; //temperature
    const double me = 1.0;
    const double vt_e = sqrt(Te / me);
    const double vmax = 10;//10.0 * sqrt(temperature / m);
    const double e = -1.0;
    const double n = 1.0;
    const double w_pe = sqrt(n*e*e / me);
    const double l_e = sqrt(Te / n*e*e);

    //special parameters
    const double u1 = 1.2; //drift speed of stream 1
    const double u2 = -1.2; //drift speed of stream 2
    const double d = 0.5;
    const double hx = 5.0; //hole length in x
    const double hv = 0.1; //hole length in v

    //definition of simulation constant
    static const int nx = 2000;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 1000;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.05;
    const int max_steps = 100;
    const double dt_max = min(dx / vmax, dv * me * k / abs(e * d));
    //
    //data saveing
    const string data_path = "./data/";
    const int data_steps = max_steps;
    const int data_num = max_steps / data_steps + 1;

    double A = CalNorm();

    double CalNorm()
    {
        double sum_f = 0.0;
        for (int i = 0; i < nx; i++)
            for (int j = 0; j < nv; j++)
                sum_f += ElecDistrib(i * dx, -vmax + j * dv);
        return sum_f * dx * dv / L;
    }
    double ElecDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double xp = x - L / 2.0;
        double p = d * pow(cosh(xp / hx), -2);
        double r1 = sqrt(1.0 / (2 * M_PI * Te)) * exp(- pow(v - u1, 2) / (2 * Te));
        double r2 = sqrt(1.0 / (2 * M_PI * Te)) * exp(- pow(v - u2, 2) / (2 * Te));
        //double fd = d * exp(-v * v / 2 / hv) * exp(-xp * xp / 2 / hx);
        //double fd = 1.0 + 4.0 * p / hx / hx - 6.0 * p * p / d / hx / hx;
        double fd = d * exp(-v * v / 2 / hv) * (-4.0 * p / hx / hx + 6.0 * p * p / d / hx / hx);

        return 0.5 * (r1 + r2) * (1 - fd);

        //return 0.5 * (r1 + r2) * (1 + d * cos(k * x));
        //double r = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v, 2) / (2 * T));
        //return r;
    }

    double GetElecInitDistrib(double x, double v)
    {
        return ElecDistrib(x, v);
        //return ElecDistrib(x, v) / A;
    }
    double GetElecFreeDistrib(double x, double v)
    {
        double r1 = sqrt(1.0 / (2 * M_PI * Te)) * exp(- pow(v - u1, 2) / (2 * Te));
        double r2 = sqrt(1.0 / (2 * M_PI * Te)) * exp(- pow(v - u2, 2) / (2 * Te));

        return 0.5 * (r1 + r2);
    }
    double GetIonInitDistrib(double x, double v, double t)
    {
        return 1.0;
    }
    double GetIonFreeDistrib(double x, double v, double t)
    {
        return 1.0;
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
