//To Use this input by changing the filename to input.h
/***********************************
 * electron acoustic wave: plateau-Maxwellian PDF
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
    const string title = "electron acoustic wave: plateau-Maxwellian PDF";

    //general parameters
    const double k = 0.4;
    const double L = 2 * M_PI / k; //simulaiton length
    const double T = 1.0; //T for cold e
    const double m = 1.0;
    const double v_th = sqrt(2 * T / m);
    const double vmax = 5;
    const double e = -1.0;
    const double lambda = sqrt(T / e / e);

    //special parameters
    const double d = 0.001;
    const double v_0 = 1.56036;
    const double dv_p = 0.3;
    const double dv_p_min = sqrt(4 * d / k / k);
    const int n_p = 10;
    const double Z = GetNormalization();

    //simulation constant
    static const int nx = 201;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 2001;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.1;
    const double dt_max = dv * m * k / abs(e) / d;
    const int max_steps = 2000;

    //data recording
    const string data_path = "./data/";
    const int data_steps = 200;
    const int data_num = max_steps / data_steps + 1;


    double f_p_unnorm(double v)
    {
        double fm = exp(-v * v / v_th / v_th) / sqrt(M_PI) / v_th;
        double fm0 = exp(-v_0 * v_0 / v_th / v_th) / sqrt(M_PI) / v_th;
        double D = 0.0;
        if (v >= 0)
            D = 1 + pow( ((v - v_0) / dv_p), n_p);
        else
            D = 1 + pow( ((v + v_0) / dv_p), n_p);
        double N = fm - fm0;
        double r = fm - N / D;
        return r;
    }

    double f_m(double v)
    {
        double fm = exp(-v * v / v_th / v_th) / sqrt(M_PI) / v_th;
        return fm;
    }

    double GetNormalization()
    {
        double n = 0.0;
        int N = 100000;
        double v_lim = 100;
        double dv_lim = 2 * v_lim / N;
        for (int i = 0; i < N; i++)
        {
            n += f_p_unnorm(-v_lim + i * dv_lim) ;
        }
        n *= dv_lim;
        return n;
    }

    double GetElecInitDistrib(double x, double v)
    {
        double wave = 1.0 + d * cos(k * x);
        //return wave * f_p_unnorm(v) / Z;
        return wave * f_m(v);
    }

    double GetIonInitDistrib(double x, double v, double t)
    {
        return 1.0;
    }

    double E_External(double x, double t)
    {
        //double T = 2 * M_PI / k / v_0;
        //double dtau = 10 * T;
        //double tau = 20 * T;
        //double g = 1.0 / (
                       //1 + pow((t - tau) / dtau, 10)
                   //);
        //double r = d * g * sin(k * x - k * v_0 * t);
        return 0.0;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "       Z = " << setw(8) << Z
             << "    v_th = " << setw(8) << v_th
             << "     v_0 = " << setw(8) << v_0 << endl;
        cout << "    dv_p = " << setw(8) << dv_p
             << "     n_p = " << setw(8) << n_p << endl;
        cout << "       d = " << setw(8) << d << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "  dt_max = " << setw(8) << dt_max << endl;
        cout << "dv_p_min = " << setw(8) << dv_p_min << endl;
        cout << "************************************" << endl;
    }
};
#endif
