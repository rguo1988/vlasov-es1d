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
    const double vmax = 10;
    const double e = -1.0;
    const double n = 1.0;
    const double w_p = sqrt(n*e*e / m);
    const double l_D = sqrt(T / n*e*e);
    const bool if_E_External = false;

    //special parameters
    const double u = 0.5;
    const double b = -4.805;
    const double d = 0.08743;//initial disturbance
    const double del = 4.5589;
    //const double Z = GetNormalizationSD();
    const double Z = 1.0;
    //const double Z = GetNormalizationDD();

    //simulation constant
    static const int nx = 1000;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 1000;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.02;
    const int max_steps = 10;
    const double dt_max = min(dx / vmax, dv * m * k / abs(e * d));


    //data recording
    const string data_path = "./data/";
    const int data_steps = max_steps;
    const int data_num = max_steps / data_steps + 1;

    double SchamelDistribution(double x, double v)
    {
        double xp = (x - L / 2.0) / del;
        //double xp = (x) / del;
        double ph = d * pow(cosh(xp), -4);
        //double ph = d * exp(-xp * xp);
        double v_waveframe = v;
        double w = pow(v_waveframe, 2) - 2.0 * ph;
        double r = 0.0;
        if (v_waveframe <= -sqrt(2.0 * ph))
            r = exp(- 0.5 * pow(-sqrt(w) + u, 2));
        else if (v_waveframe > sqrt(2.0 * ph))
            r = exp(- 0.5 * pow( sqrt(w) + u, 2));
        else
            r = exp(-b * w / 2.0 - u * u / 2.0);
        return r / sqrt(2.0 * M_PI);
    }
    double GetNormalizationDD()
    {
        double n = 0.0;
        int N = 10000;
        double x_lim = L;
        double dx_lim = x_lim / N;
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i < N; i++)
            n += DeltaPotentialDistrib(i * dx_lim);
        n *= dx_lim;
        return L / n;
    }
    double GetNormalizationSD()
    {
        double n = 0.0;
        int N = 20000;
        int M = 20000;
        double v_lim = 20;
        double x_lim = L;
        double dv_lim = 2.0 * v_lim / N;
        double dx_lim = x_lim / M;
        #pragma omp parallel for schedule(guided)
        for (int i = 0; i < N; i++)
            for(int j = 0; j < M; j++)
                n += SchamelDistribution(j * dx_lim, -v_lim + i * dv_lim) ;
        n *= dv_lim * dx_lim;
        return L / n;
    }
    double DeltaPotentialDistrib(double x)
    {
        double xp = (x - L / 2) / 4;
        double rh = 0.1 / 16.0 * (16 * pow(cosh(xp), -4) - 20 * pow(cosh(xp), -6));
        double rx = 1.0 + rh;
        return rx;
    }

    double GetElecInitDistrib(double x, double v)
    {
        //double rv = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v - u, 2) / (2 * T));
        //return DeltaPotentialDistrib(x) * rv;
        return Z * SchamelDistribution(x, v);
    }

    double GetIonInitDistrib(double x, double v, double t)
    {
        double r = Z;
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
        cout << "     beta = " << setw(6) << b
             << "        u = " << setw(6) << u << endl;
        cout << "        d = " << setw(6) << d
             << "      del = " << setw(6) << del << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(6) << dt_max << endl;
        cout << "        Z = " << setw(6) << Z << endl;
        cout << "************************************" << endl;
    }
};
#endif
