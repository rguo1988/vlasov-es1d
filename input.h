//To Use this input by changing the filename to input.h
/***********************************
 * waves in polytropic plasma
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
    const string title = "waves in polytropic plasma";
    //general parameters
    const double k = 0.05;
    const double L = 2 * M_PI / k; //simulaiton length
    const double T = 1.0; //temperature
    const double m = 1.0;
    const double vmax = 5;
    const double e = -1.0;
    const double lambda = sqrt(T / e / e);

    //definition of simulation constant
    static const int nx = 1001;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 201;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.01;
    const int max_steps = 10000;

    //special parameters
    //d must smaller than k^2*kappa
    const double d = 0.001;
    const double d_wave = 0.1;
    const double k_wave = 0.1;
    const double kappa = 5.0;
    //data recording
    const string data_path = "./data/";
    const int data_steps = 10000;
    const int data_num = max_steps / data_steps + 1;

    //calculating the normalization
    double B = d / (k * k * lambda * lambda * kappa );
    double A = L * 1.0 / GetNormalization();

    double GetElecInitDistrib(double x, double v)
    {
        //f is distribution function normalized to 1, i.e. n=1
        double fx = A * pow(1.0 - B * cos(k * x), -kappa - 1);
        //double fx = 1.0;
        double wave = 1.0 + d_wave * cos(k_wave * x);
        double Tx = T * pow(fx, -1.0 / (kappa + 1.0));
        //double Tx = T;
        double fv = sqrt(1.0 / (2 * M_PI * Tx * kappa)) * tgamma(kappa + 1.5) / tgamma(kappa + 1.0) * pow(1 + pow(v, 2) / (2 * Tx * kappa), -kappa - 1.5);
        //double fv = sqrt(1.0 / (2 * M_PI * T)) * exp(- pow(v, 2) / (2 * T));
        return wave * fx * fv;
    }
    double GetIonInitDistrib(double x, double v, double t)
    {
        double fx = A * pow(1.0 - B * cos(k * x), -kappa - 1) + d * cos(k * x);
        return fx;
    }
    double GetNormalization()
    {
        double sum = 0.0;
        for (int i = 0; i < nx_grids; i++)
        {
            sum += pow(1.0 - B * cos(k * i * dx), -kappa - 1);
        }
        return sum * dx;
    }
    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "   kappa = " << setw(8) << kappa
             << "       d = " << setw(8) << d << endl;
        cout << "  k_wave = " << setw(8) << k_wave
             << "  d_wave = " << setw(8) << d_wave << endl;
        cout << "************************************" << endl;
    }
};
#endif
