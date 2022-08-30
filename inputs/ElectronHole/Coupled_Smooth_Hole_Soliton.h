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
    const string title = "coupled smooth electron hole soliton";

    //general parameters
    const double L = 100;
    const double k = 2.0 * M_PI / L; //simulaiton length
    const double vmax = 10.0;
    const double e = -1.0;
    const double n = 1.0;
#define _ions_motion true

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
    const double u = 1.24;
    const double b = -1.28;
    const double psi = 0.3475;//initial disturbance
    const double del = 8.655;

    //simulation constant
    static const int nx = 500;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 500;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.02;
    const int max_steps = 1000;
    const double dt_max = min(dx / vmax, dv * me * k / abs(e * psi));

    VectorXd phi_sc;

    //data recording
    const string data_path = "./data/";
    const int data_steps = max_steps;
    const int data_num = max_steps / data_steps + 1;

    void CalculatePotentialSC()
    {
        phi_sc.resize(nx);
        for(int i = 0; i < nx; i++)
        {
            double xp = (i * dx - L / 2.0) / del;
            phi_sc(i) = 0.7 * psi * pow(cosh(xp), -2);// smooth EH
            //phi_sc(i) = 0.7 * psi * pow(cosh(xp), -4);  // Schamel sech4 EH
        }
        phi_sc(0) = 0.0;
        phi_sc(nx - 1) = 0.0;
        //construct laplace op
        //Dirichlet BC
        MatrixXd laplace(nx - 2, nx - 2);
        laplace.diagonal() = VectorXd::Constant(nx - 2, -2);
        laplace.diagonal(1) =  VectorXd::Constant(nx - 3, 1);
        laplace.diagonal(-1) =  VectorXd::Constant(nx - 3, 1);
        double dx2 = dx * dx;

        for(int i = 0; i < 10; i++)
        {
            VectorXd n;
            n.resize(nx);
            n.setZero();
            for (int j = 0; j < nx; j++)
            {
                double temp = 0.0;
                for(int k = 0; k < nv; k++)
                {
                    temp += SmoothDistribution(phi_sc(j), -vmax + k * dv); // smooth EH
                    //temp += SchamelDistribution(phi_sc(j), -vmax + k * dv);// Schamel sech4 EH
                }
                n(j) = temp * dv;
            }
            VectorXd r = (VectorXd::Ones(nx - 2) - n.segment(1, nx - 2)) * dx2 + laplace * phi_sc.segment(1, nx - 2);
            VectorXd m = VectorXd::Ones(nx - 2).array() * 2.0 + (n.segment(2, nx - 2) - n.segment(0, nx - 2)).array() / (phi_sc.segment(2, nx - 2) - phi_sc.segment(0, nx - 2)).array() * dx2;
            MatrixXd nr_mat(nx - 2, nx - 2);
            nr_mat.setZero();
            nr_mat.diagonal() = m;
            nr_mat.diagonal(1) = VectorXd::Constant(nx - 3, -1);
            nr_mat.diagonal(-1) = VectorXd::Constant(nx - 3, -1);
            VectorXd dphi = nr_mat.fullPivLu().solve(r);
            phi_sc.segment(1, nx - 2) += dphi;

            double err = ((dphi.array() / phi_sc.segment(1, nx - 2).array()) * (dphi.array() / phi_sc.segment(1, nx - 2).array())).sum();
            if(err < 1e-5)
            {
                cout << "Self consistent potential is solved! Iteration counts = "  << i << endl;
                break;
            }
        }
    }

    //smooth EH
    double SmoothDistribution(double phi, double v)
    {
        double w = v * v / 2.0 - phi;
        double r = 0.0;
        if (w > 0)
            r = ( exp(- 0.5 * pow(sqrt(2.0 * w) + u, 2))
                + exp(- 0.5 * pow(sqrt(2.0 * w) - u, 2)) ) / 2.0;
        else
            r = exp(-w - u * u / 2.0) * cos(sqrt(-2.0 * w) * u); //smooth EH
            //r = exp(-b * w - u * u / 2.0); //Schamel EH
        return r / sqrt(2.0 * M_PI);
    }
    //
    //Schamel sech4 EH
    //double SchamelDistribution(double phi, double v)
    //{
        //double v_waveframe = v;
        //double w = v * v / 2.0 - phi;
        //double r = 0.0;
        //if (v_waveframe <= -sqrt(2.0 * phi))
            //r = exp(- 0.5 * pow(-sqrt(2 * w) + u, 2));
        //else if (v_waveframe > sqrt(2.0 * phi))
            //r = exp(- 0.5 * pow( sqrt(2 * w) + u, 2));
        //else
            //r = exp(-b * w - u * u / 2.0);
        //return r / sqrt(2.0 * M_PI);
    //}

    double GetElecInitDistrib(double x, double v)
    {
        //return SchamelDistribution(phi_sc(x / dx), v); //Schamel sech4 EH
        return SmoothDistribution(phi_sc(x / dx), v);//smooth EH
    }

    double GetIonInitDistrib(double x, double v)
    {
        double vt = 0.2;
        double r = 0.0;
        //if (v >= -vt - u && v <= vt - u)
        //{
        //r = 0.5 / vt;
        //}
        r = exp(-0.5 * pow((v + u) / vt_i, 2)) / sqrt(2.0 * M_PI) / vt_i;
        return r;
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "        u = " << setw(6) << u << endl;
        cout << "      psi = " << setw(6) << psi
             << "      del = " << setw(6) << del << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(6) << dt_max << endl;
        cout << "************************************" << endl;
    }
};
#endif
