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
#include"diagnose.h"

using namespace std;
using namespace Eigen;

class Input
{
  protected:
    //title
    const string title = "electron trapping in ion acoustic soliton";

    //general parameters
    const double L = 60;
    const double k = 2.0 * M_PI / L; //simulaiton length
    const double vmax = 5.0;
    const double e = -1.0;
    const double n = 1.0;
#define _ions_motion true
#define _boundary_condition Debye //available bc: Dirichlet, Debye, Periodic

    const double Te = 1.0; //temperature
    const double me = 1.0;
    const double vt_e = sqrt(Te / me);
    const double w_pe = sqrt(n*e*e / me);
    const double l_e = sqrt(Te / n*e*e); //Debye length

    const double Ti = 0.1; //temperature
    const double mi = 100.0;
    const double vt_i = sqrt(Ti / mi);
    const double w_pi = sqrt(n*e*e / mi);
    const double l_i = sqrt(Ti / n*e*e);

    //special parameters
    const double cs = sqrt(Te / mi);
    const double u = 1.1; //dimensionless wave speed; unit in cs
    const double psi = 3 * (u - 1.0) * Te;
    const double del = sqrt(2.0 / (u - 1.0));

    //simulation constant
    static const int nx = 500;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 500;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.05;
    const int max_steps = 2000;
    const double dt_max = min(dx / vmax, dv * me * k / abs(e * psi));

    //data recording
    const string data_path = "./data/";
    const int data_steps = max_steps;
    const int data_num = max_steps / data_steps + 1;

    VectorXd phi_sc;
    VectorXd CalculateRho(VectorXd phi)
    {
        VectorXd rho(nx);
        for (int j = 0; j < nx; j++)
        {
            rho[j] = - exp(phi_sc[j]) + u / sqrt(u * u - 2 * phi_sc[j]);
        }
        return rho;
    }
    void CalculatePotentialSC()
    {
        //dimensionless phi; unit in Te/e; using nonlinear iteration method
        //initial guess
        phi_sc.resize(nx);
        for(int i = 0; i < nx; i++)
        {
            double xp = (i * dx - L / 2.0) / del;
            phi_sc(i) = 0.7 * psi * pow(cosh(xp), -2);
        }
        //phi_sc(0) = psi * pow(cosh(-L / 2.0 / del), -2);
        phi_sc(0) = 0.0;
        phi_sc(nx - 1) = phi_sc(0);

        //construct laplace op
        //Dirichlet BC
        MatrixXd laplace(nx - 2, nx - 2);
        laplace.diagonal() = VectorXd::Constant(nx - 2, -2);
        laplace.diagonal(1) =  VectorXd::Constant(nx - 3, 1);
        laplace.diagonal(-1) =  VectorXd::Constant(nx - 3, 1);
        double dx2 = dx * dx;

        //solve by iteration
        for(int i = 0; i < 10; i++)
        {
            VectorXd rho(nx);
            rho = CalculateRho(phi_sc);
            VectorXd r = rho.segment(1, nx - 2) * dx2 + laplace * phi_sc.segment(1, nx - 2);
            VectorXd m = VectorXd::Ones(nx - 2).array() * 2.0 - (rho.segment(2, nx - 2) - rho.segment(0, nx - 2)).array() / (phi_sc.segment(2, nx - 2) - phi_sc.segment(0, nx - 2)).array() * dx2;
            MatrixXd nr_mat(nx - 2, nx - 2);
            nr_mat.setZero();
            nr_mat.diagonal() = m;
            nr_mat.diagonal(1) = VectorXd::Constant(nx - 3, -1);
            nr_mat.diagonal(-1) = VectorXd::Constant(nx - 3, -1);
            VectorXd dphi = nr_mat.fullPivLu().solve(r);
            phi_sc.segment(1, nx - 2) += dphi;

            double err = sqrt((dphi.array() * dphi.array()).sum() / nx);
            if(err < psi * 1e-4)
            {
                cout << "Self consistent potential is solved! Iteration counts = "  << i << endl;
                break;
            }
        }
        string file = data_path + "phi_sc";
        OutputMatrix(file, phi_sc);
    }

    double GetElecInitDistrib(double x, double v)
    {
        //double r = sqrt(me / (2 * M_PI * Te)) * exp(- me * pow(v - u, 2) / (2 * Te));
        double r = sqrt(me / 2 / M_PI / Te) * exp(- me * pow(v + u * cs, 2) / (2 * Te)) * exp(phi_sc[x / dx] * Te);
        //double phi = phi_sc[x / dx];
        //double w = v * v / 2.0 - phi;
        //double r = 0.0;

        //if (v <= -sqrt(2.0 * phi))
        //r = exp(- pow(-sqrt(2 * w) + u * cs, 2) / 2 / Te);
        //else if (v > sqrt(2.0 * phi))
        //r = exp(- pow( sqrt(2 * w) + u * cs, 2) / 2 / Te);
        //else
        //r = 1.0;
        //return r * sqrt(me / 2.0 / M_PI / Te);
        return r;
    }

    double GetElecFreeDistrib(double x, double v)
    {
        double r = sqrt(mi / (2 * M_PI * Te)) * exp(- mi * pow(v + u * cs, 2) / (2 * Te));
        return r;
    }

    double GetIonInitDistrib(double x, double v)
    {
        double rx = n * u / sqrt(u * u - 2.0 * phi_sc[x / dx]);
        double rv = sqrt(mi / (2 * M_PI * Ti)) * exp(- mi * pow(v + u * cs, 2) / (2 * Ti));
        return rx * rv;
    }

    double GetIonFreeDistrib(double x, double v)
    {
        return GetIonInitDistrib(0.0, v);
    }

    void PrintSpecialParameters()
    {
        cout << "************************************" << endl;
        cout << " Special Parameters: " << endl;
        cout << "          u = " << setw(6) << u << endl;
        cout << " psi_approx = " << setw(6) << psi
             << " del_approx = " << setw(6) << del << endl;
        cout << "************************************" << endl;
        cout << " Parameters Max/Min: " << endl;
        cout << "   dt_max = " << setw(6) << dt_max << endl;
        cout << "************************************" << endl;
    }
};
#endif
