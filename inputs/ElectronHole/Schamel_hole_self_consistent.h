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
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/QR>
#include<eigen3/Eigen/Dense>

using namespace std;
using namespace Eigen;

class Input
{
  protected:
    //title
    const string title = "self-consistent Schamel electron hole";

    //general parameters
    const double L = 60;
    const double k = 2.0 * M_PI / L; //simulaiton length
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
    const double b = -2.0;
    const double d = 0.277835;//initial disturbance
    //const double alpha = -25.0;
    const double del = 4.5589;
    const double Z = 1.0;

    //simulation constant
    static const int nx = 1000;//grid num is nx-1; grid point num is nx
    static const int nx_grids = nx - 1;
    static const int nv = 1000;
    static const int nv_grids = nv - 1;
    const double dx = L / nx_grids;
    const double dv = 2 * vmax / nv_grids;
    const double dt = 0.02;
    const int max_steps = 5000;
    const double dt_max = min(dx / vmax, dv * m * k / abs(e * d));

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
            phi_sc(i) = 2.5 * d * pow(cosh(xp), -4);
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
            VectorXd n(nx);
            for (int j = 0; j < nx; j++)
            {
                double temp = 0.0;
                for(int k = 0; k < nv; k++)
                {
                    temp += SchamelDistribution(phi_sc(j), -vmax + k * dv);
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

            double err = ((dphi.array() / phi_sc.segment(1, nx - 2).array())*(dphi.array() / phi_sc.segment(1, nx - 2).array())).sum();
            if(err < 1e-5)
            {
                cout << "Self consistent potential is solved! Iteration counts = "  << i << endl;
                break;
            }
        }
    }

    double SchamelDistribution(double phi, double v)
    {
        double v_waveframe = v;
        double w = v*v/2.0 - phi;
        double r = 0.0;
        if (v_waveframe <= -sqrt(2.0 * phi))
            r = exp(- 0.5 * pow(-sqrt(2*w) + u, 2));
        else if (v_waveframe > sqrt(2.0 * phi))
            r = exp(- 0.5 * pow( sqrt(2*w) + u, 2));
        else
            r = exp(-b * w - u * u / 2.0);
        return r / sqrt(2.0 * M_PI);
    }
    //double Sech2Distrib(double x, double v)
    //{
        //double xp = (x - L / 2.0) / del;
        //double ph = d * pow(cosh(xp), -2);
        //double v_waveframe = v;
        //double w = pow(v_waveframe, 2) - 2.0 * ph;
        //double r = 0.0;
        //if (w >= 0)
            //r = exp(- 0.5 * w);
        //else
            //r = 1.0 - w + alpha * w * w;
        //return r / sqrt(2.0 * M_PI);
    //}

    double GetElecInitDistrib(double x, double v)
    {
        return SchamelDistribution(phi_sc(x / dx), v);
    }

    double GetIonInitDistrib(double x, double v, double t)
    {
        double r = Z;
        return r;
    }

    double E_External(double x, double t)
    {
        //double r = 0.0;
        //if (if_E_External)
        //{
            //double w = 1.0;
            //double tau1 = 3 * dtau;
            //double tau2 = 8 * dtau;
            //double g1 = 1.0 / ( 1 + exp(-5 * (t - tau1) / dtau));
            //double g2 = 1.0 / ( 1 + exp( 5 * (t - tau2) / dtau));
        //}
        return 0.0;
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
        cout << "************************************" << endl;
    }
};
#endif
