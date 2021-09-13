#include<iostream>
#include<iomanip>

#include<eigen3/Eigen/Dense>
#include"cubic_spline_1d.h"
#include"input.h"
#include"simulation.h"
#include"diagnose.h"
#include"poisson_solver.h"

using namespace Eigen;
using namespace std;

Simulation::Simulation()
{
    x_samples.setLinSpaced(nx, 0, L);
    v_samples.setLinSpaced(nv, -vmax, vmax);
    //initialize distribution
    f.resize(nx, nv);
    f.setZero();
    for (int i = 0; i < nx; i++)
    {
        for(int j = 0; j < nv; j++)
        {
            f(i, j) = GetElecInitDistrib(x_samples[i], v_samples[j]);
        }
    }
    for (int i = 0; i < nx; i++)
    {
        rho_i[i] =  GetIonInitDistrib(x_samples[i], 0.0, 0.0);
    }
    Ek.clear();
    Ep.clear();
    Et.clear();
}
void Simulation::Run()
{
    double x_shift = 0.0;
    double v_shift = 0.0;
    MatrixXd f_xshift(nx, nv);
    f_xshift.setZero();
    MatrixXd f_vshift(nx, nv);
    f_vshift.setZero();

    //show information
    cout << "************************************" << endl;
    cout << " Vlasov Simulation: " << title << endl;
    cout.setf(ios::left);
    cout << "************************************" << endl;
    cout << " Plasma Parameters: " << endl;
    cout << "        k = " << setw(8) << setprecision(3) << k
         << "        T = " << setw(8) << setprecision(6) << T
         << "      w_p = " << setw(8) << setprecision(6) << w_p
         << "      l_D = " << setw(8) << setprecision(6) << l_D << endl;
    string ifEEx = if_E_External ? "On" : "Off";
    cout << "     E_Ex = " << setw(8) << setprecision(6) << ifEEx << endl;
    cout << "************************************" << endl;
    cout << " Simulation Parameters: " << endl;
    cout << "        L = " << setw(8) << setprecision(6) << L
         << "       nx = " << setw(8) << setprecision(6) << nx
         << "       dx = " << setw(8) << setprecision(6) << dx << endl;
    cout << "     vmax = " << setw(8) << setprecision(6) << vmax
         << "       nv = " << setw(8) << setprecision(6) << nv
         << "       dv = " << setw(8) << setprecision(6) << dv << endl;
    cout << "       wT = " << setw(8) << setprecision(6) << max_steps*dt
         << "    steps = " << setw(8) << setprecision(6) << max_steps
         << "       dt = " << setw(8) << setprecision(6) << dt << endl;
    cout << "  dataNum = " << setw(8) << setprecision(6) << data_num << endl;

    PrintSpecialParameters();

    for(int n = 0; n < max_steps + 1; n++)
    {
        //print running process
        int percent = 100 * n / (max_steps - 1);
        if(percent % 5 == 0)
        {
            cout << "\r" << " Process: " << percent << "%" << flush;
        }
        //diagnose
        if(n % data_steps == 0)
        {
            int nn = n / data_steps;
            string filename = data_path + "data" + to_string(nn);
            OutputMatrix(filename, f);
        }

        //x shift dt/2
        //#pragma omp target teams distribute parallel for map(from:f_xshift)
        #pragma omp parallel for schedule(guided)
        for(int j = 0; j < nv; j++)
        {
            Matrix<double, nx, 1> f_fixed_v_samples = f.col(j);
            CubicSplineInterp1D cubic_spline_interp_x(x_samples, f_fixed_v_samples, CubicSplineInterp1D::periodic, CubicSplineInterp1D::equal_interval);
            for(int i = 0; i < nx; i++)
            {
                x_shift = x_samples(i) - 0.5 * dt * v_samples(j);
                ShiftAPeriod(x_shift, L);
                f_xshift(i, j) = cubic_spline_interp_x.CalcVal(x_shift);
            }
        }

        //solve E
        //calculate \rho
        rho_e = 0.5 * ( f_xshift.block(0, 0, nx, nv - 1).rowwise().sum() +
                        f_xshift.block(0, 1, nx, nv - 1).rowwise().sum() ) * dv;

        for (int i = 0; i < nx; i++)
        {
            rho_i[i] =  GetIonInitDistrib(x_samples[i], 0.0, n * dt);
        }
        rho = rho_i - rho_e;
        //PoissonSolverDirichletBC poisson_solver(rho, dx, 0.0, 0.0);
        PoissonSolverPeriodicBC poisson_solver(rho, dx);

        if(n % data_steps == 0)
        {
            int nn = n / data_steps;
            string filename = data_path + "phi" + to_string(nn);
            OutputMatrix(filename, poisson_solver.phi);
        }
        //v shift dt
        //#pragma omp target teams distribute parallel for map(from:f_vshift)
        #pragma omp parallel for schedule(guided)
        for(int i = 0; i < nx; i++)
        {
            Matrix<double, nv, 1> f_fixed_x_samples = f_xshift.row(i);
            CubicSplineInterp1D cubic_spline_interp_v(v_samples, f_fixed_x_samples, CubicSplineInterp1D::fp_zero, CubicSplineInterp1D::equal_interval);
            for(int j = 1; j < nv; j++)
            {
                v_shift = v_samples(j) - (e / m) * (poisson_solver.GetEVal(i) + E_External(i * dx, n * dt) ) * dt;
                //when v is out of range, f =0
                f_vshift(i, j) =  0.0;
                if(v_shift < vmax || v_shift > -vmax)
                {
                    f_vshift(i, j) =  cubic_spline_interp_v.CalcVal(v_shift);
                }
            }
        }

        //2nd x shift dt/2
        //#pragma omp target teams distribute parallel for map(from:f)
        #pragma omp parallel for schedule(guided)
        for(int j = 0; j < nv; j++)
        {
            Matrix<double, nx, 1> f_fixed_v_samples2 = f_vshift.col(j);
            CubicSplineInterp1D cubic_spline_interp_x2(x_samples, f_fixed_v_samples2, CubicSplineInterp1D::periodic, CubicSplineInterp1D::equal_interval);
            for(int i = 0; i < nx; i++)
            {
                x_shift = x_samples(i) - 0.5 * dt * v_samples(j);
                ShiftAPeriod(x_shift, L);
                f(i, j) = cubic_spline_interp_x2.CalcVal(x_shift);
            }
        }
        //diagnose energy
        double E_sum = 0.0;
        double Ek_sum = 0.0;
        for(int i = 0; i < nx - 1; i++)
        {
            double E_temp = poisson_solver.GetEVal(i) + E_External(i * dx, n * dt);
            E_sum +=  E_temp * E_temp;
            for(int j = 0; j < nv - 1; j++)
            {
                Ek_sum += m * pow(-vmax + j * dv, 2) * f(i, j);
            }
        }
        //average energy per particle
        double Ep_temp = 0.5 * E_sum * dx / L;
        double Ek_temp = 0.5 * Ek_sum * dx * dv / L;
        Ep.push_back(Ep_temp);
        Ek.push_back(Ek_temp);
        Et.push_back(Ek_temp + Ep_temp);
    }
    OutputVector(data_path + "kin_energy", Ek);
    OutputVector(data_path + "pot_energy", Ep);
    OutputVector(data_path + "tot_energy", Et);
    cout << endl << "Simulation Finish!" << endl;
}

void Simulation::ShiftAPeriod(double& x, double length)
{
    if(x < 0) x += length;
    else if(x >= length)  x -= length;
}
