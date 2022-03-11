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
    fe.resize(nx, nv);
    fe.setZero();
    fi.resize(nx, nv);
    fi.setZero();

    //CalculatePotentialSC();

    //electrons
    for (int i = 0; i < nx; i++)
    {
        for(int j = 0; j < nv; j++)
        {
            fe(i, j) = GetElecInitDistrib(x_samples[i], v_samples[j]);
        }
    }
    //ions
#if _ions_motion
    for (int i = 0; i < nx; i++)
    {
        for(int j = 0; j < nv; j++)
        {
            fi(i, j) = GetIonInitDistrib(x_samples[i], v_samples[j]);
        }
    }
#else
    for (int i = 0; i < nx; i++)
    {
        rho_i[i] = 1.0;
    }
#endif
    Ek.clear();
    Ep.clear();
    Et.clear();
}
void Simulation::Run()
{
//#if _ions_motion //mobile ions
    //set variables for electrons
    double xe_shift = 0.0;
    double ve_shift = 0.0;
    MatrixXd fe_xshift(nx, nv);
    fe_xshift.setZero();
    MatrixXd fe_vshift(nx, nv);
    fe_vshift.setZero();

#if _ions_motion //mobile ions
    //set variables for ions
    double xi_shift = 0.0;
    double vi_shift = 0.0;
    MatrixXd fi_xshift(nx, nv);
    fi_xshift.setZero();
    MatrixXd fi_vshift(nx, nv);
    fi_vshift.setZero();
#endif

    //show information
    cout << "************************************" << endl;
    cout << " Vlasov Simulation: " << title << endl;
    cout.setf(ios::left);
    cout << "************************************" << endl;
    cout << " Plasma Parameters: " << endl;
    cout << "       me = " << setw(8) << setprecision(6) << me
         << "       Te = " << setw(8) << setprecision(6) << Te
         << "     vt_e = " << setw(8) << setprecision(6) << vt_e << endl;
    cout << "     w_pe = " << setw(8) << setprecision(6) << w_pe
         << "     l_De = " << setw(8) << setprecision(6) << l_e << endl;
#if _ions_motion
    cout << "       mi = " << setw(8) << setprecision(6) << mi
         << "       Ti = " << setw(8) << setprecision(6) << Ti
         << "     vt_i = " << setw(8) << setprecision(6) << vt_i << endl;
    cout << "     w_pi = " << setw(8) << setprecision(6) << w_pi
         << "     l_Di = " << setw(8) << setprecision(6) << l_i << endl;
#else
    cout << "       mi = " << setw(8) << setprecision(6) << "inf"
         << "       Ti = " << setw(8) << setprecision(6) << "-"
         << "     vt_i = " << setw(8) << setprecision(6) << "-" << endl;
    cout << "     w_pi = " << setw(8) << setprecision(6) << "-"
         << "     l_Di = " << setw(8) << setprecision(6) << "-" << endl;
    cout << "Noting: Ions is Immobile!" << endl;
#endif
    cout << "************************************" << endl;
    cout << " Simulation Parameters: " << endl;
    cout << "        L = " << setw(8) << setprecision(6) << L
         << "       nx = " << setw(8) << setprecision(6) << nx
         << "       dx = " << setw(8) << setprecision(6) << dx << endl;
    cout << "        k = " << setw(8) << setprecision(6) << k << endl;
    cout << "     vmax = " << setw(8) << setprecision(6) << vmax
         << "       nv = " << setw(8) << setprecision(6) << nv
         << "       dv = " << setw(8) << setprecision(6) << dv << endl;
    cout << "       wT = " << setw(8) << setprecision(6) << max_steps*dt
         << "    steps = " << setw(8) << setprecision(6) << max_steps
         << "       dt = " << setw(8) << setprecision(6) << dt << endl;
    cout << "  dataNum = " << setw(8) << setprecision(6) << data_num << endl;

    PrintSpecialParameters();

//#if _ions_motion
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
            string filename_e = data_path + "fe" + to_string(nn);
            OutputMatrix(filename_e, fe);
            string filename_i = data_path + "fi" + to_string(nn);
            OutputMatrix(filename_i, fi);
        }

        //x shift dt/2
        #pragma omp parallel for schedule(guided)
        for(int j = 0; j < nv; j++)
        {
            //electrons
            Matrix<double, nx, 1> fe_fixed_v_samples = fe.col(j);
            CubicSplineInterp1D cubic_spline_interp_xe(x_samples, fe_fixed_v_samples, CubicSplineInterp1D::periodic, CubicSplineInterp1D::equal_interval);
            for(int i = 0; i < nx; i++)
            {
                xe_shift = x_samples(i) - 0.5 * dt * v_samples(j);
                //periodic bc
                ShiftAPeriod(xe_shift, L);
                //open bc
                //if(xe_shift > L || xe_shift < 0.0)
                //fe_xshift(i, j) = GetElecFreeDistrib(i * dx, -vmax + j * dv);
                //else
                fe_xshift(i, j) = cubic_spline_interp_xe.CalcVal(xe_shift);
            }
#if _ions_motion
            //ions
            Matrix<double, nx, 1> fi_fixed_v_samples = fi.col(j);
            CubicSplineInterp1D cubic_spline_interp_xi(x_samples, fi_fixed_v_samples, CubicSplineInterp1D::periodic, CubicSplineInterp1D::equal_interval);
            for(int i = 0; i < nx; i++)
            {
                xi_shift = x_samples(i) - 0.5 * dt * v_samples(j);
                //periodic bc
                ShiftAPeriod(xi_shift, L);
                //open bc
                //if(xi_shift > L || xi_shift < 0.0)
                //fi_xshift(i, j) = GetIonFreeDistrib(i * dx, -vmax + j * dv);
                //else
                fi_xshift(i, j) = cubic_spline_interp_xi.CalcVal(xi_shift);
            }
#endif
        }

        //solve E
#if _ions_motion
        MatrixXd df = fi_xshift - fe_xshift;
        rho = 0.5 * ( df.block(0, 0, nx, nv - 1).rowwise().sum() + df.block(0, 1, nx, nv - 1).rowwise().sum() ) * dv;
#else
        rho_e = 0.5 * ( fe_xshift.block(0, 0, nx, nv - 1).rowwise().sum() + fe_xshift.block(0, 1, nx, nv - 1).rowwise().sum() ) * dv;
        rho = rho_i - rho_e;
#endif
        PoissonSolverDirichletBC poisson_solver(rho, dx, 0.0, 0.0);
        //PoissonSolverPeriodicBC poisson_solver(rho, dx);

        if(n % data_steps == 0)
        {
            int nn = n / data_steps;
            string filename = data_path + "phi" + to_string(nn);
            OutputMatrix(filename, poisson_solver.phi);
        }
        //v shift dt
        #pragma omp parallel for schedule(guided)
        for(int i = 0; i < nx; i++)
        {
            //electrons
            Matrix<double, nv, 1> fe_fixed_x_samples = fe_xshift.row(i);
            CubicSplineInterp1D cubic_spline_interp_ve(v_samples, fe_fixed_x_samples, CubicSplineInterp1D::fp_zero, CubicSplineInterp1D::equal_interval);
            for(int j = 1; j < nv; j++)
            {
                ve_shift = v_samples(j) - (e / me) * poisson_solver.E(i) * dt;
                //when v is out of range, f =0
                fe_vshift(i, j) =  0.0;
                if(ve_shift < vmax || ve_shift > -vmax)
                {
                    fe_vshift(i, j) =  cubic_spline_interp_ve.CalcVal(ve_shift);
                }
            }
#if _ions_motion
            //ions
            Matrix<double, nv, 1> fi_fixed_x_samples = fi_xshift.row(i);
            CubicSplineInterp1D cubic_spline_interp_vi(v_samples, fi_fixed_x_samples, CubicSplineInterp1D::fp_zero, CubicSplineInterp1D::equal_interval);
            for(int j = 1; j < nv; j++)
            {
                vi_shift = v_samples(j) + (e / mi) * poisson_solver.E(i) * dt;
                //when v is out of range, f =0
                fi_vshift(i, j) =  0.0;
                if(vi_shift < vmax || vi_shift > -vmax)
                {
                    fi_vshift(i, j) =  cubic_spline_interp_vi.CalcVal(vi_shift);
                }
            }
#endif
        }

        //2nd x shift dt/2
        #pragma omp parallel for schedule(guided)
        for(int j = 0; j < nv; j++)
        {
            //electrons
            Matrix<double, nx, 1> fe_fixed_v_samples2 = fe_vshift.col(j);
            CubicSplineInterp1D cubic_spline_interp_xe2(x_samples, fe_fixed_v_samples2, CubicSplineInterp1D::periodic, CubicSplineInterp1D::equal_interval);
            for(int i = 0; i < nx; i++)
            {
                xe_shift = x_samples(i) - 0.5 * dt * v_samples(j);
                //periodic bc
                ShiftAPeriod(xe_shift, L);
                //open bc
                //if(xe_shift > L || xe_shift < 0.0)
                //fe(i, j) = GetElecFreeDistrib(i * dx, -vmax + j * dv);
                //else
                fe(i, j) = cubic_spline_interp_xe2.CalcVal(xe_shift);
            }
#if _ions_motion
            //ions
            Matrix<double, nx, 1> fi_fixed_v_samples2 = fi_vshift.col(j);
            CubicSplineInterp1D cubic_spline_interp_xi2(x_samples, fi_fixed_v_samples2, CubicSplineInterp1D::periodic, CubicSplineInterp1D::equal_interval);
            for(int i = 0; i < nx; i++)
            {
                xi_shift = x_samples(i) - 0.5 * dt * v_samples(j);
                //periodic bc
                ShiftAPeriod(xi_shift, L);
                //open bc
                //if(xi_shift > L || xi_shift < 0.0)
                //fi(i, j) = GetIonFreeDistrib(i * dx, -vmax + j * dv);
                //else
                fi(i, j) = cubic_spline_interp_xi2.CalcVal(xi_shift);
            }
#endif
        }
        //diagnose energy
        double E_sum = 0.0;
        double Ek_sum = 0.0;
        for(int i = 0; i < nx - 1; i++)
        {
            double E_temp = poisson_solver.E(i);
            E_sum +=  E_temp * E_temp;
            for(int j = 0; j < nv - 1; j++)
            {
#if _ions_motion
                Ek_sum += pow(-vmax + j * dv, 2) * (me * fe(i, j) + mi * fi(i, j));
#else
                Ek_sum += pow(-vmax + j * dv, 2) * me * fe(i, j);
#endif
            }
        }
        //average energy per particle
        double Ep_temp = 0.5 * E_sum * dx / L;
        double Ek_temp = 0.5 * Ek_sum * dx * dv / L;
        Ep.push_back(Ep_temp);
        Ek.push_back(Ek_temp);
        Et.push_back(Ek_temp + Ep_temp);
    }

//#else  //immobile ions
    //for(int n = 0; n < max_steps + 1; n++)
    //{
    ////print running process
    //int percent = 100 * n / (max_steps - 1);
    //if(percent % 5 == 0)
    //{
    //cout << "\r" << " Process: " << percent << "%" << flush;
    //}
    ////diagnose
    //if(n % data_steps == 0)
    //{
    //int nn = n / data_steps;
    //string filename_e = data_path + "fe" + to_string(nn);
    //OutputMatrix(filename_e, fe);
    //string filename_i = data_path + "fi" + to_string(nn);
    //OutputMatrix(filename_i, fi);
    //}

    ////x shift dt/2
    ////#pragma omp target teams distribute parallel for map(from:f_xshift)
    //#pragma omp parallel for schedule(guided)
    //for(int j = 0; j < nv; j++)
    //{
    ////electrons
    //Matrix<double, nx, 1> fe_fixed_v_samples = fe.col(j);
    //CubicSplineInterp1D cubic_spline_interp_xe(x_samples, fe_fixed_v_samples, CubicSplineInterp1D::periodic, CubicSplineInterp1D::equal_interval);
    //for(int i = 0; i < nx; i++)
    //{
    //xe_shift = x_samples(i) - 0.5 * dt * v_samples(j);
    ////if(xe_shift > L || xe_shift < 0.0)
    ////periodic bc
    //ShiftAPeriod(xe_shift, L);
    ////open bc
    ////fe_xshift(i, j) = GetElecFreeDistrib(i * dx, -vmax + j * dv);
    ////else
    //fe_xshift(i, j) = cubic_spline_interp_xe.CalcVal(xe_shift);
    //}
    //}

    ////solve E
    ////calculate \rho
    //rho_e = 0.5 * ( fe_xshift.block(0, 0, nx, nv - 1).rowwise().sum() + fe_xshift.block(0, 1, nx, nv - 1).rowwise().sum() ) * dv;
    //rho = rho_i - rho_e;
    ////PoissonSolverDirichletBC poisson_solver(rho, dx, 0.0, 0.0);
    //PoissonSolverPeriodicBC poisson_solver(rho, dx);

    //if(n % data_steps == 0)
    //{
    //int nn = n / data_steps;
    //string filename = data_path + "phi" + to_string(nn);
    //OutputMatrix(filename, poisson_solver.phi);
    //}
    ////v shift dt
    ////#pragma omp target teams distribute parallel for map(from:f_vshift)
    //#pragma omp parallel for schedule(guided)
    //for(int i = 0; i < nx; i++)
    //{
    ////electrons
    //Matrix<double, nv, 1> fe_fixed_x_samples = fe_xshift.row(i);
    //CubicSplineInterp1D cubic_spline_interp_ve(v_samples, fe_fixed_x_samples, CubicSplineInterp1D::fp_zero, CubicSplineInterp1D::equal_interval);
    //for(int j = 1; j < nv; j++)
    //{
    //ve_shift = v_samples(j) - (e / me) * poisson_solver.E(i) * dt;
    ////when v is out of range, f =0
    //fe_vshift(i, j) =  0.0;
    //if(ve_shift < vmax || ve_shift > -vmax)
    //{
    //fe_vshift(i, j) =  cubic_spline_interp_ve.CalcVal(ve_shift);
    //}
    //}
    //}

    ////2nd x shift dt/2
    ////#pragma omp target teams distribute parallel for map(from:f)
    //#pragma omp parallel for schedule(guided)
    //for(int j = 0; j < nv; j++)
    //{
    ////electrons
    //Matrix<double, nx, 1> fe_fixed_v_samples2 = fe_vshift.col(j);
    //CubicSplineInterp1D cubic_spline_interp_xe2(x_samples, fe_fixed_v_samples2, CubicSplineInterp1D::periodic, CubicSplineInterp1D::equal_interval);
    //for(int i = 0; i < nx; i++)
    //{
    //xe_shift = x_samples(i) - 0.5 * dt * v_samples(j);
    ////if(xe_shift > L || xe_shift < 0.0)
    ////periodic bc
    //ShiftAPeriod(xe_shift, L);
    ////open bc
    ////fe(i, j) = GetElecFreeDistrib(i * dx, -vmax + j * dv);
    ////else
    //fe(i, j) = cubic_spline_interp_xe2.CalcVal(xe_shift);
    //}
    //}
    ////diagnose energy
    //double E_sum = 0.0;
    //double Ek_sum = 0.0;
    //for(int i = 0; i < nx - 1; i++)
    //{
    //double E_temp = poisson_solver.E(i);
    //E_sum +=  E_temp * E_temp;
    //for(int j = 0; j < nv - 1; j++)
    //{
    //Ek_sum += pow(-vmax + j * dv, 2) * (me * fe(i, j));
    //}
    //}
    ////average energy per particle
    //double Ep_temp = 0.5 * E_sum * dx / L;
    //double Ek_temp = 0.5 * Ek_sum * dx * dv / L;
    //Ep.push_back(Ep_temp);
    //Ek.push_back(Ek_temp);
    //Et.push_back(Ek_temp + Ep_temp);
    //}
//#endif
    OutputVector(data_path + "kin_energy", Ek);
    OutputVector(data_path + "pot_energy", Ep);
    OutputVector(data_path + "tot_energy", Et);
    cout << endl << " Simulation Finish!" << endl;
}

void Simulation::ShiftAPeriod(double& x, double length)
{
    if(x < 0) x += length;
    else if(x >= length)  x -= length;
}
