/***************************
 * Solver of Poisson equation
 * Author: rguo
 * Data:   2022-8-31
 * Note:   Shooting method solving Poisson equation with Periodic, Dirichlet and Debye BC
 *         Periodic BC requires neutral total charges
 ****************************/
#include<eigen3/Eigen/Core>
#include"poisson_solver.h"
#include<iostream>

using namespace Eigen;

PoissonSolver::PoissonSolver(VectorXd _rho, double _dx): nx(_rho.size()), dx(_dx)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
}

PoissonSolverPeriodicBC::PoissonSolverPeriodicBC(VectorXd _rho, double _dx): PoissonSolver(_rho, _dx)
{
    Calculate();
}

//solve poisson equation with periodic condition
void PoissonSolverPeriodicBC::Calculate()
{
    //construct laplace opretor
    //here laplace is a reduced laplace oprator
    //e.g. a block(1,1)-(3,3)
    //-2  1  0  1
    //1 -2  1  0
    //0  1 -2  1
    //1  0  1 -2

    double dx2 = dx * dx;
    double r[nx - 1];
    double b[nx - 1];
    r[0] = b[0] = 0.0;
    r[1] = b[1] = 0.0;

    for(int i = 1; i <= nx - 3; i++)
    {
        r[i + 1] = -1.0 / (r[i] - 2);
        b[i + 1] = (-rho[i] * dx2 - b[i]) / (r[i] - 2);
    }

    phi(nx - 2) = (-rho(nx - 2) * dx2 - b[nx - 2]) / (r[nx - 2] - 2);

    for(int j = nx - 2; j >= 2; j--)
    {
        phi[j - 1] = r[j] * phi[j] + b[j];
    }

    phi[0] = 0.0;
    phi[nx - 1] = phi(0);
    //zero potential can be defined as \int_L \phi \dd x = 0
    //VectorXd phi_sum = VectorXd::Constant(nx, phi.sum() / (nx - 1));
    //phi -= phi_sum;

    //calculate E
    for(int i = 0; i <= nx - 2; i++)
    {
        int i_minus = (i == 0) ? (nx - 2) : (i - 1);
        int i_plus = (i == nx - 2) ? (0) : (i + 1);
        E[i] = -(phi[i_plus] - phi[i_minus]) / (2 * dx);
    }
    E[nx - 1] = E[0];
}

PoissonSolverDirichletBC::PoissonSolverDirichletBC(VectorXd _rho, double _dx, double c1, double c2): PoissonSolver(_rho, _dx)
{
    phi[0] = c1;
    phi[nx - 1] = c2;
    Calculate();
}

//solve poisson equation with Dirichlet condition
void PoissonSolverDirichletBC::Calculate()
{
    //ode: y''=f(x,y,y'); y(a) = c1; y(b) = c2;
    //solve: y''=f(x,y,y'); y(a) = c1; y'(a) = t1; obtain: y(b) = y(b, t1);
    //do k times:
    //tk = tk-1 - (y(b,t_k-1)-c2) (t_k-1-t_k-2)/ (y(b,t_k-1)-y(b,t_k-2))
    double t2 = 0.0;//set t2
    double t1 = 0.0006; //set t1
    double t0 = 0.6 * t1; //set t0
    double dx2 = dx * dx;
    double phi00, phi01, phi02 = 0.0;
    double phi10, phi11, phi12 = 0.0;
    for(int j = 1; j < 100; j++)
    {
        phi00 = phi10 = phi[0];
        phi01 = phi[0] + t0 * dx;
        phi11 = phi[0] + t1 * dx;
        for(int i = 1; i <= nx - 2; i++)
        {
            phi02 = -rho[i] * dx2 + 2.0 * phi01 - phi00;
            phi00 = phi01;
            phi01 = phi02;
            phi12 = -rho[i] * dx2 + 2.0 * phi11 - phi10;
            phi10 = phi11;
            phi11 = phi12;
        }
        if(abs(phi12 - phi[nx - 1]) < 1e-6)
            break;
        t2 = t1 - (phi12 - phi[nx - 1]) * (t1 - t0) / (phi12 - phi02);
        t0 = t1;
        t1 = t2;
    }

    phi[1] = phi[0] + t1 * dx;
    for(int i = 1; i <= nx - 2; i++)
        phi[i + 1] = -rho[i] * dx2 + 2.0 * phi[i] - phi[i - 1];

    //calculate E
    for(int i = 0; i <= nx - 1; i++)
    {
        double phi_plus = (i == nx - 1) ? phi[nx - 1] : phi[i + 1];
        double phi_minus = (i == 0) ? phi[0] : phi[i - 1];
        E[i] = -(phi_plus - phi_minus) / (2 * dx);
    }
}

PoissonSolverDebyeBC::PoissonSolverDebyeBC(VectorXd _rho, double _dx, double _l_D): PoissonSolver(_rho, _dx), l_D(_l_D)
{
    Calculate();
}

//solve poisson equation with Debye condition: phi' = phi/lambda_D and phi' = -phi/lambda_D
void PoissonSolverDebyeBC::Calculate()
{
    //ode: y''=f(x,y,y'); y'(a) = y(a)/l_D; y'(b) = -y(b)/l_D;
    //solve: y''=f(x,y,y'); y(a) = t1; y'(a) = t1/l_D; obtain: e(t1) = y'(b,t1) + y(b,t1)/l_D;
    //do k times:
    //tk = tk-1 - e(t_k-1)*(t_k-1-t_k-2)/ (e(t_k-1)-e(t_k-2))

    double t2 = 0.0;//set t2
    double t1 = 0.1; //set t1
    double t0 = 0; //set t0
    double dx2 = dx * dx;
    double phi00, phi01, phi02 = 0.0;
    double phi10, phi11, phi12 = 0.0;
    for(int j = 1; j < 100; j++)
    {
        phi00 = t0;
        phi10 = t1;
        phi01 = t0 + 2.0 * t0 / l_D * dx;
        phi11 = t1 + 2.0 * t1 / l_D * dx;
        for(int i = 1; i <= nx - 2; i++)
        {
            phi02 = -rho[i] * dx2 + 2.0 * phi01 - phi00;
            phi00 = phi01;
            phi01 = phi02;
            phi12 = -rho[i] * dx2 + 2.0 * phi11 - phi10;
            phi10 = phi11;
            phi11 = phi12;
        }
        double e0 = (phi01 - phi00) / 2.0 / dx + phi01 / l_D;
        double e1 = (phi11 - phi10) / 2.0 / dx + phi11 / l_D;
        if(abs(e1) < 1e-6)
            break;
        t2 = t1 - e1 * (t1 - t0) / (e1 - e0);
        t0 = t1;
        t1 = t2;
    }

    phi[0] = t1;
    phi[1] = t1 + 2.0 * t1 / l_D * dx;
    for(int i = 1; i <= nx - 2; i++)
        phi[i + 1] = -rho[i] * dx2 + 2.0 * phi[i] - phi[i - 1];

    //calculate E
    for(int i = 0; i <= nx - 1; i++)
    {
        double phi_plus = (i == nx - 1) ? phi[nx - 1] : phi[i + 1];
        double phi_minus = (i == 0) ? phi[0] : phi[i - 1];
        E[i] = -(phi_plus - phi_minus) / (2 * dx);
    }
}
double PoissonSolver::GetPhiVal(int x_index)
{
    return phi[x_index];
}
double PoissonSolver::GetEVal(int x_index)
{
    return E[x_index];
}
