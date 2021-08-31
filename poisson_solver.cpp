/***************************
 * Solver of Poisson equation
 * Author: rguo
 * Data:   2021-8-26
 * Note:   A direct Thomas algorithm with Periodic BC
 *         Periodic BC requires neutral total charges
 ****************************/
#include<eigen3/Eigen/Core>
#include"poisson_solver.h"
#include<iostream>

using namespace Eigen;

PoissonSolverPeriodicBC::PoissonSolverPeriodicBC(VectorXd _rho, double _dx): nx(_rho.size()), dx(_dx)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
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

PoissonSolverDirichletBC::PoissonSolverDirichletBC(VectorXd _rho, double _dx, double c1, double c2): nx(_rho.size()), dx(_dx)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
    phi[0] = c1;
    phi[nx - 1] = c2;
    Calculate();
}

//solve poisson equation with Dirichlet condition
void PoissonSolverDirichletBC::Calculate()
{
    //construct laplace opretor
    //-2  1  0  0
    //1 -2  1  0
    //0  1 -2  1
    //0  0  1 -2

    double dx2 = dx * dx;
    double r[nx - 1];
    double b[nx - 1];
    r[0] = b[0] = 0.0;
    r[1] = b[1] = 0.0;
    rho[1] += phi[0] / dx2;
    rho[nx - 2] += phi[nx - 1] / dx2;
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

    //calculate E
    for(int i = 0; i <= nx - 1; i++)
    {
        double phi_plus = (i == nx - 1) ? phi[nx - 1] : phi[i + 1];
        double phi_minus = (i == 0) ? phi[0] : phi[i - 1];
        E[i] = -(phi_plus - phi_minus) / (2 * dx);
    }
}

PoissonSolverNaturalBC::PoissonSolverNaturalBC(VectorXd _rho, double _dx, double _d1, double _d2): nx(_rho.size()), dx(_dx), d1(_d1), d2(_d2)
{
    rho = _rho;
    rho[1] -= d1 / dx;
    rho[nx - 1] += d2 / dx;
    phi.resize(nx);
    E.resize(nx);
    Calculate();
}

//solve poisson equation with Natural condition
void PoissonSolverNaturalBC::Calculate()
{
    //construct laplace opretor
    //-1  1  0  0
    // 1 -2  1  0
    // 0  1 -2  1
    // 0  0  1 -1

    double dx2 = dx * dx;
    double r[nx];
    double b[nx];
    r[0] = b[0] = 0.0;
    r[1] = b[1] = 0.0;
    r[2] = 1.0;
    b[2] = rho[1] * dx2;

    for(int i = 2; i <= nx - 2; i++)
    {
        r[i + 1] = -1.0 / (r[i] - 2);
        b[i + 1] = (-rho[i] * dx2 - b[i]) / (r[i] - 2);
    }

    std::cout<<r[nx-2]<<std::endl;
    phi(nx - 1) = (-rho(nx - 1) * dx2 - b[nx - 1]) / (r[nx - 1] - 1.0);
    for(int j = nx - 1; j >= 2; j--)
    {
        phi[j - 1] = r[j] * phi[j] + b[j];
    }
    phi[0] = -d1 * dx + phi[1];

    //calculate E
    for(int i = 1; i <= nx - 2; i++)
    {
        E[i] = -(phi[i + 1] - phi[i - 1]) / (2 * dx);
    }
    E[0] = -d1;
    E[nx - 1] = -d2;
}
double PoissonSolver::GetPhiVal(int x_index)
{
    return phi[x_index];
}
double PoissonSolver::GetEVal(int x_index)
{
    return E[x_index];
}
