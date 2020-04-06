/***************************
 * Solver of Poisson equation
 * Author: rguo
 * Data:   2020-3-10
 * Note:   A direct Thomas algorithm with periodic BC
 ****************************/
#include<eigen3/Eigen/Core>
#include"poisson_solver.h"

using namespace Eigen;

PoissonSolver::PoissonSolver(VectorXd _rho, double _dx): nx(_rho.size()), dx(_dx)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
    Calculate();
}

//solve poisson equation with periodic condition
void PoissonSolver::Calculate()
{
    //construct laplace opretor
    //here laplace is a reduced laplace oprator
    //e.g. a block(1,1)-(3,3)
    //-2  1  0  1
    //1 -2  1  0
    //0  1 -2  1
    //1  0  1 -2

    double dx2 = dx*dx;
    double r[nx-1];
    double b[nx-1];
    //VectorXd r(nx - 1);
    //VectorXd b(nx - 1);
    r[0] = b[0] = 0.0;
    r[1] = b[1] = 0.0;
    for(int i = 1; i <= nx - 3; i++)
    {
        r[i + 1] = -1.0 / (r[i] - 2);
        b[i + 1] = (-rho[i]*dx2 - b[i]) / (r[i] - 2);
    }
    phi(nx - 2) = (-rho(nx - 2)*dx2 - b[nx - 2]) / (r[nx - 2] - 2);
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
    for(int i = 0; i < nx - 1; i++)
    {
        int i_minus = (i == 0) ? (nx - 2) : (i - 1);
        int i_plus = (i == nx - 2) ? (0) : (i + 1);
        E[i] = -(phi[i_plus] - phi[i_minus]) / (2 * dx);
    }
    E[nx - 1] = E[0];
}
double PoissonSolver::GetPhiVal(int x_index)
{
    return phi[x_index];
}
double PoissonSolver::GetEVal(int x_index)
{
    return E[x_index];
}
