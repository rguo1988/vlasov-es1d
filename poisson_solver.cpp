/***************************
 * Solver of Poisson equation
 * Author: rguo
 * Data:   2021-8-26
 * Note:   solve laplace matrix equation by LU method
 ****************************/
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>
#include<eigen3/Eigen/QR>
#include"poisson_solver.h"

using namespace Eigen;

PoissonSolverPeriodicBC::PoissonSolverPeriodicBC(VectorXd _rho, double _dx): nx(_rho.size()), dx(_dx)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
    Calculate();
}

void PoissonSolverPeriodicBC::Calculate()
{
    //construct laplace opretor
    //-2  1  0  1
    //1 -2  1  0
    //0  1 -2  1
    //1  0  1 -2

    MatrixXd laplace(nx - 1, nx - 1);
    VectorXd diag_element = VectorXd::Constant(nx - 1, -2);
    VectorXd sub_diag_element = VectorXd::Constant(nx - 2, 1);
    laplace.diagonal() = diag_element;
    laplace.diagonal(1) = sub_diag_element;
    laplace.diagonal(-1) = sub_diag_element;
    laplace(0, nx - 2) = 1;
    laplace(nx - 2, 0) = 1;
    laplace /= pow(dx, 2);

    VectorXd phi_seg = laplace.fullPivLu().solve(-1.0 * rho.segment(0, nx - 1));
    phi.segment(0, nx - 1) = phi_seg;
    //phi(0) = 0.0;
    phi(nx - 1) = phi(0);
    //zero potential is defined as \int_L \phi \dd x = 0
    //VectorXd phi_sum = VectorXd::Constant(nx, phi.sum() / (nx - 1));
    //phi -= phi_sum;

    //calculate E
    for(int i = 0; i < nx - 1; i++)
    {
        int i_minus = (i == 0) ? (nx - 2) : (i - 1);
        int i_plus = (i == nx - 2) ? (0) : (i + 1);
        E(i) = -(phi(i_plus) - phi(i_minus)) / (2 * dx);
    }
    E(nx - 1) = E(0);
}

PoissonSolverDirichletBC::PoissonSolverDirichletBC(VectorXd _rho, double _dx, double _c1, double _c2): nx(_rho.size()), dx(_dx)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
    phi[0] = _c1;
    phi[nx - 1] = _c2;
    Calculate();
}

void PoissonSolverDirichletBC::Calculate()
{
    //construct laplace opretor
    //-2  1  0  0
    //1 -2  1  0
    //0  1 -2  1
    //0  0  1 -2

    double dx2 = dx * dx;
    rho[1] += phi[0] / dx2;
    rho[nx - 2] += phi[nx - 1] / dx2;

    MatrixXd laplace(nx - 2, nx - 2);
    VectorXd diag_element = VectorXd::Constant(nx - 2, -2);
    VectorXd sub_diag_element = VectorXd::Constant(nx - 3, 1);
    laplace.diagonal() = diag_element;
    laplace.diagonal(1) = sub_diag_element;
    laplace.diagonal(-1) = sub_diag_element;
    laplace /= dx2;

    VectorXd phi_seg = laplace.fullPivLu().solve(-1.0 * rho.segment(1, nx - 2));
    phi.segment(1, nx - 2) = phi_seg;

    //calculate E
    for(int i = 0; i < nx - 1; i++)
    {
        int i_minus = (i == 0) ? (nx - 2) : (i - 1);
        int i_plus = (i == nx - 2) ? (0) : (i + 1);
        E(i) = -(phi(i_plus) - phi(i_minus)) / (2 * dx);
    }
    E(nx - 1) = E(0);
}

PoissonSolverNaturalBC::PoissonSolverNaturalBC(VectorXd _rho, double _dx, double _d1, double _d2): nx(_rho.size()), dx(_dx), d1(_d1), d2(_d2)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
    //2nd order
    rho[0] -= 2.0 * d1 / dx;
    rho[nx - 1] += 2.0 * d2 / dx;
    //1st order
    //rho[1] -= d1 / dx;
    //rho[nx - 2] += d2 / dx;
    Calculate();
}

//2nd order
void PoissonSolverNaturalBC::Calculate()
{
    //construct laplace opretor
    //2nd order
    //-2  2  0  0
    // 1 -2  1  0
    // 0  1 -2  1
    // 0  0  2 -2

    MatrixXd laplace(nx, nx);
    VectorXd diag_element = VectorXd::Constant(nx, -2.0);
    VectorXd sub_diag_element = VectorXd::Constant(nx - 1, 1.0);
    laplace.diagonal() = diag_element;
    laplace.diagonal(1) = sub_diag_element;
    laplace.diagonal(-1) = sub_diag_element;
    laplace(0, 1) = 2.0;
    laplace(nx - 1, nx - 2) = 2.0;
    laplace /= pow(dx, 2);

    VectorXd phi_seg = laplace.fullPivLu().solve(-1.0 * rho);
    phi = phi_seg;

    //1st order
    //-1  1  0  0
    // 1 -2  1  0
    // 0  1 -2  1
    // 0  0  1 -1

    //MatrixXd laplace(nx - 2, nx - 2);
    //VectorXd diag_element = VectorXd::Constant(nx - 2, -2.0);
    //VectorXd sub_diag_element = VectorXd::Constant(nx - 3, 1.0);
    //laplace.diagonal() = diag_element;
    //laplace.diagonal(1) = sub_diag_element;
    //laplace.diagonal(-1) = sub_diag_element;
    //laplace(0, 0) = -1.0;
    //laplace(nx - 3, nx - 3) = -1.0;
    //laplace /= pow(dx, 2);

    //VectorXd phi_seg = laplace.fullPivLu().solve(-1.0 * rho.segment(1, nx - 2));
    //phi.segment(1, nx - 2) = phi_seg;
    //phi[0] = phi[1] - d1 * dx;
    //phi[nx - 1] = phi[nx - 2] + d2 * dx;

    //calculate E
    for(int i = 1; i <= nx - 2; i++)
    {
        E[i] = -(phi[i + 1] - phi[i - 1]) / (2.0 * dx);
    }
    E[0] = -d1;
    E[nx - 1] = -d2;
}

PoissonSolverRobinBC::PoissonSolverRobinBC(VectorXd _rho, double _dx, double _c1, double _d2): nx(_rho.size()), dx(_dx), c1(_c1), d2(_d2)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
    phi[0] = c1;
    rho[1] += phi[0] / dx / dx;
    rho[nx - 1] += d2 / dx;
    Calculate();
}

void PoissonSolverRobinBC::Calculate()
{
    //construct laplace opretor
    //-2  1  0  0
    //1 -2  1  0
    //0  1 -2  1
    //0  0  1 -1

    MatrixXd laplace(nx - 1, nx - 1);
    VectorXd diag_element = VectorXd::Constant(nx - 1, -2);
    VectorXd sub_diag_element = VectorXd::Constant(nx - 2, 1);
    laplace.diagonal() = diag_element;
    laplace.diagonal(1) = sub_diag_element;
    laplace.diagonal(-1) = sub_diag_element;
    laplace(nx - 2, nx - 2) = -1;
    laplace /= pow(dx, 2);

    VectorXd phi_seg = laplace.fullPivLu().solve(-1.0 * rho.segment(1, nx - 1));
    phi.segment(1, nx - 1) = phi_seg;

    //calculate E
    for(int i = 1; i <= nx - 2; i++)
    {
        E[i] = -(phi[i + 1] - phi[i - 1]) / (2 * dx);
    }
    E[0] = -(phi[1] - phi[0] / (2 * dx));
    E[nx - 1] = -d2;
}

PoissonSolverTwiceIntegral::PoissonSolverTwiceIntegral(VectorXd _rho, double _dx, double _c1, double _d2): nx(_rho.size()), dx(_dx), c1(_c1), d2(_d2)
{
    rho = _rho;
    phi.resize(nx);
    E.resize(nx);
    Calculate();
}
void PoissonSolverTwiceIntegral::Calculate()
{
    for(int i = 0; i <= nx - 1; i++)
    {
        double E_sum = -1.0 * d2;
        for(int j = 0; j <= i - 1; j++)
        {
            E_sum += 0.5 * (rho[j] + rho[j + 1]);
            //E_sum += rho[j];
        }
        E[i] = E_sum * dx;
    }
    for(int k = 0; k <= nx - 1; k++)
    {
        double phi_sum = c1;
        for(int l = 0; l <= k - 1; l++)
        {
            phi_sum += -0.5 * (E[l] + E[l + 1]);
            //phi_sum += -1.0 * E[l];
        }
        phi[k] = phi_sum * dx;
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
