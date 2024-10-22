#ifndef _poisson_solver_h
#define _poisson_solver_h
#include<eigen3/Eigen/Core>
using namespace Eigen;

class PoissonSolver
{
  public:
    const int nx;
    const double dx;
    VectorXd rho;
    VectorXd phi;
    VectorXd E;

    double GetPhiVal(int x_index);
    double GetEVal(int x_index);
    PoissonSolver(VectorXd _rho, double _dx);
};

class PoissonSolverPeriodicBC: public PoissonSolver
{
    void Calculate();
  public:
    PoissonSolverPeriodicBC(VectorXd _rho, double _dx);
};

class PoissonSolverDirichletBC: public PoissonSolver
{
    void Calculate();
  public:
    PoissonSolverDirichletBC(VectorXd _rho, double _dx, double c1, double c2);
};

class PoissonSolverDebyeBC: public PoissonSolver
{
    const double l_D;
    void Calculate();
  public:
    PoissonSolverDebyeBC(VectorXd _rho, double _dx, double _l_D);
};
#endif
