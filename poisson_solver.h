#ifndef _poisson_solver_h
#define _poisson_solver_h
#include<eigen3/Eigen/Core>
using namespace Eigen;
class PoissonSolver
{
  public:
    VectorXd rho;
    VectorXd phi;
    VectorXd E;
    double GetPhiVal(int x_index);
    double GetEVal(int x_index);
};
class PoissonSolverPeriodicBC: public PoissonSolver
{
    const int nx;
    const double dx;
    void Calculate();
  public:
    PoissonSolverPeriodicBC(VectorXd _rho, double _dx);
};
class PoissonSolverDirichletBC: public PoissonSolver
{
    const int nx;
    const double dx;
    void Calculate();
  public:
    PoissonSolverDirichletBC(VectorXd _rho, double _dx, double _c1, double _c2);
};
class PoissonSolverNaturalBC: public PoissonSolver
{
    const int nx;
    const double dx;
    const double d1;
    const double d2;
    void Calculate();
  public:
    PoissonSolverNaturalBC(VectorXd _rho, double _dx, double _d1, double _d2);
};
class PoissonSolverRobinBC: public PoissonSolver
{
    const int nx;
    const double dx;
    const double c1;
    const double d2;
    void Calculate();
  public:
    PoissonSolverRobinBC(VectorXd _rho, double _dx, double _c1, double _d2);
};
class PoissonSolverTwiceIntegral: public PoissonSolver
{
    const int nx;
    const double dx;
    const double c1;
    const double d2;
    void Calculate();
  public:
    PoissonSolverTwiceIntegral(VectorXd _rho, double _dx, double _c1, double _d2);
};
#endif
