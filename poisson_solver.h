#ifndef _poisson_solver_h
#define _poisson_solver_h
#include<eigen3/Eigen/Core>
using namespace Eigen;
class PoissonSolver
{
    const int nx;
    const double dx;
    VectorXd rho;
    //VectorXd E;

    void Calculate();
  public:
    PoissonSolver(VectorXd _rho, double _dx);
    double GetPhiVal(int x_index);
    double GetEVal(int x_index);
    VectorXd phi;
    VectorXd E;
};
#endif
