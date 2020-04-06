#ifndef _cubic_spline_1d_h
#define _cubic_spline_1d_h
#include<eigen3/Eigen/Core>
using namespace Eigen;
class CubicSplineInterp1D
{
    const int samples;

    VectorXd x_samples;
    VectorXd y_samples;

    VectorXd m;
    VectorXd h;
    VectorXd s_a;
    VectorXd s_b;
    VectorXd s_c;
    VectorXd s_d;

    void CalcCoefficientsPeriodic();
    void CalcCoefficientsPeriodicEI();
    void CalcCoefficientsNatrual();
    void CalcCoefficientsNatrualEI();
    void CalcCoefficientsFpZero();
    void CalcCoefficientsFpZeroEI();

  public:
    enum BoundaryCondition {periodic, natural, fp_zero};
    enum Interval {equal_interval, inequal_interval};
    CubicSplineInterp1D(VectorXd _x, VectorXd _y, BoundaryCondition bc, Interval interval);
    double CalcVal(double x);
};
#endif
