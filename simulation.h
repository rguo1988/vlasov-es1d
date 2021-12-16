#ifndef _simulation_h
#define _simulation_h
#include"input.h"
#include<eigen3/Eigen/Core>
#include<vector>
class Simulation: public Input
{
    Eigen::Matrix<double, nx, 1> x_samples;
    Eigen::Matrix<double, nv, 1> v_samples;
    Eigen::Matrix<double, nx, 1> rho_i;
    Eigen::Matrix<double, nx, 1> rho_e;
    Eigen::Matrix<double, nx, 1> rho;
    vector<double> Ek;
    vector<double> Ep;
    vector<double> Et;

  public:
    //Eigen::Matrix<double, nx, nv> f;
    Eigen::MatrixXd fe;
    Eigen::MatrixXd fi;
    Simulation();
    void Run();
    void ShiftAPeriod(double& x,double L);
};
#endif
