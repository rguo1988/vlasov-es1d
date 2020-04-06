/*********************************
 * Cubic Spline Interpolation for 1D
 * Author: rguo
 * Data:   2020-3-10
 * Note:  1. Suitable for samples with non equal interval,
 *        2. Optimized for equal interval samples
 *        3. Periodic bc and natrual bc
 *        4. Tridiagonal matrix algorithm(TDMA) is used to solve matrix equation for M
**********************************/
#include"cubic_spline_1d.h"
#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Dense>

//cubic spline interpolation with periodic condition
CubicSplineInterp1D::CubicSplineInterp1D(VectorXd _x, VectorXd _y, BoundaryCondition boundary_condition, Interval interval):
    samples(_x.size())
{
    x_samples = _x;
    y_samples = _y;

    m.resize(samples);
    h.resize(samples);
    s_a.resize(samples - 1);
    s_b.resize(samples - 1);
    s_c.resize(samples - 1);
    s_d.resize(samples - 1);

    switch(boundary_condition)
    {
    case periodic:
        switch(interval)
        {
        case equal_interval:
            CalcCoefficientsPeriodicEI();
            break;
        case inequal_interval:
            CalcCoefficientsPeriodic();
            break;
        }
        break;
    case natural:
        switch(interval)
        {
        case equal_interval:
            CalcCoefficientsNatrualEI();
            break;
        case inequal_interval:
            CalcCoefficientsNatrual();
            break;
        }
        break;
    case fp_zero:
        switch(interval)
        {
        case equal_interval:
            CalcCoefficientsFpZeroEI();
            break;
        case inequal_interval:
            CalcCoefficientsFpZero();
            break;
        }
        break;
    default:
        CalcCoefficientsPeriodic();
    }
}
void CubicSplineInterp1D::CalcCoefficientsPeriodicEI()
{
    //solve tridiagonal matrix by thomas algorithm (TDMA)
    //calculating h
    double h = x_samples[1] - x_samples[0];

    //calculating d_i
    double d[samples];
    for(int i = 1; i < samples - 1; i++)
    {
        d[i] = 6 * ( y_samples[i + 1] - 2 * y_samples[i] + y_samples[i - 1] )  / h;
    }

    //periodic condition
    d[samples - 1] = 6 * ( y_samples[1] - 2 * y_samples[0] + y_samples[samples - 2] ) / h;
    d[0] = d[samples - 1];

    double gamma[samples];
    double beta[samples];
    double s[samples];
    gamma[0] = beta[0] = s[0] = 0.0;
    gamma[1] = beta[1] = 0.0;
    s[1] = 1.0;

    for(int i = 1; i < samples - 1; i++)
    {
        double denominator = gamma[i] + 4.0;
        gamma[i + 1] = -1.0 / denominator;
        beta[i + 1]  = (d[i] / h - beta[i]) / denominator;
        s[i + 1]     = -s[i] / denominator;
    }

    double v[samples];
    double t[samples];
    v[0] = t[0] = 0.0;
    v[samples - 1] = 1.0;
    t[samples - 1] = 0.0;
    for(int i = samples - 1; i >= 2; i--)
    {
        v[i - 1] = gamma[i] * v[i] + s[i];
        t[i - 1] = gamma[i] * t[i] + beta[i];
    }
    m[samples - 1] = ( d[samples - 1] / h - t[1] - t[samples - 2] ) /
                     ( v[1] + v[samples - 2] + 4 );
    for(int i = samples - 2; i >= 1; i--)
    {
        m[i] = v[i] * m[samples - 1] + t[i];
    }
    m[0] = m[samples - 1];

    //calculate  s_a,s_b,s_c,s_d
    for(int i = 0; i < samples - 1; i++)
    {
        s_a(i) = m(i + 1) / h / 6.0;
        s_b(i) = m(i) / h / 6.0;
        s_c(i) = y_samples(i + 1) / h - m(i + 1) * h / 6.0 ;
        s_d(i) = y_samples(i) / h - m(i) * h / 6.0 ;
    }
}
void CubicSplineInterp1D::CalcCoefficientsPeriodic()
{
    //solve tridiagonal matrix by varied thomas algorithm (TDMA)
    //calculating h_i
    for (int i = 0; i < samples - 1; i++)
    {
        h[i] = x_samples[i + 1] - x_samples[i];
    }
    h[samples - 1] = h[0]; //periodic condition

    //calculating d_i
    double d[samples];
    for(int i = 1; i < samples - 1; i++)
    {
        d[i] = 6 * ( (y_samples[i + 1] - y_samples[i]) / h[i] - (y_samples[i] - y_samples[i - 1]) / h[i - 1] );
    }

    //periodic condition
    d[samples - 1] = 6 * ( (y_samples[1] - y_samples[0]) / h[0] - (y_samples[samples - 1] - y_samples[samples - 2]) / h[samples - 2] );
    d[0] = d[samples - 1];

    double a[samples];
    double b[samples];
    double c[samples];
    double gamma[samples];
    double beta[samples];
    double s[samples];
    a[0] = b[0] = c[0] = 0.0;
    gamma[0] = beta[0] = s[0] = 0.0;
    gamma[1] = beta[1] = 0.0;
    s[1] = 1.0;
    for(int i = 1; i < samples; i++)
    {
        a[i] = h[i - 1];
        b[i] = 2 * (h[i - 1] + h[i]);
        c[i] = h[i];
    }
    for(int i = 1; i < samples - 1; i++)
    {
        double denominator = a[i] * gamma[i] + b[i];
        gamma[i + 1] = -c[i] / denominator;
        beta[i + 1]  = (d[i] - a[i] * beta[i]) / denominator;
        s[i + 1]     = -a[i] * s[i] / denominator;
    }
    double v[samples];
    double t[samples];
    v[0] = t[0] = 0.0;
    v[samples - 1] = 1.0;
    t[samples - 1] = 0.0;
    for(int i = samples - 1; i >= 2; i--)
    {
        v[i - 1] = gamma[i] * v[i] + s[i];
        t[i - 1] = gamma[i] * t[i] + beta[i];
    }
    m[samples - 1] = ( d[samples - 1] - c[samples - 1] * t[1] - a[samples - 1] * t[samples - 2] ) /
                     ( c[samples - 1] * v[1] + a[samples - 1] * v[samples - 2] + b[samples - 1] );
    for(int i = samples - 2; i >= 1; i--)
    {
        m[i] = v[i] * m[samples - 1] + t[i];
    }
    m[0] = m[samples - 1];

    //calculate  s_a,s_b,s_c,s_d
    for(int i = 0; i < samples - 1; i++)
    {
        s_a(i) = m(i + 1) / h(i) / 6.0;
        s_b(i) = m(i) / h(i) / 6.0;
        s_c(i) = y_samples(i + 1) / h(i) - m(i + 1) * h(i) / 6.0 ;
        s_d(i) = y_samples(i) / h(i) - m(i) * h(i) / 6.0 ;
    }
}

void CubicSplineInterp1D::CalcCoefficientsNatrual()
{
    //solve tridiagonal matrix by varied thomas algorithm (TDMA)
    //calculating h_i
    for (int i = 0; i < samples - 1; i++)
    {
        h[i] = x_samples[i + 1] - x_samples[i];
    }
    h[samples - 1] = 0.0; //useless in natrual boundary condition

    //calculating d_i
    double d[samples];
    for(int i = 1; i < samples - 1; i++)
    {
        d[i] = 6 * ( (y_samples[i + 1] - y_samples[i]) / h[i] - (y_samples[i] - y_samples[i - 1]) / h[i - 1] );
    }

    //natrual boundary condition
    d[0] = d[samples - 1] = 0.0;

    double a[samples];
    double b[samples];
    double c[samples];
    double gamma[samples];
    double beta[samples];

    a[0] = c[0] = 0.0;
    b[0] = 1;
    for(int i = 1; i <= samples - 2; i++)
    {
        a[i] = h[i - 1];
        b[i] = 2 * (h[i - 1] + h[i]);
        c[i] = h[i];
    }
    a[samples - 1] = 0.0;
    b[samples - 1] = 1.0;

    gamma[0] = beta[0] = 0.0;
    for(int i = 0; i < samples - 1; i++)
    {
        double denominator = a[i] * gamma[i] + b[i];
        gamma[i + 1] = -c[i] / denominator;
        beta[i + 1]  = (d[i] - a[i] * beta[i]) / denominator;
    }

    m[0] = m[samples - 1] = 0.0;
    for(int i = samples - 1; i >= 2; i--)
    {
        m[i - 1] = gamma[i] * m[i] + beta[i];
    }

    //calculate  s_a,s_b,s_c,s_d
    for(int i = 0; i < samples - 1; i++)
    {
        s_a(i) = m(i + 1) / h(i) / 6.0;
        s_b(i) = m(i) / h(i) / 6.0;
        s_c(i) = y_samples(i + 1) / h(i) - m(i + 1) * h(i) / 6.0 ;
        s_d(i) = y_samples(i) / h(i) - m(i) * h(i) / 6.0 ;
    }
}
void CubicSplineInterp1D::CalcCoefficientsNatrualEI()
{
    //solve tridiagonal matrix by varied thomas algorithm (TDMA)
    //calculating h_i
    double h = x_samples[1] - x_samples[0];
    //calculating d_i
    double d[samples];
    for(int i = 1; i < samples - 1; i++)
    {
        d[i] = 6 * ( y_samples[i + 1] - 2 * y_samples[i] + y_samples[i - 1]  ) / h;
    }
    //natrual boundary condition
    d[0] = d[samples - 1] = 0.0;

    double gamma[samples];
    double beta[samples];
    gamma[0] = beta[0] = 0.0;
    gamma[1] = beta[1] = 0.0;
    for(int i = 1; i < samples - 1; i++)
    {
        double denominator =  gamma[i] + 4.0;
        gamma[i + 1] = -1.0 / denominator;
        beta[i + 1]  = (d[i] / h - beta[i]) / denominator;
    }
    m[0] = m[samples - 1] = 0.0;
    for(int i = samples - 1; i >= 2; i--)
    {
        m[i - 1] = gamma[i] * m[i] + beta[i];
    }

    //calculate  s_a,s_b,s_c,s_d
    for(int i = 0; i < samples - 1; i++)
    {
        s_a(i) = m(i + 1) / h / 6.0;
        s_b(i) = m(i) / h / 6.0;
        s_c(i) = y_samples(i + 1) / h - m(i + 1) * h / 6.0 ;
        s_d(i) = y_samples(i) / h - m(i) * h / 6.0 ;
    }
}
void CubicSplineInterp1D::CalcCoefficientsFpZero()
{
    //solve tridiagonal matrix by varied thomas algorithm (TDMA)
    //calculating h_i
    for (int i = 0; i < samples - 1; i++)
    {
        h[i] = x_samples[i + 1] - x_samples[i];
    }
    h[samples - 1] = 0.0; //useless in FpZero boundary condition

    //calculating d_i
    double d[samples];
    for(int i = 1; i <= samples - 2; i++)
    {
        d[i] = 6 * ( (y_samples[i + 1] - y_samples[i]) / h[i] - (y_samples[i] - y_samples[i - 1]) / h[i - 1] );
    }

    //FpZero boundary condition
    d[0]           = 6 * (y_samples[1] - y_samples[0]) / h[0];
    d[samples - 1] = 6 * (y_samples[samples - 2] - y_samples[samples - 1]) / h[samples - 2];

    double a[samples];
    double b[samples];
    double c[samples];
    double gamma[samples];
    double beta[samples];

    a[0] = 0.0;
    b[0] = 2 * h[0];
    c[0] = h[0];
    for(int i = 1; i <= samples - 2; i++)
    {
        a[i] = h[i - 1];
        b[i] = 2 * (h[i - 1] + h[i]);
        c[i] = h[i];
    }
    a[samples - 1] = h[samples - 2];
    b[samples - 1] = 2 * h[samples - 2];
    c[samples - 1] = 0.0;

    gamma[0] = beta[0] = 0.0;
    for(int i = 0; i < samples - 1; i++)
    {
        double denominator = a[i] * gamma[i] + b[i];
        gamma[i + 1] = -c[i] / denominator;
        beta[i + 1]  = (d[i] - a[i] * beta[i]) / denominator;
    }

    m[samples - 1] = (d[samples - 1] - a[samples - 1] * beta[samples - 1]) /
                     (a[samples - 1] * gamma[samples - 1] + b[samples - 1]);
    for(int i = samples - 1; i >= 1; i--)
    {
        m[i - 1] = gamma[i] * m[i] + beta[i];
    }

    //calculate  s_a,s_b,s_c,s_d
    for(int i = 0; i < samples - 1; i++)
    {
        s_a(i) = m(i + 1) / h(i) / 6.0;
        s_b(i) = m(i) / h(i) / 6.0;
        s_c(i) = y_samples(i + 1) / h(i) - m(i + 1) * h(i) / 6.0 ;
        s_d(i) = y_samples(i) / h(i) - m(i) * h(i) / 6.0 ;
    }
}
void CubicSplineInterp1D::CalcCoefficientsFpZeroEI()
{
    //solve tridiagonal matrix by varied thomas algorithm (TDMA)
    //calculating h_i
    double h = x_samples[1] - x_samples[0];
    //calculating d_i
    double d[samples];
    for(int i = 1; i <= samples - 2; i++)
    {
        d[i] = 6 * ( y_samples[i + 1] - 2 * y_samples[i] + y_samples[i - 1] ) / h;
    }

    //FpZero boundary condition
    d[0]           = 6 * (y_samples[1] - y_samples[0]) / h;
    d[samples - 1] = 6 * (y_samples[samples - 2] - y_samples[samples - 1]) / h;

    double gamma[samples];
    double beta[samples];

    gamma[0] = beta[0] = 0.0;
    gamma[1] = -0.5;
    beta[1] = d[0] / (2 * h);
    for(int i = 1; i < samples - 1; i++)
    {
        double denominator = gamma[i] + 4.0;
        gamma[i + 1] = -1.0 / denominator;
        beta[i + 1]  = (d[i] / h - beta[i]) / denominator;
    }

    m[samples - 1] = (d[samples - 1] / h - beta[samples - 1]) /
                     (gamma[samples - 1] + 2.0);
    for(int i = samples - 1; i >= 1; i--)
    {
        m[i - 1] = gamma[i] * m[i] + beta[i];
    }

    //calculate  s_a,s_b,s_c,s_d
    for(int i = 0; i < samples - 1; i++)
    {
        s_a(i) = m(i + 1) / h / 6.0;
        s_b(i) = m(i) / h / 6.0;
        s_c(i) = y_samples(i + 1) / h - m(i + 1) * h / 6.0 ;
        s_d(i) = y_samples(i) / h - m(i) * h / 6.0 ;
    }
}
double CubicSplineInterp1D::CalcVal(double x)
{
    int i = 0;
    double r = 0.0;

    //determine x in range (i,i+1)
    for(int j = 0; j < samples - 1; j++)
    {
        if(x >= x_samples(j))
        {
            i = j;
        }
        else
        {
            break;
        }
    }
    r = s_a(i)   * pow( x - x_samples(i),     3)
        + s_b(i) * pow( x_samples(i + 1) - x, 3)
        + s_c(i) * ( x - x_samples(i)     )
        + s_d(i) * ( x_samples(i + 1) - x );
    return r;
}
