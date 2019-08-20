//
// Created by huchao on 19-8-20.
//
#include <ImuStateParameterization.h>
#include <Eigen/Dense>

using namespace Eigen;

namespace hc{

    bool ImuStateParameterization::Plus(const double *x, const double *delta, double *x_plus_delta) const {

/*    cout <<"delta : " <<delta[0]<<","<<delta[1]<<","<<delta[2]<<","<<delta[3]<<","<<delta[4]<<","<<delta[5]<<","<<delta[6]
         <<","<<delta[7]<<","<<delta[8]<<endl; */

        const double norm_delta =
                sqrt(delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2]);
        if (norm_delta > 0.0) {
            const double sin_delta_by_delta = (sin(norm_delta) / norm_delta);
            double q_delta[4];
            q_delta[0] = cos(norm_delta);
            q_delta[1] = sin_delta_by_delta * delta[0];
            q_delta[2] = sin_delta_by_delta * delta[1];
            q_delta[3] = sin_delta_by_delta * delta[2];
            QuaternionProduct(q_delta, x, x_plus_delta);



        }else {
            for (int i = 0; i < 4; ++i)
                x_plus_delta[i] = x[i];
        }

        for (int i = 4; i < 10; ++i)
            x_plus_delta[i] = x[i] + delta[ i-1 ];

        return true;

    }

    bool ImuStateParameterization::ComputeJacobian(const double *x, double *jacobian) const {

        Eigen::Map<MatrixXd> Jac(jacobian,9,10);
        Jac.setZero();
        Jac.block<3,3>(0,0) = Matrix3d::Identity();
        Jac.block<3,3>(3,4) = Matrix3d::Identity();
        Jac.block<3,3>(6,7) = Matrix3d::Identity();

        return true;
    }

}
