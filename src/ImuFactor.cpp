//
// Created by huchao on 19-8-20.
//

#include <ImuFactor.h>
#include <Eigen/Dense>
#include <sophus/so3.hpp>
#include <math_utils.h>
using namespace Eigen;
using namespace Sophus;
namespace hc{

     bool ImuFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const {

        Quaterniond Ri(parameters[0][0],parameters[0][1],parameters[0][2],parameters[0][3]);
        Vector3d Pi(parameters[0][4],parameters[0][5],parameters[0][6]);
        Vector3d Vi(parameters[0][7],parameters[0][8],parameters[0][9]);

        Quaterniond Rj(parameters[1][0],parameters[1][1],parameters[1][2],parameters[1][3]);
        Vector3d Pj(parameters[1][4],parameters[1][5],parameters[1][6]);
        Vector3d Vj(parameters[1][7],parameters[1][8],parameters[1][9]);

        Vector3d ba_i(parameters[2][0],parameters[2][1],parameters[2][2]);
        Vector3d bg_i(parameters[2][3],parameters[2][4],parameters[2][5]);

        Map<Vector3d> delta_fai(residuals,3);
        Map<Vector3d> delta_p(residuals+3,3);
        Map<Vector3d> delta_v(residuals+6,3);

        Eigen::Matrix<double,9,9> sqrt_info = Eigen::LLT<Eigen::Matrix<double,9,9>>(PreIntegrateionResult->COVariance().inverse()).matrixL().transpose();

        Matrix3d delta_rot = PreIntegrateionResult->DeltaRij().toRotationMatrix().transpose() *
                             Ri.toRotationMatrix().transpose() * Rj.toRotationMatrix();
        SO3d deltaFai(delta_rot);
        delta_fai = deltaFai.log();

        Vector3d delta_bg = bg_i - PreIntegrateionResult->Bias_g();
        Vector3d delta_ba = ba_i - PreIntegrateionResult->Bias_a();

        PreIntegrateionResult->corrertPreIntegration(delta_ba,delta_bg);

        double tij = PreIntegrateionResult->DeltaTij();

        delta_v = Ri.toRotationMatrix().transpose() * (Vj - Vi - g*tij) - PreIntegrateionResult->DeltaVij();
        delta_p = Ri.toRotationMatrix().transpose() * (Pj - Pi - Vi*tij - 0.5*g*tij*tij) -
                    PreIntegrateionResult->DeltaPij();

        Eigen::Map<Eigen::Matrix<double,9,1> > residual(residuals);
        residual = sqrt_info * residual;

//         cout<<"delta_fai: "<<delta_fai<<endl;
//         cout<<"delta_p: "<<delta_p<<endl;
//         cout<<"delta_v: "<<delta_v<<endl;



        if(jacobians) {

            Matrix3d Jr,Jr2;
            Jr = computerJr(delta_fai);

            Vector3d D_r_bg_big = PreIntegrateionResult->DRBG()*delta_bg;
            Jr2 = computerJr(D_r_bg_big);

            if (jacobians[0]) {  //Res: fai,p,v Param: Ri,Pi,Vi
                Eigen::Map<Eigen::Matrix<double,9,10,Eigen::RowMajor> > jacobian0(jacobians[0]);
                jacobian0.setZero();
                jacobian0.block<3, 3>(0,0) = -Jr.inverse() * Rj.toRotationMatrix().transpose() *
                                               Ri.toRotationMatrix();
                jacobian0.block<3, 3>(3,0) = hat(Ri.toRotationMatrix().transpose() * (Pj - Pi - Vi * tij - 0.5 * g * tij * tij));
                jacobian0.block<3, 3>(3,4) = -Matrix3d::Identity();
                jacobian0.block<3, 3>(3,7) = -Ri.toRotationMatrix().transpose() * tij;
                jacobian0.block<3, 3>(6,0) = hat(Ri.toRotationMatrix().transpose() * (Vj - Vi - g * tij));
                jacobian0.block<3, 3>(6,7) = -Ri.toRotationMatrix().transpose();
                jacobian0 = sqrt_info * jacobian0;
            }

            if (jacobians[1]) {
                Eigen::Map<Eigen::Matrix<double,9,10,Eigen::RowMajor> > jacobian1(jacobians[1]);
                jacobian1.setZero();
                jacobian1.block<3, 3>(0,0) = Jr.inverse();
                jacobian1.block<3, 3>(3,4) = Ri.toRotationMatrix().transpose() * Rj.toRotationMatrix();
                jacobian1.block<3, 3>(6,7) = Ri.toRotationMatrix().transpose();
                jacobian1 = sqrt_info * jacobian1;
            }
            if (jacobians[2]) {
                Eigen::Map<Eigen::Matrix<double,9,6,Eigen::RowMajor> > jacobian2(jacobians[2]);
                jacobian2.setZero();
                jacobian2.block<3, 3>(0,3) = -Jr.inverse() * delta_rot.transpose() * Jr2 *
                                               PreIntegrateionResult->DRBG();
                jacobian2.block<3, 3>(3,0) = -PreIntegrateionResult->DPBA();
                jacobian2.block<3, 3>(3,3) = -PreIntegrateionResult->DPBG();
                jacobian2.block<3, 3>(6,0) = -PreIntegrateionResult->DVBA();
                jacobian2.block<3, 3>(6,3) = -PreIntegrateionResult->DVBG();
                jacobian2 = sqrt_info * jacobian2;
            }
        }

    }
}