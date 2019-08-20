//
// Created by huchao on 19-8-18.
//

#include <ManifoldPreIntegration.h>
#include <sophus/so3.hpp>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <math_utils.h>
#include <math.h>
using namespace Eigen;
using namespace Sophus;

namespace hc{

    void ManifoldPreIntegration::update(const Vector3d &acc, const Vector3d &gyr, const double &dt) {

        Eigen::Matrix<double,9,9> A(Eigen::Matrix<double,9,9>::Zero());

        Eigen::Matrix<double,9,6> B(Eigen::Matrix<double,9,6>::Zero());

        Vector3d bias_a = ProcessState->bias_a;
        Vector3d bias_g = ProcessState->bias_g;

        SO3d SO3_deltakk_1 = SO3d::exp((gyr-bias_g)*dt);
        Matrix3d  rot_deltakk_1 = SO3_deltakk_1.matrix();

        A.block<3,3>(0,0) = (rot_deltakk_1).transpose();
        A.block<3,3>(3,0) = -(ProcessState->deltaRij).toRotationMatrix() * hat(acc-bias_a) * dt;
        A.block<3,3>(6,0) = -0.5 * ProcessState->deltaRij.toRotationMatrix() * hat(acc-bias_a) * dt * dt;
        A.block<3,3>(3,3) = Matrix3d::Identity();
        A.block<3,3>(6,3) = Matrix3d::Identity() * dt;
        A.block<3,3>(6,6) = Matrix3d::Identity();

        Vector3d fai = gyr - bias_g;
        double fai_norm = fai.norm();
        Matrix3d Jr;
        if(fai_norm==0){
            Jr = Matrix3d::Identity() - 0.5 * hat(fai) + 1.0 / 6.0 * hat(fai) * hat(fai);
        }
        else {
            Jr = Matrix3d::Identity() - cos(fai_norm) / pow(fai_norm, 2) * hat(fai) +
                          (fai_norm - sin(fai_norm)) / pow(fai_norm, 3) * hat(fai) * hat(fai);
        }

        B.block<3,3>(0,0) = Jr * dt;
        B.block<3,3>(3,3) = ProcessState->deltaRij.toRotationMatrix() * dt;
        B.block<3,3>(6,3) = ProcessState->deltaRij.toRotationMatrix() * dt * dt;

        Matrix6d C(Matrix6d::Zero());
        C.block<3,3>(0,0) = p_->getGyrCov() / dt;
        C.block<3,3>(3,3) = p_->getAccCov() / dt;

        Covariance = A * Covariance * A.transpose() + B * C * B.transpose();

        //save old Rij, used in the bias jacobian update
        Quaterniond oldRij = ProcessState->deltaRij;

        ProcessState->deltaPij = ProcessState->deltaPij + ProcessState->deltaVij * dt +
                                + 0.5 * ProcessState->deltaRij.toRotationMatrix() * (acc-bias_a) * dt * dt;
        ProcessState->deltaVij = ProcessState->deltaVij + ProcessState->deltaRij * (acc-bias_a) * dt;
        ProcessState->deltaRij = ProcessState->deltaRij * Quaterniond(rot_deltakk_1);
        deltaTij += dt;

        //bias jacobian update
        D_p_ba += -1.5 * oldRij.toRotationMatrix() * dt * dt;
        D_p_bg += -1.5 * oldRij.toRotationMatrix() * hat(acc-bias_a) * D_r_bg * dt * dt;
        D_v_ba += -oldRij.toRotationMatrix() * dt;
        D_v_bg += -oldRij.toRotationMatrix() * hat(acc-bias_a) * D_r_bg * dt;
        D_r_bg = rot_deltakk_1.transpose() * D_r_bg - Jr * dt;
    }

    void ManifoldPreIntegration::corrertPreIntegration(const Vector3d &delta_bias_a, const Vector3d &delta_bias_g) {

        ProcessState->deltaPij = ProcessState->deltaPij + D_p_ba * delta_bias_a + D_p_bg * delta_bias_g;
        ProcessState->deltaVij = ProcessState->deltaVij + D_v_ba * delta_bias_a + D_v_bg * delta_bias_g;
        SO3d SO3_Rij_bias = SO3d::exp(D_r_bg*delta_bias_g);
        ProcessState->deltaRij = ProcessState->deltaRij * Quaterniond(SO3_Rij_bias.matrix());
    }
}
