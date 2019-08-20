//
// Created by huchao on 19-8-18.
//

#ifndef PREINTEGRATION_PREINTEGRATIONPARAMS_H
#define PREINTEGRATION_PREINTEGRATIONPARAMS_H
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
namespace hc {

    class PreintegrationParams {
    public:
        PreintegrationParams(){}
        ~PreintegrationParams(){}

        inline void setAccCov(Matrix3d& _accCov){
            accelerometerCovariance = _accCov;
        }
        inline Matrix3d getAccCov(){
            return accelerometerCovariance;
        }

        inline void setGyrCov(Matrix3d& _gyrCov){
            gyroscopeCovariance = _gyrCov;
        }
        inline Matrix3d getGyrCov(){
            return gyroscopeCovariance;
        }

        inline void setGravity(Vector3d& _gravity){
            w_gravity = _gravity;
        }
        inline Vector3d getGravity(){
            return w_gravity;
        }

    private:
        Matrix3d accelerometerCovariance; //< continuous-time "Covariance" of accelerometer
        Matrix3d gyroscopeCovariance; //< continuous-time "Covariance" of accelerometer
        Vector3d w_gravity;   //< gravity vector in world frame

    }; //end class PreintegrationParams

} //end namespace

#endif //PREINTEGRATION_PREINTEGRATIONPARAMS_H
