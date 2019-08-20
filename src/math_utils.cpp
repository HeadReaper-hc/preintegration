//
// Created by huchao on 19-8-20.
//
#include <math_utils.h>

namespace hc{
    Matrix3d hat(const Vector3d& x){
        Matrix3d result;
        result << 0, -x(2), x(1),
                x(2), 0, -x(0),
                -x(1), x(0), 0;
        return result;
    }

    Matrix3d computerJr(const Vector3d& fai){
        double fai_norm = fai.norm();
        Matrix3d jr;
        if(fai_norm==0) {
            jr = Matrix3d::Identity() - 0.5 * hat(fai) + 1.0 / 6.0 * hat(fai) * hat(fai);
        }
        else{
            jr = Matrix3d::Identity() - (1-cos(fai_norm))/pow(fai_norm,2) * hat(fai) +
                 (fai_norm-sin(fai_norm))/pow(fai_norm,3) * hat(fai) * hat(fai);
        }
        return jr;

    }
}
