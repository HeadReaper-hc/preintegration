//
// Created by huchao on 19-8-19.
//

#ifndef PREINTEGRATION_MATH_UTILS_H
#define PREINTEGRATION_MATH_UTILS_H
#include <Eigen/Dense>
using namespace Eigen;

namespace hc{

    Matrix3d hat(const Vector3d& x);


    Matrix3d computerJr(const Vector3d& fai);


}

#endif //PREINTEGRATION_MATH_UTILS_H
