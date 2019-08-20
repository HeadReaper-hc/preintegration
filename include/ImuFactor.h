//
// Created by huchao on 19-8-20.
//

#ifndef PREINTEGRATION_IMUFACTOR_H
#define PREINTEGRATION_IMUFACTOR_H

#include <ceres/ceres.h>
#include <ManifoldPreIntegration.h>
namespace hc {

class ImuFactor : public ceres::SizedCostFunction<9,10,10,6>{   //{9:residuals 10:Pi,Vi,Ri 10:Pj,Vj,Rj 6:ba_i,bg_i
public:
    ImuFactor(ManifoldPreIntegration* _PreIntegrationResult ,Vector3d _g):PreIntegrateionResult(_PreIntegrationResult),g(_g){}
    virtual bool Evaluate(double const *const *parameters, double *residuals,
                          double **jacobians) const;

private:
    ManifoldPreIntegration* PreIntegrateionResult;
    Vector3d g;
};

};

#endif //PREINTEGRATION_IMUFACTOR_H
