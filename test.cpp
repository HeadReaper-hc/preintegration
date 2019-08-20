//
// Created by huchao on 19-8-19.
//
#include <Eigen/Dense>
#include <iostream>
#include <common.h>
#include <boost/shared_ptr.hpp>
using namespace Eigen;
using namespace std;
using namespace hc;

int main(){
    boost::shared_ptr<PreintegrationParams> param(new PreintegrationParams());

    Matrix3d GyrCov;
    GyrCov << 0.1,0,0,
              0,0.1,0,
              0,0,0.1;
    Matrix3d AccCov;
    AccCov << 0.1,0,0,
              0,0.1,0,
              0,0,0.1;
    param->setGyrCov(GyrCov);
    param->setAccCov(AccCov);

    Vector3d bias_a(0.1,0,0);
    Vector3d bias_g(0,0,0);

    ManifoldPreIntegration temp(bias_a,bias_g,param);
    for(int i=0;i<20;i++){
        temp.update(Vector3d(1,0,0),Vector3d(0,0,0),0.01);
    }

    temp.corrertPreIntegration(Vector3d(-0,0,0),Vector3d(0,0,0));

    cout<<"DeltaRij : "<<temp.DeltaRij().toRotationMatrix()<<endl;

    cout<<"DeltaVij : "<<temp.DeltaVij()<<endl;

    cout<<"DeltaPij : "<<temp.DeltaPij()<<endl;

    cout<<"Covriance : "<<temp.COVariance()<<endl;

    cout<<"DeltaTij : "<<temp.DeltaTij()<<endl;

}

