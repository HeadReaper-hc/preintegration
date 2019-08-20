#include <Eigen/Dense>
#include <iostream>
#include <common.h>
#include <boost/shared_ptr.hpp>
#include <ceres/ceres.h>
#include "glog/logging.h"
using namespace Eigen;
using namespace std;
using namespace hc;
using namespace ceres;
int main(int argc,char** argv){
    //google::InitGoogleLogging(argv[0]);

    boost::shared_ptr<PreintegrationParams> param(new PreintegrationParams());

    Matrix3d GyrCov;
    GyrCov << 0.001,0,0,
            0,0.001,0,
            0,0,0.001;
    Matrix3d AccCov;
    AccCov << 0.001,0,0,
            0,0.001,0,
            0,0,0.001;
    param->setGyrCov(GyrCov);
    param->setAccCov(AccCov);

    Vector3d bias_a(0,0,0);
    Vector3d bias_g(0,0,0);

    ManifoldPreIntegration temp(bias_a,bias_g,param);
    for(int i=0;i<20;i++){
        temp.update(Vector3d(1,0,9.8),Vector3d(0,0,0),0.01);
    }

//    Quaterniond Ri(1,0,0,0);
//    Vector3d Pi(0,0,0);
//    Vector3d Vi(0,0,0);
//
//    Quaterniond Rj(1,0,0,0);
//    Vector3d Pj(0.10,0.005,-0.003);
//    Vector3d Vj(5,0.01,-0.01);

    double param0[10] = {1,0,0,0,0,0,0,0,0,0};
    double param1[10] = {1,0,0,0,0.024,0.005,-0.003,0.22,0.01,-0.01};
    double param2[6] = {0,0,0,0,0,0};

    Problem problem;
    ceres::LocalParameterization* local_parameterization = new ImuStateParameterization();
    problem.AddParameterBlock(param0,10,local_parameterization);

    ceres::LocalParameterization* local_parameterization1 = new ImuStateParameterization();
    problem.AddParameterBlock(param1,10,local_parameterization1);

    problem.AddParameterBlock(param2,6);

    CostFunction* cost_function = new ImuFactor(&temp,Vector3d(0,0,-9.8));
    problem.AddResidualBlock(cost_function,NULL,param0,param1,param2);

    Solver::Options options;
    options.minimizer_progress_to_stdout = true;
    Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    std::cout <<summary.FullReport()<<"\n";


    cout<<"DeltaRij : "<<temp.DeltaRij().toRotationMatrix()<<endl;

    cout<<"DeltaVij : "<<temp.DeltaVij()<<endl;

    cout<<"DeltaPij : "<<temp.DeltaPij()<<endl;

    cout<<"Covriance : "<<temp.COVariance()<<endl;

    cout<<"DeltaTij : "<<temp.DeltaTij()<<endl;

    cout<<"param0:"<<endl;
    for(int i=0;i<10;i++) {
        cout<<param0[i]<<" ";
    }
    cout<<endl;
    cout<<"param1:"<<endl;
    for(int i=0;i<10;i++) {
        cout<<param1[i]<<" ";
    }
    cout<<endl;
    cout<<"bias:"<<endl;
    for(int i=0;i<6;i++){
        cout<<param2[i]<<" ";
    }


}
