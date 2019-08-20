//
// Created by huchao on 19-8-18. ref:GTSAM
//

#ifndef PREINTEGRATION_PREINTEGRATIONBASE_H
#define PREINTEGRATION_PREINTEGRATIONBASE_H

#include <PreintegrationParams.h>
#include <boost/shared_ptr.hpp>
#include <Eigen/Dense>
using namespace Eigen;
namespace hc{

    class PreintegrationBase{
    public:
        typedef PreintegrationParams Params;

    protected:
        boost::shared_ptr<Params> p_;

        //Time interval from i to j
        double deltaTij;

        //jacobian: pij/ba pij/bg vij/ba vij/bg rij/bg
        Matrix3d D_p_ba;
        Matrix3d D_p_bg;
        Matrix3d D_v_ba;
        Matrix3d D_v_bg;
        Matrix3d D_r_bg;

        //order: R(lie algebra), V, P
        Matrix<double,9,9> Covariance;

        PreintegrationBase() = delete;

        PreintegrationBase(const boost::shared_ptr<Params>& p):p_(p),deltaTij(0),
                                                               Covariance(Matrix<double,9,9>::Zero()){
            D_p_ba.setZero();
            D_p_bg.setZero();
            D_v_ba.setZero();
            D_v_bg.setZero();
            D_r_bg.setZero();
        }

        virtual ~PreintegrationBase(){}

    public:
        //Re-initialize PreintegratedMeasurements
        virtual void resetIntegration() =0;

        //return deltaPij
        virtual Vector3d DeltaPij() const =0;

        //return deltaRij
        virtual Quaterniond DeltaRij() const =0;

        //return deltaVij
        virtual Vector3d DeltaVij() const =0;

        //return bias_a
        virtual Vector3d Bias_a() const =0;

        //return bias_g
        virtual Vector3d Bias_g() const =0;

        //set bias_a
        virtual void setBias_a(const Vector3d& newBias_a) =0;

        //set bias_g;
        virtual void setBias_g(const Vector3d& newBias_g) =0;

        inline Matrix3d DPBA(){ return D_p_ba; }

        inline Matrix3d DPBG(){ return D_p_bg; }

        inline Matrix3d DVBA(){ return D_v_ba; }

        inline Matrix3d DVBG(){ return D_v_bg; }

        inline Matrix3d DRBG(){ return D_r_bg; }

        //return params
        inline boost::shared_ptr<Params>& params() {
            return p_;
        }

        //return deltaTij
        inline double DeltaTij() const { return deltaTij; }

        //return Covariance
        inline Matrix<double,9,9> COVariance() const { return Covariance; }

        //update preintegrated measurements
        virtual void update(const Vector3d& acc, const Vector3d& gyr, const double& dt) = 0;

        //corrert the preIntegration by bias
        virtual void corrertPreIntegration(const Vector3d& delta_bias_a, const Vector3d& delta_bias_g) =0;
    };  //end class PreIntegrationBase
} //end namespace hc


#endif //PREINTEGRATION_PREINTEGRATIONBASE_H
