//
// Created by huchao on 19-8-18.
//

#ifndef PREINTEGRATION_MANIFOLDPREINTEGRATION_H
#define PREINTEGRATION_MANIFOLDPREINTEGRATION_H

#include <PreintegrationBase.h>
namespace hc{

    class ManifoldPreIntegration : public PreintegrationBase{
    public:
        typedef PreintegrationParams Params;
        struct DeltaImuState{
            DeltaImuState(const Vector3d& _bias_a,const Vector3d& _bias_g):
                            deltaPij(Vector3d(0,0,0)),
                            deltaVij(Vector3d(0,0,0)),
                            deltaRij(Quaterniond(1,0,0,0)){
                bias_a = _bias_a;
                bias_g = _bias_g;
            }

            Vector3d deltaPij;
            Vector3d deltaVij;
            Quaterniond deltaRij;
            Vector3d bias_a;
            Vector3d bias_g;
            void reset(){
                deltaPij = Vector3d(0,0,0);
                deltaVij = Vector3d(0,0,0);
                deltaRij = Quaterniond(1,0,0,0);
            }
        };

        ManifoldPreIntegration(const Vector3d& bias_a,const Vector3d& bias_g,
                               const boost::shared_ptr<Params>& p):PreintegrationBase(p){
            ProcessState = new DeltaImuState(bias_a,bias_g);
        }

        ~ManifoldPreIntegration(){ delete ProcessState; }

        virtual Vector3d DeltaPij() const override {
            return ProcessState->deltaPij;
        }

        virtual Vector3d DeltaVij() const override {
            return ProcessState->deltaVij;
        }

        virtual Quaterniond DeltaRij() const override {
            return ProcessState->deltaRij;
        }

        virtual Vector3d Bias_a() const override {
            return ProcessState->bias_a;
        }

        virtual Vector3d Bias_g() const override {
            return ProcessState->bias_g;
        }

        virtual void setBias_a(const Vector3d& newBias_a){
            ProcessState->bias_a = newBias_a;
        }

        virtual void setBias_g(const Vector3d& newBias_g){
            ProcessState->bias_g = newBias_g;
        }

        virtual void resetIntegration() override {
            ProcessState->reset();
            Covariance = Matrix<double,9,9>::Zero();
            deltaTij = 0;
            D_p_ba.setZero();
            D_p_bg.setZero();
            D_v_ba.setZero();
            D_v_bg.setZero();
            D_r_bg.setZero();
        }

        virtual void update(const Vector3d& acc, const Vector3d& gyr, const double& dt) override ;

        virtual void corrertPreIntegration(const Vector3d& delta_bias_a, const Vector3d& delta_bias_g) override ;

    private:
        DeltaImuState* ProcessState;
    };

    typedef struct ManifoldPreIntegration::DeltaImuState ImuState;
}

#endif //PREINTEGRATION_MANIFOLDPREINTEGRATION_H
