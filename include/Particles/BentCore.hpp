#ifndef BENT_CORE_HPP_
#define BENT_CORE_HPP_

#include "BaseParticle.hpp"


class BentCore: public BaseParticle
{
public:
    BentCore();
    
    Eigen::Matrix3Xd Backbone;
    void Build(int) override;
    
protected:
    double D_HARD_;
    
private:
    uint   N_S_;
    uint   N_RES_;
    
    double L_CTR_;
    double R_BTC_;
    double GAMMA_;
    double CHI_;
    double L_X_;
    double L_Y_;
    double L_Z_;
};

#endif
