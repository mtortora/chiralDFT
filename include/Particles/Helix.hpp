#ifndef HELIX_HPP_
#define HELIX_HPP_

#include "BaseParticle.hpp"


class Helix: public BaseParticle
{
public:
    Helix();
    
    Eigen::Matrix3Xd Backbone;
    void Build(int) override;
    
protected:
    double D_HARD_;
    
private:
    uint   N_S_;
    uint   N_RES_;
    
    double L_CTR_;
    double R_HLX_;
    double P_HLX_;
    double L_Z_;
};

#endif
