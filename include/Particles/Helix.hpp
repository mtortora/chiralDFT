#ifndef HELIX_HPP_
#define HELIX_HPP_

#include "BaseParticle.hpp"


class Helix: public BaseParticle
{
private:
    uint   N_S_;
    uint   N_RES_;
    
    double L_CTR_;
    double R_HLX_;
    double P_HLX_;
    double L_X_;
    double L_Y_;
    double L_Z_;
    
protected:
    double R_HARD_;
        
public:
    Helix();
    
    Eigen::Matrix3Xd Backbone;
    void Build(int) override;
    
    virtual ~Helix() {}
};

#endif
