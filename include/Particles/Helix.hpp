#ifndef HELIX_HPP_
#define HELIX_HPP_

#include "BaseParticle.hpp"


template<typename number>
class Helix: public BaseParticle<number>
{
public:
    Helix();
    
    Matrix3X<number> Backbone;
    void Build(int) override;
    
protected:
    number D_HARD_;
    
private:
    uint   N_S_;
    uint   N_RES_;
    
    number L_CTR_;
    number R_HLX_;
    number P_HLX_;
    number L_Z_;
};

#endif
