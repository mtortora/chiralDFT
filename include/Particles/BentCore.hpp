#ifndef BENT_CORE_HPP_
#define BENT_CORE_HPP_

#include "BaseParticle.hpp"


template<typename number>
class BentCore: public BaseParticle<number>
{
public:
    BentCore();
    
    Matrix3X<number> Backbone;
    void Build(int) override;
    
protected:
    number D_HARD_;
    
private:
    uint   N_S_;
    uint   N_RES_;
    
    number L_CTR_;
    number R_BTC_;
    number GAMMA_;
    number CHI_;
    number L_X_;
    number L_Y_;
    number L_Z_;
};

#endif
