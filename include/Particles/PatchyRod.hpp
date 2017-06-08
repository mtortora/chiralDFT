#ifndef PATCHY_ROD_HPP_
#define PATCHY_ROD_HPP_

#include "BaseParticle.hpp"


template<typename number>
class PatchyRod: public BaseParticle<number>
{
public:
    PatchyRod();
    
    Matrix3X<number> Backbone;
    Matrix3X<number> Patches;
    
    void Build(int) override;
    void Parse(std::mt19937_64&) override {}
    
protected:
    number EPSILON_WCA_;
    number E_CUT_;
    number DH_PREFACTOR_;
    
    number R_WCA_;
    number R_CUT_;
    number MINUS_KAPPA_;
    
private:
    uint   N_PATCH_;
    uint   N_BACK_;
    uint   N_RES_;
    
    number L_Z_;
    number P_PATCH_;
};

#endif
