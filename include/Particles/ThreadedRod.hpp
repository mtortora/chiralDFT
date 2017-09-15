#ifndef THREADED_ROD_HPP_
#define THREADED_ROD_HPP_

#include "BaseParticle.hpp"


template<typename number>
class ThreadedRod: public BaseParticle<number>
{
public:
    ThreadedRod();
    
    Matrix3X<number> Patches;

    void Build(int) override;    
    void Parse(std::mt19937_64&) override {}
    
protected:
    number E_CUT_;
    number DH_PREFACTOR_;

    number D_HARD_;
    number L_Z_;
    number R_CUT_;
    number MINUS_KAPPA_;
    
private:
    uint   N_PATCH_;
    uint   N_RES_;
    
    number P_PATCH_;
    number R_PATCH_;
};

#endif
