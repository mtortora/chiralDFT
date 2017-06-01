#ifndef THREADED_ROD_HPP_
#define THREADED_ROD_HPP_

#include "BaseParticle.hpp"


class ThreadedRod: public BaseParticle
{
public:
    ThreadedRod();
    
    void Build(int) override;
    
#if (!USE_DH)
    void Parse(std::mt19937_64&) override {}
#endif
    
protected:
    double D_HARD_;
    double L_Z_;
    double R_CUT_;
    double E_CUT_;
    double DH_PREFACTOR_;
    double MINUS_KAPPA_;
    
private:
    uint   N_PATCH_;
    uint   N_RES_;
    
    double P_PATCH_;
    double R_PATCH_;
};

#endif
