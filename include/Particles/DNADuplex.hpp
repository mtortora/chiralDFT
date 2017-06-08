#ifndef DNA_DUPLEX_HPP_
#define DNA_DUPLEX_HPP_

#include "BaseParticle.hpp"


template<typename number>
class DNADuplex: public BaseParticle<number>
{
public:
    DNADuplex();
    
    // Load random configuration from BHierarchy->Forest
    void Parse(std::mt19937_64& rng_engine) override
    {
        uint idx_conf = rng_engine() % N_CONF_;
        this->Hull    = &(this->BVH).Forest[idx_conf];
    }
    
    void Build(int) override;
  
protected:
    uint   N_STAR_;

    number EXCL_EPS_;
    number EXCL_B1_;
    number BETA_R_;
    number E_CUT_;
    
    number DELTA_SCREEN_;
    number COULOMB_FACTOR_;
    number EPSILON_DNA_;
    number EPSILON_WATER_;
    number PREFACTOR_;
    number DH_PREFACTOR_;
    number OFFSET_;
    number DH_OFFSET_;
    number B_CUT_;
    
    number EXCL_S1_;
    number EXCL_R1_;
    number EXCL_RC1_;
    number LAMBDA_;
    number MINUS_KAPPA_;
    number R_CUT_;
    number R_STAR_;
    number DELTA_R_;
    
private:
    uint   N_CONF_;

    number V_CYT_;
    number V_GUA_;
    number V_ADE_;
    number V_THY_;
    number V_BCK_;
};

#endif
