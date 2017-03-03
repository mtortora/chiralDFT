#ifndef DNA_DUPLEX_HPP_
#define DNA_DUPLEX_HPP_

#include "BaseParticle.hpp"


class DNADuplex: public BaseParticle
{
private:
    uint   N_CONF_;
    uint   idx_conf_;

    double V_CYT_;
    double V_GUA_;
    double V_ADE_;
    double V_THY_;
    double V_BCK_;
    
protected:
    uint   N_STAR_;
    
    double EXCL_EPS_;
    double EXCL_S1_;
    double EXCL_R1_;
    double EXCL_B1_;
    double EXCL_RC1_;
    double BETA_R_;
    double E_CUT_;
    double DELTA_SCREEN_;
    double COULOMB_FACTOR_;
    double EPSILON_DNA_;
    double EPSILON_WATER_;
    double LAMBDA_;
    double MINUS_KAPPA_;
    double PREFACTOR_;
    double DH_PREFACTOR_;
    double OFFSET_;
    double DH_OFFSET_;
    double B_CUT_;
    double R_CUT_;
    double R_STAR_;
    double DELTA_R_;
    
public:
    DNADuplex();

    // Load random configuration from BHierarchy->Forest
    void Parse(std::mt19937_64& rng_engine) override
    {
        idx_conf_ = rng_engine() % N_CONF_;
        BHull     = &BHierarchy->Forest[idx_conf_];
    }
        
    void Build(int) override;
    
    virtual ~DNADuplex() {}
};

#endif
