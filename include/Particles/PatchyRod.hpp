#ifndef PATCHY_ROD_HPP_
#define PATCHY_ROD_HPP_

#include "BaseParticle.hpp"


class PatchyRod: public BaseParticle
{
private:
    uint   N_PATCH_;
    uint   N_BACK_;
    uint   N_RES_;
    
    double L_Z_;
    double P_PATCH_;
    
protected:
    double R_WCA_;
    double EPSILON_WCA_;
    double R_CUT_;
    double E_CUT_;
    double DH_PREFACTOR_;
    double MINUS_KAPPA_;
    
public:
    PatchyRod();
    
    Eigen::Matrix3Xd Backbone;
    Eigen::Matrix3Xd Patches;
    
    void Build(int) override;
    void Parse(std::mt19937_64&) override {BHull = BHierarchy;}
    
    virtual ~PatchyRod() {}
};

#endif
