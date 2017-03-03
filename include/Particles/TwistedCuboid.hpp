#ifndef TWISTED_CUBOID_HPP_
#define TWISTED_CUBOID_HPP_

#include "rapid2/RAPID.H"
#include "BaseParticle.hpp"


class TwistedCuboid: public BaseParticle
{
private:
    uint   N_X_;
    uint   N_Y_;
    uint   N_Z_;

    double L_X_;
    double L_Y_;
    double L_Z_;
    double TWIST_;
    double R_BCK_;
    double P_BCK_;
    
    void SaveWireframe(const Eigen::Matrix3Xd&);
    
protected:
    double R_THRESHOLD_;

public:
    TwistedCuboid();

    RAPID_model* Mesh;

    void Build(int) override;
    void Tesselate(const Eigen::Matrix3Xd&, uint);

#if (USE_RAPID)
    void Parse(std::mt19937_64&) override {}
#endif
    
    virtual ~TwistedCuboid() {delete Mesh;}
};

#endif
