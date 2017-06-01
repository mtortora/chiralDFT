#ifndef BASE_PARTICLE_HPP
#define BASE_PARTICLE_HPP

// ===================================================================
/**
 * Header-only BaseParticle class, from which all particles types must
 * inherit. Implements the random axis/orientation generators.
 */
// ===================================================================
/*
 * BaseParticle.hpp:  Version 1.0
 * Created 22/06/2016 by Maxime Tortora
 */
// ===================================================================

#include <random>
#include <string>

#include "BHierarchy.hpp"


class BaseParticle
{
public:
    BaseParticle()
    {
        static int id = 0;
        id_           = id++;
        
        Axis          = Eigen::Vector3d::UnitZ();
        Orientation   = Eigen::Matrix3d::Identity();
    }

    uint   N_DELTA_L;
    
    double SIGMA_R = 1.;
    double R_INTEG;
    double V_INTEG;
    double V_EFF;
    double V0;

    Eigen::Vector3d Axis;
    Eigen::Matrix3d Orientation;
    
    // Bounding hull and hierarchy
    BTree*     Hull;
    BHierarchy BVH;
	
    double GetTheta() {return theta_;}

    // Random configuration generators
    inline void SetRandomAxis(std::mt19937_64& rng_engine, std::uniform_real_distribution<double>& rng_distrib)
    {
        phi_   = rng_distrib(rng_engine) * 2.*PI;
        theta_ = rng_distrib(rng_engine) * PI;
        
        Axis  << sin(theta_)*cos(phi_),
                 sin(theta_)*sin(phi_),
                 cos(theta_);
		
        // Assume the main particle axis of all base configurations is initially borne by ez
        Hull->Axis = Axis;
    }
    
    inline void SetRandomOrientation(std::mt19937_64& rng_engine, std::uniform_real_distribution<double>& rng_distrib)
    {
        alpha_ = rng_distrib(rng_engine) * 2.*PI;
        
        // Euler rotation matrix - ZYZ convention
        Orientation << cos(alpha_)*cos(theta_)*cos(phi_) - sin(alpha_)*sin(phi_),
                      -sin(alpha_)*cos(theta_)*cos(phi_) - cos(alpha_)*sin(phi_),
                       sin(theta_)*cos(phi_),
                       cos(alpha_)*cos(theta_)*sin(phi_) + sin(alpha_)*cos(phi_),
                      -sin(alpha_)*cos(theta_)*sin(phi_) + cos(alpha_)*cos(phi_),
                       sin(theta_)*sin(phi_),
                      -cos(alpha_)*sin(theta_),
                       sin(alpha_)*sin(theta_),
                       cos(theta_);
		
        Hull->Center      = Orientation * Hull->Center_p;
        Hull->Orientation = Orientation * Hull->Orientation_p;
    }
	
    virtual void Build(int) = 0;
    virtual void Parse(std::mt19937_64&) {Hull = &BVH.Forest[0];}
    
protected:
    int id_;
    
private:
    double phi_;
    double alpha_;
    double theta_;
};

#endif
