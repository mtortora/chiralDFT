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


template<typename number>
class BaseParticle
{
public:
    BaseParticle()
    {
        static int id = 0;
        id_           = id++;
        
        Axis          = Vector3<number>::UnitZ();
        Orientation   = Matrix33<number>::Identity();
    }

    uint   N_DELTA_L;
    
    number SIGMA_R = 1.;
    number R_INTEG;
    number V_INTEG;
    number V_EFF;
    number V0;

    Vector3<number> Axis;
    Matrix33<number> Orientation;
    
    // Bounding hull and hierarchy
    BTree<number>*     Hull;
    BHierarchy<number> BVH;
	
    number GetTheta() {return theta_;}

    // Random configuration generators
    inline void SetRandomAxis(std::mt19937_64& rng_engine, std::uniform_real_distribution<number>& rng_distrib)
    {
        phi_   = rng_distrib(rng_engine) * 2.*PI;
        theta_ = rng_distrib(rng_engine) * PI;
        
        Axis  << sin(theta_)*cos(phi_),
                 sin(theta_)*sin(phi_),
                 cos(theta_);
		
        // Assume the main particle axis of all base configurations is initially borne by ez
        Hull->Axis = Axis;
    }
    
    inline void SetRandomOrientation(std::mt19937_64& rng_engine, std::uniform_real_distribution<number>& rng_distrib)
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
    number phi_;
    number alpha_;
    number theta_;
};

#endif
