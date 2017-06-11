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
    
    number phi;
    number alpha;
    number theta;
    
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
	
    // Random configuration generators
    inline void SetRandomAxis(std::mt19937_64& rng_engine, std::uniform_real_distribution<number>& rng_distrib)
    {
        phi   = rng_distrib(rng_engine) * 2.*PI;
        theta = rng_distrib(rng_engine) * PI;
        
        Axis  << sin(theta)*cos(phi),
                 sin(theta)*sin(phi),
                 cos(theta);
		
        // Assume the main particle axis of all base configurations is initially borne by ez
        Hull->Axis = Axis;
    }
    
    inline void SetRandomOrientation(std::mt19937_64& rng_engine, std::uniform_real_distribution<number>& rng_distrib)
    {
        alpha = rng_distrib(rng_engine) * 2.*PI;
        
        // Euler rotation matrix - ZYZ convention
        Orientation << cos(alpha)*cos(theta)*cos(phi) - sin(alpha)*sin(phi),
                      -sin(alpha)*cos(theta)*cos(phi) - cos(alpha)*sin(phi),
                       sin(theta)*cos(phi),
                       cos(alpha)*cos(theta)*sin(phi) + sin(alpha)*cos(phi),
                      -sin(alpha)*cos(theta)*sin(phi) + cos(alpha)*cos(phi),
                       sin(theta)*sin(phi),
                      -cos(alpha)*sin(theta),
                       sin(alpha)*sin(theta),
                       cos(theta);
		
        Hull->Center      = Orientation * Hull->Center_p;
        Hull->Orientation = Orientation * Hull->Orientation_p;
    }
	
    virtual void Build(int) = 0;
    virtual void Parse(std::mt19937_64&) {Hull = &BVH.Forest[0];}
    
protected:
    int id_;
};

#endif
