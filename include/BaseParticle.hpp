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
        
        U             = Vector3<number>::UnitZ();
        V             = Vector3<number>::UnitX();
        W             = Vector3<number>::UnitY();
        
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

    Vector3<number>  U;
    Vector3<number>  V;
    Vector3<number>  W;

    Matrix33<number> Orientation;
    
    // Bounding hull and hierarchy
    BTree     <number>* Hull;
    BHierarchy<number>  BVH;
    
    // Axis setters - U, V, W are the respective long, medium, short axes
    inline void SetU(number theta_, number phi_)
    {
        U <<  sin(theta_)*cos(phi_),
              sin(theta_)*sin(phi_),
              cos(theta_);
    }
    
    inline void SetV(number alpha_, number theta_, number phi_)
    {
        V <<  cos(alpha_)*cos(theta_)*cos(phi_) - sin(alpha_)*sin(phi_),
              cos(alpha_)*cos(theta_)*sin(phi_) + sin(alpha_)*cos(phi_),
             -cos(alpha_)*sin(theta_);
    }
    
    inline void SetW(number alpha_, number theta_, number phi_)
    {
        W << -sin(alpha_)*cos(theta_)*cos(phi_) - cos(alpha_)*sin(phi_),
             -sin(alpha_)*cos(theta_)*sin(phi_) + cos(alpha_)*cos(phi_),
              sin(alpha_)*sin(theta_);
    }
	
    // Random configuration generators
    inline void SetRandomAxis(std::mt19937_64& rng_engine, std::uniform_real_distribution<number>& rng_distrib)
    {
        phi   = rng_distrib(rng_engine) * 2.*PI;
        theta = rng_distrib(rng_engine) * PI;
        
        SetU(theta, phi);
		
        // Assume the long particle axis of all base configurations is initially borne by z
        Hull->Axis = U;
    }
    
    inline void SetRandomOrientation(std::mt19937_64& rng_engine, std::uniform_real_distribution<number>& rng_distrib)
    {
        alpha = rng_distrib(rng_engine) * 2.*PI;
        
        SetV(alpha, theta, phi);
        SetW(alpha, theta, phi);

        // Euler rotation matrix - ZYZ convention
        Orientation.col(0) = V;
        Orientation.col(1) = W;
        Orientation.col(2) = U;
		
        Hull->Center       = Orientation * Hull->Center_p;
        Hull->Orientation  = Orientation * Hull->Orientation_p;
    }
	
    virtual void Build(int) = 0;
    virtual void Parse(std::mt19937_64&) {Hull = &BVH.Forest[0];}
    
protected:
    int id_;
};

#endif
