#ifndef BASE_PARTICLE_HPP
#define BASE_PARTICLE_HPP

// ===================================================================
/**
 * Contains the BaseParticle class, from which all particles types
 * must inherit. Also implements random axis/orientation generators
 */
// ===================================================================
/*
 * BaseParticle.hpp:  Version 1.0
 * Created 22/06/2016 by Maxime Tortora
 */
// ===================================================================

#include <random>
#include <string>

#include "BTree.hpp"


class BaseParticle
{
private:
    double phi_;
    double alpha_;
    double theta_;
	
protected:
    int id_;
	
public:
    BaseParticle()
    {
        static int id = 0;

        BHierarchy    = new BTree();
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
    
    // Bounding hull + tree hierarchy
    BTree* BHull;
    BTree* BHierarchy;
	
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
        BHull->Axis = Axis;
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
		
        BHull->Center      = Orientation * BHull->Center_p;
        BHull->Orientation = Orientation * BHull->Orientation_p;
    }
	
    virtual void Build(int) = 0;
    virtual void Parse(std::mt19937_64&) {BHull = &BHierarchy->Forest[0];}
	
    virtual ~BaseParticle() {delete BHierarchy;}
};

#endif
