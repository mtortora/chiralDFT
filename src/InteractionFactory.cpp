// ===================================================================
/**
 * Interaction Factory class.
 * Provides full template specialisation for all particle types
 */
// ===================================================================
/*
 * InteractionFactory.cpp: Version 1.0
 * Created 26/07/2016 by Maxime Tortora
 */
// ===================================================================

#include "InteractionFactory.hpp"

using namespace Eigen;


/* ************************* */
/* BentCore */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
double InteractionFactory<BentCore>::MayerInteraction(const Vector3d& R_cm,
                                                      BentCore* Particle1, BentCore* Particle2)
{
    double energy(0.);
    double mayer_interaction(0.);
    
    // In this case, energy counts the number of overlapping beads
    RecursiveInteraction(R_cm, Particle1->BHull, Particle2->BHull, &energy, 0.);
    
    if ( energy > 0. ) mayer_interaction = 1.;
    
    return mayer_interaction;
}


// ============================


/* ************************* */
/* DNADuplex */
/* ************************* */

// ============================
/* Repulsive Lennard-Jones (WCA) potential */
// ============================
double InteractionFactory<DNADuplex>::RepulsiveLJ_(double r, double sigma, double r_star, double b, double r_cut)
{
    double energy(0.);
    
    // Short-range LJ potential with quadratic cutoff for r_star < r < r_cut
    if( r < r_cut )
    {
        if( r < r_star )
        {
            double lj_part = SQR(CUB(sigma/r));
            energy         = 4. * EXCL_EPS_ * (SQR(lj_part) - lj_part);
        }
        
        else energy = EXCL_EPS_ * b * SQR(r - r_cut);
    }
    
    return energy;
}

// ============================
/* Debye-Huckel interaction potential */
// ============================
double InteractionFactory<DNADuplex>::DebyeHuckel_(double r)
{
    double energy(0.);
    
    if ( MODE_DH == DH_OXDNA )
    {
        // Debye-Huckel with quadratic cutoff for R_STAR_ < r < R_CUT_
        if ( r < R_CUT_ )
        {
            if ( r < R_STAR_ ) energy = exp(r * MINUS_KAPPA_) * (DH_PREFACTOR_ / r);
            else               energy = B_CUT_ * SQR(r - R_CUT_);
        }
    }
    
    if ( MODE_DH == DH_FERRARINI )
    {
        // Short-range Coulomb + long-range DH with linear interpolation for EXCL_S1_ < r < R_STAR_
        if ( r < R_CUT_ )
        {
            if      ( r < EXCL_S1_ ) energy = (PREFACTOR_ / r);
            else if ( r < R_STAR_ )  energy = (DH_OFFSET_ - OFFSET_) * (r - EXCL_S1_)/(R_STAR_ - EXCL_S1_) + OFFSET_;
            else                     energy = (DH_PREFACTOR_ / r) * exp(r * MINUS_KAPPA_);
        }
    }
    
    return energy;
}

// ============================
/* Mayer interaction function */
// ============================
double InteractionFactory<DNADuplex>::MayerInteraction(const Vector3d& R_cm,
                                                       DNADuplex* Particle1, DNADuplex* Particle2)
{
    double energy(0.);
    
    // Work out pairwise interaction energies recursively
    RecursiveInteraction(R_cm, Particle1->BHull, Particle2->BHull, &energy, E_CUT_);
    
    return 1. - exp(-BETA_R_*energy);
}


// ============================


/* ************************* */
/* Helix */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
double InteractionFactory<Helix>::MayerInteraction(const Vector3d& R_cm,
                                                   Helix* Particle1, Helix* Particle2)
{
    double energy(0.);
    double mayer_interaction(0.);
    
    // In this case, energy counts the number of overlapping beads
    RecursiveInteraction(R_cm, Particle1->BHull, Particle2->BHull, &energy, 0.);
    
    if ( energy > 0. ) mayer_interaction = 1.;
    
    return mayer_interaction;
}


// ============================


/* ************************* */
/* PatchyRod */
/* ************************* */

// ============================
/* WCA interaction potential */
// ============================
double InteractionFactory<PatchyRod>::RepulsiveWCA_(double r)
{
    double energy(0.);
    
    if ( r < R_WCA_ )
    {
        double lj_part = SQR(CUB(SIGMA_R/r));
        energy         = 4. * EPSILON_WCA_ * (SQR(lj_part) - lj_part) + EPSILON_WCA_;
    }
    
    return energy;
}

// ============================
/* Debye-Huckel interaction potential */
// ============================
double InteractionFactory<PatchyRod>::DebyeHuckel_(double r)
{
    double energy(0.);
    
    if ( r < R_CUT_ ) energy = exp(r * MINUS_KAPPA_) * (DH_PREFACTOR_ / r);
    
    return energy;
}

// ============================
/* Mayer interaction function */
// ============================
double InteractionFactory<PatchyRod>::MayerInteraction(const Vector3d& R_cm,
                                                       PatchyRod* Particle1, PatchyRod* Particle2)
{
    double energy(0.);
    
    // Backbone-backbone WCA repulsion
    MatrixXd Backbone1 = Particle1->Orientation * Particle1->Backbone;
    MatrixXd Backbone2 = Particle2->Orientation * Particle2->Backbone;
    
    for ( uint idx_vtx1 = 0; idx_vtx1 < Backbone1.cols() && energy < E_CUT_; ++idx_vtx1 )
    {
        for ( uint idx_vtx2 = 0; idx_vtx2 < Backbone2.cols() && energy < E_CUT_; ++idx_vtx2 )
        {
            Vector3d R_sep = R_cm + Backbone2.col(idx_vtx2) - Backbone1.col(idx_vtx1);
            double   norm  = R_sep.norm();
            
            energy        += RepulsiveWCA_(norm);
        }
    }
    
    if ( energy < E_CUT_ )
    {
        MatrixXd Patches1 = Particle1->Orientation * Particle1->Patches;
        MatrixXd Patches2 = Particle2->Orientation * Particle2->Patches;
        
        if ( USE_DH )
        {
            // Assume pairwise-additivity for site-site electrostatic interactions
            for ( uint idx_vtx1 = 0; idx_vtx1 < Patches1.cols() && energy < E_CUT_; ++idx_vtx1 )
            {
                for ( uint idx_vtx2 = 0; idx_vtx2 < Patches2.cols() && energy < E_CUT_; ++idx_vtx2 )
                {
                    Vector3d R_sep = R_cm + Patches2.col(idx_vtx2) - Patches1.col(idx_vtx1);
                    double   norm  = R_sep.norm();
                    
                    energy        += DebyeHuckel_(norm);
                }
            }
        }
        
        else
        {
            // Hard spherical thread of radius SIGMA_R
            for ( uint idx_vtx1 = 0; idx_vtx1 < Patches1.cols() && energy < E_CUT_; ++idx_vtx1 )
            {
                for ( uint idx_vtx2 = 0; idx_vtx2 < Patches2.cols() && energy < E_CUT_; ++idx_vtx2 )
                {
                    Vector3d R_sep = R_cm + Patches2.col(idx_vtx2) - Patches1.col(idx_vtx1);
                    double   norm  = R_sep.norm();
                    
                    if ( norm < SIGMA_R ) return 1.;
                }
            }
        }
    }
    
    return 1. - exp(-energy);
}


// ============================


/* ************************* */
/* ThreadedRod */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
double InteractionFactory<ThreadedRod>::MayerInteraction(const Vector3d& R_cm,
                                                         ThreadedRod* Particle1, ThreadedRod* Particle2)
{
    if ( !USE_DH ) return 1.;
    
    double energy(0.);
    RecursiveInteraction(R_cm, Particle1->BHull, Particle2->BHull, &energy, E_CUT_);
    
    return 1. - exp(-energy);
}


// ============================


/* ************************* */
/* TwistedCuboid */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
double InteractionFactory<TwistedCuboid>::MayerInteraction(const Vector3d& R_cm,
                                                           TwistedCuboid* Particle1, TwistedCuboid* Particle2)
{
    double energy(0.);
    double mayer_interaction(0.);
    
    // In this case, energy counts the number of overlapping vertices
    RecursiveInteraction(R_cm, Particle1->BHull, Particle2->BHull, &energy, 0.);

    if ( energy > 0. ) mayer_interaction = 1.;
    
    return mayer_interaction;
}


// ============================


/* ************************* */
/* TwistedPentagon */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
double InteractionFactory<TwistedPentagon>::MayerInteraction(const Vector3d& R_cm,
                                                             TwistedPentagon* Particle1, TwistedPentagon* Particle2)
{
    double energy(0.);
    double mayer_interaction(0.);
    
    // In this case, energy counts the number of overlapping vertices
    RecursiveInteraction(R_cm, Particle1->BHull, Particle2->BHull, &energy, 0.);
    
    if ( energy > 0. ) mayer_interaction = 1.;
    
    return mayer_interaction;
}
