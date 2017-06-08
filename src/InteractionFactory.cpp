// ===================================================================
/**
 * Interaction Factory class.
 * Provides full template specialisation for all particle types.
 */
// ===================================================================
/*
 * InteractionFactory.cpp: Version 1.0
 * Created 26/07/2016 by Maxime Tortora
 */
// ===================================================================

#include "InteractionFactory.hpp"


/* ************************* */
/* BentCore */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<BentCore<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                      BentCore<number>* Particle1, BentCore<number>* Particle2)
{
    number energy(0.);
    number mayer_interaction(0.);
    
    // In this case, energy counts the number of overlapping beads
    this->RecursiveInteraction(R_cm, Particle1->Hull, Particle2->Hull, &energy, 0.);
    
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
template<typename number>
number InteractionFactory<DNADuplex<number>, number>::RepulsiveLJ_(number r, number sigma, number r_star, number b, number r_cut)
{
    number energy(0.);
    
    // Short-range LJ potential with quadratic cutoff for r_star < r < r_cut
    if( r < r_cut )
    {
        if( r < r_star )
        {
            number lj_part = SQR(CUB(sigma/r));
            energy         = 4. * this->EXCL_EPS_ * (SQR(lj_part) - lj_part);
        }
        
        else energy = this->EXCL_EPS_ * b * SQR(r - r_cut);
    }
    
    return energy;
}

// ============================
/* Debye-Huckel interaction potential */
// ============================
template<typename number>
number InteractionFactory<DNADuplex<number>, number>::DebyeHuckel_(number r)
{
    number energy(0.);
    
    if ( MODE_DH == DH_OXDNA )
    {
        // Debye-Huckel with quadratic cutoff for R_STAR_ < r < R_CUT_
        if ( r < this->R_CUT_ )
        {
            if ( r < this->R_STAR_ ) energy = exp(r * this->MINUS_KAPPA_) * (this->DH_PREFACTOR_ / r);
            else                     energy = this->B_CUT_ * SQR(r - this->R_CUT_);
        }
    }
    
    if ( MODE_DH == DH_FERRARINI )
    {
        // Short-range Coulomb + long-range DH with linear interpolation for EXCL_S1_ < r < R_STAR_
        if ( r < this->R_CUT_ )
        {
            if      ( r < this->EXCL_S1_ ) energy = (this->PREFACTOR_ / r);
            else if ( r < this->R_STAR_ )  energy = (this->DH_OFFSET_ - this->OFFSET_) * (r - this->EXCL_S1_)/(this->R_STAR_ - this->EXCL_S1_) + this->OFFSET_;
            else                           energy = (this->DH_PREFACTOR_ / r) * exp(r * this->MINUS_KAPPA_);
        }
    }
    
    return energy;
}

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<DNADuplex<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                       DNADuplex<number>* Particle1, DNADuplex<number>* Particle2)
{
    number energy(0.);
    
    // Work out pairwise interaction energies recursively
    this->RecursiveInteraction(R_cm, Particle1->Hull, Particle2->Hull, &energy, this->E_CUT_);
    
    return 1. - exp(-this->BETA_R_*energy);
}


// ============================


/* ************************* */
/* Helix */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<Helix<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                   Helix<number>* Particle1, Helix<number>* Particle2)
{
    number energy(0.);
    number mayer_interaction(0.);
    
    // In this case, energy counts the number of overlapping beads
    this->RecursiveInteraction(R_cm, Particle1->Hull, Particle2->Hull, &energy, 0.);
    
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
template<typename number>
number InteractionFactory<PatchyRod<number>, number>::RepulsiveWCA_(number r)
{
    number energy(0.);
    
    if ( r < this->R_WCA_ )
    {
        number lj_part = SQR(CUB(this->SIGMA_R/r));
        energy         = 4. * this->EPSILON_WCA_ * (SQR(lj_part) - lj_part) + this->EPSILON_WCA_;
    }
    
    return energy;
}

// ============================
/* Debye-Huckel interaction potential */
// ============================
template<typename number>
number InteractionFactory<PatchyRod<number>, number>::DebyeHuckel_(number r)
{
    number energy(0.);
    
    if ( r < this->R_CUT_ ) energy = exp(r * this->MINUS_KAPPA_) * (this->DH_PREFACTOR_ / r);
    
    return energy;
}

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<PatchyRod<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                       PatchyRod<number>* Particle1, PatchyRod<number>* Particle2)
{
    number energy(0.);
    
    // Backbone-backbone WCA repulsion
    Matrix3X<number> Backbone1 = Particle1->Orientation * Particle1->Backbone;
    Matrix3X<number> Backbone2 = Particle2->Orientation * Particle2->Backbone;
    
    for ( uint idx_vtx1 = 0; idx_vtx1 < Backbone1.cols() && energy < this->E_CUT_; ++idx_vtx1 )
    {
        for ( uint idx_vtx2 = 0; idx_vtx2 < Backbone2.cols() && energy < this->E_CUT_; ++idx_vtx2 )
        {
            Vector3<number> R_sep = R_cm + Backbone2.col(idx_vtx2) - Backbone1.col(idx_vtx1);
            number          norm  = R_sep.norm();
            
            energy               += RepulsiveWCA_(norm);
        }
    }
    
    if ( energy < this->E_CUT_ )
    {
        Matrix3X<number> Patches1 = Particle1->Orientation * Particle1->Patches;
        Matrix3X<number> Patches2 = Particle2->Orientation * Particle2->Patches;
        
        if ( USE_DH )
        {
            // Assume pairwise-additivity for site-site electrostatic interactions
            for ( uint idx_vtx1 = 0; idx_vtx1 < Patches1.cols() && energy < this->E_CUT_; ++idx_vtx1 )
            {
                for ( uint idx_vtx2 = 0; idx_vtx2 < Patches2.cols() && energy < this->E_CUT_; ++idx_vtx2 )
                {
                    Vector3<number> R_sep = R_cm + Patches2.col(idx_vtx2) - Patches1.col(idx_vtx1);
                    number          norm  = R_sep.norm();
                    
                    energy               += DebyeHuckel_(norm);
                }
            }
        }
        
        else
        {
            // Hard spherical thread of radius SIGMA_R
            for ( uint idx_vtx1 = 0; idx_vtx1 < Patches1.cols() && energy < this->E_CUT_; ++idx_vtx1 )
            {
                for ( uint idx_vtx2 = 0; idx_vtx2 < Patches2.cols() && energy < this->E_CUT_; ++idx_vtx2 )
                {
                    Vector3<number> R_sep = R_cm + Patches2.col(idx_vtx2) - Patches1.col(idx_vtx1);
                    number          norm  = R_sep.norm();
                    
                    if ( norm < this->SIGMA_R ) return 1.;
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
template<typename number>
number InteractionFactory<ThreadedRod<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                         ThreadedRod<number>* Particle1, ThreadedRod<number>* Particle2)
{
    if ( !USE_DH ) return 1.;
    else
    {
        BNode<number> Node1 = *Particle1->Hull;
        BNode<number> Node2 = *Particle2->Hull;

        Node1.l_ch = this->L_Z_    / 2.;
        Node1.l_cr = this->D_HARD_ / 2.;
        
        Node2.l_ch = this->L_Z_    / 2.;
        Node2.l_cr = this->D_HARD_ / 2.;
        
        // Hard SC backbone
        if ( Utils<number>::OverlapBoundSC(R_cm, &Node1, &Node2) ) return 1.;
        
        number energy(0.);
        this->RecursiveInteraction(R_cm, Particle1->Hull, Particle2->Hull, &energy, this->E_CUT_);
        
        return 1. - exp(-energy);
    }
}


// ============================


/* ************************* */
/* TwistedCuboid */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<TwistedCuboid<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                           TwistedCuboid<number>* Particle1, TwistedCuboid<number>* Particle2)
{
    number energy(0.);
    number mayer_interaction(0.);
    
    // In this case, energy counts the number of overlapping vertices
    this->RecursiveInteraction(R_cm, Particle1->Hull, Particle2->Hull, &energy, 0.);

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
template<typename number>
number InteractionFactory<TwistedPentagon<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                             TwistedPentagon<number>* Particle1, TwistedPentagon<number>* Particle2)
{
    number energy(0.);
    number mayer_interaction(0.);
    
    // In this case, energy counts the number of overlapping vertices
    this->RecursiveInteraction(R_cm, Particle1->Hull, Particle2->Hull, &energy, 0.);
    
    if ( energy > 0. ) mayer_interaction = 1.;
    
    return mayer_interaction;
}


template struct InteractionFactory<MESOGEN<float>, float>;
template struct InteractionFactory<MESOGEN<double>, double>;
