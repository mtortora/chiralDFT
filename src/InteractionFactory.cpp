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
            number lj_part = pow(sigma/r, 6);
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
/* FlexibleChain */
/* ************************* */

// ============================
/* WCA with abrupt cutoff */
// ============================
template<typename number>
number InteractionFactory<FlexibleChain<number>, number>::RepulsiveWCA_(number r)
{
    number energy(0.);
        
    if ( r < this->RC_WCA_ )
    {
        number lj_part = SQR(CUB(this->SIGMA_R/r));
        energy         = 4. * this->EPSILON_ * (SQR(lj_part) - lj_part) + this->EPSILON_;
    }
        
    return energy;
}

// ============================
/* Debye-Huckel interaction potential */
// ============================
template<typename number>
number InteractionFactory<FlexibleChain<number>, number>::DebyeHuckel_(number r)
{
    number energy(0.);

    if ( r > this->SIGMA_R )
    {
        number h = r - this->SIGMA_R;
        
        if ( h < this->RC_DH_ )
        {
            number y2h = exp(-this->MINUS_KAPPA_*h) * SQR(atanh(exp(this->MINUS_KAPPA_*h/2.)*this->TYS_));
            energy     = this->DH_PREFACTOR_/r * y2h * log(1.+exp(this->MINUS_KAPPA_*h));
            //energy = this->DH_PREFACTOR_/r * (2.*atanh(exp(this->MINUS_KAPPA_*h))+log(1.-exp(2.*this->MINUS_KAPPA_*h)));
        }
    }
    
    return energy;
}

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<FlexibleChain<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                           FlexibleChain<number>* Particle1, FlexibleChain<number>* Particle2)
{
    number energy(0.);
    
    // Work out pairwise interaction energies recursively
    this->RecursiveInteraction(R_cm, Particle1->Hull, Particle2->Hull, &energy, this->E_CUT_);
    
    return 1. - exp(-energy);
}


// ============================


/* ************************* */
/* FlexibleHelix */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<FlexibleHelix<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                           FlexibleHelix<number>* Particle1, FlexibleHelix<number>* Particle2)
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
/* FlexiblePatchyRod */
/* ************************* */

// ============================
/* LJ interaction potential */
// ============================
template<typename number>
number InteractionFactory<FlexiblePatchyRod<number>, number>::RepulsiveLJ_(number r, number d, number range)
{
    number energy(0.);
    
    if ( r < range )
    {
        number lj_part = pow(d/r, 6);
        energy         = 4. * this->EPSILON_WCA_ * (SQR(lj_part) - lj_part) + this->EPSILON_WCA_;
    }
    
    return energy;
}

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<FlexiblePatchyRod<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                               FlexiblePatchyRod<number>* Particle1, FlexiblePatchyRod<number>* Particle2)
{
    number energy(0.);
    
    uint idx_conf1 = Particle1->idx_conf;
    uint idx_conf2 = Particle2->idx_conf;

    Matrix3X<number> Backbone1 = Matrix3X<number>::Map(Particle1->Backbones.data() + 3 * idx_conf1*Particle1->N_BCK, 3, Particle1->N_BCK);
    Matrix3X<number> Backbone2 = Matrix3X<number>::Map(Particle2->Backbones.data() + 3 * idx_conf2*Particle2->N_BCK, 3, Particle2->N_BCK);

    Backbone1 = Particle1->Orientation * Backbone1;
    Backbone2 = Particle2->Orientation * Backbone2;

    // Backbone-backbone WCA repulsion
    for ( uint idx_vtx1 = 0; idx_vtx1 < Particle1->N_BCK && energy < this->E_CUT_; ++idx_vtx1 )
    {
        for ( uint idx_vtx2 = 0; idx_vtx2 < Particle2->N_BCK && energy < this->E_CUT_; ++idx_vtx2 )
        {
            Vector3<number> R_sep = R_cm + Backbone2.col(idx_vtx2) - Backbone1.col(idx_vtx1);
            number          norm  = R_sep.norm();
            
            energy               += RepulsiveLJ_(norm, this->D_BACK_, this->R_BACK_);
        }
    }
    
    if ( energy < this->E_CUT_ )
    {
        Matrix3X<number> Patch1 = Matrix3X<number>::Map(Particle1->Patches.data() + 3 * idx_conf1*Particle1->N_BCK, 3, Particle1->N_BCK);
        Matrix3X<number> Patch2 = Matrix3X<number>::Map(Particle2->Patches.data() + 3 * idx_conf2*Particle2->N_BCK, 3, Particle2->N_BCK);
        
        Patch1 = Particle1->Orientation * Patch1;
        Patch2 = Particle2->Orientation * Patch2;
        
        // Lorentz-Berthelot patch-backbone repulsion
        for ( uint idx_vtx1 = 0; idx_vtx1 < Particle1->N_BCK && energy < this->E_CUT_; ++idx_vtx1 )
        {
            for ( uint idx_vtx2 = 0; idx_vtx2 < Particle2->N_BCK && energy < this->E_CUT_; ++idx_vtx2 )
            {
                Vector3<number> R_sep = R_cm + Backbone2.col(idx_vtx2) - Patch1.col(idx_vtx1);
                number          norm  = R_sep.norm();
                
                energy               += RepulsiveLJ_(norm, this->D_LB_, this->R_LB_);
            }
        }
        
        for ( uint idx_vtx1 = 0; idx_vtx1 < Particle1->N_BCK && energy < this->E_CUT_; ++idx_vtx1 )
        {
            for ( uint idx_vtx2 = 0; idx_vtx2 < Particle2->N_BCK && energy < this->E_CUT_; ++idx_vtx2 )
            {
                Vector3<number> R_sep = R_cm + Patch2.col(idx_vtx2) - Backbone1.col(idx_vtx1);
                number          norm  = R_sep.norm();
                
                energy               += RepulsiveLJ_(norm, this->D_LB_, this->R_LB_);
            }
        }
        
        // Patch-patch LJ repulsion
        for ( uint idx_vtx1 = 0; idx_vtx1 < Particle1->N_BCK && energy < this->E_CUT_; ++idx_vtx1 )
        {
            for ( uint idx_vtx2 = 0; idx_vtx2 < Particle2->N_BCK && energy < this->E_CUT_; ++idx_vtx2 )
            {
                Vector3<number> R_sep = R_cm + Patch2.col(idx_vtx2) - Patch1.col(idx_vtx1);
                number          norm  = R_sep.norm();
                
                energy               += RepulsiveLJ_(norm, this->D_PATCH_, this->R_PATCH_);
            }
        }
    }
    
    return 1. - exp(-energy);
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
        number lj_part = pow(this->SIGMA_R/r, 6);
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
    return (( r < this->R_CUT_ ) ? exp(r * this->MINUS_KAPPA_) * (this->DH_PREFACTOR_ / r) : 0.);
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
/* Debye-Huckel interaction potential */
// ============================
template<typename number>
number InteractionFactory<ThreadedRod<number>, number>::DebyeHuckel_(number r)
{
    return (( r < this->R_CUT_ ) ? exp(r * this->MINUS_KAPPA_) * (this->DH_PREFACTOR_ / r) : 0.);
}

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<ThreadedRod<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                         ThreadedRod<number>* Particle1, ThreadedRod<number>* Particle2)
{
    // Hard backbone checks are handled in MCIntegrator.cpp in the absence of DH interactions
    if ( !USE_DH ) return 1.;
    else
    {
        number energy(0.);

        Matrix3X<number> Patches1 = Particle1->Orientation * Particle1->Patches;
        Matrix3X<number> Patches2 = Particle2->Orientation * Particle2->Patches;
        
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
/* TwistedHexagon */
/* ************************* */

// ============================
/* Mayer interaction function */
// ============================
template<typename number>
number InteractionFactory<TwistedHexagon<number>, number>::MayerInteraction(const Vector3<number>& R_cm,
                                                                            TwistedHexagon<number>* Particle1, TwistedHexagon<number>* Particle2)
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
