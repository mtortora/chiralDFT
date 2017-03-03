#ifndef INTERACTION_FACTORY_HPP_
#define INTERACTION_FACTORY_HPP_

#include "BaseInteraction.hpp"
#include "Particles/BentCore.hpp"
#include "Particles/DNADuplex.hpp"
#include "Particles/Helix.hpp"
#include "Particles/PatchyRod.hpp"
#include "Particles/ThreadedRod.hpp"
#include "Particles/TwistedCuboid.hpp"
#include "Particles/TwistedPentagon.hpp"


// Full template specialisation
template<class ParticleType>
struct InteractionFactory {};


// BentCore
template<>
struct InteractionFactory<BentCore>: public BaseInteraction<InteractionFactory<BentCore> >, BentCore
{
    double MayerInteraction(const Eigen::Vector3d&, BentCore*, BentCore*);
    
    inline double PairInteraction(double r) {return (r < R_HARD_);}
};


// DNADuplex
template<>
struct InteractionFactory<DNADuplex>: public BaseInteraction<InteractionFactory<DNADuplex> >, DNADuplex
{
    double MayerInteraction(const Eigen::Vector3d&, DNADuplex*, DNADuplex*);
    
    inline double PairInteraction(double r)
    {
        if ( USE_DH )
        {
            return DebyeHuckel_(r) + RepulsiveLJ_(r, EXCL_S1_, EXCL_R1_, EXCL_B1_, EXCL_RC1_);
        }
        
        return RepulsiveLJ_(r, EXCL_S1_, EXCL_R1_, EXCL_B1_, EXCL_RC1_);
    }
    
private:
    double RepulsiveLJ_(double, double , double, double, double);
    double DebyeHuckel_(double);
};


// Helix
template<>
struct InteractionFactory<Helix>: public BaseInteraction<InteractionFactory<Helix> >, Helix
{
    double MayerInteraction(const Eigen::Vector3d&, Helix*, Helix*);
    
    inline double PairInteraction(double r) {return (r < R_HARD_);}
};


// PatchyRod
template<>
struct InteractionFactory<PatchyRod>: public PatchyRod
{
    double MayerInteraction(const Eigen::Vector3d&, PatchyRod*, PatchyRod*);
    
private:
    double RepulsiveWCA_(double);
    double DebyeHuckel_(double);
};


// ThreadedRod
template<>
struct InteractionFactory<ThreadedRod>: public BaseInteraction<InteractionFactory<ThreadedRod> >, ThreadedRod
{
    double MayerInteraction(const Eigen::Vector3d&, ThreadedRod*, ThreadedRod*);
    
    inline double PairInteraction(double r)
    {
        // Debye-Huckel with abrupt cutoff at r = R_CUT_
        return (( r < R_CUT_ ) ? exp(r * MINUS_KAPPA_) * (DH_PREFACTOR_ / r) : 0.);
    }
};


// TwistedCuboid
template<>
struct InteractionFactory<TwistedCuboid>: public BaseInteraction<InteractionFactory<TwistedCuboid> >, TwistedCuboid
{
    double MayerInteraction(const Eigen::Vector3d&, TwistedCuboid*, TwistedCuboid*);
    
    inline double PairInteraction(double r) {return (r < R_THRESHOLD_);}
};


// TwistedPentagon
template<>
struct InteractionFactory<TwistedPentagon>: public BaseInteraction<InteractionFactory<TwistedPentagon> >, TwistedPentagon
{
    double MayerInteraction(const Eigen::Vector3d&, TwistedPentagon*, TwistedPentagon*);
    
    inline double PairInteraction(double r) {return (r < R_THRESHOLD_);}
};

#endif
