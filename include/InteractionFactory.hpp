#ifndef INTERACTION_FACTORY_HPP_
#define INTERACTION_FACTORY_HPP_

#include "BaseInteraction.hpp"
#include "Particles/BentCore.hpp"
#include "Particles/DNADuplex.hpp"
#include "Particles/FlexibleHelix.hpp"
#include "Particles/Helix.hpp"
#include "Particles/PatchyRod.hpp"
#include "Particles/ThreadedRod.hpp"
#include "Particles/TriangularPrism.hpp"
#include "Particles/TwistedCuboid.hpp"
#include "Particles/TwistedPentagon.hpp"


// Requires full template specialisation
template<class ParticleType, typename number>
struct InteractionFactory {};


// BentCore
template<typename number>
struct InteractionFactory<BentCore<number>, number>: public BaseInteraction<InteractionFactory<BentCore<number>, number>, number>, BentCore<number>
{
    number MayerInteraction(const Vector3<number>&, BentCore<number>*, BentCore<number>*);
    
    inline number PairInteraction(number r) {return (r < this->D_HARD_);}
};


// DNADuplex
template<typename number>
struct InteractionFactory<DNADuplex<number>, number>: public BaseInteraction<InteractionFactory<DNADuplex<number>, number>, number>, DNADuplex<number>
{
    number MayerInteraction(const Vector3<number>&, DNADuplex<number>*, DNADuplex<number>*);
    
    inline number PairInteraction(number r)
    {
        if ( USE_DH )
        {
            return DebyeHuckel_(r) + RepulsiveLJ_(r, this->EXCL_S1_, this->EXCL_R1_, this->EXCL_B1_, this->EXCL_RC1_);
        }
        
        return RepulsiveLJ_(r, this->EXCL_S1_, this->EXCL_R1_, this->EXCL_B1_, this->EXCL_RC1_);
    }
    
private:
    number RepulsiveLJ_(number, number , number, number, number);
    number DebyeHuckel_(number);
};


// FlexibleHelix
template<typename number>
struct InteractionFactory<FlexibleHelix<number>, number>: public BaseInteraction<InteractionFactory<FlexibleHelix<number>, number>, number>, FlexibleHelix<number>
{
    number MayerInteraction(const Vector3<number>&, FlexibleHelix<number>*, FlexibleHelix<number>*);
    
    inline number PairInteraction(number r) {return (r < this->D_HARD_);}
};


// Helix
template<typename number>
struct InteractionFactory<Helix<number>, number>: public BaseInteraction<InteractionFactory<Helix<number>, number>, number>, Helix<number>
{
    number MayerInteraction(const Vector3<number>&, Helix<number>*, Helix<number>*);
    
    inline number PairInteraction(number r) {return (r < this->D_HARD_);}
};


// PatchyRod
template<typename number>
struct InteractionFactory<PatchyRod<number>, number>: public PatchyRod<number>
{
    number MayerInteraction(const Vector3<number>&, PatchyRod<number>*, PatchyRod<number>*);
    
private:
    number RepulsiveWCA_(number);
    number DebyeHuckel_(number);
};


// ThreadedRod
template<typename number>
struct InteractionFactory<ThreadedRod<number>, number>: public BaseInteraction<InteractionFactory<ThreadedRod<number>, number>, number>, ThreadedRod<number>
{
    number MayerInteraction(const Vector3<number>&, ThreadedRod<number>*, ThreadedRod<number>*);
    
    inline number PairInteraction(number r)
    {
        // Debye-Huckel with abrupt cutoff at r = R_CUT_
        return (( r < this->R_CUT_ ) ? exp(r * this->MINUS_KAPPA_) * (this->DH_PREFACTOR_ / r) : 0.);
    }
};


// TriangularPrism
template<typename number>
struct InteractionFactory<TriangularPrism<number>, number>: public TriangularPrism<number> {};


// TwistedCuboid
template<typename number>
struct InteractionFactory<TwistedCuboid<number>, number>: public BaseInteraction<InteractionFactory<TwistedCuboid<number>, number>, number>, TwistedCuboid<number>
{
    number MayerInteraction(const Vector3<number>&, TwistedCuboid<number>*, TwistedCuboid<number>*);
    
    inline number PairInteraction(number r) {return (r < this->R_THRESHOLD_);}
};


// TwistedPentagon
template<typename number>
struct InteractionFactory<TwistedPentagon<number>, number>: public BaseInteraction<InteractionFactory<TwistedPentagon<number>, number>, number>, TwistedPentagon<number>
{
    number MayerInteraction(const Vector3<number>&, TwistedPentagon<number>*, TwistedPentagon<number>*);
    
    inline number PairInteraction(number r) {return (r < this->R_THRESHOLD_);}
};

#endif
