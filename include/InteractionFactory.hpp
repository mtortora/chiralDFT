#ifndef INTERACTION_FACTORY_HPP_
#define INTERACTION_FACTORY_HPP_

#include "BaseInteraction.hpp"
#include "Particles/BentCore.hpp"
#include "Particles/DNADuplex.hpp"
#include "Particles/FlexibleChain.hpp"
#include "Particles/FDGromos.hpp"
#include "Particles/FlexibleHelix.hpp"
#include "Particles/FlexiblePatchyRod.hpp"
#include "Particles/Helix.hpp"
#include "Particles/PatchyRod.hpp"
#include "Particles/ThreadedRod.hpp"
#include "Particles/TriangularPrism.hpp"
#include "Particles/TwistedCuboid.hpp"
#include "Particles/TwistedHexagon.hpp"
#include "Particles/TwistedPentagon.hpp"


// Requires full template specialisation
template<class ParticleType, typename number>
struct InteractionFactory {};


// BentCore
template<typename number>
struct InteractionFactory<BentCore<number>, number>: public BaseInteraction<InteractionFactory<BentCore<number>, number>, number>, BentCore<number>
{
    number MayerInteraction(const Vector3<number>&, BentCore<number>*, BentCore<number>*);
    
    inline number PairInteraction(number r, number, number, uint, uint) {return (r < this->D_HARD_);}
};


// DNADuplex
template<typename number>
struct InteractionFactory<DNADuplex<number>, number>: public BaseInteraction<InteractionFactory<DNADuplex<number>, number>, number>, DNADuplex<number>
{
    number MayerInteraction(const Vector3<number>&, DNADuplex<number>*, DNADuplex<number>*);
    
    inline number PairInteraction(number r, number, number, uint, uint)
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


// FDGromos
template<typename number>
struct InteractionFactory<FDGromos<number>, number>: public BaseInteraction<InteractionFactory<FDGromos<number>, number>, number>, FDGromos<number>
{
    number MayerInteraction(const Vector3<number>&, FDGromos<number>*, FDGromos<number>*);
    
    inline number PairInteraction(number r, number c1, number c2, uint i1, uint i2)
    {
        number energy(0.);
        
        if ( (i1 < 12) && (i2 < 12) )
            energy += LJ126_(r, this->C6[i1][i2], this->C12[i1][i2]);
        
        if ( USE_DH )
        {
            if ( (c1 != 0.) && (c2 != 0.) )
                energy += ReactionField_(r, c1, c2);
        }
        
        return energy;
    }
    
private:
    number LJ126_(number, number, number);
    number ReactionField_(number, number, number);
};


// FlexibleChain
template<typename number>
struct InteractionFactory<FlexibleChain<number>, number>: public BaseInteraction<InteractionFactory<FlexibleChain<number>, number>, number>, FlexibleChain<number>
{
    number MayerInteraction(const Vector3<number>&, FlexibleChain<number>*, FlexibleChain<number>*);
    
    inline number PairInteraction(number r, number, number, uint, uint)
    {
        if ( USE_DH )
        {
            return DebyeHuckel_(r) + RepulsiveWCA_(r);
        }
        
        return RepulsiveWCA_(r);
    }
    
private:
    number RepulsiveWCA_(number);
    number DebyeHuckel_(number);
};


// FlexibleHelix
template<typename number>
struct InteractionFactory<FlexibleHelix<number>, number>: public BaseInteraction<InteractionFactory<FlexibleHelix<number>, number>, number>, FlexibleHelix<number>
{
    number MayerInteraction(const Vector3<number>&, FlexibleHelix<number>*, FlexibleHelix<number>*);
    
    inline number PairInteraction(number r, number, number, uint, uint) {return (r < this->D_HARD_);}
};


// Helix
template<typename number>
struct InteractionFactory<Helix<number>, number>: public BaseInteraction<InteractionFactory<Helix<number>, number>, number>, Helix<number>
{
    number MayerInteraction(const Vector3<number>&, Helix<number>*, Helix<number>*);
    
    inline number PairInteraction(number r, number, number, uint, uint) {return (r < this->D_HARD_);}
};


// FlexiblePatchyRod
template<typename number>
struct InteractionFactory<FlexiblePatchyRod<number>, number>: public FlexiblePatchyRod<number>
{
    number MayerInteraction(const Vector3<number>&, FlexiblePatchyRod<number>*, FlexiblePatchyRod<number>*);
    
private:
    number RepulsiveLJ_(number, number, number);
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
struct InteractionFactory<ThreadedRod<number>, number>: public ThreadedRod<number>
{
    number MayerInteraction(const Vector3<number>&, ThreadedRod<number>*, ThreadedRod<number>*);
    
private:
    number DebyeHuckel_(number);
};


// TriangularPrism
template<typename number>
struct InteractionFactory<TriangularPrism<number>, number>: public TriangularPrism<number> {};


// TwistedCuboid
template<typename number>
struct InteractionFactory<TwistedCuboid<number>, number>: public BaseInteraction<InteractionFactory<TwistedCuboid<number>, number>, number>, TwistedCuboid<number>
{
    number MayerInteraction(const Vector3<number>&, TwistedCuboid<number>*, TwistedCuboid<number>*);
    
    inline number PairInteraction(number r, number, number, uint, uint) {return (r < this->R_THRESHOLD_);}
};


// TwistedHexagon
template<typename number>
struct InteractionFactory<TwistedHexagon<number>, number>: public BaseInteraction<InteractionFactory<TwistedHexagon<number>, number>, number>, TwistedHexagon<number>
{
    number MayerInteraction(const Vector3<number>&, TwistedHexagon<number>*, TwistedHexagon<number>*);
    
    inline number PairInteraction(number r, number, number, uint, uint) {return (r < this->R_THRESHOLD_);}
};


// TwistedPentagon
template<typename number>
struct InteractionFactory<TwistedPentagon<number>, number>: public BaseInteraction<InteractionFactory<TwistedPentagon<number>, number>, number>, TwistedPentagon<number>
{
    number MayerInteraction(const Vector3<number>&, TwistedPentagon<number>*, TwistedPentagon<number>*);
    
    inline number PairInteraction(number r, number, number, uint, uint) {return (r < this->R_THRESHOLD_);}
};

#endif
