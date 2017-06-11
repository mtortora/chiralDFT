#ifndef ODF_OPTIMISER_HPP_
#define ODF_OPTIMISER_HPP_

#include "MCIntegrator.hpp"


template<template<typename number> class ParticleType, typename number>
class OdfOptimiser: public MCIntegrator<ParticleType, number>
{
public:
    OdfOptimiser();
    
    ArrayX<number> Alpha_grid;
    ArrayX<number> Theta_grid;
    ArrayX<number> Phi_grid;

    void BinodalAnalysis(const ArrayX<number>&);
    void EnergyGrid(const ArrayX<number>&, ArrayX<number>*);
    void ODFGrid(const ArrayX<number>&,
                 ArrayX<number>*, ArrayX<number>*, ArrayX<number>*, ArrayX<number>*,
                 ArrayXX<number>*);
    
private:
    ArrayX<number> Psi_iso_;
    
    number RotationalEnt(const ArrayX<number>&);
    number OrderParam   (const ArrayX<number>&);
    number VirialCoeff  (const ArrayX<number>&, const ArrayX<number>&);
    
    number FreeEnergy(number, const ArrayX<number>&, const ArrayX<number>&);
    
    ArrayX<number> LegendreCoeffs(const ArrayX<number>&);
    ArrayX<number> SequentialOptimiser(number, const ArrayX<number>&);

    Vector2<number> ODFThermo(number, const ArrayX<number>&, const ArrayX<number>&);
};

#endif
