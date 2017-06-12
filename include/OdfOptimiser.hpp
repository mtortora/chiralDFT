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

    void BinodalAnalysis(const ArrayXX<number>&, int, int);
    void EnergyGrid(const ArrayXX<number>&, ArrayX<number>*, int, int);
    void ODFGrid(const ArrayXX<number>&,
                 ArrayX<number>*, ArrayX<number>*, ArrayX<number>*, ArrayX<number>*,
                 ArrayXX<number>*, int, int);
    
private:
    ArrayX<number> Psi_iso_;
    
    number RotationalEnt(const ArrayX<number>&);
    number OrderParam   (const ArrayX<number>&);
    number VirialCoeff  (const ArrayX<number>&, const ArrayXX<number>&);
    
    number FreeEnergy(number, const ArrayX<number>&, const ArrayXX<number>&);
    
    ArrayX<number> LegendreCoeffs(const ArrayX<number>&);
    ArrayX<number> SequentialOptimiser(number, const ArrayXX<number>&, int, int);

    Vector2<number> ODFThermo(number, const ArrayX<number>&, const ArrayXX<number>&);
};

#endif
