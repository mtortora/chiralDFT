#ifndef ODF_OPTIMISER_HPP_
#define ODF_OPTIMISER_HPP_

#include "MCIntegrator.hpp"


template<template<typename number> class ParticleType, typename number>
class OdfOptimiser: public MCIntegrator<ParticleType, number>
{
public:
    OdfOptimiser();
    
    ArrayX<number> Theta_grid;
    
    void BinodalAnalysis(const MatrixXX<number>&);
    void EnergyGrid(const MatrixXX<number>&, ArrayX<number>*);
    void ODFGrid(const MatrixXX<number>&,
                 ArrayX<number>*, ArrayX<number>*, ArrayX<number>*, ArrayX<number>*,
                 ArrayXX<number>*);
    
private:
    ArrayX<number> Psi_iso_;
    
    number VirialCoeff(const ArrayX<number>&, const MatrixXX<number>&);
    number FreeEnergy(number, const ArrayX<number>&, const MatrixXX<number>&);
    
    ArrayX<number> LegendreCoeffs(const ArrayX<number>&);
    ArrayX<number> SequentialOptimiser(number, const MatrixXX<number>&);
    
    Vector2<number> ODFThermo(number, const ArrayX<number>&, const MatrixXX<number>&);
};

#endif
