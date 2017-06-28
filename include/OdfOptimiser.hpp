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
                 ArrayX<number>*, ArrayX<number>*, ArrayX<number>*,
                 ArrayXX<std::complex<number> >*, ArrayXX<number>*, int, int);
    
private:
    ArrayX<double> Psi_iso_;
    
    double RotationalEnt(const ArrayX<double>&);
    double VirialCoeff  (const ArrayX<double>&, const ArrayXX<number>&);
    
    double FreeEnergy(number, const ArrayX<double>&, const ArrayXX<number>&);

    ArrayX<double> LegendreCoeffs(const ArrayX<double>&);
    ArrayX<double> SequentialOptimiser(number, const ArrayXX<number>&, int, int);

    void OrderParams(const ArrayX<double>&, Matrix33<double>*, ArrayX<std::complex<double> >*);
    Vector2<double> ODFThermo(number, const ArrayX<double>&, const ArrayXX<number>&);
};

#endif
