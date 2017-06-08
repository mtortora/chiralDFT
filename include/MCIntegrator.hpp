#ifndef MC_INTEGRATOR_HPP_
#define MC_INTEGRATOR_HPP_

#include <random>

#include "InteractionFactory.hpp"


template<template<typename number> class ParticleType, typename number>
class MCIntegrator
{
public:
    MCIntegrator();
    
    number ExcludedIntegrator(number);
    
    InteractionFactory<ParticleType<number>, number> IManager;
    
    ArrayX<number> Q_grid;
    ArrayX<number> Eta_grid;
    
    void MCInit(int, int, int);
    void FullIntegrator    (MatrixXX<number>*, ArrayX<number>*, ArrayX<number>*);
    void LegendreIntegrator(MatrixXX<number>*, ArrayX<number>*, ArrayX<number>*);
    void LegendreIntegrator(MatrixXX<number>*, number);
    void FrankIntegrator(const ArrayXX<number>&,
                         ArrayX<number>*, ArrayX<number>*, ArrayX<number>*, ArrayX<number>*);
    
private:
    ParticleType<number> Particle1_;
    ParticleType<number> Particle2_;
    
    ullint N_PER_PROC_;

    ullint ctr_mc_;
    ullint ctr_ov_;
    
    number mayer_interaction_;
    number t_start_;
    number t_end_;
    number t_elapsed_;

    Vector3<number> R_cm_;
    
    ArrayX<number>  X_grid_;
    ArrayX<number>  Y_grid_;
    ArrayX<number>  Z_grid_;

    ArrayX<number>  Theta_grid_;
    
    ArrayX<uint>    Exc_grid_;

    std::mt19937_64 rng_engine_;
    std::uniform_real_distribution<number> rng_distrib_{0., 1.};
	
    void MCReset();
    void SyncCheck(bool);
    void PruneGrid(number);
    void ConfigGenerator();
};

#endif
