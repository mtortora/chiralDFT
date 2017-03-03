#ifndef MC_INTEGRATOR_HPP_
#define MC_INTEGRATOR_HPP_

#include <random>

#include "InteractionFactory.hpp"


template<class ParticleType>
class MCIntegrator
{
private:
    ParticleType    Particle1_;
    ParticleType    Particle2_;
    
    ullint          N_PER_PROC_;

    ullint          ctr_mc_;
    ullint          ctr_ov_;
    
    double          mayer_interaction_;
    double          t_start_;
    double          t_end_;
    double          t_elapsed_;

    Eigen::Vector3d R_cm_;
    
    Eigen::ArrayXd  X_grid_;
    Eigen::ArrayXd  Y_grid_;
    Eigen::ArrayXd  Z_grid_;

    Eigen::ArrayXi  Exc_grid_;

    std::mt19937_64 rng_engine_;
    std::uniform_real_distribution<double> rng_distrib_{0., 1.};
	
    void MCReset();
    void SyncCheck(bool);
    void PruneGrid(double);
    void ConfigGenerator();
    
public:
    MCIntegrator();

    double ExcludedIntegrator(double);
    
    InteractionFactory<ParticleType> IManager;

    Eigen::ArrayXd Q_grid;
    Eigen::ArrayXd Eta_grid;

    void MCInit(int, int, int);
    void FullIntegrator(Eigen::MatrixXd*);
    void LegendreIntegrator(double, Eigen::MatrixXd*);
    void FrankIntegrator(const Eigen::MatrixXd&,
                         Eigen::ArrayXd*, Eigen::ArrayXd*, Eigen::ArrayXd*, Eigen::ArrayXd*, Eigen::ArrayXd*,
                         Eigen::ArrayXd*);
};

#endif
