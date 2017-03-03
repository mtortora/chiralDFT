// ===================================================================
/**
 * Monte-Carlo integrator templated base class.
 * All particle headers must be included in InteractionFactory.hpp
 */
// ===================================================================
/*
 * MCIntegrator.cpp: Version 2.5
 * Created 16/09/2015 by Maxime Tortora
 */
// ===================================================================

#include <mpi.h>

#include "Legendre.hpp"
#include "MCIntegrator.hpp"

using namespace Eigen;


// ============================
/* Class constructor */
// ============================
template<class ParticleType>
MCIntegrator<ParticleType>::MCIntegrator()
{
    Q_grid   = ArrayXd::LinSpaced(N_STEPS_Q, Q_MIN, Q_MAX);
    Eta_grid = ArrayXd::LinSpaced(N_STEPS_ETA, ETA_MIN, ETA_MAX);
}

// ============================
/* Initialise MPI parameters and seed RNGs */
// ============================
template<class ParticleType>
void MCIntegrator<ParticleType>::MCInit(int seed, int mpi_rank, int mpi_size)
{
    // Build particle models
    Particle1_.Build(mpi_rank);
    Particle2_.Build(mpi_rank);
    
    // Number of MC steps to be performed by each thread
    N_PER_PROC_      = N_MC / mpi_size;
    
    // Assign extra constants to interaction manager
    IManager.R_INTEG = Particle1_.R_INTEG;
    IManager.V_INTEG = Particle1_.V_INTEG;
    IManager.V0      = Particle1_.V0;
    IManager.V_EFF   = Particle1_.V_EFF;

    // Generate independant seeds for each thread
    std::mt19937_64::result_type engine_seed = seed + mpi_rank;
    
    // Seed 64-bit Mersenne twister RNG engine
    rng_engine_.seed(engine_seed);
}

// ============================
/* Reset MC counters */
// ============================
template<class ParticleType>
void MCIntegrator<ParticleType>::MCReset()
{
    MPI_Barrier(MPI_COMM_WORLD);
    
    t_start_ = MPI_Wtime();
    
    ctr_mc_  = 0;
    ctr_ov_  = 0;
}

// ============================
/* Thread sync check every (1/N_SYNC)% */
// ============================
template<class ParticleType>
void MCIntegrator<ParticleType>::SyncCheck(bool prune)
{
    if ( ctr_mc_ % (N_PER_PROC_ / N_SYNC) == 0 )
    {
        ullint ctr_mc;
        
        // Aggregate simulation statistics
        MPI_Reduce(&ctr_mc_, &ctr_mc, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
        
        t_end_     = MPI_Wtime();
        t_elapsed_ = t_end_ - t_start_;
        t_start_   = t_end_;
        
        LogTxt("-------");
        LogCya("Completed %llu out of %llu MC steps", ctr_mc, (ullint)N_MC);
        LogInf("Time: %f ws (%.1f configurations/sec)", t_elapsed_, N_MC / (N_SYNC*t_elapsed_));
        
        if ( prune )
        {
            MPI_Allreduce(MPI_IN_PLACE, Exc_grid_.data(), Exc_grid_.size(), MPI_INT, MPI_BAND, MPI_COMM_WORLD);
            
            ullint ctr_ov = NX*NY*NZ - Exc_grid_.sum();
            
            double frac_p = (ctr_ov - ctr_ov_) * 100. / ((double)NX*NY*NZ);
            double frac_t =  ctr_ov            * 100. / ((double)NX*NY*NZ);
            
            LogInf("Pruned %.3f%% of the grid (total: %.3f%%)", frac_p, frac_t);
            
            ctr_ov_ = ctr_ov;
        }
        
        else
        {
            ullint ctr_ov;
            
            MPI_Reduce(&ctr_ov_, &ctr_ov, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
            
            LogInf("%.3f%% configurations interacted", ctr_ov * 100.*N_SYNC/N_MC);
            
            ctr_ov_ = 0;
        }
    }
}

#if (MODE == MODE_EXC)

// ============================
/* Prune grid points outside of excluded volume manifold */
// ============================
template<class ParticleType>
void MCIntegrator<ParticleType>::PruneGrid(double r_hard)
{
    Vector3d Grid_point;
    
    MatrixXd Backbone   = Particle2_.Orientation * Particle2_.Backbone;
    Backbone.colwise() += R_cm_;
    
    // Iterate over grid points
    for ( uint idx_x = 0; idx_x < NX; ++idx_x )
    {
        for ( uint idx_y = 0; idx_y < NY; ++idx_y )
        {
            for ( uint idx_z = 0; idx_z < NZ; ++idx_z )
            {
                uint idx_point = idx_x*NY*NZ + idx_y*NZ + idx_z;
                
                if ( Exc_grid_(idx_point) )
                {
                    Grid_point << X_grid_(idx_x), Y_grid_(idx_y), Z_grid_(idx_z);
                    
                    // Prune grid points located in non-overlapping particles pairs
                    for ( uint idx_vtx = 0; idx_vtx < Backbone.cols(); ++idx_vtx )
                    {
                        Vector3d R_sep = Backbone.col(idx_vtx) - Grid_point;
                        
                        if ( R_sep.norm() < r_hard )
                        {
                            Exc_grid_(idx_point) = 0;
                            break;
                        }
                    }
                }
            }
        }
    }
}

// ============================
/* MC integration of the effective particle volume */
// ============================
template<class ParticleType>
double MCIntegrator<ParticleType>::ExcludedIntegrator(double r_hard)
{
    MCReset();
    
    LogTxt("------------");
    LogTxt("Integrating effective excluded volume...");
    
    // Box spatial grid
    X_grid_   = ArrayXd::LinSpaced(NX, -Particle1_.BHull->l_xh, Particle1_.BHull->l_xh);
    Y_grid_   = ArrayXd::LinSpaced(NY, -Particle1_.BHull->l_yh, Particle1_.BHull->l_yh);
    Z_grid_   = ArrayXd::LinSpaced(NZ, -Particle1_.BHull->l_zh, Particle1_.BHull->l_zh);
    
    Exc_grid_ = ArrayXi::Constant(NX*NY*NZ, 1);

    while ( ctr_mc_ < N_PER_PROC_ )
    {
        ctr_mc_++;
        
        SyncCheck(true);
        
        // Fetch random configurations if relevant
        Particle1_.Parse(rng_engine_);
        Particle2_.Parse(rng_engine_);
        
        // Random center-of-mass to center-of-mass separation vector
        R_cm_ = Vector3d::NullaryExpr([&](int) {return rng_distrib_(rng_engine_)-0.5;}) * 2.*IManager.R_INTEG;
        
        // Bounding sphere overlap test
        if ( R_cm_.norm() < IManager.R_INTEG )
        {
            Particle2_.SetRandomAxis(rng_engine_, rng_distrib_);
            
            // Bounding spherocylinder (SC) overlap test
            if ( Utils::OverlapBoundSC(R_cm_, Particle1_.BHull, Particle2_.BHull) )
            {
                Particle2_.SetRandomOrientation(rng_engine_, rng_distrib_);
                
                // Oriented Bounding Box (OBB) overlap test
                if ( Utils::OverlapBoundOB(R_cm_, Particle1_.BHull, Particle2_.BHull) )
                {
                    mayer_interaction_ = IManager.MayerInteraction(R_cm_, &Particle1_, &Particle2_);
                    
                    if ( mayer_interaction_ == 0. ) PruneGrid(r_hard);
                }
            }
        }
    }
    
    double f_exc = Exc_grid_.sum() / ((double)NX*NY*NZ);
    double v_box = 8. * Particle1_.BHull->l_xh*Particle1_.BHull->l_yh*Particle1_.BHull->l_zh;
    
    return (v_box * f_exc);
}

#endif

// ============================
/* Random configuration generator */
// ============================
template<class ParticleType>
void MCIntegrator<ParticleType>::ConfigGenerator()
{
    mayer_interaction_ = 0.;
    
    while ( (mayer_interaction_ == 0.) && (ctr_mc_ < N_PER_PROC_) )
    {
        ctr_mc_++;
        
        SyncCheck(false);
        
        // Fetch random configurations if relevant
        Particle1_.Parse(rng_engine_);
        Particle2_.Parse(rng_engine_);
        
        // Random center-of-mass to center-of-mass separation vector
        R_cm_ = Vector3d::NullaryExpr([&](int) {return rng_distrib_(rng_engine_)-0.5;}) * 2.*IManager.R_INTEG;
        
        // Bounding sphere overlap test
        if ( R_cm_.norm() < IManager.R_INTEG )
        {
            Particle1_.SetRandomAxis(rng_engine_, rng_distrib_);
            Particle2_.SetRandomAxis(rng_engine_, rng_distrib_);
            
            // Bounding spherocylinder (SC) overlap test
            if ( Utils::OverlapBoundSC(R_cm_, Particle1_.BHull, Particle2_.BHull) )
            {
                Particle1_.SetRandomOrientation(rng_engine_, rng_distrib_);
                Particle2_.SetRandomOrientation(rng_engine_, rng_distrib_);
                
                #if USE_RAPID
                    // Use RAPID library for collision detection if implemented
                    double pos1[3];
                    double pos2[3];
                    double rot1[3][3];
                    double rot2[3][3];
                    
                    // Convert Eigen containers to C arrays
                    for ( uint idx1 = 0; idx1 < 3; ++idx1 )
                    {
                        pos1[idx1] = 0.;
                        pos2[idx1] = R_cm_(idx1);
                        
                        for ( uint idx2 = 0; idx2 < 3; ++idx2 )
                        {
                            rot1[idx1][idx2] = Particle1_.Orientation(idx1, idx2);
                            rot2[idx1][idx2] = Particle2_.Orientation(idx1, idx2);
                        }
                    }
                    
                    // RAPID collision query
                    RAPID_Collide(rot1, pos1, Particle1_.Mesh, rot2, pos2, Particle2_.Mesh, RAPID_FIRST_CONTACT);
                    
                    if ( RAPID_num_contacts > 0 ) mayer_interaction_ = 1.;
                
                #else
                    // Oriented Bounding Box (OBB) overlap test
                    if ( Utils::OverlapBoundOB(R_cm_, Particle1_.BHull, Particle2_.BHull) )
                    {
                        // Bounding Volume Hierarchy-accelerated energy computations
                        mayer_interaction_ = IManager.MayerInteraction(R_cm_, &Particle1_, &Particle2_);
                    }
                
                #endif
            }
        }
    }
}

// ============================
/* MC integration of the full angle-dependant excluded volume */
// ============================
template<class ParticleType>
void MCIntegrator<ParticleType>::FullIntegrator(MatrixXd* E_out)
{
    MCReset();

    LogTxt("------------");
    LogTxt("Integrating full Onsager angle-dependant second-virial coefficient");
    
    ArrayXd  Theta_grid = ArrayXd::LinSpaced(N_STEPS_THETA+1, 0, PI);
    MatrixXd E_loc      = MatrixXd::Zero(N_STEPS_THETA, N_STEPS_THETA);
    
    while ( ctr_mc_ < N_PER_PROC_ )
    {
        ConfigGenerator();
        
        if ( mayer_interaction_ != 0. )
        {
            ctr_ov_++;
            
            uint   idx_theta1(0);
            uint   idx_theta2(0);
            
            double theta1 = Particle1_.GetTheta();
            double theta2 = Particle2_.GetTheta();
            
            for ( uint idx_ang = 0; idx_ang < N_STEPS_THETA; ++idx_ang )
            {
                if ( (theta1 > Theta_grid(idx_ang)) && (theta1 <= Theta_grid(idx_ang+1)) ) idx_theta1 = idx_ang;
                if ( (theta2 > Theta_grid(idx_ang)) && (theta2 <= Theta_grid(idx_ang+1)) ) idx_theta2 = idx_ang;
            }
            
            E_loc(idx_theta1, idx_theta2) += mayer_interaction_;
        }
    }
    
    *E_out = E_loc * IManager.V_INTEG/N_PER_PROC_ / SQR(D_THETA);
}

// ============================
/* MC integration of the Legendre coefficients */
// ============================
template<class ParticleType>
void MCIntegrator<ParticleType>::LegendreIntegrator(double q_macro, MatrixXd* E_out)
{
    MCReset();

    LogTxt("------------");
    LogTxt("Integrating Legendre-projected second-virial coefficient - q=%.5f", q_macro);
    
    Vector3d N_q;
    MatrixXd E_loc = MatrixXd::Zero(N_L, N_L);
    
    while ( ctr_mc_ < N_PER_PROC_ )
    {
        ConfigGenerator();
        
        if ( mayer_interaction_ != 0. )
        {
            ctr_ov_++;
            
            double theta1 = Particle1_.GetTheta();
            double theta2 = Particle2_.GetTheta();
            
            N_q << sin(q_macro*R_cm_(1)), 0., cos(q_macro*R_cm_(1));
            
            for ( uint l1 = 0; l1 < N_L; l1 += IManager.N_DELTA_L )
            {
                for ( uint l2 = 0; l2 < N_L; l2 += IManager.N_DELTA_L )
                {
                    E_loc(l1,l2) += sin(theta1)*sin(theta2) * mayer_interaction_
                                  * sqrt((2.*(double)l1 + 1.)*(2.*(double)l2 + 1.)/4.)
                                  * Legendre::Pn(l1, Particle2_.Axis.dot(N_q))
                                  * Legendre::Pn(l2, cos(theta1));
                }
            }
        }
    }
    
    *E_out = E_loc * IManager.V_INTEG/N_PER_PROC_;
}

// ============================
/* Work out preliminary chiral properties */
// ============================
template<class ParticleType>
void MCIntegrator<ParticleType>::FrankIntegrator(const MatrixXd& Psi_in, ArrayXd* Kt_out,
                                                 ArrayXd* K1_out, ArrayXd* K2_out, ArrayXd* K3_out,
                                                 ArrayXd* V_out, ArrayXd* F_out)
{
    MCReset();

    LogTxt("------------");
    LogTxt("Starting preliminary chiral run...");
    
    MatrixXd Psi_dot(N_STEPS_THETA, N_STEPS_ETA);

    ArrayXd  Eff_grid   = Eta_grid * this->IManager.V_EFF/this->IManager.V0;
    
    ArrayXd  G_grid     = (1. - 3/4.*Eff_grid) / SQR(1. - Eff_grid);
    ArrayXd  N_grid     = Eta_grid * CUB(IManager.SIGMA_R)/IManager.V0;
    
    ArrayXd  Theta_grid = ArrayXd::LinSpaced(N_STEPS_THETA+1, 0, PI);
    
    ArrayXd  Kt         = ArrayXd::Zero(N_STEPS_ETA);
    ArrayXd  K1         = ArrayXd::Zero(N_STEPS_ETA);
    ArrayXd  K2         = ArrayXd::Zero(N_STEPS_ETA);
    ArrayXd  K3         = ArrayXd::Zero(N_STEPS_ETA);

    ArrayXd  F_loc      = ArrayXd::Zero(N_STEPS_ETA);
    
    ArrayXd  V_r        = ArrayXd::Zero(N_STEPS_THETA);
    ArrayXd  V_l        = ArrayXd::Zero(N_STEPS_THETA);
    
    // Work out Psi differentials
    for ( uint idx_ang = 0; idx_ang < N_STEPS_THETA; ++idx_ang )
    {
        if ( idx_ang < N_STEPS_THETA-1 ) Psi_dot.row(idx_ang) = Psi_in.row(idx_ang+1) - Psi_in.row(idx_ang);
        else                             Psi_dot.row(idx_ang) = -Psi_dot.row(0);
    }
    
    Psi_dot /= D_THETA;
    
    while ( ctr_mc_ < N_PER_PROC_ )
    {
        ConfigGenerator();
        
        if ( mayer_interaction_ != 0. )
        {
            ctr_ov_++;
            
            uint   idx_theta (0);
            uint   idx_theta1(0);
            uint   idx_theta2(0);
            
            double theta1 = Particle1_.GetTheta();
            double theta2 = Particle2_.GetTheta();
            
            double theta  = acos((Particle1_.Axis).dot(Particle2_.Axis));
            double deter  = R_cm_.dot((Particle1_.Axis).cross(Particle2_.Axis));
            
            // Work out configuration angles
            for ( uint idx_ang = 0; idx_ang < N_STEPS_THETA; ++idx_ang )
            {
                if ( (theta1 > Theta_grid(idx_ang)) && (theta1 <= Theta_grid(idx_ang+1)) ) idx_theta1 = idx_ang;
                if ( (theta2 > Theta_grid(idx_ang)) && (theta2 <= Theta_grid(idx_ang+1)) ) idx_theta2 = idx_ang;
                if ( (theta  > Theta_grid(idx_ang)) && (theta  <= Theta_grid(idx_ang+1)) ) idx_theta  = idx_ang;
            }
            
            // Update twist modulus integral & Frank elastic constants
            Kt -= mayer_interaction_ * sin(theta1)
                * Psi_in. row(idx_theta1).adjoint().array()
                * Psi_dot.row(idx_theta2).adjoint().array()
                * Particle2_.Axis(1)
                * R_cm_(0);
            
            K1 -= mayer_interaction_
                * Psi_dot.row(idx_theta1).adjoint().array()
                * Psi_dot.row(idx_theta2).adjoint().array()
                * Particle1_.Axis(0) * Particle2_.Axis(0)
                * SQR(R_cm_(0));
			
            K2 -= mayer_interaction_
                * Psi_dot.row(idx_theta1).adjoint().array()
                * Psi_dot.row(idx_theta2).adjoint().array()
                * Particle1_.Axis(1) * Particle2_.Axis(1)
                * SQR(R_cm_(0));
            
            K3 -= mayer_interaction_
                * Psi_dot.row(idx_theta1).adjoint().array()
                * Psi_dot.row(idx_theta2).adjoint().array()
                * Particle1_.Axis(0) * Particle2_.Axis(0)
                * SQR(R_cm_(2));
            
            // Update right- and left-handed interaction integrals
            if ( deter > 0. )
            {
                F_loc += mayer_interaction_ * sin(theta1)*sin(theta2)
                       * Psi_in.row(idx_theta1).adjoint().array()
                       * Psi_in.row(idx_theta2).adjoint().array();
                
                V_r(idx_theta) += mayer_interaction_;
            }
            
            else
            {
                F_loc -= mayer_interaction_ * sin(theta1)*sin(theta2)
                       * Psi_in.row(idx_theta1).adjoint().array()
                       * Psi_in.row(idx_theta2).adjoint().array();
                
                V_l(idx_theta) += mayer_interaction_;
            }
        }
    }
    
    *V_out  = (V_r - V_l) / (V_r + V_l);
    
    *F_out  = F_loc * IManager.V_INTEG/N_PER_PROC_;

    *Kt_out = Kt * SQR(N_grid)/2. * G_grid * IManager.V_INTEG/N_PER_PROC_ / pow(IManager.SIGMA_R, 4);
    
    *K1_out = K1 * SQR(N_grid)/2. * G_grid * IManager.V_INTEG/N_PER_PROC_ / pow(IManager.SIGMA_R, 5);
    *K2_out = K2 * SQR(N_grid)/2. * G_grid * IManager.V_INTEG/N_PER_PROC_ / pow(IManager.SIGMA_R, 5);
    *K3_out = K3 * SQR(N_grid)/2. * G_grid * IManager.V_INTEG/N_PER_PROC_ / pow(IManager.SIGMA_R, 5);
}

template class MCIntegrator<MESOGEN>;
