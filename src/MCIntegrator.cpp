// ===================================================================
/**
 * Monte-Carlo integrator templated base class.
 * All particle headers must be included in InteractionFactory.hpp.
 */
// ===================================================================
/*
 * MCIntegrator.cpp: Version 2.5
 * Created 16/09/2015 by Maxime Tortora
 */
// ===================================================================

#include "Legendre.hpp"
#include "MCIntegrator.hpp"


// ============================
/* Class constructor */
// ============================
template<template<typename number> class ParticleType, typename number>
MCIntegrator<ParticleType, number>::MCIntegrator()
{
    Q_grid   = ArrayX<number>::LinSpaced(N_STEPS_Q, Q_MIN, Q_MAX);
    Eta_grid = ArrayX<number>::LinSpaced(N_STEPS_ETA, ETA_MIN, ETA_MAX);
}

// ============================
/* Initialise MPI parameters and seed RNGs */
// ============================
template<template<typename number> class ParticleType, typename number>
void MCIntegrator<ParticleType, number>::MCInit(int seed, int mpi_rank, int mpi_size)
{
    // Generate independant seeds for each thread
    std::mt19937_64::result_type engine_seed = seed + mpi_rank;
    
    // Seed 64-bit Mersenne twister RNG engine
    rng_engine_.seed(engine_seed);
    
    // Build particle models
    Particle1_.Build(mpi_rank);
    Particle2_.Build(mpi_rank);
    
    // Build bootstrap index maps
    Bootstrap_map1_  = Particle1_.BootstrapMap(rng_engine_);
    Bootstrap_map2_  = Particle2_.BootstrapMap(rng_engine_);
    
    // Number of MC steps to be performed by each thread
    N_PER_PROC_      = N_MC / mpi_size;
    
    // Assign extra constants to interaction manager
    IManager.R_INTEG = Particle1_.R_INTEG;
    IManager.V_INTEG = Particle1_.V_INTEG;
    
    IManager.V0      = Particle1_.V0;
    IManager.V_EFF   = Particle1_.V_EFF;
}

// ============================
/* Reset MC counters */
// ============================
template<template<typename number> class ParticleType, typename number>
void MCIntegrator<ParticleType, number>::MCReset()
{
    MPI_Barrier(MPI_COMM_WORLD);
    
    t_start_ = MPI_Wtime();
    
    ctr_mc_  = 0;
    ctr_ov_  = 0;
}

// ============================
/* Thread sync check every (1/N_SYNC)% */
// ============================
template<template<typename number> class ParticleType, typename number>
void MCIntegrator<ParticleType, number>::SyncCheck(bool prune)
{
    if ( ctr_mc_ % (N_PER_PROC_ / N_SYNC) == 0 )
    {
        ullint ctr_mc;
        
        // Aggregate simulation statistics
        MPI_Reduce(&ctr_mc_, &ctr_mc, 1, Utils<ullint>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
        
        t_end_     = MPI_Wtime();
        t_elapsed_ = t_end_ - t_start_;
        t_start_   = t_end_;
        
        LogTxt("-------");
        LogCya("Completed %llu out of %llu MC steps", ctr_mc, (ullint)N_MC);
        LogInf("Time: %f ws (%.1f configurations/sec)", t_elapsed_, N_MC / (N_SYNC*t_elapsed_));
        
        if ( prune )
        {
            // Cull all points pruned by any thread
            MPI_Allreduce(MPI_IN_PLACE, Exc_grid_.data(), Exc_grid_.size(), Utils<uint>().MPI_type, MPI_BAND, MPI_COMM_WORLD);
            
            ullint ctr_ov = NX*NY*NZ - Exc_grid_.sum();
            
            number frac_p = (ctr_ov - ctr_ov_) * 100. / ((number)NX*NY*NZ);
            number frac_t =  ctr_ov            * 100. / ((number)NX*NY*NZ);
            
            LogInf("Pruned %.3f%% of the grid (total: %.3f%%)", frac_p, frac_t);
            
            ctr_ov_ = ctr_ov;
        }
        
        else
        {
            ullint ctr_ov;
            
            MPI_Reduce(&ctr_ov_, &ctr_ov, 1, Utils<ullint>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
            
            LogInf("%.3f%% configurations interacted", ctr_ov * 100.*N_SYNC/N_MC);
            
            ctr_ov_ = 0;
        }
    }
}

#if (MODE_SIM == MODE_EXC)

// ============================
/* Prune grid points outside of excluded volume manifold */
// ============================
template<template<typename number> class ParticleType, typename number>
void MCIntegrator<ParticleType, number>::PruneGrid(number r_hard)
{
    Vector3<number> Grid_point;
    
    Matrix3X<number> Backbone = Particle2_.Orientation * Particle2_.Backbone;
    Backbone.colwise()       += R_cm_;
    
    // Traverse grid points
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
                    
                    // Prune points located within non-overlapping bead pairs
                    for ( uint idx_vtx = 0; idx_vtx < Backbone.cols(); ++idx_vtx )
                    {
                        Vector3<number> R_sep = Backbone.col(idx_vtx) - Grid_point;
                        
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
template<template<typename number> class ParticleType, typename number>
number MCIntegrator<ParticleType, number>::ExcludedIntegrator(number r_hard)
{
    // Box spatial grid
    X_grid_   = ArrayX<number>::LinSpaced(NX, -Particle1_.Hull->l_xh, Particle1_.Hull->l_xh);
    Y_grid_   = ArrayX<number>::LinSpaced(NY, -Particle1_.Hull->l_yh, Particle1_.Hull->l_yh);
    Z_grid_   = ArrayX<number>::LinSpaced(NZ, -Particle1_.Hull->l_zh, Particle1_.Hull->l_zh);
    
    Exc_grid_ = ArrayXi::Constant(NX*NY*NZ, 1);

    LogTxt("------------");
    LogTxt("Integrating effective excluded volume...");
    
    MCReset();

    while ( ctr_mc_ < N_PER_PROC_ )
    {
        ctr_mc_++;
        
        SyncCheck(true);
        
        // Fetch random configurations if relevant
        Particle1_.Parse(rng_engine_);
        Particle2_.Parse(rng_engine_);
        
        // Random center-of-mass to center-of-mass separation vector
        R_cm_ = Vector3<number>::NullaryExpr([&](int) {return rng_distrib_(rng_engine_)-0.5;}) * 2.*IManager.R_INTEG;
        
        // Bounding sphere overlap test
        if ( R_cm_.norm() < IManager.R_INTEG )
        {
            Particle2_.SetRandomAxis(rng_engine_, rng_distrib_);
            
            // Bounding spherocylinder (SC) overlap test
            if ( Utils<number>::OverlapBoundSC(R_cm_, Particle1_.Hull, Particle2_.Hull) )
            {
                Particle2_.SetRandomOrientation(rng_engine_, rng_distrib_);
                
                // Oriented Bounding Box (OBB) overlap test
                if ( Utils<number>::OverlapBoundOB(R_cm_, Particle1_.Hull, Particle2_.Hull) )
                {
                    mayer_interaction_ = IManager.MayerInteraction(R_cm_, &Particle1_, &Particle2_);
                    
                    if ( mayer_interaction_ == 0. ) PruneGrid(r_hard);
                }
            }
        }
    }
    
    number f_exc = Exc_grid_.sum() / ((number)NX*NY*NZ);
    number v_box = 8. * Particle1_.Hull->l_xh*Particle1_.Hull->l_yh*Particle1_.Hull->l_zh;
    
    return (v_box * f_exc);
}

#endif

// ============================
/* Random configuration generator */
// ============================
template<template<typename number> class ParticleType, typename number>
void MCIntegrator<ParticleType, number>::ConfigGenerator()
{
    mayer_interaction_ = 0.;
    
    while ( (mayer_interaction_ == 0.) && (ctr_mc_ < N_PER_PROC_) )
    {
        ctr_mc_++;
        
        SyncCheck(false);
        
        // Fetch random configurations if relevant
        Particle1_.Parse(rng_engine_, Bootstrap_map1_);
        Particle2_.Parse(rng_engine_, Bootstrap_map2_);
        
        // Random center-of-mass to center-of-mass separation vector
        R_cm_ = Vector3<number>::NullaryExpr([&](int) {return rng_distrib_(rng_engine_)-0.5;}) * 2.*IManager.R_INTEG;
        
        // Bounding sphere overlap test
        if ( R_cm_.norm() < IManager.R_INTEG )
        {
            Particle1_.SetRandomAxis(rng_engine_, rng_distrib_);
            Particle2_.SetRandomAxis(rng_engine_, rng_distrib_);
            
            // Bounding spherocylinder (SC) overlap test
            if ( Utils<number>::OverlapBoundSC(R_cm_, Particle1_.Hull, Particle2_.Hull) )
            {
                Particle1_.SetRandomOrientation(rng_engine_, rng_distrib_);
                Particle2_.SetRandomOrientation(rng_engine_, rng_distrib_);
                
                #if USE_RAPID
                    // Use RAPID library for collision detection if implemented
                    double pos1[3];
                    double pos2[3];
                    double rot1[3][3];
                    double rot2[3][3];
                    
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
                    if ( Utils<number>::OverlapBoundOB(R_cm_, Particle1_.Hull, Particle2_.Hull) )
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
/* MC integration of virial coefficients */
// ============================
template<template<typename number> class ParticleType, typename number>
void MCIntegrator<ParticleType, number>::VirialIntegrator(ArrayXX<number>* E_out, ArrayXX<number>* V_b,
                                                          ArrayX<number>* V_r, ArrayX<number>* V_l)
{
    V_b->setZero(N_THETA, N_THETA);
    
    V_r->setZero(N_THETA);
    V_l->setZero(N_THETA);
    
    E_out->setZero(E_DIM, E_DIM);

    LogTxt("------------");
    LogTxt("Integrating angle-dependant second-virial coefficients");
    
    MCReset();

    while ( ctr_mc_ < N_PER_PROC_ )
    {
        ConfigGenerator();
        
        if ( mayer_interaction_ != 0. )
        {
            ctr_ov_++;
            
            Vector3<number> U1 = Particle1_.U;
            Vector3<number> U2 = Particle2_.U;
            
            Vector3<number> V1 = Particle1_.V;
            Vector3<number> V2 = Particle2_.V;
            
            // Project V2 on the normal plane of U1
            number nu;
            number v2Eu1    = V2.dot(U1);
            
            V2             -= v2Eu1 * U1;
            number norm2    = V2.norm();
            
            if ( norm2 < TOL_SC ) nu = PI/2.;
            else                  nu = acos(V1.dot(V2) / norm2);
            
            // Work out configuration handedness
            number theta    = acos(U1.dot(U2));
            number deter    = R_cm_.dot(U1.cross(U2));
            
            uint idx_nu     = floor(nu    /PI * (number)N_THETA);
            uint idx_theta  = floor(theta /PI * (number)N_THETA);
            
            idx_nu          = fmin(idx_nu,    N_THETA-1);
            idx_theta       = fmin(idx_theta, N_THETA-1);
            
            // Fetch configuration angles
            number alpha1   = Particle1_.alpha;
            number alpha2   = Particle2_.alpha;
            
            number theta1   = Particle1_.theta;
            number theta2   = Particle2_.theta;
            
            number phi1     = Particle1_.phi;
            number phi2     = Particle2_.phi;
            
            uint idx_alpha1 = floor(alpha1/(2.*PI) * (number)N_ALPHA);
            uint idx_alpha2 = floor(alpha2/(2.*PI) * (number)N_ALPHA);
            
            uint idx_theta1 = floor(theta1/PI      * (number)N_THETA);
            uint idx_theta2 = floor(theta2/PI      * (number)N_THETA);
            
            uint idx_phi1   = floor(phi1/(2.*PI)   * (number)N_PHI);
            uint idx_phi2   = floor(phi2/(2.*PI)   * (number)N_PHI);
            
            // Avoid overflow for single-precision floats
            idx_alpha1      = fmin(idx_alpha1, N_ALPHA-1);
            idx_alpha2      = fmin(idx_alpha2, N_ALPHA-1);
            
            idx_theta1      = fmin(idx_theta1, N_THETA-1);
            idx_theta2      = fmin(idx_theta2, N_THETA-1);
            
            idx_phi1        = fmin(idx_phi1,   N_PHI-1);
            idx_phi2        = fmin(idx_phi2,   N_PHI-1);
            
            // Angular excluded volume integrals
            (*V_b)(idx_theta, idx_nu) += mayer_interaction_;
            
            if ( deter > 0. ) (*V_r)(idx_theta) += mayer_interaction_*sin(theta1)*sin(theta2);
            else              (*V_l)(idx_theta) += mayer_interaction_*sin(theta1)*sin(theta2);
            
            if ( ODF_TYPE == ODF_FULL ) (*E_out).at(idx_alpha1, idx_theta1, idx_phi1,
                                                    idx_alpha2, idx_theta2, idx_phi2) += mayer_interaction_;
            
            if ( ODF_TYPE == ODF_LEGENDRE )
            {
                for ( uint l1 = 0; l1 < N_L; l1 += IManager.N_DELTA_L )
                {
                    for ( uint l2 = 0; l2 < N_L; l2 += IManager.N_DELTA_L )
                    {
                        (*E_out)(l1,l2) += sin(theta1)*sin(theta2) * mayer_interaction_
                                         * sqrt((2.*(number)l1 + 1.)*(2.*(number)l2 + 1.)/4.)
                                         * Legendre<number>::Pn(l1, cos(theta2))
                                         * Legendre<number>::Pn(l2, cos(theta1));
                    }
                }
            }
        }
    }
    
    number prefactor = ((ODF_TYPE == ODF_FULL) ? 1./SQR(D_ALPHA*D_THETA*D_PHI) : 1.);
    
    *V_b   *= IManager.V_INTEG/N_PER_PROC_ / SQR(D_THETA);
    
    *V_r   /= N_PER_PROC_ * D_THETA;
    *V_l   /= N_PER_PROC_ * D_THETA;
    
    *E_out *= IManager.V_INTEG/N_PER_PROC_ * prefactor;
}

// ============================
/* MC integration of the cholesteric Legendre coefficients */
// ============================
template<template<typename number> class ParticleType, typename number>
void MCIntegrator<ParticleType, number>::LegendreIntegrator(ArrayXX<number>* E_out, number q_macro)
{
    Vector3<number> N_q;
    
    E_out->setZero(E_DIM, E_DIM);
    
    LogTxt("------------");
    LogTxt("Integrating Legendre-projected second-virial coefficients - q=%.5f", q_macro);
    
    MCReset();

    while ( ctr_mc_ < N_PER_PROC_ )
    {
        ConfigGenerator();
        
        if ( mayer_interaction_ != 0. )
        {
            ctr_ov_++;
            
            number theta1 = Particle1_.theta;
            number theta2 = Particle2_.theta;
            
            N_q << sin(q_macro*R_cm_(1)), 0., cos(q_macro*R_cm_(1));
            
            for ( uint l1 = 0; l1 < N_L; l1 += IManager.N_DELTA_L )
            {
                for ( uint l2 = 0; l2 < N_L; l2 += IManager.N_DELTA_L )
                {
                    (*E_out)(l1,l2) += sin(theta1)*sin(theta2) * mayer_interaction_
                                     * sqrt((2.*(number)l1 + 1.)*(2.*(number)l2 + 1.)/4.)
                                     * Legendre<number>::Pn(l1, Particle2_.U.dot(N_q))
                                     * Legendre<number>::Pn(l2, cos(theta1));
                }
            }
        }
    }
    
    *E_out *= IManager.V_INTEG/N_PER_PROC_;
}

// ============================
/* Work out preliminary chiral properties */
// ============================
template<template<typename number> class ParticleType, typename number>
void MCIntegrator<ParticleType, number>::FrankIntegrator(const ArrayXX<number>& Psi_in,
                                                         ArrayX<number>* K1_out, ArrayX<number>* K2_out, ArrayX<number>* K3_out,
                                                         ArrayX<number>* Kt_out)
{
    ArrayXX<number> Psi_dot(Psi_in.rows(), Psi_in.cols());

    ArrayX<number>  Eff_grid = Eta_grid * this->IManager.V_EFF/this->IManager.V0;
    
    ArrayX<number>  G_grid   = (1. - 3/4.*Eff_grid) / SQR(1. - Eff_grid);
    ArrayX<number>  N_grid   = Eta_grid * CUB(IManager.SIGMA_R)/IManager.V0;
    
    K1_out->setZero(N_STEPS_ETA);
    K2_out->setZero(N_STEPS_ETA);
    K3_out->setZero(N_STEPS_ETA);
    Kt_out->setZero(N_STEPS_ETA);
    
    // Work out Psi differentials
    for ( uint idx_alpha = 0; idx_alpha < N_ALPHA; ++idx_alpha )
    {
        for ( uint idx_theta = 0; idx_theta < N_THETA; ++idx_theta )
        {
            for ( uint idx_phi = 0; idx_phi < N_PHI; ++idx_phi )
            {
                if ( idx_theta < N_THETA-1 )
                {
                    Psi_dot.col_at(idx_alpha, idx_theta, idx_phi) = Psi_in.col_at(idx_alpha, idx_theta+1, idx_phi)
                                                                  - Psi_in.col_at(idx_alpha, idx_theta,   idx_phi);
                }
                
                else Psi_dot.col_at(idx_alpha, idx_theta, idx_phi) = -Psi_dot.col_at(idx_alpha, 0, idx_phi);
            }
        }
    }
    
    Psi_dot /= D_THETA;
    
    LogTxt("------------");
    LogTxt("Starting preliminary chiral run...");
    
    MCReset();

    while ( ctr_mc_ < N_PER_PROC_ )
    {
        ConfigGenerator();
        
        if ( mayer_interaction_ != 0. )
        {
            ctr_ov_++;

            // Fetch configuration angles
            number alpha1   = Particle1_.alpha;
            number alpha2   = Particle2_.alpha;
            
            number theta1   = Particle1_.theta;
            number theta2   = Particle2_.theta;
            
            number phi1     = Particle1_.phi;
            number phi2     = Particle2_.phi;
            
            uint idx_alpha1 = floor(alpha1/(2.*PI) * (number)N_ALPHA);
            uint idx_alpha2 = floor(alpha2/(2.*PI) * (number)N_ALPHA);
            
            uint idx_theta1 = floor(theta1/PI      * (number)N_THETA);
            uint idx_theta2 = floor(theta2/PI      * (number)N_THETA);
            
            uint idx_phi1   = floor(phi1/(2.*PI)   * (number)N_PHI);
            uint idx_phi2   = floor(phi2/(2.*PI)   * (number)N_PHI);
            
            // Avoid overflow for single-precision floats
            idx_alpha1      = fmin(idx_alpha1, N_ALPHA-1);
            idx_alpha2      = fmin(idx_alpha2, N_ALPHA-1);
            
            idx_theta1      = fmin(idx_theta1, N_THETA-1);
            idx_theta2      = fmin(idx_theta2, N_THETA-1);
            
            idx_phi1        = fmin(idx_phi1,   N_PHI-1);
            idx_phi2        = fmin(idx_phi2,   N_PHI-1);
            
            // Update Frank elastic moduli
            *K1_out -= mayer_interaction_
                     * Psi_dot.col_at(idx_alpha1, idx_theta1, idx_phi1)
                     * Psi_dot.col_at(idx_alpha2, idx_theta2, idx_phi2)
                     * Particle1_.U(0) * Particle2_.U(0)
                     * SQR(R_cm_(0));
			
            *K2_out -= mayer_interaction_
                     * Psi_dot.col_at(idx_alpha1, idx_theta1, idx_phi1)
                     * Psi_dot.col_at(idx_alpha2, idx_theta2, idx_phi2)
                     * Particle1_.U(1) * Particle2_.U(1)
                     * SQR(R_cm_(0));
            
            *K3_out -= mayer_interaction_
                     * Psi_dot.col_at(idx_alpha1, idx_theta1, idx_phi1)
                     * Psi_dot.col_at(idx_alpha2, idx_theta2, idx_phi2)
                     * Particle1_.U(0) * Particle2_.U(0)
                     * SQR(R_cm_(2));
            
            // Update chiral strength
            *Kt_out -= mayer_interaction_ * sin(theta1)
                     * Psi_in. col_at(idx_alpha1, idx_theta1, idx_phi1)
                     * Psi_dot.col_at(idx_alpha2, idx_theta2, idx_phi2)
                     * Particle2_.U(1)
                     * R_cm_(0);
        }
    }
    
    *K1_out *= SQR(N_grid)/2. * G_grid * IManager.V_INTEG/N_PER_PROC_ / pow(IManager.SIGMA_R, 5);
    *K2_out *= SQR(N_grid)/2. * G_grid * IManager.V_INTEG/N_PER_PROC_ / pow(IManager.SIGMA_R, 5);
    *K3_out *= SQR(N_grid)/2. * G_grid * IManager.V_INTEG/N_PER_PROC_ / pow(IManager.SIGMA_R, 5);
    
    *Kt_out *= SQR(N_grid)/2. * G_grid * IManager.V_INTEG/N_PER_PROC_ / pow(IManager.SIGMA_R, 4);
}

template class MCIntegrator<MESOGEN, float>;
template class MCIntegrator<MESOGEN, double>;
