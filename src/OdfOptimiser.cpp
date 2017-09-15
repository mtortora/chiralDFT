// ===================================================================
/**
 * Orientation distribution function convenience class - always uses
 * double-precision ODF calculations for numerical stability.
 */
// ===================================================================
/*
 * OdfOptimiser.cpp: Version 2.5
 * Created 16/09/2015 by Maxime Tortora
 */
// ===================================================================

#include "Legendre.hpp"
#include "OdfOptimiser.hpp"
#include "SphericalFunctions/WignerDMatrices.hpp"


// ============================
/* Class constructor */
// ============================
template<template<typename number> class ParticleType, typename number>
OdfOptimiser<ParticleType, number>::OdfOptimiser()
{
    Alpha_grid = ArrayX<number>::LinSpaced(N_ALPHA, D_ALPHA/2., 2.*PI-D_ALPHA/2.);
    Theta_grid = ArrayX<number>::LinSpaced(N_THETA, D_THETA/2., 1.*PI-D_THETA/2.);
    Phi_grid   = ArrayX<number>::LinSpaced(N_PHI,   D_PHI/2.,   2.*PI-D_PHI/2.);

    Psi_iso_   = ArrayX<double>::Constant(N_ALPHA*N_THETA*N_PHI, 1./(8.*SQR(PI)));
}

// ============================
/* Work out the ODF simulation range through binodal decomposition analysis */
// ============================
template<template<typename number> class ParticleType, typename number>
void OdfOptimiser<ParticleType, number>::BinodalAnalysis(const ArrayXX<number>& E_in, int mpi_rank, int mpi_size)
{
    LogTxt("------------");
    LogGre("Working out the isotropic/nematic binodals...");
    
    uint   ctr_nr(0);

    bool   is_stable(true);
    bool   converged(false);

    double eta_min(ETA_MIN);
    double eta_max(ETA_MAX);
    
    double binodal_iso;
    double binodal_nem;
    
    Vector2<double> Binodals;
    Vector2<double> D_Binodals;
    
    Binodals << eta_min, eta_max;
    
    // Solve the coexistence equations using a relaxed Newton-Raphson (NR) method
    while ( is_stable && (!converged) && (ctr_nr < N_STEPS_NR) )
    {
        try
        {
            double eta_iso1 = Binodals(0);
            double eta_nem1 = Binodals(1);
            
            double eta_iso2 = eta_iso1 + TOL_BIN;
            double eta_nem2 = eta_nem1 - TOL_BIN;
            
            // Enforce density bounds for NR convergence
            if ( (Binodals.array() < 0.).any() || (Binodals.array() > eta_max).any() )
            {
                Binodals(0) = eta_min;
                
                throw std::runtime_error("Newton-Raphson method went out of bounds");
            }
            
            ArrayX<double> Psi_nem1 = SequentialOptimiser(eta_nem1, E_in, mpi_rank, mpi_size);
            ArrayX<double> Psi_nem2 = SequentialOptimiser(eta_nem2, E_in, mpi_rank, mpi_size);
            
            ctr_nr++;

            // Revert to higher density if upper binodal is not in the nematic range
            if ( (Psi_nem1 - Psi_iso_).abs().maxCoeff() < TOL_NEM )
            {
                if ( Binodals(1) == eta_max ) is_stable = false;
                else Binodals(1) -= GAMMA_NR * D_Binodals(1);
            }
            
            else
            {
                // Thermo_Xi is the vector [p, mu] for phase X at density i
                Vector2<double> Thermo_iso1 = ODFThermo(eta_iso1, Psi_iso_, E_in);
                Vector2<double> Thermo_iso2 = ODFThermo(eta_iso2, Psi_iso_, E_in);

                Vector2<double> Thermo_nem1 = ODFThermo(eta_nem1, Psi_nem1, E_in);
                Vector2<double> Thermo_nem2 = ODFThermo(eta_nem2, Psi_nem2, E_in);
                
                // Compute NR Jacobian through pressure/potential differentials
                Matrix22<double> Jcb;
                Vector2 <double> D_Thermo;
                
                Jcb.col(0)   = (Thermo_iso2 - Thermo_iso1) / TOL_BIN;
                Jcb.col(1)   = (Thermo_nem2 - Thermo_nem1) / TOL_BIN;

                // Solve direct NR equation using Householder QR decomposition
                D_Thermo     = Thermo_nem1 - Thermo_iso1;
                D_Binodals   = Jcb.colPivHouseholderQr().solve(D_Thermo);
                
                double max_b = D_Binodals.cwiseAbs().maxCoeff();
                double max_t = D_Thermo  .cwiseAbs().maxCoeff();
                
                // Relaxed NR iteration - classical NR can be recovered by setting GAMMA_NR to 1
                Binodals    += GAMMA_NR * D_Binodals;
                converged    = ((max_b < TOL_BIN) && (max_t < TOL_BIN)) ;
            }
        }
        
        // Handle ODF convergence failures
        catch ( std::exception& e )
        {
            eta_max    -= (ETA_MAX - ETA_MIN) / (N_STEPS_ETA - 1.);
            Binodals(1) = eta_max;
            
            if ( Binodals(1) <= 0. ) throw std::runtime_error("Couldn't work out ODFs");
        }
    }
    
    if ( !is_stable ) throw std::runtime_error("No stable nematic phase found in this density range");

    LogTxt("------------");

    // Print binodal analysis summary
    if ( converged )
    {
        binodal_iso = Binodals(0);
        binodal_nem = Binodals(1);
        
        LogCya("Relaxed Newton-Raphson method converged in %u steps", ctr_nr);
        LogTxt("I/N coexistence range: [%.2f%%, %.2f%%]", 100.*binodal_iso, 100.*binodal_nem);
        
        eta_min = binodal_nem;
    }
    
    else {LogRed("Unable to work out I/N spinodals");}
    LogTxt("Simulation range: [%.2f%%, %.2f%%]", 100.*eta_min, 100.*eta_max);
    
    this->Eta_grid = ArrayX<number>::LinSpaced(N_STEPS_ETA, eta_min, eta_max);
}

// ============================
/* ODF sequential solver - requires double precision to converge */
// ============================
template<template<typename number> class ParticleType, typename number>
ArrayX<double> OdfOptimiser<ParticleType, number>::SequentialOptimiser(number eta, const ArrayXX<number>& E_in,
                                                                       int mpi_rank, int mpi_size)
{
    bool   converged(false);
    
    double n_dens  = eta / this->IManager.V0;
    double eta_eff = eta * this->IManager.V_EFF/this->IManager.V0;
    double g_pl    = (1. - 3/4.*eta_eff) / SQR(1. - eta_eff);
    
    // Gaussian initial guess
    ArrayX<double> Psi (N_ALPHA*N_THETA*N_PHI);
    ArrayX<double> SinT(N_ALPHA*N_THETA*N_PHI);
    
    for ( uint idx_alpha = 0; idx_alpha < N_ALPHA; ++idx_alpha )
    {
        for ( uint idx_theta = 0; idx_theta < N_THETA; ++idx_theta )
        {
            for ( uint idx_phi = 0; idx_phi < N_PHI; ++idx_phi )
            {
                double theta = Theta_grid(idx_theta);

                Psi .at(idx_alpha, idx_theta, idx_phi) = exp(-SQR(theta));
                SinT.at(idx_alpha, idx_theta, idx_phi) = sin(theta);
            }
        }
    }
    
    for ( uint ctr_odf = 0; ctr_odf < N_STEPS_ODF; ++ctr_odf )
    {
        ArrayX<double> Psi_dummy = Psi;
        
        // Self-consistency equation iteration
        if ( ODF_TYPE == ODF_FULL )
        {
            Psi.setZero();
            
            // Grid-parallelised ODF optimisation
            int idx(0);
            
            for  ( uint idx_alpha1 = 0; idx_alpha1 < N_ALPHA; ++idx_alpha1 )
            {
                for ( uint idx_theta1 = 0; idx_theta1 < N_THETA; ++idx_theta1 )
                {
                    for ( uint idx_phi1 = 0; idx_phi1 < N_PHI; ++idx_phi1 )
                    {
                        if ( idx % mpi_size == mpi_rank )
                        {
                            // Compact inner sum for Eigen vectorisation
                            double onsager_sum = (E_in.col_at(idx_alpha1, idx_theta1, idx_phi1).template cast<double>()
                                                  * Psi_dummy * SinT).sum()
                                               * D_ALPHA*D_THETA*D_PHI;
                
                            Psi(idx) = exp(-g_pl * n_dens * onsager_sum);
                        }
                        
                        ++idx;
                    }
                }
            }
            
            MPI_Allreduce(MPI_IN_PLACE, Psi.data(), Psi.size(), Utils<double>().MPI_type, MPI_SUM, MPI_COMM_WORLD);
        }
        
        else if ( ODF_TYPE == ODF_LEGENDRE )
        {
            ArrayX<double> Psi_l = LegendreCoeffs(Psi_dummy);
            
            for ( uint idx_theta = 0; idx_theta < N_THETA; ++idx_theta )
            {
                double theta        = Theta_grid(idx_theta);
                double legendre_sum = 0.;
                
                for ( uint l1 = 0; l1 < N_L; l1 += this->IManager.N_DELTA_L )
                {
                    for ( uint l2 = 0; l2 < N_L; l2 += this->IManager.N_DELTA_L )
                    {
                        legendre_sum += (E_in(l1, l2) + E_in(l2, l1))/2.
                                      * sqrt((2.*(double)l1 + 1.)/2.) * Legendre<double>::Pn(l1, cos(theta))
                                      * Psi_l(l2);
                    }
                }
                
                Psi(idx_theta) = exp(-g_pl * n_dens/SQR(2.*PI) * legendre_sum);
            }
        }
        
        // Renormalisation & convergence check, enforcing head-tail symmetry if relevant
        if ( !IS_BIAXIAL ) Psi = (Psi + Psi.reverse()) / 2.;
        
        double norm = (Psi * SinT).sum() * D_ALPHA*D_THETA*D_PHI;
        
        Psi        /= norm;
        double err  = abs(Psi - Psi_dummy).maxCoeff();
        
        if ( err < TOL_ODF )
        {
            LogInf("ODF converged in %d steps - eta=%.4f", ctr_odf, eta);
            converged = true;
            
            break;
        }
    }
    
    if ( !converged )
    {
        LogErr("Reached maximum number of iterations - eta=%.4f", eta);
        throw std::runtime_error("ODF convergence failed");
    }
    
    return Psi;
}

// ============================
/* Legendre Psi_l coefficients */
// ============================
template<template<typename number> class ParticleType, typename number>
ArrayX<double> OdfOptimiser<ParticleType, number>::LegendreCoeffs(const ArrayX<double>& Psi_in)
{
    ArrayX<double> Psi_l(N_L);
    
    for ( uint l = 0; l < N_L; l += this->IManager.N_DELTA_L )
    {
        double legendre_sum = 0.;
        
        for ( uint idx_theta = 0; idx_theta < N_THETA; ++idx_theta )
        {
            double theta  = Theta_grid(idx_theta);
            
            legendre_sum += sqrt((2.*(double)l + 1.)/2.) * Legendre<double>::Pn(l, cos(theta))
                          * Psi_in(idx_theta) * sin(theta) * D_THETA;
        }
        
        Psi_l(l) = legendre_sum;
    }
    
    return Psi_l;
}

// ============================
/* Rotational entropy */
// ============================
template<template<typename number> class ParticleType, typename number>
double OdfOptimiser<ParticleType, number>::RotationalEnt(const ArrayX<double>& Psi_in)
{
    double s_rot(0.);
    
    // Compute Shannon entropy of the ODF
    for ( uint idx_alpha = 0; idx_alpha < N_ALPHA; ++idx_alpha )
    {
        for ( uint idx_theta = 0; idx_theta < N_THETA; ++idx_theta )
        {
            for ( uint idx_phi = 0; idx_phi < N_PHI; ++idx_phi )
            {
                double theta = Theta_grid(idx_theta);

                s_rot += Psi_in.at(idx_alpha, idx_theta, idx_phi)
                       * log(Psi_in.at(idx_alpha, idx_theta, idx_phi))
                       * sin(theta) * D_ALPHA*D_THETA*D_PHI;
            }
        }
    }
    
    return s_rot;
}

// ============================
/* Nematic order parameters */
// ============================
template<template<typename number> class ParticleType, typename number>
void OdfOptimiser<ParticleType, number>::OrderParams(const ArrayX<double>& Psi_in,
                                                     Matrix33<double>* Frame, ArrayX<std::complex<double> >* S,
                                                     int mpi_rank, int mpi_size)
{
    S->setZero(SIZE_L);
    
    Matrix33<double> Qu = Matrix33<double>::Zero();

    SphericalFunctions::WignerDMatrix WMat;

    // Averaged long-axis order parameter tensor Q
    for ( uint idx_alpha = 0; idx_alpha < N_ALPHA; ++idx_alpha )
    {
        for ( uint idx_theta = 0; idx_theta < N_THETA; ++idx_theta )
        {
            for ( uint idx_phi = 0; idx_phi < N_PHI; ++idx_phi )
            {
                double alpha = Alpha_grid(idx_alpha);
                double theta = Theta_grid(idx_theta);
                double phi   = Phi_grid  (idx_phi);

                this->IManager.SetU(theta, phi);
                Vector3<double> U = this->IManager.U.template cast<double>();
                
                // Averaged Q-tensor
                Qu(0,0) += U(0)*U(0) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                Qu(1,0) += U(1)*U(0) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                Qu(2,0) += U(2)*U(0) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                Qu(0,1) += U(0)*U(1) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                Qu(1,1) += U(1)*U(1) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                Qu(2,1) += U(2)*U(1) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                Qu(0,2) += U(0)*U(2) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                Qu(1,2) += U(1)*U(2) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                Qu(2,2) += U(2)*U(2) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                
                // Order parameters
                int idx(0);
                
                WMat.SetRotation(phi, theta, alpha);

                for ( int idx_l = 0; idx_l < N_L+1; ++idx_l )
                {
                    for ( int idx_mp = -idx_l; idx_mp < idx_l+1; ++idx_mp )
                    {
                        for ( int idx_m = -idx_l; idx_m < idx_l+1; ++idx_m )
                        {
                            if ( idx % mpi_size == mpi_rank )
                            {
                                (*S)(idx) += WMat(idx_l,idx_mp,idx_m) * Psi_in.at(idx_alpha,idx_theta,idx_phi) * sin(theta) * D_ALPHA*D_THETA*D_PHI;
                            }
                            
                            ++idx;
                        }
                        
                    }
                }
            }
        }
    }
    
    MPI_Allreduce(MPI_IN_PLACE, S->data(), S->size(), MPI_CXX_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

    // Normalise Q tensor
    Qu = 3/2.*Qu - 1/2.*Matrix33<double>::Identity();
    
    // Q spectral analysis
    Eigen::SelfAdjointEigenSolver<Matrix33<double> > Su(Qu);
    
    // Work out nematic frame
    *Frame = Su.eigenvectors().rowwise().reverse();
    
    if ( Frame->determinant() < 0. ) Frame->col(0) *= -1.;
}

// ============================
/* Second-virial coefficient */
// ============================
template<template<typename number> class ParticleType, typename number>
double OdfOptimiser<ParticleType, number>::VirialCoeff(const ArrayX<double>& Psi_in, const ArrayXX<number>& E_in)
{
    double b2(0.);
    
    // Orientation-independant second virial coefficient
    if ( ODF_TYPE == ODF_FULL )
    {
        for  ( uint idx_alpha1 = 0; idx_alpha1 < N_ALPHA; ++idx_alpha1 )
        {
            for ( uint idx_theta1 = 0; idx_theta1 < N_THETA; ++idx_theta1 )
            {
                for ( uint idx_phi1 = 0; idx_phi1 < N_PHI; ++idx_phi1 )
                {
                    for  ( uint idx_alpha2 = 0; idx_alpha2 < N_ALPHA; ++idx_alpha2 )
                    {
                        for ( uint idx_theta2 = 0; idx_theta2 < N_THETA; ++idx_theta2 )
                        {
                            for ( uint idx_phi2 = 0; idx_phi2 < N_PHI; ++idx_phi2 )
                            {
                                double theta1 = Theta_grid(idx_theta1);
                                double theta2 = Theta_grid(idx_theta2);
                                
                                b2 += E_in.at(idx_alpha1, idx_theta1, idx_phi1,
                                              idx_alpha2, idx_theta2, idx_phi2)/2.
                                    * Psi_in.at(idx_alpha1, idx_theta1, idx_phi1)
                                    * Psi_in.at(idx_alpha2, idx_theta2, idx_phi2)
                                    * sin(theta1)*sin(theta2) * SQR(D_ALPHA*D_THETA*D_PHI);
                            }
                        }
                    }
                }
            }
        }
    }
    
    // Legendre-projected second virial coefficient
    if ( ODF_TYPE == ODF_LEGENDRE )
    {
        ArrayX<double> Psi_l = LegendreCoeffs(Psi_in);
        
        for ( uint l1 = 0; l1 < N_L; l1 += this->IManager.N_DELTA_L )
        {
            for ( uint l2 = 0; l2 < N_L; l2 += this->IManager.N_DELTA_L )
            {
                b2 += Psi_l(l1)*Psi_l(l2) * E_in(l1, l2)/2.;
            }
        }
    }
    
    return b2;
}

// ============================
/* Onsager-PL free energy */
// ============================
template<template<typename number> class ParticleType, typename number>
double OdfOptimiser<ParticleType, number>::FreeEnergy(number eta, const ArrayX<double>& Psi_in, const ArrayXX<number>& E_in)
{
    double b2      = VirialCoeff(Psi_in, E_in);
    double eta_eff = eta * this->IManager.V_EFF/this->IManager.V0;
    
    double g_pl    = (1. - 3/4.*eta_eff) / SQR(1. - eta_eff);
    double n_resc  = eta * CUB(this->IManager.SIGMA_R)/this->IManager.V0;
    
    double s_rot   = RotationalEnt(Psi_in);

    // Ideal and excess free energies
    double f_id    = n_resc * s_rot;
    double f_exc   = b2/this->IManager.V0 * n_resc * eta*g_pl;
    double f_tot   = f_id + f_exc;
    
    return f_tot;
}

// ============================
/* Work out osmotic pressure and chemical potential */
// ============================
template<template<typename number> class ParticleType, typename number>
Vector2<double> OdfOptimiser<ParticleType, number>::ODFThermo(number eta, const ArrayX<double>& Psi_in,
                                                              const ArrayXX<number>& E_in)
{
    Vector2<double> Thermo;

    double b2      = VirialCoeff(Psi_in, E_in);
    double eta_eff = eta * this->IManager.V_EFF/this->IManager.V0;

    double s_rot   = RotationalEnt(Psi_in);

    // Decompose second-virial thermodynamics into ideal and excess contributions
    double p_id    = eta;
    double mu_id   = log(eta) + s_rot;
    
    double p_exc   = b2/this->IManager.V0 * eta * (eta_eff - SQR(eta_eff)/2.) / CUB(1.-eta_eff);
    double mu_exc  = b2/this->IManager.V0 * eta * (8. - 9.*eta_eff + 3.*SQR(eta_eff)) / (4.*CUB(1.-eta_eff));
    
    // Adimensioned total osmotic pressure & chemical potential
    Thermo << p_id+p_exc, mu_id+mu_exc;
    
    return Thermo;
}

// ============================
/* Energy grid iterator */
// ============================
template<template<typename number> class ParticleType, typename number>
void OdfOptimiser<ParticleType, number>::EnergyGrid(const ArrayXX<number>& E_in, ArrayX<number>* F_out,
                                                    int mpi_rank, int mpi_size)
{
    LogTxt("------------");
    LogGre("Computing free energies...");
    
    F_out->setZero(N_STEPS_ETA);
    
    for ( uint idx_eta = 0; idx_eta < N_STEPS_ETA; ++idx_eta )
    {
        number eta = this->Eta_grid(idx_eta);
        
        try
        {
            ArrayX<double> Psi = SequentialOptimiser(eta, E_in, mpi_rank, mpi_size);
            (*F_out)(idx_eta)  = FreeEnergy(eta, Psi, E_in);
        }
        
        catch ( std::exception& e ) {LogErr("%s - discard higher density data", e.what());}
    }
}

// ============================
/* ODF grid iterator */
// ============================
template<template<typename number> class ParticleType, typename number>
void OdfOptimiser<ParticleType, number>::ODFGrid(const ArrayXX<number>& E_in,
                                                 ArrayX<number>* P_out, ArrayX<number>* Mu_out, ArrayX<number>* F_out,
                                                 ArrayXX<std::complex<number> >* S_out, ArrayXX<number>* Psi_out,
                                                 int mpi_rank, int mpi_size)
{
    LogTxt("------------");
    LogGre("Computing equilibrium ODFs and %d Wigner order parameters...", SIZE_L);

    *P_out   = ArrayX<number>(N_STEPS_ETA);
    *Mu_out  = ArrayX<number>(N_STEPS_ETA);
    *F_out   = ArrayX<number>(N_STEPS_ETA);
    
    *S_out   = ArrayXX<number>(N_STEPS_ETA, SIZE_L);
    *Psi_out = ArrayXX<number>(N_STEPS_ETA, N_ALPHA*N_THETA*N_PHI);

    for ( uint idx_eta = 0; idx_eta < N_STEPS_ETA; ++idx_eta )
    {
        number eta             = this->Eta_grid(idx_eta);
        
        // Equilibrium ODF
        ArrayX <double> Psi    = SequentialOptimiser(eta, E_in, mpi_rank, mpi_size);
        Vector2<double> Thermo = ODFThermo(eta, Psi, E_in);
        
        Matrix33<double> Frame;
        ArrayX<std::complex<double> > S;
        
        OrderParams(Psi, &Frame, &S, mpi_rank, mpi_size);
        
        // Equilibrium pressure, chemical potential & free energy
        (*P_out) (idx_eta)    = Thermo(0);
        (*Mu_out)(idx_eta)    = Thermo(1);
        (*F_out) (idx_eta)    = FreeEnergy(eta, Psi, E_in);
        
        // Nematic order parameter & density-dependent ODF
        S_out  ->row(idx_eta) = S  .cast<std::complex<number> >();
        Psi_out->row(idx_eta) = Psi.cast<number>();
    }
}

template class OdfOptimiser<MESOGEN, float>;
template class OdfOptimiser<MESOGEN, double>;
