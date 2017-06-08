// ===================================================================
/**
 * Orientation distribution function derived handler class - objects
 * and methods preceded by "this->" are inherited from MCIntegrator.
 */
// ===================================================================
/*
 * OdfOptimiser.cpp: Version 2.5
 * Created 16/09/2015 by Maxime Tortora
 */
// ===================================================================

#include "Legendre.hpp"
#include "OdfOptimiser.hpp"


// ============================
/* Class constructor */
// ============================
template<template<typename number> class ParticleType, typename number>
OdfOptimiser<ParticleType, number>::OdfOptimiser()
{
    Theta_grid = ArrayX<number>::LinSpaced(N_STEPS_THETA, 0., PI);
    Psi_iso_   = ArrayX<number>::Constant (N_STEPS_THETA, 1./(8.*SQR(PI)));
}

// ============================
/* Work out the ODF simulation range through binodal decomposition analysis */
// ============================
template<template<typename number> class ParticleType, typename number>
void OdfOptimiser<ParticleType, number>::BinodalAnalysis(const MatrixXX<number>& E_in)
{
    LogTxt("------------");
    LogGre("Working out the isotropic/nematic binodals...");
    
    uint     ctr_nr(0);

    bool     is_stable(true);
    bool     converged(false);

    number   eta_min(ETA_MIN);
    number   eta_max(ETA_MAX);
    
    number   binodal_iso;
    number   binodal_nem;
    
    Vector2<number> Binodals;
    Vector2<number> D_Binodals;
    
    Binodals << eta_min, eta_max;
    
    // Solve the coexistence equations using a relaxed Newton-Raphson (NR) method
    while ( is_stable && (!converged) && (ctr_nr < N_STEPS_NR) )
    {
        try
        {
            number  eta_iso1 = Binodals(0);
            number  eta_nem1 = Binodals(1);
            
            number  eta_iso2 = eta_iso1 + TOL_BIN;
            number  eta_nem2 = eta_nem1 - TOL_BIN;
            
            ArrayX<number> Psi_nem1 = SequentialOptimiser(eta_nem1, E_in);
            ArrayX<number> Psi_nem2 = SequentialOptimiser(eta_nem2, E_in);
            
            // Revert to higher density if upper binodal is not in the nematic range
            if ( (Psi_nem1 - Psi_iso_).abs().maxCoeff() < TOL_NEM )
            {
                if ( Binodals(1) == eta_max ) is_stable = false;
                else Binodals(1) -= GAMMA_NR * D_Binodals(1);
            }
            
            else
            {
                // Enforce density bounds for NR convergence
                if ( (Binodals.array() < 0.).any() || (Binodals.array() > eta_max).any() )
                {
                    Binodals(0) = eta_min;
                    
                    throw std::runtime_error("Newton-Raphson method went out of bounds");
                }

                // Thermo_Xi is the vector [p, mu] for phase X at density i
                Vector2<number> Thermo_iso1 = ODFThermo(eta_iso1, Psi_iso_, E_in);
                Vector2<number> Thermo_iso2 = ODFThermo(eta_iso2, Psi_iso_, E_in);

                Vector2<number> Thermo_nem1 = ODFThermo(eta_nem1, Psi_nem1, E_in);
                Vector2<number> Thermo_nem2 = ODFThermo(eta_nem2, Psi_nem2, E_in);
                
                // Compute NR Jacobian through pressure/potential differentials
                Matrix22<number> Jcb;

                Jcb.col(0)           = (Thermo_iso2 - Thermo_iso1) / TOL_BIN;
                Jcb.col(1)           = (Thermo_nem2 - Thermo_nem1) / TOL_BIN;

                // Solve direct NR equation using Householder QR decomposition
                Vector2<number> D_Thermo    = Thermo_nem1 - Thermo_iso1;
                D_Binodals           = Jcb.colPivHouseholderQr().solve(D_Thermo);
                
                number max_b         = D_Binodals.cwiseAbs().maxCoeff();
                number max_t         = D_Thermo  .cwiseAbs().maxCoeff();
                
                // Relaxed NR iteration - classical NR can be recovered by setting GAMMA_NR to 1
                Binodals            += GAMMA_NR * D_Binodals;
                converged            = ((max_b < TOL_BIN) && (max_t < TOL_BIN)) ;
                
                ctr_nr++;
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
/* ODF sequential solver */
// ============================
template<template<typename number> class ParticleType, typename number>
ArrayX<number> OdfOptimiser<ParticleType, number>::SequentialOptimiser(number eta, const MatrixXX<number>& E_in)
{
    bool   converged(false);

    number n_dens  = eta / this->IManager.V0;
    number eta_eff = eta * this->IManager.V_EFF/this->IManager.V0;
    number g_pl    = (1. - 3/4.*eta_eff) / SQR(1. - eta_eff);
    
    // Gaussian initial guess
    ArrayX<number> Psi = exp(-Theta_grid*Theta_grid);
    
    for ( uint ctr_odf = 0; ctr_odf < N_STEPS_ODF; ++ctr_odf )
    {
        ArrayX<number> Psi_dummy = Psi;
        
        // Self-consistency equation iteration
        if ( ODF_TYPE == ODF_FULL )
        {
            for ( uint idx_theta1 = 0; idx_theta1 < N_STEPS_THETA; ++idx_theta1 )
            {
                number onsager_sum = 0.;
                
                for ( uint idx_theta2 = 0; idx_theta2 < N_STEPS_THETA; ++idx_theta2 )
                {
                    number theta2 = Theta_grid(idx_theta2);
                    
                    onsager_sum  += (E_in(idx_theta1, idx_theta2) + E_in(idx_theta2, idx_theta1))/2.
                                  * sin(theta2) * D_THETA
                                  * Psi_dummy(idx_theta2);
                }
                
                Psi(idx_theta1) = exp(-g_pl * n_dens/SQR(2.*PI) * onsager_sum);
            }
        }
        
        if ( ODF_TYPE == ODF_LEGENDRE )
        {
            ArrayX<number> Psi_l = LegendreCoeffs(Psi);
            
            for ( uint idx_theta = 0; idx_theta < N_STEPS_THETA; ++idx_theta )
            {
                number theta        = Theta_grid(idx_theta);
                number legendre_sum = 0.;
            
                for ( uint l1 = 0; l1 < N_L; l1 += this->IManager.N_DELTA_L )
                {
                    for ( uint l2 = 0; l2 < N_L; l2 += this->IManager.N_DELTA_L )
                    {
                        legendre_sum += (E_in(l1, l2) + E_in(l2, l1))/2.
                                      * sqrt((2.*(number)l1 + 1.)/2.) * Legendre<number>::Pn(l1, cos(theta))
                                      * Psi_l(l2);
                    }
                }
            
                Psi(idx_theta) = exp(-g_pl * n_dens/SQR(2.*PI) * legendre_sum);
            }
        }
        
        // Renormalisation & convergence check, enforcing head-tail symmetry
        Psi           = (Psi + Psi.reverse()) / 2.;
        
        ArrayX<number> Dummy = Psi * sin(Theta_grid);
        Psi          /= SQR(2.*PI)*D_THETA * Dummy.sum();
        
        number err    = abs(Psi - Psi_dummy).maxCoeff();
        
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
ArrayX<number> OdfOptimiser<ParticleType, number>::LegendreCoeffs(const ArrayX<number>& Psi)
{
    ArrayX<number> Psi_l(N_L);
    
    for ( uint l = 0; l < N_L; l += this->IManager.N_DELTA_L )
    {
        number legendre_sum = 0.;
        
        for ( uint idx_theta = 0; idx_theta < N_STEPS_THETA; ++idx_theta )
        {
            number theta  = Theta_grid(idx_theta);
            legendre_sum += sqrt((2.*(number)l + 1.)/2.) * Legendre<number>::Pn(l, cos(theta))
                          * sin(theta) * Psi(idx_theta) * D_THETA;
        }
        
        Psi_l(l) = legendre_sum;
    }
    
    return Psi_l;
}

// ============================
/* Second-virial coefficient */
// ============================
template<template<typename number> class ParticleType, typename number>
number OdfOptimiser<ParticleType, number>::VirialCoeff(const ArrayX<number>& Psi_in, const MatrixXX<number>& E_in)
{
    number b2(0.);
    
    // Orientation-independant second virial coefficient
    if ( ODF_TYPE == ODF_FULL )
    {
        for ( uint idx_theta1 = 0; idx_theta1 < N_STEPS_THETA; ++idx_theta1 )
        {
            for ( uint idx_theta2 = 0; idx_theta2 < N_STEPS_THETA; ++idx_theta2 )
            {
                number theta1 = Theta_grid(idx_theta1);
                number theta2 = Theta_grid(idx_theta2);
                
                b2 += sin(theta1)*sin(theta2) * SQR(D_THETA)
                    * Psi_in(idx_theta1)*Psi_in(idx_theta2)
                    * E_in(idx_theta1, idx_theta2)/2.;
            }
        }
    }
    
    // Legendre-projected second virial coefficient
    if ( ODF_TYPE == ODF_LEGENDRE )
    {
        ArrayX<number> Psi_l = LegendreCoeffs(Psi_in);
        
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
number OdfOptimiser<ParticleType, number>::FreeEnergy(number eta, const ArrayX<number>& Psi_in, const MatrixXX<number>& E_in)
{
    number b2      = VirialCoeff(Psi_in, E_in);
    number eta_eff = eta * this->IManager.V_EFF/this->IManager.V0;
    
    number g_pl    = (1. - 3/4.*eta_eff) / SQR(1. - eta_eff);
    number n_resc  = eta * CUB(this->IManager.SIGMA_R)/this->IManager.V0;
    
    // Ideal and excess free energies
    number f_id    = n_resc * SQR(2.*PI)*D_THETA * (sin(Theta_grid) * Psi_in * log(Psi_in)).sum();
    number f_exc   = b2/this->IManager.V0 * n_resc * eta*g_pl;
    number f_tot   = f_id + f_exc;
    
    return f_tot;
}

// ============================
/* Work out osmotic pressure and chemical potential */
// ============================
template<template<typename number> class ParticleType, typename number>
Vector2<number> OdfOptimiser<ParticleType, number>::ODFThermo(number eta, const ArrayX<number>& Psi_in, const MatrixXX<number>& E_in)
{
    Vector2<number> Thermo;

    number b2      = VirialCoeff(Psi_in, E_in);
    number eta_eff = eta * this->IManager.V_EFF/this->IManager.V0;

    // Decompose second-virial thermodynamics into ideal and excess contributions
    number p_id    = eta;
    number mu_id   = log(eta) + SQR(2.*PI)*D_THETA * (sin(Theta_grid) * Psi_in * log(Psi_in)).sum();
    
    number p_exc   = b2/this->IManager.V0 * eta * (eta_eff - SQR(eta_eff)/2.) / CUB(1.-eta_eff);
    number mu_exc  = b2/this->IManager.V0 * eta * (8. - 9.*eta_eff + 3.*SQR(eta_eff)) / (4.*CUB(1.-eta_eff));
    
    // Adimensioned total osmotic pressure & chemical potential
    Thermo << p_id+p_exc, mu_id+mu_exc;
    
    return Thermo;
}

// ============================
/* Energy grid iterator */
// ============================
template<template<typename number> class ParticleType, typename number>
void OdfOptimiser<ParticleType, number>::EnergyGrid(const MatrixXX<number>& E_in, ArrayX<number>* F_out)
{
    LogTxt("------------");
    LogGre("Computing free energies...");
    
    ArrayX<number> F_loc = ArrayX<number>::Zero(N_STEPS_ETA);
    
    for ( uint idx_eta = 0; idx_eta < N_STEPS_ETA; ++idx_eta )
    {
        number eta = this->Eta_grid(idx_eta);
        
        try
        {
            ArrayX<number> Psi    = SequentialOptimiser(eta, E_in);
            F_loc(idx_eta) = FreeEnergy(eta, Psi, E_in);
        }
        
        catch ( std::exception& e ) {LogErr("%s - discard higher density data", e.what());}
    }
    
    *F_out = F_loc;
}

// ============================
/* ODF grid iterator */
// ============================
template<template<typename number> class ParticleType, typename number>
void OdfOptimiser<ParticleType, number>::ODFGrid(const MatrixXX<number>& E_in, ArrayX<number>* P_out, ArrayX<number>* Mu_out,
                                                 ArrayX<number>* F_out, ArrayX<number>* S_out, ArrayXX<number>* Psi_out)
{
    LogTxt("------------");
    LogGre("Computing equilibrium ODFs and nematic order parameters...");

    ArrayX<number>  P_loc (N_STEPS_ETA);
    ArrayX<number>  Mu_loc(N_STEPS_ETA);
    ArrayX<number>  F_loc (N_STEPS_ETA);
    ArrayX<number>  S_loc (N_STEPS_ETA);
    
    ArrayXX<number> Psi_grd(N_STEPS_ETA, N_STEPS_THETA);
    
    for ( uint idx_eta = 0; idx_eta < N_STEPS_ETA; ++idx_eta )
    {
        number   eta         = this->Eta_grid(idx_eta);
        
        ArrayX<number>  Psi         = SequentialOptimiser(eta, E_in);
        Vector2<number> Thermo      = ODFThermo(eta, Psi, E_in);
        
        // Equilibrium pressure, chemical potential & free energy
        P_loc(idx_eta)       = Thermo(0);
        Mu_loc(idx_eta)      = Thermo(1);
        F_loc(idx_eta)       = FreeEnergy(eta, Psi, E_in);
        
        // Nematic order parameter & density-dependent ODF
        ArrayX<number> S_dummy      = (3.*SQR(cos(Theta_grid)) - 1.) / 2.;
        S_dummy             *= SQR(2.*PI)*D_THETA * sin(Theta_grid) * Psi;
        
        S_loc(idx_eta)       = S_dummy.sum();
        Psi_grd.row(idx_eta) = Psi;
    }
    
    *P_out   = P_loc;
    *Mu_out  = Mu_loc;
    *F_out   = F_loc;
    *S_out   = S_loc;
    *Psi_out = Psi_grd;
}

template class OdfOptimiser<MESOGEN, float>;
template class OdfOptimiser<MESOGEN, double>;
