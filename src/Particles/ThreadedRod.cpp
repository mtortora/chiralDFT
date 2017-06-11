// ===================================================================
/**
 * Threaded rod derived particle class.
 * Hard spherocylinder potentially decorated with helical charges.
 */
// ===================================================================
/*
 * ThreadedRod.cpp: Version 1.0
 * Created 12/05/2016 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/ThreadedRod.hpp"


template<typename number>
ThreadedRod<number>::ThreadedRod()
{
    // Bounding leaf parameter
    this->BVH.SetLeafParameter(3);

    this->N_DELTA_L = 2;

    // Rod parameters - if the Lagerwall Debye-Huckel potential is used, model lengths are reported in nanometers
    N_PATCH_        = 17;
    N_RES_          = 15;
    
    D_HARD_         = 1. * this->SIGMA_R;

    L_Z_            = 12. * D_HARD_;
    P_PATCH_        = 10. * D_HARD_;
    R_PATCH_        = 0.5 * D_HARD_;
    
    this->V0        = PI/4. * SQR(D_HARD_) * L_Z_ + PI/6. * CUB(D_HARD_);
    this->V_EFF     = this->V0;
    
    // Debye-Huckel parameters
    if ( USE_DH )
    {
        if ( MODE_DH == DH_WENSINK )
        {
            number c_dens = L_Z_ / (number)N_PATCH_;
            
            DH_PREFACTOR_ = 1. * SQR(c_dens);
            MINUS_KAPPA_  = -20. / L_Z_;
        }
        
        else if ( MODE_DH == DH_LAGERWALL )
        {
            number z_s    = 15. * (number)N_PATCH_;
            number eta_c  = 0.19;
            number l_bjer = 0.07 * D_HARD_;
            
            number rho_c  = eta_c / this->V0;
            number c_dens = z_s / (number)N_PATCH_;
            
            DH_PREFACTOR_ = l_bjer * SQR(c_dens);
            MINUS_KAPPA_  = -sqrt(4*PI * l_bjer * (z_s*rho_c + 2*C_SALT));
        }
        
        else throw std::runtime_error("Unsupported electrostatics model for threaded rods");
        
        // Needs to fully include electrostatics + hard cylindrical backbone
        R_CUT_ = fmax(5. / -MINUS_KAPPA_, D_HARD_+2.*R_PATCH_);
        E_CUT_ = 20.;
    }
    
    else R_CUT_ = 0.;
    
    this->R_INTEG = L_Z_ + R_CUT_ + D_HARD_;
    this->V_INTEG = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void ThreadedRod<number>::Build(int mpi_rank)
{
    Matrix3X<number> Patches  (3, N_PATCH_);
    Matrix3X<number> Wireframe(3, SQR(N_RES_));
    
    ArrayX<number> Z_grid     = ArrayX<number>::LinSpaced(N_PATCH_, 0., L_Z_);
    ArrayX<number> Phi_grid   = ArrayX<number>::LinSpaced(N_RES_, 0., 2.*PI);
    ArrayX<number> Theta_grid = ArrayX<number>::LinSpaced(N_RES_, 0., PI);

    std::string DATA_PATH;
    
    if ( mpi_rank == MPI_MASTER ) DATA_PATH = __DATA_PATH;
    else                          DATA_PATH = "/dev/null";
    
    std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");
    
    Patches.row(0) = R_PATCH_ * cos(2.*PI/P_PATCH_ * Z_grid);
    Patches.row(1) = R_PATCH_ * sin(2.*PI/P_PATCH_ * Z_grid);
    Patches.row(2) = Z_grid;
    
    Vector3<number> Center;
    Vector3<number> Center_of_mass = Patches.rowwise().mean();
    
    // Draw spherical beads in lieu of interaction patches
    for ( uint idx_center = 0; idx_center < N_PATCH_; ++idx_center )
    {
        Patches.col(idx_center) -= Center_of_mass;
        Center                   = Patches.col(idx_center);
        
        for ( uint idx_theta = 0; idx_theta < N_RES_; ++idx_theta )
        {
            number theta_ = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                number phi_        = Phi_grid(idx_phi);
                uint   idx         = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) = Center;
                
                Wireframe(0, idx) += R_CUT_/2. * sin(theta_)*cos(phi_);
                Wireframe(1, idx) += R_CUT_/2. * sin(theta_)*sin(phi_);
                Wireframe(2, idx) += R_CUT_/2. * cos(theta_);
                
                file_wireframe    << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    // Draw cylinder
    Wireframe.resize(3, N_RES_);

    ArrayX<number> Ax_grid = ArrayX<number>::LinSpaced(N_RES_, -L_Z_/2., L_Z_/2.);
    
    for ( uint idx_center = 0; idx_center < N_RES_; ++idx_center )
    {
        Center << 0., 0., Ax_grid(idx_center);
        
        for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
        {
            number phi_            = Phi_grid(idx_phi);
            
            Wireframe.col(idx_phi) = Center;
            
            Wireframe(0, idx_phi) += D_HARD_/2. * cos(phi_);
            Wireframe(1, idx_phi) += D_HARD_/2. * sin(phi_);
            
            file_wireframe        << Wireframe.col(idx_phi).adjoint() << std::endl;
        }

        file_wireframe << std::endl;
    }
    
    // Draw spherical caps
    Wireframe.resize(3, SQR(N_RES_)/2);

    for ( uint idx_cap = 0; idx_cap < 2; ++idx_cap )
    {
        Theta_grid = ArrayX<number>::LinSpaced(N_RES_/2, (idx_cap+1) * PI/2., idx_cap * PI/2.);
        
        for ( uint idx_theta = 0; idx_theta < N_RES_/2; ++idx_theta )
        {
            number theta_ = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                number phi_         = Phi_grid(idx_phi);
                uint   idx          = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) << 0., 0., L_Z_/2. - idx_cap*L_Z_;
                
                Wireframe(0, idx)  += D_HARD_/2. * sin(theta_)*cos(phi_);
                Wireframe(1, idx)  += D_HARD_/2. * sin(theta_)*sin(phi_);
                Wireframe(2, idx)  += D_HARD_/2. * cos(theta_);
                
                file_wireframe     << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    file_wireframe.close();
    
    if ( USE_DH )
    {
        // Build bounding volume hierarchy
        this->BVH.Build(Patches, R_CUT_);
        
        if ( this->id_ == 1 ) this->BVH.PrintBuildInfo();
    }
    
    else
    {
        this->Hull       = &this->BVH;
        
        this->Hull->l_xh =  D_HARD_ / 2.;
        this->Hull->l_yh =  D_HARD_ / 2.;
        this->Hull->l_zh = (L_Z_+D_HARD_) / 2.;
        
        this->Hull->l_ch = L_Z_    / 2.;
        this->Hull->l_cr = D_HARD_ / 2.;
    }
}

template class ThreadedRod<float>;
template class ThreadedRod<double>;
