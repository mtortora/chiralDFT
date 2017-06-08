// ===================================================================
/**
 * Patchy rod derived particle class.
 * Linear Yukawa spheres potentially decorated with helical charges.
 */
// ===================================================================
/*
 * PatchyRod.cpp: Version 1.0
 * Created 25/11/2016 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/PatchyRod.hpp"


template<typename number>
PatchyRod<number>::PatchyRod()
{
    this->N_DELTA_L = 2;
    this->Hull      = &this->BVH;
    
    // Rod geometric parameters
    N_PATCH_        = 21;
    N_BACK_         = 21;
    N_RES_          = 15;
    
    L_Z_            = 10. * this->SIGMA_R;
    P_PATCH_        = 5.  * L_Z_;
    
    this->V0        = PI/4. * SQR(this->SIGMA_R) * L_Z_ + PI/6. * CUB(this->SIGMA_R);
    this->V_EFF     = this->V0;
    
    // WCA parameters
    R_WCA_          = pow(2., 1./6) * this->SIGMA_R;
    EPSILON_WCA_    = 1.;
    
    E_CUT_          = 20.;
    
    // Debye-Huckel parameters
    if ( USE_DH )
    {
        DH_PREFACTOR_    = 250./ (number)SQR(N_PATCH_);
        MINUS_KAPPA_     = -1. / (this->SIGMA_R);
        R_CUT_           = 5.  * this->SIGMA_R;
        
        // Bounding volume parameters
        this->Hull->l_xh = (this->SIGMA_R + R_CUT_) / 2.;
        this->Hull->l_yh = (this->SIGMA_R + R_CUT_) / 2.;
        this->Hull->l_zh = (L_Z_    + R_CUT_) / 2.;
        
        this->Hull->l_ch = (this->SIGMA_R + L_Z_)   / 2.;
        this->Hull->l_cr = (this->SIGMA_R + R_CUT_) / 2.;
    }
    
    else
    {
        R_CUT_           =  R_WCA_;
        
        // Bounding volume parameters
        this->Hull->l_xh = (R_WCA_ + this->SIGMA_R) / 2.;
        this->Hull->l_yh = (R_WCA_ + this->SIGMA_R) / 2.;
        this->Hull->l_zh = (R_WCA_ + L_Z_)    / 2.;
        
        this->Hull->l_ch =  L_Z_              / 2.;
        this->Hull->l_cr = (R_WCA_ + this->SIGMA_R) / 2.;
    }
    
    this->R_INTEG = 2. * (this->Hull->l_ch + this->Hull->l_cr);
    this->V_INTEG = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void PatchyRod<number>::Build(int mpi_rank)
{
    number r_patch = USE_DH ? -1./(2.*MINUS_KAPPA_) : this->SIGMA_R/2.;
    
    Matrix3X<number> Patches_ (3, N_PATCH_);
    Matrix3X<number> Backbone_(3, N_BACK_);
    Matrix3X<number> Wireframe(3, SQR(N_RES_));
    
    ArrayX<number> Z_grid     = ArrayX<number>::LinSpaced(N_PATCH_, -L_Z_/2., L_Z_/2.);
    ArrayX<number> Phi_grid   = ArrayX<number>::LinSpaced(N_RES_, 0., 2.*PI);
    ArrayX<number> Theta_grid = ArrayX<number>::LinSpaced(N_RES_, 0., PI);
    
    std::string DATA_PATH;
    
    if ( mpi_rank == MPI_MASTER ) DATA_PATH = __DATA_PATH;
    else                          DATA_PATH = "/dev/null";
    
    std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");
    
    Backbone_.row(0) = ArrayX<number>::Zero(N_BACK_);
    Backbone_.row(1) = ArrayX<number>::Zero(N_BACK_);
    Backbone_.row(2) = ArrayX<number>::LinSpaced(N_BACK_, -L_Z_/2., L_Z_/2.);

    Patches_ .row(0) = this->SIGMA_R/2. * cos(2.*PI/P_PATCH_ * Z_grid);
    Patches_ .row(1) = this->SIGMA_R/2. * sin(2.*PI/P_PATCH_ * Z_grid);
    Patches_ .row(2) = Z_grid;
    
    // Draw spherical beads in lieu of interaction sites
    for ( uint idx_center = 0; idx_center < N_BACK_; ++idx_center )
    {
        Vector3<number> Center = Backbone_.col(idx_center);
        
        for ( uint idx_theta = 0; idx_theta < N_RES_; ++idx_theta )
        {
            number theta = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                number phi         = Phi_grid(idx_phi);
                uint   idx         = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) = Center;
                
                Wireframe(0, idx) += this->SIGMA_R/2. * sin(theta)*cos(phi);
                Wireframe(1, idx) += this->SIGMA_R/2. * sin(theta)*sin(phi);
                Wireframe(2, idx) += this->SIGMA_R/2. * cos(theta);
                
                file_wireframe    << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    for ( uint idx_center = 0; idx_center < N_PATCH_; ++idx_center )
    {
        Vector3<number> Center = Patches_.col(idx_center);
        
        for ( uint idx_theta = 0; idx_theta < N_RES_; ++idx_theta )
        {
            number theta = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                number phi         = Phi_grid(idx_phi);
                uint   idx         = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) = Center;
                
                Wireframe(0, idx) += r_patch * sin(theta)*cos(phi);
                Wireframe(1, idx) += r_patch * sin(theta)*sin(phi);
                Wireframe(2, idx) += r_patch * cos(theta);
                
                file_wireframe    << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    Patches  = Patches_;
    Backbone = Backbone_;
}

template class PatchyRod<float>;
template class PatchyRod<double>;
