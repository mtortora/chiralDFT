// ===================================================================
/**
 * Patchy rod derived particle class
 */
// ===================================================================
/*
 * PatchyRod.cpp: Version 1.0
 * Created 25/11/2016 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/PatchyRod.hpp"

using namespace Eigen;


PatchyRod::PatchyRod()
{
    N_DELTA_L    = 2;
    
    // Rod geometric parameters
    N_PATCH_     = 21;
    N_BACK_      = 21;
    N_RES_       = 15;
    
    L_Z_         = 10. * SIGMA_R;
    P_PATCH_     = 5.  * L_Z_;
    
    V0           = PI/4. * SQR(SIGMA_R) * L_Z_ + PI/6. * CUB(SIGMA_R);
    V_EFF        = V0;
    
    // WCA parameters
    R_WCA_       = pow(2., 1./6) * SIGMA_R;
    EPSILON_WCA_ = 1.;
    
    E_CUT_       = 20.;
    
    // Debye-Huckel parameters
    if ( USE_DH )
    {
        DH_PREFACTOR_    = 250./ (double)SQR(N_PATCH_);
        MINUS_KAPPA_     = -1. / (SIGMA_R);
        R_CUT_           = 5.  * SIGMA_R;
        
        // Bounding volume parameters
        BHierarchy->l_xh = (SIGMA_R + R_CUT_) / 2.;
        BHierarchy->l_yh = (SIGMA_R + R_CUT_) / 2.;
        BHierarchy->l_zh = (L_Z_    + R_CUT_) / 2.;
        
        BHierarchy->l_ch = (SIGMA_R + L_Z_)   / 2.;
        BHierarchy->l_rc = (SIGMA_R + R_CUT_) / 2.;
    }
    
    else
    {
        R_CUT_           =  R_WCA_;
        
        // Bounding volume parameters
        BHierarchy->l_xh = (R_WCA_ + SIGMA_R) / 2.;
        BHierarchy->l_yh = (R_WCA_ + SIGMA_R) / 2.;
        BHierarchy->l_zh = (R_WCA_ + L_Z_)    / 2.;
        
        BHierarchy->l_ch =  L_Z_              / 2.;
        BHierarchy->l_rc = (R_WCA_ + SIGMA_R) / 2.;
    }
    
    R_INTEG = 2. * (BHierarchy->l_ch + BHierarchy->l_rc);
    V_INTEG = CUB(2.*R_INTEG) * 16.*pow(PI, 6);
}

// ============================
/* Build particle model */
// ============================
void PatchyRod::Build(int mpi_rank)
{
    double    r_patch = USE_DH ? -1./(2.*MINUS_KAPPA_) : SIGMA_R/2.;
    
    Matrix3Xd Patches_ (3, N_PATCH_);
    Matrix3Xd Backbone_(3, N_BACK_);
    Matrix3Xd Wireframe(3, SQR(N_RES_));
    
    ArrayXd   Z_grid     = ArrayXd::LinSpaced(N_PATCH_, -L_Z_/2., L_Z_/2.);
    ArrayXd   Phi_grid   = ArrayXd::LinSpaced(N_RES_, 0., 2.*PI);
    ArrayXd   Theta_grid = ArrayXd::LinSpaced(N_RES_, 0., PI);
    
    std::string DATA_PATH;
    
    if ( mpi_rank == MPI_MASTER ) DATA_PATH = __DATA_PATH;
    else                          DATA_PATH = "/dev/null";
    
    std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");
    
    Backbone_.row(0) = ArrayXd::Zero(N_BACK_);
    Backbone_.row(1) = ArrayXd::Zero(N_BACK_);
    Backbone_.row(2) = ArrayXd::LinSpaced(N_BACK_, -L_Z_/2., L_Z_/2.);

    Patches_ .row(0) = SIGMA_R/2. * cos(2.*PI/P_PATCH_ * Z_grid);
    Patches_ .row(1) = SIGMA_R/2. * sin(2.*PI/P_PATCH_ * Z_grid);
    Patches_ .row(2) = Z_grid;
    
    // Draw spherical beads in lieu of interaction sites
    for ( uint idx_center = 0; idx_center < N_BACK_; ++idx_center )
    {
        Vector3d Center = Backbone_.col(idx_center);
        
        for ( uint idx_theta = 0; idx_theta < N_RES_; ++idx_theta )
        {
            double theta = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                double phi         = Phi_grid(idx_phi);
                uint   idx         = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) = Center;
                
                Wireframe(0, idx) += SIGMA_R/2. * sin(theta)*cos(phi);
                Wireframe(1, idx) += SIGMA_R/2. * sin(theta)*sin(phi);
                Wireframe(2, idx) += SIGMA_R/2. * cos(theta);
                
                file_wireframe    << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    for ( uint idx_center = 0; idx_center < N_PATCH_; ++idx_center )
    {
        Vector3d Center = Patches_.col(idx_center);
        
        for ( uint idx_theta = 0; idx_theta < N_RES_; ++idx_theta )
        {
            double theta = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                double phi         = Phi_grid(idx_phi);
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
