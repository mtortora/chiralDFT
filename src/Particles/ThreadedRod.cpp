// ===================================================================
/**
 * Threaded rod derived particle class
 */
// ===================================================================
/*
 * ThreadedRod.cpp: Version 1.0
 * Created 12/05/2016 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/ThreadedRod.hpp"

using namespace Eigen;


ThreadedRod::ThreadedRod()
{
    // Bounding tree properties
    BHierarchy->SetTreeProperties(0);

    N_DELTA_L = 2;

    // Rod parameters - if the Lagerwall Debye-Huckel potential is used, model lengths are reported in nanometers
    N_PATCH_  = 17;
    N_RES_    = 15;
    
    D_HARD_   = 10. * SIGMA_R;

    L_Z_      = 20. * D_HARD_;
    P_PATCH_  = 20. * D_HARD_;
    
    V0        = PI/4. * SQR(D_HARD_) * L_Z_ + PI/6. * CUB(D_HARD_);
    V_EFF     = V0;

    R_CUT_    = 0.;
    
    // Debye-Huckel parameters
    if ( USE_DH )
    {
        if ( MODE_DH == DH_WENSINK )
        {
            double c_dens = L_Z_ / (double)N_PATCH_;
            
            DH_PREFACTOR_ = 1. * SQR(c_dens);
            MINUS_KAPPA_  = -20. / L_Z_;
        }
        
        else if ( MODE_DH == DH_LAGERWALL )
        {
            double z_s    = 1000.;
            double eta_c  = 0.05;
            double l_bjer = 0.7;
            
            double rho_c  = eta_c / V0;
            double c_dens = z_s / (double)N_PATCH_;
            
            DH_PREFACTOR_ = l_bjer * SQR(c_dens);
            MINUS_KAPPA_  = -sqrt(4*PI * l_bjer * (z_s*rho_c + 2*C_SALT));
        }
        
        else throw std::runtime_error("Unsupported electrostatics model for threaded rods");
        
        R_CUT_ = 10. / -MINUS_KAPPA_;
        E_CUT_ = 20.;
    }
    
    R_INTEG = L_Z_ + R_CUT_ + D_HARD_;
    V_INTEG = CUB(2.*R_INTEG) * 16.*pow(PI, 6);
}

// ============================
/* Build particle model */
// ============================
void ThreadedRod::Build(int mpi_rank)
{
    Matrix3Xd Patches  (3, N_PATCH_);
    Matrix3Xd Wireframe(3, SQR(N_RES_));
    
    ArrayXd   Z_grid     = ArrayXd::LinSpaced(N_PATCH_, 0., L_Z_);
    ArrayXd   Phi_grid   = ArrayXd::LinSpaced(N_RES_, 0., 2.*PI);
    ArrayXd   Theta_grid = ArrayXd::LinSpaced(N_RES_, 0., PI);

    std::string DATA_PATH;
    
    if ( mpi_rank == MPI_MASTER ) DATA_PATH = __DATA_PATH;
    else                          DATA_PATH = "/dev/null";
    
    std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");
    
    Patches.row(0) = D_HARD_/2. * cos(2.*PI/P_PATCH_ * Z_grid);
    Patches.row(1) = D_HARD_/2. * sin(2.*PI/P_PATCH_ * Z_grid);
    Patches.row(2) = Z_grid;
    
    Vector3d Center;
    Vector3d Center_of_mass = Patches.rowwise().mean();
    
    // Draw spherical beads in lieu of interaction patches
    for ( uint idx_center = 0; idx_center < N_PATCH_; ++idx_center )
    {
        Patches.col(idx_center) -= Center_of_mass;
        Center                     = Patches.col(idx_center);
        
        for ( uint idx_theta = 0; idx_theta < N_RES_; ++idx_theta )
        {
            double theta = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                double phi = Phi_grid(idx_phi);
                uint   idx = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) = Center;
                
                Wireframe(0, idx) += R_CUT_/2. * sin(theta)*cos(phi);
                Wireframe(1, idx) += R_CUT_/2. * sin(theta)*sin(phi);
                Wireframe(2, idx) += R_CUT_/2. * cos(theta);
                
                file_wireframe << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    // Draw cylinder
    Wireframe.resize(3, N_RES_);

    ArrayXd Ax_grid = ArrayXd::LinSpaced(N_RES_, -L_Z_/2., L_Z_/2.);
    
    for ( uint idx_center = 0; idx_center < N_RES_; ++idx_center )
    {
        Center << 0., 0., Ax_grid(idx_center);
        
        for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
        {
            double phi = Phi_grid(idx_phi);
            
            Wireframe.col(idx_phi) = Center;
            
            Wireframe(0, idx_phi) += D_HARD_/2. * cos(phi);
            Wireframe(1, idx_phi) += D_HARD_/2. * sin(phi);
            
            file_wireframe << Wireframe.col(idx_phi).adjoint() << std::endl;
        }

        file_wireframe << std::endl;
    }
    
    // Draw spherical caps
    Wireframe.resize(3, SQR(N_RES_)/2);

    for ( uint idx_cap = 0; idx_cap < 2; ++idx_cap )
    {
        Theta_grid = ArrayXd::LinSpaced(N_RES_/2, (idx_cap+1) * PI/2., idx_cap * PI/2.);
        
        for ( uint idx_theta = 0; idx_theta < N_RES_/2; ++idx_theta )
        {
            double theta = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                double phi = Phi_grid(idx_phi);
                uint   idx = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) << 0., 0., L_Z_/2. - idx_cap*L_Z_;
                
                Wireframe(0, idx) += D_HARD_/2. * sin(theta)*cos(phi);
                Wireframe(1, idx) += D_HARD_/2. * sin(theta)*sin(phi);
                Wireframe(2, idx) += D_HARD_/2. * cos(theta);
                
                file_wireframe << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    file_wireframe.close();
    
    if ( USE_DH )
    {
        // Build bounding volume hierarchy
        BHierarchy->AllocateForest(1);
        BTree* Tree = &BHierarchy->Forest[0];
        
        BHierarchy->RecursiveBuild(Tree, Patches, R_CUT_);
        
        if ( id_ == 1 ) BHierarchy->PrintBuildInfo();
    }
    
    else
    {
        BHierarchy->l_xh = D_HARD_ / 2.;
        BHierarchy->l_yh = D_HARD_ / 2.;
        BHierarchy->l_zh = L_Z_    / 2.;
        
        BHierarchy->l_ch = L_Z_    / 2.;
        BHierarchy->l_rc = D_HARD_ / 2.;
    }
}
