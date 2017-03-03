// ===================================================================
/**
 * Helix derived particle class
 */
// ===================================================================
/*
 * Helix.cpp: Version 2.1
 * Created 20/09/2015 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/Helix.hpp"

using namespace Eigen;


Helix::Helix()
{
    // Bounding tree properties
    BHierarchy->SetTreeProperties(8);

    N_DELTA_L   = 2;

    // Helix parameters
    N_S_        = 500;
    N_RES_      = 20;
    
    R_HARD_     = 3.3   * SIGMA_R;
    L_CTR_      = 880.  * SIGMA_R;
    
    R_HLX_      = 146.  * SIGMA_R;
    P_HLX_      = 2800. * SIGMA_R;
    
    L_X_        = 2.*R_HLX_;
    L_Y_        = 2.*R_HLX_;
    L_Z_        = L_CTR_ / sqrt(1. + SQR(2.*PI * R_HLX_/P_HLX_));
    
    double d_cc = sqrt(SQR(2.*R_HLX_*sin(PI/P_HLX_ * L_Z_/(N_S_-1.))) + SQR(L_Z_/(N_S_-1.)));
    
    R_INTEG     = sqrt(SQR(L_Z_+R_HARD_) + SQR(2.*R_HLX_+R_HARD_));
    V_INTEG     = CUB(2.*R_INTEG) * 16.*pow(PI, 6);
    
    V0          = PI/6.*CUB(R_HARD_) * (1. + (N_S_-1.)/2. * (3.*d_cc/R_HARD_ - CUB(d_cc/R_HARD_)));
    V_EFF       = PI/6.*CUB(R_HARD_) * (1. + (N_S_-1.) * (3.*d_cc/R_HARD_ - CUB(d_cc/R_HARD_)/2. - 3.*sqrt(1.-SQR(d_cc/(2.*R_HARD_))) * asin(d_cc/(2.*R_HARD_))));
}

// ============================
/* Build particle model */
// ============================
void Helix::Build(int mpi_rank)
{
    Matrix3Xd Backbone_(3, N_S_);
    Matrix3Xd Wireframe(3, SQR(N_RES_));
    
    ArrayXd   Z_grid     = ArrayXd::LinSpaced(N_S_, 0., L_Z_);
    ArrayXd   Theta_grid = ArrayXd::LinSpaced(N_RES_, 0., PI);
    ArrayXd   Phi_grid   = ArrayXd::LinSpaced(N_RES_, 0., 2.*PI);

    std::string DATA_PATH;
    
    // Redirect slave output to /dev/null
    if ( mpi_rank == MPI_MASTER ) DATA_PATH = __DATA_PATH;
    else                          DATA_PATH = "/dev/null";
    
    std::ofstream file_backbone (DATA_PATH + "/backbone.out");
    std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");

    Backbone_.row(0)        = R_HLX_ * cos(2.*PI/P_HLX_ * Z_grid);
    Backbone_.row(1)        = R_HLX_ * sin(2.*PI/P_HLX_ * Z_grid);
    Backbone_.row(2)        = Z_grid;
    
    // Set center of mass to the origin and main axis to e_z
    Vector3d Center_of_mass = Backbone_.rowwise().mean();
    Backbone_               = Backbone_.colwise() - Center_of_mass;
    
    Matrix3d Rot            = Utils::PCA(Backbone_);
    Backbone_               = Rot.transpose() * Backbone_;
    
	for ( uint idx_center = 0; idx_center < N_S_; ++idx_center )
    {
        Vector3d Center = Backbone_.col(idx_center);
        
        // Draw spherical beads
        for ( uint idx_theta = 0; idx_theta < N_RES_; ++idx_theta )
        {
            double theta = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                double phi = Phi_grid(idx_phi);
                uint   idx = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) = Center;
                
                Wireframe(0, idx) += R_HARD_/2. * sin(theta)*cos(phi);
                Wireframe(1, idx) += R_HARD_/2. * sin(theta)*sin(phi);
                Wireframe(2, idx) += R_HARD_/2. * cos(theta);
                
                file_wireframe << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    file_backbone << Backbone_.transpose();
    
    file_backbone .close();
    file_wireframe.close();
    
    // Build bounding volume hierarchy
    BHierarchy->AllocateForest(1);
    
    BTree* Tree = &BHierarchy->Forest[0];
    BHierarchy->RecursiveBuild(Tree, Backbone_, R_HARD_);
    
    if ( id_ == 1 ) BHierarchy->PrintBuildInfo();
    
    Backbone = Backbone_;
}
