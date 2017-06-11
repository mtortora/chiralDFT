// ===================================================================
/**
 * Helix derived particle class.
 */
// ===================================================================
/*
 * Helix.cpp: Version 2.1
 * Created 20/09/2015 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/Helix.hpp"


template<typename number>
Helix<number>::Helix()
{
    // Bounding leaf parameter
    this->BVH.SetLeafParameter(3);

    this->N_DELTA_L = 2;

    // Helix parameters
    N_S_            = 15;
    N_RES_          = 20;
    
    D_HARD_         = 1.   * this->SIGMA_R;
    L_CTR_          = 10.  * this->SIGMA_R;
    
    R_HLX_          = 0.4  * this->SIGMA_R;
    P_HLX_          = 8. * this->SIGMA_R;
    
    L_Z_            = L_CTR_ / sqrt(1. + SQR(2.*PI * R_HLX_/P_HLX_));
    
    number d_cc     = sqrt(SQR(2.*R_HLX_*sin(PI/P_HLX_ * L_Z_/(N_S_-1.))) + SQR(L_Z_/(N_S_-1.)));
    
    this->R_INTEG   = sqrt(SQR(L_Z_+D_HARD_) + SQR(2.*R_HLX_+D_HARD_));
    this->V_INTEG   = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
    
    this->V0        = PI/6.*CUB(D_HARD_) * (1. + (N_S_-1.)/2. * (3.*d_cc/D_HARD_ - CUB(d_cc/D_HARD_)));
    this->V_EFF     = PI/6.*CUB(D_HARD_) * (1. + (N_S_-1.) * (3.*d_cc/D_HARD_ - CUB(d_cc/D_HARD_)/2. - 3.*sqrt(1.-SQR(d_cc/(2.*D_HARD_))) * asin(d_cc/(2.*D_HARD_))));
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void Helix<number>::Build(int mpi_rank)
{
    Matrix3X<number> Backbone_(3, N_S_);
    Matrix3X<number> Wireframe(3, SQR(N_RES_));
    
    ArrayX<number> Z_grid     = ArrayX<number>::LinSpaced(N_S_, 0., L_Z_);
    ArrayX<number> Theta_grid = ArrayX<number>::LinSpaced(N_RES_, 0., PI);
    ArrayX<number> Phi_grid   = ArrayX<number>::LinSpaced(N_RES_, 0., 2.*PI);

    std::string DATA_PATH;
    
    // Redirect slave output to /dev/null
    if ( mpi_rank == MPI_MASTER ) DATA_PATH = __DATA_PATH;
    else                          DATA_PATH = "/dev/null";
    
    std::ofstream file_backbone (DATA_PATH + "/backbone.out");
    std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");

    Backbone_.row(0) = R_HLX_ * cos(2.*PI/P_HLX_ * Z_grid);
    Backbone_.row(1) = R_HLX_ * sin(2.*PI/P_HLX_ * Z_grid);
    Backbone_.row(2) = Z_grid;
    
    // Set center of mass to the origin and main axis to e_z
    Vector3<number> Center_of_mass = Backbone_.rowwise().mean();
    Backbone_            = Backbone_.colwise() - Center_of_mass;
    
    Matrix33<number> Rot = Utils<number>::PCA(Backbone_);
    Backbone_            = Rot.transpose() * Backbone_;
    
	for ( uint idx_center = 0; idx_center < N_S_; ++idx_center )
    {
        Vector3<number> Center = Backbone_.col(idx_center);
        
        // Draw spherical beads
        for ( uint idx_theta = 0; idx_theta < N_RES_; ++idx_theta )
        {
            number theta_ = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                number phi_        = Phi_grid(idx_phi);
                uint   idx         = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) = Center;
                
                Wireframe(0, idx) += D_HARD_/2. * sin(theta_)*cos(phi_);
                Wireframe(1, idx) += D_HARD_/2. * sin(theta_)*sin(phi_);
                Wireframe(2, idx) += D_HARD_/2. * cos(theta_);
                
                file_wireframe    << Wireframe.col(idx).adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    file_backbone << Backbone_.transpose();
    
    file_backbone .close();
    file_wireframe.close();
    
    Backbone = Backbone_;

    // Build bounding volume hierarchy
    this->BVH.Build(Backbone, D_HARD_);
    
    if ( this->id_ == 1 ) this->BVH.PrintBuildInfo();
}

template class Helix<float>;
template class Helix<double>;
