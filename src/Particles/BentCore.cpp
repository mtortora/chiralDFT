// ===================================================================
/**
 * Bent-core mesogen derived particle class.
 */
// ===================================================================
/*
 * BentCore.cpp: Version 2.1
 * Created 20/09/2015 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/BentCore.hpp"


template<typename number>
BentCore<number>::BentCore()
{
    // Bounding leaf parameter
    this->BVH.SetLeafParameter(3);
    
    this->N_DELTA_L = 2;

    // Bent-core parameters
    N_S_            = 11;
    N_RES_          = 20;
    
    D_HARD_         = 1.   * this->SIGMA_R;
    L_CTR_          = 10.  * this->SIGMA_R;
    CHI_            = 130. * PI/180.;
    
    R_BTC_          = L_CTR_ / (PI - CHI_);
    GAMMA_          = (PI - CHI_) / 2.;
    
    number d_cc     = 2.*R_BTC_ * sin(GAMMA_ / (N_S_-1.));
    
    L_X_            = R_BTC_ * (1. - cos(GAMMA_) + std::abs(1. + cos(GAMMA_) - 2.*sin(GAMMA_)/GAMMA_));
    L_Y_            = 0.;
    L_Z_            = 2.*R_BTC_ * sin(GAMMA_);
    
    this->R_INTEG   = 2.*R_BTC_ * sqrt(SQR(sin(GAMMA_)/GAMMA_) - sin(2.*GAMMA_)/GAMMA_ + 1.) + D_HARD_;
    this->V_INTEG   = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
    
    this->V0        = PI/6.*CUB(D_HARD_) * (1. + (N_S_-1.)/2. * (3.*d_cc/D_HARD_ - CUB(d_cc/D_HARD_)));
    this->V_EFF     = PI/6.*CUB(D_HARD_) * (1. + (N_S_-1.) * (3.*d_cc/D_HARD_ - CUB(d_cc/D_HARD_)/2. - 3.*sqrt(1.-SQR(d_cc/(2.*D_HARD_))) * asin(d_cc/(2.*D_HARD_))));
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void BentCore<number>::Build(int mpi_rank)
{
    Matrix3X<number> Backbone_(3, N_S_);
    Matrix3X<number> Wireframe(3, SQR(N_RES_));
    
    ArrayX<number>   Gamma_grid = ArrayX<number>::LinSpaced(N_S_, -GAMMA_, GAMMA_);
    ArrayX<number>   Theta_grid = ArrayX<number>::LinSpaced(N_RES_, 0., PI);
    ArrayX<number>   Phi_grid   = ArrayX<number>::LinSpaced(N_RES_, 0., 2.*PI);
 
    std::string DATA_PATH;
    
    // Redirect slave output to /dev/null
    if ( mpi_rank == MPI_MASTER ) DATA_PATH = __DATA_PATH;
    else                          DATA_PATH = "/dev/null";
    
    std::ofstream file_backbone (DATA_PATH + "/backbone.out");
    std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");

    Backbone_.row(0) = R_BTC_ * (1. - cos(Gamma_grid));
    Backbone_.row(1) = ArrayX<number>::Zero(N_S_);
    Backbone_.row(2) = R_BTC_ * sin(Gamma_grid);
    
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
            number theta = Theta_grid(idx_theta);
            
            for ( uint idx_phi = 0; idx_phi < N_RES_; ++idx_phi )
            {
                number phi         = Phi_grid(idx_phi);
                uint   idx         = idx_theta*N_RES_ + idx_phi;
                
                Wireframe.col(idx) = Center;
                
                Wireframe(0, idx) += D_HARD_/2. * sin(theta)*cos(phi);
                Wireframe(1, idx) += D_HARD_/2. * sin(theta)*sin(phi);
                Wireframe(2, idx) += D_HARD_/2. * cos(theta);
                
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

template class BentCore<float>;
template class BentCore<double>;
