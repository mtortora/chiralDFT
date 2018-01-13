// ===================================================================
/**
 * Flexible patchy rod derived particle class.
 * Particle model from http://dx.doi.org/10.1039/c7sm02077e
 */
// ===================================================================
/*
 * FlexiblePatchyRod.cpp: Version 1.0
 * Created 18/12/2017 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/FlexiblePatchyRod.hpp"


template<typename number>
FlexiblePatchyRod<number>::FlexiblePatchyRod()
{
    this->N_DELTA_L = 2;
		
    // WCA/LJ parameters
	EPSILON_WCA_    = 1.;
	E_CUT_          = 20.;
	
    D_BACK_         = 1. * this->SIGMA_R;
    R_BACK_         = pow(2., 1./6) * this->SIGMA_R;

    D_PATCH_        = 0.1 * D_BACK_;
    R_PATCH_        = 2.5 * this->SIGMA_R;
	
    D_LB_           = (D_BACK_+D_PATCH_)/2.;
    R_LB_           = (R_BACK_+R_PATCH_)/2.;
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void FlexiblePatchyRod<number>::Build(int mpi_rank)
{
    uint N_TOT;
    
    // Load configurations from trajectory files on master thread
    if ( mpi_rank == MPI_MASTER )
    {
        ArrayX<uint> Sizes_bck;
        ArrayX<uint> Sizes_ptc;
        
        std::string DATA_PATH = __DATA_PATH;
        
        std::string file_bck  = DATA_PATH + "/bck.in";
        std::string file_ptc  = DATA_PATH + "/ptc.in";

        Utils<number>::Load(file_bck, &Backbones, &Sizes_bck);
        Utils<number>::Load(file_ptc, &Patches,   &Sizes_ptc);

        if ( (Backbones.size() == 0) || (Patches.size() == 0) ) throw std::runtime_error("Unreadable input file(s)");
        if ( (Sizes_bck != Sizes_ptc).all() ) throw std::runtime_error("Incompatible backbone/patch input files");

        N_BCK   = Sizes_bck(0);
        
        N_TOT   = Sizes_bck.sum();
        N_CONF_ = Sizes_bck.size();

        if ( (Sizes_bck != N_BCK).any() ) throw std::runtime_error("Found configurations of multiple sizes in backbone file");
    }
		
    // Broadcast data to slave threads
    MPI_Bcast(&N_BCK,   1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&N_TOT,   1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&N_CONF_, 1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
    
    if ( mpi_rank != MPI_MASTER )
    {
        Backbones.resize(3, N_TOT);
        Patches  .resize(3, N_TOT);
    }
    
    MPI_Bcast(Backbones.data(), Backbones.size(), Utils<number>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(Patches  .data(), Patches  .size(), Utils<number>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
    
    (this->BVH).Forest.resize(N_CONF_);

    for ( uint idx_conf_ = 0; idx_conf_ < N_CONF_; ++idx_conf_ )
    {
        // Set center of masses to the origin and main axes to e_z
        Matrix3X<number> Backbone      = Matrix3X<number>::Map(Backbones.data() + 3 * idx_conf_*N_BCK, 3, N_BCK);
        Matrix3X<number> Patch         = Matrix3X<number>::Map(Patches  .data() + 3 * idx_conf_*N_BCK, 3, N_BCK);

        Vector3<number> Center_of_mass = Backbone.rowwise().mean();
        
        Backbone                       = Backbone.colwise() - Center_of_mass;
        Patch                          = Patch   .colwise() - Center_of_mass;

        Matrix33<number> Rot           = Utils<number>::PCA(Backbone);
        
        Backbone                       = Rot.transpose() * Backbone;
        Patch                          = Rot.transpose() * Patch;
        
        // Build root bounding volumes
        number z_inf                   = Patch.row(2).minCoeff();
        number z_sup                   = Patch.row(2).maxCoeff();

        number r_max                   = Patch.block(0, 0, 2, N_BCK).colwise().norm().maxCoeff();
        number z_max                   = fmax(std::abs(z_inf),std::abs(z_sup));
                
        this->Hull                     = &(this->BVH).Forest[idx_conf_];
        
        this->Hull->l_xh               = r_max + R_PATCH_/2.;
        this->Hull->l_yh               = r_max + R_PATCH_/2.;
        this->Hull->l_zh               = z_max + R_PATCH_/2.;
        
        this->Hull->l_cr               = r_max + R_PATCH_/2.;
        this->Hull->l_ch               = z_max;
        
        Backbones.block(0, N_BCK*idx_conf_, 3, N_BCK) = Backbone;
        Patches  .block(0, N_BCK*idx_conf_, 3, N_BCK) = Patch;
    }
    
    if ( this->id_ == 1 ) LogTxt("Loaded particle trajectory file: %d configurations, %d interaction sites", N_CONF_, 2*N_TOT);
    
    this->R_INTEG = 2*Patches.colwise().norm().maxCoeff() + R_PATCH_;
    this->V_INTEG = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
    
    this->V0      = 11.76167;
    this->V_EFF   = 11.76167;
}

template class FlexiblePatchyRod<float>;
template class FlexiblePatchyRod<double>;
