// ===================================================================
/**
 * FlexibleHelix derived particle class.
 */
// ===================================================================
/*
 * FlexibleHelix.cpp: Version 1.0
 * Created 19/07/2017 by Maxime Tortora
 */
// ===================================================================

#include "Utils.hpp"
#include "Particles/FlexibleHelix.hpp"


template<typename number>
FlexibleHelix<number>::FlexibleHelix()
{
	// Bounding leaf parameter
	this->BVH.SetLeafParameter(3);
	
	this->N_DELTA_L = 2;
	
	// Helix parameters
	N_S_            = 15;
	
	D_HARD_         = 1.  * this->SIGMA_R;
	L_CTR_          = 10. * this->SIGMA_R;
	
	R_HLX_          = 0.4 * this->SIGMA_R;
	P_HLX_          = 8.  * this->SIGMA_R;
	
	L_Z_            = L_CTR_ / sqrt(1. + SQR(2.*PI * R_HLX_/P_HLX_));
	
	number d_cc     = sqrt(SQR(2.*R_HLX_*sin(PI/P_HLX_ * L_Z_/(N_S_-1.))) + SQR(L_Z_/(N_S_-1.)));
	
	this->V0        = PI/6.*CUB(D_HARD_) * (1. + (N_S_-1.)/2. * (3.*d_cc/D_HARD_ - CUB(d_cc/D_HARD_)));
	this->V_EFF     = PI/6.*CUB(D_HARD_) * (1. + (N_S_-1.) * (3.*d_cc/D_HARD_ - CUB(d_cc/D_HARD_)/2. - 3.*sqrt(1.-SQR(d_cc/(2.*D_HARD_))) * asin(d_cc/(2.*D_HARD_))));
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void FlexibleHelix<number>::Build(int mpi_rank)
{
	uint N_CONF;
	uint N_TOT;

	ArrayX<uint>     Sizes;
	Matrix3X<number> Backbones;
	
	// Load configurations from trajectory files on master thread
	if ( mpi_rank == MPI_MASTER )
	{
		std::string DATA_PATH   = __DATA_PATH;
		std::string filename_in = DATA_PATH + "/trajectory.in";
		
		Utils<number>::Load(filename_in, &Backbones, &Sizes);
		
		if ( Backbones.size() == 0 ) throw std::runtime_error("Unreadable input trajectory file");
		
		N_CONF = Sizes.size();
		N_TOT  = Sizes.sum();
		
		if ( N_TOT != N_CONF*N_S_ ) throw std::runtime_error("Trajectory file inconsistent with input parameters");
	}
	
	// Broadcast data to slave threads
	MPI_Bcast(&N_CONF, 1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&N_TOT,  1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	
	if ( mpi_rank != MPI_MASTER )
	{
		Sizes.resize(N_CONF);
		Backbones.resize(3, N_TOT);
	}
	
	MPI_Bcast(Sizes.data(),     Sizes.size(),     Utils<uint>()  .MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	MPI_Bcast(Backbones.data(), Backbones.size(), Utils<number>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	
	// Set center of masses to the origin and main axes to e_z
	for ( uint idx_conf = 0; idx_conf < N_CONF; ++idx_conf )
	{
		Matrix3X<number> Backbone      = Matrix3X<number>::Map(Backbones.data() + 3 * idx_conf*N_S_, 3, N_S_);
		
		Vector3<number> Center_of_mass = Backbone.rowwise().mean();
		Backbone                       = Backbone.colwise() - Center_of_mass;
		
		Matrix33<number> Rot           = Utils<number>::PCA(Backbone);
		Backbone                       = Rot.transpose() * Backbone;
		
		Backbones.block(0, N_S_*idx_conf, 3, N_S_) = Backbone;
	}
	
	// Build bounding volume hierarchy
	this->BVH.Build(Backbones, D_HARD_, Sizes);
	
	// Print simulation parameters
	if ( this->id_ == 1 )
	{
		LogTxt("Loaded particle trajectory file: %d configurations, %d interaction sites", N_CONF, N_TOT);
		
		this->BVH.PrintBuildInfo();
	}
	
	N_CONF_ = N_CONF;
	
	this->R_INTEG = 2*Backbones.colwise().norm().maxCoeff() + D_HARD_;
	this->V_INTEG = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
}

template class FlexibleHelix<float>;
template class FlexibleHelix<double>;
