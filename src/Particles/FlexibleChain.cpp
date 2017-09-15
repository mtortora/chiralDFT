// ===================================================================
/**
 * FlexibleChain derived particle class.
 */
// ===================================================================
/*
 * FlexibleChain.cpp: Version 1.0
 * Created 19/07/2017 by Maxime Tortora
 */
// ===================================================================

#include "Utils.hpp"
#include "Particles/FlexibleChain.hpp"


template<typename number>
FlexibleChain<number>::FlexibleChain()
{
	// Bounding leaf parameter
	this->BVH.SetLeafParameter(3);
	
	this->N_DELTA_L = 2;
	
	// Interaction parameters
	EPSILON_        = 1.;
	
	E_CUT_          = 100.;
	R_CUT_          = pow(2., 1./6) * this->SIGMA_R;
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void FlexibleChain<number>::Build(int mpi_rank)
{
	uint N_S;
	uint N_TOT;
	uint N_CONF;

	ArrayX<uint>     Sizes;
	Matrix3X<number> Backbones;
	
	// Load configurations from trajectory files on master thread
	if ( mpi_rank == MPI_MASTER )
	{
		std::string DATA_PATH   = __DATA_PATH;
		std::string filename_in = DATA_PATH + "/trajectory.in";
		
		Utils<number>::Load(filename_in, &Backbones, &Sizes);
		
		if ( Backbones.size() == 0 ) throw std::runtime_error("Unreadable input trajectory file");
		
		N_S    = Sizes(0);

		N_TOT  = Sizes.sum();
		N_CONF = Sizes.size();
		
		if ( (Sizes != N_S).any() ) throw std::runtime_error("Found configurations of multiple sizes in trajectory file");
	}
	
	// Broadcast data to slave threads
	MPI_Bcast(&N_S,    1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&N_TOT,  1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	MPI_Bcast(&N_CONF, 1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	
	if ( mpi_rank != MPI_MASTER )
	{
		Sizes.resize(N_CONF);
		Backbones.resize(3, N_TOT);
	}
	
	MPI_Bcast(Sizes.data(),     Sizes.size(),     Utils<uint>()  .MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	MPI_Bcast(Backbones.data(), Backbones.size(), Utils<number>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	
	number d_cc = 0.;
	
	// Set center of masses to the origin and main axes to e_z
	for ( uint idx_conf = 0; idx_conf < N_CONF; ++idx_conf )
	{
		Matrix3X<number> Backbone      = Matrix3X<number>::Map(Backbones.data() + 3 * idx_conf*N_S, 3, N_S);
		
		Vector3<number> Center_of_mass = Backbone.rowwise().mean();
		Backbone                       = Backbone.colwise() - Center_of_mass;
		
		Matrix33<number> Rot           = Utils<number>::PCA(Backbone);
		Backbone                       = Rot.transpose() * Backbone;
		
		// Work-out average bead separation distance
		for ( uint idx_s = 1; idx_s < N_S; ++idx_s )
		{
			Vector3<number> V_cc = Backbone.col(idx_s) - Backbone.col(idx_s-1);
			d_cc += V_cc.norm();
		}
		
		Backbones.block(0, N_S*idx_conf, 3, N_S) = Backbone;
	}
	
	d_cc /= (number)(N_TOT-N_CONF);
	
	// Build bounding volume hierarchy
	this->BVH.Build(Backbones, R_CUT_, Sizes);
	
	// Print simulation parameters
	if ( this->id_ == 1 )
	{
		LogTxt("Loaded particle trajectory file: %d configurations, %d interaction sites", N_CONF, N_TOT);
		
		this->BVH.PrintBuildInfo();
	}
	
	N_CONF_ = N_CONF;
	
	this->R_INTEG = 2*Backbones.colwise().norm().maxCoeff() + R_CUT_;
	this->V_INTEG = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
	
	this->V0      = PI/6.*CUB(this->SIGMA_R) * (1. + (N_S-1.)/2. * (3.*d_cc/this->SIGMA_R - CUB(d_cc/this->SIGMA_R)));
	this->V_EFF   = PI/6.*CUB(this->SIGMA_R) * (1. + (N_S-1.) * (3.*d_cc/this->SIGMA_R - CUB(d_cc/this->SIGMA_R)/2. - 3.*sqrt(1.-SQR(d_cc/(2.*this->SIGMA_R))) * asin(d_cc/(2.*this->SIGMA_R))));
}

template class FlexibleChain<float>;
template class FlexibleChain<double>;
