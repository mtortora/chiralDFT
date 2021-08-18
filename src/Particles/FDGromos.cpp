// ===================================================================
/**
 * FDGromos derived particle class.
 */
// ===================================================================
/*
 * FDGromos.cpp: Version 1.0
 * Created 08/06/2021 by Maxime Tortora
 */
// ===================================================================

#include "Utils.hpp"
#include "Particles/FDGromos.hpp"


template<typename number>
FDGromos<number>::FDGromos()
{
	// Bounding leaf parameter
	this->BVH.SetLeafParameter(17);
	
	this->N_DELTA_L = 2;
	
	// Cutoff distance (in Angstroms)
	R_CUT_  = 35.;
	
	// Cutoff energy (in units kT)
	E_CUT_  = 20.;

	// Solvent dielectric constant
	number epsilon = 80.;
	
	number kappa   = 50.29037937112641 * sqrt(C_SALT / (epsilon*T_ABS));
	number e_conv  = 1.2027239485976831e8 / T_ABS;
		
	DH_PREFACTOR_  = 167101.00207275996 / T_ABS;
	C_RF_          = (2*(1-epsilon) * (1+kappa*R_CUT_) - epsilon * SQR(kappa*R_CUT_)) / ((1+2*epsilon) * (1+kappa*R_CUT_) + epsilon * SQR(kappa*R_CUT_));
	
	/* GROMOS 53A6 force field parameters */
	// From Oostenbrink et al., J. Comput. Chem. (2004)
	
	// C6
	C6[0][0]  = 0.07790;
	C6[0][1]  = 0.08205;
	C6[0][2]  = 0.08740;
	C6[0][3]  = 0.08168;
	C6[0][4]  = 0.07605;
	C6[0][5]  = 0.06201;
	C6[0][6]  = 0.06201;
	C6[0][7]  = 0.06201;
	C6[0][8]  = 0.06087;
	C6[0][9]  = 0.06087;
	C6[0][10] = 0.06087;
	C6[0][11] = 0.08823;

	C6[1][1]  = 0.08642;
	C6[1][2]  = 0.09205;
	C6[1][3]  = 0.08603;
	C6[1][4]  = 0.08010;
	C6[1][5]  = 0.06531;
	C6[1][6]  = 0.06531;
	C6[1][7]  = 0.06531;
	C6[1][8]  = 0.06411;
	C6[1][9]  = 0.06411;
	C6[1][10] = 0.06411;
	C6[1][11] = 0.09293;
	
	C6[2][2]  = 0.09805;
	C6[2][3]  = 0.09164;
	C6[2][4]  = 0.08532;
	C6[2][5]  = 0.06957;
	C6[2][6]  = 0.06957;
	C6[2][7]  = 0.06957;
	C6[2][8]  = 0.06828;
	C6[2][9]  = 0.06828;
	C6[2][10] = 0.06828;
	C6[2][11] = 0.09898;
	
	C6[3][3]  = 0.08564;
	C6[3][4]  = 0.07974;
	C6[3][5]  = 0.06502;
	C6[3][6]  = 0.06502;
	C6[3][7]  = 0.06502;
	C6[3][8]  = 0.06382;
	C6[3][9]  = 0.06382;
	C6[3][10] = 0.06382;
	C6[3][11] = 0.09250;

	C6[4][4]  = 0.07425;
	C6[4][5]  = 0.06054;
	C6[4][6]  = 0.06054;
	C6[4][7]  = 0.06054;
	C6[4][8]  = 0.05942;
	C6[4][9]  = 0.05942;
	C6[4][10] = 0.05942;
	C6[4][11] = 0.08613;
	
	C6[5][5]  = 0.04936;
	C6[5][6]  = 0.04936;
	C6[5][7]  = 0.04936;
	C6[5][8]  = 0.04845;
	C6[5][9]  = 0.04845;
	C6[5][10] = 0.04845;
	C6[5][11] = 0.07023;
	
	C6[6][6]  = 0.04936;
	C6[6][7]  = 0.04936;
	C6[6][8]  = 0.04845;
	C6[6][9]  = 0.04845;
	C6[6][10] = 0.04845;
	C6[6][11] = 0.07023;

	C6[7][7]  = 0.04936;
	C6[7][8]  = 0.04845;
	C6[7][9]  = 0.04845;
	C6[7][10] = 0.04845;
	C6[7][11] = 0.07023;
	
	C6[8][8]  = 0.04756;
	C6[8][9]  = 0.04756;
	C6[8][10] = 0.04756;
	C6[8][11] = 0.06894;
	
	C6[9][9]  = 0.04756;
	C6[9][10] = 0.04756;
	C6[9][11] = 0.06894;
	
	C6[10][10] = 0.04756;
	C6[10][11] = 0.06894;
	
	C6[11][11] = 0.09992;

	// C12
	C12[0][0]  = 9.850;
	C12[0][1]  = 7.577;
	C12[0][2]  = 7.131;
	C12[0][3]  = 7.223;
	C12[0][4]  = 6.188;
	C12[0][5]  = 3.873;
	C12[0][6]  = 3.873;
	C12[0][7]  = 3.873;
	C12[0][8]  = 3.138;
	C12[0][9]  = 3.292;
	C12[0][10] = 2.912;
	C12[0][11] = 5.968;

	C12[1][1]  = 5.828;
	C12[1][2]  = 5.485;
	C12[1][3]  = 5.556;
	C12[1][4]  = 4.760;
	C12[1][5]  = 2.979;
	C12[1][6]  = 2.979;
	C12[1][7]  = 2.979;
	C12[1][8]  = 2.414;
	C12[1][9]  = 2.532;
	C12[1][10] = 2.240;
	C12[1][11] = 4.591;
	
	C12[2][2]  = 5.162;
	C12[2][3]  = 5.229;
	C12[2][4]  = 4.480;
	C12[2][5]  = 2.804;
	C12[2][6]  = 2.804;
	C12[2][7]  = 2.804;
	C12[2][8]  = 2.272;
	C12[2][9]  = 2.383;
	C12[2][10] = 2.108;
	C12[2][11] = 4.320;
	
	C12[3][3]  = 5.297;
	C12[3][4]  = 4.538;
	C12[3][5]  = 2.840;
	C12[3][6]  = 2.840;
	C12[3][7]  = 2.840;
	C12[3][8]  = 2.302;
	C12[3][9]  = 2.414;
	C12[3][10] = 2.136;
	C12[3][11] = 4.377;

	C12[4][4]  = 3.888;
	C12[4][5]  = 2.433;
	C12[4][6]  = 2.433;
	C12[4][7]  = 2.433;
	C12[4][8]  = 1.972;
	C12[4][9]  = 2.068;
	C12[4][10] = 1.830;
	C12[4][11] = 3.750;
	
	C12[5][5]  = 1.523;
	C12[5][6]  = 1.891;
	C12[5][7]  = 2.091;
	C12[5][8]  = 1.482;
	C12[5][9]  = 1.544;
	C12[5][10] = 1.891;
	C12[5][11] = 2.347;
	
	C12[6][6]  = 1.841;
	C12[6][7]  = 2.035;
	C12[6][8]  = 1.442;
	C12[6][9]  = 1.503;
	C12[6][10] = 1.841;
	C12[6][11] = 2.347;

	C12[7][7]  = 2.250;
	C12[7][8]  = 1.595;
	C12[7][9]  = 1.662;
	C12[7][10] = 2.035;
	C12[7][11] = 2.347;
	
	C12[8][8]  = 1.000;
	C12[8][9]  = 1.178;
	C12[8][10] = 0.928;
	C12[8][11] = 1.902;
	
	C12[9][9]  = 1.227;
	C12[9][10] = 1.503;
	C12[9][11] = 1.994;
	
	C12[10][10] = 0.8611;
	C12[10][11] = 1.765;
	
	C12[11][11] = 3.616;

	for ( int i = 0; i < 12; ++i )
	{
		for ( int j = 0; j < i; ++j )
		{
			C6[i][j] = C6[j][i];
			C12[i][j] = C12[j][i];
		}
	}
		
	for ( int i = 0; i < 12; ++i )
	{
		for ( int j = 0; j < 12; ++j )
		{
			C6[i][j] = SQR(C6[i][j]) * e_conv;
			C12[i][j] = SQR(C12[i][j]) * e_conv;
		}
	}
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void FDGromos<number>::Build(int mpi_rank)
{
	uint N_S;
	uint N_TOT;
	uint N_CONF;

	ArrayX<uint>     Sizes, Types;
	ArrayX<number>   Charges;
	
	Matrix3X<number> Backbones;
	
	// Load configurations from trajectory files on master thread
	if ( mpi_rank == MPI_MASTER )
	{
		std::string DATA_PATH   = __DATA_PATH;
		std::string filename_in = DATA_PATH + "/trajectory.in";
		
		Utils<number>::Load(filename_in, &Backbones, &Charges, &Types, &Sizes);
		
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
		Types.resize(N_TOT);
		Charges.resize(N_TOT);
		Backbones.resize(3, N_TOT);

		Sizes.resize(N_CONF);
	}
	
	MPI_Bcast(Types.data(),     Types.size(),     Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	MPI_Bcast(Sizes.data(),     Sizes.size(),     Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	
	MPI_Bcast(Charges.data(), Charges.size(), Utils<number>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	MPI_Bcast(Backbones.data(), Backbones.size(), Utils<number>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	
	number v0 = 0.;
	
	// Set center of masses to the origin and main axes to e_z
	for ( uint idx_conf = 0; idx_conf < N_CONF; ++idx_conf )
	{
		Matrix3X<number> Backbone      = Matrix3X<number>::Map(Backbones.data() + 3 * idx_conf*N_S, 3, N_S);
		
		Vector3<number> Center_of_mass = Backbone.rowwise().mean();
		Backbone                       = Backbone.colwise() - Center_of_mass;
		
		Matrix33<number> Rot           = Utils<number>::PCA(Backbone);
		Backbone                       = Rot.transpose() * Backbone;
		
		ArrayX<number> Vertex_z        = Backbone.row(2).array();
		
		number r_max = Backbone.block(0, 0, 2, Backbone.cols()).colwise().norm().maxCoeff();
		number lz = Vertex_z.maxCoeff() - Vertex_z.minCoeff();
		
		v0 += PI * SQR(r_max) * lz;
		
		Backbones.block(0, N_S*idx_conf, 3, N_S) = Backbone;
	}
	
	v0 /= (number)(N_CONF);
	
	// Build bounding volume hierarchy
	this->BVH.Build(Backbones, Charges, Types, R_CUT_, Sizes);
	
	// Print simulation parameters
	if ( this->id_ == 1 )
	{
		LogTxt("Loaded particle trajectory file: %d configurations, %d interaction sites", N_CONF, N_TOT);
		
		this->BVH.PrintBuildInfo();
	}
	
	N_CONF_       = N_CONF;
	
	this->R_INTEG = 2*Backbones.colwise().norm().maxCoeff() + R_CUT_;
	this->V_INTEG = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
	
	this->V0      = v0;
	this->V_EFF   = v0;
}

template class FDGromos<float>;
template class FDGromos<double>;
