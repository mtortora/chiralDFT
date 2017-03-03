// ===================================================================
/**
 * DNA duplex derived particle class
 */
// ===================================================================
/*
 * DNADuplex.cpp: Version 1.2
 * Created 22/02/2016 by Maxime Tortora
 */
// ===================================================================

#include <mpi.h>

#include "Utils.hpp"
#include "Particles/DNADuplex.hpp"

using namespace Eigen;


DNADuplex::DNADuplex()
{
    // Bounding tree properties
    BHierarchy->SetTreeProperties(4);
    
    N_DELTA_L       = 2;
    
    // Angstroms and inverse temperature to OxDNA units conversion factors
    SIGMA_R         = 0.11739845;
    BETA_R_         = 3000. / T_ABS;

    // B-DNA nucleotide volumes in OxDNA units - from Nadassy et al., Nucl. Acids Res. (2001)
    V_ADE_          = 136.2 * CUB(SIGMA_R);
    V_GUA_          = 143.8 * CUB(SIGMA_R);
    V_CYT_          = 113.6 * CUB(SIGMA_R);
    V_THY_          = 132.4 * CUB(SIGMA_R);
    V_BCK_          = 174.8 * CUB(SIGMA_R);
    
    // oxDNA interaction parameters
    N_STAR_         = 3;
    
    EXCL_EPS_       = 2.;
    EXCL_S1_        = 0.7;
    EXCL_R1_        = 0.675;
    EXCL_B1_        = 892.016223343;
    EXCL_RC1_       = 0.711879214356;

    // Cutoff energy
    E_CUT_          = 20. / BETA_R_;
    
    // Ferrarini interaction parameters
    DELTA_SCREEN_   = 0.24;
    COULOMB_FACTOR_ = 6.5390469482;
    
    // Dielectric permittivities
    EPSILON_DNA_    = 2.;
    EPSILON_WATER_  = 80.;
    
    // Debye-Huckel parameters
    if ( USE_DH )
    {
        // Debye length
        LAMBDA_      = 0.002334412 * sqrt(T_ABS) * sqrt(EPSILON_WATER_) / sqrt(C_SALT);
        MINUS_KAPPA_ = -1. / LAMBDA_;
        
        if ( MODE_DH == DH_OXDNA )
        {
            R_STAR_       = N_STAR_ * LAMBDA_;
            DH_PREFACTOR_ = 0.0543;
            
            B_CUT_        = DH_PREFACTOR_*exp(MINUS_KAPPA_*R_STAR_) * SQR(R_STAR_+LAMBDA_)/(4.*CUB(R_STAR_)*SQR(LAMBDA_));
            R_CUT_        = 2. * LAMBDA_ * SQR(N_STAR_) / (N_STAR_ + 1.);
        }
        
        // DNA electrostatics model of Tombolato et al., JCP (2005)
        else if ( MODE_DH == DH_FERRARINI )
        {
            R_STAR_       = 8. * SIGMA_R;
            
            PREFACTOR_    = COULOMB_FACTOR_ * SQR(DELTA_SCREEN_) / EPSILON_DNA_;
            DH_PREFACTOR_ = COULOMB_FACTOR_ * SQR(DELTA_SCREEN_) / EPSILON_WATER_;
            
            OFFSET_       = (PREFACTOR_    / EXCL_S1_);
            DH_OFFSET_    = (DH_PREFACTOR_ / R_STAR_) * exp(MINUS_KAPPA_*R_STAR_);
            
            R_CUT_        = 10. * LAMBDA_;
        }
        
        else throw std::runtime_error("Unsupported electrostatics model for DNA");
    }
    
    else R_CUT_ = EXCL_RC1_;
    
    DELTA_R_    = (R_CUT_ - EXCL_RC1_) / 2.;
}

// ============================
/* Build particle model */
// ============================
void DNADuplex::Build(int mpi_rank)
{
    uint      N_CONF;
    uint      N_NUCL;
    
    Matrix3Xd Backbones;

    // Load configurations from trajectory files on master thread
    if ( mpi_rank == MPI_MASTER )
    {
        uint n_vertices;
        
        std::string DATA_PATH   = __DATA_PATH;
        std::string filename_in = DATA_PATH + "/trajectory.in";
        
        Utils::Load(filename_in, &Backbones, &n_vertices);
        
        if ( n_vertices == 0 ) throw std::runtime_error("Unreadable DNA input file");
        
        N_CONF = Backbones.cols() / n_vertices;
        N_NUCL = n_vertices;
    }
    
    // Broadcast data to slave threads
    MPI_Bcast(&N_CONF, 1, MPI_UNSIGNED, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&N_NUCL, 1, MPI_UNSIGNED, MPI_MASTER, MPI_COMM_WORLD);

    if ( mpi_rank != MPI_MASTER ) Backbones.resize(3, N_CONF*N_NUCL);
    
    MPI_Bcast(Backbones.data(), Backbones.size(), MPI_DOUBLE, MPI_MASTER, MPI_COMM_WORLD);
    
    // Build bounding volume hierarchy
    BHierarchy->AllocateForest(N_CONF);
    
    for ( uint idx_conf = 0; idx_conf < N_CONF; ++idx_conf )
    {
        BTree* Tree        = &BHierarchy->Forest[idx_conf];
        Matrix3Xd Vertices = Matrix3Xd::Map(Backbones.data() + 3 * idx_conf*N_NUCL, 3, N_NUCL);
        
        BHierarchy->RecursiveBuild(Tree, Vertices, R_CUT_);
    }
    
    // Print simulation parameters
    if ( id_ == 1 )
    {
        LogTxt("Loaded DNA trajectory file: %d configurations, %d nucleotides", N_CONF, N_NUCL);
        LogInf("Running temperature: %f K", T_ABS);
        
        if ( USE_DH )
        {
            LogInf("Running salt concentration: %f mol/L", (float)C_SALT);
            LogInf("Corresponding Debye length: %f nm", LAMBDA_ / (10.*SIGMA_R));
            
            if ( MODE_DH == DH_OXDNA )     LogInf("Using oxDNA-parametrised Debye-Huckel potential");
            if ( MODE_DH == DH_FERRARINI ) LogInf("Using Ferrarini's Debye-Huckel model potential");
        }
        
        BHierarchy->PrintBuildInfo();
    }
    
    // Set particle properties
    N_CONF_ = N_CONF;
    
    R_INTEG = 2*Backbones.colwise().norm().maxCoeff() + R_CUT_;
    V_INTEG = CUB(2.*R_INTEG) * 16.*pow(PI, 6);
    
    V0      = N_NUCL * (V_BCK_ + (V_CYT_+V_GUA_)/2.);
    V_EFF   = V0;
}
