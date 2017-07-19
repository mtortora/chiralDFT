// ===================================================================
/**
 * DNA duplex derived particle class.
 */
// ===================================================================
/*
 * DNADuplex.cpp: Version 1.2
 * Created 22/02/2016 by Maxime Tortora
 */
// ===================================================================

#include "Utils.hpp"
#include "Particles/DNADuplex.hpp"


template<typename number>
DNADuplex<number>::DNADuplex()
{
    // Bounding leaf parameter
    this->BVH.SetLeafParameter(10);
    
    this->N_DELTA_L = 2;
    
    // Angstroms and inverse temperature to OxDNA units conversion factors
    this->SIGMA_R   = 0.11739845;
    BETA_R_         = 3000. / T_ABS;

    // B-DNA nucleotide volumes in OxDNA units - from Nadassy et al., Nucl. Acids Res. (2001)
    V_ADE_          = 136.2 * CUB(this->SIGMA_R);
    V_GUA_          = 143.8 * CUB(this->SIGMA_R);
    V_CYT_          = 113.6 * CUB(this->SIGMA_R);
    V_THY_          = 132.4 * CUB(this->SIGMA_R);
    V_BCK_          = 174.8 * CUB(this->SIGMA_R);
    
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
            R_STAR_       = 8. * this->SIGMA_R;
            
            PREFACTOR_    = COULOMB_FACTOR_ * SQR(DELTA_SCREEN_) / EPSILON_DNA_;
            DH_PREFACTOR_ = COULOMB_FACTOR_ * SQR(DELTA_SCREEN_) / EPSILON_WATER_;
            
            OFFSET_       = (PREFACTOR_    / EXCL_S1_);
            DH_OFFSET_    = (DH_PREFACTOR_ / R_STAR_) * exp(MINUS_KAPPA_*R_STAR_);
            
            R_CUT_        = 5. * LAMBDA_;
        }
        
        else throw std::runtime_error("Unsupported electrostatics model for DNA");
    }
    
    else R_CUT_ = EXCL_RC1_;
    
    DELTA_R_    = (R_CUT_ - EXCL_RC1_) / 2.;
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void DNADuplex<number>::Build(int mpi_rank)
{
    uint N_CONF;
    uint N_NUCL;
    
    ArrayX<uint>     Sizes;
    Matrix3X<number> Backbones;

    // Load configurations from trajectory files on master thread
    if ( mpi_rank == MPI_MASTER )
    {
        std::string DATA_PATH   = __DATA_PATH;
        std::string filename_in = DATA_PATH + "/trajectory.in";
        
        Utils<number>::Load(filename_in, &Backbones, &Sizes);
        
        if ( Backbones.size() == 0 ) throw std::runtime_error("Unreadable DNA input file");
        
        N_CONF = Sizes.size();
        N_NUCL = Sizes.sum();        
    }
    
    // Broadcast data to slave threads
    MPI_Bcast(&N_CONF, 1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(&N_NUCL, 1, Utils<uint>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);

    if ( mpi_rank != MPI_MASTER )
    {
        Sizes.resize(N_CONF);
        Backbones.resize(3, N_NUCL);
    }
    
    MPI_Bcast(Sizes.data(),     Sizes.size(),     Utils<uint>()  .MPI_type, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Bcast(Backbones.data(), Backbones.size(), Utils<number>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
    
    // Build bounding volume hierarchy
    this->BVH.Build(Backbones, R_CUT_, Sizes);

    // Print simulation parameters
    if ( this->id_ == 1 )
    {
        LogTxt("Loaded DNA trajectory file: %d configurations, %d nucleotides", N_CONF, N_NUCL);
        LogInf("Running temperature: %f K", T_ABS);
        
        if ( USE_DH )
        {
            LogInf("Running salt concentration: %f mol/L", (float)C_SALT);
            LogInf("Corresponding Debye length: %f nm", LAMBDA_ / (10.*this->SIGMA_R));
            
            if ( MODE_DH == DH_OXDNA )     LogInf("Using oxDNA-parametrised Debye-Huckel potential");
            if ( MODE_DH == DH_FERRARINI ) LogInf("Using Ferrarini's Debye-Huckel model potential");
        }
        
        this->BVH.PrintBuildInfo();
    }
    
    N_CONF_ = N_CONF;
    
    this->R_INTEG = 2*Backbones.colwise().norm().maxCoeff() + R_CUT_;
    this->V_INTEG = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
    
    this->V0      = ((float)N_NUCL) / ((float)N_CONF) * (V_BCK_ + (V_CYT_+V_GUA_)/2.);
    this->V_EFF   = this->V0;
}

template class DNADuplex<float>;
template class DNADuplex<double>;
