// ===================================================================
/**
 * Simulation manager class - handles computations calls and data
 * management.
 */
// ===================================================================
/*
 * SimManager.cpp: Version 2.5
 * Created 16/09/2015 by Maxime Tortora
 */
// ===================================================================

#include <string>
#include <fstream>

#include "SimManager.hpp"


// ============================
/* Class constructor */
// ============================
template<typename number>
SimManager<number>::SimManager(int mpi_rank, int mpi_size)
{
    data_path_ = __DATA_PATH;
	
    mpi_rank_  = mpi_rank;
    mpi_size_  = mpi_size;
    
    t_start_   = std::chrono::high_resolution_clock::now();
    
    // Disable standard output on slave processes
    if ( mpi_rank_ == MPI_MASTER ) setbuf(stdout, NULL);
    else                           fclose(stdout);
    
    LogTxt("*****************");
    LogPur("DFT calculations for chiral nematic LCs");
    LogTxt("Running compiled build %d.%d", __VERSION_MAJOR__, __VERSION_MINOR__);
    LogTxt("Using Eigen vsn. %d.%d.%d", EIGEN_WORLD_VERSION, EIGEN_MAJOR_VERSION, EIGEN_MINOR_VERSION);
    
    if      ( MODE_SIM == MODE_FULL ) LogTxt("Set-up for full cholesteric run");
    else if ( MODE_SIM == MODE_PERT ) LogTxt("Set-up for perturbative run");
    else if ( MODE_SIM == MODE_EXC )  LogTxt("Set-up for excluded volume calculations only");
    
    else throw std::runtime_error("Unsupported simulation scheme");
}

// ============================
/* Check MPI setup and seed RNGs */
// ============================
template<typename number>
void SimManager<number>::MPIInit()
{
    LogTxt("************");
    LogRed("Running with %u core(s)", mpi_size_);
    
    if ( mpi_rank_ == MPI_MASTER )
    {
        // Generate unique seed on master thread
        int   seed;
        FILE* tmp = fopen("/dev/urandom", "rb");
        
        if ( (tmp != NULL) && (fread((void*) &seed, sizeof(seed), 1, tmp) != 0) )
        {
            seed_ = seed;
            LogInf("Using entropy-harvested random seed: %d", seed_);
        }
        
        else
        {
            seed_ = time(NULL);
            LogInf("Using system time as RNG seed: %d", seed_);
        }
        
        fclose(tmp);
    }

    // Broadcast to slave threads
    MPI_Bcast(&seed_, 1, Utils<int>().MPI_type, MPI_MASTER, MPI_COMM_WORLD);
	
    LogBlu("Loading the %s particle template...", __MESOGEN);
    SimHandler.MCInit(seed_, mpi_rank_, mpi_size_);
}

// ============================
/* Reference run */
// ============================
template<typename number>
void SimManager<number>::InitRun()
{
    if ( MODE_SIM == MODE_EXC )
    {
        LogTxt("************");
        LogBlu("Effective excluded volume run");
        
        number v_eff = SimHandler.ExcludedIntegrator(SimHandler.IManager.SIGMA_R/2.);
        
        LogPur("Effective excluded volume of the %s particle: %f", __MESOGEN, v_eff);
    }
    
    else
    {
        LogTxt("************");
        LogBlu("Reference perturbative run");
        
        // Work-out second-virial coefficients by distributed Monte-Carlo integration
        if      ( ODF_TYPE == ODF_FULL  )    SimHandler.FullIntegrator    (&E_ref, &V_r, &V_l);
        else if ( ODF_TYPE == ODF_LEGENDRE ) SimHandler.LegendreIntegrator(&E_ref, &V_r, &V_l);
        
        else throw std::runtime_error("Unrecognised ODF computation scheme");
        
        // MPI average virial-type coefficients
        MatrixXX<number> E_ref_ = E_ref / mpi_size_;
        MPI_Allreduce(E_ref_.data(), E_ref.data(), E_ref.size(), Utils<number>().MPI_type, MPI_SUM, MPI_COMM_WORLD);

        // Work-out ODFs and Frank elastic constants
        SimHandler.BinodalAnalysis(E_ref);
        SimHandler.ODFGrid(E_ref, &P_res, &Mu_res, &F_ref, &S_res, &Psi_grd);
        SimHandler.FrankIntegrator(Psi_grd, &K1, &K2, &K3, &Kt);
        
        // Thread-safe data save in binary format
        MPI_Status  status;

        MPI_File    file_k1;
        MPI_File    file_k2;
        MPI_File    file_k3;
        MPI_File    file_kt;
        
        std::string filename_k1 = data_path_ + "/k1_threaded.out";
        std::string filename_k2 = data_path_ + "/k2_threaded.out";
        std::string filename_k3 = data_path_ + "/k3_threaded.out";
        std::string filename_kt = data_path_ + "/kt_threaded.out";
        
        // Set thread write offsets
        MPI_Offset  offset_k1   = sizeof(number)*K1.size() * mpi_rank_;
        MPI_Offset  offset_k2   = sizeof(number)*K2.size() * mpi_rank_;
        MPI_Offset  offset_k3   = sizeof(number)*K3.size() * mpi_rank_;
        MPI_Offset  offset_kt   = sizeof(number)*Kt.size() * mpi_rank_;
        
        MPI_File_open(MPI_COMM_WORLD, filename_k1.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file_k1);
        MPI_File_open(MPI_COMM_WORLD, filename_k2.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file_k2);
        MPI_File_open(MPI_COMM_WORLD, filename_k3.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file_k3);
        MPI_File_open(MPI_COMM_WORLD, filename_kt.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file_kt);

        MPI_File_seek(file_k1, offset_k1, MPI_SEEK_SET);
        MPI_File_seek(file_k2, offset_k2, MPI_SEEK_SET);
        MPI_File_seek(file_k3, offset_k3, MPI_SEEK_SET);
        MPI_File_seek(file_kt, offset_kt, MPI_SEEK_SET);

        MPI_File_write(file_k1, K1.data(), K1.size(), Utils<number>().MPI_type, &status);
        MPI_File_write(file_k2, K2.data(), K2.size(), Utils<number>().MPI_type, &status);
        MPI_File_write(file_k3, K3.data(), K3.size(), Utils<number>().MPI_type, &status);
        MPI_File_write(file_kt, Kt.data(), Kt.size(), Utils<number>().MPI_type, &status);

        MPI_File_close(&file_k1);
        MPI_File_close(&file_k2);
        MPI_File_close(&file_k3);
        MPI_File_close(&file_kt);
    
        // Equilibrium pitch from perturbation theory
        Qp_min = Kt / K2;
        
        Qp_inf = Qp_min;
        Qp_sup = Qp_min;

        Kt_inf = Kt;
        Kt_sup = Kt;
        
        K1_inf = K1;
        K1_sup = K1;
        
        K2_inf = K2;
        K2_sup = K2;
        
        K3_inf = K3;
        K3_sup = K3;
    }
}

// ============================
/* Sweeping cholesteric run */
// ============================
template<typename number>
void SimManager<number>::LandscapeRun()
{
    MatrixXX<number>   F_lnd_(N_STEPS_ETA, N_STEPS_Q);
    
    MPI_File   file_leg;
    MPI_Status status;
	
    // Iterate over q-grid
    for ( uint idx_q = 0; idx_q < N_STEPS_Q; ++idx_q )
    {
        LogTxt("************");
        LogBlu("Iteration %u out of %u", idx_q+1, N_STEPS_Q);
        
        ArrayX<number>  F_grd;
        MatrixXX<number> E_loc;
        
        // q_resc is the chiral wavevector in simulation units
        number q_macro = SimHandler.Q_grid(idx_q);
        number q_resc  = q_macro / SimHandler.IManager.SIGMA_R;
        
        SimHandler.LegendreIntegrator(&E_loc, q_resc);
        SimHandler.EnergyGrid(E_loc, &F_grd);
        
        // Thread-safe data save in binary format - decoded by script resources/processing/post_process.py
        std::string filename = data_path_ + "/legendre_matrix_" + std::to_string(q_macro) + ".out";
        
        MPI_Offset  offset   = sizeof(number)*E_loc.size() * mpi_rank_;

        MPI_File_open (MPI_COMM_WORLD, filename.c_str(), MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file_leg);
        MPI_File_seek (file_leg, offset, MPI_SEEK_SET);
        MPI_File_write(file_leg, E_loc.data(), E_loc.size(), Utils<number>().MPI_type, &status);
        MPI_File_close(&file_leg);
        
        F_lnd_.col(idx_q) = F_grd;
    }
    
    F_lnd = F_lnd_;
}

// ============================
/* Locate potential surface minima */
// ============================
template<typename number>
void SimManager<number>::MinSurf()
{
    IndexX<number> min_q;
    ArrayX<number> Q_min_(N_STEPS_ETA);
    
    for ( uint idx_eta = 0; idx_eta < N_STEPS_ETA; ++idx_eta )
    {
        F_lnd.row(idx_eta).minCoeff(&min_q);
        Q_min_(idx_eta) = SimHandler.Q_grid(min_q);
    }
        
    Q_min = Q_min_;
    Q_inf = Q_min_;
    Q_sup = Q_min_;
}

// ============================
/* MPI averages over all threads */
// ============================
template<typename number>
void SimManager<number>::Gather()
{
    // Reduced variables for MPI sum
    ArrayX<number>  P_res_  = P_res  / mpi_size_;
    ArrayX<number>  Mu_res_ = Mu_res / mpi_size_;
    
    ArrayX<number>  F_ref_  = F_ref  / mpi_size_;
    ArrayX<number>  S_res_  = S_res  / mpi_size_;
    
    ArrayX<number>  V_r_    = V_r    / mpi_size_;
    ArrayX<number>  V_l_    = V_l    / mpi_size_;
    
    ArrayX<number>  K1_     = K1     / mpi_size_;
    ArrayX<number>  K2_     = K2     / mpi_size_;
    ArrayX<number>  K3_     = K3     / mpi_size_;
    ArrayX<number>  Kt_     = Kt     / mpi_size_;

    ArrayX<number>  Qp_min_ = Qp_min / mpi_size_;

    ArrayX<number>  Qp_inf_ = Qp_inf;
    ArrayX<number>  Qp_sup_ = Qp_sup;
    
    ArrayX<number>  K1_inf_ = K1_inf;
    ArrayX<number>  K1_sup_ = K1_sup;
    
    ArrayX<number>  K2_inf_ = K2_inf;
    ArrayX<number>  K2_sup_ = K2_sup;
    
    ArrayX<number>  K3_inf_ = K3_inf;
    ArrayX<number>  K3_sup_ = K3_sup;
    
    ArrayX<number>  Kt_inf_ = Kt_inf;
    ArrayX<number>  Kt_sup_ = Kt_sup;
    
    if ( MODE_SIM == MODE_FULL )
    {
        MatrixXX<number> F_lnd_ = F_lnd / mpi_size_;

        ArrayX<number>  Q_min_ = Q_min / mpi_size_;
        ArrayX<number>  Q_inf_ = Q_inf;
        ArrayX<number>  Q_sup_ = Q_sup;
        
        MPI_Reduce(F_lnd_.data(), F_lnd.data(), F_lnd.size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
        
        // Equilibrium pitch from full self-consistent run
        MPI_Reduce(Q_min_.data(), Q_min.data(), Q_min.size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
        MPI_Reduce(Q_inf_.data(), Q_inf.data(), Q_inf.size(), Utils<number>().MPI_type, MPI_MIN, MPI_MASTER, MPI_COMM_WORLD);
        MPI_Reduce(Q_sup_.data(), Q_sup.data(), Q_sup.size(), Utils<number>().MPI_type, MPI_MAX, MPI_MASTER, MPI_COMM_WORLD);
    }

    MPI_Reduce(V_r_   .data(), V_r   .data(), V_r   .size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(V_l_   .data(), V_l   .data(), V_l   .size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    
    // Thermodynamic variables
    MPI_Reduce(P_res_ .data(), P_res .data(), P_res .size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(Mu_res_.data(), Mu_res.data(), Mu_res.size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(F_ref_ .data(), F_ref .data(), F_ref .size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(S_res_ .data(), S_res .data(), S_res .size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);

    // Equilibrium pitch from perturbation theory
    MPI_Reduce(Qp_min_.data(), Qp_min.data(), Qp_min.size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(Qp_inf_.data(), Qp_inf.data(), Qp_inf.size(), Utils<number>().MPI_type, MPI_MIN, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(Qp_sup_.data(), Qp_sup.data(), Qp_sup.size(), Utils<number>().MPI_type, MPI_MAX, MPI_MASTER, MPI_COMM_WORLD);

    // Frank elastic constants
    MPI_Reduce(K1_    .data(), K1    .data(), K1.    size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(K1_inf_.data(), K1_inf.data(), K1_inf.size(), Utils<number>().MPI_type, MPI_MIN, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(K1_sup_.data(), K1_sup.data(), K1_sup.size(), Utils<number>().MPI_type, MPI_MAX, MPI_MASTER, MPI_COMM_WORLD);
    
    MPI_Reduce(K2_    .data(), K2    .data(), K2.    size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(K2_inf_.data(), K2_inf.data(), K2_inf.size(), Utils<number>().MPI_type, MPI_MIN, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(K2_sup_.data(), K2_sup.data(), K2_sup.size(), Utils<number>().MPI_type, MPI_MAX, MPI_MASTER, MPI_COMM_WORLD);
    
    MPI_Reduce(K3_    .data(), K3    .data(), K3.    size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(K3_inf_.data(), K3_inf.data(), K3_inf.size(), Utils<number>().MPI_type, MPI_MIN, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(K3_sup_.data(), K3_sup.data(), K3_sup.size(), Utils<number>().MPI_type, MPI_MAX, MPI_MASTER, MPI_COMM_WORLD);
    
    MPI_Reduce(Kt_    .data(), Kt    .data(), Kt    .size(), Utils<number>().MPI_type, MPI_SUM, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(Kt_inf_.data(), Kt_inf.data(), Kt_inf.size(), Utils<number>().MPI_type, MPI_MIN, MPI_MASTER, MPI_COMM_WORLD);
    MPI_Reduce(Kt_sup_.data(), Kt_sup.data(), Kt_sup.size(), Utils<number>().MPI_type, MPI_MAX, MPI_MASTER, MPI_COMM_WORLD);
}

// ============================
/* Save aggregated data on master thread */
// ============================
template<typename number>
void SimManager<number>::Save()
{
    if ( mpi_rank_ == MPI_MASTER )
    {
        // Output files
        std::ofstream file_exc(data_path_ + "/excluded_matrix.out");
        std::ofstream file_ord(data_path_ + "/order_param.out");
        std::ofstream file_dv (data_path_ + "/delta_v.out");
        std::ofstream file_df (data_path_ + "/delta_f.out");
        std::ofstream file_per(data_path_ + "/q_pert.out");
        std::ofstream file_psi(data_path_ + "/psi.out");
        std::ofstream file_k1 (data_path_ + "/k1.out");
        std::ofstream file_k2 (data_path_ + "/k2.out");
        std::ofstream file_k3 (data_path_ + "/k3.out");
        std::ofstream file_kt (data_path_ + "/kt.out");
        std::ofstream file_mu (data_path_ + "/mu.out");
        std::ofstream file_p  (data_path_ + "/p.out");
        
        // Legendre reference matrix
        file_exc << "Particle volume: "  << SimHandler.IManager.V0    << std::endl;
        file_exc << "Effective volume: " << SimHandler.IManager.V_EFF << std::endl;

        file_exc << E_ref;
        
        // Angle-dependant excluded volume, assuming head-tail particle symmetry
        ArrayX<number> V_chi  = V_r - V_l;
        ArrayX<number> V_ave  = (V_chi - V_chi.reverse()) / 2.;

        ArrayX<number> V_nrm  = V_ave / (V_r + V_l);
        V_ave.tail(N_STEPS_THETA/2) = V_ave.head(N_STEPS_THETA/2).reverse();
        
        // Thermodynamically-averaged excluded volume
        ArrayX<number> N_grid = SimHandler.Eta_grid * CUB(SimHandler.IManager.SIGMA_R)/SimHandler.IManager.V0;
        ArrayX<number> F_chi  = (Psi_grd.rowwise() * (V_ave * sin(SimHandler.Theta_grid)).transpose()).rowwise().sum();
        
        F_chi         *= -N_grid * 4.*SQR(PI) * D_THETA;

        for ( uint idx_theta = 0; idx_theta < N_STEPS_THETA; ++idx_theta )
        {
            number theta = SimHandler.Theta_grid(idx_theta);
            
            file_dv << theta << ' ' << V_nrm(idx_theta) << ' ' << V_chi(idx_theta) << std::endl;
        }
        
        for ( uint idx_eta = 0; idx_eta < N_STEPS_ETA; ++idx_eta )
        {
            number eta = SimHandler.Eta_grid(idx_eta);
            
            file_df  << eta << ' ' << F_chi (idx_eta) << std::endl;
            
            // Torque field, elastic constants and equilibrium pitch from perturbation theory
            file_k1  << eta << ' ' << K1    (idx_eta) << ' ' << K1_inf(idx_eta) << ' ' << K1_sup(idx_eta) << std::endl;
            file_k2  << eta << ' ' << K2    (idx_eta) << ' ' << K2_inf(idx_eta) << ' ' << K2_sup(idx_eta) << std::endl;
            file_k3  << eta << ' ' << K3    (idx_eta) << ' ' << K3_inf(idx_eta) << ' ' << K3_sup(idx_eta) << std::endl;
            file_kt  << eta << ' ' << Kt    (idx_eta) << ' ' << Kt_inf(idx_eta) << ' ' << Kt_sup(idx_eta) << std::endl;

            file_per << eta << ' ' << Qp_min(idx_eta) << ' ' << Qp_inf(idx_eta) << ' ' << Qp_sup(idx_eta) << std::endl;
            
            // Nematic order parameter, pressure and chemical potential
            file_ord << eta << ' ' << S_res (idx_eta) << std::endl;
            file_mu  << eta << ' ' << Mu_res(idx_eta) << std::endl;
            file_p   << eta << ' ' << P_res (idx_eta) << std::endl;
            
            // Save ODFs
            for ( uint idx_theta = 0; idx_theta < N_STEPS_THETA; ++idx_theta )
            {
                number theta = SimHandler.Theta_grid(idx_theta);
                
                file_psi << eta << ' ' << theta << ' ' << Psi_grd(idx_eta, idx_theta) << std::endl;
            }
            
            file_psi << std::endl;
        }
        
        if ( MODE_SIM == MODE_FULL )
        {
            std::ofstream file_lnd(data_path_ + "/energy_landscape.out");
            std::ofstream file_min(data_path_ + "/q_full.out");

            // Energy landscape and full-run equilibrium pitch
            for ( uint idx_eta = 0; idx_eta < N_STEPS_ETA; ++idx_eta )
            {
                number eta = SimHandler.Eta_grid(idx_eta);
                
                file_min  << eta << ' ' << Q_min(idx_eta) << ' ' << Q_inf(idx_eta) << ' '<< Q_sup(idx_eta) << std::endl;

                for ( uint idx_q = 0; idx_q < N_STEPS_Q; ++idx_q )
                {
                    number q_macro = SimHandler.Q_grid(idx_q);
                    number delta   = F_lnd(idx_eta, idx_q) - F_ref(idx_eta);
                    
                    file_lnd << eta << ' ' << q_macro << ' ' << delta << std::endl;
                }
                
                file_lnd << std::endl;
            }
            
            file_lnd.close();
            file_min.close();
        }
        
        file_exc.close();
        file_ord.close();
        file_dv .close();
        file_df .close();
        file_per.close();
        file_psi.close();
        file_k1 .close();
        file_k2 .close();
        file_k3 .close();
        file_kt .close();
        file_mu .close();
        file_p  .close();
    }
}

// ============================
/* Class destructor */
// ============================
template<typename number>
SimManager<number>::~SimManager()
{
    t_end_     = std::chrono::high_resolution_clock::now();
    t_elapsed_ = t_end_ - t_start_;
    
    LogTxt("*****************");
    LogBlu("Total runtime: %fh", t_elapsed_.count() / 3600.);
}

template class SimManager<float>;
template class SimManager<double>;
