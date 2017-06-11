#ifndef SIM_MANAGER_HPP_
#define SIM_MANAGER_HPP_

#include <chrono>

#include "OdfOptimiser.hpp"


template<typename number>
class SimManager
{
public:
    SimManager(int, int);
    ~SimManager();

    void MPIInit();
    void InitRun();
    void LandscapeRun();
    void MinSurf();
    void Gather();
    void Save();
    
private:
    int seed_;

    int mpi_rank_;
    int mpi_size_;
    
    void MPISave(const std::string&, const ArrayX<number>&);
    
    OdfOptimiser<MESOGEN, number> SimHandler;

    std::string data_path_;
    
    std::chrono::high_resolution_clock::time_point t_start_;
    std::chrono::high_resolution_clock::time_point t_end_;
    
    std::chrono::duration<double> t_elapsed_;
    
    ArrayX<number>  V_r;
    ArrayX<number>  V_l;
    
    ArrayX<number>  Qp_min;
    ArrayX<number>  Qp_inf;
    ArrayX<number>  Qp_sup;
    
    ArrayX<number>  Q_min;
    ArrayX<number>  Q_inf;
    ArrayX<number>  Q_sup;
    
    ArrayX<number>  F_ref;
    
    ArrayX<number>  P_res;
    ArrayX<number>  Mu_res;
    ArrayX<number>  S_res;
    
    ArrayX<number>  K1;
    ArrayX<number>  K1_inf;
    ArrayX<number>  K1_sup;
    
    ArrayX<number>  K2;
    ArrayX<number>  K2_inf;
    ArrayX<number>  K2_sup;
    
    ArrayX<number>  K3;
    ArrayX<number>  K3_inf;
    ArrayX<number>  K3_sup;
    
    ArrayX<number>  Kt;
    ArrayX<number>  Kt_inf;
    ArrayX<number>  Kt_sup;
    
    ArrayXX<number> Psi_grd;
    
    ArrayX<number>  E_ref;
    ArrayXX<number> F_lnd;
};

#endif
