#ifndef SIM_MANAGER_HPP_
#define SIM_MANAGER_HPP_

#include <chrono>

#include "OdfOptimiser.hpp"


class SimManager
{
private:
    int seed_;

    int mpi_rank_;
    int mpi_size_;
    
    std::string data_path_;
    
    std::chrono::high_resolution_clock::time_point t_start_;
    std::chrono::high_resolution_clock::time_point t_end_;
    
    std::chrono::duration<double> t_elapsed_;
    
public:
    SimManager(int, int);
    
    Eigen::ArrayXd  Qp_min;
    Eigen::ArrayXd  Qp_inf;
    Eigen::ArrayXd  Qp_sup;
    
    Eigen::ArrayXd  Q_min;
    Eigen::ArrayXd  Q_inf;
    Eigen::ArrayXd  Q_sup;
    
    Eigen::ArrayXd  F_ref;
    
    Eigen::ArrayXd  P_res;
    Eigen::ArrayXd  Mu_res;
    Eigen::ArrayXd  S_res;
    Eigen::ArrayXd  V_res;
    
    Eigen::ArrayXd  Kt;
    Eigen::ArrayXd  Kt_inf;
    Eigen::ArrayXd  Kt_sup;
    
    Eigen::ArrayXd  K1;
    Eigen::ArrayXd  K1_inf;
    Eigen::ArrayXd  K1_sup;
    
    Eigen::ArrayXd  K2;
    Eigen::ArrayXd  K2_inf;
    Eigen::ArrayXd  K2_sup;
    
    Eigen::ArrayXd  K3;
    Eigen::ArrayXd  K3_inf;
    Eigen::ArrayXd  K3_sup;
    
    Eigen::ArrayXd  F_res;

    Eigen::MatrixXd E_ref;
    Eigen::MatrixXd F_lnd;
    
    OdfOptimiser<MESOGEN> SimHandler;
    
    void MPIInit();
    void InitRun(int);
    void LandscapeRun();
    void MinSurf();
    void Gather();
    void Save();
    
    ~SimManager();
};

#endif
