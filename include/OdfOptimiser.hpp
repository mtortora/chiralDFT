#ifndef ODF_OPTIMISER_HPP_
#define ODF_OPTIMISER_HPP_

#include "MCIntegrator.hpp"


template<class ParticleType>
class OdfOptimiser: public MCIntegrator<ParticleType>
{    
private:
    Eigen::ArrayXd Psi_iso_;
    
    double VirialCoeff(const Eigen::ArrayXd&, const Eigen::MatrixXd&, int);
    double FreeEnergy(double, const Eigen::ArrayXd&, const Eigen::MatrixXd&, int);
    
    Eigen::ArrayXd LegendreCoeffs(const Eigen::ArrayXd&);
    Eigen::ArrayXd SequentialOptimiser(double, const Eigen::MatrixXd&, int);
    
    Eigen::Vector2d ODFThermo(double, const Eigen::ArrayXd&, const Eigen::MatrixXd&, int);

public:
    OdfOptimiser();

    Eigen::ArrayXd Theta_grid;

    void BinodalAnalysis(const Eigen::MatrixXd&, int);
    void EnergyGrid(const Eigen::MatrixXd&, Eigen::ArrayXd*, int);
    void ODFGrid(const Eigen::MatrixXd&,
                 Eigen::ArrayXd*, Eigen::ArrayXd*, Eigen::ArrayXd*, Eigen::ArrayXd*,
                 Eigen::MatrixXd*, int);
};

#endif
