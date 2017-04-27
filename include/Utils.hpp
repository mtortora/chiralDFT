#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <vector>

#include "BNode.hpp"


namespace Utils
{
    bool OverlapBoundSC(const Eigen::Vector3d&, BNode*, BNode*);
    bool OverlapBoundOB(const Eigen::Vector3d&, BNode*, BNode*);
    
    Eigen::Matrix3d PCA(const Eigen::Matrix3Xd&);

    void Load(const std::string&, Eigen::Matrix3Xd*, Eigen::ArrayXi*);
}

#endif
