#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <mpi.h>
#include <vector>

#include "BNode.hpp"


template<typename number>
struct Utils
{
    Utils();
    
    MPI_Datatype MPI_type;

    static bool OverlapBoundSC(const Vector3<number>&, BNode<number>*, BNode<number>*);
    static bool OverlapBoundOB(const Vector3<number>&, BNode<number>*, BNode<number>*);
    
    static Matrix33<number> PCA(const Matrix3X<number>&);

    static void Load(const std::string&, Matrix3X<number>*, ArrayX<uint>*);
};

#endif
