// ===================================================================
/**
 * Implements robust overlap methods for convex bounding volumes.
 * Also includes parser for processed trajectory files.
 */
// ===================================================================
/*
 * Utils.cpp: Version 1.0
 * Created 26/07/2016 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Utils.hpp"


// ============================
/* Wrappers for MPI datatypes */
// ============================
template<> Utils<float> ::Utils(): MPI_type(MPI_FLOAT)              {}
template<> Utils<double>::Utils(): MPI_type(MPI_DOUBLE)             {}

template<> Utils<int>   ::Utils(): MPI_type(MPI_INT)                {}
template<> Utils<uint>  ::Utils(): MPI_type(MPI_UNSIGNED)           {}
template<> Utils<ullint>::Utils(): MPI_type(MPI_UNSIGNED_LONG_LONG) {}

// ============================
/* Bounding spherocylinders overlap test */
// ============================
template<typename number>
bool Utils<number>::OverlapBoundSC(const Vector3<number>& R_cm, BNode<number>* Node1, BNode<number>* Node2)
{
    number xlambda;
    number xmu;
    number aux1;
    number aux2;
    number rm2;

    Vector3<number> R_sep;
    
    // Root spherocylinders centers are forced to the particle center of mass 
    if      ( Node1->is_root && Node2->is_root ) R_sep = R_cm;
    else if ( Node1->is_root )                   R_sep = R_cm + Node2->Center;
    else if ( Node2->is_root )                   R_sep = R_cm - Node1->Center;
    else                                         R_sep = R_cm + Node2->Center - Node1->Center;

    number rEo1  = Node1->Axis.dot(R_sep);
    number rEo2  = Node2->Axis.dot(R_sep);;
    number o1Eo2 = Node1->Axis.dot(Node2->Axis);
    
    number rsqr  = SQR(Node1->l_cr + Node2->l_cr);
    number cc    = 1. - SQR(o1Eo2);
    
    /* VEGA ET AL., A FAST ALGORITHM TO EVALUATE THE SHORTEST DISTANCE BETWEEN RODS */
    // Case Axis1, Axis2 colinear
    if ( cc < TOL_SC )
    {
        if ( std::abs(rEo1) > TOL_SC )
        {
            xlambda = copysign(Node1->l_ch, rEo1);
            xmu     = xlambda*o1Eo2 - rEo2;
            
            if ( std::abs(xmu) > Node2->l_ch ) xmu = copysign(Node2->l_ch, xmu);
        }
        
        else
        {
            xlambda = 0.;
            xmu     = 0.;
        }
    }
    
    // General case
    else
    {
        xlambda = ( rEo1 - o1Eo2*rEo2) / cc;
        xmu     = (-rEo2 + o1Eo2*rEo1) / cc;
        
        if ( std::abs(xlambda) > Node1->l_ch || std::abs(xmu) > Node2->l_ch )
        {
            aux1 = std::abs(xlambda) - Node1->l_ch;
            aux2 = std::abs(xmu)     - Node2->l_ch;
            
            if ( aux1 > aux2 )
            {
                xlambda = copysign(Node1->l_ch, xlambda);
                xmu     = xlambda*o1Eo2 - rEo2;
                
                if ( std::abs(xmu) > Node2->l_ch ) xmu = copysign(Node2->l_ch, xmu);
            }
            
            else
            {
                xmu     = copysign(Node2->l_ch, xmu);
                xlambda = xmu*o1Eo2 + rEo1;
                
                if ( std::abs(xlambda) > Node1->l_ch ) xlambda = copysign(Node1->l_ch, xlambda);
            }
        }
    }
    
    // Minimum line-to-line distance
    rm2 = R_sep.squaredNorm()  + SQR(xlambda)    + SQR(xmu)
        - 2.*xlambda*xmu*o1Eo2 - 2.*xlambda*rEo1 + 2.*xmu*rEo2;
    
    return (rm2 < rsqr);
}

// ============================
/* OBB overlap test */
// ============================
template<typename number>
bool Utils<number>::OverlapBoundOB(const Vector3<number>& R_cm, BNode<number>* Node1, BNode<number>* Node2)
{
    number r1;
    number r2;
    
    Vector3<number>  E1;
    Vector3<number>  E2;
    Vector3<number>  T;
    
    Matrix33<number> Rot;
    Matrix33<number> Abs_rot;
    
    Vector3<number>  R_sep = R_cm + Node2->Center - Node1->Center;
    
    /* ADAPTED FROM C. ERICSON, REAL-TIME COLLISION DETECTION */
    // Node half-dimensions
    E1 << Node1->l_xh, Node1->l_yh, Node1->l_zh;
    E2 << Node2->l_xh, Node2->l_yh, Node2->l_zh;
    
    // Rotation matrix projecting the frame of Node2 onto Node1
    Rot     = Node1->Orientation.transpose() * Node2->Orientation;
    Abs_rot = Rot.array().abs() + TOL_OB;
    
    // Project separation vector in the frame of Node1
    T = Node1->Orientation.transpose() * R_sep;
    
    // Test axes L = Node1.x, L = Node1.y, L = Node1.z
    for ( uint i = 0; i < 3; ++i )
    {
        r1 = E1(i);
        r2 = E2.dot(Abs_rot.row(i));
        
        if ( std::abs(T(i)) > r1 + r2 ) return false;
    }
    
    // Test axes L = Node2.x, L = Node2.y, L = Node2.z
    for ( uint i = 0; i < 3; ++i )
    {
        r1 = E1.dot(Abs_rot.col(i));
        r2 = E2(i);
        
        if ( std::abs(T.dot(Rot.col(i))) > r1 + r2 ) return false;
    }
    
    // Test axis L = (Node1.x).cross(Node2.x)
    r1 = E1(1)*Abs_rot(2,0) + E1(2)*Abs_rot(1,0);
    r2 = E2(1)*Abs_rot(0,2) + E2(2)*Abs_rot(0,1);
    
    if ( std::abs(T(2)*Rot(1,0) - T(1)*Rot(2,0)) > r1 + r2 ) return false;
    
    // Test axis L = (Node1.x).cross(Node2.y)
    r1 = E1(1)*Abs_rot(2,1) + E1(2)*Abs_rot(1,1);
    r2 = E2(0)*Abs_rot(0,2) + E2(2)*Abs_rot(0,0);
    
    if ( std::abs(T(2)*Rot(1,1) - T(1)*Rot(2,1)) > r1 + r2 ) return false;
    
    // Test axis L = (Node1.x).cross(Node2.z)
    r1 = E1(1)*Abs_rot(2,2) + E1(2)*Abs_rot(1,2);
    r2 = E2(0)*Abs_rot(0,1) + E2(1)*Abs_rot(0,0);
    
    if ( std::abs(T(2)*Rot(1,2) - T(1)*Rot(2,2)) > r1 + r2 ) return false;
    
    // Test axis L = (Node1.y).cross(Node2.x)
    r1 = E1(0)*Abs_rot(2,0) + E1(2)*Abs_rot(0,0);
    r2 = E2(1)*Abs_rot(1,2) + E2(2)*Abs_rot(1,1);
    
    if ( std::abs(T(0)*Rot(2,0) - T(2)*Rot(0,0)) > r1 + r2 ) return false;
    
    // Test axis L = (Node1.y).cross(Node2.y)
    r1 = E1(0)*Abs_rot(2,1) + E1(2)*Abs_rot(0,1);
    r2 = E2(0)*Abs_rot(1,2) + E2(2)*Abs_rot(1,0);
    
    if ( std::abs(T(0)*Rot(2,1) - T(2)*Rot(0,1)) > r1 + r2 ) return false;
    
    // Test axis L = (Node1.y).cross(Node2.z)
    r1 = E1(0)*Abs_rot(2,2) + E1(2)*Abs_rot(0,2);
    r2 = E2(0)*Abs_rot(1,1) + E2(1)*Abs_rot(1,0);
    
    if ( std::abs(T(0)*Rot(2,2) - T(2)*Rot(0,2)) > r1 + r2 ) return false;
    
    // Test axis L = (Node1.z).cross(Node2.x)
    r1 = E1(0)*Abs_rot(1,0) + E1(1)*Abs_rot(0,0);
    r2 = E2(1)*Abs_rot(2,2) + E2(2)*Abs_rot(2,1);
    
    if ( std::abs(T(1)*Rot(0,0) - T(0)*Rot(1,0)) > r1 + r2 ) return false;
    
    // Test axis L = (Node1.z).cross(Node2.y)
    r1 = E1(0)*Abs_rot(1,1) + E1(1)*Abs_rot(0,1);
    r2 = E2(0)*Abs_rot(2,2) + E2(2)*Abs_rot(2,0);
    
    if ( std::abs(T(1)*Rot(0,1) - T(0)*Rot(1,1)) > r1 + r2 ) return false;
    
    // Test axis L = (Node1.z).cross(Node2.z)
    r1 = E1(0)*Abs_rot(1,2) + E1(1)*Abs_rot(0,2);
    r2 = E2(0)*Abs_rot(2,1) + E2(1)*Abs_rot(2,0);
    
    if ( std::abs(T(1)*Rot(0,2) - T(0)*Rot(1,2)) > r1 + r2 ) return false;
    
    return true;
}

// ============================
/* Vertex Principal Component Analysis (PCA) by Singular Value Decomposition (SVD) */
// ============================
template<typename number>
Matrix33<number> Utils<number>::PCA(const Matrix3X<number>& Vertices_in)
{
    // Set vertex center of mass to the origin if needed
    Vector3<number>  Center_of_mass = Vertices_in.rowwise().mean();
    Matrix3X<number> Vertices_cm    = Vertices_in.colwise() - Center_of_mass;
    
    Eigen::JacobiSVD<Matrix3X<number> > SVD(Vertices_cm, Eigen::ComputeThinU);
    Matrix33<number> Rot = SVD.matrixU();
    
    // Swap axes so that Rot.y and Rot.z correspond to the directions of minimal and maximal spread, respectively
    Rot.col(0).swap(Rot.col(1));
    Rot.col(1).swap(Rot.col(2));
    
    // Enforce right-handedness of node frame
    if ( Rot.determinant() < 0. ) Rot.col(2) *= -1.;
    
    return Rot;
}

// ============================
/* Configuration file parser with Eigen format conversion */
// ============================
template<typename number>
void Utils<number>::Load(const std::string& filename, Matrix3X<number>* Vertices, ArrayX<number>* Charges, ArrayX<uint>* Types, ArrayX<uint>* Sizes)
{
    number coeff;
    
    bool   is_first_line(false);
    
    // Size counters
    uint   rows  (0);
    uint   cols  (0);
    uint   rows_t(0);
    
    std::string   line;
    std::ifstream input_file(filename);
    
    if ( !input_file.good() ) throw std::runtime_error("Couldn't open input file " + filename);
    
    // Store unkown number of configurations into std::vector containers
    std::vector<number> data_buffer;
    std::vector<uint>   size_buffer;
    
    while ( std::getline(input_file, line) )
    {
        if ( !line.empty() )
        {
            is_first_line = true;
            
            std::istringstream stream(line);
            rows += 1;
            
            while ( stream >> coeff )
            {
                if ( rows == 1 ) cols += 1;
                data_buffer.push_back(coeff);
            }
        }
        
        else
        {
            if ( is_first_line )
            {
                is_first_line = false;
                size_buffer.push_back(rows - rows_t);
                
                rows_t = rows;
            }
        }
    }
    
    input_file.close();
    
    if ( size_buffer.size() < 1 ) size_buffer.push_back(rows);
    
    // Map std::vector to Matrix3X<number> in column-major default order    
    Matrix5X<number> Data = Matrix5X<number>::Map(&data_buffer[0], cols, rows);
    
    *Vertices = Data.block(1, 0, 3, rows);
    
    *Charges  = Data.block(4, 0, 1, rows).transpose();
    *Types    = Data.block(0, 0, 1, rows).transpose().template cast<uint>();
    
    *Sizes    = ArrayX<uint>::Map(&size_buffer[0], size_buffer.size());
}

template struct Utils<float>;
template struct Utils<double>;
