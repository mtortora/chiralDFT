// ===================================================================
/**
 * Implements robust overlap methods for convex bounding volumes
 */
// ===================================================================
/*
 * Utils.cpp: Version 1.0
 * Created 26/07/2016 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Utils.hpp"

using namespace Eigen;

    
// ============================
/* Bounding spherocylinders overlap test */
// ============================
bool Utils::OverlapBoundSC(const Vector3d& R_cm, BNode* Node1, BNode* Node2)
{
    double   xlambda;
    double   xmu;
    double   aux1;
    double   aux2;
    double   rm2;

    Vector3d R_sep;
    
    // Root spherocylinders centers are forced to the particle center of mass 
    if      ( Node1->is_root && Node2->is_root ) R_sep = R_cm;
    else if ( Node1->is_root )                   R_sep = R_cm + Node2->Center;
    else if ( Node2->is_root )                   R_sep = R_cm - Node1->Center;
    else                                         R_sep = R_cm + Node2->Center - Node1->Center;

    double rEo1  = Node1->Axis.dot(R_sep);
    double rEo2  = Node2->Axis.dot(R_sep);;
    double o1Eo2 = Node1->Axis.dot(Node2->Axis);
    
    double rsqr  = SQR(Node1->l_cr + Node2->l_cr);
    double cc    = 1. - SQR(o1Eo2);
    
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
bool Utils::OverlapBoundOB(const Vector3d& R_cm, BNode* Node1, BNode* Node2)
{
    double   r1;
    double   r2;
    
    Vector3d E1;
    Vector3d E2;
    Vector3d T;
    
    Matrix3d Rot;
    Matrix3d Abs_rot;
    
    Vector3d R_sep = R_cm + Node2->Center - Node1->Center;
    
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
Matrix3d Utils::PCA(const Matrix3Xd& Vertices_in)
{
    // Set vertex center of mass to the origin if needed
    Vector3d  Center_of_mass = Vertices_in.rowwise().mean();
    Matrix3Xd Vertices_cm    = Vertices_in.colwise() - Center_of_mass;
    
    JacobiSVD<Matrix3Xd> SVD(Vertices_cm, ComputeThinU);
    Matrix3d Rot = SVD.matrixU();
    
    // Swap axes so that Rot.y and Rot.z correspond to the directions of minimal and maximal spread, respectively
    Rot.col(0).swap(Rot.col(1));
    Rot.col(1).swap(Rot.col(2));
    
    // Enforce right-handedness of node frame
    if ( Rot.determinant() < 0. ) Rot.col(2) *= -1.;
    
    return Rot;
}

// ============================
/* Configuration file reader with Eigen format conversion */
// ============================
void Utils::Load(const std::string& filename, Matrix3Xd* Vertices, uint* particle_size)
{
    double coeff;
    
    bool   is_first_conf(true);
    
    // Size counters
    uint   rows(0);
    uint   cols(0);
    
    std::string   line;
    std::ifstream input_file(filename);
    
    if ( !input_file.good() ) throw std::runtime_error("Couldn't open input file " + filename);
    
    // Store unkown number of coefficients into a std::vector container
    std::vector<double> buffer;
    
    while ( std::getline(input_file, line) )
    {
        if ( !line.empty() )
        {
            std::istringstream stream(line);
            rows += 1;
            
            while ( stream >> coeff )
            {
                if ( rows == 1 ) cols += 1;
                buffer.push_back(coeff);
            }
        }
        
        else
        {
            if ( is_first_conf )
            {
                *particle_size = rows;
                is_first_conf  = false;
            }
        }
    }
    
    input_file.close();
    
    if ( is_first_conf ) *particle_size = rows;
    
    // Map std::vector to Eigen::Matrix3Xd in column-major default order
    *Vertices = Matrix3Xd::Map(&buffer[0], cols, rows);
}
