// ===================================================================
/**
 * Specialised root node for the bounding volume hierarchy tree.
 * Forest properties only need to be set for flexible particles
 */
// ===================================================================
/*
 * BTree.cpp: Version 1.1
 * Created 20/07/2016 by Maxime Tortora
 */
// ===================================================================

#include "BTree.hpp"

using namespace Eigen;


BTree::BTree()
{
    is_root       = true;
    
    idx_depth     = 0;
    max_depth     = 0;
    
    nodes_built   = 0;
    leaves_built  = 0;
    
    nodes_alloced = 0;
    trees_alloced = 0;
    
    Axis          = Vector3d::UnitZ();
    
    Center        = Vector3d::Zero();
    Center_p      = Vector3d::Zero();
    
    Orientation   = Matrix3d::Identity();
    Orientation_p = Matrix3d::Identity();
}

// ============================
/* Recursive memory allocation */
// ============================
void BTree::RecursiveAllocate(BNode* Node)
{
    if ( (Node == this) && (max_depth == 0) ) nodes_alloced++;
    else
    {
        Node->NodeI = &Tree[nodes_alloced++];
        Node->NodeS = &Tree[nodes_alloced++];

        Split(Node);

        if ( (Node->NodeI->idx_depth < max_depth) && (Node->NodeS->idx_depth < max_depth) )
        {
            RecursiveAllocate(Node->NodeI);
            RecursiveAllocate(Node->NodeS);
        }
    }
}

// ============================
/* Recursive memory deallocation */
// ============================
void BTree::RecursiveDeallocate(BNode* Node)
{
    if ( Node->is_leaf ) delete Node->Vertices;
    else
    {
        if ( Node->NodeI->is_leaf ) delete Node->NodeI->Vertices;
        else RecursiveDeallocate(Node->NodeI);

        if ( Node->NodeS->is_leaf ) delete Node->NodeS->Vertices;
        else RecursiveDeallocate(Node->NodeS);
    }
}

// ============================
// Vertex memory allocation & assignment
// ============================
void BTree::AllocateLeaf(BNode* Node, const Eigen::Matrix3Xd& Vertices_)
{
    Node->is_leaf  = true;
    Node->Vertices = new(std::nothrow) Eigen::Matrix3Xd;
    
    if ( !Node->Vertices ) throw std::runtime_error("Vertex memory allocation failed");
    
    *Node->Vertices = Vertices_;
    
    leaves_built++;
}

// ============================
/* Single tree allocator */
// ============================
void BTree::AllocateTree()
{        
    // Allocate children node memory in std::vector Tree
    if ( max_depth > 0 ) Tree.resize(n_nodes);
    
    RecursiveAllocate(this);
    
    trees_alloced++;
}

// ============================
/* Allocates a forest of size num_trees */
// ============================
void BTree::AllocateForest(uint num_trees)
{
    // Initialise root nodes via the default constructor, called in the resize() method
    if ( num_trees > 0 ) Forest.resize(num_trees);
    
    for ( uint idx_tree = 0; idx_tree < num_trees; ++idx_tree )
    {
        BTree* Tree_ = &Forest[idx_tree];
        *Tree_       = *this;
        
        Tree_->AllocateTree();

        nodes_alloced += Tree_->nodes_alloced;
        trees_alloced += Tree_->trees_alloced;
    }
}

// ============================
/* Build bounding volume hierarchy */
// ============================
void BTree::RecursiveBuild(BNode* Node, const Matrix3Xd& Vertices_in, double range)
{
    // Vertices are expressed in the parent frame
    Vector3d  Center_of_mass = Vertices_in.rowwise().mean();
    Matrix3Xd Vertices_cm    = Vertices_in.colwise() - Center_of_mass;
    
    // Principal Component Analysis (PCA) fails for 3xn matrices if n < 3 so set transformation to identity
    if ( Vertices_in.cols() < 3 ) Node->Orientation_p = Matrix3d::Identity();
    else
    {
        // Orientation_p_ is expressed in the parent frame
        Node->Orientation_p = Utils::PCA(Vertices_cm);
        
        // Project Vertices_cm in child frame
        Vertices_cm         = Node->Orientation_p.transpose() * Vertices_cm;
    }
    
    // Work out bounding box dimensions
    ArrayXd Vertex_x = Vertices_cm.row(0).array();
    ArrayXd Vertex_y = Vertices_cm.row(1).array();
    ArrayXd Vertex_z = Vertices_cm.row(2).array();

    double  x_min    = Vertex_x.minCoeff();
    double  y_min    = Vertex_y.minCoeff();
    double  z_min    = Vertex_z.minCoeff();

    double  x_max    = Vertex_x.maxCoeff();
    double  y_max    = Vertex_y.maxCoeff();
    double  z_max    = Vertex_z.maxCoeff();

    Node->l_xh       = (x_max-x_min + range) / 2.;
    Node->l_yh       = (y_max-y_min + range) / 2.;
    Node->l_zh       = (z_max-z_min + range) / 2.;

    // Force root spherocylinder center to the vertex center of mass
    if ( !Node->is_root )
    {
        Vertices_cm.row(0).array() -= (x_max+x_min) / 2.;
        Vertices_cm.row(1).array() -= (y_max+y_min) / 2.;
    }

    // Bounding spherocylinder dimensions
    double r_max          = Vertices_cm.block(0, 0, 2, Vertices_cm.cols()).colwise().norm().maxCoeff();
    
    Node->l_rc            = r_max + range/2.;
    Node->l_ch            = Node->is_root ? fmax(std::abs(z_min),std::abs(z_max)) : (z_max-z_min) / 2.;

    // Box center expressed in parent frame
    Node->Center_p       << (x_max+x_min)/2., (y_max+y_min)/2., (z_max+z_min)/2.;
    Node->Center_p        = Node->Orientation_p * Node->Center_p + Center_of_mass;
    
    // Vertices expressed in child frame
    MatrixXd Vertices_new = Node->Orientation_p.transpose() * Vertices_in;

    nodes_built++;

    // Leaf nodes are allocated if maximum depth is reached, or if the number of enclosed vertices is < n_vert_max
    if ( (Node->idx_depth == max_depth) || (Vertices_in.cols() < n_vert_max) )
    {
        AllocateLeaf(Node, Vertices_new);
        return;
    }

    // Split box in 2 along the axis of maximal dispersion
    ArrayXb Vertex_inf = (Vertex_z <= (z_max+z_min)/2.);
    ArrayXb Vertex_sup = (Vertex_z >  (z_max+z_min)/2.);
    
    uint    num_inf    = Vertex_inf.cast<uint>().sum();
    uint    num_sup    = Vertex_sup.cast<uint>().sum();
    
    uint    ctr_inf    = 0;
    uint    ctr_sup    = 0;

    Matrix3Xd Vertices_inf(3, num_inf);
    Matrix3Xd Vertices_sup(3, num_sup);

    for ( uint idx_vtx = 0; idx_vtx < Vertices_in.cols(); ++idx_vtx )
    {
        if ( Vertex_inf(idx_vtx) ) Vertices_inf.col(ctr_inf++) = Vertices_new.col(idx_vtx);
        if ( Vertex_sup(idx_vtx) ) Vertices_sup.col(ctr_sup++) = Vertices_new.col(idx_vtx);
    }
    
    RecursiveBuild(Node->NodeI, Vertices_inf, range);
    RecursiveBuild(Node->NodeS, Vertices_sup, range);
    
    return;
}

// Forest nodes are allocated solely through std::vector's and do not need to be freed manually
BTree::~BTree() {if ( (nodes_alloced > 0) && (Forest.size() == 0) ) RecursiveDeallocate(this);}
