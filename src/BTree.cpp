// ===================================================================
/**
 * Specialised root node for a bounding volume tree.
 * Children nodes are allocated in the std::vector Tree.
 */
// ===================================================================
/*
 * BTree.cpp: Version 1.1
 * Created 20/07/2016 by Maxime Tortora
 */
// ===================================================================

#include "BTree.hpp"


template<typename number>
BTree<number>::BTree()
{
    // Root node parameters
    this->is_root       = true;
    
    this->idx_depth     = 0;
    this->max_depth     = 0;
    
    this->nodes_built   = 0;
    this->leaves_built  = 0;
    
    this->vert_alloced  = 0;
    this->nodes_alloced = 0;
    
    // Set the principal axes of all initial configurations parallel to the reference frame
    this->Axis          = Vector3<number>::UnitZ();
    
    this->Center        = Vector3<number>::Zero();
    this->Center_p      = Vector3<number>::Zero();
    
    this->Orientation   = Matrix33<number>::Identity();
    this->Orientation_p = Matrix33<number>::Identity();
}

// ============================
/* Recursive memory allocation */
// ============================
template<typename number>
void BTree<number>::RecursiveAllocate(BNode<number>* Node)
{
    if ( (Node == this) && (max_depth == 0) ) nodes_alloced++;
    else
    {
        Node->NodeI = &Tree[nodes_alloced++];
        Node->NodeS = &Tree[nodes_alloced++];

        Node->Split();

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
template<typename number>
void BTree<number>::RecursiveDeallocate(BNode<number>* Node)
{
    if ( Node->is_leaf )
    {
        delete Node->Vertices;
        delete Node->Vtypes;
    }
    
    else
    {
        if ( Node->NodeI->is_leaf )
        {
            delete Node->NodeI->Vertices;
            delete Node->NodeI->Vtypes;
        }
        
        else RecursiveDeallocate(Node->NodeI);

        if ( Node->NodeS->is_leaf )
        {
            delete Node->NodeS->Vertices;
            delete Node->NodeS->Vtypes;
        }
        
        else RecursiveDeallocate(Node->NodeS);
    }
}

// ============================
/* Single-tree allocator */
// ============================
template<typename number>
void BTree<number>::Allocate(uint num_vert)
{    
    // Binary tree size parameters
    uint tree_depth = floor(log((float)num_vert/(float)m) / log(2.) + 0.5);
    uint n_leaves   = pow(2, tree_depth);
    
    uint n_nodes    = 2 * (n_leaves-1);
    max_depth       = tree_depth;
    
    // Node constructor is called by the resize() method
    if ( max_depth > 0 ) Tree.resize(n_nodes);
    
    RecursiveAllocate(this);
}

// ============================
/* Single-tree constructor */
// ============================
template<typename number>
void BTree<number>::Build(const Matrix3X<number>& Vertices_, const ArrayX<number>& Charges_, const ArrayX<uint>& Types_, number range, uint m_)
{
    this->m       = m_;
    uint num_vert = Vertices_.cols();

    Allocate(num_vert);
    RecursiveBuild(this, Vertices_, Charges_, Types_, range);
}

// ============================
/*  Specialised leaf constructor */
// ============================
template<typename number>
void BTree<number>::BuildLeaf(BNode<number>* Node, const Matrix3X<number>& Vertices_, const ArrayX<number>& Charges_, const ArrayX<uint>& Types_)
{
    Node->is_leaf  = true;
    
    Node->Vertices = new(std::nothrow) Matrix3X<number>;
    
    Node->Vcharges = new(std::nothrow) ArrayX<number>;
    Node->Vtypes = new(std::nothrow) ArrayX<uint>;
    
    if ( !Node->Vertices ) throw std::runtime_error("Vertex memory allocation failed");
    if ( !Node->Vtypes ) throw std::runtime_error("Type memory allocation failed");

    *Node->Vertices = Vertices_;
    
    *Node->Vcharges = Charges_;
    *Node->Vtypes = Types_;
    
    vert_alloced   += Vertices_.cols();
    
    leaves_built++;
}

// ============================
/* Build bounding volume hierarchy */
// ============================
template<typename number>
void BTree<number>::RecursiveBuild(BNode<number>* Node, const Matrix3X<number>& Vertices_in, const ArrayX<number>& Charges_in, const ArrayX<uint>& Types_in, number range)
{
    // Vertices are expressed in the parent frame
    Vector3<number>  Center_of_mass = Vertices_in.rowwise().mean();
    Matrix3X<number> Vertices_cm    = Vertices_in.colwise() - Center_of_mass;
    
    // Principal Component Analysis (PCA) fails for 3xn matrices if n < 3 so set transformation to identity
    if ( Vertices_in.cols() < 3 ) Node->Orientation_p = Matrix33<number>::Identity();
    else
    {
        // Orientation_p_ is expressed in the parent frame
        Node->Orientation_p = Utils<number>::PCA(Vertices_cm);
        
        // Project Vertices_cm in child frame
        Vertices_cm         = Node->Orientation_p.transpose() * Vertices_cm;
    }
    
    // Work out bounding box dimensions
    ArrayX<number> Vertex_x = Vertices_cm.row(0).array();
    ArrayX<number> Vertex_y = Vertices_cm.row(1).array();
    ArrayX<number> Vertex_z = Vertices_cm.row(2).array();

    number  x_min           = Vertex_x.minCoeff();
    number  y_min           = Vertex_y.minCoeff();
    number  z_min           = Vertex_z.minCoeff();

    number  x_max           = Vertex_x.maxCoeff();
    number  y_max           = Vertex_y.maxCoeff();
    number  z_max           = Vertex_z.maxCoeff();

    Node->l_xh              = (x_max-x_min + range) / 2.;
    Node->l_yh              = (y_max-y_min + range) / 2.;
    Node->l_zh              = (z_max-z_min + range) / 2.;

    // Force root spherocylinder center to the vertex center of mass
    if ( !Node->is_root )
    {
        Vertices_cm.row(0).array() -= (x_max+x_min) / 2.;
        Vertices_cm.row(1).array() -= (y_max+y_min) / 2.;
    }

    // Bounding spherocylinder dimensions
    number r_max   = Vertices_cm.block(0, 0, 2, Vertices_cm.cols()).colwise().norm().maxCoeff();
    
    Node->l_cr     = r_max + range/2.;
    Node->l_ch     = Node->is_root ? fmax(std::abs(z_min),std::abs(z_max)) : (z_max-z_min)/2.;

    // Box center expressed in parent frame
    Node->Center_p << (x_max+x_min)/2., (y_max+y_min)/2., (z_max+z_min)/2.;
    Node->Center_p =  Node->Orientation_p * Node->Center_p + Center_of_mass;
    
    // Vertices expressed in child frame
    Matrix3X<number> Vertices_new = Node->Orientation_p.transpose() * Vertices_in;

    nodes_built++;

    // Leaf nodes are allocated if maximum depth is reached, or if the number of enclosed vertices is < m
    if ( (Node->idx_depth == max_depth) || (Vertices_in.cols() < m) )
    {
        BuildLeaf(Node, Vertices_new, Charges_in, Types_in);
        return;
    }

    // Split box in 2 along the axis of maximal dispersion
    ArrayX<bool> Vertex_inf = (Vertex_z <= (z_max+z_min)/2.);
    ArrayX<bool> Vertex_sup = (Vertex_z >  (z_max+z_min)/2.);
    
    uint num_inf = Vertex_inf.cast<uint>().sum();
    uint num_sup = Vertex_sup.cast<uint>().sum();
    
    uint ctr_inf = 0;
    uint ctr_sup = 0;

    Matrix3X<number> Vertices_inf(3, num_inf);
    Matrix3X<number> Vertices_sup(3, num_sup);

    ArrayX<number> Charges_inf(num_inf);
    ArrayX<number> Charges_sup(num_sup);
    
    ArrayX<uint> Types_inf(num_inf);
    ArrayX<uint> Types_sup(num_sup);
    
    for ( uint idx_vtx = 0; idx_vtx < Vertices_in.cols(); ++idx_vtx )
    {
        if ( Vertex_inf(idx_vtx) )
        {
            Vertices_inf.col(ctr_inf) = Vertices_new.col(idx_vtx);
            
            Charges_inf(ctr_inf) = Charges_in(idx_vtx);
            Types_inf(ctr_inf) = Types_in(idx_vtx);

            ++ctr_inf;
        }
        
        if ( Vertex_sup(idx_vtx) )
        {
            Vertices_sup.col(ctr_sup) = Vertices_new.col(idx_vtx);
            
            Charges_sup(ctr_sup) = Charges_in(idx_vtx);
            Types_sup(ctr_sup) = Types_in(idx_vtx);

            ++ctr_sup;
        }
    }
    
    RecursiveBuild(Node->NodeI, Vertices_inf, Charges_inf, Types_inf, range);
    RecursiveBuild(Node->NodeS, Vertices_sup, Charges_sup, Types_sup, range);
    
    return;
}

template class BTree<float>;
template class BTree<double>;
