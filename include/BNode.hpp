#ifndef BNODE_HPP
#define BNODE_HPP

// ===================================================================
/**
 * Base bounding node structure for broad overlap tests. Includes both
 * an Oriented-Bounding Box (OBB) and a Sphero-Cylinder (SC) template.
 */
// ===================================================================
/*
 * BNode.hpp: Version 1.1
 * Created 22/07/2016 by Maxime Tortora
 */
// ===================================================================

#include "globals.hpp"


template<typename number>
struct BNode
{
    BNode()
    {
        is_root = false;
        is_leaf = false;
    }
    
    virtual ~BNode() {}
    
    // Pointers to children nodes
    BNode<number>* NodeI;
    BNode<number>* NodeS;
    
    // Tree locators
    bool   is_root;
    bool   is_leaf;
    
    uint   idx_depth;

    // OBB dimensions
    number l_xh;
    number l_yh;
    number l_zh;

    // SC dimensions
    number l_ch;
    number l_cr;

    // 3d properties - Center_p_ and Orientation_p_ are expressed in the parent frame
    Vector3<number> Axis;
    
    Vector3<number> Center;
    Vector3<number> Center_p;

    Matrix33<number> Orientation;
    Matrix33<number> Orientation_p;
    
    // Enclosed vertices - only allocated for leaf nodes
    Matrix3X<number>* Vertices;
    
    ArrayX<number>* Vcharges;
    ArrayX<uint>* Vtypes;

    BNode<number>& operator=(const BNode<number>& Node)
    {
        l_xh        = Node.l_xh;
        l_yh        = Node.l_yh;
        l_zh        = Node.l_zh;
        
        l_ch        = Node.l_ch;
        l_cr        = Node.l_cr;
        
        Axis        = Node.Axis;
        Center      = Node.Center;
        Orientation = Node.Orientation;

        return *this;
    }
    
    // Initialise child nodes
    inline void Split()
    {
        *NodeI           = *this;
        *NodeS           = *this;
        
        NodeI->idx_depth = idx_depth + 1;
        NodeS->idx_depth = idx_depth + 1;
    }
    
    // Align child frame with parent
    inline void RotateFrame(const BNode* Parent)
    {
        Center      = Parent->Orientation * Center_p;
        Orientation = Parent->Orientation * Orientation_p;
        
        if ( MODE_TREE == TREE_SC ) Axis = Orientation.col(2);
    }
};

#endif
