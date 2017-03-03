#ifndef BNODE_HPP
#define BNODE_HPP

// ===================================================================
/**
 * Base bounding node structure for broad overlap tests. Includes both
 * an Oriented-Bounding Box (OBB) and a Sphero-Cylinder (SC) template
 */
// ===================================================================
/*
 * BNode.hpp:  Version 1.1
 * Created 22/07/2016 by Maxime Tortora
 */
// ===================================================================

#include "globals.hpp"


struct BNode
{
    BNode()
    {
        is_root = false;
        is_leaf = false;
    }
    
    // Pointers to children nodes
    BNode* NodeI;
    BNode* NodeS;
    
    // Tree locators
    bool   is_root;
    bool   is_leaf;
    
    uint   idx_depth;

    // OBB dimensions
    double l_xh;
    double l_yh;
    double l_zh;

    // SC dimensions
    double l_ch;
    double l_rc;

    // 3d properties - Center_p_ and Orientation_p_ are expressed in the parent frame
    Eigen::Vector3d Axis;
    
    Eigen::Vector3d Center;
    Eigen::Vector3d Center_p;

    Eigen::Matrix3d Orientation;
    Eigen::Matrix3d Orientation_p;
    
    // Enclosed vertices - only allocated for leaf nodes
    Eigen::Matrix3Xd* Vertices;

    // Overloaded assignment operator for tree traversal
    BNode& operator=(const BNode& Node)
    {
        l_xh        = Node.l_xh;
        l_yh        = Node.l_yh;
        l_zh        = Node.l_zh;
        
        l_ch        = Node.l_ch;
        l_rc        = Node.l_rc;
        
        Center      = Node.Center;
        Orientation = Node.Orientation;

        return *this;
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
