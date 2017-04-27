#ifndef BTREE_HPP_
#define BTREE_HPP_

#include "Utils.hpp"


struct BTree: public BNode
{
    BTree();
    ~BTree();

    // Tree properties
    uint n_nodes;
    uint max_depth;
    uint n_vert_max;
    
    // Build counters
    uint nodes_built;
    uint leaves_built;
    
    uint vert_alloced;
    uint nodes_alloced;
    uint trees_alloced;

    // Memory storage to be assigned for the bounding volume hierarchy
    std::vector<BNode> Tree;
    std::vector<BTree> Forest;

    // Recursive constructors
	void RecursiveAllocate(BNode*);
	void RecursiveDeallocate(BNode*);
    void RecursiveBuild(BNode*, const Eigen::Matrix3Xd&, double);

    // Allocators
    void AllocateForest(uint);
    void AllocateTree();
    void AllocateLeaf(BNode*, const Eigen::Matrix3Xd&);

    // Overloaded assignment operator for forest traversal
    BTree& operator=(const BTree& Tree_)
    {
        max_depth   = Tree_.max_depth;
        n_nodes     = Tree_.n_nodes;
        n_vert_max  = Tree_.n_vert_max;
        
        return *this;
    }
    
    // Performance can be optimised by fine-tuning max_depth and n_vert_max
    inline void SetTreeProperties(uint max_depth_, uint n_vert_max_=3)
    {
        uint n_leaves = pow(2, max_depth_);
        
        n_nodes       = 2 * (n_leaves-1);
        n_vert_max    = fmax(3, n_vert_max_);
        
        max_depth     = max_depth_;
    }
    
    // Build log
    inline void PrintBuildInfo()
    {
        LogTxt("Root nodes successfully allocated: %d", trees_alloced);
        LogTxt("Maximum children nodes allocated: %d", nodes_alloced);
        LogGre("Built %d total bounding structures, inc. %d leaf nodes", nodes_built, leaves_built);
        LogGre("Total number of vertices enclosed: %d", vert_alloced);
    }
};

#endif
