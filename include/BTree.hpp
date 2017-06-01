#ifndef BTREE_HPP_
#define BTREE_HPP_

#include "Utils.hpp"


class BTree: public BNode
{
public:
    BTree();
    
    virtual ~BTree() {if ( nodes_alloced > 0 ) RecursiveDeallocate(this);}
    
    uint m;
    
    // Build counters
    uint nodes_built;
    uint leaves_built;
    
    uint vert_alloced;
    uint nodes_alloced;
    
    void Build(const Eigen::Matrix3Xd&, double, uint);
    
private:
    uint max_depth;

    // Memory storage to be assigned for the bounding volume hierarchy
    std::vector<BNode> Tree;

    void Allocate(uint);

    // Recursive allocators and constructors
	void RecursiveAllocate(BNode*);
	void RecursiveDeallocate(BNode*);
    void RecursiveBuild(BNode*, const Eigen::Matrix3Xd&, double);

    void BuildLeaf(BNode*, const Eigen::Matrix3Xd&);
};

#endif
