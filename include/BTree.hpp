#ifndef BTREE_HPP_
#define BTREE_HPP_

#include "Utils.hpp"


template<typename number>
class BTree: public BNode<number>
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
    
    void Build(const Matrix3X<number>&, const ArrayX<number>&, const ArrayX<uint>&, number, uint);
    
private:
    uint max_depth;

    // Memory storage to be assigned for the bounding volume hierarchy
    std::vector<BNode<number> > Tree;

    void Allocate(uint);

    // Recursive allocators and constructors
	void RecursiveAllocate(BNode<number>*);
	void RecursiveDeallocate(BNode<number>*);
    void RecursiveBuild(BNode<number>*, const Matrix3X<number>&, const ArrayX<number>&, const ArrayX<uint>&, number);

    void BuildLeaf(BNode<number>*, const Matrix3X<number>&, const ArrayX<number>&, const ArrayX<uint>&);
};

#endif
