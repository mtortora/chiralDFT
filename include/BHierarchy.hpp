#ifndef BHIERARCHY_HPP_
#define BHIERARCHY_HPP_

#include "BTree.hpp"


template<typename number>
class BHierarchy: public BTree<number>
{
public:
	BHierarchy();
	
	// Memory storage to be assigned for the bounding volume hierarchy
	std::vector<BTree<number> > Forest;
	
	// Hierarchy constructors
	void Build(const Matrix3X<number>&, number);
	void Build(const Matrix3X<number>&, number, const ArrayX<uint>&);
	
	// Performance can be optimised by fine-tuning m
	inline void SetLeafParameter(uint m_=3) {this->m = fmax(3, m_);}
	
	// Build log
	inline void PrintBuildInfo()
	{
		LogTxt("Tree hierarchies successfully built: %d", trees_alloced);
		LogTxt("Maximum children nodes allocated: %d", total_alloced);
		LogGre("Built %d total bounding structures, inc. %d leaf nodes", total_built, total_leaves);
		LogGre("Total number of vertices enclosed: %d", total_vert);
	}
	
private:
	// Forest properties
	uint trees_alloced;
	
	uint total_built;
	uint total_alloced;
	uint total_leaves;
	uint total_vert;
};

#endif
