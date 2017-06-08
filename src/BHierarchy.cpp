// ===================================================================
/**
 * Bounding volume hierarchy class.
 */
// ===================================================================
/*
 * BHierarchy.cpp: Version 1.0
 * Created 31/06/2017 by Maxime Tortora
 */
// ===================================================================


#include "BHierarchy.hpp"


template<typename number>
BHierarchy<number>::BHierarchy()
{
	this->m       = 3;
	trees_alloced = 0;
	
	total_built   = 0;
	total_alloced = 0;
	
	total_vert    = 0;
	total_leaves  = 0;
}

// ============================
/* Build single tree hierarchy */
// ============================
template<typename number>
void BHierarchy<number>::Build(const Matrix3X<number>& Vertices_, number range)
{
	if ( Forest.size() == 0 ) Forest.resize(1);
	
	BTree<number>* Tree_ = &Forest[trees_alloced++];
	
	Tree_->Build(Vertices_, range, this->m);
	
	total_built   += Tree_->nodes_built;
	total_alloced += Tree_->nodes_alloced;

	total_vert    += Tree_->vert_alloced;
	total_leaves  += Tree_->leaves_built;
}

// ============================
/* Overloaded forest constructor */
// ============================
template<typename number>
void BHierarchy<number>::Build(const Matrix3X<number>& Vertices_trees, number range, const ArrayX<uint>& Size_trees)
{
	uint num_trees = Size_trees.size();
	
	if ( num_trees > 0 ) Forest.resize(num_trees);
	else throw std::runtime_error("Found empty size array for forest allocation");

	for ( uint idx_tree = 0; idx_tree < num_trees; ++idx_tree )
	{
		uint n_vert = Size_trees(idx_tree);
		Matrix3X<number> Vertices_ = Matrix3X<number>::Map(Vertices_trees.data() + 3 * total_vert, 3, n_vert);
		
		Build(Vertices_, range);
	}
}

template class BHierarchy<float>;
template class BHierarchy<double>;
