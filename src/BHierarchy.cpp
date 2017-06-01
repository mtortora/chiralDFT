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

using namespace Eigen;


BHierarchy::BHierarchy()
{
	m             = 3;
	trees_alloced = 0;
	
	total_built   = 0;
	total_alloced = 0;
	
	total_vert    = 0;
	total_leaves  = 0;
}

// ============================
/* Build single tree hierarchy */
// ============================
void BHierarchy::Build(const Matrix3Xd& Vertices_, double range)
{
	if ( Forest.size() == 0 ) Forest.resize(1);
	
	BTree* Tree_   = &Forest[trees_alloced++];
	
	Tree_->Build(Vertices_, range, m);
	
	total_built   += Tree_->nodes_built;
	total_alloced += Tree_->nodes_alloced;

	total_vert    += Tree_->vert_alloced;
	total_leaves  += Tree_->leaves_built;
}

// ============================
/* Overloaded forest constructor */
// ============================
void BHierarchy::Build(const Matrix3Xd& Vertices_trees, double range, const ArrayXi& Size_trees)
{
	uint num_trees = Size_trees.size();
	
	if ( num_trees > 0 ) Forest.resize(num_trees);
	else throw std::runtime_error("Found empty size array for forest allocation");

	for ( uint idx_tree = 0; idx_tree < num_trees; ++idx_tree )
	{
		uint      n_vert    = Size_trees(idx_tree);
		Matrix3Xd Vertices_ = Matrix3Xd::Map(Vertices_trees.data() + 3 * total_vert, 3, n_vert);
		
		Build(Vertices_, range);
	}
}
