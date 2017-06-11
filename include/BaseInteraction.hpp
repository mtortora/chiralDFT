#ifndef BASE_INTERACTION_HPP
#define BASE_INTERACTION_HPP

// ===================================================================
/**
 * Header-only BaseInteraction class
 * Makes use of the Curiously Recurring Template Pattern (CRTP):
 * https://en.wikipedia.org/wiki/Curiously_recurring_template_pattern
 * --
 * Speeds calculations up by roughly 30% over virtual function calls.
 */
// ===================================================================
/*
 * BaseParticle.hpp: Version 1.1
 * Created 29/06/2016 by Maxime Tortora
 */
// ===================================================================

#include "Utils.hpp"


template<class Interaction, typename number>
struct BaseInteraction
{
    // ============================
    /* Recursive BNode collision detection and energy computations */
    // ============================
    void RecursiveInteraction(const Vector3<number>& R_cm, BNode<number>* Node1, BNode<number>* Node2,
                              number* energy, number cutoff)
    {
        bool overlap(true);

        // Early escape if cutoff energy is already reached
        if ( *energy > cutoff ) return;
        
        // Avoid overlap checks on root node, as these were already handled by the MC integrator
        if ( (!Node1->is_root) || (!Node2->is_root) )
        {
            if ( MODE_TREE == TREE_OB ) overlap = Utils<number>::OverlapBoundOB(R_cm, Node1, Node2);
            if ( MODE_TREE == TREE_SC ) overlap = Utils<number>::OverlapBoundSC(R_cm, Node1, Node2);
        }
        
        // Stop tree traversal if boxes are disjoint
        if ( !overlap ) return;
        
        // Work out leaf interaction energies
        if ( (Node1->is_leaf) && (Node2->is_leaf) )
        {
            *energy += LeafInteraction(R_cm, Node1, Node2, *energy, cutoff);
            
            return;
        }
        
        // Descend into the children of Node1
        if ( (Node2->is_leaf) || (!Node1->is_leaf && (Node1->l_zh > Node2->l_zh)) )
        {
            Node1->NodeI->RotateFrame(Node1);
            Node1->NodeS->RotateFrame(Node1);
            
            RecursiveInteraction(R_cm, Node1->NodeI, Node2, energy, cutoff);
            RecursiveInteraction(R_cm, Node1->NodeS, Node2, energy, cutoff);
        }
        
        // Descend into the children of Node2
        else
        {
            Node2->NodeI->RotateFrame(Node2);
            Node2->NodeS->RotateFrame(Node2);
            
            RecursiveInteraction(R_cm, Node1, Node2->NodeI, energy, cutoff);
            RecursiveInteraction(R_cm, Node1, Node2->NodeS, energy, cutoff);
        }
        
        return;
    }
	
    // ============================
    /* Work out leaf-leaf interaction energy */
    // ============================
    number LeafInteraction(const Vector3<number>& R_cm, BNode<number>* Node1, BNode<number>* Node2,
                           number energy, number cutoff)
    {
        number energy_(0.);
        
        Matrix3X<number> Voxel1 = Node1->Orientation * (*Node1->Vertices);
        Matrix3X<number> Voxel2 = Node2->Orientation * (*Node2->Vertices);
        
        for ( uint idx_vtx1 = 0; idx_vtx1 < Voxel1.cols() && energy_ <= cutoff-energy; ++idx_vtx1 )
        {
            for ( uint idx_vtx2 = 0; idx_vtx2 < Voxel2.cols() && energy_ <= cutoff-energy; ++idx_vtx2 )
            {
                Vector3<number> R_sep = R_cm + Voxel2.col(idx_vtx2) - Voxel1.col(idx_vtx1);
                number norm           = R_sep.norm();
                
                // Significantly faster than a direct virtual function call
                energy_ += static_cast<Interaction*>(this)->PairInteraction(norm);
            }
        }
        
        return energy_;
    }
};

#endif
