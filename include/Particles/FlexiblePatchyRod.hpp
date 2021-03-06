#ifndef FLEXIBLE_PATCHY_ROD_HPP_
#define FLEXIBLE_PATCHY_ROD_HPP_

#include "BaseParticle.hpp"


template<typename number>
class FlexiblePatchyRod: public BaseParticle<number>
{
public:
    FlexiblePatchyRod();
    
    Matrix3X<number> Backbones;
    Matrix3X<number> Patches;
	
	uint N_BCK;
	uint idx_conf;

	// Build bootstrap index map
	ArrayX<uint> BootstrapMap(std::mt19937_64& rng_engine) override
	{
		return ArrayX<uint>::NullaryExpr(N_CONF_, [&]() {return rng_engine() % N_CONF_;});
	}
	
	// Load random configuration from BHierarchy->Forest
	void Parse(std::mt19937_64& rng_engine, ArrayX<uint>& BMap) override
	{
		idx_conf = rng_engine() % N_CONF_;
		idx_conf = BMap[idx_conf];
		
		this->Hull    = &(this->BVH).Forest[idx_conf];
	}
	
	void Build(int) override;
	
protected:
    number EPSILON_WCA_;
    number E_CUT_;
	
    number D_BACK_;
	number R_BACK_;
	
	number D_PATCH_;
	number R_PATCH_;
	
	number D_LB_;
	number R_LB_;

private:
    uint   N_CONF_;
};

#endif
