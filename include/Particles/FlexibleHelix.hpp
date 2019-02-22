#ifndef FLEXIBLE_HELIX_HPP_
#define FLEXIBLE_HELIX_HPP_

#include "BaseParticle.hpp"


template<typename number>
class FlexibleHelix: public BaseParticle<number>
{
public:
	FlexibleHelix();
	
	// Build bootstrap index map
	ArrayX<uint> BootstrapMap(std::mt19937_64& rng_engine) override
	{
		return ArrayX<uint>::NullaryExpr(N_CONF_, [&]() {return rng_engine() % N_CONF_;});
	}
	
	// Load random configuration from BHierarchy->Forest
	void Parse(std::mt19937_64& rng_engine, ArrayX<uint>& BMap) override
	{
		uint idx_conf = rng_engine() % N_CONF_;
		uint idx_boot = BMap[idx_conf];
		
		this->Hull    = &(this->BVH).Forest[idx_boot];
	}
	
	void Build(int) override;
	
protected:
	number D_HARD_;
	
private:
	uint   N_CONF_;
};

#endif
