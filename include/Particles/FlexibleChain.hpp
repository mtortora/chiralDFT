#ifndef FLEXIBLE_CHAIN_HPP_
#define FLEXIBLE_CHAIN_HPP_

#include "BaseParticle.hpp"


template<typename number>
class FlexibleChain: public BaseParticle<number>
{
public:
	FlexibleChain();
	
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
	number E_CUT_;
	number R_CUT_;
	
	number RC_WCA_;
	number RC_DH_;
	
	number DH_PREFACTOR_;
	number MINUS_KAPPA_;

	number TYS_;
	number EPSILON_;
	
private:
	uint   N_CONF_;
};

#endif
