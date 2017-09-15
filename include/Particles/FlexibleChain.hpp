#ifndef FLEXIBLE_CHAIN_HPP_
#define FLEXIBLE_CHAIN_HPP_

#include "BaseParticle.hpp"


template<typename number>
class FlexibleChain: public BaseParticle<number>
{
public:
	FlexibleChain();
	
	// Load random configuration from BHierarchy->Forest
	void Parse(std::mt19937_64& rng_engine) override
	{
		uint idx_conf = rng_engine() % N_CONF_;
		this->Hull    = &(this->BVH).Forest[idx_conf];
	}
	
	void Build(int) override;
	
protected:
	number E_CUT_;
	number R_CUT_;
	number EPSILON_;
	
private:
	uint   N_CONF_;
};

#endif
