#ifndef FLEXIBLE_HELIX_HPP_
#define FLEXIBLE_HELIX_HPP_

#include "BaseParticle.hpp"


template<typename number>
class FlexibleHelix: public BaseParticle<number>
{
public:
	FlexibleHelix();
	
	// Load random configuration from BHierarchy->Forest
	void Parse(std::mt19937_64& rng_engine) override
	{
		uint idx_conf = rng_engine() % N_CONF_;
		this->Hull    = &(this->BVH).Forest[idx_conf];
	}
	
	void Build(int) override;
	
protected:
	number D_HARD_;
	
private:
	uint   N_CONF_;
};

#endif
