#ifndef FD_GROMOS_HPP_
#define FD_GROMOS_HPP_

#include "BaseParticle.hpp"


template<typename number>
class FDGromos: public BaseParticle<number>
{
public:
	FDGromos();
	
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
	
	number DH_PREFACTOR_;
	number MINUS_KAPPA_;
	number C_RF_;
	
	number C6[12][12];
	number C12[12][12];
	
private:
	uint   N_CONF_;
};

#endif
