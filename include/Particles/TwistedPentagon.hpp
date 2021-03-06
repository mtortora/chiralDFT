#ifndef TWISTED_PENTAGON_HPP_
#define TWISTED_PENTAGON_HPP_

#include "rapid2/RAPID.H"
#include "BaseParticle.hpp"


template<typename number>
class TwistedPentagon: public BaseParticle<number>
{
public:
	TwistedPentagon();
	~TwistedPentagon() {delete Mesh;}
	
	RAPID_model* Mesh;
	
	void Build(int) override;
	void Tesselate(const Matrix3X<number>&, uint*);
	
#if (USE_RAPID)
	void Parse(std::mt19937_64&, ArrayX<uint>&) override {}
#endif
	
protected:
	number R_THRESHOLD_;
	
private:
	uint   N_R_;
	uint   N_Z_;
	
	number R_PNT_;
	number L_PNT_;
	number L_Z_;
	number TWIST_;
	number R_BCK_;
	number P_BCK_;
	
	void SaveWireframe(const Matrix3X<number>&);
	void SaveMesh     (const Matrix3X<number>&, uint);
};

#endif
