#ifndef TWISTED_HEXAGON_HPP_
#define TWISTED_HEXAGON_HPP_

#include "rapid2/RAPID.H"
#include "BaseParticle.hpp"


template<typename number>
class TwistedHexagon: public BaseParticle<number>
{
public:
	TwistedHexagon();
	~TwistedHexagon() {delete Mesh;}
	
	RAPID_model* Mesh;
	
	void Build(int) override;
	void Tesselate(const Matrix3X<number>&, uint*);
	
#if (USE_RAPID)
	void Parse(std::mt19937_64&) override {}
#endif
	
protected:
	number R_THRESHOLD_;
	
private:
	uint   N_R_;
	uint   N_Z_;
	
	number R_HEX_;
	number L_Z_;
	number TWIST_;
	number R_BCK_;
	number P_BCK_;
	
	void SaveWireframe(const Matrix3X<number>&);
	void SaveMesh     (const Matrix3X<number>&, uint);
};

#endif
