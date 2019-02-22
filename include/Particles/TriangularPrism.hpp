#ifndef TRIANGULAR_PRISM_HPP_
#define TRIANGULAR_PRISM_HPP_

#include "rapid2/RAPID.H"
#include "BaseParticle.hpp"


template<typename number>
class TriangularPrism: public BaseParticle<number>
{
public:
	TriangularPrism();
	~TriangularPrism() {delete Mesh;}
	
	RAPID_model* Mesh;
	
	void Build(int) override;
	void Tesselate(const Matrix3X<number>&);
	
	void Parse(std::mt19937_64&, ArrayX<uint>&) override {}
	
private:
	number L_X_;
	number L_Y_;
	number L_Z_;
	number TWIST_;
	number GAMMA_;
	
	void SaveMesh(const Matrix3X<number>&);
};

#endif
