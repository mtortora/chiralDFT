#ifndef TWISTED_PENTAGON_HPP_
#define TWISTED_PENTAGON_HPP_

#include "rapid2/RAPID.H"
#include "BaseParticle.hpp"


class TwistedPentagon: public BaseParticle
{
private:
	uint   N_R_;
	uint   N_Z_;
	
	double R_PNT_;
	double L_PNT_;
	double L_Z_;
	double TWIST_;
	double R_BCK_;
	double P_BCK_;
	
	void SaveWireframe(const Eigen::Matrix3Xd&);
	void SaveMesh     (const Eigen::Matrix3Xd&, uint);
	
protected:
	double R_THRESHOLD_;
	
public:
	TwistedPentagon();
	
	RAPID_model* Mesh;
	
	void Build(int) override;
	void Tesselate(const Eigen::Matrix3Xd&, uint*);
	
#if (USE_RAPID)
	void Parse(std::mt19937_64&) override {}
#endif
	
	virtual ~TwistedPentagon() {delete Mesh;}
};

#endif
