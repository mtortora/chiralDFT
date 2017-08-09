// ===================================================================
/**
 * Triangular prism derived particle class.
 */
// ===================================================================
/*
 * TwistedCuboid.cpp: Version 2.1
 * Created 31/07/2017 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/TriangularPrism.hpp"


template<typename number>
TriangularPrism<number>::TriangularPrism()
{
	this->N_DELTA_L = 2;
	
	// Prism parameters
	GAMMA_          = 60. * PI/180.;
	TWIST_          = 0.  * PI/180.;
	
	L_X_            = 1.   * this->SIGMA_R;
	
	number w        = L_X_/PI * (1. + 1./cos(GAMMA_));

	L_Z_            = 6.   * w;
	L_Y_            = L_X_ * tan(GAMMA_)/2.;
	
	this->V0        = L_X_*L_Y_*L_Z_/2.;
	this->V_EFF     = this->V0;
	
	this->R_INTEG   = sqrt(SQR(L_X_) + SQR(L_Z_) + SQR(2.*L_Y_/3.));
	this->V_INTEG   = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
	
	// Allocate RAPID mesh
	Mesh            = new RAPID_model;
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void TriangularPrism<number>::Build(int mpi_rank)
{
	Matrix3X<number> Wireframe(3, 6);
	
	Vector3<number>  Z_axis = Vector3<number>::UnitZ();

	Wireframe.col(0) << -L_X_/2., -L_Y_/3., -L_Z_/2.;
	Wireframe.col(1) <<  L_X_/2., -L_Y_/3., -L_Z_/2.;
	Wireframe.col(2) <<  0.,    2.*L_Y_/3., -L_Z_/2.;

	Wireframe.col(3) << -L_X_/2., -L_Y_/3.,  L_Z_/2.;
	Wireframe.col(4) <<  L_X_/2., -L_Y_/3.,  L_Z_/2.;
	Wireframe.col(5) <<  0.,    2.*L_Y_/3.,  L_Z_/2.;

	// Rotate upper face by angle TWIST_ around e_z
	Matrix3X<number> Face_s = Matrix3X<number>::Map(Wireframe.data()+9, 3, 3);
	Face_s = Eigen::AngleAxis<number>(TWIST_, Z_axis).toRotationMatrix() * Face_s;
	
	Wireframe.block(0, 3, 3, 3)    = Face_s;
	
	// Set center of mass to the origin and main axis to e_z
	Vector3<number> Center_of_mass = Wireframe.rowwise().mean();
	Wireframe                      = Wireframe.colwise() - Center_of_mass;
	
	Matrix33<number> Rot           = Utils<number>::PCA(Wireframe);
	Wireframe                      = Rot.transpose() * Wireframe;
	
	// Build the RAPID_model mesh for RAPID collision detection
	Tesselate(Wireframe);

	if ( mpi_rank == MPI_MASTER ) SaveMesh(Wireframe);

	this->Hull       = &this->BVH;
	
	// Bounding volume parameters
	this->Hull->l_xh = L_X_/2.;
	this->Hull->l_yh = std::max(L_Y_/2., L_X_*sin(TWIST_)/2.);
	this->Hull->l_zh = L_Z_/2.;
	
	this->Hull->l_ch = this->Hull->l_zh;
	this->Hull->l_cr = L_X_/2.;
}

// ============================
/* Wireframe triangulation for RAPID interference tests */
// ============================
template<typename number>
void TriangularPrism<number>::Tesselate(const Matrix3X<number>& Wireframe)
{
	uint   ctr_tri(0);
	
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];
	double p5[3];
	double p6[3];
	
	Mesh->BeginModel();
	
	// Allocate vertices
	p1[0] = Wireframe(0, 0);
	p1[1] = Wireframe(1, 0);
	p1[2] = Wireframe(2, 0);
	
	p2[0] = Wireframe(0, 1);
	p2[1] = Wireframe(1, 1);
	p2[2] = Wireframe(2, 1);
	
	p3[0] = Wireframe(0, 2);
	p3[1] = Wireframe(1, 2);
	p3[2] = Wireframe(2, 2);
	
	p4[0] = Wireframe(0, 3);
	p4[1] = Wireframe(1, 3);
	p4[2] = Wireframe(2, 3);
	
	p5[0] = Wireframe(0, 4);
	p5[1] = Wireframe(1, 4);
	p5[2] = Wireframe(2, 4);
	
	p6[0] = Wireframe(0, 5);
	p6[1] = Wireframe(1, 5);
	p6[2] = Wireframe(2, 5);
	
	// Build mesh
	Mesh->AddTri(p1, p2, p3, ctr_tri++);
	Mesh->AddTri(p1, p2, p4, ctr_tri++);
	Mesh->AddTri(p1, p3, p4, ctr_tri++);
	Mesh->AddTri(p2, p3, p5, ctr_tri++);
	Mesh->AddTri(p2, p4, p5, ctr_tri++);
	Mesh->AddTri(p3, p4, p6, ctr_tri++);
	Mesh->AddTri(p3, p5, p6, ctr_tri++);
	Mesh->AddTri(p4, p5, p6, ctr_tri++);
	
	Mesh->EndModel();
	
	if ( this->id_ == 1 ) LogTxt("Running with %d-triangle tesselated mesh", ctr_tri);
}

// ============================
/* Save mesh to PLY file */
// ============================
template<typename number>
void TriangularPrism<number>::SaveMesh(const Matrix3X<number>& Wireframe)
{
	std::string   DATA_PATH = __DATA_PATH;
	std::ofstream file_mesh(DATA_PATH + "/mesh.ply");
	
	file_mesh << "ply" << std::endl;
	file_mesh << "format ascii 1.0" << std::endl;
	
	// Vertex and face formatting
	file_mesh << "element vertex "  << 6 << std::endl;
	file_mesh << "property float x" << std::endl;
	file_mesh << "property float y" << std::endl;
	file_mesh << "property float z" << std::endl;
	
	file_mesh << "element face "  << 8 << std::endl;
	file_mesh << "property list uchar int vertex_index" << std::endl;
	
	file_mesh << "end_header" << std::endl;
	
	file_mesh << Wireframe.transpose() << std::endl;
	
	file_mesh << 3 << " " << 0 << " " << 1 << " " << 2 << std::endl;
	file_mesh << 3 << " " << 0 << " " << 1 << " " << 3 << std::endl;
	file_mesh << 3 << " " << 0 << " " << 2 << " " << 3 << std::endl;
	file_mesh << 3 << " " << 1 << " " << 2 << " " << 4 << std::endl;
	file_mesh << 3 << " " << 1 << " " << 3 << " " << 4 << std::endl;
	file_mesh << 3 << " " << 2 << " " << 3 << " " << 5 << std::endl;
	file_mesh << 3 << " " << 2 << " " << 4 << " " << 5 << std::endl;
	file_mesh << 3 << " " << 3 << " " << 4 << " " << 5 << std::endl;

	
	file_mesh.close();
}

template class TriangularPrism<float>;
template class TriangularPrism<double>;
