// ===================================================================
/**
 * Twisted pentagon derived particle class.
 */
// ===================================================================
/*
 * TwistedPentagon.cpp: Version 2.1
 * Created 22/02/2017 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/TwistedPentagon.hpp"


template<typename number>
TwistedPentagon<number>::TwistedPentagon()
{
	// Bounding leaf parameter
	this->BVH.SetLeafParameter(5);
	
	this->N_DELTA_L = 2;
	
	// Twisted pentagon parameters
	N_R_            = 10;
	N_Z_            = 2000;
	
	R_PNT_          = 3.3  * this->SIGMA_R;
	L_PNT_          = R_PNT_ * 2.*sin(PI/5.);
	
	// Helical backbone parameters
	number pitch    = 3.3  * this->SIGMA_R;
	number l_ctr    = 880. * this->SIGMA_R;

	R_BCK_          = 146. * this->SIGMA_R;
	P_BCK_          = 2800. * this->SIGMA_R;
	
	L_Z_            = l_ctr / sqrt(1. + SQR(2.*PI * R_BCK_/P_BCK_));
	TWIST_          = 2.*PI/5. * L_Z_/pitch;
	
	R_THRESHOLD_    = sqrt(SQR(L_PNT_/(N_R_-1.)) + SQR(L_Z_/(N_Z_-1.)))/2.;
	
	this->R_INTEG   = sqrt(SQR(2.*R_PNT_+2.*R_BCK_) + SQR(L_Z_)) + R_THRESHOLD_;
	this->V_INTEG   = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
	
	this->V0        = SQR(L_PNT_)*L_Z_ * 5./4. * tan(54. * PI/180.);
	this->V_EFF     = this->V0;
	
	// Allocate RAPID mesh
	Mesh            = new RAPID_model;
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void TwistedPentagon<number>::Build(int mpi_rank)
{
	number l_h = L_PNT_ / tan(36. * PI/180.);
	
	Matrix3X<number> Edge(3, N_R_);
	
	Matrix3X<number> Backbone (3, N_Z_);
	Matrix3X<number> Wireframe(3, 5*N_R_*N_Z_);
	Matrix3X<number> Face_r   (3, N_R_*N_Z_);
	
	Vector3<number>  Z_axis     = Vector3<number>::UnitZ();
	
	ArrayX<number>   R_grid     = ArrayX<number>::LinSpaced(N_R_, -L_PNT_/2., L_PNT_/2.);
	ArrayX<number>   Z_grid     = ArrayX<number>::LinSpaced(N_Z_,  0.,      L_Z_);
	
	ArrayX<number>   Alpha_grid = ArrayX<number>::LinSpaced(N_Z_,  0.,      TWIST_);
	
	// Generate helical backbone
	Backbone.row(0)             = R_BCK_ * cos(2.*PI/P_BCK_ * Z_grid);
	Backbone.row(1)             = R_BCK_ * sin(2.*PI/P_BCK_ * Z_grid);
	Backbone.row(2)             = Z_grid;
	
	// Generate first edges
	Edge.row(0)                 = R_grid;
	Edge.row(1)                 = VectorX<number>::Constant(N_R_, -l_h/2.);
	Edge.row(2)                 = VectorX<number>::Constant(N_R_, -L_Z_/2.);
	
	// Generate transversal faces
	for ( uint idx_f = 0; idx_f < 5; ++idx_f )
	{
		for ( uint idx_z = 0; idx_z < N_Z_; ++idx_z )
		{
			number alpha_     = Alpha_grid(idx_z);
			
			Matrix3X<number> R_edge = Eigen::AngleAxis<number>(alpha_, Z_axis).toRotationMatrix() * Edge;
			R_edge.colwise() += Backbone.col(idx_z);
			
			Face_r.block(0, idx_z * N_R_, 3, N_R_) = R_edge;
		}
		
		Wireframe.block(0, idx_f * N_R_*N_Z_, 3, N_R_*N_Z_) = Face_r;
		
		Edge = Eigen::AngleAxis<number>(2.*PI/5., Z_axis).toRotationMatrix() * Edge;
	}
	
	// Set center of mass to the origin and main axis to e_z
	Vector3<number> Center_of_mass = Wireframe.rowwise().mean();
	Wireframe                      = Wireframe.colwise() - Center_of_mass;
	
	Matrix33<number> Rot           = Utils<number>::PCA(Wireframe);
	Wireframe                      = Rot.transpose() * Wireframe;
	
	// Save to file on master thread
	if ( mpi_rank == MPI_MASTER ) SaveWireframe(Wireframe);
	
	// Build the RAPID_model mesh for RAPID collision detection
	if ( USE_RAPID )
	{
		uint num_tri;
		
		this->Hull       = &this->BVH;
		
		// Bounding volume parameters
		this->Hull->l_xh = R_BCK_ + R_PNT_;
		this->Hull->l_yh = R_BCK_ + R_PNT_;
		this->Hull->l_zh = L_Z_ / 2.;
		
		this->Hull->l_ch = this->Hull->l_zh;
		this->Hull->l_cr = R_BCK_ + R_PNT_;
		
		Tesselate(Wireframe, &num_tri);
		
		if ( mpi_rank == MPI_MASTER ) SaveMesh(Wireframe, num_tri);
	}
	
	// Build bounding volume hierarchy
	else
	{
		// Build bounding volume hierarchy
		this->BVH.Build(Wireframe, R_THRESHOLD_);
		
		if ( this->id_ == 1 ) this->BVH.PrintBuildInfo();
	}
}


// ============================
/* Wireframe triangulation for RAPID interference tests */
// ============================
template<typename number>
void TwistedPentagon<number>::Tesselate(const Matrix3X<number>& Wireframe, uint* num_tri)
{
	uint   ctr_tri(0);
	
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];
	
	Mesh->BeginModel();
	
	for ( uint idx_f = 0; idx_f < 5; ++idx_f )
	{
		for ( uint idx_z = 0; idx_z < N_Z_-1; ++idx_z )
		{
			for ( uint idx_r = 0; idx_r < N_R_-1; ++idx_r )
			{
				p1[0] = Wireframe(0, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r);
				p1[1] = Wireframe(1, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r);
				p1[2] = Wireframe(2, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r);
				
				p2[0] = Wireframe(0, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + 1);
				p2[1] = Wireframe(1, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + 1);
				p2[2] = Wireframe(2, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + 1);
				
				p3[0] = Wireframe(0, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + N_R_);
				p3[1] = Wireframe(1, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + N_R_);
				p3[2] = Wireframe(2, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + N_R_);
				
				p4[0] = Wireframe(0, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + N_R_ + 1);
				p4[1] = Wireframe(1, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + N_R_ + 1);
				p4[2] = Wireframe(2, idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + N_R_ + 1);
				
				Mesh->AddTri(p1, p2, p3, ctr_tri);
				Mesh->AddTri(p4, p2, p3, ctr_tri+1);
				
				ctr_tri += 2;
			}
		}
	}
	
	Mesh->EndModel();
	
	if ( this->id_ == 1 ) LogTxt("Running with %d-triangle tesselated mesh", ctr_tri);
	
	*num_tri = ctr_tri;
}

// ============================
/* Save wireframe to gnuplot-readable file */
// ============================
template<typename number>
void TwistedPentagon<number>::SaveWireframe(const Matrix3X<number>& Wireframe)
{
	std::string   DATA_PATH = __DATA_PATH;
	std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");
	
	for ( uint idx_f = 0; idx_f < 5; ++idx_f )
	{
		for ( uint idx_z = 0; idx_z < N_Z_; ++idx_z )
		{
			for ( uint idx_r = 0; idx_r < N_R_; ++idx_r )
			{
				uint idx = idx_f*N_R_*N_Z_ + idx_z*N_R_ + idx_r;
				Vector3<number> Vertex = Wireframe.col(idx);
				
				file_wireframe << Vertex.adjoint() << std::endl;
			}
			
			file_wireframe << std::endl;
		}
		
		file_wireframe << std::endl;
	}
	
	file_wireframe.close();
}

// ============================
/* Save mesh to PLY file */
// ============================
template<typename number>
void TwistedPentagon<number>::SaveMesh(const Matrix3X<number>& Wireframe, uint num_tri)
{
	uint ctr_tri(0);
	
	std::string   DATA_PATH = __DATA_PATH;
	std::ofstream file_mesh(DATA_PATH + "/mesh.ply");
	
	file_mesh << "ply" << std::endl;
	file_mesh << "format ascii 1.0" << std::endl;
	
	// Vertex and face formatting
	file_mesh << "element vertex "  << 2*num_tri << std::endl;
	file_mesh << "property float x" << std::endl;
	file_mesh << "property float y" << std::endl;
	file_mesh << "property float z" << std::endl;

	file_mesh << "element face "  << num_tri << std::endl;
	file_mesh << "property list uchar int vertex_index" << std::endl;
	
	file_mesh << "end_header" << std::endl;
	
	// Build vertex list
	for ( uint idx_f = 0; idx_f < 5; ++idx_f )
	{
		for ( uint idx_z = 0; idx_z < N_Z_-1; ++idx_z )
		{
			for ( uint idx_r = 0; idx_r < N_R_-1; ++idx_r )
			{
				Vector3<number> Vtx1 = Wireframe.col(idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r);
				Vector3<number> Vtx2 = Wireframe.col(idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + 1);
				Vector3<number> Vtx3 = Wireframe.col(idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + N_R_);
				Vector3<number> Vtx4 = Wireframe.col(idx_f * N_R_*N_Z_ + idx_z*N_R_ + idx_r + N_R_ + 1);
				
				file_mesh << Vtx1(0) << " " << Vtx1(1) << " " << Vtx1(2) << std::endl;
				file_mesh << Vtx2(0) << " " << Vtx2(1) << " " << Vtx2(2) << std::endl;
				file_mesh << Vtx3(0) << " " << Vtx3(1) << " " << Vtx3(2) << std::endl;
				file_mesh << Vtx4(0) << " " << Vtx4(1) << " " << Vtx4(2) << std::endl;
			}
		}
	}
	
	// Build face list
	while ( ctr_tri < num_tri )
	{
		file_mesh << 3 << " " << 2*ctr_tri   << " " << 2*ctr_tri+1 << " " << 2*ctr_tri+2 << std::endl;
		file_mesh << 3 << " " << 2*ctr_tri+3 << " " << 2*ctr_tri+1 << " " << 2*ctr_tri+2 << std::endl;

		ctr_tri += 2;
	}
	
	file_mesh.close();
}

template class TwistedPentagon<float>;
template class TwistedPentagon<double>;
