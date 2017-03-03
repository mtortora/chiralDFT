// ===================================================================
/**
 * Twisted pentagon derived particle class
 */
// ===================================================================
/*
 * TwistedPentagon.cpp: Version 2.1
 * Created 22/02/2017 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/TwistedPentagon.hpp"

using namespace Eigen;


TwistedPentagon::TwistedPentagon()
{
	// Bounding tree properties
	BHierarchy->SetTreeProperties(14, 5);
	
	N_DELTA_L    = 2;
	
	// Twisted pentagon parameters
	N_R_         = 10;
	N_Z_         = 2000;
	
	R_PNT_       = 3.3  * SIGMA_R;
	L_PNT_       = R_PNT_ * 2.*sin(PI/5.);

	L_Z_         = 880. * SIGMA_R;
	
	double pitch = 3.3  * SIGMA_R;
	TWIST_       = 2.*PI/5. * L_Z_/pitch;

	// Helical backbone parameters
	R_BCK_       = 0.  * SIGMA_R;
	P_BCK_       = 20. * SIGMA_R;
	
	R_THRESHOLD_ = sqrt(SQR(L_PNT_/(N_R_-1.)) + SQR(L_Z_/(N_Z_-1.)))/2.;
	R_INTEG      = sqrt(SQR(2.*R_PNT_+2.*R_BCK_) + SQR(L_Z_)) + R_THRESHOLD_;
	V_INTEG      = CUB(2.*R_INTEG) * 16.*pow(PI, 6);
	V0           = SQR(L_PNT_)*L_Z_ * 5./4. * tan(54. * PI/180.);
	V_EFF        = V0;
	
	Mesh         = new RAPID_model;
}

// ============================
/* Build particle model */
// ============================
void TwistedPentagon::Build(int mpi_rank)
{
	double    l_h = L_PNT_ / tan(36. * PI/180.);
	
	Matrix3Xd Edge(3, N_R_);
	
	Matrix3Xd Backbone (3, N_Z_);
	Matrix3Xd Wireframe(3, 5*N_R_*N_Z_);
	Matrix3Xd Face_r   (3, N_R_*N_Z_);
	
	Vector3d  Z_axis     = Vector3d::UnitZ();
	
	ArrayXd   R_grid     = VectorXd::LinSpaced(N_R_, -L_PNT_/2., L_PNT_/2.);
	ArrayXd   Z_grid     = VectorXd::LinSpaced(N_Z_,  0.,      L_Z_);
	
	ArrayXd   Alpha_grid = VectorXd::LinSpaced(N_Z_,  0.,      TWIST_);
	
	// Generate helical backbone
	Backbone.row(0) = R_BCK_ * cos(2.*PI/P_BCK_ * Z_grid);
	Backbone.row(1) = R_BCK_ * sin(2.*PI/P_BCK_ * Z_grid);
	Backbone.row(2) = Z_grid;
	
	// Generate first edges
	Edge.row(0) = R_grid;
	Edge.row(1) = VectorXd::Constant(N_R_, -l_h/2.);
	Edge.row(2) = VectorXd::Constant(N_R_, -L_Z_/2.);
	
	// Generate transversal faces
	for ( uint idx_f = 0; idx_f < 5; ++idx_f )
	{
		for ( uint idx_z = 0; idx_z < N_Z_; ++idx_z )
		{
			double alpha      = Alpha_grid(idx_z);
			
			Matrix3Xd R_edge  = AngleAxisd(alpha, Z_axis).toRotationMatrix() * Edge;
			R_edge.colwise() += Backbone.col(idx_z);
			
			Face_r.block(0, idx_z * N_R_, 3, N_R_) = R_edge;
		}
		
		Wireframe.block(0, idx_f * N_R_*N_Z_, 3, N_R_*N_Z_) = Face_r;
		
		Edge = AngleAxisd(2.*PI/5., Z_axis).toRotationMatrix() * Edge;
	}
	
	// Set center of mass to the origin and main axis to e_z
	Vector3d Center_of_mass = Wireframe.rowwise().mean();
	Wireframe               = Wireframe.colwise() - Center_of_mass;
	
	Matrix3d Rot            = Utils::PCA(Wireframe);
	Wireframe               = Rot.transpose() * Wireframe;
	
	// Save to file on master thread
	if ( mpi_rank == MPI_MASTER ) SaveWireframe(Wireframe);
	
	// Build the RAPID_model mesh for RAPID collision detection
	if ( USE_RAPID )
	{
		BHull       = BHierarchy;
		
		// Bounding volume parameters
		BHull->l_xh = R_BCK_ + R_PNT_;
		BHull->l_yh = R_BCK_ + R_PNT_;
		BHull->l_zh = L_Z_ / 2.;
		
		BHull->l_ch = BHull->l_zh;
		BHull->l_rc = R_BCK_ + R_PNT_;
		
		Tesselate(Wireframe, mpi_rank);
	}
	
	// Build bounding volume hierarchy
	else
	{
		BHierarchy->AllocateForest(1);
		
		BTree* Tree = &BHierarchy->Forest[0];
		BHierarchy->RecursiveBuild(Tree, Wireframe, R_THRESHOLD_);
		
		if ( id_ == 1 ) BHierarchy->PrintBuildInfo();
	}
}


// ============================
/* Wireframe triangulation for RAPID interference tests */
// ============================
void TwistedPentagon::Tesselate(const Matrix3Xd& Wireframe, uint mpi_rank)
{
	uint   ctr_tri(0);
	
	double p1[3];
	double p2[3];
	double p3[3];
	double p4[3];
	
	std::string DATA_PATH;
	
	// Redirect slave thread output to /dev/null
	if ( mpi_rank == MPI_MASTER ) DATA_PATH = __DATA_PATH;
	else                          DATA_PATH = "/dev/null";
	
	// Tesselated mesh file can be displayed by the resources/utils/display_mesh script
	std::ofstream file_mesh(DATA_PATH + "/mesh.out");
	
	Mesh->BeginModel();
	
	// Tesselate transversal faces
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
				
				file_mesh << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
				file_mesh << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
				file_mesh << p3[0] << " " << p3[1] << " " << p3[2] << std::endl;
				file_mesh << p1[0] << " " << p1[1] << " " << p1[2] << std::endl;
				file_mesh << std::endl << std::endl;
				
				file_mesh << p4[0] << " " << p4[1] << " " << p4[2] << std::endl;
				file_mesh << p2[0] << " " << p2[1] << " " << p2[2] << std::endl;
				file_mesh << p3[0] << " " << p3[1] << " " << p3[2] << std::endl;
				file_mesh << p4[0] << " " << p4[1] << " " << p4[2] << std::endl;
				file_mesh << std::endl << std::endl;
				
				Mesh->AddTri(p1, p2, p3, ctr_tri);
				Mesh->AddTri(p4, p2, p3, ctr_tri+1);
				
				ctr_tri += 2;
			}
		}
	}
	
	Mesh->EndModel();
	file_mesh.close();
	
	if ( id_ == 1 ) LogTxt("Running with %d-triangle tesselated mesh", ctr_tri);
}

// ============================
/* Save wireframe to gnuplot-readable file */
// ============================
void TwistedPentagon::SaveWireframe(const Matrix3Xd& Wireframe)
{
	std::string   DATA_PATH = __DATA_PATH;
	std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");
	
	// Save transversal faces
	for ( uint idx_f = 0; idx_f < 5; ++idx_f )
	{
		for ( uint idx_z = 0; idx_z < N_Z_; ++idx_z )
		{
			for ( uint idx_r = 0; idx_r < N_R_; ++idx_r )
			{
				uint idx        = idx_f*N_R_*N_Z_ + idx_z*N_R_ + idx_r;
				Vector3d Vertex = Wireframe.col(idx);
				
				file_wireframe << Vertex.adjoint() << std::endl;
			}
			
			file_wireframe << std::endl;
		}
		
		file_wireframe << std::endl;
	}
	
	file_wireframe.close();
}
