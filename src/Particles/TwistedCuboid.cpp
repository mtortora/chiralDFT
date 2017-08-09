// ===================================================================
/**
 * Twisted cuboid derived particle class.
 */
// ===================================================================
/*
 * TwistedCuboid.cpp: Version 2.1
 * Created 20/09/2015 by Maxime Tortora
 */
// ===================================================================

#include <fstream>

#include "Particles/TwistedCuboid.hpp"


template<typename number>
TwistedCuboid<number>::TwistedCuboid()
{
    // Bounding leaf parameter
    this->BVH.SetLeafParameter(4);

    this->N_DELTA_L = 2;

    // Cuboid parameters
    N_X_            = 10;
    N_Y_            = 10;
    N_Z_            = 1000;

    L_X_            = 1.   * this->SIGMA_R;
    L_Y_            = 1.   * this->SIGMA_R;
    L_Z_            = 100. * this->SIGMA_R;
    
    number alpha_   = 0.5;
    
    // Helical backbone parameters
    R_BCK_          = 0.  * this->SIGMA_R;
    P_BCK_          = 40. * this->SIGMA_R;
    
    // Thread angle in radians
    number mu       = 80. * PI/180.;
    TWIST_          = 2. * L_Z_ / (sqrt(SQR(L_X_)+SQR(L_Y_)) * tan(mu));

    this->V0        = L_X_*L_Y_*L_Z_;
    this->V_EFF     = this->V0;
    
    R_THRESHOLD_    = (1.+alpha_) * fmin(L_X_/((number)N_X_), fmin(L_Y_/((number)N_Y_), L_Z_/((number)N_Z_)));

    #if (!USE_RAPID)
    // Rescale cuboid dimensions to account for finite HS radii
    L_X_           *= (1. - (1.+alpha_)/((number)N_X_));
    L_Y_           *= (1. - (1.+alpha_)/((number)N_Y_));
    L_Z_           *= (1. - (1.+alpha_)/((number)N_Z_));
    #endif
    
    this->R_INTEG   = sqrt(SQR(L_X_+R_BCK_) + SQR(L_Y_+R_BCK_) + SQR(L_Z_)) + R_THRESHOLD_;
    this->V_INTEG   = CUB(2.*this->R_INTEG) * 16.*pow(PI, 6);
    
    // Allocate RAPID mesh
    Mesh            = new RAPID_model;
}

// ============================
/* Build particle model */
// ============================
template<typename number>
void TwistedCuboid<number>::Build(int mpi_rank)
{
    Matrix3X<number> Face_xz(3, N_X_*N_Z_), Edge_x(3, N_X_);
    Matrix3X<number> Face_yz(3, N_Y_*N_Z_), Edge_y(3, N_Y_);
    Matrix3X<number> Face_xy(3, N_X_*N_Y_);
    
    Matrix3X<number> Backbone (3, N_Z_);
    Matrix3X<number> Wireframe(3, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + 2*N_X_*N_Y_);

    Vector3<number>  Z_axis     = Vector3<number>::UnitZ();
    
    ArrayX<number>   X_grid     = ArrayX<number>::LinSpaced(N_X_, -L_X_/2., L_X_/2.);
    ArrayX<number>   Y_grid     = ArrayX<number>::LinSpaced(N_Y_, -L_Y_/2., L_Y_/2.);
    ArrayX<number>   Z_grid     = ArrayX<number>::LinSpaced(N_Z_,  0.,      L_Z_);
    ArrayX<number>   Alpha_grid = ArrayX<number>::LinSpaced(N_Z_,  0.,      TWIST_);

    // Generate helical backbone
    Backbone.row(0)             = R_BCK_ * cos(2.*PI/P_BCK_ * Z_grid);
    Backbone.row(1)             = R_BCK_ * sin(2.*PI/P_BCK_ * Z_grid);
    Backbone.row(2)             = Z_grid;
    
    // Generate first edges
    Edge_x.row(0)               = X_grid;
    Edge_x.row(1)               = VectorX<number>::Constant(N_X_, -L_Y_/2.);
    Edge_x.row(2)               = VectorX<number>::Constant(N_X_, -L_Z_/2.);

    Edge_y.row(0)               = VectorX<number>::Constant(N_Y_, -L_X_/2.);
    Edge_y.row(1)               = Y_grid;
    Edge_y.row(2)               = VectorX<number>::Constant(N_Y_, -L_Z_/2.);

    // Generate transversal faces
    for ( uint idx_f = 0; idx_f < 2; ++idx_f )
    {
        Edge_x = Eigen::AngleAxis<number>(idx_f*PI, Z_axis).toRotationMatrix() * Edge_x;
        Edge_y = Eigen::AngleAxis<number>(idx_f*PI, Z_axis).toRotationMatrix() * Edge_y;

        for ( uint idx_z = 0; idx_z < N_Z_; ++idx_z )
        {
            number alpha_             = Alpha_grid(idx_z);
            
            Matrix3X<number> R_edge_x = Eigen::AngleAxis<number>(alpha_, Z_axis).toRotationMatrix() * Edge_x;
            Matrix3X<number> R_edge_y = Eigen::AngleAxis<number>(alpha_, Z_axis).toRotationMatrix() * Edge_y;

            R_edge_x.colwise() += Backbone.col(idx_z);
            R_edge_y.colwise() += Backbone.col(idx_z);
            
            Face_xz.block(0, idx_z * N_X_, 3, N_X_) = R_edge_x;
            Face_yz.block(0, idx_z * N_Y_, 3, N_Y_) = R_edge_y;
        }
        
        Wireframe.block(0, idx_f * N_X_*N_Z_, 3, N_X_*N_Z_)               = Face_xz;
        Wireframe.block(0, 2*N_X_*N_Z_ + idx_f * N_Y_*N_Z_, 3, N_Y_*N_Z_) = Face_yz;
    }
    
    // Generate normal faces
    for ( uint idx_f = 0; idx_f < 2; ++idx_f )
    {
        for ( uint idx_x = 0; idx_x < N_X_; ++idx_x )
        {
            for ( uint idx_y = 0; idx_y < N_Y_; ++idx_y )
            {
                uint idx        = idx_x*N_Y_ + idx_y;
                
                Face_xy(0, idx) = X_grid(idx_x);
                Face_xy(1, idx) = Y_grid(idx_y);
                Face_xy(2, idx) = -L_Z_/2.;
            }
        }
        
        Face_xy = Eigen::AngleAxis<number>(idx_f*TWIST_, Z_axis).toRotationMatrix() * Face_xy;
        Face_xy.colwise() += Backbone.col(idx_f * (N_Z_-1));
        
        Wireframe.block(0, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_f * N_X_*N_Y_, 3, N_X_*N_Y_) = Face_xy;
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
        this->Hull->l_xh = R_BCK_ + sqrt(SQR(L_X_) + SQR(L_Y_)) / 2.;
        this->Hull->l_yh = R_BCK_ + sqrt(SQR(L_X_) + SQR(L_Y_)) / 2.;
        this->Hull->l_zh = L_Z_ / 2.;
        
        this->Hull->l_ch = this->Hull->l_zh;
        this->Hull->l_cr = R_BCK_ + sqrt(SQR(L_X_) + SQR(L_Y_)) / 2.;
        
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
void TwistedCuboid<number>::Tesselate(const Matrix3X<number>& Wireframe, uint* num_tri)
{
    uint   ctr_tri(0);

    double p1[3];
    double p2[3];
    double p3[3];
    double p4[3];
    
    Mesh->BeginModel();

    // Tesselate (XZ) faces
    for ( uint idx_fy = 0; idx_fy < 2; ++idx_fy )
    {
        for ( uint idx_z = 0; idx_z < N_Z_-1; ++idx_z )
        {
            for ( uint idx_x = 0; idx_x < N_X_-1; ++idx_x )
            {
                p1[0] = Wireframe(0, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x);
                p1[1] = Wireframe(1, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x);
                p1[2] = Wireframe(2, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x);
                
                p2[0] = Wireframe(0, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + 1);
                p2[1] = Wireframe(1, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + 1);
                p2[2] = Wireframe(2, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + 1);
            
                p3[0] = Wireframe(0, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + N_X_);
                p3[1] = Wireframe(1, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + N_X_);
                p3[2] = Wireframe(2, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + N_X_);
                
                p4[0] = Wireframe(0, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + N_X_ + 1);
                p4[1] = Wireframe(1, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + N_X_ + 1);
                p4[2] = Wireframe(2, idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + N_X_ + 1);
                
                Mesh->AddTri(p1, p2, p3, ctr_tri);
                Mesh->AddTri(p4, p2, p3, ctr_tri+1);
                
                ctr_tri += 2;
            }
        }
    }

    // Tesselate (YZ) faces
    for ( uint idx_fx = 0; idx_fx < 2; ++idx_fx )
    {
        for ( uint idx_z = 0; idx_z < N_Z_-1; ++idx_z )
        {
            for ( uint idx_y = 0; idx_y < N_Y_-1; ++idx_y )
            {
                p1[0] = Wireframe(0, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y);
                p1[1] = Wireframe(1, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y);
                p1[2] = Wireframe(2, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y);
                
                p2[0] = Wireframe(0, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + 1);
                p2[1] = Wireframe(1, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + 1);
                p2[2] = Wireframe(2, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + 1);
                
                p3[0] = Wireframe(0, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + N_Y_);
                p3[1] = Wireframe(1, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + N_Y_);
                p3[2] = Wireframe(2, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + N_Y_);
                
                p4[0] = Wireframe(0, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + N_Y_ + 1);
                p4[1] = Wireframe(1, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + N_Y_ + 1);
                p4[2] = Wireframe(2, 2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + N_Y_ + 1);
                
                Mesh->AddTri(p1, p2, p3, ctr_tri);
                Mesh->AddTri(p4, p2, p3, ctr_tri+1);
                
                ctr_tri += 2;
            }
        }
    }

    // Tesselate (XY) faces
    for ( uint idx_fz = 0; idx_fz < 2; ++idx_fz )
    {
        for ( uint idx_x = 0; idx_x < N_X_-1; ++idx_x )
        {
            for ( uint idx_y = 0; idx_y < N_Y_-1; ++idx_y )
            {
                p1[0] = Wireframe(0, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y);
                p1[1] = Wireframe(1, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y);
                p1[2] = Wireframe(2, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y);
                
                p2[0] = Wireframe(0, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + 1);
                p2[1] = Wireframe(1, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + 1);
                p2[2] = Wireframe(2, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + 1);
                
                p3[0] = Wireframe(0, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + N_Y_);
                p3[1] = Wireframe(1, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + N_Y_);
                p3[2] = Wireframe(2, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + N_Y_);
                
                p4[0] = Wireframe(0, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + N_Y_ + 1);
                p4[1] = Wireframe(1, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + N_Y_ + 1);
                p4[2] = Wireframe(2, 2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + N_Y_ + 1);
                
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
void TwistedCuboid<number>::SaveWireframe(const Matrix3X<number>& Wireframe)
{
    std::string   DATA_PATH = __DATA_PATH;
    std::ofstream file_wireframe(DATA_PATH + "/wireframe.out");
    
    // Save (XZ) faces
    for ( uint idx_fy = 0; idx_fy < 2; ++idx_fy )
    {
        for ( uint idx_z = 0; idx_z < N_Z_; ++idx_z )
        {
            for ( uint idx_x = 0; idx_x < N_X_; ++idx_x )
            {
                uint idx = idx_fy*N_X_*N_Z_ + idx_z*N_X_ + idx_x;
                Vector3<number> Vertex = Wireframe.col(idx);
                
                file_wireframe << Vertex.adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    // Save (YZ) faces
    for ( uint idx_fx = 0; idx_fx < 2; ++idx_fx )
    {
        for ( uint idx_z = 0; idx_z < N_Z_; ++idx_z )
        {
            for ( uint idx_y = 0; idx_y < N_Y_; ++idx_y )
            {
                uint idx = 2*N_X_*N_Z_ + idx_fx*N_Y_*N_Z_ + idx_z*N_Y_ + idx_y;
                Vector3<number> Vertex = Wireframe.col(idx);
                
                file_wireframe << Vertex.adjoint() << std::endl;
            }
            
            file_wireframe << std::endl;
        }
        
        file_wireframe << std::endl;
    }
    
    // Save (XY) faces
    for ( uint idx_fz = 0; idx_fz < 2; ++idx_fz )
    {
        for ( uint idx_x = 0; idx_x < N_X_; ++idx_x )
        {
            for ( uint idx_y = 0; idx_y < N_Y_; ++idx_y )
            {
                uint idx = 2*N_X_*N_Z_ + 2*N_Y_*N_Z_+ idx_fz*N_X_*N_Y_ + idx_x*N_Y_ + idx_y;
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
void TwistedCuboid<number>::SaveMesh(const Matrix3X<number>& Wireframe, uint num_tri)
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
    
    // Build (XZ) vertex list
    for ( uint idx_fy = 0; idx_fy < 2; ++idx_fy )
    {
        for ( uint idx_z = 0; idx_z < N_Z_-1; ++idx_z )
        {
            for ( uint idx_x = 0; idx_x < N_X_-1; ++idx_x )
            {
                Vector3<number> Vtx1 = Wireframe.col(idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x);
                Vector3<number> Vtx2 = Wireframe.col(idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + 1);
                Vector3<number> Vtx3 = Wireframe.col(idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + N_X_);
                Vector3<number> Vtx4 = Wireframe.col(idx_fy * N_X_*N_Z_ + idx_z*N_X_ + idx_x + N_X_ + 1);

                file_mesh << Vtx1(0) << " " << Vtx1(1) << " " << Vtx1(2) << std::endl;
                file_mesh << Vtx2(0) << " " << Vtx2(1) << " " << Vtx2(2) << std::endl;
                file_mesh << Vtx3(0) << " " << Vtx3(1) << " " << Vtx3(2) << std::endl;
                file_mesh << Vtx4(0) << " " << Vtx4(1) << " " << Vtx4(2) << std::endl;
            }
        }
    }

    // Build (YZ) vertex list
    for ( uint idx_fx = 0; idx_fx < 2; ++idx_fx )
    {
        for ( uint idx_z = 0; idx_z < N_Z_-1; ++idx_z )
        {
            for ( uint idx_y = 0; idx_y < N_Y_-1; ++idx_y )
            {
                Vector3<number> Vtx1 = Wireframe.col(2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y);
                Vector3<number> Vtx2 = Wireframe.col(2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + 1);
                Vector3<number> Vtx3 = Wireframe.col(2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + N_Y_);
                Vector3<number> Vtx4 = Wireframe.col(2*N_X_*N_Z_ + idx_fx * N_Y_*N_Z_ + idx_z*N_Y_ + idx_y + N_Y_ + 1);
                
                file_mesh << Vtx1(0) << " " << Vtx1(1) << " " << Vtx1(2) << std::endl;
                file_mesh << Vtx2(0) << " " << Vtx2(1) << " " << Vtx2(2) << std::endl;
                file_mesh << Vtx3(0) << " " << Vtx3(1) << " " << Vtx3(2) << std::endl;
                file_mesh << Vtx4(0) << " " << Vtx4(1) << " " << Vtx4(2) << std::endl;
            }
        }
    }
    
    // Build (XY) vertex list
    for ( uint idx_fz = 0; idx_fz < 2; ++idx_fz )
    {
        for ( uint idx_x = 0; idx_x < N_X_-1; ++idx_x )
        {
            for ( uint idx_y = 0; idx_y < N_Y_-1; ++idx_y )
            {
                Vector3<number> Vtx1 = Wireframe.col(2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y);
                Vector3<number> Vtx2 = Wireframe.col(2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + 1);
                Vector3<number> Vtx3 = Wireframe.col(2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + N_Y_);
                Vector3<number> Vtx4 = Wireframe.col(2*N_X_*N_Z_ + 2*N_Y_*N_Z_ + idx_fz * N_X_*N_Y_ + idx_x*N_Y_ + idx_y + N_Y_ + 1);
                
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

template class TwistedCuboid<float>;
template class TwistedCuboid<double>;
