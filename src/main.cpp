// ===================================================================
/**
 * Monte-Carlo based DFT simulations of liquid-crystalline chiral
 * nematic phases with object-oriented MPI implementation.
 * -------
 * Requires working Eigen/OpenMPI installs (see ReadMe.md)
 */
// ===================================================================
/*
 * main.cpp: Version 2.5
 * Created 16/09/2015 by Maxime Tortora
 */
// ===================================================================

#include <mpi.h>

#include "SimManager.hpp"


int main(int argc, char* argv[])
{
    int mpi_rank;
    int mpi_size;
    
    // Setup MPI
    MPI_Init(&argc, &argv);
    
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    
    try
    {
        // Simulation constructor
        SimManager MCsim(mpi_rank, mpi_size);
        
        // Seed MPI MC integrator and load particle templates
        MCsim.MPIInit();
        
        // Initial perturbative run
        MCsim.InitRun(MC_TYPE);
        
        if ( MODE == MODE_FULL )
        {
            // Main chiral sweeping run
            MCsim.LandscapeRun();
        
            // Minimise chiral free energy surface
            MCsim.MinSurf();
        }
        
        if ( MODE != MODE_EXC )
        {
            // Broadcast results to master thread
            MCsim.Gather();
        
            // Save aggregated data
            MCsim.Save();
        }
    }
    
    // Treat exceptions as fatal
    catch ( std::exception& e )
    {
        LogErr("%s on thread %u - aborting", e.what(), mpi_rank);
        
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        exit     (EXIT_FAILURE);
    }
    
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}
