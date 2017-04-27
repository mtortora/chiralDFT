#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_

#include <eigen3/Eigen/Dense>


// ============================
/* Simulation options */
// ============================

// Particle type to be used - implemented: BentCore, DNADuplex, Helix, PatchyRod, ThreadedRod, TwistedCuboid, TwistedPentagon
#define MESOGEN       DNADuplex

// Full run switch - set to 0 for full run, 1 for perturbative run or 2 for excluded volume only (BentCore, Helix)
#define MODE          1

// Set to 0 for full functional minimisation or 1 for Legendre-projected run (only for preliminary run)
#define MC_TYPE       0

// Use RAPID library for collision queries - only implemented for the TwistedCuboid & TwistedPentagon particle templates
#define USE_RAPID     0

// Set to 1 to enable Debye-Huckel interactions
#define USE_DH        1

// Type of soft interaction - 0: oxDNA-parametrised electrostatics, 1: Ferrarini, 2: Wensink, 3: Lagerwall
#define MODE_DH       0

// Bounding tree hierarchy mode - set to 0 for Oriented Bounding Box or 1 for SpheroCylinder
#define MODE_TREE     0


// ============================
/* Model parameters */
// ============================

// Maximum number of Monte-Carlo steps for excluded volume integration
#define N_MC          1E14

// Soft interaction parameters - only relevant for the DNADuplex template
#define T_ABS         293.16
#define C_SALT        0.26

// Set maximum range of concentrations to be simulated
#define ETA_MIN       0.01
#define ETA_MAX       0.27

// Set inverse pitch grid range
#define Q_MIN        -0.0002
#define Q_MAX         0.0005

// Gamma coefficient for the relaxed Newton-Raphson method
#define GAMMA_NR      0.5

// Parameter grid sizes
#define N_L           20
#define N_SYNC        50
#define N_STEPS_Q     15
#define N_STEPS_ETA   100
#define N_STEPS_THETA 250
#define N_STEPS_ODF   2500
#define N_STEPS_NR    5000

// Numerical tolerance coefficients
#define TOL_NEM       1e-2
#define TOL_BIN       1e-4
#define TOL_ODF       1e-6
#define TOL_SC        1e-6
#define TOL_OB        1e-6

// Excluded volume grid dimensions
#define NX            150
#define NY            150
#define NZ            1500


// ============================
/* Custom symbols */
// ============================

#define PI            3.1415926535897932

#define D_THETA       (PI / (N_STEPS_THETA-1.))

#define MPI_MASTER    0

// Model options
#define MODE_FULL     0
#define MODE_PERT     1
#define MODE_EXC      2

#define ODF_FULL      0
#define ODF_LEGENDRE  1

#define DH_OXDNA      0
#define DH_FERRARINI  1
#define DH_WENSINK    2
#define DH_LAGERWALL  3

#define TREE_OB       0
#define TREE_SC       1

#define SQR(x)        ((x)*(x))
#define CUB(x)        ((x)*(x)*(x))


// ============================
/* Colored logs */
// ============================

#define COLORS_ESCAPE "\033["

// ANSI escape sequences
#ifdef __XCODE__

#define COLORS_RESET COLORS_ESCAPE ";"
#define COLORS_RED "fg255,0,0;"
#define COLORS_GRE "fg0,255,0;"
#define COLORS_BLU "fg0,0,255;"
#define COLORS_PUR "fg255,0,255;"
#define COLORS_CYA "fg20,150,255;"
#define COLORS_TXT "fg0,0,0;"

#define COLORS_INF "fg0,0,0;"
#define COLORS_ERR "fg255,0,0;"

#else

#define COLORS_RESET COLORS_ESCAPE "0m"
#define COLORS_RED "1;31m"
#define COLORS_GRE "1;32m"
#define COLORS_BLU "1;34m"
#define COLORS_PUR "1;35m"
#define COLORS_CYA "1;36m"
#define COLORS_TXT "1;37m"

#define COLORS_INF "0;37m"
#define COLORS_ERR "0;31m"

#endif

// Loggers through variadic macros
#define LogRed(frmt, ...) fprintf(stdout, (COLORS_ESCAPE COLORS_RED frmt COLORS_RESET "\n"), ##__VA_ARGS__)
#define LogGre(frmt, ...) fprintf(stdout, (COLORS_ESCAPE COLORS_GRE frmt COLORS_RESET "\n"), ##__VA_ARGS__)
#define LogBlu(frmt, ...) fprintf(stdout, (COLORS_ESCAPE COLORS_BLU frmt COLORS_RESET "\n"), ##__VA_ARGS__)
#define LogPur(frmt, ...) fprintf(stdout, (COLORS_ESCAPE COLORS_PUR frmt COLORS_RESET "\n"), ##__VA_ARGS__)
#define LogCya(frmt, ...) fprintf(stdout, (COLORS_ESCAPE COLORS_CYA frmt COLORS_RESET "\n"), ##__VA_ARGS__)
#define LogTxt(frmt, ...) fprintf(stdout, (COLORS_ESCAPE COLORS_TXT frmt COLORS_RESET "\n"), ##__VA_ARGS__)

#define LogInf(frmt, ...) fprintf(stdout, (COLORS_ESCAPE COLORS_INF frmt COLORS_RESET "\n"), ##__VA_ARGS__)
#define LogErr(frmt, ...) fprintf(stderr, (COLORS_ESCAPE COLORS_ERR frmt COLORS_RESET "\n"), ##__VA_ARGS__)


// ============================
/* Custom typedefs */
// ============================

typedef unsigned int           uint;
typedef unsigned long long int ullint;

typedef Eigen::MatrixXd::Index                              IndexXd;
typedef Eigen::Array <bool, Eigen::Dynamic, 1>              ArrayXb;
typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXb;


// ============================
/* Path & string handling */
// ============================

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define __DATA_PATH TOSTRING(__DPATH__)
#define __MESOGEN TOSTRING(MESOGEN)

#endif
