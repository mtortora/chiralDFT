#ifndef PARAMS_HPP_
#define PARAMS_HPP_


// ============================
/* Simulation options */
// ============================

// Particle type to be used - implemented: BentCore, DNADuplex, Helix, PatchyRod, ThreadedRod, TwistedCuboid, TwistedPentagon
#define MESOGEN     DNADuplex

// Full run switch - MODE_PERT for perturbative run, MODE_FULL for full run or MODE_EXC for excluded volume (BentCore, Helix)
#define MODE_SIM    MODE_PERT

// Set to ODF_FULL for full functional minimisation or ODF_LEGENDRE for Legendre-projected run (only for preliminary run)
#define ODF_TYPE    ODF_FULL

// Type of soft interaction if relevant - 0: oxDNA-parametrised electrostatics, 1: Ferrarini, 2: Wensink, 3: Lagerwall
#define MODE_DH     DH_OXDNA

// Bounding tree hierarchy mode - set to TREE_OB for Oriented Bounding Box or TREE_SC for SpheroCylinder
#define MODE_TREE   TREE_OB

// Use RAPID library for collision queries - only implemented for the TwistedCuboid & TwistedPentagon particle templates
#define USE_RAPID   0

// Set to 1 to enable Debye-Huckel interactions
#define USE_DH      1


// ============================
/* Model parameters */
// ============================

// Maximum number of Monte-Carlo steps for excluded volume integration
#define N_MC        1E14

// Soft interaction parameters - only relevant for the DNADuplex template
#define T_ABS       293.16
#define C_SALT      0.26

// Set maximum range of concentrations to be simulated
#define ETA_MIN     0.01
#define ETA_MAX     0.27

// Set inverse pitch grid range
#define Q_MIN      -0.0002
#define Q_MAX       0.0005

// Gamma coefficient for the relaxed Newton-Raphson method
#define GAMMA_NR    0.5

// Parameter grid sizes
#define N_L         20
#define N_SYNC      50
#define N_STEPS_Q   15
#define N_STEPS_ETA 100
#define N_STEPS_ODF 10000
#define N_STEPS_NR  100

#define N_ALPHA     1
#define N_THETA     250
#define N_PHI       1

// Numerical tolerance coefficients
#define TOL_NEM     1e-2
#define TOL_BIN     1e-4
#define TOL_ODF     1e-6
#define TOL_SC      1e-6
#define TOL_OB      1e-6

// Excluded volume grid dimensions
#define NX          150
#define NY          150
#define NZ          1500

#endif
