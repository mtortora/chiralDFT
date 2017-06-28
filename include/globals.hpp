#ifndef GLOBALS_HPP_
#define GLOBALS_HPP_

#include <iostream>

#include "params.hpp"


// ============================
/* Path & string handling */
// ============================

#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)

#define __DATA_PATH TOSTRING(__DPATH__)
#define __MESOGEN TOSTRING(MESOGEN)


// ============================
/* Custom symbols */
// ============================

#define MPI_MASTER    0

#define MODE_PERT     0
#define MODE_FULL     1
#define MODE_EXC      2

#define ODF_FULL      0
#define ODF_LEGENDRE  1

#define DH_OXDNA      0
#define DH_FERRARINI  1
#define DH_WENSINK    2
#define DH_LAGERWALL  3

#define TREE_OB       0
#define TREE_SC       1

#define PI            3.1415926535897932

#define D_ALPHA       (2.*PI / N_ALPHA)
#define D_THETA       (1.*PI / N_THETA)
#define D_PHI         (2.*PI / N_PHI)

#define IS_BIAXIAL    ((N_ALPHA > 1) || (N_PHI > 1))
#define E_DIM         ((ODF_TYPE == ODF_FULL) ? N_ALPHA*N_THETA*N_PHI : N_L)

#define SIZE_L        ((N_L+1)*(2*N_L+1) + 2*N_L * (N_L+1)*(2*N_L+1)/3)


#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))


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

#define __EIGDENSE TOSTRING(__EIGDENSE__)
#define EIGEN_DENSEBASE_PLUGIN __EIGDENSE

#include <eigen3/Eigen/Dense>


typedef unsigned int           uint;
typedef unsigned long long int ullint;

template<typename T> using Matrix22 = Eigen::Matrix<T, 2, 2>;
template<typename T> using Matrix33 = Eigen::Matrix<T, 3, 3>;
template<typename T> using Matrix3X = Eigen::Matrix<T, 3, Eigen::Dynamic>;
template<typename T> using MatrixXX = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

template<typename T> using Vector2  = Eigen::Matrix<T, 2, 1>;
template<typename T> using Vector3  = Eigen::Matrix<T, 3, 1>;
template<typename T> using VectorX  = Eigen::Matrix<T, Eigen::Dynamic, 1>;

template<typename T> using ArrayX   = Eigen::Array<T, Eigen::Dynamic, 1>;
template<typename T> using ArrayXX  = Eigen::Array<T, Eigen::Dynamic, Eigen::Dynamic>;

#endif
