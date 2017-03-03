#ifndef LEGENDRE_HPP_
#define LEGENDRE_HPP_

// ===================================================================
/**
 * Legendre polynomials through a fast recursive iteration algorithm
 * Adapted from James A. Chappell:
 * https://github.com/jachappell/Legendre-polynomials/blob/master/legendre.h
 */
// ===================================================================
/*
 * Legendre.hpp:  Version 1.0
 * Created 12/09/2015 by Maxime Tortora
 */
// ===================================================================


namespace Legendre
{
    // n = 0
    static inline double P0(double)   {return 1.;}
    
    // n = 1
    static inline double P1(double x) {return x;}
    
    // n = 2
    static inline double P2(double x) {return ((3. * x*x) - 1.) * 0.5;}
    
    // n
    static inline double Pn(unsigned int n, double x)
    {
        if      ( n == 0 )          return P0(x);
        else if ( n == 1 )          return P1(x);
        else if ( n == 2 )          return P2(x);
        
        if ( x ==  1. )             return 1.;
        if ( x == -1. )             return ((n % 2 == 0) ? 1. : -1.);
        if ( (x == 0.) && (n % 2) ) return 0.;
        
        double pnm1(P2(x));
        double pnm2(P1(x));
        double pn(pnm1);
        
        for ( unsigned int l = 3; l <= n; ++l )
        {
            pn = (((2. * (double)l) - 1.) * x * pnm1 - (((double)l - 1.) * pnm2)) / (double)l;
            pnm2 = pnm1;
            pnm1 = pn;
        }
        
        return pn;
    }
}

#endif
