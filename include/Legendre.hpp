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
 * Legendre.hpp: Version 1.0
 * Created 12/09/2015 by Maxime Tortora
 */
// ===================================================================


template<typename number>
struct Legendre
{
    // n = 0
    static inline number P0(number)   {return 1.;}
    
    // n = 1
    static inline number P1(number x) {return x;}
    
    // n = 2
    static inline number P2(number x) {return ((3. * x*x) - 1.) * 0.5;}
    
    // n
    static inline number Pn(unsigned int n, number x)
    {
        if      ( n == 0 )          return P0(x);
        else if ( n == 1 )          return P1(x);
        else if ( n == 2 )          return P2(x);
        
        if ( x ==  1. )             return 1.;
        if ( x == -1. )             return ((n % 2 == 0) ? 1. : -1.);
        if ( (x == 0.) && (n % 2) ) return 0.;
        
        number pnm1(P2(x));
        number pnm2(P1(x));
        number pn(pnm1);
        
        for ( unsigned int l = 3; l <= n; ++l )
        {
            pn = (((2. * (number)l) - 1.) * x * pnm1 - (((number)l - 1.) * pnm2)) / (number)l;
            pnm2 = pnm1;
            pnm1 = pn;
        }
        
        return pn;
    }
};

#endif
