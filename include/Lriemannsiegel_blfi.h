#ifndef Lriemannsiegel_blfi_H
#define Lriemannsiegel_blfi_H

#include "L.h"

using namespace std;

//inline Double sinc(Double u) {
//// Computes sinc(u)=sin(u)/u
//
//        // return 1.;
//        if(abs(u) < 1.e-5) {
//                return  1/120.0 * pow(u, 4) - 1.0/6.0 * pow(u, 2) + 1;
//        }
//        else {
//                return sin(u) / u;
//        }
//}


// Computes the kernel function in BLFI
// XXXXXX Needs to be made multiprecision
inline Double blfi_kernel(Double x) {

        // return 1.;

        Double one=Double(1);
        if( x < 0 ) return Double(0);
        else if(x < 1.e-5) {
                //return one/39916800*pow(x, 5) + one/362880*pow(x, 4) + one/5040*pow(x, 3) + one/120 *pow(x,2) + one/6 * x + 1;
                return ((((one/39916800*x + one/362880)*x + one/5040)*x + one/120) *x + one/6) * x + 1;
        }
        else {
                Double sqrt_x = sqrt(x);
                return sinh(sqrt_x)/sqrt_x;
        }
}

#endif
