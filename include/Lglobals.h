/*

   Copyright (C) 2001,2002,2003,2004 Michael Rubinstein

   This file is part of the L-function package L.

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   Check the License for details. You should have received a copy of it, along
   with the package; see the file 'COPYING'. If not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

*/


#ifndef Lglobals_H
#define Lglobals_H

using namespace std;

#include <valarray>

//#ifdef INCLUDE_PARI
    //#include <pari/pari.h>
//#endif

//set Double according to the choice specified in the Makefile --------------
#ifdef USE_DOUBLE
    typedef double Double;
#endif

#ifdef USE_LONG_DOUBLE
    typedef long double Double;
#endif

#ifdef USE_BAILEY_DD
    #include <qd/dd_real.h>
    //#include "dd_real.h"
    typedef dd_real Double;
#endif

#ifdef USE_BAILEY_QD
    #include <qd/qd_real.h>
    //#include "qd_real.h"
    typedef qd_real Double;
#endif

#ifdef USE_MPFR
    #include "Lgmpfrxx.h"
    typedef mpfr_class Double;
#endif

#ifdef USE_MPFRCPP
    #include "mpreal.h"
    //#include "Lmpreal.h"
    using namespace mpfr;
    typedef mpreal Double;

    template<class T> inline void reset(T& t) {
        mpfr_clear(t.get_mpfr_t());
        mpfr_init(t.get_mpfr_t());
    }


/*
mpfr_srcptr get_mpfr_t() const { return mp; }
mpfr_ptr get_mpfr_t() { return mp; }
*/


#endif
//---------------------------------------------------------------------------


#include "Lcomplex.h"     //for complex numbers
//#include <complex>     //for complex numbers
typedef complex<Double> Complex;

//#ifdef USE_LONG_DOUBLE
//    #include "Lcommon_ld.h"
//#else
//    #include "Lcommon.h"
//#endif
#include "Lcommon.h"

#include<limits>
#include<iostream>

//#include "Lint_complex.h"
#include <math.h>

//--------Incomplete gamma function global variables--------

extern Complex last_z;         // the last z to be considered in inc_GAMMA
extern Complex last_w;         // the last w to be considered in inc_GAMMA
extern Complex last_comp_inc_GAMMA; // g(last_z,last_w)

extern Complex last_z_GAMMA;  //the last z to be considered in GAMMA
extern Complex last_log_G;    //the last log(GAMMA(z));

extern Double temme_a[1002],temme_g[501];

//used in Temme's asymptotic expansion of the
//incomplete gamma function
//XXXX might need more terms if I go to higher precision

//----------variables related to my cosine function --------------

extern Double *cos_taylor; //table of taylor coefficients for cosine
extern int cos_taylor_arraysize;
extern int number_cos_taylor_terms; //the number of taylor terms per series. Should be even.

extern Double one_over_cos_taylor_arraysize;
extern Double twoPi_over_cos_taylor_arraysize;




//-----Constants----------------------------------------------
extern Double Pi;
extern Double twoPi;
extern Double one_over_twoPi;
extern Double log_2Pi;
extern Complex I;

extern bool only_use_dirichlet_series; //whether to compute just using the Dirichlet series
extern int N_use_dirichlet_series; //if so, how many terms in the Dirichlet series to use.


//-----Global variables----------------------------------------
extern int my_verbose;       // verbosity level: 0 means no verbose

extern int DIGITS, DIGITS2; // precision and sacrifice
extern int DIGITS3; // how many digits to output
extern int DIGITS_xxx; // how many digits to output as determined by the explicit formula
extern Double xxx_max_DIFF; //maximum difference recorded comparing lhs to rhs of the explicit formula

extern Double tolerance;
extern Double tolerance_sqrd;
extern Double tolerance2;
extern Double tolerance3;

extern int global_derivative;  //used to specify which derivative to compute

extern int max_n; //the largest n used in a dirichlet series while computing a value

extern Double A; //controls the 'support' of g(w) in Riemann sum method
extern Double incr; //the increment in the Riemann sum method
extern Double tweak; //used in value_via_Riemann_sum to play with the angle

extern Double *LG;         // lookup table for log(n)
extern Double *two_inverse_SQUARE_ROOT;         // lookup table for sqrt(n)
extern int number_sqrts;    // how many sqrt(n)'s to store
extern int number_logs;    // how many log(n)'s to store

extern Double *bernoulli;  // lookup table for bernoulli numbers
extern Double rs_remainder[40][72]; //taylor coefficients for Riemann Siegel correction terms

extern Double hermite_norm[201]; // stores 1/sqrt(2^n*n!*sqrt(Pi)). Used in explicit formula, Hermite test function

//extern bool kronecker_table_available;
//extern int **kronecker_table; //lookup table for the kronecker function
//extern int kronecker_bound; //how large a table to make

extern bool print_warning;

extern Long my_LLONG_MAX;

extern int *prime_table;
extern int number_primes;

extern const Double sin_cof[];//={1.,-1./6.,1./120.,-1./5040.,1./362880.,-1./39916800.};
extern const Double sin_tol;

// Riemann Siegel band limited interpolation ----------------------------

extern bool try_use_blfi; //whether to compute dirichlet series using blf interpolation
extern Double blfi_interval_length;
extern Complex* block_value;

//================= older blfi variables to remove once new blfi is working
extern bool do_blfi;
extern const Double sinh_mult_fac;
extern const int sin_terms;

extern const Double blfi_block_growth; // how fast blfi blocks grow as we traverse the main sum, keep as is for now
extern const Double beta_fac_mult;  // controls the density of blfi sampling and the number of blfi terms needed
extern const Double blfi_fac;  // efficiency of the blfi interpolation sum relative to an RS sum of same length
extern const Double pts_array_fac;

extern const int rs_blfi_N;

extern Double *klog0; //log(k) at the beginning
extern Double *klog2; //log(k) at the end if needed
extern Double *ksqrt0; // 1/sqrt(k) at the beginning
extern Double *ksqrt2;// 1/sqrt(k) at the end if needed
extern int *num_blocks; // number of blocks
extern int *size_blocks;// size of blocks
extern Double *trig; // stores correction terms

extern Double **klog1; //log(k) in the middle if needed
extern Double **ksqrt1; // 1/sqrt(k) in the middle if needed
extern Double **klog_blfi; //initial term
extern Double **qlog_blfi; //band-width
extern Double **piv_org; //original pivot
extern Double **bbeta; //beta
extern Double **blambda; //lambda
extern Double **bepsilon; //epsilon
extern Double **arg_blfi; //arg_blfi
extern Double **inv_arg_blfi; //inv_arg_blfi

extern Double ***qlog_blfi_dense; // log(1+k/v) terms
extern Double ***qsqrt_blfi_dense; // 1/sqrt(1+k/v)
extern int ***blfi_done_left; //block done or not
extern int ***blfi_done_right; //block done or not
extern Double ***blfi_val_re_left; //real value of block
extern Double ***blfi_val_re_right; //real value of block
extern Double ***blfi_val_im_left; //imag value of block
extern Double ***blfi_val_im_right; //imag value of block

extern int length_org; // length of the main sum
extern int length_split; // length of the portion of the main sum to be evaluated directly
extern int lgdiv; // number of divisions of the main sum into intervals of the form [N,2N)
extern int max_pts; // max number of interpolation points allowed
extern int range; // number of blfi interpolation points needed
extern int blfi_block_size_org; // starting length of the blfi block
extern int total_blocks;

extern Double bc;
extern Double bc2;
extern Double kernel_fac;
extern Double ler;
extern Double mult_fac;
extern Double approx_blfi_mean_spacing;
extern Double interval_length;
extern Double error_tolerance;
extern Double input_mean_spacing;
extern Double input_mean_spacing_given;


//-----intializing and cleaning up routines----------------------

void initialize_globals(int n=200);
void delete_globals();

void extend_LG_table(int m);
void extend_sqrt_table(int m);

void extend_prime_table(int m);
int get_prime(int j);

void initialize_cos_array();

void initialize_rs_remainder1();
void initialize_rs_remainder2();
void initialize_rs_remainder3();
void initialize_rs_remainder4();
void initialize_rs_remainder5();
void initialize_rs_remainder6();
void initialize_rs_remainder7();
void initialize_rs_remainder8();

//----- used in one of the gamma routines. put it here since it is called
//----- during initialize_globals

Double dfac(int i);



inline Double my_norm(Complex z)
{
    return(real(z*conj(z)));
}

#endif
