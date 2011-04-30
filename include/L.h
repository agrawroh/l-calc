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


#ifndef L_H
#define L_H

#ifdef _OPENMP
#include <omp.h>
#endif

#include <iomanip>          //for manipulating output such as setprecision
#include <fstream>          //for file input output
#include <string.h>         //for string functions such as strcpy
#include <string>
#include <iostream>         //for ostrstream
//#include <strstream>      //depreceated

#include <math.h>
#include<vector>

#include "Lglobals.h"            //for global variables
#include "Lmisc.h"               //things like sn or LOG
#include "Lgamma.h"              //incomplete gamma function code
//#include "Lprecomp.h"
#include "Lriemannsiegel.h"      //Riemann Siegel formula
#include "Lriemannsiegel_blfi.h" //Hiary's Riemann Siegel formula using band limited interpolation
#include "Lnumbertheory.h"       //basic number theory routines
#include "Lelliptic.h"           //for initializing an elliptic curve L-function, if we have pari



//-----THE L Function Class-----------------------------------
// For the sake of uniformity, I assume that 
// L(s) = \sum_1^\infty b(n)/n^s satisfies a functional equation
// of the form:
// Let Lambda(s) = Q^s product_{j=1}^a Gamma(gamma_j s + lambda_j) L(s)
// functional equation:
// Lambda(s) = OMEGA \conjugate{\Lambda(1-\conjugate{s})}
// We also allow for simple poles of Lambda(s) to take into account
// the Riemann zeta function.


template <class ttype>
class L_function
{
public:


    char *name;                             // the name of the L_function

    int what_type_L;                        // -1 for zeta
                                            // 0 for unknown
                                            // 1 for periodic, i.e an L(s,chi), Dirichlet L-function
                                            // -2 for L(s,chi) but where the
                                            // number of coeffs computed is < period.
                                            // 2 for cusp form (in S_K(Gamma_0(N))
                                            // 3 for Maass form for SL_2(Z)
                                            // 4 for dedekind zeta function
                                            // 5 etc for future types

    int number_of_dirichlet_coefficients;   // the number of dirichlet coefficients

    ttype *dirichlet_coefficient;           // the dirichlet coefficients, i.e. a(n) above.

    long long period;                       //stores the period.
                                            //0 if not periodic.

    Double Q;                               // from the Q^s factor in the functional equation
    Complex OMEGA;                          // the root number 
                                            // we assume |OMEGA| = 1 for the following:
                                            // rotating L(s) so as to be real on the
                                            // critical line, and in check_funct_equation.

    int a;                                  // quasi degree of the L-function
    Double *gamma;                          // for the GAMMA factors. Either 1/2 or 1
    Complex *lambda;                        // for the GAMMA factors

                                            // for notational convenience we label 
                                            // them 1..a

    int number_of_poles;                    // the number of poles
    Complex *pole;                          // poles of the L-function
    Complex *residue;                       // residues at the poles

    //Complex S_0;                            // taylor series about S_0
    //Complex *taylor_series;                 // taylor coefficients
    //precomputation2 *local_series_a;        // precomputed taylor expansions for f1 integrand
    //precomputation2 *local_series_b;        // precomputed taylor expansions for f2 integrand

    //-----Constructor: default initilization is to the Riemann zeta function-------

    L_function ()
    {
        if(my_verbose>1)
            cout << "zeta constructor called\n";
        name = new char[5];
        strcpy(name,"zeta");
        what_type_L=-1;  // this is how I know it is zeta
        number_of_dirichlet_coefficients=0;
        dirichlet_coefficient = new ttype[1];
        period=0;
        Q=1/sqrt(Pi);
        OMEGA=1.;
        a=1;
        gamma = new Double[2];
        lambda = new Complex[2];
        gamma[1]=.5;
        lambda[1]=0;
        number_of_poles=2;
        pole = new Complex[3];
        residue = new Complex[3];
        pole[1]=1;
        residue[1]=1;
        pole[2]=0;
        residue[2]=-1;


        //S_0=.5;
        //taylor_series= new Complex[number_taylor_coeffs+1];

        //local_series_a=new precomputation2[number_local_series];
        //local_series_b=new precomputation2[number_local_series];
    }


    //-----Constructor: initialize the L-function from given data------------------
    L_function (const char *NAME, int what_type, int N, ttype *coeff, long long Period,
    Double q,  Complex w, int A, Double *g, Complex *l,
    int n_poles, Complex *p, Complex *r)
    {
        if(my_verbose>1)
            cout << "constructor called\n";
        int k;
        bool use_legendre_duplication = false; //at present there's no need to call legendre.

        name = new char[strlen(NAME)+1];
        strcpy(name,NAME);
        what_type_L=what_type;
        number_of_dirichlet_coefficients=N;
        dirichlet_coefficient = new ttype[N+1];
        for (k=1;k<=N;k++)
        {
            dirichlet_coefficient[k]= coeff[k];
            if(my_verbose>1&&k<=10)
            cout << "setting dirichlet coefficient" << k << " "
            << coeff[k]<< " "
            << dirichlet_coefficient[k]<< endl;
        }
        period=Period;
        Q=q;
        OMEGA=w;
        a=A;
        if(use_legendre_duplication){
            for (k=1;k<=A;k++) if (1.1-g[k]<.2&&A>1) a++; //i.e. if g[k]=1 as opposed to 1/2
        }
        gamma = new Double[a+1];
        lambda = new Complex[a+1];
        int j=A+1;
        for (k=1;k<=A;k++)
        {
            if (use_legendre_duplication&&1.1-g[k]<.2&&A>1)
            {
                gamma[k]=g[k]*.5;
                gamma[j]=g[k]*.5;
                lambda[k]=l[k]*.5;
                lambda[j]=l[k]*.5+.5;
                Q=2*Q;
                j++;
            }
            else
            {
                gamma[k]=g[k];
                lambda[k]=l[k];
            }
        }
        number_of_poles=n_poles;
        pole = new Complex[n_poles+1];
        residue = new Complex[n_poles+1];
        for (k=1;k<=n_poles;k++)
        {
            pole[k]=p[k];
            residue[k]=r[k];
        }

        if(my_verbose>2){
            cout << "    what_type_L: " << what_type_L << endl;
            cout << "    number_of_dirichlet_coefficients: " << number_of_dirichlet_coefficients << endl;
            cout << "    Period: " << period << endl;
        }


    }


    //-----Constructor: initialize the L-function from given data------------------
    L_function (const char *NAME, int what_type, int N, ttype *coeff, long long Period,
    Double q,  Complex w, int A, Double *g, Complex *l) //this one assumes no poles
    {
        if(my_verbose>1)
            cout << "constructor called\n";
        int k;
        bool use_legendre_duplication = false;

        name = new char[strlen(NAME)+1];
        strcpy(name,NAME);
        what_type_L=what_type;
        number_of_dirichlet_coefficients=N;
        dirichlet_coefficient = new ttype[N+1];
        for (k=1;k<=N;k++)
        {
            dirichlet_coefficient[k]= coeff[k];
            if(my_verbose>1&&k<=10)
            cout << "setting dirichlet coefficient" << k << " "
            << coeff[k]<< " "
            << dirichlet_coefficient[k]<< endl;
        }
        period=Period;
        Q=q;
        OMEGA=w;
        a=A;
        if(use_legendre_duplication){
            for (k=1;k<=A;k++) if (1.1-g[k]<.2&&A>1) a++; //i.e. if g[k]=1 as opposed to 1/2
        }
        gamma = new Double[a+1];
        lambda = new Complex[a+1];
        int j=A+1;
        for (k=1;k<=A;k++)
        {
            if (use_legendre_duplication&&1.1-g[k]<.2&&A>1)
            {
                gamma[k]=g[k]*.5;
                gamma[j]=g[k]*.5;
                lambda[k]=l[k]*.5;
                lambda[j]=l[k]*.5+.5;
                Q=2*Q;
                j++;
            }
            else
            {
                gamma[k]=g[k];
                lambda[k]=l[k];
            }
        }
        number_of_poles=0;
        pole = new Complex[1];
        residue = new Complex[1];

    }

    //-----Constructor: initialize the L-function for an elliptic curve -----------------
#ifdef INCLUDE_PARI
    L_function (char *a1, char *a2, char *a3, char *a4, char *a6, int N_terms)
    {

        if(my_verbose>1)
            cout << "L elliptic constructor called\n";
        name = new char[15];
        strcpy(name,"elliptic curve");
        what_type_L=2;
        number_of_dirichlet_coefficients=N_terms;
        dirichlet_coefficient = new ttype[N_terms+1];

        period=0;

        a=1;

        //GAMMA factor is GAMMA(s+1/2)
        gamma=new Double[2];
        lambda=new Complex[2];
        gamma[1]=1.;
        lambda[1]=.5;

        data_E(a1,a2,a3,a4,a6,N_terms,dirichlet_coefficient);
        //coeff[n], if n > 1, is the nth dirichlet coefficient normalized by sqrt(n).
        //coeff[0] comes back with the sign of the functional equation which we pass to w.
        //coeff[1] comes back with the conductor of E. We then set it to 1.


        Q=sqrt(dirichlet_coefficient[1])/(2*Pi);
        dirichlet_coefficient[1]=1.;

        OMEGA=dirichlet_coefficient[0];


        number_of_poles=0;
        pole = new Complex[1];
        residue = new Complex[1];



    }

#endif //ifdef INCLUDE_PARI

    //-----Copy constructor-------------------------
    L_function (const L_function &L)
    {

        if(my_verbose>1)
            cout << "copy called\n";

        int k;
        name = new char[strlen(L.name)+1];
        strcpy(name,L.name);
        what_type_L=L.what_type_L;
        number_of_dirichlet_coefficients=L.number_of_dirichlet_coefficients;;
        dirichlet_coefficient = new ttype[number_of_dirichlet_coefficients+1];
        for (k=1;k<=number_of_dirichlet_coefficients;k++)
        {
            dirichlet_coefficient[k]= L.dirichlet_coefficient[k];
            if(my_verbose>1&&k<=10)
            cout << "setting dirichlet coefficient" << k << " "
            << L.dirichlet_coefficient[k] << " "
            << dirichlet_coefficient[k]<< endl;
        }
        period=L.period;
        #if defined(USE_MPFR) || defined(USE_MPFRCPP)
            reset(Q);
            reset(OMEGA);
        #endif
        Q=L.Q;
        OMEGA=L.OMEGA;
        a=L.a;
        gamma = new Double[a+1];
        lambda = new Complex[a+1];
        for (k=1;k<=a;k++)
        {
            gamma[k]=L.gamma[k];
            lambda[k]=L.lambda[k];
        }
        number_of_poles=L.number_of_poles;
        pole = new Complex[number_of_poles+1];
        residue = new Complex[number_of_poles+1];
        for (k=1;k<=number_of_poles;k++)
        {
            pole[k]=L.pole[k];
            residue[k]=L.residue[k];
        }

        //S_0=L.S_0;
        //taylor_series= new Complex[number_taylor_coeffs+1];
        //for(k=0;k<=number_taylor_coeffs;k++)taylor_series[k]=L.taylor_series[k];

        //local_series_a=new precomputation2[number_local_series];
        //local_series_b=new precomputation2[number_local_series];
        //for(k=0;k<number_local_series;k++){
            //local_series_a[k]=L.local_series_a[k];
            //local_series_b[k]=L.local_series_b[k];
        //}

    }

    //-----Assignment operator--------------------
    L_function & operator = (const L_function &L)
    {
        int k;
        if(my_verbose>1)
            cout << "assignment called\n";
        if (this != &L)
        {
            delete [] name;
            name=new char[strlen(L.name)+1];
            strcpy(name,L.name);
            what_type_L=L.what_type_L;
            number_of_dirichlet_coefficients=L.number_of_dirichlet_coefficients;
            delete [] dirichlet_coefficient;
            dirichlet_coefficient = new ttype[number_of_dirichlet_coefficients+1];
            for (k=1;k<=number_of_dirichlet_coefficients;k++)
                dirichlet_coefficient[k]= L.dirichlet_coefficient[k];
            period=L.period;
            #ifdef USE_MPFR
                reset(Q);
                reset(OMEGA);
            #endif
            Q=L.Q;
            OMEGA=L.OMEGA;
            a=L.a;
            delete [] gamma;
            gamma = new Double[a+1];
            delete [] lambda;
            lambda = new Complex[a+1];
            for (k=1;k<=a;k++)
            {
                gamma[k]=L.gamma[k];
                lambda[k]=L.lambda[k];
            }
            number_of_poles=L.number_of_poles;
            delete [] pole;
            pole = new Complex[number_of_poles+1];
            delete [] residue;
            residue = new Complex[number_of_poles+1];
            for (k=1;k<=number_of_poles;k++)
            {
                pole[k]=L.pole[k];
                residue[k]=L.residue[k];
            }


            //S_0=L.S_0;
            //delete [] taylor_series;

            //taylor_series= new Complex[number_taylor_coeffs+1];
            //for(k=0;k<=number_taylor_coeffs;k++)taylor_series[k]=L.taylor_series[k];

            //local_series_a=new precomputation2[number_local_series];
            //local_series_b=new precomputation2[number_local_series];
            //for(k=0;k<number_local_series;k++){
                //local_series_a[k]=L.local_series_a[k];
                //local_series_b[k]=L.local_series_b[k];
            //}


        }
        return *this;
    }

    //-----Destructor: free allocated memory------------------
    ~L_function ()
    {
        if(my_verbose>1)
            cout << "destructor called\n";
        delete [] name;
        delete [] dirichlet_coefficient;
        delete [] gamma;
        delete [] lambda;
        delete [] pole;
        delete [] residue;
        //delete [] taylor_series;
        //delete [] local_series_a;
        //delete [] local_series_b;
    }


    //-----addition operator--------------------
    //returns the 'L-function' whose basic data (funct eqn, dirichlet coefficients) 
    //is that of this added to that of L
    //not particularly useful- was used for a crazy experiment to deform L-functions
    L_function operator + (const L_function &L)
    {
        L_function L2;
        int k;
        if(my_verbose>1)
            cout << "addition called\n";
        L2.name=new char[1];
        strcpy(L2.name,"");
        L2.what_type_L=L.what_type_L;
        L2.number_of_dirichlet_coefficients=L.number_of_dirichlet_coefficients;
        L2.dirichlet_coefficient = new ttype[number_of_dirichlet_coefficients+1];
        for (k=1;k<=number_of_dirichlet_coefficients;k++)
            L2.dirichlet_coefficient[k]= dirichlet_coefficient[k] +L.dirichlet_coefficient[k];
        L2.period=L.period;
        #ifdef USE_MPFR //XXXXXXXXXXXX don't think this is needed since the assignment call does same
            reset(L2.Q);
            reset(L2.OMEGA);
        #endif
        L2.Q=Q+L.Q;
        L2.OMEGA=OMEGA+L.OMEGA;
        L2.a=L.a;
        L2.gamma = new Double[a+1];
        L2.lambda = new Complex[a+1];
        for (k=1;k<=a;k++)
        {
            L2.gamma[k]=gamma[k]+L.gamma[k];
            L2.lambda[k]=lambda[k]+L.lambda[k];
        }
        L2.number_of_poles=number_of_poles;
        L2.pole = new Complex[number_of_poles+1];
        L2.residue = new Complex[number_of_poles+1];
        for (k=1;k<=number_of_poles;k++)
        {
            L2.pole[k]=pole[k]+L.pole[k];
            L2.residue[k]=residue[k]+L.residue[k];
        }

        return L2;

    }
    //-----multiplication operator--------------------
    //returns the L-function whose basic data (funct eqn, dirichlet coefficients)
    //is that of this times the scalar t
    //not particularly useful- was used for a crazy experiment
    L_function operator * (Double t)
    {
        L_function L2;
        int k;
        if(my_verbose>1)
            cout << "addition called\n";
        L2.name=new char[1];
        strcpy(L2.name,"");
        L2.what_type_L=what_type_L;
        L2.number_of_dirichlet_coefficients=number_of_dirichlet_coefficients;
        L2.dirichlet_coefficient = new ttype[number_of_dirichlet_coefficients+1];
        for (k=1;k<=number_of_dirichlet_coefficients;k++)
            L2.dirichlet_coefficient[k]= dirichlet_coefficient[k]*t;
        L2.period=period;
        #ifdef USE_MPFR //XXXXXXXXXXXX don't think this is needed since the assignment call does same
            reset(L2.Q);
            reset(L2.OMEGA);
        #endif
        L2.Q=Q*t;
        L2.OMEGA=OMEGA*t;
        L2.a=a;
        L2.gamma = new Double[a+1];
        L2.lambda = new Complex[a+1];
        for (k=1;k<=a;k++)
        {
            L2.gamma[k]=gamma[k]*t;
            L2.lambda[k]=lambda[k]*t;
        }
        L2.number_of_poles=number_of_poles;
        L2.pole = new Complex[number_of_poles+1];
        L2.residue = new Complex[number_of_poles+1];
        for (k=1;k<=number_of_poles;k++)
        {
            L2.pole[k]=pole[k]*t;
            L2.residue[k]=residue[k]*t;
        }

        return L2;
    }


    //#include "Lprint.h"         //printing routine
    void print_data_L(int N=10); //prints basic data for an L-function

    //#include "Lnumberzeros.h"   //for N(T) and S(T)
    Double Nmain(Double T); //main term in N(T), i.e. N(T)-S(T)
    //Computes S(T) and N(T). Not rigorous, but practical.
    //assumes that T is not a zero of L(1/2+iT) or L(1/2-iT).
    Double S(Double T);
    Double N(Double T);
    Double Nmain_inverse(Long n); //finds T such that Nmain(T)=n
    Double density_zeros(Double T); //returns local density of zeros


    //#include "Lgram.h"          //for finding gram points
    Double initialize_gram(Double t);
    Double next_gram(Double t);
    Double nth_gram(Long n);


    //#include "Ldirichlet_series.h" //for computing Dirichlet series
    Complex partial_dirichlet_series(Complex s, long long N1, long long N2);
    Complex dirichlet_series(Complex s, long long N);

    Complex partial_dirichlet_series_via_blfi(Complex s, long long N1, long long N, long long K, Double epsilon);
    Complex dirichlet_series_via_blfi(Complex s, long long N, Double blfi_interval_length, Double epsilon);
    Complex dirichlet_series_block_blfi(Complex s, long long v, long long K, int c, int c0,
            long long center_location, long long index0, Double tau, Double beta, Double lambda, Double eps);
    Complex get_block_value_directly(Complex s, long long v, long long K);

    //#include "Ltaylor_series.h" //for computing taylor series for Dirichlet series
    //void compute_taylor_series(int N, int K, Complex s_0, Complex *series);
    //void compute_local_contribution(Complex *local_series,Complex z1,Complex z2);
    //void compute_taylor();
    //Complex value_via_taylor_series(Complex s,int n=0,char *return_type="pure");

    //#include "Lvalue.h"         //value via Riemann sum, via gamma sum, various options for value
    Complex find_delta(Complex s,Double g);
    Complex value_via_Riemann_sum(Complex s, const char *return_type="pure");
    Complex value_via_gamma_sum(Complex s, const char *return_type="pure");
    Complex value(Complex s, int derivative = 0, const char *return_type="pure", const char*method="default");

    //#include "Lfind_zeros.h"    //finding zeros routine
    Double zeros_zoom_brent(Double L1, Double L2, Double u, Double v);
    void find_zeros(Double t1, Double t2, Double step_size, const char* filename="cout", const char* message_stamp="");
    void  find_zeros_v(Double t1, Double t2, Double step_size,  vector<Double> &result);//This is the same as above function
    //void find_zeros_via_gram(Double t1, Long count=0,Double max_refine=1025, const char* filename="cout", const char* message_stamp="");
    int compute_rank(bool print_rank=false);
    void verify_rank(int rank);
    //void find_zeros_via_N(Long count=0,bool do_negative=true,Double max_refine=1025, int rank=-1, bool test_explicit_formula=false, const char* filename="cout", const char* message_stamp="");
    //void find_zeros_via_N_v(Long count,bool do_negative,Double max_refine, int rank, bool do_test_explicit_formula, vector<Double>  &result);
    int find_zeros(Long count, Long start_N=0,Double max_refine=1025, int rank=-1, const char* message_stamp="");
    //int find_zeros(Long count, Long start_N=0,Double max_refine=1025, int rank=-1, bool do_test_explicit_formula=false, const char* message_stamp="");
    bool is_complex();

    //#include "Ldokchitser.h"    //dokchitser algorithm for what he calls phi(t), i.e. inverse
    void phi_series(int precision);

    //#include "Lexplicit_formula.h"
    int dirichlet_coeffs_log_diff(int num_coeff, Complex *c);
    int test_explicit_formula(Double x_0, Double alpha, Double *zero_table, int number_zeros, Complex *c, int num_coeffs, const char *method, int n, vector<Double> polynomial);
    Double rhs_explicit_formula(Double x_0, Double alpha, Complex *c, int num_coeffs, const char *method, int n, vector<Double> polynomial);
    int plot_explicit_formula(Double alpha, Double x, Double x2, Double step_size, const char *method, int num_coeffs, Double *rhs_array);


};



//templated class code should be kept in .h files
#include "Ldirichlet_series.h" //for computing Dirichlet series
#include "Lprint.h"         //printing routine
#include "Lnumberzeros.h"   //for N(T) and S(T)
#include "Lgram.h"          //for finding gram points
//#include "Ltaylor_series.h" //for computing taylor series for Dirichlet series
#include "Lvalue.h"         //value via Riemann sum, via gamma sum, various options for value
#include "Lfind_zeros.h"    //finding zeros routine
#include "Ldokchitser.h"    //dokchitser algorithm for what he calls phi(t), i.e. inverse
                            //mellin transform

#include "Lexplicit_formula.h" //for testing zeros with the explicit formula

#endif
