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

#include "Lmisc.h"

#ifndef Ldirichlet_series_H
#define Ldirichlet_series_H


template <class ttype>
Complex L_function <ttype>::
partial_dirichlet_series(Complex s, long long N1, long long N2)
{
    Complex z=0.;
    long long m,n;

    if(what_type_L==-1)   //i.e. if the Riemann zeta function
        for(n=N1;n<=N2;n++) z+=lcalc_exp(-s*LOG(n));
        //for(n=N1;n<=N2;n++) z+=exp(-s*LOG(n));
    else if(what_type_L!=1) //if not periodic
        for(n=N1;n<=N2;n++) z+=dirichlet_coefficient[n]*lcalc_exp(-s*LOG(n));
        //for(n=N1;n<=N2;n++) z+=dirichlet_coefficient[n]*exp(-s*LOG(n));
    else //if periodic
        for(n=N1;n<=N2;n++)
        {
            m=n%period; if(m==0)m=period;
            z+=dirichlet_coefficient[m]*lcalc_exp(-s*LOG(n));
        }
    return z;
}

template <class ttype>
Complex L_function <ttype>::
dirichlet_series(Complex s, long long N=-1)
{
    Complex z=0.;
    long long m,n;
     if(N==-1) N=number_of_dirichlet_coefficients;
    if(N>number_of_dirichlet_coefficients&&what_type_L!=-1&&what_type_L!=1)
    {

        if(print_warning){
            print_warning=false;
            cerr << "# WARNING from dirichlet series- we don't have enough Dirichlet coefficients." << endl;
            cerr << "# Will use the maximum possible, though the output ";
            cerr << "# will not necessarily be accurate." << endl;
        }
        N=number_of_dirichlet_coefficients;
    }
    if(what_type_L==-1)   //i.e. if the Riemann zeta function
        for(n=1;n<=N;n++) z+=lcalc_exp(-s*LOG(n));
    else if(what_type_L!=1) //if not periodic
        for(n=1;n<=N;n++) z+=dirichlet_coefficient[n]*lcalc_exp(-s*LOG(n));
    else //if periodic
        for(n=1;n<=N;n++)
        {
            m=n%period; if(m==0)m=period;
            z+=dirichlet_coefficient[m]*lcalc_exp(-s*LOG(n));
        }
    return z;

}



// ======================================================================================
// ======== Some functions useful for experiments and debugging the blfi routine ========
// ======================================================================================

template <class ttype>
Complex L_function <ttype>::
partial_dirichlet_series_via_blfi(Complex s, long long N1=1, long long N=-1, long long K = 1, Double epsilon = 1.e-14) {
//this function computes sum_{n=N1}^N a_n*n^(-s) in blocks of fixed length K,
//for testing purposes

    Complex z=0.;
    long long m,n;
    long long length_of_remainder_series = (N - N1 + 1) % K;
    long long number_of_blocks = (N - N1 - length_of_remainder_series + 1) / K;

    //compute the remainder series
    z = z + partial_dirichlet_series(s, N - length_of_remainder_series + 1, N);
    //cout << "z: " << z << endl;

    //compute the rest of the sum, which has length divisible by K
    for(n = 0; n < number_of_blocks; n++) {
            z = z + dirichlet_series_block_blfi(s, N1 + n*K, K, epsilon/number_of_blocks);
            //cout << n << " z: " << z << endl;
    }
    return z;
}

// ======================================================================================
// ======================== The blfi routine ============================================
// ======================================================================================

template <class ttype>
Complex L_function <ttype>::
get_block_value_directly(Complex s, long long v, long long K) {
// Function to compute sum_{k=0}^{K-1} a_{v+k} * exp(i*t*log(1 + k/v)) / (v+k)^sigma
// directly, where t=imag(s) and sigma=real(s)

    long long m,k;
    Complex z = 0;
    //Complex power = 0;
    Double LG_v_k;
    Double LG_v=LOG(v);
    //long long v_k;
    Double x=-real(s),y=-imag(s);


    if(what_type_L==-1) {   //i.e. if the Riemann zeta function
        for(k = 0; k < K; k++) {
            //power = real(s)*LOG(v + k) + I*imag(s)*log(1 + (Double)k / (Double)v); 
            LG_v_k=LOG(v+k);
            //power = Complex(x*LG_v_k,y*(LG_v_k-LG_v));
            //z+=exp(power);
            z+= exp(x*LG_v_k)*lcalc_expIt(y*(LG_v_k-LG_v));
        }
    }
    else if(what_type_L!=1) {//if not periodic
        for(k = 0; k < K; k++) {
            //power = real(s)*LOG(v + k) + I*imag(s)*log(1 + (Double)k / (Double)v);
            //power = Complex(x*LG_v_k,y*(LG_v_k-LG_v));
            //z+= dirichlet_coefficient[v + k] * exp(power);
            LG_v_k=LOG(v+k);
            z+= dirichlet_coefficient[v + k]*exp(x*LG_v_k)*lcalc_expIt(y*(LG_v_k-LG_v));
        }
    }
    else { //if periodic
        for(k = 0; k < K; k++) {
            //power = real(s)*LOG(v + k) + I*imag(s)*log(1 + (Double)k/ (Double)v);
            //power = Complex(x*LG_v_k,y*(LG_v_k-LG_v));
            m = (v + k) % period;
            if(m == 0) m = period;
            //z+= dirichlet_coefficient[m] * exp(power);
            LG_v_k=LOG(v+k);
            z+= dirichlet_coefficient[m]*exp(x*LG_v_k)*lcalc_expIt(y*(LG_v_k-LG_v));
        }
    }

    return z;
}


template <class ttype>
Complex L_function <ttype>::
dirichlet_series_via_blfi(Complex s, long long N= -1,Double blfi_interval_length = 1000, Double epsilon = 1.e-9) {
// Function to compute sum_{n=1}^N a_n*n^(-s) to within epsilon using precomputed data.
// blfi_interval_length (optional) is the approximate length of the interval where
// the sum is to be evaluated

    if(N == -1) N = number_of_dirichlet_coefficients;
    if(N > number_of_dirichlet_coefficients && what_type_L!=-1 && what_type_L!=1) N = number_of_dirichlet_coefficients;
    if(N < 1000) return dirichlet_series(s, N);

    //cout << "number_of_dir_coeffs: " << number_of_dirichlet_coefficients << endl;
    //cout << "N: " << N << endl;

    // Some fine-tuning parameters:

    // 1. Kmin is the minimum (or starting) block size;
    // can be adjusted according to one's needs;
    // we chose 50 for Kmin because, typically,
    // the interpolation formula involves that many terms, and
    // so one may as well compute the block directly.
    const int Kmin = 50;

    // 2. v0_factor helps control the starting point of blfi interplation;
    // its precise value can make a noticeable difference in the running time, 
    // but to find an optimal value requires knowledge of the density of points 
    // where the user wishes to evaluate the dirichlet series...  
    // A loose rule of thumb: as v0_factor decreases, the size of the blocks
    // increases (which is good), but also the frequency of redoing the 
    // precomputation increases (which is bad).
    const Double v0_factor = 0.05;

    // 3. beta = beta_factor*tau, where beta determines the density
    // of grid points where blocks in the dirichlet series are to be
    // precomputed. In general, the algorithm is not 
    // very sensetive to beta as long as its value is 
    // neither too small nor too large. 
    // A loose rule of thumb: as beta_factor increases, the interpolation
    // formula requires fewer terms (which is good), but also
    // the frequency of redoing the precomputation increases
    // (which is bad).
    const int beta_factor = 3;

    //XXXXXXXXXXX casting needs to be fixed, example v0= below
    // Precomputation parameters
    static int initialized = 0; // whether precomputation is already done?
    static long long N0; // length of dirichlet series for which the precomputation is done
    static Complex s0; // the s used in in the precomputation
    static Double v_over_K0; // this is approx. the highest frequency in each block; it's used to calculate tau0 below 
    static Double tau0; // parameter used in blfi interpolation
    static Double beta0; //parameter used in blfi interpolation (controls density of grid points where dirichlet series is precomputed)
    static Double lambda0; // parameter used in blfi interpolation (in our version of blfi, lambda0 is determined by tau0 and beta0)
    static Double eps0; // parameter used in blfi interpolation (in our version of blfi, eps0 is determined by tau0 and beta0)
    static long long index0; //index based on which precomputation was done
    static int c0; // parameter used in blfi interpolation (controls the accuracy with which dirichlet series is eval. using blfi) 
    static long long v0; // the starting point of blfi interpolation (the sum from 1 to v0 is eval. directly)
    //static long long length0; // length of the sum to be eval. using blfi (which is the sum from v0+1 up to N0)
    static Double s_range; //this is approx. the range where interpolation can be used  
    static Double c_range; // 2*c_range is approx. the number of grid points where each block is precomputed 

    // A few other useful things
    int c = Int(ceil(-log(epsilon) + (1 - real(s))*LOG(N)));
    //cout << "c = " << c << endl;
    //cout << "s = " << s << endl;
    //cout << "N = " << N << endl;
    //cout << "epsilon = " << epsilon << endl;
    long long K; // to keep track of block size
    long long v; // to keep track of starting point of the block
    long long length; // to keep track of length of remaining sum

    // to keep track of where precomputed data for each block 
    // is stored (useful for retrieving precomputed data later).
    long long center_location = c0;

    // We (re)do the precomputation if one of the following holds:
    //  1. precomputation is not already done
    //  2. length of dirichlet series changes (this condition can be vastly relaxed)
    //  3. real(s) changes (this condition can be relaxed)
    //  4. there'is a big change in c (which controls accuracy)
    //  5. current imag(s) is too far from imag(s0)
   if(my_verbose>1) cout << "#         Entering dirichlet_series_via_blfi, s: "
                         << s << " N: "
                         << N << " blfi_interval_length: "
                         << K << " epsilon: "
                         << epsilon << " c0: "
                         << c0 << " center_location: "
                         << center_location << endl;


   if(initialized == 0 || N != N0 || abs(real(s) - real(s0)) > 1.e-14 || c > c_range || abs(imag(s) - imag(s0)) > s_range) {

       size_t size = 10000000;
       if(initialized == 0) block_value = new Complex[size];
       if(initialized == 1) {
          delete [] block_value;
          block_value = new Complex[size];
       }

       initialized = 1;
       N0 = N;
       s0 = s;

       v_over_K0 = pow(Double(1)*N0, Double(1)/2);
       v_over_K0 = v0_factor * my_min(v_over_K0, 1.1*blfi_interval_length); // the 0.125 can be adjusted according to one's needs

       // Note: the 0.5 in the definition of tau0 below is there
       // because we later center the blocks (centering means
       // that the range of frequencies in each block is in
       // [-v_over_K0/2, v_over_K0/2], rather than [0, v_over_K0],
       // which is useful since it leads to savings during the precomputation.)
       tau0 = 0.5 / v_over_K0;

       beta0 = beta_factor * tau0; // blfi requires beta > tau, so beta_factor must be > 1
       lambda0 = (beta_factor + 1) * tau0 / 2; // lambda = (beta + tau)/2
       eps0 = (beta_factor - 1) * tau0 / 2; // eps = (beta - tau)/2

       index0 = Long(floor(imag(s0) * beta0 / Pi));

       v0 = Long(Kmin*ceil(v_over_K0));
       v0 = my_min(N0, v0);

       // There's a lot of flexibility in choosing the next
       // three parameters (though their choices are not
       // completely independent of each other).
       c0 = Int(my_max(Double(1), c * beta0 / (Pi * (beta0 - tau0) / 2) + 6));
       c_range = c + 3;
       s_range = 5*Pi / beta0;

       K = 1;
       v = v0;
       length = N0 - v0;
       center_location = c0;

       // The precomputation loop
       while(length > 0) {
          K = Long(floor(v / v_over_K0));
          K = my_min(K, length);
          for(int counter = -c0; counter < c0; counter++) {
             block_value[center_location + counter] = get_block_value_directly(real(s0) + Double(1)*(index0 + counter)*I*Pi/ beta0, v + 1, K);
          }
          v = v + K;
          length = length - K;
          center_location += 2 * c0;
       }

       return dirichlet_series_via_blfi(s, N, blfi_interval_length,epsilon);

    }

    Complex z = 0;

    // First, we compute initial sum: sum_{n=1}^v a_n*n^(-s)
    z = partial_dirichlet_series(s, 1, v0);

    // Then we compute remaining sum: sum_{n=v+1}^N a_n*n^(-s), in blocks, via blfi
    K = 1;
    v = v0;
    length = N0 - v0;
    center_location = c0;

    while(length > 0) {
          K = Long(floor(v / v_over_K0)); // block size
          K = my_min(K, length);
          z = z + dirichlet_series_block_blfi(s, v + 1, K, c, c0, center_location, index0, tau0, beta0, lambda0, eps0);
          v = v + K;
          length = length - K;
          center_location += 2 * c0;
          //cout << length << " z: " << z << endl;
    }

    //cout << endl << "Initial sum = " << (float)(100 * (Double)v0 / (Double) N) << "%, block total = " << floor(center_location/(2*c0)) << ", c = " << c << ", c0 = " << c0 << endl;
    return z;
}


template <class ttype>
Complex L_function <ttype>::
dirichlet_series_block_blfi(Complex s, long long v, long long K, int c, int c0, long long center_location, long long index0, Double tau, Double beta, Double lambda, Double eps) {
// Function to compute the block sum_{k=0}^{K-1} a_{v+k}*(v+k)^(-s) via blfi

   if(my_verbose>1) cout << "#         Entering dirichlet_series_block_blfi with: "
                         << s << " "
                         << v << " "
                         << K << " "
                         << c << " "
                         << c0 << " "
                         << center_location << " "
                         << index0 << " "
                         << tau << " "
                         << beta << " "
                         << lambda << " "
                         << eps << endl;

   if(v < 100 || K  < 50) return partial_dirichlet_series(s, v, v + K -1);

   Complex z = 0;
   Double Pi_over_beta = Pi/beta;
   long long index1 = Long(floor(imag(s)/Pi_over_beta));
   int difference = index1 - index0;
   Double u0 = index1 * Pi_over_beta - imag(s);

   // For efficiency reasons, the normalization factor below
   // is computed here (rather than in blfi_kernel(..),
   // which is where it belongs on paper).
   Double normalization = c / sinh(Double(c));

   // The next two while loops are the blfi interplation formula:

   Double u = u0;
   Complex exp_factor = lcalc_expIt(tau * u);
   Complex exp_inc = lcalc_expIt(tau * Pi_over_beta);
   int counter = 0;

   long long center_location_plus_difference=center_location + difference;

   Double c_over_eps = c/eps;

   while(u < c_over_eps) {
      // z = z + block_value[center_location + difference + counter] * exp(I * tau * u) * sinc(lambda * u) * blfi_kernel(u, c, eps);
      //z+=block_value[center_location_plus_difference + counter] * exp_factor * sinc(lambda*u) * blfi_kernel((c - eps*u)*(c+eps*u));
      z+=block_value[center_location_plus_difference + counter] * exp_factor * sinc(Double(lambda*u)) * blfi_kernel(c*c - eps*eps*u*u);
      u+=Pi_over_beta;
      exp_factor *= exp_inc; //XXXXXXXX errors accumulate- careful
      counter++;
   }

   if(counter > c0 - 1) {
      cout << "blfi error! block_value array is being called out of range" << endl;
      exit(1);
   }

   u = u0 - Pi_over_beta;
   exp_factor = lcalc_expIt(tau * u);
   exp_inc = conj(exp_inc);
   counter = -1;


   while(u > -c_over_eps) {
      // z = z + block_value[center_location + difference + counter] * exp(I * tau * u) * sinc(lambda * u) * blfi_kernel(u, c, eps);
      z+=block_value[center_location_plus_difference + counter] * exp_factor * sinc(Double(lambda*u)) * blfi_kernel(c*c - eps*eps*u*u);
      //z+=block_value[center_location_plus_difference + counter] * exp_factor * sinc(lambda*u) * blfi_kernel((c - eps*u)*(c+eps*u));
      u-=Pi_over_beta;
      exp_factor *= exp_inc; //XXXXXXXX errors accumulate- careful
      counter--;
   }

   if(counter < -c0 + 1) {
      cout << "blfi error! block_value array is being called out of range" << endl;
      exit(1);
   }

   // Last, we multiply z by the normalization factors in the blfi formula
   //z = z * lambda / beta * exp(-I * imag(s) * LOG(v));
   z*=normalization*lambda / beta * lcalc_expIt(-imag(s) * LOG(v));
   return z;
}










#endif
