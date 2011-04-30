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


#include "L.h"
#include "Lriemannsiegel.h"



// ZETA FUNCTION
Complex Zeta(Complex s, const char *return_type) {

  Complex L_value;

  L_value= siegel(s);


  //this returns Zeta(s)
  if (!strcmp(return_type,"pure"))
      return L_value;

  //returns Zeta(s) rotated to be real on critical line
  else if (!strcmp(return_type,"rotated pure"))
      return L_value*exp(I*(imag(log_GAMMA(s/2)) - (imag(s)/2)*log(Pi)));

  return L_value;

}


// RIEMANN-SIEGEL ZETA EVALUATION
Complex siegel(Complex s) {

    Double theta, result = 0;
    Double errorTerms;
    Double tau,p, t;
    int N, n;

    Double z;


    tau = sqrt(abs(imag(s))/(2*Pi));
    N = Int(floor(tau));
    p = tau - N;

    t = imag(s);

    theta = imag(log_GAMMA(s/2)) - (imag(s)/2)*log(Pi);

//cout <<setprecision(16);
//cout << t << " rs theta: " << exp(-I*theta) << endl;


    if (my_verbose>1) cout <<"Main sum is " << N << " terms long" << endl;
    two_inverse_sqrt(N); LOG(N); //call to extend two_inverse_SQUARE_ROOT and LG tables
                                 //ahead of time as far as needed. More efficient this way

    if(N>7000){ //parallelization only makes sense if the length is substantial enough to make up for overhead
        #pragma omp parallel for reduction(+: result)
        //for (n = 1; n <=N; n++){
        for (n = N; n >0; n--){
          //result += two_inverse_SQUARE_ROOT[n]*cos(t*LG[n]-theta);
          #ifndef USE_DOUBLE
              __builtin_prefetch(&two_inverse_SQUARE_ROOT[n-1],0,0);
              __builtin_prefetch(&LG[n-1],0,0);
          #endif
          result += two_inverse_SQUARE_ROOT[n]*lcalc_cos(theta-t*LG[n]);
        }
    }
    else{
        //for (n = 1; n <=N; n++){
        for (register int n = N; n >0; n--){ //checking n>0 is faster than checking n<=N
          #ifndef USE_DOUBLE
              __builtin_prefetch(&two_inverse_SQUARE_ROOT[n-1],0,0);
              __builtin_prefetch(&LG[n-1],0,0);
          #endif
          result += two_inverse_SQUARE_ROOT[n]*lcalc_cos(theta-t*LG[n]);
        }
/*

        __builtin_prefetch(&result,1,3);
        for (n = 1; n<=N-4; n+=4){ //checking n>0 is faster than checking n<=N
          __builtin_prefetch(&two_inverse_SQUARE_ROOT[n+4],0,0);
          __builtin_prefetch(&LG[n+4],0,0);
          result += two_inverse_SQUARE_ROOT[n]*lcalc_cos(theta-t*LG[n]);
          result += two_inverse_SQUARE_ROOT[n+1]*lcalc_cos(theta-t*LG[n+1]);
          result += two_inverse_SQUARE_ROOT[n+2]*lcalc_cos(theta-t*LG[n+2]);
          result += two_inverse_SQUARE_ROOT[n+3]*lcalc_cos(theta-t*LG[n+3]);
        }
        for (; n<=N; n++)
            result += two_inverse_SQUARE_ROOT[n]*lcalc_cos(theta-t*LG[n]);
*/
    }


    z=p-.5;

    max_n=N;


/*
    Double remainder=0;
    Double one_over_tau=1./tau,r=1.,C;

    Double zz[136];
    zz[0]=1;
    for(int n=1;n<=135;n++)
        zz[n]=zz[n-1]*z; //stores powers of z;

    int j=0;
    do{
        C=0.;
        int shift = j%2;
        for(int n=0;n<=71;n++){
            C+=rs_remainder[j][n]*zz[2*n+shift];
            //cout << j << " " << n << " " << rs_remainder[j][n] << " " << zz[2*n+shift] << " " << C <<endl;
        }
        remainder+=C*r;
        cout << "with rs terms:" << j << " " << (result+pow(tau,-Double(1)/2)*(1-2*((N-1)%2))*remainder)*exp(-I*theta) <<endl;
        r*=one_over_tau;
        j++;
    }while(1e10*r>tolerance&&j<=39);

*/

    errorTerms=rs_remainder_terms(z,tau);
    errorTerms*=pow(tau,-Double(1)/2)*(1-2*((N-1)%2));


    //cout << setprecision(DIGITS3);
    //cout << "rs main term: " << result*exp(-I*theta) << endl;
    //cout << "rs remainder term: " << errorTerms * exp(-I*theta) << endl;

    result += errorTerms;


    // result is Z.  Now we rotate to get zeta
    return (result* exp(-I*theta));



}

/*
Complex rs_main_sum(int N,Double theta){

    cout << "calling rs_main_sum with N,theta: " << N << " " << theta << endl;
    if(N==0) return 0;

    int N0 = N/210;
    N0=210*N0;
    Complex tail=0;
    for(int n = N0+1;n<=N;n++)
        tail+= two_inverse_SQUARE_ROOT[n]*lcalc_cos(theta-t*LG[n]);

    if(N0==0) return tail;
    //have a lookup table here (i.e. if N<=210)

    //sieve out by 2,3,5,7
    Complex s2,s3,s5,s7,s23,s25,s27,s35,s37,s57,s235,s237,s257,s357,s2357;

    s2=

    return (tail+rs_main_sum()*);

}

*/

Double rs_remainder_terms(Double z, Double tau){

    Double remainder=0;
    Double one_over_tau=1./tau,r=1.,C;


    Double zz[144];
    zz[0]=1;
    for(int n=1;n<=143;n++)
        zz[n]=zz[n-1]*z; //stores powers of z;

    int j=0,n;
    do{
        C=0.;
        int shift = j%2;

        /*
        for(int n=0;n<=71;n++){
            C+=rs_remainder[j][n]*zz[2*n+shift];
        }
        Break this into 4 blocks and check size of terms after each block
        Compromise between checking size after each term and not doing any checking
        */
        for(n=0;n<=20;n++){
            C+=rs_remainder[j][n]*zz[2*n+shift];
        }
        if(my_norm(rs_remainder[j][n]*zz[2*n+shift])>tolerance_sqrd){
            for(n=21;n<=35;n++){
                C+=rs_remainder[j][n]*zz[2*n+shift];
            }
            if(my_norm(rs_remainder[j][n]*zz[2*n+shift])>tolerance_sqrd){
                for(n=36;n<=54;n++){
                    C+=rs_remainder[j][n]*zz[2*n+shift];
                }
                if(my_norm(rs_remainder[j][n]*zz[2*n+shift])>tolerance_sqrd){
                    for(n=55;n<=71;n++){
                        C+=rs_remainder[j][n]*zz[2*n+shift];
                    }
                }
            }
        }

        remainder+=C*r;

        //cout << "rs remainder:" << j << " " << C << " " << C*r <<endl;
        r*=one_over_tau;
        j++;
    }while(r>tolerance&&j<=39);

    return remainder;
}

