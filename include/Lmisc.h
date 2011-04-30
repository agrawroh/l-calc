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


#ifndef Lmisc_H
#define Lmisc_H

#include "Lglobals.h"
#include<vector>

vector<string> &split(const string &s, char delim, vector<string> &elems);
vector<string> split(const string &s, char delim);

vector<Double> split_Double(const string &s, char delim);
vector<Double> &split_Double(const string &s, char delim, vector<Double> &elems);


template<class T> int sn(T x);

inline Double LOG(int n)
{
    if(n>number_logs) extend_LG_table(n);
    return LG[n];
}



//inline Double two_inverse_sqrt(double n) // Used for long
inline Double two_inverse_sqrt(int n)
{
    //int m=Int(n);
    //if(m>number_sqrts) extend_sqrt_table(m);
    //return two_inverse_SQUARE_ROOT[m];

    if(n>number_sqrts) extend_sqrt_table(n);
    return two_inverse_SQUARE_ROOT[n];
}

template <class type1, class type2>
inline type1 my_max(type1 a, type2 b)
{
     type1 r,b_converted;
     r=a;
     b_converted = (type1) b;
     if(b_converted > a) r= b_converted;
     return r;
}

template <class type1, class type2>
inline type1 my_min(type1 a, type2 b)
{
     type1 r,b_converted;
     r=a;
     b_converted = (type1) b;
     if(b_converted < a) r= b_converted;
     return r;
}

inline Double lcalc_cos(Double x)
{


    Double t2=x*one_over_twoPi;
    Double t=t2-floor(t2);


    int n=Int(t*cos_taylor_arraysize);


    //Double Y = 2*Pi*(t-(n+.5)/cos_taylor_arraysize);
    Double Y = twoPi*t-twoPi_over_cos_taylor_arraysize*(n+.5);


    //Double debug_cos;

    Double *p;
    p=&cos_taylor[number_cos_taylor_terms*n];

    if(DIGITS<17){

        //debug_cos=abs(((p[3]*Y + p[2])*Y + p[1])*Y + p[0]-cos(x));
        //if(debug_cos>1e-10){
            //cout <<"center on: " << twoPi_over_cos_taylor_arraysize*(n+.5) << endl;
            //cout <<"Y: " << Y << endl;
            //cout <<"n: " << n << endl;
            //cout << "lcalc_cos(" << x << ") -cos= " << ((p[3]*Y + p[2])*Y + p[1])*Y + p[0]-cos(x) << endl;
            //cout << "lcalc_cos term: " << 0 << " " << p[0] << endl;
            //cout << "lcalc_cos term: " << 1 << " " << p[1] << endl;
            //cout << "lcalc_cos term: " << 2 << " " << p[2] << endl;
            //cout << "lcalc_cos term: " << 3 << " " << p[3] << endl;
        //}


        return (((p[3]*Y + p[2])*Y + p[1])*Y + p[0]);

    }

    n=number_cos_taylor_terms-1;
    t=p[n];
    do{
        //cout << "lcalc_cos term: " << n << " " << p[n] << endl;
        n--;
        t=t*Y+p[n];
    }while(n>0);

    //Double debug_cos=t-cos(x);
    //if(debug_cos>1e-10){
        //cout <<"center on: " << twoPi_over_cos_taylor_arraysize*(n+.5) << endl;
        //cout <<"center on: " << twoPi_over_cos_taylor_arraysize*(n+Double(1)/2) << endl;
        //cout <<"Y: " << Y << endl;
        //cout <<"n: " << n << endl;
        //cout << "lcalc_cos(" << x << ") -cos= " << ((p[3]*Y + p[2])*Y + p[1])*Y + p[0]-cos(x) << endl;
        //cout << "lcalc_cos term: " << 0 << " " << p[0] << endl;
        //cout << "lcalc_cos term: " << 1 << " " << p[1] << endl;
        //cout << "lcalc_cos term: " << 2 << " " << p[2] << endl;
        //cout << "lcalc_cos term: " << 3 << " " << p[3] << endl;
        //cout << "lcalc_cos term: " << 4 << " " << p[4] << endl;
        //cout << "lcalc_cos term: " << 5 << " " << p[5] << endl;
    //}

    return t;

}

#define lcalc_sin(x) lcalc_cos(x-Pi/2)

inline void lcalc_cis(Double& c, Double& s, Double t){
    c=lcalc_cos(t);
    s=lcalc_sin(t);
}


inline void lcalc_cis_b(Double& c, Double& s, Double t){
    long long m=Long(floor(t/Pi));
    c=lcalc_cos(t);
    s=(1.-2*(m%2))*sqrt(1.-c*c);
}

inline Complex lcalc_expIt(Double t){
    Double c=lcalc_cos(t), s=lcalc_sin(t);
    return Complex(c,s);
}

#ifdef USE_MPFR
    #define lcalc_exp(x) exp(x)
#else
inline Double lcalc_exp(Double z){
    return exp(z);
}

inline Complex lcalc_exp(Complex z){
    Double t=imag(z),r=exp(real(z));
    Double c=r*lcalc_cos(t), s=r*lcalc_sin(t);
    return Complex(c,s);
}
#endif

inline Complex lcalc_expIt_b(Double t){
    Double c=lcalc_cos(t), s;
    long long m=Long(floor(t/Pi));
    s=(1.-2*(m%2))*sqrt(1.-c*c);
    return Complex(c,s);
}

inline int my_sgn(Double t){
   if(t<0) return (-1);
   return 1;
}

//compute sinc(x)=sin(x)/x


//XXXX needs to be made multiprecision
inline Complex sinc(Complex u){
    Complex ans=1;
    if(abs(u) > sin_tol) ans=sin(u)/u;
    //if(abs(u)<= sin_tol){
    else{
        Complex u2=u*u,temp=u2;
        for(int j=1;j<sin_terms;j++){
            ans+=sin_cof[j]*temp;
            temp*=u2;
        }
    }

    return ans;
}

//XXXX needs to be made multiprecision
inline Double sinc(Double u){
    Double ans=1;
    if(abs(u) > sin_tol) ans=lcalc_sin(u)/u;
    //if(abs(u)<= sin_tol){
    else{
        Double u2=u*u,temp=u2;
        for(int j=1;j<sin_terms;j++){
            ans+=sin_cof[j]*temp;
            temp*=u2;
        }
    }

    return ans;
}

// template <class ttype>
// inline ttype sinc(ttype u){
//     ttype ans=1;
//     if(abs(u) > sin_tol) ans=sin(u)/u;
//     //if(abs(u)<= sin_tol){
//     else{
//         ttype u2=u*u,temp=u2;
//         for(int j=1;j<sin_terms;j++){
//             ans+=sin_cof[j]*temp;
//             temp*=u2;
//         }
//     }
// 
//     return ans;
// }

template <class ttype>
inline ttype Horner(vector<Double> b, ttype x )
{

    //cout << "# enter horner with: " << x << endl;
    vector<Double>::size_type deg=b.size()-1;

    ttype result = b[deg];
    for(int j=deg-1; j >= 0 ; --j){
        //cout << "HORNER: " << j << " " << b[j] << endl;
        result = result * x + b[j];
    }
    return result;
}


#endif
