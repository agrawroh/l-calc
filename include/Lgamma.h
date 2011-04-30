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



#ifndef Lgamma_H
#define Lgamma_H

#include "Lglobals.h"
#include "Lmisc.h"
#include <iomanip>          //for manipulating output such as setprecision
#include <iostream>         //for ostrstream
#include <cstring>
//using namespace std;





template <class ttype>
precise(ttype) log_GAMMA (ttype z,int n=0); //return ttype since I also want to
                                            //use it for N(T) which is real.
template <class ttype>
Complex GAMMA (ttype z);

template <class ttype,class ttype2>
Complex GAMMA (ttype z,ttype2 delta);

//computes normalized incomplete gamma function G(z,w)=w^{-z} Gamma(z,w)  =====
template <class ttype>
Complex inc_GAMMA (ttype z,ttype w, const char *method="temme");  //computes G(z,w)

template <class ttype>
ttype asympt_GAMMA (ttype z,ttype w);  //computes G(z,w) via asymptotic series

template <class ttype>
ttype cfrac_GAMMA (ttype z,ttype w);  //computes G(z,w) via continued fraction

template <class ttype>
ttype comp_inc_GAMMA (ttype z,ttype w);  //computes g(z,w)

//=====================================================================

inline Double real(Double z){
    return z;
}
inline Double imag(Double z){
    return 0;
}

//used in Temme's asymptotics of the incomplete gamma function
Complex Q(Complex z, Complex w);
Complex erfc(Complex z);
Complex erfc2(Complex z);

template <class ttype>
Complex gamma_sum(Complex s, int what_type,ttype *coeff, int N, Double g, Complex l, Double Q, Long Period, Complex delta=1, const char *method="temme");
//computes a sum of G(z,w)'s (as in (3.3.20)).

//compute the nth derivative of log(GAMMA(z))
//n=0 gives log_GAMMA(z).
//Only bothers to handle Re(z)<0 if n=0.
template <class ttype>
precise(ttype) log_GAMMA (ttype z,int n=0)
{
    int M;
    precise(ttype) log_G,r,r2,y;

    Double xx=real((Complex)z),yy=imag((Complex)z);
    //if(xx<0) xx=-xx;

    //check if z is zero. If so return nan
    if(my_norm(z)<tolerance_sqrd){
        return numeric_limits<float>::quiet_NaN();
    }

    //else if Re(z)<0 return log(Gamma(z+1)/z)
    if(xx<0){
       if(n==0){
         return(log_GAMMA(z+1)-log(z));
       }
       else{
           cout << "error in log_GAMMA: derivative called with Re(z)<0" << endl;
       }
    }


    Double x;
    int i,m;


    //assume the remainder stopping at the mth term is roughly bernoulli(2m)/(z+M)^(2m).
    //Now bernoulli(2m) is roughly (2m)!/(2Pi)^(2m). So remainder is more or less
    //(2m/(2ePi(z+M))^(2m). Setting 2m = Digits, we see that we need m/(ePi(z+M))
    //to be roughly 1/10 in size, i.e. |z+M| should be roughly 10*m/(ePi)=10*Digits/(2ePi).
    //squaring, gives us how large M should be.

    //n==0 leads to |z+M| >10*2m/(2*Pi*e) with 2m=Digits. However, when we differentiate
    //we end up with extra powers of (z+M) in the denominators, but extra (2m-1)...(2m+n-2)
    //in the numerators. So assume further that |z+M| > 2m+n = Digits+n

    if(n==0){
        //.343 = 100./(4*Pi*Pi*exp(2.))
        if((xx*xx+yy*yy)> .343*DIGITS*DIGITS) M=0;
        else M=Int(ceil(sqrt((DIGITS*DIGITS*.343-yy*yy))-xx+1));
    }
    else{
        if((xx*xx+yy*yy)> .343*(DIGITS+n)*(DIGITS+n)) M=0;
        else M=Int(ceil(sqrt(((DIGITS+n)*(DIGITS+n)*.343-yy*yy))-xx+1));
    }

    if(n==0){
       log_G=log_2Pi/2+(z+M-.5)*log(z+M)-(z+M);
    }
    else if(n==1)
       log_G=log(z+M)-.5/(z+M);
    else{
       r=Double(1);
       for(m=1;m<=n-1;m++){
          r=-r*m/(z+M);
       }
       log_G=log_G-r/(n-1)-.5*r/(z+M);
    }

    r=Double(1);
    for(m=1;m<=n;m++){
       r=-r*m/(z+M);
    }
    r=r/(z+M);

    r2=1/((z+M)*(z+M));
    m=2;
    x=my_norm(r);
    do{
        y=bernoulli[m]*r;
        log_G+=y/(m*(m-1));

        r*=(m+n-1)*(m+n)*r2/((m-1)*m);
        m+=2;
    }while(m<=DIGITS&&tolerance_sqrd<x*my_norm(y));
    //}while(m<=DIGITS);

    for(m=0;m<=M-1;m++){
       if(n==0){
           log_G-=log(z+m); //XXX might be faster to multiply and take one log,
                                 //but careful about overflow errors

       }
       else{
           r=Double(1);
           for(i=1;i<=n;i++){
              r=-r*i/(z+m);
           }
           log_G+=r/n;
       }
    }

    return log_G;
}

// computes GAMMA(z) delta^(-z)
// assumes the ttype is at least real (i.e., Double or Complex)
// I compute this as exp(log(GAMMA(z))
// since there are 1/2 as many terms and the error analysis
//is easier.
//Temme's book, chapter 3 gives the error analysis for
//log(GAMMA(z)), z>0, and mentions Spira's estimate for complex values of z.
//He also suggests using the reflection formula for negative values.

//unlike log_GAMMA, I assume that Re(z)>=0. WHile I could simplify a bit
// and call log_GAMMA, I prefer to keep this since it is slightly shorter
template <class ttype,class ttype2>
Complex GAMMA (ttype z1, ttype2 delta)
{

    precise(ttype) z=z1;
    Complex log_G;
    Double xx=real((Complex)z),yy=imag((Complex)z);
    if(xx<0) xx=-xx;

    if(z==last_z_GAMMA) log_G=last_log_G;

    else{
        int M;
        precise(ttype) r,r2,y;

        Double x;
        int m;

        //assume the remainder stopping at the mth term is roughly bernoulli(2m)/(z+M)^(2m).
        //Now bernoulli(2m) is roughly (2m)!/(2Pi)^(2m). So remainder is more or less
        //(2m/(2ePi(z+M))^(2m). Setting 2m = Digits, we see that we need m/(ePi(z+M))
        //to be roughly 1/10 in size, i.e. |z+M| should be roughly 10*m/(ePi)=10*Digits/(2ePi).
        //squaring, gives us how large M should be.

        //.343 = 100./(4*Pi*Pi*exp(2.))
        if((xx*xx+yy*yy)> .343*DIGITS*DIGITS) M=0;
        else M=Int(ceil(sqrt((DIGITS*DIGITS*.343-yy*yy))-xx+1));

//cout << " GAMMA M : " << M << endl;
        log_G=log_2Pi/2+(z+M-.5)*log(z+M)-(z+M);

        r2=(z+Double(M));
        r=r2*r2;
        m=2;
        x=my_norm(r2);
        do{
            y=bernoulli[m]/(m*(m-1)*r2);
            log_G+=y;
            r2*=r;
            m+=2;
        }while(m<=DIGITS&&tolerance_sqrd<x*my_norm(y));

        // could combine into a single log, but careful about overflow
        //for(m=0;m<=M-1;m++) log_G-=log(z+m);

        //... so to control for overflow, break up 0 to M-1 into blocks of length 10,
        //i.e. do up to ten multiplications, and then a log. Our z's are not going to be
        //bigger than 10^30 or so, so 10 multiplications won't overflow.

        bool do_log=false;
        r=Double(1);
        for(m=0;m<=M-1;m++){
            do_log=true;
            r*=(z+m);
            if((m+1)%10==0){
                log_G-=log(r);
                r=Double(1);
                do_log=false;
            }
        }
        if(do_log)log_G-=log(r);


        last_log_G=log_G;
        last_z_GAMMA=z;
    }




    log_G-=z*log(delta);


    //if(my_verbose>2) cout << "GAMMA function("<<z<<","<<delta<<")= " << exp(log_G) << endl; //XXXXX
    //return exp(log_G);
    return lcalc_exp(log_G);

}

template <class ttype>
Complex GAMMA (ttype z){
    return GAMMA(z,Double(1));
}


//XXXXXXXXXX
//causes problems if z is an integer <=0

//not meant to be used for very negative values of Re(z)
//since a recursion is used to get z positive.
//computes G(z,w), so there's an extra w^(-z) factor.
template <class ttype>
Complex inc_GAMMA (ttype z,ttype w, const char *method="temme")
{

    Complex G;
    Double h,x,y;
    ttype tmp=z+1;


 //cout << setprecision(16);

    if(my_verbose>5) cout << "#                        inc_GAMMA called. G(" << z<< " , " << w << ")" << endl;

    //inc_GAMMA(z,w) = GAMMA(z)-comp_inc_Gamma(z,w), and at z=0, the rhs has a NAN - NAN, and
    //and cancellation problems nearby. So we avoid this by calling the continued fraction
    //for inc_GAMMA(z,w) near z=0.
    if(my_norm(z)<.01){
        if(my_verbose>5) cout << "#                        calling cfrac_GAMMA from here" << endl;
        return cfrac_GAMMA(z,w);
    }

    if(my_norm(z-1)<tolerance_sqrd){
        return(lcalc_exp(-w)/w);
    }


    if(real((Complex)z)<=0){
        //XXXX tmp equals z+1. compiler complains when I try inc_GAMMA(z+1,w)
        //when mpfr is used. Regards z+1 here as a binary expression rather than a
        //Complex or Double
            return ((w*inc_GAMMA(tmp,w)-lcalc_exp(-w))/z);
        }

        y=my_norm(w);
        Double tmp_norm = my_norm(z);

        //if((my_norm(z)>100&&y>my_norm(z*1.01))||!strcmp(method,"continued fraction"))
        if((tmp_norm>100&&y>tmp_norm*1.02)||!strcmp(method,"continued fraction"))
        {
            if(my_verbose>5) cout << "#                        calling cfrac_GAMMA from this spot" << endl;
            return cfrac_GAMMA(z,w);
        }

        x=imag((Complex)z);
        if(x<0)x=-x;
        //x=x-sqrt(x);
        x*=.99; //temme works better here than x-sqrt(x), at least asymptotically

        if(y<1600||y<x*x)  //XXXX roughly if |w|<|imag(z)| , or if |w| is small
                       //XXXX NOTE:
                       //XXXX one apparent problem with using g(z,w)
                       //XXXX when Re w is large, is that g(z,w) is
                       //XXXX close to GAMMA(z)*w^(-z).
                       //XXXX So subtracting to get G(z,w) loses me precision in
                       //XXXX the 'non leading zero' digits.
                       //XXXX On the other hand... if we count the leading zeros,
                       //XXXX then we have our full precision! And we should count them!
                       //XXXX When we sum the terms in (33) of my Computation methods paper,
                       //XXXX the initial ones start
                       //XXXX off comparatively large... so the leading zeros count in
                       //XXXX the later terms (i.e. even if we had more precision for
                       //XXXX the later G(z,w) terms it would be lost when we summed
                       //XXXX against the earlier terms in (33)).
                       //XXXX note... we can use this to *slightly* speed up
                       //XXXX the computation in that we don't always need to ask g(z,w)
                       //XXXX for too many digits since only the first few will be useful.
                       //XXXX i haven't implemented this, though.
        {
            last_z=z;
            last_w=w;

            last_comp_inc_GAMMA=comp_inc_GAMMA(z,w);
        G=GAMMA(z,w)-last_comp_inc_GAMMA;
        if(my_verbose>5){
            cout << "#                        last_comp_inc_GAMMA("<<z<<","<<w<<")= " << last_comp_inc_GAMMA  << endl; //XXXXX
            cout << "#                        series incGAMMA("<<z<<","<<w<<")= " << G << endl; //XXXXX
        }
        return G;
    }

    //XXXXX condition here and above should depend on precision
    //if(y>=10&&y<1600)
    //{
        //return cfrac_GAMMA(z,w);
    //}

    //use temme's uniform asymptotics
    if(!strcmp(method,"temme")&&y<my_norm(z*1.2)) 
    {
        G=Q(z,w)*GAMMA(z,w);
        if(my_verbose>5) cout << "#                        temme GAMMA("<<z<<","<<w<<")= " << G << endl; //XXXXX
        return G;
    }

    //asymptotic series... should never be called, but doesn't hurt to have it
    x=abs(z);
    h=(DIGITS+2)*2.3026+1; //XXXX+2 is to be safe
    x=x+h+sqrt(h*h+4*h*x); //XXXX (3.3.48) of my thesis generalized for DIGITS
    if(y>x*x)                //XXXX if abs(w) > (3.3.48)
    {
        return asympt_GAMMA(z,w);
    }


    // this is called if abs(z)<10 and w is around size 40
    //return ((w*inc_GAMMA(z+1,w)-exp(-w))/z);

    return ((w*inc_GAMMA(tmp,w)-lcalc_exp(-w))/z);


}


template <class ttype>
ttype cfrac_GAMMA (ttype z,ttype w)  //computes G(z,w) via continued fraction
{

        ttype G;

        //if(my_verbose>3) cout << "called cfrac_GAMMA("<<z<<","<<w<<")"<< endl;
        /*
        //old code commented out. Old code used backward recursion iterated until
        //convergence to within tolerance was satisfied. This requires extra steps
        //compared to forward recursion.

        //Complex tmp=1;
        int m;
        int M=2;
        bool escape=false;


        for(int m=M;m>=1;m--)
            //tmp=w+(m-z)/(1+m/tmp);
            tmp=w+tmp*(m-z)/(tmp+m);

        do{
            M=2*M;
            G=1;

            for(int m=M;m>=1;m--)
                G=w+G*(m-z)/(G+m);
            if(my_norm(1-tmp/G)<tolerance_sqrd)escape=true;
            else tmp=G;
            cout << "cfrac with number terms: "<<M << endl;
        }while(!escape&&M<100000000);
        */

        // newer version uses forward recursion. We also avoid the
        // expense of dividing, keeping track of the numerators
        // and denominators, dividing after convergence is achieved.

        int n;
        //ttype P1=1.,P2=w,Q1=0.,Q2=1., P3,Q3;
        ttype P1=Double(1),P2=w,Q1=0.,Q2=Double(1), c1=-z;

        n=0;
        do{

            /* this is forward recursion. however, the implementation below is faster
            n++;
            P3=P2+(n-z)*P1;
            Q3=Q2+(n-z)*Q1;
            P1=P2;P2=P3;
            Q1=Q2;Q2=Q3;

            P3=w*P2+n*P1;
            Q3=w*Q2+n*Q1;
            P1=P2;P2=P3;
            Q1=Q2;Q2=Q3;
            */


            //each block here is equivalent to the above. I unroll and repeat four times.
            n++; c1+=Double(1);
            P1=c1*P1+P2, Q1=c1*Q1+Q2;
            P2=w*P1+n*P2; Q2=w*Q1+n*Q2;

            n++; c1+=Double(1);
            P1=c1*P1+P2, Q1=c1*Q1+Q2;
            P2=w*P1+n*P2; Q2=w*Q1+n*Q2;

            n++; c1+=Double(1);
            P1=c1*P1+P2, Q1=c1*Q1+Q2;
            P2=w*P1+n*P2; Q2=w*Q1+n*Q2;

            n++; c1+=Double(1);
            P1=c1*P1+P2, Q1=c1*Q1+Q2;
            P2=w*P1+n*P2; Q2=w*Q1+n*Q2;


            //to prevent overflow, XXX check this
            if(n%8==0&&(real(P2)>1.e40||real(P2)<-1.e40||imag(P2)>1.e40||imag(P2)<-1.e40)){
                P1*=1.e-40;
                P2*=1.e-40;
                Q1*=1.e-40;
                Q2*=1.e-40;

            }

            //cout << "cfrac1 : " << n<< " " <<my_norm(Q2*P1-P2*Q1) << " " << my_norm(Q2*P1*tolerance)<<endl;
            //cout << "cfrac called n equals: " << n << " P2/Q2 " << P2/Q2 << endl;
        } while(my_norm(1-P2*Q1/(Q2*P1))>tolerance_sqrd&&n<1000000);
        //if (my_verbose>3) cout << "cfrac called n equals: " << n << " " << abs(z/w) << endl;

        G=P2/Q2;


        if(n>999999){
             cout << "Continued fraction for G(z,w) failed to converge. z = "
             << z << "  w = " << w << endl;
             exit(1);
        }


        return lcalc_exp(-w)/G; //XXXXX mpfr should I precise(ttype) w ?

}

template <class ttype>
ttype asympt_GAMMA (ttype z,ttype w)  //computes G(z,w) via asymptotic series
{

        if(my_verbose>5) cout << "#                        called asympt_GAMMA("<<z<<","<<w<<")"<< endl;
        ttype G=0;
        ttype r=Double(1);
        int j=0;
        do
        {
            G=G+r;
            r=r*(z-1-j)/w;
            //cout << j << " " << G << " " << r <<" " << tolerance <<endl; //XXXXX
            j++;
        }while(my_norm(r)>tolerance_sqrd);
        G=G*lcalc_exp(-w)/w;
        //cout << "asymptotics GAMMA("<<z<<","<<w<<")= " << G << endl; //XXXXX
        return G;
}


template <class ttype>
ttype comp_inc_GAMMA (ttype z,ttype w)  //computes g(z,w)
{

    ttype g;
    Double t;

    //Complex u=w*w;
    //u=u*u;

    //if(my_verbose>3) cout << "called comp_inc_GAMMA("<<z<<","<<w<<")"<< endl;

    t=my_norm(w/z);


    if(t>.9801 || my_norm(w)<.36){

        //double minus_Re_z = -real((Complex)z);

        //below and do loop is just a fancy way to add 1+w/u+w^2/(u*(u+1))+w^3/(u*(u+1)*(u+2))+....
        //designed to cut back on divisions and to do several terms at once.
        ttype r=Double(1), u=z+1,w2=w*w, w3=w2*w,w4=w2*w2,t1=1/(u*(u+1)*(u+2)*(u+3)*(u+4)*(u+5)),w6=w3*w3;
        ttype c1=10 + w, c2= 10 * w + w2 + 35, c3= 35 * w + w3 + 9 * w2 + 50, c4= 26 * w2 + 50 * w + w4 + 7 * w3 + 24,
        c5= 24 * w + 4 * w4 + 24 * w2 + w4 * w + 12 * w3;
        g=0.;
        do{

            g+=r*t1*(u+5)*(((((u + c1)*u + c2) * u + c3) * u + c4 )* u + c5);
            r*=t1*w6;
            u+=6;
            t1=1/(u*(u+1)*(u+2)*(u+3)*(u+4)*(u+5));

            //g+=r*t1*(u+2)*(u*(u+1)+w*(u+1)+w*w);
            //r*=t1*w3;
            //u+=3;
            //t1=1/(u*(u+1)*(u+2));
            //g+=r*t1*(u+2)*(u*(u+1)+w*(u+1)+w*w);
            //r*=t1*w3;
            //u+=3;
            //t1=1/(u*(u+1)*(u+2));

            //g+=r;
            //r*=w/u;
            //g+=r;
            //r*=w/(u+1);
            //g+=r;
            //r*=w/(u+2);
            //u+=3;



        //}while(my_norm(r)>tolerance_sqrd||real((Complex)z)<=-m);
        //}while(my_norm(r)>tolerance_sqrd||m<minus_Re_z); //isn't the second part vacuous XXXXXXXXXX
        }while(my_norm(r)>tolerance_sqrd);
        g*=lcalc_exp(-w)/z;
    }

    else{

        /*-----------------------------------------------------------
        // old iterated backward recursion replaced by forward recursion

        int M=2;
        Complex tmp=1;
        bool escape=false;


        for(int m=M;m>=1;m--){
            if(m%2==0){
               //tmp= z+m-1+(m/2)*w/tmp;
               tmp= z+m-1+m*.5*w/tmp;
            }
            else{
                tmp= z+m-1-(z+(m-1)*.5)*w/tmp;
            }
        }

        do{
            M=2*M;
            g=1;
            for(int m=M;m>=1;m--){
                if(m%2==0){
                    g= z+m-1+m*.5*w/g;
                }
                else{
                    g= z+m-1-(z+(m-1)*.5)*w/g;

                }
            }
            if(abs(1-tmp/g)<tolerance)escape=true;
            else tmp=g;

        }while(!escape&&M<10000000);

        //XXXXXX maybe put an error message + exit here if M is big
        if(M==10000000){
            cout << "continued fraction for complimentary incomplete gamma failed to converge.";
            cout<< endl << "z = " << z << "  w = " << w << endl;
            exit(1);
        }


        */

        //cout <<setprecision(20) << endl;
        //cout << "g1: " << g << endl;

        int n=0;
        //ttype P1=1.,P2=z,P3,Q1=0.,Q2=1.,Q3;
        //ttype P1=1.,P2=z,Q1=0.,Q2=1.,u=.5*w,t1=z,t2=z*w,t3=0;
        ttype P1=Double(1),P2=z,Q1=Double(0),Q2=Double(1),t1=z,t2=z*w,t3=Double(0);

        do{
            //forward recursion.
            //n++;
            //P3=(z+n)*P2-(z+(n-1)*.5)*w*P1;
            //Q3=(z+n)*Q2-(z+(n-1)*.5)*w*Q1;
            //P1=P2;P2=P3;
            //Q1=Q2;Q2=Q3;

            //n++;
            //P3=(z+n)*P2+n*u*P1;
            //Q3=(z+n)*Q2+n*u*Q1;

            //P1=P2;P2=P3;
            //Q1=Q2;Q2=Q3;


            //faster implementation of the forward recursion. Each block below is
            //equivalent to the above.
            t1+=Double(1);
            P1=t1*P2-t2*P1; Q1=t1*Q2-t2*Q1;
            t3+=w; t1+=Double(1);t2+=w;
            P2=t1*P1+t3*P2; Q2=t1*Q1+t3*Q2;

            t1+=Double(1);
            P1=t1*P2-t2*P1; Q1=t1*Q2-t2*Q1;
            t3+=w; t1+=Double(1);t2+=w;
            P2=t1*P1+t3*P2; Q2=t1*Q1+t3*Q2;

            t1+=Double(1);
            P1=t1*P2-t2*P1; Q1=t1*Q2-t2*Q1;
            t3+=w; t1+=Double(1);t2+=w;
            P2=t1*P1+t3*P2; Q2=t1*Q1+t3*Q2;

            t1+=Double(1);
            P1=t1*P2-t2*P1; Q1=t1*Q2-t2*Q1;
            t3+=w; t1+=Double(1);t2+=w;
            P2=t1*P1+t3*P2; Q2=t1*Q1+t3*Q2;

            n++;

//cout << P2/Q2 << " " << norm(Q2*P1-P2*Q1) / norm(Q2*P1*tolerance) <<endl;

            //to prevent overlow
            if(real(P2)>1.e50||real(P2)<-1.e50||imag(P2)>1.e50||imag(P2)<-1.e50){
                P1*=1.e-50;
                P2*=1.e-50;
                Q1*=1.e-50;
                Q2*=1.e-50;

            }

        //} while(n<3||(my_norm(Q2*P1-P2*Q1)>my_norm(Q2*P1*tolerance)&&n<1000000));
        } while(my_norm(1-P2*Q1/(Q2*P1))>tolerance_sqrd&&n<100000);

        g=P2/Q2;

        //cout<< "using cfrac for comp inc " << t << " " << n << endl;

        if(n>99999){
             cout << "Mofu. Continued fraction for g(z,w) failed to converge. z = "
             << z << "  w = " << w << endl;
             exit(1);
        }



        g=lcalc_exp(-w)/g;

    }

    //if(my_verbose>3) cout << "exit comp_inc_GAMMA("<<z<<","<<w<<")"<< endl;
    return g;

}

template <class ttype>
Complex gamma_sum(Complex s, int what_type, ttype *coeff, int N, Double g, Complex l, Double Q, Long Period, Complex delta=Double(1), const char *method="temme")
{
    Complex SUM=0;

    Complex z,w;
    Complex G;
    Complex r;
    Complex u;


    int n=1;
    int n2=1;

    Double x,MAX=0;
    bool escape=false;
    bool is_z_real=false;
    bool is_w_real=false;


    z=g*s+l;

    if(my_norm(imag(z))<tolerance_sqrd) is_z_real=true;
    if(my_norm(imag(delta))<tolerance_sqrd) is_w_real=true;

    w=delta/Q;
    if(g<.6) w=w*w;     //i.e. if g=1/2


    //y=abs(z)+abs(real(z))+1;


    if(what_type==-1)   //i.e. if the Riemann zeta function
    do{
        w=Pi*n*n*delta*delta;
        G=inc_GAMMA(z,w,method);
        SUM+=G;
        n++;

        x=my_norm(SUM);
        if(x>MAX) MAX=x;

        if(real(w-z)>10) //10 is kind of arbitrary . G(z,w) will decay like e(-Re(w)), so 
                         //we'll escape once Re(w) is around log(10)*DIGITS. So, 10 is okay
                         //as a place to start checking
        {
            if (my_norm(G)<MAX*tolerance_sqrd) escape=true;
        }
    }while(!escape);

    //XXXX checking abs(G)>tolerance is not so smart... since two of the ways 
    //XXXX I compute G, series and nielsen, only gives us 10E-15, say, regardless
    //XXXX how small it really is (since we subtract g(z,w) from GAMMA(z)w^(-z)
    //XXXX and these are nearly equal when w is big enough.
    //XXXX Best to compare estimate for largest of terms yet to be added
    //XXXX and see if this is <MAX*tolerance, where MAX is the maximum
    //XXXX attained by the partial sums.

    else
    do
    {


        w=n*delta/Q;
        //i.e. if g=1/2
        if(g<.6) {
            w=w*w;
        }


        if(l==0){
            u=Double(1);
        }
        else{
            u=lcalc_exp(LOG(n)*l/g);  //XXX this can be stored if it is called repeatedly
        }


        if(coeff[n2]!=0)
        {
            //if both are real, we should send as Doubles. That way we avoid complex arithmetic
            //which is more time consuming
            if(is_z_real&&is_w_real){
                G=inc_GAMMA(Double(real(z)),Double(real(w)),method);
                if (my_verbose>5) cout << "#                        GAMMA SUM with doubles, n^(l/g) b(n) G(" << Double(real(z)) << ", " << Double(real(w)) << ") = " 
                << G*u*coeff[n2] << " SUM = " << SUM << endl;
                //cout<<"both z and w are real\n";
            }
            else{
                G=inc_GAMMA(z,w,method);
                if (my_verbose>5) cout << "#                        GAMMA SUM, G = " << G << endl;
                //cout<<"none are real\n";
            }
            SUM+=G*u*coeff[n2];
        }


        n++; n2++;

        x=my_norm(SUM);
        if(x>MAX) MAX=x;

        if(real(w-z)>10) //10 is kind of arbitrary . G(z,w) will decay like e(-Re(w)), so 
                         //we'll escape once Re(w) is around log(10)*DIGITS. So, 10 is okay
                         //as a place to start checking
        {
            //y3=4*exp(-real(w))/(y2+1);
            //if (my_norm(u)*y3*y3*n*n<MAX*tolerance_sqrd) escape=true;
            if (my_norm(u*G)*n*n<MAX*tolerance_sqrd) escape=true;
        }

        //if(n2>Period&&Period>1) n2=(n2-Period);
        if(n2>Period&&Period>1) n2-=Period;

    }while(n2<=N&&!escape);


    //XXXXnote, we copy the tolerance feature for zeta, but make sure
    //XXXXto also include b(n)<sigma_0(n)< n and n^(Re l/g) factors.
    //XXXX(that is why we have the n*u in the while).
    if(n2>N&&what_type!=-1)
    {

        if(print_warning){
            print_warning=false;
            cout << "WARNING from gamma sum- we don't have enough Dirichlet coefficients." << endl;
            cout << "Will use the maximum possible, though the output ";
            cout << "will not necessarily be accurate." << endl;
        }
        //exit(1);

    }


    max_n = n;
    if(my_verbose>5) cout << "#                        s = " << s << "gamma_sum called, number terms used: " << n << endl;

    return SUM;

}

#endif
