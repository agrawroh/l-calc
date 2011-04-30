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


//This file contains functions for computing N(T).
//N(T) is the number of zeros in the critical strip with |t|<T. notice that this counts
//zeros both below, on, and above the real axis. This is the approach that I decided was
//simplest and most general, rather than distinguish the two cases L is self dual, or not.
//The distinction between self dual or not is taken into account in the Lfind_zeros.h file.

//Nmain(T) computes N(T) without S(T), Nmain(T)=N(T)-S(T)
//(N(T) is the number of zeros in the critical strip with |t|<T)
template <class ttype>
Double L_function <ttype>::
Nmain(Double T)
{
    Double x;
    int i;

    x=T*log(Q)*2/Pi; //contribution from the conductor

    //contribution from the gamma factors
    for(i=1;i<=a;i++){
        x=x+imag(log_GAMMA( (.5+I*T)*gamma[i]+lambda[i] ) - log_GAMMA( (.5-I*T)*gamma[i]+lambda[i] )) /Pi;
        //Note: the imag() takes care of the contribution from the contragredient L-function that
        //arises on applying the functional equation,
        // along with the 1/Pi instead 1/(2Pi).
        // the log_GAMMA(*+IT) - log_GAMMA(*-IT) takes care of a contribution from the one Gamma factor,
        //at T and -T.
    }


    // contribution from the pole at s=1 if there is one
    // (actually considers poles on the line Re_s=1)
    for(i=1;i<=number_of_poles;i++)
      if(real(pole[i])>-10E-7&&real(pole[i])<1+10E-7) x=x+1;

    return x;

}


//Computes S(T) by estimating via a contour integral. Not rigorous, but practical.
//assumes that T is not a zero of L(sigma+iT) or L(sigma-iT), sigma >= 1/2.
//one should avoid calling very close to a zero. As a consequence, computing
//the Nth zero, in the situation that the neighbouring zeros
//are nearby, becomes challenging. Can improve by doing a Turing test for nth zero
//but, for now, below is okay.

template <class ttype>
Double L_function <ttype>::
S(Double T)
{
    Double delta_theta=0; //keeps track of the change in arg L(s)

    Complex s=.5+I*T;

    Double theta;
    Complex z=value(s);
    Double theta_0=imag(log(z));
    Double theta_1=theta_0;

    Double x;
    long long n=0;

    Double deg=0; //degree of the L-function
    for(int j=1;j<=a;j++) deg+=gamma[j];
    deg*=2;


    //sigma is used to determine how far to the right to go to be sure that L(s) is done wrapping
    //around 0. We go far enough so that L(s) stays close to 1, i.e. in same halfplane, so that
    //the contribution of the vertical integral is easy.
    Double sigma; // comparing the Euler product to the Euler product of zeta, and using
    //the ramanujan bound, we have, for sigma > 1:
    //a(n) <= d_{deg}(n), i.e. nth dirichlet coeff bounded by {deg}-divisor numbers, coeffs of zeta(s)^{deg}.
    //Thus L(s) is close to 1 (in same halfplane) when d_{deg}(2)/2^sigma+d_{deg}(3)/3^sigma+...<1,
    //i.e. when zeta(sigma)^{deg} < 2. Now zeta(sigma)< 1 + 1/2^sigma+...+ 1/M^sigma+M^{1-sigma}/(sigma-1)) (compare the terms M+1,...
    //to the integral of t^{-sigma} from M to infinity.

    int M=100;
    sigma=1.;
    do{
        sigma+=.1;
        x=1;
        for(int m=2;m<=M;m++) x+=pow(Double(m),-sigma);
        x+=pow(Double(M),1-sigma)/(sigma-1);
        x=pow(x,deg);
    }while(x>2);
    cout<<setprecision(DIGITS);
    if(my_verbose>3) cout << "S(T): rectangle: 1/2-iT, sigma-iT,sigma+iT,1/2+iT, T= " << T << " ,sigma= " << sigma << endl;



    bool is_complex=false;

    if(what_type_L!=-1) //if not the zeta function
    do{
        n++;
        if(abs(imag(Complex(dirichlet_coefficient[n])))>.00001) is_complex=true;
    }while(n<number_of_dirichlet_coefficients&&!is_complex);


    cout<<setprecision(DIGITS);
    if(my_verbose>4)
        cout << "S(T): theta_0=" << theta_0 << endl;

    do{
        s+=.1*abs(z)/abs(value(s,1)); //increment by one tenth over the log derivative, thus L incr should be about one tenth
        if(real(s)>sigma) s=sigma+I*T;
        //s+=.1;
        z=value(s);
        theta=imag(log(z));
        cout<<setprecision(DIGITS);
        x=theta-theta_1;
        delta_theta-=x; //keep track of the change in the arg
        if(my_verbose>4)
            cout << "S(T): theta  =" << theta << " delta_theta= " << delta_theta << " delta delta=" << x << " s= " << s << endl;

        //if we've crossed branch cut of the logarithm add or subtract 2pi.
        //4 is an approximation to 2pi which should be sufficiently small to take into account increment size.
        if(-x<-4){
            delta_theta+=twoPi;
            if(my_verbose>4)
                cout << "S(T): theta_jump  =" << delta_theta << " s= " << s << endl;
        }
        if(-x>4){
            delta_theta-=twoPi;
            if(my_verbose>4)
                cout << "S(T): theta_jump  =" << delta_theta << " s= " << s << endl;
        }
        theta_1=theta;

    }while(real(s)<sigma);

    //if is self dual the horizontal contribution below the axis is the same as above the axis.
    if(!is_complex){
         //the +theta is to take into account the change in arg along the vertical
         if (my_verbose>3) cout << "SSS not complex S(" << T << ")=" << 2*(delta_theta+theta)/Pi << endl;
         return (2*(delta_theta+theta)/Pi);
    }

    cout<<setprecision(DIGITS);
    if(my_verbose>3)
        cout << "S(T): L is complex. Evaluate delta arg L below real axis..." << endl;

    delta_theta+=theta - imag(log(value(conj(s)))); //the vertical contribution
    s=.5-I*T;
    z=value(s);
    theta_0=imag(log(z));
    theta_1=theta_0;
    do{
        s+=.1*abs(z)/abs(value(s,1)); //increment by one tenth over the log derivative, thus L incr should be about one tenth
        //s+=.1;
        if(real(s)>sigma) s=sigma-I*T;
        z=value(s);
        theta=imag(log(z));
        cout<<setprecision(DIGITS);
        x=theta-theta_1;
        delta_theta+=x; //keep track of the change in the arg
        if(my_verbose>4)
            cout << "S(T): theta  =" << theta << " delta_theta= " << delta_theta << " s= " << s << endl;

        //if we've crossed branch cut of the logarithm add or subtract 2pi.
        //4 is an approximation to 2pi which should be sufficiently small to take into account increment size.
        if(x<-4){
            delta_theta+=twoPi;
            if(my_verbose>4)
                cout << "S(T): theta_jump  =" << delta_theta << " s= " << s << endl;
        }
        if(x>4){
            delta_theta-=twoPi;
            if(my_verbose>4)
                cout << "S(T): theta_jump  =" << delta_theta << " s= " << s << endl;
        }
        theta_1=theta;

    }while(real(s)<sigma);


    if (my_verbose>3) cout << "SSS is complex S(" << T << ")=" << delta_theta/Pi << endl;
    return (delta_theta/Pi);

}


//Computes N(T)
//N(T) is the number of zeros in the critical strip with |t|<T. notice that this counts
//zeros both below, on, and above the real axis.
//Note: the arg principle is used with Lambda(s) which only has non-trivial zeros
//so we can take as large a rectangle as we want to the left and right. Davenport is a
//bit weird in his treatment of L(s,chi) zero at 0 or -1. Not sure why he thinks that's
//relevant.
template <class ttype>
Double L_function <ttype>::
N(Double T)
{
    if(T<tolerance) return 0; //N should not be called with T=0. But return 0 if it is.
    Double x1=Nmain(T),x2=S(T);
    cout <<setprecision(DIGITS);
    if(my_verbose>2) cout << "N(T): T = " << T << " Nmain(T)= " << x1 << ", S(T) = " << x2 << endl;
    return (x1+x2);

}


//find the nth gram point t such that Nmain(t)=n
//similar to function nth_gram, but without the factor of 1/2. This allows us to capture n odd,
//(for example, in the case of an odd rank elliptic curve).
template <class ttype>
Double L_function <ttype>::
Nmain_inverse(Long n){

     if(n==0) return 0.;

     Double t=1000.,r,x,b=(n*1.)*(n*1.)*tolerance_sqrd;
     do{
         x=(double(n)-Nmain(t));
         r=(Nmain(t-.05)-Nmain(t+.05))*10; //this is essentially Newton's method, but with a crude approximation
                                    //for the derivative. Works reasonably well. Reason for approx: in case
                                    //n, hence t, is large, we won't have much precision left over for the derivative.
         t-=x/r;
         if(my_verbose>3)
             cout << "Nmain_inverse("<< n << "): Nmain(" << t << " )=" << Nmain(t) << "  , difference:" << x << " vs " << b << endl;
     }while(x*x>b);

     return t;
}

//rough approximation for the local density of zeros
//the interval is of length .1, but N counts zeros above and below,
//so we divide by .2 rather than by .1. Notice this is a rough approx,
//so we don't care about .2 and .05 being interpreted by the compiler 
//as doubles rather than Doubles

template <class ttype>
Double L_function <ttype>:: density_zeros(Double T){
    return ((Nmain(T+.05)-Nmain(T-.05))/.2);
}
