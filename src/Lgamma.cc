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


#include "Lgamma.h"


//temme's uniform asymptotics coded by Kyle Wichert
//based on Temmes paper: The asymptotic expansion of the incomplete
//gamma functions, SIAM J. Math Anal, vol 10, number 4, July 1978

Complex Q(Complex z, Complex w)
{

  if(my_verbose>5) cout << "entering temme Q z= " << z << " " << " w= " << w << endl;

  Complex b=1.,sum=.5,eta,c[501],X,r=0.,R,v;
  int k=0,n=2;



  eta=(w-z);
  eta=eta/z;
 //cout << eta<< endl;
  do{
    //b=b*(-1)*( ((Double) n)/((Double)(n+1)))*eta;
    b=-b*eta;
    v=b/(n+1);
    sum=sum+v;
    n++;
  }while(abs(v)>tolerance);
  //cout << "my eta = " << exp(log(2*(eta-log(1+eta)))/2) << endl;
  eta=eta*sqrt(2*sum);

  //cout << "eta ratio = " <<eta/sqrt(2*(eta-log(1+eta)) << endl;

//must do one run for k==0, then the general stuff below begins to work
  n=0;
  c[0]=-1.;
  c[0]=c[0]/3; //if we do c[0]=-1./3 we get problems with long double. And if we do c[0]=-1/3.L we might get problems with multiprecision
  Complex tmp=1.;
  Complex tmp2;
  while(abs(((n+2)*temme_a[n+2]*tmp)/c[0]) >= tolerance)
  {
    n++;
    tmp=tmp*eta;
    c[0]+=(n+2)*temme_a[n+2]*tmp;
  }
  r=c[0];

  do
  {
    k++;
    // (3.9)
    n=0;
    c[k]=0;
    tmp=1.;
    do
    {
      X= n==0? -temme_g[k]/3:(Double)0; Double df=1;
      for (int i= n==0? 1:0;i<=k;i++){
          if(i>0) df*=(n+2*i);
          X+=temme_g[k-i]*df*(n+2*i+2)*temme_a[n+2*i+2];
          //X+=temme_g[k-i]*dfac(n+2*i)/dfac(n)*(n+2*i+2)*temme_a[n+2*i+2];
      }
      c[k]+=X*tmp;
      n++; tmp=tmp*eta;
    } while(my_norm(X*tmp/c[k]) >= my_norm(tolerance*eta));
    // End (3.9)

    //tmp=1.;
    //for(int i=1;i<=k;i++) tmp=tmp/z;
    //r+=c[k]*tmp;
    tmp2=c[k]*pow(1./z,k);
    r+=tmp2;


} while(my_norm(tmp2/r) >= tolerance_sqrd);

  R=1.0/sqrt(2*Pi*z) * exp(-.5*z*eta*eta) * r;

  //erfc depends on eta*sqrt(z/2) so eta*=sqrt(z/2)
  eta=eta*sqrt(z/2);



  if(my_verbose>5) cout << "exit temme" << endl;
  return .5*erfc2(eta)+R;

}

//XXXXX Kyle's formward recursion in erfc is a bit funny.
//below works okay for 16 digits... but something will
//have to replace this for higher precision.
//can use sqrt(pi) erfc(z) = incomplete gamma(1/2,z^2) =z G(1/2,z^2)
//and is now implemeneted in erfc2 below. So erfc can be deleted.
Complex erfc(Complex z)
{

  if(my_verbose>5) cout << "entering erfc\n";

  Complex x;
  int n;

  if(real(z)<0) return (2-erfc(-z));

  if(my_norm(z) > 64) {

    if(abs(z) > 70) n=10;
    else if(abs(z) > 14) n=20;
    else n=30;
    x=2*z;
    for (int i=n;i>0;i-=2)
      x=2*z+i/x;
    x=2*exp((-1)*z*z)/(sqrt(Pi)*x);
    return x;
  }

  else {
    if(abs(z) < 1.4) n=20;
    else if(abs(z) < 4.2) n=55;
    else n=125;
    x=2*n+1;
    for (int i=2*n-1;i>0;i-=2)
      x=i+(1-2*(((i+1)/2)%2))*(i+1)*z*z/x;
    x=1-2*z*exp(-z*z)/(sqrt(Pi)*x);
    return x;
  }

}


Complex erfc2(Complex z)
{
  if(my_verbose>5) cout << "entering erfc2\n";

  Complex x;

  if(real(z)<0) return (2-erfc2(-z));
  //inc_GAMMA requires two variables to be of same type
  if(my_norm(z) > .5) return (cfrac_GAMMA(.5+0*I,z*z)*z/sqrt(Pi)); 
  else return ((GAMMA(.5,z*z)-comp_inc_GAMMA(.5+0*I,z*z))*z/sqrt(Pi));


}

