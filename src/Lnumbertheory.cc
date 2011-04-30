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


#include "Lnumbertheory.h"

using namespace std;


bool isprime(long long N) 
{


  long long a;
  double maxrand,b;

//  From http://www.utm.edu/research/primes/prove/prove2_3.html
//  If n < 1,373,653 is a both 2 and 3-SPRP, then n is prime. 
//  If n < 25,326,001 is a 2, 3 and 5-SPRP, then n is prime. 
//  If n < 25,000,000,000 is a 2, 3, 5 and 7-SPRP, then either n = 3,215,031,751 or n is prime. (This is actually true for n < 118,670,087,467.) 
//  If n < 2,152,302,898,747 is a 2, 3, 5, 7 and 11-SPRP, then n is prime. 
//  If n < 3,474,749,660,383 is a 2, 3, 5, 7, 11 and 13-SPRP, then n is prime. 
//  If n < 341,550,071,728,321 is a 2, 3, 5, 7, 11, 13 and 17-SPRP, then n is prime. 

  if(N==2 || N==3 || N==5 || N==7 || N==11 || N== 13 || N== 17) return true;
  if(!N%2 || !N%3 || !N%5 || !N%7 || !N%11 || !N%13 || !N%17) return false;

  if(N<1373653) { if(RM(2,N) && RM(3,N)) return true;
                  else return false; }
  if(N<25326001) { if(RM(2,N) && RM(3,N) && RM(5,N)) return true;
                   else return false; }
  if(N<strtoll("25000000000",0,10)) { if(N!=strtoll("3215031751",0,10) && RM(2,N) && RM(3,N) && RM(5,N) && RM(7,N)) return true;
                      else return false; }
  if(N<strtoll("2152302898747",0,10)) { if(RM(2,N) && RM(3,N) && RM(5,N) && RM(7,N) && RM(11,N)) return true;
                        else return false; }
  if(N<strtoll("3474749660383",0,10)) { if(RM(2,N) && RM(3,N) && RM(5,N) && RM(7,N) && RM(11,N) && RM(13,N)) return true;
                        else return false; }
  if(N<strtoll("341550071728321",0,10)) { if(RM(2,N) && RM(3,N) && RM(5,N) && RM(7,N) && RM(11,N) && RM(13,N) && RM(17,N)) return true;
                          else return false; }

  maxrand=RAND_MAX;

  if(N>strtoll("341550071728321",0,10)) {
    srand((unsigned int)time((time_t *)NULL)); 
    for (int i=1;i<=20;i++) {
      b=(rand())/(maxrand+1);
      a=(long long)(N*b+1);
      if(!RM(a,N) && !(a==N)) return false;
    }
  }
  return true;

}


bool RM(long long a, long long N)
{
  long long q,t,b,e;

//factoring N-1 to (2^t)*q
  q=N-1;
  t=0;
  while(!(q&1)) {
    q >>= 1;
    t++;
  }

// Set b=a^q Mod N
  b=1;
  while(q>0)  {
    if (1 & q) b=multmodN(b,a,N);
    a=multmodN(a,a,N);
    q >>= 1;
  }
  e=0;

  if(b==1)  return true;
  while(b!=1 && b!=(N-1) && e<=(t-2)) {
    b=multmodN(b,b,N);
    e++;
  }
  if(b==(N-1)) return true;
  else return false;


}


long long multmodN(long long a, long long b, long long N)
{
  long long s;

  if((a<1073741823) && (b<1073741823)) return (a*b)%N;

  s=0;
  while(a > 0) {
    if(a & 1) s=(s+b)%N;
    b=(2*b)%N;
    a>>=1;
  }
  return s;
}


long long nextprime(long long n)
{
    long long m=n;

    do{
        m++;
    }while(!isprime(m));

    return m;
}

//assumes n is positive
bool issquarefree(long long n)
{

    long long i=1;

    for(i=2;(i*i)<=n&&n%(i*i)!=0;i++);

    if(n%(i*i)==0) return false;
    else return true;

}

// checks if n is a fundamental discriminant
bool isfunddiscriminant(long long n)
{
    long long m=n;
    if(n<0) m=-n;

    if(m%2==1)
    {
        if((n<-2&&m%4==3)||(n>2&&m%4==1) )
            return(issquarefree(m));
        else return false;

    }
    else if(m%8==4)
    { 
        if((n<-2&&(m/4)%4==1)||(n>2&&(m/4)%4==3) )
            return(issquarefree(m/4));
        else return false;

    }

    else if(m%16==8) return(issquarefree(m/8));

    else return false;
}

long long nextfunddiscriminant(long long n)
{
    long long m=n;

    do{
        m++;
    }while(!isfunddiscriminant(m));

    return m;
}


//=================== Kronecker routines =========================

//takes any integers n and m and computes the kronecker symbol (n|m).
//integer airthmetic such as n%4 is much faster than long long arithmetic.
//so one should use my_kronecker(int n,int m) rather than
//my_kronecker(long long n,long long m) whenever possible

//XXX below I think we should replace the %2 %4 and %8 operations 
//with something that looks just at the last 1,2,3 bits. Should be much faster

int my_kronecker(int n,int m)
{

    int i,s,t;
    bool p;


    //if(n%2==0&&m%2==0) return(0);
    if(!(n & 1) && !(m & 1)) return(0);

    if(m==0) return(0);
    s=1;

    //if m is negative
    if(m<0){
        m=-m;
        if(n<0) s = -s;
    }

    //if m is even
    p=false;
    while(!(m & 1)){ //while m is even
        //m=m/2;
        m=m>>1; //divide m by 2 by bit shifting
        p=!p;
    }
    if(p){ //i.e. was m divisable by an odd power of 2
        //t=n%8;
        t=(n&7);
        if(t<0)t=t+8; //is this step needed? yes! t might be negative
        if(t==3||t==5) s=-s;
    }

    //if n is negative
    if(n<0){
        n=-n;
        if((m&3)==3)s=-s;
        //if(m%4==3)s=-s;
    }

    if(m==1) return(s);

    //if(kronecker_table_available && n<kronecker_bound && m<kronecker_bound) return (s*kronecker_table[n][m]);

    n=n%m;

    if(n==0)return(0);

    else{
        p=false;
        while(!(n & 1)){ //if n is even
                n = n>>1;
                p = !p;
        }
        if(p){
                ////t = m%8;
                short int t = m&7;
                if(t==3||t==5) s = -s;
                //if(m&1&&(((m>>1)&1)^((m>>2)&1))) s = -s;
        }
        if(n!=1){
            if((m&3)==3&&(n&3)==3)s = -s; //uses quadratic reciprocity
            i=n;
            n = m%n; //so n -> smaller than n

            if(n==0)return(0);
            m=i; // m -> previous n, so we now have n < m

            //at this point n and m are both positive and m is odd and n < m
            //so we can carry out a simplified routine
            do{
                //cout << m << " " << n << endl;

                //if(kronecker_table_available && m<kronecker_bound) return (s*kronecker_table[n][m]);

                bool p=false;
                while(!(n & 1)){ //if n is even  //can be sped up by just counting the low zeros? nope.
                        n = n>>1; //divide by 2
                        p = !p;
                }
                if(p){
                //int power_2; //power of 2 that divides n
                //for (power_2 = 0 ; power_2<32 ;power_2++) //this assumes that ints are 32 bits. Will that stay the same forever.
                    //if (n&(1<<power_2)) break;
                //n=n>>power_2;
                //if(power_2&1){
                        short int t = m&7;
                        if(t==3||t==5) s = -s;
                        //if(m&1&&(((m>>1)&1)^((m>>2)&1))) s = -s;
                }
                if(n!=1){
                   if((m&3)==3&&(n&3)==3)s = -s; //uses quadratic reciprocity
                   //if(((m&3)&(n&3))==3)s = -s; //uses quadratic reciprocity
                   int i=n;
                   n = m%n;
                   m=i;
                   if(n==0)return(0);
                }
                else m=n;
            }while(m!=1);

        }
    }
    return(s);

}


//takes any integers n and m and computes the kronecker symbol (n|m)
int my_kronecker(long long n,long long m)
{

    long long i;
    int p,s,t;

    if(n%2==0&&m%2==0) return(0);

    if(m==0) return(0);
    s=1;

    //if m is negative
    if(m<0){
        m=-m;
        if(n<0) s = -s;
    }

    //if m is even
    p=0;
    while(m%2==0){
        m=m/2;
        p=1-p;
    }
    if(p==1){ //i.e. m divisable by an odd power of 2
        t=n%8;
        if(t<0)t=t+8; //is this step needed? yes! t might be negative
        if(t==3||t==5) s=-s;
    }

    //if n is negative
    if(n<0){
        n=-n;
        if(m%4==3)s=-s;
    }

    if(m==1) return(s);


    n=n%m;

    if(n==0)return(0);

    else{
        p=0;
        while(n%2 == 0){ //if n is even
                n = n/2;
                p = 1-p;
        }
        if(p==1){
                t = m%8;
                if(t==3) s = -s;
                if(t==5) s = -s;
        }
        if(n!=1){
           if(m%4==3){
                   if(n%4==3)s = -s; //uses quadratic reciprocity
           }
           i=n;
           n = m%n; //so n -> smaller than n
           if(n==0)return(0);
           m=i; // m -> previous n, so we now have n < m

           //at this point n and m are both positive and m is odd and n < m
           //so we can call a simplified routine
           if(m < INT_MAX) //faster to use ints if longs aren't needed
           {
               //s=s*simplified_jacobi(int(n),int(m));
               int nn,mm,ii;
               nn=Int(n);
               mm=Int(m);
               do{
                   p=0;
                   while(nn%2 == 0){ //if n is even
                           nn = nn/2;
                           p = 1-p;
                   }
                   if(p==1){
                           t = mm%8;
                           if(t==3) s = -s;
                           if(t==5) s = -s;
                   }
                   if(nn!=1){
                      if(mm%4==3){
                              if(nn%4==3)s = -s; //uses quadratic reciprocity
                      }
                      ii=nn;
                      nn = mm%nn;
                      if(nn==0)return(0);
                      mm=ii;
                   }
                   else mm=nn;
               }while(mm!=1);
           }
           else s=s*simplified_jacobi(n,m);

        }
    }
    return(s);

}

//assumes n and m are positive and m is odd and n < m
int simplified_jacobi(int n,int m)
{

    int p,s,t,i;


    s=1;
    do{
        p=0;
        while(n%2 == 0){ //if n is even
                n = n/2;
                p = 1-p;
        }
        if(p==1){
                t = m%8;
                if(t==3) s = -s;
                if(t==5) s = -s;
        }
        if(n!=1){
           if(m%4==3){
                   if(n%4==3)s = -s; //uses quadratic reciprocity
           }
           i=n;
           n = m%n;
           if(n==0)return(0);
           m=i;
        }
        else m=n;
    }while(m!=1);
    return(s);
}

//assumes n and m are positive and m is odd and n < m
int simplified_jacobi(long long n,long long m)
{

    long long i;
    int p,s,t;

    s=1;
    p=0;
    while(n%2 == 0){ //if n is even
            n = n/2;
            p = 1-p;
    }
    if(p==1){
            t = m%8;
            if(t==3) s = -s;
            if(t==5) s = -s;
    }
    if(n!=1){
       if(m%4==3){
               if(n%4==3)s = -s; //uses quadratic reciprocity
       }
       i=n;
       n = m%n;
       if(n==0)return(0);
       m=i;
       //at this point n and m are both positive and m is odd
       //so we can call a simplified routine
       if(m < INT_MAX)
       s=s*simplified_jacobi(Int(n),Int(m));
       else
       s=s*simplified_jacobi(n,m);
    }
    return(s);

}

//===================end Kronecker routines ======================

//d a fund. discriminant
//returns L(1,chi_d) where chi_d(n) is (d|n), Kronecker's symbol.
//Can be used to check L-function computation of L(1,chi),
//and as alternative method (though less efficient) for computing
//the class number of Q(srt(d)).

Double L_1_chi(int d)
{

    int j=0;
    Double L=0;

    if(isfunddiscriminant(d))
    {
        if(d<0){
            for(j=1;j<=-d;j++) L=L+j*my_kronecker(d,j);
            L=-L*Pi*exp(-Double(3)/2*log(-Double(1)*d));
        }
        else{
            for(j=1;j<=d-1;j++){
                L=L+my_kronecker(d,j)*log(sin(j*Pi/d));
            }
            L=-L*exp(-log(Double(1)*d)/2);
        }

        return L;

    }
    else return 0;


}

//returns the class number of the quadratic field q(sqrt(D))
//as described by (15) and (16) of Davenport.

int class_number(int d)
{
    int h,w=2;

    Double t,x,u;
    Double epsilon;

    if(isfunddiscriminant(d)){
        if(d<0){
            if(d==-4)w=4;
            if(d==-3)w=6;
            h=Int(w*exp(log(-Double(1)*d)/2)*L_1_chi(d)/(2*Pi)+.000001);
        }
        else{
            u=0.;
            do{
                //xxx here I should replace with
                //a continued fraction algorithm
                u=u+1;
                x=sqrt(4+Double(d)*u*u)+.000001;
                t=floor(x);
            }while(x-t>.00001);
            cout << "t = " << t << "  u = " << u << endl;
            epsilon=(t+u*sqrt(Double(d)))/2;
            h=Int(exp(log(Double(1)*d)/2)*L_1_chi(d)/log(epsilon)+.000001);
        }
        return h;
    }
    else return 0;

}


//returns co[n]=tau(n)/n^(11/2), n=0...N.
//these are the (normalized) coefficients of the weight 12 modular form
void ramanujan_tau(Double *c0, int N_terms)
{

    Double *c;
    Double *c1;
    Double x=(-1+sqrt(1+8*Double(N_terms)))/2+.0000001;
    int i=Int(x);
    int k,m,n;

    c= new Double[i+2];
    c1= new Double[N_terms+2];

    for(i=0;i<=N_terms;i++)c0[i]=0.;
    c0[0]=1;

    for(n=0;n<=x;n++){c[n]=(2*n+1);if(n%2==1)c[n]=-c[n];}

    for(k=1;k<=8;k++){
       for(m=0;m<=N_terms;m++){
          c1[m]=0.;
          x=(-1+sqrt(1+8*Double(m)))/2+.00000001;
          for(n=0;n<=x;n++){
               c1[m]+=c[n]*c0[m-(n*(n+1))/2];
          }
       }
       for(m=0;m<=N_terms;m++)c0[m]=c1[m];
    }

    c0[0]=0.;
    c0[1]=1;
    for(n=0;n<=N_terms-1;n++){
      x=1./(Double)(n+1);
      x=x*x*x;
      c0[n+1]=c1[n]*x*x*sqrt((Double)(n+1)); //i.e. divided by (n+1)^5.5
    }

    delete [] c;
    delete [] c1;


}

long long gcd(long long a,long long b)
{
    long long r;

    if(a>b){r=a;a=b;b=r;}

    if(a==0)return b;

    do{
        r=(b%a);b=a;a=r;
    }while(r>0);
    return b;
}


//computes a^k mod q using binary expansion of k
long long power_mod_q(long long a, long long k,long long q)
{
    long long r,m;

    if(k==0) return 1;

    //cout << a <<"^"<<k << " mod " << q ;
    m=a; r=1;
    do{
        if(k%2==0) k=k/2;
        else{
            k=(k-1)/2;
            r=(r*m)%q;
        }
        m=(m*m)%q;
    }while(k>0);

    //cout <<  " = " << r << endl;
    return(r);
}


//returns a primitive root for the cyclic
//group of (Z/p^alpha Z). Assumes p is an odd prime.
//first find a generator r mod p, then checks if r or r+p
//is a generator mod p^alpha.
//algorithm is from Henri Cohen's book
int prim_root(long long p, int alpha)
{
    int a,j,k;
    long long **factors;
    bool is_primitive; //is 'a' a primitive root or not

    factors = new long long *[30];
    for(k=0;k<=29;k++) factors[k]= new long long[3];
    for(k=0;k<=29;k++)
    for(j=0;j<=2;j++)factors[k][j]=0;

    factor(p-1,factors);

    //find a generator mod p
    a=1;
    do{
        is_primitive=true;
        a++;
        if(factors[0][1]>0){
            if(my_kronecker((long long)a,p)!=-1) is_primitive=false;
            //cout << 2 << " there " << my_kronecker(a,p)<< " " << is_primitive << endl;
        }
        k=1;
        if(is_primitive&&factors[0][2]>0)
        do{
            if(power_mod_q(a,(p-1)/factors[k][0],p)==1)
                is_primitive=false;
            //cout << (p-1)/factors[k][0]  << " here " << is_primitive << endl;
            k++;
        }while(k<=factors[0][2]&&is_primitive);

    }while(!is_primitive&&a<p-1);

    for(k=0;k<=29;k++) delete[](factors[k]);
    delete[] (factors);

    if(alpha==1) return a;

    if(power_mod_q(a,p-1,p*p)!=1) return a;
    else return(a+p);
}


//factor returns a vector containing:
//factors[k][0]=p, prime dividing q
//factors[k][1]=alpha, power of that prime dividing q
//factors[k][2]=generator for the cyclic group mod p^alpha
//factors[0][0]=2, factors[0][1]=power of 2 dividing q
//factors[0][2]=k, number of odd p's diving q

void factor(long long q, long long **factors)
{
     long long m;
     int k=0;
     long long Q=q;

    factors[0][0]=2;
    factors[0][1]=0;

    while(Q%2==0){
        factors[0][1]++;
        Q=Q/2;
    }

    for(m=3;m*m<=q;m++)
    {
         while(Q%m==0&&isprime(m)){
             k++;
             factors[k][0]=m;
             factors[k][1]=0;
             while(Q%m==0){
                 factors[k][1]++;
                 Q=Q/m;
             }
         }


    }
    if(Q*Q>q){
        k++;
        factors[k][0]=Q;
        factors[k][1]=1;
    }

    factors[0][2]=k;
    for(m=1;m<=k;m++)factors[m][2]=prim_root(factors[m][0],factors[m][1]);
}

//characters are implemented as described in Davenport, especially
//(1) of Chapter 4, pg 29.
int characters()
{

    long long q;   //conductor
    long long q_1; //conductor of the inducing character chi_1
    int n; //as in chi[n]
    int K; //number of distinct odd primes dividing q

    int j,k; //loop variables
    int loop1,loop2; //looping bounds
    int a,b; //various ints

    Complex u,w,z;  //used to generate the character

    Complex *chi;   //the character
    Complex *chi_1; //the primitive character inducing chi
    Complex *sum;   //to check that sum over chi is zero
    Complex SUM;   //to check that sum over n is zero
    Double total=0;  //sum of abs(SUM) over all chi
    long long **factors;  //stores the factors of q

    int *m, *v;    //as in (1) chapter 4 of Davenport
    int m_prime, v_prime; 
    int *phi; //stores the values of phi(p^alpha)
    int *power_p;// stores p_j^alpha_j

    int PHI; //is set to phi(q/2^alpha_0)
    int M; //for looping 0..PHI-1
    int dial; //for creating an "analog dial" looping through the m[j]'s

    int r; //generator for the cyclic group of (Z/p^alpha Z)

    bool is_primitive; //for storing whether chi is primitive or not


    factors = new long long *[30];
    for(k=0;k<=29;k++) factors[k]= new long long[3];
    for(k=0;k<=29;k++)
    for(j=0;j<=2;j++)factors[k][j]=0;

    cin >> q;
    factor(q,factors);

    K=factors[0][2];

    for(k=0;k<=K;k++)
        cout << q <<" " << factors[k][0] << " " << factors[k][1] <<" " << factors[k][2] << endl;

    m = new int [K+1];
    v = new int [K+1];
    phi = new int [K+1];
    power_p = new int [K+1];
    chi = new Complex[q+1];
    chi_1 = new Complex[q+1];
    sum = new Complex[q+1];
    for(n=0;n<=q;n++)sum[n]=0;

    PHI=1;
    //shortcut for evaluating phi(p_j^alpha_j)
    for(j=1;j<=K;j++){
        phi[j]=(factors[j][0]-1); //(p-1)
        power_p[j]=power_mod_q(factors[j][0],factors[j][1],q+1);//p^alpha

        if(factors[j][1]>1)
        phi[j]=phi[j]*power_p[j]/factors[j][0];

        PHI=PHI*phi[j];
    }


    if(factors[0][1]>2) loop1=power_mod_q(2,factors[0][1]-2,q);
    //shortcut for evaluating 2^(alpha_0-2) (<q so mod q doesn't change result)
    else loop1=1;

    if(factors[0][1]<=1) loop2=1;
    else loop2=2;

   for(j=1;j<=K;j++)
        m[j]=0;

    for(m_prime=0;m_prime<loop1;m_prime++)
    for(m[0]=0;m[0]<loop2;m[0]++)
    for(M=0;M<PHI;M++)
    {

        is_primitive=true;


        for(j=1;j<=K;j++)
            if(gcd(m[j],factors[j][0])>1) is_primitive=false;

        // m[j]'s,m_prime are set. now compute chi(n)


        for(n=1;n<=q-1;n++){
            chi[n]=1;
            chi_1[n]=1;
        }
        chi[0]=0;
        chi_1[0]=0;
        chi[q]=0;
        chi_1[q]=0;

        q_1=q;

        for(k=1;k<=K;k++)
        {

            b=power_mod_q(factors[k][0],factors[k][1],q+1);//=p^alpha
            //q+1 is needed so that p<q+1 in the event q is prime

            //if m[k]==0 drop that character from chi_1
            q_1=q_1/gcd(m[k],b);

            n=1;
            r=factors[k][2]; //generator mod p^alpha
            z=exp(2.*m[k]*I*Pi/(1.*phi[k])); u=1.;
            for(v[k]=1;v[k]<=phi[k];v[k]++)
            {
                u=u*z;
                n=(n*r)%b;
                a=n;
                do{
                    chi[a]=chi[a]*u;
                    if(m[k]!=0)chi_1[a]=chi_1[a]*u;
                    a=a+b;
                }while(a<q);
            }//for v[k]

        }//for k

        //do powers of 2
        if(factors[0][1]==1){
            is_primitive=false;
            q_1=q_1/2;
        }
        if(factors[0][1]==2&&m[0]==1){
            for(n=3;n<=q;n=n+4){
                chi[n]=-chi[n];
                chi_1[n]=-chi_1[n];
            }
        }
        if(factors[0][1]==2&&m[0]==0){
            is_primitive=false;
            q_1=q_1/4;
        }
        if(factors[0][1]>2){

            if(m_prime%2==0){
                is_primitive=false;
                if(m[0]==0&&m_prime==0)q_1=q_1/(4*loop1); //strip powers of 2
                else q_1=q_1/gcd(m_prime,loop1);
            }

            n=1;
            if(m[0]==0) w=1.; else w=-1.;
            z=exp(2.*m_prime*I*Pi/(1.*loop1));
            u=1.;
            b=4*loop1;//=2^alpha_0
            for(v[0]=1;v[0]<=2;v[0]++){
                n=(-n)%b;
                if(n<0)n=n+b;
                u=u*w;
                for(v_prime=1;v_prime<=loop1;v_prime++){
                    n=(5*n)%b;
                    u=u*z;
                    a=n;
                    do{
                        chi[a]=chi[a]*u;
                        if(m[0]!=0||m_prime!=0)chi_1[a]=chi_1[a]*u;
                        a=a+b;
                    }while(a<q);

                }
            }
        }

        for(n=0;n<=q;n++){
            if(gcd(n,q)>1)chi[n]=0;
            if(gcd(n,q_1)>1)chi_1[n]=0;
        }

        cout << "inducing q_1 = " << q_1 << endl;
        cout << "M = " << M << " ;" ;
        cout << m_prime << " ";
        for(n=0;n<=K;n++) cout << m[n] << " ";
        if(is_primitive) cout << "is primitive" << endl;
        else cout << "NOT primitive" << endl;

        //end compute chi(n)
        // check that sum over characters is 0;
        SUM=0;
        for(n=0;n<=q;n++){
            sum[n]=sum[n]+chi[n];
            SUM=SUM+chi[n];
            cout << q << " " << n << " " <<  chi[n] << " " << chi_1[n] << endl;
        }

        cout <<"SUM over n = " << abs(SUM) <<endl;
        total=total+abs(SUM);

        dial=0;
        if(PHI>1)
        do{
            dial++;
            m[dial]=(m[dial]+1)%phi[dial];
        }while(m[dial]==0&&dial<K);

        cout<< "--------------------------------------------\n";
    }//for m_prime. Done looping through all characters.

    for(n=0;n<=q;n++){
        cout << q << " sum[" << n << "] = " <<  sum[n] << endl;
    }

    //cout << "total = " << total << endl;

    for(k=0;k<=29;k++) delete(factors[k]);
    delete (factors);

    delete [] chi;
    delete [] chi_1;
    delete [] phi;
    delete [] power_p;
    delete [] m;
    delete [] v;
    delete [] sum;

    return 0;
}


/*
//returns an array of Ramanujan's tau_n/n^(11/2) for n <= N_terms
void tau(Double *coeff,int N_terms){

    Double x;
    Double *c,*c0;
    int i,n,k,m;


    c0=new Double[N_terms+1];
    c=new Double[N_terms+1];


    for(i=0;i<=N_terms;i++)c0[i]=0;
    c0[0]=1;

    x=rint(-1+sqrt(1+8.*N_terms))/2;
    for(n=0;n<=x;n++){
        c[n]=(2*n+1);
        if(n%2==1)c[n]=-c[n];
    }

    for(k=1;k<=8;k++){
       for(m=0;m<=N_terms;m++){
          coeff[m]=0;
          x=rint(-1+sqrt(1+8.*m))/2;
          for(n=0;n<=x;n++){
               coeff[m]=coeff[m]+c[n]*c0[m-(n*(n+1))/2];
          }
       }
       if(k<8) for(m=0;m<=N_terms;m++) c0[m]=coeff[m];
    }

    for(m=N_terms;m>=1;m--) coeff[m] = coeff[m-1]*pow(1.*m,-5.5);

    //cout <<setprecision(16);
    //for(m=1;m<=N_terms;m++) cout << m << " " << coeff[m] << endl;

    delete [] c0;
    delete [] c;

}
*/
