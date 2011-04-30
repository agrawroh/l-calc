// This file partially due to Kevin McGown, with modifications by Michael Rubinstein

#ifndef Lexplicit_formula_H
#define Lexplicit_formula_H


//#include "Lglobals.h"
//#include "Lmisc.h"


/***************************************************************************************************/

/*
//// the phi in the explicit formula
//inline Complex xxx_phi(Double alpha, Double x_0, Complex x)
inline Complex xxx_phi(Complex x, Double alpha = 1, Double x_0=0, int n=0)
{
    return exp(-alpha*(x-x_0)*(x-x_0));
}

//its fourier transform
//inline Complex xxx_phi_hat(Double alpha, Double x_0, Complex x)
inline Complex xxx_phi_hat(Complex x, Double alpha = 1, Double x_0=0, int n=0)
{
    return sqrt(Pi/alpha) * exp(-2*Pi*I*x_0*x - Pi*Pi*x*x/alpha);
}

*/


template <class ttype>
inline ttype hermite_H(ttype y,int n=0){
    if(n==0) return 1;
    else if(n==1) return 2*y;
    else{
        ttype H_a = 1, H_b = 2*y, H_c;
        int count = 1;
        do{
            H_c= 2*(y*H_b-count*H_a);
            H_a=H_b;
            H_b=H_c;
            count++;
        }while(count!=n);
        return H_c;
    }
}

inline Double Heaviside(Double t){

    if(t>0) return Double(1);
    else if(t<0) return Double(0);
    return Double(1)/2;
}

// the phi in the explicit formula.
// We allow Complex x because poles from zeta, for example, might enter, as
// xxx_phi(alpha, x_0, (pole[j]-0.5)/I)
template <class ttype>
ttype xxx_phi(ttype x, Double alpha = Double(1), Double x_0=Double(0), const char *method="hermite", int n=0,vector<Double> polynomial=vector<Double>())
{

    ttype phi;
    //nth hermite function, scaled
    if(!strcmp(method,"hermite")){
        ttype y=(x-x_0)*sqrt(2*alpha);  //center on x_0 and scale by alpha. Because of the y*y/2 we use sqrt(2*alpha) so as to get a factor of exp(-alpha(x-x_0)^2)
        phi=hermite_norm[n]*hermite_H(y,n)*exp(-y*y/2); //the nth Hermite (physicist) function
    }
    //polynomial times power sinc
    else if(!strcmp(method,"sinc")&&!polynomial.empty()){
        ttype y=(x-x_0)*alpha;  //center on x_0 and scale by alpha
        phi=pow(sinc(y),n)*Horner(polynomial,y);
    }
    //nth sinc power, scaled
    else if(!strcmp(method,"sinc")){
        ttype y=(x-x_0)*alpha;  //center on x_0 and scale by alpha
        phi=pow(sinc(y),n);
    }

    else if(!strcmp(method,"alice")){
        ttype y=(x-x_0);  //center on x_0. scale of alpha = 12 is built in below
        phi= -exp(-alpha*y*y)*sqrt(alpha*alpha*alpha)*(128.0*pow(y,10.0)*alpha*alpha*alpha*
        alpha*alpha*alpha*alpha*alpha*alpha+576.0*y*y*y*y*y*y*y*y*alpha*alpha*alpha*alpha*alpha*alpha*
        alpha*pow(log(7.0),2.0)+840.0*y*y*y*y*y*y*alpha*alpha*alpha*alpha*alpha*pow(log(7.0),
        4.0)+464.0*y*y*y*y*alpha*alpha*alpha*pow(log(7.0),6.0)+72.0*y*y*alpha*pow(log(7.0
        ),8.0)-3780.0*alpha*alpha*alpha*alpha-2880.0*y*y*y*y*y*y*y*y*alpha*alpha*alpha*alpha*
        alpha*alpha*alpha*alpha-1392.0*y*y*alpha*alpha*pow(log(7.0),6.0)-36.0*pow(log(7.0),8.0)
        -8064.0*y*y*y*y*y*y*alpha*alpha*alpha*alpha*alpha*alpha*pow(log(7.0),2.0)-6300.0*y*
        y*y*y*alpha*alpha*alpha*alpha*pow(log(7.0),4.0)+37800.0*alpha*alpha*alpha*alpha*alpha*y*y+
        3780.0*alpha*alpha*alpha*pow(log(7.0),2.0)-50400.0*alpha*alpha*alpha*alpha*alpha*alpha*y*y*y*
        y-30240.0*alpha*alpha*alpha*alpha*y*y*pow(log(7.0),2.0)-1575.0*alpha*alpha*pow(log(7.0),
        4.0)+20160.0*alpha*alpha*alpha*alpha*alpha*alpha*alpha*y*y*y*y*y*y+9450.0*alpha*alpha*alpha*
        y*y*pow(log(7.0),4.0)+30240.0*alpha*alpha*alpha*alpha*alpha*y*y*y*y*pow(log(7.0),
        2.0)+348.0*alpha*pow(log(7.0),6.0))/sqrt(pow(0.3141592653589793E1,21.0))/128.0;
        //phi:=simplify(int(exp(-Pi^2*u*u/alpha)*exp(2*Pi*I*u*x)*u*u*(u+log(7)/(2*Pi))*(u-log(7)/(2*Pi))*(u-2*log(7)/(2*Pi)
        //*(u+2*log(7)/(2*Pi))*(u-2*log(7)/(2*Pi))*(u+2*log(7)/(2*Pi))*(u-3*log(7)/(2*Pi))*(u+3*log(7)/(2*Pi)),u=-infinity..infinity));

        Double phi_of_0=-exp(0.0)*sqrt(alpha*alpha*alpha)*(-3780.0*alpha*alpha*alpha*alpha-36.0*pow(log(
        7.0),8.0)+3780.0*alpha*alpha*alpha*pow(log(7.0),2.0)-1575.0*alpha*alpha*pow(log(7.0),4.0)
        +348.0*alpha*pow(log(7.0),6.0))/sqrt(pow(0.3141592653589793E1,21.0))/128.0;


        phi=phi/phi_of_0;


    }
    return phi;
}

//its fourier transform
inline Complex xxx_phi_hat(Double u, Double alpha = Double(1), Double x_0=Double(0), const char *method="hermite", int n=0,vector<Double> polynomial=vector<Double>())
{

    Complex phi_hat;

    if(!strcmp(method,"hermite")){
        Double y=twoPi*u/sqrt(2*alpha);
        Double phi_0= hermite_norm[n]*hermite_H(y,n)*exp(-y*y/2);

        Complex I_factor;
        switch(n%4){
            case 0:
                I_factor=1;
                break;
            case 1:
                I_factor=-I;
                break;
            case 2:
                I_factor=-1;
                break;
            case 3:
                I_factor=I;
                break;
        }
        phi_hat=I_factor*sqrt(Pi/alpha)*phi_0*lcalc_expIt(-twoPi*x_0*u);
    }

    else if(!strcmp(method,"sinc")&&!polynomial.empty()){
        Double w=twoPi*u/alpha;
        Complex r=Double(1);
        Double r2=Double(1);
        for(int j = 1;j<=n-1;j++) r2*=2*j;
        r2=Pi/r2;

        //vector<Double>::size_type deg=polynomial.size()-1;
        int deg=int(polynomial.size()-1);

        Complex s=Double(0);

        //for(vector<Double>::size_type i=0;i<=deg;i++){
        for(int i=0;i<=deg;i++){
            Complex I_factor;
            switch(i%4){
                case 0:
                    I_factor=1;
                    break;
                case 1:
                    I_factor=I;
                    break;
                case 2:
                    I_factor=-1;
                    break;
                case 3:
                    I_factor=-I;
                    break;
            }
            r=r2*I_factor*polynomial[i];
            for(int j=0;j<=n;j++){
                s+=r*pow(w+n-2*j,n-1-i)*Heaviside(w+n-2*j);
                //cout << "            sinc_hat sum: " << j << " " << pow(w+n-2*j,n-1)*Heaviside(w+n-2*j)  << " :r: " << r << " " << n-j << endl;
                r=-(n-j)*r/(j+1);
            }
            r2*=(n-i-1);
        }


        phi_hat=s*lcalc_expIt(-twoPi*x_0*u)/alpha;
    }

    else if(!strcmp(method,"sinc")){
        Double w=twoPi*u/alpha;
        Double r=Double(1);
        for(int j = 1;j<=n-1;j++) r*=2*j;
        r=Pi/r;
        Double s=Double(0);
        for(int j=0;j<=n;j++){
            s+=r*pow(w+n-2*j,n-1)*Heaviside(w+n-2*j);
//cout << "            sinc_hat sum: " << j << " " << pow(w+n-2*j,n-1)*Heaviside(w+n-2*j)  << " :r: " << r << " " << n-j << endl;
            r=-(n-j)*r/(j+1);
        }
        phi_hat=s*lcalc_expIt(-twoPi*x_0*u)/alpha;
    }

    else if(!strcmp(method,"alice")){
        /*
        Double w=twoPi*u/alpha;
        Double r=Double(1);
        for(int j = 1;j<=n-1;j++) r*=2*j;
        r=Pi/r;
        Double s=Double(0);
        for(int j=0;j<=n;j++){
            s+=r*pow(w+n-2*j,n-1)*Heaviside(w+n-2*j);
//cout << "            sinc_hat sum: " << j << " " << pow(w+n-2*j,n-1)*Heaviside(w+n-2*j)  << " :r: " << r << " " << n-j << endl;
            r=-(n-j)*r/(j+1);
        }
        phi_hat=I*w*(u-log(7.)/(twoPi))*(u+log(7.)/(twoPi))*s*lcalc_expIt(-twoPi*x_0*u)/alpha;
        */

        /*
        Double y=twoPi*u/sqrt(2*alpha);
        Double phi_0= hermite_norm[n]*hermite_H(y,n)*exp(-y*y/2);

        Complex I_factor;
        switch(n%4){
            case 0:
                I_factor=1;
                break;
            case 1:
                I_factor=-I;
                break;
            case 2:
                I_factor=-1;
                break;
            case 3:
                I_factor=-I;
                break;
        }
        phi_hat=I_factor*sqrt(Pi/alpha)*phi_0*lcalc_expIt(-twoPi*x_0*u)*(u-log(7.)/(twoPi))*(u+log(7.)/(twoPi))*(u-2*log(7.)/(twoPi))*(u+2*log(7.)/(twoPi));
       */
       phi_hat=lcalc_expIt(-twoPi*x_0*u)*exp(-0.3141592653589793E1*0.3141592653589793E1*u*u/alpha)*u*u*pow(
       0.3141592653589793E1*u-log(7.0),2.0)*pow(0.3141592653589793E1*u+log(7.0),2.0)/(
       0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1*
       0.3141592653589793E1*0.3141592653589793E1*0.3141592653589793E1*
       0.3141592653589793E1*0.3141592653589793E1)*(4.0*0.3141592653589793E1*
       0.3141592653589793E1*u*u-pow(log(7.0),2.0))*(4.0*0.3141592653589793E1*
       0.3141592653589793E1*u*u-9.0*pow(log(7.0),2.0))/16.0;

       Double phi_of_0=-exp(0.0)*sqrt(alpha*alpha*alpha)*(-3780.0*alpha*alpha*alpha*alpha-36.0*pow(log(
       7.0),8.0)+3780.0*alpha*alpha*alpha*pow(log(7.0),2.0)-1575.0*alpha*alpha*pow(log(7.0),4.0)
       +348.0*alpha*pow(log(7.0),6.0))/sqrt(pow(0.3141592653589793E1,21.0))/128.0;

       phi_hat=phi_hat/phi_of_0;

       //phi_hat is simplify((exp(-Pi^2*u*u/alpha)*u*u*(u+log(7)/(2*Pi))*(u-log(7)/(2*Pi))*(u-2*log(7)/(2*Pi))*(u+2*log(7)/(2*Pi))
       //*(u-2*log(7)/(2*Pi))*(u+2*log(7)/(2*Pi))*(u-3*log(7)/(2*Pi))*(u+3*log(7)/(2*Pi)))); TIMES lcalc_expIt(-twoPi*x_0*u)

    }


    return phi_hat;
}
/***************************************************************************************************

With the above choice of phi, and the choice of alpha used in Lfind_zeros.h

Prime power side (i.e. RHS):
---------------------------

In phi_hat we call it with x=log(n)/(2 Pi). To get 10^(-DIGITS), roughly exp(-2.3*DIGITS),
we need Pi^2 (log(n)/(2Pi))^2 / alpha to be roughly 2.3 DIGITS. Solving for n, this requires n as big
as exp(sqrt(2.3 DIGITS 4alpha)). So alpha can't be too big, otherwise we get slaughtered on the number
of n needed. In Lfind_zeros.h I take alpha =2.3/DIGITS so that n needs to be as big as 100, very reasonable.

On the zeros side (LHS):
-----------------------

We want to see how many zeros we need to take so that alpha(x-x_0)^2 is >= 2.3 DIGITS.
If alpha = 2.3/DIGITS, this restricts |x-x_0|  < DIGITS.
In Lfind_zeros.h I take the number of zeros to be a bit more than 2*DIGITS*(density of zeros at largest height being considered).
density is Nmain(T+1)-Nmain(T), which is for an interval of length 2. That's okay because I store positive and
negative zeros in the array of zeros that is passed to the explicit formula. The extra 2 on the outside
is to account for the fact that we want zeros on either side of x_0, i.e. x-x_0 can be positive or negative, so
we need to double the number of zeros considered. Typically, though, half the zeros give an exponentially small amount,
and their contribution is not evaluated (those that are in the opposite halfplane). So the number of exp's
we need to evaluate, with the choice of alpha=2.3/DIGITS, is DIGITS*density. DIGITS is typically something like 15,30, or 60
and density is a handful, a number that one can count on fingers, so somewhere near one or two hundred terms on the LHS
typically get evaluated.


***************************************************************************************************/

//this should be (is) called just once, not each time test is called
//We could also speed this up by restricting n in c(n) to prime powers n,
//but in practice, I use this for small (around 100) num_coeffs, so it's not a big deal
//to overdo it a bit.

template <class ttype>
int L_function <ttype>::
dirichlet_coeffs_log_diff(int num_coeffs, Complex *c)
{

  if(my_verbose>1) cout << "#        Computing " << num_coeffs << " Dirichlet coefficients of the logarithmic derivative" << endl;

  Complex *b;
  b = new Complex[num_coeffs+1];
  int j, n, d1, ind;
  Complex total, total2, temp;


  if (what_type_L != 1 && what_type_L != -1
      && num_coeffs > number_of_dirichlet_coefficients)
  {
      cout << "Don't have enough Dirichlet coefficients." << endl;
      return 1;
  }

  //make table of prime powers from table of primes (i.e. from get_prime)
  int *isprimepower;
  isprimepower = new int[num_coeffs+1];
  for(n=0;n<=num_coeffs;n++) isprimepower[n]=0;

  j = 0;
  int p = get_prime(j);

  while (p <= num_coeffs)
  {
    n = 1;
    do
    {
      n *= p;
      isprimepower[n]=p;
    }while(double(n)*p <= num_coeffs);
    j++;
    p = get_prime(j);
  }



  b[1] = 1;

  for (n=2;n<=num_coeffs;n++) if(isprimepower[n]!=0)
  {
      total = 0.;
      total2 = 0.;

      j=1;
      p=isprimepower[n];
      do
      {
          d1 = n/j;
          if (what_type_L == -1)
              temp = b[j];
          else if (what_type_L == 1)
          {
            ind = d1 % period;
            if (ind == 0)
                ind = period;
            temp = dirichlet_coefficient[ind]*b[j];
          }
          else
              temp = dirichlet_coefficient[d1]*b[j];
          total -= temp;
          total2 += temp*LOG(d1);
          //cout << j << " : " << n << endl;
          j*=p;
       }while(j<n);
      /*
      for (j=1;j<=n/2;j++) if (n % j == 0)
      {
          d1 = n/j;
          if (what_type_L == -1)
              temp = b[j];
          else if (what_type_L == 1)
          {
            ind = d1 % period;
            if (ind == 0)
                ind = period;
            temp = dirichlet_coefficient[ind]*b[j];
          }
          else
              temp = dirichlet_coefficient[d1]*b[j];
          total -= temp;
          total2 += temp*LOG(d1);
        }
       */
       b[n] = total;
       c[n] = total2;
       if(my_verbose>4) cout << "                c[" << n << "] = " << c[n] << endl;
  }
  else{
      b[n]=Double(0);
      c[n]=Double(0);
  }

  delete [] b;
  delete [] isprimepower;
  return 0;
}

/************************************************************************************************/

//XXXXXXXXXXXXXXXXXX this can (should) move to .cc file. Does not require any L_function data.
//all moved to rhs_explicit_formula
template <class ttype>
int L_function <ttype>::
test_explicit_formula(Double x_0, Double alpha, Double *zero_table, int number_zeros, Complex *c, int num_coeffs, const char *method="hermite",int n=0,vector<Double> polynomial=vector<Double>())
{

  if(my_verbose>1) cout << "#        Entering test_explicit_formula" << endl;
  Double LHS, RHS;
  int flag;

  // COMPUTE RHS
  RHS = rhs_explicit_formula(x_0,alpha,c,num_coeffs,method,n,polynomial);

  /*** COMPUTE LHS ***/
  LHS = 0.;
  Double u=Double(0);
  //we will want to truncate this automatically
  for (int j=0;j<=number_zeros-1;j++)
  {

    u=zero_table[j]-x_0;
    u=u*u*alpha;
    if(u<2.3*DIGITS+2){
        LHS += real(xxx_phi(zero_table[j],alpha,x_0,method,n,polynomial));
        //cout << "LHS, u: " << u << " , contribution: " << exp(-u) << endl;
    }
  }



  Double DIFF=abs(LHS-RHS);
  Double RELATIVE_DIFF=DIFF/(abs(x_0)+1); // denominator to take into account
  // the fact that our test function localizes at x_0, i.e is a function of x_0-x, so
  // precision is lost due to this subtraction.

  if(RELATIVE_DIFF>xxx_max_DIFF) xxx_max_DIFF=RELATIVE_DIFF;

  if(my_verbose!=0){
      cout << "alpha = " << alpha << ", x_0 = " << x_0 << ", ";
      //cout << " x_0=" << x_0 << ",";
      cout << "DIFF = " << DIFF << ", RELATIVE_DIFF = " << RELATIVE_DIFF << ", ";
      cout << "max_relative_DIFF = " << xxx_max_DIFF << ", ";
  }


  if(xxx_max_DIFF<tolerance) DIGITS_xxx=DIGITS;
  else DIGITS_xxx=Int(-log(xxx_max_DIFF)/2.3+1); //adjust output precision for the zeros.
  if(DIGITS_xxx<2) DIGITS_xxx=2; //to prevent cout errors, never set output precision less than 2 digits.

  if (DIFF < 1e-2) //if results agree to at least two places after the decimal
  {
    flag = 0;
    if(my_verbose!=0) cout << "PASS." << endl;
  }
  else
  {
    flag = 1;
    if(my_verbose!=0) cout << "FAIL!" << endl;
  }

  return flag;
}

template <class ttype>
Double L_function <ttype>::
rhs_explicit_formula(Double x_0, Double alpha, Complex *c, int num_coeffs, const char *method="hermite",int n=0,vector<Double> polynomial=vector<Double>())
{

  if(my_verbose>1){
      cout << "#        Entering rhs_explicit_formula with method:" << method << "_" << n << endl;
      cout << "#        polynomial:";
      for(vector<string>::size_type j = 0; j != polynomial.size(); j++) {
          cout << " " << polynomial[j];
      }
      cout << endl;
  }

  Double t, t_begin, t_end, t_step;
  Double D;
  Double total;
  Double term1, term2, term3;
  int p, m, x, j;
  Double temp;
  Double RHS;
  //Complex *c;
  //int num_coeffs;

  //num_coeffs = 150; //XXXXXXXXXX should depend on test required
  //c = new Complex[num_coeffs+1]; //XXXXXXX move to L.h
  //dirichlet_coeffs_log_diff(num_coeffs, c);

  //compute the possible contribution from poles

  term1 = 0.;
  for (j=1;j<=number_of_poles;j++){
    term1 += real(xxx_phi((pole[j]-0.5)/I,alpha,x_0,method,n,polynomial));
    //cout << "pole: " << (pole[j]-.5)/I << " XXXXXX  term1: " << term1 << endl;
  }
  // compute the contribution from the Gamma factors (integral of the log diff)

  //t_step = .25*Pi/sqrt(2.3*DIGITS*alpha); //bogus formula
  t_step = .01;




  total = Double(0);

  //for our choice of phi hermite_0, this gives DIGITS precision
  if(!strcmp(method,"hermite")||!strcmp(method,"alice")){
      D = ceil(sqrt(DIGITS*2.3/alpha)/t_step) * t_step;
 D =D*4;
      t_begin = x_0 - D;
      t_end = x_0 + D;

      for (t=t_begin;t<=t_end;t=t+t_step)
      {

        temp = Double(0);

        if(!strcmp(method,"alice"))
        temp +=  12*real(log_GAMMA(gamma[j]/2 + I*t*gamma[j]+lambda[j], 1)*gamma[j]); //'1' tells log_GAMMA to compute the logarithmic derivative
        else
        for (j=1;j<=this->a;j++)
        {
          temp +=  2*real(log_GAMMA(gamma[j]/2 + I*t*gamma[j]+lambda[j], 1)*gamma[j]); //'1' tells log_GAMMA to compute the logarithmic derivative
          //temp +=  log_GAMMA(gamma[j]/2 - I*t*gamma[j]+conj(lambda[j]), 1)*gamma[j];
        }


        total = total + real(xxx_phi(t,alpha, x_0, method,n,polynomial)) * temp; //xxx_phi in this line is real, but the value is returned as a Complex,
                                                     //so we take it's real part since total is a Double
        //cout << t << "term2 " << total*t_step << endl;
      }
  }

  //use tanh sinh, or should I say, sinh sinh method
  if(!strcmp(method,"sinc")){
      for(int m=-400;m<=400;m++)
      {

        temp = Double(0);
        Double sinh_tm=sinh(t_step*m);
        t=sinh(sinh_tm)+x_0;
        for (j=1;j<=this->a;j++)
        {
          temp +=  2*real(log_GAMMA(gamma[j]/2 + I*t*gamma[j]+lambda[j], 1)*gamma[j]); //'1' tells log_GAMMA to compute the logarithmic derivative
          //temp +=  log_GAMMA(gamma[j]/2 - I*t*gamma[j]+conj(lambda[j]), 1)*gamma[j];
        }
        total = total + real(xxx_phi(t,alpha, x_0, method,n,polynomial)) * temp * cosh(sinh_tm) * cosh(t_step*m); //xxx_phi in this line is real, but the value is returned as a Complex,
                                                                                                     //so we take it's real part since total is a Double
        //cout << t << "term2 " << total*t_step << endl;
      }
  }

  term2 = t_step*total + 2*log(Q)*real(xxx_phi_hat(Double(0),alpha,x_0,method,n,polynomial));
  term2 = 1/(2*Pi) * term2;

  //compute the contribution from the Dirichlet coefficients

  x = num_coeffs;
  //extend_prime_table(x);

  j = 0;
  p = get_prime(j);
  term3 = Double(0);

  while (p <= x)
  {
    //m = p;
    m = 1;
    do
    {
      m *= p;
      //temp = 2*real(c[m]*xxx_phi_hat(alpha, x_0, LOG(m)/(2*Pi)));// + conj(c[m])*xxx_phi_hat(alpha, x_0, -LOG(m)/(2*Pi));
      temp = 2*real(c[m]*xxx_phi_hat(LOG(m)/(2*Pi),alpha,x_0,method,n,polynomial));// + conj(c[m])*xxx_phi_hat(alpha, x_0, -LOG(m)/(2*Pi));
      //cout << m << "           alice phi_hat: " << xxx_phi_hat(LOG(m)/(2*Pi),alpha,x_0,method,n) << endl;
      term3 += two_inverse_sqrt(m)*temp;
    }while(double(m)*p <= x);
    j++;
    p = get_prime(j);
  }
  term3 = 1/(4*Pi) * term3;

  /*** COMPUTE RHS ***/

  RHS = term1 + term2 - term3;


  /*** Display Results ***/

  if(my_verbose > 2)
  {
    cout << endl << endl;
    cout << "#            *** RHS Explicit Formula for L ***" << endl;
    cout << "#            alpha = " << alpha << endl;
    cout << "#            x_0 = " << x_0 << endl;
    cout << "#            D = " << D << endl;
    cout << "#            TERM 1:  " << term1 << endl;
    cout << "#            TERM 2:  " << term2 << endl;
    cout << "#            TERM 3:  " << term3 << endl;
    cout << "#            RHS:  " << RHS << endl;
  }

  return RHS;
}

template <class ttype>
int L_function <ttype>::
plot_explicit_formula(Double alpha, Double x=Double(0), Double x2=Double(100), Double step_size=.01, const char *xxx_phi_method="hermite", int num_coeffs=-1, Double *rhs_store=NULL){


    if(my_verbose>0) cout << "#    Enter plot_explicit_formula with:"
        << alpha << " "
        << x << " "
        << x2 << " "
        << step_size << " "
        << xxx_phi_method << " "
        << num_coeffs << endl;

    //cout <<"#   STEPSIZE = " << step_size << endl;
    //cout <<"#   x = " << x << endl;
    //cout <<"#   x2 = " << x2 << endl;

    int number_of_log_diff_coeff;
    int n=0; // as in hermite_n, or sinc_n = sin(alpha x)^n/(alpha x)^n. hermite_0 is const exp(-alpha x^2).


    //========================================================================================
    // ========  parse the xxx_phi_method string to determine which method to use  ===========
    //========================================================================================

    //const char *ptr = xxx_phi_method;
    char method[100];
    vector<Double> polynomial;

    vector<string> method_parsed = split(xxx_phi_method, '_');

    //if(my_verbose>4) for(vector<string>::size_type j = 0; j != method_parsed.size(); j++) {
    //    cout << "#                    explicit formula method, parsed: " << method_parsed[j] << endl;
    //}

    strcpy(method,method_parsed[0].c_str());
    if(method_parsed.size()>1) n=atoi(method_parsed[1].c_str());
    if(my_verbose>0) cout << "#    plot explicit formula, phi method = "<< method << "_" << n << endl;
    if(method_parsed.size()>2){
        //parse the polynomial portion of the string. For example, (sin(x)/x)^4(a+bx+cx^2)  would be notated as sinc_4_a,b,c
        polynomial= split_Double(method_parsed[2].c_str(),',');
        if(my_verbose>0){
            cout << "#    polynomial:";
            for(vector<string>::size_type j = 0; j != polynomial.size(); j++) {
                cout << " " << polynomial[j];
            }
            cout << endl;
        }
    }

    //========================================================================================
    // ========  decide how many terms of the Dirichlet series of L'/L to use  ===============
    //========================================================================================

    double N_as_double; //as double because we might exceed INT_MAX
    if(!strcmp(method,"hermite")){
        N_as_double = lcalc_to_double(exp(sqrt(2.3*DIGITS*4*alpha)));
        if(my_verbose>0) cout << "#    plot explicit formula, N_as_double = " << N_as_double << endl;
    }

    //if(!strcmp(method,"sinc")||!strcmp(method,"alice")){
    if(!strcmp(method,"sinc")){
        N_as_double = lcalc_to_double(exp(n*alpha));
        //if(my_verbose>1) cout << "#        plot explicit formula, N_as_double = " << N_as_double << endl;
    }

    if(!strcmp(method,"alice")){
        N_as_double = 1e7;
    }

    if(num_coeffs>-1) N_as_double=double(num_coeffs); // if num_coeffs is specified, limit to that

    if (what_type_L != 1 && what_type_L != -1 && N_as_double> number_of_dirichlet_coefficients){
        cout << "# WARNING Don't have enough (" << N_as_double << ") Dirichlet coefficients. " << endl;
        cout << "# WARNING Will use the maximum, "<< number_of_dirichlet_coefficients <<" , available." << endl;
        number_of_log_diff_coeff = number_of_dirichlet_coefficients;


        Double tmp= log(Double(number_of_log_diff_coeff))/(2*Pi);
        //cout << "WARNING These can only give the rhs of the explicit formula to within roughly: " << exp(-Pi*Pi/alpha*tmp*tmp) << endl;
        cout << "# WARNING These can only give the rhs of the explicit formula to within roughly: " << abs(xxx_phi_hat(tmp,alpha,Double(0),method,n,polynomial)) << endl;
    }

    else if ((what_type_L == 1 || what_type_L == -1)&&N_as_double>100000000){
        cout << "# WARNING Will use p^k < 100000000" << endl;
        number_of_log_diff_coeff=100000000; //for zeta or Dirichlet L-functions I've set the max number of terms at 1000000
        Double tmp= log(Double(number_of_log_diff_coeff))/(twoPi);
        //cout << "These can only give compute the rhs of the explicit formula to within: " << exp(-Pi*Pi/alpha*tmp*tmp) << endl;
        //cout << "#tmp: " << tmp;
        //cout << "#alpha: " << alpha;
        //XXXXXXXXXXXXXXX  needs to be fixed to take into account the slow rate of decay of the prime power sum, i.e one term does not dominate the tail
        cout << "# WARNING These can only give the rhs of the explicit formula to within roughly: " << abs(xxx_phi_hat(tmp,alpha,Double(0),method,n,polynomial)) << endl;

    }

    else number_of_log_diff_coeff=Int(N_as_double);
    //else number_of_log_diff_coeff=Int(exp(sqrt(2.3*DIGITS*4*alpha)));

    if(my_verbose>1){
        cout << "#        plot_explicit_formula called with: alpha = "<< alpha << ", " << x << " <= x0 <= " << x2 << endl;
        cout << "#        Will use: " << number_of_log_diff_coeff << " terms of the logarithmic derivative" << endl;
    }

    //========================================================================================
    // ========  compute the coefficients of the Dirichlet series of -L'/L  ==================
    //========================================================================================

    Complex *xxx_log_diff_coeffs;
    xxx_log_diff_coeffs = new Complex[number_of_log_diff_coeff+1];

    dirichlet_coeffs_log_diff(number_of_log_diff_coeff, xxx_log_diff_coeffs);

    Double t;


    //========================================================================================
    // ========  increment t from x to x2, and evaluate sum_gamma phi(gamma-t) via the rhs ===
    // ========  of the explicit formula, i.e. from the prime powers side                  ===
    //========================================================================================

    int count=0;
    t=x;
    do{


        Double rhs=rhs_explicit_formula(t,alpha,xxx_log_diff_coeffs,number_of_log_diff_coeff,method,n,polynomial);

        if(rhs_store){
            rhs_store[count]=rhs;
            //cout << method<< "_" << n << " ::::::::: " << t << " " << rhs_store[count] << endl;
            count++;
            t+=step_size;
        }
        else{
            cout << t;
            cout << " " << rhs << endl;
            t+=step_size;
        }


    }while(t<=x2);


    //delete [] rhs_store;
    delete [] xxx_log_diff_coeffs;
    return 0;
}


#endif
