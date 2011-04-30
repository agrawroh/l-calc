//#include <Lfunction/L.h>
#include "L.h"

//a simple program illustrating a few features of the L_function class.

int main (int argc, char *argv[])
{


    initialize_globals(); //initialize global variables. This *must* be called.

    Double x,y;
    Complex s;

    L_function<int> zeta; //default L-function is the Riemann zeta function

    L_function<int> L4; //will be assigned below to be L(s,chi_{-4}), 
                        //chi{-4} being the real quadratic character mod 4

    L_function<Complex> L5; //will be assigned below to be L(s,chi), 
                            //with chi being a complex character mod 5


//==================== Initialize the L-functions ==================================
    //one drawback of arrays- the index starts at 0 rather than 1
    //so each array below is declared to be one entry larger than it is. 
    //I prefer this so that referring to the array elements is more 
    //straightforward, for example coeff[1] refers
    //to the first Dirichlet cefficient rather than coeff[0]. But to make up
    //for this, one needs to have a bogus entry (such as 0) at the start of each array.

    int coeff_L4[5] = {0,1,0,-1,0}; //the Dirichlet coefficients, periodic of period 4.
    Double gamma_L4[2] = {0,.5}; //the gamma factor Gamma(gamma s + lambda) is Gamma(s/2+1/2)
    Complex lambda_L4[2] = {0,.5}; //the lambda
    Complex pole_L4[1] = {0}; //no pole
    Complex residue_L4[1] = {0}; //no residue


    L4=L_function<int>("L4",1,4,coeff_L4,4,sqrt(4/Pi),1,1,gamma_L4,lambda_L4,0,pole_L4,residue_L4);

    // "L4" is the name of the L-function 
    //  1 - what_type, 1 stands for periodic Dirichlet coefficients
    //  4 - N_terms, number of Dirichlet coefficients given
    //  coeff_L4  - array of Dirichlet coefficients
    //  4 - period (0 if coeffs are not periodic)
    //  sqrt(4/Pi) - the Q^s that appears in the functional equation 
    //  1 - sign of the functional equation
    //  1 - number of gamma factors of the form Gamma(gamma s + lambda), gamma = .5 or 1
    //  gamma_L4  - array of gamma's (each gamma is .5 or 1)
    //  lambda_L4  - array of lambda's (given as complex numbers)
    //  0 - number of poles. Typically there won't be any poles.
    //  pole_L4 - array of poles, in this case none
    //  residue_L4 - array of residues, in this case none

    //  Note: one can call the constructor without the last three arguements when number of poles = 0
    //  as in:
    L4 = L_function<int>("L4",1,4,coeff_L4,4,sqrt(4/Pi),1,1,gamma_L4,lambda_L4);


    Complex coeff_L5[6] = {0,1,I,-I,-1,0}; 

    Complex gauss_sum=0.;
    for(int n=1;n<=4; n++) gauss_sum=gauss_sum+coeff_L5[n]*exp(n*2*I*Pi/5);


    L5=L_function<Complex>("L5",1,5,coeff_L5,5,sqrt(5/Pi),gauss_sum/(I*sqrt(5.)),1,gamma_L4,lambda_L4);
    // "L5" is the name of the L-function 
    //  1 - what_type, 1 stands for periodic Dirichlet coefficients
    //  5 - N_terms, number of Dirichlet coefficients given
    //  coeff_L5  - array of Dirichlet coefficients
    //  5 - period (0 if coeffs are not periodic)
    //  sqrt(5/Pi), the Q^s that appears in the functional equation 
    //  gauss_sum/sqrt(5) - omega of the functional equation
    //  1 - number of gamma factors of the form Gamma(gamma s + lambda), gamma = .5 or 1
    //  gamma_L4  - L5 has same gamma factor as L4
    //  lambda_L4  - ditto


    //=========== print basic data for the L-function ========================================
    zeta.print_data_L();
    L4.print_data_L();
    L5.print_data_L();



    //=========== print some L-values ========================================

    x= .5; y =0.;
    cout << "zeta" << x+I*y << " = " << zeta.value(x+I*y) << endl;
    cout << "L4" << x+I*y << " = " << L4.value(x+I*y) << endl;
    cout << "L5" << x+I*y << " = " << L5.value(x+I*y) << endl;

    x= 1; y =0.;
    cout << "L4" << x+I*y << " = " << L4.value(x+I*y) << endl;
    cout << "L5" << x+I*y << " = " << L5.value(x+I*y) << endl;

    //=========== find and print some zeros ========================================
    //find zeros of zeta up to height 100 taking steps of size .1, looking for sign
    //changes on the critical line. Some zeros can be missed in this fashion.
    //First column gives the imaginary part of the zeros.
    //Second column outputted is related to S(T) and should be small on average
    //(larger values means zeros were missed).
    zeta.find_zeros(Double(0),Double(100),Double(.1));

    //find the first 100 zeros of zeta. This also verifies RH and does not omit zeros.
    zeta.find_zeros(100); //will *not* look for zeros below the 
                                      //real axis as they come in conjugate pairs


    //do same for L4 and L5
    L4.find_zeros(Double(0),Double(100),Double(.1));
    L4.find_zeros(100);

    L5.find_zeros(Double(0),Double(100),Double(.1));
    L5.find_zeros(100); //will look for zeros above and below the real axis since is not self-dual



}
