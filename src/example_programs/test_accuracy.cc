//#include <Lfunction/L.h>
#include "L.h"
#include <time.h>
//#ifdef INCLUDE_PARI
//    #include <pari/pari.h>
 //   #undef init
//#endif


void wait(double x);
char *time_stamp();
void start_timer();
void end_timer();
void end_timer(char *message);
clock_t start,end;      //for timing
time_t now;             //for printing the current time

template <class ttype, class ttype2>
int compare_accuracy(ttype a, ttype2 b){
    Double x=abs(a-b)/(abs(a)+abs(b));
    if(x<tolerance) return DIGITS;
    return -Int(log(x)* 0.4342944819032518); // 0.4342944819032518 = 1/log(10)
}

template <class ttype>
int compare_to_zero(ttype a){
    Double x=abs(a);
    return -Int(log(x)* 0.4342944819032518); // 0.4342944819032518 = 1/log(10)
}

//a simple program illustrating a few features of the L_function class.

int main (int argc, char *argv[])
{

    initialize_globals(); //initialize global variables. This *must* be called.


    Double x;
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
    //for this, we insert a bogus entry (such as 0) at the start of each array.

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

    L5=L_function<Complex>("L5",1,5,coeff_L5,5,sqrt(Double(5)/Pi),gauss_sum/(I*sqrt(Double(5))),1,gamma_L4,lambda_L4);
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


    //cout << setprecision(DIGITS);
    //x=Double(2); cout << "2: "<< x << endl;
    //x=2.; cout << "2: "<< x << endl;
    //x=0.; cout << "0: "<< x << endl;
    //x=.5; cout << ".5: "<< x << endl;
    //int nnn=3;
    //x=nnn+.5; cout << "n+.5: "<< x << endl;
    //x=.25; cout << ".25: "<< x << endl;
    //x=Pi; cout << "Pi: "<< x << endl;

    x=zeta.initialize_gram(0.);
    for(int n=1;n<=100;n++){
        cout << x << " " << zeta.value(.5+I*x,0,"rotated pure") << endl;
        x=zeta.next_gram(x);
    }

    x=L5.initialize_gram(0.);
    for(int n=1;n<=100;n++){
        cout << x << " " << L5.value(.5+I*x,0,"rotated pure") << endl;
        x=L5.next_gram(x);
    }

    cout << "======== Testing basic functions and constants ==================================" << endl;

cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";
cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX\n";

    cout << "  lcalc_cos(Pi) agrees with -1 to: " << compare_accuracy(lcalc_cos(Pi),-1) << " digits." <<endl;
    cout << "  lcalc_cos(2)^2 +lcalc_sin(2)^2 agrees with 1 to: " << compare_accuracy(lcalc_cos(2)*lcalc_cos(2)+lcalc_sin(2)*lcalc_sin(2),1) << " digits." <<endl;
    cout << "  GAMMA(2) agrees with 1 to: " << compare_accuracy(GAMMA(Double(2)),1) << " digits." <<endl;

    cout << "    --------------- timing tests basic functions" << endl;

    Double tmp=0.;

    x=twoPi/1000000;
    cout << "      plain cos: time to do 'for(int j=1;j<=1000000;j++)tmp+=cos(j*x);': ";
    tmp=0.;
    start_timer();
    for(int j=1;j<=1000000;j++) tmp+=cos(j*x);
    end_timer();
    cout << "      Accuracy: " << compare_to_zero(tmp) << " digits."  << endl;
    cout << "      REMARK: above we expect to lose an extra 6 digits just for magnifying error in x by j=10^6" << endl;


    cout << "      let x=twoPi/1000000 ";
    cout << "      lcalc_cos: time to do 'for(int j=1;j<=1000000;j++) tmp+=lcalc_cos(j*x);': ";
    start_timer();
    for(int j=1;j<=1000000;j++) tmp+=lcalc_cos(j*x);
    end_timer();
    cout << "      Accuracy: " << compare_to_zero(tmp) << " digits."  << endl;







    cout << "======== Testing Incomplete Gamma function routines==================================" << endl;
    cout << "    --------------- comparing continued fraction vs Temme's asymptotics near the transition zone" << endl;
    Complex z,w;
    z=1+100*I; w = z+9; cout << "      cfrac_GAMMA(1+100*I, 10+100*I) agrees with temme to: " << compare_accuracy(cfrac_GAMMA(z,w) ,Q(z,w)*GAMMA(z,w)) << " digits." << endl;
    z=1+1000*I; w = z+9; cout << "      cfrac_GAMMA(1+1000*I, 10+1000*I) agrees with temme to: " << compare_accuracy(cfrac_GAMMA(z,w) ,Q(z,w)*GAMMA(z,w)) << " digits." << endl;
    z=1+10000*I; w = z+9; cout << "      cfrac_GAMMA(1+10000*I, 10+10000*I) agrees with temme to: " << compare_accuracy(cfrac_GAMMA(z,w) ,Q(z,w)*GAMMA(z,w)) << " digits." << endl;
    z=1+100000*I; w = z+9; cout << "      cfrac_GAMMA(1+100000*I, 10+10000*I) agrees with temme to: " << compare_accuracy(cfrac_GAMMA(z,w) ,Q(z,w)*GAMMA(z,w)) << " digits." << endl;


    cout << "======== Testing Dirichlet L-functions ==================================" << endl;
    cout << "  L(1,chi_{-4})  agrees with Pi/4 to: " << compare_accuracy(L4.value((Double)1),Pi/4) <<  " digits." << endl;


    x=str_to_Double("-1.46035450880958681288949915251529801246722933101258149054288608782553052947450062527641937546335681951449637467986952958389234371035889426181923283975");

#ifdef INCLUDE_PARI
    cout << "======== Testing an elliptic curve L-function ==================================" << endl;

    char a1[2]; strcpy(a1,"0");
    char a2[2]; strcpy(a2,"0");
    char a3[2]; strcpy(a3,"0");
    char a4[2]; strcpy(a4,"4");
    char a6[2]; strcpy(a6,"0");

    pari_init_opts(400000000,2,INIT_DFTm); // the last option is to prevent
    //pari from giving its interrupt signal when its elliptic curve a_p
    //algorithm is called and interrupted with ctrl-c.
    L_function<Double> L32(a1,a2,a3,a4,a6,10000);
    //cout << "  L_{32A}(1/2)  agrees with GAMMA(1/2)GAMMA(1/4)/GAMMA(3/4)/8 to: " << compare_accuracy(L32.value(Double(1)/Double(2)),GAMMA(Double(1)/Double(2))*GAMMA(Double(1)/Double(4))/GAMMA(Double(3)/Double(4))/Double(8)) <<  " digits." << endl;
    //the formula given comes from computing the real period and applying bsd.
    //see Dave Husemoller elliptic curves, prop 6.2, page 185.
#endif

    cout << "======== Testing zeta ==================================" << endl;
    x=str_to_Double("-1.46035450880958681288949915251529801246722933101258149054288608782553052947450062527641937546335681951449637467986952958389234371035889426181923283975");
    cout << "  zeta(.5) agrees with maple value ";
    if(DIGITS-DIGITS2>150) cout << "(stored 150 digits) ";
    cout << "to: " << compare_accuracy(zeta.value(.5),x) << " digits." << endl;

    Complex L1, L2;
    Double t=1;
    cout << "  Comparing gamma sum method vs riemann sum method for zeta(s):" << endl;
    do{
        s=-3.5+I*t;
        do{
            s+=.5;
            L1 = zeta.value(s,0,"pure","Gamma sum");
            L2 = zeta.value(s,0,"pure","Riemann sum");
            cout << setprecision(10);
            cout << "s: " << s << " Agreement: " << compare_accuracy(L1,L2) <<  " digits." << endl;
        }while(real(s)<3.5);
        t*=10;
    }while(t<1e3);


    t=10000;
    do_blfi=false;
    cout << "  Comparing Riemann Siegel formula vs gamma sum method for zeta(1/2+it):" << endl;
    do{
        s=Complex(.5,t);
        L1 = zeta.value(s,0,"pure","Gamma sum");
        L2 = zeta.value(s,0,"pure");
        cout << "t: " << t << " Agreement: " << compare_accuracy(L1,L2) <<  " digits." << endl;
        t*=10;
    }while(t<1e11);



    t=1e6; Double t2;
    input_mean_spacing=.1;
    cout << "  Comparing Riemann Siegel formula vs Band Limited Interpolation for zeta(1/2+it):" << endl;
    do{

        t2=t;
        L1=L2=0;
        do_blfi=false;
        start_timer();
        for(int j=1;j<=10000;j++){
            t2+=input_mean_spacing;
            L1 += zeta.value(Complex(.5,t2),0,"pure");
        }
        end_timer();
        t2=t;
        do_blfi=true;
        start_timer();
        for(int j=1;j<=10000;j++){
            t2+=input_mean_spacing;
            L2 += zeta.value(Complex(.5,t2),0,"pure");
        }
        end_timer();
        cout << "t: " << t << " " << L1 << " " << " " << L2 << " Agreement: " << compare_accuracy(L1,L2) <<  " digits." << endl;
        t*=10;
    }while(t<1e17);


}

//-----Timing functions--------------------------------------------------------

void wait(double x)
{
    start=clock();
    do
    {
        end=clock();
    }while(end-start<x*CLOCKS_PER_SEC);
}

char *time_stamp()
{
    time(&now);
    return ctime(&now);
}

void start_timer()
{
    //cout << "Starting timer at " << time_stamp();
    start=clock();
}

void end_timer()
{
    end=clock();
    //cout << "Ending timer at " << time_stamp();
    //cout << "  Took " << (end-start)*1./CLOCKS_PER_SEC << " seconds." << endl;
    printf("%4lf seconds\n",lcalc_to_double((end-start)*1./CLOCKS_PER_SEC));
    //cout << lcalc_to_double((end-start)*1./CLOCKS_PER_SEC) << " seconds" <<endl;
}

