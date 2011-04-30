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

#include "Lglobals.h"


//-----Global variables---------------------------------------


int my_verbose=0;       // verbosity level: 0 means no verbose

int DIGITS, DIGITS2; // precision and sacrifice
int DIGITS3; // how many digits to output
int DIGITS_xxx=0; // how many digits to output as determined by the explicit formula
Double xxx_max_DIFF=0.; //maximum difference recorded comparing lhs to rhs of the explicit formula

Double tolerance;
Double tolerance_sqrd;
Double tolerance2;
Double tolerance3;

int global_derivative; //which derivative to compute
int max_n; // largest n used in a dirichlet series while computing an L-value.

Double A; //controls the 'support' of g(w) in the Riemann sum method
Double incr; //the increment in the Riemann sum method
Double tweak; //used in value_via_Riemann_sum to play with the angle

Double *LG;             // lookup table for log(n)
Double *two_inverse_SQUARE_ROOT;    // lookup table for sqrt(n)
int number_sqrts=0;     // how many sqrt(n)'s to store
int number_logs=0;      // how many log(n)'s to store

int *prime_table;
int number_primes=0;

Double *bernoulli;             // lookup table for bernoulli numbers

Double rs_remainder[40][72]; //taylor coefficients for Riemann Siegel correction terms

Double hermite_norm[201]; // stores 1/sqrt(2^n*n!*sqrt(Pi)). Used in explicit formula, Hermite test function


//bool kronecker_table_available=false;
//int **kronecker_table; //lookup table for the kronecker function
//int kronecker_bound; //how large a table to make


bool print_warning=true;

Long my_LLONG_MAX; //equals LLONG_MAX defined in limits.h. For
//reasons doing with compilation bugs, I determine it myself

Double Pi;
Double twoPi;
Double one_over_twoPi;
Double log_2Pi;
Complex I;

bool only_use_dirichlet_series=false; //whether to compute just using the Dirichlet series
int N_use_dirichlet_series; //if so, how many terms in the Dirichlet series to use.

//------Incomplete gamma function global variables--------------

Complex last_z;         // the last z to be considered in inc_GAMMA
Complex last_w;         // the last w to be considered in inc_GAMMA
Complex last_comp_inc_GAMMA; // g(last_z,last_w)

Complex last_z_GAMMA;  //the last z to be considered in GAMMA
Complex last_log_G;    //the last log(GAMMA(z));

Double temme_a[1002],temme_g[501];
//used in Temme's asymptotic expansion of the
//incomplete gamma function
//XXXX might need more terms if I go to higher precision

//----------variables related to my cosine function --------------

Double *cos_taylor; //table of taylor coefficients for cosine
int cos_taylor_arraysize;
int number_cos_taylor_terms; //the number of taylor terms per series. Should be even.

Double one_over_cos_taylor_arraysize;
Double twoPi_over_cos_taylor_arraysize;


const Double sin_cof[]={1,-Double(1)/6,Double(1)/120,-Double(1)/5040,Double(1)/362880,-Double(1)/39916800};
const Double sin_tol=1.e-5;
// Riemann Siegel band limited interpolation arrays and variables ----------

bool try_use_blfi=false; //whether to compute dirichlet series using blf interpolation
Double blfi_interval_length;
Complex* block_value;


//========================= older blfi variables to remove once the new one works
bool do_blfi=false;
const Double sinh_mult_fac=Double(1)/2;
const int sin_terms=2;


const Double blfi_block_growth=2; // how fast blfi blocks grow as we traverse the main sum, keep as is for now
const Double beta_fac_mult=6;  // controls the density of blfi sampling and the number of blfi terms needed
const Double blfi_fac=.085;  // efficiency of the blfi interpolation sum relative to an RS sum of same length
const Double pts_array_fac=10000;

const int rs_blfi_N=3;

Double *klog0; //log(k) at the beginning
Double *klog2; //log(k) at the end if needed
Double *ksqrt0; // 1/sqrt(k) at the beginning
Double *ksqrt2;// 1/sqrt(k) at the end if needed
int *num_blocks; // number of blocks
int *size_blocks;// size of blocks
Double *trig; // stores correction terms
Double *zz; // store powers of fractional part

Double **klog1; //log(k) in the middle if needed
Double **ksqrt1; // 1/sqrt(k) in the middle if needed
Double **klog_blfi; //initial term
Double **qlog_blfi; //band-width
Double **piv_org; //original pivot
Double **bbeta; //beta
Double **blambda; //lambda
Double **bepsilon; //epsilon
Double **arg_blfi; //arg_blfi
Double **inv_arg_blfi; //inv_arg_blfi

Double ***qlog_blfi_dense; // log(1+k/v) terms
Double ***qsqrt_blfi_dense; // 1/sqrt(1+k/v)
int ***blfi_done_left; //block done or not
int ***blfi_done_right; //block done or not
Double ***blfi_val_re_left; //real value of block
Double ***blfi_val_re_right; //real value of block
Double ***blfi_val_im_left; //imag value of block
Double ***blfi_val_im_right; //imag value of block

int length_org=0; // length of the main sum
int length_split=0; // length of the portion of the main sum to be evaluated directly
int lgdiv=0; // number of divisions of the main sum into intervals of the form [N,2N)
int max_pts=0; // max number of interpolation points allowed
int range=0; // number of blfi interpolation points needed
int blfi_block_size_org=0; // starting length of the blfi block
int total_blocks=0;

Double bc=0;
Double bc2=0;
Double kernel_fac=0;
Double ler=0;
Double mult_fac=0;
Double approx_blfi_mean_spacing=0;
Double interval_length=0;
Double error_tolerance=0;
Double input_mean_spacing=0;
Double input_mean_spacing_given=.01; //rough anticipated sampling distance between t's for zeta(1/2+It).
//Should be set if known.

//------------------------------------------------------------------



//-----intializing and cleaning up routines----------------------

//n is the number of log(n)'s to precompute
//call this function after changing the DIGITS global variable
void initialize_globals(int n){

//#ifdef INCLUDE_PARI

    //pari_init_opts(400000000,2,INIT_DFTm); // the last option is to prevent
   ////pari from giving its interrupt signal when its elliptic curve a_p
   ////algorithm is called and interrupted with ctrl-c.

//#endif


#if defined(USE_MPFR) || defined(USE_MPFRCPP)

    if(!DIGITS) DIGITS=100;
    if(!DIGITS2) DIGITS2=3;
    DIGITS3=DIGITS-DIGITS2;
    //__mpfr_default_fp_bit_precision=Int(DIGITS*log(10.)/log(2.));
    //mpfr_set_default_prec(__mpfr_default_fp_bit_precision);
    #ifdef USE_MPFR
        mpfr_set_default_prec(Int(DIGITS*log(10.)/log(2.)));
    #endif
    #ifdef USE_MPFRCPP
cout << "setting output precision to: " << Int(DIGITS*log(10.)/log(2.)) << endl;
        Double::set_default_prec(Int(DIGITS*log(10.)/log(2.)));
    #endif
    reset(Pi);
    reset(twoPi);
    reset(one_over_twoPi);
    reset(log_2Pi);
    reset(twoPi_over_cos_taylor_arraysize);
    reset(one_over_cos_taylor_arraysize);
    reset(I);
    reset(tolerance);
    reset(last_z);
    reset(last_w);
    reset(last_comp_inc_GAMMA);
    reset(last_z_GAMMA);
    reset(last_log_G);
    loop(i,0,1002) reset(temme_a[i]);
    loop(i,0,501) reset(temme_g[i]);
    loop(i,0,39) loop(j,0,71) reset(rs_remainder[i][j]);
    mpfr_const_pi(Pi.get_mpfr_t(), __gmp_default_rounding_mode);
    //twoPi=2*Pi;
    //one_over_twoPi=Double(1)/twoPi;
    //log_2Pi=log(twoPi);
#else
    if(!DIGITS){
        typedef std::numeric_limits< Double > tmp_Double;
        DIGITS=tmp_Double::digits10;
        if(my_verbose>0) cout << "#    DIGITS set to " << DIGITS << endl;
    }
    if(!DIGITS2) DIGITS2=3;
    DIGITS3=DIGITS-DIGITS2;


//unsigned int oldcw; //Bailey has some comment about an x86 fix, but I'm not sure if it is needed
//fpu_fix_start(&oldcw);

    //Pi= 4*atan((Double)1);
    Pi= 4*atan(Double(1));
    //twoPi=2*Pi;
    //log_2Pi= log((Double)2*Pi);
    //log_2Pi= log(twoPi);

#endif

    twoPi=2*Pi;
    one_over_twoPi=Double(1)/twoPi;
    log_2Pi=log(twoPi);

    cout.precision(DIGITS);
    if(my_verbose>0){
        cout << "#    Pi set to " << Pi << endl;
        cout << "#    log_2Pi set to " << log_2Pi << endl;
    }

    I=Complex(0,1);
    if(my_verbose>0) cout << "#    I set to " << I << endl;

    tolerance=pow(Double(.1),DIGITS);
    tolerance_sqrd=tolerance*tolerance;
    tolerance2=pow(Double(.1),(DIGITS-DIGITS2+1));
    tolerance3=pow(Double(.1),(DIGITS3+1));

    A=1./(16*Pi*Pi*23*DIGITS/10.);
    incr=2*Pi*.5/(log(10.)*DIGITS); //.5 here is current v in Riemann sum method 
    tweak=1;

    last_z=1;
    last_w=0;
    last_comp_inc_GAMMA=0;
    last_z_GAMMA=1;
    last_log_G=0;

    global_derivative=0;
    max_n=1;

    numeric_limits<long long> ll;
    my_LLONG_MAX = ll.max(); //used once, but in the elliptic curve module Lellipitic.cc
    //pari defines max #define max(a,b) ((a)>(b)?(a):(b)) and this causes problems with ll.max
    //so I put it here as a global.
    //also look at: http://www.google.com/answers/threadview?id=334912
    //Re: g++ problems with strtoll (LLONG_MIN, LLONG_MAX, ULLONG_MAX not defined)
    //for an explanation of why we just don't call LLONG_MAX (compiler inconsistencies).

    int j,k;
    Double r,x,dsum;

    number_logs=n;
    if(LG) delete[] LG;
    LG = new Double[n+1];
    for(k=1;k<=n;k++) LG[k]=log((Double)k);

    number_sqrts=n;
    if(two_inverse_SQUARE_ROOT) delete[] two_inverse_SQUARE_ROOT;
    two_inverse_SQUARE_ROOT = new Double[n+1];
    for(k=1;k<=n;k++) two_inverse_SQUARE_ROOT[k]=2/sqrt((Double)k);

    //initilialize the cosine taylor series array
    initialize_cos_array();

    //I break this up into 8 routines because g++ crashes when trying to optimize what
    //would otherwise be a very long (though simple!) routine
    initialize_rs_remainder1();
    initialize_rs_remainder2();
    initialize_rs_remainder3();
    initialize_rs_remainder4();
    initialize_rs_remainder5();
    initialize_rs_remainder6();
    initialize_rs_remainder7();
    initialize_rs_remainder8();

    if(bernoulli) delete[] bernoulli;
    bernoulli= new Double[DIGITS+1];


    bernoulli[0]=1;
    for(k=1;k<=DIGITS;k++){
        r=k+1; x=0.;
        for(j=1;j<=k;j++){
            r=r*(k+1-j)*1./(j+1);
            x=x-bernoulli[k-j]*r;
        }
        bernoulli[k]=x/(k+1);
    }

    hermite_norm[0]=1/sqrt(sqrt(Pi));
    for (int n=1; n<=200;n++)
        hermite_norm[n]=hermite_norm[n-1]/sqrt(Double(2*n));


    //XXXXX this should depend on DIGITS
    temme_a[1]=1;

    temme_a[2]=(Double) 3;
    temme_a[2]=1/temme_a[2]; //doing temme_a[2]=1./3 causes problems with long doubles or multiprecision.
    temme_a[3]=(Double)36;
    temme_a[3]=1/temme_a[3];

    for (int i=4;i<=281;i++) {
        dsum=0.;
        for (int j=2;j<=(i-1);j++)
            dsum+=j*temme_a[j]*temme_a[i-j+1];
        temme_a[i]=(temme_a[i-1]-dsum)/((Double) (i+1));
    }

    Double df=1;
    for(int i=1;i<=140;i++){
        df*=(2*i+1);
        temme_g[i]=(1-2*(i%2))*df*temme_a[2*i+1]; //gives overflow error when i=150;
        //temme_g[i]=(1-2*(i%2))*dfac(2*i+1)*temme_a[2*i+1]; //gives overflow error when i=150;
    }
    temme_g[0]=1;

    extend_prime_table(100);


}

void delete_globals(){

    if(LG) delete [] LG;
    if(two_inverse_SQUARE_ROOT) delete [] two_inverse_SQUARE_ROOT;
    if(bernoulli) delete [] bernoulli;
    if(prime_table) delete [] prime_table;
    if(cos_taylor) delete [] cos_taylor;
    //if(hermite_norm) delete [] hermite_norm;

}

void extend_LG_table(int m){


    int n;

    Double *tmp_LG; //used to copy over the old values
    tmp_LG = new Double[number_logs+1];
    for(n=1;n<=number_logs;n++) tmp_LG[n]=LG[n];

    delete [] LG;
    int new_number_logs=(int)(1.5*m); // extend table by 50 percent

    LG = new Double[new_number_logs+1];
    for(n=1;n<=number_logs;n++) LG[n]=tmp_LG[n];
    for(n=number_logs+1;n<=new_number_logs;n++) LG[n]=log((Double)n);
    number_logs=new_number_logs;

    if(my_verbose>0)
    cout << "#    extended log table to: " << number_logs << endl;

    delete [] tmp_LG;

}

void extend_sqrt_table(int m){


    int n;

    Double *tmp_sqrt; //used to copy over the old values
    tmp_sqrt = new Double[number_sqrts+1];
    for(n=1;n<=number_sqrts;n++) tmp_sqrt[n]=two_inverse_SQUARE_ROOT[n];

    delete [] two_inverse_SQUARE_ROOT;
    int new_number_sqrts=(int)(1.5*m); // extend table by 50 percent

    two_inverse_SQUARE_ROOT = new Double[new_number_sqrts+1];
    for(n=1;n<=number_sqrts;n++) two_inverse_SQUARE_ROOT[n]=tmp_sqrt[n];
    for(n=number_sqrts+1;n<=new_number_sqrts;n++) two_inverse_SQUARE_ROOT[n]=2/sqrt((Double)n);
    number_sqrts=new_number_sqrts;

    if(my_verbose>0)
    cout << "#    extended sqrt table to: " << number_sqrts << endl;

    delete [] tmp_sqrt;

}


// a lot of Q time is spend here. This can be precomputed XXXXX
Double dfac(int i)
{
    Double t=1;
    int n=i;
    if(i==1 || i==0) return 1;
    while(n>0) {
        t*=n;
        n-=2;
    }
    return t;
}

void extend_prime_table(int x)
{

  //if(x < prime_table[number_primes-1]) return; //if prime table is already big enough

  bool *sieve;
  int sqrtx;
  int n, p;
  int number_primes_guess;

  sieve = new bool[x+1];

  if (prime_table) delete [] prime_table;
  number_primes = 0;

  // Use explicit prime number theorem to allocate enough memory for primes
  number_primes_guess = max((int)ceil(x/(log((double)x)-4))+1, 100);

  prime_table = new int[number_primes_guess];

  if (my_verbose > 0)
  {
    cout << "#    extending prime table to: " << x << endl;
    //cout << "    guessed " << number_primes_guess << " primes; ";
  }
  for (n=0;n<=x;n++)
    sieve[n] = 1;

  p = 2;
  //sqrtx = (int) sqrt((double)x);
  sqrtx = Int(sqrt(Double(x)));
  while (p <= sqrtx)
  {
    n = 2*p;
    while (n <= x)
    {
      sieve[n] = 0;
      n = n + p;
    }
    do
      p++;
    while (sieve[p] == 0);
  }
  for (n=2;n<=x;n++)
  {
    if (sieve[n] == 1)
    {
      prime_table[number_primes] = n;
      number_primes++;
    }
  }
  delete [] sieve;

  if (my_verbose > 0) cout << "#    found " << number_primes << " primes." << endl;

  return;
}

int get_prime(int j)
{
  while (j >= number_primes)
  {
    extend_prime_table(prime_table[number_primes-1]*2);
  }

  return prime_table[j];
}


/*-------------------------------------------------------------------------------
Initialize the cosine Taylor series array
---------------------------------------------------------------------------------*/
void initialize_cos_array(){

    if(DIGITS<17)
        cos_taylor_arraysize=3000; //number of Taylor series to store
    else //else we should have a denser collection of series so as to reduce the number of terms needed
        cos_taylor_arraysize=100000; //number of Taylor series to store

    number_cos_taylor_terms=0; //how many Taylor coefficients to compute

    double r=1.;
    int n=0;
    if(DIGITS<17) number_cos_taylor_terms=4;
    else{
        do{
            n+=2;
            r*=.25e-10/((n-1)*n); // r is (1/2 * 1/100000)^n/n! to reflect max size of nth Taylor coeff
        }while(r>tolerance);
        number_cos_taylor_terms=n;
    }

    if(my_verbose>0)
        cout << "#    Will compute " << number_cos_taylor_terms << " terms of each cos taylor series." << endl;

    one_over_cos_taylor_arraysize=Double(1)/cos_taylor_arraysize;
    twoPi_over_cos_taylor_arraysize=2*Pi/cos_taylor_arraysize;

    if(cos_taylor) delete [] cos_taylor;
    cos_taylor = new Double[(cos_taylor_arraysize+1)*number_cos_taylor_terms];


    Double x,f;
    int count=0;
    for(int m = 0; m <= cos_taylor_arraysize; m++)
    {
        x=2*Pi*(m+.5)/cos_taylor_arraysize;

        Double cos_x=cos(x), sin_x=sin(x);
        f=1.;
        for(int j=0;j!=number_cos_taylor_terms;j=j+2){
            cos_taylor[count]=cos_x*f;
            f/=(j+1); //Bailey's package doesn't allow for r/=Double(j+1), and r/=(j+1) casts j+1 as a double
            count++;
            cos_taylor[count]=-sin_x*f;
            f/=(j+2);
            count++;
            f=-f;
        }
    }


}

/*-----------------------------------------------------------------------------------------
initialize_rs_remainder initializes the taylor series about z=0, for the correction terms in
the Riemann Siegel formula. The first few correction terms, with P[j]=Psi^(j)(2(z+1/2))
P[0]=cos(Pi*(u^2/2-u-1/8))/cos(Pi*u):

These are the coefficients of z^(2n) or z^(2n+1) depending on whether we are looking at
rs_remainder[even] or rs_remainder[odd] respectively.

The number, 71, of coefficients given is what is needed to compute about 64
digits precision of all the correction terms listed. This was determined taking
into account the size of the coefficients and the powers of |z| <= 1/2.

------------------------------------------------------------------------------------------*/


void initialize_rs_remainder1(){

    //============================= n = 0 ================================================
    rs_remainder[0][0]=str_to_Double(".3826834323650897717284599840303988667613445624856270414338006356");
    rs_remainder[0][1]=str_to_Double("1.748961872310081797441185869485327948292166004169452127465514443");
    rs_remainder[0][2]=str_to_Double("2.118025207685496373184564278264168887589167928938778762639083009");
    rs_remainder[0][3]=str_to_Double("-.8707216670511480739189240773823940902852396975500125749287602886");
    rs_remainder[0][4]=str_to_Double("-3.473311224346516707306411669375796783623576822241419998509642672");
    rs_remainder[0][5]=str_to_Double("-1.662694730899932449643136301192857319075022559408689464728555836");
    rs_remainder[0][6]=str_to_Double("1.216731288919232134476893528044169908502254933813501815960480030");
    rs_remainder[0][7]=str_to_Double("1.301430416100797577300605380997857517360085651554830737525560034");
    rs_remainder[0][8]=str_to_Double(".3051102182736167242108987123980541829057180919498709488873851053e-1");
    rs_remainder[0][9]=str_to_Double("-.3755803051545095242798193212293423794084270648139628309988198050");
    rs_remainder[0][10]=str_to_Double("-.1085784416564065974354697590132864187244451846829672423172812207");
    rs_remainder[0][11]=str_to_Double(".5183290299954962337576051067322495156659652970039163191501317935e-1");
    rs_remainder[0][12]=str_to_Double(".2999948061990227592040084956911813060279448783088444499485510180e-1");
    rs_remainder[0][13]=str_to_Double("-.2275939670612564226019948510205498195977410217168016633552534643e-2");
    rs_remainder[0][14]=str_to_Double("-.4382647416580338305940070135851452276667106360520495258064816116e-2");
    rs_remainder[0][15]=str_to_Double("-.4064230183729846993072327211603183384158540507475599675397975319e-3");
    rs_remainder[0][16]=str_to_Double(".4006097785422113927891031460770685785807272493435465423563616889e-3");
    rs_remainder[0][17]=str_to_Double(".8971057991388841297834181953786891780312293618586957566459432886e-4");
    rs_remainder[0][18]=str_to_Double("-.2302565002723910711610294525738991854154840094916001285826516636e-4");
    rs_remainder[0][19]=str_to_Double("-.9380006601906792484719729401274735908193460377157883523565550115e-5");
    rs_remainder[0][20]=str_to_Double(".6323514947609107504249861239594295018383639220840328596632428537e-6");
    rs_remainder[0][21]=str_to_Double(".6551022819231501666212231231334105099853135960394307090094723019e-6");
    rs_remainder[0][22]=str_to_Double(".2210523745552697258660868903818763926839469833511570635200483589e-7");
    rs_remainder[0][23]=str_to_Double("-.3322316176445628835031335170176244552138467244689322365188454675e-7");
    rs_remainder[0][24]=str_to_Double("-.3734910989933656081764604760152216088162375339767051103604655047e-8");
    rs_remainder[0][25]=str_to_Double(".1244506706079773919515100002493662033911835955288378516911042123e-8");
    rs_remainder[0][26]=str_to_Double(".2476820537650219184251435483288424602185504196320793903112029936e-9");
    rs_remainder[0][27]=str_to_Double("-.3284272816891627194458969154064130192953139674775330616686265696e-10");
    rs_remainder[0][28]=str_to_Double("-.1130540685229840367788134533541012727211926225115744806746082945e-10");
    rs_remainder[0][29]=str_to_Double(".4565463979588693927593011016423466252865616217039763529476023304e-12");
    rs_remainder[0][30]=str_to_Double(".3959848094524921519585471190127944112380583726017443959325651733e-12");
    rs_remainder[0][31]=str_to_Double(".7849566221259617317141828461287821447964091640589424933731432712e-14");
    rs_remainder[0][32]=str_to_Double("-.1105904315099123319372344571747309322893399621754990551092223760e-13");
    rs_remainder[0][33]=str_to_Double("-.7738543987641508317060001628750396777823417852899949946095867167e-15");
    rs_remainder[0][34]=str_to_Double(".2485775555027137218484004673471218482480566334404120628760275479e-15");
    rs_remainder[0][35]=str_to_Double(".3051479718882721790971875899825612649459073783286182795478366016e-16");
    rs_remainder[0][36]=str_to_Double("-.4414297887793302845228285032740064416390158138269782203523927546e-17");
    rs_remainder[0][37]=str_to_Double("-.8631388878188414739303188085218835728845711510658488892976082562e-18");
    rs_remainder[0][38]=str_to_Double(".5701292196842975217572887724477185708269470336908862987409860339e-19");
    rs_remainder[0][39]=str_to_Double(".1952964016419934107653433477234329708725454373631441870349199857e-19");
    rs_remainder[0][40]=str_to_Double("-.3370766713534960218133595818688514788624085206762193356358473641e-21");
    rs_remainder[0][41]=str_to_Double("-.3679459871576221269046597696095914442234506560731310625782020193e-21");
    rs_remainder[0][42]=str_to_Double("-.7311865182444788001771694062355879145998697877668104758987711862e-23");
    rs_remainder[0][43]=str_to_Double(".5869094638676538817500150747027581719946076649786144678312984633e-23");
    rs_remainder[0][44]=str_to_Double(".3130759211365692456865742857348097669373380175940887327504255828e-24");
    rs_remainder[0][45]=str_to_Double("-.7947839566038058637226172211111547660815460372684777700543035068e-25");
    rs_remainder[0][46]=str_to_Double("-.7084209480928370484910261300898518950323663886857773160562843635e-26");
    rs_remainder[0][47]=str_to_Double(".9022840660124028291282305834073381848912841389857162152974282452e-27");
    rs_remainder[0][48]=str_to_Double(".1224410947927371892746059097886467990567919755542546652398447088e-27");
    rs_remainder[0][49]=str_to_Double("-.8218858109260206931969615869478086402389595136113365075747891492e-29");
    rs_remainder[0][50]=str_to_Double("-.1761392489687557943153965229854080331724879786596025787578272129e-29");
    rs_remainder[0][51]=str_to_Double(".5136687660743046661443512914349478774464151046287771087609017091e-31");
    rs_remainder[0][52]=str_to_Double(".2181932145256603096324670258546611035588273976666275085740502657e-31");
    rs_remainder[0][53]=str_to_Double("-.1658478523782396526460223467087688542025568822467604169920921936e-34");
    rs_remainder[0][54]=str_to_Double("-.2364464947075049478311070992062523165713863688449818205811482134e-33");
    rs_remainder[0][55]=str_to_Double("-.5570484477199485894027672593223193030649037672712930763891463003e-35");
    rs_remainder[0][56]=str_to_Double(".2255485933729007125784223110408832656442480664948692844803966492e-35");
    rs_remainder[0][57]=str_to_Double(".1064235399184662951041323493934114021130851545087465169966826535e-36");
    rs_remainder[0][58]=str_to_Double("-.1891121787521712141529437286759380472641744838115078873221579860e-37");
    rs_remainder[0][59]=str_to_Double("-.1388297220262483580247299617028760235625406043768336009408265213e-38");
    rs_remainder[0][60]=str_to_Double(".1377546581031590690632065463233382714214703858124249070516350230e-39");
    rs_remainder[0][61]=str_to_Double(".1479755095507785485137406071803213915509004092109922800605576699e-40");
    rs_remainder[0][62]=str_to_Double("-.8440239487797857712758312893209576176610376448060127647327691299e-42");
    rs_remainder[0][63]=str_to_Double("-.1363029128165632378718801657715797269829339558167659601250267910e-42");
    rs_remainder[0][64]=str_to_Double(".3952185926800835118526308802240734693344257085185263761026165167e-44");
    rs_remainder[0][65]=str_to_Double(".1112181856478376656704571061597801339345109958020899315847926019e-44");
    rs_remainder[0][66]=str_to_Double("-.8441003145417726264431880594984599475584364419001824366226682452e-47");
    rs_remainder[0][67]=str_to_Double("-.8139403608089633378545982760949371764318728704639810807953186667e-47");
    rs_remainder[0][68]=str_to_Double("-.8658059439497875938975979419498216110785486360122997778295175466e-49");
    rs_remainder[0][69]=str_to_Double(".5373865326463434022149365916316237939767249405521302032506388331e-49");
    rs_remainder[0][70]=str_to_Double(".1492965137935395555494865149269851797895335590512157309217780204e-50");
    rs_remainder[0][71]=str_to_Double("-.3203895770351591499389099698930529595901436713541334026177156438e-51");
    //============================= n = 1 ================================================
    rs_remainder[1][0]=str_to_Double("-.5365020525675069405998280791133349931854094486128644149955350563e-1");
    rs_remainder[1][1]=str_to_Double(".1102781874108148243989636207191692989275858047802601063557263910");
    rs_remainder[1][2]=str_to_Double("1.231720015431522631319565291622059740675161839824290073064915250");
    rs_remainder[1][3]=str_to_Double("1.263496486279945788417554821912133208678939468373697075193327193");
    rs_remainder[1][4]=str_to_Double("-1.695108997559503018449447399065957275034671030332175857033976893");
    rs_remainder[1][5]=str_to_Double("-2.999871196765010088955487358941413331102813741642153765571236484");
    rs_remainder[1][6]=str_to_Double("-.1081994495989920864269225778743784884165987156472166204500396911");
    rs_remainder[1][7]=str_to_Double("1.940766294621271268793876325397157775541178365431267499318739991");
    rs_remainder[1][8]=str_to_Double(".7838423561500686532884345748869425947391937969131127571308010394");
    rs_remainder[1][9]=str_to_Double("-.5054829667900365918790214132617270722931110255931258455021902826");
    rs_remainder[1][10]=str_to_Double("-.3845072349605797405134273885311015390964058075107949074977645536");
    rs_remainder[1][11]=str_to_Double(".3747264646531532067594447494022893716348312493754551415921687630e-1");
    rs_remainder[1][12]=str_to_Double(".9092026610973176317258142450575652865159861194046748976303154858e-1");
    rs_remainder[1][13]=str_to_Double(".1044923755006450921820113972658784258869072792274009009949183093e-1");
    rs_remainder[1][14]=str_to_Double("-.1258297965158341649747892224592184374339262437755601621304228301e-1");
    rs_remainder[1][15]=str_to_Double("-.3399503721151274085058948861372571436810272995565060041851080837e-2");
    rs_remainder[1][16]=str_to_Double(".1041095053771489126829542406554543241783057979522741680508448288e-2");
    rs_remainder[1][17]=str_to_Double(".5010949051118486860355652672743046866628161804683724486046036172e-3");
    rs_remainder[1][18]=str_to_Double("-.3956359669003181559547118556963366701092080477448761396224467605e-4");
    rs_remainder[1][19]=str_to_Double("-.4762459245357189638654098302680350067152074032296393922787948549e-4");
    rs_remainder[1][20]=str_to_Double("-.1853935533808513227343490645691168880364883755693923107091910294e-5");
    rs_remainder[1][21]=str_to_Double(".3193691808006897204046635393432675467575125346290484899516982711e-5");
    rs_remainder[1][22]=str_to_Double(".4090780760850606632650894536770176800340649760323484507745940534e-6");
    rs_remainder[1][23]=str_to_Double("-.1544662433257663212843757232731043276970446466653502900257529683e-6");
    rs_remainder[1][24]=str_to_Double("-.3466307491769133172225594059340728816492418636839426453719970358e-7");
    rs_remainder[1][25]=str_to_Double(".5158711258806154784708610822160350430917870446344875918171846119e-8");
    rs_remainder[1][26]=str_to_Double(".1984539255640794420228408155024203027231284992805072501990723482e-8");
    rs_remainder[1][27]=str_to_Double("-.8920820862551490847995104313036936836896951633200655796655741832e-10");
    rs_remainder[1][28]=str_to_Double("-.8581017807796222642200248697440284940668548770771906299066496921e-10");
    rs_remainder[1][29]=str_to_Double("-.1879955001383284952091505117756490946347641216713160133427567968e-11");
    rs_remainder[1][30]=str_to_Double(".2917821950594353786495376346779228616084995421569854809284642925e-11");
    rs_remainder[1][31]=str_to_Double(".2242464328378943531236480385491720050653277006398917057622406871e-12");
    rs_remainder[1][32]=str_to_Double("-.7888938771825912503789827138301722977156873722942032659255932917e-13");
    rs_remainder[1][33]=str_to_Double("-.1057780490885248615717739010307786089269476961111912877648531614e-13");
    rs_remainder[1][34]=str_to_Double(".1667168683572909141230354845174074085971721441403969602841325051e-14");
    rs_remainder[1][35]=str_to_Double(".3543209091148631748047711103433592916545224534161707530197544734e-15");
    rs_remainder[1][36]=str_to_Double("-.2538100978709331063520804208206192515208887478353350424917034514e-16");
    rs_remainder[1][37]=str_to_Double("-.9408538863065022553744148351305851607660085456506536020786888477e-17");
    rs_remainder[1][38]=str_to_Double(".1753756925869445233040896031340920471304866563584622301354242462e-18");
    rs_remainder[1][39]=str_to_Double(".2063488014466568312581210570838094226725448336467330832028832412e-18");
    rs_remainder[1][40]=str_to_Double(".4411928113648420582213991918316532303532740517285137567286017006e-20");
    rs_remainder[1][41]=str_to_Double("-.3803617151918602561643312674252832035308052762428172622215378654e-20");
    rs_remainder[1][42]=str_to_Double("-.2175600072560313868834821329762752228408490524556768784576479505e-21");
    rs_remainder[1][43]=str_to_Double(".5912801531051567197011905610917804959701945749434592634214584441e-22");
    rs_remainder[1][44]=str_to_Double(".5633679619144557647576236873879203702403414424709177382321018482e-23");
    rs_remainder[1][45]=str_to_Double("-.7658980063274947684213182668909980152060512105040827125639682583e-24");
    rs_remainder[1][46]=str_to_Double("-.1107844784922138392642822981317749399831154713797705756179522942e-24");
    rs_remainder[1][47]=str_to_Double(".7916068568868281039969716721614354812753595587169014971076992205e-26");
    rs_remainder[1][48]=str_to_Double(".1803625771154524408702802176825639979668691351911176577122323939e-26");
    rs_remainder[1][49]=str_to_Double("-.5585135656831858814978042741646634162249813384595625260006295205e-28");
    rs_remainder[1][50]=str_to_Double("-.2516175897434856683865157065968131168568922965407905889755718071e-28");
    rs_remainder[1][51]=str_to_Double(".2026130327807167588430946014072365383755288789700774277465245651e-31");
    rs_remainder[1][52]=str_to_Double(".3056858247034789733103896874969721148886225639687825174988644399e-30");
    rs_remainder[1][53]=str_to_Double(".7613159498422358695237845780448877035242455759525044529283107112e-32");
    rs_remainder[1][54]=str_to_Double("-.3255388546517976033162161462212407484981270569219792719430305513e-32");
    rs_remainder[1][55]=str_to_Double("-.1620570912911252714423660013689806489443099166593337839450269654e-33");
    rs_remainder[1][56]=str_to_Double(".3035352968480857707385284761125837961782382114460614986345535726e-34");
    rs_remainder[1][57]=str_to_Double(".2346590639455370467393147654643735576565152042167545634303277797e-35");
    rs_remainder[1][58]=str_to_Double("-.2449884046100996033493915270679145550890953396661909210687545739e-36");
    rs_remainder[1][59]=str_to_Double("-.2766593248341706529251328063100758311245005317335025323794627335e-37");
    rs_remainder[1][60]=str_to_Double(".1657564429979149980937974751298444076540748877083999340629066863e-38");
    rs_remainder[1][61]=str_to_Double(".2809547647427386122974638213095790666998904870635342878811834154e-39");
    rs_remainder[1][62]=str_to_Double("-.8543791170081281630723442153713930031445605791340579281459304403e-41");
    rs_remainder[1][63]=str_to_Double("-.2519694336290890372652350190339048304195249421480125016524660413e-41");
    rs_remainder[1][64]=str_to_Double(".2002680470982924708928914679402041268341391538793405265352335231e-43");
    rs_remainder[1][65]=str_to_Double(".2020933283254223213319362691652012862599645936821392723483114570e-43");
    rs_remainder[1][66]=str_to_Double(".2248159189250192233975359010272303659519706751317436616437074237e-45");
    rs_remainder[1][67]=str_to_Double("-.1458325137686845832992211208208641080493949936778943554558195129e-45");
    rs_remainder[1][68]=str_to_Double("-.4231561046429640475100709951304927768913015307828674480691313910e-47");
    rs_remainder[1][69]=str_to_Double(".9478571564046474716378087655414378016919497792612989906882479686e-48");
    rs_remainder[1][70]=str_to_Double(".4510975164437871041082476036614058789725916066100683815297094031e-49");
    rs_remainder[1][71]=str_to_Double("-.5526148131408408156854272262032530507538460210068722893031550465e-50");
    //============================= n = 2 ================================================
    rs_remainder[2][0]=str_to_Double(".5188542830293168493784581519230959565968684337910516563725522453e-2");
    rs_remainder[2][1]=str_to_Double(".123786335522538984133826974438351529467800123179513757157434873e-2");
    rs_remainder[2][2]=str_to_Double("-.1813750572516699741149189640941362145639915982441846550217386872");
    rs_remainder[2][3]=str_to_Double(".1429149274853212654116560337651435620542974718825139004882758770");
    rs_remainder[2][4]=str_to_Double("1.330339176668756532509933299985457123460053019550360441500129172");
    rs_remainder[2][5]=str_to_Double(".3522472353403733677532765550583613028780635138198853104043319929");
    rs_remainder[2][6]=str_to_Double("-2.421001595891950723781530543340501105290598649926956109950809918");
    rs_remainder[2][7]=str_to_Double("-1.676078702253810885333461814923720644025343553285547698599466478");
    rs_remainder[2][8]=str_to_Double("1.368941672332837218423491538070760058466225708102037768187777931");
    rs_remainder[2][9]=str_to_Double("1.553901943022298322145639526559347726284816347939935962719058333");
    rs_remainder[2][10]=str_to_Double("-.1722164273472998051958258699891772541754273973596261998684118949");
    rs_remainder[2][11]=str_to_Double("-.6359068055045430988970490235584462364640473159972306193885511651");
    rs_remainder[2][12]=str_to_Double("-.9911649873041208105423564341369990594844217098597685233543465367e-1");
    rs_remainder[2][13]=str_to_Double(".1403348006738700895073825489831586845516559897399154564334105620");
    rs_remainder[2][14]=str_to_Double(".4782352019827292236438803506512481834959874758665861603701292800e-1");
    rs_remainder[2][15]=str_to_Double("-.1735604064147978079795864709223282676958084522417238134753381982e-1");
    rs_remainder[2][16]=str_to_Double("-.1022501253402859184447660413126300828489699854485411724418751122e-1");
    rs_remainder[2][17]=str_to_Double(".9274149159794887899427001437137195009241932389235423857716463922e-3");
    rs_remainder[2][18]=str_to_Double(".1357219437237338534525336199577230127124914915844177600486918535e-2");
    rs_remainder[2][19]=str_to_Double(".6413690120293880089962387363945326797172680505533411309712024540e-4");
    rs_remainder[2][20]=str_to_Double("-.1230080569819662988334232293655311402280396551206588402380682124e-3");
    rs_remainder[2][21]=str_to_Double("-.1831350740478920255476755439796205141772519046833586489637629776e-4");
    rs_remainder[2][22]=str_to_Double(".7821628604322627308501399384618720544215010952141294955072844881e-5");
    rs_remainder[2][23]=str_to_Double(".2008754248475994550349852939191568214658932787821910191580723550e-5");
    rs_remainder[2][24]=str_to_Double("-.3353276539318571373727497272414531442842343783888858788470755517e-6");
    rs_remainder[2][25]=str_to_Double("-.1461602091741823092645100971227603722249074918815445596076851204e-6");
    rs_remainder[2][26]=str_to_Double(".7261497384040072462492982995434329661188846726079154435292332630e-8");
    rs_remainder[2][27]=str_to_Double(".7894805679006706236084907605818256889581374370239231168349805093e-8");
    rs_remainder[2][28]=str_to_Double(".1958985823464410454317580804192010254335541455617833977089122370e-9");
    rs_remainder[2][29]=str_to_Double("-.3302802050431173020571795860805970713803754608013272113733947522e-9");
    rs_remainder[2][30]=str_to_Double("-.2814894537376278741861907725299490592433181536751411318242799992e-10");
    rs_remainder[2][31]=str_to_Double(".1084007931514484249791451039697496807389120518929504436689161628e-10");
    rs_remainder[2][32]=str_to_Double(".1599196020009304118050184603894982648958542823320220320854446259e-11");
    rs_remainder[2][33]=str_to_Double("-.2747810753378103212034303107642090144553821479447927143862905424e-12");
    rs_remainder[2][34]=str_to_Double("-.6388781373997463403763185062550446407540785966976864145644069464e-13");
    rs_remainder[2][35]=str_to_Double(".4962566999974764002256006236544031318448644875810288924118613802e-14");
    rs_remainder[2][36]=str_to_Double(".2004197125300311226023004552980093897519480495592919131886378507e-14");
    rs_remainder[2][37]=str_to_Double("-.4010496701027515651656466633168224001349500913393396369974124636e-16");
    rs_remainder[2][38]=str_to_Double("-.5148132304571040165881993125620101574142275952653379800450895118e-16");
    rs_remainder[2][39]=str_to_Double("-.1195244426717772932380474043773942073755042240935339474801243240e-17");
    rs_remainder[2][40]=str_to_Double(".1102685403966278855251475971021765714464277319778563363604454228e-17");
    rs_remainder[2][41]=str_to_Double(".6795774139729971823155656406959849234114039699054070704521200177e-19");
    rs_remainder[2][42]=str_to_Double("-.1977606633239386282831932356947581611393910149068921970832895106e-19");
    rs_remainder[2][43]=str_to_Double("-.2021942546058184397562097142109100000175443881208631806826397018e-20");
    rs_remainder[2][44]=str_to_Double(".2935996867116108350510423984728159196690561941213445122947831671e-21");
    rs_remainder[2][45]=str_to_Double(".4542483274947033670191671390364098035512820956915378066795739206e-22");
    rs_remainder[2][46]=str_to_Double("-.3456821516533061096436395196183796559342151287556858152344723988e-23");
    rs_remainder[2][47]=str_to_Double("-.8402320445256338890443881022902716454143124568705212451743524883e-24");
    rs_remainder[2][48]=str_to_Double(".2761371188047659029030127894443544642362570394453722483025445564e-25");
    rs_remainder[2][49]=str_to_Double(".1324929066999056949387818109040867356635638334410877081741693355e-25");
    rs_remainder[2][50]=str_to_Double("-.1050776970686853300589259902992068638509066800794046016781908778e-28");
    rs_remainder[2][51]=str_to_Double("-.1810704553381788899211119574092082008532411483027973272559175585e-27");
    rs_remainder[2][52]=str_to_Double("-.4784862049560126534665946219009029687872785439940885402400357990e-29");
    rs_remainder[2][53]=str_to_Double(".2159567786232025343287959436916685939853223031377111517018751164e-29");
    rs_remainder[2][54]=str_to_Double(".1137118102421619321357862590028576417251691904426181063216769491e-30");
    rs_remainder[2][55]=str_to_Double("-.2245794777288454820977450654627683907949212710601628164129762297e-31");
    rs_remainder[2][56]=str_to_Double("-.1832328658704393547248690414893160303227098906592775457159827026e-32");
    rs_remainder[2][57]=str_to_Double(".2013843883809104262503120530869166472813491815484739636664920461e-33");
    rs_remainder[2][58]=str_to_Double(".2395573014978892315609041832200335164480647339653139685119975405e-34");
    rs_remainder[2][59]=str_to_Double("-.1508298284205927013320407993746939427184954484276036427002021898e-35");
    rs_remainder[2][60]=str_to_Double("-.2688607987932643399427419600341578464704205939312913642561813354e-36");
    rs_remainder[2][61]=str_to_Double(".8575423237585002118417249924366817742084135102941323225353421964e-38");
    rs_remainder[2][62]=str_to_Double(".2656292237730002720203625456519774633665949699637172575868362246e-38");
    rs_remainder[2][63]=str_to_Double("-.2206055045036210101995327310978220837119660139520293720981361363e-40");
    rs_remainder[2][64]=str_to_Double("-.2339955909317179683558085048899964100985021577655558600722811781e-40");
    rs_remainder[2][65]=str_to_Double("-.2730804473883891779168319473394591170390777018622464639453357831e-42");
    rs_remainder[2][66]=str_to_Double(".1849288237007360000607543435828474589399049800986458495722618335e-42");
    rs_remainder[2][67]=str_to_Double(".5614334108441998543271759911147538331813936849172293356285499999e-44");
    rs_remainder[2][68]=str_to_Double("-.1312878788426408099877922616136057508109694463369659632744544384e-44");
    rs_remainder[2][69]=str_to_Double("-.6527178400453218926673118820008423419455138825602322472446911710e-46");
    rs_remainder[2][70]=str_to_Double(".8339432317534500418546503973478976528809634863496996077338294016e-47");
    rs_remainder[2][71]=str_to_Double(".6016652956692163453127980444072369443837961742310867160119377791e-48");
    //============================= n = 3 ================================================
    rs_remainder[3][0]=str_to_Double("-.2679432181438913808539671459890456247712707906335677313857985118e-2");
    rs_remainder[3][1]=str_to_Double(".2995372109103514963731329491569917265027452033956800816390676319e-1");
    rs_remainder[3][2]=str_to_Double("-.425701725418286979850193511168771363167068238751826558828839419e-1");
    rs_remainder[3][3]=str_to_Double("-.2899796577980388750689320947866888815887297721455598771928237875");
    rs_remainder[3][4]=str_to_Double(".4888831999235445972537474640716858037813574483608842133316313471");
    rs_remainder[3][5]=str_to_Double("1.230855876395746081193125043362941191563787449863696246736851778");
    rs_remainder[3][6]=str_to_Double("-.829756070852740870417969104329757351019575809355355388794472459");
    rs_remainder[3][7]=str_to_Double("-2.249763536666566866520450126599034809328793200254604172941573256");
    rs_remainder[3][8]=str_to_Double(".784513996100547137936547362018357998522630710168227301340620972e-1");
    rs_remainder[3][9]=str_to_Double("1.746749280086889400391986666452186788901871710076338108110046796");
    rs_remainder[3][10]=str_to_Double(".4596808097974993510923730617316921983800453595353511641167104562");
    rs_remainder[3][11]=str_to_Double("-.6619353471039774946433904000898284322110681939392204256589640137");
    rs_remainder[3][12]=str_to_Double("-.3159044103617363457897963297331576998293888302893566020648435956");
    rs_remainder[3][13]=str_to_Double(".1284479254520749598851184747620936581235148743065505409736145651");
    rs_remainder[3][14]=str_to_Double(".1007338271662615230096945020751295986900045444657351483326066935");
    rs_remainder[3][15]=str_to_Double("-.9530183848825267759504659842297419233227549892358832451269198468e-2");
    rs_remainder[3][16]=str_to_Double("-.1926442168751408889840098069714443332236390669436686767610994963e-1");
    rs_remainder[3][17]=str_to_Double("-.1246463715876929171247907164583021249997061086328213380941630655e-2");
    rs_remainder[3][18]=str_to_Double(".2424396964110308573972152458406465097941364324670569449815574393e-2");
    rs_remainder[3][19]=str_to_Double(".4376476977418570182756129039558366912216413635432052537203211066e-3");
    rs_remainder[3][20]=str_to_Double("-.2071403268700179127591307830406688286287589059508390609497533219e-3");
    rs_remainder[3][21]=str_to_Double("-.6274344504186515560526109580298036054689301090914734216007345058e-4");
    rs_remainder[3][22]=str_to_Double(".1157534381459566934837892089893160300949712084669639203604928599e-4");
    rs_remainder[3][23]=str_to_Double(".5883854924540379783885978856970775611374486310354052173336527910e-5");
    rs_remainder[3][24]=str_to_Double("-.3124677400696336220869614490760334642292560368465762425780391793e-6");
    rs_remainder[3][25]=str_to_Double("-.4024065775498959500979814931371302363202232653378608595627599258e-6");
    rs_remainder[3][26]=str_to_Double("-.1199110779489632960574695390403727489410829824397694493877994125e-7");
    rs_remainder[3][27]=str_to_Double(".2096375406393870831831451427047264343353148469652996999242324695e-7");
    rs_remainder[3][28]=str_to_Double(".2020356022540215377863668667015120281871034580343968610739441865e-8");
    rs_remainder[3][29]=str_to_Double("-.8440146463909390057210900639795926834173393338929557655183172241e-9");
    rs_remainder[3][30]=str_to_Double("-.1388884542004012860683518982881711564708584076229068752872807913e-9");
    rs_remainder[3][31]=str_to_Double(".2588490692171973475445882063278485048500790529560507539281713470e-10");
    rs_remainder[3][32]=str_to_Double(".6664830790556666110923464935176428511664371023463559874884268671e-11");
    rs_remainder[3][33]=str_to_Double("-.5577569833891270686294858594819087420399727426154212273070758926e-12");
    rs_remainder[3][34]=str_to_Double("-.2487835961168490237006826345380064416286838376423413301025166521e-12");
    rs_remainder[3][35]=str_to_Double(".5220722370970689138374683797285494931069896548498206606119202476e-14");
    rs_remainder[3][36]=str_to_Double(".7534873081244042107799095498887804739983358688288527655859170939e-14");
    rs_remainder[3][37]=str_to_Double(".1939902746711014610211859578370698106882691237827836293922935896e-15");
    rs_remainder[3][38]=str_to_Double("-.1886892371609056930616238509868980790178541237920031538180917032e-15");
    rs_remainder[3][39]=str_to_Double("-.1264648274431615753029201218976088904011166245496494288932336626e-16");
    rs_remainder[3][40]=str_to_Double(".3925715193537363628220821482952320067339667591007380668913785082e-17");
    rs_remainder[3][41]=str_to_Double(".4336765160753024127108265354322415839247262303028638483424728802e-18");
    rs_remainder[3][42]=str_to_Double("-.6712223979407535352141207908587198686668979393504288427028262169e-19");
    rs_remainder[3][43]=str_to_Double("-.1117617655812504315047214564206255534610356247643023224294215900e-19");
    rs_remainder[3][44]=str_to_Double(".9037721692447813963650827747877129994039250836834884689687076031e-21");
    rs_remainder[3][45]=str_to_Double(".2358420992377783032537762433314377854533566287188115872951590411e-21");
    rs_remainder[3][46]=str_to_Double("-.8186656220775764697454458642055239169949831422164076937650930018e-23");
    rs_remainder[3][47]=str_to_Double("-.4220063889328535987081346884658665420324043508898444680071182678e-23");
    rs_remainder[3][48]=str_to_Double(".2501227571969917714633802750263392734925558369967503088026016526e-26");
    rs_remainder[3][49]=str_to_Double(".6511680744527076486987537821863833518130122028463335858236301682e-25");
    rs_remainder[3][50]=str_to_Double(".1839698886935773602206369768948396936599065361854325991336374616e-26");
    rs_remainder[3][51]=str_to_Double("-.8727433864941407738169935208881192974961179188650678314794713790e-27");
    rs_remainder[3][52]=str_to_Double("-.4883030058951234812293308360846894772669257165056428433685290870e-28");
    rs_remainder[3][53]=str_to_Double(".1015425520269571153620606882556981522349940250488063674739205764e-28");
    rs_remainder[3][54]=str_to_Double(".8775787114724583039731962267514419003634008867102590656442188584e-30");
    rs_remainder[3][55]=str_to_Double("-.1014498160995071485512863055846453800793178695092950706466086388e-30");
    rs_remainder[3][56]=str_to_Double("-.1275507652142258453153002777874888608514841062385783423520631668e-31");
    rs_remainder[3][57]=str_to_Double(".8431439269975897906157770840742867219604818237583823155367037607e-33");
    rs_remainder[3][58]=str_to_Double(".1586054976031077851323966817442329012155794660479230997436087377e-33");
    rs_remainder[3][59]=str_to_Double("-.5296049637341945715943131637451395284563116190687810920698634707e-35");
    rs_remainder[3][60]=str_to_Double("-.1730430977261074878244986997062617712712268554785579585557577753e-35");
    rs_remainder[3][61]=str_to_Double(".1489206163973747041405588065097828146432842403158232108860654901e-37");
    rs_remainder[3][62]=str_to_Double(".1678080730918902715403424382115207676854118647780847713827591643e-37");
    rs_remainder[3][63]=str_to_Double(".2068873658586969031800836247160167061808775168131963776618492703e-39");
    rs_remainder[3][64]=str_to_Double("-.1455617799521557674736064597911522755681581178284852093372616385e-39");
    rs_remainder[3][65]=str_to_Double("-.4638695324787393398732594857508508447190776145537142896146911600e-41");
    rs_remainder[3][66]=str_to_Double(".1131055616938559978550863517411184968371121790194158523632492210e-41");
    rs_remainder[3][67]=str_to_Double(".5888537634624262251619790160814897972895651314979256396035594213e-43");
    rs_remainder[3][68]=str_to_Double("-.7842462198145411808357168734989140926769179678700024298106463372e-44");
    rs_remainder[3][69]=str_to_Double("-.5915472795540936257172733836918002350041589695351662384183402727e-45");
    rs_remainder[3][70]=str_to_Double(".4794250178748227387774200622977038527803734156716932155752715805e-46");
    rs_remainder[3][71]=str_to_Double(".5099225830100147988637500242725074627077812771865756163972601046e-47");
    //============================= n = 4 ================================================
    rs_remainder[4][0]=str_to_Double(".4648338936176338185363046255956724354485860691074530376452807266e-3");
    rs_remainder[4][1]=str_to_Double("-.4022642946136188303911539891451814630430097902705504598541582648e-2");
    rs_remainder[4][2]=str_to_Double(".3847177051796126883591306852717719532474217625903384345558459671e-2");
    rs_remainder[4][3]=str_to_Double(".6581175135809486002088309200741040383272395460605692178095043786e-1");
    rs_remainder[4][4]=str_to_Double("-.1960412434369444911769552844820478128181618815936251459890780685");
    rs_remainder[4][5]=str_to_Double("-.2085405368635885324440001279449438203234580420195613692109034102");
    rs_remainder[4][6]=str_to_Double(".9507754185141750945847757415105821521830035725416265875592100136");
    rs_remainder[4][7]=str_to_Double(".5341535312914873976051759245989361225324868753185661228592636984");
    rs_remainder[4][8]=str_to_Double("-1.676349441176340079591164482034037707003135957209431766401918656");
    rs_remainder[4][9]=str_to_Double("-1.076747157875128992787846635104317574103797674675059661871399845");
    rs_remainder[4][10]=str_to_Double("1.235339301656596985287883611892512115661353117755840936062685838");
    rs_remainder[4][11]=str_to_Double("1.025782534005727577183489495779140805889096082787530259516143222");
    rs_remainder[4][12]=str_to_Double("-.4012409579398854437872813752331312237303908069928824172891920426");
    rs_remainder[4][13]=str_to_Double("-.5036663995108303447958525759160364872472544843477390579371587341");
    rs_remainder[4][14]=str_to_Double(".3573487795502744985807080163387400744158012297415665630718782597e-1");
    rs_remainder[4][15]=str_to_Double(".1443176308678541662428523949584380865967285552611645757553758834");
    rs_remainder[4][16]=str_to_Double(".1509152741790346941712677290432031023801463213357815713096432652e-1");
    rs_remainder[4][17]=str_to_Double("-.2609887477919436131761773965447995033486781180826275762756856902e-1");
    rs_remainder[4][18]=str_to_Double("-.6126628379519261749049099089484084884310358690602785827170594982e-2");
    rs_remainder[4][19]=str_to_Double(".3077503129870841184767877821668240661402108730901688788895006903e-2");
    rs_remainder[4][20]=str_to_Double(".1156247893408875231612012042200893670513948398270706378729907477e-2");
    rs_remainder[4][21]=str_to_Double("-.2277596675847212747280773395343114508470543461686500838632842945e-3");
    rs_remainder[4][22]=str_to_Double("-.1418963711818144443268157989356342335647556491822723224852159341e-3");
    rs_remainder[4][23]=str_to_Double(".7464860307955919453122409844503131688218383899984026844265286205e-5");
    rs_remainder[4][24]=str_to_Double(".1247970164540911661744499888468706566148553134996569362559610657e-4");
    rs_remainder[4][25]=str_to_Double(".4863945184002094619079980847461804462152503463470913417839327607e-6");
    rs_remainder[4][26]=str_to_Double("-.8210237414123167233865939358809409002494144133157268831160227734e-6");
    rs_remainder[4][27]=str_to_Double("-.9223258397495269288607204081584781580088836148996588405093703994e-7");
    rs_remainder[4][28]=str_to_Double(".4103687848816232529632371802825885010876013394861621696346476127e-7");
    rs_remainder[4][29]=str_to_Double(".7693057400376143855280143293567183647516801467160054844879027788e-8");
    rs_remainder[4][30]=str_to_Double("-.1537147568328222211023312384545635374025853250104522752222525986e-8");
    rs_remainder[4][31]=str_to_Double("-.4466468429045027206377692204417846294428946015836642524345793882e-9");
    rs_remainder[4][32]=str_to_Double(".3970731966941122353746662947759771827234520837523736961598222843e-10");
    rs_remainder[4][33]=str_to_Double(".1999346572482580040523377510284601313961006790465337240344112574e-10");
    rs_remainder[4][34]=str_to_Double("-.4197063294194712427002733055999701516441602618430863591948920490e-12");
    rs_remainder[4][35]=str_to_Double("-.7193102535061939087703220212364117334514686364484405463891678520e-12");
    rs_remainder[4][36]=str_to_Double("-.2121104297127854704619101316316375903988272954006908784795614398e-13");
    rs_remainder[4][37]=str_to_Double(".2120381631691065169683165159581948945230839450709766300344189646e-13");
    rs_remainder[4][38]=str_to_Double(".1567294225005494689039459701767655838046966511740643478722432911e-14");
    rs_remainder[4][39]=str_to_Double("-.5148771682757868078568874743182918070245889225295216911833053471e-15");
    rs_remainder[4][40]=str_to_Double("-.6208386562023204241449806821228673631354938405607489728639572381e-16");
    rs_remainder[4][41]=str_to_Double(".1019131351286850584823594475163904872146334445405952092345899959e-16");
    rs_remainder[4][42]=str_to_Double(".1842732566741699181839487649250654664436733665340695320945622750e-17");
    rs_remainder[4][43]=str_to_Double("-.1575250047669901779960776554680229989420549140501067749914488024e-18");
    rs_remainder[4][44]=str_to_Double("-.4455454438251618493289259464345667218612347060577014223023427903e-19");
    rs_remainder[4][45]=str_to_Double(".1617972888601144027962903317301425648053406347896177216725123720e-20");
    rs_remainder[4][46]=str_to_Double(".9085049972738517360450768238831513163021037043583753199348248459e-21");
    rs_remainder[4][47]=str_to_Double("-.2851800665099397146486113131928392952523601281851695247798405378e-25");
    rs_remainder[4][48]=str_to_Double("-.1589095970500260126671710620260708584956202450500485498755221303e-22");
    rs_remainder[4][49]=str_to_Double("-.4855789087358807891560933071209457562407748428982016463125698495e-24");
    rs_remainder[4][50]=str_to_Double(".2402263458416748505984603301785530575831418052452900061579662036e-24");
    rs_remainder[4][51]=str_to_Double(".1437805467127560752982844046965402060742051189337373886190506696e-25");
    rs_remainder[4][52]=str_to_Double("-.3137622561968839653145487484454059528935259670532164736160030139e-26");
    rs_remainder[4][53]=str_to_Double("-.2887588652852895931257051337179553883152110351465917631659589840e-27");
    rs_remainder[4][54]=str_to_Double(".3503015805998071776805228674507699088521181649355472506036648849e-28");
    rs_remainder[4][55]=str_to_Double(".4677515322129049626375258623476206560179814504428448844379279949e-29");
    rs_remainder[4][56]=str_to_Double("-.3238363631479144654658530221050527606679526117047133523253420456e-30");
    rs_remainder[4][57]=str_to_Double("-.6461014168871971612254522478254635069746250645126544021373362603e-31");
    rs_remainder[4][58]=str_to_Double(".2249765334791500565630855810929371602377359543682190595575065198e-32");
    rs_remainder[4][59]=str_to_Double(".7804301349199836053705802001405756362276742205743117714200259861e-33");
    rs_remainder[4][60]=str_to_Double("-.6854861432990437371053587399702554447822290144970456696717610123e-35");
    rs_remainder[4][61]=str_to_Double("-.8351818373456204482619786725669344630477206852169412144677508211e-35");
    rs_remainder[4][62]=str_to_Double("-.1099571491050093679436022093430887251053795059376283307978729490e-36");
    rs_remainder[4][63]=str_to_Double(".7969964794108458245823274359014548055190841267135175304238848201e-37");
    rs_remainder[4][64]=str_to_Double(".2678706527046017138870079447482713712479258603265490753039529808e-38");
    rs_remainder[4][65]=str_to_Double("-.6792791458823540311434723624690444394515895304771475038284063566e-39");
    rs_remainder[4][66]=str_to_Double("-.3716070157811945691201916748225526866605165523575876283208400034e-40");
    rs_remainder[4][67]=str_to_Double(".5151577446385547086214044608773864871642223931055798200988725145e-41");
    rs_remainder[4][68]=str_to_Double(".4074456510803250864003273048589351415640800570763178679748138000e-42");
    rs_remainder[4][69]=str_to_Double("-.3435057604464650032419249507337626796208900410279563478463118123e-43");
    rs_remainder[4][70]=str_to_Double("-.3825594661195704322403965103641265552151095378666080003851076403e-44");
    rs_remainder[4][71]=str_to_Double(".1955435891480970278986002129400507819228727958483191314706583124e-45");
}

void initialize_rs_remainder2(){

    //============================= n = 5 ================================================
    rs_remainder[5][0]=str_to_Double(".2268681184573736317655795724547260210455363682623334900669596536e-3");
    rs_remainder[5][1]=str_to_Double(".110812468537183880897586725283788730363583255267509658128569738e-2");
    rs_remainder[5][2]=str_to_Double("-.1621857925555009106408484258686492831450236048360289119307015894e-1");
    rs_remainder[5][3]=str_to_Double(".527650340539874166272412666564860258836519964959117942525573335e-1");
    rs_remainder[5][4]=str_to_Double(".257088020090332399929001011109522381250583682397945016926838668e-1");
    rs_remainder[5][5]=str_to_Double("-.3805866044080639726443599184814563489241886009667079439382142314");
    rs_remainder[5][6]=str_to_Double(".2253198789264231532297692698983807286036841717261828395641525476");
    rs_remainder[5][7]=str_to_Double("1.034457331649522172113044996573887395688477415370301953152880234");
    rs_remainder[5][8]=str_to_Double("-.552825769705081378988884752967346186992602637953434753897178080");
    rs_remainder[5][9]=str_to_Double("-1.528771264107807299627365717141686780707417230431288465105495018");
    rs_remainder[5][10]=str_to_Double(".3282836642771958367203166939405941916743653763952277881627628156");
    rs_remainder[5][11]=str_to_Double("1.229110218540087062384250012396773759266820075643503402194989110");
    rs_remainder[5][12]=str_to_Double(".409369393831152983068928979090236736516465734108250477648738873e-1");
    rs_remainder[5][13]=str_to_Double("-.5586040472642019344273587677564412197909161368941257096202408745");
    rs_remainder[5][14]=str_to_Double("-.1124197636805911539678843978960873452852771901007544894484920030");
    rs_remainder[5][15]=str_to_Double(".1521267771179559182929594014480911341729511939344523291000325154");
    rs_remainder[5][16]=str_to_Double(".5173718845528038784023625510664489108383008846701324059173062168e-1");
    rs_remainder[5][17]=str_to_Double("-.2561227689700728294043343196049761335148654471052070423314194961e-1");
    rs_remainder[5][18]=str_to_Double("-.1296367251404617794428713962276982428546663350262095795910638942e-1");
    rs_remainder[5][19]=str_to_Double(".2545557481861163278061927441881923854541391028667670480670925048e-2");
    rs_remainder[5][20]=str_to_Double(".2119331951087777528850732134144290697104365325973636685647821547e-2");
    rs_remainder[5][21]=str_to_Double("-.919139194515677754051761292342158746748012565228746726008202390e-4");
    rs_remainder[5][22]=str_to_Double("-.2441346653385527265704955250928993826325905348522116083486654489e-3");
    rs_remainder[5][23]=str_to_Double("-.1369798269228338712377224163763806177794253043726382002567986883e-4");
    rs_remainder[5][24]=str_to_Double(".2062078503328423806426076694852059952894030142771333713542236005e-4");
    rs_remainder[5][25]=str_to_Double(".2817724196357407322489162817540872607359396128181385112871488638e-5");
    rs_remainder[5][26]=str_to_Double("-.1297434492844366892788263656446682046863355347251652923744601059e-5");
    rs_remainder[5][27]=str_to_Double("-.2855985999713360746253679530199828297787916260926629212627561241e-6");
    rs_remainder[5][28]=str_to_Double(".5996997144661264597007764789892657941213829334590929374099272044e-7");
    rs_remainder[5][29]=str_to_Double(".2021348781334734777703179524299182106100008976692958768295549687e-7");
    rs_remainder[5][30]=str_to_Double("-.1862487313117858697546423357131803394713530729064442250777739639e-8");
    rs_remainder[5][31]=str_to_Double("-.1094398866128935069379459248278640402118159318666804151562244175e-8");
    rs_remainder[5][32]=str_to_Double(".2097177143058943481490331635090271638727097767277815068185182062e-10");
    rs_remainder[5][33]=str_to_Double(".4716061490690834725326687908589104778556664079204462285238923027e-10");
    rs_remainder[5][34]=str_to_Double(".1660848125399721185577729739569608257984716378048390230870592218e-11");
    rs_remainder[5][35]=str_to_Double("-.1648934964877150791326314781634704659279443751546512065027590421e-11");
    rs_remainder[5][36]=str_to_Double("-.1370048598826852761793798904433521625593266096361000404512010497e-12");
    rs_remainder[5][37]=str_to_Double(".4704583253855874122419598994231168136843721590046948356771437678e-13");
    rs_remainder[5][38]=str_to_Double(".6278234664212507900954230558352656378592012888192752505729235781e-14");
    rs_remainder[5][39]=str_to_Double("-.1084019160311783064284715994462769454143543070278268662743342373e-14");
    rs_remainder[5][40]=str_to_Double("-.2154967254218182570800585975897301417438958906874026297245256116e-15");
    rs_remainder[5][41]=str_to_Double(".1930486830439871756306385408742024773103767408295489779791121655e-16");
    rs_remainder[5][42]=str_to_Double(".5997429222622066113666307058165348595049087860661892556546072710e-17");
    rs_remainder[5][43]=str_to_Double("-.2243115269799842199123485798824839891673653462967363858475648161e-18");
    rs_remainder[5][44]=str_to_Double("-.1399960448814351064103665114014619211581494361521577553045693524e-18");
    rs_remainder[5][45]=str_to_Double("-.1573328024432791737709616113099065837009974058352280262947149302e-21");
    rs_remainder[5][46]=str_to_Double(".2787807175667093170381390351984794306901551826244171590111480379e-20");
    rs_remainder[5][47]=str_to_Double(".9357632174576192753892888577385646406223791663047908443875957536e-22");
    rs_remainder[5][48]=str_to_Double("-.4772495089145947151762674742391260157960251655880704145680791266e-22");
    rs_remainder[5][49]=str_to_Double("-.3083117830565907986165246371614013195127398136752116634119946924e-23");
    rs_remainder[5][50]=str_to_Double(".7022826029624433026117658919774867867914968721212057626467043761e-24");
    rs_remainder[5][51]=str_to_Double(".6930786438671730599400032968509700643301773009880529848721035390e-25");
    rs_remainder[5][52]=str_to_Double("-.8789181023477559386086828335765965382448544732654489901676967324e-26");
    rs_remainder[5][53]=str_to_Double("-.1254417835885284531308591700188383115878428764035128397559756414e-26");
    rs_remainder[5][54]=str_to_Double(".9059268824299502037634782532574550741819020496640796180153352658e-28");
    rs_remainder[5][55]=str_to_Double(".1930011825643519774997795736281543853466001561423743562921712459e-28");
    rs_remainder[5][56]=str_to_Double("-.6964403302695697121647242479038030237658945583303185129687080740e-30");
    rs_remainder[5][57]=str_to_Double("-.2588001648141067119188746234917477842616803133035223347431479623e-30");
    rs_remainder[5][58]=str_to_Double(".2265034962924886974324950280823532915235781966789222479945611751e-32");
    rs_remainder[5][59]=str_to_Double(".3064286857860032598707231022984075009897000160538619570925626586e-32");
    rs_remainder[5][60]=str_to_Double(".4370580087904805862346243600085832347525100682317932289506986663e-34");
    rs_remainder[5][61]=str_to_Double("-.3224880701747766128532439360454992032206377973011201545910340297e-34");
    rs_remainder[5][62]=str_to_Double("-.1150385610069575017758592152022344131395183603604066382575009773e-35");
    rs_remainder[5][63]=str_to_Double(".3021738083317245590294358706235866037858460483976999994485219433e-36");
    rs_remainder[5][64]=str_to_Double(".1744833828214870483187728248910485021536914138168665958468576613e-37");
    rs_remainder[5][65]=str_to_Double("-.2511795139281475206625186084486680377537986954150010755199392439e-38");
    rs_remainder[5][66]=str_to_Double("-.2091040798312279857356439648180256939618508529592165829933169256e-39");
    rs_remainder[5][67]=str_to_Double(".1830240263076620162602634276403427479718285642541744724225397052e-40");
    rs_remainder[5][68]=str_to_Double(".2142070573196337307664134855455450755715103842642726051635520040e-41");
    rs_remainder[5][69]=str_to_Double("-.1134909769625905495662059744216350960035474159116599959766344093e-42");
    rs_remainder[5][70]=str_to_Double("-.1937406004565893102607560453283145408740559472336979152853843892e-43");
    rs_remainder[5][71]=str_to_Double(".5516036048474587480469824943264955580700529132765099496780734509e-45");
    //============================= n = 6 ================================================
    rs_remainder[6][0]=str_to_Double(".3369099840108093441362067822587549375066729534975297568975124893e-4");
    rs_remainder[6][1]=str_to_Double("-.4873038727737406480295032472233028432563632085005725216145266475e-3");
    rs_remainder[6][2]=str_to_Double(".349130411512094926417402699315129393670294393856515674741323471e-2");
    rs_remainder[6][3]=str_to_Double("-.1063618141082453509111878788291748882930066705104248284514335167e-1");
    rs_remainder[6][4]=str_to_Double("-.7962052861482918178126557852123920026799507312595925465677646540e-2");
    rs_remainder[6][5]=str_to_Double(".1237587562368654036956689016651519869792596687711062662300882420");
    rs_remainder[6][6]=str_to_Double("-.1849404122581204893704871913488802579844861136191330271075214910");
    rs_remainder[6][7]=str_to_Double("-.3039358023967954827404643682409931865632086817321452347044462009");
    rs_remainder[6][8]=str_to_Double(".7612833126395631593858738695847551136161118485823422639212486192");
    rs_remainder[6][9]=str_to_Double(".4067440568556812181503317961426565759576454095859639868278072570");
    rs_remainder[6][10]=str_to_Double("-1.230172180854170894524558619247602788934931978244710434314527213");
    rs_remainder[6][11]=str_to_Double("-.511764085569652212484019879135123414096589013315228506913328479");
    rs_remainder[6][12]=str_to_Double(".9962463615472546763398167390372930181848687971836693499065588187");
    rs_remainder[6][13]=str_to_Double(".4705671616186110701138140986022424692156551718511813611064183666");
    rs_remainder[6][14]=str_to_Double("-.4414445866526113678663497828277119590003617242870310357564017074");
    rs_remainder[6][15]=str_to_Double("-.2591849331053527531268263958546430019384073445613348299485560122");
    rs_remainder[6][16]=str_to_Double(".1111768899354234391658542954393085189354025221539071211183481031");
    rs_remainder[6][17]=str_to_Double(".8794868546608422363440501450752219755242093088753561478003880587e-1");
    rs_remainder[6][18]=str_to_Double("-.1480327188610340568392786361596650532654937177121919635717442354e-1");
    rs_remainder[6][19]=str_to_Double("-.1961041509857540989375996120205385942162165137767343493246919658e-1");
    rs_remainder[6][20]=str_to_Double(".3165009964103191689966684960715715971912528221062961947635971537e-3");
    rs_remainder[6][21]=str_to_Double(".3027401420822915542180871970066622504418900234482105094536991271e-2");
    rs_remainder[6][22]=str_to_Double(".2675580524826069148186449642704906293010115129525888211662781586e-3");
    rs_remainder[6][23]=str_to_Double("-.3349691662856057683766762689886514150581375109521905256161205961e-3");
    rs_remainder[6][24]=str_to_Double("-.5891092180278285448464916604326564743977718510985204134783673410e-4");
    rs_remainder[6][25]=str_to_Double(".2696637124579705432341111363195226755009379204013550864774914954e-4");
    rs_remainder[6][26]=str_to_Double(".7272656793940831261384020726473889082352202777885477625190386987e-5");
    rs_remainder[6][27]=str_to_Double("-.1554913482011107290164014792092043360483791926064948845866653128e-5");
    rs_remainder[6][28]=str_to_Double("-.6327824247833729259487234225727375215923746466910161761986684644e-6");
    rs_remainder[6][29]=str_to_Double(".5795829140308803407619064803234753373235208251458548688669783223e-7");
    rs_remainder[6][30]=str_to_Double(".4183986878632716358512315042641233265914971306069977165865493309e-7");
    rs_remainder[6][31]=str_to_Double("-.5919718812989810592372697697669629323937964442708188216148235381e-9");
    rs_remainder[6][32]=str_to_Double("-.2179917361243006977534793805068682675317724607837657165835877212e-8");
    rs_remainder[6][33]=str_to_Double("-.9616778969107069123165787279127995830391267806544261056090873273e-10");
    rs_remainder[6][34]=str_to_Double(".9118051994604895199514405775033642332764181665780583461868406468e-10");
    rs_remainder[6][35]=str_to_Double(".8727140632665969006881597826257348821532456071911835670474668940e-11");
    rs_remainder[6][36]=str_to_Double("-.3079465314837680607818338002809228472492269593086574809821537210e-11");
    rs_remainder[6][37]=str_to_Double("-.4631465693789377240871472950163143104247115525053484479740628129e-12");
    rs_remainder[6][38]=str_to_Double(".8309120714986330773009072453496190739383094745645986902556754389e-13");
    rs_remainder[6][39]=str_to_Double(".1846131870254684616144460260203846225604405224846354274154864766e-13");
    rs_remainder[6][40]=str_to_Double("-.1710444988418924939251014536177866207357329400747634828305986284e-14");
    rs_remainder[6][41]=str_to_Double("-.5943186199715446270133309466852221928950862801890696145477251031e-15");
    rs_remainder[6][42]=str_to_Double(".2234686781892984202553971807599516355499590192616738752698802079e-16");
    rs_remainder[6][43]=str_to_Double(".1595996125012342222020015818522130736854063993697651513558972698e-16");
    rs_remainder[6][44]=str_to_Double(".5095142324320878914905531798153290836345740411930184939507802857e-19");
    rs_remainder[6][45]=str_to_Double("-.3635333622938522891142084273840452846720724805427394858541595731e-18");
    rs_remainder[6][46]=str_to_Double("-.1366445990800952052044889203304050952111517155736274758397111301e-19");
    rs_remainder[6][47]=str_to_Double(".7078193679409975268806739749200259345566792857046871353046570767e-20");
    rs_remainder[6][48]=str_to_Double(".4991815709595988529196514651188904623920511029811512701359224488e-21");
    rs_remainder[6][49]=str_to_Double("-.1178041071520382429135426832834065334628704961413101733891910816e-21");
    rs_remainder[6][50]=str_to_Double("-.1257858085990709190400371153174270236499255385581562060621299485e-22");
    rs_remainder[6][51]=str_to_Double(".1658120287971612571588844964422223148595786657229993812419403780e-23");
    rs_remainder[6][52]=str_to_Double(".2550258889093815884447042940696640184817758533195955352736736333e-24");
    rs_remainder[6][53]=str_to_Double("-.1909921414436297570359760739243996940552145699143001289219985737e-25");
    rs_remainder[6][54]=str_to_Double("-.4383020889237005617785280597522604188883369608270221857869787486e-26");
    rs_remainder[6][55]=str_to_Double(".1624227654860267017846256480520145417364308435276122803978865872e-27");
    rs_remainder[6][56]=str_to_Double(".6543273717865298202572588824679517278678693514588473113756391562e-28");
    rs_remainder[6][57]=str_to_Double("-.5498318829780754197907860611152250207323256662676830678645713595e-30");
    rs_remainder[6][58]=str_to_Double("-.8595705350898431647186361503577041975902715211254815191669456901e-30");
    rs_remainder[6][59]=str_to_Double("-.1351632159976319176182909118619126945746362246733663438547973751e-31");
    rs_remainder[6][60]=str_to_Double(".1000260113304237776986895784338352013230132608680782861975609750e-31");
    rs_remainder[6][61]=str_to_Double(".3817096982395822119363234283590050429304603369275629288170783107e-33");
    rs_remainder[6][62]=str_to_Double("-.1032916157081309317448138118355417373917442999592399261510696635e-33");
    rs_remainder[6][63]=str_to_Double("-.6330803096933377249089783664972227947950412452285849655851737971e-35");
    rs_remainder[6][64]=str_to_Double(".9431649388172835745739282985479887039210428381201630609473077951e-36");
    rs_remainder[6][65]=str_to_Double(".8303883112307444921924782841747355646449900519217092196421182523e-37");
    rs_remainder[6][66]=str_to_Double("-.7524381234587142820011878299254934633010832170529172687456627525e-38");
    rs_remainder[6][67]=str_to_Double("-.9296822226775147640122256958407706234620782224268291709312248294e-39");
    rs_remainder[6][68]=str_to_Double(".5089568861907786778842429809074938576238222558133138606965476794e-40");
    rs_remainder[6][69]=str_to_Double(".9171031188778719888156699584022612771137458997780191857929450267e-41");
    rs_remainder[6][70]=str_to_Double("-.2683636790413776342290546468519923831167859849541271668983782016e-42");
    rs_remainder[6][71]=str_to_Double("-.8102120001770112472856656017558753670657464110296286720754380710e-43");
    //============================= n = 7 ================================================
    rs_remainder[7][0]=str_to_Double(".661247991827990458049746373433313070922531163215730472605721038e-4");
    rs_remainder[7][1]=str_to_Double("-.4467040957733873337506070694246162214508175497041172973358109342e-3");
    rs_remainder[7][2]=str_to_Double(".108402320680893125596527781997613448978493176413599955371941299e-2");
    rs_remainder[7][3]=str_to_Double(".5028543891765805977995962556174682509319071651609414977609598934e-2");
    rs_remainder[7][4]=str_to_Double("-.388614855153086405858180918510636012316864325587678966995342461e-1");
    rs_remainder[7][5]=str_to_Double(".77079567414100727304423627138374079685435615842216077450349588e-1");
    rs_remainder[7][6]=str_to_Double(".635596974406339688749284441132411910117352645452372160669351512e-1");
    rs_remainder[7][7]=str_to_Double("-.407459627303950790998848963594477476324482254043503971063186953");
    rs_remainder[7][8]=str_to_Double(".180337521119586400633904729347255753625837878219708919669982063");
    rs_remainder[7][9]=str_to_Double(".806430248560645287342294611351272286031115332452363879977275353");
    rs_remainder[7][10]=str_to_Double("-.517835801831444421078377800420096392932087606536289958079751334");
    rs_remainder[7][11]=str_to_Double("-.948227179582071540024994573250021236058309141361867958860643474");
    rs_remainder[7][12]=str_to_Double(".4719558561190385010296183967332796054845211767053407875801053398");
    rs_remainder[7][13]=str_to_Double(".7154101522815643706225164007518739467486364626873060968554799166");
    rs_remainder[7][14]=str_to_Double("-.1929666274743222194464450367746942580888543084062762843681687864");
    rs_remainder[7][15]=str_to_Double("-.3455616621306974043043210220577291427713887414339271290847147194");
    rs_remainder[7][16]=str_to_Double(".2925214982106184809592537830130667159633231954367349237733353701e-1");
    rs_remainder[7][17]=str_to_Double(".1090930109969402379235325117175732359789889140506327013714848863");
    rs_remainder[7][18]=str_to_Double(".530846569235287248925265612273940634351897823852723909669909366e-2");
    rs_remainder[7][19]=str_to_Double("-.2329344363057829379108824158412218458159063978798067098170069653e-1");
    rs_remainder[7][20]=str_to_Double("-.3516663416181302692165542468873629173446991249506306832321780017e-2");
    rs_remainder[7][21]=str_to_Double(".3464823737000550636845925098337627427690074264125476822591333269e-2");
    rs_remainder[7][22]=str_to_Double(".8444826314873811786456431059324788267392654369201950998259015491e-3");
    rs_remainder[7][23]=str_to_Double("-.3639543345208182880784780716589371852714507467138615029676040438e-3");
    rs_remainder[7][24]=str_to_Double("-.1277945452163360895114915037644285286929907615164417241699911593e-3");
    rs_remainder[7][25]=str_to_Double(".2645611034419118536266006019080363699190134214443607236231236850e-4");
    rs_remainder[7][26]=str_to_Double(".1381224527283913540588745909761865798902654096277882245132824019e-4");
    rs_remainder[7][27]=str_to_Double("-.1166448033035099648195465005747185266648303288982226588896607128e-5");
    rs_remainder[7][28]=str_to_Double("-.1128008010585149068611335660373361554341365757137896326202493991e-5");
    rs_remainder[7][29]=str_to_Double(".4291445682960736814794383938905543141795370996482948806004869734e-8");
    rs_remainder[7][30]=str_to_Double(".7183011442560939382266899752607954916216148850031251630812677043e-7");
    rs_remainder[7][31]=str_to_Double(".4170671957474236004438985010584367457995136867753402025690356774e-8");
    rs_remainder[7][32]=str_to_Double("-.3629407210769605673301310030361376241435561663865467778050791511e-8");
    rs_remainder[7][33]=str_to_Double("-.4124827698855936141100797363179030230492985827101553699337615458e-9");
    rs_remainder[7][34]=str_to_Double(".1463063425008355236824236842428696515885816486097719106294346996e-9");
    rs_remainder[7][35]=str_to_Double(".2538568050776675706244088745754913210137098178543441167734552112e-10");
    rs_remainder[7][36]=str_to_Double("-.4651815505981773057619837363284631161205841502605093485185122187e-11");
    rs_remainder[7][37]=str_to_Double("-.1180447181823237917552708093352868022725821554124398632473600007e-11");
    rs_remainder[7][38]=str_to_Double(".1109491490124069838522889043593547500875861081310677026428011795e-12");
    rs_remainder[7][39]=str_to_Double(".4419451281952939507262857430455806517309432081024162832544084064e-13");
    rs_remainder[7][40]=str_to_Double("-.1608478314692817457646430025618066008135823942882319203495958542e-14");
    rs_remainder[7][41]=str_to_Double("-.1372760528978133677452787761459005446446416788675322289175116208e-14");
    rs_remainder[7][42]=str_to_Double("-.9086382245607191715080190600732218298369306606849440408660847856e-17");
    rs_remainder[7][43]=str_to_Double(".3595073949423357697679907365079175890156113516918640735826046780e-16");
    rs_remainder[7][44]=str_to_Double(".1547586226971373490568001508409332163980733934446219575084401535e-17");
    rs_remainder[7][45]=str_to_Double("-.7998943503611832961665467683731967890021913230670371753491578006e-18");
    rs_remainder[7][46]=str_to_Double("-.6244382266774407959348653002567984915458722353982349594816358119e-19");
    rs_remainder[7][47]=str_to_Double(".1512013439917889744955571490256804547073079196872566022266841282e-19");
    rs_remainder[7][48]=str_to_Double(".1766128759196134967813271099346415404674576847808886053742902456e-20");
    rs_remainder[7][49]=str_to_Double("-.2401693347993949933432615682475046566512063332229057692759212722e-21");
    rs_remainder[7][50]=str_to_Double("-.4021771026997475598702188609888496811858628489699491268026067427e-22");
    rs_remainder[7][51]=str_to_Double(".3098079996533695589787589508459210254039032680910636482744525405e-23");
    rs_remainder[7][52]=str_to_Double(".7744398621026775072727826138550526439562492507818493248643448038e-24");
    rs_remainder[7][53]=str_to_Double("-.2909981194257874154768724205152330917959192046062653839481371171e-25");
    rs_remainder[7][54]=str_to_Double("-.1291082624999983239957836263818205874761241512882644513533973139e-25");
    rs_remainder[7][55]=str_to_Double(".9807515174106098654941507139834393551955798156670047905205143018e-28");
    rs_remainder[7][56]=str_to_Double(".1887333137378436886849534951012611488633314120007043892834029053e-27");
    rs_remainder[7][57]=str_to_Double(".3337379339311476664796512071359233556187688462326621455642103500e-29");
    rs_remainder[7][58]=str_to_Double("-.2435245390390282198995175473065294503038275552521377958661641041e-29");
    rs_remainder[7][59]=str_to_Double("-.1003631088133617771724002997698913112591514205472443897221805723e-30");
    rs_remainder[7][60]=str_to_Double(".2778613665867153558389058374142338003778961175427033723641899610e-31");
    rs_remainder[7][61]=str_to_Double(".1819986993653430338051062755397912920873024893625180618767476379e-32");
    rs_remainder[7][62]=str_to_Double("-.2793546435221929707687870745098653683232150185355166326198647002e-33");
    rs_remainder[7][63]=str_to_Double("-.2616240434791790041615866882782344184840413326730224422325601643e-34");
    rs_remainder[7][64]=str_to_Double(".2444838365551365240276705768221768651450480298531634742050818354e-35");
    rs_remainder[7][65]=str_to_Double(".3206790733180182291418478975399294534607353269711249227467579576e-36");
    rs_remainder[7][66]=str_to_Double("-.1806328497603823869039976371457812899309898598577289164487492181e-37");
    rs_remainder[7][67]=str_to_Double("-.3456706892212677392599500933314400898500362039488357324485606698e-38");
    rs_remainder[7][68]=str_to_Double(".1032974534105303305961166465592694123245457866832884822012046492e-39");
    rs_remainder[7][69]=str_to_Double(".3329686648686949868859199302541080735033863180051876334138450276e-40");
    rs_remainder[7][70]=str_to_Double("-.3000829886017220073170851514030268387630975596272358629192991445e-42");
    rs_remainder[7][71]=str_to_Double("-.2892036076727618558540539424567653118105208633761207162638876874e-42");
    //============================= n = 8 ================================================
    rs_remainder[8][0]=str_to_Double(".2419753613611796494579307999781989797564314445520583825307023191e-5");
    rs_remainder[8][1]=str_to_Double("-.161135227707040515320034104919326519389534204238072485359860555e-4");
    rs_remainder[8][2]=str_to_Double(".217180825329948501527320529679229707677062766592229533037367006e-3");
    rs_remainder[8][3]=str_to_Double("-.234415555034883341984196720617466494306936554633722546410667893e-2");
    rs_remainder[8][4]=str_to_Double(".115526317963676488921750533227471284924000688951598776921835867e-1");
    rs_remainder[8][5]=str_to_Double("-.2392447916109697179933463939763139246141596336331153299192225029e-1");
    rs_remainder[8][6]=str_to_Double("-.1553080439636881371492114696816252106253682408268673858375961942e-1");
    rs_remainder[8][7]=str_to_Double(".168054572159558936954740597170917457769895061144591064528064365");
    rs_remainder[8][8]=str_to_Double("-.207678931024271260913231432594353874524366994924598599636596762");
    rs_remainder[8][9]=str_to_Double("-.27051056273432962066986475859254594164779091002682615928607293");
    rs_remainder[8][10]=str_to_Double(".6603245174239842309412423925938244169438123150833972560003950701");
    rs_remainder[8][11]=str_to_Double(".1748462936027383322939933814326744590207756402710286068959825818");
    rs_remainder[8][12]=str_to_Double("-.9023846153142662280326148216826598372079768481974082016035446012");
    rs_remainder[8][13]=str_to_Double("-.9498217340106973193382779269377694984778099422737592994465980997e-1");
    rs_remainder[8][14]=str_to_Double(".7024752578716655015874175759138485362789031813535900012593122101");
    rs_remainder[8][15]=str_to_Double(".92687598380457370647763391117313096900804296389131910123489664e-1");
    rs_remainder[8][16]=str_to_Double("-.3377367994776856130030180421675088338677498137351726536666414152");
    rs_remainder[8][17]=str_to_Double("-.684030352687347166519736821592708602223062916477107946470542401e-1");
    rs_remainder[8][18]=str_to_Double(".1049550768889425174212015202010015528620341303321821917042918630");
    rs_remainder[8][19]=str_to_Double(".2965526457990261642520012481651756646282285476522058285463865458e-1");
    rs_remainder[8][20]=str_to_Double("-.2178017896559150545661822080342468889718757691109041585922248832e-1");
    rs_remainder[8][21]=str_to_Double("-.8154606200067587184415758535378944049446464985985335622912862668e-2");
    rs_remainder[8][22]=str_to_Double(".3059089916511288332613685836077729327156741902523070377001926696e-2");
    rs_remainder[8][23]=str_to_Double(".1535812550622837006269242853624889540024452899990316430179211476e-2");
    rs_remainder[8][24]=str_to_Double("-.2820164500854979909831584311188410400800938564790490836810061536e-3");
    rs_remainder[8][25]=str_to_Double("-.2091780782247145789697429974333680293048337808438283552071705772e-3");
    rs_remainder[8][26]=str_to_Double(".1385348732230897416735703428519229321124262393795603893232324020e-4");
    rs_remainder[8][27]=str_to_Double(".2139724177292542128222894020990172918930690751621997701335285497e-4");
    rs_remainder[8][28]=str_to_Double(".3300672135711366070046201082798348022732273453537439571381929197e-6");
    rs_remainder[8][29]=str_to_Double("-.1686934540219032128952523504213973667544744151471101322569874289e-5");
    rs_remainder[8][30]=str_to_Double("-.1352846282130496583509276461806437235224168418842873216159038337e-6");
    rs_remainder[8][31]=str_to_Double(".1041373739421406092425798648009587087258428109577929523707181522e-6");
    rs_remainder[8][32]=str_to_Double(".1457894343123659628662712498052625283957262873998805286327325164e-7");
    rs_remainder[8][33]=str_to_Double("-.5057361584708640309782462512030552672347484220533712091381715377e-8");
    rs_remainder[8][34]=str_to_Double("-.1043268891915953221269035381874564148813580976969853159979260837e-8");
    rs_remainder[8][35]=str_to_Double(".1906936980846097103664431782114090167963269295160619917834465832e-9");
    rs_remainder[8][36]=str_to_Double(".5689148160741357046045305461553582868048458645687125919321097513e-10");
    rs_remainder[8][37]=str_to_Double("-.5271318630772878794929230608740181866314466685349110762226037954e-11");
    rs_remainder[8][38]=str_to_Double("-.2491989495032678996520555126949694488209273367148923230945341567e-11");
    rs_remainder[8][39]=str_to_Double(".8238228532680538616464567618144168140493705556489575675219964814e-13");
    rs_remainder[8][40]=str_to_Double(".9007401495214419721106117047967990995068559817680842714350860422e-13");
    rs_remainder[8][41]=str_to_Double(".1083873681739085921992291228596999332253907242098286527983871514e-14");
    rs_remainder[8][42]=str_to_Double("-.2727604743096212702376042527747674846083512483013430885546483449e-14");
    rs_remainder[8][43]=str_to_Double("-.1378732838834482965271885025322789082306885549527045152777108868e-15");
    rs_remainder[8][44]=str_to_Double(".6971051630397714654465783960173385533160110051634629287127040610e-16");
    rs_remainder[8][45]=str_to_Double(".6124616772218422294134297053865581907870708806071215718384610261e-17");
    rs_remainder[8][46]=str_to_Double("-.1503346615206087725552567320922803341318013238510753025880924390e-17");
    rs_remainder[8][47]=str_to_Double("-.1947105325947110714186941973395169267155043159290380223417175505e-18");
    rs_remainder[8][48]=str_to_Double(".2704359958252626821487200411886801167099382100887528640016412511e-19");
    rs_remainder[8][49]=str_to_Double(".4993964532829073201021088181564376354626717661487122031186525547e-20");
    rs_remainder[8][50]=str_to_Double("-.3913733717110084808750436451236420762463474499098289423480414342e-21");
    rs_remainder[8][51]=str_to_Double("-.1080927120432324227203707776769199314589534605723092091599462218e-21");
    rs_remainder[8][52]=str_to_Double(".4045812478658757497223171871298707541029104365373671573075442551e-23");
    rs_remainder[8][53]=str_to_Double(".2018979016416460018772775156887709072435640399514585307089775391e-23");
    rs_remainder[8][54]=str_to_Double("-.1241837099059039145237836783771899912114026493656852787212433261e-25");
    rs_remainder[8][55]=str_to_Double("-.3294693460906857832243851880904833337315846524088631789025523584e-25");
    rs_remainder[8][56]=str_to_Double("-.6691652356564758862530667643133317026540594306221394546022797434e-27");
    rs_remainder[8][57]=str_to_Double(".4727954223851406170834075365557376968155178274243955808465721734e-27");
    rs_remainder[8][58]=str_to_Double(".2127726816768317861513618283372793264480785849042440925345745272e-28");
    rs_remainder[8][59]=str_to_Double("-.5977076931705627866574062814114799976115523757634419835752268406e-29");
    rs_remainder[8][60]=str_to_Double("-.4217685450958975135550222172469847259827862186025898019152261778e-30");
    rs_remainder[8][61]=str_to_Double(".6632601056328305560579767956441299893727857459091867283034149358e-31");
    rs_remainder[8][62]=str_to_Double(".6653349415010250735329750692474082402059271262668644020863913956e-32");
    rs_remainder[8][63]=str_to_Double("-.6380564610268546494341234069101173629074423807493547953257650424e-33");
    rs_remainder[8][64]=str_to_Double("-.8944499105208653964844515619975310525027949450854280792611044969e-34");
    rs_remainder[8][65]=str_to_Double(".5155455668869606936908244365963146370079352055295601520452262668e-35");
    rs_remainder[8][66]=str_to_Double(".1055616036494807930891933018049471005024815707948181697840655161e-35");
    rs_remainder[8][67]=str_to_Double("-.3194438182225804971720772806255501282243654527349569742591192387e-37");
    rs_remainder[8][68]=str_to_Double("-.1110872393075272748382008548138442407214727241882240111123955592e-37");
    rs_remainder[8][69]=str_to_Double(".9560947781982692852840109791834595536774183072644647770410290674e-40");
    rs_remainder[8][70]=str_to_Double(".1051683994493524478694709828209870603526746253338371639509238847e-39");
    rs_remainder[8][71]=str_to_Double(".1054989626160314762952970062377498817344097792927229813696590797e-41");
    //============================= n = 9 ================================================
    rs_remainder[9][0]=str_to_Double(".1376824100605469009292335676394516239888861367521744944722496262e-4");
    rs_remainder[9][1]=str_to_Double("-.108364270244188672629688513496529594439311519833524032467488744e-3");
    rs_remainder[9][2]=str_to_Double(".69612874087076726169990683663026654888436215664693745671877686e-3");
    rs_remainder[9][3]=str_to_Double("-.281532866123007466276180755344147262153999786510880531230243647e-2");
    rs_remainder[9][4]=str_to_Double(".5211281628124410228192697472653099980997356157495377808818186088e-2");
    rs_remainder[9][5]=str_to_Double(".78567631757874179365894532788622301586864409284439037616441030e-2");
    rs_remainder[9][6]=str_to_Double("-.64877602152968781227523025658710969316513257696629688352369781e-1");
    rs_remainder[9][7]=str_to_Double(".115185965479761063118114592985233127894675601948208938660286640");
    rs_remainder[9][8]=str_to_Double(".489527561061987324005191037501187984512355581522185218899875882e-1");
    rs_remainder[9][9]=str_to_Double("-.398461953562061221652052487637042192216962943609214552355476511");
    rs_remainder[9][10]=str_to_Double(".238884035084224009773504834958609098621476128970168436108191737");
    rs_remainder[9][11]=str_to_Double(".554851452505780727587580620470739470017863769833974503414736627");
    rs_remainder[9][12]=str_to_Double("-.5214374827560175018429225363178229362055500552561288022978670073");
    rs_remainder[9][13]=str_to_Double("-.486302248016181898655679208925925416694911618274638579176511555");
    rs_remainder[9][14]=str_to_Double(".4731537582420426817840706387616517989002381783186655188007774479");
    rs_remainder[9][15]=str_to_Double(".3150572871288756587959013579044756471126933988833807236671834268");
    rs_remainder[9][16]=str_to_Double("-.2397340129568678155898446659852318800070883460952231754196686160");
    rs_remainder[9][17]=str_to_Double("-.1501053083905259975796775868251289749386291059870070975457248196");
    rs_remainder[9][18]=str_to_Double(".7418116204376772293468409683330984352811593706131756416542396037e-1");
    rs_remainder[9][19]=str_to_Double(".5092186596602571200602606441703345371189010849337906546124767327e-1");
    rs_remainder[9][20]=str_to_Double("-.1438850581640066164057834165356595295045527962709924234004705169e-1");
    rs_remainder[9][21]=str_to_Double("-.1229766806834080213617534373206101762152132276200855722078937484e-1");
    rs_remainder[9][22]=str_to_Double(".1648126474815363314901282073311734681615822708242178349631603320e-2");
    rs_remainder[9][23]=str_to_Double(".2157082055276682867558326352115977635360068895873063467141127601e-2");
    rs_remainder[9][24]=str_to_Double("-.6290960205561819658436423525258380097271759931359724807667212051e-4");
    rs_remainder[9][25]=str_to_Double("-.2815677670395356128881643987063682569177933319013642786648135792e-3");
    rs_remainder[9][26]=str_to_Double("-.1404322368507607505690421123069111596730570834255120611760815941e-4");
    rs_remainder[9][27]=str_to_Double(".2792943537597993942591524855771476755260539318003087377450572712e-4");
    rs_remainder[9][28]=str_to_Double(".3243566185480260313546597570613943505173731559689835626055123482e-5");
    rs_remainder[9][29]=str_to_Double("-.2135094804475155589440420261485192875665265506866251004680990752e-5");
    rs_remainder[9][30]=str_to_Double("-.3850204796740483147734145865018022633610330675017917779236015697e-6");
    rs_remainder[9][31]=str_to_Double(".1262321946494359714948568480504453358170412345158328588980915030e-6");
    rs_remainder[9][32]=str_to_Double(".3220457404531872942975615368436846596858712472580421623286092808e-7");
    rs_remainder[9][33]=str_to_Double("-.5678387745876475048044755187432310046425415527370502574122709391e-8");
    rs_remainder[9][34]=str_to_Double("-.2072956031473247965550826190180574054350367949571372039778288862e-8");
    rs_remainder[9][35]=str_to_Double(".1810952142102159831073943755455042752110834917888368766202362597e-9");
    rs_remainder[9][36]=str_to_Double(".1069763684880553597652090473077403056780586039675762944405421972e-9");
    rs_remainder[9][37]=str_to_Double("-.2840745400507622301092809947144508190083660398600264222261788912e-11");
    rs_remainder[9][38]=str_to_Double("-.4530282452974083317998490605899976250655160073990906283552443834e-11");
    rs_remainder[9][39]=str_to_Double("-.9217631705481312968588159440510193783003425677982619153357132349e-13");
    rs_remainder[9][40]=str_to_Double(".1596370000250060551005105579516467354530147127913694604402303230e-12");
    rs_remainder[9][41]=str_to_Double(".9732296766819740330968115413679316677569917851712719769239451314e-14");
    rs_remainder[9][42]=str_to_Double("-.4713204718412881774187374277089942288161708820149043344984249877e-14");
    rs_remainder[9][43]=str_to_Double("-.4752160023605214767639760749828577424632124440242646972672224866e-15");
    rs_remainder[9][44]=str_to_Double(".1165242824659787975643776097322552175876309152502475892051705874e-15");
    rs_remainder[9][45]=str_to_Double(".1701201508355189107182246634570690241260554751793801386291021139e-16");
    rs_remainder[9][46]=str_to_Double("-.2382434239389920307915207012787691693263299817253114854239896814e-17");
    rs_remainder[9][47]=str_to_Double("-.4929586782976705453979777280353248473863374534424510928824978691e-18");
    rs_remainder[9][48]=str_to_Double(".3872523475978005687352036567529375007881759640398224158904173731e-19");
    rs_remainder[9][49]=str_to_Double(".1203539251254132102045735024534160500669825526987887906836740990e-19");
    rs_remainder[9][50]=str_to_Double("-.4374561572024989786637978966861274209634820377116682224002821664e-21");
    rs_remainder[9][51]=str_to_Double("-.2527582373149387813746644366584512932243939918968873487507987811e-21");
    rs_remainder[9][52]=str_to_Double(".967946239481972177147849265251849680770961772225612311274547252e-24");
    rs_remainder[9][53]=str_to_Double(".4620262267655266682025646465098969062801358560549653032601568582e-23");
    rs_remainder[9][54]=str_to_Double(".1101092245101671508038021766977079574443192255029845135301517082e-24");
    rs_remainder[9][55]=str_to_Double("-.7397532219261315734952410561039929065884428290523520127098899988e-25");
    rs_remainder[9][56]=str_to_Double("-.3681236161279897199657740806695464428125877473068413404053339093e-26");
    rs_remainder[9][57]=str_to_Double(".1039222173805369377059845511761777968066119350440925917832825592e-26");
    rs_remainder[9][58]=str_to_Double(".7974917714183369230166283745922806812452636098250456977658259616e-28");
    rs_remainder[9][59]=str_to_Double("-.1276103493805737661280516610905634620927081036299817745701146607e-28");
    rs_remainder[9][60]=str_to_Double("-.1382415981710439957450224132979427980646269670058083186172403647e-29");
    rs_remainder[9][61]=str_to_Double(".1352137017073816565480185301969592289177646096360895797940153767e-30");
    rs_remainder[9][62]=str_to_Double(".2042194342739036853370062381649934007935000942982755957829663713e-31");
    rs_remainder[9][63]=str_to_Double("-.1195986498970525953924789837568432210189545939654611520438345722e-32");
    rs_remainder[9][64]=str_to_Double("-.2644238004507043826176214824754839912082243752765180310422344921e-33");
    rs_remainder[9][65]=str_to_Double(".8013921228289212738061259760229178112258837385234612671129320464e-35");
    rs_remainder[9][66]=str_to_Double(".3046365532143794971423427454439098865121783489450713987978984141e-35");
    rs_remainder[9][67]=str_to_Double("-.2399459435225171428237338620117825978733799915463900433522831375e-37");
    rs_remainder[9][68]=str_to_Double("-.3149984606211840968590435531713400823832065776911351074103569831e-37");
    rs_remainder[9][69]=str_to_Double("-.3614163533711639354211069269893247280598871887067492351197528107e-39");
    rs_remainder[9][70]=str_to_Double(".2937283365044645134464868937521564679376516136588658865760464954e-39");
    rs_remainder[9][71]=str_to_Double(".8677336575533969969884560178354170784055361856253746172052836303e-41");
}

void initialize_rs_remainder3(){

    //============================= n = 10 ================================================
    rs_remainder[10][0]=str_to_Double("-.200010251733325102715374940882139315136190707141381121364212981e-6");
    rs_remainder[10][1]=str_to_Double(".109915017824018848094463468639833101973945193417857621283750681e-4");
    rs_remainder[10][2]=str_to_Double("-.101662428071697269285977946002664140944912759196747738635708191e-3");
    rs_remainder[10][3]=str_to_Double(".39098972580504788716075893715441097933956302052801722323973164e-3");
    rs_remainder[10][4]=str_to_Double("-.2673636266785313170826751501191872019695118081462935379886635e-3");
    rs_remainder[10][5]=str_to_Double("-.48379817156048473476016690749796717568810433979144107907805178e-2");
    rs_remainder[10][6]=str_to_Double(".2472958071818887415651879906086388532647352146387746327246728095e-1");
    rs_remainder[10][7]=str_to_Double("-.49296166385310937277469980362586951336310951296581514798168420e-1");
    rs_remainder[10][8]=str_to_Double("-.113949855227625323512777892015695208339683915691552519326902e-3");
    rs_remainder[10][9]=str_to_Double(".184701867717972746394085681708334948038488100996531551645535725");
    rs_remainder[10][10]=str_to_Double("-.25652183258970756679842869250040678085230601433944366642876626");
    rs_remainder[10][11]=str_to_Double("-.1537136523491344716789118450166422360574956407745594678894168021");
    rs_remainder[10][12]=str_to_Double(".5595524945513237940253969451758250565370829650944609065867556081");
    rs_remainder[10][13]=str_to_Double("-.5196110807521286395885843155822692660020246679224218765603963253e-1");
    rs_remainder[10][14]=str_to_Double("-.6005581196365656400182026498926514782930391618233048146349808237");
    rs_remainder[10][15]=str_to_Double(".1447901603090104752531594220428840303057390154230600982715063563");
    rs_remainder[10][16]=str_to_Double(".409077663263046482611545614498970082538917526581099826289633185");
    rs_remainder[10][17]=str_to_Double("-.88226894891095460986660178828648269088238045241585185017616135e-1");
    rs_remainder[10][18]=str_to_Double("-.1905302290152115054318272010908567336233547120941922978433699616");
    rs_remainder[10][19]=str_to_Double(".2479786856070810205337174143377731659640535758461250747679361124e-1");
    rs_remainder[10][20]=str_to_Double(".622859398931108331245640011802062558141511273007901405499633822e-1");
    rs_remainder[10][21]=str_to_Double("-.2394281857139954170720651709315519763006730527670054736903368574e-2");
    rs_remainder[10][22]=str_to_Double("-.1457863465975384080461590199040836082323715212636740913871239911e-1");
    rs_remainder[10][23]=str_to_Double("-.6176319463115360612225718524659588546078763228714235710028329441e-3");
    rs_remainder[10][24]=str_to_Double(".2494689886531393690107451916768315218379772563077489055566982858e-2");
    rs_remainder[10][25]=str_to_Double(".2800494678003457673734935524421432877546666667468460088419293981e-3");
    rs_remainder[10][26]=str_to_Double("-.3180982936548568630770414676905002321571498577314709218886230114e-3");
    rs_remainder[10][27]=str_to_Double("-.5634491871741518997024490035201732277080427684568932550415348681e-4");
    rs_remainder[10][28]=str_to_Double(".3062431690880365343700527956693953818237480587902591138628233994e-4");
    rs_remainder[10][29]=str_to_Double(".7528918412754401007992949572336405140910150006969852293415074012e-5");
    rs_remainder[10][30]=str_to_Double("-.2230237520897740475532574098716783494657552075464066422652203963e-5");
    rs_remainder[10][31]=str_to_Double("-.7429880360858513566558141324525080388466342845825521021947802280e-6");
    rs_remainder[10][32]=str_to_Double(".1201540201322338150102931085920210694229161780194381252516405042e-6");
    rs_remainder[10][33]=str_to_Double(".5693779860296284444932269205556093958597345980359688452350521715e-7");
    rs_remainder[10][34]=str_to_Double("-.4343982540597992795703031043395894444368667896275990444442741789e-8");
    rs_remainder[10][35]=str_to_Double("-.3490824878496632779997950817764277597489871700371112151205892455e-8");
    rs_remainder[10][36]=str_to_Double(".5362166673179062641879783206255229425733732983522578620877874441e-10");
    rs_remainder[10][37]=str_to_Double(".1745730307526982170014959683245304178874260097695675726505456084e-9");
    rs_remainder[10][38]=str_to_Double(".5747766032462888050204814881038293399283755881313085539198498393e-11");
    rs_remainder[10][39]=str_to_Double("-.7210071865734769384926806262198387593855945283264528062413443060e-11");
    rs_remainder[10][40]=str_to_Double("-.5452600618879835801078620189164563809872691770026905834539338744e-12");
    rs_remainder[10][41]=str_to_Double(".2474771874353507340886944292036320098403424334993606148387266490e-12");
    rs_remainder[10][42]=str_to_Double(".2929558775066477615851447581479350239791661487231643106311958610e-13");
    rs_remainder[10][43]=str_to_Double("-.7050629734215047884691521858267403799365468957000821562350024061e-14");
    rs_remainder[10][44]=str_to_Double("-.1183891200998738338866300385208670750463941263894134232757182578e-14");
    rs_remainder[10][45]=str_to_Double(".1644132536945956232299718861645194645231287388597567296873927645e-15");
    rs_remainder[10][46]=str_to_Double(".3889490100578073098548618086939324205197978806359465126575846816e-16");
    rs_remainder[10][47]=str_to_Double("-.3001248054431501521183824193793396655151674006766789671692407932e-17");
    rs_remainder[10][48]=str_to_Double("-.1075273284658564082378717263914081955573857688168413497878155519e-17");
    rs_remainder[10][49]=str_to_Double(".3653513564093827118282632929905360254603032221506506596362228287e-19");
    rs_remainder[10][50]=str_to_Double(".2548930403832848232743801180914552678816871267063121313590490532e-19");
    rs_remainder[10][51]=str_to_Double("-.81399853521447662111721950310747240724560527484077642728698898e-23");
    rs_remainder[10][52]=str_to_Double("-.5238599674828963393791558992818695910971412326748250146630830390e-21");
    rs_remainder[10][53]=str_to_Double("-.1495468338228999156166040896985658427897481705923010192746580707e-22");
    rs_remainder[10][54]=str_to_Double(".9390859276698796548395184138713790268291570036079774279323797663e-23");
    rs_remainder[10][55]=str_to_Double(".5240187565675688380339192326028855985156739331855888203856608138e-24");
    rs_remainder[10][56]=str_to_Double("-.1470600765911961933502995614886608145011637812546729789635628646e-24");
    rs_remainder[10][57]=str_to_Double("-.1240708663766456111843019339698890234258014275610242696106576929e-25");
    rs_remainder[10][58]=str_to_Double(".2003648584110135517712205657007365991193441319270161715188620816e-26");
    rs_remainder[10][59]=str_to_Double(".2366845959303346202846382487890727679621291510489916712170958183e-27");
    rs_remainder[10][60]=str_to_Double("-.2343075786547928548779916345011978292151783737975459635978823605e-28");
    rs_remainder[10][61]=str_to_Double("-.3849848627987263366911454084557128031823715588376810825958212031e-29");
    rs_remainder[10][62]=str_to_Double(".2270282252656306817116277357931043804939859083778879264318867862e-30");
    rs_remainder[10][63]=str_to_Double(".5480903329623832148190376120069041973393015913998787503808990175e-31");
    rs_remainder[10][64]=str_to_Double("-.1639380956343126364130389665888766190128627695063334142691154603e-32");
    rs_remainder[10][65]=str_to_Double("-.6928178643933973850925558631110725891359200598319665128063139279e-33");
    rs_remainder[10][66]=str_to_Double(".4658396421575405556883456247622233455496140292393081487029780184e-35");
    rs_remainder[10][67]=str_to_Double(".7841323079407138480061875538415373438136915285671728220056486984e-35");
    rs_remainder[10][68]=str_to_Double(".1046794084650898913903334863470432255835971002659967888805987809e-36");
    rs_remainder[10][69]=str_to_Double("-.7983260958428511830798046319455052866386799234078092532549570194e-37");
    rs_remainder[10][70]=str_to_Double("-.2566554026508714423927865219694942503868032061996890597547926600e-38");
    rs_remainder[10][71]=str_to_Double(".7322848046062478928494747155208361188850852663035252176081717418e-39");
    //============================= n = 11 ================================================
    rs_remainder[11][0]=str_to_Double(".21165343104016635137549519094363868539951860633078857488439109e-5");
    rs_remainder[11][1]=str_to_Double("-.611323345941696862992786502551499330269926543050509173698489669e-5");
    rs_remainder[11][2]=str_to_Double(".41527172154005453326037823463713344868962603755628951286617560e-4");
    rs_remainder[11][3]=str_to_Double("-.3868293393145709271403552429496153225818508874845636725828674542e-3");
    rs_remainder[11][4]=str_to_Double(".23730480240985707620631263674471708404945559084269904382307986e-2");
    rs_remainder[11][5]=str_to_Double("-.8852454011086692157099234552299179964732165421955553292704796584e-2");
    rs_remainder[11][6]=str_to_Double(".17413458638826571495247447712673856180241873061530936761052797e-1");
    rs_remainder[11][7]=str_to_Double("-.1271291600046961659528223276116329107606981618795457129835412741e-2");
    rs_remainder[11][8]=str_to_Double("-.7959730639378094834104612999516332520507357839188864962811333236e-1");
    rs_remainder[11][9]=str_to_Double(".16172806405795674034593950352932421802864814550666142201266305");
    rs_remainder[11][10]=str_to_Double("-.1942362458094586210342998958823967450212103848717432903176727e-1");
    rs_remainder[11][11]=str_to_Double("-.335953689520006483230136670756292061007673316096167237044177394");
    rs_remainder[11][12]=str_to_Double(".310594395387957567317177559633681794935246925609119513340660255");
    rs_remainder[11][13]=str_to_Double(".2929882624665191362008845047617326347700151248932018782785597983");
    rs_remainder[11][14]=str_to_Double("-.466712051558845099345098021261729240124204705223153730243897228");
    rs_remainder[11][15]=str_to_Double("-.1510551092580185978051970086084093177775940930073186147409771394");
    rs_remainder[11][16]=str_to_Double(".364666339737236362609864332935301735523070653154010490719675575");
    rs_remainder[11][17]=str_to_Double(".65105344040678932742000398953336595350162749012295255670675133e-1");
    rs_remainder[11][18]=str_to_Double("-.179660971556704541039400344698145091083204184875840662323238492");
    rs_remainder[11][19]=str_to_Double("-.2900643199150295389589825999934704828496853139825612252717650470e-1");
    rs_remainder[11][20]=str_to_Double(".599036571244580403075200134252934516769726143981305407261458985e-1");
    rs_remainder[11][21]=str_to_Double(".1131592835421748363973709828211329778222608013254007025579190430e-1");
    rs_remainder[11][22]=str_to_Double("-.1404734474141129503357148409730665963757359070010866298046218456e-1");
    rs_remainder[11][23]=str_to_Double("-.328375923302889685571779828183676856306047624978310440622479613e-2");
    rs_remainder[11][24]=str_to_Double(".2376938716065603686149606514617921783718571547288554688108638925e-2");
    rs_remainder[11][25]=str_to_Double(".6894697578080202611824283282023922446258891126768336085118855154e-3");
    rs_remainder[11][26]=str_to_Double("-.2947034080790354134269615175136545243687684924480591017062979819e-3");
    rs_remainder[11][27]=str_to_Double("-.1070535892593512696335106023954265224670536839303126730100638330e-3");
    rs_remainder[11][28]=str_to_Double(".267707258239237602584392420388893758165727000684505975311470904e-4");
    rs_remainder[11][29]=str_to_Double(".1265356938485185330851600626079086680170192074868796383922948631e-4");
    rs_remainder[11][30]=str_to_Double("-.1720736566460480470188739347300630270726197495842094315061616040e-5");
    rs_remainder[11][31]=str_to_Double("-.1168125155008184038263798403223465974179173097562647589290956316e-5");
    rs_remainder[11][32]=str_to_Double(".6629701445841325526650439937147264070532064172489474754241764754e-7");
    rs_remainder[11][33]=str_to_Double(".8597849983555671071153958914110228889938993256210873663034705790e-7");
    rs_remainder[11][34]=str_to_Double(".2470117013946503377872051978371390299648727876033949837353818272e-9");
    rs_remainder[11][35]=str_to_Double("-.5125558151094924321153212516985128889282010456144398065963599154e-8");
    rs_remainder[11][36]=str_to_Double("-.2660561132567809498645082404187023547429843265342347421028364742e-9");
    rs_remainder[11][37]=str_to_Double(".2501930865677914949406309107246465716288118852778277920762168325e-9");
    rs_remainder[11][38]=str_to_Double(".2417429432052166486607405595684247915490367938321894475898861173e-10");
    rs_remainder[11][39]=str_to_Double("-.1005452740771584439238783537123101159465235027135887210576425257e-10");
    rs_remainder[11][40]=str_to_Double("-.1435383152605079704227952959496032699133883290584671149982826276e-11");
    rs_remainder[11][41]=str_to_Double(".3319224601443056545035917212861920857932939786172825025350780964e-12");
    rs_remainder[11][42]=str_to_Double(".6572219032894278441302233284100273325992740344014155337745532342e-13");
    rs_remainder[11][43]=str_to_Double("-.8854315734716001895230528417307040477569471342604547470247449425e-14");
    rs_remainder[11][44]=str_to_Double("-.2458278216046970013854643491837685821044353257499190616015402957e-14");
    rs_remainder[11][45]=str_to_Double(".1810292649759519053424916451811804522439146711755050316279564321e-15");
    rs_remainder[11][46]=str_to_Double(".7729071598795982190422016759430116644857402228336233866018842459e-16");
    rs_remainder[11][47]=str_to_Double("-.2307443609207449495883460946513190721327201555508462770011842282e-17");
    rs_remainder[11][48]=str_to_Double("-.2076950864010546452303140387831702904196822969671235296281079769e-17");
    rs_remainder[11][49]=str_to_Double("-.9970144112277938496669281823697060263996965030184244515337276706e-20");
    rs_remainder[11][50]=str_to_Double(".4818939924905500107145576535482409381409535733226077953455560104e-19");
    rs_remainder[11][51]=str_to_Double(".1680019574374630907712106632045106405274119932226424082733395179e-20");
    rs_remainder[11][52]=str_to_Double("-.9708367800563526718558165921319045062926220040631615435887197514e-21");
    rs_remainder[11][53]=str_to_Double("-.6168506190638026842304800069016327681725020033988337593016323903e-22");
    rs_remainder[11][54]=str_to_Double(".1700431480139217025130095105364747252374716125159057646661854135e-22");
    rs_remainder[11][55]=str_to_Double(".1597059780289260197770657798814885160956242920570126696691659555e-23");
    rs_remainder[11][56]=str_to_Double("-.2577793364932476123390994291316696563833790949533373153647418058e-24");
    rs_remainder[11][57]=str_to_Double("-.3358484228486994760563337832476173401448348781018821295619987224e-25");
    rs_remainder[11][58]=str_to_Double(".3333304652653837050489565736087638339193363094626237549555783237e-26");
    rs_remainder[11][59]=str_to_Double(".6028083180842705758524525413392829874940030461262325703527082365e-27");
    rs_remainder[11][60]=str_to_Double("-.3538373039695676684435026950537949702543019579282605322922783535e-28");
    rs_remainder[11][61]=str_to_Double("-.9458315052688700216397260328918369065671257344302579519302424617e-29");
    rs_remainder[11][62]=str_to_Double(".2737008730350189074464039207504689795502083540456722628803960173e-30");
    rs_remainder[11][63]=str_to_Double(".1314904833507041921168653475263682840641939519319331173102371011e-30");
    rs_remainder[11][64]=str_to_Double("-.6610194282652660567896046256325876233095396330527464695250145006e-33");
    rs_remainder[11][65]=str_to_Double("-.1632706499807649420375254041917498135603875555897175983648331930e-32");
    rs_remainder[11][66]=str_to_Double("-.2575941234695098776172193403470145450441075813543431807987498170e-34");
    rs_remainder[11][67]=str_to_Double(".1818871629199218170515212704808351374788482350030463455018514365e-34");
    rs_remainder[11][68]=str_to_Double(".6428878790348618896033886590505515337485202357915681482256781624e-36");
    rs_remainder[11][69]=str_to_Double("-.1820633386868972593865333035443217538295013710113676222938353301e-36");
    rs_remainder[11][70]=str_to_Double("-.1007277805983749916988809401070044183584551127501243668772522259e-37");
    rs_remainder[11][71]=str_to_Double(".1634388722051081243541490318413254334175766174354475514492210062e-38");
    //============================= n = 12 ================================================
    rs_remainder[12][0]=str_to_Double("-.1508368668385969365261101701115484328803366208461720473455238975e-6");
    rs_remainder[12][1]=str_to_Double(".345152382706460117152765663592458026557684467501190689712565667e-5");
    rs_remainder[12][2]=str_to_Double("-.3163079970682537769033502841979296579690321071859990602923100828e-4");
    rs_remainder[12][3]=str_to_Double(".18477881939289544824372303332698653707716029541089389984030020e-3");
    rs_remainder[12][4]=str_to_Double("-.8212214504547483808956763929153827395208669082854166350647656e-3");
    rs_remainder[12][5]=str_to_Double(".26207728653465900521279963811134543293140232506215929906335097e-2");
    rs_remainder[12][6]=str_to_Double("-.4634609776800809563921295622231960177732942524437155009549200205e-2");
    rs_remainder[12][7]=str_to_Double("-.212609479329723124857933955616823699887921713116868491291690e-2");
    rs_remainder[12][8]=str_to_Double(".35350140836666518946515943021761013746598405234341261378205951e-1");
    rs_remainder[12][9]=str_to_Double("-.8330205498109722710226483427242331846744117616862788036563853422e-1");
    rs_remainder[12][10]=str_to_Double(".4769406016343326915089056595446703512144204436851113213912489495e-1");
    rs_remainder[12][11]=str_to_Double(".15220239667306243472894576346090334114924634750401642787889859");
    rs_remainder[12][12]=str_to_Double("-.29040977180861098482284614063481761322589754293496894969384671");
    rs_remainder[12][13]=str_to_Double(".6517184743767042644420031837330542056023853112381857798715400056e-3");
    rs_remainder[12][14]=str_to_Double(".41261621920563120474309157053885229307209929619237533962229093");
    rs_remainder[12][15]=str_to_Double("-.200056556918061040160349440589081206247483014875495140574353070");
    rs_remainder[12][16]=str_to_Double("-.3219527992056398144738301365566190779003299217236293437361181785");
    rs_remainder[12][17]=str_to_Double(".2209708796175579758108755336495269892706772225477668149170008566");
    rs_remainder[12][18]=str_to_Double(".1742653094845546642764614318954581591228407428674463728677836736");
    rs_remainder[12][19]=str_to_Double("-.1248936070733797457403411743330484173411017869635853498896344609");
    rs_remainder[12][20]=str_to_Double("-.7214496415900399695320416065116352016116831256252559072395068153e-1");
    rs_remainder[12][21]=str_to_Double(".4418908774924856448968024507279067683093304101937454980592341046e-1");
    rs_remainder[12][22]=str_to_Double(".231806846743802046651167734630583239120209742951494364143250067e-1");
    rs_remainder[12][23]=str_to_Double("-.105197194029303320316111123293269318698004655603427888450672606e-1");
    rs_remainder[12][24]=str_to_Double("-.5711108774671099848955186260229378280546652927621619962997863104e-2");
    rs_remainder[12][25]=str_to_Double(".173782072426392687129035851780963086793039396943696172348488003e-2");
    rs_remainder[12][26]=str_to_Double(".1074611664894627474963199529406444543305697426027101623997659183e-2");
    rs_remainder[12][27]=str_to_Double("-.198828141693152015417699185791900148492492882582624827327287811e-3");
    rs_remainder[12][28]=str_to_Double("-.1556608642299057655209033166399719799604539843793494041180000500e-3");
    rs_remainder[12][29]=str_to_Double(".1466931116407916765766869439503828136026894373819496322989860395e-4");
    rs_remainder[12][30]=str_to_Double(".1760073805465567197170382617065729446738912339274848112691959576e-4");
    rs_remainder[12][31]=str_to_Double("-.4280564574492022194778831583854638099088659927242574601249008249e-6");
    rs_remainder[12][32]=str_to_Double("-.1576655329537554036563197722036734415427796662703629441879256759e-5");
    rs_remainder[12][33]=str_to_Double("-.5053126924230625742983416984988015930465095523135127101840840875e-7");
    rs_remainder[12][34]=str_to_Double(".1133570903094462098947436752495674790744376166561465733733935014e-6");
    rs_remainder[12][35]=str_to_Double(".9168890421537814752597980066504092860817074070415606429002528240e-8");
    rs_remainder[12][36]=str_to_Double("-.6604549395252872443952485663541677015429857122764511751907632693e-8");
    rs_remainder[12][37]=str_to_Double("-.8421034534685374556761105029775852055527735230340562838594568791e-9");
    rs_remainder[12][38]=str_to_Double(".3132573664174607454129023990833689331878884940288741311346540512e-9");
    rs_remainder[12][39]=str_to_Double(".5570795860529623953560156546893149421531655009307822650658882468e-10");
    rs_remainder[12][40]=str_to_Double("-.1205086330206021787372095677462153809418922293562918617680573415e-10");
    rs_remainder[12][41]=str_to_Double("-.2905295877543638154159589097681835085956043085081194441321727209e-11");
    rs_remainder[12][42]=str_to_Double(".3683900909692115233894939012772635445208627550182924210238936899e-12");
    rs_remainder[12][43]=str_to_Double(".1243560795233411699039304718665125518692260722597420618804925110e-12");
    rs_remainder[12][44]=str_to_Double("-.8366822230029150478036838153268133728694845477602641315091404094e-14");
    rs_remainder[12][45]=str_to_Double("-.4469098219261017154562206608313486055335980022110286546734327230e-14");
    rs_remainder[12][46]=str_to_Double(".1043575983424116826620410605232985077784338612251970513486149163e-15");
    rs_remainder[12][47]=str_to_Double(".1367970296694513688166111406257136426493982701549668943071899949e-15");
    rs_remainder[12][48]=str_to_Double(".1656241290364946538959134903593915699702293814210641944350148345e-17");
    rs_remainder[12][49]=str_to_Double("-.3599406046049861643401029920017240798348002899970699667758924856e-17");
    rs_remainder[12][50]=str_to_Double("-.1560119210227633590246369509481771047075483048899340027724700468e-18");
    rs_remainder[12][51]=str_to_Double(".8183000259200118973240180241056723768044751775310952358139619714e-19");
    rs_remainder[12][52]=str_to_Double(".6019799262454336677635572017104756371101610739027155187764329429e-20");
    rs_remainder[12][53]=str_to_Double("-.1608810146836092338300298791620370012436271699955540060795654460e-20");
    rs_remainder[12][54]=str_to_Double("-.1706440289279893543300234778171123282063562695639794778384911972e-21");
    rs_remainder[12][55]=str_to_Double(".2721456888825307776340826883569365122942699996499650020591504887e-22");
    rs_remainder[12][56]=str_to_Double(".3963815957359582846395570528708032192500449740692576778647289119e-23");
    rs_remainder[12][57]=str_to_Double("-.3897899774423225597541266762960558540937752982816121603083258762e-24");
    rs_remainder[12][58]=str_to_Double("-.7869742636657603367726941868552952678574807536470126918389723652e-25");
    rs_remainder[12][59]=str_to_Double(".4529260370546856481562164123027446641071199051109592043462268233e-26");
    rs_remainder[12][60]=str_to_Double(".1364361106777330598452115990304564712767197923575131011537243039e-26");
    rs_remainder[12][61]=str_to_Double("-.3713912572719327714178209872314611577707831032677721801642543812e-28");
    rs_remainder[12][62]=str_to_Double("-.2091365532715337713592304490343478334003863253711437882106831996e-28");
    rs_remainder[12][63]=str_to_Double(".5507797562941319909030132235816520655386921974961309393507147053e-31");
    rs_remainder[12][64]=str_to_Double(".2855992618028541265083459333113573582804579429932089955531096662e-30");
    rs_remainder[12][65]=str_to_Double(".5399061562323227445061764534935868158303994839709319514191346949e-32");
    rs_remainder[12][66]=str_to_Double("-.3489508378492824204320637463070703446525430047929929672922305557e-32");
    rs_remainder[12][67]=str_to_Double("-.1370808405685465496593694956491893685902998378670901078575713764e-33");
    rs_remainder[12][68]=str_to_Double(".3819764194433909132530329491176408506223058113710780129185258550e-34");
    rs_remainder[12][69]=str_to_Double(".2300238162639776331563663095432004997606834883390315130974556207e-35");
    rs_remainder[12][70]=str_to_Double("-.3738387021035816976653963691858418830930777552155856036424150794e-36");
    rs_remainder[12][71]=str_to_Double("-.3132814189745311174375308706825675380284368145520483590633075020e-37");
    //============================= n = 13 ================================================
    rs_remainder[13][0]=str_to_Double(".301005458693798194582407315993381554197638724730300479381722863e-6");
    rs_remainder[13][1]=str_to_Double(".3585549411877917859036018845277689509287462003563722501922965e-5");
    rs_remainder[13][2]=str_to_Double("-.335307220852230665597159434915735752025857471767390503637178438e-4");
    rs_remainder[13][3]=str_to_Double(".1604144446717479749805474056101558028019692368678567542239030e-3");
    rs_remainder[13][4]=str_to_Double("-.477424412019625694987430826900770780828994798044648766864320e-3");
    rs_remainder[13][5]=str_to_Double(".4904594326842085978203173875963363236673479615199180263840500e-3");
    rs_remainder[13][6]=str_to_Double(".2882678293006571999461422374072091093789580611609656137083619e-2");
    rs_remainder[13][7]=str_to_Double("-.1628504031868478275935059417097396358645489513275304968820141341e-1");
    rs_remainder[13][8]=str_to_Double(".38383076843178547798786025584421752006242931509920854375010250e-1");
    rs_remainder[13][9]=str_to_Double("-.3195075621995552270263845368180978813748382286749083021206417044e-1");
    rs_remainder[13][10]=str_to_Double("-.6156222557857874015068671713632738029969441807980771315598023547e-1");
    rs_remainder[13][11]=str_to_Double(".1882279062468161510557575813925930412196050690005345534086294361");
    rs_remainder[13][12]=str_to_Double("-.1140940003016310787067414150934634859352083373828600535630917");
    rs_remainder[13][13]=str_to_Double("-.21335419278725302598787167831932020282007917566166017573893441");
    rs_remainder[13][14]=str_to_Double(".32650671240362276333288390148882778619857597139132620905375995");
    rs_remainder[13][15]=str_to_Double(".64771463214201300458981317085473728957347707555687780975941185e-1");
    rs_remainder[13][16]=str_to_Double("-.3402078786942617733049041137211334542074536109677754959361445200");
    rs_remainder[13][17]=str_to_Double(".43168443965751925736187372884706108283367422475964485355945715e-1");
    rs_remainder[13][18]=str_to_Double(".214094334879350578408273350390472358395748478870216969036258300");
    rs_remainder[13][19]=str_to_Double("-.46569210895565534896893896276298298724502994149779988566763747e-1");
    rs_remainder[13][20]=str_to_Double("-.9371243257043141144055948270983598882227583430376453663830026508e-1");
    rs_remainder[13][21]=str_to_Double(".195212644642370539733204743278025223835733932484922648329513528e-1");
    rs_remainder[13][22]=str_to_Double(".301107056898841253565851518369406417020454824451193097785515095e-1");
    rs_remainder[13][23]=str_to_Double("-.4666802979786240355985985015886228895748093735405273339882604302e-2");
    rs_remainder[13][24]=str_to_Double("-.7263830392552472838342253502926087743120300170105083717319507492e-2");
    rs_remainder[13][25]=str_to_Double(".6501385392587072383643969412608330902893152274040476903663918850e-3");
    rs_remainder[13][26]=str_to_Double(".1334711125403111430328103907574211060803423441104119481030478005e-2");
    rs_remainder[13][27]=str_to_Double("-.356058330903165438904811432778010871625871116278455032246310768e-4");
    rs_remainder[13][28]=str_to_Double("-.1893915600228281063887513434510545587641342661613338535438162449e-3");
    rs_remainder[13][29]=str_to_Double("-.55553022147161660391906509935131257701217445573182356070274192e-5");
    rs_remainder[13][30]=str_to_Double(".2103841559345797927584840507305546804454005092241240332394603748e-4");
    rs_remainder[13][31]=str_to_Double(".1670405725027517967595708834537229328693427409582619826373424574e-5");
    rs_remainder[13][32]=str_to_Double("-.1851960580329576588315921220175200874804939355222657942976545321e-5");
    rs_remainder[13][33]=str_to_Double("-.2342410311506625166245511885212276337878217594135223861810341647e-6");
    rs_remainder[13][34]=str_to_Double(".1303581647683373553512178072245901898270788043097696152251256195e-6");
    rs_remainder[13][35]=str_to_Double(".2279237447328949718934596745187816819601279721724937105080570233e-7");
    rs_remainder[13][36]=str_to_Double("-.7364291523208950086859746252238147769538621497650152308737870803e-8");
    rs_remainder[13][37]=str_to_Double("-.1700243405876827054849999808438268352389631064309188668758761121e-8");
    rs_remainder[13][38]=str_to_Double(".3319001627934553992057240947729621304117702899397046650564552539e-9");
    rs_remainder[13][39]=str_to_Double(".1017498735735324204598548577406541432324248742410517655124000055e-9");
    rs_remainder[13][40]=str_to_Double("-.1160983522982928413331607885574294839600787853894651866826484439e-10");
    rs_remainder[13][41]=str_to_Double("-.5015926812914263485711205801981708355294560405940322173364751244e-11");
    rs_remainder[13][42]=str_to_Double(".2865666532068730511192808727956460928317180246371155058323227159e-12");
    rs_remainder[13][43]=str_to_Double(".2072732982715559500980281249415822698524883632305400876608697838e-12");
    rs_remainder[13][44]=str_to_Double("-.281847416071920738287834435444190012727089802109477880013651183e-14");
    rs_remainder[13][45]=str_to_Double("-.7267048693591685321338034230324663539493871016412025242906448872e-14");
    rs_remainder[13][46]=str_to_Double("-.1629203855855047077295794172927551018441779082543931674701164090e-15");
    rs_remainder[13][47]=str_to_Double(".2179397391386059238891392726980861926147322252830108423285015844e-15");
    rs_remainder[13][48]=str_to_Double(".1194287951203185628099327204495586458375692595227620396358900809e-16");
    rs_remainder[13][49]=str_to_Double("-.5616438321281450916485967244559606987914834320938448265699271319e-17");
    rs_remainder[13][50]=str_to_Double("-.4871471050031139928447174506557913870815449204399454931259098610e-18");
    rs_remainder[13][51]=str_to_Double(".1244206254529418894604362498652516309626581721383800734327987870e-18");
    rs_remainder[13][52]=str_to_Double(".1515357479207126977650783224667218208211801249646514299707689288e-19");
    rs_remainder[13][53]=str_to_Double("-.2355296338162821723644633188111641062987487522453118523641396410e-20");
    rs_remainder[13][54]=str_to_Double("-.3897799603884405450778479888183399989504681816969455730874455979e-21");
    rs_remainder[13][55]=str_to_Double(".3741337777060713838611449948798126392220856305069530106396208186e-22");
    rs_remainder[13][56]=str_to_Double(".8583477152934329524703280890485231380506887873341471072107294221e-23");
    rs_remainder[13][57]=str_to_Double("-.4746713254600655354787135826971370566849736605288046358426907657e-24");
    rs_remainder[13][58]=str_to_Double("-.1648861777886659197536013938877771677394967736024679229929030030e-24");
    rs_remainder[13][59]=str_to_Double(".4048384669668228160081576567343409838913242216788958981947103387e-26");
    rs_remainder[13][60]=str_to_Double(".2794508349634659747220664211010393219982271241046456961847279956e-26");
    rs_remainder[13][61]=str_to_Double(".1818354408623866346609164305419693121117533458277970454273131080e-29");
    rs_remainder[13][62]=str_to_Double("-.4208205656613361892492277295880439685410803160446723641812700457e-28");
    rs_remainder[13][63]=str_to_Double("-.9645953361042346598094964136893217124895328661573491956795593255e-30");
    rs_remainder[13][64]=str_to_Double(".5653221102859016366452700523850876361464162585221603949043641444e-30");
    rs_remainder[13][65]=str_to_Double(".2496561611797540297169749884029370237785872729665175040692127952e-31");
    rs_remainder[13][66]=str_to_Double("-.6782725886135230104556275689313273864529052667306230280754563660e-32");
    rs_remainder[13][67]=str_to_Double("-.4487381739989302947895955991297408570172679824917795833474260643e-33");
    rs_remainder[13][68]=str_to_Double(".7251474731046680924394032124130654739129337604641587843665221370e-34");
    rs_remainder[13][69]=str_to_Double(".6610744683554431432982245459447150543598679412921420856638003743e-35");
    rs_remainder[13][70]=str_to_Double("-.6855503863223420342605783056478319386691212647910824905209454845e-36");
    rs_remainder[13][71]=str_to_Double("-.8425316415900249787155358128980557893403460059190477747328615519e-37");
    //============================= n = 14 ================================================
    rs_remainder[14][0]=str_to_Double("-.634876965749489559891104278350630941963454835330814654153775840e-7");
    rs_remainder[14][1]=str_to_Double(".59233662516197255693325149591217524992622851053719525610515890e-6");
    rs_remainder[14][2]=str_to_Double("-.219345125806225604826896256228929729529844418786099576814217e-5");
    rs_remainder[14][3]=str_to_Double(".52509656959864466812145736365556557890362796241549440565159e-5");
    rs_remainder[14][4]=str_to_Double("-.352730700322151181201340533688675564297229264562876495451923e-4");
    rs_remainder[14][5]=str_to_Double(".34261153534401559466084662360204783713257875270774349100480352e-3");
    rs_remainder[14][6]=str_to_Double("-.19935412597142656504528042296875086978188890460725011175725244e-2");
    rs_remainder[14][7]=str_to_Double(".72219153014732530960608712652568274830854654268933117255964980e-2");
    rs_remainder[14][8]=str_to_Double("-.158974827097533100372213630248505337788197368317105218307286e-1");
    rs_remainder[14][9]=str_to_Double(".1443274512919071753370466461148723111095588486773250739970863e-1");
    rs_remainder[14][10]=str_to_Double(".2674005306252058139709689683140544725052773397142743599590211e-1");
    rs_remainder[14][11]=str_to_Double("-.10470853517253971906671858058828277620750633323494278054546064");
    rs_remainder[14][12]=str_to_Double(".1121861592288677524191202034626954388069062039777199281729608595");
    rs_remainder[14][13]=str_to_Double(".661097624889320949645518788099962613891107035804962711752389e-1");
    rs_remainder[14][14]=str_to_Double("-.2663433346978368904675055090761909811822399223752337279675974929");
    rs_remainder[14][15]=str_to_Double(".1309626670195893609634358052687762936471675990144019778272006105");
    rs_remainder[14][16]=str_to_Double(".23016934313732815216982566631466066634617036530190824491482824");
    rs_remainder[14][17]=str_to_Double("-.237844031993620753932791403182949669100451889285972396226169924");
    rs_remainder[14][18]=str_to_Double("-.1054166822978624279645074916183477261406434230588707079299873086");
    rs_remainder[14][19]=str_to_Double(".1853399489311877693086847100163220193377788732952645605650074134");
    rs_remainder[14][20]=str_to_Double(".3035663992260268424523533761198290986952623896725555178384999e-1");
    rs_remainder[14][21]=str_to_Double("-.900798431752424370015936257732200435602246094281557344388578663e-1");
    rs_remainder[14][22]=str_to_Double("-.7357606767703472903310348111017918036090648223000343140025016e-2");
    rs_remainder[14][23]=str_to_Double(".304301720797921815057398029265862264272538579359876539669614457e-1");
    rs_remainder[14][24]=str_to_Double(".2177278034980529952883316825148338322217139027779048802534279317e-2");
    rs_remainder[14][25]=str_to_Double("-.7500911821308715917282775118535566326938417649404297950611096121e-2");
    rs_remainder[14][26]=str_to_Double("-.6842137072243453879396819676765575886798066207276922372486835106e-3");
    rs_remainder[14][27]=str_to_Double(".138768589146594698873303974368813567275737644375073204634708362e-2");
    rs_remainder[14][28]=str_to_Double(".1708360231362019607628385184988855630737279296656904510925995516e-3");
    rs_remainder[14][29]=str_to_Double("-.1965354344609274446818687538996731246034380127703412363978461152e-3");
    rs_remainder[14][30]=str_to_Double("-.3156448858850178939881517000868166313067417886456074137688343995e-4");
    rs_remainder[14][31]=str_to_Double(".2162854185504323040808778860340797558878841996413172904839305169e-4");
    rs_remainder[14][32]=str_to_Double(".4386224422929653102078212695485294668961605165716121470146704967e-5");
    rs_remainder[14][33]=str_to_Double("-.1868035633816018671281609990221237797616092375544677371177283307e-5");
    rs_remainder[14][34]=str_to_Double("-.4718974749798863510944243668288236432372352779973631011373201507e-6");
    rs_remainder[14][35]=str_to_Double(".1269986584228152799723788186858626684077502604017957248674543342e-6");
    rs_remainder[14][36]=str_to_Double(".4034046297438722661281329358029776584854782023962914239264597357e-7");
    rs_remainder[14][37]=str_to_Double("-.672870269333983775057767991968270791559776814351609234909060622e-8");
    rs_remainder[14][38]=str_to_Double("-.2798396154443932940040483518883683269896458947301563166854442540e-8");
    rs_remainder[14][39]=str_to_Double(".2664900929657003179617291181340737635744041231556194711239027306e-9");
    rs_remainder[14][40]=str_to_Double(".1601820834337158794947573807762338700234579483864822403619536016e-9");
    rs_remainder[14][41]=str_to_Double("-.670089783879033520445728526901232397629961192877016565187197082e-11");
    rs_remainder[14][42]=str_to_Double("-.7665932183588291831942450945922655011722779575464004995688854575e-11");
    rs_remainder[14][43]=str_to_Double("-.5137770282701548714308039296823630688084836726363263975116602e-14");
    rs_remainder[14][44]=str_to_Double(".3098385121615456533121547043840889324199728625695553977349812769e-12");
    rs_remainder[14][45]=str_to_Double(".1144585285992595453275618191242170954469527396352128134643549928e-13");
    rs_remainder[14][46]=str_to_Double("-.1065194965058705746087000506104424125106189523843708923365069662e-13");
    rs_remainder[14][47]=str_to_Double("-.7502697742028400130004458375216042345516770353388768535423298675e-15");
    rs_remainder[14][48]=str_to_Double(".3127205362102781980673161000690156179231919105714223377684801909e-15");
    rs_remainder[14][49]=str_to_Double(".3263009521103028929635447180708860893094249601769898822686210346e-16");
    rs_remainder[14][50]=str_to_Double("-.7837650370144907396516006579635424836313949197816377335392104117e-17");
    rs_remainder[14][51]=str_to_Double("-.1117685485466045706722035012325844617711130916714811359234056549e-17");
    rs_remainder[14][52]=str_to_Double(".1664893454646717650188278927297624899915720360000052375981371913e-18");
    rs_remainder[14][53]=str_to_Double(".3193436079615326070429291293712804315738285834802093211882738021e-19");
    rs_remainder[14][54]=str_to_Double("-.2934209909597740923820317030685996000109355539860689234528198577e-20");
    rs_remainder[14][55]=str_to_Double("-.7824624985499790299537941495521611355864102424930569899261371509e-21");
    rs_remainder[14][56]=str_to_Double(".4041690442083199862012365733275033728925739283645798497519616488e-22");
    rs_remainder[14][57]=str_to_Double(".1670678252401687069732859857589872247238365860032048000765898759e-22");
    rs_remainder[14][58]=str_to_Double("-.3453148575569658653952102134479465390376271518284700642910528953e-24");
    rs_remainder[14][59]=str_to_Double("-.3140138712997385039920716873812582129400225573348705927050271017e-24");
    rs_remainder[14][60]=str_to_Double("-.1593826897450697674522008135548060303119583230919525075848178390e-26");
    rs_remainder[14][61]=str_to_Double(".5229388558463073770099963738876580006970374327801888369317087608e-26");
    rs_remainder[14][62]=str_to_Double(".1468294220364702275364030685820510138853881758856070502687258193e-27");
    rs_remainder[14][63]=str_to_Double("-.7744571912038852265488444408209240192714614768468477188642298066e-28");
    rs_remainder[14][64]=str_to_Double("-.3890844843123859636548107852283739386304178271704796603996734988e-29");
    rs_remainder[14][65]=str_to_Double(".1020897724905672362724738346058109698211163779195257658780325174e-29");
    rs_remainder[14][66]=str_to_Double(".7496914018154406318452456584937241028933112135893802156399349804e-31");
    rs_remainder[14][67]=str_to_Double("-.1194721693913612424313182439293399207134247984193617573979137156e-31");
    rs_remainder[14][68]=str_to_Double("-.1196086250006163685722224245240891088312466541484383029372307214e-32");
    rs_remainder[14][69]=str_to_Double(".1230973029880040590413255724935008850700542304631433565820183528e-33");
    rs_remainder[14][70]=str_to_Double(".1655140441725888850586967541046523695982172888865240035318994265e-34");
    rs_remainder[14][71]=str_to_Double("-.1094966840838214586379202046741222633732963924660291549527952993e-35");
}

void initialize_rs_remainder4(){

    //============================= n = 15 ================================================
    rs_remainder[15][0]=str_to_Double(".434804329620086519789267724488997588729035911407664398405421e-7");
    rs_remainder[15][1]=str_to_Double(".1525605194860272908508420340414011740212939846406514468410131e-5");
    rs_remainder[15][2]=str_to_Double("-.1145688040703370298526165451693932930388300719188426987068933e-4");
    rs_remainder[15][3]=str_to_Double(".5519802938879232045283052195921220834421934717880196527715855e-4");
    rs_remainder[15][4]=str_to_Double("-.241437458445229585183902133535226102478490398074175963318028e-3");
    rs_remainder[15][5]=str_to_Double(".9220770310156539481220858741448358358036328887329105944607090013e-3");
    rs_remainder[15][6]=str_to_Double("-.2752404474028618330281121494823806540471363973092809718638675898e-2");
    rs_remainder[15][7]=str_to_Double(".552942357347273973234517128764720829020496334492095003931333319e-2");
    rs_remainder[15][8]=str_to_Double("-.415467375778968232117161411886929776280798445476377058139969e-2");
    rs_remainder[15][9]=str_to_Double("-.1427178896292178943204333225886802311898002616368208876248635414e-1");
    rs_remainder[15][10]=str_to_Double(".5498048370557003620983485170925955881973771349401914608940079e-1");
    rs_remainder[15][11]=str_to_Double("-.7637565339909961069203746245438850724583525718987955992942311813e-1");
    rs_remainder[15][12]=str_to_Double("-.2597141179636499642469646150711814662753799083861370639866797471e-2");
    rs_remainder[15][13]=str_to_Double(".16221168039324495177882513799413234504262371881070873526041993");
    rs_remainder[15][14]=str_to_Double("-.183601111470682525802410413087398521388394881776981506628520099");
    rs_remainder[15][15]=str_to_Double("-.623991948046159971487434423597822698634978865651153321393709996e-1");
    rs_remainder[15][16]=str_to_Double(".2622218357184940982625329729144209420284164352030455625207950100");
    rs_remainder[15][17]=str_to_Double("-.8079024114493556672459250072915239781998723355590797476228559e-1");
    rs_remainder[15][18]=str_to_Double("-.1845158837894114097826219049951757212477956220350493617084429119");
    rs_remainder[15][19]=str_to_Double(".109837653895598771761631621149344791390786582782522602993633919");
    rs_remainder[15][20]=str_to_Double(".8491059504531552355904886682335659707461819936941410771285200374e-1");
    rs_remainder[15][21]=str_to_Double("-.6521539396487817738487420349552171565985304679495860990655024186e-1");
    rs_remainder[15][22]=str_to_Double("-.2956069155770168643531431066170974575378785403975596581893958218e-1");
    rs_remainder[15][23]=str_to_Double(".24312899018572550236270583413852439436462770453382770736379253e-1");
    rs_remainder[15][24]=str_to_Double(".8410050308476713467844084910845550311658324065108316429244386665e-2");
    rs_remainder[15][25]=str_to_Double("-.6308121047873095878336156310262269997222558383396520203761270241e-2");
    rs_remainder[15][26]=str_to_Double("-.198402450436794218096966216874566331595805282175283861103078606e-2");
    rs_remainder[15][27]=str_to_Double(".1195115609821401698332109508731799695477350790634850542784794297e-2");
    rs_remainder[15][28]=str_to_Double(".381708069218145751701644435069431049859792686262128262848768617e-3");
    rs_remainder[15][29]=str_to_Double("-.1699041700630159562269034378446164689161628848188364931614401338e-3");
    rs_remainder[15][30]=str_to_Double("-.5901941077549909114883056751608833821302361680508597948130663163e-4");
    rs_remainder[15][31]=str_to_Double(".1840007451687004684436037618314013058031549195424342789519366494e-4");
    rs_remainder[15][32]=str_to_Double(".7315764688133816395485608106443396680181302921379616660760471547e-5");
    rs_remainder[15][33]=str_to_Double("-.152207674889485519365261197899935892912178665762169402980414415e-5");
    rs_remainder[15][34]=str_to_Double("-.731601151586106474849549021549096001996409578774215814819927609e-6");
    rs_remainder[15][35]=str_to_Double(".9439147392893215769889182234736455188066126487176173942404139794e-7");
    rs_remainder[15][36]=str_to_Double(".5963207566687846799155105786033279836879085894217059070839190634e-7");
    rs_remainder[15][37]=str_to_Double("-.4054236568525632387249334564081894376653927448711212466605383749e-8");
    rs_remainder[15][38]=str_to_Double("-.4005482959157948090377778933063165916839331866497546139860752637e-8");
    rs_remainder[15][39]=str_to_Double(".7729150580229152519815568419326416184730145835564142412940631366e-10");
    rs_remainder[15][40]=str_to_Double(".223985251941058583189443032270468751258730448640227178839915317e-9");
    rs_remainder[15][41]=str_to_Double(".479507284870432741169106432017394139773696623680800334415147192e-11");
    rs_remainder[15][42]=str_to_Double("-.1051662865641408241708155992993939708034794527686929622001629689e-10");
    rs_remainder[15][43]=str_to_Double("-.6052471843853957448746069842546124918016471282219795950454185955e-12");
    rs_remainder[15][44]=str_to_Double(".417211548642944373953586786190378470098642539442384875761464675e-12");
    rs_remainder[15][45]=str_to_Double(".3844528437233678595612181749982288087588456539620631585444552734e-13");
    rs_remainder[15][46]=str_to_Double("-.1403086962699868719364173840451128222697193866575566703809891430e-13");
    rs_remainder[15][47]=str_to_Double("-.1802143454394559472549054047557009430536070214948036199034870445e-14");
    rs_remainder[15][48]=str_to_Double(".3995247304918818167183758644784431932362222189046681117068710649e-15");
    rs_remainder[15][49]=str_to_Double(".6830284140156203630697714662770519576966668146089837383408792505e-16");
    rs_remainder[15][50]=str_to_Double("-.954540016094821443314761354037573474091969053209157136334347106e-17");
    rs_remainder[15][51]=str_to_Double("-.2176243518567399716354382768153737679407020624747557253000884353e-17");
    rs_remainder[15][52]=str_to_Double(".1863928980108733853232674643121024526446262405279035164466427190e-18");
    rs_remainder[15][53]=str_to_Double(".5954835265107553686945580340539341105696394575488325245771337521e-19");
    rs_remainder[15][54]=str_to_Double("-.2754056061643906755408226329834601453207164717511185478475573961e-20");
    rs_remainder[15][55]=str_to_Double("-.1418198857338740562199190351618909073208654590160955429956585798e-20");
    rs_remainder[15][56]=str_to_Double(".2159875989576581640467167619408925127018806983227228349720214162e-22");
    rs_remainder[15][57]=str_to_Double(".2966067530663158149388074221114862831506615995599337521417625867e-22");
    rs_remainder[15][58]=str_to_Double(".3253760756812655179842025583484952363854802079928868402470254346e-24");
    rs_remainder[15][59]=str_to_Double("-.5479637009371575994152997420749448635029437311853643897691772790e-24");
    rs_remainder[15][60]=str_to_Double("-.1901572918107549483547510157159288815767605433828754292763855193e-25");
    rs_remainder[15][61]=str_to_Double(".8971983097147223080730558911536399003448540718917679621250057996e-26");
    rs_remainder[15][62]=str_to_Double(".5191814980901962999448305886714645529511421224019018640090041901e-27");
    rs_remainder[15][63]=str_to_Double("-.1302714227417273786730502803999488502387434125976667385831074916e-27");
    rs_remainder[15][64]=str_to_Double("-.1073964571141773499425516220841972258538657971236522211581536153e-28");
    rs_remainder[15][65]=str_to_Double(".1672263313015638171763379733308572378000641106271738026981306449e-29");
    rs_remainder[15][66]=str_to_Double(".1858476597049012868266578107555232212722325389041365549703772238e-30");
    rs_remainder[15][67]=str_to_Double("-.1880504130238477414176617171108059424214176063979738830881896952e-31");
    rs_remainder[15][68]=str_to_Double("-.2797163823759525181007671953659174141793775235311574238854612174e-32");
    rs_remainder[15][69]=str_to_Double(".1812618409383961973072058237943455965808133986295227959905218033e-33");
    rs_remainder[15][70]=str_to_Double(".3735034742458605541914848015593711913013090269280748081963109177e-34");
    rs_remainder[15][71]=str_to_Double("-.1416921604309863397348479573500657042929816101082800133298855596e-35");
    //============================= n = 16 ================================================
    rs_remainder[16][0]=str_to_Double("-.25232821113668129668111639654714956334225022914255076576940814e-7");
    rs_remainder[16][1]=str_to_Double(".10338138665522437340811187925386356773370379302006894631820e-7");
    rs_remainder[16][2]=str_to_Double(".16361969973186584529446435664598611825327888699550466661808191e-5");
    rs_remainder[16][3]=str_to_Double("-.1338100801852404520207153396427405900387421176398010297142242987e-4");
    rs_remainder[16][4]=str_to_Double(".70183748501949013320671330429184348918669095301527123112456e-4");
    rs_remainder[16][5]=str_to_Double("-.2691924767609373615077373294291935138283117440930550643127328e-3");
    rs_remainder[16][6]=str_to_Double(".75158196927974707869745185626714203621703762132897524374134832e-3");
    rs_remainder[16][7]=str_to_Double("-.128669022750222419821900559022583866550831725036981541349371547e-2");
    rs_remainder[16][8]=str_to_Double("-.13900422618858683220948224746893218025995911416492950321070e-3");
    rs_remainder[16][9]=str_to_Double(".8748914618817418690560405610882723080059370138289042773686975343e-2");
    rs_remainder[16][10]=str_to_Double("-.2823666041148973174472509210830852334268732215907506387653372520e-1");
    rs_remainder[16][11]=str_to_Double(".435734008226501415250284123138571454589961441365352500732541e-1");
    rs_remainder[16][12]=str_to_Double("-.118321758743102782016959263240724186922109122911876905835477e-1");
    rs_remainder[16][13]=str_to_Double("-.857814690965989425748059098084105623199066902515760571948876e-1");
    rs_remainder[16][14]=str_to_Double(".15294900227518691896505077437292605914969250574250636647239296");
    rs_remainder[16][15]=str_to_Double("-.415335514293653616504260985852436996607931879780938702473244e-1");
    rs_remainder[16][16]=str_to_Double("-.1756188782948516902493359072481805863641478906509248850266250");
    rs_remainder[16][17]=str_to_Double(".185870310229733468662898295266146611435747307330484883573383093");
    rs_remainder[16][18]=str_to_Double(".6273604273732335693877310713761938262168390058783773558992680e-1");
    rs_remainder[16][19]=str_to_Double("-.1853808772299809589125897028739246272427381740950321986881973064");
    rs_remainder[16][20]=str_to_Double(".1952982097757173870441522959579507807317625557066747005509066e-1");
    rs_remainder[16][21]=str_to_Double(".105983602917737303895576067331213663971214568749909814805690367");
    rs_remainder[16][22]=str_to_Double("-.2854079588057870554349518792780854540941529284837277322755096974e-1");
    rs_remainder[16][23]=str_to_Double("-.42088094198849210146464058605156540805247819831890883235867836e-1");
    rs_remainder[16][24]=str_to_Double(".136194715395091387216847274910560368672804581758532522659703114e-1");
    rs_remainder[16][25]=str_to_Double(".1261900212607064789996303469706299030266209189133143897988824588e-1");
    rs_remainder[16][26]=str_to_Double("-.39308358969684493200213684044093443466880684223046094915337512e-2");
    rs_remainder[16][27]=str_to_Double("-.2962751064877650012395937461867466353465603053221240063240991713e-2");
    rs_remainder[16][28]=str_to_Double(".7756583104175614288512717923386441733463955531238940101213565377e-3");
    rs_remainder[16][29]=str_to_Double(".5524836672954665183443647439941939422683580720297345597587320105e-3");
    rs_remainder[16][30]=str_to_Double("-.109361933940654208062972346782689326639094931789231804622899072e-3");
    rs_remainder[16][31]=str_to_Double("-.8239999522811155878473757249018614021659247632846321863702587197e-4");
    rs_remainder[16][32]=str_to_Double(".1108096063271608579702748721556185402661060463695020125975774217e-4");
    rs_remainder[16][33]=str_to_Double(".9896872399200362120422174545718906231853493891274101315837906571e-5");
    rs_remainder[16][34]=str_to_Double("-.7680920766330038577852729784232995589377743423275243729552458089e-6");
    rs_remainder[16][35]=str_to_Double("-.9652555139925845484016583102084830929256818573877369933716568302e-6");
    rs_remainder[16][36]=str_to_Double(".2750265915375656847660344918011709604241044820593270710366949994e-7");
    rs_remainder[16][37]=str_to_Double(".7714732626367980307727348531207959202736207964758982736677675105e-7");
    rs_remainder[16][38]=str_to_Double(".108133405987197745826359395873763530239815829339920340546800388e-8");
    rs_remainder[16][39]=str_to_Double("-.5097879892054933264875120405916197709668443917855228648913912290e-8");
    rs_remainder[16][40]=str_to_Double("-.2645173128139683841029179177779085905839091211639462477997506053e-9");
    rs_remainder[16][41]=str_to_Double(".2806728006605569283711692851033292643014600324317285431652228927e-9");
    rs_remainder[16][42]=str_to_Double(".2452857038713191487393910001067485621376945965757968183722497883e-10");
    rs_remainder[16][43]=str_to_Double("-.1294893738177784113404064986308046781816278343771128218456276330e-10");
    rs_remainder[16][44]=str_to_Double("-.1594470414734308147501026147504940117663251072882589888998319166e-11");
    rs_remainder[16][45]=str_to_Double(".5019169216078299965487135318529197032830977776838394808457544133e-12");
    rs_remainder[16][46]=str_to_Double(".8158618766256252245578002348066312895677933838369787950706258925e-13");
    rs_remainder[16][47]=str_to_Double("-.1630583112071892485045682316932289281332422697067808138489290310e-13");
    rs_remainder[16][48]=str_to_Double("-.3443377787626976289340644824863103203800230843178737234456350776e-14");
    rs_remainder[16][49]=str_to_Double(".438697208785196068669137300337891837994162077495926577842308028e-15");
    rs_remainder[16][50]=str_to_Double(".1229428670539534397605645679100455492527372283512251388520948490e-15");
    rs_remainder[16][51]=str_to_Double("-.9442858580986280286312807796039439601311940265617926566398030798e-17");
    rs_remainder[16][52]=str_to_Double("-.3773201682188492144589015249878020255525779080938108532327084838e-17");
    rs_remainder[16][53]=str_to_Double(".1457048911794188830046817067252213020552656532436521644841218683e-18");
    rs_remainder[16][54]=str_to_Double(".1006383063385393191239497712055513703608822933501193136177444153e-18");
    rs_remainder[16][55]=str_to_Double("-.7831093369776714936337272944502950506448765085037985643676664642e-21");
    rs_remainder[16][56]=str_to_Double("-.2350813285696783818544938297427931701636100498823835855181253506e-20");
    rs_remainder[16][57]=str_to_Double("-.4414065754936126108812125401785070577439529317635270873706270668e-22");
    rs_remainder[16][58]=str_to_Double(".4834466103411788220269796823035254952897939275428546167276730900e-22");
    rs_remainder[16][59]=str_to_Double(".2090752930014892547759887125521687635910423757218413110866998830e-23");
    rs_remainder[16][60]=str_to_Double("-.8778415259192988671977976020483502507376484448671914781161126930e-24");
    rs_remainder[16][61]=str_to_Double("-.5927731733240214783412150155312736775371447493043477207540393156e-25");
    rs_remainder[16][62]=str_to_Double(".1407702805313513910028027378814276077023365601897592164511795466e-25");
    rs_remainder[16][63]=str_to_Double(".1319424184839528737159988120272227433077731487504666062820549279e-26");
    rs_remainder[16][64]=str_to_Double("-.1986325971940390153520151003177032465181901479814946889454699170e-27");
    rs_remainder[16][65]=str_to_Double("-.2481274153996121548820752083063328904808116655864423541042794915e-28");
    rs_remainder[16][66]=str_to_Double(".2440829578199159003128685555697604452688884516500798644496291934e-29");
    rs_remainder[16][67]=str_to_Double(".4069901578941669346767125002374343426564748289984023002419403470e-30");
    rs_remainder[16][68]=str_to_Double("-.2548163332444505124833677186824830867861469558670110039955431952e-31");
    rs_remainder[16][69]=str_to_Double("-.5923496051715518432951652657467586773541102597162114501036770468e-32");
    rs_remainder[16][70]=str_to_Double(".2117507147866338000663300913338027350316610085908943988507561940e-33");
    rs_remainder[16][71]=str_to_Double(".7730533865258583197740419968560902164789069706349187456078450316e-34");
    //============================= n = 17 ================================================
    rs_remainder[17][0]=str_to_Double(".1508733845352269450942701338238778405898169467177673702107962e-7");
    rs_remainder[17][1]=str_to_Double(".289903600011706500104466425071344636418279840901315548146661e-6");
    rs_remainder[17][2]=str_to_Double("-.8909425694286618512342871339630282530545523770042539733327e-6");
    rs_remainder[17][3]=str_to_Double("-.2947955935035210447875438955835852088450012068999799055285e-5");
    rs_remainder[17][4]=str_to_Double(".24406426725603480204720377573531746525718790444348046472909e-4");
    rs_remainder[17][5]=str_to_Double("-.616573376342053384078552250839961114398929985970269672966658e-4");
    rs_remainder[17][6]=str_to_Double("-.5753626414457513964742302424085463566719032608687174831528e-4");
    rs_remainder[17][7]=str_to_Double(".11543521382127133779462258312410802950705654746915404692184498e-2");
    rs_remainder[17][8]=str_to_Double("-.5233953693988480886590159081774134804269512768116954921087574e-2");
    rs_remainder[17][9]=str_to_Double(".1396251929435170281221417545138909318603967071062353544375745e-1");
    rs_remainder[17][10]=str_to_Double("-.217930013175789903201870521292460104054279365020114870785274e-1");
    rs_remainder[17][11]=str_to_Double(".8648714391217946927541232727941323643906912979207543598115749215e-2");
    rs_remainder[17][12]=str_to_Double(".44225746593606265651856465572506670133171368220977646499231402e-1");
    rs_remainder[17][13]=str_to_Double("-.10362798895080332921848137635099202823890070616173480293503337");
    rs_remainder[17][14]=str_to_Double(".7244890891036394897055829234326081769927068986357962771591427590e-1");
    rs_remainder[17][15]=str_to_Double(".8057783587339997718573718639965127975096283250255874773081452e-1");
    rs_remainder[17][16]=str_to_Double("-.1853599834369384918536744119901482068137255727092815850066360295");
    rs_remainder[17][17]=str_to_Double(".594281359562321938273174989517577852325863217378178923029280e-1");
    rs_remainder[17][18]=str_to_Double(".14662599769974310322985856703609016234320240035847800943967358");
    rs_remainder[17][19]=str_to_Double("-.1265550003372316472770950866232606403155517463545428566176387");
    rs_remainder[17][20]=str_to_Double("-.5563444416128160315779319457940054547622375009340195272777868e-1");
    rs_remainder[17][21]=str_to_Double(".9404771265488205253418102221432918254640551247939815255983938e-1");
    rs_remainder[17][22]=str_to_Double(".8592271296952544192476600562509598712705619695284849213981350e-2");
    rs_remainder[17][23]=str_to_Double("-.4304994874461247118145727066982793951025661394746484136673973e-1");
    rs_remainder[17][24]=str_to_Double(".12484312389529498028806592887409114320174806113314415623270381e-2");
    rs_remainder[17][25]=str_to_Double(".13957476716550914394002535495254274906739853001626150814305371e-1");
    rs_remainder[17][26]=str_to_Double("-.9026416508210596329349322849798278848727686036373760366038386209e-3");
    rs_remainder[17][27]=str_to_Double("-.3403950399430311882911313025771920565867237943977000318024023342e-2");
    rs_remainder[17][28]=str_to_Double(".2028326995167502151567574055733860920949384827337849290640523635e-3");
    rs_remainder[17][29]=str_to_Double(".6439334412727463254398015746918297344754776991500475132470195813e-3");
    rs_remainder[17][30]=str_to_Double("-.2262800057323340075897862781451803854317333808575190159532940666e-4");
    rs_remainder[17][31]=str_to_Double("-.9624413263671015249796321788190207062180538126875926885872125906e-4");
    rs_remainder[17][32]=str_to_Double(".2873326391577896301807188262262158699929080859418489860217196160e-6");
    rs_remainder[17][33]=str_to_Double(".1151889407272672605473905787472336171766368875696124268090496668e-4");
    rs_remainder[17][34]=str_to_Double(".3580550169546709859772277855908679762935475653501492658299866357e-6");
    rs_remainder[17][35]=str_to_Double("-.1116417260775391002063452830714480366872232580543752009119425311e-5");
    rs_remainder[17][36]=str_to_Double("-.72408031955972583287573143564356076848834111396241511485649549e-7");
    rs_remainder[17][37]=str_to_Double(".8848000206490392706108608644686367731242281342609968308892600493e-7");
    rs_remainder[17][38]=str_to_Double(".8688885079090160389977118474562914179877569932562439408436684120e-8");
    rs_remainder[17][39]=str_to_Double("-.5780255359359196300835086384815131825100985358745859930250726651e-8");
    rs_remainder[17][40]=str_to_Double("-.764536285716359448113858539013791466858767806870869524351557446e-9");
    rs_remainder[17][41]=str_to_Double(".3130223288205494768098134510326685434595611172901477963073117002e-9");
    rs_remainder[17][42]=str_to_Double(".5300159792017686406684854737727223870558152678324417938757709067e-10");
    rs_remainder[17][43]=str_to_Double("-.1407956076014579511464306274804124460349322136829852334847680381e-10");
    rs_remainder[17][44]=str_to_Double("-.3003098354027087389381514311436680745289895662697948052135033238e-11");
    rs_remainder[17][45]=str_to_Double(".5237629188204275143905783327447997583570921383827271829006059898e-12");
    rs_remainder[17][46]=str_to_Double(".1422861504261275124331725721839163647092038494771082298119896328e-12");
    rs_remainder[17][47]=str_to_Double("-.158433066247086038068810674745904200152313667126168485806607901e-13");
    rs_remainder[17][48]=str_to_Double("-.5727510845218945586180640344340015063235349495316101712526808929e-14");
    rs_remainder[17][49]=str_to_Double(".370639516947096995489881886264755491688165389682503490240997604e-15");
    rs_remainder[17][50]=str_to_Double(".1981853167834444133263501110555145602402683759847903429182000668e-15");
    rs_remainder[17][51]=str_to_Double("-.5571468792418714449724879187727158792082422953045569092902827700e-17");
    rs_remainder[17][52]=str_to_Double("-.5947471011110637311635665886113673857779372498857565934365046015e-17");
    rs_remainder[17][53]=str_to_Double("-.1367436956344527495556934927755435092281296166474717284749141357e-19");
    rs_remainder[17][54]=str_to_Double(".1558226301638915428337389214133946697692578245720509983339083291e-18");
    rs_remainder[17][55]=str_to_Double(".4536519736156092810562901235821977866275962281281983073883500324e-20");
    rs_remainder[17][56]=str_to_Double("-.3580714250136485072100335193985577105727721240710374622080156511e-20");
    rs_remainder[17][57]=str_to_Double("-.1945945493389012780364215075567680398535259925929063608223856686e-21");
    rs_remainder[17][58]=str_to_Double(".7234668692317907410016442652122816573673640541402637785703472749e-22");
    rs_remainder[17][59]=str_to_Double(".5780862740107711856292372648023640314736026697107960658090361461e-23");
    rs_remainder[17][60]=str_to_Double("-.1284871336906479751157190322882503568987504413659247748724861355e-23");
    rs_remainder[17][61]=str_to_Double("-.1388929354811616390753888248337071698055913631537899745789083434e-24");
    rs_remainder[17][62]=str_to_Double(".1996963243611553414291107468520319761550661869426805732973680925e-25");
    rs_remainder[17][63]=str_to_Double(".2845344337099495026998205765689509983412750161030127590128290659e-26");
    rs_remainder[17][64]=str_to_Double("-.2683692034590708135094271650868087899201866060258659469774382799e-27");
    rs_remainder[17][65]=str_to_Double("-.5097826562580504542916307189961626598838048129963155444565515730e-28");
    rs_remainder[17][66]=str_to_Double(".3029269246263507560026533008425187909691711097949175926427232572e-29");
    rs_remainder[17][67]=str_to_Double(".8105371988474197858300523291212639639680852386621436515612733811e-30");
    rs_remainder[17][68]=str_to_Double("-.2650702202077498899219779258126937041094093644283419965388094302e-31");
    rs_remainder[17][69]=str_to_Double("-.1154330992104625730290479965904467585725526598562144299518975514e-31");
    rs_remainder[17][70]=str_to_Double(".1249898250235006853399257069642073413087883996363968026287688301e-33");
    rs_remainder[17][71]=str_to_Double(".1481601282346930791540816430781736956858775626821619169279912180e-33");
    //============================= n = 18 ================================================
    rs_remainder[18][0]=str_to_Double("-.11214875444432193401482746314280340648692577466924778520996729e-7");
    rs_remainder[18][1]=str_to_Double("-.3317199243780702390322953897269008152983852559810168464169350e-7");
    rs_remainder[18][2]=str_to_Double(".844522340843154643006521601851836279547094990399279553687696e-6");
    rs_remainder[18][3]=str_to_Double("-.42295989907524188635533341081000026534118876405965380251725e-5");
    rs_remainder[18][4]=str_to_Double(".16172914518711025208184806193619150150799336894676682365226661e-4");
    rs_remainder[18][5]=str_to_Double("-.61426269371421028316948908584705578856549726658060218613899e-4");
    rs_remainder[18][6]=str_to_Double(".240294386935883215264808997179746468840079173802635077958701e-3");
    rs_remainder[18][7]=str_to_Double("-.8826392518752885436997866055325807823157596289727586194020e-3");
    rs_remainder[18][8]=str_to_Double(".27081731729009218941530930805645551321071594031971682335599e-2");
    rs_remainder[18][9]=str_to_Double("-.6270892645329922506708796377528065868794297050193135455249134348e-2");
    rs_remainder[18][10]=str_to_Double(".93030301007948184933417610483001474723086893844967323714152e-2");
    rs_remainder[18][11]=str_to_Double("-.2897593653702664867785280357429465868573256907396213404964865870e-2");
    rs_remainder[18][12]=str_to_Double("-.244306330325411807080913367469913629269016053055773098435226e-1");
    rs_remainder[18][13]=str_to_Double(".63727438787351106264295013657017285667124191542116062154640e-1");
    rs_remainder[18][14]=str_to_Double("-.64279890701527816426614213732089330629799668617153466540748e-1");
    rs_remainder[18][15]=str_to_Double("-.22172444834865865855379867043979081793122530899184705077172e-1");
    rs_remainder[18][16]=str_to_Double(".134968114527140280278287939425873760681639818445161518311041");
    rs_remainder[18][17]=str_to_Double("-.1162144666089393762137702781852379885012624604326697974301701281");
    rs_remainder[18][18]=str_to_Double("-.5615083733688657915162869318419762134108630953416873167749647815e-1");
    rs_remainder[18][19]=str_to_Double(".15722978567712162372402137451312617958205005998452307872951213");
    rs_remainder[18][20]=str_to_Double("-.41454133336952961047123718466309355324755640282768083094948644e-1");
    rs_remainder[18][21]=str_to_Double("-.9664236704399644697779434255815278513316942047902493592414482758e-1");
    rs_remainder[18][22]=str_to_Double(".58945628003113274002270842947456103107270313029290645093142177e-1");
    rs_remainder[18][23]=str_to_Double(".3646462831567729924472193088601354799711607684176081579815860e-1");
    rs_remainder[18][24]=str_to_Double("-.3381585318567707515267415174518909349061869599074406020454683e-1");
    rs_remainder[18][25]=str_to_Double("-.974666868178964686562048390737170850881920022043802555111217e-2");
    rs_remainder[18][26]=str_to_Double(".12333294443940063965797573245965678542910244918146539501396793e-1");
    rs_remainder[18][27]=str_to_Double(".21029503768595147192879208039857577649085565440384102630352567e-2");
    rs_remainder[18][28]=str_to_Double("-.32162169903703682032927192482466498739903428069551497745612068e-2");
    rs_remainder[18][29]=str_to_Double("-.4120672152022592335848279569830051347541922578530318947674733e-3");
    rs_remainder[18][30]=str_to_Double(".6321687480122072463469561427052023688891481010577805464069105474e-3");
    rs_remainder[18][31]=str_to_Double(".7595304184702296732569217697947565014100954701729497072677589e-4");
    rs_remainder[18][32]=str_to_Double("-.964894035075706831684325752703929764341560113608194861372860552e-4");
    rs_remainder[18][33]=str_to_Double("-.1249704652118926416455138130210251489282824230793051349079676999e-4");
    rs_remainder[18][34]=str_to_Double(".1166210233991812550506844448333438685650355156412326517931870889e-4");
    rs_remainder[18][35]=str_to_Double(".1730439456441528573692288604118203028245850013474894306919304559e-5");
    rs_remainder[18][36]=str_to_Double("-.113196863881744243540337605270709879609574218901210589934326767e-5");
    rs_remainder[18][37]=str_to_Double("-.196406281367024044923544596770177768216007234815408521100351876e-6");
    rs_remainder[18][38]=str_to_Double(".8914295866905872730063966809835642295238633261439230047239709973e-7");
    rs_remainder[18][39]=str_to_Double(".1822920786725796649006860821492806687641530822168983456956864645e-7");
    rs_remainder[18][40]=str_to_Double("-.5732104051601784719569618778323734255862018026474122339281249373e-8");
    rs_remainder[18][41]=str_to_Double("-.1394580030570260861994296653032241720890041043255262351373720008e-8");
    rs_remainder[18][42]=str_to_Double(".301352040102534939972802869413493287710963078702013735352275199e-9");
    rs_remainder[18][43]=str_to_Double(".8891113545551410839931676649512522386886033704761666672839504233e-10");
    rs_remainder[18][44]=str_to_Double("-.1285596079681621408785835279289708364484253337551960502180649007e-10");
    rs_remainder[18][45]=str_to_Double("-.4776465924689651801753694999540932943676354594855457151159123424e-11");
    rs_remainder[18][46]=str_to_Double(".4332336093775687057524278152423527958055259169244378125865684049e-12");
    rs_remainder[18][47]=str_to_Double(".2184015338760217268197274430519177667395961747359056880787204208e-12");
    rs_remainder[18][48]=str_to_Double("-.1057783826980566304585416506670016291285601174260632447426940000e-13");
    rs_remainder[18][49]=str_to_Double("-.8573860739482407598689458354845961053053818729316348922365381010e-14");
    rs_remainder[18][50]=str_to_Double(".1184986055658689556652661394030109109232247915728744179019803314e-15");
    rs_remainder[18][51]=str_to_Double(".2910799026431033512378761313013665065313965241799485377835377900e-15");
    rs_remainder[18][52]=str_to_Double(".463084651031315395730631622180337834993493235298207315574264467e-17");
    rs_remainder[18][53]=str_to_Double("-.8595007111257437275629456387492532088232472575900648002151701173e-17");
    rs_remainder[18][54]=str_to_Double("-.3685096722377018788843977955261698044092425125392957676004890053e-18");
    rs_remainder[18][55]=str_to_Double(".221631854320781152464555477358789680639764076852442077283564223e-18");
    rs_remainder[18][56]=str_to_Double(".1527640900080784683399463564333523438878123902296867166806853846e-19");
    rs_remainder[18][57]=str_to_Double("-.5000722181606194532865467162371544459127151745881483398768392194e-20");
    rs_remainder[18][58]=str_to_Double("-.4801872360981991991535654168926589656363688728893862320437958727e-21");
    rs_remainder[18][59]=str_to_Double(".9864368553888870397037326339557726472733173273557776373501524732e-22");
    rs_remainder[18][60]=str_to_Double(".1250383778147782885822302839153031432184320733585086908434191428e-22");
    rs_remainder[18][61]=str_to_Double("-.1691811320240345434469365662336614139633494178179400981770153721e-23");
    rs_remainder[18][62]=str_to_Double("-.2798606914553475188185270902532123296082936851128833502123127917e-24");
    rs_remainder[18][63]=str_to_Double(".2486719577457661390257858802053997980597449526376118390739505964e-25");
    rs_remainder[18][64]=str_to_Double(".5491438632704468288608117615856050782099204623075911607878571207e-26");
    rs_remainder[18][65]=str_to_Double("-.3023194654382244973696834198194647559328989643048212310071957386e-27");
    rs_remainder[18][66]=str_to_Double("-.9562312041268983745231256307206796121211862552543333547730491986e-28");
    rs_remainder[18][67]=str_to_Double(".2736545936236159358468072282282618927395311875421082589149988489e-29");
    rs_remainder[18][68]=str_to_Double(".1489678930986206165959335800616893561166424341446568018837634136e-29");
    rs_remainder[18][69]=str_to_Double("-.98836438909808713952325607256088716741888423770439846500331187e-32");
    rs_remainder[18][70]=str_to_Double("-.2087756286635254817372652269690704011848528089808595964774338509e-31");
    rs_remainder[18][71]=str_to_Double("-.2658273563096959947491104823869406527100234435854900896663006059e-33");
    //============================= n = 19 ================================================
    rs_remainder[19][0]=str_to_Double(".893219301398407836136000594145699357440270071363454300500523e-8");
    rs_remainder[19][1]=str_to_Double("-.2352992755280715375309625593904885190664546294853716165924507e-7");
    rs_remainder[19][2]=str_to_Double(".6648769833825603522331639267808394419974874296690907267741750605e-6");
    rs_remainder[19][3]=str_to_Double("-.6100989591484968953283145811430119115963316695617530039551e-5");
    rs_remainder[19][4]=str_to_Double(".3146192923576133694972619314469282657148314967496542505937497891e-4");
    rs_remainder[19][5]=str_to_Double("-.12227193231562109555842090968655776787762340770484708398193e-3");
    rs_remainder[19][6]=str_to_Double(".39913073313964033815167191491711582185449367676308588985967e-3");
    rs_remainder[19][7]=str_to_Double("-.108913596147007744647440566451910114707419799664650430938702e-2");
    rs_remainder[19][8]=str_to_Double(".229190507206248266267109401892472700490689765800952794144383e-2");
    rs_remainder[19][9]=str_to_Double("-.2956789627959746165711910274817945782022487929938380071341e-2");
    rs_remainder[19][10]=str_to_Double("-.80317117239361146862239726497295438826763531557850266536661e-3");
    rs_remainder[19][11]=str_to_Double(".14825819599812946600237404911965805410140191306582443779453377e-1");
    rs_remainder[19][12]=str_to_Double("-.3741177278538565017196545702348355534334098643642112496176095358e-1");
    rs_remainder[19][13]=str_to_Double(".44937284314627207405472919563422840912705249931489231257446e-1");
    rs_remainder[19][14]=str_to_Double("-.1616228249372123532386224243557123642211588484447200489124530789e-2");
    rs_remainder[19][15]=str_to_Double("-.832933447041546491266275207783226450582673267654513853745829e-1");
    rs_remainder[19][16]=str_to_Double(".1162051514498264593565733259183057177765245613681771280096904463");
    rs_remainder[19][17]=str_to_Double("-.173017569412618894412450019815452162044679887223193932158287e-1");
    rs_remainder[19][18]=str_to_Double("-.120420672031887204955919033102486160931340140509777096505314");
    rs_remainder[19][19]=str_to_Double(".109721089901054021234828596996201778090964184586769079999966158");
    rs_remainder[19][20]=str_to_Double(".3604710008891768940976411745710202238014737741458111474894392e-1");
    rs_remainder[19][21]=str_to_Double("-.989706866574561632008432274853468667911774334615294049498620e-1");
    rs_remainder[19][22]=str_to_Double(".1522397841403196845641875441683258233735148378159997817295525e-1");
    rs_remainder[19][23]=str_to_Double(".494849889482267643613804575768397028053339870076078318535236e-1");
    rs_remainder[19][24]=str_to_Double("-.1820375699822167683578191111269691681173883203555714034842129913e-1");
    rs_remainder[19][25]=str_to_Double("-.1698525422739174925405440663211815143416143209134199238007548e-1");
    rs_remainder[19][26]=str_to_Double(".837373665342686419576298890570138869987301621458258655776859272e-2");
    rs_remainder[19][27]=str_to_Double(".444985647167817297344430508028152803598349406876602438091548228e-2");
    rs_remainder[19][28]=str_to_Double("-.2453157348750148449181206266824472704222543462553573124169797504e-2");
    rs_remainder[19][29]=str_to_Double("-.9468601892906575993160114999222675046814402810552940024035305e-3");
    rs_remainder[19][30]=str_to_Double(".5149316259195852753511413097068524780005572043365003298606352395e-3");
    rs_remainder[19][31]=str_to_Double(".1683982108560343778584715249116584530203051410240936215766781843e-3");
    rs_remainder[19][32]=str_to_Double("-.815602898768529173469880774730939064968510970656814822493497607e-4");
    rs_remainder[19][33]=str_to_Double("-.2516249958833280052778555818308503328219615204390026586680006747e-4");
    rs_remainder[19][34]=str_to_Double(".1003236939159816786733475331337081983230679040568943880246262278e-4");
    rs_remainder[19][35]=str_to_Double(".314612482312430510627664913961778474442243614497362883887425479e-5");
    rs_remainder[19][36]=str_to_Double("-.9751142820332806778977589347470042921062179514923282267072132481e-6");
    rs_remainder[19][37]=str_to_Double("-.3281248948923273315895685005171367460468419961551139212521510262e-6");
    rs_remainder[19][38]=str_to_Double(".7558719452191458838430997129895780872616285425693345296287771437e-7");
    rs_remainder[19][39]=str_to_Double(".2857466405967806922858716674995013241341086242318970744846000706e-7");
    rs_remainder[19][40]=str_to_Double("-.4674398100571735190797794459266470970948918024404710185769617186e-8");
    rs_remainder[19][41]=str_to_Double("-.2087188471905070387230263817413049232370110514850883128359075831e-8");
    rs_remainder[19][42]=str_to_Double(".2272257183589369002203801631990418555641701772992621158183723735e-9");
    rs_remainder[19][43]=str_to_Double(".1287007552232827090068421947377851371964952568069874903441326950e-9");
    rs_remainder[19][44]=str_to_Double("-.82339719250182236061479690127334348524891875695256011968905520e-11");
    rs_remainder[19][45]=str_to_Double("-.6747105762434667669028873853152101514571543683175536030532442438e-11");
    rs_remainder[19][46]=str_to_Double(".1789216993437365003775598465516206724306048453152607098742454086e-12");
    rs_remainder[19][47]=str_to_Double(".3028030594287818433570842062532734514866163446439692554891235448e-12");
    rs_remainder[19][48]=str_to_Double(".1713478886983649247516937625573226686645671445739787652098892720e-14");
    rs_remainder[19][49]=str_to_Double("-.1170526225307306080264851362530802079203325647730645138718653293e-13");
    rs_remainder[19][50]=str_to_Double("-.4027970555409631315919744567433939952647586760213896214400271020e-15");
    rs_remainder[19][51]=str_to_Double(".3917148058650304435694162213655391715027139901814853517428984983e-15");
    rs_remainder[19][52]=str_to_Double(".2406176601792288695130376589424858251945062807215112266055766997e-16");
    rs_remainder[19][53]=str_to_Double("-.1138850436107528770139586364275219826021947216749069804230564712e-16");
    rs_remainder[19][54]=str_to_Double("-.1006998280907373396132834770190999470626147365517393302829556922e-17");
    rs_remainder[19][55]=str_to_Double(".288079895693258733230792091043855749368431260963062180306504446e-18");
    rs_remainder[19][56]=str_to_Double(".3384195218497852247284063963624280905524401681761291848309456386e-19");
    rs_remainder[19][57]=str_to_Double("-.6329822533294203615511978156437846519470118129313364388932662002e-20");
    rs_remainder[19][58]=str_to_Double("-.9598017053616922822457293793812784702449602207131537229860649059e-21");
    rs_remainder[19][59]=str_to_Double(".1199620556150856269861670910858657189622349187371820550035630648e-21");
    rs_remainder[19][60]=str_to_Double(".2355454631405162555447897494085756105695351911371573556959927653e-22");
    rs_remainder[19][61]=str_to_Double("-.192596839359071217531885456976775912998223861156661678392539399e-23");
    rs_remainder[19][62]=str_to_Double("-.5077574912818985789017902415933429176189235922539706468094588709e-24");
    rs_remainder[19][63]=str_to_Double(".2501625241419679889228086600796526049304880677540703402072487839e-25");
    rs_remainder[19][64]=str_to_Double(".9710905319576994603154409595871330053205210167241365056377985173e-26");
    rs_remainder[19][65]=str_to_Double("-.2258799008287111966221925865899617489480655127603954165639713364e-27");
    rs_remainder[19][66]=str_to_Double("-.1659281028906037107686636422689673392627956167718637722824943188e-27");
    rs_remainder[19][67]=str_to_Double(".2069118901207932426807165268808470041978906672938810378255098985e-30");
    rs_remainder[19][68]=str_to_Double(".2545509089536622051266256045145645193310269202035005936336625667e-29");
    rs_remainder[19][69]=str_to_Double(".469268543116234001747909829544313304796934812752799593282677215e-31");
    rs_remainder[19][70]=str_to_Double("-.3517482815042533087348690790241659826939486063386451317771097369e-31");
    rs_remainder[19][71]=str_to_Double("-.1303808835340158797237268684506669190657297386188059143535716673e-32");
}

void initialize_rs_remainder5(){

    //============================= n = 20 ================================================
    rs_remainder[20][0]=str_to_Double("-.559302842663312450752457473747594409990839151017179012967210e-8");
    rs_remainder[20][1]=str_to_Double("-.12001667116897246037574423172226482013582892702056476542745e-7");
    rs_remainder[20][2]=str_to_Double(".2463276787487751715805293747649575130108965345234973954295553232e-6");
    rs_remainder[20][3]=str_to_Double("-.153008043787720061856987356778766149480631038320467060041820e-6");
    rs_remainder[20][4]=str_to_Double("-.416342558387336146913664128919706344172832431285017822742e-5");
    rs_remainder[20][5]=str_to_Double(".25723399331161154622359070922070610264683315165563382717410e-4");
    rs_remainder[20][6]=str_to_Double("-.9736319936396126111211248318394378050036841213099636354286e-4");
    rs_remainder[20][7]=str_to_Double(".2688570840606231347814122711023983976921156263034735441094e-3");
    rs_remainder[20][8]=str_to_Double("-.499797514326645662307957121722369871706323263302364080228e-3");
    rs_remainder[20][9]=str_to_Double(".26403947894794554153191845157236153393095519613019757266524828e-3");
    rs_remainder[20][10]=str_to_Double(".2124268070600467170258932019680955721118841609890905072185454685e-2");
    rs_remainder[20][11]=str_to_Double("-.9321830671098570223614704234209813507876793962118872077909e-2");
    rs_remainder[20][12]=str_to_Double(".21181126771688620949354418138100085097737259559641589718144e-1");
    rs_remainder[20][13]=str_to_Double("-.2717262083989055815702647865135199636551996805885749258186e-1");
    rs_remainder[20][14]=str_to_Double(".621865842305152018436275814230075721102914477330276578982e-2");
    rs_remainder[20][15]=str_to_Double(".4807518662707914244836036715401955668522358447613733973592e-1");
    rs_remainder[20][16]=str_to_Double("-.91293456235616242884855301244949703974987175081010097680244e-1");
    rs_remainder[20][17]=str_to_Double(".5133645453722076779730386072859425919918217307404061094311118284e-1");
    rs_remainder[20][18]=str_to_Double(".6384906418065581826516806999617155084007358634888228932931462298e-1");
    rs_remainder[20][19]=str_to_Double("-.122696797648925451415675097245367984913190693311563005109341");
    rs_remainder[20][20]=str_to_Double(".354449774402490964442593716948083489431388456432205377961723e-1");
    rs_remainder[20][21]=str_to_Double(".8271862902750968877075915516271564933115960316255405306805158223e-1");
    rs_remainder[20][22]=str_to_Double("-.7080807432060148872593085441464862131144253413140415907901863122e-1");
    rs_remainder[20][23]=str_to_Double("-.2304311484552353099788484694567921337910676001833098481789835941e-1");
    rs_remainder[20][24]=str_to_Double(".4705371488748176183519854635511774185018353053615904016564952452e-1");
    rs_remainder[20][25]=str_to_Double("-.1095595042228500943440156372421498992370786972871274669710180161e-2");
    rs_remainder[20][26]=str_to_Double("-.1930162807064304106893106336922344555927795416982772490620151108e-1");
    rs_remainder[20][27]=str_to_Double(".3200644516180608309981732428699851817475666637632741362200953e-2");
    rs_remainder[20][28]=str_to_Double(".56942986967043530247104573012904422267696208407876050765424053e-2");
    rs_remainder[20][29]=str_to_Double("-.12949943535577906457080876522768355148183182086103598857114111e-2");
    rs_remainder[20][30]=str_to_Double("-.12975977232970669694547809089129441792984677948043791632565493e-2");
    rs_remainder[20][31]=str_to_Double(".3124944359778612269597137118807554074712442666367065077272635224e-3");
    rs_remainder[20][32]=str_to_Double(".23706977947092637222450730806865135114857822300525582604515162e-3");
    rs_remainder[20][33]=str_to_Double("-.5290698939702912104708180416137269928079663209182838509811096126e-4");
    rs_remainder[20][34]=str_to_Double("-.3538784489301323493324358149027231130801536802998077084949777737e-4");
    rs_remainder[20][35]=str_to_Double(".6655957342243245089585791117593305991420682450798527337398486e-5");
    rs_remainder[20][36]=str_to_Double(".435854229012926449489706633609982803589202146551512246801226329e-5");
    rs_remainder[20][37]=str_to_Double("-.636563449300243775759985186941532772962949357345028125314774611e-6");
    rs_remainder[20][38]=str_to_Double("-.4458324200375894998381027968253973016789602259620243093893906490e-6");
    rs_remainder[20][39]=str_to_Double(".46215409076922612175283753339490399215723463157746827372508833e-7");
    rs_remainder[20][40]=str_to_Double(".381025456476420257615565759385081490591035771271037057696748475e-7");
    rs_remainder[20][41]=str_to_Double("-.2444199599041122599693563944157882422491880190336356996123127114e-8");
    rs_remainder[20][42]=str_to_Double("-.2737963349371651711552213391619198140166373923241790689495167225e-8");
    rs_remainder[20][43]=str_to_Double(".7805498355972045218119587642574879029059103911218413227466489090e-10");
    rs_remainder[20][44]=str_to_Double(".1665018348331920787912382595185826828500715037537653614337278068e-9");
    rs_remainder[20][45]=str_to_Double(".5792976319943755394833096386019811651948139255537610772708862606e-12");
    rs_remainder[20][46]=str_to_Double("-.8623303100688499503535309682909941944461354792052189195287352580e-11");
    rs_remainder[20][47]=str_to_Double("-.2813278607935765193938476650550979836841051675748389922242157089e-12");
    rs_remainder[20][48]=str_to_Double(".382536858755860871782254192024615471127597447725171774116566971e-12");
    rs_remainder[20][49]=str_to_Double(".2296090957897114158446056569392659774139431337810941186619297697e-13");
    rs_remainder[20][50]=str_to_Double("-.1460341899598904883697953591750178133195797713004630773376313187e-13");
    rs_remainder[20][51]=str_to_Double("-.1270357376212275819549572937974599424799822557343743515726355890e-14");
    rs_remainder[20][52]=str_to_Double(".481277371827133233952418761207284941297548961535753287699258468e-15");
    rs_remainder[20][53]=str_to_Double(".5543181310101470886247763386865725237372356166999921377935717687e-16");
    rs_remainder[20][54]=str_to_Double("-.1370542956252104879160135478383034584212169592686339061237238901e-16");
    rs_remainder[20][55]=str_to_Double("-.2013454999133333843978365222831218634630272352613975043725042356e-17");
    rs_remainder[20][56]=str_to_Double(".3363269649302911704580337679742249466730234210835071130381455484e-18");
    rs_remainder[20][57]=str_to_Double(".6256214968664215780978719925984976303454194929873443265275458533e-19");
    rs_remainder[20][58]=str_to_Double("-.7045758211618840007411681957980479965726535272132048987740703931e-20");
    rs_remainder[20][59]=str_to_Double("-.1690658059505709181240185070387164869867450377068249406154007373e-20");
    rs_remainder[20][60]=str_to_Double(".1230336056129295197193061446598605217176443342602571893853706005e-21");
    rs_remainder[20][61]=str_to_Double(".4018300199495667110610667095733778125870913868384716584772341992e-22");
    rs_remainder[20][62]=str_to_Double("-.167843226719880182611527381668236139363976426663073501947866323e-23");
    rs_remainder[20][63]=str_to_Double("-.8467986405479907550034493354198588920548551604910496726957805130e-24");
    rs_remainder[20][64]=str_to_Double(".1382724473036608320847827055558452367916917101532766606757801609e-25");
    rs_remainder[20][65]=str_to_Double(".1591701246510204919020327154424765222556984095343597519340616184e-25");
    rs_remainder[20][66]=str_to_Double(".8973056769669143769840084663105365162457600555192032460077145795e-28");
    rs_remainder[20][67]=str_to_Double("-.2680212871739970151018363707388520342675953903742659871813148210e-27");
    rs_remainder[20][68]=str_to_Double("-.6865875749719477717795630837358422984085091967085229759098589696e-29");
    rs_remainder[20][69]=str_to_Double(".4054588195699398418537656534696950654521864107222740337357040886e-29");
    rs_remainder[20][70]=str_to_Double(".1817883171681163721773076349652627658341737372570910257135934376e-30");
    rs_remainder[20][71]=str_to_Double("-.5517713378067080795248255916294684628179192415895362614978836824e-31");
    //============================= n = 21 ================================================
    rs_remainder[21][0]=str_to_Double(".5214988790425134738396328983589494860150043784157477080235683e-8");
    rs_remainder[21][1]=str_to_Double("-.56807547663837929463483045372438429573201350527872060596905e-7");
    rs_remainder[21][2]=str_to_Double(".3621940790907443402866181596674591295807166639776396817135e-6");
    rs_remainder[21][3]=str_to_Double("-.1425165181281159612671644870716295517307991480699532226088e-5");
    rs_remainder[21][4]=str_to_Double(".359742435555611466099331272371999698871569676755801826195027088e-5");
    rs_remainder[21][5]=str_to_Double("-.6050669926069108415374895410320694742588501955175946966148582e-5");
    rs_remainder[21][6]=str_to_Double(".14074019024125357534519595428763762610436962954973579158882633e-4");
    rs_remainder[21][7]=str_to_Double("-.880675400176230882070232219092724565808651134516133144417e-4");
    rs_remainder[21][8]=str_to_Double(".494563391313192411774131204091628096233027524093820331323124e-3");
    rs_remainder[21][9]=str_to_Double("-.195210389489795148254967477140400673892387373054437967725016e-2");
    rs_remainder[21][10]=str_to_Double(".5559398253375692703496408598600247682909968237497740617308328e-2");
    rs_remainder[21][11]=str_to_Double("-.112793506879329337687389445729476177885257727921017115000238e-1");
    rs_remainder[21][12]=str_to_Double(".1437463052895683433150072641284581254943534801270776464426e-1");
    rs_remainder[21][13]=str_to_Double("-.3723464502586179066226409077966146085654693015379228412049e-2");
    rs_remainder[21][14]=str_to_Double("-.283227816284062422685143250414668666002958802164933902317622e-1");
    rs_remainder[21][15]=str_to_Double(".637889633964918565707535438336706795479674919704658751259708e-1");
    rs_remainder[21][16]=str_to_Double("-.554615149155916765105625113484404089811467699340481373596211e-1");
    rs_remainder[21][17]=str_to_Double("-.20382468559977280477758287366629059065222607775774831692352643e-1");
    rs_remainder[21][18]=str_to_Double(".99359048472316694129876606615270060001260783084260667045585e-1");
    rs_remainder[21][19]=str_to_Double("-.7796390039178850241710217733393648786868048899778210974777780641e-1");
    rs_remainder[21][20]=str_to_Double("-.3299865579309309212195191030824694836807108742563612654085960034e-1");
    rs_remainder[21][21]=str_to_Double(".9091585107370445280622358892818525364109030637981899287803307964e-1");
    rs_remainder[21][22]=str_to_Double("-.2771718153998941872116106923341779648993102239760894887347231294e-1");
    rs_remainder[21][23]=str_to_Double("-.4624927164135199210817951660962964870536250585394585209944743316e-1");
    rs_remainder[21][24]=str_to_Double(".32837756711308315198632479690605102498471918007148330055464612e-1");
    rs_remainder[21][25]=str_to_Double(".130282571033176932167015127027517299937153903469164384196794671e-1");
    rs_remainder[21][26]=str_to_Double("-.1688773975965796289341866765880925368514215424123816925900925e-1");
    rs_remainder[21][27]=str_to_Double("-.19004879965399284406105546188915047567665985545598434670953365e-2");
    rs_remainder[21][28]=str_to_Double(".5677426780469825894984862649068927948516277801453887984826098e-2");
    rs_remainder[21][29]=str_to_Double(".6668798429212052020112832079750701654724732879439718593864e-6");
    rs_remainder[21][30]=str_to_Double("-.14021936302070796472639988483600847986591602922793529704353195e-2");
    rs_remainder[21][31]=str_to_Double(".63055661613933278367668903870696353719551732048159013746846692e-4");
    rs_remainder[21][32]=str_to_Double(".2686399080984684191529337882480254842604896994019152045359931e-3");
    rs_remainder[21][33]=str_to_Double("-.1484577726836869775125598537285189013602674367894334534510253925e-4");
    rs_remainder[21][34]=str_to_Double("-.4113202033591240209330883599270254776690396196711918348411571438e-4");
    rs_remainder[21][35]=str_to_Double(".1907918662026141649971133876993471617807317219631367523770615007e-5");
    rs_remainder[21][36]=str_to_Double(".5124895898641966400802300512136908843488547682842604499950580135e-5");
    rs_remainder[21][37]=str_to_Double("-.1409219987375291111330058045917970425279829105535522398021936139e-6");
    rs_remainder[21][38]=str_to_Double("-.5260846947166151442823803019668694244493607550142809888385128383e-6");
    rs_remainder[21][39]=str_to_Double(".2140626226029965056882424378236729269498764701670827093670885553e-8");
    rs_remainder[21][40]=str_to_Double(".4492272719344014320232504741632308977970800108901887638396049176e-7");
    rs_remainder[21][41]=str_to_Double(".9392974701501557470813314439933604723915682893194575865582698e-9");
    rs_remainder[21][42]=str_to_Double("-.3217054809079186391518359183735904995685232267826704577846907376e-8");
    rs_remainder[21][43]=str_to_Double("-.14829074452136002072117171067514389115942019794406817750362330e-9");
    rs_remainder[21][44]=str_to_Double(".194592661447236135379622033165405979832105671579088059984675353e-9");
    rs_remainder[21][45]=str_to_Double(".13856515100379785296670417589030887721046664325116353653649005e-10");
    rs_remainder[21][46]=str_to_Double("-.1000251187918157130672553954492933381597218213578180851572624287e-10");
    rs_remainder[21][47]=str_to_Double("-.9669731072006973832022779713469517025339407817181902592839328093e-12");
    rs_remainder[21][48]=str_to_Double(".4390122210153106605863354636346333243797967731848422899874266793e-12");
    rs_remainder[21][49]=str_to_Double(".5421613324543399507209385926215184133598933789595045209414924565e-13");
    rs_remainder[21][50]=str_to_Double("-.1650110452376429356270030479266574800309508027826308291757007434e-13");
    rs_remainder[21][51]=str_to_Double("-.2530771435568298769063531193640253027491785681737903017117226876e-14");
    rs_remainder[21][52]=str_to_Double(".5312727437905039645213598328764152796153771717605102540687768839e-15");
    rs_remainder[21][53]=str_to_Double(".1004887996555911234918015781718507489880009935473749419900127655e-15");
    rs_remainder[21][54]=str_to_Double("-.1458807944537135516258189787223571828041403768129790282939798545e-16");
    rs_remainder[21][55]=str_to_Double("-.3444287390262935749176898979720668667108309588552198802718509471e-17");
    rs_remainder[21][56]=str_to_Double(".337124407311232805927929840137169291448044525446379787938468293e-18");
    rs_remainder[21][57]=str_to_Double(".1030137407634612035763006392622085261022373544871195040955129494e-18");
    rs_remainder[21][58]=str_to_Double("-.633449832500145028788148660691097788626465647322601181702782194e-20");
    rs_remainder[21][59]=str_to_Double("-.2710928477046665651052271832665218512478887510200227540380392202e-20");
    rs_remainder[21][60]=str_to_Double(".8717311392996137407106777968797967875305364632788128941938665869e-22");
    rs_remainder[21][61]=str_to_Double(".6318308263095272333178445676495985737287150524555344310421738463e-22");
    rs_remainder[21][62]=str_to_Double("-.4687024410008022443735822757563284289583580549496886236377044273e-24");
    rs_remainder[21][63]=str_to_Double("-.1310845783899113846415503938423741710596834519197821767223139820e-23");
    rs_remainder[21][64]=str_to_Double("-.189115985906916152715098623921704736539691183315556622145523434e-25");
    rs_remainder[21][65]=str_to_Double(".2430071248373067041617569798798008437537156436754271117199342020e-25");
    rs_remainder[21][66]=str_to_Double(".8431972603337747377108923034899697586393812134271801496768295613e-27");
    rs_remainder[21][67]=str_to_Double("-.4035367142072070069027577165896507281396832753167620227931533183e-27");
    rs_remainder[21][68]=str_to_Double("-.2202306133595501851548417530506666531833361435126036679249756405e-28");
    rs_remainder[21][69]=str_to_Double(".6008434738539668762335196699515232404917230860670090558421387173e-29");
    rs_remainder[21][70]=str_to_Double(".4525566763918074161185803374594291433995430231069661415070733191e-30");
    rs_remainder[21][71]=str_to_Double("-.8012760131124511637940906635656260494312730294438760866965206398e-31");
    //============================= n = 22 ================================================
    rs_remainder[22][0]=str_to_Double("-.309997000615056935668075726553925618062997167979569089237068485e-8");
    rs_remainder[22][1]=str_to_Double(".1843928111707119917369347502576927610267927663889550384782e-8");
    rs_remainder[22][2]=str_to_Double(".464156716380096037490959580079233905107020598078042452746e-7");
    rs_remainder[22][3]=str_to_Double(".316396654828262253218373637458036008281372740360240131949e-6");
    rs_remainder[22][4]=str_to_Double("-.315039343816098419654711671405949188519611527702721176921e-5");
    rs_remainder[22][5]=str_to_Double(".136847039897873591120673766752574654255222109893168523275e-4");
    rs_remainder[22][6]=str_to_Double("-.4592953336320214165115287773075810795562794630289169849303e-4");
    rs_remainder[22][7]=str_to_Double(".14411201103208211600190464452109472816586340133394553571851e-3");
    rs_remainder[22][8]=str_to_Double("-.437323319664516814080437923889408698812330185167518800261e-3");
    rs_remainder[22][9]=str_to_Double(".1217635326965792199401145635217069378936444152768226586403e-2");
    rs_remainder[22][10]=str_to_Double("-.28860995603967786685634186822030906244163257231558910571621652e-2");
    rs_remainder[22][11]=str_to_Double(".532123503523754742879809636871814385777221300490099643817e-2");
    rs_remainder[22][12]=str_to_Double("-.629762425251891344130755644641840006230810205808034186821e-2");
    rs_remainder[22][13]=str_to_Double(".269846415545838920358693452357052454627416919298857211339215e-3");
    rs_remainder[22][14]=str_to_Double(".1802403641470595421043894442881904267341926103499872112783e-1");
    rs_remainder[22][15]=str_to_Double("-.4211498921546480925129116538165377870229135664239336821591155740e-1");
    rs_remainder[22][16]=str_to_Double(".455704715180513531929298635613763048182720183292382420211e-1");
    rs_remainder[22][17]=str_to_Double("-.218413386530613223977075846361909877632926979217901107216e-2");
    rs_remainder[22][18]=str_to_Double("-.66958738869935741850231112127810652025093812342561840668364614e-1");
    rs_remainder[22][19]=str_to_Double(".8648141262256154466531746456284849253481747904755356030859e-1");
    rs_remainder[22][20]=str_to_Double("-.146228934673930730767712447052156793714198228204680985432114e-1");
    rs_remainder[22][21]=str_to_Double("-.7323390524548785696186995195033622367968070625403427801187406329e-1");
    rs_remainder[22][22]=str_to_Double(".67266602633867151988901600280962299426749970572683863160818e-1");
    rs_remainder[22][23]=str_to_Double(".13774572070622475247878554472292728211626186157531548910183e-1");
    rs_remainder[22][24]=str_to_Double("-.5110930806701509321817953521657264896679806209306744644619264650e-1");
    rs_remainder[22][25]=str_to_Double(".12968113826030021424505759096418750023827845998910761746061e-1");
    rs_remainder[22][26]=str_to_Double(".2120632680540720995889504117182062509383574065531286813224689232e-1");
    rs_remainder[22][27]=str_to_Double("-.110554349446851778437020721089427424607371802718273251908194e-1");
    rs_remainder[22][28]=str_to_Double("-.5815066347805108962245841892249882245237267279567931008668231798e-2");
    rs_remainder[22][29]=str_to_Double(".453750839077071143253851904261243379828418877762849169476583e-2");
    rs_remainder[22][30]=str_to_Double(".1163717006161331430089993879539645326920083183477340741950075e-2");
    rs_remainder[22][31]=str_to_Double("-.1252843414967218215119183092671953626630819923603469074570626e-2");
    rs_remainder[22][32]=str_to_Double("-.1860259034905761030006485508500163722025028917669224815343917560e-3");
    rs_remainder[22][33]=str_to_Double(".257052124092013102932675028690422152255182926529902154790514184e-3");
    rs_remainder[22][34]=str_to_Double(".2626822972139994091816376617616377792867951893611598745709027097e-4");
    rs_remainder[22][35]=str_to_Double("-.4108602385056982817913913354387323924290295053283238285384297254e-4");
    rs_remainder[22][36]=str_to_Double("-.3536577622223708553402910351512665434265364946440284089721925846e-5");
    rs_remainder[22][37]=str_to_Double(".5257364359121964404192823239198024829788449395360239885356179778e-5");
    rs_remainder[22][38]=str_to_Double(".4537265415399493354214680739375397760407014950810142617937980e-6");
    rs_remainder[22][39]=str_to_Double("-.5482936077899370024715259399052582733748740833239195348129220563e-6");
    rs_remainder[22][40]=str_to_Double("-.5252785622473343281277456165504759388115956327702136481692611e-7");
    rs_remainder[22][41]=str_to_Double(".4720917668241875550226048534488811171901704937004723327060024035e-7");
    rs_remainder[22][42]=str_to_Double(".5234553974960330226442780214031305084433932353650122255694355207e-8");
    rs_remainder[22][43]=str_to_Double("-.3389109241035626356514809843642094893054534121096036521366214215e-8");
    rs_remainder[22][44]=str_to_Double("-.4399520305827151149289937751404796452270507054726659546985446996e-9");
    rs_remainder[22][45]=str_to_Double(".2044059135683583754891550743920849129263986761875552224884749368e-9");
    rs_remainder[22][46]=str_to_Double(".3111120450494229539007045036276054989745545892401650242420060026e-10");
    rs_remainder[22][47]=str_to_Double("-.1041431720994465211166364508610639839299692601603604107651723125e-10");
    rs_remainder[22][48]=str_to_Double("-.1861452485580167755500591805973522004202008233568789283953373987e-11");
    rs_remainder[22][49]=str_to_Double(".44954668045018439813404681850237338009239715777134548398716696e-12");
    rs_remainder[22][50]=str_to_Double(".9501870541230635973181085948020038947225348902139423000540647392e-13");
    rs_remainder[22][51]=str_to_Double("-.1642843930469619261834117030859916295522040304545263717999674635e-13");
    rs_remainder[22][52]=str_to_Double("-.417387291498550803623790902248234616799625576847295815546545590e-14");
    rs_remainder[22][53]=str_to_Double(".50463451181592050200982555801086562453022489500897624573251703e-15");
    rs_remainder[22][54]=str_to_Double(".1590578660021272310897909216541463999559207043956681452230815325e-15");
    rs_remainder[22][55]=str_to_Double("-.1276236559934777716525016652339707541919936892590083596455751825e-16");
    rs_remainder[22][56]=str_to_Double("-.5296608503845608549235513324316292217325792814314705524229268830e-17");
    rs_remainder[22][57]=str_to_Double(".2508907075247513780447833437608151765077537236183242483527736544e-18");
    rs_remainder[22][58]=str_to_Double(".1550959631982582332604733374302284944825365059637960815856890552e-18");
    rs_remainder[22][59]=str_to_Double("-.3080565980659035864321201462219732900032811601027534278978262811e-20");
    rs_remainder[22][60]=str_to_Double("-.4014965536035691409968022870476421771483125476811695410798882148e-20");
    rs_remainder[22][61]=str_to_Double("-.1621411565537426070798447908524445082023975249303702339903245262e-22");
    rs_remainder[22][62]=str_to_Double(".9228646102599192504952156276248507261245380220411083042154233634e-22");
    rs_remainder[22][63]=str_to_Double(".237130224241252518855826145607504403025334663412573919633793284e-23");
    rs_remainder[22][64]=str_to_Double("-.1889787379178682682529686027610478741106206698830005291928089259e-23");
    rs_remainder[22][65]=str_to_Double("-.874474402858401733889465157189141251387180236724902372116620088e-25");
    rs_remainder[22][66]=str_to_Double(".345483739519389304044865565196178562208775384530605912053310375e-25");
    rs_remainder[22][67]=str_to_Double(".2312152100377092640681237079027971521855627404715258124929771241e-26");
    rs_remainder[22][68]=str_to_Double("-.5642016398113051422818842630355387857047671333795940973509798654e-27");
    rs_remainder[22][69]=str_to_Double("-.5021273934855637490973556909145204871717553906436501447361275024e-28");
    rs_remainder[22][70]=str_to_Double(".8217259505012851540295499716276774611551536609580526211680703118e-29");
    rs_remainder[22][71]=str_to_Double(".9399443629796439004170523419198730318625903648751605706132251884e-30");
    //============================= n = 23 ================================================
    rs_remainder[23][0]=str_to_Double(".249124013604233893726847799293720371732178871236422792795698485e-8");
    rs_remainder[23][1]=str_to_Double("-.3963765519245909635656528671585905944142855298579106129674610535e-7");
    rs_remainder[23][2]=str_to_Double(".119864967447921451139223261169227246905833511984517170341218e-6");
    rs_remainder[23][3]=str_to_Double(".2685603376434480145920831169993053213391008505550927119040e-6");
    rs_remainder[23][4]=str_to_Double("-.347640492470504451212839604734251874747048440319487409415e-5");
    rs_remainder[23][5]=str_to_Double(".17283355899821056242313684497599366897361274529730543309e-4");
    rs_remainder[23][6]=str_to_Double("-.64153982178704594095135797884777307356098753364955460664e-4");
    rs_remainder[23][7]=str_to_Double(".19918011436144857520721362925028287538693363112128086818209e-3");
    rs_remainder[23][8]=str_to_Double("-.528260926876479318056540473335990492964492996075325344527e-3");
    rs_remainder[23][9]=str_to_Double(".11656650816618410574191674172810779246051828319575840881108e-2");
    rs_remainder[23][10]=str_to_Double("-.1959612602644998761438888294866843473177326394894785273905e-2");
    rs_remainder[23][11]=str_to_Double(".18137380881566997337499107386274388327516510988171916547168997e-2");
    rs_remainder[23][12]=str_to_Double(".19500405509147916978536774280007130191937353347247704344034e-2");
    rs_remainder[23][13]=str_to_Double("-.1228105528343310548930970612891060973068089368659497386018e-1");
    rs_remainder[23][14]=str_to_Double(".2696899530583600735145736462717835369463389666070443299251600e-1");
    rs_remainder[23][15]=str_to_Double("-.323633672097585078532714169061227836217862238388132998115723997e-1");
    rs_remainder[23][16]=str_to_Double(".92154344791551651877077764800530839565694809276101859042e-2");
    rs_remainder[23][17]=str_to_Double(".4091595692844045478259028299259965353757221208510496451162e-1");
    rs_remainder[23][18]=str_to_Double("-.7481187490347544886043704056736859533055830833060706879704041e-1");
    rs_remainder[23][19]=str_to_Double(".4201504689135814765776002365892272415074475656103969842874e-1");
    rs_remainder[23][20]=str_to_Double(".39300018675334510205459138869522130149066755181374709065163306e-1");
    rs_remainder[23][21]=str_to_Double("-.7815409896181117152560167082068564701507223673214862053543510e-1");
    rs_remainder[23][22]=str_to_Double(".270796258963986222505701833033478732758271653417514208420912e-1");
    rs_remainder[23][23]=str_to_Double(".4141798670447060948584721870458883411792005585864430748745662327e-1");
    rs_remainder[23][24]=str_to_Double("-.4027224215881692176766791308578271558271872159263548340646226e-1");
    rs_remainder[23][25]=str_to_Double("-.6112565300369720783414850927944906416156183493880882462428652614e-2");
    rs_remainder[23][26]=str_to_Double(".2251314441332011247944543156153481214580611559028345271807735836e-1");
    rs_remainder[23][27]=str_to_Double("-.374056351451428106765925962130667477027560719535053581339354e-2");
    rs_remainder[23][28]=str_to_Double("-.784852846314578406112596185582385132214730190889430075300823e-2");
    rs_remainder[23][29]=str_to_Double(".2647328420755578268875081722790163211901937959119220973210753e-2");
    rs_remainder[23][30]=str_to_Double(".197422945006136882262477859053843353063547800840361779442885510e-2");
    rs_remainder[23][31]=str_to_Double("-.8948780955040608042203107381847142949620691696247742794710045197e-3");
    rs_remainder[23][32]=str_to_Double("-.388823915670746784591830342487881127441072757207104782672911500e-3");
    rs_remainder[23][33]=str_to_Double(".2047854508566268202107764984816483300070589904118885142431913452e-3");
    rs_remainder[23][34]=str_to_Double(".633379602091950994220176047423107522741586085333944556084601e-4");
    rs_remainder[23][35]=str_to_Double("-.3497742513705118605906316378824437184942627629416892171884088331e-4");
    rs_remainder[23][36]=str_to_Double("-.88199645438288212589917743680252982923399254397865021122802297e-5");
    rs_remainder[23][37]=str_to_Double(".4666444232229650318076695509675709239799246900919244583386286139e-5");
    rs_remainder[23][38]=str_to_Double(".1063554781863184674574366299094318258587308018745761627817603871e-5");
    rs_remainder[23][39]=str_to_Double("-.4992645845534463852682470666789312355576795265091248425173888534e-6");
    rs_remainder[23][40]=str_to_Double("-.1109994320528679140222865802944382243813249734502735512494428142e-6");
    rs_remainder[23][41]=str_to_Double(".4357583907317462169413108196574689729196596617060438635381787457e-7");
    rs_remainder[23][42]=str_to_Double(".9977049912014657834194016402625902822457876764325639735144140028e-8");
    rs_remainder[23][43]=str_to_Double("-.3138840281088605931001024666683717988227623591870193256397458213e-8");
    rs_remainder[23][44]=str_to_Double("-.7696686788513788484820202160867107001126273971591625583706273252e-9");
    rs_remainder[23][45]=str_to_Double(".1879904387544788297441173712505379842052639924802275914196498662e-9");
    rs_remainder[23][46]=str_to_Double(".509513897253534504332149580736605500541234106596136571774536099e-10");
    rs_remainder[23][47]=str_to_Double("-.9391780426469726938299988915507305240022404646295524149282766217e-11");
    rs_remainder[23][48]=str_to_Double("-.2902076228084179304369017936346118844026032057668236385944244127e-11");
    rs_remainder[23][49]=str_to_Double(".3903741365660486626832628740319175250490899001497094008127663075e-12");
    rs_remainder[23][50]=str_to_Double(".1428441344938101001776205606409710024331974742881902498632509349e-12");
    rs_remainder[23][51]=str_to_Double("-.1332344667467724071227840254001561212135689474983908681254099966e-13");
    rs_remainder[23][52]=str_to_Double("-.6107458468921150481361991896002771827224800506392058843155730715e-14");
    rs_remainder[23][53]=str_to_Double(".3592915506375497102773978840426967340083024309718374959056764507e-15");
    rs_remainder[23][54]=str_to_Double(".2280403260263845622052894074363625999352311085273714673909192627e-15");
    rs_remainder[23][55]=str_to_Double("-.6739402173290784625567019888221076033559053984676322527583460182e-17");
    rs_remainder[23][56]=str_to_Double("-.7473323855559305831867062458961057047001556109852696151569244367e-17");
    rs_remainder[23][57]=str_to_Double(".3027774188812580797039606869681330517054349591344742813373889819e-19");
    rs_remainder[23][58]=str_to_Double(".2159445468678434506249856217550286101579549226103405527763208848e-18");
    rs_remainder[23][59]=str_to_Double(".40752029412517590284521927482962283341587730791628638002877218e-20");
    rs_remainder[23][60]=str_to_Double("-.5522879111972254720000444741914890590565198204525562541193757880e-20");
    rs_remainder[23][61]=str_to_Double("-.2224674264148637591689548998909591400348730357328187502825757045e-21");
    rs_remainder[23][62]=str_to_Double(".1253892841999962125713106066656605756642989322564827057601965267e-21");
    rs_remainder[23][63]=str_to_Double(".7675993049361667749553386309411398554239409730532344822084246994e-23");
    rs_remainder[23][64]=str_to_Double("-.253159530888647797411915869498940701582018857272786448871173347e-23");
    rs_remainder[23][65]=str_to_Double("-.2097056419066623303406298613906788060447549610734736070195714927e-24");
    rs_remainder[23][66]=str_to_Double(".454599027376049000306995262908429943053055401731040682977771655e-25");
    rs_remainder[23][67]=str_to_Double(".4847451797006564272916999513651727355436149999840978706204099776e-26");
    rs_remainder[23][68]=str_to_Double("-.7243534285004073612188205799346368358415886495993183468847657017e-27");
    rs_remainder[23][69]=str_to_Double("-.9769798841478138512745204082192766433924192712662766827196663428e-28");
    rs_remainder[23][70]=str_to_Double(".1017455033293290481712107837104802540102158695596985984315140847e-28");
    rs_remainder[23][71]=str_to_Double(".1745971187356643764352039186215056691402730788219422767025471416e-29");
    //============================= n = 24 ================================================
    rs_remainder[24][0]=str_to_Double("-.18485708755558750754148376868272284972413561032781583942764e-8");
    rs_remainder[24][1]=str_to_Double(".666232209499619835906562702946817737892876872442961116828194e-8");
    rs_remainder[24][2]=str_to_Double("-.24879395218847093413888909096533362795937419818085311917e-8");
    rs_remainder[24][3]=str_to_Double(".954443242977865326958005908554366641263065833403037379384e-7");
    rs_remainder[24][4]=str_to_Double("-.40806819988813343605636239279398093746546229839573315217e-6");
    rs_remainder[24][5]=str_to_Double("-.31500227816251957818634152763650776040495698063168151411e-6");
    rs_remainder[24][6]=str_to_Double(".8381490755242933873354701900080814358233238661584115789e-5");
    rs_remainder[24][7]=str_to_Double("-.40199459426900707863955655663796211532338648110862016913e-4");
    rs_remainder[24][8]=str_to_Double(".122765286005918777996138408204171724116575167883513281779e-3");
    rs_remainder[24][9]=str_to_Double("-.26516136653951939218242903651146819183675929082803630962e-3");
    rs_remainder[24][10]=str_to_Double(".336560634946909235550606119153274956494753418875251983585575634e-3");
    rs_remainder[24][11]=str_to_Double(".21701117639475019922994291396235718104840007314684723751e-3");
    rs_remainder[24][12]=str_to_Double("-.265403701115887776084975742426924698438542311239226211662e-2");
    rs_remainder[24][13]=str_to_Double(".842550468518621952908315916497098268671609604481236723354e-2");
    rs_remainder[24][14]=str_to_Double("-.1669283954747522305345196742180612786405423246371931387825e-1");
    rs_remainder[24][15]=str_to_Double(".206165808135091894016233335663531092307809414236982594608722e-1");
    rs_remainder[24][16]=str_to_Double("-.82444178975194279919186593778084053364131759856361237129e-2");
    rs_remainder[24][17]=str_to_Double("-.246929622110421614088470102938041207764868280606935780765e-1");
    rs_remainder[24][18]=str_to_Double(".57130919539022243166520059731906678447073307533724334775633e-1");
    rs_remainder[24][19]=str_to_Double("-.4954819728840221414863720397378342134617733612019362571711950799e-1");
    rs_remainder[24][20]=str_to_Double("-.91717131211760527884368806048480154602000198100827269594e-2");
    rs_remainder[24][21]=str_to_Double(".6616800428978750200879997864273919412251378998978972518369820495e-1");
    rs_remainder[24][22]=str_to_Double("-.5451083409633044793715035486580716800859012949665112897226e-1");
    rs_remainder[24][23]=str_to_Double("-.12666019358301822981732664792716411394578576702526900581252e-1");
    rs_remainder[24][24]=str_to_Double(".499758260399684387764533085787310986029155275922362549588090e-1");
    rs_remainder[24][25]=str_to_Double("-.2029041421166157799764337515389684192852410619720263791832989e-1");
    rs_remainder[24][26]=str_to_Double("-.19470026604051705645172985140129702066308894329963689098510e-1");
    rs_remainder[24][27]=str_to_Double(".1794021854871665180274068228159639929999849031827885554116462e-1");
    rs_remainder[24][28]=str_to_Double(".3077645518073521651544567869241236315169355624606247614705336683e-2");
    rs_remainder[24][29]=str_to_Double("-.783693630368717858015039625016677228896478875368246083815805e-2");
    rs_remainder[24][30]=str_to_Double(".5014625711058663652926532501389207183106473776026379956823983421e-3");
    rs_remainder[24][31]=str_to_Double(".2307571453163250784560209739773370407459521045512604175920012370e-2");
    rs_remainder[24][32]=str_to_Double("-.4133527633414169794126884231925632588771908795490412011347526119e-3");
    rs_remainder[24][33]=str_to_Double("-.5106238504653671941994486899688369451336417795847232314719947e-3");
    rs_remainder[24][34]=str_to_Double(".1221784431564856187912688836475551708667883218428773524625423162e-3");
    rs_remainder[24][35]=str_to_Double(".899830127045900866418054775361906130916906098565115493244687e-4");
    rs_remainder[24][36]=str_to_Double("-.2365586861346374218991191034876847309074593367769314401004110e-4");
    rs_remainder[24][37]=str_to_Double("-.13064322341458720448844944060361408126740460603281429444604436e-4");
    rs_remainder[24][38]=str_to_Double(".338622402789405277202409966801905460528181007011319377091388605e-5");
    rs_remainder[24][39]=str_to_Double(".1592844265836129246430289631985524141871908628363878512030571e-5");
    rs_remainder[24][40]=str_to_Double("-.3766880084081752098137380702151023484382405772409180623022663268e-6");
    rs_remainder[24][41]=str_to_Double("-.1647184216093678593882856425170058568112468749079305659379551698e-6");
    rs_remainder[24][42]=str_to_Double(".3341998807917489877226196453760834187328279226211862044395946e-7");
    rs_remainder[24][43]=str_to_Double(".1452470960083086081984601564867699547487286710490600118834548500e-7");
    rs_remainder[24][44]=str_to_Double("-.2397625541876055427077696886452726172238516774063676427501834e-8");
    rs_remainder[24][45]=str_to_Double("-.10962331877402629164135818735215413251501278339611705826565375e-8");
    rs_remainder[24][46]=str_to_Double(".139709096422287232625804642367614930624246186107666137575587209e-9");
    rs_remainder[24][47]=str_to_Double(".7108094535487093229989628368857760910443497900444984350469220483e-10");
    rs_remainder[24][48]=str_to_Double("-.65631470375761050558257196571450057651210868248101247187291531e-11");
    rs_remainder[24][49]=str_to_Double("-.397630421103337536507679331490485914314425629003137528356299963e-11");
    rs_remainder[24][50]=str_to_Double(".24091222040930193355000399291592373802786534239226513436545579e-12");
    rs_remainder[24][51]=str_to_Double(".1927838908336445140470567576848346981135079011343226206404329120e-12");
    rs_remainder[24][52]=str_to_Double("-.62054455542777861426743856619665502332015516043010683638898472e-14");
    rs_remainder[24][53]=str_to_Double("-.8139005177875631929024524420697107089546068749227438551497781176e-14");
    rs_remainder[24][54]=str_to_Double(".5438787680180475434057246390711971490688923019083130302684524e-16");
    rs_remainder[24][55]=str_to_Double(".3005815401381318338612735999832681128396978355620231836282105998e-15");
    rs_remainder[24][56]=str_to_Double(".4985652006062300187631866739521546330072050822617668932697970e-17");
    rs_remainder[24][57]=str_to_Double("-.9751086944183827640342458354214724628818047791268904708218602499e-17");
    rs_remainder[24][58]=str_to_Double("-.37411512436879616914539323195665312181677398290596290588632413e-18");
    rs_remainder[24][59]=str_to_Double(".2788564808175114452462950494615411167458706664865618772348696945e-18");
    rs_remainder[24][60]=str_to_Double(".1658683308791252420966274751787780603071641858739863571842749486e-19");
    rs_remainder[24][61]=str_to_Double("-.7048542780497318579831593668799759532954025143379733150327475019e-20");
    rs_remainder[24][62]=str_to_Double("-.570128330758043118534799991117886833533961727849178814741741750e-21");
    rs_remainder[24][63]=str_to_Double(".157697958688926180646279572360546171139395208972778339467918550e-21");
    rs_remainder[24][64]=str_to_Double(".1636838874976325782431034372362463620666434611359702278052973576e-22");
    rs_remainder[24][65]=str_to_Double("-.312168394797053109287338435441645514763483453354979990702657957e-23");
    rs_remainder[24][66]=str_to_Double("-.4057039751041095341834494966795799883225268899331231704469763051e-24");
    rs_remainder[24][67]=str_to_Double(".5449443163844526246512024067815125717355547614370185675465459595e-25");
    rs_remainder[24][68]=str_to_Double(".8842824350743067734157026506438203805039335257158650021013541737e-26");
    rs_remainder[24][69]=str_to_Double("-.8317859045290941826643473718330432967456627693169375294254104429e-27");
    rs_remainder[24][70]=str_to_Double("-.1715075026138022570989084345203687630764169921942605338764049026e-27");
    rs_remainder[24][71]=str_to_Double(".1088759008650021236972729208473922489526069918747553090463625784e-28");
}

void initialize_rs_remainder6(){

    //============================= n = 25 ================================================
    rs_remainder[25][0]=str_to_Double(".76657016671196908302009721734151318078635115956582371421634e-9");
    rs_remainder[25][1]=str_to_Double("-.229727845144909888640847629869883678977771711859752973322931e-7");
    rs_remainder[25][2]=str_to_Double(".49265072349833148451598154449009215829226441076933006030777e-7");
    rs_remainder[25][3]=str_to_Double(".27182310936594360277883259532499670051080119010153509711772e-6");
    rs_remainder[25][4]=str_to_Double("-.17610265986250322593031389767168773535395311952446463935e-5");
    rs_remainder[25][5]=str_to_Double(".518417860470620083859264589577016493780247462270192844115e-5");
    rs_remainder[25][6]=str_to_Double("-.1090206819758161236733657079695476156786852784462515872271e-4");
    rs_remainder[25][7]=str_to_Double(".2181652665796434035014385944690367504915395354534371604107274e-4");
    rs_remainder[25][8]=str_to_Double("-.5940132876086994010261650930734347831806856031697079825e-4");
    rs_remainder[25][9]=str_to_Double(".2167461820440636935337427184154215553121675348650178429880201168e-3");
    rs_remainder[25][10]=str_to_Double("-.7718199923434104712215905292063985333618254349333534361e-3");
    rs_remainder[25][11]=str_to_Double(".2293327940720560662678974121474059838442590513738915051697998653e-2");
    rs_remainder[25][12]=str_to_Double("-.5411246153435761245300470720186962691789349488503179697828733e-2");
    rs_remainder[25][13]=str_to_Double(".9689747213451778064232627646801654724472497752836518689725e-2");
    rs_remainder[25][14]=str_to_Double("-.116811633875915947465853464650289511099549505934195498786e-1");
    rs_remainder[25][15]=str_to_Double(".459011976129258824260043800070530596185206036675865737163e-2");
    rs_remainder[25][16]=str_to_Double(".1604294239450097640741298502021249158005265776777211992535e-1");
    rs_remainder[25][17]=str_to_Double("-.411120294947795839095186772293701447062054283484161812920768145e-1");
    rs_remainder[25][18]=str_to_Double(".451004661732390385356045578018586990441411299623905509457854869e-1");
    rs_remainder[25][19]=str_to_Double("-.89162872929973879770075669001198443908831382609099062792512460e-2");
    rs_remainder[25][20]=str_to_Double("-.4562848136414713339401102766844145424773038436638319237691e-1");
    rs_remainder[25][21]=str_to_Double(".6244943016098243148959155958267372673297645261322230579066856398e-1");
    rs_remainder[25][22]=str_to_Double("-.1719261046216144220131444619071661065300785195114223688788054699e-1");
    rs_remainder[25][23]=str_to_Double("-.39158467457798534974455323151041878725959056010680974128762e-1");
    rs_remainder[25][24]=str_to_Double(".4119753172078465025498729982912213985119899184847772294656788064e-1");
    rs_remainder[25][25]=str_to_Double(".1102568638295974728396062644921072022930223114712438412589199288e-2");
    rs_remainder[25][26]=str_to_Double("-.247193415463639282120779617184755077314009325617169329378025e-1");
    rs_remainder[25][27]=str_to_Double(".97891996781841964783462691683149445620348207903113755401776e-2");
    rs_remainder[25][28]=str_to_Double(".787483994415847034747956109674486396872652883425098615884811190e-2");
    rs_remainder[25][29]=str_to_Double("-.61105865351434157392400530313186404899482965065402524660464998e-2");
    rs_remainder[25][30]=str_to_Double("-.1412162237897989311944712095607332498603488009337682346194644493e-2");
    rs_remainder[25][31]=str_to_Double(".215274253780013371823102025041100550403575211856124133841609540e-2");
    rs_remainder[25][32]=str_to_Double(".8839102324797855126348931834903708497838145270924397026983482751e-4");
    rs_remainder[25][33]=str_to_Double("-.5346417332447029819902753064002476210424837219719431555693079556e-3");
    rs_remainder[25][34]=str_to_Double(".250527971821886344467513473593344553365026878593598971768897e-4");
    rs_remainder[25][35]=str_to_Double(".1018239534214040764996501822810905683868481645214952229243330e-3");
    rs_remainder[25][36]=str_to_Double("-.894072676641162940820212773020448869531261704268809621752476e-5");
    rs_remainder[25][37]=str_to_Double("-.1554253536738975395341522119710761703215208632725144248413819768e-4");
    rs_remainder[25][38]=str_to_Double(".1575307059773149229774666023470186478194712362152968381232962e-5");
    rs_remainder[25][39]=str_to_Double(".1951300551134358113617728674749921746043920296592381650934003119e-5");
    rs_remainder[25][40]=str_to_Double("-.1907426371906409422747277063672336812690933796863109569295872427e-6");
    rs_remainder[25][41]=str_to_Double("-.2047454151572484798786527345658425467305037495709906542634163030e-6");
    rs_remainder[25][42]=str_to_Double(".1715540921730322900761688090617036558467711097148282811917899517e-7");
    rs_remainder[25][43]=str_to_Double(".1814631142310870888599034483562677342888624716180461105963047314e-7");
    rs_remainder[25][44]=str_to_Double("-.11656325601867486767913197295699902049895231535611125646211927e-8");
    rs_remainder[25][45]=str_to_Double("-.13690325183100593163549794013237019150216524032781738529847030e-8");
    rs_remainder[25][46]=str_to_Double(".5786099481044009944556538826299090152627360147973008398378644e-10");
    rs_remainder[25][47]=str_to_Double(".88482161807543831797212751885595265373649390881711493378250499e-10");
    rs_remainder[25][48]=str_to_Double("-.1749521213780614426300571982950722493828420052629330725554008e-11");
    rs_remainder[25][49]=str_to_Double("-.49269393824802922235933865323476066182187996140041640442680618e-11");
    rs_remainder[25][50]=str_to_Double("-.112520200909405024466196727062658537925435339979062357883838e-13");
    rs_remainder[25][51]=str_to_Double(".2375990760691316384034088755760254045682286702058198290892655994e-12");
    rs_remainder[25][52]=str_to_Double(".5604776231741525799446594368108207211082730553995444136762775613e-14");
    rs_remainder[25][53]=str_to_Double("-.9970857184138244231361151127264029555956540231407147897116789305e-14");
    rs_remainder[25][54]=str_to_Double("-.44118094226031959772564584886216801098804975526240097543258171e-15");
    rs_remainder[25][55]=str_to_Double(".3656582915429640454516410974242063312158632857748183329658694985e-15");
    rs_remainder[25][56]=str_to_Double(".2363964474134216671145191355868332323801765707359576640077220414e-16");
    rs_remainder[25][57]=str_to_Double("-.117595333957705145846629799286330229261523550368254680006937623e-16");
    rs_remainder[25][58]=str_to_Double("-.1004663929433059863303188822002952154553853844313322214034435732e-17");
    rs_remainder[25][59]=str_to_Double(".3324740997495497583313815453173510132785319225971960110141334116e-18");
    rs_remainder[25][60]=str_to_Double(".3574971255485572597253219498885525752934979209151022931782450419e-19");
    rs_remainder[25][61]=str_to_Double("-.827256392152008113844236537705958578026052370475735351328489143e-20");
    rs_remainder[25][62]=str_to_Double("-.109448893753371208302142483105911170470268642640233538291340260e-20");
    rs_remainder[25][63]=str_to_Double(".1809481152288189158993480117382050580146945966305552572119069462e-21");
    rs_remainder[25][64]=str_to_Double(".2930942653258994659582202217890824207441532221394440393560536479e-22");
    rs_remainder[25][65]=str_to_Double("-.3462852203461362147590590820024010221409663322576633707100040270e-23");
    rs_remainder[25][66]=str_to_Double("-.694301428097124599246707123683054855074807139869790147676179337e-24");
    rs_remainder[25][67]=str_to_Double(".5730702084814392855907592601344195230616906575637345832741170713e-25");
    rs_remainder[25][68]=str_to_Double(".1466969069525409221569938417965213920147586504104236526343173124e-25");
    rs_remainder[25][69]=str_to_Double("-.79804559994647681759817115549964261089420912405183672888092697e-27");
    rs_remainder[25][70]=str_to_Double("-.2782152666555381914097855497001821998339628642955447173736308697e-27");
    rs_remainder[25][71]=str_to_Double(".8692076048541801267565243176974202150312782222607077289003399068e-29");
    //============================= n = 26 ================================================
    rs_remainder[26][0]=str_to_Double("-.11595880039454603974728216735562518021802314099441985379198e-8");
    rs_remainder[26][1]=str_to_Double(".73787446285845688772617738401136762771474607380378472744084957e-8");
    rs_remainder[26][2]=str_to_Double("-.1460072539016022178050177651921493557979818242287559473646e-7");
    rs_remainder[26][3]=str_to_Double("-.1720088748867783995881997776220799373824312181998872985e-7");
    rs_remainder[26][4]=str_to_Double(".400217468019059187362275396661237702388192388184187198188e-6");
    rs_remainder[26][5]=str_to_Double("-.236292823335592117536681644455930911377885053253772615586e-5");
    rs_remainder[26][6]=str_to_Double(".9526966221326168563300441030574303888895799315756357367e-5");
    rs_remainder[26][7]=str_to_Double("-.3103922013133336995791209299968175689619419479273820612265e-4");
    rs_remainder[26][8]=str_to_Double(".8958532079663728251560659775168562748124187149793313777e-4");
    rs_remainder[26][9]=str_to_Double("-.2428092646906380316320439055304615263692307582435282636e-3");
    rs_remainder[26][10]=str_to_Double(".6247666843558932968282280431632616665513126725101097967027e-3");
    rs_remainder[26][11]=str_to_Double("-.147603909691202630692898549369690418524920056208451820924e-2");
    rs_remainder[26][12]=str_to_Double(".3025111050766487186928028639360146649400295003082927958e-2");
    rs_remainder[26][13]=str_to_Double("-.4970907710604819625750778944577847727995649786680963361252e-2");
    rs_remainder[26][14]=str_to_Double(".5529861037367971078808802685170570489879330682038207197368e-2");
    rs_remainder[26][15]=str_to_Double("-.1028176203578875988485866571477420457151050117919721639527774065e-2");
    rs_remainder[26][16]=str_to_Double("-.11673681464139515403441145076395367264518925210049016223e-1");
    rs_remainder[26][17]=str_to_Double(".2893761750689957200490393164417024340131591678305333358187620258e-1");
    rs_remainder[26][18]=str_to_Double("-.3594369491947415427990047268974591336007451521600407975482894611e-1");
    rs_remainder[26][19]=str_to_Double(".1576368429852552804420065630034078457825041275929040794500542e-1");
    rs_remainder[26][20]=str_to_Double(".2721366848638277852957758391179815154387209173489696659980674456e-1");
    rs_remainder[26][21]=str_to_Double("-.565753823317167870982054107893576668260574601401723917604e-1");
    rs_remainder[26][22]=str_to_Double(".3649356895769360606156970449543901040441431052511178522536232383e-1");
    rs_remainder[26][23]=str_to_Double(".1779225133111790785212610712178218332259347136606453478249301882e-1");
    rs_remainder[26][24]=str_to_Double("-.46851192093438505512471629639871323321945532534558715809941e-1");
    rs_remainder[26][25]=str_to_Double(".2168776219057982546755955501231866451477438271909761144453e-1");
    rs_remainder[26][26]=str_to_Double(".1746178658131926799954243797495188837578003220064288469905277731e-1");
    rs_remainder[26][27]=str_to_Double("-.22154129213803025692716314217676047498431005295717517117962e-1");
    rs_remainder[26][28]=str_to_Double(".8316790913458902482563431368346880026939168213900075524595e-3");
    rs_remainder[26][29]=str_to_Double(".9881296827243122056136201447289855211455185131854553501394688781e-2");
    rs_remainder[26][30]=str_to_Double("-.33418535803791354565947858232913308069119650191604592906343e-2");
    rs_remainder[26][31]=str_to_Double("-.272271143269626062007804039556652892250814092708163084328807e-2");
    rs_remainder[26][32]=str_to_Double(".1602510698766405713844752588719828486184166476286099284602982e-2");
    rs_remainder[26][33]=str_to_Double(".5124073602607563950263543367640423773766326677399515851592410663e-3");
    rs_remainder[26][34]=str_to_Double("-.4643526477295939732209730620766975571846713607133809345882304476e-3");
    rs_remainder[26][35]=str_to_Double("-.687357723105155177613710057511962424940739322249793966121276e-4");
    rs_remainder[26][36]=str_to_Double(".9731129779291620263988102897838954704873569980231438665781269891e-4");
    rs_remainder[26][37]=str_to_Double(".6732608216982505659821498522270897332928325152199710195482029702e-5");
    rs_remainder[26][38]=str_to_Double("-.1582142236726166870907912431191849467682810542910183330573907200e-4");
    rs_remainder[26][39]=str_to_Double("-.5074884073181683941429800831932647378506143674457112798457401565e-6");
    rs_remainder[26][40]=str_to_Double(".2071307623076027933733316583521291358173900674078969129364465370e-5");
    rs_remainder[26][41]=str_to_Double(".3665592532063188484776512343711332955012223806719144737258369e-7");
    rs_remainder[26][42]=str_to_Double("-.22335488362362190640324833064890922881220659208082678033132985e-6");
    rs_remainder[26][43]=str_to_Double("-.3725160539308286897540408453605557535477630358045965522962538020e-8");
    rs_remainder[26][44]=str_to_Double(".2014005805425303805527602850499553816263732317240600145367562615e-7");
    rs_remainder[26][45]=str_to_Double(".475768596323902682984146407104168183574054440658005723952751e-9");
    rs_remainder[26][46]=str_to_Double("-.15352647401139046278985047971076420525860835084321327828236486e-8");
    rs_remainder[26][47]=str_to_Double("-.541202533926139664650053788866743744587590708015850762187129e-10");
    rs_remainder[26][48]=str_to_Double(".9978013847363903889501623632355187298964439635435832477913555686e-10");
    rs_remainder[26][49]=str_to_Double(".4950205252439213468073292247940770163013881641932728097670409233e-11");
    rs_remainder[26][50]=str_to_Double("-.5567416218284825018907402055438972217922511746471411996758913207e-11");
    rs_remainder[26][51]=str_to_Double("-.3654900515079984436989574700418692418981218152494943956280965260e-12");
    rs_remainder[26][52]=str_to_Double(".2682466328042151059102683452675562460958681026476014755910265159e-12");
    rs_remainder[26][53]=str_to_Double(".2224891269212856641772220585418067930933674082673978202819083955e-13");
    rs_remainder[26][54]=str_to_Double("-.112144233876943536786012037970541134057613395314902648231064446e-13");
    rs_remainder[26][55]=str_to_Double("-.11388704804340734502031552739264710453748471016952570878502841e-14");
    rs_remainder[26][56]=str_to_Double(".4083320996425250337518349808027484098986877051078369325586489224e-15");
    rs_remainder[26][57]=str_to_Double(".4980635708676719960509994413115903096372325802395673423073047641e-16");
    rs_remainder[26][58]=str_to_Double("-.1298101346084868467341486570438072770242275576873800688098950108e-16");
    rs_remainder[26][59]=str_to_Double("-.1884727868527172695769293516639782172363915348119945924478308154e-17");
    rs_remainder[26][60]=str_to_Double(".3605273140188690699774359340242773494523327410620927023103500500e-18");
    rs_remainder[26][61]=str_to_Double(".6234633711186325063736754482296501966809800411747834689016795607e-19");
    rs_remainder[26][62]=str_to_Double("-.8728759416638805234440761500364452913512865554592985139450093071e-20");
    rs_remainder[26][63]=str_to_Double("-.1818103949068455508824818635956370733451373045546151870443199526e-20");
    rs_remainder[26][64]=str_to_Double(".1829147921171879672771663871475612246025830085507128426553757187e-21");
    rs_remainder[26][65]=str_to_Double(".4706594074534784642164761724994862932078044121474012929874663371e-22");
    rs_remainder[26][66]=str_to_Double("-.3260825942707134560560457898062935277141417194409623213444074076e-23");
    rs_remainder[26][67]=str_to_Double("-.1087985748846153759989619673009410345801850284512386624111599446e-23");
    rs_remainder[26][68]=str_to_Double(".4738126700444252636970686322585963309179080423114850357849045427e-25");
    rs_remainder[26][69]=str_to_Double(".2256901040441656094692091367325822285209547722285691480856266961e-25");
    rs_remainder[26][70]=str_to_Double("-.4901838518819286141082937377457960079656162363746708917166980546e-27");
    rs_remainder[26][71]=str_to_Double("-.4218583534660784516800292144264575340610842943200692254898421815e-27");
    //============================= n = 27 ================================================
    rs_remainder[27][0]=str_to_Double("-.253596966794264274572006927184481368945187746648402640537e-9");
    rs_remainder[27][1]=str_to_Double("-.11720232202495111792694776493598257847860007473668468770139e-7");
    rs_remainder[27][2]=str_to_Double(".36050322711252613064015965690286483737825400563512696178e-7");
    rs_remainder[27][3]=str_to_Double(".303106590377923966511194622157962813796296474471362730e-7");
    rs_remainder[27][4]=str_to_Double("-.123980443094728351641477193500283819165213053726796621e-6");
    rs_remainder[27][5]=str_to_Double("-.1268217607684121783195364978268940502024042658454574328e-5");
    rs_remainder[27][6]=str_to_Double(".955381956478501672215870744481029594372916704698895894e-5");
    rs_remainder[27][7]=str_to_Double("-.381279541890753870012090447358657694529582596572385650e-4");
    rs_remainder[27][8]=str_to_Double(".11650503187915357301300858336292024270795246593576609014e-3");
    rs_remainder[27][9]=str_to_Double("-.303273584051968652958318435903475024741771287754602262866010e-3");
    rs_remainder[27][10]=str_to_Double(".688143494634774364919434469300552012562899943732537355e-3");
    rs_remainder[27][11]=str_to_Double("-.13249100268412690070275767113011323845777516733199521199e-2");
    rs_remainder[27][12]=str_to_Double(".19900015254935870716170239026797852609854587436472067260264e-2");
    rs_remainder[27][13]=str_to_Double("-.173699720100809446610971915447202337602966433624464050643211136e-2");
    rs_remainder[27][14]=str_to_Double("-.13545182718537965812862860079256676330531316149976992144e-2");
    rs_remainder[27][15]=str_to_Double(".9137941671025989160296648469680547224618142184301127187500166408e-2");
    rs_remainder[27][16]=str_to_Double("-.20166248111948050457672301569495253486904951112709960154e-1");
    rs_remainder[27][17]=str_to_Double(".262853760448589159316668836830698289531436919551432169115e-1");
    rs_remainder[27][18]=str_to_Double("-.1544533919949585400011905457549166697476413483567131687807188041e-1");
    rs_remainder[27][19]=str_to_Double("-.14999569755212808202953288427608708856321577538054914641e-1");
    rs_remainder[27][20]=str_to_Double(".4519889079052693706461090941661308989240183553892168278512e-1");
    rs_remainder[27][21]=str_to_Double("-.4337218842975804849859967791467266702334389294529317805220884e-1");
    rs_remainder[27][22]=str_to_Double(".268486254142280321077083386911253414360482010131918732611415e-2");
    rs_remainder[27][23]=str_to_Double(".389710418556781122950679746144938291267131425481068577132040e-1");
    rs_remainder[27][24]=str_to_Double("-.3751362323443661626378319665654067693596891205467390098047539898e-1");
    rs_remainder[27][25]=str_to_Double("-.1952800740763542821514275391808842376625120325418956484513393409e-3");
    rs_remainder[27][26]=str_to_Double(".2519081252762531718839990505366024210987787603504588114512206903e-1");
    rs_remainder[27][27]=str_to_Double("-.1411586735120332994498813613239957196908268572536164113595841e-1");
    rs_remainder[27][28]=str_to_Double("-.6495566870841559868772337693075805444177528030306839649183e-2");
    rs_remainder[27][29]=str_to_Double(".9140046214359222442219743156414192099918455914838113663249466525e-2");
    rs_remainder[27][30]=str_to_Double("-.333288185100477627545161868244116249217413473854220773576413129e-3");
    rs_remainder[27][31]=str_to_Double("-.3244482250525438137741721648870108541882309313099349217309197558e-2");
    rs_remainder[27][32]=str_to_Double(".8210170673111222715998785805768657239216326268101969058754239126e-3");
    rs_remainder[27][33]=str_to_Double(".7863633866551573897684181881899859345409203698801003846283968013e-3");
    rs_remainder[27][34]=str_to_Double("-.321252944677937918453314212673903176431265204882838218947221e-3");
    rs_remainder[27][35]=str_to_Double("-.1429859601733961468721246798873779090131449610616098201157703082e-3");
    rs_remainder[27][36]=str_to_Double(".779366916696756437844645158700657647296712842444973080611647e-4");
    rs_remainder[27][37]=str_to_Double(".2078700402059628258450672091253460910976610618380183913126157992e-4");
    rs_remainder[27][38]=str_to_Double("-.1384646643080257199941289383822006394497177351669880306052806e-4");
    rs_remainder[27][39]=str_to_Double("-.25479575029365722889817651821947708437803441414710406516334618e-5");
    rs_remainder[27][40]=str_to_Double(".192152883346553670970178938962522579071796076481655367320019e-5");
    rs_remainder[27][41]=str_to_Double(".2749374698658085677620308023061640986728295142322060843011482014e-6");
    rs_remainder[27][42]=str_to_Double("-.21553843184684489209764146361421085380037765321006446713044112e-6");
    rs_remainder[27][43]=str_to_Double("-.2678719978060674737231282365449360676759274623479981250002443450e-7");
    rs_remainder[27][44]=str_to_Double(".1996179418907781044934312699731646734696245591455997145074780e-7");
    rs_remainder[27][45]=str_to_Double(".2365106131912666700697747550806106941089571127616467958611867949e-8");
    rs_remainder[27][46]=str_to_Double("-.15487954905301137766872792627490324829335328224910661607308461e-8");
    rs_remainder[27][47]=str_to_Double("-.1872055380711823506210023721982990560922482523419165281357489015e-9");
    rs_remainder[27][48]=str_to_Double(".1017514869035844689374402041409447510574836585122000755417761177e-9");
    rs_remainder[27][49]=str_to_Double(".1310683466075991203590875741832574728273440992051000956856963084e-10");
    rs_remainder[27][50]=str_to_Double("-.57064955681989027949096810899387760066533992116383548721488335e-11");
    rs_remainder[27][51]=str_to_Double("-.8043750569333337720220678041882759161648193587078479077029230756e-12");
    rs_remainder[27][52]=str_to_Double(".274902564094002207348703156032136689993778582196124765905284240e-12");
    rs_remainder[27][53]=str_to_Double(".43121861680315358615189354760219615710771620825090580835284471e-13");
    rs_remainder[27][54]=str_to_Double("-.1142653085715370388030641086162444285003512143072456650534532609e-13");
    rs_remainder[27][55]=str_to_Double("-.2020770564929415922471164817364070634415745248435859869250402600e-14");
    rs_remainder[27][56]=str_to_Double(".4108550757631167939837885471572564283568855701980947379118379e-15");
    rs_remainder[27][57]=str_to_Double(".8303789359191505263179557924847046249049509006157687354131853174e-16");
    rs_remainder[27][58]=str_to_Double("-.1277907569993752087009973991489662844891902321305630894649599759e-16");
    rs_remainder[27][59]=str_to_Double("-.3004994015556739666594931321515134002371992225894713971529343801e-17");
    rs_remainder[27][60]=str_to_Double(".3424351474780383352493048429892114568229145024122210888093036777e-18");
    rs_remainder[27][61]=str_to_Double(".96218911269226076521378485388198322490535628020565843515036727e-19");
    rs_remainder[27][62]=str_to_Double("-.7813593617980869535237248340776313099489957891044732030657540339e-20");
    rs_remainder[27][63]=str_to_Double("-.2738751755647819987422054692050929361786939596560659488427637319e-20");
    rs_remainder[27][64]=str_to_Double(".1474752701364973419918651084324504951908826767884111325994251294e-21");
    rs_remainder[27][65]=str_to_Double(".6960296497377661122970822401117491726287136370139132820584464054e-22");
    rs_remainder[27][66]=str_to_Double("-.2121858850730519434249408529309628060925576029168922809796166785e-23");
    rs_remainder[27][67]=str_to_Double("-.1585689160293468786797872175646787842018070403766714652636565576e-23");
    rs_remainder[27][68]=str_to_Double(".15901904139350604119448988580483455527801228835581873584695282e-25");
    rs_remainder[27][69]=str_to_Double(".3249717289823212046168109295429437380808346772047038707339565319e-25");
    rs_remainder[27][70]=str_to_Double(".2731074153192675573307638343819314750906583958947626883378477295e-27");
    rs_remainder[27][71]=str_to_Double("-.600885328948682388545106511237623603193819604932529890342160739e-27");
    //============================= n = 28 ================================================
    rs_remainder[28][0]=str_to_Double("-.749196820477556248222630117136560946501196050980930311776987e-9");
    rs_remainder[28][1]=str_to_Double(".6558101631810910602561968216295646137232297006443406866878957e-8");
    rs_remainder[28][2]=str_to_Double("-.190065315060603907961722046792158244647659035961204451797e-7");
    rs_remainder[28][3]=str_to_Double("-.154524465456899782172215284492710067315021690870535784e-7");
    rs_remainder[28][4]=str_to_Double(".30232080473534091063677454713082082133281687210665237325755e-6");
    rs_remainder[28][5]=str_to_Double("-.9525450692453036336680395948499215152245741553321078345414e-6");
    rs_remainder[28][6]=str_to_Double(".1075952086356034344855216809631473433764977485331136228898e-5");
    rs_remainder[28][7]=str_to_Double(".30044013994638833332460494485516610957637421906363168e-5");
    rs_remainder[28][8]=str_to_Double("-.208941176872633309672533280813778183017309127792046433183523e-4");
    rs_remainder[28][9]=str_to_Double(".714965639233761696248733420045602263443900128307126808284610443e-4");
    rs_remainder[28][10]=str_to_Double("-.1779709536657417628374892810779092673602626723623602610e-3");
    rs_remainder[28][11]=str_to_Double(".329369938327507389154248979623326393754590802725404309e-3");
    rs_remainder[28][12]=str_to_Double("-.3671265472283117730894358013511224505794919784550779276161318826e-3");
    rs_remainder[28][13]=str_to_Double("-.21471971832537991382886492364866281888699676000333822e-3");
    rs_remainder[28][14]=str_to_Double(".2387334112528439474011099362389253112843313598066825353941851637e-2");
    rs_remainder[28][15]=str_to_Double("-.711620953639042980054830407480418338394967310164638312003e-2");
    rs_remainder[28][16]=str_to_Double(".137424347121126274495906663360153892246431045451294735476e-1");
    rs_remainder[28][17]=str_to_Double("-.1785208060788390391511470374961743698617716953460271381234041680e-1");
    rs_remainder[28][18]=str_to_Double(".1177432839285199461294059276135593765006294809552677839448624e-1");
    rs_remainder[28][19]=str_to_Double(".8647786060661272749086254087324177603937771628675878580e-2");
    rs_remainder[28][20]=str_to_Double("-.33996788269929641021108769321298907779445997631019567930e-1");
    rs_remainder[28][21]=str_to_Double(".4153265387279030617147553859108535983031021304111658202923905494e-1");
    rs_remainder[28][22]=str_to_Double("-.1615595259926164858209613120826746399669375924824093388842002573e-1");
    rs_remainder[28][23]=str_to_Double("-.2507227666049781439637879441015267365063566345484958264335618845e-1");
    rs_remainder[28][24]=str_to_Double(".4217501609110943861811942741792739204539720299127886647074897392e-1");
    rs_remainder[28][25]=str_to_Double("-.17959966307969139219470783292843842202605903280074866258127e-1");
    rs_remainder[28][26]=str_to_Double("-.1725377805525526429299890031113143280651938608116528702576307352e-1");
    rs_remainder[28][27]=str_to_Double(".2394488388731438199086479039479827354291987882092322735773824241e-1");
    rs_remainder[28][28]=str_to_Double("-.3935068765836802221968727873786875172999103754339324599158883975e-2");
    rs_remainder[28][29]=str_to_Double("-.106533123139555136420819058305601678039512037266465209805595e-1");
    rs_remainder[28][30]=str_to_Double(".63511841120714686542519046655165379695519147994137190621850e-2");
    rs_remainder[28][31]=str_to_Double(".21919859439737149295229437262490298636727474711075384208687633e-2");
    rs_remainder[28][32]=str_to_Double("-.299221213664768190193859699223492957456020783398363376729453e-2");
    rs_remainder[28][33]=str_to_Double("-.7637842135482931586365622484329143288816547393195056793861387590e-6");
    rs_remainder[28][34]=str_to_Double(".8760419047883342996038421342593340895249218741943972572420326636e-3");
    rs_remainder[28][35]=str_to_Double("-.138957112806171371182539167571624087051766805530125275280273680e-3");
    rs_remainder[28][36]=str_to_Double("-.1862542477377321986562316882269207639760627411433890598250628501e-3");
    rs_remainder[28][37]=str_to_Double(".478120779551785487371140015419409325469711007027316506507017019e-4");
    rs_remainder[28][38]=str_to_Double(".30994223757343316864432442071756161146628921451545505097355e-4");
    rs_remainder[28][39]=str_to_Double("-.9951219886564098325270102017040784044675686285753498512234112331e-5");
    rs_remainder[28][40]=str_to_Double("-.4232755503855732827316117321081124181939904886149306532773827067e-5");
    rs_remainder[28][41]=str_to_Double(".1513993597026867581499950085560337524120584366624957516089017356e-5");
    rs_remainder[28][42]=str_to_Double(".4897710675501390232399839207795906912883069328303647284731830015e-6");
    rs_remainder[28][43]=str_to_Double("-.1801073296496054014407297471062196360969735526564768193664736747e-6");
    rs_remainder[28][44]=str_to_Double("-.4898032916277485460944521306245479263245655458825989046700609e-7");
    rs_remainder[28][45]=str_to_Double(".1733857705343617954378078348462096207694416625695623429841961027e-7");
    rs_remainder[28][46]=str_to_Double(".4276828992839483599535241399504778419740906456726967657757339729e-8");
    rs_remainder[28][47]=str_to_Double("-.1379154107305475068100006577133473141583313167810454604249226451e-8");
    rs_remainder[28][48]=str_to_Double("-.3272531378452496989200781908195051353554265099000848065393231181e-9");
    rs_remainder[28][49]=str_to_Double(".9190430018053806247715435602742739416243289934659134496936386099e-10");
    rs_remainder[28][50]=str_to_Double(".2195732706467141430712762503571199911455387546068343559891922390e-10");
    rs_remainder[28][51]=str_to_Double("-.5179623351653217955293629076434500952764347394218875592486887413e-11");
    rs_remainder[28][52]=str_to_Double("-.1292031586469230680126929212547832910916761349556010443255316746e-11");
    rs_remainder[28][53]=str_to_Double(".2484110615229854003620502751138527145686683387957702129206784430e-12");
    rs_remainder[28][54]=str_to_Double(".6673724496106995885287501667654655149426943857185923729020928408e-13");
    rs_remainder[28][55]=str_to_Double("-.1016740878911516823381650078645182509919427043151621869850637816e-13");
    rs_remainder[28][56]=str_to_Double("-.30320239575224085265987685672386092438912509577140341614747447e-14");
    rs_remainder[28][57]=str_to_Double(".3546978918371700276530357129243130130244193682007772709626266249e-15");
    rs_remainder[28][58]=str_to_Double(".1215107407777588859077819771073802647395665835623975618737792941e-15");
    rs_remainder[28][59]=str_to_Double("-.1046089419773684802098925748028314399262542746319882779731701790e-16");
    rs_remainder[28][60]=str_to_Double("-.4310221209945927332799983725157768987554212627554076808614154049e-17");
    rs_remainder[28][61]=str_to_Double(".2550054821314864689676380840205124079992452579336075286327588548e-18");
    rs_remainder[28][62]=str_to_Double(".1358231599711011440399315992746098298520340633116073437677371538e-18");
    rs_remainder[28][63]=str_to_Double("-.4827241058583366835801948190153766658590356454969338292245272226e-20");
    rs_remainder[28][64]=str_to_Double("-.3816094007030157088833979957821002302698925462615937848788021723e-20");
    rs_remainder[28][65]=str_to_Double(".5549040166564424969327919428510067508136632278508035560215251e-22");
    rs_remainder[28][66]=str_to_Double(".9592523157362015383157080900292234252249493629163745747424857109e-22");
    rs_remainder[28][67]=str_to_Double(".430473367682833104978714242083339777957561320662283712829323e-24");
    rs_remainder[28][68]=str_to_Double("-.2164052427693684629479871857304872811714897572454695417799082928e-23");
    rs_remainder[28][69]=str_to_Double("-.4796193166638951977651804006149079511167237107671242486777627213e-25");
    rs_remainder[28][70]=str_to_Double(".4393162482547674261303628247031620277427774873084433311156360982e-25");
    rs_remainder[28][71]=str_to_Double(".17172912063421951727017883397084792509849148130564646542748312e-26");
    //============================= n = 29 ================================================
    rs_remainder[29][0]=str_to_Double("-.823824226421150634141409185505449366814055024843790078873e-9");
    rs_remainder[29][1]=str_to_Double("-.422059008570911648636501263991231208246553885345897518313e-8");
    rs_remainder[29][2]=str_to_Double(".275248395458094343478797152222232379853847004231870701528773e-7");
    rs_remainder[29][3]=str_to_Double("-.83583035891856681075212261181314009220631722484100394304e-7");
    rs_remainder[29][4]=str_to_Double(".31564608198700959803723209162211331843129801070801417297e-6");
    rs_remainder[29][5]=str_to_Double("-.130606191557994729575670766239286231988712017791675375e-5");
    rs_remainder[29][6]=str_to_Double(".4167899035050484201326022258165592388119233824110418039773e-5");
    rs_remainder[29][7]=str_to_Double("-.98467337635562622194163445828600467236162766014759244e-5");
    rs_remainder[29][8]=str_to_Double(".185835910335491579641861246368879323056943302676613014e-4");
    rs_remainder[29][9]=str_to_Double("-.34466099475394960851336578344844184195543186816761175e-4");
    rs_remainder[29][10]=str_to_Double(".85988878463089488994814256305918458121339374823665475e-4");
    rs_remainder[29][11]=str_to_Double("-.2781369645424096821651752342142789821163135337097858884e-3");
    rs_remainder[29][12]=str_to_Double(".876894217249270755337429968518472201185492338982051676e-3");
    rs_remainder[29][13]=str_to_Double("-.23510360930467231047646032311851562968780882873563826910e-2");
    rs_remainder[29][14]=str_to_Double(".514122294584946837996316855318235202491887424433584371e-2");
    rs_remainder[29][15]=str_to_Double("-.8858590464282491236671825208080197927940474421099944043234274550e-2");
    rs_remainder[29][16]=str_to_Double(".110878883927664927646943236277561103231007036875626632759e-1");
    rs_remainder[29][17]=str_to_Double("-.7246016673347254097184044129408540888999602082023775509569e-2");
    rs_remainder[29][18]=str_to_Double("-.6186602253205208940872612429630567409728146376180263657802604207e-2");
    rs_remainder[29][19]=str_to_Double(".25243151560057238469011486161495296660328193210332014782460101e-1");
    rs_remainder[29][20]=str_to_Double("-.35516564500935319601674625201055276172134020991930905579e-1");
    rs_remainder[29][21]=str_to_Double(".2193595742775025163122923656924782555517199000758324019910530557e-1");
    rs_remainder[29][22]=str_to_Double(".11920036880930129193402792088351230348540645200815912426691e-1");
    rs_remainder[29][23]=str_to_Double("-.3812436961555753161353945050049087002867809372471444584448e-1");
    rs_remainder[29][24]=str_to_Double(".302261962515979044279295068010436286612964162023347188425318432e-1");
    rs_remainder[29][25]=str_to_Double(".3272309095148262058960615541182795246689727066175637405181581666e-2");
    rs_remainder[29][26]=str_to_Double("-.253492817087109763145511033509789873219726678508250768780e-1");
    rs_remainder[29][27]=str_to_Double(".160189748001625779791958475014851361539513604356938020836e-1");
    rs_remainder[29][28]=str_to_Double(".517585726707544356293272319373611445703058216291296330021e-2");
    rs_remainder[29][29]=str_to_Double("-.1121965944453955468272391574236291544614269386882446279430e-1");
    rs_remainder[29][30]=str_to_Double(".258044780131536344347126983254012235084824956847532579393949e-2");
    rs_remainder[29][31]=str_to_Double(".3746493884896237553167310381453060495315593171779067979910e-2");
    rs_remainder[29][32]=str_to_Double("-.21477889910728752694968375365553211747817100865854758771446e-2");
    rs_remainder[29][33]=str_to_Double("-.6898184884379290475105620205778779580952104629311067679941e-3");
    rs_remainder[29][34]=str_to_Double(".7881790391678080104223463127875302105142618779175530018577568616e-3");
    rs_remainder[29][35]=str_to_Double(".447220744142181866121683749287689127649840109796985347372e-4");
    rs_remainder[29][36]=str_to_Double("-.19374870404826837775293213922229380808740535260155773223283989e-3");
    rs_remainder[29][37]=str_to_Double(".1274306247461084050874653990249333860382624413584083207502e-4");
    rs_remainder[29][38]=str_to_Double(".3587917580698515514519021931756543737025941078869599787109738514e-4");
    rs_remainder[29][39]=str_to_Double("-.47911306904014397130125948051671915159749529869275094210377e-5");
    rs_remainder[29][40]=str_to_Double("-.53056541674189691329843262764479570802778219271317774754631e-5");
    rs_remainder[29][41]=str_to_Double(".90316070755520814081364179090609510074430308618768549561547e-6");
    rs_remainder[29][42]=str_to_Double(".6488642645445940582707311651650864527603904909373685363378094441e-6");
    rs_remainder[29][43]=str_to_Double("-.1201581298350592742004139515818962748739785851495725274663657e-6");
    rs_remainder[29][44]=str_to_Double("-.6710301118938010351192672009217975322598065837329741793588885266e-7");
    rs_remainder[29][45]=str_to_Double(".1234444167710468179329197454924481712212341388665293963089686259e-7");
    rs_remainder[29][46]=str_to_Double(".5950783631869236333139433616491094852762200121835719633901545281e-8");
    rs_remainder[29][47]=str_to_Double("-.101897463154561198262578089717417692547933543203127729497694397e-8");
    rs_remainder[29][48]=str_to_Double("-.4564802349910672052591130617655485720503583596325277605123372595e-9");
    rs_remainder[29][49]=str_to_Double(".690314205702185335569867139893869792776073382659456527252335e-10");
    rs_remainder[29][50]=str_to_Double(".3045982335120740035758886256108799863595556160577414036712638149e-10");
    rs_remainder[29][51]=str_to_Double("-.388238714985342749447954939244258779585429735412019575180585e-11");
    rs_remainder[29][52]=str_to_Double("-.1775365080876277810020645995181245692641907730821118726912298795e-11");
    rs_remainder[29][53]=str_to_Double(".1819884910971724101939400856254928539692407534910489473156855296e-12");
    rs_remainder[29][54]=str_to_Double(".9071235829721497627788819040323406999672062021915827249777414729e-13");
    rs_remainder[29][55]=str_to_Double("-.7075353429047733487667860565777888308830175037777804929664767915e-14");
    rs_remainder[29][56]=str_to_Double("-.4077448811579711203169594466663890264073819756596240710321954966e-14");
    rs_remainder[29][57]=str_to_Double(".2234062236327633983687034605084105987493463810906828229843969116e-15");
    rs_remainder[29][58]=str_to_Double(".16180707311537446526543567025583556880953415981732417918417312e-15");
    rs_remainder[29][59]=str_to_Double("-.5373746792665102947883761018441846774636292310694477248400495074e-17");
    rs_remainder[29][60]=str_to_Double("-.56890870966174326973019212201664295681914444705751873748840029e-17");
    rs_remainder[29][61]=str_to_Double(".7533799689474658620746028461731594015938224128167307211525861e-19");
    rs_remainder[29][62]=str_to_Double(".17784265472751863656917060774268573967818127759116425433940179e-18");
    rs_remainder[29][63]=str_to_Double(".9479987729044547035170848474637646329529655485075207574706592609e-21");
    rs_remainder[29][64]=str_to_Double("-.4959085813708501013907773949227035139686585225921621181935691070e-20");
    rs_remainder[29][65]=str_to_Double("-.1132805197483122522877323171508634951061636291999628428195601824e-21");
    rs_remainder[29][66]=str_to_Double(".1237150077562487080665890323961876541792656462042425814377014734e-21");
    rs_remainder[29][67]=str_to_Double(".4912715863320108799381307000141285139682343739647212446556550e-23");
    rs_remainder[29][68]=str_to_Double("-.2768101568036405385387288133801770799325970684402873704263263797e-23");
    rs_remainder[29][69]=str_to_Double("-.1561900067933986715741998844139510073237315310513187516041867046e-24");
    rs_remainder[29][70]=str_to_Double(".5565380441442192820356373663917214855725533868387073079848604445e-25");
    rs_remainder[29][71]=str_to_Double(".4095692056580320525535694106637019336903841364098362647514008058e-26");
}

void initialize_rs_remainder7(){

    //============================= n = 30 ================================================
    rs_remainder[30][0]=str_to_Double("-.489261152572672696663979474601699414750819035038062393718e-9");
    rs_remainder[30][1]=str_to_Double(".524983328121722736096612047226550553537111485239071471538e-8");
    rs_remainder[30][2]=str_to_Double("-.20243433012341022729597878248191234072545873974876425576e-7");
    rs_remainder[30][3]=str_to_Double(".1766675308052405408702758573303377410935429008295571549e-7");
    rs_remainder[30][4]=str_to_Double(".87989370687685367158882389482661518251417502738984931e-7");
    rs_remainder[30][5]=str_to_Double("-.936774868540269454616044293155964960400313531016196939391e-7");
    rs_remainder[30][6]=str_to_Double("-.12665241236041028429228266745995911951006427540843011345e-5");
    rs_remainder[30][7]=str_to_Double(".70521416053031705681775748909639028673174819310082347701e-5");
    rs_remainder[30][8]=str_to_Double("-.22974413512869367063918594308739005775026112893591289161e-4");
    rs_remainder[30][9]=str_to_Double(".60602760502509647271367864656133473397550728274233543e-4");
    rs_remainder[30][10]=str_to_Double("-.1473798932403361736146902954902377317383763570089404798e-3");
    rs_remainder[30][11]=str_to_Double(".348573722434372486182959139795034896420521729832292446e-3");
    rs_remainder[30][12]=str_to_Double("-.798896562296801493723607269823254046438623564333580635540664745e-3");
    rs_remainder[30][13]=str_to_Double(".1708779998843349339190917290618869353773136668976900570e-2");
    rs_remainder[30][14]=str_to_Double("-.3241695104223534858627192721853082498330086698349545309e-2");
    rs_remainder[30][15]=str_to_Double(".5117130649383651997593179792985095187603108900968546983e-2");
    rs_remainder[30][16]=str_to_Double("-.5981335107687833594640074799270779154736065964039092757e-2");
    rs_remainder[30][17]=str_to_Double(".31984775175406990972958739721181254024604430540727178e-2");
    rs_remainder[30][18]=str_to_Double(".564639212409623608547738472187882698996875687897706694e-2");
    rs_remainder[30][19]=str_to_Double("-.190311191234340709123470580638538131539722440600162561347524e-1");
    rs_remainder[30][20]=str_to_Double(".284433762728436630008468534003712185960376551400631771461e-1");
    rs_remainder[30][21]=str_to_Double("-.2206915150624458578962292939429218259208467800214031389e-1");
    rs_remainder[30][22]=str_to_Double("-.28537797689683357254239810427934540074943663603874461855508e-2");
    rs_remainder[30][23]=str_to_Double(".30185145659747836795916391054585668233916992953869270015598e-1");
    rs_remainder[30][24]=str_to_Double("-.348910666911264243326170156529724511261959101112622633221e-1");
    rs_remainder[30][25]=str_to_Double(".10397304454301645792927898324316232001918558811833463111453899e-1");
    rs_remainder[30][26]=str_to_Double(".190080578043180200381176732802739404738880510319995159076562312e-1");
    rs_remainder[30][27]=str_to_Double("-.2395681268106882460606096728648789380809034715150720273707e-1");
    rs_remainder[30][28]=str_to_Double(".512031185304778407871531188106628488873858780442486366897e-2");
    rs_remainder[30][29]=str_to_Double(".1096602510075151364600869098913773366671752319597389201359774e-1");
    rs_remainder[30][30]=str_to_Double("-.8777074658807157631096899342738966110017307815044964110551188e-2");
    rs_remainder[30][31]=str_to_Double("-.1069553202415177799877886907647860068282732891966736416427e-2");
    rs_remainder[30][32]=str_to_Double(".4152171700880297064845174211275739052313749282488325787268327858e-2");
    rs_remainder[30][33]=str_to_Double("-.9936128381276857738196788012406316447528980108266698843790344570e-3");
    rs_remainder[30][34]=str_to_Double("-.1125848982888321767194754015837534681131353444831229297860197904e-2");
    rs_remainder[30][35]=str_to_Double(".56415211210638438822264768366213508855989531237811075308990181e-3");
    rs_remainder[30][36]=str_to_Double(".1957819763582515207300580848014499070293912032140893130032493017e-3");
    rs_remainder[30][37]=str_to_Double("-.1677119309762634910238474935330098428648085982826045720585667023e-3");
    rs_remainder[30][38]=str_to_Double("-.20976619360579804329875361615275641203505902449909477228067e-4");
    rs_remainder[30][39]=str_to_Double(".3495920590052310136815556406053395719446882883703606158747745755e-4");
    rs_remainder[30][40]=str_to_Double(".7825363253395685552713340822144253422869322816417445881733696128e-6");
    rs_remainder[30][41]=str_to_Double("-.5612641008382723483403461133987704043031112390172983500900508684e-5");
    rs_remainder[30][42]=str_to_Double(".177103520018921694531397164799851264379931224767079761504117e-6");
    rs_remainder[30][43]=str_to_Double(".727565530875682610628307908747291178495766197723652907957206e-6");
    rs_remainder[30][44]=str_to_Double("-.4278627022925643074449990282832595212466415379255202407569892672e-7");
    rs_remainder[30][45]=str_to_Double("-.7832596323339395555664531892060806912028323416659807717928666006e-7");
    rs_remainder[30][46]=str_to_Double(".5421195075349787855402079259984386916765019908147410548648663685e-8");
    rs_remainder[30][47]=str_to_Double(".7130908340172162691577566555040177448003547431047065474349530e-8");
    rs_remainder[30][48]=str_to_Double("-.488322856907198194754560590237083143304396008864939115177616e-9");
    rs_remainder[30][49]=str_to_Double("-.5557803654168309810775739963218296490737426466822266674238124726e-9");
    rs_remainder[30][50]=str_to_Double(".3363106477923630974560598220031066234373524487964747691695307950e-10");
    rs_remainder[30][51]=str_to_Double(".3740692499881033740577572659596895674814186611650090705475843561e-10");
    rs_remainder[30][52]=str_to_Double("-.180032629747367253607637361916998539209183248076459797354420e-11");
    rs_remainder[30][53]=str_to_Double("-.2188534075758802446366138484550694603880935391693368236616571301e-11");
    rs_remainder[30][54]=str_to_Double(".7292707799789006345927919141715024826573098880515124019419702149e-13");
    rs_remainder[30][55]=str_to_Double(".11190411048795646463193078135572549143064189684530198272878836e-12");
    rs_remainder[30][56]=str_to_Double("-.194872048366255728182092428904760487965063841761397437459124e-14");
    rs_remainder[30][57]=str_to_Double("-.5024156324786383755870079252995686625961172837595818130174445e-14");
    rs_remainder[30][58]=str_to_Double(".6006045394999547650230894217832968453215527502356564300220157507e-17");
    rs_remainder[30][59]=str_to_Double(".1989013525872186979651000597839019378053287635912826612060183431e-15");
    rs_remainder[30][60]=str_to_Double(".297093102297563540631158467945153643211157991595617946626681e-17");
    rs_remainder[30][61]=str_to_Double("-.6970132616851495380644383810266780103172086527841227905527696780e-17");
    rs_remainder[30][62]=str_to_Double("-.2151636398141589713453818365759227633409925931897416066934313648e-18");
    rs_remainder[30][63]=str_to_Double(".2169602627857947345574816384848288265848291195771456036292232911e-18");
    rs_remainder[30][64]=str_to_Double(".1013138056703166509932895032714650646955961778580484757139912e-19");
    rs_remainder[30][65]=str_to_Double("-.6016872885501412756146901068156696273450445382336379126257142e-20");
    rs_remainder[30][66]=str_to_Double("-.3771319516860677648297150660612802480912950600475528463075503911e-21");
    rs_remainder[30][67]=str_to_Double(".1490351102668121848241791520658557933686457863147929168144268358e-21");
    rs_remainder[30][68]=str_to_Double(".1180883648833381306527419101596632559726145905138458068507682427e-22");
    rs_remainder[30][69]=str_to_Double("-.3302888733317140593321149549512359653080674264462393220813999460e-23");
    rs_remainder[30][70]=str_to_Double("-.3203218243324538939101624613848188311822810582718968845741204138e-24");
    rs_remainder[30][71]=str_to_Double(".655409083867126280423002769080788874373884206651899960081751233e-25");
    //============================= n = 31 ================================================
    rs_remainder[31][0]=str_to_Double("-.112366878866443341732533548916560583445006773451806685467e-8");
    rs_remainder[31][1]=str_to_Double(".9215299033854126920475162200020787456575136881882895169230e-9");
    rs_remainder[31][2]=str_to_Double(".1390149258566388923739877071180050809465169409714470279e-7");
    rs_remainder[31][3]=str_to_Double("-.903924189665126429002825003981974837397560861313698768e-7");
    rs_remainder[31][4]=str_to_Double(".280402632309097497089917547347554502916878840861066409e-6");
    rs_remainder[31][5]=str_to_Double("-.363650890209864959556303815427177186918452603873980055390981028e-6");
    rs_remainder[31][6]=str_to_Double("-.841130278825544792784743599865487786274471063097540699430229e-6");
    rs_remainder[31][7]=str_to_Double(".68428787332221378290344844775012030447502792098320788e-5");
    rs_remainder[31][8]=str_to_Double("-.261928121703109294459787132898081174613046732586734291e-4");
    rs_remainder[31][9]=str_to_Double(".7754043152269101736738834798196450424381372806231609311e-4");
    rs_remainder[31][10]=str_to_Double("-.197775695656667024471720632893366658303478608983060171e-3");
    rs_remainder[31][11]=str_to_Double(".45072135825020850658767216521487878113641661891963715e-3");
    rs_remainder[31][12]=str_to_Double("-.9199870496871416909505585191358471954800927703276562e-3");
    rs_remainder[31][13]=str_to_Double(".16399103354725288583048011409858349609342667744232485e-2");
    rs_remainder[31][14]=str_to_Double("-.2398068352273634060025093777706532127297690839788815501517e-2");
    rs_remainder[31][15]=str_to_Double(".2420107574218220496727039391885057793906299661056432989e-2");
    rs_remainder[31][16]=str_to_Double("-.2311919045535261704987900655221974383324472014619704964173882813e-3");
    rs_remainder[31][17]=str_to_Double("-.5635555639073707434189352866780371728814200786244201544763e-2");
    rs_remainder[31][18]=str_to_Double(".14601180014694841130570023250878231267557394720798975003692e-1");
    rs_remainder[31][19]=str_to_Double("-.21792588001704207876925120729995510623021102810804144406665e-1");
    rs_remainder[31][20]=str_to_Double(".190171372171404606345017795704618763961964187133367677971267e-1");
    rs_remainder[31][21]=str_to_Double("-.1719465116312909876474864879660914309810393291158525824001309029e-2");
    rs_remainder[31][22]=str_to_Double("-.2224020835437325267354723073936762824881397689972070070413e-1");
    rs_remainder[31][23]=str_to_Double(".33789704174852092259139504015537929995058925025265272606e-1");
    rs_remainder[31][24]=str_to_Double("-.1983223942036546768739842494781204352470704888108034367405e-1");
    rs_remainder[31][25]=str_to_Double("-.89937576473994134279906477782338766061896893425468565996725e-2");
    rs_remainder[31][26]=str_to_Double(".2538274942780326674024698605828487046535586856086362819551519105e-1");
    rs_remainder[31][27]=str_to_Double("-.15519609437294319003933700624315044779921437018464874064e-1");
    rs_remainder[31][28]=str_to_Double("-.4990339405443508653196972540992278697331154293440942518933322e-2");
    rs_remainder[31][29]=str_to_Double(".1255260136655518938531941408103319820926242730934838808367e-1");
    rs_remainder[31][30]=str_to_Double("-.45624644029603651564267488992644928348886548544551400612e-2");
    rs_remainder[31][31]=str_to_Double("-.3736067436073313439745923651353676380145480730319217976269044e-2");
    rs_remainder[31][32]=str_to_Double(".352798221602809797768034183287621892661579378598671314279e-2");
    rs_remainder[31][33]=str_to_Double(".1665765019223610138890181620888318357593685590828360880377e-3");
    rs_remainder[31][34]=str_to_Double("-.1259509343174128381385647890805498617006555934060270405375e-2");
    rs_remainder[31][35]=str_to_Double(".2670492100193697849048492426424778695138715195983761069048487682e-3");
    rs_remainder[31][36]=str_to_Double(".2903944733097282166975326397251318280079398450937289618973970604e-3");
    rs_remainder[31][37]=str_to_Double("-.11634995641615849119664247465191364120922322262864926259766e-3");
    rs_remainder[31][38]=str_to_Double("-.4778088697986586591622019942303492756835817737732261348889e-4");
    rs_remainder[31][39]=str_to_Double(".2877946553954117820692541685214409204783349122059561421906749305e-4");
    rs_remainder[31][40]=str_to_Double(".5881322058139710115999439305039958225912826551655017940371518981e-5");
    rs_remainder[31][41]=str_to_Double("-.5127299903073003905866127524439351337530100028749381329354142799e-5");
    rs_remainder[31][42]=str_to_Double("-.55745083938577189175600651058753187288621649393028425808680e-6");
    rs_remainder[31][43]=str_to_Double(".71320248191533360039715481148981004581911696143292639566674e-6");
    rs_remainder[31][44]=str_to_Double(".4192144878827799492999775908365219236373952823798609381728758111e-7");
    rs_remainder[31][45]=str_to_Double("-.80679714511172268002154607172373060519134607204212592708722e-7");
    rs_remainder[31][46]=str_to_Double("-.26632009314115194581921551515251241036153321978111194927382085e-8");
    rs_remainder[31][47]=str_to_Double(".760682996178938266272972464771499521085709716666718345370546295e-8");
    rs_remainder[31][48]=str_to_Double(".16519113418990471071444033922063017874662388001029563385781077e-9");
    rs_remainder[31][49]=str_to_Double("-.6075679945860251405710685131522672982059024650336307784427713697e-9");
    rs_remainder[31][50]=str_to_Double("-.1189894964289033245444598821874468200985796983412409705444803276e-10");
    rs_remainder[31][51]=str_to_Double(".4158588222976657464819597436279240638100456563630036175427429e-10");
    rs_remainder[31][52]=str_to_Double(".9708608024189314875488695393190092341204691647493588213080433826e-12");
    rs_remainder[31][53]=str_to_Double("-.2460481481364122833136014213434145473811152123068991800078919123e-11");
    rs_remainder[31][54]=str_to_Double("-.7600548427245661663489763843138312720292460733942018819211263114e-13");
    rs_remainder[31][55]=str_to_Double(".1267057624528961475827865623854002182346276665381401417113530028e-12");
    rs_remainder[31][56]=str_to_Double(".5177590007927757588657295895695072026793018237571793479064330179e-14");
    rs_remainder[31][57]=str_to_Double("-.5711265153573086744287084738006339727241008807759799768974634528e-14");
    rs_remainder[31][58]=str_to_Double("-.2994977087135543250818270222368163280874764172866567250128194437e-15");
    rs_remainder[31][59]=str_to_Double(".2264158770851663181264899227140357941076815073277043437632331147e-15");
    rs_remainder[31][60]=str_to_Double(".1475950735986085214829365233768562994177960577830406948244246734e-16");
    rs_remainder[31][61]=str_to_Double("-.7926582789715641839107617070655604781004659382206961416614380216e-17");
    rs_remainder[31][62]=str_to_Double("-.6260060402092189439507462863711894037285344871974084072009948337e-18");
    rs_remainder[31][63]=str_to_Double(".245887110797589988572013793744688577426125097491040049106610509e-18");
    rs_remainder[31][64]=str_to_Double(".2309599310659422275989222572753196457500733197373779726689612615e-19");
    rs_remainder[31][65]=str_to_Double("-.6776295885473871070128096029758488181881130451816365600248016092e-20");
    rs_remainder[31][66]=str_to_Double("-.7483113929728690290627076681135950099177488449547789570131223923e-21");
    rs_remainder[31][67]=str_to_Double(".1661823441389295815144946949876149805488299435233451208087279567e-21");
    rs_remainder[31][68]=str_to_Double(".2146699836211093095705865834347329777425658686989021025610503687e-22");
    rs_remainder[31][69]=str_to_Double("-.3628172907916144497181819774857827365199566351943788862723066109e-23");
    rs_remainder[31][70]=str_to_Double("-.5490835204960151132847116636714340562010515924347230092331259308e-24");
    rs_remainder[31][71]=str_to_Double(".7041198581092128549467258096920569726321943298541216155223611968e-25");
    //============================= n = 32 ================================================
    rs_remainder[32][0]=str_to_Double("-.3152145323174493799425038664812579923533988926192163278039e-9");
    rs_remainder[32][1]=str_to_Double(".383202503386141049734094620023916118964357608889993821280e-8");
    rs_remainder[32][2]=str_to_Double("-.1856495021880269373317560578909116032708179575120600271e-7");
    rs_remainder[32][3]=str_to_Double(".427338877176076513913443101925966733794222193389436587e-7");
    rs_remainder[32][4]=str_to_Double("-.486599893649170182697524368308811407407592026119160485e-7");
    rs_remainder[32][5]=str_to_Double(".10607293753187827214927861006009788739155696532702668e-6");
    rs_remainder[32][6]=str_to_Double("-.4715893927951378681230153754134245056610494132949604687e-6");
    rs_remainder[32][7]=str_to_Double(".9762280263666987679066853260678412353228049202458465392023716782e-6");
    rs_remainder[32][8]=str_to_Double(".871680625632131787784459910824105210923694477357683301614318e-6");
    rs_remainder[32][9]=str_to_Double("-.128495919632132832034379339199656522551312980769223948009791029e-4");
    rs_remainder[32][10]=str_to_Double(".51003933736738897925945788271692677907224998813742344457023628e-4");
    rs_remainder[32][11]=str_to_Double("-.140308710054727308881114818747635151527703440936210712331e-3");
    rs_remainder[32][12]=str_to_Double(".3052278724817377467767395272203292473509056137141212e-3");
    rs_remainder[32][13]=str_to_Double("-.53009410887193642844788640555349441223884339496673782314e-3");
    rs_remainder[32][14]=str_to_Double(".66083424551088896733787311329322525846777904302769571e-3");
    rs_remainder[32][15]=str_to_Double("-.2537978392905218926078360986798724837893767822523202e-3");
    rs_remainder[32][16]=str_to_Double("-.1492014900891883005999068315830735150099930585127391674923805427e-2");
    rs_remainder[32][17]=str_to_Double(".53844129763624062376053458411818494721835023993808651e-2");
    rs_remainder[32][18]=str_to_Double("-.1114876079923779184297677975325458606270301882942500870556673948e-1");
    rs_remainder[32][19]=str_to_Double(".1599694036625622050128071691529168936040034404328253571581e-1");
    rs_remainder[32][20]=str_to_Double("-.1463217299343638241182498624021842955576223670599724091e-1");
    rs_remainder[32][21]=str_to_Double(".2854691163786062083408209180866474108030867330450256992246e-2");
    rs_remainder[32][22]=str_to_Double(".1619037102462348856803705411635772787816406376583502753e-1");
    rs_remainder[32][23]=str_to_Double("-.2974225030713989647917031079663114402980684751552699159e-1");
    rs_remainder[32][24]=str_to_Double(".2421016142144789227486492719203459205636920304018568578058674195e-1");
    rs_remainder[32][25]=str_to_Double("-.58836118569388250234567128016291800747829969109965409e-3");
    rs_remainder[32][26]=str_to_Double("-.2154874186843754100938159231158154625306340972027661749452043195e-1");
    rs_remainder[32][27]=str_to_Double(".2231794676470309526508642523795323828371455747637283710125810591e-1");
    rs_remainder[32][28]=str_to_Double("-.4069153953548120080854959590634987708201115826548450555350723566e-2");
    rs_remainder[32][29]=str_to_Double("-.1157611221261540716908789706964449809914799886715982224584152e-1");
    rs_remainder[32][30]=str_to_Double(".1037309195018704956000798668116119687433384621523153100065e-1");
    rs_remainder[32][31]=str_to_Double("-.186380505988200800419126705534510593233429021803598204e-4");
    rs_remainder[32][32]=str_to_Double("-.4961666789431291475012282745713267326340115051198408751329e-2");
    rs_remainder[32][33]=str_to_Double(".2209423884194725618522718867408538059368038338097179428690854650e-2");
    rs_remainder[32][34]=str_to_Double(".10789950149750784122437681732497863859813100597831290125505e-2");
    rs_remainder[32][35]=str_to_Double("-.1112417962257655475546015914172552567583447364588570930453126585e-2");
    rs_remainder[32][36]=str_to_Double("-.3419337218127137785146560925061746965041233174840184821718e-4");
    rs_remainder[32][37]=str_to_Double(".3184575132827459839479607674412890463333683258832841739031e-3");
    rs_remainder[32][38]=str_to_Double("-.516710517915097047884322698033227862170645066064945120532e-4");
    rs_remainder[32][39]=str_to_Double("-.63809394430278441286130763453350793788829770311252593146155e-4");
    rs_remainder[32][40]=str_to_Double(".1875034823540121369487633091697254965472378088995615013002e-4");
    rs_remainder[32][41]=str_to_Double(".973818838131934148474016451756302324952292119843055392948794e-5");
    rs_remainder[32][42]=str_to_Double("-.395218698694916040746000195302907910246847025845013422621614610e-5");
    rs_remainder[32][43]=str_to_Double("-.1193312989167783582935602586204353880157404994855319812837647144e-5");
    rs_remainder[32][44]=str_to_Double(".6073402957299170461384947895628282410339894320812514008838256934e-6");
    rs_remainder[32][45]=str_to_Double(".1224445514872259970160106044824738674639026037757173779336213323e-6");
    rs_remainder[32][46]=str_to_Double("-.73415901735413244889413240052575586674875502946960687850192e-7");
    rs_remainder[32][47]=str_to_Double("-.109145670963344020474742347793330078629924639762103882035095e-7");
    rs_remainder[32][48]=str_to_Double(".7251111527890871331064883331201677887218930840813858278023789196e-8");
    rs_remainder[32][49]=str_to_Double(".8709116062412869829446172157679317541176020493513270392294945535e-9");
    rs_remainder[32][50]=str_to_Double("-.5987030838251272501683353652386919039091212405578573416691160698e-9");
    rs_remainder[32][51]=str_to_Double("-.6335558328253288839473418112083843765807765991637327610245129933e-10");
    rs_remainder[32][52]=str_to_Double(".41964236929639848730597424525192122450416203701770327800132422e-10");
    rs_remainder[32][53]=str_to_Double(".4219761410414471162561837094237783722643113147760563165271817251e-11");
    rs_remainder[32][54]=str_to_Double("-.2524750725659850590206308184364404379573288735095813916668752178e-11");
    rs_remainder[32][55]=str_to_Double("-.2557725676099747095039533386781691372936030670794641203445025367e-12");
    rs_remainder[32][56]=str_to_Double(".1314867901055554201676385301115853344985361803731792118563146728e-12");
    rs_remainder[32][57]=str_to_Double(".1397569616448280609350994304028335886148399935439854884844658782e-13");
    rs_remainder[32][58]=str_to_Double("-.5966768969002699526827505923270672324148480018764173199454025726e-14");
    rs_remainder[32][59]=str_to_Double("-.6831133117621255001901551021417112082206766648409127292809950168e-15");
    rs_remainder[32][60]=str_to_Double(".2371763938270738115467525092921438188960045881053104286469762826e-15");
    rs_remainder[32][61]=str_to_Double(".2974408231539370440427585264467833003347497546207617543179942466e-16");
    rs_remainder[32][62]=str_to_Double("-.8292014979401536940587607014117665869275086896905717445145549729e-17");
    rs_remainder[32][63]=str_to_Double("-.1152630607522122184091640882519331687143249692081218011563257208e-17");
    rs_remainder[32][64]=str_to_Double(".2557350651219766921513480941673485292957183310263759086699664235e-18");
    rs_remainder[32][65]=str_to_Double(".3980425508852350389115608677476698910636758121376177910808025109e-19");
    rs_remainder[32][66]=str_to_Double("-.69690326316277625540215971454248794343491019805739192871120578e-20");
    rs_remainder[32][67]=str_to_Double("-.1228174351488758729429805516554669904205241233591727595227789871e-20");
    rs_remainder[32][68]=str_to_Double(".1677825191293790384147342978680173929948186724309539301585989803e-21");
    rs_remainder[32][69]=str_to_Double(".3397105832022883666030606372448403500983700170350953363916820926e-22");
    rs_remainder[32][70]=str_to_Double("-.3558545505986153241490217104385267768198535936459938631671213186e-23");
    rs_remainder[32][71]=str_to_Double("-.8452925715718129945873737672634615430495319974755126470666598e-24");
    //============================= n = 33 ================================================
    rs_remainder[33][0]=str_to_Double("-.125926719114624740658416628978415423231488731977999779126940e-8");
    rs_remainder[33][1]=str_to_Double(".436548012790143619477923190772865084035426551818806478578235338e-8");
    rs_remainder[33][2]=str_to_Double("-.201985653862133132187744789980311687770962827958164230e-8");
    rs_remainder[33][3]=str_to_Double("-.50638594760121512300983086087052577779197294097000524e-7");
    rs_remainder[33][4]=str_to_Double(".176679352530001779318044696415385479329140010744179411666506074e-6");
    rs_remainder[33][5]=str_to_Double("-.852024513363722833500509733684468427563769921062014e-7");
    rs_remainder[33][6]=str_to_Double("-.10570480203117726041625714455200176237523645534167109e-5");
    rs_remainder[33][7]=str_to_Double(".41971014329312824363623006358508126769208113017910357e-5");
    rs_remainder[33][8]=str_to_Double("-.9253091279390813159478076877748913927223419268435814e-5");
    rs_remainder[33][9]=str_to_Double(".1411591381517039044592608472044417728733091858211142562004276342e-4");
    rs_remainder[33][10]=str_to_Double("-.1551139408781567709034319011049023798466599599704437e-4");
    rs_remainder[33][11]=str_to_Double(".157344171792508948734245371329117654823715793227844666e-4");
    rs_remainder[33][12]=str_to_Double("-.4422583878049932772541470229913990812652836703973854945437156e-4");
    rs_remainder[33][13]=str_to_Double(".205470589203909252852097853201790970616327142747989620937e-3");
    rs_remainder[33][14]=str_to_Double("-.7478342886142753537979425758601324030290873241669457e-3");
    rs_remainder[33][15]=str_to_Double(".20873932056724450042470694061340823124924217799631934e-2");
    rs_remainder[33][16]=str_to_Double("-.4622490402446878297888307244808088411189438118214781e-2");
    rs_remainder[33][17]=str_to_Double(".81563350717870472984640987130874158477714932980343410e-2");
    rs_remainder[33][18]=str_to_Double("-.110582117749695351473515469003189287590415415602368921e-1");
    rs_remainder[33][19]=str_to_Double(".100397882091743968296817690204098314573135454016169849487052e-1");
    rs_remainder[33][20]=str_to_Double("-.19352007332245909527325916832278682545458923905374328261309825e-2");
    rs_remainder[33][21]=str_to_Double("-.123468155056394721438350805712457585865725353751112991621e-1");
    rs_remainder[33][22]=str_to_Double(".249369772272217118552021632165670574857780088386308574e-1");
    rs_remainder[33][23]=str_to_Double("-.245924356841294924176522866829949559748526475253128047812266e-1");
    rs_remainder[33][24]=str_to_Double(".73975641902389964676107449545768247977553553741424940e-2");
    rs_remainder[33][25]=str_to_Double(".1531776435172920414736482421633788906698275902586522577438e-1");
    rs_remainder[33][26]=str_to_Double("-.2446268670743109363984628872173921778340486478828711835036250682e-1");
    rs_remainder[33][27]=str_to_Double(".1271837131601461573387288435018359095891610399907517412837e-1");
    rs_remainder[33][28]=str_to_Double(".63173053112818589483507926585622911128069263687590737206214e-2");
    rs_remainder[33][29]=str_to_Double("-.13479211898741183917721535747508961926331530134133144844e-1");
    rs_remainder[33][30]=str_to_Double(".5792498152439406009986006824017904905074109239135152324914657949e-2");
    rs_remainder[33][31]=str_to_Double(".36170581009307295009483157579558324266337631251522569945e-2");
    rs_remainder[33][32]=str_to_Double("-.4728534046928303007160926305075278461537926350189249682985774776e-2");
    rs_remainder[33][33]=str_to_Double(".631282903380971717328429520036817160204394710832766587946418e-3");
    rs_remainder[33][34]=str_to_Double(".1588055610880092176074362123842061178515712721225080396199320809e-2");
    rs_remainder[33][35]=str_to_Double("-.7616663217025650756277674165554581797545416242119356330837106075e-3");
    rs_remainder[33][36]=str_to_Double("-.27974972526222387685035994666793415946976882254904035894757e-3");
    rs_remainder[33][37]=str_to_Double(".2839188662255712570557720951667463436428305091440012436881864737e-3");
    rs_remainder[33][38]=str_to_Double(".132084097038973463766613351291637863885662592638768596189e-4");
    rs_remainder[33][39]=str_to_Double("-.6747199678011714049068401018453757993189472492703817956524262895e-4");
    rs_remainder[33][40]=str_to_Double(".6833007444468299332003543246184922842164702625642683158643188898e-5");
    rs_remainder[33][41]=str_to_Double(".118355123062889428427225674653316023705610688175467063552676e-4");
    rs_remainder[33][42]=str_to_Double("-.2296853101270335918434934126232592403719755159610341158272324712e-5");
    rs_remainder[33][43]=str_to_Double("-.1640789053953247137625609004125298491025927097691313122178888784e-5");
    rs_remainder[33][44]=str_to_Double(".4255549828775411418202976323997178814196713786487711465802e-6");
    rs_remainder[33][45]=str_to_Double(".18766035934290890749521038078851100421097642398673783079713e-6");
    rs_remainder[33][46]=str_to_Double("-.5717386790155846627767161159066003434047654652854279862713982327e-7");
    rs_remainder[33][47]=str_to_Double("-.182556437020518451870148556452632659032372865467499572470119e-7");
    rs_remainder[33][48]=str_to_Double(".6045326290059829772527190646858904863095363014849053266135771628e-8");
    rs_remainder[33][49]=str_to_Double(".15438238177579383763112899618867428792752685722291874033657e-8");
    rs_remainder[33][50]=str_to_Double("-.5229962512180735605756295995713262427270295476826754822206477182e-9");
    rs_remainder[33][51]=str_to_Double("-.1151296485567821317545678275804832767287039918483619272473395468e-9");
    rs_remainder[33][52]=str_to_Double(".3787237160144407536843327368365386675371048402788428678244575e-10");
    rs_remainder[33][53]=str_to_Double(".7631429792251888685939017945357155050660993785787751431110235957e-11");
    rs_remainder[33][54]=str_to_Double("-.2330287760837134258943115931886260098020335565174910956328209726e-11");
    rs_remainder[33][55]=str_to_Double("-.451124892506581198905429969963395869698256910651100957877271e-12");
    rs_remainder[33][56]=str_to_Double(".1231294489747481908739528427758090381936906263717466909916838918e-12");
    rs_remainder[33][57]=str_to_Double(".2380218827798766814216220156448905145101589922777509431212161e-13");
    rs_remainder[33][58]=str_to_Double("-.5630351472579213587510776646707322295535429156538796655371898e-14");
    rs_remainder[33][59]=str_to_Double("-.11209768474407806671196824874787682132746973538496304202165453e-14");
    rs_remainder[33][60]=str_to_Double(".2240530980410989174539692157969057515742876145851197172335729629e-15");
    rs_remainder[33][61]=str_to_Double(".4714053921711639827452043572335509033506099626272929678946549751e-16");
    rs_remainder[33][62]=str_to_Double("-.7787459150376290424390142783324407649043026486574846449899858208e-17");
    rs_remainder[33][63]=str_to_Double("-.17720398809575411110270004460972438480881724546073387711296840e-17");
    rs_remainder[33][64]=str_to_Double(".2367839953510128157944299069342620260283290300092387234191456689e-18");
    rs_remainder[33][65]=str_to_Double(".5964912073739579652822030936997832231094404973814479943948507037e-19");
    rs_remainder[33][66]=str_to_Double("-.6290718418280956504112403276102540560042409832333613839446968433e-20");
    rs_remainder[33][67]=str_to_Double("-.18021522173753032077317142209957831994411371812675047874720470e-20");
    rs_remainder[33][68]=str_to_Double(".1452089066998944426619676911337627432597500528078239877795678580e-21");
    rs_remainder[33][69]=str_to_Double(".4899881405839185848523463363529035671079850757046812021730880757e-22");
    rs_remainder[33][70]=str_to_Double("-.2871403392572686740389304934942983703793054863747855631494668e-23");
    rs_remainder[33][71]=str_to_Double("-.1202264010876669331848807758975893953720726538900954607908779299e-23");
    //============================= n = 34 ================================================
    rs_remainder[34][0]=str_to_Double("-.192507231769022901607049398067425369275785684792639292291e-9");
    rs_remainder[34][1]=str_to_Double(".24477378772098052767905152841660309184171338310919376716e-8");
    rs_remainder[34][2]=str_to_Double("-.1430942884844988674677514955718778438212066029697588395e-7");
    rs_remainder[34][3]=str_to_Double(".503317916278879590473494613231065944052787499078577353302307654e-7");
    rs_remainder[34][4]=str_to_Double("-.11395926757466201073168073875838696035596987320046377059995e-6");
    rs_remainder[34][5]=str_to_Double(".149235406389371645055626869381629692830895423934743033578e-6");
    rs_remainder[34][6]=str_to_Double(".1332475575371975406789999605275517298152173617361361e-6");
    rs_remainder[34][7]=str_to_Double("-.16162359999108245325379310462232301213865461116686402e-5");
    rs_remainder[34][8]=str_to_Double(".6280904541679258682318625809950710659598233621705273390e-5");
    rs_remainder[34][9]=str_to_Double("-.17843441155752344925874897474702885004341575246831503874e-4");
    rs_remainder[34][10]=str_to_Double(".428852214269339068488921972741796100989984239548050e-4");
    rs_remainder[34][11]=str_to_Double("-.9430244521525331698169267338698786639571328788062787e-4");
    rs_remainder[34][12]=str_to_Double(".20141539499258192712932960864005203861802188951180057001e-3");
    rs_remainder[34][13]=str_to_Double("-.4314548372835344655986924241303291127589758357287449533411194e-3");
    rs_remainder[34][14]=str_to_Double(".917899524293692212613978127038714593647425251344715355365e-3");
    rs_remainder[34][15]=str_to_Double("-.18623911524150193601522089464650242712845530553514995935e-2");
    rs_remainder[34][16]=str_to_Double(".3430311626250643151937034549668067756555812948763654776678e-2");
    rs_remainder[34][17]=str_to_Double("-.544477394323610920622931802897989694409527329991702204e-2");
    rs_remainder[34][18]=str_to_Double(".69211506087727742063989726428003944591844881054178643e-2");
    rs_remainder[34][19]=str_to_Double("-.58656235855722837512239960457468226010152685028423823e-2");
    rs_remainder[34][20]=str_to_Double(".1417675787972644039095392748512302820101346832551852e-3");
    rs_remainder[34][21]=str_to_Double(".10193188523582010407877163293028275429366327581477573364e-1");
    rs_remainder[34][22]=str_to_Double("-.20498618086387646898428217037022889781808624467205577666e-1");
    rs_remainder[34][23]=str_to_Double(".22531183065536442257040973193431910714428498555948760305e-1");
    rs_remainder[34][24]=str_to_Double("-.109009008249132834563419498931196561140456524599194040344e-1");
    rs_remainder[34][25]=str_to_Double("-.92519326691696350078133838158826740923654262306468228247936023e-2");
    rs_remainder[34][26]=str_to_Double(".230151438500352916346851324976032185336928944331119683e-1");
    rs_remainder[34][27]=str_to_Double("-.186677749077141959304432253652125648043657491196023523099e-1");
    rs_remainder[34][28]=str_to_Double(".9281654759707502145098026770181538929737835310955246613386211525e-3");
    rs_remainder[34][29]=str_to_Double(".12733684225166562713716190706202261758862629204973174856e-1");
    rs_remainder[34][30]=str_to_Double("-.11112758445051059536125928921975649965260780935314927426e-1");
    rs_remainder[34][31]=str_to_Double(".543841222164105258023864867428958198888244842906313982e-3");
    rs_remainder[34][32]=str_to_Double(".55885715840458566140908040055832936609820928342791717114e-2");
    rs_remainder[34][33]=str_to_Double("-.3378265124027403083958302960789389423213305573277741450454666212e-2");
    rs_remainder[34][34]=str_to_Double("-.791836574435879185063987745295845994190170284995699291828e-3");
    rs_remainder[34][35]=str_to_Double(".1654322179038962354564745117654728189319905956593345231837997736e-2");
    rs_remainder[34][36]=str_to_Double("-.31464349365753171555086351198012559279828968206549646330e-3");
    rs_remainder[34][37]=str_to_Double("-.4298618687863616819742673295356822312864014190876699550531393775e-3");
    rs_remainder[34][38]=str_to_Double(".2022045329741114603037224051669371671722721831235219715153768130e-3");
    rs_remainder[34][39]=str_to_Double(".6668380493519698652706672085385676686290191722521027775398700213e-4");
    rs_remainder[34][40]=str_to_Double("-.5952580270341832763863095544931711263689387672717799505905e-4");
    rs_remainder[34][41]=str_to_Double("-.4854082899513782937775025319644008956143547054328933456752518435e-5");
    rs_remainder[34][42]=str_to_Double(".11979931097972635080840775263159763288374594971600354549641e-4");
    rs_remainder[34][43]=str_to_Double("-.4376874766882802803481783840913934047013381844657671107583e-6");
    rs_remainder[34][44]=str_to_Double("-.1842105882520892407820561780839657859619289472204434531859242346e-5");
    rs_remainder[34][45]=str_to_Double(".19479515388070376139708792041573128756809109186699212400340e-6");
    rs_remainder[34][46]=str_to_Double(".2286226419234801776200449363605169328110532677628714383376210097e-6");
    rs_remainder[34][47]=str_to_Double("-.33919369382252271501715090659435375501903023354519816285228e-7");
    rs_remainder[34][48]=str_to_Double("-.2368897253318819090487763044230155968137055193383448877535628e-7");
    rs_remainder[34][49]=str_to_Double(".40898690484900789529169654620514095317963390430772737207549e-8");
    rs_remainder[34][50]=str_to_Double(".209622601116613859704529921553133809074326837704131919970155e-8");
    rs_remainder[34][51]=str_to_Double("-.3828329649900837927971334009676885701147956228236595862202679160e-9");
    rs_remainder[34][52]=str_to_Double("-.16087294806647996811406324498828926830626770121047347462836e-9");
    rs_remainder[34][53]=str_to_Double(".291519376102995984533901669198440611286244622370108046169825e-10");
    rs_remainder[34][54]=str_to_Double(".108176434433545748307625169502462902315920267871987859718608e-10");
    rs_remainder[34][55]=str_to_Double("-.1851578386921345655658992190424978878694725418040905938446405e-11");
    rs_remainder[34][56]=str_to_Double("-.64163901588595803880760830412054531622018183293734744057954e-12");
    rs_remainder[34][57]=str_to_Double(".995908968269835988044583503605314231897555794396694593995128e-13");
    rs_remainder[34][58]=str_to_Double(".3371998401717026317441482670753450971532532567817615653114718343e-13");
    rs_remainder[34][59]=str_to_Double("-.4579666352861368239297840627055051145984210708636484957972966e-14");
    rs_remainder[34][60]=str_to_Double("-.1575163037949826401617480037426739140281202608702951349520149e-14");
    rs_remainder[34][61]=str_to_Double(".181019295441879060125600598820983674420132038639058068277280e-15");
    rs_remainder[34][62]=str_to_Double(".6558166436978693341345104490813791881645168666124205568247271234e-16");
    rs_remainder[34][63]=str_to_Double("-.6158678079978684059777199620699593784513684497502956681748496712e-17");
    rs_remainder[34][64]=str_to_Double("-.2439935723394217701325159810614270833351878972044739827098973898e-17");
    rs_remainder[34][65]=str_to_Double(".1796412954279067685001687833232372154577752377272419541105518323e-18");
    rs_remainder[34][66]=str_to_Double(".8132958012717572145345794119295605249831050657485440298973811e-19");
    rs_remainder[34][67]=str_to_Double("-.4433050234386129920776785743479895859153604758872883411261426045e-20");
    rs_remainder[34][68]=str_to_Double("-.243531776989835789842034989569723906680670735346816863705697912e-20");
    rs_remainder[34][69]=str_to_Double(".8937839357220338344161212592118276947283215082219928987844814e-22");
    rs_remainder[34][70]=str_to_Double(".6568599547935580946617644954194771703400262222241176994965728e-22");
    rs_remainder[34][71]=str_to_Double("-.1324720386287066568480808575805574526476601119499694988595978e-23");
}

void initialize_rs_remainder8(){

    //============================= n = 35 ================================================
    rs_remainder[35][0]=str_to_Double("-.1292502087742916993549366499803021229889837230066613826119091950e-8");
    rs_remainder[35][1]=str_to_Double(".641730009743498734426475028435026465095989428653253899e-8");
    rs_remainder[35][2]=str_to_Double("-.1633426586542647315763820544748318400187010031279587494e-7");
    rs_remainder[35][3]=str_to_Double(".202692095306419587635855679112792347762851752215235e-8");
    rs_remainder[35][4]=str_to_Double(".744963768004172559301061256522084403605608353550567e-7");
    rs_remainder[35][5]=str_to_Double("-.13417173743636431051261286509917646120720801639119830e-6");
    rs_remainder[35][6]=str_to_Double(".74610596651129710341889152001059887988490112183066e-7");
    rs_remainder[35][7]=str_to_Double("-.57649007524731443441373954012146147777682555004194e-6");
    rs_remainder[35][8]=str_to_Double(".46292204087644875054842063055474081278632933001616e-5");
    rs_remainder[35][9]=str_to_Double("-.192495265079354906150708567531313097250097062272590e-4");
    rs_remainder[35][10]=str_to_Double(".57505421063189676890194098599288571074748864242409e-4");
    rs_remainder[35][11]=str_to_Double("-.1426037895408896739916747887664974331406711985082444331584931e-3");
    rs_remainder[35][12]=str_to_Double(".315948500114834211289991435465108709851838488197043288e-3");
    rs_remainder[35][13]=str_to_Double("-.6460073701070750206243977670209347357541739959901163772e-3");
    rs_remainder[35][14]=str_to_Double(".12210202826014557089323225851838101548013849346873738561e-2");
    rs_remainder[35][15]=str_to_Double("-.20872421859586621019955740739311945093283911684916866e-2");
    rs_remainder[35][16]=str_to_Double(".3086499093975314981736446778074860384736759785448263e-2");
    rs_remainder[35][17]=str_to_Double("-.3608828149986586501509989820845952577560384764430636898880586932e-2");
    rs_remainder[35][18]=str_to_Double(".2465318913337226443730833729293854760466031916310859035e-2");
    rs_remainder[35][19]=str_to_Double(".167970513365194932184325153047032434912097218352396497167e-2");
    rs_remainder[35][20]=str_to_Double("-.89770480923311418037677756052365755293874298967793503e-2");
    rs_remainder[35][21]=str_to_Double(".167202243755623846247125141373817356464835430130905113402359e-1");
    rs_remainder[35][22]=str_to_Double("-.192784395167259479691625170789547488794138356389282058707262e-1");
    rs_remainder[35][23]=str_to_Double(".11625215523009870855776420947407062130778899339797295311702e-1");
    rs_remainder[35][24]=str_to_Double(".478908675043073780987666450207053668093479529894993024e-2");
    rs_remainder[35][25]=str_to_Double("-.19794369967689850999545196852015142904068481202628584973025237e-1");
    rs_remainder[35][26]=str_to_Double(".213847504685650559538969256329408920509424399611790670e-1");
    rs_remainder[35][27]=str_to_Double("-.77475750585301775734461326745801089715251195218831987716504236e-2");
    rs_remainder[35][28]=str_to_Double("-.891652844874373139782171213092587597100533422150965253e-2");
    rs_remainder[35][29]=str_to_Double(".1405412086118678163606215754866749912253099353963077685330579e-1");
    rs_remainder[35][30]=str_to_Double("-.60143003981288290963365834528547993587455071489816259589422213e-2");
    rs_remainder[35][31]=str_to_Double("-.38130909479702195402093035118895611630481567982780809349e-2");
    rs_remainder[35][32]=str_to_Double(".571817792393477905421812807946247131402042099148140015365e-2");
    rs_remainder[35][33]=str_to_Double("-.1459063245527250104383149400284811534072648158213346728363388614e-2");
    rs_remainder[35][34]=str_to_Double("-.17584662087303921038260838784343098798378624373165007532121201e-2");
    rs_remainder[35][35]=str_to_Double(".13432149957734456512339027787089728031820071265793981322119400e-2");
    rs_remainder[35][36]=str_to_Double(".118741225165073420279244867559437407675839076744399198059656664e-3");
    rs_remainder[35][37]=str_to_Double("-.470022566995775718191316130498761349161327382173532824698e-3");
    rs_remainder[35][38]=str_to_Double(".95658830350395788323145357104020923678815756931447793486e-4");
    rs_remainder[35][39]=str_to_Double(".1006864676191916945835948771064943191063056641177585764566621672e-3");
    rs_remainder[35][40]=str_to_Double("-.426977971174597218774376140012111318019611132227870657536e-4");
    rs_remainder[35][41]=str_to_Double("-.144251224841879419065807371963080918981924214195235847930e-4");
    rs_remainder[35][42]=str_to_Double(".1031353663526552678941299905721835734398565065875151157748325e-4");
    rs_remainder[35][43]=str_to_Double(".133203567793218023200911939941287594624197058760335736889786e-5");
    rs_remainder[35][44]=str_to_Double("-.177981158136396934691291851078895568667507203266525220320499e-5");
    rs_remainder[35][45]=str_to_Double("-.51145485873971101451312878364468572847383754064546629914658e-7");
    rs_remainder[35][46]=str_to_Double(".2399742385339462399369802448119642441136327239207841288618999239e-6");
    rs_remainder[35][47]=str_to_Double("-.6663587253159332459371924907833334054838755196113733980011460301e-8");
    rs_remainder[35][48]=str_to_Double("-.2645205743644837761474826968390134196852400586974627433433e-7");
    rs_remainder[35][49]=str_to_Double(".15951428773659853638772500146666950015054384282283270987935e-8");
    rs_remainder[35][50]=str_to_Double(".2450595740273498139792841426900154493471281981125895268823481675e-8");
    rs_remainder[35][51]=str_to_Double("-.1895798053379443084817212624629417531988675634611186752372291706e-9");
    rs_remainder[35][52]=str_to_Double("-.1943471430835000129348553724282492778062381671156458254256079905e-9");
    rs_remainder[35][53]=str_to_Double(".1624207408915209363249551950724071271622337970533800482640780280e-10");
    rs_remainder[35][54]=str_to_Double(".1336249276563213557146364255981653898293898829552593370401084125e-10");
    rs_remainder[35][55]=str_to_Double("-.109637519867373687726378209383011396311618263008530402437222e-11");
    rs_remainder[35][56]=str_to_Double("-.8037212807993510293062146431120207241883628624847659089150166e-12");
    rs_remainder[35][57]=str_to_Double(".604041464058744187588393744516364690080492135042856659150946e-13");
    rs_remainder[35][58]=str_to_Double(".42569358823412796298809001902455371984149854480890352126225e-13");
    rs_remainder[35][59]=str_to_Double("-.2756541860149444160399060407495685151345759206150065226027935353e-14");
    rs_remainder[35][60]=str_to_Double("-.1995595480647666115131921947555896647842799143577653318333162808e-14");
    rs_remainder[35][61]=str_to_Double(".1043081456020475947255080690003202882187168247549714275805155862e-15");
    rs_remainder[35][62]=str_to_Double(".831468195874329923415086922338578548745791511552356547711856e-16");
    rs_remainder[35][63]=str_to_Double("-.322033009574519564481582260634868853669474608578100469735242e-17");
    rs_remainder[35][64]=str_to_Double("-.3090265757215819063242892682873070803342667946255868946865556e-17");
    rs_remainder[35][65]=str_to_Double(".7672351531707803872705468878782004551392191552498293679088892045e-19");
    rs_remainder[35][66]=str_to_Double(".1027917472460753508621360573765907954931060215443440823017104059e-18");
    rs_remainder[35][67]=str_to_Double("-.11279708329669477811216386128826928262835057977261525784719025e-20");
    rs_remainder[35][68]=str_to_Double("-.30694909669154152593008584542713115588983335000617481478727345e-20");
    rs_remainder[35][69]=str_to_Double("-.8066824588457113721997681009173594089094676953203515725514923525e-23");
    rs_remainder[35][70]=str_to_Double(".8252093717299973591468771278614848767912940061263587810035150e-22");
    rs_remainder[35][71]=str_to_Double(".131439919961467376662644627221558029836192041273811092620793e-23");
    //============================= n = 36 ================================================
    rs_remainder[36][0]=str_to_Double("-.10132462000009898037596299661976315388937480449333272948e-9");
    rs_remainder[36][1]=str_to_Double(".1155511556086084093221073631661089513791901509592081719e-8");
    rs_remainder[36][2]=str_to_Double("-.8162304715914063394191258932232414155675294538096250e-8");
    rs_remainder[36][3]=str_to_Double(".41425972642484227572670376548151514026608482506317614e-7");
    rs_remainder[36][4]=str_to_Double("-.128157995089880013421338384828568545121346399389265182e-6");
    rs_remainder[36][5]=str_to_Double(".2066110176084640495170569978690348940547280016618204e-6");
    rs_remainder[36][6]=str_to_Double(".2390301483078705143288995890270752309690351220346519289489511800e-7");
    rs_remainder[36][7]=str_to_Double("-.92791672796734346601631463567044692768882437616018e-6");
    rs_remainder[36][8]=str_to_Double(".203418691561390623438109870428678379394091240093420e-5");
    rs_remainder[36][9]=str_to_Double("-.544715802429087044204619591061380553631188878683e-6");
    rs_remainder[36][10]=str_to_Double("-.10427963626827935425635881943593479286705705416287e-4");
    rs_remainder[36][11]=str_to_Double(".440898454949894900395597070787407446224332246550196468e-4");
    rs_remainder[36][12]=str_to_Double("-.1227409226006056750591540239029251811098694658318311843307813331e-3");
    rs_remainder[36][13]=str_to_Double(".278619311205175646811693755377815535195979956711964e-3");
    rs_remainder[36][14]=str_to_Double("-.544181536434008147969196597768934156922431678569640e-3");
    rs_remainder[36][15]=str_to_Double(".914234086797791096959546863232973429669097833004352353e-3");
    rs_remainder[36][16]=str_to_Double("-.1259506250667003695701765026472231852341359191104721e-2");
    rs_remainder[36][17]=str_to_Double(".1200384191527786472138618149423216353547407011930250872e-2");
    rs_remainder[36][18]=str_to_Double("-.4215509158228597156495311901595706408421290910601e-4");
    rs_remainder[36][19]=str_to_Double("-.29963461257903609230736357315455561541152340900084605545e-2");
    rs_remainder[36][20]=str_to_Double(".803491825547707714560960875149131281641660202466804403e-2");
    rs_remainder[36][21]=str_to_Double("-.1346363389998895904745768148370846579102718977254057168967787767e-1");
    rs_remainder[36][22]=str_to_Double(".1558943042696484276290848325525839945983216749559900780058553663e-1");
    rs_remainder[36][23]=str_to_Double("-.10434059085646931656994548625061422253740436124965834541e-1");
    rs_remainder[36][24]=str_to_Double("-.2308802509758455944745043793904978628407401825551642517e-2");
    rs_remainder[36][25]=str_to_Double(".1634036990338352247135929793007628357268868109057536695039362070e-1");
    rs_remainder[36][26]=str_to_Double("-.2153515816622867748553018336528611913938878178481977636827262e-1");
    rs_remainder[36][27]=str_to_Double(".126788526670781885198299685910567928074690623218434868e-1");
    rs_remainder[36][28]=str_to_Double(".3761333917739750247592153275873959649666529395004798333319604e-2");
    rs_remainder[36][29]=str_to_Double("-.14116264725582314079161614377638746232585530579523838154e-1");
    rs_remainder[36][30]=str_to_Double(".1089123381700647673040053342311777538455151866161465841878e-1");
    rs_remainder[36][31]=str_to_Double("-.1800368203988721083194823081319135992349081436823269043054482142e-3");
    rs_remainder[36][32]=str_to_Double("-.6274295161092814605598782707312916476104044073477520810616123221e-2");
    rs_remainder[36][33]=str_to_Double(".4336209428773065319193796939881958549286654282486302829e-2");
    rs_remainder[36][34]=str_to_Double(".4480447694097702425336929476924791444918859931818085311802018e-3");
    rs_remainder[36][35]=str_to_Double("-.2129617611465719130526964167471283439059215799646202942692467846e-2");
    rs_remainder[36][36]=str_to_Double(".7920797705217850827409187422565155136915943602000080389182443511e-3");
    rs_remainder[36][37]=str_to_Double(".4499115462653618963562758870929330025480025169117453690205777e-3");
    rs_remainder[36][38]=str_to_Double("-.4101275106205806640427451301356484164892141857131577422107382362e-3");
    rs_remainder[36][39]=str_to_Double("-.11797717657652659921242742912916191269283691306003016007e-4");
    rs_remainder[36][40]=str_to_Double(".1117766161588730536601081323147128749513815851471133532504053063e-3");
    rs_remainder[36][41]=str_to_Double("-.2097418986904333132473996267208415161152135915973259161828e-4");
    rs_remainder[36][42]=str_to_Double("-.20539379383538177319551235641631487026465414978715490041951e-4");
    rs_remainder[36][43]=str_to_Double(".72640556694661459780955020877493743899679437959960289207913e-5");
    rs_remainder[36][44]=str_to_Double(".275218793921852638143532970999357213545284731903562242532e-5");
    rs_remainder[36][45]=str_to_Double("-.14779058342611324745396661348110189356048350718509080586596e-5");
    rs_remainder[36][46]=str_to_Double("-.27719183841324634940413647588586301508171330762958959100220e-6");
    rs_remainder[36][47]=str_to_Double(".2207526981458772200417105718208890640351484647741300541520e-6");
    rs_remainder[36][48]=str_to_Double(".2098995798128551151589123735974009536893668705066043289015084160e-7");
    rs_remainder[36][49]=str_to_Double("-.2613839072911635480824821625529037351965810672653072015178271252e-7");
    rs_remainder[36][50]=str_to_Double("-.1145032357667973115332993361293854715575604949700062187109197581e-8");
    rs_remainder[36][51]=str_to_Double(".2551577528790225391162145208724253668591685345435220477998e-8");
    rs_remainder[36][52]=str_to_Double(".3737303328882582148569813064725188418351455170341350010486e-10");
    rs_remainder[36][53]=str_to_Double("-.2103307277239177188608076869201281505693247913697114349137462174e-9");
    rs_remainder[36][54]=str_to_Double(".183485466470964668004753581114432279996702222016772385157e-12");
    rs_remainder[36][55]=str_to_Double(".1487769053052376275668941786671444116399785050851292632396025169e-10");
    rs_remainder[36][56]=str_to_Double("-.1094707375947689140604358679595009722660512247532051304284e-12");
    rs_remainder[36][57]=str_to_Double("-.9133964343053010744949281029457470821251645354781416392642565310e-12");
    rs_remainder[36][58]=str_to_Double(".7028607436384514416892956098391875889777822108009631073554457359e-14");
    rs_remainder[36][59]=str_to_Double(".4908550627598651964479466254329147853555618363742354581201412145e-13");
    rs_remainder[36][60]=str_to_Double("-.1889548967063625440649047464324355698685556296582396468254e-15");
    rs_remainder[36][61]=str_to_Double("-.2324188536254485765109135148786291385002255234874136295358148404e-14");
    rs_remainder[36][62]=str_to_Double("-.649321213105872558251814142672937580267461444666941818393712e-17");
    rs_remainder[36][63]=str_to_Double(".9748174885615290459224842504808118803979860719201897805021222518e-16");
    rs_remainder[36][64]=str_to_Double(".1098211789754239074143137656624918849899693069556204967725081568e-17");
    rs_remainder[36][65]=str_to_Double("-.3637902192199298029834954827111179327357567200650833836324449679e-17");
    rs_remainder[36][66]=str_to_Double("-.7608254886864097422069375707266964560970253674753642101198398120e-19");
    rs_remainder[36][67]=str_to_Double(".1212648213277055704748332185625332584648969040223693466424186394e-18");
    rs_remainder[36][68]=str_to_Double(".3799263258479334012806459523368945593287909773023591540310798891e-20");
    rs_remainder[36][69]=str_to_Double("-.3622868008480932249779993955625857378189624136554523228661537796e-20");
    rs_remainder[36][70]=str_to_Double("-.1532464804320015065216433416012747171205997963235191033989705903e-21");
    rs_remainder[36][71]=str_to_Double(".9729837965447823053330486492222492393247368086861860139750345e-22");
    //============================= n = 37 ================================================
    rs_remainder[37][0]=str_to_Double("-.125851385698293493154239831934711924089106965048647193665e-8");
    rs_remainder[37][1]=str_to_Double(".723413209513262243031562961448571467837924696097763827e-8");
    rs_remainder[37][2]=str_to_Double("-.26217002755728744694559032214534998998265934033622887e-7");
    rs_remainder[37][3]=str_to_Double(".5067305116009696252520078016164459076272182709591897e-7");
    rs_remainder[37][4]=str_to_Double("-.3685136470868087991275987824019378542688333710150866937e-7");
    rs_remainder[37][5]=str_to_Double("-.10178907108550096489497582416109392272997940478780641011e-6");
    rs_remainder[37][6]=str_to_Double(".56528866367656664387721833323987922911278966027989315e-6");
    rs_remainder[37][7]=str_to_Double("-.17553561017349953381161073752470898612813758301689617e-5");
    rs_remainder[37][8]=str_to_Double(".41102788484809374361786380883990662248056838879808e-5");
    rs_remainder[37][9]=str_to_Double("-.72838549664714485510641734433777336942584122553346e-5");
    rs_remainder[37][10]=str_to_Double(".8485138175260253638033755614788104336867602429108e-5");
    rs_remainder[37][11]=str_to_Double("-.75264299780813131150401142856107574476221026692546290e-6");
    rs_remainder[37][12]=str_to_Double("-.263226792318731777504934505809705588360790059895817e-4");
    rs_remainder[37][13]=str_to_Double(".7853303009577642558295422347492589899596644689866e-4");
    rs_remainder[37][14]=str_to_Double("-.13395382038335912363433043363882478307654418644918e-3");
    rs_remainder[37][15]=str_to_Double(".9716023864825293157561288352294640642141243395209e-4");
    rs_remainder[37][16]=str_to_Double(".2624873125403860504395313350787931116718069888966810021757663e-3");
    rs_remainder[37][17]=str_to_Double("-.133673575291817025308205778152607805027541683865266e-2");
    rs_remainder[37][18]=str_to_Double(".35482563325180041693442327244197263466540291502984613774e-2");
    rs_remainder[37][19]=str_to_Double("-.6931287755567078281511219926995997360736575867759958e-2");
    rs_remainder[37][20]=str_to_Double(".1048658349753150297218090593474868153049932175191249489438967087e-1");
    rs_remainder[37][21]=str_to_Double("-.118418645470644860716725610007831460749767052129164076e-1");
    rs_remainder[37][22]=str_to_Double(".8129822462333736967966580061359811500655976492012616108486338438e-2");
    rs_remainder[37][23]=str_to_Double(".1540836685326333153816949438716797035550399515625425859308791483e-2");
    rs_remainder[37][24]=str_to_Double("-.135167229033874024934978525666852673965058890714290512e-1");
    rs_remainder[37][25]=str_to_Double(".201661142025764651173210634639318073140549856626795742e-1");
    rs_remainder[37][26]=str_to_Double("-.153455125622248249586680333098491610800741418667572174809245e-1");
    rs_remainder[37][27]=str_to_Double(".1104648081811323824963049745416832026416964849249789e-2");
    rs_remainder[37][28]=str_to_Double(".120120618721044809175326449236879977541161672353296235e-1");
    rs_remainder[37][29]=str_to_Double("-.139394603503197651463348908372607519896395323129467588933e-1");
    rs_remainder[37][30]=str_to_Double(".505748241046162368692750587837849424692601286826073206e-2");
    rs_remainder[37][31]=str_to_Double(".458861525086563389579832217217447768182655587061864896e-2");
    rs_remainder[37][32]=str_to_Double("-.6546990393871301217081134636563432792127509002095002772970620501e-2");
    rs_remainder[37][33]=str_to_Double(".209535501384878704459966550045852030530220699553347963002285111e-2");
    rs_remainder[37][34]=str_to_Double(".1877319487488733446767504257137854709491004147144946952889528888e-2");
    rs_remainder[37][35]=str_to_Double("-.1930969049118912553211487426878819010038173706576739813e-2");
    rs_remainder[37][36]=str_to_Double(".1668680297932751148120524385802003565824449781305200986e-3");
    rs_remainder[37][37]=str_to_Double(".62696700493374599239748721473947405115483657288558114092e-3");
    rs_remainder[37][38]=str_to_Double("-.2786822345098771423874000039102521385673885710737849688564838505e-3");
    rs_remainder[37][39]=str_to_Double("-.9958854262235954952243459397729022943843920607596458470e-4");
    rs_remainder[37][40]=str_to_Double(".101201837927929439292967323811658700716051419556295505338e-3");
    rs_remainder[37][41]=str_to_Double(".1269143848558429443932630377735727982593301306190958077369389901e-5");
    rs_remainder[37][42]=str_to_Double("-.2258418552184992063601313088004731701330408456048437253517e-4");
    rs_remainder[37][43]=str_to_Double(".34481258752632552506596282329086968671034863603451621563e-5");
    rs_remainder[37][44]=str_to_Double(".36368252295762864173648930260177899401884456745466805343268749e-5");
    rs_remainder[37][45]=str_to_Double("-.9958655420531707135690393552328850937525188838658797253815e-6");
    rs_remainder[37][46]=str_to_Double("-.4525794620067569340030107475083139987541972756112628089595e-6");
    rs_remainder[37][47]=str_to_Double(".1744628579177816123871384789190352002343117039361763123016337586e-6");
    rs_remainder[37][48]=str_to_Double(".45393035147962766970833761343938681293390930787351665154448e-7");
    rs_remainder[37][49]=str_to_Double("-.2276379187436989764276082275954038592117431639857185353162e-7");
    rs_remainder[37][50]=str_to_Double("-.3793853202899683704801073349302419411898078903612672113400e-8");
    rs_remainder[37][51]=str_to_Double(".2375100140294395382861835062638951022235455714882657371827e-8");
    rs_remainder[37][52]=str_to_Double(".2728544110815449296986873557397475192634006262409665987996e-9");
    rs_remainder[37][53]=str_to_Double("-.20544295309580414482980952993754314687084975617692116935970e-9");
    rs_remainder[37][54]=str_to_Double("-.17467367294025766968905552852958331782033084681249972668975e-10");
    rs_remainder[37][55]=str_to_Double(".1505985604599013632964187450288514564530685395053521203918649e-10");
    rs_remainder[37][56]=str_to_Double(".102853831085985797577033067884709667081147357325724983726124e-11");
    rs_remainder[37][57]=str_to_Double("-.9495491258428853160933637541525411989710320121813171808744432572e-12");
    rs_remainder[37][58]=str_to_Double("-.5702656079116909594646703784673005587619234061347543499244322240e-13");
    rs_remainder[37][59]=str_to_Double(".5205198807346364438667498760492174321053536905828253563972277706e-13");
    rs_remainder[37][60]=str_to_Double(".2995381642866958053626159379200222223415319890038410449567891908e-14");
    rs_remainder[37][61]=str_to_Double("-.2501014394017916847468235032893151681488149855762453607552905e-14");
    rs_remainder[37][62]=str_to_Double("-.147594164404227102027373176728540008075902154152053682396504e-15");
    rs_remainder[37][63]=str_to_Double(".10601034573212118837610657761596670405276056743645442094778022e-15");
    rs_remainder[37][64]=str_to_Double(".671205983015869271120750891771775730772269139702527357964791e-17");
    rs_remainder[37][65]=str_to_Double("-.398484343704824551085561702308502139028444229604615816796184e-17");
    rs_remainder[37][66]=str_to_Double("-.277787661836161139554607506958932531372730918406077098065980e-18");
    rs_remainder[37][67]=str_to_Double(".1334137179378878836912240345804279904881150100802365470747858925e-18");
    rs_remainder[37][68]=str_to_Double(".1037369126154267642766288792066639345449250886296724921585713425e-19");
    rs_remainder[37][69]=str_to_Double("-.3993041015873985830930103035836306727028419532452359783692260e-20");
    rs_remainder[37][70]=str_to_Double("-.3483245150196653566732190492238741459157636526300114245555738610e-21");
    rs_remainder[37][71]=str_to_Double(".10716058685466784477874045354534536574714782541345109122559706e-21");
    //============================= n = 38 ================================================
    rs_remainder[38][0]=str_to_Double("-.296996907574416393195689857563448230170778900700911078e-10");
    rs_remainder[38][1]=str_to_Double("-.112745753647832484272955476447929239298178059773502e-10");
    rs_remainder[38][2]=str_to_Double("-.98225887232324539938696098130459941720140154558318966750e-9");
    rs_remainder[38][3]=str_to_Double(".1981445573911599511604386856752084095014996412891824e-7");
    rs_remainder[38][4]=str_to_Double("-.97804315474444568383372743325048557144254615250087863e-7");
    rs_remainder[38][5]=str_to_Double(".23435145181532491725191309999641240088572760642327e-6");
    rs_remainder[38][6]=str_to_Double("-.303727061847196054398225419093342246740998705103590e-6");
    rs_remainder[38][7]=str_to_Double(".2822360670176317737375114078056400699519930587176e-6");
    rs_remainder[38][8]=str_to_Double("-.937003739745800098632175141163975383171460060093e-6");
    rs_remainder[38][9]=str_to_Double(".439690039322264447983074052494468383010497151633404526145e-5");
    rs_remainder[38][10]=str_to_Double("-.1377607703404132267314841578649863570701038371104250368e-4");
    rs_remainder[38][11]=str_to_Double(".320898658009562524385120090609455901603702866042064565e-4");
    rs_remainder[38][12]=str_to_Double("-.624368316283337044788653251241310136415150887420663e-4");
    rs_remainder[38][13]=str_to_Double(".112626562941701420479323519997729619545057374671045128e-3");
    rs_remainder[38][14]=str_to_Double("-.2078162118612907004255110057110797963217973694645124e-3");
    rs_remainder[38][15]=str_to_Double(".41348484292810300570163294846525208771279435190792e-3");
    rs_remainder[38][16]=str_to_Double("-.8634236514876258405563213944584446271196124336472e-3");
    rs_remainder[38][17]=str_to_Double(".17653940326866943545020077830965838850045629901219e-2");
    rs_remainder[38][18]=str_to_Double("-.332253088804441564418729744215070309603828187274468746761e-2");
    rs_remainder[38][19]=str_to_Double(".54974732744391561688692301079206601047856952167186e-2");
    rs_remainder[38][20]=str_to_Double("-.7637990379672731333403946048181516580025105439926282305e-2");
    rs_remainder[38][21]=str_to_Double(".8233317270323187673433148709072058710931245047295336e-2");
    rs_remainder[38][22]=str_to_Double("-.53398836116923567759711308169379021223835282458552753e-2");
    rs_remainder[38][23]=str_to_Double("-.19423721628489630068109159101163413843545643071231894e-2");
    rs_remainder[38][24]=str_to_Double(".11562326883992980700338284630616070156115373799653382e-1");
    rs_remainder[38][25]=str_to_Double("-.181504853279387312126930971831798937269544802054490594351037e-1");
    rs_remainder[38][26]=str_to_Double(".160626688112084654835483068571854625756612489942230628e-1");
    rs_remainder[38][27]=str_to_Double("-.469788713872277583253282594347435104963605827510674589562e-2");
    rs_remainder[38][28]=str_to_Double("-.8946663865749716621313891160603124868802644189526180164753441e-2");
    rs_remainder[38][29]=str_to_Double(".149326008811376885325179045140630344760577582584057743e-1");
    rs_remainder[38][30]=str_to_Double("-.94474261924269298505800254209853976756923446972520682914529832e-2");
    rs_remainder[38][31]=str_to_Double("-.120818323872237853920250603763379938843480797931644893e-2");
    rs_remainder[38][32]=str_to_Double(".7148237521810424226672233636888743691927733855041861868506826420e-2");
    rs_remainder[38][33]=str_to_Double("-.49841193952899480987392959050896042060106346821241657544e-2");
    rs_remainder[38][34]=str_to_Double("-.259857478124183933286485221763785375543330802474364076e-3");
    rs_remainder[38][35]=str_to_Double(".2569579281218615099118380886970572688883017748847526288508990015e-2");
    rs_remainder[38][36]=str_to_Double("-.13161673985034461070761093436728277479605134993943498960e-2");
    rs_remainder[38][37]=str_to_Double("-.3799402944304608037279623633988156696074971141823544223114929e-3");
    rs_remainder[38][38]=str_to_Double(".63957314677726041660949260796233726696770825876790682959e-3");
    rs_remainder[38][39]=str_to_Double("-.11404399102334503161229322662532659243254601778617988040e-3");
    rs_remainder[38][40]=str_to_Double("-.15403332348750878832698910397960455294666050686323238910462580e-3");
    rs_remainder[38][41]=str_to_Double(".7404143851202069791900115959707285940733349609177709194714248087e-4");
    rs_remainder[38][42]=str_to_Double(".20109319450911984560856684759914049642020812769167535405e-4");
    rs_remainder[38][43]=str_to_Double("-.20707068677310643751978682289617101918785294155831779737e-4");
    rs_remainder[38][44]=str_to_Double("-.4505602789937159911203283813759937472299979577596679295178614665e-6");
    rs_remainder[38][45]=str_to_Double(".38995074210272498020183790854925497335219365740784241647289e-5");
    rs_remainder[38][46]=str_to_Double("-.4170199199247730183235793226989805150682206356944732293564525176e-6");
    rs_remainder[38][47]=str_to_Double("-.5555080406025840889961349188598872436884210540190291978458552605e-6");
    rs_remainder[38][48]=str_to_Double(".1084240437133445013806770127635201210285095721237302891674887114e-6");
    rs_remainder[38][49]=str_to_Double(".634021500043245369794240735611638909522366467653744681647709489e-7");
    rs_remainder[38][50]=str_to_Double("-.1676699425559193473571722360439500567853293897171951856267697361e-7");
    rs_remainder[38][51]=str_to_Double("-.6019467538129650221680118355734453751756823163423392425543385458e-8");
    rs_remainder[38][52]=str_to_Double(".1933064679780747871529611386599485658002406358075298540969887470e-8");
    rs_remainder[38][53]=str_to_Double(".48914368053227757903932427414413782670489765312618206349084e-9");
    rs_remainder[38][54]=str_to_Double("-.178749626215524080280245133060696319358695506444660616664119e-9");
    rs_remainder[38][55]=str_to_Double("-.3480746665850025606132371603959327422896328998767866963869930112e-10");
    rs_remainder[38][56]=str_to_Double(".1374106276015850311842633262705326206093853842831597650886116350e-10");
    rs_remainder[38][57]=str_to_Double(".22078412438207115051665828853404133743545602272327441321726077e-11");
    rs_remainder[38][58]=str_to_Double("-.8971762474819894607209820140397491980897413011211539588093985e-12");
    rs_remainder[38][59]=str_to_Double("-.1263562660051197002752781112135482426308941477726034042848523513e-12");
    rs_remainder[38][60]=str_to_Double(".5047218844915619726801736284497973183519992927877731813011657786e-13");
    rs_remainder[38][61]=str_to_Double(".65672821025772105873250672130709044889215186943270033095701e-14");
    rs_remainder[38][62]=str_to_Double("-.24719313062433392635222907064512560000170841744278776487873633e-14");
    rs_remainder[38][63]=str_to_Double("-.3105724918872145458408291528545222675336975379900401963748720770e-15");
    rs_remainder[38][64]=str_to_Double(".1062259726754226241850369286070259149795703040447782433048439319e-15");
    rs_remainder[38][65]=str_to_Double(".1335233760008465932349410693277589690319029857333511560776335222e-16");
    rs_remainder[38][66]=str_to_Double("-.402985912742462175391941476719206595521045961593237948001527e-17");
    rs_remainder[38][67]=str_to_Double("-.520951344341450185711778430523683469418486433139837938998893e-18");
    rs_remainder[38][68]=str_to_Double(".1356172474381895237210627475501200298443596582498849312295554920e-18");
    rs_remainder[38][69]=str_to_Double(".184185455774248742241934062357560670163754919444160744316974e-19");
    rs_remainder[38][70]=str_to_Double("-.406398437959004672844768488595041792280242304759312016055421e-20");
    rs_remainder[38][71]=str_to_Double("-.58978358819780283507272193033787430031941975708539580628385940e-21");
    //============================= n = 39 ================================================
    rs_remainder[39][0]=str_to_Double("-.11778155140561865022474007360136291627105871723913099577e-8");
    rs_remainder[39][1]=str_to_Double(".692027087134460916398783687941073258193572150555823359e-8");
    rs_remainder[39][2]=str_to_Double("-.2986989108450549079739526731419396960975656914236427e-7");
    rs_remainder[39][3]=str_to_Double(".836586536194096702689354712828102738246455618073368731e-7");
    rs_remainder[39][4]=str_to_Double("-.14668080051901231583433397372678832691076373179689e-6");
    rs_remainder[39][5]=str_to_Double(".9200830329912496103568775424327306872780434697228e-7");
    rs_remainder[39][6]=str_to_Double(".307883774489220238722563762528060051893002168040470e-6");
    rs_remainder[39][7]=str_to_Double("-.9849701362125522789867573292713854778852124027680e-6");
    rs_remainder[39][8]=str_to_Double(".751871732852145746187277102358556392158073836957332e-6");
    rs_remainder[39][9]=str_to_Double(".3461414698603466964660938941920325242061962306604652895e-5");
    rs_remainder[39][10]=str_to_Double("-.170037420934183889733377728314661303723512014970020211e-4");
    rs_remainder[39][11]=str_to_Double(".481451343878592584155460119432766811525090253472e-4");
    rs_remainder[39][12]=str_to_Double("-.1103575089371650600814843541033746152030648486133e-3");
    rs_remainder[39][13]=str_to_Double(".2277897019651994372604127418642250869046499098367e-3");
    rs_remainder[39][14]=str_to_Double("-.4446696533239209900643816985137643475579535057910e-3");
    rs_remainder[39][15]=str_to_Double(".835328087222736362634334995134883785644527623023834e-3");
    rs_remainder[39][16]=str_to_Double("-.1500466904035690324718665395719229819961074427091400638797154e-2");
    rs_remainder[39][17]=str_to_Double(".2518002562430753672697653651971302480029387203731120901671303346e-2");
    rs_remainder[39][18]=str_to_Double("-.38103760942360726703834898944105692786126582958672838997657e-2");
    rs_remainder[39][19]=str_to_Double(".49344906938361640960994554349787455315673940950101e-2");
    rs_remainder[39][20]=str_to_Double("-.49317858562891215669325173468924639140181128711540e-2");
    rs_remainder[39][21]=str_to_Double(".2543416997312145334244115828378069432550342089869209181805971487e-2");
    rs_remainder[39][22]=str_to_Double(".29256858959835764451255500934323975002780584354778e-2");
    rs_remainder[39][23]=str_to_Double("-.10313087424881196100117803183532863818599587303008849080451336e-1");
    rs_remainder[39][24]=str_to_Double(".15986766580153608214613899127972374265300328437235738557003990e-1");
    rs_remainder[39][25]=str_to_Double("-.153778939206319432926443170032848033441630287912911066e-1");
    rs_remainder[39][26]=str_to_Double(".6696676764457837620011565155665195225023728762692630600e-2");
    rs_remainder[39][27]=str_to_Double(".599509732690700644159959916369874868628951571419812561e-2");
    rs_remainder[39][28]=str_to_Double("-.14378497044690106971887834310165217541799761355946389404e-1");
    rs_remainder[39][29]=str_to_Double(".12533996631268443310289550293616707058741268314760791e-1");
    rs_remainder[39][30]=str_to_Double("-.2816061601763462294371801880564965406922754594867505536464868814e-2");
    rs_remainder[39][31]=str_to_Double("-.5976619681322330031485591961472290481434275008063039187e-2");
    rs_remainder[39][32]=str_to_Double(".719904249260227264978553394449393387637277517420831040e-2");
    rs_remainder[39][33]=str_to_Double("-.2361019324518937577591048754621999846769719479464246618728632073e-2");
    rs_remainder[39][34]=str_to_Double("-.2100879342762256852842229651775148556850541285423767198882e-2");
    rs_remainder[39][35]=str_to_Double(".24954953962296603868777846831918588529685872474370215881272e-2");
    rs_remainder[39][36]=str_to_Double("-.507121457428014567310906254772925113181703369831945629e-3");
    rs_remainder[39][37]=str_to_Double("-.740901570904398281604618293591836001941727295753288784e-3");
    rs_remainder[39][38]=str_to_Double(".5140854724734694676531544306288668946827430574229123485e-3");
    rs_remainder[39][39]=str_to_Double(".4499789490093588224413899655851102958567023738510447224e-4");
    rs_remainder[39][40]=str_to_Double("-.1700068706873858661978128947334613227085960890782703927e-3");
    rs_remainder[39][41]=str_to_Double(".3772561002165384793141456971583179737132147890016043861080581280e-4");
    rs_remainder[39][42]=str_to_Double(".32733807759295544091753137473689191539769434832071041438479e-4");
    rs_remainder[39][43]=str_to_Double("-.15706827785155691437477218284238329957832603894942656476532e-4");
    rs_remainder[39][44]=str_to_Double("-.380376058444530486427524577843937150849786055605922613970204628e-5");
    rs_remainder[39][45]=str_to_Double(".3559600178504534606237860049060601272885887290571397510874626e-5");
    rs_remainder[39][46]=str_to_Double(".166343748081451267061909755196808409587488185986759679965892e-6");
    rs_remainder[39][47]=str_to_Double("-.576105985066215019547935745501731321227859550149335626807881e-6");
    rs_remainder[39][48]=str_to_Double(".3252621081316829358552781980497612050710653375961536532473e-7");
    rs_remainder[39][49]=str_to_Double(".72838969756085890943677598446381139529662419391026702779921e-7");
    rs_remainder[39][50]=str_to_Double("-.894471158012397040176941983171678819458926027984720764403e-8");
    rs_remainder[39][51]=str_to_Double("-.7542783086105660464530686915538334152353176832367389088957e-8");
    rs_remainder[39][52]=str_to_Double(".127275248911813516123097569749338400480904672688990988328798e-8");
    rs_remainder[39][53]=str_to_Double(".65972218786431097673367851255325034521011909878647860134677e-9");
    rs_remainder[39][54]=str_to_Double("-.132045008592477696175989014857845059751497544647800297673236e-9");
    rs_remainder[39][55]=str_to_Double("-.49819454779832463500086744760840305255303095699168807112280e-10");
    rs_remainder[39][56]=str_to_Double(".1092436991818715174851455374969163877835012849190688602963594837e-10");
    rs_remainder[39][57]=str_to_Double(".3300879752787209080499724498060011956364167599869040210643932026e-11");
    rs_remainder[39][58]=str_to_Double("-.7501609120501517136069537083427382647491446763924230874397799666e-12");
    rs_remainder[39][59]=str_to_Double("-.1941025548169847670647518793927538383114099147285218326075558200e-12");
    rs_remainder[39][60]=str_to_Double(".4373571342331071628990799375278702553705021809871830340364885655e-13");
    rs_remainder[39][61]=str_to_Double(".1020831348297672606758700038170116352058053064496748541956310e-13");
    rs_remainder[39][62]=str_to_Double("-.2196790858790791998117390346130956762817719838088158028439740296e-14");
    rs_remainder[39][63]=str_to_Double("-.4825094062659631974190326649090474489562921644856941078303852969e-15");
    rs_remainder[39][64]=str_to_Double(".9603521573149251674581173490452659348821967456905224237390416788e-16");
    rs_remainder[39][65]=str_to_Double(".2055690326532900394328444493427546492903346149122492719982472896e-16");
    rs_remainder[39][66]=str_to_Double("-.3680947128084808614243036284567174515531390404480921273153066490e-17");
    rs_remainder[39][67]=str_to_Double("-.7908857863603367870760944668510531748142430100966951553487348952e-18");
    rs_remainder[39][68]=str_to_Double(".1243633802298846027386934973503112788387184422811429762193369482e-18");
    rs_remainder[39][69]=str_to_Double(".27516361556212980827726595717536428006663251644786773407977e-19");
    rs_remainder[39][70]=str_to_Double("-.3717123510874455798419265274206844830432066907187659752445829e-20");
    rs_remainder[39][71]=str_to_Double("-.8669473557997118151171302540625103328996445123780531622876493384e-21");
}

