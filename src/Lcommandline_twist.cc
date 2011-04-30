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



#include "Lcommandline_twist.h"


/*===========================================================================

Dictionary of twists:

degree 1:

   L(s,chi), chi primitive of conductor q, has functional equation

   OMEGA = tau(chi)/sqrt(q) if chi(-1)= 1
   OMEGA = tau(chi)/(I*sqrt(q)) if chi(-1)= -1

   Q=sqrt(q/Pi)

   Gamma(s/2) if chi(-1)= 1
   Gamma((s+1)/2) if chi(-1)= -1

   When chi = kronecker's symbol, OMEGA is always equal to one

   reference: davenport chapter 9. For some reason he has an (Pi/q)^(-(s+a)/2)
   in(13), pg 71 rather than a (Pi/q)^(-s/2). The extra constant factor makes
   no difference, though.

degree 2:

   --------------------------------------------------------------------
   L(s,f) a cusp form of weight k and level N, i.e. in S_k(Gamma_0(N)).

   Then, (see for example Knapp, Theorem 9.8 page 270) Hecke obtained the
   functional equation for L(s,f) with

   underlying OMEGA = i^k * epsilon = plus or minus one.
   underlying Q = sqrt(N)/(2*Pi).

   epsilon is the eigenvalue of the Fricke involution = 1 or -1

   assume that the conductor of chi, q, satisfies (q,N)=1, and that
   chi is primitive.

   XXXXXXXXXXXXXXXXX
   What happens for imprimitive chi's when we twist
   a cusp form? What happens when (q,N) > 1.
   XXXXXXXXXXXXXXXXx

   Then the functional for L(s,f,chi) has

   Q=q*sqrt(N)/(2Pi), i.e. underlying Q (of L(s,f)) gets multiplied by q.
   OMEGA = (underlying OMEGA) chi(N) tau(chi)^2 / q

   when chi = kronecker's sympbol,
   OMEGA = (underlying OMEGA) chi(-N)


   references: Shimura, theorem 3.66, Knapp
   Bump, Koblitz.

   --------------------------------------------------------------------

   L-functions of Maass forms for SL_2(Z)

   underlying Q = 1/Pi
   OMEGA = 1

   twisting by primitive chi gives

   Q = q/Pi
   OMEGA = 1.

   XXXXXX
   What happens for higher level and for imprimitive chi?



   --------------------------------------------------------------------

higher degree:

   see what's known.



=============================================================================*/

//if what_to_do="values" then count and step_size are not used
//and we compute L(s=x+I*y, chi_d).
//if what_to_do="zeros" then count and step_size are used
//and we compute zeros of L(1/2+I*t, chi_d) by looking
//for sign changes in steps of step_size. We start at
//t1=x and go up to t2=y, or, if count >0 we look for
//the first count zeros above x.


//int quadratic_twists(Long D1, Long D2,Double x,Double y,Long count,Long start_N,Double step_size,const char *what_to_do,bool do_only_even_twists,bool test_explicit_formula,int desired_rank)
int quadratic_twists(Long D1, Long D2,Double x,Double y,Long count,Long start_N,Double step_size,const char *what_to_do,bool do_only_even_twists,int desired_rank)
{


    int k,n;

    char message_stamp[300];
    ostringstream os2;


    Complex s,u;

    Long  d;
    char tmp_name[300];
    ostringstream os3;



    int *coeff_chi;

    Double y2=0; //used in determining how many twisted coefficients are needed


    int tmp_a; //the quasi degree of the twisted (and also underlying) L-function

    if(current_L_type==1){
        tmp_a=int_L.a;
    }
    if(current_L_type==2){
        tmp_a=Double_L.a;
    }
    if(current_L_type==3){
        tmp_a=Complex_L.a;
    }

    Double *tmp_g;
    Complex *tmp_l;
    tmp_g=new Double[tmp_a+1];
    tmp_l=new Complex[tmp_a+1];

    Complex underlying_OMEGA; // OMEGA of the underlying L-function
    Complex tmp_OMEGA; // OMEGA of the twisted L-function
    int what_type; // what_type of the underlying L-function
    int tmp_what_type; // what_type of the twisted L-function

    int underlying_N_terms; //number of dirichlet coefficients
                            //of the underlying L-function

    Double underlying_Q; //Q for the underlying L-function

    Long underlying_cond;

    //first copy the GAMMA factors of the underlying L-function
    //afterwards we modify... first according to a dictionary
    //of twists and if this doesn't work then by searching...
    //latter still needs to be figured out... non-linear system to
    //consider.
    if(current_L_type==1){
        //tmp_a=int_L.a;
        for(k=1;k<=tmp_a;k++){
             tmp_g[k]=int_L.gamma[k];
             tmp_l[k]=int_L.lambda[k];
             y2=y2+real(tmp_l[k]);
        }
        underlying_OMEGA=int_L.OMEGA;
        what_type=int_L.what_type_L;
        underlying_N_terms=int_L.number_of_dirichlet_coefficients;
        underlying_Q=int_L.Q;
    }
    if(current_L_type==2){
        //tmp_a=Double_L.a;
        for(k=1;k<=tmp_a;k++){
             tmp_g[k]=Double_L.gamma[k];
             tmp_l[k]=Double_L.lambda[k];
             y2=y2+real(tmp_l[k]);
        }
        underlying_OMEGA=Double_L.OMEGA;
        what_type=Double_L.what_type_L;
        underlying_N_terms=Double_L.number_of_dirichlet_coefficients;
        underlying_Q=Double_L.Q;
    }
    if(current_L_type==3){
        //tmp_a=Complex_L.a;
        for(k=1;k<=tmp_a;k++){
             tmp_g[k]=Complex_L.gamma[k];
             tmp_l[k]=Complex_L.lambda[k];
             y2=y2+real(tmp_l[k]);
        }
        underlying_OMEGA=Complex_L.OMEGA;
        what_type=Complex_L.what_type_L;
        underlying_N_terms=Complex_L.number_of_dirichlet_coefficients;
        underlying_Q=Complex_L.Q;
    }

    Complex *tmp_cmplx; //used to set poles and residues to none
                        //since twists never have poles.
    tmp_cmplx=new Complex[1]; tmp_cmplx[0]=0;


    Long r; // r = |d|
    Long tmp_period;

    int N_terms;// number of chi(n)'s to compute
    int tmp_N_terms;// number of b(n)*chi(n)'s to compute

    Double tmp_Q;

    int *int_twisted_coeffs;
    Double *Double_twisted_coeffs;
    Complex *Complex_twisted_coeffs;

    if(what_type==2)
    underlying_cond=Long(rint(double(underlying_Q*2*Pi*underlying_Q*2*Pi)));

    else underlying_cond=1;

/*
    //precompute a table of the kronecker symbol up to 1000
    if(!kronecker_table_available&&abs((D1-D2)*sqrt(abs(D1)+abs(D2)))>1e6){
        kronecker_bound=10000;
        if(my_verbose>1) cout << "precomputing a table of the kronecker symbol (m|n), n< " << kronecker_bound << endl;
        kronecker_table = new int *[kronecker_bound];
        for(int k=0;k<kronecker_bound;k++) kronecker_table[k]= new int[kronecker_bound];

        for(int j=0;j<kronecker_bound;j++)
        for(int k=0;k<kronecker_bound;k++) kronecker_table[j][k]=my_kronecker(j,k);
        kronecker_table_available=true;

    }

*/

    d=nextfunddiscriminant(D1-1);
    for(; d<=D2; d=nextfunddiscriminant(d))
    if(what_type!=2||(what_type==2&&gcd(underlying_cond,sn(d)*d)==1))
    {

         tmp_what_type=what_type;
         if(tmp_what_type==-1)tmp_what_type=-2; //this gets updated to 1 if N_terms=r;



        tmp_OMEGA=underlying_OMEGA;

        tmp_period=0;
        if (d<0) r=-d; else r=d;

        os3.str("");
        //os3.seekp(0);
        if(d<0)
        //os3 <<"L_chi_-"<< r << ends;
        os3 <<"L_chi_-"<< r;
        else
        //os3 <<"L_chi_"<< r << ends;
        os3 <<"L_chi_"<< r;
        strcpy(tmp_name,os3.str().c_str());

        //Double C; // as in the 2*C/Pi

        if(what_type==-1) // -1 is for zeta
        {
            //if(my_kronecker(d,d-1)==1) tmp_l[1]=0.;
            //this is -1 for d<0 and 1 for d>0, so ... no need for the call to kronecker
            if(d>0) tmp_l[1]=0;
            else tmp_l[1]=(Double)1/2;
            tmp_Q=sqrt(double(r)/Pi); //XXXXXXXXXXXXX Might this cause trouble in higher precision.
                                      //I did this as a quick way to fix overload not defined for long long / mpreal

            //C = 2*DIGITS2*log(10.);
            Double T;

            if(!strcmp(what_to_do,"zeros"))
                T = abs(max(abs(x),abs(y)));
            else //if what_to_do=="values"
                T = abs(y);

            if(count>0)
                T=T+(count+100.)*2*Pi/log(T+3);

            Complex delta=int_L.find_delta(1+I*T+tmp_l[1],.5);
            tmp_N_terms = Int(2.3 * DIGITS*tmp_Q/real(delta));
            do{
                tmp_N_terms=(int)(tmp_N_terms*1.3);
            }while(tmp_N_terms*real(delta)/tmp_Q-log((1.+real(tmp_l[1]))*tmp_N_terms)<2.3*DIGITS);
            //}while(tmp_N_terms*real(delta)/tmp_Q-log((1.+real(tmp_l[1]))*tmp_N_terms)<2.3*DIGITS);

            tmp_N_terms=(int)(tmp_N_terms*1.3);

//tmp_N_terms/=10; //XXXXXXXXXXXXXXXXXX krtek

            N_terms=tmp_N_terms; // for Dirichlet L-functions
            if(r<=1000)N_terms=r;

            if(my_verbose>1){
                cout << "delta = " << delta << " tmp_Q = " << tmp_Q << endl;
                cout << "N_terms for quad twist: " << N_terms << endl;
            }

        }


        // If is a degree > 1 L-function
        if(what_type!=-1)
        {


           //for now only cusp forms are implemented
           if(what_type==0||what_type>3||what_type==1||what_type==-2){
               cout << endl << "At present only know how to twist zeta, ";
               cout << "cusp form (in S_k(Gamma_0(N)) L-functions," << endl;
               cout << "and Maass forms for the full modular group" << endl;
               delete [] tmp_g;
               delete [] tmp_l;
               delete [] tmp_cmplx;
               return 0;
           }

           // if a twist of a cusp form
           if(what_type==2||what_type==3){

               tmp_Q=underlying_Q*double(r); //XXXXXXXXXXXXXXXXXXXXXX another overload mpreal and long long fix XXXXX check in higher precision

               if(what_type==2){
                   tmp_OMEGA =tmp_OMEGA*my_kronecker(d,underlying_cond);
                   if(d<0) tmp_OMEGA=-tmp_OMEGA;
               }

               //C = DIGITS2*log(10.);
               Double T;


               if(!strcmp(what_to_do,"zeros"))
                   T = abs(max(abs(x),abs(y)));
               else //if what_to_do=="values"
                   T = abs(y);

               if(count>0)
               T=T+(count+100.)*Pi/log(T+3);

               Complex delta=int_L.find_delta(1+I*T+tmp_l[1],1);
               tmp_N_terms = Int(2.3 * DIGITS*tmp_Q/real(delta));
               do{
                   tmp_N_terms=(int)(tmp_N_terms*1.3);
               }while(tmp_N_terms*real(delta)/tmp_Q-log((1.+real(tmp_l[1]))*tmp_N_terms)<2.3*DIGITS);
               tmp_N_terms=(int)(tmp_N_terms*1.3);


           }

           if(tmp_N_terms>underlying_N_terms)
               tmp_N_terms=underlying_N_terms;

           N_terms=tmp_N_terms;

        } //if is a degree >1 L-function

        if (N_terms>=r)
        {
            N_terms=r;
            if(what_type==-1){
                tmp_N_terms=r;
                tmp_what_type=1; //(i.e. is a Dirichlet series with all coeffs)
                tmp_period=r;
            }
        }


    //if do_only_even_twists and is an even twist or !do_only_even_twists
    if((do_only_even_twists&&abs(tmp_OMEGA-1)<1.e-12)||!do_only_even_twists)
    {

            coeff_chi=new int[N_terms+1];
            coeff_chi[0]=0;
            for(n=1;n<=N_terms;n++){
                if(d<INT_MAX&&d>-INT_MAX)
                    coeff_chi[n]=my_kronecker(Int(d),n);
                else
                    coeff_chi[n]=my_kronecker(d,(Long) n);
            }


            if(current_L_type==1){
                if(what_type==-1){ //i.e. if twists of zeta
                    if (my_verbose>1) cout << "calling L_function for zeta twist with " << tmp_N_terms << " coefficients" << endl;
                    int_L2= L_function<int>(tmp_name,tmp_what_type,tmp_N_terms,coeff_chi,tmp_period,tmp_Q,tmp_OMEGA,tmp_a,tmp_g,tmp_l,0,tmp_cmplx,tmp_cmplx);
                }
                else{
                    int_twisted_coeffs= new int[tmp_N_terms+1]; //xxxxxxxxxxxx moved here
                    for(n=1;n<=tmp_N_terms;n++) int_twisted_coeffs[n]=
                        int_L.dirichlet_coefficient[n]*coeff_chi[n%r];
                    if (my_verbose>1) cout << "calling L_function for twist with " << tmp_N_terms << " coefficients" << endl;
                    int_L2= L_function<int>(tmp_name,tmp_what_type,tmp_N_terms,int_twisted_coeffs,tmp_period,tmp_Q,tmp_OMEGA,tmp_a,tmp_g,tmp_l,0,tmp_cmplx,tmp_cmplx);
                }
            }
            if(current_L_type==2){
                Double_twisted_coeffs= new Double[tmp_N_terms+1];
                for(n=1;n<=tmp_N_terms;n++) Double_twisted_coeffs[n]=
                    Double_L.dirichlet_coefficient[n]*coeff_chi[n%r];
                if (my_verbose>1) cout << "calling L_function for twist with " << tmp_N_terms << " coefficients" << endl;
                Double_L2= L_function<Double>(tmp_name,tmp_what_type,tmp_N_terms,Double_twisted_coeffs,tmp_period,tmp_Q,tmp_OMEGA,tmp_a,tmp_g,tmp_l,0,tmp_cmplx,tmp_cmplx);
            }
            if(current_L_type==3){
                Complex_twisted_coeffs= new Complex[tmp_N_terms+1];
                for(n=1;n<=tmp_N_terms;n++) Complex_twisted_coeffs[n]=
                    Complex_L.dirichlet_coefficient[n]*coeff_chi[n%r];
                if (my_verbose>1) cout << "calling L_function for twist with " << tmp_N_terms << " coefficients" << endl;
                Complex_L2= L_function<Complex>(tmp_name,tmp_what_type,tmp_N_terms,Complex_twisted_coeffs,tmp_period,tmp_Q,tmp_OMEGA,tmp_a,tmp_g,tmp_l,0,tmp_cmplx,tmp_cmplx);
            }

            os2.str("");
            //os2.seekp(0);
            if(d<0)
            os2 <<"-"<< r ;
            else
            os2 << r ;

            strcpy(message_stamp,os2.str().c_str());

            //cout << "stamp="<<message_stamp<<endl;

    //---------------------find zeros or compute values...
            if(!strcmp(what_to_do,"zeros")||!strcmp(what_to_do,"zeros and ranks"))
            {
                bool do_zeros=true;
                if(!strcmp(what_to_do,"zeros and ranks")){
                    int analytic_rank=0;
                    if(abs(u)<.000001){
                        switch(current_L_type)
                                {
                                    case 1:
                                        analytic_rank=int_L2.compute_rank();
                                        break;
                                    case 2:
                                        analytic_rank=Double_L2.compute_rank();
                                        break;
                                    case 3:
                                        analytic_rank=Complex_L2.compute_rank();
                                        break;
                                }
                    }
                    if(analytic_rank!=desired_rank){
                        do_zeros=false;
                    }


                }
                if(do_zeros)
                {
                    switch(current_L_type)
                    {
                        case 1:
                            if(count==0)
                            int_L2.find_zeros(x,y,step_size,"cout",message_stamp);
                            else
                            int_L2.find_zeros(count,start_N,step_size,-1,message_stamp);
                            //int_L2.find_zeros(count,start_N,step_size,-1,test_explicit_formula,message_stamp);
                            break;
                        case 2:
                            if(count==0)
                            Double_L2.find_zeros(x,y,step_size,"cout",message_stamp);
                            else
                            Double_L2.find_zeros(count,start_N,step_size,-1,message_stamp);
                            //Double_L2.find_zeros(count,start_N,step_size,-1,test_explicit_formula,message_stamp);
                            break;
                        case 3:
                            if(count==0)
                            Complex_L2.find_zeros(x,y,step_size,"cout",message_stamp);
                            else
                            Complex_L2.find_zeros(count,start_N,step_size,-1,message_stamp);
                            //Complex_L2.find_zeros(count,start_N,step_size,-1,test_explicit_formula,message_stamp);
                            break;
                    }


                }
            }
            if(!strcmp(what_to_do,"values")||!strcmp(what_to_do,"values and ranks"))
            {
                s=x+I*y;
                switch(current_L_type)
                {
                    case 1:
                        u=int_L2.value(s,global_derivative,"pure");
                        break;
                    case 2:
                        u=Double_L2.value(s,global_derivative,"pure");
                        break;
                    case 3:
                        u=Complex_L2.value(s,global_derivative,"pure");
                        break;
                }


                if(desired_rank==-1){ //default... -1 means not to restrict the rank
                    //cout << setprecision(DIGITS3);
                    cout << message_stamp << " " << real(u) << " " << imag(u);
                }

                if(!strcmp(what_to_do,"values and ranks")){
                    int analytic_rank;
                    if(abs(u)>.000001 && desired_rank==-1) cout << " 0";
                    else{
                        switch(current_L_type)
                        {
                            case 1:
                                analytic_rank=int_L2.compute_rank();
                                break;
                            case 2:
                                analytic_rank=Double_L2.compute_rank();
                                break;
                            case 3:
                                analytic_rank=Complex_L2.compute_rank();
                                break;
                        }
                        if(analytic_rank==desired_rank){ //restrict output to those with specified rank
                            //cout << setprecision(DIGITS3);
                            cout << message_stamp << " " << real(u) << " " << imag(u) << endl;
                        }
                        else if(desired_rank==-1) cout << " " << analytic_rank;
                    }


                }
                if(desired_rank==-1) cout << endl;
                //cout << message_stamp << " " << abs(u) << endl;

                //cout << d << " " << real(u) << " vs finite sum formula" << endl<< L_1_chi(d) <<" " << L_1_chi(d) * sqrt(abs(Double(d)))/(2.L*Pi) << endl;
            }
    //--------------------------------------------------
            delete [] coeff_chi;
            if(current_L_type==1){
                if(what_type!=-1)delete [] int_twisted_coeffs;
            }
            if(current_L_type==2){
                delete [] Double_twisted_coeffs;
            }
            if(current_L_type==3){
                delete [] Complex_twisted_coeffs;
            }

        }//if do_only_even_twists

    } // for d
    delete [] tmp_g;
    delete [] tmp_l;
    delete [] tmp_cmplx;


    return 0;

}

//note: the twists in all_twists are all treated as Complex L-functions,
//including the twists that are real.
//for elliptic curves we check if (q,N)=1.

//twist_type: -1 just one complex primitive twist, -2 all non-real primitive twists, 0 all primitive twists, 1 all twists with
//conjugate pairs only appearing once, 2 all twists with conjugate pairs appearing twice
//int all_twists(Long D1, Long D2,Double x,Double y,Long count, Long start_N,Double step_size,const char *what_to_do,int twist_type,int print_character,bool test_explicit_formula)
int all_twists(Long D1, Long D2,Double x,Double y,Long count, Long start_N,Double step_size,const char *what_to_do,int twist_type,int print_character)
{

    //variables for twisting
    Long q;   //conductor
    Long q_1; //conductor of the inducing character chi_1
    int n; //as in chi[n]
    int K; //number of distinct odd primes dividing q

    int i,j,k; //loop variables
    int loop1,loop2; //looping bounds
    int a,b; //various ints

    Complex u,w,z;  //used to generate the character

    Complex *chi;   //the character
    Complex *chi_1; //the primitive character inducing chi
    Complex *sum;   //to check that sum over chi is zero
    Complex SUM;    //to check that sum over n is zero
    //Double total=0;  //sum of abs(SUM) over all chi
    Long **factors;  //stores the factors of q

    Long *m, *v;    //as in (1) chapter 4 of Davenport
    Long m_prime, v_prime;
    Long *phi; //stores the values of phi(p^alpha)
    Long *power_p;// stores p_j^alpha_j

    Long PHI; //is set to phi(q/2^alpha_0)
    Long M; //for looping 0..PHI-1
    Long dial; //for creating an "analog dial" looping through the m[j]'s

    Long gr; //generator for the cyclic group of (Z/p^alpha Z)

    bool is_primitive; //for storing whether chi is primitive or not
    bool only_do_primitive; // only do primitive twists or not
    bool is_fresh; //false if we already did the conjugate character 
                   //and current_L_type is not Complex

    factors = new Long *[30];
    for(k=0;k<=29;k++) factors[k]= new Long[3];

    //----------------------------------------

    char message_stamp[300];
    ostringstream os2;


    Complex s;

    Long d;
    const char *tmp_name="L_chi";

    i=0;


    int tmp_a; //the quasi degree of the twisted (and also underlying) L-function
    Double y2=0; //used in determining how many twisted coefficients are needed

    if(current_L_type==1){
        tmp_a=int_L.a;
    }
    if(current_L_type==2){
        tmp_a=Double_L.a;
    }
    if(current_L_type==3){
        tmp_a=Complex_L.a;
    }

    Double *tmp_g;
    Complex *tmp_l;
    tmp_g=new Double[tmp_a+1];
    tmp_l=new Complex[tmp_a+1];

    Complex underlying_OMEGA; // OMEGA of the underlying L-function
    Complex tmp_OMEGA; // OMEGA of the twisted L-function
    int what_type; // what_type of the underlying L-function
    int tmp_what_type; // what_type of the twisted L-function

    int underlying_N_terms; //number of dirichlet coefficients
                            //of the underlying L-function

    Double underlying_Q; //Q for the underlying L-function

    //first copy the GAMMA factors of the underlying L-function
    //afterwards we modify... first according to a dictionary
    //of twists and if this doesn't work then by searching...
    //latter still needs to be figured out... non-linear system to
    //consider.
    if(current_L_type==1){
        //tmp_a=int_L.a;
        for(k=1;k<=tmp_a;k++){
             tmp_g[k]=int_L.gamma[k];
             tmp_l[k]=int_L.lambda[k];
             y2=y2+real(tmp_l[k]);
        }
        underlying_OMEGA=int_L.OMEGA;
        what_type=int_L.what_type_L;
        underlying_N_terms=int_L.number_of_dirichlet_coefficients;
        underlying_Q=int_L.Q;
    }
    if(current_L_type==2){
        //tmp_a=Double_L.a;
        for(k=1;k<=tmp_a;k++){
             tmp_g[k]=Double_L.gamma[k];
             tmp_l[k]=Double_L.lambda[k];
             y2=y2+real(tmp_l[k]);
        }
        underlying_OMEGA=Double_L.OMEGA;
        what_type=Double_L.what_type_L;
        underlying_N_terms=Double_L.number_of_dirichlet_coefficients;
        underlying_Q=Double_L.Q;
    }
    if(current_L_type==3){
        //tmp_a=Complex_L.a;
        for(k=1;k<=tmp_a;k++){
             tmp_g[k]=Complex_L.gamma[k];
             tmp_l[k]=Complex_L.lambda[k];
             y2=y2+real(tmp_l[k]);
        }
        underlying_OMEGA=Complex_L.OMEGA;
        what_type=Complex_L.what_type_L;
        underlying_N_terms=Complex_L.number_of_dirichlet_coefficients;
        underlying_Q=Complex_L.Q;
    }

    tmp_what_type=what_type;

    Complex *tmp_cmplx; //used to set poles and residues to none
                        //since twists never have poles.
                        //except for the trivial twist of zeta !
                        //but this special case is dealt with seperately
    tmp_cmplx=new Complex[1]; tmp_cmplx[0]=0;


    Long r; // r = |d|
    Long tmp_period;

    int tmp_N_terms;// number of b(n)*chi(n)'s to compute

    Double tmp_Q;
    int number_chi; //counts characters

    Complex *Complex_twisted_coeffs;

    Long underlying_cond; //is set to N, conductor of underlying curve

    if(what_type==2)
        underlying_cond = Long(rint(double(underlying_Q*2*Pi*underlying_Q*2*Pi)));
    else
        underlying_cond=1;

    if(twist_type==1||twist_type==2)
        only_do_primitive=false; // do all twists
    else
        only_do_primitive=true; // only do primitive twists


    for(d=D1; d<=D2; d++)
    if(what_type!=2||(what_type==2&&gcd(underlying_cond,d)==1))
    {

        number_chi=0;
        q=d;
        for(k=0;k<=29;k++)
        for(j=0;j<=2;j++)factors[k][j]=0;
        factor(q,factors);

        K=factors[0][2];

        //for(k=0;k<=K;k++)
            //cout << q <<" " << factors[k][0] << " " << factors[k][1] <<" " << factors[k][2] << endl;

        m = new Long [K+1];
        v = new Long [K+1];
        phi = new Long [K+1];
        power_p = new Long [K+1];
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


        r=d; // a bit redundant since we are using q,r,d for the 
             // same value. But this allows me to mouse Quadratic_twists


// xxxxxxxxxxxx trivial character do not compute...
//since we don't have any poles
// UPDATE APRIL 4 2000 I think this has been taken care of!

//-----------begin character loop over m's
    for(m_prime=0;m_prime<loop1;m_prime++)
    for(m[0]=0;m[0]<loop2;m[0]++)
    for(M=0;M<PHI;M++)
    {

        tmp_OMEGA=underlying_OMEGA;
        tmp_period=0;

        //total=0;
        is_primitive=true;


        for(j=1;j<=K;j++)
            if(gcd(m[j],factors[j][0])>1) is_primitive=false;

        if(factors[0][1]==1) is_primitive=false;
        if(factors[0][1]==2&&m[0]==0) is_primitive=false;
        if(factors[0][1]>2) if(m_prime%2==0) is_primitive=false;

        number_chi++;

 
        is_fresh = true;
        if(current_L_type<=2&&twist_type!=2){
            //i.e. if the dirichlet coeffs are real and we don't insist on
            //conjugate characters appearing twice
            if(K>0){
                j=0;
                do{
                    j++;
                }while((m[j]==phi[j]/2||m[j]==0)&&j<K);
                if(m[j]>phi[j]/2) is_fresh=false;
                else if(j==K&&(m[j]==phi[j]/2||m[j]==0)&&m_prime>(loop1)/2) is_fresh=false;
            }
            else if(m_prime>(loop1)/2) is_fresh=false;

            //if(is_fresh) cout << "FRESH"<< endl; else cout << "is conjugate"<< endl;
            //cout << 0 << " ------  " << loop2 << " " << m[0] << " " << endl;
            //cout << 0 << " ------  " << loop1 << " " << m_prime << " " << endl;
            //for(j=1;j<=K;j++) cout << j << " ------  " << phi[j] << " " << m[j] << endl;
            //false means we already saw it's conjugate
        }


//-----if primitive, or if: not do only primitive (i.e. if do all). We also
//-----make sure we haven't already done the conjugate character when
//-----current_L_type is int or Double

if(is_fresh&&(is_primitive || !only_do_primitive))
{

        // m[j]'s,m_prime are set. now compute chi(n)

        for(n=1;n<=q-1;n++){
            chi[n]=1.;
            chi_1[n]=1.;
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
            gr=factors[k][2]; //generator mod p^alpha
            z=exp(2*(Double)((double)m[k])*I*Pi/(Double)((double)phi[k])); u=1.;
            for(v[k]=1;v[k]<=phi[k];v[k]++)
            {
                u=u*z; //XXXXXXXXXXXX this loses precision if phi[k] is large
                n=(n*gr)%b;
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
            //is_primitive=false;
            q_1=q_1/2;
        }
        if(factors[0][1]==2&&m[0]==1){
            for(n=3;n<=q;n=n+4){
                chi[n]=-chi[n];
                chi_1[n]=-chi_1[n];
            }
        }
        if(factors[0][1]==2&&m[0]==0){
            //is_primitive=false;
            q_1=q_1/4;
        }
        if(factors[0][1]>2){

            if(m_prime%2==0){
                //is_primitive=false;
                if(m[0]==0&&m_prime==0)q_1=q_1/(4*loop1); //strip powers of 2
                else q_1=q_1/gcd(m_prime,loop1);
            }

            n=1;
            if(m[0]==0) w=1.; else w=-1.;
            z=exp(2*Double((double)m_prime)*I*Pi/Double(loop1));
            u=1.;
            b=4*loop1;//=2^alpha_0
            for(v[0]=1;v[0]<=2;v[0]++){
                n=(-n)%b;
                if(n<0)n=n+b;
                u=u*w;
                for(v_prime=1;v_prime<=loop1;v_prime++){
                    n=(5*n)%b;
                    u=u*z;  //XXXXXXXXXXXXXXXXXXXX loses precision if loop1 is large
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

        //cout << "inducing q_1 = " << q_1 << endl;
        //cout << "M = " << M << " ;" ;
        //cout << m_prime << " ";

        //for(n=0;n<=K;n++) cout << m[n] << " ";
        //if(is_primitive) cout << "is primitive" << endl;
        //else cout << "NOT primitive" << endl;

        //end compute chi(n)
        // check that sum over characters is 0;

        //cout << setprecision(DIGITS3);

        SUM=0;
        //if(is_fresh) cout << "fresh==============================\n " ;
        //if(!is_fresh) cout << "old*******************************\n " ;
        for(n=0;n<=q;n++){
            sum[n]=sum[n]+chi[n];
            SUM=SUM+chi[n];

            if(print_character==1&&abs(chi[n])>.1){
                cout << q << " " <<number_chi<<" "<< q_1 << " " << n << " " <<  real(chi[n]) << " " <<imag(chi[n])<< endl;
            }

            //if(print_character==1&&abs(chi[n])>.1&&(is_primitive||number_chi==1))
            //cout << q << " " << 0 << " " <<number_chi<<" "<< n << " " <<  real(chi[n]) << " " <<imag(chi[n])<< endl; //i.e. if primitive or if the trivial character
            //else if(print_character==1&&abs(chi[n])>.1)
            //cout << q << " " << 1 << " " <<number_chi<<" "<< n << " " <<  real(chi[n]) << " " <<imag(chi[n])<< endl;
            //else if(print_character==1&&abs(chi[n])<.1&&abs(chi_1[n])>.1)
            //cout << q << " " << q_1 << " " <<number_chi<<" "<< n << " " <<  real(chi_1[n]) << " " <<imag(chi_1[n])<< endl;
        }

        if(print_character==2&&is_primitive)
            cout << q << " chi_primitive " <<number_chi<<" "<< q-1 << " " <<  real(chi[q-1]) << " " <<imag(chi[q-1])<< endl;
        else if(print_character==2)
            cout << q << " chi " <<number_chi<<" "<< q-1 << " " <<  real(chi[q-1]) << " " <<imag(chi[q-1])<< endl;

        //cout <<"SUM over n = " << abs(SUM) <<endl;
        //total=total+abs(SUM);


        //aha for all_twists I do need to compute the Gauss sum anyways
        //so is an O(q) computation for which I need all q chi(n)'s
        //so O(sqrt(q)) is out for now Dirichlet L-functions...
        //but we could use the fft to do all the gauss sums in O(log(q)) time. 
        //Easiest when //q is prime. Might implement this in the future.

        Complex tau_chi=gauss_sum(chi_1,q_1); // Gauss sum
        //cout << "Gauss sum: " << tau_chi << endl;

        if(what_type==-1) // -1 is for zeta
        {
            if(abs(chi_1[q_1-1]-1)<.1){
                tmp_l[1]=0;
                tmp_OMEGA=tau_chi/sqrt((Double)((double)q_1));
            }
            else{
                tmp_l[1]=(Double)1/2;
                tmp_OMEGA=tau_chi/(I*sqrt((Double)((double)q_1)));

            }
            tmp_Q=sqrt(double(q_1)/Pi); //XXXXXXXXXXXXXXXXXXXXXX another overload mpreal and long long fix XXXXX check in higher precision
        }


        // If is a degree > 1 L-function
        if(what_type!=-1)
        {


           //for now only cusp forms are implemented
           if(what_type==0||what_type>3||what_type==1||what_type==-2){
               cout << endl << "At present only know how to twist zeta, ";
               cout << "cusp form (in S_k(Gamma_0(N)) L-functions," << endl;
               cout << "and Maass forms for the full modular group" << endl;

               delete [] chi;
               delete [] chi_1;
               delete [] phi;
               delete [] power_p;
               delete [] m;
               delete [] v;
               delete [] sum;
               for(k=0;k<=29;k++) delete(factors[k]);
               delete (factors);

               delete [] tmp_g;
               delete [] tmp_l;
               delete [] tmp_cmplx;

               return 0;

           }


           if(what_type==2||what_type==3){

               tmp_Q=(double(q_1)*underlying_Q);
               //XXXXXXXXXXXXXXXXXXXXXX another overload mpreal and long long fix XXXXX check in higher precision

               tmp_OMEGA=tmp_OMEGA*chi_1[underlying_cond%q_1]*tau_chi*tau_chi/Double((double)q_1);

               //tmp_Q=underlying_Q*r;

               //Double C; // as in the 2*C/Pi

               //C = DIGITS2*log(10.);
               Double T;


               if(!strcmp(what_to_do,"zeros"))
                   T = abs(max(abs(x),abs(y)));
               else //if what_to_do=="values"
                   T = abs(y);

               if(count>0)
               T=T+(count+100.)*Pi/log(T+3);

               Complex delta=int_L.find_delta(1+I*T+tmp_l[1],1);
               tmp_N_terms = Int(2.3 * DIGITS*tmp_Q/real(delta));
               do{
                   tmp_N_terms=(int)(tmp_N_terms*1.3);
               }while(tmp_N_terms*real(delta)/tmp_Q-log((1.+real(tmp_l[1]))*tmp_N_terms)<2.3*DIGITS);
               tmp_N_terms=(int)(tmp_N_terms*1.3);


           }

           if(tmp_N_terms>underlying_N_terms)
               tmp_N_terms=underlying_N_terms;


        }

        //tmp_N_terms=underlying_N_terms;


        if(what_type==-1){
            tmp_N_terms=q_1;
            tmp_what_type=1; //(i.e. is a Dirichlet series with all coeffs)
            tmp_period=q_1;
        }




        if(what_type==-1){ //i.e. if twists of zeta
            if(abs(SUM)>.5) //the trivial character so set to zeta
                Complex_L2= L_function<Complex>();
            else{
                if (my_verbose>1) cout << "calling L_function for twist with " << tmp_N_terms << " coefficients" << endl;
                Complex_L2= L_function<Complex>(tmp_name,tmp_what_type,tmp_N_terms,chi_1,tmp_period,tmp_Q,tmp_OMEGA,tmp_a,tmp_g,tmp_l,0,tmp_cmplx,tmp_cmplx);
            }
        }
        else{
            Complex_twisted_coeffs= new Complex[tmp_N_terms+1];
            for(n=1;n<=tmp_N_terms;n++)
            {
              switch(current_L_type){
               case 1:
                Complex_twisted_coeffs[n]=
                int_L.dirichlet_coefficient[n]*chi_1[n%q_1];
                break;
               case 2:
                Complex_twisted_coeffs[n]=
                Double_L.dirichlet_coefficient[n]*chi_1[n%q_1];
                break;
               case 3:
                Complex_twisted_coeffs[n]=
                Complex_L.dirichlet_coefficient[n]*chi_1[n%q_1];
                break;
              }
            }
            if (my_verbose>1) cout << "calling L_function for twist with " << tmp_N_terms << " coefficients" << endl;
            Complex_L2= L_function<Complex>(tmp_name,tmp_what_type,tmp_N_terms,Complex_twisted_coeffs,tmp_period,tmp_Q,tmp_OMEGA,tmp_a,tmp_g,tmp_l,0,tmp_cmplx,tmp_cmplx);
            delete [] Complex_twisted_coeffs;
        }

        //s=x+I*y;
        //os2 << setprecision(11);
        //os2.seekp(0);
        //os2 << q << " " << Complex_L2.dirichlet_series(s,100000)<< ends;
        //os2 << q << " " <<number_chi<< ends;

        os2.str("");
        //os2.seekp(0);
        os2 << q << " " <<number_chi;

        strcpy(message_stamp,os2.str().c_str());
        //cout << "stamp="<<message_stamp<<endl;

//---------------------find zeros or compute values...
        Double max_imag_chi=0;
        for(n=1;n<q_1;n++){
            if(abs(imag(chi_1[n]))>max_imag_chi)
                 max_imag_chi=abs(imag(chi_1[n]));
        }

        if(!strcmp(what_to_do,"zeros")&&
       ((twist_type<0&&max_imag_chi>.1)||twist_type>=0)
          )
        {

                if(count==0)
                Complex_L2.find_zeros(x,y,step_size,"cout",message_stamp);

                if(current_L_type==3||max_imag_chi>.1)
                Complex_L2.find_zeros(count,start_N,step_size,-1,message_stamp);
                //Complex_L2.find_zeros(count,start_N,step_size,-1,test_explicit_formula,message_stamp);
                else Complex_L2.find_zeros(count,start_N,step_size,-1,message_stamp);
                //else Complex_L2.find_zeros(count,start_N,step_size,-1,test_explicit_formula,message_stamp);
                if(twist_type==-1){M=PHI;m[0]=loop2;m_prime=loop1;}

        }



        if((!strcmp(what_to_do,"values")||!strcmp(what_to_do,"values and ranks"))&&
           ((twist_type<0&&max_imag_chi>.1)||twist_type>=0)
          )
        {
            s=x+I*y;
            u=Complex_L2.value(s,global_derivative,"pure");

            //cout << setprecision(DIGITS3);
            //cout << message_stamp << " " << real(u) << endl;

            cout << message_stamp << " " << real(u) << " " << imag(u);
        if(!strcmp(what_to_do,"values and ranks")){
                if(abs(u)>.000001) cout << " 0";
            else cout << " " << Complex_L2.compute_rank();
        }
        cout << endl;

            //else cout << d << " " << real(u) << endl;

            //else cout << d << " " << real(u) << " " << L_1_chi(d) <<endl;
        if(twist_type==-1){M=PHI;m[0]=loop2;m_prime=loop1;}
        }
//--------------------------------------------------


}
//-----end if primtive, or if not do only primitive (i.e. if do all)


        dial=0;
        if(PHI>1)
        do{
            dial++;
            m[dial]=(m[dial]+1)%phi[dial];
        }while(m[dial]==0&&dial<K);

        //cout<< "--------------------------------------------\n";
    }//for m_prime. Done looping through all characters.
    // } is not aligned since I combined 2 files with the mouse


        //Complex_L2=Complex_L; //XXXXXXX what is this good for 

        delete [] chi;
        delete [] chi_1;
        delete [] phi;
        delete [] power_p;
        delete [] m;
        delete [] v;
        delete [] sum;



    } // for d
    for(k=0;k<=29;k++) delete[](factors[k]);
    delete [] (factors);
    delete [] tmp_g;
    delete [] tmp_l;
    delete [] tmp_cmplx;
                  


    return 0;

}


