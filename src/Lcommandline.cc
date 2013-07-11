/*

Copyright (C) Michael Rubinstein

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


//---------------------------------------------------
//
// Command line interface to the L function package
// By Michael Rubinstein
//
//---------------------------------------------------

#include "Lcommandline.h"
#include "cmdline.h"

void pari_err_recover_nop(long errnum) {return;}

int main (int argc, char *argv[])
{


    Long n;
    bool do_zeros=false,do_values=false,file_to_open=false,do_twists=false,file_of_s_values=false;
    bool file2_to_open=false,do_print=false,do_interpolate=false;
    char data_filename[1000]; //filename of file containing data for L-function.
    char data_filename2[1000]; //filename of file containing data for L-function.
    char xxx_phi_method[100];
    int number_coeff_print=30;
    bool do_hardy=false;
    bool do_rhs_explicit=false; //whether to evaluate the rhs of the explicit formula
    Double alpha; //the alpha in phi(x)=exp(-alpha x^2)

    //used for elliptic curves
    // y^2 + a1 xy + a3 y = x^3 + a2 x^2 + a4 x + a6
    bool do_elliptic_curve=false;
    char a1[200];
    char a2[200];
    char a3[200];
    char a4[200];
    char a6[200];
    int rank=-1;

    bool do_tau=false;
    bool check_rank=false;
    bool do_compute_rank=false;
    bool limit_the_rank=false;
    bool do_only_even_twists=false;


    char s_file_name[1000];  //file of s values

    Double x=0.,y=0.,step_size=10000000;
    Double x2=0.,y2=0.;
    Long count=0;
    Long start_N=0.; //for finding zeros above the Nth zero.

    Double T2=0.; //for computing N(T) at T=T2
    bool do_NT=false;

    char twist_type; //type of twist to do: quadratic, primitive, all, nth
    Long D1=0,D2=0; // for twists

    int print_character=0; //1 is all chi(n) , 2 is just chi(-1)

    bool do_exit=false;
    gengetopt_args_info args_info;
    /* call the cmdline parser */
    if (cmdline_parser (argc, argv, &args_info) != 0){
        exit(1) ;
    }


    if(args_info.zeros_given){ /* Whether zeros was given.  */
        do_zeros=true;
        count=atoll(args_info.zeros_arg);
    }
    if(args_info.N_given){ /* for starting zero search above Nth zero */
        //start_N= str_to_Double(args_info.N_arg);
        start_N= atoll(args_info.N_arg);
    }

    if(args_info.NT_given){ /* for computing N(T) at T=T2 */
        do_NT=true;
        #ifdef USE_LONG_DOUBLE
            sscanf(args_info.NT_arg,"%Lg",&T2);
        #elif USE_MPFR
            T2=args_info.NT_arg;
        #elif USE_MPFRCPP
            T2=args_info.NT_arg;
        #elif USE_BAILEY_DD
            T2=args_info.NT_arg;
        #elif USE_BAILEY_QD
            T2=args_info.NT_arg;
        #else
            sscanf(args_info.NT_arg,"%lg",&T2);
        #endif
    }

    if(args_info.zeros_interval_given){ /* Whether zeros-interval was given.  */
        do_zeros=true;
        if(!args_info.x_given||!args_info.y_given){
            cout << "x and y options must be used" << endl;
            exit(1);
        }
        if(!args_info.stepsize_given){
            cout << "stepsize option must be used" << endl;
            exit(1);
        }
        #ifdef USE_LONG_DOUBLE
            sscanf(args_info.x_arg,"%Lg",&x);
            sscanf(args_info.y_arg,"%Lg",&y);
        #elif USE_MPFR
            x=args_info.x_arg;
            y=args_info.y_arg;
        #elif USE_MPFRCPP
            x=args_info.x_arg;
            y=args_info.y_arg;
        #elif USE_BAILEY_DD
            x=args_info.x_arg;
            y=args_info.y_arg;
        #elif USE_BAILEY_QD
            x=args_info.x_arg;
            y=args_info.y_arg;
        #else
            sscanf(args_info.x_arg,"%lg",&x);
            sscanf(args_info.y_arg,"%lg",&y);
        #endif
        if(y<x){
            cout << "x should be < than y" << endl;
            exit(1);
        }
        #ifdef USE_LONG_DOUBLE
            sscanf(args_info.stepsize_arg,"%Lg",&step_size);
        #elif USE_MPFR
            step_size=args_info.stepsize_arg;
        #elif USE_MPFRCPP
            step_size=args_info.stepsize_arg;
        #elif USE_BAILEY_DD
            step_size=args_info.stepsize_arg;
        #elif USE_BAILEY_QD
            step_size=args_info.stepsize_arg;
        #else
            sscanf(args_info.stepsize_arg,"%lg",&step_size);
        #endif
        //step_size=atof(args_info.stepsize_arg);
    }

    /* Whether to check the explicit formula  */
    if(args_info.rhs_explicit_formula_given){
        do_rhs_explicit=true;
        strcpy(xxx_phi_method,args_info.rhs_explicit_formula_arg);

        if(!args_info.x_given||!args_info.X_given){
            cout << "x and X options must be used" << endl;
            exit(1);
        }
        if(!args_info.stepsize_given){
            cout << "stepsize option must be used" << endl;
            exit(1);
        }
        #ifdef USE_LONG_DOUBLE
            sscanf(args_info.x_arg,"%Lg",&x);
            sscanf(args_info.X_arg,"%Lg",&x2);
            sscanf(args_info.stepsize_arg,"%Lg",&step_size);
            sscanf(args_info.alpha_arg,"%Lg",&alpha);
        #elif USE_MPFR
            x=args_info.x_arg;
            x2=args_info.X_arg;
            step_size=args_info.stepsize_arg;
            alpha=args_info.alpha_arg;
        #elif USE_MPFRCPP
            x=args_info.x_arg;
            x2=args_info.X_arg;
            step_size=args_info.stepsize_arg;
            alpha=args_info.alpha_arg;
        #elif USE_BAILEY_DD
            x=args_info.x_arg;
            x2=args_info.X_arg;
            step_size=args_info.stepsize_arg;
            alpha=args_info.alpha_arg;
        #elif USE_BAILEY_QD
            x=args_info.x_arg;
            x2=args_info.X_arg;
            step_size=args_info.stepsize_arg;
            alpha=args_info.alpha_arg;
        #else
            sscanf(args_info.x_arg,"%lg",&x);
            sscanf(args_info.X_arg,"%lg",&x2);
            sscanf(args_info.stepsize_arg,"%lg",&step_size);
            sscanf(args_info.alpha_arg,"%lg",&alpha);
        #endif
        if(x2<x){
            cout << "x should be < than X" << endl;
            exit(1);
        }
        if(alpha<=0){
            cout << "alpha should be >0" << endl;
            exit(1);
        }
    }

    if(args_info.value_given){ /* Whether value was given.  */

        do_values=true;
        //if(!args_info.x_given||!args_info.y_given){
            //cout << "x and y options must be used" << endl;
            //exit(1);
        //}
        #ifdef USE_LONG_DOUBLE
            sscanf(args_info.x_arg,"%Lg",&x);
            sscanf(args_info.y_arg,"%Lg",&y);
        #elif USE_MPFR
            x=args_info.x_arg;
            y=args_info.y_arg;
        #elif USE_MPFRCPP
            x=args_info.x_arg;
            y=args_info.y_arg;
        #elif USE_BAILEY_DD
            x=args_info.x_arg;
            y=args_info.y_arg;
        #elif USE_BAILEY_QD
            x=args_info.x_arg;
            y=args_info.y_arg;
        #else
            sscanf(args_info.x_arg,"%lg",&x);
            sscanf(args_info.y_arg,"%lg",&y);
        #endif
    }
    if(args_info.value_file_given){ /* Whether value-file was given.  */
        do_values=true;
        file_of_s_values=true;
        strcpy(s_file_name,args_info.value_file_arg);
    }
    if(args_info.value_line_segment_given){/* Whether value-line-segment was given.  */
        do_values=true;
        if(!args_info.x_given||!args_info.y_given||
        !args_info.X_given||!args_info.Y_given){
            cout << "x,y,X,Y options must be used" << endl;
            do_exit=true;
        }
        if(!args_info.number_samples_given){
            cout << "--number-samples option must be used" << endl;
            do_exit=true;
        }
        if(do_exit) exit(1);
        #ifdef USE_LONG_DOUBLE
            sscanf(args_info.x_arg,"%Lg",&x);
            sscanf(args_info.y_arg,"%Lg",&y);
            sscanf(args_info.X_arg,"%Lg",&x2);
            sscanf(args_info.Y_arg,"%Lg",&y2);
        #elif USE_MPFR
            x=args_info.x_arg;
            y=args_info.y_arg;
            x2=args_info.X_arg;
            y2=args_info.Y_arg;
        #elif USE_MPFRCPP
            x=args_info.x_arg;
            y=args_info.y_arg;
            x2=args_info.X_arg;
            y2=args_info.Y_arg;
        #elif USE_BAILEY_DD
            x=args_info.x_arg;
            y=args_info.y_arg;
            x2=args_info.X_arg;
            y2=args_info.Y_arg;
        #elif USE_BAILEY_QD
            x=args_info.x_arg;
            y=args_info.y_arg;
            x2=args_info.X_arg;
            y2=args_info.Y_arg;
        #else
            sscanf(args_info.x_arg,"%lg",&x);
            sscanf(args_info.y_arg,"%lg",&y);
            sscanf(args_info.X_arg,"%lg",&x2);
            sscanf(args_info.Y_arg,"%lg",&y2);
        #endif
    }
    if(args_info.number_samples_given){ /* Whether derivative was given.  */
        count=args_info.number_samples_arg;
    }
    if(args_info.hardy_given){ /* Whether Hardy Z option was given.  */
        do_hardy=true;
    }

    if(args_info.elliptic_curve_given){ /* Whether elliptic-curve was given.  */
        do_elliptic_curve=true;
        if(!args_info.a1_given){ cout << "must specify a1 with --a1" << endl; do_exit=true;}
        if(!args_info.a2_given){ cout << "must specify a2 with --a2" << endl; do_exit=true;}
        if(!args_info.a3_given){ cout << "must specify a3 with --a3" << endl; do_exit=true;}
        if(!args_info.a4_given){ cout << "must specify a4 with --a4" << endl; do_exit=true;}
        if(!args_info.a6_given){ cout << "must specify a6 with --a6" << endl; do_exit=true;}
        if(do_exit) exit(1);
        strcpy(a1,args_info.a1_arg);
        strcpy(a2,args_info.a2_arg);
        strcpy(a3,args_info.a3_arg);
        strcpy(a4,args_info.a4_arg);
        strcpy(a6,args_info.a6_arg);
    }
    if(args_info.file_input_given){ /* Whether file-input was given.  */
        file_to_open=true;
        strcpy(data_filename,args_info.file_input_arg);
    }
    if(args_info.url_given){ /* Whether file-input was given.  */
        char str1[]="wget -O temporary_url_file_lcalc "; 
        if(system(strcat(str1,args_info.url_arg))!=0){
             cout << "Error retrieving file: " << args_info.url_arg << endl;
             do_exit=true;
        }
        if(do_exit) exit(1);
        file_to_open=true;
        strcpy(data_filename,"temporary_url_file_lcalc");
    }
    if(args_info.interpolate_given){ /* Whether interpolate was given.  */
        do_interpolate=true;
        file2_to_open=true;
        strcpy(data_filename2,args_info.interpolate_arg);
    }
    if(args_info.output_character_given){ /* Whether output-character was given.  */
        print_character = atoi(args_info.output_character_arg);
    }
    if(args_info.output_data_given){ /* Whether output-data was given.  */
        do_print=true;
        number_coeff_print=args_info.output_data_arg;
    }
    if(args_info.precision_given){ /* Whether precision was given.  */
         DIGITS=args_info.precision_arg;
         DIGITS3=DIGITS-DIGITS2;

    }
    if(args_info.sacrifice_given){ /* Whether sacrifice was given.  */
         DIGITS2=args_info.sacrifice_arg;
         DIGITS3=DIGITS-DIGITS2;

    }
    if(args_info.rank_compute_given){ /* Whether rank-compute was given.  */
        do_compute_rank=true;
    }
    if(args_info.rank_verify_given){ /* Whether rank-verify was given.  */
        check_rank=true;
        rank=args_info.rank_verify_arg;
    }
    if(args_info.rank_limit_given){ /* Whether rank-limit was given.  */
        limit_the_rank=true;
        rank=args_info.rank_limit_arg;
    }
    if(args_info.tau_given){ /* Whether tau was given.  */
        do_tau=true;
    }
    if(args_info.twist_quadratic_given){ /* Whether twist-quadratic was given.  */
        do_twists=true;
        twist_type='q';
    }
    if(args_info.twist_quadratic_even_given){ /* Whether twist-quadratic was given.  */
        do_twists=true;
        do_only_even_twists=true;
        twist_type='q';
    }
    if(args_info.twist_primitive_given){ /* Whether twist-primitive was given.  */
        do_twists=true;
        twist_type='p';
    }
    if(args_info.twist_all_given){ /* Whether twist-all was given.  */
        do_twists=true;
        twist_type='A';
    }
    if(args_info.twist_all_no_conj_pairs_given){ /* Whether twist-all-no-conj-pairs was given.  */
        do_twists=true;
        twist_type='a';
    }
    if(args_info.twist_complex_no_conj_pairs_given){ /* Whether twist-complex-no-conj-pairs was given.  */
        do_twists=true;
        twist_type='c';
    }
    if(args_info.twist_generic_given){ /* Whether twist-complex-no-conj-pairs was given.  */
        do_twists=true;
        twist_type='g';
    }
    if(args_info.degree_given){ /* Whether degree was given.  */
        do_twists=true;
        twist_type='n';
        n=args_info.degree_arg;
    }
    if(do_twists){
        if(!args_info.start_given){ 
            cout << "must specify starting discriminant/conductor (--finish or -f)" << endl; 
            do_exit=true;
        }
        if(!args_info.finish_given){ 
            cout << "must specify finishing discriminant/conductor (with --start or -s)" << endl;
            do_exit=true;
        }
        if(do_exit) exit(1);
        D1=strtoll(args_info.start_arg,0,10); D2=strtoll(args_info.finish_arg,0,10);
    }

    #ifdef _OPENMP
        if(args_info.openmp_given) /* whether a number of threads was specified */
        {
            /* get the total number of CPUs/cores available for OpenMP */
            int    NCPU = omp_get_num_procs();
            if (args_info.openmp_arg>NCPU){
                cout << "ERROR: You only have " << NCPU << " processors available." << endl;
                cout << "Specify fewer processors" << endl;
                exit(1);
            }
            omp_set_num_threads(args_info.openmp_arg);
        }
        else //default is 1
            omp_set_num_threads(1);
    #endif
    #ifndef _OPENMP
        if(args_info.openmp_given){ /* whether a number of threads was specified */
            cout << endl;
            cout << "ERROR: lcalc has not been compiled with openmp enabled." << endl;
            cout << "You need to uncomment the line: " << endl;
            cout << "#OPENMP_FLAG = -fopenmp" << endl;
            cout << "in the Makefile, and recompile by typing:"<< endl;
            cout << endl;
            cout << "make clean" << endl;
            cout << "make" << endl;
            cout << endl;
            exit(1);
        }
    #endif

    if(args_info.verbose_given){ /* Whether vebosity level was specified.  */
        my_verbose=args_info.verbose_arg;
        if(my_verbose==-33) do_blfi=true;
    }

    //placed here rather than at top since
    //precision depends on user input.
    initialize_globals();
    cout << setprecision(DIGITS3);


    //the L-functions used throughout are global and
    //need to be re-initialized now that we have initialized
    //Pi in initialize_globals.
    initialize_commandline_globals();
    if(args_info.derivative_given){ /* Whether derivative was given.  */
        global_derivative = args_info.derivative_arg;
        //DIGITS3=(int)(DIGITS3/pow(2.,global_derivative)); // now taken care of in Lvalue.h
    }

    if(do_zeros){
        input_mean_spacing_given = (6.28/log((count+start_N)+100.))*2./DIGITS; //only used
        //for the zeta function, and then in the band limited interpolation scheme. Is a
        //rough estimate for the average distance between sample points.
        //Dividing by DIGITS/2 more than accounts for searching for and then zooming in
        //on zeros.
        //cout << " input_mean_spacing_given " << input_mean_spacing_given << endl;
        //if(input_mean_spacing_given<.3) do_blfi=true;
    }

    if(args_info.use_blfi_given){ /* Whether to compute Dirichlet series using blf interpolation */
        try_use_blfi=true;
        blfi_interval_length=args_info.use_blfi_arg;
        if (blfi_interval_length<=0) blfi_interval_length = 1;
    }

    if(args_info.use_dirichlet_series_given){ /* Whether to compute just using the Dirichlet series */
        only_use_dirichlet_series=true;
        if(!args_info.number_terms_given){
            cout << "--number-terms option must be used" << endl;
            exit(1);
        }
        //N_use_dirichlet_series=args_info.use_dirichlet_series_arg; //how many terms to use
        N_use_dirichlet_series=args_info.number_terms_arg; //how many terms to use
    }




    A=1./(16*Pi*Pi*2.3*DIGITS); //controls the 'support' of g(w)

    if(file2_to_open){
        initialize_new_L(data_filename2);
        switch(current_L_type)
        {
            case 1:
                int_L2=int_L;
                break;
            case 2:
                Double_L2=Double_L;
                break;
            case 3:
                Complex_L2=Complex_L;
                break;
        }

    }

    if(file_to_open){
        initialize_new_L(data_filename);
        if(args_info.url_given){ /* Whether file-input was given.  */
            system("rm temporary_url_file_lcalc"); 
        }
    }
    else current_L_type=1; //default is zeta, and type 1, i.e. int, is simplest

    if(do_interpolate) //compute the zeros on the critical line
    //of the L-series that interpolate between f and f2. Must be used in combination with 
    //with the option: -zi x y increment 
    //By interpolate I mean take t*(basic data of L) and (1-t)*(basic data of L2).
    {
         L_interpolate(x,y,step_size);
         return 0;
    }


    Double *coeff;

    if(do_elliptic_curve||do_tau){ //determine how many dirichlet coefficients are needed
                                   //and initialize the curve


        int N_terms; //number of dirichlet coefficients
        Double T;
        Double tmp_Q; //the maximum normalized conductor (including twists)

        Double t2;


        if(do_elliptic_curve){
#ifdef INCLUDE_PARI

            t2=.5; //t2=.5 because of the GAMMA(s+1/2)

            pari_init_opts(400000000,2,INIT_DFTm); // the last option is to prevent
            //pari from giving its interrupt signal when its elliptic curve a_p
            //algorithm is called and interrupted with ctrl-c.

            coeff = new Double[3];
            //compute the conductor which is copied to coeff[1]
            data_E(a1,a2,a3,a4,a6,2,coeff);

            tmp_Q=sqrt(coeff[1])/(2*Pi);
            delete [] coeff;
#else
            cout << "You need to uncomment the line: PARI_DEFINE = -DINCLUDE_PARI" <<endl;
            cout << "in the Makefile and do: 'make clean', then 'make' if you wish to use" <<endl;
            cout << "elliptic curve L-functions. Requires that you already have pari installed" <<endl;
            cout << "on your machine." <<endl;
            exit(1);
#endif //ifdef INCLUDE_PARI
        }


        if(do_tau){
            t2=5.5;     //t2=5.5 because of the GAMMA(s+11/2)
            tmp_Q=1./(2*Pi);
        }



        if(do_twists) tmp_Q=tmp_Q*max(abs(1.*D1),abs(1.*D2)); // take into account twists


        if(do_zeros)
            T = max(abs(x),abs(y));
        else
            T = max(abs(y),abs(y2));


        if(count>0&&!args_info.number_samples_given)
            T=T+(start_N+count+100.)*Pi/log(T+3);

        if(my_verbose>0) cout << "#    T = " << T << endl;

        //based on the number of terms used in the gamma_sum routine. possible updates to gamma_sum condition should 
        //thus be reflected here
        Complex delta=int_L.find_delta(1+I*T+t2,1);
        if(args_info.number_terms_given) N_terms=args_info.number_terms_arg;
        else{
            N_terms = Int(2.3 * DIGITS*tmp_Q/real(delta));
            do{
                N_terms=(int)(N_terms*1.3);
                if(my_verbose>3) cout << "N_terms to precompute = " << N_terms << endl;
            }while(N_terms*real(delta)/tmp_Q-log((1.+t2)*N_terms)<2.3*DIGITS);
            N_terms=(int)(N_terms*1.3+40);
        }
        if(my_verbose>0) cout << "#    N_terms to precompute = " << N_terms << endl;



#ifdef INCLUDE_PARI
        if(do_elliptic_curve){
            void (*saved_err_recover)(long) = cb_pari_err_recover;
            cb_pari_err_recover = pari_err_recover_nop;
            allocatemem(N_terms*16 + 1000000); //XXXXXXXXX this should depend on whether we're double or long double or mpfr double
            cb_pari_err_recover = saved_err_recover;


             if (my_verbose>0) cout << "#    Will precompute " << N_terms << " elliptic L-function dirichlet coefficients..." << endl;
             current_L_type=2;

             Double_L=L_function<Double>(a1,a2,a3,a4,a6,N_terms);
        }
#endif //ifdef INCLUDE_PARI

        if(do_tau){

            tmp_Q=1./(2*Pi);

            coeff =  new Double[N_terms+1];
            ramanujan_tau(coeff,N_terms);

            Double *g;
            Complex *l;
            Complex *p;
            Complex *r;

            g=new Double[2];
            l=new Complex[2];
            g[1]=1.;
            l[1]=5.5;

            p = new Complex[1];
            r = new Complex[1];



            if (my_verbose>0) cout << "#    Will precompute " << N_terms << " Ramanujan tau(n) dirichlet coefficients..." << endl;
            current_L_type=2; //the normalized dirichlet coeffs are real
            Double_L=L_function<Double>("Ramanujan Tau",2,N_terms,coeff,0,tmp_Q,1,1,g,l,0,p,r);

            delete [] g;
            delete [] l;
            delete [] p;
            delete [] r;
            delete [] coeff;

        }



    }


    if(do_print) print_L(number_coeff_print);

    if(check_rank&&rank>=0) verify_rank(rank);

    if(do_compute_rank&&!do_twists) compute_rank();

    if(do_NT){
        switch(current_L_type)
        {
            case 1:
                cout << "N("<< T2 << ") = " << int_L.N(T2) << endl;
                break;
            case 2:
                cout << "N("<< T2 << ") = " << Double_L.N(T2) << endl;
                break;
            case 3:
                cout << "N("<< T2 << ") = " << Complex_L.N(T2) << endl;
                break;
        }

    }

    if(do_values){
        if(do_twists){
            switch(twist_type){
                case 'q':
                    if(do_compute_rank)
                        quadratic_twists(D1,D2,x,y,0,0,0,"values and ranks",do_only_even_twists);
                        //quadratic_twists(D1,D2,x,y,0,0,0,"values and ranks",do_only_even_twists,do_explicit);
                    else if(limit_the_rank)
                        //quadratic_twists(D1,D2,x,y,0,0,0,"values and ranks",do_only_even_twists,rank,do_explicit);
                        quadratic_twists(D1,D2,x,y,0,0,0,"values and ranks",do_only_even_twists,rank);
                    else
                        quadratic_twists(D1,D2,x,y,0,0,0,"values",do_only_even_twists);
                        //quadratic_twists(D1,D2,x,y,0,0,0,"values",do_only_even_twists,do_explicit);
                    break;
                case 'p':
                    if(do_compute_rank)
                        all_twists(D1,D2,x,y,0,0,0,"values and ranks",0,print_character);
                        //all_twists(D1,D2,x,y,0,0,0,"values and ranks",0,print_character,do_explicit);
                    else
                        all_twists(D1,D2,x,y,0,0,0,"values",0,print_character);
                    break;
                case 'a':
                    if(do_compute_rank)
                        all_twists(D1,D2,x,y,0,0,0,"values and ranks",1,print_character);
                    else
                        all_twists(D1,D2,x,y,0,0,0,"values",1,print_character);
                    break;
                case 'A':
                    if(do_compute_rank)
                        all_twists(D1,D2,x,y,0,0,0,"values and ranks",2,print_character);
                    else
                        all_twists(D1,D2,x,y,0,0,0,"values",2,print_character);
                    break;
                case 'g':
                    if(do_compute_rank)
                        all_twists(D1,D2,x,y,0,0,0,"values and ranks",-1,print_character);
                    else
                        all_twists(D1,D2,x,y,0,0,0,"values",-1,print_character);
                        //all_twists(D1,D2,x,y,0,0,0,"values",-1,print_character,do_explicit);
                    break;
                case 'c':
                    if(do_compute_rank)
                        all_twists(D1,D2,x,y,0,0,0,"values and ranks",-2,print_character);
                        //all_twists(D1,D2,x,y,0,0,0,"values and ranks",-2,print_character,do_explicit);
                    else
                        all_twists(D1,D2,x,y,0,0,0,"values",-2,print_character);
                        //all_twists(D1,D2,x,y,0,0,0,"values",-2,print_character,do_explicit);
                    break;
            }
        }
        else if(file_of_s_values){
             if(!do_hardy) compute_values(x,y,"pure",s_file_name);
             else compute_values(x,y,"rotated pure",s_file_name);
        }
        else{

            input_mean_spacing_given=(y2-y)/double(count); //used in Riemann Siegel band limited routine
            //if(input_mean_spacing_given<.3) do_blfi=true;

            //cout << " input_mean_spacing_given " << input_mean_spacing_given << endl;

            if(!do_hardy) compute_values(x,y,"pure","",x2,y2,count);
            else compute_values(x,y,"rotated pure","",x2,y2,count);
        }

    }
    else if(do_zeros){
        if(do_twists){
            switch(twist_type){
                case 'q':
                    if(limit_the_rank)
                    quadratic_twists(D1,D2,x,y,count,start_N,step_size,"zeros and ranks",do_only_even_twists,rank);
                    //quadratic_twists(D1,D2,x,y,count,start_N,step_size,"zeros and ranks",do_only_even_twists,do_explicit,rank);
                    else
                    quadratic_twists(D1,D2,x,y,count,start_N,step_size,"zeros",do_only_even_twists);
                    break;
                case 'p':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"zeros",0,print_character);
                    break;
                case 'a':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"zeros",1,print_character);
                    break;
                case 'A':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"zeros",2,print_character);
                    break;
                case 'g':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"zeros",-1,print_character);
                    break;
                case 'c':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"zeros",-2,print_character);
                    break;

            }
        }
        else{
            compute_zeros(x,y,step_size,count,start_N,rank);
        }
    }
    else if(do_rhs_explicit){

        int N_terms=-1;
        if(args_info.number_terms_given) N_terms=args_info.number_terms_arg;
        switch(current_L_type)
        {
            case 1:
                int_L.plot_explicit_formula(alpha, x, x2, step_size, xxx_phi_method,N_terms);
                break;
            case 2:
                Double_L.plot_explicit_formula(alpha, x, x2, step_size, xxx_phi_method,N_terms);
                break;
            case 3:
                Complex_L.plot_explicit_formula(alpha, x, x2, step_size, xxx_phi_method,N_terms);
                break;
        }
    }
    else{ //else allow for printing of characters. "print character" is a dummy char str,
          //so all_twists will print out the character without doing zeros or values
        if(do_twists){
            switch(twist_type){
                case 'p':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"print character",0,print_character);
                    break;
                case 'a':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"print character",1,print_character);
                    break;
                case 'A':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"print character",2,print_character);
                    break;
                case 'g':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"print character",-1,print_character);
                    break;
                case 'c':
                    all_twists(D1,D2,x,y,count,start_N,step_size,"print character",-2,print_character);
                    break;

            }
        }

    }


    delete_globals();
    cmdline_parser_free (&args_info); /* release allocated memory */

    return 0;
}
