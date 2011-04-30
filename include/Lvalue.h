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

    template <class ttype>
    Complex L_function <ttype>::
    find_delta (Complex z,Double g)
    {

        //cout << " find delta z g " << z << " " << g<< endl;
        Double sigma=real(z);
        Double t=imag(z); if(t<0) t=-t;
        Double r=abs(z);
        Double theta=atan(t/sigma);
        Double epsilon;

        Double a=-theta;
        Double b=0.;
        Double c,f1,f3;
        Double local_tolerance=.01/(t+100); if (local_tolerance<tolerance) local_tolerance=tolerance;

        f1=sigma*log(sigma/(r*lcalc_cos(theta+a))) - t*a;

        if(f1<=DIGITS2*2.3) epsilon=-theta;
        else{
            do{
                c=(a+b)/2;
                f3=sigma*log(sigma/(r*lcalc_cos(theta+c))) - t*c;
                if(f3>DIGITS2*2.3)a=c;
                else b=c;
                //cout<< "theta+epsilon: " << a+theta<< endl;
            } while(b-a>local_tolerance);
            epsilon=a;
        }
        //if(imag(z)>=0) cout << " returning delta: " << exp(I*(theta+epsilon)*g) << endl;
        //else cout << " returning delta: " << exp(-I*(theta+epsilon)*g) << endl;

        if(imag(z)>=0) return exp(I*(theta+epsilon)*g);
        else return exp(-I*(theta+epsilon)*g);
   }


    //computes (3.2.5) as a Riemann sum using
    //g(w) = exp(A*(w-s)^2) * delta^(-w)

    template <class ttype>
    Complex L_function <ttype>::
    value_via_Riemann_sum(Complex s, const char *return_type)
    {
        int j,k,m,mm,n;
        Complex r,z;
        Complex SUM=0.;

        Double tmp;

        Complex L_value=0.;
        Double theta;

        Complex *DELTA; // variant on (3.3.10), without the theta

        Double t_0;
        Double c;    //controls speed of convergence but at the expense
                     //of loss of precision.

        //Double v=.5; // the v in (3.2.5)
        Double v=1-real(s); // the v in (3.2.5)
        if(v<.5)v=.5;
        //Double incr; // incrememnt size in the Riemann sum- is now a global variable

        Complex dirichletseries;
        Complex * __restrict__ dirichlet_vector; //used to compute Dirichlet series incrementally
        Complex * __restrict__ dirichlet_vector_copy; //used to compute Dirichlet series incrementally
        Complex * __restrict__ dirichlet_multiplier; // stores powers of n^{-I incr}
        Complex * __restrict__ dirichlet_multiplier_conj; // stores powers of n^{I incr}. Prefer to store since
                                                      // this will make vectorization easier.

        //int M; //number of terms to take in the Riemann sum
        //M is no longer used. escape is determined numerically. Is fine to
        //escape this way. Gamma and exp behave predictably

        int N; //the number of terms to take in the Dirichlet series
        Double r1=0,r2=0,r3=0,r4=0,r5=0,mynorm_r,local_average, max_integrand; //used to decide when to truncate the Riemann sum


        //cout << "a= " << a << endl;
        theta=0.;
        for(j=1;j<=a;j++)
        {
            theta = theta+gamma[j];
            //cout << "theta = " << theta << endl;
        }

        c=DIGITS2*log(10.); //i.e exp(-c)=10^(-DIGITS2)
                            //this sacrifices roughly at most DIGITS2
                            //out of DIGITS precision.


        //incr is now a global variable
        //incr=2*Pi*v/(log(10.)*DIGITS);

        //M=Int(sqrt(log(10.)*DIGITS/A)/incr);



        DELTA = new Complex[a+1];

        //if(abs(t_0)<=c/(2*sqrt(A))) tmp= 2*sqrt(A);
        //else tmp=c/t_0;


        Double c1=0.;
        for(j=1;j<=a;j++){

            t_0=imag(gamma[j]*s+lambda[j]);

        //cout << "tmp_" << j << " =   " << tmp  << endl;
        //cout << "t_0" << " =   " << t_0  << endl;

            if(abs(t_0)<=2*c/(Pi*a)) tmp= Pi/2;
            else tmp=abs(c/(t_0*a));

            if(t_0>=0)r=1; else r=-1;

            DELTA[j]= exp(I*r*(Pi/2-tweak*tmp)); //tweak defaults to 1. It allows me to globally set a slightly different angle
                                                 //for the purpose of testing precision or looking for L-functions

            c1=c1+gamma[j]*tmp;

            //DELTA[j]=find_delta(s*gamma[j]+lambda[j],gamma[j]);


        }


        for(k=1;k<=number_of_poles;k++){
           z =A*(pole[k]-s)*(pole[k]-s);
           for(j=1;j<=a;j++) z=z-log(DELTA[j])*(gamma[j]*pole[k]+lambda[j]);
           //the 5 below is for kicks. 2.3 would have been fine.
           if(real(z)>-5*DIGITS)
               L_value=L_value+residue[k]*exp(z)/(s-pole[k]);
        }
        //cout << "poles contribute: " << L_value << endl;

        //the rough estimate: G(z,(N DELTA/Q)^2) is, in size,
        //roughly exp(-Re((N*DELTA/Q)^(1/theta))) and we want this
        //to be > 2.3 DIGITS XXXXXXXXX check this

        double N_as_double = lcalc_to_double(Q*exp(log(2.3*DIGITS*theta/c1)*theta)+10); //as double because we might exceed INT_MAX

        if(N_as_double>number_of_dirichlet_coefficients&&what_type_L!=-1&&what_type_L!=1)
        {

            if(print_warning){
                print_warning=false;
                cout << "WARNING from Riemann sum- we don't have enough Dirichlet coefficients." << endl;
                cout << "Will use the maximum possible, though the output ";
                cout << "will not necessarily be accurate." << endl;
            }
            N=number_of_dirichlet_coefficients;
        }
        else N=Int(Q*exp(log(2.3*DIGITS*theta/c1)*theta)+10);

        if(N>number_logs) extend_LG_table(N);


        dirichlet_vector= new Complex[N+1]; //initially stores a(n)/n^{s+v} (or 1-s instead of s)
        dirichlet_vector_copy= new Complex[N+1]; // used to do negative m
        dirichlet_multiplier= new Complex[N+1];
        dirichlet_multiplier_conj= new Complex[N+1];
        #pragma omp parallel for //shared(N,I,incr,dirichlet_multiplier) private(n)
        for(n=1;n<=N;n++){
            dirichlet_multiplier[n]=exp(-I*LOG(n)*incr); //mutiplying by this repeatedly will accumulate error in one direction, so might want to improve
            dirichlet_multiplier_conj[n]=conj(dirichlet_multiplier[n]);

            if(what_type_L==-1)   //i.e. if the Riemann zeta function
                dirichlet_vector[n]=exp(-(s+v)*LOG(n));
            else if(what_type_L!=1) //if not periodic
                dirichlet_vector[n]=dirichlet_coefficient[n]*exp(-(s+v)*LOG(n));
            else //if periodic
            {
                m=n%period; if(m==0) m=period;
                dirichlet_vector[n]=dirichlet_coefficient[m]*exp(-(s+v)*LOG(n));
            }
            dirichlet_vector_copy[n]=dirichlet_vector[n];
        }


        max_n=N;

        if(my_verbose>1){
            cout << "#        s =  " << s << "  Will use  " << N << " terms of the Dirichlet series" << endl;
            cout << "#        entering Riemann sum. incr =  " << incr << " A = " << A<<endl;
        }

        Double log_Q=log(Q);

/* old riemann sum. escape was fixed ahead of time and was not efficient.
        for(m=-M;m<=M;m++){
            r=exp(A*(v+I*incr*m)*(v+I*incr*m)+log_Q*(s+v+I*incr*m));
            for(j=1;j<=a;j++){
                //cout << "gamma[j]*(s+v+I*incr*m)+lambda[j]= " << gamma[j]*(s+v+I*incr*m)+lambda[j] << endl;
                //cout << "DELTA[j] =   " << DELTA[j]  << endl;
                r=r*GAMMA(gamma[j]*(s+v+I*incr*m)+lambda[j],DELTA[j]);
            }
            //cout << "r= " << r << endl;

            //r=r*this->dirichlet_series(s+v+I*incr*m,N,false);

            r=r*this->dirichlet_series(s+v+I*incr*m,N);
            //r=r*this->dirichlet_series(s+v+I*incr*m,N);

            //dirichlet series part needs to be more precise
            SUM=SUM+r/(v+I*incr*m);
            //if(m%100==0)
            cout << "m= " << m << "  SUM1 = " << SUM << endl;
        }
        SUM=SUM*incr/(2*Pi);
        //cout << "m= " << m << "  SUM1 = " << SUM << endl;
*/

        max_integrand=0.; //the max of the integrand, without the dirichlet series factor
        mm=0;
        //first do the terms m >=0
        do{
            for(m=mm;m<=mm+99;m++){
                Complex zz=v+I*incr*m;;
                //r=exp(A*(v+I*incr*m)*(v+I*incr*m)+log_Q*(s+v+I*incr*m));
                r=exp(A*zz*zz+log_Q*(s+zz));
                for(j=1;j<=a;j++){
                    //r=r*GAMMA(gamma[j]*(s+v+I*incr*m)+lambda[j],DELTA[j]);
                    r*=GAMMA(gamma[j]*(s+zz)+lambda[j],DELTA[j]);
                }
                r/=zz;

                //treat these as doubles, to save on multiprecision arithmetic
                mynorm_r=my_norm(r);
                if(mynorm_r>max_integrand) max_integrand = mynorm_r;

                r1=r2;r2=r3;r3=r4;r4=r5;r5=mynorm_r;
                local_average=(r1+r2+r3+r4+r5)/5;

                //XXXXXX replaced by vectorized version
                //r=r*this->dirichlet_series(s+v+I*incr*m,N);

                if(try_use_blfi&&N>999){
                    //cout << "# go: dirichlet_series_via_blfi" << endl;
                    dirichletseries= this->dirichlet_series_via_blfi(s+zz, N,blfi_interval_length);
                }
                else{
                    dirichletseries=0.;
                    for(j=1;j<=N;j++){ //XXXXXXXXXX might be good candidate for unrolling and also decrementing rather than incrementing
                        dirichletseries+=dirichlet_vector[j];
                        dirichlet_vector[j]*=dirichlet_multiplier[j];
                    }
                }

                r*=dirichletseries;
                //cout << "1 dirichletseries: " << dirichletseries << endl;
                //cout << "1 thischletseries: " << this->dirichlet_series(s+v+I*incr*m,N) << endl;


                SUM=SUM+r;
                if(my_verbose>2) cout << "#            m= " << m << "  SUM = " << SUM << "  r= " << r  << "  local average = " << local_average << " max=  "<< max_integrand<< endl;

            }
            mm=m;
        }while(local_average>max_integrand*tolerance_sqrd);

        m=-1;
        //then do the terms negative m
        do{
            Complex zz=v+I*incr*m;;
            r=exp(A*zz*zz+log_Q*(s+zz));
            for(j=1;j<=a;j++){
                r*=GAMMA(gamma[j]*(s+zz)+lambda[j],DELTA[j]);
            }
            r/=zz;

            mynorm_r=my_norm(r);
            if(mynorm_r>max_integrand) max_integrand = mynorm_r;
            r1=r2;r2=r3;r3=r4;r4=r5;r5=mynorm_r;
            local_average=(r1+r2+r3+r4+r5)/5;

            //r=r*this->dirichlet_series(s+v+I*incr*m,N);

            if(try_use_blfi&&N>999){
                dirichletseries= this->dirichlet_series_via_blfi(s+zz, N,blfi_interval_length);
            }
            else{
                dirichletseries=0.;
                for(j=1;j<=N;j++){ //XXXXXXXXXX might be good candidate for unrolling and also decrementing rather than incrementing
                    dirichlet_vector_copy[j]*=dirichlet_multiplier_conj[j];
                    dirichletseries+=dirichlet_vector_copy[j];
                }
            }

            //dirichletseries=0.;
            //for(j=1;j<=N;j++){
                //dirichlet_vector_copy[j]*=dirichlet_multiplier_conj[j];
                //dirichletseries+=dirichlet_vector_copy[j];
            //}

            r=r*dirichletseries;
            //cout << "2 dirichletseries: " << dirichletseries << endl;
            //cout << "2 thischletseries: " << this->dirichlet_series(s+v+I*incr*m,N) << endl;


            SUM=SUM+r;
            if(my_verbose>2) cout << "#            m= " << m << "  SUM = " << SUM << "  r= " << r  << "  local average = " << local_average << " max=  "<< max_integrand<< endl;
            m--;
        }while(m>-100||local_average>max_integrand*tolerance_sqrd);

        SUM*=incr/(2*Pi);


        //no longer needed
        //r=0;
        //for(j=1;j<=a;j++)r=r+lambda[j];
        //SUM=SUM*exp(log(DELTA)*r);

        //cout << "m= " << m << "  SUM1 = " << SUM << endl;

        L_value=L_value+SUM;


        if(real(s)!=.5){ //do the second sum i.e. for f_2

            v=real(s);

            if(what_type_L==-1)   //i.e. if the Riemann zeta function
            {
                #pragma omp parallel for shared(N,dirichlet_vector,s,v) private(n)
                for(n=1;n<=N;n++) dirichlet_vector[n]=exp(-conj(1-s+v)*LOG(n));
            }
            else if(what_type_L!=1) //if not periodic
            {
                #pragma omp parallel for shared(N,dirichlet_vector,s,v) private(n)
                for(n=1;n<=N;n++) dirichlet_vector[n]=dirichlet_coefficient[n]*exp(-conj(1-s+v)*LOG(n));
            }
            else //if periodic
            {
                for(n=1;n<=N;n++)
                {
                    m=n%period; if(m==0)m=period;
                    dirichlet_vector[n]=dirichlet_coefficient[m]*exp(-conj(1-s+v)*LOG(n));
                }
            }

            for(n=1;n<=N;n++) dirichlet_vector_copy[n]=dirichlet_vector[n];

            SUM=0.;


/*
            for(m=-M;m<=M;m++){
                r=exp(A*(v+I*incr*m)*(v+I*incr*m)+log_Q*(1-s+v+I*incr*m));
                for(j=1;j<=a;j++)
                        r=r*GAMMA(gamma[j]*(1-s+v+I*incr*m)+conj(lambda[j]),1/DELTA[j]);
                //r=r*conj(this->dirichlet_series(conj(1-s+v+I*incr*m),N,false));

                r=r*conj(this->dirichlet_series(conj(1-s+v+I*incr*m),N));

                //dirichlet series part needs to be more precise
                SUM=SUM+r/(v+I*incr*m);
                //if(m%100==0)
                //cout << "m= " << m << "  SUM2 = " << SUM << endl;
            }
            SUM=SUM*incr/(2*Pi);
*/

            max_integrand=0.; //the max of the integrand, without the dirichlet series factor
            m=0;
            //first do the terms m >=0
            do{
                Complex zz=v+I*incr*m;
                r=exp(A*zz*zz+log_Q*(1-s+zz));
                for(j=1;j<=a;j++){
                    r*=GAMMA(gamma[j]*(1-s+zz)+conj(lambda[j]),1/DELTA[j]);
                }
                r/=zz;


                mynorm_r=my_norm(r);
                if(mynorm_r>max_integrand) max_integrand = mynorm_r;
                r1=r2;r2=r3;r3=r4;r4=r5;r5=mynorm_r;
                local_average=(r1+r2+r3+r4+r5)/5;

                //r=r*conj(this->dirichlet_series(conj(1-s+v+I*incr*m),N));
                if(try_use_blfi&&N>999){
                    dirichletseries= this->dirichlet_series_via_blfi(1-s+zz, N,blfi_interval_length);
                }
                else{
                    dirichletseries=0.;
                    for(j=1;j<=N;j++){ //XXXXXXXXXX might be good candidate for unrolling and also decrementing rather than incrementing
                        dirichletseries+=dirichlet_vector[j];
                        dirichlet_vector[j]*=dirichlet_multiplier_conj[j];
                    }
                }
                //dirichletseries=0.;
                //for(j=1;j<=N;j++){
                    //dirichletseries+=dirichlet_vector[j];
                    //dirichlet_vector[j]*=dirichlet_multiplier_conj[j];
                //}
                r*=conj(dirichletseries);
                //cout << "3 dirichletseries: " << dirichletseries << endl;
                //cout << "3 thischletseries: " << this->dirichlet_series(conj(1-s+v+I*incr*m),N) << endl;

                SUM+=r;
                if(my_verbose>2) cout << "#            m= " << m << "  SUM = " << SUM << "  r= " << r  << "  local average = " << local_average << " max=  "<< max_integrand<< endl;
                m++;
            }while(m<100||local_average>max_integrand*tolerance_sqrd);

            m=-1;
            //then do the terms negative m
            do{
                Complex zz=v+I*incr*m;
                r=exp(A*zz*zz+log_Q*(1-s+zz));
                for(j=1;j<=a;j++){
                    r*=GAMMA(gamma[j]*(1-s+zz)+conj(lambda[j]),1/DELTA[j]);
                }
                r/=zz;

                mynorm_r=my_norm(r);
                if(mynorm_r>max_integrand) max_integrand = mynorm_r;
                r1=r2;r2=r3;r3=r4;r4=r5;r5=mynorm_r;
                local_average=(r1+r2+r3+r4+r5)/5;

                //r=r*conj(this->dirichlet_series(conj(1-s+v+I*incr*m),N));
                if(try_use_blfi&&N>999){
                    dirichletseries= this->dirichlet_series_via_blfi(1-s+zz, N,blfi_interval_length);
                }
                else{
                    dirichletseries=0.;
                    for(j=1;j<=N;j++){ //XXXXXXXXXX might be good candidate for unrolling and also decrementing rather than incrementing
                        dirichlet_vector_copy[j]*=dirichlet_multiplier[j];
                        dirichletseries+=dirichlet_vector_copy[j];
                    }
                }
                //dirichletseries=0.;
                //for(j=1;j<=N;j++){
                    //dirichlet_vector_copy[j]*=dirichlet_multiplier[j];
                    //dirichletseries+=dirichlet_vector_copy[j];
                //}
                r*=conj(dirichletseries);
                //cout << "4 dirichletseries: " << dirichletseries << endl;
                //cout << "4 thischletseries: " << this->dirichlet_series(conj(1-s+v+I*incr*m),N) << endl;

                SUM+=r;
                if(my_verbose>2) cout << "#            m= " << m << "  SUM = " << SUM << "  r= " << r  << "  local average = " << local_average << " max=  "<< max_integrand<< endl;
                m--;
            }while(m>-100||local_average>max_integrand*tolerance_sqrd);

            SUM*=incr/(2*Pi);



         }
        else SUM =conj(SUM);

        for(j=1;j<=a;j++){
            r=-gamma[j]-2*real(lambda[j]);
            SUM*=exp(log(DELTA[j])*r);
        }

        //cout << "m= " << m << "  SUM2 = " << SUM << endl;
        //cout << "r= " << r <<  endl;

        L_value=L_value+OMEGA*SUM;

        delete [] dirichlet_vector;
        delete [] dirichlet_vector_copy;
        delete [] dirichlet_multiplier;
        delete [] dirichlet_multiplier_conj;

        //this returns L(s)
        if (!strcmp(return_type,"pure"))
        {
            z=1;
            for(j=1;j<=a;j++)
                z*=GAMMA(gamma[j]*s+lambda[j],DELTA[j]);
            //cout << "pure  " << L_value*exp(-log(Q)*s)/z << endl;
            delete [] DELTA;
            return L_value*exp(-log(Q)*s)/z;
        }

        //returns L(s) rotated to be real on critical line
        //assumes |OMEGA|=1. Valid assumption since
        //LAMBDA(1/2+it) = OMEGA conj(LAMBDA(1/2+it))
        else if (!strcmp(return_type,"rotated pure"))
        {
            r=1;
            for(j=1;j<=a;j++)
                r*=GAMMA(gamma[j]*s+lambda[j],DELTA[j]);
            z=0.;
            for(j=1;j<=a;j++)
                z+=log(DELTA[j])*real(gamma[j]*s+lambda[j]);
            //cout << "rotated pure  " <<  L_value*exp(-log(Q)*real(s)-.5*log(OMEGA))*exp(z)/abs(r) << endl;
            delete [] DELTA;
            return L_value*exp(-log(Q)*real(s)-.5*log(OMEGA))*exp(z)/abs(r);
        }

        //else return Lambda(s) OMEGA^(-1/2) delta^(Re(s)).
        //This returns a real number (though, as a Complex)
        //on the critical line assuming |OMEGA|=1. Valid assumption
        //since LAMBDA(1/2+it) = OMEGA conj(LAMBDA(1/2+it))
        else if(!strcmp(return_type,"normalized and real"))
        {
            z=0.;
            for(j=1;j<=a;j++)
                z+=log(DELTA[j])*real(gamma[j]*s+lambda[j]);
            //cout << "normalized and real  " << L_value*exp(z-.5*log(OMEGA)) << endl;
            delete [] DELTA;
            return L_value*exp(z-.5*log(OMEGA));
        }

        return L_value*exp(-log(Q)*s)/z;
    }


    // implements (3.3.20) with no precomputations.
    // DIGITS is how much precision we would like.
    // DIGITS2 is how much precision (out of DIGITS)
    // we are willing to sacrifice for the sake of
    template <class ttype>
    Complex L_function <ttype>::
    value_via_gamma_sum(Complex s, const char *return_type)
    {
        Complex L_value=0.;
        Double theta=gamma[1];  // equals gamma_1
        Double c;    //controls speed of convergence but at the expense
                     //of loss of precision.
        Complex DELTA;  //(3.3.10)

        Complex u;
        int k;

        c=DIGITS2*log(10.)/theta; //i.e exp(-c theta)=10^(-DIGITS2)
                                  //this sacrifices roughly at most DIGITS2
                                  //out of DIGITS precision.

        //if(abs(t_0)<=2*c/Pi) DELTA=1;
        //else if(t_0>=0) DELTA = exp(I*theta*(Pi/2-c/t_0));
        //else            DELTA = exp(I*theta*(-Pi/2-c/t_0));
        DELTA=find_delta(s*gamma[1]+lambda[1],gamma[1])*exp(2*Pi*I*tweak);


        u=log(DELTA);
        for(k=1;k<=number_of_poles;k++)
           L_value+=residue[k]*exp(-u*pole[k])/(s-pole[k]);



        u=gamma_sum(s, what_type_L, dirichlet_coefficient,
            number_of_dirichlet_coefficients, gamma[1], lambda[1],
            Q, period, DELTA);


        L_value+=exp(log(DELTA/Q)*lambda[1]/gamma[1])*u;


        if(real(s)!=.5)
            u=gamma_sum(1-conj(s), what_type_L, dirichlet_coefficient,
            number_of_dirichlet_coefficients,gamma[1],lambda[1],Q,period,DELTA);
        u=conj(u);


        L_value+=(OMEGA/DELTA)*exp(-log(DELTA*Q)*conj(lambda[1])/gamma[1])*u;

        //this returns L(s)
        if (!strcmp(return_type,"pure"))
        {
            u=log(DELTA/Q)/gamma[1];
//cout << "returning " << L_value << " divided by  " << (GAMMA(gamma[1]*s+lambda[1],exp(u))*exp(u*lambda[1])) << endl;//XXXXXXXXXXXXXXXXXX
            return L_value/(GAMMA(gamma[1]*s+lambda[1],exp(u))*exp(u*lambda[1]));
        }

        //returns L(s) rotated to be real on critical line
        //assumes |OMEGA|=1.
        else if (!strcmp(return_type,"rotated pure"))
        {
            u=log(DELTA/Q)/gamma[1];
            u=abs(GAMMA(gamma[1]*s+lambda[1],exp(u))*exp(u*lambda[1]));
            //u=GAMMA(gamma s + lambda)*(delta/Q)^(-s)
            return L_value*exp(log(DELTA)*real(s)-.5*log(OMEGA))/u;
        }

        //else return Lambda(s) OMEGA^(-1/2) delta^(Re(s)).
        //This returns a real number (though, as a Complex) 
        //on the critical line assuming |OMEGA|=1
        else if(!strcmp(return_type,"normalized and real"))
            return L_value*exp(log(DELTA)*real(s)-.5*log(OMEGA));

        else // return L(s)
        {
            u=log(DELTA/Q)/gamma[1];
            return L_value/(GAMMA(gamma[1]*s+lambda[1],exp(u))*exp(u*lambda[1]));
        }

    }


    template <class ttype>
    Complex L_function <ttype>::
    value(Complex s, int derivative, const char *return_type, const char *method)
    {
      Complex L_value;

      //apply functional equation if Re(s) is < .4 (it's only as we approach 0 that we start to
      //run into problems. So while we could have <.5, no need to do so, and it avoids complications
      //due to things such as .49999999999999
      if(real(s)<.4&& derivative==0){

          //Complex z; Double distance;

          if(my_verbose>2)
              cout << "#            applying functional equation. L(s) -> conj(L(1-conj(s)))." << endl;
          Complex gamma_product=0;
          for(int j=1;j<=this->a;j++){
              //check for trivial zeros
              //z=this->gamma[j]*s+(this->lambda[j]); distance=mynorm(z);

              //else use functional equation and multipy/divide out by Gamma factors
              gamma_product+=log_GAMMA((this->gamma[j])*(1-s)+conj(this->lambda[j]));
              gamma_product-=log_GAMMA((this->gamma[j])*s+(this->lambda[j]));
          }
          if(my_verbose>3)
              cout << "#                gamma product:" << gamma_product << endl;
          //if is a nan we're either at a trivial zero, or not, if a zero produced
          //by the Gamma factors cancels against a pole produced from L(s) and functional equation.
          //Note: here we're in Re(s)<.4, so the only relevant pole is essentially s=0 corresponding
          //to having zeta as a factor. But may as well leave it more general to allow, for example for
          //translations in the t aspect.
          //if(isnan(real(gamma_product))){
          Double rg= Double(real(gamma_product));
          if(rg!=rg){ //i.e. if is nan   NOTE: compiling with --fast-math enabled breaks this comparison, unless one also uses -fno-finite-math-only
               //if(my_verbose>4) cout << "#                    nan detected." << endl;
               if((this->number_of_poles)>0){ //need to check that we're not at pole of Lambda
                   for(int j=1;j<=number_of_poles;j++){
                       if(my_norm(s-this->pole[j])<tolerance_sqrd){ //if we're at a pole of Lambda
                           Complex r=0.,r2,z;
                           for(int i=1;i<=this->a;i++){
                               z=(this->gamma[i])*s+(this->lambda[i]);
                               if(my_norm(z)>tolerance_sqrd) r-=log_GAMMA(z);
                               else {r-=s*log(this->Q); r2=(this->gamma[i]);}
                           }
                           //cout << " j: " << j << endl;
                           //cout << " gamma factor: " << this->gamma[j] << endl;
                           //cout << " residue: " << this->residue[j] << endl;
                           //cout << " exp(r): " << exp(r) << endl;
                           r=exp(r)*r2*(this->residue[j]);
                           return r;
                       }
                   }
               }
               return 0; //else we're at a trivial zero
          }
          //at this point, we're neither a trivial zero or a pole of Lambda(s), and we apply functional equation
          gamma_product=exp(gamma_product);
          gamma_product*=(this->OMEGA)*exp(log(this->Q)*(1-2*s));
          return (gamma_product*conj(this->value(1-conj(s),0,return_type,method))); //XXXXXXX need to do not pure case too
      }

      if(derivative==0){


           //if(only_use_dirichlet_series){
           //    L_value= this->dirichlet_series(s,N_use_dirichlet_series);
           //    return L_value;
           //}
           if(only_use_dirichlet_series){
               if(try_use_blfi){
                   L_value= this->dirichlet_series_via_blfi(s, N_use_dirichlet_series,blfi_interval_length);
                   return L_value;
               }
               else{
                   L_value= this->dirichlet_series(s,N_use_dirichlet_series);
                   return L_value;
               }
           }







          //if(what_type_L==-1&&real(s)==.5&&abs(imag(s))>500) return Zeta(s,return_type);

          //uses Riemann Siegel. This is good only up to limited precision.
          //last condition in the if takes into account that Riemann Sigel is an asymptotic expansion
          //The first 40 remainder terms are coded to an accuracy up to 60 or so (64) digits.
          // This gives roughly O(t^{-20}) accuracy.



          double rs_start_at;
          if(DIGITS<32) rs_start_at = 100.;
          else rs_start_at = 3000.; //for 60 digits, I'll want t>3000 or so

          //Double alpha= sqrt(abs(imag(s))/twoPi);

          //L_value = this->value_via_Riemann_sum(s,return_type);


          //if(what_type_L==-1&&real(s)==.5&&abs(imag(s))>rs_start_at&&DIGITS<66||!strcmp(method,"Riemann Siegel"))){
          if(what_type_L==-1&&real(s)==.5&&abs(imag(s))>rs_start_at&&DIGITS<66&&!strcmp(method,"default")){

               //int success;
               ////Double error_tol=1E-14;



// OLD BLFI... these commented lines need to be replaced with a call to a new rs as below
// For now, I've placed that rs routine in the old_blfi_etc directory
//
//               if(do_blfi&&abs(imag(s))>2e6){
//                   //cout << " rs blfi " ;
//                   L_value= rs(imag(s),tolerance,input_mean_spacing_given,success,return_type);
//                   //output_detail(imag(s));
//
//                   if(success!=1){
//                       cout << "blfi routine failed, calling plain Riemann Siegel." << endl;
//                       do_blfi=false;
//                   }
//
///* XXXXXXXXXXXXX
//
//Complex  blfi_rs = L_value-Zeta(s,return_type);
////Complex  blfi_rs2 = L_value-this->value(s,0,"pure","Gamma sum");
//if(abs(blfi_rs)>1e-10){
//    cout << "LARGE: ";
//    cout<<setprecision(10);
//    cout << imag(s) << " difference:blfi-rs= " << abs(blfi_rs) << endl;
//    //cout << imag(s) << " difference:blfi-gamma sum= " << abs(blfi_rs2) << endl;
//}
//
//XXXXXXXXXXXXX */
//
////blfi_rs = L_value-this->value(s,0,"pure","Gamma sum");
////cout << imag(s) << " difference:blfi-riemann= " << real(blfi_rs) << " " << imag(blfi_rs) << " " << abs(blfi_rs) << endl;
//
//               }
//             else
//             if(!do_blfi||abs(imag(s))<=2e6 ||alpha-floor(alpha)>.9999){
                   //cout << " rs " ;
                   L_value = Zeta(s,return_type);
//             }

               //1.7725 is Pi^(.5), to account for the Q^\pm s in the approximate functional eqn
               DIGITS3=Int((DIGITS-log(log(1.*max_n*1.7725+3)*abs(imag(s))/6.28+3)/2.3)*(4./(4+abs(global_derivative))));
               cout << setprecision(DIGITS3);
               if (my_verbose>5) cout << "#                        Setting output precision to: " << DIGITS3 << endl;
               tolerance3=pow(Double(.1),(DIGITS3+1));
               //cout << setprecision(DIGITS);
               //cout << s << " riemann siegel (" << DIGITS3 << " ): " <<  L_value << endl;
               if (my_verbose>1) cout << "#        L_value via Riemann Siegel: " << L_value << endl;

               return L_value;
          }



         if(a==1&&strcmp(method,"Riemann sum")){ //and the method isn't Riemann sum
              //cout << "gamma sum: " << method <<  endl;


              L_value = this->value_via_gamma_sum(s,return_type);
         }
         else if(a>1||!strcmp(method,"Riemann sum")){
              //cout << "riemann sum: " << method << endl;
              L_value = this->value_via_Riemann_sum(s,return_type);
         }

         DIGITS3=Int( (DIGITS-DIGITS2-log(log(1.*max_n*Q+3)*abs(imag(s))/6.28+3)/2.3));
         cout << setprecision(DIGITS3);

         tolerance3=pow(Double(.1),(DIGITS3+1));

         if(my_verbose>1) cout << "#        calling L:  " << s << " = " << L_value << endl;
         if(my_verbose>3) cout << "#                output precision set to  " << DIGITS3 << " DIGITS" << endl;


         return L_value;
      }
      else if(derivative==-1){ //simple way to use existing framework. -1 stands for logarithmic derivative

          L_value=this->value(s,0,return_type,method);
          return(this->value(s,1,return_type,method)/L_value); //order, i.e. value then derivative, is important since
                                                        //derivative sets output to lower precision so should be called 2nd

      }
      // approximates derivatives using central differences with O(h^4)
      // by taking linear combinations of f(s+mh) = sum f[k](s)/k! (mh)^k, i.e.
      // Taylor expansion around s, f[k] is kth derivative of f.
      // We take the linear combination, with some choices of m,
      // of this that gives f[n](s) + O(h^4).
      else if(derivative>0&&derivative<26){

          Complex der=0,z,h;
          int DIGITStmp=1;
          int DIGITS2tmp=DIGITS2;
          DIGITS2=1;



          DIGITS3=Int( (DIGITS-DIGITS2-log(log(1.*max_n*Q+3)*abs(imag(s))/6.28+3)/2.3))+2;

          switch(derivative)
          {
             case 1:
                 //max of the partial sums of the coeffs is       0.6
                 DIGITStmp=DIGITS3+1;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-2*h,0,return_type,method); der+=1*z/12;
                 z=this->value(s-h,0,return_type,method); der+=-2*z/3;
                 z=this->value(s+h,0,return_type,method); der+=2*z/3;
                 z=this->value(s+2*h,0,return_type,method); der+=-1*z/12;
                 break;
             case 2:
                 //max of the partial sums of the coeffs is       1.2
                 DIGITStmp=DIGITS3+1;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-2*h,0,return_type,method); der+=-1*z/12;
                 z=this->value(s-h,0,return_type,method); der+=4*z/3;
                 z=this->value(s,0,return_type,method); der+=-5*z/2;
                 z=this->value(s+h,0,return_type,method); der+=4*z/3;
                 z=this->value(s+2*h,0,return_type,method); der+=-1*z/12;
                 break;
             case 3:
                 //max of the partial sums of the coeffs is       0.9
                 DIGITStmp=DIGITS3+1;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-3*h,0,return_type,method); der+=1*z/8;
                 z=this->value(s-2*h,0,return_type,method); der+=-1*z;
                 z=this->value(s-h,0,return_type,method); der+=13*z/8;
                 z=this->value(s+h,0,return_type,method); der+=-13*z/8;
                 z=this->value(s+2*h,0,return_type,method); der+=1*z;
                 z=this->value(s+3*h,0,return_type,method); der+=-1*z/8;
                 break;
             case 4:
                 //max of the partial sums of the coeffs is       4.7
                 DIGITStmp=DIGITS3+1;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-3*h,0,return_type,method); der+=-1*z/6;
                 z=this->value(s-2*h,0,return_type,method); der+=2*z;
                 z=this->value(s-h,0,return_type,method); der+=-13*z/2;
                 z=this->value(s,0,return_type,method); der+=28*z/3;
                 z=this->value(s+h,0,return_type,method); der+=-13*z/2;
                 z=this->value(s+2*h,0,return_type,method); der+=2*z;
                 z=this->value(s+3*h,0,return_type,method); der+=-1*z/6;
                 break;
             case 5:
                 //max of the partial sums of the coeffs is       3.0
                 DIGITStmp=DIGITS3+1;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-4*h,0,return_type,method); der+=1*z/6;
                 z=this->value(s-3*h,0,return_type,method); der+=-3*z/2;
                 z=this->value(s-2*h,0,return_type,method); der+=13*z/3;
                 z=this->value(s-h,0,return_type,method); der+=-29*z/6;
                 z=this->value(s+h,0,return_type,method); der+=29*z/6;
                 z=this->value(s+2*h,0,return_type,method); der+=-13*z/3;
                 z=this->value(s+3*h,0,return_type,method); der+=3*z/2;
                 z=this->value(s+4*h,0,return_type,method); der+=-1*z/6;
                 break;
             case 6:
                 //max of the partial sums of the coeffs is      18.8
                 DIGITStmp=DIGITS3;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-4*h,0,return_type,method); der+=-1*z/4;
                 z=this->value(s-3*h,0,return_type,method); der+=3*z;
                 z=this->value(s-2*h,0,return_type,method); der+=-13*z;
                 z=this->value(s-h,0,return_type,method); der+=29*z;
                 z=this->value(s,0,return_type,method); der+=-75*z/2;
                 z=this->value(s+h,0,return_type,method); der+=29*z;
                 z=this->value(s+2*h,0,return_type,method); der+=-13*z;
                 z=this->value(s+3*h,0,return_type,method); der+=3*z;
                 z=this->value(s+4*h,0,return_type,method); der+=-1*z/4;
                 break;
             case 7:
                 //max of the partial sums of the coeffs is      10.3
                 DIGITStmp=DIGITS3;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-5*h,0,return_type,method); der+=5*z/24;
                 z=this->value(s-4*h,0,return_type,method); der+=-13*z/6;
                 z=this->value(s-3*h,0,return_type,method); der+=69*z/8;
                 z=this->value(s-2*h,0,return_type,method); der+=-17*z;
                 z=this->value(s-h,0,return_type,method); der+=63*z/4;
                 z=this->value(s+h,0,return_type,method); der+=-63*z/4;
                 z=this->value(s+2*h,0,return_type,method); der+=17*z;
                 z=this->value(s+3*h,0,return_type,method); der+=-69*z/8;
                 z=this->value(s+4*h,0,return_type,method); der+=13*z/6;
                 z=this->value(s+5*h,0,return_type,method); der+=-5*z/24;
                 break;
             case 8:
                 //max of the partial sums of the coeffs is      77.0
                 DIGITStmp=DIGITS3;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-5*h,0,return_type,method); der+=-1*z/3;
                 z=this->value(s-4*h,0,return_type,method); der+=13*z/3;
                 z=this->value(s-3*h,0,return_type,method); der+=-23*z;
                 z=this->value(s-2*h,0,return_type,method); der+=68*z;
                 z=this->value(s-h,0,return_type,method); der+=-126*z;
                 z=this->value(s,0,return_type,method); der+=154*z;
                 z=this->value(s+h,0,return_type,method); der+=-126*z;
                 z=this->value(s+2*h,0,return_type,method); der+=68*z;
                 z=this->value(s+3*h,0,return_type,method); der+=-23*z;
                 z=this->value(s+4*h,0,return_type,method); der+=13*z/3;
                 z=this->value(s+5*h,0,return_type,method); der+=-1*z/3;
                 break;
             case 9:
                 //max of the partial sums of the coeffs is      36.5
                 DIGITStmp=DIGITS3;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-6*h,0,return_type,method); der+=1*z/4;
                 z=this->value(s-5*h,0,return_type,method); der+=-3*z;
                 z=this->value(s-4*h,0,return_type,method); der+=15*z;
                 z=this->value(s-3*h,0,return_type,method); der+=-41*z;
                 z=this->value(s-2*h,0,return_type,method); der+=261*z/4;
                 z=this->value(s-h,0,return_type,method); der+=-54*z;
                 z=this->value(s+h,0,return_type,method); der+=54*z;
                 z=this->value(s+2*h,0,return_type,method); der+=-261*z/4;
                 z=this->value(s+3*h,0,return_type,method); der+=41*z;
                 z=this->value(s+4*h,0,return_type,method); der+=-15*z;
                 z=this->value(s+5*h,0,return_type,method); der+=3*z;
                 z=this->value(s+6*h,0,return_type,method); der+=-1*z/4;
                 break;
             case 10:
                 //max of the partial sums of the coeffs is     318.5
                 DIGITStmp=DIGITS3-1;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-6*h,0,return_type,method); der+=-5*z/12;
                 z=this->value(s-5*h,0,return_type,method); der+=6*z;
                 z=this->value(s-4*h,0,return_type,method); der+=-75*z/2;
                 z=this->value(s-3*h,0,return_type,method); der+=410*z/3;
                 z=this->value(s-2*h,0,return_type,method); der+=-1305*z/4;
                 z=this->value(s-h,0,return_type,method); der+=540*z;
                 z=this->value(s,0,return_type,method); der+=-637*z;
                 z=this->value(s+h,0,return_type,method); der+=540*z;
                 z=this->value(s+2*h,0,return_type,method); der+=-1305*z/4;
                 z=this->value(s+3*h,0,return_type,method); der+=410*z/3;
                 z=this->value(s+4*h,0,return_type,method); der+=-75*z/2;
                 z=this->value(s+5*h,0,return_type,method); der+=6*z;
                 z=this->value(s+6*h,0,return_type,method); der+=-5*z/12;
                 break;
             case 11:
                 //max of the partial sums of the coeffs is     131.6
                 DIGITStmp=DIGITS3-1;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-7*h,0,return_type,method); der+=7*z/24;
                 z=this->value(s-6*h,0,return_type,method); der+=-4*z;
                 z=this->value(s-5*h,0,return_type,method); der+=575*z/24;
                 z=this->value(s-4*h,0,return_type,method); der+=-248*z/3;
                 z=this->value(s-3*h,0,return_type,method); der+=1441*z/8;
                 z=this->value(s-2*h,0,return_type,method); der+=-748*z/3;
                 z=this->value(s-h,0,return_type,method); der+=1529*z/8;
                 z=this->value(s+h,0,return_type,method); der+=-1529*z/8;
                 z=this->value(s+2*h,0,return_type,method); der+=748*z/3;
                 z=this->value(s+3*h,0,return_type,method); der+=-1441*z/8;
                 z=this->value(s+4*h,0,return_type,method); der+=248*z/3;
                 z=this->value(s+5*h,0,return_type,method); der+=-575*z/24;
                 z=this->value(s+6*h,0,return_type,method); der+=4*z;
                 z=this->value(s+7*h,0,return_type,method); der+=-7*z/24;
                 break;
             case 12:
                 //max of the partial sums of the coeffs is    1320.0
                 DIGITStmp=DIGITS3-2;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-7*h,0,return_type,method); der+=-1*z/2;
                 z=this->value(s-6*h,0,return_type,method); der+=8*z;
                 z=this->value(s-5*h,0,return_type,method); der+=-115*z/2;
                 z=this->value(s-4*h,0,return_type,method); der+=248*z;
                 z=this->value(s-3*h,0,return_type,method); der+=-1441*z/2;
                 z=this->value(s-2*h,0,return_type,method); der+=1496*z;
                 z=this->value(s-h,0,return_type,method); der+=-4587*z/2;
                 z=this->value(s,0,return_type,method); der+=2640*z;
                 z=this->value(s+h,0,return_type,method); der+=-4587*z/2;
                 z=this->value(s+2*h,0,return_type,method); der+=1496*z;
                 z=this->value(s+3*h,0,return_type,method); der+=-1441*z/2;
                 z=this->value(s+4*h,0,return_type,method); der+=248*z;
                 z=this->value(s+5*h,0,return_type,method); der+=-115*z/2;
                 z=this->value(s+6*h,0,return_type,method); der+=8*z;
                 z=this->value(s+7*h,0,return_type,method); der+=-1*z/2;
                 break;
             case 13:
                 //max of the partial sums of the coeffs is     482.2
                 DIGITStmp=DIGITS3-1;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-8*h,0,return_type,method); der+=1*z/3;
                 z=this->value(s-7*h,0,return_type,method); der+=-31*z/6;
                 z=this->value(s-6*h,0,return_type,method); der+=36*z;
                 z=this->value(s-5*h,0,return_type,method); der+=-895*z/6;
                 z=this->value(s-4*h,0,return_type,method); der+=1222*z/3;
                 z=this->value(s-3*h,0,return_type,method); der+=-1521*z/2;
                 z=this->value(s-2*h,0,return_type,method); der+=2860*z/3;
                 z=this->value(s-h,0,return_type,method); der+=-4147*z/6;
                 z=this->value(s+h,0,return_type,method); der+=4147*z/6;
                 z=this->value(s+2*h,0,return_type,method); der+=-2860*z/3;
                 z=this->value(s+3*h,0,return_type,method); der+=1521*z/2;
                 z=this->value(s+4*h,0,return_type,method); der+=-1222*z/3;
                 z=this->value(s+5*h,0,return_type,method); der+=895*z/6;
                 z=this->value(s+6*h,0,return_type,method); der+=-36*z;
                 z=this->value(s+7*h,0,return_type,method); der+=31*z/6;
                 z=this->value(s+8*h,0,return_type,method); der+=-1*z/3;
                 break;
             case 14:
                 //max of the partial sums of the coeffs is    5469.8
                 DIGITStmp=DIGITS3-2;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-8*h,0,return_type,method); der+=-7*z/12;
                 z=this->value(s-7*h,0,return_type,method); der+=31*z/3;
                 z=this->value(s-6*h,0,return_type,method); der+=-84*z;
                 z=this->value(s-5*h,0,return_type,method); der+=1253*z/3;
                 z=this->value(s-4*h,0,return_type,method); der+=-4277*z/3;
                 z=this->value(s-3*h,0,return_type,method); der+=3549*z;
                 z=this->value(s-2*h,0,return_type,method); der+=-20020*z/3;
                 z=this->value(s-h,0,return_type,method); der+=29029*z/3;
                 z=this->value(s,0,return_type,method); der+=-21879*z/2;
                 z=this->value(s+h,0,return_type,method); der+=29029*z/3;
                 z=this->value(s+2*h,0,return_type,method); der+=-20020*z/3;
                 z=this->value(s+3*h,0,return_type,method); der+=3549*z;
                 z=this->value(s+4*h,0,return_type,method); der+=-4277*z/3;
                 z=this->value(s+5*h,0,return_type,method); der+=1253*z/3;
                 z=this->value(s+6*h,0,return_type,method); der+=-84*z;
                 z=this->value(s+7*h,0,return_type,method); der+=31*z/3;
                 z=this->value(s+8*h,0,return_type,method); der+=-7*z/12;
                 break;
             case 15:
                 //max of the partial sums of the coeffs is    1865.5
                 DIGITStmp=DIGITS3-2;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-9*h,0,return_type,method); der+=3*z/8;
                 z=this->value(s-8*h,0,return_type,method); der+=-13*z/2;
                 z=this->value(s-7*h,0,return_type,method); der+=413*z/8;
                 z=this->value(s-6*h,0,return_type,method); der+=-249*z;
                 z=this->value(s-5*h,0,return_type,method); der+=1625*z/2;
                 z=this->value(s-4*h,0,return_type,method); der+=-1883*z;
                 z=this->value(s-3*h,0,return_type,method); der+=6279*z/2;
                 z=this->value(s-2*h,0,return_type,method); der+=-3653*z;
                 z=this->value(s-h,0,return_type,method); der+=10153*z/4;
                 z=this->value(s+h,0,return_type,method); der+=-10153*z/4;
                 z=this->value(s+2*h,0,return_type,method); der+=3653*z;
                 z=this->value(s+3*h,0,return_type,method); der+=-6279*z/2;
                 z=this->value(s+4*h,0,return_type,method); der+=1883*z;
                 z=this->value(s+5*h,0,return_type,method); der+=-1625*z/2;
                 z=this->value(s+6*h,0,return_type,method); der+=249*z;
                 z=this->value(s+7*h,0,return_type,method); der+=-413*z/8;
                 z=this->value(s+8*h,0,return_type,method); der+=13*z/2;
                 z=this->value(s+9*h,0,return_type,method); der+=-3*z/8;
                 break;
             case 16:
                 //max of the partial sums of the coeffs is   22641.7
                 DIGITStmp=DIGITS3-3;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-9*h,0,return_type,method); der+=-2*z/3;
                 z=this->value(s-8*h,0,return_type,method); der+=13*z;
                 z=this->value(s-7*h,0,return_type,method); der+=-118*z;
                 z=this->value(s-6*h,0,return_type,method); der+=664*z;
                 z=this->value(s-5*h,0,return_type,method); der+=-2600*z;
                 z=this->value(s-4*h,0,return_type,method); der+=7532*z;
                 z=this->value(s-3*h,0,return_type,method); der+=-16744*z;
                 z=this->value(s-2*h,0,return_type,method); der+=29224*z;
                 z=this->value(s-h,0,return_type,method); der+=-40612*z;
                 z=this->value(s,0,return_type,method); der+=135850*z/3;
                 z=this->value(s+h,0,return_type,method); der+=-40612*z;
                 z=this->value(s+2*h,0,return_type,method); der+=29224*z;
                 z=this->value(s+3*h,0,return_type,method); der+=-16744*z;
                 z=this->value(s+4*h,0,return_type,method); der+=7532*z;
                 z=this->value(s+5*h,0,return_type,method); der+=-2600*z;
                 z=this->value(s+6*h,0,return_type,method); der+=664*z;
                 z=this->value(s+7*h,0,return_type,method); der+=-118*z;
                 z=this->value(s+8*h,0,return_type,method); der+=13*z;
                 z=this->value(s+9*h,0,return_type,method); der+=-2*z/3;
                 break;
             case 17:
                 //max of the partial sums of the coeffs is    7345.0
                 DIGITStmp=DIGITS3-2;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-10*h,0,return_type,method); der+=5*z/12;
                 z=this->value(s-9*h,0,return_type,method); der+=-8*z;
                 z=this->value(s-8*h,0,return_type,method); der+=214*z/3;
                 z=this->value(s-7*h,0,return_type,method); der+=-392*z;
                 z=this->value(s-6*h,0,return_type,method); der+=5933*z/4;
                 z=this->value(s-5*h,0,return_type,method); der+=-4080*z;
                 z=this->value(s-4*h,0,return_type,method); der+=8364*z;
                 z=this->value(s-3*h,0,return_type,method); der+=-12784*z;
                 z=this->value(s-2*h,0,return_type,method); der+=28067*z/2;
                 z=this->value(s-h,0,return_type,method); der+=-28288*z/3;
                 z=this->value(s+h,0,return_type,method); der+=28288*z/3;
                 z=this->value(s+2*h,0,return_type,method); der+=-28067*z/2;
                 z=this->value(s+3*h,0,return_type,method); der+=12784*z;
                 z=this->value(s+4*h,0,return_type,method); der+=-8364*z;
                 z=this->value(s+5*h,0,return_type,method); der+=4080*z;
                 z=this->value(s+6*h,0,return_type,method); der+=-5933*z/4;
                 z=this->value(s+7*h,0,return_type,method); der+=392*z;
                 z=this->value(s+8*h,0,return_type,method); der+=-214*z/3;
                 z=this->value(s+9*h,0,return_type,method); der+=8*z;
                 z=this->value(s+10*h,0,return_type,method); der+=-5*z/12;
                 break;
             case 18:
                 //max of the partial sums of the coeffs is   93593.5
                 DIGITStmp=DIGITS3-3;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-10*h,0,return_type,method); der+=-3*z/4;
                 z=this->value(s-9*h,0,return_type,method); der+=16*z;
                 z=this->value(s-8*h,0,return_type,method); der+=-321*z/2;
                 z=this->value(s-7*h,0,return_type,method); der+=1008*z;
                 z=this->value(s-6*h,0,return_type,method); der+=-17799*z/4;
                 z=this->value(s-5*h,0,return_type,method); der+=14688*z;
                 z=this->value(s-4*h,0,return_type,method); der+=-37638*z;
                 z=this->value(s-3*h,0,return_type,method); der+=76704*z;
                 z=this->value(s-2*h,0,return_type,method); der+=-252603*z/2;
                 z=this->value(s-h,0,return_type,method); der+=169728*z;
                 z=this->value(s,0,return_type,method); der+=-187187*z;
                 z=this->value(s+h,0,return_type,method); der+=169728*z;
                 z=this->value(s+2*h,0,return_type,method); der+=-252603*z/2;
                 z=this->value(s+3*h,0,return_type,method); der+=76704*z;
                 z=this->value(s+4*h,0,return_type,method); der+=-37638*z;
                 z=this->value(s+5*h,0,return_type,method); der+=14688*z;
                 z=this->value(s+6*h,0,return_type,method); der+=-17799*z/4;
                 z=this->value(s+7*h,0,return_type,method); der+=1008*z;
                 z=this->value(s+8*h,0,return_type,method); der+=-321*z/2;
                 z=this->value(s+9*h,0,return_type,method); der+=16*z;
                 z=this->value(s+10*h,0,return_type,method); der+=-3*z/4;
                 break;
             case 19:
                 //max of the partial sums of the coeffs is   28836.2
                 DIGITStmp=DIGITS3-3;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-11*h,0,return_type,method); der+=11*z/24;
                 z=this->value(s-10*h,0,return_type,method); der+=-29*z/3;
                 z=this->value(s-9*h,0,return_type,method); der+=765*z/8;
                 z=this->value(s-8*h,0,return_type,method); der+=-1768*z/3;
                 z=this->value(s-7*h,0,return_type,method); der+=60781*z/24;
                 z=this->value(s-6*h,0,return_type,method); der+=-8037*z;
                 z=this->value(s-5*h,0,return_type,method); der+=155363*z/8;
                 z=this->value(s-4*h,0,return_type,method); der+=-36176*z;
                 z=this->value(s-3*h,0,return_type,method); der+=206397*z/4;
                 z=this->value(s-2*h,0,return_type,method); der+=-162146*z/3;
                 z=this->value(s-h,0,return_type,method); der+=424099*z/12;
                 z=this->value(s+h,0,return_type,method); der+=-424099*z/12;
                 z=this->value(s+2*h,0,return_type,method); der+=162146*z/3;
                 z=this->value(s+3*h,0,return_type,method); der+=-206397*z/4;
                 z=this->value(s+4*h,0,return_type,method); der+=36176*z;
                 z=this->value(s+5*h,0,return_type,method); der+=-155363*z/8;
                 z=this->value(s+6*h,0,return_type,method); der+=8037*z;
                 z=this->value(s+7*h,0,return_type,method); der+=-60781*z/24;
                 z=this->value(s+8*h,0,return_type,method); der+=1768*z/3;
                 z=this->value(s+9*h,0,return_type,method); der+=-765*z/8;
                 z=this->value(s+10*h,0,return_type,method); der+=29*z/3;
                 z=this->value(s+11*h,0,return_type,method); der+=-11*z/24;
                 break;
             case 20:
                 //max of the partial sums of the coeffs is  386308.0
                 DIGITStmp=DIGITS3-4;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-11*h,0,return_type,method); der+=-5*z/6;
                 z=this->value(s-10*h,0,return_type,method); der+=58*z/3;
                 z=this->value(s-9*h,0,return_type,method); der+=-425*z/2;
                 z=this->value(s-8*h,0,return_type,method); der+=4420*z/3;
                 z=this->value(s-7*h,0,return_type,method); der+=-43415*z/6;
                 z=this->value(s-6*h,0,return_type,method); der+=26790*z;
                 z=this->value(s-5*h,0,return_type,method); der+=-155363*z/2;
                 z=this->value(s-4*h,0,return_type,method); der+=180880*z;
                 z=this->value(s-3*h,0,return_type,method); der+=-343995*z;
                 z=this->value(s-2*h,0,return_type,method); der+=1621460*z/3;
                 z=this->value(s-h,0,return_type,method); der+=-2120495*z/3;
                 z=this->value(s,0,return_type,method); der+=772616*z;
                 z=this->value(s+h,0,return_type,method); der+=-2120495*z/3;
                 z=this->value(s+2*h,0,return_type,method); der+=1621460*z/3;
                 z=this->value(s+3*h,0,return_type,method); der+=-343995*z;
                 z=this->value(s+4*h,0,return_type,method); der+=180880*z;
                 z=this->value(s+5*h,0,return_type,method); der+=-155363*z/2;
                 z=this->value(s+6*h,0,return_type,method); der+=26790*z;
                 z=this->value(s+7*h,0,return_type,method); der+=-43415*z/6;
                 z=this->value(s+8*h,0,return_type,method); der+=4420*z/3;
                 z=this->value(s+9*h,0,return_type,method); der+=-425*z/2;
                 z=this->value(s+10*h,0,return_type,method); der+=58*z/3;
                 z=this->value(s+11*h,0,return_type,method); der+=-5*z/6;
                 break;
             case 21:
                 //max of the partial sums of the coeffs is  113050.0
                 DIGITStmp=DIGITS3-4;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-12*h,0,return_type,method); der+=1*z/2;
                 z=this->value(s-11*h,0,return_type,method); der+=-23*z/2;
                 z=this->value(s-10*h,0,return_type,method); der+=125*z;
                 z=this->value(s-9*h,0,return_type,method); der+=-1707*z/2;
                 z=this->value(s-8*h,0,return_type,method); der+=4102*z;
                 z=this->value(s-7*h,0,return_type,method); der+=-29449*z/2;
                 z=this->value(s-6*h,0,return_type,method); der+=40831*z;
                 z=this->value(s-5*h,0,return_type,method); der+=-178125*z/2;
                 z=this->value(s-4*h,0,return_type,method); der+=307173*z/2;
                 z=this->value(s-3*h,0,return_type,method); der+=-207043*z;
                 z=this->value(s-2*h,0,return_type,method); der+=208658*z;
                 z=this->value(s-h,0,return_type,method); der+=-133399*z;
                 z=this->value(s+h,0,return_type,method); der+=133399*z;
                 z=this->value(s+2*h,0,return_type,method); der+=-208658*z;
                 z=this->value(s+3*h,0,return_type,method); der+=207043*z;
                 z=this->value(s+4*h,0,return_type,method); der+=-307173*z/2;
                 z=this->value(s+5*h,0,return_type,method); der+=178125*z/2;
                 z=this->value(s+6*h,0,return_type,method); der+=-40831*z;
                 z=this->value(s+7*h,0,return_type,method); der+=29449*z/2;
                 z=this->value(s+8*h,0,return_type,method); der+=-4102*z;
                 z=this->value(s+9*h,0,return_type,method); der+=1707*z/2;
                 z=this->value(s+10*h,0,return_type,method); der+=-125*z;
                 z=this->value(s+11*h,0,return_type,method); der+=23*z/2;
                 z=this->value(s+12*h,0,return_type,method); der+=-1*z/2;
                 break;
             case 22:
                 //max of the partial sums of the coeffs is 1592120.8
                 DIGITStmp=DIGITS3-5;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-12*h,0,return_type,method); der+=-11*z/12;
                 z=this->value(s-11*h,0,return_type,method); der+=23*z;
                 z=this->value(s-10*h,0,return_type,method); der+=-275*z;
                 z=this->value(s-9*h,0,return_type,method); der+=6259*z/3;
                 z=this->value(s-8*h,0,return_type,method); der+=-22561*z/2;
                 z=this->value(s-7*h,0,return_type,method); der+=46277*z;
                 z=this->value(s-6*h,0,return_type,method); der+=-449141*z/3;
                 z=this->value(s-5*h,0,return_type,method); der+=391875*z;
                 z=this->value(s-4*h,0,return_type,method); der+=-3378903*z/4;
                 z=this->value(s-3*h,0,return_type,method); der+=4554946*z/3;
                 z=this->value(s-2*h,0,return_type,method); der+=-2295238*z;
                 z=this->value(s-h,0,return_type,method); der+=2934778*z;
                 z=this->value(s,0,return_type,method); der+=-9552725*z/3;
                 z=this->value(s+h,0,return_type,method); der+=2934778*z;
                 z=this->value(s+2*h,0,return_type,method); der+=-2295238*z;
                 z=this->value(s+3*h,0,return_type,method); der+=4554946*z/3;
                 z=this->value(s+4*h,0,return_type,method); der+=-3378903*z/4;
                 z=this->value(s+5*h,0,return_type,method); der+=391875*z;
                 z=this->value(s+6*h,0,return_type,method); der+=-449141*z/3;
                 z=this->value(s+7*h,0,return_type,method); der+=46277*z;
                 z=this->value(s+8*h,0,return_type,method); der+=-22561*z/2;
                 z=this->value(s+9*h,0,return_type,method); der+=6259*z/3;
                 z=this->value(s+10*h,0,return_type,method); der+=-275*z;
                 z=this->value(s+11*h,0,return_type,method); der+=23*z;
                 z=this->value(s+12*h,0,return_type,method); der+=-11*z/12;
                 break;
             case 23:
                 //max of the partial sums of the coeffs is  442940.7
                 DIGITStmp=DIGITS3-4;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-13*h,0,return_type,method); der+=13*z/24;
                 z=this->value(s-12*h,0,return_type,method); der+=-27*z/2;
                 z=this->value(s-11*h,0,return_type,method); der+=3839*z/24;
                 z=this->value(s-10*h,0,return_type,method); der+=-3595*z/3;
                 z=this->value(s-9*h,0,return_type,method); der+=25461*z/4;
                 z=this->value(s-8*h,0,return_type,method); der+=-76406*z/3;
                 z=this->value(s-7*h,0,return_type,method); der+=954569*z/12;
                 z=this->value(s-6*h,0,return_type,method); der+=-198099*z;
                 z=this->value(s-5*h,0,return_type,method); der+=9541895*z/24;
                 z=this->value(s-4*h,0,return_type,method); der+=-3860021*z/6;
                 z=this->value(s-3*h,0,return_type,method); der+=6619239*z/8;
                 z=this->value(s-2*h,0,return_type,method); der+=-2421854*z/3;
                 z=this->value(s-h,0,return_type,method); der+=3038461*z/6;
                 z=this->value(s+h,0,return_type,method); der+=-3038461*z/6;
                 z=this->value(s+2*h,0,return_type,method); der+=2421854*z/3;
                 z=this->value(s+3*h,0,return_type,method); der+=-6619239*z/8;
                 z=this->value(s+4*h,0,return_type,method); der+=3860021*z/6;
                 z=this->value(s+5*h,0,return_type,method); der+=-9541895*z/24;
                 z=this->value(s+6*h,0,return_type,method); der+=198099*z;
                 z=this->value(s+7*h,0,return_type,method); der+=-954569*z/12;
                 z=this->value(s+8*h,0,return_type,method); der+=76406*z/3;
                 z=this->value(s+9*h,0,return_type,method); der+=-25461*z/4;
                 z=this->value(s+10*h,0,return_type,method); der+=3595*z/3;
                 z=this->value(s+11*h,0,return_type,method); der+=-3839*z/24;
                 z=this->value(s+12*h,0,return_type,method); der+=27*z/2;
                 z=this->value(s+13*h,0,return_type,method); der+=-13*z/24;
                 break;
             case 24:
                 //max of the partial sums of the coeffs is 6552378.0
                 DIGITStmp=DIGITS3-5;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-13*h,0,return_type,method); der+=-1*z;
                 z=this->value(s-12*h,0,return_type,method); der+=27*z;
                 z=this->value(s-11*h,0,return_type,method); der+=-349*z;
                 z=this->value(s-10*h,0,return_type,method); der+=2876*z;
                 z=this->value(s-9*h,0,return_type,method); der+=-16974*z;
                 z=this->value(s-8*h,0,return_type,method); der+=76406*z;
                 z=this->value(s-7*h,0,return_type,method); der+=-272734*z;
                 z=this->value(s-6*h,0,return_type,method); der+=792396*z;
                 z=this->value(s-5*h,0,return_type,method); der+=-1908379*z;
                 z=this->value(s-4*h,0,return_type,method); der+=3860021*z;
                 z=this->value(s-3*h,0,return_type,method); der+=-6619239*z;
                 z=this->value(s-2*h,0,return_type,method); der+=9687416*z;
                 z=this->value(s-h,0,return_type,method); der+=-12153844*z;
                 z=this->value(s,0,return_type,method); der+=13104756*z;
                 z=this->value(s+h,0,return_type,method); der+=-12153844*z;
                 z=this->value(s+2*h,0,return_type,method); der+=9687416*z;
                 z=this->value(s+3*h,0,return_type,method); der+=-6619239*z;
                 z=this->value(s+4*h,0,return_type,method); der+=3860021*z;
                 z=this->value(s+5*h,0,return_type,method); der+=-1908379*z;
                 z=this->value(s+6*h,0,return_type,method); der+=792396*z;
                 z=this->value(s+7*h,0,return_type,method); der+=-272734*z;
                 z=this->value(s+8*h,0,return_type,method); der+=76406*z;
                 z=this->value(s+9*h,0,return_type,method); der+=-16974*z;
                 z=this->value(s+10*h,0,return_type,method); der+=2876*z;
                 z=this->value(s+11*h,0,return_type,method); der+=-349*z;
                 z=this->value(s+12*h,0,return_type,method); der+=27*z;
                 z=this->value(s+13*h,0,return_type,method); der+=-1*z;
                 break;
             case 25:
                 //max of the partial sums of the coeffs is 1735290.6
                 DIGITStmp=DIGITS3-5;
                 if(DIGITStmp<1)DIGITStmp=1;
                 h=I*pow(Double(.1),DIGITStmp/(derivative+4.));
                 z=this->value(s-14*h,0,return_type,method); der+=7*z/12;
                 z=this->value(s-13*h,0,return_type,method); der+=-47*z/3;
                 z=this->value(s-12*h,0,return_type,method); der+=201*z;
                 z=this->value(s-11*h,0,return_type,method); der+=-1639*z;
                 z=this->value(s-10*h,0,return_type,method); der+=38125*z/4;
                 z=this->value(s-9*h,0,return_type,method); der+=-42030*z;
                 z=this->value(s-8*h,0,return_type,method); der+=145820*z;
                 z=this->value(s-7*h,0,return_type,method); der+=-407330*z;
                 z=this->value(s-6*h,0,return_type,method); der+=3715305*z/4;
                 z=this->value(s-5*h,0,return_type,method); der+=-1739375*z;
                 z=this->value(s-4*h,0,return_type,method); der+=2667885*z;
                 z=this->value(s-3*h,0,return_type,method); der+=-3297165*z;
                 z=this->value(s-2*h,0,return_type,method); der+=12517865*z/4;
                 z=this->value(s-h,0,return_type,method); der+=-1931540*z;
                 z=this->value(s+h,0,return_type,method); der+=1931540*z;
                 z=this->value(s+2*h,0,return_type,method); der+=-12517865*z/4;
                 z=this->value(s+3*h,0,return_type,method); der+=3297165*z;
                 z=this->value(s+4*h,0,return_type,method); der+=-2667885*z;
                 z=this->value(s+5*h,0,return_type,method); der+=1739375*z;
                 z=this->value(s+6*h,0,return_type,method); der+=-3715305*z/4;
                 z=this->value(s+7*h,0,return_type,method); der+=407330*z;
                 z=this->value(s+8*h,0,return_type,method); der+=-145820*z;
                 z=this->value(s+9*h,0,return_type,method); der+=42030*z;
                 z=this->value(s+10*h,0,return_type,method); der+=-38125*z/4;
                 z=this->value(s+11*h,0,return_type,method); der+=1639*z;
                 z=this->value(s+12*h,0,return_type,method); der+=-201*z;
                 z=this->value(s+13*h,0,return_type,method); der+=47*z/3;
                 z=this->value(s+14*h,0,return_type,method); der+=-7*z/12;
                 break;
          }

          //the factor 4./(4+derivative) is to account for loss of precision in our scheme
          //for numerical differentiation. DIGITS4 accounts for the maximum of the partial sums
          //of the coefficients
          DIGITS3= (DIGITStmp*4)/(4+derivative)-1;
          if(derivative>1) DIGITS3-=1;
          if(DIGITS3<1) DIGITS3=1;
          cout << setprecision(DIGITS3);
          if (my_verbose>1) cout << "#        Setting L derivative output precision to: " << DIGITS3 << endl;

          //tolerance3=pow(1./10,(DIGITS3+1));

          DIGITS2=DIGITS2tmp;

          return(der/pow(h,derivative));
      }
      else{
          cout << "Error. Specified derivative must be >= -1 and <= 25" << endl;
          exit(1);
      }

    }
