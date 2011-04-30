#ifndef Ldokchitser_H
#define Ldokchitser_H

//finding the explicit taylor series for \phi(t) using Dokchitser algo

#define MYDIGITS 5 // estimate of precision ... will be set using the precision variable

void mult_poly_taylor(Complex *, Complex *, Complex *, int );

template <class ttype>
void L_function <ttype>::
phi_series(int precision)
{
	cout << "-----------------------------------------------"<< endl << endl;
	cout << "phi series for " << name << " L_function" << endl << endl;
	
	int j,k;
		
	// constructing the equivalence classes Lambda[k] for k = 1 to N
	
	int pordtmp[a+1];
	Complex diff;
	Complex *lambda_k;
	int *l;
	
	for (j=1;j<=a;j++)
    	pordtmp[j]= 1;
    	
    for (j=1;j<=a;j++)
        for (k=1;k<=a;k++)
            if (j != k)
            {    
                diff = 2*(lambda[j] - lambda[k]);
                if((imag(diff)==0) && (fmod(real(diff),Double(2)) == 0) && (real(diff)<=0))
                {
                    pordtmp[j]+=pordtmp[k];
                    pordtmp[k]=0;
                }
            }
        	
    Complex temp_lambda_k[a];
    int temp_l[a];
    j=1;
    for (k=1;k<=a;k++)
        if (pordtmp[k]!=0)
        {
            temp_lambda_k[j]= lambda[k];
            temp_l[j]=pordtmp[k];
            j++;
        }
    
    int N = j-1;   
    lambda_k = new Complex[N+1];
    l = new int[N+1];
    for (j=1;j<=N;j++)
    {
    	lambda_k[j] = temp_lambda_k[j];
    	l[j] = temp_l[j];
    }
	
	cout << "-----------------------------------------------"<< endl << endl;
	cout << "There are "<< N << " equivalence classes Lambda[j]"<<endl;
    cout<< "The equivalence classes Lambda[j] for poles are represented by"<<endl;
    for (j=1;j<=N;j++)
    {
    	cout << "lambda_k["<<j<<"] = "<< lambda_k[j] << " with order l["<<j<<"] = "<<l[j]<<endl;
    }
    cout<<endl;
    
    // compute the values m[j] for the respective lambda_k[j]
    
    Complex m[N+1];
    for (j=1;j<=N;j++)
    	m[j] = -2*lambda_k[j] + 2;
	
	
	// compute sum_{k=1}^a log Gamma((s+m[j]+2*lambda[k])/2) for each j
	
	int n,fact_n;
	Complex log_Gamma[N+1][a+1][MYDIGITS+1];
	Complex sum_log_Gamma[N+1][MYDIGITS+1];
	
	for (j=1;j<=N;j++)
	for (n=0;n<=MYDIGITS;n++)
			sum_log_Gamma[j][n] = 0;
	
	for (j=1;j<=N;j++)
	{
		for (k=1;k<=a;k++)
		{
			fact_n = 1;
			for (n=0;n<=MYDIGITS;n++)
			{
				if (n!=0)
					fact_n = fact_n*n;	
			
				log_Gamma[j][k][n] = pow(Double(1)/2,n)*(log_GAMMA((m[j]/2 + lambda[k]),n))/fact_n;
				sum_log_Gamma[j][n] = sum_log_Gamma[j][n] + log_Gamma[j][k][n];
			}
		}
	}	
	cout << endl;
	
	// compute the exponential taylor series for gamma = exp(sum_log_Gamma)
	
	Complex exp_sum_log_Gamma[N+1][MYDIGITS+1][MYDIGITS+1]; // symmetric functions
	Complex gamma[N+1][MYDIGITS+1]; // gamma(s+m[j]) for j = 1 to N
	Complex temp_gamma[MYDIGITS+1];
	Complex temp_mult_gamma[MYDIGITS+1];
	Complex temp_exp_sum_log_Gamma[MYDIGITS+1];
	int fact_n_k;
		
	for (j=1;j<=N;j++)
	{
		for (n=0;n<=MYDIGITS;n++)
			exp_sum_log_Gamma[j][0][n] = 0;			
		exp_sum_log_Gamma[j][0][0] = exp(sum_log_Gamma[j][0]);
		
		for (k=1;k<=MYDIGITS;k++)
		{
			fact_n_k = 1;
						
			for (n=0;n<=MYDIGITS;n++)
			{
				if(n%k == 0)
				{
					if (n/k != 0)
						fact_n_k = fact_n_k*(n/k);
						
					exp_sum_log_Gamma[j][k][n] = pow(sum_log_Gamma[j][k],n/k)/fact_n_k;
				}
				else
					exp_sum_log_Gamma[j][k][n] = 0;
			}
		}		
	}
	
	for (j=1;j<=N;j++)
	{
		
		for (n=0;n<=MYDIGITS;n++)
			temp_mult_gamma[n] = exp_sum_log_Gamma[j][0][n];
						
		for (k=1;k<=MYDIGITS;k++)
		{
			for (n=0;n<=MYDIGITS;n++)
			{
				temp_exp_sum_log_Gamma[n] = exp_sum_log_Gamma[j][k][n];
				temp_gamma[n] = temp_mult_gamma[n];
			}
				
			mult_poly_taylor(temp_gamma,temp_exp_sum_log_Gamma,temp_mult_gamma,MYDIGITS);
		}
		
		for (n=0;n<=MYDIGITS;n++)
			gamma[j][n] = temp_mult_gamma[n];		
					
	}
	
	cout << "-----------------------------------------------"<< endl;
    cout << "The gamma(s+m[j]) coefficients are as follows"<< endl<<endl;
    
    for(j=1;j<=N;j++)
    {
		cout<<"gamma(s+m["<<j<<"]) = "<<"gamma(s+"<<m[j]<<") = "<<gamma[j][0];
    	for (n=1;n<=MYDIGITS;n++)
    		cout << " + " << gamma[j][n] <<" s^"<<n;
    		
    	cout<<endl;
    }
    
    cout << "-----------------------------------------------"<< endl;
    
}

#endif
