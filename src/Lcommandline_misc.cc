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
   GNU General Public License for more details

   Check the License for details. You should have received a copy of it, along
   with the package; see the file 'COPYING'. If not, write to the Free Software
   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA._

*/

#include "Lcommandline_misc.h"

//returns 0 if initialization goes smoothly, 1 if there's an error
void initialize_new_L(char *file_name)
{

    int i,j;

    Double x,y;

    fstream file;

    // basic data for the L-function (see the class L_function for full comments)
    int N_terms;
    int what_type;
    Long Period;
    Double q;
    Complex w;
    int A;
    Double *g;
    Complex *l;
    int n_poles;
    Complex *p;
    Complex *r;

    //int missing_data=0;   // XXXX specifies if any of the basic data is missing
                            // NOT YET IMPLEMENTED... so far I assume
                            // that all data has been specified.


    file.open(file_name, ios::in);
    if(!file.is_open())
    {
        cout << endl << "FAILED TO FIND THE FILE: " << file_name<<endl;
        exit(1);
    }

    file >> current_L_type; //1 for int, 2 for Double, 3 for Complex

    file >> what_type; //0 for unknown 1 for periodic, etc

    file >> N_terms; //number of dirichlet coefficients to use
    file >> Period; if(Period<N_terms&&Period>0) N_terms=Period; //Period should be 0 if not periodic

    file >> A; //quasi degree
    g=new Double[A+1];   //gamma's of the GAMMA factors
    l=new Complex[A+1];  //lambda's of the GAMMA factors
    for(i=1;i<=A;i++){
        file >>x;
        g[i]=x;
        file >>x;  file >>y;
        l[i]=Complex(x,y);
    }

    file >> q; // example: Q=1/sqrt(Pi) for Zeta
    file >> x;  file >> y; w=Complex(x,y);

    file >> n_poles;
    p = new Complex[n_poles+1];
    r = new Complex[n_poles+1];
    for (i=1;i<=n_poles;i++)
    {
        file >>x;  file >>y;
        p[i]=Complex(x,y);

        file >>x;  file >>y;
        r[i]=Complex(x,y);
    }


    int *coeff_int;
    Double *coeff_Double;
    Complex *coeff_Complex;


    //XXXXXXXXXXXXXX make more efficient... so as to assume a_n0's don't necessarily appear
    //in the file
    switch(current_L_type)
    {
        case 1:  //i.e. integer coeffs
            coeff_int = new int[N_terms+1];
            for(i=1;i<=N_terms;i++)
            {
                file >> j;
                coeff_int[i]=j;
            }
            break;
        case 2:  //i.e. real coeffs
            coeff_Double = new Double[N_terms+1];
            for(i=1;i<=N_terms;i++)
            {
                file >> x;
                coeff_Double[i]=x;
            }
            break;
        case 3:  //i.e. complex coeffs
            coeff_Complex = new Complex[N_terms+1];
            for(i=1;i<=N_terms;i++)
            {
                file >> x; file >> y;
                coeff_Complex[i]=Complex(x,y);
            }
            break;
    }
    file.close();


    switch(current_L_type)
    {
        case 1:
            int_L=L_function<int>(file_name,what_type,N_terms,coeff_int,Period,q,w,A,g,l,n_poles,p,r);
            delete [] coeff_int;
            break;
        case 2:
            Double_L=L_function<Double>(file_name,what_type,N_terms,coeff_Double,Period,q,w,A,g,l,n_poles,p,r);
            delete [] coeff_Double;
            break;
        case 3:
            Complex_L=L_function<Complex>(file_name,what_type,N_terms,coeff_Complex,Period,q,w,A,g,l,n_poles,p,r);
            delete [] coeff_Complex;
            break;
    }

    delete [] g;
    delete [] l;
    delete [] p;
    delete [] r;

}

void print_L(int N_terms){
    switch(current_L_type)
    {
        case 1:
            int_L.print_data_L(N_terms);
            break;
        case 2:
            Double_L.print_data_L(N_terms);
            break;
        case 3:
            Complex_L.print_data_L(N_terms);
            break;
    }

}

