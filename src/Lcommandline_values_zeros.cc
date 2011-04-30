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



#include "L.h"
#include "Lcommandline_values_zeros.h"


//-----functions--------------------------------------------------------------
void compute_values(Double x,Double y,const char *return_type,const char *file_name,Double x3,Double y3,Long count)
{

    fstream file;

    Complex s,u;
    Double x2,y2;

    //cout << setprecision(DIGITS3);

    if(strcmp(file_name,"")) //i.e. if a file_name has been provided...
    {
        file.open(file_name, ios::in);
        if(!file.is_open())
        {
            cout << endl << "FAILED TO FIND THE FILE: " << file_name<< endl;
            exit(1);
        }

        while (!file.eof())
        {
            file >> x2;
            file >> y2;
            s=Complex(x2,y2);
            if(!file.eof())
            switch(current_L_type)
            {
                case 1:
                    u=int_L.value(s,global_derivative,return_type); 
                    //global_derivative specifies which derivative to compute. Default is 0 which
                    //is just the plain L-value
                    cout << real(u) << " " << imag(u) << endl;
                    //cout << imag(s) << " " << abs(u)<< endl;
                    break;
                case 2:
                    u=Double_L.value(s,global_derivative,return_type);
                    cout << real(u) << " " << imag(u) << endl;
                    break;
                case 3:
                    u=Complex_L.value(s,global_derivative,return_type);
                    cout << real(u) << " " << imag(u) << endl;
                    break;
            }
        }
        file.close();
    }

    else  // else use the x and y provided
    {
        Long n=0;
        s=Complex(x,y);
        //#pragma omp sections
        {
        do{
            n++;
            cout << setprecision(DIGITS);
            if(count>0)cout << real(s) << " " << imag(s) <<" ";
            cout << setprecision(DIGITS3);
            switch(current_L_type)
            {
                case 1:
                    u=int_L.value(s,global_derivative,return_type);
                    cout << real(u) << " " << imag(u) << endl;
                    //cout << imag(s) << " " << abs(u)<< endl;
                    break;
                case 2:
                    u=Double_L.value(s,global_derivative,return_type);
                    cout << real(u) << " " << imag(u) << endl;
                    break;
                case 3:
                    u=Complex_L.value(s,global_derivative,return_type);
                    cout << real(u) << " " << imag(u) << endl;
                    break;
            }
            s=s+Complex(x3-x,y3-y)/(double)count;
        } while(n<count);
        }
    }


}


//------------------------------------------------
//  Zero finding code

//looks for zeros on the critical line starting at 1/2+It1 by
//advancing in steps of size step_size (which can be negative)
//and looking for sign changes. If count==0 we look for
//zeros between t1=x and t2=y. If count>0, we look for the first
//count zeros (starting at t1 not yet implemented) using steps of no smaller
//than step_size.
//if count >0 then t2 is not used, so it's value is irrelevant.

void compute_zeros(Double x, Double y,Double step_size, Long count,Long start_N,int rank, bool test_explicit_formula)
{

    if(current_L_type==1){
        if(count==0)
        int_L.find_zeros(x,y,step_size);
        else int_L.find_zeros(count,start_N,step_size,rank);
        //else int_L.find_zeros(count,start_N,step_size,rank,test_explicit_formula);
    }
    if(current_L_type==2){
        if(count==0)
        Double_L.find_zeros(x,y,step_size);
        else Double_L.find_zeros(count,start_N,step_size,rank);//,test_explicit_formula);
    }
    if(current_L_type==3){
        if(count==0)
        Complex_L.find_zeros(x,y,step_size);
        Complex_L.find_zeros(count,start_N,step_size,rank);//,test_explicit_formula);
    }

}


//find the zeros of the function obtained by interpolating L and L2. Assumes
//they are of the same quasi degree and have same number of dirichlet 
//coefficients. The interpolated L-function has all its basic data as t times the
//data of the first L-function and (1-t) times the data of the second L-function.
//does not work with integer L_functions since interpolation gives at least doubles
void L_interpolate(Double x, Double y,Double step_size, int n)
{

    double t;
    char message[300];
    ostringstream os3;

    t=0;

    do{

        os3.str("");

        os3 << t << endl;
        strcpy(message, os3.str().c_str());

        if(current_L_type==2){
            Double_L3=Double_L*(1-t)+Double_L2*t;
            //Double_L3.find_zeros(x,y,step_size,"cout",message_stamp);
            Double_L3.find_zeros(x,y,step_size,"cout",message);
        }
        if(current_L_type==3){
            Complex_L3=Complex_L*(1-t)+Complex_L2*t;
            //Complex_L3.find_zeros(x,y,step_size,"cout",message_stamp);
            Complex_L3.find_zeros(x,y,step_size,"cout",message);
        }
        t=t+1./n;
    }while(t<=1);

}

void compute_rank(){
    switch(current_L_type)
    {
        case 1:
            int_L.compute_rank(true);
            break;
        case 2:
            Double_L.compute_rank(true);
            break;
        case 3:
            Complex_L.compute_rank(true);
            break;
    }
}

void verify_rank(int rank){
    switch(current_L_type)
    {
        case 1:
            int_L.verify_rank(rank);
            break;
        case 2:
            Double_L.verify_rank(rank);
            break;
        case 3:
            Complex_L.verify_rank(rank);
            break;
    }
}


