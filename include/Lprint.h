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


    //-----------prints the basic data of the L-function----------------

    template <class ttype>
    void L_function <ttype>::
    print_data_L (int N)   // prints the basic data for the L-function
    {

        int k;
        cout << "DIGITSi3 SET TO " << DIGITS3 << endl;
        cout << setprecision (DIGITS3);

        cout << "-----------------------------------------------" << endl << endl;
        cout << "Name of L_function: " << name << endl;
        if(number_of_dirichlet_coefficients>0)
        cout << "number of dirichlet coefficients = " << number_of_dirichlet_coefficients << endl;
        else cout << "All coefficients are equal to 1" << endl;
        if (what_type_L==1) cout << "coefficients are periodic" << endl;

        k=0;
        if (number_of_dirichlet_coefficients>0) do
        {
           k++;
           cout << "b[" << k << "] = " << dirichlet_coefficient[k] << endl;
        } while (k<N&&k<number_of_dirichlet_coefficients);
        cout << endl;

        cout << "Q = " << Q << endl;
        cout << "OMEGA = " << OMEGA << endl;
        cout << "a = " << a << " (the quasi degree)" << endl;
        for (k=1;k<=a;k++)
        {
            cout << "gamma[" << k << "] =" << gamma[k] << "    ";
            cout << "lambda[" << k << "] =" << lambda[k] << endl;
        }
        cout << endl << endl;

        cout << "number of poles (of the completed L function) = ";
        cout << number_of_poles << endl;

        for (k=1;k<=number_of_poles;k++)
        {
            cout << "pole[" << k << "] =" << pole[k] << "    ";
            cout << "residue[" << k << "] =" << residue[k] << endl;
        }
        cout << "-----------------------------------------------" << endl << endl;


    }

