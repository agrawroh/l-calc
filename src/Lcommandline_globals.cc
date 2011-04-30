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

#include "Lcommandline_globals.h"

//-----Global variables--------------------------------

int current_L_type;     //1,2,3:For the 3 diff types of L-functions: int, Double, Complex

L_function<int> int_L,int_L2,int_L3;
L_function<Double> Double_L,Double_L2,Double_L3;
L_function<Complex> Complex_L,Complex_L2,Complex_L3;

// below is needed since default initialization (i.e. default constructor)
// includes things like Pi. Pi is computed
// in Lgobals.cc: Linitialize_globals() according to required precision.
// So we need to reintialize after first calling Linitialize_globals
// in Lcommandline.cc
void initialize_commandline_globals(){

    L_function<int> tmp_int_L,tmp_int_L2;
    L_function<Double> tmp_Double_L,tmp_Double_L2;
    L_function<Complex> tmp_Complex_L,tmp_Complex_L2;


    int_L = tmp_int_L;
    int_L2 = tmp_int_L2;

    Double_L = tmp_Double_L;
    Double_L2 = tmp_Double_L2;

    Complex_L = tmp_Complex_L;
    Complex_L2 = tmp_Complex_L2;



}


