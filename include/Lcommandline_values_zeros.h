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


#ifndef Lcommandline_values_zeros_H
#define Lcommandline_values_zeros_H

#include "L.h"              //for the L_function class
//#include "Lcommandline_numbertheory.h"  //for number theory functions
#include "Lcommandline_globals.h"      //command line global variables


//-----functions--------------------------------------------------------------
void compute_values(Double x,Double y,const char *return_type="pure",const char *file_name="",Double x3=0,Double y3=0,Long count=0);
void compute_zeros(Double x, Double y,Double step_size, Long count=0,Long start_N=0,int rank=-1,bool test_explicit_formula=false);
void L_interpolate(Double x, Double y,Double step_size, int n=1000);

void verify_rank(int rank);
void compute_rank();


#endif

