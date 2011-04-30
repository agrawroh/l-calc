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


#ifndef Lcommandline_twist_H
#define Lcommandline_twist_H

#include "L.h"             //for the L_function class
//#include "Lcommandline_numbertheory.h"  //for number theory functions
#include "Lcommandline_globals.h"      //command line global variables
//#include "Lelliptic.h"

//===============functions================================================

//int quadratic_twists(Long D1, Long D2,Double x,Double y,Long count,Long start_N,Double step_size, const char *what_to_do,bool do_only_even_twists,bool test_explicit_formula,int desired_rank=-1);
int quadratic_twists(Long D1, Long D2,Double x,Double y,Long count,Long start_N,Double step_size, const char *what_to_do,bool do_only_even_twists,int desired_rank=-1);
//int all_twists(Long D1, Long D2,Double x,Double y,Long count,Long start_N,Double step_size, const char *what_to_do,int twist_type,int print_character,bool test_explicit_formula);
int all_twists(Long D1, Long D2,Double x,Double y,Long count,Long start_N,Double step_size, const char *what_to_do,int twist_type,int print_character);

#endif
