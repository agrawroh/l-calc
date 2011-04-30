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


#ifndef Lcommandline_globals_H
#define Lcommandline_globals_H

#include "L.h"

//-----Global variables--------------------------------

extern int current_L_type;     //1,2,3:For the 3 diff types of L-functions: int, Double, Complex

extern L_function<int> int_L,int_L2,int_L3;
extern L_function<Double> Double_L,Double_L2,Double_L3;
extern L_function<Complex> Complex_L,Complex_L2,Complex_L3;

void initialize_commandline_globals();

#endif
