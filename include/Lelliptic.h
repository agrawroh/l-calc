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


#ifndef Lelliptic_H
#define Lelliptic_H


#include "L.h"             //for the L_function class
//#include "Lcommandline_numbertheory.h"  //for number theory functions
//#include "Lcommandline_globals.h"      //command line global variables
#include <iostream>


#ifdef INCLUDE_PARI
#include <pari/pari.h>     //for pari's elliptic curve functions 
#undef init                //pari has a '#define init pari_init' which
                           //causes trouble with the stream.h init.
                           //pari also causes trouble with things like abs.
                           //we place the pari include first since otherwise it
                           //messes up.

void data_E(char *a1, char *a2, char *a3, char *a4, char *a6,int N_terms,Double * coeff);
#endif //ifdef INCLUDE_PARI



#endif
