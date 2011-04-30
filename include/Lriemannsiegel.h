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


//Reiemann Siegel formula.
//http://www.dtc.umn.edu/~odlyzko/doc/arch/fast.zeta.eval.pdf
#ifndef Lriemannsiegel_H
#define Lriemannsiegel_H

//#include "Lglobals.h"           //for global variables
//#include "Lgamma.h"           //for global variables
//#include "Lmisc.h"           //for global variables

Complex Zeta(Complex s, const char *return_type="pure");
Double rs_remainder_terms(Double z, Double tau);
Complex siegel(Complex s);

#endif
