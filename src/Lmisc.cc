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


#include "Lmisc.h"

//inline Double LOG(double n) // Used for long
//{
    //if(n<=number_logs) return LG[Int(n)];
    //else return log(Double(n));
//}


//Double SQRT(double n) // Used for long
//{
    //if(n<=number_sqrts) return SQUARE_ROOT[Int(n)];
    //else return sqrt(Double(n));
//}


//For parsing c++ strings 
vector<string> &split(const string &s, char delim, vector<string> &elems) {
    stringstream ss(s);
    string item;
    while(getline(ss, item, delim)) {
            elems.push_back(item);
    }
    return elems;
}


vector<string> split(const string &s, char delim) {
    vector<string> elems;
    return split(s, delim, elems);
}

vector<Double> &split_Double(const string &s, char delim, vector<Double> &elems) {
    stringstream ss(s);
    string item;
    Double x;
    while(getline(ss, item, delim)) {
        #ifdef USE_LONG_DOUBLE
            sscanf(item.c_str(),"%Lg",&x);
        #elif USE_MPFR
            x=item.c_str();
        #elif USE_MPFRCPP
            x=item.c_str();
        #elif USE_BAILEY_DD
            x=item.c_str();
        #elif USE_BAILEY_QD
            x=item.c_str();
        #else
            sscanf(item.c_str(),"%lg",&x);
        #endif
            elems.push_back(x);
    }
    return elems;
}


vector<Double> split_Double(const string &s, char delim) {
    vector<Double> elems;
    return split_Double(s, delim, elems);
}
