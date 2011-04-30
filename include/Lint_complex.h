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


#ifndef Lint_complex_H
#define Lint_complex_H

inline Complex operator + (Complex x, int y)
{
  return Complex (real (x) + y, imag (x));
}

inline Complex operator + (int x, Complex y)
{
  return Complex (x + real (y), imag (y));
}

inline Complex operator - (Complex x, int y)
{
  return Complex (real (x) - y, imag (x));
}

inline Complex operator - (int x, Complex y)
{
  return Complex (x - real (y), - imag (y));
}

inline Complex operator * (Complex x, int y)
{
  return Complex (real (x) * y, imag (x) * y);
}

inline Complex operator * (int x, Complex y)
{
  return Complex (x * real (y), x * imag (y));
}

inline Complex operator / (Complex x, int y)
//Complex operator / (Complex x, int y)
{
  return Complex (real (x) / y, imag (x) / y);
}

inline Complex operator / (int x, Complex y)
//Complex operator / (int x, Complex y)
{
 return (((Double)x)/y);
}


inline bool operator == (Complex x, int y)
{
  return real (x) == y && imag (x) == 0;
}

inline bool operator == (int x, Complex y)
{
  return x == real (y) && imag (y) == 0;
}

inline bool operator != (Complex x, int y)
{
  return real (x) != y || imag (x) != 0;
}

inline bool operator != (int x, Complex y)
{
  return x != real (y) || imag (y) != 0;
}

inline Complex operator + (Complex x, Long y)
{
  return Complex (real (x) + y, imag (x));
}

inline Complex operator + (Long x, Complex y)
{
  return Complex (x + real (y), imag (y));
}

inline Complex operator - (Complex x, Long y)
{
  return Complex (real (x) - y, imag (x));
}

inline Complex operator - (Long x, Complex y)
{
  return Complex (x - real (y), - imag (y));
}


inline Complex operator * (Complex x, Long y)
{
  return Complex (real (x) * y, imag (x) * y);
}

inline Complex operator * (Long x, Complex y)
{
  return Complex (x * real (y), x * imag (y));
}

inline Complex operator / (Complex x, Long y)
//Complex operator / (Complex x, Long y)
{
  return Complex (real (x) / y, imag (x) / y);
}

inline Complex operator / (Long x, Complex y)
//Complex operator / (Long x, Complex y)
{
  return (((Double)x)/y);
}


inline bool operator == (Complex x, Long y)
{
  return real (x) == y && imag (x) == 0;
}

inline bool operator == (Long x, Complex y)
{
  return x == real (y) && imag (y) == 0;
}

inline bool operator != (Complex x, Long y)
{
  return real (x) != y || imag (x) != 0;
}

inline bool operator != (Long x, Complex y)
{
  return x != real (y) || imag (y) != 0;
}

#endif
