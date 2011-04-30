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

/*
  Michael Rubinstein's addition to mpreal.h
*/

#ifndef Lmpreal_H
#define Lmpreal_H

#include "mpreal.h"

#define cpplongop(op) \
inline mpreal operator op(const mpreal& x, const long long& y) { return x op (double)y; } \
inline mpreal operator op(const long long& y, const mpreal& x) { return (double)y op x; }

cpplongop(+) cpplongop(*) cpplongop(-) cpplongop(/)

template<class T> inline void reset(T& t) {
    mpfr_clear(t.get_mpfr_t());
    mpfr_init(t.get_mpfr_t());
}

#define cppop(op) \
template<typename T> inline mpreal op(const T& x) {\
        mpreal ret,y=x;\
        mpfr_##op(ret.get_mpfr_t(), y.get_mpfr_t(), __gmp_default_rounding_mode);\
        return ret;\
}
#define cpp_two_op(op) \
template<typename T> inline mpreal op(const T& x) {\
        mpreal ret,y=x;\
        mpfr_##op(ret.get_mpfr_t(), y.get_mpfr_t());\
        return ret;\
}

cppop(abs)
cppop(rint)
cpp_two_op(trunc)
cpp_two_op(floor)
cpp_two_op(ceil)
cppop(sqrt)
cppop(log)
cppop(log2)
cppop(log10)
cppop(exp)
cppop(exp2)
cppop(exp10)
cppop(cos)
cppop(sin)
cppop(tan)
cppop(sec)
cppop(csc)
cppop(cot)
cppop(acos)
cppop(asin)
cppop(atan)
cppop(cosh)
cppop(sinh)
cppop(tanh)
cppop(sech)
cppop(csch)
cppop(coth)
cppop(acosh)
cppop(asinh)
cppop(atanh)
cppop(log1p)
cppop(expm1)
cppop(eint)
cppop(gamma)
cppop(lngamma)
cppop(zeta)
cppop(erf)
cppop(erfc)

template<typename T, typename U> inline mpreal pow(const T& a, const U& b) {
    return exp(log(a)*b);
}

template<typename T> inline mpreal atan2(const T& y, const T& x) {
    mpreal ret,Y=y,X=x;\
    mpfr_atan2(ret.get_mpfr_t(), Y.get_mpfr_t(), X.get_mpfr_t(),__gmp_default_rounding_mode);\
    return ret;\
}

/*
__GMP_DEFINE_UNARY_FUNCTION(abs, __gmp_abs_function)
__GMP_DEFINE_UNARY_FUNCTION(rint, __gmp_rint_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(trunc, __gmp_trunc_function)
__GMP_DEFINE_UNARY_FUNCTION(floor, __gmp_floor_function)
__GMP_DEFINE_UNARY_FUNCTION(ceil, __gmp_ceil_function)
__GMP_DEFINE_UNARY_FUNCTION(sqrt, __gmp_sqrt_function)
__GMP_DEFINE_UNARY_FUNCTION(log, __gmp_log_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(log2, __gmp_log2_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(log10, __gmp_log10_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(exp, __gmp_exp_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(exp2, __gmp_exp2_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(exp10, __gmp_exp10_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(cos, __gmp_cos_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(sin, __gmp_sin_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(tan, __gmp_tan_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(sec, __gmp_sec_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(csc, __gmp_csc_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(cot, __gmp_cot_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(acos, __gmp_acos_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(asin, __gmp_asin_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(atan, __gmp_atan_function) // not in gmpxx.h
__GMP_DEFINE_BINARY_FUNCTION(atan2, __gmp_atan2_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(cosh, __gmp_cosh_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(sinh, __gmp_sinh_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(tanh, __gmp_tanh_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(sech, __gmp_sech_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(csch, __gmp_csch_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(coth, __gmp_coth_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(acosh, __gmp_acosh_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(asinh, __gmp_asinh_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(atanh, __gmp_atanh_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION_UI(mpfr_t, fac_ui, __gmp_fac_ui_function) // not gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(log1p, __gmp_log1p_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(expm1, __gmp_expm1_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(eint, __gmp_eint_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(gamma, __gmp_gamma_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(lngamma, __gmp_lngamma_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(zeta, __gmp_zeta_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(erf, __gmp_erf_function) // not in gmpxx.h
__GMP_DEFINE_UNARY_FUNCTION(erfc, __gmp_erfc_function) // not in gmpxx.h
*/

//inline const mpreal fmod (const mpreal& x, const mpreal& y, mp_rnd_t rnd_mode)
inline const mpreal fmod (const mpreal& x, const mpreal& y)
{
    mpreal a;
    mp_prec_t yp, xp;

    yp = y.get_prec();
    xp = x.get_prec();

    a.set_prec(yp>xp?yp:xp);

    //mpfr_fmod(a.mp, x.mp, y.mp, rnd_mode);
    //mpfr_fmod(a.mp, x.mp, y.mp, __gmp_default_rounding_mode);
    mpfr_fmod(a.get_mpfr_t(), x.get_mpfr_t(), y.get_mpfr_t(), __gmp_default_rounding_mode);

return a;
}


#endif
