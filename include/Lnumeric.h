
/*
#define cxxop(op) \
template<typename T> inline mpfr_class op(const T& x) {\
	mpfr_class ret,y=x;\
	mpfr_##op(ret.get_mpfr_t(), y.get_mpfr_t(), __gmp_default_rounding_mode);\
	return ret;\
}

cxxop(log) cxxop(exp) cxxop(cos) cxxop(sin) cxxop(tan)
cxxop(atan) cxxop(cosh) cxxop(sinh) cxxop(rint)

*/

#define cxxlongop(op) \
inline mpfr_class operator op(const mpfr_class& x, const long long& y) { return x op (double)y; } \
inline mpfr_class operator op(const long long& y, const mpfr_class& x) { return (double)y op x; }

cxxlongop(+) cxxlongop(*) cxxlongop(-) cxxlongop(/)

#define cpplongop(op) \
inline mpreal operator op(const mpreal& x, const long long& y) { return x op (double)y; } \
inline mpreal operator op(const long long& y, const mpreal& x) { return (double)y op x; }

cpplongop(+) cpplongop(*) cpplongop(-) cpplongop(/)

/*
template<typename T> inline mpfr_class atan2(const T& y, const T& x) {
	if(x==0)
		if(y>0) return Pi/2;
		else if(y<0) return -Pi/2;
		else return 0;
	else if(x>0) return atan(y/x);
	else
		if(y>0) return atan(y/x)+Pi;
		else return atan(y/x)-Pi;
}

*/


/*
template<typename T, typename U> inline mpfr_class pow(const T& a, const U& b) {
	return exp(log(a)*b);
}
*/


// Called by the init function when the precision is changed
template<class T> inline void reset(T& t) {
    mpfr_clear(t.get_mpfr_t());
    mpfr_init(t.get_mpfr_t());
}
