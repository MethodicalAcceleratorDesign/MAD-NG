#ifndef MAD_TPSA_HPP
#define MAD_TPSA_HPP

/*
 o-----------------------------------------------------------------------------o
 |
 | Simple C++ wrapper to GTPSA
 |
 | Methodical Accelerator Design - Copyright (c) 2016+
 | Support: http://cern.ch/mad  - mad at cern.ch
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  Purpose:
  - Provide memory management and operator overloading for GTPSA, aka tpsa_t.

 o-----------------------------------------------------------------------------o
 */

extern "C" {
#include "mad_tpsa.h"
}

#include <cmath>
#include <memory>

// --- types ------------------------------------------------------------------o

namespace mad {

// private class to build dtor as a type and avoid extra function ptr in tpsa.
struct tpsa_del_ {
  void operator()(tpsa_t *t) { mad_tpsa_del(t); }
};

// public class to use tpsa, i.e. wrapper to tpsa_t.
using tpsa = std::unique_ptr<tpsa_t, tpsa_del_>;

// private class to manage temporaries in expressions, i.e. avoid allocation.
struct tpsa_tmp_ : tpsa {
  explicit
  tpsa_tmp_(tpsa_t    * a) : tpsa(a)       {}
  tpsa_tmp_(tpsa_tmp_&& a) : tpsa(std::move(a)) {}
};

// public class to locally wrap tpsa_t and provide assignments operators
struct ref {
  ref(tpsa_t *a) : t(*a) {}
  ref(tpsa_t &a) : t( a) {}

  ref()                 = delete;
  ref(ref&&)            = delete;
  ref(const ref&)       = delete;
  void operator=(ref&&) = delete;

  void operator =(const tpsa   &a) { mad_tpsa_copy(   &*a,&t); }
  void operator =(const tpsa_t &a) { mad_tpsa_copy(   & a,&t); }
  void operator =(const tpsa_t *a) { mad_tpsa_copy(     a,&t); }

  void operator+=(const tpsa   &a) { mad_tpsa_add (&t,&*a,&t); }
  void operator+=(const tpsa_t &a) { mad_tpsa_add (&t,& a,&t); }
  void operator+=(const tpsa_t *a) { mad_tpsa_add (&t,  a,&t); }
  void operator+=(       num_t  a) { mad_tpsa_set0(&t,  1, a); }

  void operator-=(const tpsa   &a) { mad_tpsa_sub (&t,&*a,&t); }
  void operator-=(const tpsa_t &a) { mad_tpsa_sub (&t,& a,&t); }
  void operator-=(const tpsa_t *a) { mad_tpsa_sub (&t,  a,&t); }
  void operator-=(       num_t  a) { mad_tpsa_set0(&t,  1,-a); }

  void operator*=(const tpsa   &a) { mad_tpsa_mul (&t,&*a,&t); }
  void operator*=(const tpsa_t &a) { mad_tpsa_mul (&t,& a,&t); }
  void operator*=(const tpsa_t *a) { mad_tpsa_mul (&t,  a,&t); }
  void operator*=(       num_t  a) { mad_tpsa_scl (&t,  a,&t); }

  void operator/=(const tpsa   &a) { mad_tpsa_div (&t,&*a,&t); }
  void operator/=(const tpsa_t &a) { mad_tpsa_div (&t,& a,&t); }
  void operator/=(const tpsa_t *a) { mad_tpsa_div (&t,  a,&t); }
  void operator/=(       num_t  a) { mad_tpsa_scl (&t,1/a,&t); }

private:
  tpsa_t &t;
};

} // mad

// --- operators --------------------------------------------------------------o

namespace mad {

// --- ctors ---

inline tpsa_tmp_
newt ()
{
  return tpsa_tmp_(mad_tpsa_newd(mad_desc_curr, mad_tpsa_default));
}

inline tpsa_tmp_
newt (const tpsa_t &a)
{
  return tpsa_tmp_(mad_tpsa_new(&a, mad_tpsa_default));
}

// --- unary ---

inline tpsa_tmp_
operator- (const tpsa_t &a) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_scl(&a, -1, c.get());
  return c;
}

inline tpsa_tmp_ operator-(const tpsa &a) { return -*a; }

// --- add ---

inline tpsa_tmp_
operator+ (const tpsa_t &a, const tpsa_t &b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_add(&a, &b, c.get());
  return c;
}

inline tpsa_tmp_
operator+ (const tpsa_t &a, num_t b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_copy(&a, c.get());
  mad_tpsa_set0(c.get(), 1, b);
  return c;
}

inline tpsa_tmp_
operator+ (tpsa_tmp_ &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_add(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator+ (tpsa_tmp_ &a, const tpsa_t &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_add(c.get(), &b, c.get());
  return c;
}

inline tpsa_tmp_
operator+ (const tpsa_t &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(b.release());
  mad_tpsa_add(&a, c.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator+ (tpsa_tmp_ &a, num_t b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_set0(c.get(), 1, b);
  return c;
}

inline tpsa_tmp_ operator+(const tpsa   &a, const tpsa   &b) { return *a + *b; }
inline tpsa_tmp_ operator+(const tpsa   &a, const tpsa_t &b) { return *a +  b; }
inline tpsa_tmp_ operator+(const tpsa   &a, const tpsa_t *b) { return *a + *b; }
inline tpsa_tmp_ operator+(const tpsa_t &a, const tpsa   &b) { return  a + *b; }
inline tpsa_tmp_ operator+(const tpsa_t *a, const tpsa   &b) { return *a + *b; }
inline tpsa_tmp_ operator+(const tpsa_t &a, const tpsa_t *b) { return  a + *b; }
inline tpsa_tmp_ operator+(const tpsa_t *a, const tpsa_t &b) { return *a +  b; }
inline tpsa_tmp_ operator+(const tpsa   &a,        num_t  b) { return *a +  b; }
inline tpsa_tmp_ operator+(       num_t  a, const tpsa   &b) { return *b +  a; }
inline tpsa_tmp_ operator+(       num_t  a, const tpsa_t &b) { return  b +  a; }

inline tpsa_tmp_ operator+(   tpsa_tmp_ &a, const tpsa_t *b) { return  a + *b; }
inline tpsa_tmp_ operator+(const tpsa_t *a,    tpsa_tmp_ &b) { return *a +  b; }
inline tpsa_tmp_ operator+(       num_t  a,    tpsa_tmp_ &b) { return  b +  a; }

// --- sub ---

inline tpsa_tmp_
operator- (const tpsa_t &a, const tpsa_t &b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_sub(&a, &b, c.get());
  return c;
}

inline tpsa_tmp_
operator- (const tpsa_t &a, num_t b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_copy(&a, c.get());
  mad_tpsa_set0(c.get(), 1, -b);
  return c;
}

inline tpsa_tmp_
operator- (num_t a, const tpsa_t &b) {
  tpsa_tmp_ c(newt(b));
  mad_tpsa_scl(&b, -1, c.get());
  mad_tpsa_set0(c.get(), 1, a);
  return c;
}

inline tpsa_tmp_
operator- (tpsa_tmp_ &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_sub(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator- (tpsa_tmp_ &a, const tpsa_t &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_sub(c.get(), &b, c.get());
  return c;
}

inline tpsa_tmp_
operator- (const tpsa_t &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(b.release());
  mad_tpsa_sub(&a, c.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator- (tpsa_tmp_ &a, num_t b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_set0(c.get(), -1, b);
  return c;
}

inline tpsa_tmp_
operator- (num_t a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(b.release());
  mad_tpsa_scl (c.get(),-1, c.get());
  mad_tpsa_set0(c.get(), 1, a);
  return c;
}

inline tpsa_tmp_ operator-(const tpsa   &a, const tpsa   &b) { return *a - *b; }
inline tpsa_tmp_ operator-(const tpsa   &a, const tpsa_t &b) { return *a -  b; }
inline tpsa_tmp_ operator-(const tpsa   &a, const tpsa_t *b) { return *a - *b; }
inline tpsa_tmp_ operator-(const tpsa_t &a, const tpsa   &b) { return  a - *b; }
inline tpsa_tmp_ operator-(const tpsa_t *a, const tpsa   &b) { return *a - *b; }
inline tpsa_tmp_ operator-(const tpsa_t &a, const tpsa_t *b) { return  a - *b; }
inline tpsa_tmp_ operator-(const tpsa_t *a, const tpsa_t &b) { return *a -  b; }
inline tpsa_tmp_ operator-(const tpsa   &a,        num_t  b) { return *a -  b; }
inline tpsa_tmp_ operator-(       num_t  a, const tpsa   &b) { return  a - *b; }

inline tpsa_tmp_ operator-(   tpsa_tmp_ &a, const tpsa_t *b) { return  a - *b; }
inline tpsa_tmp_ operator-(const tpsa_t *a,    tpsa_tmp_ &b) { return *a -  b; }

// --- mul ---

inline tpsa_tmp_
operator* (const tpsa_t &a, const tpsa_t &b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_mul(&a, &b, c.get());
  return c;
}

inline tpsa_tmp_
operator* (const tpsa_t &a, num_t b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_scl(&a, b, c.get());
  return c;
}

inline tpsa_tmp_
operator* (tpsa_tmp_ &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_mul(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator* (tpsa_tmp_ &a, const tpsa_t &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_mul(c.get(), &b, c.get());
  return c;
}

inline tpsa_tmp_
operator* (const tpsa_t &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(b.release());
  mad_tpsa_mul(&a, c.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator* (tpsa_tmp_ &a, num_t b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_scl(c.get(), b, c.get());
  return c;
}

inline tpsa_tmp_ operator*(const tpsa   &a, const tpsa   &b) { return *a * *b; }
inline tpsa_tmp_ operator*(const tpsa   &a, const tpsa_t &b) { return *a *  b; }
inline tpsa_tmp_ operator*(const tpsa   &a, const tpsa_t *b) { return *a * *b; }
inline tpsa_tmp_ operator*(const tpsa_t &a, const tpsa   &b) { return  a * *b; }
inline tpsa_tmp_ operator*(const tpsa_t *a, const tpsa   &b) { return *a * *b; }
inline tpsa_tmp_ operator*(const tpsa_t &a, const tpsa_t *b) { return  a * *b; }
inline tpsa_tmp_ operator*(const tpsa_t *a, const tpsa_t &b) { return *a *  b; }
inline tpsa_tmp_ operator*(const tpsa   &a,        num_t  b) { return *a *  b; }
inline tpsa_tmp_ operator*(       num_t  a, const tpsa   &b) { return *b *  a; }
inline tpsa_tmp_ operator*(       num_t  a, const tpsa_t &b) { return  b *  a; }

inline tpsa_tmp_ operator*(   tpsa_tmp_ &a, const tpsa_t *b) { return  a * *b; }
inline tpsa_tmp_ operator*(const tpsa_t *a,    tpsa_tmp_ &b) { return *a *  b; }
inline tpsa_tmp_ operator*(       num_t  a,    tpsa_tmp_ &b) { return  b *  a; }

// --- div ---

inline tpsa_tmp_
operator/ (const tpsa_t &a, const tpsa_t &b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_div(&a, &b, c.get());
  return c;
}

inline tpsa_tmp_
operator/ (const tpsa_t &a, num_t b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_scl(&a, 1/b, c.get());
  return c;
}

inline tpsa_tmp_
operator/ (num_t a, const tpsa_t &b) {
  tpsa_tmp_ c(newt(b));
  mad_tpsa_inv(&b, a, c.get());
  return c;
}

inline tpsa_tmp_
operator/ (tpsa_tmp_ &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_div(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator/ (tpsa_tmp_ &a, const tpsa_t &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_div(c.get(), &b, c.get());
  return c;
}

inline tpsa_tmp_
operator/ (const tpsa_t &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(b.release());
  mad_tpsa_div(&a, c.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator/ (tpsa_tmp_ &a, num_t b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_scl(c.get(), 1/b, c.get());
  return c;
}

inline tpsa_tmp_
operator/ (num_t a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(b.release());
  mad_tpsa_inv(c.get(), a, c.get());
  return c;
}

inline tpsa_tmp_ operator/(const tpsa   &a, const tpsa   &b) { return *a / *b; }
inline tpsa_tmp_ operator/(const tpsa   &a, const tpsa_t &b) { return *a /  b; }
inline tpsa_tmp_ operator/(const tpsa   &a, const tpsa_t *b) { return *a / *b; }
inline tpsa_tmp_ operator/(const tpsa_t &a, const tpsa   &b) { return  a / *b; }
inline tpsa_tmp_ operator/(const tpsa_t *a, const tpsa   &b) { return *a / *b; }
inline tpsa_tmp_ operator/(const tpsa_t &a, const tpsa_t *b) { return  a / *b; }
inline tpsa_tmp_ operator/(const tpsa_t *a, const tpsa_t &b) { return *a /  b; }
inline tpsa_tmp_ operator/(const tpsa   &a,        num_t  b) { return *a /  b; }
inline tpsa_tmp_ operator/(       num_t  a, const tpsa   &b) { return  a / *b; }

inline tpsa_tmp_ operator/(   tpsa_tmp_ &a, const tpsa_t *b) { return  a / *b; }
inline tpsa_tmp_ operator/(const tpsa_t *a,    tpsa_tmp_ &b) { return *a /  b; }

// --- pow ---

inline tpsa_tmp_
pow (const tpsa_t &a, const tpsa_t &b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_pow(&a, &b, c.get());
  return c;
}

inline tpsa_tmp_
pow (const tpsa_t &a, int b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_powi(&a, b, c.get());
  return c;
}

inline tpsa_tmp_
pow (const tpsa_t &a, num_t b) {
  tpsa_tmp_ c(newt(a));
  mad_tpsa_pown(&a, b, c.get());
  return c;
}

inline tpsa_tmp_
pow (num_t a, const tpsa_t &b) {
  tpsa_tmp_ c(newt(b));
  mad_tpsa_scl(&b, std::log(a), c.get());
  mad_tpsa_exp(c.get(), c.get());
  return c;
}

inline tpsa_tmp_
pow (tpsa_tmp_ &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_pow(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
pow (tpsa_tmp_ &a, const tpsa_t &b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_pow(c.get(), &b, c.get());
  return c;
}

inline tpsa_tmp_
pow (const tpsa_t &a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(b.release());
  mad_tpsa_pow(&a, c.get(), c.get());
  return c;
}

inline tpsa_tmp_
pow (tpsa_tmp_ &a, int b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_powi(c.get(), b, c.get());
  return c;
}

inline tpsa_tmp_
pow (tpsa_tmp_ &a, num_t b) {
  tpsa_tmp_ c(a.release());
  mad_tpsa_pown(c.get(), b, c.get());
  return c;
}

inline tpsa_tmp_
pow (num_t a, tpsa_tmp_ &b) {
  tpsa_tmp_ c(b.release());
  mad_tpsa_scl(c.get(), std::log(a), c.get());
  mad_tpsa_exp(c.get(), c.get());
  return c;
}

inline tpsa_tmp_ pow(const tpsa   &a, const tpsa   &b) { return pow(*a,*b); }
inline tpsa_tmp_ pow(const tpsa_t &a, const tpsa_t *b) { return pow( a,*b); }
inline tpsa_tmp_ pow(const tpsa_t *a, const tpsa_t &b) { return pow(*a, b); }
inline tpsa_tmp_ pow(const tpsa   &a, const tpsa_t &b) { return pow(*a, b); }
inline tpsa_tmp_ pow(const tpsa   &a, const tpsa_t *b) { return pow(*a,*b); }
inline tpsa_tmp_ pow(const tpsa_t &a, const tpsa   &b) { return pow( a,*b); }
inline tpsa_tmp_ pow(const tpsa_t *a, const tpsa   &b) { return pow(*a,*b); }
inline tpsa_tmp_ pow(const tpsa   &a,        int    b) { return pow(*a, b); }
inline tpsa_tmp_ pow(const tpsa   &a,        num_t  b) { return pow(*a, b); }
inline tpsa_tmp_ pow(       num_t  a, const tpsa   &b) { return pow( a,*b); }

inline tpsa_tmp_ pow(   tpsa_tmp_ &a, const tpsa_t *b) { return pow( a,*b); }
inline tpsa_tmp_ pow(const tpsa_t *a,    tpsa_tmp_ &b) { return pow(*a, b); }

// warning: the operator ^ hasn't the expected precedence and associativity...

inline tpsa_tmp_ operator^(const tpsa   &a, const tpsa   &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t &a, const tpsa_t &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t &a, const tpsa_t *b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t *a, const tpsa_t &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa   &a, const tpsa_t &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa   &a, const tpsa_t *b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t &a, const tpsa   &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t *a, const tpsa   &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t &a,        int    b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t &a,        num_t  b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa   &a,        int    b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa   &a,        num_t  b) { return pow(a,b); }
inline tpsa_tmp_ operator^(       num_t  a, const tpsa   &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(       num_t  a, const tpsa_t &b) { return pow(a,b); }

inline tpsa_tmp_ operator^(   tpsa_tmp_ &a,    tpsa_tmp_ &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(   tpsa_tmp_ &a, const tpsa_t &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t &a,    tpsa_tmp_ &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(   tpsa_tmp_ &a,        int    b) { return pow(a,b); }
inline tpsa_tmp_ operator^(   tpsa_tmp_ &a,        num_t  b) { return pow(a,b); }
inline tpsa_tmp_ operator^(       num_t  a,    tpsa_tmp_ &b) { return pow(a,b); }
inline tpsa_tmp_ operator^(   tpsa_tmp_ &a, const tpsa_t *b) { return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_t *a,    tpsa_tmp_ &b) { return pow(a,b); }

} // mad

// --- functions ---

namespace mad {

// --- unary ---

#define FUN(F) \
\
inline tpsa_tmp_ F (const tpsa_t &a) { \
  tpsa_tmp_ c(newt(a)); \
  mad_tpsa_ ## F (&a, c.get()); \
  return c; \
} \
\
inline tpsa_tmp_ F (tpsa_tmp_ &a) { \
  tpsa_tmp_ c(a.release()); \
  mad_tpsa_ ## F (c.get(), c.get()); \
  return c; \
} \
inline tpsa_tmp_ F (const tpsa &a) { return F(*a); }

FUN(abs   );
FUN(sqrt  );
FUN(exp   );
FUN(log   );
FUN(sin   );
FUN(cos   );
FUN(tan   );
FUN(cot   );
FUN(sinc  );
FUN(sinh  );
FUN(cosh  );
FUN(tanh  );
FUN(coth  );
FUN(sinhc );
FUN(asin  );
FUN(acos  );
FUN(atan  );
FUN(acot  );
FUN(asinc );
FUN(asinh );
FUN(acosh );
FUN(atanh );
FUN(acoth );
FUN(asinhc);
FUN(erf   );
FUN(erfc  );

#undef FUN

} // mad

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_HPP
