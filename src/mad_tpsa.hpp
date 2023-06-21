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
#include <cstdio>
#include <cassert>
#include <memory>

// --- debug ------------------------------------------------------------------o

#define TRC(...)

#ifndef TRC
#define TRC(...) \
  (printf("%s:%3d:%12s: ", __FILE__, __LINE__, __func__), \
   printf(__VA_ARGS__), printf("\n"));
#endif

// --- types ------------------------------------------------------------------o

namespace mad {

// private class to build dtor as a type and avoid extra function ptr in tpsa.
struct tpsa_del_ {
  void operator()(tpsa_t *t) { TRC("tpsa_t* %p", (void*)t) mad_tpsa_del(t); }
};

// public class to use tpsa, i.e. wrapper to tpsa_t with memory management.
using tpsa = std::unique_ptr<tpsa_t, tpsa_del_>;

// public class to *locally* wrap tpsa_t without memory management.
struct tpsa_ref {
  explicit tpsa_ref(tpsa_t *a) : t(*a) { TRC("tpsa_t* %p", (void*)a      ) }
  explicit tpsa_ref(tpsa   &a) : t(*a) { TRC("tpsa& %p"  , (void*)a.get()) }
           tpsa_ref(tpsa_t &a) : t( a) { TRC("tpsa_t& %p", (void*)&a     ) }

  tpsa_ref()                     = delete;
  tpsa_ref(tpsa_ref&&)           = delete;
  tpsa_ref(const tpsa_ref&)      = delete;
  tpsa_ref(std::nullptr_t)       = delete;
  void operator=(tpsa_ref&&)     = delete;
  void operator=(std::nullptr_t) = delete;

  tpsa_t* get()       const { return &t; }
  tpsa_t& operator*() const { return  t; }

  void operator= (const tpsa     &a) { TRC("tpa %p", (void*)a.get()) mad_tpsa_copy(&*a,&t); }
  void operator= (const tpsa_ref &a) { TRC("ref %p", (void*)a.get()) mad_tpsa_copy(&*a,&t); }
  void operator= (      num_t     a) { TRC("num")                    mad_tpsa_set0(&t,0,a); }

  void operator+=(const tpsa     &a) { mad_tpsa_add (&t,&*a,&t); }
  void operator+=(const tpsa_ref &a) { mad_tpsa_add (&t,&*a,&t); }
  void operator+=(      num_t     a) { mad_tpsa_set0(&t,  1, a); }

  void operator-=(const tpsa     &a) { mad_tpsa_sub (&t,&*a,&t); }
  void operator-=(const tpsa_ref &a) { mad_tpsa_sub (&t,&*a,&t); }
  void operator-=(      num_t     a) { mad_tpsa_set0(&t,  1,-a); }

  void operator*=(const tpsa     &a) { mad_tpsa_mul (&t,&*a,&t); }
  void operator*=(const tpsa_ref &a) { mad_tpsa_mul (&t,&*a,&t); }
  void operator*=(      num_t     a) { mad_tpsa_scl (&t,  a,&t); }

  void operator/=(const tpsa     &a) { mad_tpsa_div (&t,&*a,&t); }
  void operator/=(const tpsa_ref &a) { mad_tpsa_div (&t,&*a,&t); }
  void operator/=(      num_t     a) { mad_tpsa_scl (&t,1/a,&t); }

private:
  tpsa_t &t;
};

// shortcut
using ref = tpsa_ref;

} // mad

// --- operators --------------------------------------------------------------o

namespace mad {

// --- ctors ---

inline tpsa
newt () {  TRC("void")
  return tpsa(mad_tpsa_newd(mad_desc_curr, mad_tpsa_default));
}

inline tpsa
newt (int mo) {  TRC("int")
  return tpsa(mad_tpsa_newd(mad_desc_curr, mo));
}

inline tpsa
newt (const tpsa_ref &a, int mo=mad_tpsa_default) {  TRC("ref")
  return tpsa(mad_tpsa_new(a.get(), mo));
}

// --- unary ---

inline tpsa
operator- (const tpsa_ref &a) {  TRC("ref")
  tpsa c(newt(a));
  mad_tpsa_scl(a.get(), -1, c.get());
  return c;
}

inline tpsa
operator- (const tpsa &a) {  TRC("tpa")
  return -tpsa_ref(*a);
}

// --- add ---

inline tpsa
operator+ (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa c(newt(a));
  mad_tpsa_add(a.get(), b.get(), c.get());
  return c;
}

inline tpsa
operator+ (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa c(newt(a));
  mad_tpsa_copy(a.get(), c.get());
  mad_tpsa_set0(c.get(), 1, b);
  return c;
}

inline tpsa operator+(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return *a+*b; }
inline tpsa operator+(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return *a+ b; }
inline tpsa operator+(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return  a+*b; }
inline tpsa operator+(const tpsa     &a,       num_t     b) { TRC("tpa,num") return *a+ b; }
inline tpsa operator+(      num_t     a, const tpsa     &b) { TRC("num,tpa") return *b+ a; }
inline tpsa operator+(      num_t     a, const tpsa_ref &b) { TRC("num,ref") return  b+ a; }

// --- sub ---

inline tpsa
operator- (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa c(newt(a));
  mad_tpsa_sub(a.get(), b.get(), c.get());
  return c;
}

inline tpsa
operator- (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa c(newt(a));
  mad_tpsa_copy(a.get(), c.get());
  mad_tpsa_set0(c.get(), 1, -b);
  return c;
}

inline tpsa
operator- (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  tpsa c(newt(b));
  mad_tpsa_scl (b.get(),-1, c.get());
  mad_tpsa_set0(c.get(), 1, a);
  return c;
}

inline tpsa operator-(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return *a-*b; }
inline tpsa operator-(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return *a- b; }
inline tpsa operator-(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return  a-*b; }
inline tpsa operator-(const tpsa     &a,       num_t     b) { TRC("tpa,num") return *a- b; }
inline tpsa operator-(      num_t     a, const tpsa     &b) { TRC("num,tpa") return  a-*b; }

// --- mul ---

inline tpsa
operator* (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa c(newt(a));
  mad_tpsa_mul(a.get(), b.get(), c.get());
  return c;
}

inline tpsa
operator* (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa c(newt(a));
  mad_tpsa_scl(a.get(), b, c.get());
  return c;
}

inline tpsa operator*(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return *a**b; }
inline tpsa operator*(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return *a* b; }
inline tpsa operator*(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return  a**b; }
inline tpsa operator*(const tpsa     &a,       num_t     b) { TRC("tpa,num") return *a* b; }
inline tpsa operator*(      num_t     a, const tpsa     &b) { TRC("num,tpa") return *b* a; }
inline tpsa operator*(      num_t     a, const tpsa_ref &b) { TRC("num,ref") return  b* a; }

// --- div ---

inline tpsa
operator/ (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa c(newt(a));
  mad_tpsa_div(a.get(), b.get(), c.get());
  return c;
}

inline tpsa
operator/ (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa c(newt(a));
  mad_tpsa_scl(a.get(), 1/b, c.get());
  return c;
}

inline tpsa
operator/ (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  tpsa c(newt(b));
  mad_tpsa_inv(b.get(), a, c.get());
  return c;
}

inline tpsa operator/(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return *a/ *b; }
inline tpsa operator/(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return *a/  b; }
inline tpsa operator/(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return  a/ *b; }
inline tpsa operator/(const tpsa     &a,       num_t     b) { TRC("tpa,num") return *a/  b; }
inline tpsa operator/(      num_t     a, const tpsa     &b) { TRC("num,tpa") return  a/ *b; }

// --- pow ---

inline tpsa
pow (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa c(newt(a));
  mad_tpsa_pow(a.get(), b.get(), c.get());
  return c;
}

inline tpsa
pow (const tpsa_ref &a, int b) {  TRC("ref,int")
  tpsa c(newt(a));
  mad_tpsa_powi(a.get(), b, c.get());
  return c;
}

inline tpsa
pow (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa c(newt(a));
  mad_tpsa_pown(a.get(), b, c.get());
  return c;
}

inline tpsa
pow (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  tpsa c(newt(b));
  mad_tpsa_scl(b.get(), std::log(a), c.get());
  mad_tpsa_exp(c.get(), c.get());
  return c;
}

inline tpsa pow(const tpsa     &a, const tpsa     &b) { TRC("tpa,int") return pow(*a,*b); }
inline tpsa pow(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return pow(*a, b); }
inline tpsa pow(const tpsa_ref &a, const tpsa     &b) { TRC("ref,int") return pow( a,*b); }
inline tpsa pow(const tpsa     &a,       int       b) { TRC("tpa,int") return pow(*a, b); }
inline tpsa pow(const tpsa     &a,       num_t     b) { TRC("tpa,num") return pow(*a, b); }
inline tpsa pow(      num_t     a, const tpsa     &b) { TRC("num,tpa") return pow( a,*b); }

// warning: the operator ^ hasn't the expected precedence and associativity...

inline tpsa operator^(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return pow(a,b); }
inline tpsa operator^(const tpsa_ref &a, const tpsa_ref &b) { TRC("ref,ref") return pow(a,b); }
inline tpsa operator^(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return pow(a,b); }
inline tpsa operator^(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return pow(a,b); }
inline tpsa operator^(const tpsa_ref &a,       int       b) { TRC("ref,int") return pow(a,b); }
inline tpsa operator^(const tpsa_ref &a,       num_t     b) { TRC("ref,num") return pow(a,b); }
inline tpsa operator^(const tpsa     &a,       int       b) { TRC("tpa,int") return pow(a,b); }
inline tpsa operator^(const tpsa     &a,       num_t     b) { TRC("tpa,num") return pow(a,b); }
inline tpsa operator^(      num_t     a, const tpsa     &b) { TRC("num,tpa") return pow(a,b); }
inline tpsa operator^(      num_t     a, const tpsa_ref &b) { TRC("num,ref") return pow(a,b); }

} // mad

// --- functions ---

namespace mad {

inline num_t sqr(      num_t      a) { TRC("num") return a*a; }
inline tpsa  sqr(const tpsa      &a) { TRC("tpa") return a*a; }
inline tpsa  sqr(const tpsa_ref  &a) { TRC("ref") return a*a; }

// --- unary ---

#define FUN(F) \
\
inline tpsa F (const tpsa_ref &a) { \
  TRC("ref") \
  tpsa c(newt(a)); \
  mad_tpsa_ ## F (a.get(), c.get()); \
  return c; \
} \
\
inline tpsa F (const tpsa &a) { TRC("tpa") return F(*a); }

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

#undef TRC

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_HPP
