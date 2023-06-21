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
#include <memory>

// --- debug ------------------------------------------------------------------o

#if 0
#define TRC(...) \
  (printf("%s:%3d:%12s: ", __FILE__, __LINE__, __func__), \
   printf(__VA_ARGS__), printf("\n"));
#else
#define TRC(...)
#endif

// --- types ------------------------------------------------------------------o

namespace mad {

// private class to build dtor as a type and avoid extra function ptr in tpsa.
struct tpsa_del_ {
  void operator()(tpsa_t *t) { TRC("tpsa_t* %p", (void*)t) mad_tpsa_del(t); }
};

// public class to use tpsa, i.e. wrapper to tpsa_t with memory management.
using tpsa = std::unique_ptr<tpsa_t, tpsa_del_>;

// private class to manage temporaries in expressions, i.e. save allocations.
struct tpsa_tmp_ : tpsa {
  explicit
  tpsa_tmp_(tpsa_t     *a) : tpsa(a)            { TRC("tpsa_t* %p",    (void*)a      ) }
  tpsa_tmp_(tpsa_tmp_ &&a) : tpsa(std::move(a)) { TRC("tpsa_tmp&& %p", (void*)a.get()) }
};

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

inline tpsa_tmp_
newt () {  TRC("void")
  return tpsa_tmp_(mad_tpsa_newd(mad_desc_curr, mad_tpsa_default));
}

inline tpsa_tmp_
newt (int mo) {  TRC("int")
  return tpsa_tmp_(mad_tpsa_newd(mad_desc_curr, mo));
}

inline tpsa_tmp_
newt (const tpsa_ref &a, int mo=mad_tpsa_default) {  TRC("ref")
  return tpsa_tmp_(mad_tpsa_new(a.get(), mo));
}

// --- temps ---

inline tpsa_tmp_&
cct (const tpsa_tmp_ &a) {
  return const_cast<tpsa_tmp_&>(a);
}

// --- unary ---

inline tpsa_tmp_
operator- (const tpsa_ref &a) {  TRC("ref")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_scl(a.get(), -1, c.get());
  return c;
}

inline tpsa_tmp_
operator- (const tpsa_tmp_ &a) {  TRC("tmp")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_scl(c.get(), -1, c.get());
  return c;
}

inline tpsa_tmp_
operator- (const tpsa &a) {  TRC("tpa")
  return -*a;
}

// --- add ---

inline tpsa_tmp_
operator+ (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_add(a.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator+ (const tpsa_tmp_ &a, const tpsa_tmp_ &b) {  TRC("tmp,tmp")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_add(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator+ (const tpsa_tmp_ &a, const tpsa_ref &b) {  TRC("tmp,ref")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_add(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator+ (const tpsa_ref &a, const tpsa_tmp_ &b) {  TRC("ref,tmp")
  tpsa_tmp_ c(cct(b).release());
  mad_tpsa_add(a.get(), c.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator+ (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_copy(a.get(), c.get());
  mad_tpsa_set0(c.get(), 1, b);
  return c;
}

inline tpsa_tmp_
operator+ (const tpsa_tmp_ &a, num_t b) {  TRC("tmp,num")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_set0(c.get(), 1, b);
  return c;
}

inline tpsa_tmp_ operator+(const tpsa      &a, const tpsa      &b) { TRC("tpa,tpa") return *a+*b; }
inline tpsa_tmp_ operator+(const tpsa      &a, const tpsa_ref  &b) { TRC("tpa,ref") return *a+ b; }
inline tpsa_tmp_ operator+(const tpsa      &a, const tpsa_tmp_ &b) { TRC("tpa,tmp") return *a+ b; }
inline tpsa_tmp_ operator+(const tpsa_ref  &a, const tpsa      &b) { TRC("ref,tpa") return  a+*b; }
inline tpsa_tmp_ operator+(const tpsa_tmp_ &a, const tpsa      &b) { TRC("tmp,tpa") return  a+*b; }
inline tpsa_tmp_ operator+(const tpsa      &a,       num_t      b) { TRC("tpa,num") return *a+ b; }
inline tpsa_tmp_ operator+(      num_t      a, const tpsa      &b) { TRC("num,tpa") return *b+ a; }
inline tpsa_tmp_ operator+(      num_t      a, const tpsa_ref  &b) { TRC("num,ref") return  b+ a; }
inline tpsa_tmp_ operator+(      num_t      a, const tpsa_tmp_ &b) { TRC("num,tmp") return  b+ a; }

// --- sub ---

inline tpsa_tmp_
operator- (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_sub(a.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator- (const tpsa_tmp_ &a, const tpsa_tmp_ &b) {  TRC("tmp,tmp")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_sub(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator- (const tpsa_tmp_ &a, const tpsa_ref &b) {  TRC("tmp,ref")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_sub(c.get(), b.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator- (const tpsa_ref &a, const tpsa_tmp_ &b) {  TRC("ref,tmp")
  tpsa_tmp_ c(cct(b).release());
  mad_tpsa_sub(a.get(), c.get(), c.get());
  return c;
}

inline tpsa_tmp_
operator- (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_copy(a.get(), c.get());
  mad_tpsa_set0(c.get(), 1, -b);
  return c;
}

inline tpsa_tmp_
operator- (const tpsa_tmp_ &a, num_t b) {  TRC("tmp,num")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_set0(c.get(), -1, b);
  return c;
}

inline tpsa_tmp_
operator- (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  tpsa_tmp_ c(newt(b));
  mad_tpsa_scl (b.get(),-1, c.get());
  mad_tpsa_set0(c.get(), 1, a);
  return c;
}

inline tpsa_tmp_
operator- (num_t a, const tpsa_tmp_ &b) {  TRC("num,tmp")
  tpsa_tmp_ c(cct(b).release());
  mad_tpsa_scl (c.get(),-1, c.get());
  mad_tpsa_set0(c.get(), 1, a);
  return c;
}

inline tpsa_tmp_ operator-(const tpsa      &a, const tpsa      &b) { TRC("tpa,tpa") return *a-*b; }
inline tpsa_tmp_ operator-(const tpsa      &a, const tpsa_ref  &b) { TRC("tpa,ref") return *a- b; }
inline tpsa_tmp_ operator-(const tpsa      &a, const tpsa_tmp_ &b) { TRC("tpa,tmp") return *a- b; }
inline tpsa_tmp_ operator-(const tpsa_ref  &a, const tpsa      &b) { TRC("ref,tpa") return  a-*b; }
inline tpsa_tmp_ operator-(const tpsa_tmp_ &a, const tpsa      &b) { TRC("tmp,tpa") return  a-*b; }
inline tpsa_tmp_ operator-(const tpsa      &a,       num_t      b) { TRC("tpa,num") return *a- b; }
inline tpsa_tmp_ operator-(      num_t      a, const tpsa      &b) { TRC("num,tpa") return  a-*b; }

// --- mul ---

inline tpsa_tmp_
operator* (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_mul(a.get(), b.get(), c.get());
  return c;
}

//inline tpsa_tmp_
//operator* (const tpsa_tmp_ &a, const tpsa_tmp_ &b) {  TRC("tmp,tmp")
//  tpsa_tmp_ c(cct(a).release());
//  mad_tpsa_mul(c.get(), b.get(), c.get());
//  return c;
//}

//inline tpsa_tmp_
//operator* (const tpsa_tmp_ &a, const tpsa_ref &b) {  TRC("tmp,ref")
//  tpsa_tmp_ c(cct(a).release());
//  mad_tpsa_mul(c.get(), b.get(), c.get());
//  return c;
//}

//inline tpsa_tmp_
//operator* (const tpsa_ref &a, const tpsa_tmp_ &b) {  TRC("ref,tmp")
//  tpsa_tmp_ c(cct(b).release());
//  mad_tpsa_mul(a.get(), c.get(), c.get());
//  return c;
//}

inline tpsa_tmp_
operator* (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_scl(a.get(), b, c.get());
  return c;
}

inline tpsa_tmp_
operator* (const tpsa_tmp_ &a, num_t b) {  TRC("tmp,num")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_scl(c.get(), b, c.get());
  return c;
}

inline tpsa_tmp_ operator*(const tpsa      &a, const tpsa      &b) { TRC("tpa,tpa") return *a**b; }
inline tpsa_tmp_ operator*(const tpsa      &a, const tpsa_ref  &b) { TRC("tpa,ref") return *a* b; }
//inline tpsa_tmp_ operator*(const tpsa      &a, const tpsa_tmp_ &b) { TRC("tpa,tmp") return *a* b; }
inline tpsa_tmp_ operator*(const tpsa_ref  &a, const tpsa      &b) { TRC("ref,tpa") return  a**b; }
//inline tpsa_tmp_ operator*(const tpsa_tmp_ &a, const tpsa      &b) { TRC("tmp,tpa") return  a**b; }
inline tpsa_tmp_ operator*(const tpsa      &a,       num_t      b) { TRC("tpa,num") return *a* b; }
inline tpsa_tmp_ operator*(      num_t      a, const tpsa      &b) { TRC("num,tpa") return *b* a; }
inline tpsa_tmp_ operator*(      num_t      a, const tpsa_ref  &b) { TRC("num,ref") return  b* a; }
inline tpsa_tmp_ operator*(      num_t      a, const tpsa_tmp_ &b) { TRC("num,tmp") return  b* a; }


// --- div ---

inline tpsa_tmp_
operator/ (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_div(a.get(), b.get(), c.get());
  return c;
}

//inline tpsa_tmp_
//operator/ (const tpsa_tmp_ &a, const tpsa_tmp_ &b) {  TRC("tmp,tmp")
//  tpsa_tmp_ c(cct(a).release());
//  mad_tpsa_div(c.get(), b.get(), c.get());
//  return c;
//}

//inline tpsa_tmp_
//operator/ (const tpsa_tmp_ &a, const tpsa_ref &b) {  TRC("tmp,ref")
//  tpsa_tmp_ c(cct(a).release());
//  mad_tpsa_div(c.get(), b.get(), c.get());
//  return c;
//}

//inline tpsa_tmp_
//operator/ (const tpsa_ref &a, const tpsa_tmp_ &b) {  TRC("ref,tmp")
//  tpsa_tmp_ c(cct(b).release());
//  mad_tpsa_div(a.get(), c.get(), c.get());
//  return c;
//}

inline tpsa_tmp_
operator/ (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_scl(a.get(), 1/b, c.get());
  return c;
}

inline tpsa_tmp_
operator/ (const tpsa_tmp_ &a, num_t b) {  TRC("tmp,num")
  tpsa_tmp_ c(cct(a).release());
  mad_tpsa_scl(c.get(), 1/b, c.get());
  return c;
}

inline tpsa_tmp_
operator/ (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  tpsa_tmp_ c(newt(b));
  mad_tpsa_inv(b.get(), a, c.get());
  return c;
}

//inline tpsa_tmp_
//operator/ (num_t a, const tpsa_tmp_ &b) {  TRC("num,tmp")
//  tpsa_tmp_ c(cct(b).release());
//  mad_tpsa_inv(c.get(), a, c.get());
//  return c;
//}

inline tpsa_tmp_ operator/(const tpsa      &a, const tpsa      &b) { TRC("tpa,tpa") return *a/ *b; }
inline tpsa_tmp_ operator/(const tpsa      &a, const tpsa_ref  &b) { TRC("tpa,ref") return *a/  b; }
//inline tpsa_tmp_ operator/(const tpsa      &a, const tpsa_tmp_ &b) { TRC("tpa,tmp") return *a/  b; }
inline tpsa_tmp_ operator/(const tpsa_ref  &a, const tpsa      &b) { TRC("ref,tpa") return  a/ *b; }
//inline tpsa_tmp_ operator/(const tpsa_tmp_ &a, const tpsa      &b) { TRC("tmp,tpa") return  a/ *b; }
inline tpsa_tmp_ operator/(const tpsa      &a,       num_t      b) { TRC("tpa,num") return *a/  b; }
inline tpsa_tmp_ operator/(      num_t      a, const tpsa      &b) { TRC("num,tpa") return  a/ *b; }

// --- pow ---

inline tpsa_tmp_
pow (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_pow(a.get(), b.get(), c.get());
  return c;
}

//inline tpsa_tmp_
//pow (const tpsa_tmp_ &a, const tpsa_tmp_ &b) {  TRC("tmp,tmp")
//  tpsa_tmp_ c(cct(a).release());
//  mad_tpsa_pow(c.get(), b.get(), c.get());
//  return c;
//}

//inline tpsa_tmp_
//pow (const tpsa_tmp_ &a, const tpsa_ref &b) {  TRC("tmp,ref")
//  tpsa_tmp_ c(cct(a).release());
//  mad_tpsa_pow(c.get(), b.get(), c.get());
//  return c;
//}

//inline tpsa_tmp_
//pow (const tpsa_ref &a, const tpsa_tmp_ &b) {  TRC("ref,tmp")
//  tpsa_tmp_ c(cct(b).release());
//  mad_tpsa_pow(a.get(), c.get(), c.get());
//  return c;
//}

inline tpsa_tmp_
pow (const tpsa_ref &a, int b) {  TRC("ref,int")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_powi(a.get(), b, c.get());
  return c;
}

//inline tpsa_tmp_
//pow (const tpsa_tmp_ &a, int b) {  TRC("tmp,int")
//  tpsa_tmp_ c(cct(a).release());
//  mad_tpsa_powi(c.get(), b, c.get());
//  return c;
//}

inline tpsa_tmp_
pow (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  tpsa_tmp_ c(newt(a));
  mad_tpsa_pown(a.get(), b, c.get());
  return c;
}

//inline tpsa_tmp_
//pow (const tpsa_tmp_ &a, num_t b) {  TRC("tmp,num")
//  tpsa_tmp_ c(cct(a).release());
//  mad_tpsa_pown(c.get(), b, c.get());
//  return c;
//}

inline tpsa_tmp_
pow (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  tpsa_tmp_ c(newt(b));
  mad_tpsa_scl(b.get(), std::log(a), c.get());
  mad_tpsa_exp(c.get(), c.get());
  return c;
}

//inline tpsa_tmp_
//pow (num_t a, const tpsa_tmp_ &b) {  TRC("num,tmp")
//  tpsa_tmp_ c(cct(b).release());
//  mad_tpsa_scl(c.get(), std::log(a), c.get());
//  mad_tpsa_exp(c.get(), c.get());
//  return c;
//}

inline tpsa_tmp_ pow(const tpsa      &a, const tpsa      &b) { TRC("tpa,tpa") return pow(*a,*b); }
inline tpsa_tmp_ pow(const tpsa      &a, const tpsa_ref  &b) { TRC("tpa,ref") return pow(*a, b); }
//inline tpsa_tmp_ pow(const tpsa      &a, const tpsa_tmp_ &b) { TRC("tpa,tmp") return pow(*a, b); }
inline tpsa_tmp_ pow(const tpsa_ref  &a, const tpsa      &b) { TRC("ref,tpa") return pow( a,*b); }
//inline tpsa_tmp_ pow(const tpsa_tmp_ &a, const tpsa      &b) { TRC("tmp,tpa") return pow( a,*b); }
inline tpsa_tmp_ pow(const tpsa      &a,       int        b) { TRC("tpa,int") return pow(*a, b); }
inline tpsa_tmp_ pow(const tpsa      &a,       num_t      b) { TRC("tpa,num") return pow(*a, b); }
inline tpsa_tmp_ pow(      num_t      a, const tpsa      &b) { TRC("num,tpa") return pow( a,*b); }

// warning: the operator ^ hasn't the expected precedence and associativity...

inline tpsa_tmp_ operator^(const tpsa      &a, const tpsa      &b) { TRC("tpa,tpa") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa      &a, const tpsa_ref  &b) { TRC("tpa,ref") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa      &a, const tpsa_tmp_ &b) { TRC("tpa,tmp") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_ref  &a, const tpsa      &b) { TRC("ref,tpa") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_ref  &a, const tpsa_ref  &b) { TRC("ref,ref") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_ref  &a, const tpsa_tmp_ &b) { TRC("ref,tmp") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_tmp_ &a, const tpsa      &b) { TRC("tmp,tpa") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_tmp_ &a, const tpsa_ref  &b) { TRC("tmp,ref") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_tmp_ &a, const tpsa_tmp_ &b) { TRC("tmp,tmp") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa      &a,       int        b) { TRC("tpa,int") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_ref  &a,       int        b) { TRC("ref,int") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_tmp_ &a,       int        b) { TRC("tmp,int") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa      &a,       num_t      b) { TRC("tpa,num") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_ref  &a,       num_t      b) { TRC("ref,num") return pow(a,b); }
inline tpsa_tmp_ operator^(const tpsa_tmp_ &a,       num_t      b) { TRC("tmp,num") return pow(a,b); }
inline tpsa_tmp_ operator^(      num_t      a, const tpsa      &b) { TRC("num,tpa") return pow(a,b); }
inline tpsa_tmp_ operator^(      num_t      a, const tpsa_ref  &b) { TRC("num,ref") return pow(a,b); }
inline tpsa_tmp_ operator^(      num_t      a, const tpsa_tmp_ &b) { TRC("num,tmp") return pow(a,b); }

} // mad

// --- functions ---

namespace mad {

inline num_t     sqr(      num_t     a) { TRC("num") return a*a; }
inline tpsa_tmp_ sqr(const tpsa     &a) { TRC("tpa") return a*a; }
inline tpsa_tmp_ sqr(const tpsa_ref &a) { TRC("ref") return a*a; }

// --- unary ---

#define FUN(F) \
\
inline tpsa_tmp_ F (const tpsa_ref &a) {  TRC("ref") \
  tpsa_tmp_ c(newt(a)); \
  mad_tpsa_ ## F (a.get(), c.get()); \
  return c; \
} \
\
inline tpsa_tmp_ F (const tpsa &a) { TRC("tpa") return F(*a); }

/*
inline tpsa_tmp_ F (const tpsa_tmp_ &a) { TRC("tmp") \
  tpsa_tmp_ c(cct(a).release()); \
  mad_tpsa_ ## F (c.get(), c.get()); \
  return c; \
} \
*/

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
