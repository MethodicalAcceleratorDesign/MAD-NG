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

// comment to disable temporaries and traces
#define TPSA_USE_TMP 1
//#define TPSA_USE_TRC 1

// --- includes ---------------------------------------------------------------o

#include <cmath>
#include <cstdio>
#include <memory>
#include <string>
#include <vector>

extern "C" {
#include "mad_tpsa.h"
}

// --- trace ------------------------------------------------------------------o

#if TPSA_USE_TRC
#define TRC(...) \
  (printf("%s:%3d:%12s: ", __FILE__, __LINE__, __func__), \
   printf(__VA_ARGS__), printf("\n"));
#else
#define TRC(...)
#endif

// --- types ------------------------------------------------------------------o

namespace mad {

// private class to build dtor as a type and avoid extra function ptr in tpsa.
namespace mad_prv_ {
struct tpsa_del_ {
  void operator()(tpsa_t *t) { TRC("tpsa_t* %p", (void*)t) mad_tpsa_del(t); }
};}

// public class to wrap tpsa_t with memory management.
using tpsa = std::unique_ptr<tpsa_t, mad_prv_::tpsa_del_>;

// public class to *locally* wrap tpsa_t without memory management.
struct tpsa_ref {
  explicit tpsa_ref(tpsa_t *a) : t(*a) { TRC("tpsa_t* %p", (void*)a      ) }
  explicit tpsa_ref(tpsa   &a) : t(*a) { TRC("tpsa& %p"  , (void*)a.get()) }
           tpsa_ref(tpsa_t &a) : t( a) { TRC("tpsa_t& %p", (void*)&a     ) }

  tpsa_ref()                           = delete; // default ctor
  tpsa_ref(tpsa_ref&&)                 = delete; // move    ctor
  tpsa_ref(const tpsa_ref&)            = delete; // copy    ctor
  tpsa_ref(std::nullptr_t)             = delete; // nullptr ctor
  tpsa_ref& operator=(tpsa_ref&&)      = delete; // move    assign
//tpsa_ref& operator=(const tpsa_ref&) = delete; // copy    assign (see below)
  tpsa_ref& operator=(std::nullptr_t)  = delete; // nullptr assign

  tpsa_t* get()       const { return &t; }
  tpsa_t& operator*() const { return  t; }

  void operator= (const tpsa     &a) { mad_tpsa_copy(   &*a,&t); }
  void operator= (const tpsa_ref &a) { mad_tpsa_copy(   &*a,&t); }
  void operator= (      num_t     a) { mad_tpsa_set0(&t,  0, a); }

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

  void operator^=(const tpsa     &a) { mad_tpsa_pow (&t,&*a,&t); }
  void operator^=(const tpsa_ref &a) { mad_tpsa_pow (&t,&*a,&t); }
  void operator^=(      num_t     a) { mad_tpsa_pown(&t,  a,&t); }
  void operator^=(      int       a) { mad_tpsa_powi(&t,  a,&t); }

  num_t operator[](idx_t i) const { return mad_tpsa_geti(&t,  i); }
  num_t operator[](str_t s) const { return mad_tpsa_gets(&t,0,s); }

  num_t operator[](const std::string& s) const {
    return mad_tpsa_gets(&t, s.size(), s.c_str());
  }
  num_t operator[](const std::vector<ord_t>& m) const {
    return mad_tpsa_getm(&t, m.size(), m.data());
  }
  num_t operator[](const std::vector<idx_t>& m) const {
    return mad_tpsa_getsm(&t, m.size(), m.data());
  }

private:
  tpsa_t &t;
};

#if TPSA_USE_TMP

// private class to manage temporaries in expressions, i.e. save allocations.
namespace mad_prv_ {
struct tpsa_tmp_ : tpsa {
  explicit
  tpsa_tmp_(tpsa_t     *a) : tpsa(a) { TRC("tpsa_t* %p", (void*)a) }
  tpsa_tmp_(tpsa_tmp_ &&a) : tpsa(std::move(a)) {}
};
inline tpsa_tmp_& cct(const tpsa_tmp_ &a) { return const_cast<tpsa_tmp_&>(a); }
} // mad_prv_

#define T mad_prv_::tpsa_tmp_
#else
#define T tpsa
#endif // TPSA_USE_TMP

} // mad

// --- operators --------------------------------------------------------------o

namespace mad {

// --- ctors ---

inline T
newt_ (const desc_t *d=mad_desc_curr, int mo=mad_tpsa_default) {
  return T(mad_tpsa_newd(d, mo));
}

inline T
newt_ (const desc_t *d, int mo, str_t nam) {
  T c(mad_tpsa_newd(d, mo));
  mad_tpsa_setnam(c.get(), nam);
  return c;
}

inline T
newt () {  TRC("void")
  return newt_();
}

inline T
newt (str_t nam) { TRC("void")
  return newt_(mad_desc_curr, mad_tpsa_default, nam);
}

inline T
newt (int mo) {  TRC("int")
  return newt_(mad_desc_curr, mo);
}

inline T
newt (int mo, str_t nam) {  TRC("int")
  return newt_(mad_desc_curr, mo, nam);
}

inline T
newt (const tpsa_ref &a) {  TRC("ref")
  return newt_(mad_tpsa_desc(a.get()));
}

inline T
newt (const tpsa_ref &a, int mo) {  TRC("ref")
  return newt_(mad_tpsa_desc(a.get()), mo);
}

inline T
newt (const tpsa_ref &a, str_t nam) { TRC("ref")
  return newt_(mad_tpsa_desc(a.get()), mad_tpsa_default, nam);
}

inline T
newt (const tpsa_ref &a, int mo, str_t nam) {  TRC("ref")
  return newt_(mad_tpsa_desc(a.get()), mo, nam);
}

inline T
newt (const tpsa &a) {  TRC("tpa")
  return newt_(mad_tpsa_desc(a.get()));
}

inline T
newt (const tpsa &a, int mo) {  TRC("tpa")
  return newt_(mad_tpsa_desc(a.get()), mo);
}

inline T
newt (const tpsa &a, str_t nam) { TRC("tpa")
  return newt_(mad_tpsa_desc(a.get()), mad_tpsa_default, nam);
}

inline T
newt (const tpsa &a, int mo, str_t nam) {  TRC("tpa")
  return newt_(mad_tpsa_desc(a.get()), mo, nam);
}

// --- unary ---

inline T
operator- (const tpsa_ref &a) {  TRC("ref")
  T c(newt(a));
  mad_tpsa_scl(a.get(), -1, c.get());
  return c;
}

inline T operator-(const tpsa &a) { TRC("tpa") return -*a; }

#if TPSA_USE_TMP

inline T
operator- (const T &a) {  TRC("tmp")
  T c(cct(a).release());
  mad_tpsa_scl(c.get(), -1, c.get());
  return c;
}

#endif // TPSA_USE_TMP

// --- add ---

inline T
operator+ (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  T c(newt(a));
  mad_tpsa_add(a.get(), b.get(), c.get());
  return c;
}

inline T
operator+ (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  T c(newt(a));
  mad_tpsa_copy(a.get(), c.get());
  mad_tpsa_set0(c.get(), 1, b);
  return c;
}

inline T operator+(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return *a+*b; }
inline T operator+(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return *a+ b; }
inline T operator+(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return  a+*b; }
inline T operator+(const tpsa     &a,       num_t     b) { TRC("tpa,num") return *a+ b; }
inline T operator+(      num_t     a, const tpsa     &b) { TRC("num,tpa") return *b+ a; }
inline T operator+(      num_t     a, const tpsa_ref &b) { TRC("num,ref") return  b+ a; }

#if TPSA_USE_TMP

inline T
operator+ (const T &a, const T &b) {  TRC("tmp,tmp")
  T c(cct(a).release());
  mad_tpsa_add(c.get(), b.get(), c.get());
  return c;
}

inline T
operator+ (const T &a, const tpsa_ref &b) {  TRC("tmp,ref")
  T c(cct(a).release());
  mad_tpsa_add(c.get(), b.get(), c.get());
  return c;
}

inline T
operator+ (const tpsa_ref &a, const T &b) {  TRC("ref,tmp")
  T c(cct(b).release());
  mad_tpsa_add(a.get(), c.get(), c.get());
  return c;
}

inline T
operator+ (const T &a, num_t b) {  TRC("tmp,num")
  T c(cct(a).release());
  mad_tpsa_set0(c.get(), 1, b);
  return c;
}

inline T operator+(const tpsa &a, const T    &b) { TRC("tpa,tmp") return *a+ b; }
inline T operator+(const T    &a, const tpsa &b) { TRC("tmp,tpa") return  a+*b; }
inline T operator+(      num_t a, const T    &b) { TRC("num,tmp") return  b+ a; }

#endif // TPSA_USE_TMP

// --- sub ---

inline T
operator- (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  T c(newt(a));
  mad_tpsa_sub(a.get(), b.get(), c.get());
  return c;
}

inline T
operator- (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  T c(newt(a));
  mad_tpsa_copy(a.get(), c.get());
  mad_tpsa_set0(c.get(), 1, -b);
  return c;
}

inline T
operator- (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  T c(newt(b));
  mad_tpsa_scl (b.get(),-1, c.get());
  mad_tpsa_set0(c.get(), 1, a);
  return c;
}

inline T operator-(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return *a-*b; }
inline T operator-(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return *a- b; }
inline T operator-(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return  a-*b; }
inline T operator-(const tpsa     &a,       num_t     b) { TRC("tpa,num") return *a- b; }
inline T operator-(      num_t     a, const tpsa     &b) { TRC("num,tpa") return  a-*b; }

#if TPSA_USE_TMP

inline T
operator- (const T &a, const T &b) {  TRC("tmp,tmp")
  T c(cct(a).release());
  mad_tpsa_sub(c.get(), b.get(), c.get());
  return c;
}

inline T
operator- (const T &a, const tpsa_ref &b) {  TRC("tmp,ref")
  T c(cct(a).release());
  mad_tpsa_sub(c.get(), b.get(), c.get());
  return c;
}

inline T
operator- (const tpsa_ref &a, const T &b) {  TRC("ref,tmp")
  T c(cct(b).release());
  mad_tpsa_sub(a.get(), c.get(), c.get());
  return c;
}

inline T
operator- (const T &a, num_t b) {  TRC("tmp,num")
  T c(cct(a).release());
  mad_tpsa_set0(c.get(), -1, b);
  return c;
}

inline T
operator- (num_t a, const T &b) {  TRC("num,tmp")
  T c(cct(b).release());
  mad_tpsa_scl (c.get(),-1, c.get());
  mad_tpsa_set0(c.get(), 1, a);
  return c;
}

inline T operator-(const tpsa &a, const T    &b) { TRC("tpa,tmp") return *a- b; }
inline T operator-(const T    &a, const tpsa &b) { TRC("tmp,tpa") return  a-*b; }

#endif // TPSA_USE_TMP

// --- mul ---

inline T
operator* (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  T c(newt(a));
  mad_tpsa_mul(a.get(), b.get(), c.get());
  return c;
}

inline T
operator* (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  T c(newt(a));
  mad_tpsa_scl(a.get(), b, c.get());
  return c;
}

inline T operator*(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return *a**b; }
inline T operator*(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return *a* b; }
inline T operator*(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return  a**b; }
inline T operator*(const tpsa     &a,       num_t     b) { TRC("tpa,num") return *a* b; }
inline T operator*(      num_t     a, const tpsa     &b) { TRC("num,tpa") return *b* a; }
inline T operator*(      num_t     a, const tpsa_ref &b) { TRC("num,ref") return  b* a; }

#if TPSA_USE_TMP

//inline T
//operator* (const T &a, const T &b) {  TRC("tmp,tmp")
//  T c(cct(a).release());
//  mad_tpsa_mul(c.get(), b.get(), c.get());
//  return c;
//}

//inline T
//operator* (const T &a, const tpsa_ref &b) {  TRC("tmp,ref")
//  T c(cct(a).release());
//  mad_tpsa_mul(c.get(), b.get(), c.get());
//  return c;
//}

//inline T
//operator* (const tpsa_ref &a, const T &b) {  TRC("ref,tmp")
//  T c(cct(b).release());
//  mad_tpsa_mul(a.get(), c.get(), c.get());
//  return c;
//}

inline T
operator* (const T &a, num_t b) {  TRC("tmp,num")
  T c(cct(a).release());
  mad_tpsa_scl(c.get(), b, c.get());
  return c;
}

inline T operator*(num_t a, const T &b) { TRC("num,tmp") return b*a; }

//inline T operator*(const tpsa &a, const T    &b) { TRC("tpa,tmp") return *a* b; }
//inline T operator*(const T    &a, const tpsa &b) { TRC("tmp,tpa") return  a**b; }

#endif // TPSA_USE_TMP

// --- div ---

inline T
operator/ (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  T c(newt(a));
  mad_tpsa_div(a.get(), b.get(), c.get());
  return c;
}

inline T
operator/ (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  T c(newt(a));
  mad_tpsa_scl(a.get(), 1/b, c.get());
  return c;
}

inline T
operator/ (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  T c(newt(b));
  mad_tpsa_inv(b.get(), a, c.get());
  return c;
}

inline T operator/(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return *a/ *b; }
inline T operator/(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return *a/  b; }
inline T operator/(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return  a/ *b; }
inline T operator/(const tpsa     &a,       num_t     b) { TRC("tpa,num") return *a/  b; }
inline T operator/(      num_t     a, const tpsa     &b) { TRC("num,tpa") return  a/ *b; }

#if TPSA_USE_TMP

//inline T
//operator/ (const T &a, const T &b) {  TRC("tmp,tmp")
//  T c(cct(a).release());
//  mad_tpsa_div(c.get(), b.get(), c.get());
//  return c;
//}

//inline T
//operator/ (const T &a, const tpsa_ref &b) {  TRC("tmp,ref")
//  T c(cct(a).release());
//  mad_tpsa_div(c.get(), b.get(), c.get());
//  return c;
//}

//inline T
//operator/ (const tpsa_ref &a, const T &b) {  TRC("ref,tmp")
//  T c(cct(b).release());
//  mad_tpsa_div(a.get(), c.get(), c.get());
//  return c;
//}

inline T
operator/ (const T &a, num_t b) {  TRC("tmp,num")
  T c(cct(a).release());
  mad_tpsa_scl(c.get(), 1/b, c.get());
  return c;
}

inline T
operator/ (num_t a, const T &b) {  TRC("num,tmp")
  T c(cct(b).release());
  mad_tpsa_inv(c.get(), a, c.get());
  return c;
}

//inline T operator/(const tpsa &a, const T    &b) { TRC("tpa,tmp") return *a/  b; }
//inline T operator/(const T    &a, const tpsa &b) { TRC("tmp,tpa") return  a/ *b; }

#endif // TPSA_USE_TMP

// --- pow ---

inline T
pow (const tpsa_ref &a, const tpsa_ref &b) {  TRC("ref,ref")
  T c(newt(a));
  mad_tpsa_pow(a.get(), b.get(), c.get());
  return c;
}

inline T
pow (const tpsa_ref &a, int b) {  TRC("ref,int")
  T c(newt(a));
  mad_tpsa_powi(a.get(), b, c.get());
  return c;
}

inline T
pow (const tpsa_ref &a, num_t b) {  TRC("ref,num")
  T c(newt(a));
  mad_tpsa_pown(a.get(), b, c.get());
  return c;
}

inline T
pow (num_t a, const tpsa_ref &b) {  TRC("num,ref")
  T c(newt(b));
  mad_tpsa_scl(b.get(), std::log(a), c.get());
  mad_tpsa_exp(c.get(), c.get());
  return c;
}

inline T pow(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return pow(*a,*b); }
inline T pow(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return pow(*a, b); }
inline T pow(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return pow( a,*b); }
inline T pow(const tpsa     &a,       int       b) { TRC("tpa,int") return pow(*a, b); }
inline T pow(const tpsa     &a,       num_t     b) { TRC("tpa,num") return pow(*a, b); }
inline T pow(      num_t     a, const tpsa     &b) { TRC("num,tpa") return pow( a,*b); }

#if TPSA_USE_TMP

//inline T
//pow (const T &a, const T &b) {  TRC("tmp,tmp")
//  T c(cct(a).release());
//  mad_tpsa_pow(c.get(), b.get(), c.get());
//  return c;
//}

//inline T
//pow (const T &a, const tpsa_ref &b) {  TRC("tmp,ref")
//  T c(cct(a).release());
//  mad_tpsa_pow(c.get(), b.get(), c.get());
//  return c;
//}

//inline T
//pow (const tpsa_ref &a, const T &b) {  TRC("ref,tmp")
//  T c(cct(b).release());
//  mad_tpsa_pow(a.get(), c.get(), c.get());
//  return c;
//}

inline T
pow (const T &a, int b) {  TRC("tmp,int")
  T c(cct(a).release());
  mad_tpsa_powi(c.get(), b, c.get());
  return c;
}

inline T
pow (const T &a, num_t b) {  TRC("tmp,num")
  T c(cct(a).release());
  mad_tpsa_pown(c.get(), b, c.get());
  return c;
}

inline T
pow (num_t a, const T &b) {  TRC("num,tmp")
  T c(cct(b).release());
  mad_tpsa_scl(c.get(), std::log(a), c.get());
  mad_tpsa_exp(c.get(), c.get());
  return c;
}

//inline T pow(const tpsa &a, const T    &b) { TRC("tpa,tmp") return pow(*a, b); }
//inline T pow(const T    &a, const tpsa &b) { TRC("tmp,tpa") return pow( a,*b); }

#endif // TPSA_USE_TMP

// warning: the operator ^ hasn't the expected precedence and associativity...

inline T operator^(const tpsa     &a, const tpsa     &b) { TRC("tpa,tpa") return pow(a,b); }
inline T operator^(const tpsa     &a, const tpsa_ref &b) { TRC("tpa,ref") return pow(a,b); }
inline T operator^(const tpsa_ref &a, const tpsa     &b) { TRC("ref,tpa") return pow(a,b); }
inline T operator^(const tpsa_ref &a, const tpsa_ref &b) { TRC("ref,ref") return pow(a,b); }
inline T operator^(const tpsa     &a,       int       b) { TRC("tpa,int") return pow(a,b); }
inline T operator^(const tpsa_ref &a,       int       b) { TRC("ref,int") return pow(a,b); }
inline T operator^(const tpsa     &a,       num_t     b) { TRC("tpa,num") return pow(a,b); }
inline T operator^(const tpsa_ref &a,       num_t     b) { TRC("ref,num") return pow(a,b); }
inline T operator^(      num_t     a, const tpsa     &b) { TRC("num,tpa") return pow(a,b); }
inline T operator^(      num_t     a, const tpsa_ref &b) { TRC("num,ref") return pow(a,b); }

#if TPSA_USE_TMP

inline T operator^(const T &a, int      b) { TRC("tmp,int") return pow(a,b); }
inline T operator^(const T &a, num_t    b) { TRC("tmp,num") return pow(a,b); }
inline T operator^(num_t    a, const T &b) { TRC("num,tmp") return pow(a,b); }

//inline T operator^(const tpsa_ref &a, const T        &b) { TRC("ref,tmp") return pow(a,b); }
//inline T operator^(const T        &a, const tpsa_ref &b) { TRC("tmp,ref") return pow(a,b); }
//inline T operator^(const T        &a, const T        &b) { TRC("tmp,tmp") return pow(a,b); }
//inline T operator^(const tpsa     &a, const T        &b) { TRC("tpa,tmp") return pow(a,b); }
//inline T operator^(const T        &a, const tpsa     &b) { TRC("tmp,tpa") return pow(a,b); }

#endif // TPSA_USE_TMP

// --- I/O ---

inline std::FILE* operator<<(std::FILE *out, const tpsa_ref &a) {
  mad_tpsa_print(&*a, 0,0,0, out);
  return out;
}

inline std::FILE* operator<<(std::FILE *out, const tpsa &a) {
  mad_tpsa_print(&*a, 0,0,0, out);
  return out;
}

// --- functions ---

inline num_t sqr(      num_t     a) { TRC("num") return a*a; }
inline T     sqr(const tpsa     &a) { TRC("tpa") return a*a; }
inline T     sqr(const tpsa_ref &a) { TRC("ref") return a*a; }

inline num_t inv(      num_t     a) { TRC("num") return 1/a; }
inline T     inv(const tpsa     &a) { TRC("tpa") return 1/a; }
inline T     inv(const tpsa_ref &a) { TRC("ref") return 1/a; }

inline num_t nrm(      num_t     a) { TRC("num") return std::abs(a);           }
inline num_t nrm(const tpsa     &a) { TRC("tpa") return mad_tpsa_nrm(a.get()); }
inline num_t nrm(const tpsa_ref &a) { TRC("ref") return mad_tpsa_nrm(a.get()); }

// --- unary ---

#define FUN(F) \
inline T F (const tpsa_ref &a) {  TRC("ref") \
  T c(newt(a)); \
  mad_tpsa_ ## F (a.get(), c.get()); \
  return c; \
} \
inline T F (const tpsa &a) { TRC("tpa") return F(*a); } \
FUN_TMP(F)

#if TPSA_USE_TMP

#define FUN_TMP(F) \
inline T F (const T &a) { TRC("tmp") \
  T c(cct(a).release()); \
  mad_tpsa_ ## F (c.get(), c.get()); \
  return c; \
}

#else
#define FUN_TMP(F)
#endif // TPSA_USE_TMP

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

} // mad

#undef T
#undef TRC
#undef FUN
#undef FUN_TMP

#undef TPSA_USE_TMP
#undef TPSA_USE_TRC

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_HPP
