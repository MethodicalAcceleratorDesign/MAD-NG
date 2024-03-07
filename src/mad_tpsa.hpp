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
#include <cassert>
#include <memory>
#include <string>
#include <vector>

extern "C" {
#include "mad_num.h"
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

// constants
const int dflt = mad_tpsa_default;
const int same = mad_tpsa_same;

// forward decl
struct tpsa;
struct tpsa_ref;
namespace mad_prv_ { struct tpsa_tmp_; }

// public abstract class implementing interface for tpsa and tpsa_ref.
// use static polymorphism, i.e. CRTP + ADL + CRT for efficiency.
template <class D>
struct tpsa_base {
  // getters
  tpsa_t*  ptr () const { return static_cast<const D*>( this)->ptr(); }
  tpsa_t&  ref () const { return static_cast<const D*>( this)->ref(); }
  const D& self() const { return static_cast<const D&>(*this);        }
        D& self()       { return static_cast<      D&>(*this);        }

  // set name
  D& set(str_t              s)         { mad_tpsa_setnam(ptr(), s);          return self(); }
  D& set(const std::string &s)         { mad_tpsa_setnam(ptr(), s.c_str());  return self(); }

  // set value and variables
  D& set (num_t a)                     { mad_tpsa_setval(ptr(), a      );    return self(); }
  D& set (num_t a, idx_t v)            { mad_tpsa_setvar(ptr(), a, v, 0);    return self(); }
  D& setp(num_t a, idx_t v)            { mad_tpsa_setprm(ptr(), a, v   );    return self(); }

  // debug
  log_t isvalid() const { return mad_ctpsa_isvalid(ptr()); }

  template <class A>
  D& operator+=(const tpsa_base<A> &a) { TRC("baz,baz") mad_tpsa_add (ptr(),a.ptr(),ptr()); return self(); }
  D& operator+=(      num_t         a) { TRC("baz,num") mad_tpsa_set0(ptr(),      1,    a); return self(); }

  template <class A>
  D& operator-=(const tpsa_base<A> &a) { TRC("baz,baz") mad_tpsa_sub (ptr(),a.ptr(),ptr()); return self(); }
  D& operator-=(      num_t         a) { TRC("baz,num") mad_tpsa_set0(ptr(),      1,   -a); return self(); }

  template <class A>
  D& operator*=(const tpsa_base<A> &a) { TRC("baz,baz") mad_tpsa_mul (ptr(),a.ptr(),ptr()); return self(); }
  D& operator*=(      num_t         a) { TRC("baz,num") mad_tpsa_scl (ptr(),      a,ptr()); return self(); }

  template <class A>
  D& operator/=(const tpsa_base<A> &a) { TRC("baz,baz") mad_tpsa_div (ptr(),a.ptr(),ptr()); return self(); }
  D& operator/=(      num_t         a) { TRC("baz,num") mad_tpsa_scl (ptr(),    1/a,ptr()); return self(); }

  template <class A>
  D& operator^=(const tpsa_base<A> &a) { TRC("baz,baz") mad_tpsa_pow (ptr(),a.ptr(),ptr()); return self(); }
  D& operator^=(      num_t         a) { TRC("baz,num") mad_tpsa_pown(ptr(),      a,ptr()); return self(); }
  D& operator^=(      int           a) { TRC("baz,int") mad_tpsa_powi(ptr(),      a,ptr()); return self(); }

  // indexing by index, monomial as string (literal), and (sparse) monomial as vector.
  num_t operator[](idx_t i) const { return i ? mad_tpsa_geti(ptr(),  i) : mad_tpsa_get0(ptr()); }
  num_t operator[](str_t s) const { return     mad_tpsa_gets(ptr(),0,s);                        }
  num_t operator[](const std::string& s) const {
    return mad_tpsa_gets(ptr(), s.size(), s.c_str());
  }
  num_t operator[](const std::vector<ord_t>& m) const {
    return mad_tpsa_getm(ptr(), m.size(), m.data());
  }
  num_t operator[](const std::vector<idx_t>& m) const {
    return mad_tpsa_getsm(ptr(), m.size(), m.data());
  }

protected:
  tpsa_base()                              = default; // abstract class
  tpsa_base(tpsa_base<D>&&)                = delete;  // move    ctor
  tpsa_base(const tpsa_base<D>&)           = delete;  // copy    ctor
  tpsa_base(std::nullptr_t)                = delete;  // nullptr ctor
  tpsa_base operator=(tpsa_base<D>&&)      = delete;  // move    assign
  tpsa_base operator=(const tpsa_base<D>&) = delete;  // copy    assign
  tpsa_base operator=(std::nullptr_t)      = delete;  // nullptr assign
};

// public class to wrap tpsa_t* _without_ memory management.
struct tpsa_ref : tpsa_base<tpsa_ref> {
  tpsa_ref(tpsa_t *a) : t(a) { TRC("tpsa_t* %p", (void*)a) }

  tpsa_t* ptr () const { return  t; }
  tpsa_t& ref () const { return *t; }
  void    swp (tpsa_ref &a) { std::swap(t,a.t); }

  template <class A>
  tpsa_ref& operator=(const tpsa_base<A> &a) { TRC("ref=baz") mad_tpsa_copy(a.ptr(),ptr()); return *this; }
  tpsa_ref& operator=(const tpsa_ref     &a) { TRC("ref=ref") mad_tpsa_copy(a.ptr(),ptr()); return *this; }
  tpsa_ref& operator=(      tpsa_ref    &&a) { TRC("ref<ref") mad_tpsa_copy(a.ptr(),ptr()); return *this; }
  tpsa_ref& operator=(      num_t         a) { TRC("ref=num") mad_tpsa_setval(ptr(), a   ); return *this; }

private:
  tpsa_ref()                               = delete;  // final   class
  tpsa_ref(tpsa_ref&&)                     = delete;  // move    ctor
  tpsa_ref(const tpsa_ref&)                = delete;  // copy    ctor
  tpsa_ref(std::nullptr_t)                 = delete;  // nullptr ctor
//tpsa_ref& operator=(tpsa_ref&&)          = delete;  // move    assign
//tpsa_ref& operator=(const tpsa_ref&)     = delete;  // copy    assign
  tpsa_ref& operator=(std::nullptr_t)      = delete;  // nullptr assign

private:
  tpsa_t *t;
};

// public class to wrap tpsa_t*[] _without_ memory management.
struct tpsa_refs {
  explicit tpsa_refs(tpsa_t *a[], int na) : n(na), t(a) { TRC("tpsa_t** %p", (void*) a) }

  const tpsa_ref operator[](int i) const { TRC("ref[]") assert(i>=0 && i<n); return tpsa_ref(t[i]); }

private:
  tpsa_refs()                              = delete;  // final   class
  tpsa_refs(tpsa_refs&&)                   = delete;  // move    ctor
  tpsa_refs(const tpsa_refs&)              = delete;  // copy    ctor
  tpsa_refs(std::nullptr_t)                = delete;  // nullptr ctor
  tpsa_refs& operator=(tpsa_refs&&)        = delete;  // move    assign
  tpsa_refs& operator=(const tpsa_refs&)   = delete;  // copy    assign
  tpsa_refs& operator=(std::nullptr_t)     = delete;  // nullptr assign

private:
  int      n;
  tpsa_t **t;
};

// public class handle tpsa_t* _with_ memory management.
struct tpsa : tpsa_base<tpsa> {
  explicit tpsa()                : t(mad_tpsa_newd(mad_desc_curr, dflt)) { TRC("dft! %p", (void*)t.get()) }
  explicit tpsa(int mo)          : t(mad_tpsa_newd(mad_desc_curr, mo  )) { TRC("int! %p", (void*)t.get()) }
  explicit tpsa(str_t s)         : t(mad_tpsa_newd(mad_desc_curr, dflt)) { TRC("dft,str! %p", (void*)t.get()) this->set(s); }
  explicit tpsa(str_t s, int mo) : t(mad_tpsa_newd(mad_desc_curr, mo  )) { TRC("int,str! %p", (void*)t.get()) this->set(s); }

  template <class A>
  explicit tpsa(const tpsa_base<A> &a)         : t(mad_tpsa_new(a.ptr(), dflt)) { TRC("&baz! %p", (void*)t.get()) }
  template <class A>
  explicit tpsa(const tpsa_base<A> &a, int mo) : t(mad_tpsa_new(a.ptr(), mo  )) { TRC("&baz,int! %p", (void*)t.get()) }

  tpsa_t* ptr () const  { return t.get(); }
  tpsa_t& ref () const  { return *t;      }
  void    swp (tpsa &a) { t.swap(a.t);    }

  template <class A>
  tpsa& operator=(const tpsa_base<A> &a) { TRC("tpa=baz") mad_tpsa_copy(a.ptr(),ptr()); return *this; }
  tpsa& operator=(const tpsa         &a) { TRC("tpa=tpa") mad_tpsa_copy(a.ptr(),ptr()); return *this; }
  tpsa& operator=(      tpsa        &&a) { TRC("tpa<tpa") mad_tpsa_copy(a.ptr(),ptr()); return *this; }
  tpsa& operator=(      num_t         a) { TRC("tpa=num") mad_tpsa_setval(ptr(), a   ); return *this; }

#if TPSA_USE_TMP // specialization for capturing temporaries
  tpsa(const mad_prv_::tpsa_tmp_&); // forward decl
#else            // replace specialization but slow...
  tpsa(const tpsa &a) : t(mad_tpsa_new(a.ptr(), dflt)) { TRC("&tpa!! %p", (void*)t.get())
    mad_tpsa_copy(a.ptr(),ptr());
  }
#endif

protected:
  explicit tpsa(tpsa_t *a) : t(a) { TRC("tpsa_t* %p", (void*)a) }

private:
//tpsa()                                   = delete;  // dflt    ctor
#if TPSA_USE_TMP
  tpsa(tpsa&&)                             = delete;  // move    ctor
  tpsa(const tpsa&)                        = delete;  // copy    ctor
#endif
  tpsa(std::nullptr_t)                     = delete;  // nullptr ctor
//tpsa& operator=(tpsa&&)                  = delete;  // move    assign
//tpsa& operator=(const tpsa&)             = delete;  // copy    assign
  tpsa& operator=(std::nullptr_t)          = delete;  // nullptr assign

  friend std::FILE* operator>>(std::FILE*, tpsa&);

private:
  struct tpsa_del_ {
    void operator()(tpsa_t *t) { TRC("~tpa! %p", (void*)t) mad_tpsa_del(t); }
  };

  std::unique_ptr<tpsa_t, tpsa_del_> t;
};

#if TPSA_USE_TMP

// private class to manage temporaries in expressions, i.e. save allocations.
namespace mad_prv_ {

struct tpsa_tmp_ : tpsa {
  explicit tpsa_tmp_() : tpsa()                 { TRC("dft! %p", (void*)ptr()) }
  tpsa_tmp_(tpsa_tmp_ &&a) : tpsa(std::move(a)) { TRC("<tmp") } // move ctor
  tpsa_tmp_(const tpsa_tmp_ &a) : tpsa(a)       { TRC("&tmp") } // copy ctor

  template <class A>
  explicit tpsa_tmp_(const tpsa_base<A> &a) : tpsa(a) { TRC("&baz") }
  explicit tpsa_tmp_(tpsa_t *a)             : tpsa(a) { TRC("*tmp") } // capture ptr (see scan)

private:
//tpsa_tmp_()                              = delete; // dflt    ctor
//tpsa_tmp_(tpsa_tmp_ &&)                  = delete; // move    ctor
//tpsa_tmp_(const tpsa_tmp_ &)             = delete; // copy    ctor
  tpsa_tmp_(std::nullptr_t)                = delete; // nullptr ctor
  tpsa_tmp_& operator=(tpsa_tmp_&&)        = delete; // move    assign
  tpsa_tmp_& operator=(const tpsa_tmp_&)   = delete; // copy    assign
  tpsa_tmp_& operator=(std::nullptr_t)     = delete; // nullptr assign
};

} // mad_prv_

#define T mad_prv_::tpsa_tmp_
inline tpsa::tpsa(const T &a) { t.swap(const_cast<T&>(a).t); TRC("&tmp") }
#else
#define T tpsa
#endif // TPSA_USE_TMP

// --- operators --------------------------------------------------------------o

// --- unary ---

template <class A>
inline T operator+ (const tpsa_base<A> &a) {  TRC("+baz")
  T c(a); return c;
}

template <class A>
inline T operator- (const tpsa_base<A> &a) {  TRC("-baz")
  T c(a); mad_tpsa_scl(a.ptr(), -1, c.ptr()); return c;
}

#if TPSA_USE_TMP

inline T operator+ (const T &a) {  TRC("+tmp")
  T c(a); return c;
}

inline T operator- (const T &a) {  TRC("-tmp")
  T c(a); mad_tpsa_scl(c.ptr(), -1, c.ptr()); return c;
}

#endif // TPSA_USE_TMP

// --- add ---

template <class A, class B>
inline T operator+ (const tpsa_base<A> &a, const tpsa_base<B> &b) {  TRC("baz+baz")
  T c(a); mad_tpsa_add(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator+ (num_t a, const tpsa_base<A> &b) {  TRC("num+baz")
  T c(b);
  mad_tpsa_copy(b.ptr(), c.ptr());
  mad_tpsa_set0(c.ptr(), 1, a);
  return c;
}

template <class A>
inline T operator+ (const tpsa_base<A> &a, num_t b) {  TRC("baz+num")
  return b+a;
}

#if TPSA_USE_TMP

template <class A>
inline T operator+ (const tpsa_base<A> &a, const T &b) {  TRC("baz+tmp")
  T c(b); mad_tpsa_add(a.ptr(), c.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator+ (const T &a, const tpsa_base<A> &b) {  TRC("tmp+baz")
  T c(a); mad_tpsa_add(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T operator+ (const T &a, const T &b) {  TRC("tmp+tmp")
  T c(a); mad_tpsa_add(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T operator+ (num_t a, const T &b) {  TRC("num+tmp")
  T c(b); mad_tpsa_set0(c.ptr(), 1, a); return c;
}

inline T operator+ (const T &a, num_t b) {  TRC("tmp+num")
  return b+a;
}

#endif // TPSA_USE_TMP

// --- sub ---

template <class A, class B>
inline T operator- (const tpsa_base<A> &a, const tpsa_base<B> &b) {  TRC("baz-baz")
  T c(a); mad_tpsa_sub(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator- (const tpsa_base<A> &a, num_t b) {  TRC("baz-num")
  T c(a);
  mad_tpsa_copy(a.ptr(), c.ptr());
  mad_tpsa_set0(c.ptr(), 1, -b);
  return c;
}

template <class A>
inline T operator- (num_t a, const tpsa_base<A> &b) {  TRC("num-baz")
  T c(b);
  mad_tpsa_scl (b.ptr(),-1, c.ptr());
  mad_tpsa_set0(c.ptr(), 1, a);
  return c;
}

#if TPSA_USE_TMP

template <class A>
inline T operator- (const tpsa_base<A> &a, const T &b) {  TRC("baz-tmp")
  T c(b); mad_tpsa_sub(a.ptr(), c.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator- (const T &a, const tpsa_base<A> &b) {  TRC("tmp-baz")
  T c(a); mad_tpsa_sub(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T operator- (const T &a, const T &b) {  TRC("tmp-tmp")
  T c(a); mad_tpsa_sub(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T operator- (const T &a, num_t b) {  TRC("tmp-num")
  T c(a); mad_tpsa_set0(c.ptr(), 1, -b); return c;
}

inline T operator- (num_t a, const T &b) {  TRC("num-tmp")
  T c(b);
  mad_tpsa_scl (c.ptr(),-1, c.ptr());
  mad_tpsa_set0(c.ptr(), 1, a);
  return c;
}

#endif // TPSA_USE_TMP

// --- mul ---

template <class A, class B>
inline T operator* (const tpsa_base<A> &a, const tpsa_base<B> &b) {  TRC("baz*baz")
  T c(a); mad_tpsa_mul(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator* (num_t a, const tpsa_base<A> &b) {  TRC("num*baz")
  T c(b); mad_tpsa_scl(b.ptr(), a, c.ptr()); return c;
}

template <class A>
inline T operator* (const tpsa_base<A> &a, num_t b) {  TRC("baz*num")
  return b*a;
}

#if TPSA_USE_TMP

template <class A>
inline T operator* (const tpsa_base<A> &a, const T &b) {  TRC("baz*tmp")
  T c(b); mad_tpsa_mul(a.ptr(), c.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator* (const T &a, const tpsa_base<A> &b) {  TRC("tmp*baz")
  T c(a); mad_tpsa_mul(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T operator* (const T &a, const T &b) {  TRC("tmp*tmp")
  T c(a); mad_tpsa_mul(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T operator*(num_t a, const T &b) { TRC("num*tmp")
  T c(b); mad_tpsa_scl(c.ptr(), a, c.ptr()); return c;
}

inline T
operator* (const T &a, num_t b) {  TRC("tmp*num")
  return b*a;
}

#endif // TPSA_USE_TMP

// --- div ---

template <class A, class B>
inline T operator/ (const tpsa_base<A> &a, const tpsa_base<B> &b) {  TRC("baz/baz")
  T c(a); mad_tpsa_div(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator/ (const tpsa_base<A> &a, num_t b) {  TRC("baz/num")
  T c(a); mad_tpsa_scl(a.ptr(), 1/b, c.ptr()); return c;
}

template <class A>
inline T operator/ (num_t a, const tpsa_base<A> &b) {  TRC("num/baz")
  T c(b); mad_tpsa_inv(b.ptr(), a, c.ptr()); return c;
}

#if TPSA_USE_TMP

template <class A>
inline T operator/ (const tpsa_base<A> &a, const T &b) {  TRC("baz/tmp")
  T c(b); mad_tpsa_div(a.ptr(), c.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator/ (const T &a, const tpsa_base<A> &b) {  TRC("tmp/baz")
  T c(a); mad_tpsa_div(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T operator/ (const T &a, const T &b) {  TRC("tmp/tmp")
  T c(a); mad_tpsa_div(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T operator/ (const T &a, num_t b) {  TRC("tmp/num")
  T c(a); mad_tpsa_scl(c.ptr(), 1/b, c.ptr()); return c;
}

inline T operator/ (num_t a, const T &b) {  TRC("num/tmp")
  T c(b); mad_tpsa_inv(c.ptr(), a, c.ptr()); return c;
}

#endif // TPSA_USE_TMP

// --- pow ---

template <class A, class B>
inline T pow (const tpsa_base<A> &a, const tpsa_base<B> &b) {  TRC("baz^baz")
  T c(a); mad_tpsa_pow(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T pow (const tpsa_base<A> &a, int b) {  TRC("baz^int")
  T c(a); mad_tpsa_powi(a.ptr(), b, c.ptr()); return c;
}

template <class A>
inline T pow (const tpsa_base<A> &a, num_t b) {  TRC("baz^num")
  T c(a); mad_tpsa_pown(a.ptr(), b, c.ptr()); return c;
}

template <class A>
inline T pow (num_t a, const tpsa_base<A> &b) {  TRC("num^baz")
  T c(b);
  mad_tpsa_scl(b.ptr(), log(a), c.ptr());
  mad_tpsa_exp(c.ptr(), c.ptr());
  return c;
}

#if TPSA_USE_TMP

template <class A>
inline T pow (const tpsa_base<A> &a, const T &b) {  TRC("baz^tmp")
  T c(b); mad_tpsa_pow(a.ptr(), c.ptr(), c.ptr()); return c;
}

template <class A>
inline T pow (const T &a, const tpsa_base<A> &b) {  TRC("tmp^baz")
  T c(a); mad_tpsa_pow(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T pow (const T &a, const T &b) {  TRC("tmp^tmp")
  T c(a); mad_tpsa_pow(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T pow (const T &a, int b) {  TRC("tmp^int")
  T c(a); mad_tpsa_powi(c.ptr(), b, c.ptr()); return c;
}

inline T pow (const T &a, num_t b) {  TRC("tmp^num")
  T c(a); mad_tpsa_pown(c.ptr(), b, c.ptr()); return c;
}

inline T pow (num_t a, const T &b) {  TRC("num^tmp")
  T c(b);
  mad_tpsa_scl(c.ptr(), log(a), c.ptr());
  mad_tpsa_exp(c.ptr(), c.ptr());
  return c;
}

#endif // TPSA_USE_TMP

// --- atan2, hypot ---

template <class A, class B>
inline T atan2 (const tpsa_base<A> &a, const tpsa_base<B> &b) {  TRC("baz,baz")
  T c(a); mad_tpsa_atan2(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A, class B>
inline T hypot (const tpsa_base<A> &a, const tpsa_base<B> &b) {  TRC("baz,baz")
  T c(a); mad_tpsa_hypot(a.ptr(), b.ptr(), c.ptr()); return c;
}

#if TPSA_USE_TMP

template <class A>
inline T atan2 (const tpsa_base<A> &a, const T &b) {  TRC("baz,tmp")
  T c(b); mad_tpsa_atan2(a.ptr(), c.ptr(), c.ptr()); return c;
}

template <class A>
inline T atan2 (const T &a, const tpsa_base<A> &b) {  TRC("tmp,baz")
  T c(a); mad_tpsa_atan2(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T atan2 (const T &a, const T &b) {  TRC("tmp,tmp")
  T c(a); mad_tpsa_atan2(c.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T hypot (const tpsa_base<A> &a, const T &b) {  TRC("baz,tmp")
  T c(b); mad_tpsa_hypot(a.ptr(), c.ptr(), c.ptr()); return c;
}

template <class A>
inline T hypot (const T &a, const tpsa_base<A> &b) {  TRC("tmp,baz")
  T c(a); mad_tpsa_hypot(c.ptr(), b.ptr(), c.ptr()); return c;
}

inline T hypot (const T &a, const T &b) {  TRC("tmp,tmp")
  T c(a); mad_tpsa_hypot(c.ptr(), b.ptr(), c.ptr()); return c;
}

#endif // TPSA_USE_TMP

// --- I/O ---

template <class A>
inline std::FILE* operator<< (std::FILE *out, const tpsa_base<A> &a) {
  mad_tpsa_print(a.ptr(), 0,0,0, out); return out;
}

inline std::FILE* operator<< (std::FILE *out, const tpsa_t &a) {
  mad_tpsa_print(&a, 0,0,0, out); return out;
}

inline std::FILE* operator>> (std::FILE *in, tpsa &a) {
  T c(mad_tpsa_scan(in)); a = c; return in;
}

// --- swap ---

inline void swap (tpsa_ref &a, tpsa_ref &b) { TRC("ref") a.swp(b); }
inline void swap (tpsa     &a, tpsa     &b) { TRC("tpa") a.swp(b); }
inline void swap (num_t    &a, num_t    &b) { TRC("num") std::swap(a,b); }

// --- functions ---

inline num_t fval   (num_t a           ) { TRC("num") return a; }
inline num_t nrm    (num_t a           ) { TRC("num") return abs(a); }
inline num_t sqr    (num_t a           ) { TRC("num") return a*a; }
inline num_t inv    (num_t a, num_t v=1) { TRC("num") return v/a; }
inline num_t invsqrt(num_t a, num_t v=1) { TRC("num") return v/sqrt(a); }
inline num_t sinc   (num_t a           ) { TRC("num") return mad_num_sinc (a); }
inline num_t sinhc  (num_t a           ) { TRC("num") return mad_num_sinhc(a); }
inline num_t asinc  (num_t a           ) { TRC("num") return mad_num_asinc(a); }

inline num_t fval(const tpsa_t *a) { TRC("tspa")
  return mad_tpsa_get0(a);
}

template <class A>
inline num_t fval (const tpsa_base<A> &a) { TRC("baz")
  return a[0];
}

inline num_t fabs (const tpsa_t *a) { TRC("tspa")
  return abs(mad_tpsa_get0(a));
}

template <class A>
inline num_t fabs (const tpsa_base<A> &a) { TRC("baz")
  return abs(a[0]);
}

template <class A>
inline num_t nrm (const tpsa_base<A> &a) { TRC("baz")
  return mad_tpsa_nrm(a.ptr());
}

template <class A>
inline T sqr (const tpsa_base<A> &a) { TRC("baz")
  T c(a); mad_tpsa_mul(a.ptr(), a.ptr(), c.ptr()); return c;
}

template <class A>
inline T inv (const tpsa_base<A> &a, num_t v=1) { TRC("baz")
  T c(a); mad_tpsa_inv(a.ptr(), v, c.ptr()); return c;
}

template <class A>
inline T invsqrt (const tpsa_base<A> &a, num_t v=1) { TRC("baz")
  T c(a); mad_tpsa_invsqrt(a.ptr(), v, c.ptr()); return c;
}

#if TPSA_USE_TMP

inline T sqr (const T &a) { TRC("tmp")
  T c(a); mad_tpsa_mul(c.ptr(), c.ptr(), c.ptr()); return c;
}

inline T inv (const T &a, num_t v=1) { TRC("tmp")
  T c(a); mad_tpsa_inv(c.ptr(), v, c.ptr()); return c;
}

inline T invsqrt (const T &a, num_t v=1) { TRC("tmp")
  T c(a); mad_tpsa_invsqrt(c.ptr(), v, c.ptr()); return c;
}

#endif // TPSA_USE_TMP

// --- unary ---

#define FUN(F) \
template <class A> \
inline T F (const tpsa_base<A> &a) {  TRC("baz") \
  T c(a); mad_tpsa_ ## F (a.ptr(), c.ptr()); return c; \
} \
FUN_TMP(F)

#if TPSA_USE_TMP

#define FUN_TMP(F) \
inline T F (const T &a) { TRC("tmp") \
  T c(a); mad_tpsa_ ## F (c.ptr(), c.ptr()); return c; \
}

#else
#define FUN_TMP(F)
#endif // TPSA_USE_TMP

FUN(abs   );
FUN(unit  );
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

//#undef TPSA_USE_TMP
//#undef TPSA_USE_TRC

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_HPP
