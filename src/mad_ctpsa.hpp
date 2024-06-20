#ifndef MAD_CTPSA_HPP
#define MAD_CTPSA_HPP

/*
 o-----------------------------------------------------------------------------o
 |
 | Simple C++ wrapper to GTPSA - complex version
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

// --- includes ---------------------------------------------------------------o

#include <complex>
#include "mad_tpsa.hpp"

extern "C" {
#include "mad_ctpsa.h"
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

#define C(...) (*(cpx_t*)&(__VA_ARGS__))
#define RE(a)  (((num_t*)&(a))[0])
#define IM(a)  (((num_t*)&(a))[1])

namespace mad {

// forward decl
struct ctpsa;
struct ctpsa_ref;
namespace mad_prv_ { struct ctpsa_tmp_; }
using CPX = std::complex<double>;

// public abstract class implementing interface for ctpsa and ctpsa_ref.
// use static polymorphism, i.e. CRTP + ADL + CRT for efficiency.
template <class D>
struct ctpsa_base {
  // getters
  ctpsa_t* ptr () const { return static_cast<const D*>( this)->ptr(); }
  ctpsa_t& ref () const { return static_cast<const D*>( this)->ref(); }
  const D& self() const { return static_cast<const D&>(*this);        }
        D& self()       { return static_cast<      D&>(*this);        }
  ord_t    mo  () const { return mad_ctpsa_ord(ptr(), false);         }

  // set name
  D& set(str_t              s)  { mad_ctpsa_nam(ptr(), s);             return self(); }
  D& set(const std::string &s)  { mad_ctpsa_nam(ptr(), s.c_str());     return self(); }

  // set value and variables
  D& set (CPX a)                 { mad_ctpsa_setval(ptr(), C(a)      ); return self(); }
  D& set (CPX a, idx_t v)        { mad_ctpsa_setvar(ptr(), C(a), v, 0); return self(); }
  D& setp(CPX a, idx_t v)        { mad_ctpsa_setprm(ptr(), C(a), v   ); return self(); }

  // debug
  log_t isvalid() const { return mad_ctpsa_isvalid(ptr()); }

  template <class A>
  D& operator+=(const ctpsa_base<A> &a) { TRC("baz,baz") mad_ctpsa_add (ptr(),a.ptr(),ptr()); return self(); }
  D& operator+=(      CPX            a) { TRC("baz,num") mad_ctpsa_seti(ptr(), 0,   1, C(a)); return self(); }

  template <class A>
  D& operator-=(const ctpsa_base<A> &a) { TRC("baz,baz") mad_ctpsa_sub (ptr(),a.ptr(),ptr()); return self(); }
  D& operator-=(      CPX            a) { TRC("baz,num") mad_ctpsa_seti(ptr(),0,1,C(a=-a,a)); return self(); }

  template <class A>
  D& operator*=(const ctpsa_base<A> &a) { TRC("baz,baz") mad_ctpsa_mul (ptr(),a.ptr(),ptr()); return self(); }
  D& operator*=(      CPX            a) { TRC("baz,num") mad_ctpsa_scl (ptr(),   C(a),ptr()); return self(); }

  template <class A>
  D& operator/=(const ctpsa_base<A> &a) { TRC("baz,baz") mad_ctpsa_div (ptr(),a.ptr(),ptr()); return self(); }
  D& operator/=(      CPX            a) { TRC("baz,num") mad_ctpsa_scl (ptr(),C(a=1./a,a),ptr()); return self(); }

  template <class A>
  D& operator^=(const ctpsa_base<A> &a) { TRC("baz,baz") mad_ctpsa_pow (ptr(),a.ptr(),ptr()); return self(); }
  D& operator^=(      CPX            a) { TRC("baz,num") mad_ctpsa_pown(ptr(),   C(a),ptr()); return self(); }
  D& operator^=(      int            a) { TRC("baz,int") mad_ctpsa_powi(ptr(),     a ,ptr()); return self(); }

  // indexing by index, monomial as string (literal), and (sparse) monomial as vector.
  CPX operator[](idx_t i) const {
    cpx_t r = mad_ctpsa_geti(ptr(),i);
    return CPX(RE(r), IM(r));
  }
  CPX operator[](str_t s) const {
    cpx_t r = mad_ctpsa_gets(ptr(),0,s);
    return CPX(RE(r), IM(r));
  }
  CPX operator[](const std::string& s) const {
    cpx_t r = mad_ctpsa_gets(ptr(), s.size(), s.c_str());
    return CPX(RE(r), IM(r));
  }
  CPX operator[](const std::vector<ord_t>& m) const {
    cpx_t r = mad_ctpsa_getm(ptr(), m.size(), m.data());
    return CPX(RE(r), IM(r));
  }
  CPX operator[](const std::vector<idx_t>& m) const {
    cpx_t r = mad_ctpsa_getsm(ptr(), m.size(), m.data());
    return CPX(RE(r), IM(r));
  }

protected:
  ctpsa_base()                               = default; // abstract class
  ctpsa_base(ctpsa_base<D>&&)                = delete;  // move    ctor
  ctpsa_base(const ctpsa_base<D>&)           = delete;  // copy    ctor
  ctpsa_base(std::nullptr_t)                 = delete;  // nullptr ctor
  ctpsa_base operator=(ctpsa_base<D>&&)      = delete;  // move    assign
  ctpsa_base operator=(const ctpsa_base<D>&) = delete;  // copy    assign
  ctpsa_base operator=(std::nullptr_t)       = delete;  // nullptr assign
};

// public class to wrap tpsa_t* _without_ memory management.
struct ctpsa_ref : ctpsa_base<ctpsa_ref> {
  ctpsa_ref(ctpsa_t *a) : t(a) { TRC("ctpsa_t* %p", (void*)a) }

  ctpsa_t* ptr () const { return  t; }
  ctpsa_t& ref () const { return *t; }
  void     swp (ctpsa_ref &a) { std::swap(t,a.t); }

  template <class A>
  ctpsa_ref& operator=(const ctpsa_base<A> &a) { TRC("ref=baz") mad_ctpsa_copy(a.ptr(),ptr()); return *this; }
  ctpsa_ref& operator=(const ctpsa_ref     &a) { TRC("ref=ref") mad_ctpsa_copy(a.ptr(),ptr()); return *this; }
  ctpsa_ref& operator=(      ctpsa_ref    &&a) { TRC("ref<ref") mad_ctpsa_copy(a.ptr(),ptr()); return *this; }
  ctpsa_ref& operator=(      CPX            a) { TRC("ref=num") mad_ctpsa_setval(ptr(), C(a)); return *this; }

private:
  ctpsa_ref()                                = delete;  // final   class
  ctpsa_ref(ctpsa_ref&&)                     = delete;  // move    ctor
  ctpsa_ref(const ctpsa_ref&)                = delete;  // copy    ctor
  ctpsa_ref(std::nullptr_t)                  = delete;  // nullptr ctor
//ctpsa_ref& operator=(ctpsa_ref&&)          = delete;  // move    assign
//ctpsa_ref& operator=(const ctpsa_ref&)     = delete;  // copy    assign
  ctpsa_ref& operator=(std::nullptr_t)       = delete;  // nullptr assign

private:
  ctpsa_t *t;
};

// public class to wrap tpsa_t*[] _without_ memory management.
struct ctpsa_refs {
  explicit ctpsa_refs(ctpsa_t *a[], int na) : n(na), t(a) { TRC("ctpsa_t** %p", (void*) a) }

  const ctpsa_ref operator[](int i) const { TRC("ref[]") assert(i>=0 && i<n); return ctpsa_ref(t[i]); }

private:
  ctpsa_refs()                               = delete;  // final   class
  ctpsa_refs(ctpsa_refs&&)                   = delete;  // move    ctor
  ctpsa_refs(const ctpsa_refs&)              = delete;  // copy    ctor
  ctpsa_refs(std::nullptr_t)                 = delete;  // nullptr ctor
  ctpsa_refs& operator=(ctpsa_refs&&)        = delete;  // move    assign
  ctpsa_refs& operator=(const ctpsa_refs&)   = delete;  // copy    assign
  ctpsa_refs& operator=(std::nullptr_t)      = delete;  // nullptr assign

private:
  int       n;
  ctpsa_t **t;
};

// public class handle tpsa_t* _with_ memory management.
struct ctpsa : ctpsa_base<ctpsa> {
#if TPSA_USE_DFT
  explicit ctpsa()         : t(mad_ctpsa_newd(mad_desc_curr,dflt)) { TRC("dft! %p", (void*)t.get()) }
#endif
  explicit ctpsa(ord_t mo) : t(mad_ctpsa_newd(mad_desc_curr,mo  )) { TRC("ord! %p", (void*)t.get()) }

  explicit ctpsa(const ctpsa &a)                   : t(mad_ctpsa_new(a.ptr(),same)) { TRC("&tpa! %p", (void*)t.get()) }
  template <class A>
  explicit ctpsa(const ctpsa_base<A> &a)           : t(mad_ctpsa_new(a.ptr(),same)) { TRC("&baz! %p", (void*)t.get()) }
  template <class A>
  explicit ctpsa(const ctpsa_base<A> &a, ord_t mo) : t(mad_ctpsa_new(a.ptr(),mo  )) { TRC("&baz,ord! %p", (void*)t.get()) }
  template <class A>
  explicit ctpsa(const ctpsa_base<A> &a,
                 const ctpsa_base<A> &b) : t(mad_ctpsa_new(a.ptr(),std::max(a.mo(), b.mo()))) { TRC("&baz,&baz! %p", (void*)t.get()) }

  template <class A>
  explicit ctpsa(const tpsa_base<A> &re) { TRC("&baz! %p", (void*)t.get()); mad_ctpsa_cplx(re.ptr(), NULL, t.get()); }
  template <class A>
  explicit ctpsa(const tpsa_base<A> &re,
                 const tpsa_base<A> &im) { TRC("&baz,&baz! %p", (void*)t.get()); mad_ctpsa_cplx(re.ptr(), im.ptr(), t.get()); }

  ctpsa_t* ptr () const { return t.get(); }
  ctpsa_t& ref () const { return *t;      }
  void     swp (ctpsa &a) { t.swap(a.t);  }

  template <class A>
  ctpsa& operator=(const ctpsa_base<A> &a) { TRC("tpa=baz") mad_ctpsa_copy(a.ptr(),ptr()); return *this; }
  ctpsa& operator=(const ctpsa         &a) { TRC("tpa=tpa") mad_ctpsa_copy(a.ptr(),ptr()); return *this; }
  ctpsa& operator=(      ctpsa        &&a) { TRC("tpa<tpa") mad_ctpsa_copy(a.ptr(),ptr()); return *this; }
  ctpsa& operator=(      CPX            a) { TRC("tpa=num") mad_ctpsa_setval(ptr(), C(a)); return *this; }

  template <class A>
  ctpsa& operator=(const  tpsa_base<A> &a) { TRC("tpa=baz") mad_ctpsa_cplx(a.ptr(),NULL,ptr()); return *this; }

#if TPSA_USE_TMP // specialization for capturing temporaries (forward decl)
  ctpsa(const mad_prv_::ctpsa_tmp_&); //
  ctpsa(const mad_prv_::ctpsa_tmp_&, const mad_prv_::ctpsa_tmp_&);
  template <class A>
  ctpsa(const mad_prv_::ctpsa_tmp_&, const ctpsa_base<A> &);
  template <class A>
  ctpsa(const ctpsa_base<A> &      , const mad_prv_::ctpsa_tmp_&);
#endif

protected:
  explicit ctpsa(ctpsa_t *a) : t(a) { TRC("ctpsa_t* %p", (void*)a) }

private:
#if !TPSA_USE_DFT
  ctpsa()                                    = delete;  // dflt    ctor
#endif
#if TPSA_USE_TMP
  ctpsa(ctpsa&&)                             = delete;  // move    ctor
#endif
  ctpsa(std::nullptr_t)                      = delete;  // nullptr ctor
//ctpsa& operator=(ctpsa&&)                  = delete;  // move    assign
//ctpsa& operator=(const ctpsa&)             = delete;  // copy    assign
  ctpsa& operator=(std::nullptr_t)           = delete;  // nullptr assign

  friend std::FILE* operator>>(std::FILE*, ctpsa&);

private:
  struct ctpsa_del_ {
    void operator()(ctpsa_t *t) { TRC("~tpa! %p", (void*)t) mad_ctpsa_del(t); }
  };

  std::unique_ptr<ctpsa_t, ctpsa_del_> t;
};

#if TPSA_USE_TMP

// private class to manage temporaries in expressions, i.e. save allocations.
namespace mad_prv_ {

struct ctpsa_tmp_ : ctpsa {
#if TPSA_USE_DFT
  explicit ctpsa_tmp_() : ctpsa()                  { TRC("dft! %p", (void*)ptr()) }
#endif
  ctpsa_tmp_(ctpsa_t *a) : ctpsa(a)                { TRC("*tmp") } // capture ptr (see scan)
  ctpsa_tmp_(ctpsa_tmp_ &&a) : ctpsa(std::move(a)) { TRC("<tmp") } // move ctor
  ctpsa_tmp_(const ctpsa_tmp_ &a) : ctpsa(a)       { TRC("&tmp") } // copy ctor

  template <class A>
  ctpsa_tmp_(const ctpsa_base<A> &a) : ctpsa(a)    { TRC("&baz") }

  ctpsa_tmp_(const ctpsa_tmp_ &a,
             const ctpsa_tmp_ &b)    : ctpsa(a,b)  { TRC("&tmp,&tmp") }
  template <class A>
  ctpsa_tmp_(const ctpsa_base<A> &a,
             const ctpsa_tmp_    &b) : ctpsa(a,b)  { TRC("&baz,&tmp") }
  template <class A>
  ctpsa_tmp_(const ctpsa_tmp_    &a,
             const ctpsa_base<A> &b) : ctpsa(a,b)  { TRC("&tmp,&baz") }
  template <class A, class B>
  ctpsa_tmp_(const ctpsa_base<A> &a,
             const ctpsa_base<B> &b) : ctpsa(a,b)  { TRC("&baz,&baz") }

private:
#if !TPSA_USE_DFT
  ctpsa_tmp_()                               = delete; // dflt    ctor
#endif
//ctpsa_tmp_(ctpsa_tmp_ &&)                  = delete; // move    ctor
//ctpsa_tmp_(const ctpsa_tmp_ &)             = delete; // copy    ctor
  ctpsa_tmp_(std::nullptr_t)                 = delete; // nullptr ctor
  ctpsa_tmp_& operator=(ctpsa_tmp_&&)        = delete; // move    assign
  ctpsa_tmp_& operator=(const ctpsa_tmp_&)   = delete; // copy    assign
  ctpsa_tmp_& operator=(std::nullptr_t)      = delete; // nullptr assign
};

} // mad_prv_

#define R mad_prv_:: tpsa_tmp_
#define T mad_prv_::ctpsa_tmp_

inline ctpsa::ctpsa(const T &a) { t.swap(const_cast<T&>(a).t); TRC("&tmp") }

inline ctpsa::ctpsa(const T &a, const T &b) {
  if (a.mo() >= b.mo()) { t.swap(const_cast<T&>(a).t); TRC("&tmp,&tmp,1") }
  else {                  t.swap(const_cast<T&>(b).t); TRC("&tmp,&tmp,2") }
}

template <class A>
inline ctpsa::ctpsa(const T &a, const ctpsa_base<A> &b) {
  if (a.mo() >= b.mo()) { t.swap(const_cast<T&>(a).t); TRC("&tmp,&baz,1") }
  else { T tb(b);         t.swap(               tb.t); TRC("&tmp,&baz,2") }
}

template <class A>
inline ctpsa::ctpsa(const ctpsa_base<A> &a, const T &b) {
  if (b.mo() >= a.mo()) { t.swap(const_cast<T&>(b).t); TRC("&baz,&tmp,1") }
  else { T ta(a);         t.swap(               ta.t); TRC("&baz,&tmp,2") }
}

#else
#define R  tpsa
#define T ctpsa
#endif // TPSA_USE_TMP

// --- operators --------------------------------------------------------------o

template <class A>
inline R real (const ctpsa_base<A> &a) {  TRC("re(baz)")
  R c(a.mo()); mad_ctpsa_real(a.ptr(),c.ptr()); return c;
}

template <class A>
inline R imag (const ctpsa_base<A> &a) {  TRC("im(baz)")
  R c(a.mo()); mad_ctpsa_imag(a.ptr(),c.ptr()); return c;
}

// --- unary ---

template <class A>
inline T operator+ (const ctpsa_base<A> &a) {  TRC("+baz")
  T c(a); return c;
}

template <class A>
inline T operator- (const ctpsa_base<A> &a) {  TRC("-baz")
  T c(a); mad_ctpsa_scl(a.ptr(), -1, c.ptr()); return c;
}

#if TPSA_USE_TMP

inline T operator+ (const T &a) {  TRC("+tmp")
  T c(a); return c;
}

inline T operator- (const T &a) {  TRC("-tmp")
  T c(a); mad_ctpsa_scl(c.ptr(), -1, c.ptr()); return c;
}

#endif // TPSA_USE_TMP

// --- add ---

template <class A, class B>
inline T operator+ (const ctpsa_base<A> &a, const ctpsa_base<B> &b) {  TRC("baz+baz")
  T c(a,b); mad_ctpsa_add(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator+ (CPX a, const ctpsa_base<A> &b) {  TRC("num+baz")
  T c(b);
  mad_ctpsa_copy(b.ptr(), c.ptr());
  mad_ctpsa_seti(c.ptr(), 0, 1, C(a));
  return c;
}

template <class A>
inline T operator+ (const ctpsa_base<A> &a, CPX b) {  TRC("baz+num")
  return b+a;
}

#if TPSA_USE_TMP

template <class A>
inline T operator+ (const ctpsa_base<A> &a, const T &b) {  TRC("baz+tmp")
  T c(a,b);
  if (a.ptr()) mad_ctpsa_add(a.ptr(), c.ptr(), c.ptr());
  else         mad_ctpsa_add(c.ptr(), b.ptr(), c.ptr());
  return c;
}

template <class A>
inline T operator+ (const T &a, const ctpsa_base<A> &b) {  TRC("tmp+baz")
  T c(a,b);
  if (a.ptr()) mad_ctpsa_add(a.ptr(), c.ptr(), c.ptr());
  else         mad_ctpsa_add(c.ptr(), b.ptr(), c.ptr());
  return c;
}

inline T operator+ (const T &a, const T &b) {  TRC("tmp+tmp")
  T c(a,b);
  if (a.ptr()) mad_ctpsa_add(a.ptr(), c.ptr(), c.ptr());
  else         mad_ctpsa_add(c.ptr(), b.ptr(), c.ptr());
  return c;
}

inline T operator+ (CPX a, const T &b) {  TRC("num+tmp")
  T c(b); mad_ctpsa_seti(c.ptr(), 0, 1, C(a)); return c;
}

inline T operator+ (const T &a, CPX b) {  TRC("tmp+num")
  return b+a;
}

#endif // TPSA_USE_TMP

// --- sub ---

template <class A, class B>
inline T operator- (const ctpsa_base<A> &a, const ctpsa_base<B> &b) {  TRC("baz-baz")
  T c(a,b); mad_ctpsa_sub(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator- (const ctpsa_base<A> &a, CPX b) {  TRC("baz-num")
  T c(a);
  mad_ctpsa_copy(a.ptr(), c.ptr());
  mad_ctpsa_seti(c.ptr(), 0, 1, C(b=-b,b));
  return c;
}

template <class A>
inline T operator- (CPX a, const ctpsa_base<A> &b) {  TRC("num-baz")
  T c(b);
  mad_ctpsa_scl (b.ptr(),-1, c.ptr());
  mad_ctpsa_seti(c.ptr(), 0, 1, C(a));
  return c;
}

#if TPSA_USE_TMP

template <class A>
inline T operator- (const ctpsa_base<A> &a, const T &b) {  TRC("baz-tmp")
  T c(a,b);
  if (a.ptr()) mad_ctpsa_sub(a.ptr(), c.ptr(), c.ptr());
  else         mad_ctpsa_sub(c.ptr(), b.ptr(), c.ptr());
  return c;
}

template <class A>
inline T operator- (const T &a, const ctpsa_base<A> &b) {  TRC("tmp-baz")
  T c(a,b);
  if (a.ptr()) mad_ctpsa_sub(a.ptr(), c.ptr(), c.ptr());
  else         mad_ctpsa_sub(c.ptr(), b.ptr(), c.ptr());
  return c;
}

inline T operator- (const T &a, const T &b) {  TRC("tmp-tmp")
  T c(a,b);
  if (a.ptr()) mad_ctpsa_sub(a.ptr(), c.ptr(), c.ptr());
  else         mad_ctpsa_sub(c.ptr(), b.ptr(), c.ptr());
  return c;
}

inline T operator- (const T &a, CPX b) {  TRC("tmp-num")
  T c(a); mad_ctpsa_seti(c.ptr(), 0, 1, C(b=-b,b)); return c;
}

inline T operator- (CPX a, const T &b) {  TRC("num-tmp")
  T c(b);
  mad_ctpsa_scl (c.ptr(),-1, c.ptr());
  mad_ctpsa_seti(c.ptr(), 0, 1, C(a));
  return c;
}

#endif // TPSA_USE_TMP

// --- mul ---

template <class A, class B>
inline T operator* (const ctpsa_base<A> &a, const ctpsa_base<B> &b) {  TRC("baz*baz")
  T c(a,b); mad_ctpsa_mul(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator* (CPX a, const ctpsa_base<A> &b) {  TRC("num*baz")
  T c(b); mad_ctpsa_scl(b.ptr(), C(a), c.ptr()); return c;
}

template <class A>
inline T operator* (const ctpsa_base<A> &a, CPX b) {  TRC("baz*num")
  return b*a;
}

#if TPSA_USE_TMP

inline T operator*(CPX a, const T &b) { TRC("num*tmp")
  T c(b); mad_ctpsa_scl(c.ptr(), C(a), c.ptr()); return c;
}

inline T operator* (const T &a, CPX b) {  TRC("tmp*num")
  return b*a;
}

#endif // TPSA_USE_TMP

// --- div ---

template <class A, class B>
inline T operator/ (const ctpsa_base<A> &a, const ctpsa_base<B> &b) {  TRC("baz/baz")
  T c(a,b); mad_ctpsa_div(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T operator/ (const ctpsa_base<A> &a, CPX b) {  TRC("baz/num")
  T c(a); mad_ctpsa_scl(a.ptr(), C(b=1./b,b), c.ptr()); return c;
}

template <class A>
inline T operator/ (CPX a, const ctpsa_base<A> &b) {  TRC("num/baz")
  T c(b); mad_ctpsa_inv(b.ptr(), C(a), c.ptr()); return c;
}

#if TPSA_USE_TMP

inline T operator/ (const T &a, CPX b) {  TRC("tmp/num")
  T c(a); mad_ctpsa_scl(c.ptr(), C(b=1./b,b), c.ptr()); return c;
}

inline T operator/ (CPX a, const T &b) {  TRC("num/tmp")
  T c(b); mad_ctpsa_inv(c.ptr(), C(a), c.ptr()); return c;
}

#endif // TPSA_USE_TMP

// --- pow ---

template <class A, class B>
inline T pow (const ctpsa_base<A> &a, const ctpsa_base<B> &b) {  TRC("baz^baz")
  T c(a,b); mad_ctpsa_pow(a.ptr(), b.ptr(), c.ptr()); return c;
}

template <class A>
inline T pow (const ctpsa_base<A> &a, int b) {  TRC("baz^int")
  T c(a); mad_ctpsa_powi(a.ptr(), b, c.ptr()); return c;
}

template <class A>
inline T pow (const ctpsa_base<A> &a, CPX b) {  TRC("baz^num")
  T c(a); mad_ctpsa_pown(a.ptr(), C(b), c.ptr()); return c;
}

template <class A>
inline T pow (CPX a, const ctpsa_base<A> &b) {  TRC("num^baz")
  T c(b);
  mad_ctpsa_scl(b.ptr(), C(a=std::log(a),a), c.ptr());
  mad_ctpsa_exp(c.ptr(), c.ptr());
  return c;
}

#if TPSA_USE_TMP

inline T pow (const T &a, int b) {  TRC("tmp^int")
  T c(a); mad_ctpsa_powi(c.ptr(), b, c.ptr()); return c;
}

inline T pow (const T &a, CPX b) {  TRC("tmp^num")
  T c(a); mad_ctpsa_pown(c.ptr(), C(b), c.ptr()); return c;
}

inline T pow (CPX a, const T &b) {  TRC("num^tmp")
  T c(b);
  mad_ctpsa_scl(c.ptr(), C(a=std::log(a),a), c.ptr());
  mad_ctpsa_exp(c.ptr(), c.ptr());
  return c;
}

#endif // TPSA_USE_TMP

// --- hypot ---

template <class A, class B>
inline T hypot (const ctpsa_base<A> &a, const ctpsa_base<B> &b) {  TRC("baz,baz")
  T c(a); mad_ctpsa_hypot(a.ptr(), b.ptr(), c.ptr()); return c;
}

// --- I/O ---

template <class A>
inline std::FILE* operator<< (std::FILE *out, const ctpsa_base<A> &a) {
  mad_ctpsa_print(a.ptr(), 0,0,0, out); return out;
}

inline std::FILE* operator<< (std::FILE *out, const ctpsa_t &a) {
  mad_ctpsa_print(&a, 0,0,0, out); return out;
}

inline std::FILE* operator>> (std::FILE *in, ctpsa &a) {
  T c(mad_ctpsa_scan(in)); a = c; return in;
}

// --- swap ---

inline void swap (ctpsa_ref &a, ctpsa_ref &b) { TRC("ref") a.swp(b); }
inline void swap (ctpsa     &a, ctpsa     &b) { TRC("tpa") a.swp(b); }
inline void swap (CPX       &a, CPX       &b) { TRC("cpx") std::swap(a,b); }

// --- functions ---

inline CPX fval   (CPX a         ) { TRC("cpx") return a; }
inline CPX nrm    (CPX a         ) { TRC("cpx") return std::abs(a); }
inline CPX sqr    (CPX a         ) { TRC("cpx") return a*a; }
inline CPX inv    (CPX a, CPX v=1) { TRC("cpx") return v/a; }
inline CPX invsqrt(CPX a, CPX v=1) { TRC("cpx") return v/std::sqrt(a); }
inline CPX sinc   (CPX a         ) { TRC("cpx") cpx_t r = mad_cpx_sinc (C(a)); return CPX(RE(r),IM(r)); }
inline CPX sinhc  (CPX a         ) { TRC("cpx") cpx_t r = mad_cpx_sinhc(C(a)); return CPX(RE(r),IM(r)); }
inline CPX asinc  (CPX a         ) { TRC("cpx") cpx_t r = mad_cpx_asinc(C(a)); return CPX(RE(r),IM(r)); }

inline CPX fval(const ctpsa_t *a) { TRC("tpsa")
  return mad_ctpsa_geti(a,0);
}

template <class A>
inline CPX fval (const ctpsa_base<A> &a) { TRC("baz")
  return a[0];
}

template <class A>
inline CPX nrm (const ctpsa_base<A> &a) { TRC("baz")
  return mad_ctpsa_nrm(a.ptr());
}

template <class A>
inline T sqr (const ctpsa_base<A> &a) { TRC("baz")
  T c(a); mad_ctpsa_mul(a.ptr(), a.ptr(), c.ptr()); return c;
}

template <class A>
inline T inv (const ctpsa_base<A> &a, CPX v=1) { TRC("baz")
  T c(a); mad_ctpsa_inv(a.ptr(), C(v), c.ptr()); return c;
}

template <class A>
inline T invsqrt (const ctpsa_base<A> &a, CPX v=1) { TRC("baz")
  T c(a); mad_ctpsa_invsqrt(a.ptr(), C(v), c.ptr()); return c;
}

#if TPSA_USE_TMP

inline T inv (const T &a, CPX v=1) { TRC("tmp")
  T c(a); mad_ctpsa_inv(c.ptr(), C(v), c.ptr()); return c;
}

inline T invsqrt (const T &a, CPX v=1) { TRC("tmp")
  T c(a); mad_ctpsa_invsqrt(c.ptr(), C(v), c.ptr()); return c;
}

#endif // TPSA_USE_TMP

// --- unary ---

#define FUN(F) \
template <class A> \
inline T F (const ctpsa_base<A> &a) {  TRC("baz") \
  T c(a); mad_ctpsa_ ## F (a.ptr(), c.ptr()); return c; \
} \
FUN_TMP(F)

#if TPSA_USE_TMP

#define FUN_TMP(F) \
inline T F (const T &a) { TRC("tmp") \
  T c(a); mad_ctpsa_ ## F (c.ptr(), c.ptr()); return c; \
}

#else
#define FUN_TMP(F)
#endif // TPSA_USE_TMP

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
#undef R
#undef C
#undef RE
#undef IM
#undef TRC
#undef FUN
#undef FUN_TMP

//#undef TPSA_USE_TMP
//#undef TPSA_USE_TRC

// --- end --------------------------------------------------------------------o

#endif // MAD_CTPSA_HPP
