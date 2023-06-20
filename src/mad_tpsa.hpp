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
  - Provide memory management and operator overloading for GTPSA.

 o-----------------------------------------------------------------------------o
 */

#include <cmath>
#include <memory>

// --- types ------------------------------------------------------------------o

namespace mad {

struct tpsa_del {
  // --- dtor ---
  void operator()(tpsa_t *t) { mad_tpsa_del(t); }
};

using tpsa = std::unique_ptr<tpsa_t, tpsa_del>;

}

// --- operators --------------------------------------------------------------o

namespace mad {

// --- ctors ---

inline tpsa
newt ()
{
  return tpsa(mad_tpsa_newd(mad_desc_curr, mad_tpsa_default));
}

inline tpsa
newt (const tpsa_t &a)
{
  return tpsa(mad_tpsa_new(&a, mad_tpsa_default));
}

// --- getter ---

inline tpsa_t&
ref (const tpsa &a)
{
  return *a.get();
}

// --- unary ---

inline tpsa
operator- (const tpsa_t &a) {
  tpsa c(newt(a));
  mad_tpsa_scl(&a, -1, c.get());
  return c;
}

inline tpsa operator-(const tpsa &a) { return -ref(a); }

// --- add ---

inline tpsa
operator+ (const tpsa_t &a, const tpsa_t &b) {
  tpsa c(newt(a));
  mad_tpsa_add(&a, &b, c.get());
  return c;
}

inline tpsa
operator+ (const tpsa_t &a, num_t b) {
  tpsa c(newt(a));
  mad_tpsa_copy(&a, c.get());
  mad_tpsa_set0(c.get(), 1, b);
  return c;
}

inline tpsa operator+(const tpsa   &a, const tpsa_t &b) { return ref(a) +     b ; }
inline tpsa operator+(const tpsa_t &a, const tpsa   &b) { return     a  + ref(b); }
inline tpsa operator+(const tpsa   &a, const tpsa   &b) { return ref(a) + ref(b); }
inline tpsa operator+(const tpsa   &a,        num_t  b) { return ref(a) +     b ; }
inline tpsa operator+(       num_t  a, const tpsa   &b) { return ref(b) +     a ; }
inline tpsa operator+(       num_t  a, const tpsa_t &b) { return     b  +     a ; }

// --- sub ---

inline tpsa
operator- (const tpsa_t &a, const tpsa_t &b) {
  tpsa c(newt(a));
  mad_tpsa_sub(&a, &b, c.get());
  return c;
}

inline tpsa
operator- (const tpsa_t &a, num_t b) {
  tpsa c(newt(a));
  mad_tpsa_copy(&a, c.get());
  mad_tpsa_set0(c.get(), 1, -b);
  return c;
}

inline tpsa
operator- (num_t a, const tpsa_t &b) {
  tpsa c(newt(b));
  mad_tpsa_scl(&b, -1, c.get());
  mad_tpsa_set0(c.get(), 1, a);
  return c;
}

inline tpsa operator-(const tpsa_t &a, const tpsa   &b) { return     a  - ref(b); }
inline tpsa operator-(const tpsa   &a, const tpsa_t &b) { return ref(a) -     b ; }
inline tpsa operator-(const tpsa   &a, const tpsa   &b) { return ref(a) - ref(b); }
inline tpsa operator-(const tpsa   &a,        num_t  b) { return ref(a) -     b ; }
inline tpsa operator-(       num_t  a, const tpsa   &b) { return     a  - ref(b); }

// --- mul ---

inline tpsa
operator* (const tpsa_t &a, const tpsa_t &b) {
  tpsa c(newt(a));
  mad_tpsa_mul(&a, &b, c.get());
  return c;
}

inline tpsa
operator* (const tpsa_t &a, num_t b) {
  tpsa c(newt(a));
  mad_tpsa_scl(&a, b, c.get());
  return c;
}

inline tpsa operator*(const tpsa_t &a, const tpsa   &b) { return     a  * ref(b); }
inline tpsa operator*(const tpsa   &a, const tpsa_t &b) { return ref(a) *     b ; }
inline tpsa operator*(const tpsa   &a, const tpsa   &b) { return ref(a) * ref(b); }
inline tpsa operator*(const tpsa   &a,        num_t  b) { return ref(a) *     b ; }
inline tpsa operator*(       num_t  a, const tpsa   &b) { return ref(b) *     a ; }
inline tpsa operator*(       num_t  a, const tpsa_t &b) { return     b  *     a ; }

// --- div ---

inline tpsa
operator/ (const tpsa_t &a, const tpsa_t &b) {
  tpsa c(newt(a));
  mad_tpsa_div(&a, &b, c.get());
  return c;
}

inline tpsa
operator/ (const tpsa_t &a, num_t b) {
  tpsa c(newt(a));
  mad_tpsa_scl(&a, 1/b, c.get());
  return c;
}

inline tpsa
operator/ (num_t a, const tpsa_t &b) {
  tpsa c(newt(b));
  mad_tpsa_inv(&b, a, c.get());
  return c;
}

inline tpsa operator/(const tpsa_t &a, const tpsa   &b) { return     a  / ref(b); }
inline tpsa operator/(const tpsa   &a, const tpsa_t &b) { return ref(a) /     b ; }
inline tpsa operator/(const tpsa   &a, const tpsa   &b) { return ref(a) / ref(b); }
inline tpsa operator/(const tpsa   &a,        num_t  b) { return ref(a) /     b ; }
inline tpsa operator/(       num_t  a, const tpsa   &b) { return     a  / ref(b); }

// --- pow ---

inline tpsa
pow (const tpsa_t &a, const tpsa_t &b) {
  tpsa c(newt(a));
  mad_tpsa_pow(&a, &b, c.get());
  return c;
}

inline tpsa
pow (const tpsa_t &a, int b) {
  tpsa c(newt(a));
  mad_tpsa_powi(&a, b, c.get());
  return c;
}

inline tpsa
pow (const tpsa_t &a, num_t b) {
  tpsa c(newt(a));
  mad_tpsa_pown(&a, b, c.get());
  return c;
}

inline tpsa
pow (num_t a, const tpsa_t &b) {
  tpsa c(newt(b));
  mad_tpsa_scl(&b, std::log(a), c.get());
  mad_tpsa_exp(c.get(), c.get());
  return c;
}

inline tpsa pow(const tpsa   &a, const tpsa_t &b) { return pow(ref(a),     b ); }
inline tpsa pow(const tpsa_t &a, const tpsa   &b) { return pow(    a , ref(b)); }
inline tpsa pow(const tpsa   &a, const tpsa   &b) { return pow(ref(a), ref(b)); }
inline tpsa pow(const tpsa   &a,        int    b) { return pow(ref(a),     b ); }
inline tpsa pow(const tpsa   &a,        num_t  b) { return pow(ref(a),     b ); }
inline tpsa pow(       num_t  a, const tpsa   &b) { return pow(    a , ref(b)); }

// warning: the operator ^ hasn't the expected precedence and associativity...

inline tpsa operator^(const tpsa_t &a, const tpsa_t &b) { return pow(a,b); }
inline tpsa operator^(const tpsa   &a, const tpsa_t &b) { return pow(a,b); }
inline tpsa operator^(const tpsa_t &a, const tpsa   &b) { return pow(a,b); }
inline tpsa operator^(const tpsa   &a, const tpsa   &b) { return pow(a,b); }
inline tpsa operator^(const tpsa_t &a,        int    b) { return pow(a,b); }
inline tpsa operator^(const tpsa_t &a,        num_t  b) { return pow(a,b); }
inline tpsa operator^(const tpsa   &a,        int    b) { return pow(a,b); }
inline tpsa operator^(const tpsa   &a,        num_t  b) { return pow(a,b); }
inline tpsa operator^(       num_t  a, const tpsa   &b) { return pow(a,b); }
inline tpsa operator^(       num_t  a, const tpsa_t &b) { return pow(a,b); }

} // mad

// --- functions ---

namespace mad {

// --- unary ---

#define FUN(F) \
\
inline tpsa F (const tpsa_t &a) { \
  tpsa c(newt(a)); \
  mad_tpsa_ ## F (&a, c.get()); \
  return c; \
} \
\
inline tpsa F (const tpsa &a) { return F(ref(a)); }

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

// --- end --------------------------------------------------------------------o

#endif // MAD_TPSA_HPP
