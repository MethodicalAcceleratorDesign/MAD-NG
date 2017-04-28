#ifndef TVAL_H
#define TVAL_H

/*
 o-----------------------------------------------------------------------------o
 |
 | Tagged values module (single header only)
 |
 | Authors: L. Deniau, laurent.deniau at cern.ch
 | Contrib: -
 |
 o-----------------------------------------------------------------------------o
 | You can redistribute this file and/or modify it under the terms of the GNU
 | General Public License GPLv3 (or later), as published by the Free Software
 | Foundation. This file is distributed in the hope that it will be useful, but
 | WITHOUT ANY WARRANTY OF ANY KIND. See http://gnu.org/licenses for details.
 o-----------------------------------------------------------------------------o

  - Small and simple module to manipulate Tagged Values (i.e. NaN): API for
    identification, boxing and unboxing TV, plus few other useful functions.

  - Provides support for few types (can be changed or extended easily):
    + nil         (special)
    + logical     (boolean, false and true)
    + integer     (46 bit signed integer)
    + number      (double precision floating point)
    + function    (generic function pointer)
    + pointer     (generic void pointer)
    + string      (constant, '\0' terminated)
    + array       (array of TV, opt. nil terminated)
    + object      (user defined objects with common header)
    + reference   (reference to TV)
    + instruction (46 bit encoded instruction)

 o-----------------------------------------------------------------------------o
*/

#include <stdbool.h>
#include <stdint.h>

// function name mangling
#ifndef FN
#define FN(f) f
#define FN_UNDEF
#endif

// --- interface --------------------------------------------------------------o

// types
typedef bool       log_t;
typedef int64_t    i64_t;
typedef uint64_t   u64_t;
typedef double     num_t;
typedef void      (fun_t)(void);
typedef void       ptr_t;
typedef const char str_t;
typedef struct arr arr_t;
typedef struct obj obj_t;
typedef union tval val_t;

// identity
static bool FN(tvisnan)  (val_t); // numerical NaN
static bool FN(tvisnul)  (val_t); // 0, +0.0, NULL
static bool FN(tvisnil)  (val_t);
static bool FN(tvisfalse)(val_t);
static bool FN(tvistrue) (val_t);

static bool FN(tvislog) (val_t);
static bool FN(tvisint) (val_t); // 46 bit signed int
static bool FN(tvisnum) (val_t); // tv num
static bool FN(tvisins) (val_t); // 46 bit instruction
static bool FN(tvisfun) (val_t); // tv fun
static bool FN(tvisptr) (val_t); // tv ptr (!log !int !num !fun)
static bool FN(tvisstr) (val_t);
static bool FN(tvisarr) (val_t);
static bool FN(tvisobj) (val_t);
static bool FN(tvisref) (val_t); // tv ref (val_t*)
static bool FN(tvisval) (val_t); // tv val (log int)

// boxing
static val_t FN(tvnan)  (void);
static val_t FN(tvnul)  (void);
static val_t FN(tvnil)  (void);
static val_t FN(tvfalse)(void);
static val_t FN(tvtrue) (void);

static val_t FN(tvlog) (log_t);
static val_t FN(tvint) (i64_t); // 46 bit signed int
static val_t FN(tvnum) (num_t);
static val_t FN(tvins) (u64_t); // 46 bit instruction
static val_t FN(tvfun) (fun_t*);
static val_t FN(tvptr) (ptr_t*);
static val_t FN(tvstr) (str_t*);
static val_t FN(tvarr) (arr_t*);
static val_t FN(tvobj) (obj_t*);
static val_t FN(tvref) (val_t*);

// unboxing
static log_t  FN(logtv) (val_t);
static i64_t  FN(inttv) (val_t); // 46 bit signed int
static num_t  FN(numtv) (val_t);
static u64_t  FN(instv) (val_t); // 46 bit instruction
static fun_t* FN(funtv) (val_t);
static ptr_t* FN(ptrtv) (val_t);
static str_t* FN(strtv) (val_t);
static arr_t* FN(arrtv) (val_t);
static obj_t* FN(objtv) (val_t);
static val_t* FN(reftv) (val_t);

// dereference
static val_t  FN(tvget) (val_t); // tv ref resolution

// typeid, typename
static int    FN(typtv) (val_t); // type id
static str_t* FN(namtv) (val_t); // type name

// debugging, display
static ptr_t* FN(hextv) (val_t); // ptr representation
static u64_t  FN(bittv) (val_t); // bit representation

// typeid
enum { TVNUM, TVNIL, TVLOG, TVINT, TVINS, TVFUN,
       TVPTR, TVSTR, TVARR, TVOBJ, TVXXX, TVREF };

// ----------------------------------------------------------------------------o
// --- implementation ---------------------------------------------------------o
// ----------------------------------------------------------------------------o

// tags (!= typeid)
#define DEF_NINF     (      0x7FF0ULL << DEF_TSHT )  // num canonical Inf
#define DEF_NNAN     (      0x7FF8ULL << DEF_TSHT )  // num canonical NaN
#define DEF_TNAN     (      0xFFF8ULL << DEF_TSHT )  // tag canonical NaN
#define DEF_TVAL     (DEF_TNAN | 0ULL << DEF_TSHT )  // value
#define DEF_TFUN     (DEF_TNAN | 1ULL << DEF_TSHT )  // function
#define DEF_TPTR     (DEF_TNAN | 2ULL << DEF_TSHT )  // pointer
#define DEF_TSTR     (DEF_TNAN | 3ULL << DEF_TSHT )  // pointer
#define DEF_TARR     (DEF_TNAN | 4ULL << DEF_TSHT )  // pointer
#define DEF_TOBJ     (DEF_TNAN | 5ULL << DEF_TSHT )  // pointer
#define DEF_TXXX     (DEF_TNAN | 6ULL << DEF_TSHT )  // pointer (unused)
#define DEF_TREF     (DEF_TNAN | 7ULL << DEF_TSHT )  // pointer
#define DEF_TMSK     (DEF_TNAN | 7ULL << DEF_TSHT )  // mask  (pointer mask)
#define DEF_TSHT     (48)                            // shift (pointer size)

#define DEF_TVNIL    (DEF_TVAL | 0ULL << DEF_TVSHT)  // value
#define DEF_TVLOG    (DEF_TVAL | 1ULL << DEF_TVSHT)  // value
#define DEF_TVINT    (DEF_TVAL | 2ULL << DEF_TVSHT)  // value
#define DEF_TVINS    (DEF_TVAL | 3ULL << DEF_TVSHT)  // value (instruction)
#define DEF_TVMSK    (DEF_TMSK | 3ULL << DEF_TVSHT)  // mask  (value mask)
#define DEF_TVSHT    (DEF_TSHT-2)                    // shift (value size)

#define DEF_TVNUL    (0)                             // value
#define DEF_TVFALSE  (DEF_TVLOG | 0)                 // value
#define DEF_TVTRUE   (DEF_TVLOG | 1)                 // value

// checks
#define TVISNAN(v)   ((v).__u == DEF_NNAN   )
#define TVISNUL(v)   ((v).__u == DEF_TVNUL  )
#define TVISNIL(v)   ((v).__u == DEF_TVNIL  )
#define TVISFALSE(v) ((v).__u == DEF_TVFALSE)
#define TVISTRUE(v)  ((v).__u == DEF_TVTRUE )

#define TVISLOG(v)   (((v).__u & DEF_TVMSK) == DEF_TVLOG)
#define TVISINT(v)   (((v).__u & DEF_TVMSK) == DEF_TVINT)
#define TVISINS(v)   (((v).__u & DEF_TVMSK) == DEF_TVINS)
#define TVISNUM(v)   (((v).__u & DEF_TNAN ) != DEF_TNAN )
#define TVISFUN(v)   (((v).__u & DEF_TMSK ) == DEF_TFUN )
#define TVISPTR(v)   (((v).__u & DEF_TMSK ) >= DEF_TPTR )
#define TVISSTR(v)   (((v).__u & DEF_TMSK ) == DEF_TSTR )
#define TVISARR(v)   (((v).__u & DEF_TMSK ) == DEF_TARR )
#define TVISOBJ(v)   (((v).__u & DEF_TMSK ) == DEF_TOBJ )
#define TVISXXX(v)   (((v).__u & DEF_TMSK ) == DEF_TXXX )
#define TVISREF(v)   (((v).__u & DEF_TMSK ) == DEF_TREF )
#define TVISVAL(v)   (((v).__u & DEF_TMSK ) == DEF_TVAL )
#define TVISHEX(v)   (((v).__u & DEF_TMSK ) >= DEF_TFUN )

// taggeg value representation
union tval {
  u64_t  __u;
  i64_t  __i;
  num_t  __d;
  fun_t* __f;
  ptr_t* __p;
  str_t* __s;
  arr_t* __a;
  obj_t* __o;
  val_t* __r;
};

// sanity checks
enum { Unsupported_architecture_IEEE_754_required =
       1/(sizeof(val_t) == sizeof(num_t)) };

// local inlined functions
#undef  inline
#define inline static inline

// --- introspection ----------------------------------------------------------o

inline bool FN(tvisnan)  (val_t v) { return TVISNAN(v);   }
inline bool FN(tvisnul)  (val_t v) { return TVISNUL(v);   }
inline bool FN(tvisnil)  (val_t v) { return TVISNIL(v);   }
inline bool FN(tvisfalse)(val_t v) { return TVISFALSE(v); }
inline bool FN(tvistrue) (val_t v) { return TVISTRUE(v);  }

inline bool FN(tvislog)  (val_t v) { return TVISLOG(v); }
inline bool FN(tvisint)  (val_t v) { return TVISINT(v); }
inline bool FN(tvisnum)  (val_t v) { return TVISNUM(v); }
inline bool FN(tvisins)  (val_t v) { return TVISINS(v); }
inline bool FN(tvisfun)  (val_t v) { return TVISFUN(v); }
inline bool FN(tvisptr)  (val_t v) { return TVISPTR(v); }
inline bool FN(tvisstr)  (val_t v) { return TVISSTR(v); }
inline bool FN(tvisarr)  (val_t v) { return TVISARR(v); }
inline bool FN(tvisobj)  (val_t v) { return TVISOBJ(v); }
inline bool FN(tvisref)  (val_t v) { return TVISREF(v); }
inline bool FN(tvisval)  (val_t v) { return TVISVAL(v); }

// --- boxing -----------------------------------------------------------------o

inline val_t FN(tvnan)  (void) { return (val_t){ DEF_NNAN    }; }
inline val_t FN(tvnul)  (void) { return (val_t){ DEF_TVNUL   }; }
inline val_t FN(tvnil)  (void) { return (val_t){ DEF_TVNIL   }; }
inline val_t FN(tvfalse)(void) { return (val_t){ DEF_TVFALSE }; }
inline val_t FN(tvtrue) (void) { return (val_t){ DEF_TVTRUE  }; }

inline val_t FN(tvlog) (log_t l)
{
  return (val_t){ DEF_TVLOG | l };
}

inline val_t FN(tvint) (i64_t i)
{
  if (i < 0)
    return (val_t){ (-i & ~DEF_TVMSK >> 1) | DEF_TVINT | 1ULL << (DEF_TVSHT-1)};
  else
    return (val_t){ ( i & ~DEF_TVMSK >> 1) | DEF_TVINT };
}

inline val_t FN(tvins) (u64_t i)
{
  return (val_t){ ( i & ~DEF_TVMSK) | DEF_TVINS };
}

inline val_t FN(tvnum) (num_t d)
{
  return d == d ? (val_t){ .__d = d } : (val_t){ DEF_NNAN };
}

inline val_t FN(tvfun) (fun_t* f)
{
  val_t v = (val_t){ .__f = f };
  v.__u = (v.__u & ~DEF_TMSK) | DEF_TFUN;
  return v;
}

inline val_t FN(tvptr) (ptr_t* p)
{
  val_t v = (val_t){ .__p = p };
  v.__u = (v.__u & ~DEF_TMSK) | DEF_TPTR;
  return v;
}

inline val_t FN(tvstr) (str_t* s)
{
  val_t v = (val_t){ .__s = s };
  v.__u = (v.__u & ~DEF_TMSK) | DEF_TSTR;
  return v;
}

inline val_t FN(tvarr) (arr_t* a)
{
  val_t v = (val_t){ .__a = a };
  v.__u = (v.__u & ~DEF_TMSK) | DEF_TARR;
  return v;
}

inline val_t FN(tvobj) (obj_t* o)
{
  val_t v = (val_t){ .__o = o };
  v.__u = (v.__u & ~DEF_TMSK) | DEF_TOBJ;
  return v;
}

inline val_t FN(tvref) (val_t* r)
{
  val_t v = (val_t){ .__r = r };
  v.__u = (v.__u & ~DEF_TMSK) | DEF_TREF;
  return v;
}

// --- unboxing ---------------------------------------------------------------o

#include <assert.h>

inline log_t FN(logtv) (val_t v)
{
  assert( TVISLOG(v) );
  return v.__u & 1;
}

inline i64_t FN(inttv) (val_t v)
{
  assert( TVISINT(v) );
  if ( v.__u & 1ULL << (DEF_TVSHT-1) )
    return -(v.__i & (~DEF_TVMSK >> 1));
  else
    return   v.__i & (~DEF_TVMSK >> 1) ;
}

inline num_t FN(numtv) (val_t v)
{
  return v.__d;
}

inline u64_t FN(instv) (val_t v)
{
  return TVISINS(v) ? v.__u & ~DEF_TVMSK : 0;
}

inline fun_t* FN(funtv) (val_t v)
{
  return TVISFUN(v) ? (val_t){ v.__u & ~DEF_TMSK } .__f : 0;
}

inline ptr_t* FN(ptrtv) (val_t v)
{
  return TVISPTR(v) ? (val_t){ v.__u & ~DEF_TMSK } .__p : 0;
}

inline str_t* FN(strtv) (val_t v)
{
  return TVISSTR(v) ? (val_t){ v.__u & ~DEF_TMSK } .__s : 0;
}

inline arr_t* FN(arrtv) (val_t v)
{
  return TVISARR(v) ? (val_t){ v.__u & ~DEF_TMSK } .__a : 0;
}

inline obj_t* FN(objtv) (val_t v)
{
  return TVISOBJ(v) ? (val_t){ v.__u & ~DEF_TMSK } .__o : 0;
}

inline val_t* FN(reftv) (val_t v)
{
  return TVISREF(v) ? (val_t){ v.__u & ~DEF_TMSK } .__r : 0;
}

// --- dereference ------------------------------------------------------------o

inline val_t FN(tvget) (val_t v)
{
  while ( TVISREF(v) ) v = *FN(reftv)(v);
  return v;
}

// --- typeid, typename -------------------------------------------------------o

inline int FN(typtv) (val_t v)
{
  return TVISNUM(v) ? 0 :                            // num: not a tv
         TVISHEX(v) ? ((v.__u >> DEF_TSHT ) & 7)+4 : // decode pointer
                      ((v.__u >> DEF_TVSHT) & 3)+1 ; // decode value
}

inline str_t* FN(namtv) (val_t v)
{
  int id = FN(typtv)(v);
  assert(0 <= id && id <= 11);
  return (str_t*[12]){
    "tvnum", "tvnil", "tvlog", "tvint", "tvins", "tvfun",
    "tvptr", "tvstr", "tvarr", "tvobj", "tvxxx", "tvref" } [id];
}

// --- debugging, display -----------------------------------------------------o

inline ptr_t* FN(hextv) (val_t v)
{
  return (val_t){ v.__u & ~DEF_TMSK } .__p;
}

inline u64_t FN(bittv) (val_t v)
{
  return v.__u;
}

// --- cleanup ----------------------------------------------------------------o

#undef inline

#ifdef FN_UNDEF
#undef FN_UNDEF
#undef FN
#endif

// --- end --------------------------------------------------------------------o

#endif // TVAL_H
