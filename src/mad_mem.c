/*
 o-----------------------------------------------------------------------------o
 |
 | Memory module implementation
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
*/

/* Note: This module adds a thread-safe front-end to the global C allocator to
   speed-up by x10+ frequent interleaved malloc and free of "small" objects,
   like e.g. in expressions evaluations. See unit test in main() below. */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include "mad_mem.h"

#define MAD_MEM_STD   0 // 1 -> use standard C allocator only.
#define MAD_MEM_CLR   0 // 1 -> replace malloc by calloc & clear pooled chunk.
#define MAD_MEM_UTEST 0 // 1 -> run standalone unit tests in main().
#define DBGMEM(P)       // P // uncomment for verbose debugging output

// --- standard allocator -----------------------------------------------------o

#if MAD_MEM_CLR == 1    // clear malloc'ed chunk
#define malloc(sz) calloc(1,sz)
#endif

#if MAD_MEM_STD == 1    // module disabled

// allocators
void*  (mad_malloc ) (size_t sz)               { return malloc (sz);        }
void*  (mad_calloc ) (size_t cnt , size_t sz ) { return calloc (cnt, sz);   }
void*  (mad_realloc) (void*  ptr_, size_t sz_) { return realloc(ptr_, sz_); }
void   (mad_free   ) (void*  ptr_)             { free(ptr_); }

// utils
size_t  mad_mcached  (void)                    { return 0; }
size_t  mad_mcollect (void)                    { return 0; }
void    mad_mdump    (FILE* fp)                { (void)fp; }

// --- generational allocator -------------------------------------------------o

#else

// --- types & constants ------------------------------------------------------o

// constant sizes
enum {
  stp_slot = sizeof(double), // memblk unit of size increment, don't touch
  max_slot = 8192,    // max 2^16-2, default 8192 slots -> max obj size is 64KB
  max_mblk = 2048,    // max 2^16-2, default 2048 slots
  max_mkch = 2097152, // max 2^32-1, default 2097152 stp_slot -> 16MB
};

// memory block
struct memblk {
  uint16_t slot; // slot index (i.e. rounded size).
  uint16_t next; // index of next memblk (i.e. linked list).
  uint32_t mark; // memory boundary marker
  union { // alignment of data
    ptrdiff_t s, *sp;
    size_t    u, *up;
    double    d, *dp;
    void         *vp;
    void        (*fp)(void);
  } data[1];
};

// memory pool
struct pool {
  uint32_t mkch;            // amount of cached memory in stp_slot unit
  uint16_t free;            // index of first free slot in mblk
  uint16_t _dum;            // for alignment
  uint16_t slot[max_slot];  // index+1 in mblk of first mbp for this size
  union {
    size_t         nxt;     // index in mblk of next free slot
    struct memblk *mbp;
  }        mblk[max_mblk];
  char     str[128];        // for debug, see pdump()
};

// macros
#define BASE(ptr) ((void*)((char*)(ptr)-stp_slot)) // ptr -> mbp
#define SIZE(idx) (((size_t)(idx)+2)*stp_slot)     // sizeof(*mbp)
#define CACHED(p) ((size_t)(p)->mkch*stp_slot)
#define IDXMAX    0xFFFF
#define SLTMAX    0xFFFFFFFF
#define MARK      0xACCEDEAD

// static sanity checks
enum {
  static_assert__max_slot_not_a_power_of_2 = 1/!(max_slot & (max_slot-1)), // not a real constraint, could be remove
  static_assert__max_mblk_not_a_power_of_2 = 1/!(max_mblk & (max_mblk-1)), // not a real constraint, could be remove
  static_assert__max_mkch_not_a_power_of_2 = 1/!(max_mkch & (max_mkch-1)), // not a real constraint, could be remove
  static_assert__max_slot_not_lt_IDXMAX    = 1/ (max_slot < IDXMAX),
  static_assert__max_mblk_not_lt_IDXMAX    = 1/ (max_mblk < IDXMAX),
  static_assert__max_mkch_not_lt_SLTMAX    = 1/ (max_mkch < SLTMAX),

  static_assert__stp_slot_not_a_power_of_2 = 1/!(stp_slot & (stp_slot-1)),
  static_assert__stp_slot_neq_sizeof       = 1/ (stp_slot == sizeof(double)),
  static_assert__stp_slot_neq_offsetof     = 1/ (stp_slot == offsetof(struct memblk,data)), // very important...
};

// --- locals -----------------------------------------------------------------o

static struct pool pool = {0,0,0,{0},{{0}},{0}};

#ifdef _OPENMP
#pragma omp threadprivate(pool)
#endif

static inline char*
pdump(struct memblk *mbp)
{
  struct pool *p = &pool;
  snprintf(p->str, sizeof(p->str), "%p {slot=%4d(%5td), next=%4d, mark=%x}%s",
         (void*)mbp, mbp->slot,
         mbp->slot == IDXMAX ? -1 : (ptrdiff_t)SIZE(mbp->slot),
         mbp->next-1, mbp->mark, mbp->mark == MARK ? "" : " (corrupted!)");
  return p->str;
}

// --- implementation ---------------------------------------------------------o

void*
(mad_malloc) (size_t size)
{
  struct pool *p = &pool;
  size_t idx = size ? (size-1) / stp_slot : 0;
  struct memblk *mbp; // mbp = mblk[slot[idx]-1]

  if (idx < max_slot && p->slot[idx]) {
    idx_t slt = p->slot[idx]-1;
    DBGMEM( printf("alloc: reuse mblk[[%d]=%d]", (int)idx, slt); )
    mbp = p->mblk[slt].mbp, p->mblk[slt].nxt = p->free;
    p->free = p->slot[idx], p->slot[idx] = mbp->next;
    p->mkch -= idx+2;
#if MAD_MEM_CLR
    memset(mbp->data, 0, size);
#endif
  } else {
    DBGMEM( printf("alloc: malloc(%2zu)", size); )
    mbp = malloc(SIZE(idx));
    if (!mbp) return warn("cannot allocate %zu bytes", size), NULL;
    mbp->slot = idx < max_slot ? idx : IDXMAX;
    mbp->mark = MARK;
    ensure((size_t)mbp > IDXMAX, "unexpected very low address"); // see collect
  }

  DBGMEM( printf(" at %s\n", pdump(mbp)); )
  return mbp->data;
}

void
(mad_free) (void *ptr)
{
  if (!ptr) return;

  struct memblk *mbp = BASE(ptr);
  idx_t idx = mbp->slot;

  ensure(mbp->mark == MARK, "invalid or corrupted allocated memory");

  if (idx == IDXMAX) {
    DBGMEM( printf("free : free mblk at %s\n", pdump(mbp)); )
    free(mbp); return;
  }

  struct pool *p = &pool;
  if (!p->free || p->mkch >= max_mkch) // no free slot (or init) or max cache
    mad_mcollect();

  idx_t slt = p->free-1;
  DBGMEM( printf("free : cache mblk[[%d]=%d]", idx, slt); )
  mbp->next = p->slot[idx], p->slot[idx] = p->free;   // update slot
  p->free = p->mblk[slt].nxt, p->mblk[slt].mbp = mbp; // store  mblk
  p->mkch += idx+2;

  DBGMEM( printf(" at %s\n", pdump(mbp)); )
}

void*
(mad_calloc) (size_t ecount, size_t esize)
{
  size_t size = ecount * esize;
  void  *ptr  = (mad_malloc)(size);
  return ptr ? memset(ptr, 0, size) : NULL;
}

void*
(mad_realloc) (void *ptr, size_t size)
{
  if (!size) return (mad_free)(ptr), NULL;
  if (!ptr ) return (mad_malloc)(size);

  DBGMEM( printf("alloc: realloc(%2zu)", size); )
  struct memblk *mbp = BASE(ptr);

  ensure(mbp->mark == MARK, "invalid or corrupted allocated memory");

  size_t idx = (size-1) / stp_slot;
  mbp = realloc(mbp, SIZE(idx));
  if (!mbp) return warn("cannot reallocate %zu bytes", size), NULL;

#if MAD_MEM_CLR
  if (mbp->slot < idx && idx < max_slot) {
    size_t off = SIZE(mbp->slot)-SIZE(idx);
    memset((char*)mbp->data+off, 0, size-off);
  }
#endif

  mbp->slot = idx < max_slot ? idx : IDXMAX;
  DBGMEM( printf(" at %s\n", pdump(mbp)); )
  return mbp->data;
}

// -- utils

size_t
mad_mcached (void)
{
  struct pool *p = &pool;
  size_t cached = CACHED(p), ccached = 0;

  for (idx_t i=0; i < max_mblk; i++)
    if (p->mblk[i].nxt > IDXMAX) // ptr
      ccached += SIZE(p->mblk[i].mbp->slot);

  ensure(ccached == cached, "corrupted cache %zu != %zu bytes", ccached, cached);

  return cached;
}

size_t
mad_mcollect (void)
{
  struct pool *p = &pool;
  size_t cached = CACHED(p);

  DBGMEM( printf("collect/clear/init cache\n"); )
  DBGMEM( printf("collecting %zu bytes\n", mad_mcached()); )

  p->mkch = 0;
  p->free = 1;

  for (idx_t i=0; i < max_slot; i++)
    p->slot[i] = 0;

  for (idx_t i=0; i < max_mblk; i++) {
    if (p->mblk[i].nxt > IDXMAX) // ptr
      free(p->mblk[i].mbp);
    p->mblk[i].nxt = i+2;
  }
  p->mblk[max_mblk-1].nxt = 0; // close linked list

  return cached;
}

// -- debug

void
mad_mdump (FILE *fp)
{
  struct pool *p = &pool;
  size_t cached = CACHED(p);

  if (!fp) fp = stdout;

  // init cache to avoid full dump of empty mblk
  if (!p->free && !cached) mad_mcollect();

  fprintf(fp, "mdump: %zu bytes\n", cached);

  // display content of slot[] when used, i.e. link to mblk[] + linked list.
  for (idx_t i=0; i < max_slot; i++) {
    idx_t slt = p->slot[i];
    if (slt) { // from slot (size) to memblk (object)
      fprintf(fp, "  slot[%4d] -> mblk[%d]", i, slt-1);
      struct memblk *mbp = p->mblk[slt-1].mbp;
      idx_t nxt = -1, j = 0;
      while (mbp->next) { // linked list
        nxt = mbp->next-1, mbp = p->mblk[nxt].mbp;
        if (++j < 8) fprintf(fp, "->[%d]", nxt);
      }
      if (j == 8) fprintf(fp,     "->[%d]\n", nxt); else
      if (j >  8) fprintf(fp, "->..->[%d]\n", nxt); else
                  fprintf(fp,           "\n");
    }
  }

  // display content of mblk[], i.e. object or not trivial link into mblk[]
  for (idx_t i=0; i < max_mblk; i++)
    if (p->mblk[i].nxt > IDXMAX)         // ptr
      fprintf(fp, "  mblk[%4d] -> %s\n", i, pdump(p->mblk[i].mbp));
    else if (i+1 == p->free)             // free
      fprintf(fp, "->mblk[%4d] -> [%d]\n", i, (int)p->mblk[i].nxt-1);
    else if (i+2 != (int)p->mblk[i].nxt) // idx
      fprintf(fp, "  mblk[%4d] -> [%d]\n", i, (int)p->mblk[i].nxt-1);
}

#endif // MAD_MEM_STD != 1

// --- unit test --------------------------------------------------------------o

#if MAD_MEM_UTEST == 1

/*
gcc -std=c99 -W -Wall -Wextra -pedantic -O3 -march=native -Wcast-align \
    -Wdisabled-optimization -Wpointer-arith -Wsign-compare -Wmissing-prototypes \
    -Wstrict-prototypes -Wunreachable-code -I. libgtpsa/mad_log.c mad_mem.c \
    -o mad_mem_ut
./mad_mem_ut > out
*/

int
main(void)
{
  enum { loop=1000000, mn = 50, n0 = 5, ni = 7, z = 6 };
  void *ptr[mn] = {0};
  size_t mcnt=0, fcnt=0;

  setvbuf(stdout, 0, _IONBF, 0);

  DBGMEM( union { size_t val; void* ptr; } p = {.ptr=ptr}; )
  DBGMEM( printf("ptr=%p, ptr as uint=0x%lx\n", p.ptr, p.val); )

  for (int k=0; k<loop; k++)
  for (int i=1, n=n0; n<mn; n+=ni, i++) {
    DBGMEM( printf("iteration %d: n=%d\n", i, n); )
    DBGMEM( printf("calling %d malloc\n", n); )
    for (int i=0; i<n; i++) ptr[i] = (++mcnt, mad_malloc((i%(13)+1)*z));
    DBGMEM( printf("status after %d malloc\n", n); mad_mdump(stdout); )
    DBGMEM( printf("calling %d free\n", n); )
    for (int i=0; i<n; i++) (++fcnt, mad_free(ptr[i]), ptr[i]=0);
    DBGMEM( printf("status after %d free\n", n); mad_mdump(stdout); )
  }
  mad_mdump(stdout);
  DBGMEM( mad_mdump(stdout); )
  printf("%9zu mallocs performed\n", mcnt);
  printf("%9zu frees   performed\n", fcnt);
  printf("%9zu bytes   cached\n", mad_mcached());
  mad_mcollect();
}

#endif

// --- end --------------------------------------------------------------------o

