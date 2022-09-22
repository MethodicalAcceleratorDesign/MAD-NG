.. index::
   Pseudo-random number generator
   PRNG

**************
Random Numbers
**************

The module :mod:`gmath` provides few Pseudo-Random Number Generators (PRNGs).The defaut implementation is the *Xoshiro256\*\** (XOR/shift/rotate) variant of the `XorShift <https://en.wikipedia.org/wiki/Xorshift>`_ PRNG familly [XORSHFT03]_, an all-purpose, rock-solid generator with a period of :math:`2^{256}-1` that supports long jumps of period :math:`2^{128}`. This PRNG is also the default implementation of recent versions of Lua (not LuaJIT, see below) and GFortran. See https://prng.di.unimi.it for details about Xoshiro/Xoroshiro PRNGs.

The module :mod:`math` of LuaJIT provides an implementation of the *Tausworthe* PRNG [TAUSWTH96]_, which has a period of :math:`2^{223}` but doesn't support long jumps, and hence uses a single global PRNG.

The module :mod:`gmath` also provides an implementation of the simple global PRNG of MAD-X for comparison.

It's worth mentionning that none of these PRNG are cryptographically secure generators, they are nevertheless superior to the commonly used *Mersenne Twister* PRNG [MERTWIS98]_, with the exception of the MAD-X PRNG.

All PRNG *functions* (except constructors) are wrappers around PRNG *methods* with the same name, and expect an optional PRNG :obj:`prng_` as first parameter. If this optional PRNG :obj:`prng_` is omitted, i.e. not provided, these functions will use the current global PRNG by default.

Contructors
===========

.. function:: randnew ()

   Return a new Xoshiro256\*\* PRNG with a period of :math:`2^{128}` that is garuanteed to not overlapp with any other Xoshiro256\*\* PRNGs, unless it is initialized with a seed.

.. function:: xrandnew ()

   Return a new MAD-X PRNG initialized with default seed 123456789. Hence, all new MAD-X PRNG will generate the same sequence until they are initialized with a user-defined seed.

Functions
=========

.. function:: randset (prng_)

   Set the current global PRNG to :obj:`prng` (if provided) and return the previous global PRNG.

.. function:: is_randgen(a)

   Return :const:`true` if :var:`a` is a PRNG, :const:`false` otherwise. This function is also available from the module :mod:`MAD.typeid`.

.. function:: is_xrandgen(a)

   Return :const:`true` if :var:`a` is a MAD-X PRNG, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

All methods are also provided as functions from the module :mod:`MAD.gmath` for convenience. If the PRNG is not provided, the current global PRNG is used instead.

.. function:: prng:randseed (seed)
              randseed ([prng_,] seed)

   Set the seed of the PRNG :obj:`prng` to :var:`seed`.

.. function:: prng:rand ()
              rand (prng_)
              
   Return a new pseudo-random number in the range :const:`[0, 1)` from the PRNG :obj:`prng`.

.. function:: prng:randi ()
              randi (prng_)

   Return a new pseudo-random number in the range of a :type:`u64_t` from the PRNG :obj:`prng` (:type:`u32_t` for the MAD-X PRNG), see C API below for details.

.. function:: prng:randn ()
              randn (prng_)

   Return a new pseudo-random gaussian number in the range :const:`[-inf, +inf]` from the PRNG :obj:`prng` by using the Box-Muller transformation (Marsaglia's polar form) to a peuso-random number in the range :const:`[0, 1)`.

.. function:: prng:randtn (cut_)
              randtn ([prng_,] cut_)

   Return a new truncated pseudo-random gaussian number in the range :const:`[-cut_, +cut_]` from the PRNG :obj:`prng` by using iteratively the method :func:`prng:randn()`. This simple algorithm is actually used for compatibility with MAD-X.
   Default: :expr:`cut_ = +inf`.

.. function:: prng:randp (lmb_)
              randp ([prng_,] lmb_)

   Return a new pseudo-random poisson number in the range :const:`[0, +inf]` from the PRNG :obj:`prng` with parameter :math:`\lambda > 0` by using the *inverse transform sampling* method on peuso-random numbers.
   Default: :expr:`lmb_ = 1`.

Iterators
=========

.. function:: ipairs(prng)
   :noindex:

   Return an :type:`ipairs` iterator suitable for generic :const:`for` loops. The generated values are those returned by :func:`prng:rand()`. 

C API
=====

.. c:type:: prng_state_t
            xrng_state_t

   The Xoshiro256\*\* and the MAD-X PRNG types.

.. c:function:: num_t mad_num_rand (prng_state_t*)

   Return a pseudo-random double precision float in the range :const:`[0, 1)`. 

.. c:function:: u64_t mad_num_randi (prng_state_t*)

   Return a pseudo-random 64 bit unsigned integer in the range :const:`[0, ULLONG_MAX]`.

.. c:function:: void mad_num_randseed (prng_state_t*, num_t seed)

   Set the seed of the PRNG.

.. c:function:: void mad_num_randjump (prng_state_t*)

   Apply a jump to the PRNG as if :math:`2^{128}` pseudo-random numbers were generated. Hence PRNGs with different number of jumps will never overlap. This function is applied to new PRNGs with an incremental number of jumps. 

.. c:function:: num_t mad_num_xrand (xrng_state_t*)

   Return a pseudo-random double precision float in the range :const:`[0, 1)` from the MAD-X PRNG.

.. c:function:: u32_t mad_num_xrandi (xrng_state_t*)

   Return a pseudo-random 32 bit unsigned integer in the range :const:`[0, UINT_MAX]` from the MAD-X PRNG.

.. c:function:: void mad_num_xrandseed (xrng_state_t*, u32_t seed)

   Set the seed of the MAD-X PRNG.

References
==========

.. [XORSHFT03] G. Marsaglia, *"Xorshift RNGs"*, Journal of Statistical Software, 8 (14), July 2003. doi:10.18637/jss.v008.i14.

.. [TAUSWTH96] P. L’Ecuyer, *“Maximally Equidistributed Combined Tausworthe Generators”*, Mathematics of Computation, 65 (213), 1996, p203–213.

.. [MERTWIS98] M. Matsumoto and T. Nishimura, *“Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator”*. ACM Trans. on Modeling and Comp. Simulation, 8 (1), Jan. 1998, p3–30.
