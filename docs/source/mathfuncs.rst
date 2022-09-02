.. index::
   constants and functions

***********************
Constants and Functions
***********************

This chapter describes some constants and functions provided by the modules :mod:`MAD.constant` and :mod:`MAD.gmath`. The module :mod:`gmath` extends the standard LUA module :mod:`math` with *generic* functions working on any types that support the methods with the same names. For example, the code :func:`gmath.sin(a)` will call :func:`math.sin(a)` if :var:`a` is a :type:`number`, otherwise it will call the method :func:`a:sin()`, i.e. delegate the invocation to :obj:`a`. This is how MAD-NG handles several types like :type:`numbers`, :type:`complex` number and :type:`TPSA` within a single *polymorphic* code that expects scalar-like behavior.

Mathematical Constants
======================

This section describes basic mathematical constants uniquely defined as macros in the C header :file:`mad_cst.h` and available from C and MAD modules as floating point double precision variables. If these mathematical constants are already provided by the system libraries, they will be used instead of their local definitions.

===================  ======================  =========================  ======================
MAD constants        C macros                C constants                Values
===================  ======================  =========================  ======================
:const:`eps`         :c:macro:`DBL_EPSILON`  :const:`mad_cst_EPS`       Smallest representable increment near one
:const:`tiny`        :c:macro:`DBL_MIN`      :const:`mad_cst_TINY`      Smallest representable number
:const:`huge`        :c:macro:`DBL_MAX`      :const:`mad_cst_HUGE`      Largest representable number
:const:`inf`         :c:macro:`INFINITY`     :const:`mad_cst_INF`       Positive infinity, :math:`1/0`
:const:`nan`         :c:macro:`NAN`          :const:`mad_cst_NAN`       Canonical NaN [#]_, :math:`0/0`
:const:`e`           :c:macro:`M_E`          :const:`mad_cst_E`         :math:`e`
:const:`log2e`       :c:macro:`M_LOG2E`      :const:`mad_cst_LOG2E`     :math:`\log_2(e)`
:const:`log10e`      :c:macro:`M_LOG10E`     :const:`mad_cst_LOG10E`    :math:`\log_{10}(e)`
:const:`ln2`         :c:macro:`M_LN2`        :const:`mad_cst_LN2`       :math:`\ln(2)`
:const:`ln10`        :c:macro:`M_LN10`       :const:`mad_cst_LN10`      :math:`\ln(10)`
:const:`lnpi`        :c:macro:`M_LNPI`       :const:`mad_cst_LNPI`      :math:`\ln(\pi)`
:const:`pi`          :c:macro:`M_PI`         :const:`mad_cst_PI`        :math:`\pi`
:const:`twopi`       :c:macro:`M_2PI`        :const:`mad_cst_2PI`       :math:`2\pi`
:const:`pi_2`        :c:macro:`M_PI_2`       :const:`mad_cst_PI_2`      :math:`\pi/2`
:const:`pi_4`        :c:macro:`M_PI_4`       :const:`mad_cst_PI_4`      :math:`\pi/4`
:const:`one_pi`      :c:macro:`M_1_PI`       :const:`mad_cst_1_PI`      :math:`1/\pi`
:const:`two_pi`      :c:macro:`M_2_PI`       :const:`mad_cst_2_PI`      :math:`2/\pi`
:const:`sqrt2`       :c:macro:`M_SQRT2`      :const:`mad_cst_SQRT2`     :math:`\sqrt 2`
:const:`sqrt3`       :c:macro:`M_SQRT3`      :const:`mad_cst_SQRT3`     :math:`\sqrt 3`
:const:`sqrtpi`      :c:macro:`M_SQRTPI`     :const:`mad_cst_SQRTPI`    :math:`\sqrt{\pi}`
:const:`sqrt1_2`     :c:macro:`M_SQRT1_2`    :const:`mad_cst_SQRT1_2`   :math:`\sqrt{1/2}`
:const:`sqrt1_3`     :c:macro:`M_SQRT1_3`    :const:`mad_cst_SQRT1_3`   :math:`\sqrt{1/3}`
:const:`one_sqrtpi`  :c:macro:`M_1_SQRTPI`   :const:`mad_cst_1_SQRTPI`  :math:`1/\sqrt{\pi}`
:const:`two_sqrtpi`  :c:macro:`M_2_SQRTPI`   :const:`mad_cst_2_SQRTPI`  :math:`2/\sqrt{\pi}`
:const:`rad2deg`     :c:macro:`M_RAD2DEG`    :const:`mad_cst_RAD2DEG`   :math:`180/\pi`
:const:`deg2rad`     :c:macro:`M_DEG2RAD`    :const:`mad_cst_DEG2RAD`   :math:`\pi/180`
===================  ======================  =========================  ======================

.. index::
   mathematical constants

Physical Constants
==================

This section describes basic physical constants uniquely defined as macros in the C header :file:`mad_cst.h` and available from C and MAD modules as floating point double precision variables.

==================  =====================  =========================  ======================
MAD constants       C macros               C constants                Values
==================  =====================  =========================  ======================
:const:`minlen`     :c:macro:`P_MINLEN`    :const:`mad_cst_MINLEN`    Minimum length tolerance, default :math:`10^{-10}` in :unit:`[m]`
:const:`minang`     :c:macro:`P_MINANG`    :const:`mad_cst_MINANG`    Minimum angle tolerance, default :math:`10^{-10}` in :unit:`[1/m]`
:const:`minstr`     :c:macro:`P_MINSTR`    :const:`mad_cst_MINSTR`    Minimum strength tolerance, default :math:`10^{-10}` in :unit:`[rad]`
==================  =====================  =========================  ======================

The following table lists some physical constants from the `CODATA 2018 <https://physics.nist.gov/cuu/pdf/wall_2018.pdf>`_ sheet.

==================  =====================  =========================  ======================
MAD constants       C macros               C constants                Values
==================  =====================  =========================  ======================
:const:`clight`     :c:macro:`P_CLIGHT`    :const:`mad_cst_CLIGHT`    Speed of light, :math:`c` in :unit:`[m/s]`
:const:`mu0`        :c:macro:`P_MU0`       :const:`mad_cst_MU0`       Permeability of vacuum, :math:`\mu_0` in :unit:`[T.m/A]`
:const:`epsilon0`   :c:macro:`P_EPSILON0`  :const:`mad_cst_EPSILON0`  Permittivity of vacuum, :math:`\epsilon_0` in :unit:`[F/m]`
:const:`qelect`     :c:macro:`P_QELECT`    :const:`mad_cst_QELECT`    Elementary electric charge, :math:`e` in :unit:`[C]`
:const:`hbar`       :c:macro:`P_HBAR`      :const:`mad_cst_HBAR`      Reduced Plack's constant, :math:`\hbar` in :unit:`[GeV.s]`
:const:`amass`      :c:macro:`P_AMASS`     :const:`mad_cst_AMASS`     Unified atomic mass, :math:`m_u\,c^2` in :unit:`[GeV]`
:const:`emass`      :c:macro:`P_EMASS`     :const:`mad_cst_EMASS`     Electron mass, :math:`m_e\,c^2` in :unit:`[GeV]`
:const:`pmass`      :c:macro:`P_PMASS`     :const:`mad_cst_PMASS`     Proton mass, :math:`m_p\,c^2` in :unit:`[GeV]`
:const:`nmass`      :c:macro:`P_NMASS`     :const:`mad_cst_NMASS`     Neutron mass, :math:`m_n\,c^2` in :unit:`[GeV]`
:const:`mumass`     :c:macro:`P_MUMASS`    :const:`mad_cst_MUMASS`    Muon mass, :math:`m_{\mu}\,c^2` in :unit:`[GeV]`
:const:`deumass`    :c:macro:`P_DEUMASS`   :const:`mad_cst_DEUMASS`   Deuteron mass, :math:`m_d\,c^2` in :unit:`[GeV]`
:const:`eradius`    :c:macro:`P_ERADIUS`   :const:`mad_cst_ERADIUS`   Classical electron radius, :math:`r_e` in :unit:`[m]`
:const:`alphaem`    :c:macro:`P_ALPHAEM`   :const:`mad_cst_ALPHAEM`   Fine-structure constant, :math:`\alpha`
==================  =====================  =========================  ======================

.. index::
   physical constants
   CODATA

Mathematical Functions
======================

Generic Real-like Functions
---------------------------

Real-like generic functions forward the call to the method of the same name from the first argument when the latter is not a :type:`number`. The C functions column lists the C implementation used when the argument is a :type:`number` and the implementation does not rely on the standard :code:`math` module.

======================  =======================================================  =============
Functions               Return values                                            C functions
======================  =======================================================  =============
:func:`abs(x)`          :math:`|x|`
:func:`acos(x)`         :math:`\cos^{-1} x`
:func:`acosh(x)`        :math:`\cosh^{-1} x`                                     :c:func:`acosh`
:func:`acot(x)`         :math:`\cot^{-1} x`
:func:`acoth(x)`        :math:`\coth^{-1} x`                                     :c:func:`atanh`
:func:`asin(x)`         :math:`\sin^{-1} x`
:func:`asinc(x)`        :math:`\frac{\sin^{-1} x}{x}`                            :c:func:`mad_num_asinc`
:func:`asinh(x)`        :math:`\sinh^{-1} x`                                     :c:func:`asinh`
:func:`asinhc(x)`       :math:`\frac{\sinh^{-1} x}{x}`                           :c:func:`mad_num_asinhc`
:func:`atan(x)`         :math:`\tan^{-1} x`
:func:`atan2(x,y)`      :math:`\tan^{-1} \frac{x}{y}`
:func:`atanh(x)`        :math:`\tanh^{-1} x`                                     :c:func:`atanh`
:func:`ceil(x)`         :math:`\operatorname{ceil}(x)`
:func:`cos(x)`          :math:`\cos x`
:func:`cosh(x)`         :math:`\cosh x`
:func:`cot(x)`          :math:`\cot x`
:func:`coth(x)`         :math:`\coth x`
:func:`exp(x)`          :math:`\exp x`
:func:`floor(x)`        :math:`\operatorname{floor}(x)`
:func:`fact(n)`         :math:`n!`                                               :c:func:`mad_num_fact` [#]_
:func:`frac(x)`         :math:`\operatorname{frac}(x)`
:func:`hypot(x,y)`      :math:`\sqrt{x^2+y^2}`                                   :c:func:`hypot`
:func:`hypot3(x,y,z)`   :math:`\sqrt{x^2+y^2+z^2}`                               :c:func:`hypot`
:func:`invsqrt(x,v_)`   :math:`\frac{v}{\sqrt x}`
:func:`invfact(n)`      :math:`\frac{1}{n!}`                                     :c:func:`mad_num_invfact`
:func:`log(x)`          :math:`\log x`
:func:`log10(x)`        :math:`\log_{10} x`
:func:`pow(x,y)`        :math:`x^y`
:func:`powi(x,n)`       :math:`x^n`                                              :c:func:`mad_num_powi`
:func:`rangle(a,r)`     :math:`a + 2\pi \operatorname{round}(\frac{r-a}{2\pi})`  :c:func:`round`
:func:`round(x)`        :math:`\operatorname{round}(x)`                          :c:func:`round`
:func:`sign(x)`         :math:`-1, 0\text{ or }1`                                :c:func:`mad_num_sign`
:func:`sign1(x)`        :math:`-1\text{ or }1`                                   :c:func:`mad_num_sign1` [#]_
:func:`sin(x)`          :math:`\sin x`
:func:`sinc(x)`         :math:`\frac{\sin x}{x}`                                 :c:func:`mad_num_sinc`
:func:`sinh(x)`         :math:`\sinh x`
:func:`sinhc(x)`        :math:`\frac{\sinh x}{x}`                                :c:func:`mad_num_sinhc`
:func:`sqrt(x)`         :math:`\sqrt{x}`
:func:`tan(x)`          :math:`\tan x`
:func:`tanh(x)`         :math:`\tanh x`
:func:`lgamma(x,tol)`   :math:`\ln|\Gamma(x)|`                                   :c:func:`lgamma`
:func:`tgamma(x,tol)`   :math:`\Gamma(x)`                                        :c:func:`tgamma`
:func:`trunc(x)`        :math:`\operatorname{trunc}(x)`
:func:`unit(x)`         :math:`\frac{x}{|x|}`
======================  =======================================================  =============

Generic Complex-like Functions
------------------------------

Complex-like generic functions forward the call to the method of the same name from the first argument when the latter is not a :type:`number`, otherwise it implements a real-like compatibility layer using the equivalent representation :math:`z=x+0i`.

====================  ==================================
Functions             Return values
====================  ==================================
:func:`cabs(z)`       :math:`|z|`
:func:`carg(z)`       :math:`\arg z`
:func:`conj(z)`       :math:`z^*`
:func:`cplx(x,y)`     :math:`x+i\,y`
:func:`imag(z)`       :math:`\Im(z)`
:func:`polar(z)`      :math:`|z|\,e^{i \arg z}`
:func:`proj(z)`       :math:`\operatorname{proj}(z)`
:func:`real(z)`       :math:`\Re(z)`
:func:`rect(z)`       :math:`\Re(z)\cos \Im(z)+i\,\Re(z)\sin \Im(z)`
:func:`reim(z)`       :math:`\Re(z), \Im(z)`
====================  ==================================

Generic Error-like Functions
----------------------------

Error-like generic functions forward the call to the method of the same name from the first argument when the latter is not a :type:`number`, otherwise it calls C wrappers to the corresponding functions from the `Faddeeva library <http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package>`_ from the MIT (see :file:`mad_num.c`).

======================  ==========================================================  ======================
Functions               Return values                                               C functions  
======================  ==========================================================  ======================
:func:`erf(z,tol_)`     :math:`\frac{2}{\sqrt\pi}\int_0^z e^{-t^2} dt`              :c:func:`mad_num_erf`      
:func:`erfc(z,tol_)`    :math:`1-\operatorname{erf}(z)`                             :c:func:`mad_num_erfc`     
:func:`erfi(z,tol_)`    :math:`-i\operatorname{erf}(i z)`                           :c:func:`mad_num_erfi`     
:func:`erfcx(z,tol_)`   :math:`e^{z^2}\operatorname{erfc}(z)`                       :c:func:`mad_num_erfcx`    
:func:`wf(z,tol_)`      :math:`e^{-z^2}\operatorname{erfc}(-i z)`                   :c:func:`mad_num_wf`       
:func:`dawson(z,tol_)`  :math:`\frac{-i\sqrt\pi}{2}e^{-z^2}\operatorname{erf}(iz)`  :c:func:`mad_num_dawson`
======================  ==========================================================  ======================

Generic MapFold-like Functions
------------------------------

MapFold-like generic functions (also known as MapReduce) forward the call to the method of the same name from the first argument when the latter is not a :type:`number`. These functions are useful when used as high-order functions passed to methods :func:`:map2()`, :func:`:foldl()` (fold left) or :func:`:foldr()` (fold right) of containers like vectors and matrices.

====================  ========================
Functions             Return values
====================  ========================
:func:`sumsqr(x,y)`   :math:`x^2 + y^2`
:func:`sumabs(x,y)`   :math:`|x| + |y|`
:func:`minabs(x,y)`   :math:`\min(|x|, |y|)`
:func:`maxabs(x,y)`   :math:`\max(|x|, |y|)`
:func:`sumysqr(x,y)`  :math:`x + y^2`
:func:`sumyabs(x,y)`  :math:`x + |y|`
:func:`minyabs(x,y)`  :math:`\min(x, |y|)`
:func:`maxyabs(x,y)`  :math:`\max(x, |y|)`
:func:`sumxsqr(x,y)`  :math:`x^2 + y`
:func:`sumxabs(x,y)`  :math:`|x| + y`
:func:`minxabs(x,y)`  :math:`\min(|x|, y)`
:func:`maxxabs(x,y)`  :math:`\max(|x|, y)`
====================  ========================

Functions for Circular Sector
-----------------------------

Basic functions for arc and cord lengths conversion rely on the following elementary relations:

.. math::

    l_{\text{arc}}  &= a r = \frac{l_{\text{cord}}}{\operatorname{sinc} \frac{a}{2}}

    l_{\text{cord}} &= 2 r \sin \frac{a}{2} = l_{\text{arc}} \operatorname{sinc} \frac{a}{2} 

where :math:`r` stands for the radius and :math:`a` for the angle of the `Circular Sector <https://en.wikipedia.org/wiki/Circular_sector>`_.

=====================  ==========================
Functions              Return values
=====================  ==========================
:func:`arc2cord(l,a)`  :math:`l_{\text{arc}} \operatorname{sinc} \frac{a}{2}`
:func:`arc2len(l,a)`   :math:`l_{\text{arc}} \operatorname{sinc} \frac{a}{2}\, \cos a`
:func:`cord2arc(l,a)`  :math:`\frac{l_{\text{cord}}}{\operatorname{sinc} \frac{a}{2}}`
:func:`cord2len(l,a)`  :math:`l_{\text{cord}} \cos a`
:func:`len2arc(l,a)`   :math:`\frac{l}{\operatorname{sinc} \frac{a}{2}\, cos a}`
:func:`len2cord(l,a)`  :math:`\frac{l}{\cos a}`
=====================  ==========================

Pseudo-Random Number Generators
===============================

The module :mod:`gmath` provides an implementation of the *Xoshiro256\*\**  variant of the `XorShift <https://en.wikipedia.org/wiki/Xorshift>`_ PRNG familly [XORSHFT03]_, an all-purpose, rock-solid generator with a period of :math:`2^{256}-1` that supports long jumps of period :math:`2^{128}`. This PRNG is also the default implementation of recent versions of Lua (not LuaJIT, see below) and GFortran. See https://prng.di.unimi.it for details about xoshiro/xoroshiro PRNGs.

The module :mod:`math` of LuaJIT provides an implementation of the *Tausworthe* PRNG [TAUSWTH96]_, which has a period of :math:`2^{223}` but doesn't support long jumps, and hence uses a single global PRNG.

The module :mod:`gmath` also provides an implementation of the simple global PRNG of MAD-X for comparison.

It's worth mentionning that none of these PRNG are cryptographically secure generators, but MAD-X PRNG excepted, they are nevertheless superior to the commonly used *Mersenne Twister* PRNG [MERTWIS98]_.

All PRNG *functions* (except constructors) are wrappers around PRNG *methods* with the same name, and expect an optional PRNG :obj:`rng_` as first parameter. If this optional PRNG :obj:`rng_` is omitted, i.e. not provided, these functions will use the current global PRNG by default.

.. function:: randnew ()

   Return a new Xoshiro256\*\* PRNG with a period of :math:`2^{128}` that is garuanteed to not overlapp with any other Xoshiro256\*\* PRNGs, unless it is initialized with a seed.

.. function:: xrandnew ()

   Return a new MAD-X PRNG initialized with default seed 123456789. Hence, all new MAD-X PRNG will generate the same sequence until they are initialized with a user-defined seed.

.. function:: randset (rng_)

   Set the current global PRNG to :obj:`rng` (if provided) and return the previous global PRNG.

.. function:: randseed (rng_, seed)
              rng:randseed (seed)

   Set the seed of the PRNG :obj:`rng` to :var:`seed`.

.. function:: rand (rng_)
              rng:rand ()

   Return a new pseudo-random number in the range ``[0, 1)`` from the PRNG :obj:`rng`.

.. function:: randi (rng_)
              rng:randi ()
              
   Return a new pseudo-random number in the range of a :type:`u64_t` from the PRNG :obj:`rng` (:type:`u32_t` for the MAD-X PRNG), see C API below for details.

.. function:: randn (rng_)
              rng:randn ()

   Return a new pseudo-random gaussian number in the range ``[-inf, +inf]`` from the PRNG :obj:`rng` by using the Box-Muller transformation (Marsaglia's polar form) to a peuso-random number in the range ``[0, 1)``.

.. function:: randtn (rng_, cut_)
              rng:randtn (cut_)

   Return a new truncated pseudo-random gaussian number in the range ``[-cut_, +cut_]`` from the PRNG :obj:`rng` by using iteratively the method :func:`rng:randn()`. This simple algorithm is actually used for compatibility with MAD-X.
   Default: :code:`cut_ = +inf`.

.. function:: randp (rng_, lmb_)
              rng:randp (lmb_)

   Return a new pseudo-random poisson number in the range ``[0, +inf]`` from the PRNG :obj:`rng` with parameter :math:`\lambda > 0` by using the *inverse transform sampling* method on peuso-random numbers.
   Default: :code:`lmb_ = 1`.

C API
-----

.. c:type:: prng_state_t
            xrng_state_t

   The Xoshiro256\*\* and the MAD-X PRNG types.

.. c:function:: num_t mad_num_rand (prng_state_t*)

   Return a pseudo-random double precision float in the range ``[0, 1)``. 

.. c:function:: u64_t mad_num_randi (prng_state_t*)

   Return a pseudo-random 64 bit unsigned integer in the range ``[0, ULLONG_MAX]``.

.. c:function:: void mad_num_randseed (prng_state_t*, num_t seed)

   Set the seed of the PRNG.

.. c:function:: void mad_num_randjump (prng_state_t*)

   Apply a jump to the PRNG like if :math:`2^{128}` pseudo-random numbers were generated. Hence PRNGs with different number of jumps will never overlapp. This function is applied to new PRNGs with an incremental number of jumps. 

.. c:function:: num_t mad_num_xrand (xrng_state_t*)

   Return a pseudo-random double precision float in the range ``[0, 1)`` from the MAD-X PRNG.

.. c:function:: u32_t mad_num_xrandi (xrng_state_t*)

   Return a pseudo-random 32 bit unsigned integer in the range ``[0, UINT_MAX]`` from the MAD-X PRNG.

.. c:function:: void mad_num_xrandseed (xrng_state_t*, u32_t seed)

   Set the seed of the MAD-X PRNG.

.. ------------------------------------------------------------

.. rubric:: Footnotes

.. [#] Canonical NaN bit patterns may differ between MAD and C for the mantissa, but both should exibit the same behavior.
.. [#] Factorial and inverse factorial support negative integers as input as it uses extended factorial definition.
.. [#] Sign and sign1 functions take care of special cases like ±0, ±inf and ±NaN.

References
==========

.. [XORSHFT03] G. Marsaglia, *"Xorshift RNGs"*, Journal of Statistical Software, 8 (14), July 2003. doi:10.18637/jss.v008.i14.

.. [TAUSWTH96] P. L’Ecuyer, *“Maximally Equidistributed Combined Tausworthe Generators”*, Mathematics of Computation, 65 (213), 1996, p203–213.

.. [MERTWIS98] M. Matsumoto and T. Nishimura, *“Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator”*. ACM Trans. on Modeling and Comp. Simulation, 8 (1), Jan. 1998, p3–30.