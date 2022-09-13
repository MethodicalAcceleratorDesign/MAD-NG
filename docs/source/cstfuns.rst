.. index::
   Constants and functions

***********************
Constants and Functions
***********************

This chapter describes some constants and functions provided by the modules :mod:`MAD.constant`, :mod:`MAD.gmath` and :mod:`MAD.gfunc`.

The module :mod:`gmath` extends the standard LUA module :mod:`math` with *generic* functions working on any types that support the methods with the same names. For example, the code :func:`gmath.sin(a)` will call :func:`math.sin(a)` if :var:`a` is a :type:`number`, otherwise it will call the method :func:`a:sin()`, i.e. delegate the invocation to :obj:`a`. This is how MAD-NG handles several types like :type:`numbers`, :type:`complex` number and :type:`TPSA` within a single *polymorphic* code that expects scalar-like behavior.

The module :mod:`gfunc` provides useful functions to help dealing with operators as functions and to manipulate functions in a `functional <https://en.wikipedia.org/wiki/Functional_programming>`_ way [#f4]_. It also provide a complete set of functions to create, combine and use *functors*, i.e. objects that behave like functions with :type:`callable` semantic.

Mathematical Constants
======================

This section describes basic mathematical constants uniquely defined as macros in the C header :file:`mad_cst.h` and available from C and MAD modules as floating point double precision variables. If these mathematical constants are already provided by the system libraries, they will be used instead of their local definitions.

===================  ======================  =========================  ======================
MAD constants        C macros                C constants                Values
===================  ======================  =========================  ======================
:const:`eps`         :c:macro:`DBL_EPSILON`  :const:`mad_cst_EPS`       Smallest representable step near one
:const:`tiny`        :c:macro:`DBL_MIN`      :const:`mad_cst_TINY`      Smallest representable number
:const:`huge`        :c:macro:`DBL_MAX`      :const:`mad_cst_HUGE`      Largest representable number
:const:`inf`         :c:macro:`INFINITY`     :const:`mad_cst_INF`       Positive infinity, :math:`1/0`
:const:`nan`         :c:macro:`NAN`          :const:`mad_cst_NAN`       Canonical NaN [#f1]_, :math:`0/0`
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

===============  ===================  =======================  ======================
MAD constants    C macros             C constants              Values
===============  ===================  =======================  ======================
:const:`minlen`  :c:macro:`P_MINLEN`  :const:`mad_cst_MINLEN`  Min length tolerance, default :math:`10^{-10}` in :unit:`[m]`
:const:`minang`  :c:macro:`P_MINANG`  :const:`mad_cst_MINANG`  Min angle tolerance, default :math:`10^{-10}` in :unit:`[1/m]`
:const:`minstr`  :c:macro:`P_MINSTR`  :const:`mad_cst_MINSTR`  Min strength tolerance, default :math:`10^{-10}` in :unit:`[rad]`
===============  ===================  =======================  ======================

The following table lists some physical constants from the `CODATA 2018 <https://physics.nist.gov/cuu/pdf/wall_2018.pdf>`_ sheet.

=================  =====================  =========================  ======================
MAD constants      C macros               C constants                Values
=================  =====================  =========================  ======================
:const:`clight`    :c:macro:`P_CLIGHT`    :const:`mad_cst_CLIGHT`    Speed of light, :math:`c` in :unit:`[m/s]`
:const:`mu0`       :c:macro:`P_MU0`       :const:`mad_cst_MU0`       Permeability of vacuum, :math:`\mu_0` in :unit:`[T.m/A]`
:const:`epsilon0`  :c:macro:`P_EPSILON0`  :const:`mad_cst_EPSILON0`  Permittivity of vacuum, :math:`\epsilon_0` in :unit:`[F/m]`
:const:`qelect`    :c:macro:`P_QELECT`    :const:`mad_cst_QELECT`    Elementary electric charge, :math:`e` in :unit:`[C]`
:const:`hbar`      :c:macro:`P_HBAR`      :const:`mad_cst_HBAR`      Reduced Plack's constant, :math:`\hbar` in :unit:`[GeV.s]`
:const:`amass`     :c:macro:`P_AMASS`     :const:`mad_cst_AMASS`     Unified atomic mass, :math:`m_u\,c^2` in :unit:`[GeV]`
:const:`emass`     :c:macro:`P_EMASS`     :const:`mad_cst_EMASS`     Electron mass, :math:`m_e\,c^2` in :unit:`[GeV]`
:const:`pmass`     :c:macro:`P_PMASS`     :const:`mad_cst_PMASS`     Proton mass, :math:`m_p\,c^2` in :unit:`[GeV]`
:const:`nmass`     :c:macro:`P_NMASS`     :const:`mad_cst_NMASS`     Neutron mass, :math:`m_n\,c^2` in :unit:`[GeV]`
:const:`mumass`    :c:macro:`P_MUMASS`    :const:`mad_cst_MUMASS`    Muon mass, :math:`m_{\mu}\,c^2` in :unit:`[GeV]`
:const:`deumass`   :c:macro:`P_DEUMASS`   :const:`mad_cst_DEUMASS`   Deuteron mass, :math:`m_d\,c^2` in :unit:`[GeV]`
:const:`eradius`   :c:macro:`P_ERADIUS`   :const:`mad_cst_ERADIUS`   Classical electron radius, :math:`r_e` in :unit:`[m]`
:const:`alphaem`   :c:macro:`P_ALPHAEM`   :const:`mad_cst_ALPHAEM`   Fine-structure constant, :math:`\alpha`
=================  =====================  =========================  ======================

.. index::
   physical constants
   CODATA

Mathematical Functions
======================

Generic Real-like Functions
---------------------------

Real-like generic functions forward the call to the method of the same name from the first argument when the latter is not a :type:`number`. The optional argument :var:`r_` represents a destination for results with reference semantic, i.e. avoiding memory allocation, which is ignored by results with value semantic. The C functions column lists the C implementation used when the argument is a :type:`number` and the implementation does not rely on the standard :code:`math` module.

===============================  =======================================================  =============
Functions                        Return values                                            C functions
===============================  =======================================================  =============
:func:`abs(x,r_)`                :math:`|x|`
:func:`acos(x,r_)`               :math:`\cos^{-1} x`
:func:`acosh(x,r_)`              :math:`\cosh^{-1} x`                                     :c:func:`acosh`
:func:`acot(x,r_)`               :math:`\cot^{-1} x`
:func:`acoth(x,r_)`              :math:`\coth^{-1} x`                                     :c:func:`atanh`
:func:`asin(x,r_)`               :math:`\sin^{-1} x`
:func:`asinc(x,r_)`              :math:`\frac{\sin^{-1} x}{x}`                            :c:func:`mad_num_asinc`
:func:`asinh(x,r_)`              :math:`\sinh^{-1} x`                                     :c:func:`asinh`
:func:`asinhc(x,r_)`             :math:`\frac{\sinh^{-1} x}{x}`                           :c:func:`mad_num_asinhc`
:func:`atan(x,r_)`               :math:`\tan^{-1} x`
:func:`atan2(x,y,r_)`            :math:`\tan^{-1} \frac{x}{y}`
:func:`atanh(x,r_)`              :math:`\tanh^{-1} x`                                     :c:func:`atanh`
:func:`ceil(x,r_)`               :math:`\operatorname{ceil}(x)`
:func:`cos(x,r_)`                :math:`\cos x`
:func:`cosh(x,r_)`               :math:`\cosh x`
:func:`cot(x,r_)`                :math:`\cot x`
:func:`coth(x,r_)`               :math:`\coth x`
:func:`exp(x,r_)`                :math:`\exp x`
:func:`floor(x,r_)`              :math:`\operatorname{floor}(x)`
:func:`frac(x,r_)`               :math:`\operatorname{frac}(x)`
:func:`hypot(x,y,r_)`            :math:`\sqrt{x^2+y^2}`                                   :c:func:`hypot`
:func:`hypot3(x,y,z,r_)`         :math:`\sqrt{x^2+y^2+z^2}`                               :c:func:`hypot`
:func:`inv(x,v_,r_)` [#f2]_      :math:`\frac{v}{x}`
:func:`invsqrt(x,v_,r_)` [#f2]_  :math:`\frac{v}{\sqrt x}`
:func:`lgamma(x,tol_,r_)`        :math:`\ln|\Gamma(x)|`                                   :c:func:`lgamma`
:func:`log(x,r_)`                :math:`\log x`
:func:`log10(x,r_)`              :math:`\log_{10} x`
:func:`pow(x,y,r_)`              :math:`x^y`
:func:`powi(x,n,r_)`             :math:`x^n`                                              :c:func:`mad_num_powi`
:func:`rangle(a,r)`              :math:`a + 2\pi \operatorname{round}(\frac{r-a}{2\pi})`  :c:func:`round`
:func:`round(x,r_)`              :math:`\operatorname{round}(x)`                          :c:func:`round`
:func:`sign(x)`                  :math:`-1, 0\text{ or }1`                                :c:func:`mad_num_sign`  [#f3]_
:func:`sign1(x)`                 :math:`-1\text{ or }1`                                   :c:func:`mad_num_sign1` [#f3]_
:func:`sin(x,r_)`                :math:`\sin x`
:func:`sinc(x,r_)`               :math:`\frac{\sin x}{x}`                                 :c:func:`mad_num_sinc`
:func:`sinh(x,r_)`               :math:`\sinh x`
:func:`sinhc(x,r_)`              :math:`\frac{\sinh x}{x}`                                :c:func:`mad_num_sinhc`
:func:`sqrt(x,r_)`               :math:`\sqrt{x}`
:func:`tan(x,r_)`                :math:`\tan x`
:func:`tanh(x,r_)`               :math:`\tanh x`
:func:`tgamma(x,tol_,r_)`        :math:`\Gamma(x)`                                        :c:func:`tgamma`
:func:`trunc(x,r_)`              :math:`\operatorname{trunc}(x)`
:func:`unit(x,r_)`               :math:`\frac{x}{|x|}`
===============================  =======================================================  =============

Generic Complex-like Functions
------------------------------

Complex-like generic functions forward the call to the method of the same name from the first argument when the latter is not a :type:`number`, otherwise it implements a real-like compatibility layer using the equivalent representation :math:`z=x+0i`. The optional argument :var:`r_` represents a destination for results with reference semantic, i.e. avoiding memory allocation, which is ignored by results with value semantic. 

=======================  ==================================
Functions                Return values
=======================  ==================================
:func:`cabs(z,r_)`       :math:`|z|`
:func:`carg(z,r_)`       :math:`\arg z`
:func:`conj(z,r_)`       :math:`z^*`
:func:`cplx(x,y,r_)`     :math:`x+i\,y`
:func:`imag(z,r_)`       :math:`\Im(z)`
:func:`polar(z,r_)`      :math:`|z|\,e^{i \arg z}`
:func:`proj(z,r_)`       :math:`\operatorname{proj}(z)`
:func:`real(z,r_)`       :math:`\Re(z)`
:func:`rect(z,r_)`       :math:`\Re(z)\cos \Im(z)+i\,\Re(z)\sin \Im(z)`
:func:`reim(z,re_,im_)`  :math:`\Re(z), \Im(z)`
=======================  ==================================

Generic Vector-like Functions
-----------------------------

Vector-like functions (also known as MapFold or MapReduce) are functions useful when used as high-order functions passed to methods like :func:`:map2()`, :func:`:foldl()` (fold left) or :func:`:foldr()` (fold right) of containers like lists, vectors and matrices.

====================  ========================
Functions             Return values
====================  ========================
:func:`sumsqr(x,y)`   :math:`x^2 + y^2`
:func:`sumabs(x,y)`   :math:`|x| + |y|`
:func:`minabs(x,y)`   :math:`\min(|x|, |y|)`
:func:`maxabs(x,y)`   :math:`\max(|x|, |y|)`
:func:`sumsqrl(x,y)`  :math:`x + y^2`
:func:`sumabsl(x,y)`  :math:`x + |y|`
:func:`minabsl(x,y)`  :math:`\min(x, |y|)`
:func:`maxabsl(x,y)`  :math:`\max(x, |y|)`
:func:`sumsqrr(x,y)`  :math:`x^2 + y`
:func:`sumabsr(x,y)`  :math:`|x| + y`
:func:`minabsr(x,y)`  :math:`\min(|x|, y)`
:func:`maxabsr(x,y)`  :math:`\max(|x|, y)`
====================  ========================

Generic Error-like Functions
----------------------------

Error-like generic functions forward the call to the method of the same name from the first argument when the latter is not a :type:`number`, otherwise it calls C wrappers to the corresponding functions from the `Faddeeva library <http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package>`_ from the MIT (see :file:`mad_num.c`). The optional argument :var:`r_` represents a destination for results with reference semantic, i.e. avoiding memory allocation, which is ignored by results with value semantic.

==========================  ==========================================================  ========================
Functions                   Return values                                               C functions  
==========================  ==========================================================  ========================
:func:`erf(z,rtol_,r_)`     :math:`\frac{2}{\sqrt\pi}\int_0^z e^{-t^2} dt`              :c:func:`mad_num_erf`      
:func:`erfc(z,rtol_,r_)`    :math:`1-\operatorname{erf}(z)`                             :c:func:`mad_num_erfc`     
:func:`erfi(z,rtol_,r_)`    :math:`-i\operatorname{erf}(i z)`                           :c:func:`mad_num_erfi`     
:func:`erfcx(z,rtol_,r_)`   :math:`e^{z^2}\operatorname{erfc}(z)`                       :c:func:`mad_num_erfcx`    
:func:`wf(z,rtol_,r_)`      :math:`e^{-z^2}\operatorname{erfc}(-i z)`                   :c:func:`mad_num_wf`       
:func:`dawson(z,rtol_,r_)`  :math:`\frac{-i\sqrt\pi}{2}e^{-z^2}\operatorname{erf}(iz)`  :c:func:`mad_num_dawson`
==========================  ==========================================================  ========================

Special Functions
-----------------

The special functions factorial and inverse factorial support negative integers as input as it uses extended factorial definition. The value are cached making the complexity of these functions in :math:`O(1)` after warmup. 

==================  ====================  =========================
Functions           Return values         C functions
==================  ====================  =========================
:func:`fact(n)`     :math:`n!`            :c:func:`mad_num_fact`
:func:`invfact(n)`  :math:`\frac{1}{n!}`  :c:func:`mad_num_invfact`
==================  ====================  =========================

Functions for Circular Sector
-----------------------------

Basic functions for arc and cord lengths conversion rely on the following elementary relations:

.. math::

    l_{\text{arc}}  &= a r = \frac{l_{\text{cord}}}{\operatorname{sinc} \frac{a}{2}}

    l_{\text{cord}} &= 2 r \sin \frac{a}{2} = l_{\text{arc}} \operatorname{sinc} \frac{a}{2} 

where :math:`r` stands for the radius and :math:`a` for the angle of the `Circular Sector <https://en.wikipedia.org/wiki/Circular_sector>`_.

=====================  =====================================
Functions              Return values
=====================  =====================================
:func:`arc2cord(l,a)`  :math:`l_{\text{arc}} \operatorname{sinc} \frac{a}{2}`
:func:`arc2len(l,a)`   :math:`l_{\text{arc}} \operatorname{sinc} \frac{a}{2}\, \cos a`
:func:`cord2arc(l,a)`  :math:`\frac{l_{\text{cord}}}{\operatorname{sinc} \frac{a}{2}}`
:func:`cord2len(l,a)`  :math:`l_{\text{cord}} \cos a`
:func:`len2arc(l,a)`   :math:`\frac{l}{\operatorname{sinc} \frac{a}{2}\, cos a}`
:func:`len2cord(l,a)`  :math:`\frac{l}{\cos a}`
=====================  =====================================

.. ----------------------------------------------

Operators as Functions
======================

The module :mod:`MAD.gfunc` provides many functions that are named version of operators and useful when operators cannot be used directly, e.g. when passed as argument or to compose together. These functions can also be retrieved from the module :mod:`MAD.gfunc.opstr` by their associated string (if available).

Math Operators
--------------

Functions for math operators are wrappers to associated mathematical operators, which themselves can be overridden by their associated metamethods.

================  =================  ===============  ===================
Functions         Return values      Operator string  Metamethods
================  =================  ===============  ===================
:func:`unm(x)`    :math:`-x`         :const:`"~"`     :func:`__unm(x,_)`
:func:`inv(x)`    :math:`1 / x`      :const:`"1/"`    :func:`__div(1,x)`
:func:`sqr(x)`    :math:`x \cdot x`  :const:`"^2"`    :func:`__mul(x,x)`
:func:`add(x,y)`  :math:`x + y`      :const:`"+"`     :func:`__add(x,y)`
:func:`sub(x,y)`  :math:`x - y`      :const:`"-"`     :func:`__sub(x,y)`
:func:`mul(x,y)`  :math:`x \cdot y`  :const:`"*"`     :func:`__mul(x,y)`
:func:`div(x,y)`  :math:`x / y`      :const:`"/"`     :func:`__div(x,y)`
:func:`mod(x,y)`  :math:`x \mod y`   :const:`"%"`     :func:`__mod(x,y)`
:func:`pow(x,y)`  :math:`x ^ y`      :const:`"^"`     :func:`__pow(x,y)`
================  =================  ===============  ===================

Vector Operators
----------------

Functions for element-wise operators [#f5]_ are wrappers to associated mathematical operators of vector-like objects, which themselves can be overridden by their associated metamethods.

=================  =====================  ===============  ====================
Functions          Return values          Operator string  Metamethods
=================  =====================  ===============  ====================
:func:`emul(x,y)`  :math:`x\,.*\,y`       :const:`".*"`    :func:`__emul(x,y)`
:func:`ediv(x,y)`  :math:`x\,./\,y`       :const:`"./"`    :func:`__ediv(x,y)`
:func:`emod(x,y)`  :math:`x\,.\%\,y`      :const:`".%"`    :func:`__emod(x,y)`
:func:`epow(x,y)`  :math:`x\,.\hat\ \ y`  :const:`".^"`    :func:`__epow(x,y)`
=================  =====================  ===============  ====================

Logical Operators
-----------------

Functions for logical operators are wrappers to associated logical operators.

=================  ====================  ===============
Functions          Return values         Operator string
=================  ====================  ===============
:func:`lfalse()`   :const:`true`                                         
:func:`ltrue()`    :const:`false`                                        
:func:`lnot(x)`    :math:`\lnot x`       :const:`"!"`                      
:func:`lbool(x)`   :math:`\lnot\lnot x`  :const:`"!!"`                       
:func:`land(x,y)`  :math:`x \land y`     :const:`"&&"`                       
:func:`lor(x,y)`   :math:`x \lor y`      :const:`"||"`                       
=================  ====================  ===============

Relational Operators
--------------------

Functions for relational operators are wrappers to associated logical operators, which themselves can be overridden by their associated metamethods. Relational ordering operators are available only for objects that are ordered.

===============  ================  ==============================  =================
Functions        Return values     Operator string                 Metamethods
===============  ================  ==============================  =================
:func:`eq(x,y)`  :math:`x = y`     :const:`"=="`                   :func:`__eq(x,y)`
:func:`ne(x,y)`  :math:`x \neq y`  :const:`"!="` or :const:`"~="`  :func:`__eq(x,y)`
:func:`lt(x,y)`  :math:`x < y`     :const:`"<"`                    :func:`__lt(x,y)`
:func:`le(x,y)`  :math:`x <= y`    :const:`"<="`                   :func:`__le(x,y)`
:func:`gt(x,y)`  :math:`x > y`     :const:`">"`                    :func:`__le(x,y)`
:func:`ge(x,y)`  :math:`x >= y`    :const:`">="`                   :func:`__lt(x,y)`
===============  ================  ==============================  =================

Object Operators
----------------

Functions for object operators are wrappers to associated object operators, which themselves can be overridden by their associated metamethods.

===================  ==============  ===============  =================
Functions            Return values   Operator string  Metamethods
===================  ==============  ===============  =================
:func:`get(x,k)`     :math:`x[k]`    :const:`"->"`    :func:`__index(x,k)`
:func:`set(x,k,v)`   :math:`x[k]=v`  :const:`"<-"`    :func:`__newindex(x,k,v)`
:func:`len(x)`       :math:`\#x`     :const:`"#"`     :func:`__len(x)`
:func:`cat(x,y)`     :math:`x .. y`  :const:`".."`    :func:`__concat(x,y)`
:func:`call(x,...)`  :math:`x(...)`  :const:`"()"`    :func:`__call(x,...)`
===================  ==============  ===============  =================

Bitwise Functions
=================

Functions for bitwise operations are those from the LuaJIT module :mod:`bit` and imported into the module :mod:`MAD.gfunc` for convenience, see http://bitop.luajit.org/api.html for details. Note that all these functions have *value semantic* and normalise their arguments to the numeric range of a 32 bit integer before use.

====================  ====================================================
Functions             Return values         
====================  ====================================================
:func:`tobit(x)`      Return the normalized value of :var:`x` to the range of a 32 bit integer      
:func:`tohex(x,n_)`   Return the hex string of :var:`x` with :var:`n` digits (:math:`n<0` use caps)    
:func:`bnot(x)`       Return the bitwise reverse of :var:`x` bits    
:func:`band(x,...)`   Return the bitwise *AND* of all arguments     
:func:`bor(x,...)`    Return the bitwise *OR* of all arguments 
:func:`bxor(x,...)`   Return the bitwise *XOR* of all arguments
:func:`lshift(x,n)`   Return the bitwise left shift of :var:`x` by :var:`n` bits with 0-bit shift-in     
:func:`rshift(x,n)`   Return the bitwise right shift of :var:`x` by :var:`n` bits with 0-bit shift-in
:func:`arshift(x,n)`  Return the bitwise right shift of :var:`x` by :var:`n` bits with sign bit shift-in
:func:`rol(x,n)`      Return the bitwise left rotation of :var:`x` by :var:`n` bits      
:func:`ror(x,n)`      Return the bitwise right rotation of :var:`x` by :var:`n` bits     
:func:`bswap(x)`      Return the swapped bytes of :var:`x`, i.e. convert big endian to/from little endian       
====================  ====================================================

Flags Functions
---------------

A flag is 32 bit unsigned integer used to store up to 32 binary states with the convention that :const:`0` means disabled/cleared and :const:`1` means enabled/set. Functions on flags are useful aliases to, or combination of, bitwise operations to manipulate their states (i.e. their bits). Flags are mainly used by the object model to keep track of hidden and user-defined states in a compact and efficient format. 

===================  ====================================================
Functions            Return values         
===================  ====================================================
:func:`bset(x,n)`    Return the flag :var:`x` with state :var:`n` enabled
:func:`bclr(x,n)`    Return the flag :var:`x` with state :var:`n` disabled   
:func:`btst(x,n)`    Return :const:`true` if state :var:`n` is enabled in :var:`x`, :const:`false` otherwise      
:func:`fbit(n)`      Return a flag with only state :var:`n` enabled    
:func:`fnot(x)`      Return the flag :var:`x` with all states flipped
:func:`fset(x,...)`  Return the flag :var:`x` with disabled states flipped if enabled in any flag passed as argument
:func:`fcut(x,...)`  Return the flag :var:`x` with enabled states flipped if disabled in any flag passed as argument 
:func:`fclr(x,f)`    Return the flag :var:`x` with enabled states flipped if enabled in :var:`f`
:func:`ftst(x,f)`    Return :const:`true` if all states enabled in :var:`f` are enabled in :var:`x`, :const:`false` otherwise 
:func:`fall(x)`      Return :const:`true` if all states are enabled in :var:`x`, :const:`false` otherwise       
:func:`fany(x)`      Return :const:`true` if any state is enabled in :var:`x`, :const:`false` otherwise    
===================  ====================================================

Special Functions
=================

The module :mod:`MAD.gfunc` provides some useful functions when passed as argument or composed with other functions.

======================  ====================================================
Functions               Return values         
======================  ====================================================
:func:`narg(...)`       Return the number of arguments      
:func:`ident(...)`      Return all arguments unchanged, i.e. functional identity    
:func:`fnil()`          Return :const:`nil`, i.e. functional nil    
:func:`ftrue()`         Return :const:`true`, i.e. functional true
:func:`ffalse()`        Return :const:`false`, i.e. functional false
:func:`fzero()`         Return :const:`0`, i.e. functional zero
:func:`fone()`          Return :const:`1`, i.e. functional one     
:func:`first(a)`        Return first argument and discard the others
:func:`second(a,b)`     Return second argument and discard the others
:func:`third(a,b,c)`    Return third argument and discard the others      
:func:`swap(a,b)`       Return first and second arguments swapped and discard the other arguments   
:func:`swapv(a,b,...)`  Return first and second arguments swapped followed by the other arguments        
:func:`echo(...)`       Return all arguments unchanged after echoing them on stdout       
======================  ====================================================

.. ---------------------------------------

Functors
========

Functors are objects that behave like functions with :type:`callable` semantic, and also like readonly arrays with :type:`indexable` semantic, where the index is translated as a unique argument into the function call. They are mainly used by the object model to distinguish them from functions which are interpreted as deferred expressions and evaluated automatically on reading, and by the Survey and Track codes to handle (user-defined) actions. 

The module :mod:`MAD.gfunc` offers few functions to expert users for creating and manipulating them.

.. function:: functor(f)

   Return a :type:`functor` that encapsulates the function (or any callable object) :var:`f`. Calling the returned functor is like calling :var:`f` itself with the same arguments. 

.. function:: compose(f, g)

   Return a :type:`functor` that encapsulates the composition of :var:`f` and :var:`g`. Calling the returned functor is like calling :math:`(f \circ g)(\dots)`. The operator :code:`f ^ g` is a shortcut for :func:`compose` if :var:`f` is a :type:`functor`.

.. function:: chain(f, g)

   Return a :type:`functor` that encapsulates the calls chain of :var:`f` and :var:`g`. Calling the returned functor is like calling :math:`f(\dots) ; g(\dots)`. The operator :code:`f .. g` is a shortcut for :func:`chain` if :var:`f` is a :type:`functor`.

.. function:: achain(f, g)

   Return a :type:`functor` that encapsulates the *AND*-ed calls chain of :var:`f` and :var:`g`. Calling the returned functor is like calling :math:`f(\dots) \land g(\dots)`.

.. function:: ochain(f, g)

   Return a :type:`functor` that encapsulates the *OR*-ed calls chain of :var:`f` and :var:`g`. Calling the returned functor is like calling :math:`f(\dots) \lor g(\dots)`.

.. function:: bind1st(f, a)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`a` as its first argument. Calling the returned functor is like calling :math:`f(a,\dots)`.

.. function:: bind2nd(f, b)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`b` as its second argument. Calling the returned functor is like calling :math:`f(a,b,\dots)` where :var:`a` may or may not be provided.

.. function:: bind3rd(f, c)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`c` as its third argument. Calling the returned functor is like calling :math:`f(a,b,c,\dots)` where :var:`a` and :var:`b` may or may not be provided.

.. function:: bind2st(f, a, b)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`a` and :var:`b` as its two first arguments. Calling the returned functor is like calling :math:`f(a,b,\dots)`.

.. function:: bind3st(f, a, b, c)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`a`, :var:`b` and :var:`c` as its three first arguments. Calling the returned functor is like calling :math:`f(a,b,c,\dots)`.

.. function:: bottom()

   Return a :type:`functor` that encapsulates the identity function :func:`ident` to define the *bottom* symbol of functors. Bottom is also available in the operator strings table :mod:`opstr` as :const:`"_|_"`.

.. function:: is_functor(a)

   Return :const:`true` if :var:`a` is a :type:`functor`, :const:`false` otherwise. This function is also available from the module :mod:`MAD.typeid`.

.. ---------------------------------------

Pseudo-Random Number Generators
===============================

The module :mod:`gmath` provides an implementation of the *Xoshiro256\*\** (XOR/shift/rotate) variant of the `XorShift <https://en.wikipedia.org/wiki/Xorshift>`_ PRNG familly [XORSHFT03]_, an all-purpose, rock-solid generator with a period of :math:`2^{256}-1` that supports long jumps of period :math:`2^{128}`. This PRNG is also the default implementation of recent versions of Lua (not LuaJIT, see below) and GFortran. See https://prng.di.unimi.it for details about Xoshiro/Xoroshiro PRNGs.

The module :mod:`math` of LuaJIT provides an implementation of the *Tausworthe* PRNG [TAUSWTH96]_, which has a period of :math:`2^{223}` but doesn't support long jumps, and hence uses a single global PRNG.

The module :mod:`gmath` also provides an implementation of the simple global PRNG of MAD-X for comparison.

It's worth mentionning that none of these PRNG are cryptographically secure generators, they are nevertheless superior to the commonly used *Mersenne Twister* PRNG [MERTWIS98]_, with the exception of the MAD-X PRNG.

All PRNG *functions* (except constructors) are wrappers around PRNG *methods* with the same name, and expect an optional PRNG :obj:`prng_` as first parameter. If this optional PRNG :obj:`prng_` is omitted, i.e. not provided, these functions will use the current global PRNG by default.

Functions and Methods
---------------------

.. function:: randnew ()

   Return a new Xoshiro256\*\* PRNG with a period of :math:`2^{128}` that is garuanteed to not overlapp with any other Xoshiro256\*\* PRNGs, unless it is initialized with a seed.

.. function:: xrandnew ()

   Return a new MAD-X PRNG initialized with default seed 123456789. Hence, all new MAD-X PRNG will generate the same sequence until they are initialized with a user-defined seed.

.. function:: randset (prng_)

   Set the current global PRNG to :obj:`prng` (if provided) and return the previous global PRNG.

.. function:: randseed (prng_, seed)
              prng:randseed (seed)

   Set the seed of the PRNG :obj:`prng` to :var:`seed`.

.. function:: rand (prng_)
              prng:rand ()

   Return a new pseudo-random number in the range ``[0, 1)`` from the PRNG :obj:`prng`.

.. function:: randi (prng_)
              prng:randi ()
              
   Return a new pseudo-random number in the range of a :type:`u64_t` from the PRNG :obj:`prng` (:type:`u32_t` for the MAD-X PRNG), see C API below for details.

.. function:: randn (prng_)
              prng:randn ()

   Return a new pseudo-random gaussian number in the range ``[-inf, +inf]`` from the PRNG :obj:`prng` by using the Box-Muller transformation (Marsaglia's polar form) to a peuso-random number in the range ``[0, 1)``.

.. function:: randtn (prng_, cut_)
              prng:randtn (cut_)

   Return a new truncated pseudo-random gaussian number in the range ``[-cut_, +cut_]`` from the PRNG :obj:`prng` by using iteratively the method :func:`prng:randn()`. This simple algorithm is actually used for compatibility with MAD-X.
   Default: :code:`cut_ = +inf`.

.. function:: randp (prng_, lmb_)
              prng:randp (lmb_)

   Return a new pseudo-random poisson number in the range ``[0, +inf]`` from the PRNG :obj:`prng` with parameter :math:`\lambda > 0` by using the *inverse transform sampling* method on peuso-random numbers.
   Default: :code:`lmb_ = 1`.

.. function:: is_randgen(a)

   Return :const:`true` if :var:`a` is a PRNG, :const:`false` otherwise. This function is also available from the module :mod:`MAD.typeid`.

.. function:: is_xrandgen(a)

   Return :const:`true` if :var:`a` is a MAD-X PRNG, :const:`false` otherwise. This function is also available from the module :mod:`MAD.typeid`.

.. function:: isa_randgen(a)

   Return :const:`true` if :var:`a` is either a PRNG or a MAD-X PRNG, :const:`false` otherwise. This function is also available from the module :mod:`MAD.typeid`.

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

   Apply a jump to the PRNG as if :math:`2^{128}` pseudo-random numbers were generated. Hence PRNGs with different number of jumps will never overlap. This function is applied to new PRNGs with an incremental number of jumps. 

.. c:function:: num_t mad_num_xrand (xrng_state_t*)

   Return a pseudo-random double precision float in the range ``[0, 1)`` from the MAD-X PRNG.

.. c:function:: u32_t mad_num_xrandi (xrng_state_t*)

   Return a pseudo-random 32 bit unsigned integer in the range ``[0, UINT_MAX]`` from the MAD-X PRNG.

.. c:function:: void mad_num_xrandseed (xrng_state_t*, u32_t seed)

   Set the seed of the MAD-X PRNG.

References
==========

.. [XORSHFT03] G. Marsaglia, *"Xorshift RNGs"*, Journal of Statistical Software, 8 (14), July 2003. doi:10.18637/jss.v008.i14.

.. [TAUSWTH96] P. L’Ecuyer, *“Maximally Equidistributed Combined Tausworthe Generators”*, Mathematics of Computation, 65 (213), 1996, p203–213.

.. [MERTWIS98] M. Matsumoto and T. Nishimura, *“Mersenne Twister: A 623-dimensionally equidistributed uniform pseudorandom number generator”*. ACM Trans. on Modeling and Comp. Simulation, 8 (1), Jan. 1998, p3–30.

.. ------------------------------------------------------------

.. rubric:: Footnotes

.. [#f4] For *true* Functional Programming, see the module :mod:`MAD.lfun`, a binding of the `LuaFun <https://github.com/luafun/luafun>`_  library adapted to the ecosystem of MAD-NG.
.. [#f1] Canonical NaN bit patterns may differ between MAD and C for the mantissa, but both should exibit the same behavior.
.. [#f2] Default: :code:`v_ = 1`. 
.. [#f3] Sign and sign1 functions take care of special cases like ±0, ±inf and ±NaN.
.. [#f5] Element-wise operators are not available directly in the programming language, here we use the Matlab-like notation for convenience.
