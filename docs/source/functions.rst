.. index::
   Functions

*********
Functions
*********

This chapter describes some functions provided by the modules :mod:`MAD.gmath` and :mod:`MAD.gfunc`.

The module :mod:`gmath` extends the standard LUA module :mod:`math` with *generic* functions working on any types that support the methods with the same names. For example, the code :func:`gmath.sin(a)` will call :func:`math.sin(a)` if :var:`a` is a :type:`number`, otherwise it will call the method :func:`a:sin()`, i.e. delegate the invocation to :obj:`a`. This is how MAD-NG handles several types like :type:`numbers`, :type:`complex` number and :type:`TPSA` within a single *polymorphic* code that expects scalar-like behavior.

The module :mod:`gfunc` provides useful functions to help dealing with operators as functions and to manipulate functions in a `functional <https://en.wikipedia.org/wiki/Functional_programming>`_ way [#f1]_.

Mathematical Functions
======================

Generic Real-like Functions
---------------------------

Real-like generic functions forward the call to the method of the same name from the first argument when the latter is not a :type:`number`. The optional argument :var:`r_` represents a destination placeholder for results with reference semantic, i.e. avoiding memory allocation, which is ignored by results with value semantic. The C functions column lists the C implementation used when the argument is a :type:`number` and the implementation does not rely on the standard :code:`math` module but on functions provided with MAD-NG or by the standard math library described in the C Programming Language Standard [ISOC99]_.

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
:func:`ceil(x,r_)`               :math:`\lceil x \rceil`                                  
:func:`cos(x,r_)`                :math:`\cos x`
:func:`cosh(x,r_)`               :math:`\cosh x`
:func:`cot(x,r_)`                :math:`\cot x`
:func:`coth(x,r_)`               :math:`\coth x`
:func:`exp(x,r_)`                :math:`\exp x`
:func:`floor(x,r_)`              :math:`\lfloor x \rfloor`
:func:`frac(x,r_)`               :math:`x - \operatorname{trunc}(x)`
:func:`hypot(x,y,r_)`            :math:`\sqrt{x^2+y^2}`                                   :c:func:`hypot`
:func:`hypot3(x,y,z,r_)`         :math:`\sqrt{x^2+y^2+z^2}`                               :c:func:`hypot`
:func:`inv(x,v_,r_)` [#f2]_      :math:`\frac{v}{x}`
:func:`invsqrt(x,v_,r_)` [#f2]_  :math:`\frac{v}{\sqrt x}`
:func:`lgamma(x,tol_,r_)`        :math:`\ln|\Gamma(x)|`                                   :c:func:`lgamma`
:func:`log(x,r_)`                :math:`\log x`
:func:`log10(x,r_)`              :math:`\log_{10} x`
:func:`powi(x,n,r_)`             :math:`x^n`                                              :c:func:`mad_num_powi`
:func:`round(x,r_)`              :math:`\lfloor x \rceil`                                 :c:func:`round`
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
:func:`trunc(x,r_)`              :math:`\lfloor x \rfloor, x\geq 0;\lceil x \rceil, x<0`
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
:func:`fabs(z,r_)`       :math:`|\Re(z)|+i\,|\Im(z)|`
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

The special function :func:`fact()` supports negative integers as input as it uses extended factorial definition, and the values are cached to make its complexity in :math:`O(1)` after warmup.

The special function :func:`rangle()` adjust the angle :var:`a` versus the *previous* right angle :var:`r`, e.g. during phase advance accumulation, to ensure proper value when passing through the :math:`\pm 2k\pi` boundaries.

===================  ================================================  =========================
Functions            Return values                                     C functions
===================  ================================================  =========================
:func:`fact(n)`      :math:`n!`                                        :c:func:`mad_num_fact`
:func:`rangle(a,r)`  :math:`a + 2\pi \lfloor \frac{r-a}{2\pi} \rceil`  :c:func:`round`
===================  ================================================  =========================

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

Element Operators
-----------------

Functions for element-wise operators [#f4]_ are wrappers to associated mathematical operators of vector-like objects, which themselves can be overridden by their associated metamethods.

====================  =====================  ===============  ====================
Functions             Return values          Operator string  Metamethods
====================  =====================  ===============  ====================
:func:`emul(x,y,r_)`  :math:`x\,.*\,y`       :const:`".*"`    :func:`__emul(x,y,r_)`
:func:`ediv(x,y,r_)`  :math:`x\,./\,y`       :const:`"./"`    :func:`__ediv(x,y,r_)`
:func:`emod(x,y,r_)`  :math:`x\,.\%\,y`      :const:`".%"`    :func:`__emod(x,y,r_)`
:func:`epow(x,y,r_)`  :math:`x\,.\hat\ \ y`  :const:`".^"`    :func:`__epow(x,y,r_)`
====================  =====================  ===============  ====================

Logical Operators
-----------------

Functions for logical operators are wrappers to associated logical operators.

=================  ===============================================================  ===============
Functions          Return values                                                    Operator string
=================  ===============================================================  ===============
:func:`lfalse()`   :const:`true`                                                    :const:`"T"`                                
:func:`ltrue()`    :const:`false`                                                   :const:`"F"`                          
:func:`lnot(x)`    :math:`\lnot x`                                                  :const:`"!"`                      
:func:`lbool(x)`   :math:`\lnot\lnot x`                                             :const:`"!!"`                       
:func:`land(x,y)`  :math:`x \land y`                                                :const:`"&&"`                       
:func:`lor(x,y)`   :math:`x \lor y`                                                 :const:`"||"`                       
:func:`lnum(x)`    :math:`\lnot x\rightarrow 0`, :math:`\lnot\lnot x\rightarrow 1`  :const:`"!#"`
=================  ===============================================================  ===============

Relational Operators
--------------------

Functions for relational operators are wrappers to associated logical operators, which themselves can be overridden by their associated metamethods. Relational ordering operators are available only for objects that are ordered.

================  =========================  ==============================  =================
Functions         Return values              Operator string                 Metamethods
================  =========================  ==============================  =================
:func:`eq(x,y)`   :math:`x = y`              :const:`"=="`                   :func:`__eq(x,y)`
:func:`ne(x,y)`   :math:`x \neq y`           :const:`"!="` or :const:`"~="`  :func:`__eq(x,y)`
:func:`lt(x,y)`   :math:`x < y`              :const:`"<"`                    :func:`__lt(x,y)`
:func:`le(x,y)`   :math:`x \leq y`           :const:`"<="`                   :func:`__le(x,y)`
:func:`gt(x,y)`   :math:`x > y`              :const:`">"`                    :func:`__le(y,x)`
:func:`ge(x,y)`   :math:`x \geq y`           :const:`">="`                   :func:`__lt(y,x)`
:func:`cmp(x,y)`  :math:`(x > y) - (x < y)`  :const:`"?="`
================  =========================  ==============================  =================

The special relational operator :func:`cmp()` returns the number :const:`1` for :math:`x<y`, :const:`-1` for :math:`x>y`, and :const:`0` otherwise.

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

C API
=====

These functions are provided for performance reason and compliance with the C API of other modules.

.. c:function:: int mad_num_sign (num_t x)

   Return an integer amongst :const:`{-1, 0, 1}` representing the sign of the :type:`number` :var:`x`.

.. c:function:: int mad_num_sign1 (num_t x)

   Return an integer amongst :const:`{-1, 1}` representing the sign of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_fact (int n)

   Return the extended factorial the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_powi (num_t x, int n)

   Return the :type:`number` :var:`x` raised to the power of the :type:`integer` :var:`n` using a fast algorithm.

.. c:function:: num_t mad_num_sinc (num_t x)

   Return the sine cardinal of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_sinhc (num_t x)

   Return the hyperbolic sine cardinal of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_asinc (num_t x)

   Return the arc sine cardinal of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_asinhc (num_t x)

   Return the hyperbolic arc sine cardinal of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_wf (num_t x, num_t relerr)

   Return the Faddeeva function of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_erf (num_t x, num_t relerr)

   Return the error function of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_erfc (num_t x, num_t relerr) 

   Return the complementary error function of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_erfcx (num_t x, num_t relerr)

   Return the scaled complementary error function of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_erfi (num_t x, num_t relerr)

   Return the imaginary error function of the :type:`number` :var:`x`.

.. c:function:: num_t mad_num_dawson (num_t x, num_t relerr)

   Return the Dawson integral for the :type:`number` :var:`x`.

.. ------------------------------------------------------------

References
==========

.. [ISOC99] ISO/IEC 9899:1999 Programming Languages - C. https://www.iso.org/standard/29237.html.

.. rubric:: Footnotes

.. [#f1] For *true* Functional Programming, see the module :mod:`MAD.lfun`, a binding of the `LuaFun <https://github.com/luafun/luafun>`_  library adapted to the ecosystem of MAD-NG.
.. [#f2] Default: :expr:`v_ = 1`. 
.. [#f3] Sign and sign1 functions take care of special cases like ±0, ±inf and ±NaN.
.. [#f4] Element-wise operators are not available directly in the programming language, here we use the Matlab-like notation for convenience.
