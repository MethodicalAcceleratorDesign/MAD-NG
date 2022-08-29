.. index::
   single: elementary mathematical functions

********************
Elementary functions
********************

This chapter describes elementary functions provided by the module :code:`gmath`. This module extends the standard module :code:`math` with *generic* functions working on any type that implements the methods with the same name. For example, the code :code:`gmath.sin(a)` will call :code:`math.sin(a)` if :code:`a` is a :type:`number`, otherwise it will calling the method :code:`a:sin()`, i.e. delegate the call to :code:`a`. This is how MAD-NG handles few types like :type:`numbers`, :type:`complex` number, :type:`matrix` and :type:`TPSA` within a single code.

Generic Operators
=================

Generic operators are named functions that rely on associated operators, which themselves can be redefined by their associated metamethods.

====================  ==============  =============
Operators             Return values   Metamethods
====================  ==============  =============
:code:`unm(x)`        :code:`-x`      __unm(x,_)
:code:`add(x,y)`      :code:`x + y`   __add(x,y)
:code:`sub(x,y)`      :code:`x - y`   __sub(x,y)
:code:`mul(x,y)`      :code:`x * y`   __mul(x,y)
:code:`div(x,y)`      :code:`x / y`   __div(x,y)
:code:`mod(x,y)`      :code:`x % y`   __mod(x,y)
:code:`pow(x,y)`      :code:`x ^ y`   __pow(x,y)
:code:`sqr(x)`        :code:`x * x`   -
:code:`inv(x)`        :code:`1 / x`   -
:code:`emul(x,y,r_)`  :code:`x .* y`  __emul(x,y,r_)
:code:`ediv(x,y,r_)`  :code:`x ./ y`  __ediv(x,y,r_)
:code:`emod(x,y,r_)`  :code:`x .% y`  __emod(x,y,r_)
:code:`epow(x,y,r_)`  :code:`x .^ y`  __epow(x,y,r_)
====================  ==============  =============

Generic Functions (real-like)
=============================

Real-like generic functions forward the call to the method of the same name from the first argument when the later is not a :type:`number`.

======================  ==================================  =============
Functions               Return values                       C functions
======================  ==================================  =============
:code:`abs    (x)`      :math:`|x|`
:code:`acos   (x)`      :math:`\cos^{-1}(x)`
:code:`acosh  (x)`      :math:`\cosh^{-1}(x)`               :func:`acosh`
:code:`acot   (x)`      :math:`\cot^{-1}(x)`
:code:`acoth  (x)`      :math:`\coth^{-1}(x)`               :func:`atanh`
:code:`asin   (x)`      :math:`\sin^{-1}(x)`
:code:`asinc  (x)`      :math:`\frac{\sin^{-1}(x)}{x}`
:code:`asinh  (x)`      :math:`\sinh^{-1}(x)`               :func:`asinh`
:code:`asinhc (x)`      :math:`\frac{\sinh^{-1}(x)}{x}`
:code:`atan   (x)`      :math:`\tan^{-1}(x)`
:code:`atan2  (x,y)`    :math:`\tan^{-1}(\frac{x}{y})`
:code:`atanh  (x)`      :math:`\tanh^{-1}(x)`               :func:`atanh`
:code:`ceil   (x)`      :math:`\operatorname{ceil}(x)`
:code:`cos    (x)`      :math:`\cos(x)`
:code:`cosh   (x)`      :math:`\cosh(x)`
:code:`cot    (x)`      :math:`\cot(x)`
:code:`coth   (x)`      :math:`\coth(x)`
:code:`deg2rad(x)`      :math:`\frac{\pi}{180} x`
:code:`exp    (x)`      :math:`\exp(x)`
:code:`floor  (x)`      :math:`\operatorname{floor}(x)`
:code:`frac   (x)`      :math:`\operatorname{frac}(x)`
:code:`hypot  (x,y)`    :math:`\sqrt{x^2+y^2}`              :func:`hypot`
:code:`hypot3 (x,y,z)`  :math:`\sqrt{x^2+y^2+z^2}`          :func:`hypot`
:code:`invsqrt(x,v_)`   :math:`\frac{v}{\sqrt x}`
:code:`log    (x)`      :math:`\log(x)`
:code:`log10  (x)`      :math:`\operatorname{log10}(x)`
:code:`pow    (x,y)`    :math:`x^y`
:code:`rad2deg(x)`      :math:`\frac{180}{pi} x`
:code:`round  (x)`      :math:`\operatorname{round}(x)`     :func:`round`
:code:`sign   (x)`      :math:`-1, 0\text{ or }1`           :func:`mad_num_sign`
:code:`sign1  (x)`      :math:`-1\text{ or }1`              :func:`mad_num_sign1`
:code:`sin    (x)`      :math:`\sin(x)`
:code:`sinc   (x)`      :math:`\frac{\sin(x)}{x}`
:code:`sinh   (x)`      :math:`\sinh(x)`
:code:`sinhc  (x)`      :math:`\frac{\sinh(x)}{x}`
:code:`sqrt   (x)`      :math:`\sqrt{x}`
:code:`tan    (x)`      :math:`\tan(x)`
:code:`tanh   (x)`      :math:`\tanh(x)`
:code:`lgamma (x,tol)`  :math:`\ln|\Gamma(x)|`              :func:`lgamma`
:code:`tgamma (x,tol)`  :math:`\Gamma(x)`                   :func:`tgamma`
:code:`trunc  (x)`      :math:`\operatorname{trunc}(x)`
:code:`unit   (x)`      :math:`\frac{x}{|x|}`
======================  ==================================  =============

Generic Functions (complex-like)
================================

Complex-like generic functions forward the call to the method of the same name from the first argument when the later is not a :type:`number`, otherwise it implements a real-like compatibility layer using the equivalent representation :math:`x+0i`.

====================  ==================================
Functions             Return values
====================  ==================================
:code:`cabs (z)`      :math:`|z|`
:code:`carg (z)`      :math:`\arg(z)`
:code:`conj (z)`      :math:`z^*`
:code:`cplx (x,y)`    :math:`x+i\,y`
:code:`imag (z)`      :math:`\Im(z)`
:code:`polar(z)`      :math:`|z|\,e^{i\arg(z)}`
:code:`proj (z)`      :math:`\operatorname{Proj}(z)`
:code:`real (z)`      :math:`\Re(z)`
:code:`rect (z)`      :math:`\Re(z)\cos(\Im(z))+i\,\Re(z)\sin(\Im(z))`
:code:`reim (z)`      :math:`(\Re(z), \Im(z))`
====================  ==================================

Generic Functions (Error-like)
==============================

Error-like generic functions forward the call to the method of the same name from the first argument when the later is not a :type:`number`, otherwise it calls a C wrapper to corresponding function from the Faddeeva library from the MIT (see :file:`mad_num.c`).

====================  ======================  =======================
Functions             C functions for reals   C functions for complex 
====================  ======================  =======================
:code:`erf  (x,tol)`  :code:`mad_num_erf`     :code:`mad_cnum_erf`   
:code:`erfc (x,tol)`  :code:`mad_num_erfc`    :code:`mad_cnum_erfc`  
:code:`erfi (x,tol)`  :code:`mad_num_erfi`    :code:`mad_cnum_erfi`  
:code:`erfcx(x,tol)`  :code:`mad_num_erfcx`   :code:`mad_cnum_erfcx` 
:code:`wf   (x,tol)`  :code:`mad_num_wf`      :code:`mad_cnum_wf`    
====================  ======================  =======================

Generic Functions (Length-Angle based)
======================================

Length-Angle base generic function relies on the following elementary relations between length and angle.

.. math::
    l_{\text{arc}}  = a r = \frac{l_{\text{cord}}}{\operatorname{sinc}(\frac{a}{2})}
    l_{\text{cord}} = 2 r \sin(\frac{a}{2}) = l_{\text{arc}} \operatorname{sinc}(\frac{a}{2}) 

=====================  ==================================
Functions              Return values
=====================  ==================================
:code:`arc2cord(l,a)`  :math:`l \operatorname{sinc}(\frac{a}{2})`
:code:`arc2len (l,a)`  :math:`l \operatorname{sinc}(\frac{a}{2}) cos(a)`
:code:`cord2arc(l,a)`  :math:`\frac{l}{\operatorname{sinc}(\frac{a}{2})}`
:code:`cord2len(l,a)`  :math:`l cos(a)`
:code:`len2arc (l,a)`  :math:`\frac{l}{\operatorname{sinc}(\frac{a}{2}) cos(a)}`
:code:`len2cord(l,a)`  :math:`\frac{l}{cos(a)}`
:code:`rangle  (a,r)`  :math:`a + 2\pi \operatorname{round}(\frac{r-a}{2\pi})`
=====================  ==================================

Generic Functions (Folding-Left based)
======================================

=====================  ==================================
Functions              Return values
=====================  ==================================
:code:`sumsqr (x,y)`   :math:`x^2 + y^2`
:code:`sumabs (x,y)`   :math:`|x| + |y|`
:code:`minabs (x,y)`   :math:`\min(|x|, |y|)`
:code:`maxabs (x,y)`   :math:`\max(|x|, |y|)`
:code:`sumysqr(x,y)`   :math:`x + y^2`
:code:`sumyabs(x,y)`   :math:`x + |y|`
:code:`minyabs(x,y)`   :math:`\min(x, |y|)`
:code:`maxyabs(x,y)`   :math:`\max(x, |y|)`
=====================  ==================================

Non-Generic Functions
=====================

===============  ==================================
Functions        C or math functions
===============  ==================================
:code:`deg`      :code:`math.deg`
:code:`fact`     :code:`mad_num_fact`, :math:`n!`
:code:`fmod`     :code:`math.fmod`
:code:`frexp`    :code:`math.frexp`
:code:`invfact`  :code:`mad_num_invfact`, :math:`1/n!`
:code:`ldexp`    :code:`math.ldexp`
:code:`max`      :code:`math.max`
:code:`min`      :code:`math.min`
:code:`modf`     :code:`math.modf`
:code:`rad`      :code:`math.rad`
===============  ==================================

Random number generators
========================

.. TODO (new chapter?)
