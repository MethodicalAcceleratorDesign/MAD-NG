.. index::
   Complex numbers

***************
Complex Numbers
***************

This chapter describes the :type:`complex` numbers as supported by MAD-NG. The module for `Complex numbers <https://en.wikipedia.org/wiki/Complex_number>`_ is not exposed, only the contructors are visible from the :mod:`MAD` environment and thus, complex numbers are handled directly by their methods or by the generic functions of the same name from the module :mod:`MAD.gmath`. Note that :type:`complex` have value semantic like a pair of :type:`number` equivalent to a C structure or an array :expr:`const num_t[2]` for direct compliance with the C API.

Types promotion
===============

The operations on complex numbers may involve other data types like real numbers leading to few combinations of types. In order to simplify the descriptions, the generic names :var:`num` and :var:`cpx` are used for real and complex numbers respectively. The following table summarizes all valid combinations of types for binary operations involving at least one :type:`complex` type:

=================  ==================  ===============
Left Operand Type  Right Operand Type  Result Type
=================  ==================  ===============
:type:`number`     :type:`complex`     :type:`complex`
:type:`complex`    :type:`number`      :type:`complex`
:type:`complex`    :type:`complex`     :type:`complex`
=================  ==================  ===============

Constructors
============

The constructors for :type:`complex` numbers are directly available from the :mod:`MAD` environment, except for the special case of the imaginary postfix, which is part of the language definition.

.. constant:: i

   The imaginary postfix that qualifies literal numbers as imaginary numbers, i.e. :const:`1i` is the imaginary unit, and :const:`1+2i` is the :type:`complex` number :math:`1+2i`.

.. function:: complex(re_, im_)

   Return the :type:`complex` number equivalent to :expr:`re + im * 1i`. Default: :expr:`re_ = 0`, :expr:`im_ = 0`.

.. function:: tocomplex(str)

   Return the :type:`complex` number decoded from the string :var:`str` containing the literal complex number :const:`"a+bi"` (with no spaces) where :var:`a` and :var:`b` are literal numbers, i.e. the strings :const:`"1"`, :const:`"2i"` and :const:`"1+2i"` will give respectively the :type:`complex` numbers :math:`1+0i`, :math:`0+2i` and :math:`1+2i`.

Attributes
==========

.. constant:: cpx.re

   The real part of the :type:`complex` number :var:`cpx`.

.. constant:: cpx.im

   The imaginary part of the :type:`complex` number :var:`cpx`.

Functions
=========

.. function:: is_complex(a)

   Return :const:`true` if :var:`a` is a :type:`complex` number, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

.. function:: is_scalar(a)

   Return :const:`true` if :var:`a` is a :type:`number` or a :type:`complex` number, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

Operator-like Methods
---------------------

=================  ===================   ===================  =============================
Functions          Return values         Metamethods          C functions                         
=================  ===================   ===================  =============================
:func:`z:unm()`    :math:`-z`            :func:`__unm(z,_)`                                
:func:`z:add(z2)`  :math:`z + z_2`       :func:`__add(z,z2)`                               
:func:`z:sub(z2)`  :math:`z - z_2`       :func:`__sub(z,z2)`                               
:func:`z:mul(z2)`  :math:`z \cdot z_2`   :func:`__mul(z,z2)`                               
:func:`z:div(z2)`  :math:`z / z_2`       :func:`__div(z,z2)`  :c:func:`mad_cpx_div_r` [#f1]_
:func:`z:mod(z2)`  :math:`z \mod z_2`    :func:`__mod(z,z2)`  :c:func:`mad_cpx_mod_r`
:func:`z:pow(z2)`  :math:`z ^ {z_2}`     :func:`__pow(z,z2)`  :c:func:`mad_cpx_pow_r`
:func:`z:eq(z2)`   :math:`z = z_2`       :func:`__eq(z,z2)`                                
=================  ===================   ===================  =============================

Real-like Methods
-----------------

=============================  ====================================================================  ============================
Functions                      Return values                                                         C functions
=============================  ====================================================================  ============================
:func:`z:abs()`                :math:`|z|`                                                           :c:func:`mad_cpx_abs_r`
:func:`z:acos()`               :math:`\cos^{-1} z`                                                   :c:func:`mad_cpx_acos_r`
:func:`z:acosh()`              :math:`\cosh^{-1} z`                                                  :c:func:`mad_cpx_acosh_r`
:func:`z:acot()`               :math:`\cot^{-1} z`                                                   :c:func:`mad_cpx_atan_r`
:func:`z:acoth()`              :math:`\coth^{-1} z`                                                  :c:func:`mad_cpx_atanh_r`
:func:`z:asin()`               :math:`\sin^{-1} z`                                                   :c:func:`mad_cpx_asin_r`
:func:`z:asinc()`              :math:`\frac{\sin^{-1} z}{z}`                                         :c:func:`mad_cpx_asinc_r`
:func:`z:asinh()`              :math:`\sinh^{-1} x`                                                  :c:func:`mad_cpx_asinh_r`
:func:`z:asinhc()`             :math:`\frac{\sinh^{-1} z}{z}`                                        :c:func:`mad_cpx_asinhc_r`
:func:`z:atan()`               :math:`\tan^{-1} z`                                                   :c:func:`mad_cpx_atan_r`
:func:`z:atanh()`              :math:`\tanh^{-1} z`                                                  :c:func:`mad_cpx_atanh_r`
:func:`z:ceil()`               :math:`\lceil\Re(z)\rceil+i\,\lceil\Im(z)\rceil`        
:func:`z:cos()`                :math:`\cos z`                                                        :c:func:`mad_cpx_cos_r`
:func:`z:cosh()`               :math:`\cosh z`                                                       :c:func:`mad_cpx_cosh_r`
:func:`z:cot()`                :math:`\cot z`                                                        :c:func:`mad_cpx_tan_r`
:func:`z:coth()`               :math:`\coth z`                                                       :c:func:`mad_cpx_tanh_r`
:func:`z:exp()`                :math:`\exp z`                                                        :c:func:`mad_cpx_exp_r`
:func:`z:floor()`              :math:`\lfloor\Re(z)\rfloor+i\,\lfloor\Im(z)\rfloor`      
:func:`z:frac()`               :math:`z - \operatorname{trunc}(z)`           
:func:`z:hypot(z2)`            :math:`\sqrt{z^2+z_2^2}`                                              [#f2]_         
:func:`z:hypot3(z2,z3)`        :math:`\sqrt{z^2+z_2^2+z_3^2}`                                        [#f2]_  
:func:`z:inv(v_)`              :math:`\frac{v}{z}`                                                   :c:func:`mad_cpx_inv_r` [#f1]_
:func:`z:invsqrt(v_)`          :math:`\frac{v}{\sqrt z}`                                             :c:func:`mad_cpx_invsqrt_r` [#f1]_
:func:`z:log()`                :math:`\log z`                                                        :c:func:`mad_cpx_log_r`
:func:`z:log10()`              :math:`\log_{10} z`                                                   :c:func:`mad_cpx_log10_r`
:func:`z:powi(n)`              :math:`z^n`                                                           :c:func:`mad_cpx_powi_r`
:func:`z:round()`              :math:`\lfloor\Re(z)\rceil+i\,\lfloor\Im(z)\rceil`     
:func:`z:sin()`                :math:`\sin z`                                                        :c:func:`mad_cpx_sin_r`
:func:`z:sinc()`               :math:`\frac{\sin z}{z}`                                              :c:func:`mad_cpx_sinc_r`
:func:`z:sinh()`               :math:`\sinh z`                                                       :c:func:`mad_cpx_sinh_r`
:func:`z:sinhc()`              :math:`\frac{\sinh z}{z}`                                             :c:func:`mad_cpx_sinhc_r`
:func:`z:sqr()`                :math:`z \cdot z`                                                                                     
:func:`z:sqrt()`               :math:`\sqrt{z}`                                                      :c:func:`mad_cpx_sqrt_r`
:func:`z:tan()`                :math:`\tan z`                                                        :c:func:`mad_cpx_tan_r`
:func:`z:tanh()`               :math:`\tanh z`                                                       :c:func:`mad_cpx_tanh_r`
:func:`z:trunc()`              :math:`\operatorname{trunc} \Re(z)+i\,\operatorname{trunc} \Im(z)`                                
:func:`z:unit()`               :math:`\frac{z}{|z|}`                                                 :c:func:`mad_cpx_unit_r`
=============================  ====================================================================  ============================

In methods :func:`inv()` and :func:`invsqrt()`, default is :expr:`v_ = 1`.

Complex-like Methods
--------------------

=================  ==============================================  ==========================
Functions          Return values                                   C functions
=================  ==============================================  ==========================
:func:`z:cabs()`   :math:`|z|`                                     :c:func:`mad_cpx_abs_r`
:func:`z:carg()`   :math:`\arg z`                                  :c:func:`mad_cpx_arg_r`
:func:`z:conj()`   :math:`z^*`                                     
:func:`z:fabs()`   :math:`|\Re(z)|+i\,|\Im(z)|`
:func:`z:imag()`   :math:`\Im(z)`                                     
:func:`z:polar()`  :math:`|z|\,e^{i \arg z}`                       :c:func:`mad_cpx_polar_r`
:func:`z:proj()`   :math:`\operatorname{proj}(z)`                  :c:func:`mad_cpx_proj_r`
:func:`z:real()`   :math:`\Re(z)`                                     
:func:`z:rect()`   :math:`\Re(z)\cos \Im(z)+i\,\Re(z)\sin \Im(z)`  :c:func:`mad_cpx_rect_r`
:func:`z:reim()`   :math:`\Re(z), \Im(z)`                                     
=================  ==============================================  ==========================

Error-like Methods
------------------

Error-like methods call C wrappers to the corresponding functions from the `Faddeeva library <http://ab-initio.mit.edu/Faddeeva>`_ from the MIT, considered as one of the most accurate and fast implementation over the complex plane [FADDEEVA]_ (see :file:`mad_num.c`).

=======================  ==========================================================  ======================
Functions                Return values                                               C functions  
=======================  ==========================================================  ======================
:func:`z:erf(rtol_)`     :math:`\frac{2}{\sqrt\pi}\int_0^z e^{-t^2} dt`              :c:func:`mad_cpx_erf_r`
:func:`z:erfc(rtol_)`    :math:`1-\operatorname{erf}(z)`                             :c:func:`mad_cpx_erfc_r`
:func:`z:erfi(rtol_)`    :math:`-i\operatorname{erf}(i z)`                           :c:func:`mad_cpx_erfi_r`
:func:`z:erfcx(rtol_)`   :math:`e^{z^2}\operatorname{erfc}(z)`                       :c:func:`mad_cpx_erfcx_r`
:func:`z:wf(rtol_)`      :math:`e^{-z^2}\operatorname{erfc}(-i z)`                   :c:func:`mad_cpx_wf_r`
:func:`z:dawson(rtol_)`  :math:`\frac{-i\sqrt\pi}{2}e^{-z^2}\operatorname{erf}(iz)`  :c:func:`mad_cpx_dawson_r`
=======================  ==========================================================  ======================

Operators
=========

The operators on :type:`complex` follow the conventional mathematical operations of `Complex numbers <https://en.wikipedia.org/wiki/Complex_number#Relations_and_operations>`__.

.. function:: -cpx

   Return a :type:`complex` resulting from the negation of the operand as computed by :func:`cpx:unm()`.

.. function:: num + cpx
              cpx + num
              cpx + cpx

   Return a :type:`complex` resulting from the sum of the left and right operands as computed by :func:`cpx:add()`.

.. function:: num - cpx
              cpx - num
              cpx - cpx

   Return a :type:`complex` resulting from the difference of the left and right operands as computed by :func:`cpx:sub()`.

.. function:: num * cpx
              cpx * num
              cpx * cpx

   Return a :type:`complex` resulting from the product of the left and right operands as computed by :func:`cpx:mul()`.

.. function:: num / cpx
              cpx / num
              cpx / cpx

   Return a :type:`complex` resulting from the division of the left and right operands as computed by :func:`cpx:div()`. If the right operand is a complex number, the division uses a robut and fast algorithm implemented in :c:func:`mad_cpx_div_r` [#f1]_.

.. function:: num % cpx
              cpx % num
              cpx % cpx

   Return a :type:`complex` resulting from the rest of the division of the left and right operands, i.e. :math:`x - y \lfloor \frac{x}{y} \rfloor`,  as computed by :func:`cpx:mod()`. If the right operand is a complex number, the division uses a robut and fast algorithm implemented in :c:func:`mad_cpx_div_r` [#f1]_.

.. function:: num ^ cpx
              cpx ^ num
              cpx ^ cpx

   Return a :type:`complex` resulting from the left operand raised to the power of the right operand as computed by :func:`cpx:pow()`.

.. function:: num == cpx
              cpx == num
              cpx == cpx

   Return :const:`false` if the real or the imaginary part differ between the left and right operands, :const:`true` otherwise. A number :var:`a` will be interpreted as :math:`a+i0` for the comparison.

C API
=====

These functions are provided for performance reason and compliance (i.e. branch cut) with the C API of other modules dealing with complex numbers like the linear and the differential algebra. For the same reason, some functions hereafter refer to the section 7.3 of the C Programming Language Standard [ISOC99CPX]_.

.. c:function:: num_t mad_cpx_abs_r (num_t x_re, num_t x_im)

   Return the modulus of the :type:`complex` :var:`x` as computed by :c:func:`cabs()`.

.. c:function:: num_t mad_cpx_arg_r (num_t x_re, num_t x_im)

   Return the argument in :math:`[-\pi, +\pi]` of the :type:`complex` :var:`x` as computed by :c:func:`carg()`.

.. c:function:: void  mad_cpx_unit_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the :type:`complex` :var:`x` divided by its modulus as computed by :c:func:`cabs()`.

.. c:function:: void  mad_cpx_proj_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the projection of the :type:`complex` :var:`x` on the Riemann sphere as computed by :c:func:`cproj()`.

.. c:function:: void  mad_cpx_rect_r (num_t  rho, num_t  ang, cpx_t *r)

   Put in :var:`r` the rectangular form of the :type:`complex` :expr:`rho * exp(i*ang)`.

.. c:function:: void  mad_cpx_polar_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the polar form of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_inv_r (num_t x_re, num_t x_im, cpx_t *r)
                cpx_t mad_cpx_inv (cpx_t x)

   Put in :var:`r` or return the inverse of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_invsqrt_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the square root of the inverse of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_div_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cpx_t *r)
                cpx_t mad_cpx_div (cpx_t x, cpx_t y)

   Put in :var:`r` or return the :type:`complex` :var:`x` divided by the :type:`complex` :var:`y`.

.. c:function:: void  mad_cpx_mod_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cpx_t *r)

   Put in :var:`r` the remainder of the division of the :type:`complex` :var:`x` by the :type:`complex` :var:`y`.

.. c:function:: void  mad_cpx_pow_r (num_t x_re, num_t x_im, num_t y_re, num_t y_im, cpx_t *r)

   Put in :var:`r` the :type:`complex` :var:`x` raised to the power of :type:`complex` :var:`y` using :c:func:`cpow()`.

.. c:function:: void  mad_cpx_powi_r (num_t x_re, num_t x_im, int n, cpx_t *r)
                cpx_t mad_cpx_powi (cpx_t x, int n)

   Put in :var:`r` or return the :type:`complex` :var:`x` raised to the power of the :type:`integer` :var:`n` using a fast algorithm.

.. c:function:: void  mad_cpx_sqrt_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the square root of the :type:`complex` :var:`x` as computed by :c:func:`csqrt()`.

.. c:function:: void  mad_cpx_exp_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the exponential of the :type:`complex` :var:`x` as computed by :c:func:`cexp()`.

.. c:function:: void  mad_cpx_log_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the natural logarithm of the :type:`complex` :var:`x` as computed by :c:func:`clog()`.

.. c:function:: void  mad_cpx_log10_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the logarithm of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_sin_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the sine of the :type:`complex` :var:`x` as computed by :c:func:`csin()`.

.. c:function:: void  mad_cpx_cos_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the cosine of the :type:`complex` :var:`x` as computed by :c:func:`ccos()`.

.. c:function:: void  mad_cpx_tan_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the tangent of the :type:`complex` :var:`x` as computed by :c:func:`ctan()`.

.. c:function:: void  mad_cpx_sinh_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the hyperbolic sine of the :type:`complex` :var:`x` as computed by :c:func:`csinh()`.

.. c:function:: void  mad_cpx_cosh_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the hyperbolic cosine of the :type:`complex` :var:`x` as computed by :c:func:`ccosh()`.

.. c:function:: void  mad_cpx_tanh_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the hyperbolic tangent of the :type:`complex` :var:`x` as computed by :c:func:`ctanh()`.

.. c:function:: void  mad_cpx_asin_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the arc sine of the :type:`complex` :var:`x` as computed by :c:func:`casin()`.

.. c:function:: void  mad_cpx_acos_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the arc cosine of the :type:`complex` :var:`x` as computed by :c:func:`cacos()`.

.. c:function:: void  mad_cpx_atan_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the arc tangent of the :type:`complex` :var:`x` as computed by :c:func:`catan()`.

.. c:function:: void  mad_cpx_asinh_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the hyperbolic arc sine of the :type:`complex` :var:`x` as computed by :c:func:`casinh()`.

.. c:function:: void  mad_cpx_acosh_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the hyperbolic arc cosine of the :type:`complex` :var:`x` as computed by :c:func:`cacosh()`.

.. c:function:: void  mad_cpx_atanh_r (num_t x_re, num_t x_im, cpx_t *r)

   Put in :var:`r` the hyperbolic arc tangent of the :type:`complex` :var:`x` as computed by :c:func:`catanh()`.

.. c:function:: void  mad_cpx_sinc_r (num_t x_re, num_t x_im, cpx_t *r)
                cpx_t mad_cpx_sinc (cpx_t x)

   Put in :var:`r` or return the sine cardinal of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_sinhc_r (num_t x_re, num_t x_im, cpx_t *r)
                cpx_t mad_cpx_sinhc (cpx_t x)

   Put in :var:`r` or return the hyperbolic sine cardinal of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_asinc_r (num_t x_re, num_t x_im, cpx_t *r)
                cpx_t mad_cpx_asinc (cpx_t x)

   Put in :var:`r` or return the arc sine cardinal of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_asinhc_r (num_t x_re, num_t x_im, cpx_t *r)
                cpx_t mad_cpx_asinhc (cpx_t x)

   Put in :var:`r` or return the hyperbolic arc sine cardinal of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_wf_r (num_t x_re, num_t x_im, num_t relerr, cpx_t *r)
                cpx_t mad_cpx_wf (cpx_t x, num_t relerr)

   Put in :var:`r` or return the Faddeeva function of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_erf_r (num_t x_re, num_t x_im, num_t relerr, cpx_t *r)
                cpx_t mad_cpx_erf (cpx_t x, num_t relerr)

   Put in :var:`r` or return the error function of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_erfc_r (num_t x_re, num_t x_im, num_t relerr, cpx_t *r)
                cpx_t mad_cpx_erfc (cpx_t x, num_t relerr)

   Put in :var:`r` or return the complementary error function of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_erfcx_r (num_t x_re, num_t x_im, num_t relerr, cpx_t *r)
                cpx_t mad_cpx_erfcx (cpx_t x, num_t relerr)

   Put in :var:`r` or return the scaled complementary error function of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_erfi_r (num_t x_re, num_t x_im, num_t relerr, cpx_t *r)
                cpx_t mad_cpx_erfi (cpx_t x, num_t relerr)

   Put in :var:`r` or return the imaginary error function of the :type:`complex` :var:`x`.

.. c:function:: void  mad_cpx_dawson_r (num_t x_re, num_t x_im, num_t relerr, cpx_t *r)
                cpx_t mad_cpx_dawson (cpx_t x, num_t relerr)

   Put in :var:`r` or return the Dawson integral for the :type:`complex` :var:`x`.

References
==========

.. [CPXDIV] R. L. Smith, *"Algorithm 116: Complex division"*, Commun. ACM, 5(8):435, 1962.

.. [CPXDIV2] M. Baudin and R. L. Smith, *"A robust complex division in Scilab"*, October 2012. http://arxiv.org/abs/1210.4539.

.. [FADDEEVA] A. Oeftiger, R. De Maria, L. Deniau et al, *"Review of CPU and GPU Faddeeva Implementations"*, IPAC2016. https://cds.cern.ch/record/2207430/files/wepoy044.pdf.

.. [ISOC99CPX] ISO/IEC 9899:1999 Programming Languages - C. https://www.iso.org/standard/29237.html.

.. ---------------------------------------

.. rubric:: Footnotes

.. [#f1] Division and inverse use a robust and fast complex division algorithm, see [CPXDIV]_ and [CPXDIV2]_ for details. 
.. [#f2] Hypot and hypot3 methods use a trivial implementation that may lead to numerical overflow/underflow.

