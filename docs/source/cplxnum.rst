.. index::
   Complex numbers

***************
Complex numbers
***************

This chapter describes the :type:`complex` numbers as supported by MAD-NG. The module for complex numbers is not exposed, only the contructors are visible. Thus, complex numbers are handled directly by their methods or by the generic functions of the same name from the module :mod:`MAD.gmath`. Note that :type:`complex` have value semantic like :type:`number`. 

Constructors
============

The constructors for :type:`complex` numbers are directly available from the :mod:`MAD` environment, except for the special case of the imaginary postfix, which is part of the language definition.

.. constant:: i

   The imaginary postfix that qualifies literal numbers as imaginary numbers, i.e. :const:`1i` is the imaginary unit, and :const:`1+2i` is the :type:`complex` number :math:`1+2i`.

.. function:: complex(re, im_)

   Return the :type:`complex` number equivalent to :code:`re + im * 1i`. Default: :code:`im_ = 0`.

.. function:: tocomplex(str)

   Return the :type:`complex` number decoded from the string :var:`str` containing the literal complex number :const:`"a+bi"` (with no spaces) where :var:`a` and :var:`b` are literal numbers, i.e. the strings :const:`"1"`, :const:`"2i"` and :const:`"1+2i"` will give respectively the :type:`complex` numbers :math:`1+0i`, :math:`0+2i` and :math:`1+2i`.

Functions
=========

.. function:: is_complex(a)

   Return :const:`true` if :var:`a` is a :type:`complex` number, :const:`false` otherwise. This function is also available from the module :mod:`MAD.typeid`.

.. function:: is_scalar(a)

   Return :const:`true` if :var:`a` is a :type:`number` or a :type:`complex` number, :const:`false` otherwise. This function is also available from the module :mod:`MAD.typeid`.

Methods
=======

Operator-like Methods
---------------------

=================  ===================   ===================  =============================
Functions          Return values         Metamethods          C functions                         
=================  ===================   ===================  =============================
:func:`z:unm()`    :math:`-z`            :func:`__unm(z,_)`                                
:func:`z:sqr()`    :math:`z \cdot z`     :func:`__mul(z,z)`                                
:func:`z:inv()`    :math:`1 / z`                              :c:func:`mad_cnum_inv_r` [#f1]_                            
:func:`z:add(z2)`  :math:`z + z_2`       :func:`__add(z,z2)`                               
:func:`z:sub(z2)`  :math:`z - z_2`       :func:`__sub(z,z2)`                               
:func:`z:mul(z2)`  :math:`z \cdot z_2`   :func:`__mul(z,z2)`                               
:func:`z:div(z2)`  :math:`z / z_2`       :func:`__div(z,z2)`  :c:func:`mad_cnum_div_r` [#f1]_                             
:func:`z:mod(z2)`  :math:`z \mod z_2`    :func:`__mod(z,z2)`  :c:func:`mad_cnum_mod_r`               
:func:`z:pow(z2)`  :math:`z ^ {z_2}`     :func:`__pow(z,z2)`  :c:func:`mad_cnum_pow_r`                               
:func:`z:eq(z2)`   :math:`z = z_2`       :func:`__eq(z,z2)`                                
=================  ===================   ===================  =============================

Real-like Methods
-----------------

=============================  ===============================  ============================
Functions                      Return values                    C functions
=============================  ===============================  ============================
:func:`z:abs()`                :math:`|z|`                      :c:func:`mad_cnum_abs_r`
:func:`z:acos()`               :math:`\cos^{-1} z`              :c:func:`mad_cnum_acos_r`
:func:`z:acosh()`              :math:`\cosh^{-1} z`             :c:func:`mad_cnum_acosh_r`
:func:`z:acot()`               :math:`\cot^{-1} z`              :c:func:`mad_cnum_atan_r`
:func:`z:acoth()`              :math:`\coth^{-1} z`             :c:func:`mad_cnum_atanh_r`
:func:`z:asin()`               :math:`\sin^{-1} z`              :c:func:`mad_cnum_asin_r`        
:func:`z:asinc()`              :math:`\frac{\sin^{-1} z}{z}`    :c:func:`mad_cnum_asinc_r`
:func:`z:asinh()`              :math:`\sinh^{-1} x`             :c:func:`mad_cnum_asinh_r`
:func:`z:asinhc()`             :math:`\frac{\sinh^{-1} z}{z}`   :c:func:`mad_cnum_asinhc_r`
:func:`z:atan()`               :math:`\tan^{-1} z`              :c:func:`mad_cnum_atan_r`        
:func:`z:atan2(z2)`            :math:`\tan^{-1} \frac{z}{z_2}`  :c:func:`mad_cnum_atan2_r`
:func:`z:atanh()`              :math:`\tanh^{-1} z`             :c:func:`mad_cnum_atanh_r`
:func:`z:ceil()`               :math:`\operatorname{ceil}(z)`   
:func:`z:cos()`                :math:`\cos z`                   :c:func:`mad_cnum_cos_r`   
:func:`z:cosh()`               :math:`\cosh z`                  :c:func:`mad_cnum_cosh_r`
:func:`z:cot()`                :math:`\cot z`                   :c:func:`mad_cnum_tan_r`
:func:`z:coth()`               :math:`\coth z`                  :c:func:`mad_cnum_tanh_r`
:func:`z:exp()`                :math:`\exp z`                   :c:func:`mad_cnum_exp_r`
:func:`z:floor()`              :math:`\operatorname{floor}(z)`     
:func:`z:frac()`               :math:`\operatorname{frac}(z)`                
:func:`z:hypot(z2)`            :math:`\sqrt{z^2+z_2^2}`         [#f2]_         
:func:`z:hypot3(z2,z3)`        :math:`\sqrt{z^2+z_2^2+z_3^2}`   [#f2]_  
:func:`z:inv(v_)` [#f3]_       :math:`\frac{v}{z}`              :c:func:`mad_cnum_inv_r` [#f1]_             
:func:`z:invsqrt(v_)` [#f3]_   :math:`\frac{v}{\sqrt z}`        :c:func:`mad_cnum_invsqrt_r` [#f1]_              
:func:`z:log()`                :math:`\log z`                   :c:func:`mad_cnum_log_r`
:func:`z:log10()`              :math:`\log_{10} z`              :c:func:`mad_cnum_log10_r`
:func:`z:pow(z2)`              :math:`z^{z_2}`                  :c:func:`mad_cnum_pow_r`  
:func:`z:powi(n)`              :math:`z^n`                      :c:func:`mad_cnum_powi_r`
:func:`z:round()`              :math:`\operatorname{round}(z)`  
:func:`z:sin()`                :math:`\sin z`                   :c:func:`mad_cnum_sin_r`   
:func:`z:sinc()`               :math:`\frac{\sin z}{z}`         :c:func:`mad_cnum_sinc_r`
:func:`z:sinh()`               :math:`\sinh z`                  :c:func:`mad_cnum_sinh_r`    
:func:`z:sinhc()`              :math:`\frac{\sinh z}{z}`        :c:func:`mad_cnum_sinhc_r`
:func:`z:sqrt()`               :math:`\sqrt{z}`                 :c:func:`mad_cnum_sqrt_r`     
:func:`z:tan()`                :math:`\tan z`                   :c:func:`mad_cnum_tan_r`
:func:`z:tanh()`               :math:`\tanh z`                  :c:func:`mad_cnum_tanh_r`
:func:`z:trunc()`              :math:`\operatorname{trunc}(z)`                      
:func:`z:unit()`               :math:`\frac{z}{|z|}`            :c:func:`mad_cnum_unit_r`
=============================  ===============================  ============================

Complex-like Methods
--------------------

=================  ==============================================  ==========================
Functions          Return values                                   C functions
=================  ==============================================  ==========================
:func:`z:cabs()`   :math:`|z|`                                     :c:func:`mad_cnum_abs_r`
:func:`z:carg()`   :math:`\arg z`                                  :c:func:`mad_cnum_arg_r`   
:func:`z:conj()`   :math:`z^*`                                     
:func:`z:imag()`   :math:`\Im(z)`                                     
:func:`z:polar()`  :math:`|z|\,e^{i \arg z}`                       :c:func:`mad_cnum_polar_r`              
:func:`z:proj()`   :math:`\operatorname{proj}(z)`                  :c:func:`mad_cnum_proj_r`                   
:func:`z:real()`   :math:`\Re(z)`                                     
:func:`z:rect()`   :math:`\Re(z)\cos \Im(z)+i\,\Re(z)\sin \Im(z)`  :c:func:`mad_cnum_rect_r`                                   
:func:`z:reim()`   :math:`\Re(z), \Im(z)`                                     
=================  ==============================================  ==========================

Error-like Methods
------------------

Error-like methods call C wrappers to the corresponding functions from the `Faddeeva library <http://ab-initio.mit.edu/wiki/index.php/Faddeeva_Package>`_ from the MIT, considered as one of the most accurate and fast implementation over the complex plane [FADDEEVA]_ (see :file:`mad_num.c`).

=======================  ==========================================================  ======================
Functions                Return values                                               C functions  
=======================  ==========================================================  ======================
:func:`z:erf(rtol_)`     :math:`\frac{2}{\sqrt\pi}\int_0^z e^{-t^2} dt`              :c:func:`mad_cnum_erf_r`      
:func:`z:erfc(rtol_)`    :math:`1-\operatorname{erf}(z)`                             :c:func:`mad_cnum_erfc_r`     
:func:`z:erfi(rtol_)`    :math:`-i\operatorname{erf}(i z)`                           :c:func:`mad_cnum_erfi_r`     
:func:`z:erfcx(rtol_)`   :math:`e^{z^2}\operatorname{erfc}(z)`                       :c:func:`mad_cnum_erfcx_r`    
:func:`z:wf(rtol_)`      :math:`e^{-z^2}\operatorname{erfc}(-i z)`                   :c:func:`mad_cnum_wf_r`       
:func:`z:dawson(rtol_)`  :math:`\frac{-i\sqrt\pi}{2}e^{-z^2}\operatorname{erf}(iz)`  :c:func:`mad_cnum_dawson_r`
=======================  ==========================================================  ======================

References
==========

.. [CPXDIV] R. L. Smith, *"Algorithm 116: Complex division"*, Commun. ACM, 5(8):435, 1962.

.. [CPXDIV2] M. Baudin and R. L. Smith, *"A robust complex division in Scilab"*, October 2012. http://arxiv.org/abs/1210.4539.

.. [FADDEEVA] A. Oeftiger, R. De Maria, L. Deniau et al, *"Review of CPU and GPU Faddeeva Implementations"*, IPAC2016. https://cds.cern.ch/record/2207430/files/wepoy044.pdf

.. ---------------------------------------

.. rubric:: Footnotes

.. [#f1] Division and inverse use a robust and fast complex division algorithm, see [CPXDIV]_ and [CPXDIV2]_ for details. 
.. [#f2] Hypot and hypot3 methods use a trivial implementation that may lead to numerical overflow/underflow.
.. [#f3] Default: :code:`v_ = 1`. 

