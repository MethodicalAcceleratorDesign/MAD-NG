***********************
Utility Functions
***********************

Operators as Functions
======================

Generic Math-like Functions
---------------------------

Generic Math-like functions are wrappers to associated operators, which themselves can be overridden by their associated metamethods.

====================  ======================  =============
Operators             Return values           Metamethods
====================  ======================  =============
:code:`unm(x)`        :math:`-x`              :func:`__unm`
:code:`add(x,y)`      :math:`x + y`           :func:`__add`
:code:`sub(x,y)`      :math:`x - y`           :func:`__sub`
:code:`mul(x,y)`      :math:`x * y`           :func:`__mul`
:code:`div(x,y)`      :math:`x / y`           :func:`__div`
:code:`mod(x,y)`      :math:`x\,\%\,y`        :func:`__mod`
:code:`pow(x,y)`      :math:`x ^ y`           :func:`__pow`
:code:`sqr(x)`        :math:`x * x`           :func:`__mul`
:code:`inv(x)`        :math:`1 / x`           :func:`__div`
:code:`emul(x,y,r_)`  :math:`x\,.*\,y`        :func:`__emul`
:code:`ediv(x,y,r_)`  :math:`x\,./\,y`        :func:`__ediv`
:code:`emod(x,y,r_)`  :math:`x\,.\%\,y`       :func:`__emod`
:code:`epow(x,y,r_)`  :math:`x\,.\hat\ \ y`   :func:`__epow`
====================  ======================  =============
