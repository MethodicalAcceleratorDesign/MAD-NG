.. index::
   functions utilities

*******************
Functions Utilities
*******************

Operators as Functions
======================

The module :mod:`MAD.gfunc` provides many functions that are named version of operators and useful when operators cannot be used directly, e.g. when passed as argument or to compose together. These functions can also be retrieved from the module :mod:`MAD.gfunc.opstr` by their associated string (if available).

Math Operators
--------------

Functions for math operators are wrappers to associated mathematical operators, which themselves can be overridden by their associated metamethods.

====================  ======================  ===============  =======================
Functions             Return values           Operator string  Metamethods
====================  ======================  ===============  =======================
:func:`unm(x)`        :math:`-x`              :const:`"~"`     :func:`__unm(x,_)`
:func:`add(x,y)`      :math:`x + y`           :const:`"+"`     :func:`__add(x,y)`
:func:`sub(x,y)`      :math:`x - y`           :const:`"-"`     :func:`__sub(x,y)`
:func:`mul(x,y)`      :math:`x * y`           :const:`"*"`     :func:`__mul(x,y)`
:func:`div(x,y)`      :math:`x / y`           :const:`"/"`     :func:`__div(x,y)`
:func:`mod(x,y)`      :math:`x\,\%\,y`        :const:`"%"`     :func:`__mod(x,y)`
:func:`pow(x,y)`      :math:`x ^ y`           :const:`"^"`     :func:`__pow(x,y)`
:func:`sqr(x)`        :math:`x * x`           :const:`"^2"`    :func:`__mul(x,x)`
:func:`inv(x)`        :math:`1 / x`           :const:`"1/"`    :func:`__div(1,x)`
:func:`emul(x,y,r_)`  :math:`x\,.*\,y`        :const:`".*"`    :func:`__emul(x,y,r_)`
:func:`ediv(x,y,r_)`  :math:`x\,./\,y`        :const:`"./"`    :func:`__ediv(x,y,r_)`
:func:`emod(x,y,r_)`  :math:`x\,.\%\,y`       :const:`".%"`    :func:`__emod(x,y,r_)`
:func:`epow(x,y,r_)`  :math:`x\,.\hat\ \ y`   :const:`".^"`    :func:`__epow(x,y,r_)`
====================  ======================  ===============  =======================

Logical Operators
-----------------

Functions for logical operators are wrappers to associated logical operators, which themselves can be overridden by their associated metamethods (if any).

====================  ======================  ==============================  =================
Functions             Return values           Operator string                 Metamethods
====================  ======================  ==============================  =================
:func:`lfalse()`      :const:`true`                                           
:func:`ltrue()`       :const:`false`                                          
:func:`lnot(x)`       :math:`\lnot x`         :const:`"!"`                      
:func:`lbool(x)`      :math:`\lnot\lnot x`    :const:`"!!"`                       
:func:`land(x,y)`     :math:`x \land y`       :const:`"&&"`                       
:func:`lor(x,y)`      :math:`x \lor y`        :const:`"||"`                       
:func:`eq(x,y)`       :math:`x == y`          :const:`"=="`                   :func:`__eq(x,y)`
:func:`ne(x,y)`       :math:`x \neq y`        :const:`"!="` or :const:`"~="`  :func:`__eq(x,y)`
:func:`lt(x,y)`       :math:`x < y`           :const:`"<"`                    :func:`__lt(x,y)`
:func:`le(x,y)`       :math:`x <= y`          :const:`"<="`                   :func:`__le(x,y)`
:func:`gt(x,y)`       :math:`x > y`           :const:`">"`                    :func:`__le(x,y)`
:func:`ge(x,y)`       :math:`x >= y`          :const:`">="`                   :func:`__lt(x,y)`
====================  ======================  ==============================  =================

Bitwise Functions
=================

Functions for bitwise operations are those from the LuaJIT module :mod:`bit` and imported into the module :mod:`MAD.gfunc` for convenience, see http://bitop.luajit.org/api.html for details. Note that all these functions have *value semantic* and normalise their arguments to the numeric range of a 32 bit integer before use.

=====================  ====================================================
Functions              Return values         
=====================  ====================================================
:func:`tobit(x)`       Return the normalized value of :var:`x` to the range of a 32 bit integer      
:func:`tohex(x,n_)`    Return the hex string of :var:`x` with :var:`n` digits (:math:`n<0` use caps)    
:func:`bnot(x)`        Return the bitwise reverse of :var:`x` bits    
:func:`band(x,...)`    Return the bitwise *AND* of all arguments     
:func:`bor(x,...)`     Return the bitwise *OR* of all arguments 
:func:`bxor(x,...)`    Return the bitwise *XOR* of all arguments
:func:`lshift(x,n)`    Return the bitwise left shift of :var:`x` by :var:`n` bits with 0-bit shift-in     
:func:`rshift(x,n)`    Return the bitwise right shift of :var:`x` by :var:`n` bits with 0-bit shift-in
:func:`arshift(x,n)`   Return the bitwise right shift of :var:`x` by :var:`n` bits with sign bit shift-in
:func:`rol(x,n)`       Return the bitwise left rotation of :var:`x` by :var:`n` bits      
:func:`ror(x,n)`       Return the bitwise right rotation of :var:`x` by :var:`n` bits     
:func:`bswap(x)`       Return the swapped bytes of :var:`x`, i.e. convert big endian to/from little endian       
=====================  ====================================================

Flags Functions
===============

A flag is 32 bit unsigned integer used to store up to 32 binary states with the convention that :const:`0` means disabled/cleared and :const:`1` means enabled/set. Functions on flags are useful aliases to -- or combinaison of -- bitwise operations to manipulate their states (i.e. their bits).

=====================  ====================================================
Functions              Return values         
=====================  ====================================================
:func:`bset(x,n)`      Return the flag :var:`x` with state :var:`n` enabled
:func:`bclr(x,n)`      Return the flag :var:`x` with state :var:`n` disabled   
:func:`btst(x,n)`      Return :const:`true` if state :var:`n` is enabled in :var:`x`, :const:`false` otherwise      
:func:`fbit(n)`        Return a flag with only state :var:`n` enabled    
:func:`fnot(x)`        Return the flag :var:`x` with all states flipped
:func:`fset(x,...)`    Return the flag :var:`x` with disabled states flipped if enabled in any flag passed as argument
:func:`fcut(x,...)`    Return the flag :var:`x` with enabled states flipped if disabled in any flag passed as argument 
:func:`fclr(x,f)`      Return the flag :var:`x` with enabled states flipped if enabled in :var:`f`
:func:`ftst(x,f)`      Return :const:`true` if all states enabled in :var:`f` are enabled in :var:`x`, :const:`false` otherwise 
:func:`fall(x)`        Return :const:`true` if all states are enabled in :var:`x`, :const:`false` otherwise       
:func:`fany(x)`        Return :const:`true` if any state is enabled in :var:`x`, :const:`false` otherwise    
=====================  ====================================================

Object Operators
================
