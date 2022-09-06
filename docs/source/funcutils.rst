.. index::
   functions utilities

********************
Functional Utilities
********************

This chapter describes useful functions provided by the module :mod:`MAD.gfunc` to help dealing with operators as functions and to manipulate functions in a `functional <https://en.wikipedia.org/wiki/Functional_programming>`_ way [#]_. It also provide a complete set of functions to create, combine and use *functors*, i.e. objects that behave like functions with :type:`callable` semantic. Functors are mainly used by the object model to distinguish from functions that are interpreted as deferred expressions and automatically evaluated on read (see Deferred Expressions), and by the tracking codes Survey and Track to deal with (user-defined) actions. 

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

Operators as Functions
======================

The module :mod:`MAD.gfunc` provides many functions that are named version of operators and useful when operators cannot be used directly, e.g. when passed as argument or to compose together. These functions can also be retrieved from the module :mod:`MAD.gfunc.opstr` by their associated string (if available).

Math Operators
--------------

Functions for math operators are wrappers to associated mathematical operators, which themselves can be overridden by their associated metamethods.

====================  =====================  ===============  =======================
Functions             Return values          Operator string  Metamethods
====================  =====================  ===============  =======================
:func:`unm(x)`        :math:`-x`             :const:`"~"`     :func:`__unm(x,_)`
:func:`inv(x)`        :math:`1 / x`          :const:`"1/"`    :func:`__div(1,x)`
:func:`sqr(x)`        :math:`x * x`          :const:`"^2"`    :func:`__mul(x,x)`
:func:`add(x,y)`      :math:`x + y`          :const:`"+"`     :func:`__add(x,y)`
:func:`sub(x,y)`      :math:`x - y`          :const:`"-"`     :func:`__sub(x,y)`
:func:`mul(x,y)`      :math:`x * y`          :const:`"*"`     :func:`__mul(x,y)`
:func:`div(x,y)`      :math:`x / y`          :const:`"/"`     :func:`__div(x,y)`
:func:`mod(x,y)`      :math:`x\,\%\,y`       :const:`"%"`     :func:`__mod(x,y)`
:func:`pow(x,y)`      :math:`x ^ y`          :const:`"^"`     :func:`__pow(x,y)`
:func:`emul(x,y,r_)`  :math:`x\,.*\,y`       :const:`".*"`    :func:`__emul(x,y,r_)` [#]_
:func:`ediv(x,y,r_)`  :math:`x\,./\,y`       :const:`"./"`    :func:`__ediv(x,y,r_)`
:func:`emod(x,y,r_)`  :math:`x\,.\%\,y`      :const:`".%"`    :func:`__emod(x,y,r_)`
:func:`epow(x,y,r_)`  :math:`x\,.\hat\ \ y`  :const:`".^"`    :func:`__epow(x,y,r_)`
====================  =====================  ===============  =======================

Logical Operators
-----------------

Functions for logical operators are wrappers to associated logical operators, which themselves can be overridden by their associated metamethods (if any).

=================  ====================  ==============================  =================
Functions          Return values         Operator string                 Metamethods
=================  ====================  ==============================  =================
:func:`lfalse()`   :const:`true`                                         
:func:`ltrue()`    :const:`false`                                        
:func:`lnot(x)`    :math:`\lnot x`       :const:`"!"`                      
:func:`lbool(x)`   :math:`\lnot\lnot x`  :const:`"!!"`                       
:func:`land(x,y)`  :math:`x \land y`     :const:`"&&"`                       
:func:`lor(x,y)`   :math:`x \lor y`      :const:`"||"`                       
:func:`eq(x,y)`    :math:`x = y`         :const:`"=="`                   :func:`__eq(x,y)`
:func:`ne(x,y)`    :math:`x \neq y`      :const:`"!="` or :const:`"~="`  :func:`__eq(x,y)`
:func:`lt(x,y)`    :math:`x < y`         :const:`"<"`                    :func:`__lt(x,y)` [#]_
:func:`le(x,y)`    :math:`x <= y`        :const:`"<="`                   :func:`__le(x,y)`
:func:`gt(x,y)`    :math:`x > y`         :const:`">"`                    :func:`__le(x,y)`
:func:`ge(x,y)`    :math:`x >= y`        :const:`">="`                   :func:`__lt(x,y)`
=================  ====================  ==============================  =================

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
===============

A flag is 32 bit unsigned integer used to store up to 32 binary states with the convention that :const:`0` means disabled/cleared and :const:`1` means enabled/set. Functions on flags are useful aliases to -- or combinaison of -- bitwise operations to manipulate their states (i.e. their bits). Flags are mainly used by the object model to keep track of hidden and user-defined states in a compact and efficient format. 

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

Functors
========

Functors are objects that behave like functions with :type:`callable` semantic, and like readonly arrays with :type:`indexable` semantic translated into function call with the index as unique argument. The module :mod:`MAD.gfunc` offers few functions to expert users for creating and manipulating them.

.. function:: functor(f)

   Return a :type:`functor` that encapsulates the function (or any callable object) :var:`f`. Calling the returned functor is like calling :var:`f` itself with the same arguments. 

.. function:: compose(f, g)

   Return a :type:`functor` that encapsulates the composition of :var:`f` and :var:`g`. Calling the returned functor is like calling :math:`(f \circ g)(\dots)`. The operator :code:`f ^ g` is a shortcut for :func:`compose` if :var:`f` is a :type:`functor`.

.. function:: chain(f, g)

   Return a :type:`functor` that encapsulates the calls chain of :var:`f` and :var:`g`. Calling the returned functor is like calling :math:`f(\dots) ; g(\dots)`. The operator :code:`f .. g` is a shortcut for :func:`chain` if :var:`f` is a :type:`functor`.

.. function:: achain(f, g)

   Return a :type:`functor` that encapsulates the *ANDed* calls chain of :var:`f` and :var:`g`. Calling the returned functor is like calling :math:`f(\dots) \land g(\dots)`.

.. function:: ochain(f, g)

   Return a :type:`functor` that encapsulates the *ORed* calls chain of :var:`f` and :var:`g`. Calling the returned functor is like calling :math:`f(\dots) \lor g(\dots)`.

.. function:: bind1st(f, a)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`a` as its first argument. Calling the returned functor is like calling :math:`f(a,\dots)`.

.. function:: bind2nd(f, b)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`b` as its second argument. Calling the returned functor is like calling :math:`f(a,b,\dots)` where :var:`a` has to be provided.

.. function:: bind3rd(f, c)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`c` as its third argument. Calling the returned functor is like calling :math:`f(a,b,c,\dots)` where :var:`a` and :var:`b` have to be provided.

.. function:: bind2st(f, a, b)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`a` and :var:`b` as its two first arguments. Calling the returned functor is like calling :math:`f(a,b,\dots)`.

.. function:: bind3st(f, a, b, c)

   Return a :type:`functor` that encapsulates :var:`f` and binds :var:`a`, :var:`b` and :var:`c` as its three first arguments. Calling the returned functor is like calling :math:`f(a,b,c,\dots)`.

.. function:: bottom()

   Return a :type:`functor` that encapsulates the identity function :func:`ident` to define the *bottom* symbol of functors. Bottom is also available in the operator strings table :mod:`opstr` as :const:`"_|_"`.

.. function:: is_functor(a)

   Return :const:`true` if :var:`a` is a :type:`functor`, :const:`false` otherwise. This function is also available from the module :mod:`MAD.typeid`.


.. ---------------------------------------

.. rubric:: Footnotes

.. [#] For *true* Functional Programming, see the module :mod:`MAD.lfun`, a binding of the `LuaFun <https://github.com/luafun/luafun>`_  library adapted to the ecosystem of MAD-NG.

.. [#] Element-wise operators are only available for vector-like containers.

.. [#] Relational ordering operators are only available for ordered objects.