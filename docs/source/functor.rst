.. index::
   Functors

********
Functors
********

This chapter describes how to create, combine and use *functors* from the :mod:`MAD` environment. Functors are objects that behave like functions with :type:`callable` semantic, and also like readonly arrays with :type:`indexable` semantic, where the index is translated as a unique argument into the function call. They are mainly used by the object model to distinguish them from functions which are interpreted as deferred expressions and evaluated automatically on reading, and by the Survey and Track tracking codes to handle (user-defined) actions. 

Constructors
============

This module provides mostly constructors to create functors from functions, functors and any objects with :type:`callable` semantic, and combine them all together.

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

Functions
=========

.. function:: is_functor(a)

   Return :const:`true` if :var:`a` is a :type:`functor`, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.
