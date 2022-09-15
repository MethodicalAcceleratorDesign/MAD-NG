.. index::
   Monomials

*********
Monomials
*********

This chapter describes *monomial* objects used by the differential algebra library to index the TPSA coefficients and available from the :mod:`MAD` environment. The module for monomials is not exposed, only the contructor is visible and thus, monomials must be handled directly by their methods.

Constructors
============

.. function:: monomial(len_, ord)

   Return a :type:`monomial` of length :var:`len` where the variable slots are set to the orders in the range of :const:`0..63` given by :var:`ord`. Default: :expr:`len = #ord`.
   
   If :var:`ord` is a :type:`number` then :var:`len` must be provided and all variable slots are set to the value of :var:`ord`.
   
   If :var:`ord` is a :type:`list` then :var:`len` can be omitted and all variable slots are set to the orders given by :var:`ord`.
   
   If :var:`ord` is a :type:`string` then :var:`len` can be omitted and all variable slots are set to the orders given by :var:`ord`, where each character in [0-9A-Za-z] is interpreted as an order in basis 62, e.g. the string :const:`"Bc"` represents a monomial of length 2 with orders 11 and 38. Note that orders > 61 are not specified when orders are encoded as strings.

Functions
=========

.. function:: is_monomial(a)

   Return :const:`true` if :var:`a` is a :type:`monomial`, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======


C API
=====