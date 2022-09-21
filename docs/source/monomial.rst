.. index::
   Monomials

*********
Monomials
*********

This chapter describes `Monomial <https://en.wikipedia.org/wiki/Monomial>`_ objects useful to encode the variables powers of `Multivariate <https://en.wikipedia.org/wiki/Multivariable_calculus>`_ `Taylor Series <https://en.wikipedia.org/wiki/Taylor_series>`_ used by the `Differential Algebra <https://en.wikipedia.org/wiki/Differential_algebra>`_ library of MAD-NG. The module for monomials is not exposed, only the contructor is visible from the :mod:`MAD` environment and thus, monomials must be handled directly by their methods. Monomial objects do not know to which variables the stored orders belong, the relationship is only through the indexes. Note that monomials are objects with reference semantic that store variable orders as 8-bit unsigned integers, thus arithmetic on variable orders occurs in the ring :math:`\mathbb{N}/2^8\mathbb{N}`. 

Constructors
============

The constructor for :type:`monomial` is directly available from the :mod:`MAD` environment.

.. function:: monomial([len_,] ord_)

   Return a :type:`monomial` of size :var:`len` with the variable orders set to the values given by :var:`ord`, as computed by :func:`mono:fill(ord_)`. If :var:`ord` is omitted then :var:`len` must be provided. Default: :expr:`len_ = #ord`, :expr:`ord_ = 0`.

Attributes
==========

.. constant:: mono.n

   The number of variable orders in :var:`mono`, i.e. its size or length.

Functions
=========

.. function:: is_monomial(a)

   Return :const:`true` if :var:`a` is a :type:`monomial`, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

The optional argument :var:`r_` represents a destination placeholder for results.

.. function:: mono:same(n_)

   Return a monomial of length :var:`n` filled with zeros. Default: :expr:`n_ = #mono`.

.. function:: mono:copy(r_)

   Return a copy of :var:`mono`.

.. function:: mono:fill(ord_)

   Return :var:`mono` with the variable orders set to the values given by :var:`ord`. Default: :expr:`ord_ = 0`.

   - If :var:`ord` is a :type:`number` then all variable orders are set to the value of :var:`ord`.
   
   - If :var:`ord` is a :type:`list` then all variable orders are set to the values given by :var:`ord`.
   
   - If :var:`ord` is a :type:`string` then all variable orders are set to the values given by :var:`ord`, where each character in the set :const:`[0-9A-Za-z]` is interpreted as a variable order in the `Basis 62 <https://en.wikipedia.org/wiki/Base62>`_, e.g. the string :const:`"Bc"` will be interpreted as a monomial with variable orders 11 and 38. Characters not in the set :const:`[0-9A-Za-z]` are not allowed and lead to an undefined behavior, meaning that orders :math:`\ge 62` cannot be safely specified through strings.

.. function:: mono:min()

   Return the minimum variable order of :var:`mono`.

.. function:: mono:max()

   Return the maximum variable order of :var:`mono`.

.. function:: mono:ord()

   Return the order of :var:`mono`, that is the sum of all the variable orders.

.. function:: mono:ordp(step_)

   Return the product of the variable orders of :var:`mono` at every :var:`step`. Default: :expr:`step_ = 1`.

.. function:: mono:ordpf(step_)

   Return the product of the factorial of the variable orders of :var:`mono` at every :var:`step`. Default: :expr:`step_ = 1`.

.. function:: mono:add(mono2, r_)

   Return the sum of the monomials :var:`mono` and :var:`mono2`, that is the sum of the all their variable orders, i.e. :math:`(o_1 + o_2) \mod 256` where :math:`o_1` and :math:`o_2` are two variable orders at the same index in :var:`mono` and :var:`mono2`.

.. function:: mono:sub(mono2, r_)

   Return the difference of the monomials :var:`mono` and :var:`mono2`, that is the subtraction of the all their variable orders, i.e. :math:`(o_1 - o_2) \mod 256` where :math:`o_1` and :math:`o_2` are two variable orders at the same index in :var:`mono` and :var:`mono2`.

.. function:: mono:concat(mono2, r_)

   Return the concatenation of the monomials :var:`mono` and :var:`mono2`.

.. function:: mono:reverse(r_)

   Return the reverse of the monomial :var:`mono`.

.. function:: mono:totable()

   Return a :type:`list` containing all the variable orders of :var:`mono`.

.. function:: mono:tostring(sep_)

   Return a :type:`string` containing all the variable orders of :var:`mono` encoded with characters in the set :const:`[0-9A-Za-z]` and separated by the :type:`string` :var:`sep`. Default: :expr:`sep_ = ''`.

Operators
=========

.. function:: #mono

   Return the number of variable orders in :var:`mono`, i.e. its length.

.. function:: mono[n]

   Return the variable order at index :var:`n` for :expr:`1 <= n <= #mono`, :const:`nil` otherwise.

.. function:: mono[n] = v

   Assign the value :var:`v` to the variable order at index :var:`n` for :expr:`1 <= n <= #mono`, otherwise raise an *"out of bounds"* error.

.. function:: mono + mono2

   Equivalent to :expr:`mono:add(mono2)`.

.. function:: mono - mono2

   Equivalent to :expr:`mono:sub(mono2)`.

.. function:: mono < mono2

   Return :const:`false` if one variable order in :var:`mono` is greater or equal to the variable order at the same index in :var:`mono2`, :const:`true` otherwise.

.. function:: mono <= mono2

   Return :const:`false` if one variable order in :var:`mono` is greater than the variable order at the same index in :var:`mono2`, :const:`true` otherwise.

.. function:: mono == mono2

   Return :const:`false` if one variable order in :var:`mono` is not equal to the variable order at the same index in :var:`mono2`, :const:`true` otherwise.

.. function:: mono .. mono2

   Equivalent to :expr:`mono:concat(mono2)`.

Iterators
=========

.. function:: ipairs(mono)
   :noindex:

   Return an :type:`ipairs` iterator suitable for generic :const:`for` loops. The generated values are those returned by :func:`mono[i]`. 

C API
=====

.. c:type:: ord_t

   The variable order type, which is an alias for 8-bit unsigned integer. In the C API, monomials are arrays of variable orders with their size :var:`n` tracked separately, i.e. :var:`a[n]`. 

.. c:function:: ssz_t mad_mono_str (ssz_t n, ord_t a[n], str_t s)

   Return the number of converted characters from the :type:`string` :var:`s` into variable orders stored to the monomial :var:`a[n]`, as decribed in the method :func:`:fill()`.

.. c:function:: str_t mad_mono_prt (ssz_t n, const ord_t a[n], char s[n+1])

   Return the :type:`string` :var:`s` filled with characters resulting from the conversion of the variable orders given in the monomial :var:`a[n]`, as decribed in the method :func:`:tostring()`.

.. c:function:: void mad_mono_fill (ssz_t n, ord_t a[n], ord_t v)

   Fill the monomial :var:`a[n]` with the variable order :var:`v`.

.. c:function:: void mad_mono_copy (ssz_t n, const ord_t a[n], ord_t r[n])

   Copy the monomial :var:`a[n]` to the monomial :var:`r[n]`.

.. c:function:: ord_t mad_mono_min  (ssz_t n, const ord_t a[n])

   Return the minimum variable order of the monomial :var:`a[n]`.

.. c:function:: ord_t mad_mono_max (ssz_t n, const ord_t a[n])

   Return the minimum variable order of the monomial :var:`a[n]`.

.. c:function:: int mad_mono_ord (ssz_t n, const ord_t a[n])

   Return the order of the monomial :var:`a[n]`.

.. c:function:: num_t mad_mono_ordp (ssz_t n, const ord_t a[n], idx_t stp)

   Return the product of the variable orders of the monomial :var:`a[n]` at every :var:`stp`.

.. c:function:: num_t mad_mono_ordpf (ssz_t n, const ord_t a[n], idx_t stp)

   Return the product of the factorial of the variable orders of the monomial :var:`a[n]` at every :var:`stp`.

.. c:function:: log_t mad_mono_eq (ssz_t n, const ord_t a[n], const ord_t b[n])

   Return :const:`FALSE` if one variable order in monomial :var:`a[n]` is not equal to the variable order at the same index in monomial :var:`b[n]`, :const:`TRUE` otherwise.

.. c:function:: log_t mad_mono_lt (ssz_t n, const ord_t a[n], const ord_t b[n])

   Return :const:`FALSE` if one variable order in monomial :var:`a[n]` is greater or equal to the variable order at the same index in monomial :var:`b[n]`, :const:`TRUE` otherwise.

.. c:function:: log_t mad_mono_le (ssz_t n, const ord_t a[n], const ord_t b[n])

   Return :const:`FALSE` if one variable order in monomial :var:`a[n]` is greater than the variable order at the same index in monomial :var:`b[n]`, :const:`TRUE` otherwise.

.. c:function:: int mad_mono_cmp (ssz_t n, const ord_t a[n], const ord_t b[n])

   Return the difference between the first variable orders that are not equal for a given index starting from the beginning in monomials :var:`a[n]` and :var:`b[n]`.

.. c:function:: int mad_mono_rcmp (ssz_t n, const ord_t a[n], const ord_t b[n])

   Return the difference between the first variable orders that are not equal for a given index starting from the end in monomials :var:`a[n]` and :var:`b[n]`.

.. c:function:: void mad_mono_add (ssz_t n, const ord_t a[n], const ord_t b[n], ord_t r[n])

   Put the sum of the monomials :var:`a[n]` and :var:`b[n]` in the monomial :var:`r[n]`.

.. c:function:: void mad_mono_sub (ssz_t n, const ord_t a[n], const ord_t b[n], ord_t r[n])

   Put the difference of the monomials :var:`a[n]` and :var:`b[n]` in the monomial :var:`r[n]`.

.. c:function:: void mad_mono_cat (ssz_t n, const ord_t a[n], ssz_t m, const ord_t b[m], ord_t r[n+m])

   Put the concatenation of the monomials :var:`a[n]` and :var:`b[m]` in the monomial :var:`r[n+m]`.

.. c:function:: void mad_mono_rev (ssz_t n, const ord_t a[n], ord_t r[n])

   Put the reverse of the monomial :var:`a[n]` in the monomial :var:`r[n]`.

.. c:function:: void mad_mono_print (ssz_t n, const ord_t a[n], FILE *fp_)

   Print the monomial :var:`a[n]` to the file :var:`fp`. Default: :expr:`fp_ = stdout`.
