.. index::
   Types

*****
Types
*****

This chapter describes some types and concepts defined in the module :mod:`MAD.typeid` and :mod:`MAD._C` (C API). The module :mod:`typeid` is extended by other module on load, e.g. :type:`is_range`, :type:`is_complex`, :type:`is_matrix`, :type:`is_tpsa``, etc...   

Typeids
=======

All the functions for type identification return only :const:`true` or :const:`false`.

Primitive Types
---------------

The following table shows the functions for identifying the primitive type of LuaJIT, i.e. using :expr:`type(a) == 'type'` 

========================  ====================================
Functions                 Return :const:`true` if :var:`a`
========================  ====================================
:func:`is_nil(a)`  is     :const:`nil`
:func:`is_boolean(a)`     is a :const:`boolean`
:func:`is_number(a)`      is a :const:`number`
:func:`is_string(a)`      is a :const:`string`
:func:`is_function(a)`    is a :const:`function`
:func:`is_table(a)`       is a :const:`table`
:func:`is_userdata(a)`    is a :const:`userdata`
:func:`is_coroutine(a)`   is a :const:`thread`
:func:`is_cdata(a)`       is a :const:`cdata`
========================  ====================================

Extended Types
--------------

The following table shows the functions for identifying the extended types, that is primitive types with some extensions or specializations. 

=========================  ====================================
Functions                  Return :const:`true` if :var:`a`
=========================  ====================================
:func:`is_nan(a)`          is   :const:`nan` (Not a Number)
:func:`is_true(a)`         is   :const:`true`
:func:`is_false(a)`        is   :const:`false`
:func:`is_logical(a)`      is a :const:`boolean` or :const:`nil`
:func:`is_finite(a)`       is a :const:`number` with :math:`|a| < \infty`
:func:`is_infinite(a)`     is a :const:`number` with :math:`|a| = \infty`
:func:`is_positive(a)`     is a :const:`number` with :math:`a > 0`
:func:`is_negative(a)`     is a :const:`number` with :math:`a < 0`
:func:`is_zpositive(a)`    is a :const:`number` with :math:`a \ge 0`
:func:`is_znegative(a)`    is a :const:`number` with :math:`a \le 0`
:func:`is_nonzero(a)`      is a :const:`number` with :math:`a \ne 0`
:func:`is_integer(a)`      is a :const:`number` with :math:`-2^{52} \le a \le 2^{52}` and no fractional part
:func:`is_natural(a)`      is an :const:`integer` with :math:`a \ge 0`
:func:`is_even(a)`         is an even :const:`integer`
:func:`is_odd(a)`          is an odd :const:`integer`
:func:`is_decimal(a)`      is not an :const:`integer`
:func:`is_emptystring(a)`  is a :const:`string` with :expr:`#a == 0`
:func:`is_identifier(a)`   is a :const:`string` with valid identifier characters, i.e. :expr:`%s*[_%a][_%w]*%s*`
:func:`is_rawtable(a)`     is a :const:`table`  with no metatable
:func:`is_emptytable(a)`   is a :const:`table`  with no element
:func:`is_file(a)`         is a :const:`userdata` with :expr:`io.type(a) ~= nil`
:func:`is_openfile(a)`     is a :const:`userdata` with :expr:`io.type(a) == "file"`
:func:`is_closedfile(a)`   is a :const:`userdata` with :expr:`io.type(a) == "closed file"`
:func:`is_emptyfile(a)`    is an :const:`openfile` with some content
=========================  ====================================

Concepts
========

==========================  ====================================
Functions                   Return :const:`true` if :var:`a`
==========================  ====================================
:func:`is_value(a)`         is :const:`nil`, a :const:`boolean`, a :const:`number` or a :const:`string`
:func:`is_lengthable(a)`    supports operation :expr:`#a`
:func:`is_iterable(a)`      supports operation :expr:`ipairs(a)`
:func:`is_mappable(a)`      supports operation :expr:`pairs(a)`
:func:`is_indexable(a)`     supports operation :expr:`a[?]`
:func:`is_extendable(a)`    supports operation :expr:`a[]=?`
:func:`is_callable(a)`      supports operation :expr:`a()`
:func:`is_equalable(a)`     supports operation :expr:`a == ?`
:func:`is_orderable(a)`     supports operation :expr:`a < ?`
:func:`is_concatenable(a)`  supports operation :expr:`a .. ?`
:func:`is_negatable(a)`     supports operation :expr:`-a`
:func:`is_addable(a)`       supports operation :expr:`a + ?`
:func:`is_subtractable(a)`  supports operation :expr:`a - ?`
:func:`is_multipliable(a)`  supports operation :expr:`a * ?`
:func:`is_dividable(a)`     supports operation :expr:`a / ?`
:func:`is_modulable(a)`     supports operation :expr:`a % ?`
:func:`is_powerable(a)`     supports operation :expr:`a ^ ?`
:func:`is_copiable(a)`      supports metamethod :expr:`__copy()`
:func:`is_sameable(a)`      supports metamethod :expr:`__same()`
:func:`is_tablable(a)`      supports metamethod :expr:`__totable()`
:func:`is_stringable(a)`    supports metamethod :expr:`__tostring()`
:func:`is_mutable(a)`       supports metamethod :expr:`__metatable()`
:func:`is_empty(a)`         supports metamethod :expr:`__pairs()` and iterator returns :const:`nil`
==========================  ====================================

=========================  ====================================
Functions                  Return :const:`true` if 
=========================  ====================================
:func:`is_same(a,b)`       :var:`a` has the same type and metatable as :var:`b`
:func:`has_member(a,b)`    operation :expr:`a.b` does not return :const:`nil`
:func:`has_method(a,b)`    operation :expr:`a:b()` does not fail
=========================  ====================================

C API
=====

.. c:type:: log_t

   The :type:`logical` type aliasing :type:`_Bool`, i.e. boolean, that holds :const:`TRUE` or :const:`FALSE`.

.. c:type:: idx_t

   The :type:`index` type aliasing :type:`int32_t`, i.e. signed 32-bit integer, that holds signed indexes in the range :math:`[-2^{31}, 2^{31}-1]`.

.. c:type:: ssz_t

   The :type:`size` type aliasing :type:`int32_t`, i.e. signed 32-bit integer, that holds signed sizes in the range :math:`[-2^{31}, 2^{31}-1]`.

.. c:type:: num_t

   The :type:`number` type aliasing :type:`double`, i.e. double precision 64-bit floating point numbers, that holds double-precision normalized number in IEC 60559 in the approximative range :math:`\{-\infty\} \cup [-\text{huge}, -\text{tiny}] \cup \{0\} \cup [\text{tiny}, \text{huge}] \cup \{\infty\}` where :math:`\text{huge} \approx 10^{308}` and :math:`\text{tiny} \approx 10^{-308}`. See :const:`MAD.constant.huge` and :const:`MAD.constant.tiny` for precise values.

.. c:type:: cnum_t

   The :type:`complex` type aliasing :type:`double _Complex`, i.e. two double precision 64-bit floating point numbers, that holds double-precision normalized number in IEC 60559.

.. c:type:: str_t

   The :type:`string` type aliasing :type:`const char*`, i.e. pointer to a readonly null-terminated array of characters.

.. c:type:: ptr_t

   The :type:`pointer` type aliasing :type:`const void*`, i.e. pointer to readonly memory of unknown/any type.

.. ------------------------------------------------------------

.. rubric:: Footnotes
