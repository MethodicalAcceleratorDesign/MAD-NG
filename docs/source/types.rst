.. index::
   Types

*****
Types
*****

This chapter describes some types identification and concepts setup defined by the module :mod:`MAD.typeid` and :mod:`MAD._C` (C API). The module :mod:`typeid` is extended by types from other modules on load like e.g. :type:`is_range`, :type:`is_complex`, :type:`is_matrix`, :type:`is_tpsa`, etc...   

Typeids
=======

All the functions listed hereafter return only :const:`true` or :const:`false` when identifying types.

Primitive Types
---------------

The following table shows the functions for identifying the primitive type of LuaJIT, i.e. using :expr:`type(a) == 'the_type'` 

========================  ====================================
Functions                 Return :const:`true` if :var:`a`
========================  ====================================
:func:`is_nil(a)`         is a :type:`nil`
:func:`is_boolean(a)`     is a :type:`boolean`
:func:`is_number(a)`      is a :type:`number`
:func:`is_string(a)`      is a :type:`string`
:func:`is_function(a)`    is a :type:`function`
:func:`is_table(a)`       is a :type:`table`
:func:`is_userdata(a)`    is a :type:`userdata`
:func:`is_coroutine(a)`   is a :type:`thread` [#f1]_
:func:`is_cdata(a)`       is a :type:`cdata`
========================  ====================================

Extended Types
--------------

The following table shows the functions for identifying the extended types, which are primitive types with some extensions, specializations or value ranges. 

=========================  ====================================
Functions                  Return :const:`true` if :var:`a`
=========================  ====================================
:func:`is_nan(a)`          is   :const:`nan` (Not a Number)
:func:`is_true(a)`         is   :const:`true`
:func:`is_false(a)`        is   :const:`false`
:func:`is_logical(a)`      is a :type:`boolean` or :const:`nil`
:func:`is_finite(a)`       is a :type:`number` with :math:`|a| < \infty`
:func:`is_infinite(a)`     is a :type:`number` with :math:`|a| = \infty`
:func:`is_positive(a)`     is a :type:`number` with :math:`a > 0`
:func:`is_negative(a)`     is a :type:`number` with :math:`a < 0`
:func:`is_zpositive(a)`    is a :type:`number` with :math:`a \ge 0`
:func:`is_znegative(a)`    is a :type:`number` with :math:`a \le 0`
:func:`is_nonzero(a)`      is a :type:`number` with :math:`a \ne 0`
:func:`is_integer(a)`      is a :type:`number` with :math:`-2^{52} \le a \le 2^{52}` and no fractional part
:func:`is_int32(a)`        is a :type:`number` with :math:`-2^{31} \le a < 2^{31}` and no fractional part
:func:`is_natural(a)`      is an :type:`integer` with :math:`a \ge 0`
:func:`is_even(a)`         is an even :type:`integer`
:func:`is_odd(a)`          is an odd :type:`integer`
:func:`is_decimal(a)`      is not an :type:`integer`
:func:`is_emptystring(a)`  is a :type:`string` with :expr:`#a == 0`
:func:`is_identifier(a)`   is a :type:`string` with valid identifier characters, i.e. :expr:`%s*[_%a][_%w]*%s*`
:func:`is_rawtable(a)`     is a :type:`table`  with no metatable
:func:`is_emptytable(a)`   is a :type:`table`  with no element
:func:`is_file(a)`         is a :type:`userdata` with :expr:`io.type(a) ~= nil`
:func:`is_openfile(a)`     is a :type:`userdata` with :expr:`io.type(a) == 'file'`
:func:`is_closedfile(a)`   is a :type:`userdata` with :expr:`io.type(a) == 'closed file'`
:func:`is_emptyfile(a)`    is an open :type:`file` with some content
=========================  ====================================

Concepts
========

Concepts are an extention of types looking at their behavior. The concepts are  more based on supported metamethods (or methods) than on the types themself and their valid range of values.

==========================  ====================================
Functions                   Return :const:`true` if :var:`a`
==========================  ====================================
:func:`is_value(a)`         is a :type:`nil`, a :type:`boolean`, a :type:`number` or a :type:`string`
:func:`is_reference(a)`     is not a :type:`value`
:func:`is_empty(a)`         is a :type:`mappable` and 1st iteration returns :const:`nil`
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
:func:`is_mutable(a)`       defines metamethod :expr:`__metatable()`
:func:`is_restricted(a)`    has metamethods for restriction, see :func:`wrestrict()`
:func:`is_protected(a)`     has metamethods for protection, see :func:`wprotect()`
:func:`is_deferred(a)`      has metamethods for deferred expressions, see :func:`deferred()`
:func:`is_same(a,b)`        has the same type and metatable as :var:`b`
==========================  ====================================

The functions in the following table are complementary to concepts and usually used to prevent an error during concepts checks.

===========================  ====================================
Functions                    Return :const:`true` if 
===========================  ====================================
:func:`has_member(a,b)`      :expr:`a[b]` is not :const:`nil`
:func:`has_method(a,f)`      :expr:`a[f]` is a :type:`callable`
:func:`has_metamethod(a,f)`  metamethod :var:`f` is defined
:func:`has_metatable(a)`     :var:`a` has a metatable
===========================  ====================================

.. function:: is_metaname(a)

   Returns :const:`true` if the :type:`string` :var:`a` is a valid metamethod name, :const:`false` otherwise.

.. function:: get_metatable(a)

   Returns the metatable of :var:`a` even if :var:`a` is a :type:`cdata`, which is not the case of :func:`getmetatable()`.

.. function:: get_metamethod(a,f)

   Returns the metamethod (or method) :var:`f` of :var:`a` even if :var:`a` is a :type:`cdata` and :var:`f` is only reachable through the metatable, or :const:`nil`.

Setting Concepts
----------------

.. data:: typeid.concept

   The :type:`table` :var:`concept` contains the lists of concepts that can be passed to the function :func:`set_concept` to prevent the use of their associated metamethods. The concepts can be combined together by adding them, e.g. :expr:`not_comparable = not_equalable + not_orderable`.

===========================  ====================================
Concepts                     Associated metamethods
===========================  ====================================
:const:`not_lengthable`      :func:`__len`
:const:`not_iterable`        :func:`__ipairs`
:const:`not_mappable`        :func:`__ipairs` and :func:`__pairs`
:const:`not_scannable`       :func:`__len`, :func:`__ipairs` and :func:`__pairs`
:const:`not_indexable`       :func:`__index`
:const:`not_extendable`      :func:`__newindex`
:const:`not_callable`        :func:`__call`
:const:`not_equalable`       :func:`__eq`
:const:`not_orderable`       :func:`__lt` and :func:`__le`
:const:`not_comparable`      :func:`__eq`, :func:`__lt` and :func:`__le`
:const:`not_concatenable`    :func:`__concat`
:const:`not_copiable`        :func:`__copy` and :func:`__same`
:const:`not_tablable`        :func:`__totable`
:const:`not_stringable`      :func:`__tostring`
:const:`not_mutable`         :func:`__metatable` and :func:`__newindex`
:const:`not_negatable`       :func:`__unm`
:const:`not_addable`         :func:`__add`
:const:`not_subtractable`    :func:`__sub`
:const:`not_additive`        :func:`__add` and :func:`__sub`
:const:`not_multipliable`    :func:`__mul`
:const:`not_dividable`       :func:`__div`
:const:`not_multiplicative`  :func:`__mul` and :func:`__div`
:const:`not_modulable`       :func:`__mod`
:const:`not_powerable`       :func:`__pow`
===========================  ====================================

.. function:: set_concept(mt, concepts, strict_)

   Return the metatable :var:`mt` after setting the metamethods associated to the combination of concepts set in :var:`concepts` to prevent their use. The concepts can be combined together by adding them, e.g. :expr:`not_comparable = not_equalable + not_orderable`. Metamethods can be overridden if :expr:`strict = false`, otherwise the overload is silently discarded. If :var:`concepts` requires :type:`iterable` but not :type:`mappable` then :func:`pairs` is equivalent to :func:`ipairs`.

.. function:: wrestrict(a)

   Return a proxy for :var:`a` which behaves like :var:`a`, except that it prevents existing indexes from being modified while allowing new ones to be created, i.e. :var:`a` is :type:`extendable`.

.. function:: wprotect(a)

   Return a proxy for :var:`a` which behaves like :var:`a`, except that it prevents existing indexes from being modified and does not allow new ones to be created, i.e. :var:`a` is :type:`readonly`.

.. function:: wunprotect(a)

   Return :var:`a` from the proxy, i.e. expect a restricted or a protected :var:`a`.

.. function:: deferred(a)

   Return a proxy for :var:`a` which behaves like :var:`a` except that elements of type :type:`function` will be considered as deferred expressions and evaluated on read, i.e. returning their results in their stead.

C Type Sizes
============

The following table lists the constants holding the size of the C types used by common :type:`cdata` like complex, matrices or TPSA. See section on `C API`_ for the description for those C types.

====================  ================  
C types sizes         C types            
====================  ================  
:const:`ctsz_log`     :c:type:`log_t`   
:const:`ctsz_idx`     :c:type:`idx_t`   
:const:`ctsz_ssz`     :c:type:`ssz_t`   
:const:`ctsz_dbl`     :c:type:`num_t`   
:const:`ctsz_cpx`     :c:type:`cpx_t`  
:const:`ctsz_str`     :c:type:`str_t`   
:const:`ctsz_ptr`     :c:type:`ptr_t`   
====================  ================  
       
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

.. c:type:: cpx_t

   The :type:`complex` type aliasing :type:`double _Complex`, i.e. two double precision 64-bit floating point numbers, that holds double-precision normalized number in IEC 60559.

.. c:type:: str_t

   The :type:`string` type aliasing :type:`const char*`, i.e. pointer to a readonly null-terminated array of characters.

.. c:type:: ptr_t

   The :type:`pointer` type aliasing :type:`const void*`, i.e. pointer to readonly memory of unknown/any type.

.. ------------------------------------------------------------

.. rubric:: Footnotes

.. [#f1] The Lua "threads" are user-level non-preemptive threads also named   coroutines.
