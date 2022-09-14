.. index::
   Numerical ranges

****************
Numerical Ranges
****************

This chapter describes *numerical ranges* directly available from the :mod:`MAD` environment. The module for numerical ranges is not exposed, only the contructors are visible and thus, numerical ranges must be handled directly by their methods. Numerical ranges are useful objects for representing numerical loops, intervals, discrete sets, (log)lines and linear spaces. Note that :type:`range` and :type:`logrange` have value semantic like :type:`number`.

Constructors
============

The constructors for :type:`range` and :type:`logrange` are directly available from the :mod:`MAD` environment, except for the special case of the concatenation operator applied to two or three numbers, which is part of the language definition as a MAD-NG extension. The :type:`logrange` behave as a the :type:`range` but they work on logarithmic scale. All constructor functions adjust the value of :var:`step` to ensure stable sizes and iterators across platforms (see the method :func:`adjust` for details).

.. constant:: start..stop
              start..stop..step

   The concatenation operator applied to two or three numbers creates a :type:`range` and does not perform any adjustment of :var:`step`. Default: :code:`step_ = 1`.

.. function:: range(start_, stop, step_)

   Return a :type:`range` object starting at :var:`start`, ending at :var:`stop` (included), with increments of size :var:`step`. Default: :code:`start_ = 0, step_ = 1`.

.. function:: nrange(start_, stop, size_)

   Return a :type:`range` object starting at :var:`start`, ending at :var:`stop` (included), with :var:`size` increments. Default: :code:`start_ = 0, size_ = 1`.

.. function:: logrange(start_, stop, step_)

   Return a :type:`logrange` object starting at :var:`start`, ending at :var:`stop` (included), with increments of size :var:`step`. Default: :code:`start_ = 0, step_ = 1`.

.. function:: nlogrange(start_, stop, size_)

   Return a :type:`logrange` object starting at :var:`start`, ending at :var:`stop` (included), with :var:`size` increments. Default: :code:`start_ = 0, size_ = 1`.

.. function:: torange(str)

   Return a :type:`range` decoded from the string :var:`str` containing a literal numerical ranges of the form :const:`"a..b"` or :const:`"a..b..c"` where :var:`a`,  :var:`b` and :var:`c` are literal numbers.

Functions
=========

.. function:: is_range(a)

   Return :const:`true` if :var:`a` is a :type:`range`, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

.. function:: is_logrange(a)

   Return :const:`true` if :var:`a` is a :type:`logrange` number, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

Unless specified, the object :var:`rng` owning the methods stands for a :type:`range` or a :type:`logrange` indifferently.

.. function:: rng:is_empty()

   Return :const:`false` if :var:`rng` contains at least one value, :const:`true` otherwise.

.. function:: rng:same()

   Return :var:`rng` itself. This method is the identity for objects with value semantic.

.. function:: rng:copy()

   Return :var:`rng` itself. This method is the identity for objects with value semantic.

.. function:: rng:size()

   Return the number of steps contained by the range.

.. function:: rng:step()

   Return the :var:`step` component of the range, which may slighlty differ from the value provided to the constructors due to adjustment. 

.. function:: rng:value(x)

   Return the interpolated value at :var:`x`, i.e. interpreting the range as a (log)line with equation :expr:`start + x * step` 

.. function:: rng:get(x)
   
   Return :func:`rng:value(x)` if the results is inside the range's bounds, :const:`nil` otherwise. 

.. function:: rng:last()

   Return the last value inside the range's bounds, :const:`nil` otherwise. 

.. function:: rng:adjust()

   Return a range with a :var:`step` adjusted.

   The internal quantity :var:`step` is adjusted if the computed size is close to an integer by :math:`±10^{-12}`. Then the following properties should hold even for rational numbers (in binary representation) given a consistent input for :var:`start`, :var:`stop`, :var:`step` and :var:`size`:

   - :expr:`range (start, stop, step):size()        == size`
   - :expr:`nrange(start, stop, size):step()        == step`
   - :expr:`range (start, stop, step):value(size-1) == stop`
   
   The maximum adjustment is :expr:`step = step * (1-eps)^2`, beyond this value it is the user reponsibility to provide better inputs.

.. function:: rng:ranges()

   Return the three numbers characterising the range :var:`rng`, namely its :var:`start`, :var:`stop` and :var:`step` in this order. 

.. function:: rng:bounds()

   Return the three numbers characterising the boundaries of the range :var:`rng`, namely its :var:`start`, :var:`last` and :var:`step` :math:`>0` in this order, :const:`nil` otherwise.

.. function:: rng:overlap(rng2)

   Return :const:`true` if :var:`rng` and :var:`rng2` overlap, i.e. have intersecting bounds, :const:`false` otherwise.
   
.. function:: rng:reverse()

   Return a range which is the reverse of the range :var:`rng`, i.e. swap :var:`start` and :var:`stop`, and reverse :var:`step`.

.. function:: rng:log()

   Return a :type:`logrange` build from the conversion of the :type:`range` :var:`rng`.

.. function:: rng:tostring()

   Return a :type:`string` encoding the range :var:`rng` into a literal numerical ranges of the form :const:`"a..b"` or :const:`"a..b..c"` where :var:`a`,  :var:`b` and :var:`c` are literal numbers.

.. function:: rng:totable()

   Return a :type:`table` filled with :func:`rng:size()` values computed by :func:`rng:value()`. Note that ranges are objects with a very small memory footprint while the generated tables can be huge.

Operators
=========

.. function:: #rng

   Return the size of the range as computed by :func:`rng:size()`.

.. function:: rng[n]

   Return the value interpolated by the range as computed by :func:`rng:get(n-1)`, i.e. assuming an index-like interpolation.

.. function:: -rng

   Return a range with all components :var:`start`, :var:`stop` and :var:`step` reversed.

.. function:: rng + num
              num + rng

   Return a range with :var:`start` and :var:`stop` shifted by :var:`num`.

.. function:: rng - num

   Return a range with :var:`start` and :var:`stop` shifted by :var:`-num`, i.e. it is equivalent to :expr:`rng + (-num)`.

.. function:: num - rng

   Return a range reversed with :var:`start` and :var:`stop` shifted by :var:`num`, i.e. it is equivalent to :expr:`num + (-rng)`.

.. function:: num * rng
              rng * num

   Return a range with :var:`start`, :var:`stop` and :var:`step` scaled by :var:`num`.

.. function:: rng / num

   Return a range with :var:`start`, :var:`stop` and :var:`step` scaled by :var:`1/num`, i.e. it is equivalent to :expr:`rng * (1/num)`.

.. function:: rng == rng2

   Return :const:`true` if :var:`rng` and :var:`rng2` are of same king, have equal :var:`start` and :var:`stop`, and their :var:`step` are within one :const:`eps` from each other, :const:`false` otherwise.

Iterators
=========

 .. function:: ipairs(rng)

   Return an :type:`ipairs` iterator suitable for generic :const:`for` loops. The generated values are those returned by :func:`rng:value(i)`. 

