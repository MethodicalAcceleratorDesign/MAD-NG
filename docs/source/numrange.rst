.. index::
   Numerical ranges

****************
Numerical Ranges
****************

This chapter describes :type:`range` and :type:`logrange` objects that are useful abstaction of numerical loops, intervals, discrete sets, (log)lines and linear spaces. The module for numerical ranges is not exposed, only the contructors are visible from the :mod:`MAD` environment and thus, numerical ranges must be handled directly by their methods. Note that :type:`range` and :type:`logrange` have value semantic like :type:`number`.

Constructors
============

The constructors for :type:`range` and :type:`logrange` are directly available from the :mod:`MAD` environment, except for the special case of the concatenation operator applied to two or three :type:`number`, which is part of the language definition as a MAD-NG extension. The :type:`logrange` behave as a the :type:`range` but they work on logarithmic scale. All constructor functions adjust the value of :var:`step` to ensure stable sizes and iterators across platforms (see the method :func:`adjust` for details).

.. constant:: start..stop
              start..stop..step

   The concatenation operator applied to two or three numbers creates a :type:`range` and does not perform any adjustment of :var:`step`. The default step for the first form is one.

.. function:: range([start_,] stop, step_)

   Return a :type:`range` object starting at :var:`start`, ending at :var:`stop` (included), with increments of size :var:`step`. Default: :expr:`start_ = 1, step_ = 1`.

.. function:: nrange([start_,] stop, size_)

   Return a :type:`range` object starting at :var:`start`, ending at :var:`stop` (included), with :var:`size` increments. Default: :expr:`start_ = 1, size_ = 100`.

.. function:: logrange([start_,] stop, step_)

   Return a :type:`logrange` object starting at :var:`start`, ending at :var:`stop` (included), with increments of size :var:`step`. Default: :expr:`start_ = 1, step_ = 1`.

.. function:: nlogrange([start_,] stop, size_)

   Return a :type:`logrange` object starting at :var:`start`, ending at :var:`stop` (included), with :var:`size` increments. Default: :expr:`start_ = 1, size_ = 100`.

.. function:: torange(str)

   Return a :type:`range` decoded from the string :var:`str` containing a literal numerical ranges of the form :const:`"a..b"` or :const:`"a..b..c"` where :var:`a`,  :var:`b` and :var:`c` are literal numbers.

Empty Ranges
^^^^^^^^^^^^

   Empty ranges of size zero can be created by fulfilling the constraints :expr:`start > stop` and :expr:`step > 0` or :expr:`start < stop` and :expr:`step < 0` in :type:`range` constructor.

Singleton Ranges
^^^^^^^^^^^^^^^^

   Singleton ranges of size one can be created by fulfilling the constraints :expr:`step > stop-start` for :expr:`start < stop` and :expr:`step < stop-start` for :expr:`stop < start` in :type:`range` constructor or :expr:`size == 1` in :type:`nrange` constructor. In this latter case, :var:`step` will be set to :expr:`step = huge * sign(stop-start)`.

Constant Ranges
^^^^^^^^^^^^^^^

   Constant ranges of infinite size can be created by fulfilling the constraints :expr:`start == stop` and :expr:`step == 0` in :type:`range` constructor or :expr:`size == inf` in :type:`nrange` constructor. The user must satify the constraint :expr:`start == stop` in both constructors to show its intention.

Attributes
==========

.. constant:: rng.start
              rng.logstart

   The component *start* of the :type:`range` and the :type:`logrange` on a linear scale. 

.. constant:: rng.stop
              rng.logstop

   The component *stop* of the :type:`range` and the :type:`logrange` on a linear scale. 

.. constant:: rng.step
              rng.logstep

   The component *step* of the :type:`range` and the :type:`logrange` on a linear scale, which may slighlty differ from the value provided to the constructors due to adjustment. 

Functions
=========

.. function:: is_range(a)
              is_logrange(a)

   Return :const:`true` if :var:`a` is respectively a :type:`range` or a :type:`logrange`, :const:`false` otherwise. These functions are only available from the module :mod:`MAD.typeid`.

.. function:: isa_range(a)

   Return :const:`true` if :var:`a` is a :type:`range` or a :type:`logrange` (i.e. is-a range), :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

Unless specified, the object :var:`rng` that owns the methods represents either a :type:`range` or a :type:`logrange`.

.. function:: rng:is_empty()

   Return :const:`false` if :var:`rng` contains at least one value, :const:`true` otherwise.

.. function:: rng:same()

   Return :var:`rng` itself. This method is the identity for objects with value semantic.

.. function:: rng:copy()

   Return :var:`rng` itself. This method is the identity for objects with value semantic.

.. function:: rng:ranges()

   Return the values of :var:`start`, :var:`stop` and :var:`step`, fully characterising the range :var:`rng`.

.. function:: rng:size()

   Return the number of values contained by the range :var:`rng`, i.e. its size that is the number of steps plus one.

.. function:: rng:value(x)

   Return the interpolated value at :var:`x`, i.e. interpreting the range  :var:`rng` as a (log)line with equation :expr:`start + x * step`.

.. function:: rng:get(x)
   
   Return :func:`rng:value(x)` if the result is inside the range's bounds, :const:`nil` otherwise. 

.. function:: rng:last()

   Return the last value inside the bounds of the range :var:`rng`, :const:`nil` otherwise. 

.. function:: rng:adjust()

   Return a range with a :var:`step` adjusted.

   The internal quantity :var:`step` is adjusted if the computed size is close to an integer by :math:`\pm10^{-12}`. Then the following properties should hold even for rational binary numbers given a consistent input for :var:`start`, :var:`stop`, :var:`step` and :var:`size`:

   - :expr:`#range(start, stop, step)               == size`
   - :expr:`nrange(start, stop, size):step()        == step`
   - :expr:`range (start, stop, step):value(size-1) == stop`
   
   The maximum adjustment is :expr:`step = step * (1-eps)^2`, beyond this value it is the user reponsibility to provide better inputs.

.. function:: rng:bounds()

   Return the values of :var:`start`, :var:`last` (as computed by :func:`rng:last()`) and :var:`step` (made positive) characterising the boundaries of the range :var:`rng`, i.e. interpreted as an interval, :const:`nil` otherwise.

.. function:: rng:overlap(rng2)

   Return :const:`true` if :var:`rng` and :var:`rng2` overlap, i.e. have intersecting bounds, :const:`false` otherwise.
   
.. function:: rng:reverse()

   Return a range which is the reverse of the range :var:`rng`, i.e. swap :var:`start` and :var:`stop`, and reverse :var:`step`.

.. function:: rng:log()

   Return a :type:`logrange` built by converting the :type:`range` :var:`rng` to logarithmic scale.

.. function:: rng:unm()

   Return a range with all components :var:`start`, :var:`stop` and :var:`step` negated.

.. function:: rng:add(num)

   Return a range with :var:`start` and :var:`stop` shifted by :expr:`num`.

.. function:: rng:sub(num)

   Return a range with :var:`start` and :var:`stop` shifted by :expr:`-num`.

.. function:: rng:mul(num)

   Return a range with :var:`stop` and :var:`step` scaled by :expr:`num`.

.. function:: rng:div(num)

   Return a range with :var:`stop` and :var:`step` scaled by :expr:`1/num`.

.. function:: rng:tostring()

   Return a :type:`string` encoding the range :var:`rng` into a literal numerical ranges of the form :const:`"a..b"` or :const:`"a..b..c"` where :var:`a`,  :var:`b` and :var:`c` are literal numbers.

.. function:: rng:totable()

   Return a :type:`table` filled with :expr:`#rng` values computed by :func:`rng:value()`. Note that ranges are objects with a very small memory footprint while the generated tables can be huge.

Operators
=========

.. function:: #rng

   Return the number of values contained by the range :var:`rng`, i.e. it is equivalent to :expr:`rng:size()`.

.. function:: rng[n]

   Return the value at index :var:`n` contained by the range :var:`rng`, i.e. it is equivalent to :expr:`rng:get(round(n-1))`.

.. function:: -rng

   Equivalent to :expr:`rng:unm()`. 

.. function:: rng + num
              num + rng

   Equivalent to :expr:`rng:add(num)`. 

.. function:: rng - num

   Equivalent to :expr:`rng:sub(num)`.

.. function:: num - rng

   Equivalent to :expr:`rng:unm():add(num)`.

.. function:: num * rng
              rng * num

   Equivalent to :expr:`rng:mul(num)`.

.. function:: rng / num

   Equivalent to :expr:`rng:div(num)`.

.. function:: rng == rng2

   Return :const:`true` if :var:`rng` and :var:`rng2` are of same king, have equal :var:`start` and :var:`stop`, and their :var:`step` are within one :const:`eps` from each other, :const:`false` otherwise.

Iterators
=========

.. function:: ipairs(rng)
   :noindex:

   Return an :type:`ipairs` iterator suitable for generic :const:`for` loops. The generated values are those returned by :func:`rng:value(i)`. 
