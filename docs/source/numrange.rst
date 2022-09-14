.. index::
   Numerical ranges

****************
Numerical Ranges
****************

This chapter describes *numerical ranges* directly available from the :mod:`MAD` environment. The module for numerical ranges is not exposed, only the contructors are visible and thus, numerical ranges must be handled directly by their methods. Numerical ranges are useful objects for representing numerical loops, intervals, discrete sets, (log)lines and linear spaces. Note that :type:`range` and :type:`logrange` have value semantic like :type:`number`.

Constructors
============

The constructors for :type:`range` and :type:`logrange` are directly available from the :mod:`MAD` environment, except for the special case of the concatenation operator applied to two or three numbers, which is part of the language definition as a MAD-NG extension. The :type:`logrange` behave as a the :type:`range` but they work on logarithmic scale. All constructor functions adjust the value of :var:`step` if :expr:`|start + size * step - stop|` :math:`< 10^{-12}` for some :var:`size` to obtain zero and ensure stable sizes and iterators across platforms.

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

   Return a :type:`range` decoded from the string :var:`str` containing the literal numerical ranges :const:`"a..b"` or :const:`"a..b..c"` where :var:`a`,  :var:`b` and :var:`c` are literal numbers.

Functions
=========

.. function:: is_range(a)

   Return :const:`true` if :var:`a` is a :type:`range`, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

.. function:: is_logrange(a)

   Return :const:`true` if :var:`a` is a :type:`logrange` number, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

