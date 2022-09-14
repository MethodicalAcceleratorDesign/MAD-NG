.. index::
   Numerical ranges

****************
Numerical Ranges
****************

This chapter describes *numerical ranges* directly available from the :mod:`MAD` environment. The module for numerical ranges is not exposed, only the contructors are visible and thus, numerical ranges must be handled directly by their methods. Numerical ranges are useful objects for representing numerical loops, intervals, discrete sets, (log)lines and linear spaces. Note that :type:`range` and :type:`logrange` have value semantic like :type:`number`.

Constructors
============

The constructors for :type:`range` and :type:`logrange` numbers are directly available from the :mod:`MAD` environment, except for the special case of the concatenation operator applied to two or three numbers, which is part of the language definition as a MAD-NG extension.

.. constant:: start..stop
              start..stop..step

   The concatenation operator applied to two or three numbers creates a :type:`range` and does not perform any adjustement. Default: :code:`step_ = 1`.

.. function:: range(start_, stop, step_)

   Return a :type:`range` object starting at :var:`start`, ending at :var:`stop` (included), with increments of size :var:`step`. The value of :var:`step` is adjusted if :expr:`|start + size * step - stop|` :math:`< 10^{-12}` for some :var:`size` to obtain zero. Default: :code:`start_ = 0, step_ = 1`.

.. function:: nrange(start_, stop, size_)

   Return a :type:`range` object starting at :var:`start`, ending at :var:`stop` (included), with :var:`size` increments. The value of :var:`step` is adjusted if :math:`|start + size * step - stop| < 10^{-12}` for some :var:`step` to obtain zero. Default: :code:`start_ = 0, size_ = 1`.

.. function:: logrange(start_, stop, step_)

   Return a :type:`logrange` object starting at :var:`start`, ending at :var:`stop` (included), with increments of size :var:`step`. The value of :var:`step` is adjusted if :math:`|start + size * step - stop| < 10^{-12}` for some :var:`size` to obtain zero. Default: :code:`start_ = 0, step_ = 1`.

.. function:: nlogrange(start_, stop, size_)

   Return a :type:`logrange` object starting at :var:`start`, ending at :var:`stop` (included), with :var:`size` increments. The value of :var:`step` is adjusted if :math:`|start + size * step - stop| < 10^{-12}` for some :var:`step` to obtain zero. Default: :code:`start_ = 0, size_ = 1`.

.. function:: torange(str)

   Return a :type:`range` decoded from the string :var:`str` containing the literal numerical ranges :const:`"a..b"` or :const:`"a..b..c"` where :var:`a`,  :var:`b` and :var:`c` are literal numbers.

Functions
=========

.. function:: is_range(a)
