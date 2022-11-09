.. index::
   Differential maps
   Taylor Series
   Truncated Power Series Algebra, TPSA

*****************
Differential Maps
*****************

This chapter describes real :type:`damap` and complex :type:`cdamap` objects as supported by MAD-NG. They are useful abstractions to represent non-linear parametric multivariate differential maps, i.e. `Diffeomorphisms <https://en.wikipedia.org/wiki/Diffeomorphism>`_, `Vector Fields <https://en.wikipedia.org/wiki/Vector_field>`_, `Exponential Maps <https://en.wikipedia.org/wiki/Exponential_map_(Lie_theory)>`_ and `Lie Derivative <https://en.wikipedia.org/wiki/Lie_algebra>`_.  The module for the differential maps is not exposed, only the contructors are visible from the :mod:`MAD` environment and thus, differential maps are handled directly by their methods or by the generic functions of the same name from the module :mod:`MAD.gmath`. Note that :type:`damap` and :type:`cdamap` are defined as C structure for direct compliance with the C API.

Introduction
============


Constructors
============

Functions
=========

Methods
=======

Operators
=========

Iterators
=========

C API
=====

