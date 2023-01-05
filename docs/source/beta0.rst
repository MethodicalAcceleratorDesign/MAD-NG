Beta0 Blocks
============
.. _ch.gen.beta0:



The ``beta0`` object is the *root object* of beta0 blocks that store information relative to the phase space at given positions, e.g. initial conditions, Poincaré section.

The ``beta0`` module extends the :doc:`typeid <types>` module with the ``is_beta0`` *function*, which returns ``true`` if its argument is a ``beta0`` object, ``false`` otherwise.

Attributes
----------

The ``beta0`` *object* provides the following attributes:

**particle**
	 A *string* specifying the name of the particle. (default: ``"positron"``).


Methods
-------

The ``beta0`` object provides the following methods:

**showdb**
	 A *method*	``([file])`` displaying the content of the particles database to ``file`` (default: ``io.stdout``).


Metamethods
-----------

The ``beta0`` object provides the following metamethods:

**__init**
	 A *metamethod*	``()`` returning \SLF after having processed the attributes with the special protect-and-update mechanism, where the order of the updates is undefined. It also creates new particle on-the-fly if the ``mass`` and the ``charge`` are defined, and then select it.




**__beta0**
	 A unique private *reference* that characterizes beta0 blocks.


Examples
--------

.. rubric:: Footnotes

