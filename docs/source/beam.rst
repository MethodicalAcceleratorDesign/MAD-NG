Beams
=====
.. _ch.gen.beam:

The ``beam`` object is the *root object* of beams that store information relative to particles and particle beams. It also provides a simple interface to the particles and nuclei database.

The ``beam`` module extends the :doc:`typeid <types>` module with the ``is_beam`` *function*, which returns ``true`` if its argument is a ``beam`` object, ``false`` otherwise.

Attributes
----------

The ``beam`` *object* provides the following attributes:

**particle**
	 A *string* specifying the name of the particle. (default: ``"positron"``).

**mass**
	 A *number* specifying the energy-mass of the particle [GeV]. (default: ``emass``).

**charge**
	 A *number* specifying the charge of the particle in [q] unit of ``qelect``. [#f1]_ (default: ``1``).

**spin**
	 A *number* specifying the spin of the particle. (default: ``0``).

**emrad**
	 A *lambda* returning the electromagnetic radius of the particle [m], 

	 :math:`\mathrm{emrad} = \mathrm{krad\_GeV}\times\mathrm{charge}^2/\mathrm{mass}` where :math:`\mathrm{krad\_GeV} = 10^{-9} (4 \pi\varepsilon_0)^{-1} q`.

**aphot**
	 A *lambda* returning the average number of photon emitted per bending unit, 

	 :math:`\mathrm{aphot} = \mathrm{kpht\_GeV}\times\mathrm{charge}^2*\mathrm{betgam}` where :math:`\mathrm{kpht\_GeV}` :math:`= \frac{5}{2\sqrt{3}}`.

**energy**
	 A *number* specifying the particle energy [GeV]. (default: ``1``).

**pc**
	 A *lambda* returning the particle momentum times the speed of light [GeV],

	 :math:`\mathrm{pc} = (\mathrm{energy}^2 - \mathrm{mass}^2)^{\frac{1}{2}}`.

**beta**
	 A *lambda* returning the particle relativistic :math:`\beta=\frac{v}{c}`,

	 :math:`\mathrm{beta} = (1 - (\mathrm{mass}/\mathrm{energy})^2)^{\frac{1}{2}}`.

**gamma**
	 A *lambda* returning the particle Lorentz factor :math:`\gamma=(1-\beta^2)^{-\frac{1}{2}}`,

	 :math:`\mathrm{gamma} = \mathrm{energy}/\mathrm{mass}`.

**betgam**
	 A *lambda* returning the product :math:`\beta\gamma`,

	 :math:`\mathrm{betgam} = (\mathrm{gamma}^2 - 1)^\frac{1}{2}`.

**pc2**
	 A *lambda* returning :math:`\mathrm{pc}^2`, avoiding the square root.

**beta2**
	 A *lambda* returning :math:`\mathrm{beta}^2`, avoiding the square root.

**betgam2**
	 A *lambda* returning :math:`\mathrm{betgam}^2`, avoiding the square root.

**brho**
	 A *lambda* returning the magnetic rigidity [T.m], 

	 :literal:`brho = GeV_c * pc/|charge|` where ``GeV_c`` = :math:`10^{9}/c`

**ex**
	 A *number* specifying the horizontal emittance :math:`\epsilon_x` [m]. (default: ``1``).

**ey**
	 A *number* specifying the vertical emittance :math:`\epsilon_y` [m]. (default: ``1``).

**et**
	 A *number* specifying the longitudinal emittance :math:`\epsilon_t` [m]. (default: :literal:`1e-3`).

**exn**
	 A *lambda* returning the normalized horizontal emittance [m], 

	 :literal:`exn = ex * betgam`.

**eyn**
	 A *lambda* returning the normalized vertical emittance [m], 

	 :literal:`eyn = ey * betgam`.

**etn**
	 A *lambda* returning the normalized longitudinal emittance [m], 

	 :literal:`etn = et * betgam`.

**nbunch**
	 A *number* specifying the number of particle bunches in the machine. (default: ``0``).

**npart**
	 A *number* specifying the number of particles per bunch. (default: ``0``).

**sigt**
	 A *number* specifying the bunch length in :math:`c \sigma_t`. (default: ``1``).

**sige**
	 A *number* specifying the relative energy spread in :math:`\sigma_E/E` [GeV]. (default: :literal:`1e-3`).


The ``beam`` *object* also implements a special protect-and-update mechanism for its attributes to ensure consistency and precedence between the physical quantities stored internally:

#.	 The following attributes are *read-only*, i.e. writing to them triggers an error:
		``mass, charge, spin, emrad, aphot``

#.	 The following attributes are *read-write*, i.e. hold values, with their accepted numerical ranges:
		``particle, energy`` :math:`>` ``mass``,
		``ex`` :math:`>0`, ``ey`` :math:`>0`, ``et`` :math:`>0`,
		``nbunch`` :math:`>0`, ``npart`` :math:`>0`, ``sigt`` :math:`>0`, ``sige`` :math:`>0`.

#.	 The following attributes are *read-update*, i.e. setting these attributes update the ``energy``, with their accepted numerical ranges:
		``pc`` :math:`>0`, :math:`0.9>` ``beta`` :math:`>0`, ``gamma`` :math:`>1`, ``betgam`` :math:`>0.1`, ``brho`` :math:`>0`,
		``pc2``, ``beta2``, ``betgam2``.
#.	 The following attributes are *read-update*, i.e. setting these attributes update the emittances ``ex``, ``ey``, and ``et`` repectively, with their accepted numerical ranges:
		``exn`` :math:`>0`, ``eyn`` :math:`>0`, ``etn`` :math:`>0`.


Methods
-------

The ``beam`` object provides the following methods:

**new_particle**
	 A *method*	``(particle, mass, charge, [spin])`` creating new particles or nuclei and store them in the particles database. The arguments specify in order the new ``particle``'s name, energy-``mass`` [GeV], ``charge`` [q], and ``spin`` (default: ``0``). These arguments can also be grouped into a *table* with same attribute names as the argument names and passed as the solely argument.

**set_variables**
	 A *method*	``(set)`` returning ``self`` with the attributes set to the pairs (*key*, *value*) contained in ``set``. This method overrides the original one to implement the special protect-and-update mechanism, but the order of the updates is undefined. It also creates new particle on-the-fly if the ``mass`` and the ``charge`` are defined, and then select it. Shortcut ``setvar``.

**showdb**
	 A *method*	``([file])`` displaying the content of the particles database to ``file`` (default: ``io.stdout``).


Metamethods
-----------

The ``beam`` object provides the following metamethods:

**__init**
	 A *metamethod*	``()`` returning ``self`` after having processed the attributes with the special protect-and-update mechanism, where the order of the updates is undefined. It also creates new particle on-the-fly if the ``mass`` and the ``charge`` are defined, and then select it.

**__newindex**
	 A *metamethod*	``(key, val)`` called by the assignment operator ``[key]=val`` to create new attributes for the pairs (*key*, *value*) or to update the underlying physical quantity of the ``beam`` objects.


The following attribute is stored with metamethods in the metatable, but has different purpose:


**__beam**
	 A unique private *reference* that characterizes beams.


Particles database
------------------

The ``beam`` *object* manages the particles database, which is shared by all ``beam`` instances. The default set of supported particles is:
		electron, positron, proton, antiproton, neutron, antineutron, ion, muon, 
		antimuon, deuteron, antideuteron, negmuon (=muon), posmuon (=antimuon).

New particles can be added to the database, either explicitly using the ``new_particle`` method, or by creating or updating a beam *object* and specifying all the attributes of a particle, i.e. ``particle``'s name, ``charge``, ``mass``, and (optional) ``spin``:

.. code-block:: lua
	
	local beam in MAD
	local nmass, pmass, mumass in MAD.constant
	
	-- create a new particle
	beam:new_particle{ particle='mymuon', mass=mumass, charge=-1, spin=1/2 }
	
	-- create a new beam and a new nucleus
	local pbbeam = beam { particle='pb208', mass=82*pmass+126*nmass, charge=82 }

The particles database can be displayed with the ``showdb`` method at any time from any beam:

.. code-block:: lua
	
	beam:showdb()  -- check that both, mymuon and pb208 are in the database.


Particle charges
----------------

The physics of \MAD is aware of particle charges. To enable the compatibility with codes like MAD-X that ignores the particle charges, the global option ``nocharge`` can be used to control the behavior of created beams as shown by the following example:

.. code-block:: lua
	
	local beam, option in MAD
	local beam1 = beam { particle="electron" } -- beam with negative charge
	print(beam1.charge, option.nocharge)       -- display: -1  false
	
	option.nocharge = true                     -- disable particle charges
	local beam2 = beam { particle="electron" } -- beam with negative charge
	print(beam2.charge, option.nocharge)       -- display:  1  true
	
	-- beam1 was created before nocharge activation...
	print(beam1.charge, option.nocharge)       -- display: -1  true

This approach ensures consistency of beams behavior during their entire lifetime. [#f2]_ 

Examples
--------



.. code-block:: lua
	
	local beam in MAD
	local lhcb1, lhcb2 in MADX
	local nmass, pmass, amass in MAD.constant
	local pbmass = 82*pmass+126*nmass
	
	-- attach a new beam with a new particle to lhcb1 and lhcb2.
	lhc1.beam = beam 'Pb208' { particle='pb208', mass=pbmass, charge=82 }
	lhc2.beam = lhc1.beam -- let sequences share the same beam...
	
	-- print Pb208 nuclei energy-mass in GeV and unified atomic mass.
	print(lhcb1.beam.mass, lhcb1.beam.mass/amass)


.. rubric:: Footnotes

.. [#f1] The ``qelect`` value is defined in the :doc:`constants` module.``
.. [#f2] The option ``rbarc`` in MAD-X is too volatile and does not ensure such consistency...
