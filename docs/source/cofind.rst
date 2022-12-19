Cofind
======
.. _ch.cmd.cofind:

The ``cofind`` command (i.e. closed orbit finder) provides a simple interface to find a closed orbit using the Newton algorithm on top of the ``track`` command.

Command synopsis
----------------



.. code-block:: lua
	:caption: Synopsis of the ``cofind`` command with default setup.
	:name: fig-cofind-synop

	mtbl, mflw = cofind} { 
		sequence=sequ,	-- sequence (required) 
		beam=nil, 	-- beam (or sequence.beam, required) 
		range=nil,  	-- range of tracking (or sequence.range) 
		dir=nil,  	-- s-direction of tracking (1 or -1) 
		s0=nil,  	-- initial s-position offset [m]
		X0=nil,  	-- initial coordinates (or damap, or beta block) 
		O0=nil,  	-- initial coordinates of reference orbit 
		deltap=nil,  	-- initial deltap(s) 
		nturn=nil,  	-- number of turns to track 
		nslice=nil,  	-- number of slices (or weights) for each element 
		mapdef=true,  	-- setup for damap (or list of, true => {}) 
		method=nil,  	-- method or order for integration (1 to 8) 
		model=nil,  	-- model for integration ('DKD' or 'TKT') 
		ptcmodel=nil,  	-- use strict PTC thick model (override option) 
		implicit=nil,  	-- slice implicit elements too (e.g. plots) 
		misalign=nil,  	-- consider misalignment 
		fringe=nil,  	-- enable fringe fields (see element.flags.fringe) 
		radiate=nil,  	-- radiate at slices 
		totalpath=nil,  -- variable 't' is the totalpath 
		save=false,  	-- create mtable and save results 
		title=nil,  	-- title of mtable (default seq.name) 
		observe=nil,  	-- save only in observed elements (every n turns) 
		savesel=nil,  	-- save selector (predicate) 
		savemap=nil,  	-- save damap in the column __map 
		atentry=nil,  	-- action called when entering an element 
		atslice=nil,  	-- action called after each element slices 
		atexit=nil,  	-- action called when exiting an element 
		ataper=nil,  	-- action called when checking for aperture 
		atsave=nil,  	-- action called when saving in mtable 
		atdebug=fnil,  	-- action called when debugging the element maps 
		codiff=1e-10,	-- finite differences step for jacobian 
		coiter=20,  	-- maximum number of iterations 
		cotol=1e-8,  	-- closed orbit tolerance (i.e.~|dX|) 
		X1=0,  		-- optional final coordinates translation 
		info=nil,  	-- information level (output on terminal) 
		debug=nil, 	-- debug information level (output on terminal) 
		usrdef=nil,  	-- user defined data attached to the mflow 
		mflow=nil,  	-- mflow, exclusive with other attributes 
	}

The ``cofind`` command format is summarized in :numref:`fig-cofind-synop`, including the default setup of the attributes. Most of these attributes are set to :const:`nil` by default, meaning that ``cofind`` relies on the ``track`` command defaults.

The ``cofind`` command supports the following attributes:

.. _cofind.attr:
	**sequence**
		The *sequence* to track. (no default, required). 
		Example: ``sequence = lhcb1``.

	**beam**
		The reference *beam* for the tracking. If no beam is provided, the command looks for a beam attached to the sequence, i.e. the attribute ``seq.beam``.
		Example: ``beam = beam 'lhcbeam' { beam-attributes }``. [#f1]_

	**range**
		A *range* specifying the span of the sequence track. If no range is provided, the command looks for a range attached to the sequence, i.e. the attribute . (default: :const:`nil`). 
		Example: ``range = "S.DS.L8.B1/E.DS.R8.B1"``.

	**dir**
		The :math:`s`-direction of the tracking: ``1`` forward, ``-1`` backward. (default: :const:`nil`). 
		Example: ``dir = - 1``.

	**s0**
		A *number* specifying the initial :math:`s`-position offset. (default: :const:`nil`). 
		Example: ``s0 = 5000``.

	**X0**
		A *mappable* (or a list of *mappable*) specifying initial coordinates ``{x,px,y,py, t,pt}``, damap, or beta block for each tracked object, i.e. particle or damap. The beta blocks are converted to damaps, while the coordinates are converted to damaps only if ``mapdef`` is specified, but both will use ``mapdef`` to setup the damap constructor. Each tracked object may also contain a ``beam`` to override the reference beam, and a *logical* ``nosave`` to discard this object from being saved in the mtable. (default: :const:`nil`). 
		Example: ``X0 = { x=1e- 3, px=- 1e- 5 }``.

	**O0**
		A *mappable* specifying initial coordinates ``{x,px,y,py,t,pt}`` of the reference orbit around which X0 definitions take place. If it has the attribute ``cofind == true``, it will be used as an initial guess to search for the reference closed orbit. (default: ``0`` ). 
		Example: ``O0 = { x=1e- 4, px=- 2e- 5, y=- 2e- 4, py=1e- 5 }``.

	**deltap**
		A *number* (or list of *number*) specifying the initial :math:`\delta_p` to convert (using the beam) and add to the ``pt`` of each tracked particle or damap. (default::const:`nil`). 
		Example: ``s0 = 5000``.

	**nturn**
		A *number* specifying the number of turn to track. (default: :const:`nil`). 
		Example: ``nturn = 2``.

	**nstep**
		A *number* specifying the number of element to track. A negative value will track all elements. (default: :const:`nil`). 
		Example: ``nstep = 1``.

	**nslice**
		A *number* specifying the number of slices or an *iterable* of increasing relative positions or a *callable* ``(elm, mflw, lw)`` returning one of the two previous kind of positions to track in the elements. The arguments of the callable are in order, the current element, the tracked map flow, and the length weight of the step. This attribute can be locally overridden by the element. (default: :const:`nil`). 
		Example: ``nslice = 5``.

	**mapdef**
		A *logical* or a *damap* specification as defined by the :doc:`DAmap <diffmap>` module to track DA maps instead of particles coordinates. A value of ``true`` is equivalent to invoke the *damap* constructor with ``{}`` as argument. A value of ``false`` or :const:`nil` disable the use of damaps and force ``cofind`` to replace each particles or damaps by seven particles to approximate their Jacobian by finite difference. (default: ``true``). 
		Example: ``mapdef = { xy=2, pt=5 }``.

	**method**
		A *number* specifying the order of integration from 1 to 8, or a *string* specifying a special method of integration. Odd orders are rounded to the next even order to select the corresponding Yoshida or Boole integration schemes. The special methods are ``simple`` (equiv. to ``DKD`` order 2), ``collim`` (equiv. to ``MKM`` order 2), and ``teapot`` (Teapot splitting order 2). (default: :const:`nil`). 
		Example: ``method = 'teapot'``.

	**model**
		A *string* specifying the integration model, either ``'DKD'`` for *Drift-Kick-Drift* thin lens integration or ``'TKT'`` for *Thick-Kick-Thick* thick lens integration. [#f2]_ (default: :const:`nil`) 
		Example: ``model = 'DKD'``.

	**ptcmodel**
		A *logical* indicating to use strict PTC model. [#f3]_ (default: :const:`nil`) 
		Example: ``ptcmodel = true``.

	**implicit**
		A *logical* indicating that implicit elements must be sliced too, e.g. for smooth plotting. (default: :const:`nil`). 
		Example: ``implicit = true``.

	**misalign**
		A *logical* indicating that misalignment must be considered. (default: :const:`nil`). 
		Example: ``misalign = true``.

	**fringe**
		A *logical* indicating that fringe fields must be considered or a *number* specifying a bit mask to apply to all elements fringe flags defined by the element module. The value ``true`` is equivalent to the bit mask , i.e. allow all elements (default) fringe fields. (default: :const:`nil`). 
		Example: ``fringe = false``.

	**radiate**
		A *logical* enabling or disabling the radiation or the *string* specifying the ``'average'`` type of radiation during the closed orbit search. The value ``true`` is equivalent to ``'average'`` and the value ``'quantum'`` is converted to ``'average'``. (default: :const:`nil`). 
		Example: ``radiate = 'average'``.

	**totalpath**
		A *logical* indicating to use the totalpath for the fifth variable ``'t'`` instead of the local path. (default: :const:`nil`). 
		Example: ``totalpath = true``.

	**save**
		A *logical* specifying to create a *mtable* and record tracking information at the observation points. The ``save`` attribute can also be a *string* specifying saving positions in the observed elements: ``"atentry"``, ``"atslice"``, ``"atexit"`` (i.e. ``true``), ``"atbound"`` (i.e. entry and exit), ``"atbody"`` (i.e. slices and exit) and ``"atall"``. (default: ``false``). 
		Example: ``save = false``.

	**title**
		A *string* specifying the title of the *mtable*. If no title is provided, the command looks for the name of the sequence, i.e. the attribute ``seq.name``. (default: :const:`nil`). 
		Example: ``title = "track around IP5"``.

	**observe**
		A *number* specifying the observation points to consider for recording the tracking information. A zero value will consider all elements, while a positive value will consider selected elements only, checked with method :meth:`:is_observed`, every :math:`>0` turns. (default: :const:`nil`). 
		Example: ``observe = 1``.

	**savesel**
		A *callable* ``(elm, mflw, lw, islc)`` acting as a predicate on selected elements for observation, i.e. the element is discarded if the predicate returns ``false``. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: :const:`nil`) 
		Example: ``savesel = \LMB e -> mylist[e.name] ~= nil``.

	**savemap**
		A *logical* indicating to save the damap in the column ``__map`` of the *mtable*. (default: :const:`nil`). 
		Example: ``savemap = true``.

	**atentry**
		A *callable* ``(elm, mflw, 0, - 1)`` invoked at element entry. The arguments are in order, the current element, the tracked map flow, zero length and the slice index ``-1``. (default: :const:`nil`). 
		Example: ``atentry = myaction``.

	**atslice**
		A *callable* ``(elm, mflw, lw, islc)`` invoked at element slice. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: :const:`nil`). 
		Example: ``atslice = myaction``.

	**atexit** 
		A *callable* ``(elm, mflw, 0, - 2)`` invoked at element exit. The arguments are in order, the current element, the tracked map flow, zero length and the slice index . (default: :const:`nil`). 
		Example: ``atexit = myaction``.

	**ataper**
		A *callable* ``(elm, mflw, lw, islc)`` invoked at element aperture checks, by default at last slice. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. If a particle or a damap hits the aperture, then its ``status~=~"lost"`` and it is removed from the list of tracked items. (default: :const:`fnil`). 
		Example: ``ataper = myaction``.

	**atsave**
		A *callable* ``(elm, mflw, lw, islc)`` invoked at element saving steps, by default at exit. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: :const:`nil`). 
		Example: ``atsave = myaction``.

	**atdebug**
		A *callable* ``(elm, mflw, lw, [msg], [...])`` invoked at the entry and exit of element maps during the integration steps, i.e. within the slices. The arguments are in order, the current element, the tracked map flow, the length weight of the integration step and a *string* specifying a debugging message, e.g. ``"map_name:0"`` for entry and ``":1"`` for exit. If the level ``debug`` :math:`\geq 4` and ``atdebug`` is not specified, the default *function* ``mdump`` is used. In some cases, extra arguments could be passed to the method. (default: :const:`fnil`). 
		Example: ``atdebug = myaction``.

	**codiff**
		A *number* specifying the finite difference step to approximate the Jacobian when damaps are disabled. If ``codiff`` is larger than :math:`100\times`\ ``cotol``, it will be adjusted to ``cotol`` :math:`/100` and a warning will be emitted. (default: ``1e- 8``). 
		Example: ``codiff = 1e- 8``.

	**coiter**
		A *number* specifying the maximum number of iteration. If this threshold is reached, all the remaining tracked objects are tagged as ``"unstable"``. (default: 20). 
		Example: ``coiter = 5``.

	**cotol**
		A *number* specifying the closed orbit tolerance. If all coordinates update of a particle or a damap are smaller than ``cotol``, then it is tagged as ``"stable"``. (default: ``1e-8``). 
		Example: ``cotol = 1e- 6``.

	**X1**
		A *mappable* specifying the coordinates ``{x,px,y,py,t,pt}`` to *subtract* to the final coordinates of the particles or the damaps. (default: ``0``). 
		Example: ``X1 = { t=100, pt=10 }``.

	**info**
		A *number* specifying the information level to control the verbosity of the output on the console. (default: :const:`nil`). 
		Example: ``info = 2``.

	**debug**
		A *number* specifying the debug level to perform extra assertions and to control the verbosity of the output on the console. (default: :const:`nil`). 
		Example: ``debug = 2``.

	**usrdef**
		Any user defined data that will be attached to the tracked map flow, which is internally passed to the elements method ``:track`` and to their underlying maps. (default: :const:`nil`). 
		Example: ``usrdef = { myvar=somevalue }``.

	**mflow** 
		A *mflow* containing the current state of a ``track`` command. If a map flow is provided, all attributes are discarded except ``nstep``, ``info`` and ``debug``, as the command was already set up upon its creation. (default: :const:`nil`). 
		Example: ``mflow = mflow0``.

The ``cofind`` command stops when all particles or damap are tagged as ``"stable"``, ``"unstable"``, ``"singular"`` or ``"lost"``. The ``cofind`` command returns the following objects in this order:

	**mtbl**
		A *mtable* corresponding to the TFS table of the ``track`` command where the ``status`` column may also contain the new values ``"stable"``, ``"unstable"`` or ``"singular"``.

	**mflw**
		A *mflow* corresponding to the map flow of the ``track`` command. The particles or damaps ``status`` are tagged and ordered by ``"stable"``, ``"unstable"``, ``"singular"``, ``"lost"`` and ``id``.

Cofind mtable
-------------
.. _sec.cofind.mtable:

The ``cofind`` command returns the ``track`` *mtable* unmodified except for the ``status`` column. The tracked objects id will appear once per iteration at the ``\$end`` marker, and other defined observation points if any, until they are removed from the list of tracked objects.

Examples
--------

TODO

.. [#f1] Initial coordinates ``X0`` may override it by providing a beam per particle or damap. 
.. [#f2] The ``TKT`` scheme (Yoshida) is automatically converted to the ``MKM`` scheme (Boole) when appropriate.
.. [#f3] In all cases, MAD-NG uses PTC setup ``time=true, exact=true``.
