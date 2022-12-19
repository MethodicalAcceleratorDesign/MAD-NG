Twiss
=====
.. _ch.cmd.twiss:

The ``twiss`` command provides a simple interface to compute the optical functions around an orbit on top of the ``track`` command, and the ``cofind`` command if the search for closed orbits is requested.

Command synopsis
----------------
.. _sec.twiss.synop:

The ``twiss`` command format is summarized in the :numref:`fig-twiss-synop`, including the default setup of the attributes. Most of these attributes are set to ``nil`` by default, meaning that ``twiss`` relies on the ``track`` and the ``cofind`` commands defaults.

.. code-block:: lua
	:name: fig-twiss-synop
	:caption: Synopsis of the ``twiss`` command with default setup.


	mtbl, mflw [, eidx] = twiss { 
		sequence=sequ,  -- sequence (required) 
		beam=nil, 	-- beam (or sequence.beam, required) 
		range=nil,  	-- range of tracking (or sequence.range) 
		dir=nil,  	-- s-direction of tracking (1 or -1) 
		s0=nil,  	-- initial s-position offset [m]
		X0=nil,  	-- initial coordinates (or damap(s), or beta block(s)) 
		O0=nil,  	-- initial coordinates of reference orbit 
		deltap=nil,  	-- initial deltap(s) 
		chrom=false,  	-- compute chromatic functions by finite difference 
		coupling=false, -- compute optical functions for non-diagonal modes 
		nturn=nil,  	-- number of turns to track 
		nstep=nil,  	-- number of elements to track 
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
		save=true,  	-- create mtable and save results 
		title=nil,  	-- title of mtable (default seq.name) 
		observe=0,  	-- save only in observed elements (every n turns) 
		savesel=nil,  	-- save selector (predicate) 
		savemap=nil,  	-- save damap in the column __map 
		atentry=nil,  	-- action called when entering an element 
		atslice=nil,  	-- action called after each element slices 
		atexit=nil,  	-- action called when exiting an element 
		ataper=nil,  	-- action called when checking for aperture 
		atsave=nil,  	-- action called when saving in mtable 
		atdebug=fnil,  	-- action called when debugging the element maps 
		codiff=nil,  	-- finite differences step for jacobian 
		coiter=nil,  	-- maximum number of iterations 
		cotol=nil,  	-- closed orbit tolerance (i.e.~|dX|) 
		X1=nil,  	-- optional final coordinates translation 
		info=nil,  	-- information level (output on terminal) 
		debug=nil, 	-- debug information level (output on terminal) 
		usrdef=nil,  	-- user defined data attached to the mflow 
		mflow=nil,  	-- mflow, exclusive with other attributes 
	}

The ``twiss`` command supports the following attributes:

.. _twiss.attr:

**sequence**
	The *sequence* to track. (no default, required). 
	Example: ``sequence = lhcb1``.

**beam**
	The reference *beam* for the tracking. If no beam is provided, the command looks for a beam attached to the sequence, i.e. the attribute ``seq.beam`` . [#f1]_ (default: ``nil``).
	Example: ``beam = beam 'lhcbeam' { beam-attributes }``.

**range**
	A *range* specifying the span of the sequence track. If no range is provided, the command looks for a range attached to the sequence, i.e. the attribute . (default: ``nil``). 
	Example: ``range = "S.DS.L8.B1/E.DS.R8.B1"``.

**dir**
	The :math:`s`-direction of the tracking: ``1`` forward, ``-1`` backward. (default: ``nil``). 
	Example: ``dir = - 1``.

**s0**
	A *number* specifying the initial :math:`s`-position offset. (default: ``nil``). 
	Example: ``s0 = 5000``.

**X0**
	A *mappable* (or a list of *mappable*) specifying initial coordinates ``{x,px,y,py, t,pt}``, damap, or beta0 block for each tracked object, i.e. particle or damap. The beta0 blocks are converted to damaps, while the coordinates are converted to damaps only if , damap, or beta0 block for each tracked object, i.e. particle or damap. The beta0 blocks are converted to damaps, while the coordinates are converted to damaps only if ``mapdef`` is specified, but both will use ``mapdef`` to setup the damap constructor. A closed orbit will be automatically searched for damaps built from coordinates. Each tracked object may also contain a ``beam`` to override the reference beam, and a *logical* ``nosave`` to discard this object from being saved in the mtable. (default: ``0``). 
	Example: ``X0 = { x=1e- 3, px=- 1e- 5 }``.

**O0** 
	A *mappable* specifying initial coordinates ``{x,px,y,py,t,pt}`` of the reference orbit around which X0 definitions take place. If it has the attribute ``cofind == true``, it will be used as an initial guess to search for the reference closed orbit. (default: ``0``). 
	Example: ``O0 = { x=1e- 4, px=- 2e- 5, y=- 2e- 4, py=1e- 5 }``.

**deltap**
	A *number* (or list of *number*) specifying the initial :math:`\delta_p` to convert (using the beam) and add to the ``pt`` of each tracked particle or damap. (default: ``nil``). 
	Example: ``s0 = 5000``.

**chrom**
	A *logical* specifying to calculate the chromatic functions by finite different using an extra :math:`\delta_p=` ``1e-6``. (default: ``false``). 
	Example: ``chrom = true``.

**coupling**
	A *logical* specifying to calculate the optical functions for coupling terms in the normalized forms. (default: ``false``). 
	Example: ``chrom = true``.

**nturn**
	A *number* specifying the number of turn to track. (default: ``nil``). 
	Example: ``nturn = 2``.

**nstep**
	A *number* specifying the number of element to track. A negative value will track all elements. (default: ``nil``). 
	Example: ``nstep = 1``.

**nslice**
	A *number* specifying the number of slices or an *iterable* of increasing relative positions or a *callable* ``(elm, mflw, lw)`` returning one of the two previous kind of positions to track in the elements. The arguments of the callable are in order, the current element, the tracked map flow, and the length weight of the step. This attribute can be locally overridden by the element. (default: ``nil``). 
	Example: ``nslice = 5``.

**mapdef** 
	A *logical* or a *damap* specification as defined by the :doc:`DAmap <diffmap>` module to track DA maps instead of particles coordinates. A value of ``true`` is equivalent to invoke the *damap* constructor with ``{}`` as argument. A value of ``false`` or ``nil`` will be internally forced to ``true`` for the tracking of the normalized forms. (default: ``true``). 
	Example: ``mapdef = { xy=2, pt=5 }``.

**method**
	A *number* specifying the order of integration from 1 to 8, or a *string* specifying a special method of integration. Odd orders are rounded to the next even order to select the corresponding Yoshida or Boole integration schemes. The special methods are ``simple`` (equiv. to ``DKD`` order 2), ``collim`` (equiv. to ``MKM`` order 2), and ``teapot`` (Teapot splitting order 2). (default: ``nil``). 
	Example: ``method = 'teapot'``.

**model**
	A *string* specifying the integration model, either ``'DKD'`` for *Drift-Kick-Drift* thin lens integration or ``'TKT'`` for *Thick-Kick-Thick* thick lens integration. [#f7]_ (default: ``nil``) 
	Example: ``model = 'DKD'``.

**ptcmodel**
	A *logical* indicating to use strict PTC model. [#f8]_ (default: ``nil``) 
	Example: ``ptcmodel = true``.

**implicit**
	A *logical* indicating that implicit elements must be sliced too, e.g. for smooth plotting. (default: ``nil``). 
	Example: ``implicit = true``.

**misalign**
	A *logical* indicating that misalignment must be considered. (default: ``nil``). 
	Example: ``misalign = true``.

**fringe**
	A *logical* indicating that fringe fields must be considered or a *number* specifying a bit mask to apply to all elements fringe flags defined by the element module. The value ``true`` is equivalent to the bit mask , i.e. allow all elements (default) fringe fields. (default: ``nil``). 
	Example: ``fringe = false``.

**radiate**
	A *logical* enabling or disabling the radiation or the *string* specifying the ``'average'`` type of radiation during the closed orbit search. The value ``true`` is equivalent to ``'average'`` and the value ``'quantum'`` is converted to ``'average'``. (default: ``nil``). 
	Example: ``radiate = 'average'``.

**totalpath**
	A *logical* indicating to use the totalpath for the fifth variable ``'t'`` instead of the local path. (default: ``nil``). 
	Example: ``totalpath = true``.

**save**
	A *logical* specifying to create a *mtable* and record tracking information at the observation points. The ``save`` attribute can also be a *string* specifying saving positions in the observed elements: ``"atentry"``, ``"atslice"``, ``"atexit"`` (i.e. ``true``), ``"atbound"`` (i.e. entry and exit), ``"atbody"`` (i.e. slices and exit) and ``"atall"``. (default: ``false``). 
	Example: ``save = false``.

**title**
	A *string* specifying the title of the *mtable*. If no title is provided, the command looks for the name of the sequence, i.e. the attribute ``seq.name``. (default: ``nil``). 
	Example: ``title = "track around IP5"``.

**observe**
	A *number* specifying the observation points to consider for recording the tracking information. A zero value will consider all elements, while a positive value will consider selected elements only, checked with method :meth:`:is_observed`, every :math:`>0` turns. (default: ``nil``). 
	Example: ``observe = 1``.

**savesel**
	A *callable* ``(elm, mflw, lw, islc)`` acting as a predicate on selected elements for observation, i.e. the element is discarded if the predicate returns ``false``. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: :const:`fnil`) 
	Example: ``savesel = \LMB e -> mylist[e.name] ~= nil``.

**savemap**
	A *logical* indicating to save the damap in the column ``__map`` of the *mtable*. (default: ``nil``). 
	Example: ``savemap = true``.

**atentry**
	 A *callable* ``(elm, mflw, 0, - 1)`` invoked at element entry. The arguments are in order, the current element, the tracked map flow, zero length and the slice index ``-1``. (default: :const:`fnil`). 
	 Example: ``atentry = myaction``.

**atslice**
	A *callable* ``(elm, mflw, lw, islc)`` invoked at element slice. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: :const:`fnil`). 
	Example: ``atslice = myaction``.

**atexit** 
	A *callable* ``(elm, mflw, 0, - 2)`` invoked at element exit. The arguments are in order, the current element, the tracked map flow, zero length and the slice index . (default: :const:`fnil`). 
	Example: ``atexit = myaction``.

**ataper**
	A *callable* ``(elm, mflw, lw, islc)`` invoked at element aperture checks, by default at last slice. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. If a particle or a damap hits the aperture, then its ``status~=~"lost"`` and it is removed from the list of tracked items. (default: :const:`fnil`). 
	Example: ``ataper = myaction``.

**atsave**
	A *callable* ``(elm, mflw, lw, islc)`` invoked at element saving steps, by default at exit. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: :const:`fnil`). 
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
	 A *number* specifying the information level to control the verbosity of the output on the console. (default: ``nil``). 
	 Example: ``info = 2``.

**debug**
	 A *number* specifying the debug level to perform extra assertions and to control the verbosity of the output on the console. (default: ``nil``). 
	 Example: ``debug = 2``.

**usrdef**
	Any user defined data that will be attached to the tracked map flow, which is internally passed to the elements method ``:track`` and to their underlying maps. (default: ``nil``). 
	Example: ``usrdef = { myvar=somevalue }``.

**mflow** 
	A *mflow* containing the current state of a ``track`` command. If a map flow is provided, all attributes are discarded except ``nstep``, ``info`` and ``debug``, as the command was already set up upon its creation. (default: ``nil``). 
	Example: ``mflow = mflow0``.


The ``twiss`` command returns the following objects in this order:

**mtbl}** A *mtable* corresponding to the augmented TFS table of the ``track`` command with the ``twiss`` command columns.

**mflw** A *mflow* corresponding to the augmented map flow of the ``track`` command with the ``twiss`` command data.

**eidx**
	 An optional *number* corresponding to the last tracked element index in the sequence when ``nstep`` was specified and stopped the command before the end of the ``range``.


Twiss mtable
------------
.. _sec.track.mtable:

The ``twiss`` command returns a *mtable* where the information described hereafter is the default list of fields written to the TFS files. [#f2]_ 

The header of the *mtable* contains the fields in the default order: [#f3]_ 

	**name**
	 The name of the command that created the ``"track"``.
	**type**
	 The type of the ``"track"``.
	**title**
	 The value of the command attribute ``title``.
	**origin**
	 The origin of the application that created the ``"MAD 1.0.0 OSX 64"``.
	**date**
	 The date of the creation of the ``"27/05/20"``.
	**time**
	 The time of the creation of the ``"19:18:36"``.
	**refcol**
	 The reference *column* for the *mtable* dictionnary, e.g. ``"name"``.
	**direction**
	 The value of the command attribute ``dir``.
	**observe**
	 The value of the command attribute ``observe``.
	**implicit**
	 The value of the command attribute ``implicit``.
	**misalign**
	 The value of the command attribute ``misalign``.
	**deltap**
	 The value of the command attribute ``deltap``.
	**lost**
	 The number of lost particle(s) or damap(s).
	**chrom**
	 The value of the command attribute ``chrom``.
	**coupling**
	 The value of the command attribute ``coupling``.
	**length**
	 The :math:`s`-length of the tracked design orbit.
	**q1**
	 The tunes of mode 1.
	**q2**
	 The tunes of mode 2.
	**q3**
	 The tunes of mode 3.
	**alfap**
	 The momentum compaction factor :math:`\alpha_p`.
	**etap**
	 The phase slip factor :math:`\eta_p`.
	**gammatr**
	 The energy gamma transition :math:`\gamma_{\text{tr}}`.
	**synch_1**
	 The first synchroton radiation integral.
	**synch_2**
	 The second synchroton radiation integral.
	**synch_3**
	 The third synchroton radiation integral.
	**synch_4**
	 The fourth synchroton radiation integral.
	**synch_5**
	 The fifth synchroton radiation integral.
	**synch_6**
	 The sixth synchroton radiation integral.
	**synch_8**
	 The eighth synchroton radiation integral.
	**range**
	 The value of the command attribute ``range``. [#f4]_ 
	**__seq**
	 The *sequence* from the command attribute ``sequence``. [#f5]_ .. _ref.twiss.mtbl1}:

The core of the *mtable* contains the columns in the default order: [#f6]_

	**name**
	 The name of the element.
	**kind**
	 The kind of the element.
	**s**
	 The :math:`s`-position at the end of the element slice.
	**l**
	 The length from the start of the element to the end of the element slice.
	**id**
	 The index of the particle or damap as provided in ``X0``.
	**x**
	 The local coordinate :math:`x` at the :math:`s`-position .
	**px**
	 The local coordinate :math:`p_x` at the :math:`s`-position.
	**y**
	 The local coordinate :math:`y` at the :math:`s`-position.
	**py**
	 The local coordinate :math:`p_y` at the :math:`s`-position.
	**t**
	 The local coordinate :math:`t` at the :math:`s`-position.
	**pt**
	 The local coordinate :math:`p_t` at the :math:`s`-position.
	**slc**
	 The slice index ranging from ``- 2`` to ``nslice``.
	**turn**
	 The turn number.
	**tdir**
	 The :math:`t`-direction of the tracking in the element.
	**eidx**
	 The index of the element in the sequence.
	**status**
	 The status of the particle or damap.
	**alfa11**
	 The optical function :math:`\alpha` of mode 1 at the :math:`s`-position.
	**beta11**
	 The optical function :math:`\beta` of mode 1 at the :math:`s`-position.
	**gama11**
	 The optical function :math:`\gamma` of mode 1 at the :math:`s`-position.
	**mu1**
	 The phase advance :math:`\mu` of mode 1 at the :math:`s`-position.
	**dx**
	 The dispersion function of :math:`x` at the :math:`s`-position.
	**dpx**
	 The dispersion function of :math:`p_x` at the :math:`s`-position.
	**alfa22**
	 The optical function :math:`\alpha` of mode 2 at the :math:`s`-position.
	**beta22**
	 The optical function :math:`\beta` of mode 2 at the :math:`s`-position.
	**gama22**
	 The optical function :math:`\gamma` of mode 2 at the :math:`s`-position.
	**mu2**
	 The phase advance :math:`\mu` of mode 2 at the :math:`s`-position.
	**dy**
	 The dispersion function of :math:`y` at the :math:`s`-position.
	**dpy**
	 The dispersion function of :math:`p_y` at the :math:`s`-position.
	**alfa33**
	 The optical function :math:`\alpha` of mode 3 at the :math:`s`-position.
	**beta33**
	 The optical function :math:`\beta` of mode 3 at the :math:`s`-position.
	**gama33**
	 The optical function :math:`\gamma` of mode 3 at the :math:`s`-position.
	**mu3**
	 The phase advance :math:`\mu` of mode 3 at the :math:`s`-position.
	**__map**
	 The damap at the :math:`s`-position. [#f5]_

The ``chrom`` attribute will add the following fields to the *mtable* header:

	**dq1**
	 The chromatic derivative of tunes of mode 1, i.e. chromaticities.
	**dq2**
	 The chromatic derivative of tunes of mode 2, i.e. chromaticities.
	**dq3**
	 The chromatic derivative of tunes of mode 3, i.e. chromaticities.

The ``chrom`` attribute will add the following columns to the *mtable*:

	**dmu1**
	 The chromatic derivative of the phase advance of mode 1 at the :math:`s`-position.
	**ddx**
	 The chromatic derivative of the dispersion function of :math:`x` at the :math:`s`-position.
	**ddpx**
	 The chromatic derivative of the dispersion function of :math:`p_x` at the :math:`s`-position.
	**wx**
	 The chromatic amplitude function of mode 1 at the :math:`s`-position.
	**phix**
	 The chromatic phase function of mode 1 at the :math:`s`-position.
	**dmu2**
	 The chromatic derivative of the phase advance of mode 2 at the :math:`s`-position.
	**ddy**
	  The chromatic derivative of the dispersion function of :math:`y` at the :math:`s`-position.
	**ddpy**
	 The chromatic derivative of the dispersion function of :math:`p_y` at the :math:`s`-position.
	**wy**
	 The chromatic amplitude function of mode 2 at the :math:`s`-position.
	**phiy**
	 The chromatic phase function of mode 2 at the :math:`s`-position.

The ``coupling`` attribute will add the following columns to the *mtable*:

	**alfa12**
	 The optical function :math:`\alpha` of coupling mode 1-2 at the :math:`s`-position.
	**beta12**
	 The optical function :math:`\beta` of coupling mode 1-2 at the :math:`s`-position.
	**gama12**
	 The optical function :math:`\gamma` of coupling mode 1-2 at the :math:`s`-position.
	**alfa13**
	 The optical function :math:`\alpha` of coupling mode 1-3 at the :math:`s`-position.
	**beta13**
	 The optical function :math:`\beta` of coupling mode 1-3 at the :math:`s`-position.
	**gama13**
	 The optical function :math:`\gamma` of coupling mode 1-3 at the :math:`s`-position.
	**alfa21**
	 The optical function :math:`\alpha` of coupling mode 2-1 at the :math:`s`-position.
	**beta21**
	 The optical function :math:`\beta` of coupling mode 2-1 at the :math:`s`-position.
	**gama21**
	 The optical function :math:`\gamma` of coupling mode 2-1 at the :math:`s`-position.
	**alfa23**
	 The optical function :math:`\alpha` of coupling mode 2-3 at the :math:`s`-position.
	**beta23**
	 The optical function :math:`\beta` of coupling mode 2-3 at the :math:`s`-position.
	**gama23**
	 The optical function :math:`\gamma` of coupling mode 2-3 at the :math:`s`-position.
	**alfa31**
	 The optical function :math:`\alpha` of coupling mode 3-1 at the :math:`s`-position.
	**beta31**
	 The optical function :math:`\beta` of coupling mode 3-1 at the :math:`s`-position.
	**gama31**
	 The optical function :math:`\gamma` of coupling mode 3-1 at the :math:`s`-position.
	**alfa32**
	 The optical function :math:`\alpha` of coupling mode 3-2 at the :math:`s`-position.
	**beta32**
	 The optical function :math:`\beta` of coupling mode 3-2 at the :math:`s`-position.
	**gama32**
	 The optical function :math:`\gamma` of coupling mode 3-2 at the :math:`s`-position.


Tracking linear normal form
---------------------------

TODO

Examples
--------

TODO


.. rubric:: Footnotes

.. [#f1] Initial coordinates ``X0`` may override it by providing a beam per particle or damap. 
.. [#f7] The ``TKT`` scheme (Yoshida) is automatically converted to the ``MKM`` scheme (Boole) when appropriate.
.. [#f8] In all cases, MAD-NG uses PTC setup ``time=true, exact=true``.
.. [#f2] The output of mtable in TFS files can be fully customized by the user.
.. [#f3] The fields from ``name`` to ``lost`` set by the ``track`` command
.. [#f4] This field is not saved in the TFS table by default.
.. [#f5] Fields and columns starting with two underscores are protected data and never saved to TFS files.
.. [#f6] The column from ``name`` to ``status`` are set by the ``track`` command.
