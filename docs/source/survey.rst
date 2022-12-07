Survey
======
.. _ch.cmd.survey:

The ``survey`` command provides a simple interface to the *geometric* tracking code. [#f1]_ The geometric tracking can be used to place the elements of a sequence in the global reference system in :numref:`fig-phy-grs`.

.. _fig.survey.synop:

.. code-block:: lua
	:caption: Synopsis of the ``survey`` command with default setup.

	mtbl, mflw [, eidx] = survey { 
		sequence=sequ,  -- sequence (required) 
		range=nil,  	-- range of tracking (or sequence.range) 
		dir=1,  	-- s-direction of tracking (1 or -1) 
		s0=0,  		-- initial s-position offset [m] 
		X0=0,  		-- initial coordinates x, y, z [m] 
		A0=0,  		-- initial angles theta, phi, psi [rad] or matrix W0 
		nturn=1,  	-- number of turns to track 
		nstep=-1,  	-- number of elements to track 
		nslice=1,  	-- number of slices (or weights) for each element 
		implicit=false,  -- slice implicit elements too (e.g.~plots) 
		misalign=false,  -- consider misalignment 
		save=true,  	-- create mtable and save results 
		title=nil,  	-- title of mtable (default seq.name) 
		observe=0,  	-- save only in observed elements (every n turns) 
		savesel=fnil,  	-- save selector (predicate) 
		savemap=false,  -- save the orientation matrix W in the column __map 
		atentry=fnil,  	-- action called when entering an element 
		atslice=fnil,  	-- action called after each element slices 
		atexit=fnil,  	-- action called when exiting an element 
		atsave=fnil,  	-- action called when saving in mtable 
		atdebug=fnil,  	-- action called when debugging the element maps 
		info=nil,  	-- information level (output on terminal) 
		debug=nil, 	-- debug information level (output on terminal) 
		usrdef=nil,  	-- user defined data attached to the mflow 
		mflow=nil,  	-- mflow, exclusive with other attributes except nstep 
	}

Command synopsis
----------------
.. _sec.survey.synop:


The ``survey`` command format is summarized in :numref:`fig.survey.synop`, including the default setup of the attributes. The ``survey`` command supports the following attributes:

.. _survey.attr:

**sequence**
	 The *sequence* to survey. (no default, required). 
 	 Example: ``sequence = lhcb1``.

**range** 
	 A *range* specifying the span of the sequence survey. If no range is provided, the command looks for a range attached to the sequence, i.e. the attribute . (default: ``nil``). 
	 Example: ``range = "S.DS.L8.B1/E.DS.R8.B1"``.

**dir** 
	 The :math:`s`-direction of the tracking: ``1`` forward, ``-1`` backward. (default: ``1``). 
	 Example: ``dir = - 1``.

**s0** 
	 A *number* specifying the initial :math:`s`-position offset. (default: ``0`` [m]). 
	 Example: ``s0 = 5000``.

**X0** 
	 A *mappable* specifying the initial coordinates ``{x,y,z}``. (default: ``0`` [m]). 
	 Example: ``X0 = { x=100, y=- 50 }``

**A0** 
	 A *mappable* specifying the initial angles ``theta``, ``phi`` and ``psi`` or an orientation *matrix* ``W0``. [#f2]_ ``W=matrix(3):rotzxy(- phi,theta,psi)}.`` (default: ``0`` [rad]). 
	 Example: ``A0 = { theta=deg2rad(30) }``

**nturn** 
	 A *number* specifying the number of turn to track. (default: ``2``). 
	 Example: ``nturn = 2``.

**nstep** 
	 A *number* specifying the number of element to track. A negative value will track all elements. (default: ``-1``). 
	 Example: ``nstep = 1``.

**nslice** 
	 A *number* specifying the number of slices or an *iterable* of increasing relative positions or a *callable* ``(elm, mflw, lw)`` returning one of the two previous kind of positions to track in the elements. The arguments of the callable are in order, the current element, the tracked map flow, and the length weight of the step. This attribute can be locally overridden by the element. (default: ``1``). 
	 Example: ``nslice = 5``.

**implicit** 
	 A *log* indicating that implicit elements must be sliced too, e.g. for smooth plotting. (default: ``false``). 
	 Example: ``implicit = true``.

**misalign** 
	 A *log* indicating that misalignment must be considered. (default: ``true``). 
	 Example: ``implicit = true``.

**save** 
	 A *log* specifying to create a *mtable* and record tracking information at the observation points. The ``save`` attribute can also be a *string* specifying saving positions in the observed elements: ``"atentry"``, ``"atslice"``, ``"atexit"`` (i.e. ``true``), ``"atbound"`` (i.e. entry and exit), ``"atbody"`` (i.e. slices and exit) and ``"atall"``. (default: ``true``). 
	 Example: ``save = false``.

**title** 
	 A *string* specifying the title of the *mtable*. If no title is provided, the command looks for the name of the sequence, i.e. the attribute ``seq.name``. (default: ``nil``). 
	 Example: ``title = "Survey around IP5"``.

**observe** 
	 A *number* specifying the observation points to consider for recording the tracking information. A zero value will consider all elements, while a positive value will consider selected elements only, checked with method ``:is_observed``, every ``observe`` :math:`>0` turns. (default: ``0``). 
	 Example: ``observe = 1``.

**savesel** 
	 A *callable* ``(elm, mflw, lw, islc)`` acting as a predicate on selected elements for observation, i.e. the element is discarded if the predicate returns ``false``. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: ``fnil``) 
	 Example: ``savesel = \LMB e -> mylist[e.name] ~= nil``.

**savemap** 
	 A *log* indicating to save the orientation matrix ``W`` in the column ``__map`` of the *mtable*. (default: ``false``). 
	 Example: ``savemap = true``.

**atentry** 
	 A *callable* ``(elm, mflw, 0, - 1)`` invoked at element entry. The arguments are in order, the current element, the tracked map flow, zero length and the slice index . (default: ``fnil``). 
	 Example: ``atentry = myaction``.

**atslice** 
	 A *callable* ``(elm, mflw, lw, islc)`` invoked at element slice. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: ``fnil``). 
	 Example: ``atslice = myaction``.

**atexit** 
	 A *callable* ``(elm, mflw, 0, - 2)`` invoked at element exit. The arguments are in order, the current element, the tracked map flow, zero length and the slice index . (default: ``fnil``). 
	 Example: ``atexit = myaction``.

**atsave** 
	 A *callable* ``(elm, mflw, lw, islc)`` invoked at element saving steps, by default at exit. The arguments are in order, the current element, the tracked map flow, the length weight of the slice and the slice index. (default: ``fnil``). 
	 Example: ``atsave = myaction``.

**atdebug** 
	 A *callable* ``(elm, mflw, lw, [msg], [...])`` invoked at the entry and exit of element maps during the integration steps, i.e. within the slices. The arguments are in order, the current element, the tracked map flow, the length weight of the integration step and a *string* specifying a debugging message, e.g. ``"map_name:0"`` for entry and ``":1"`` for exit. If the level ``debug`` :math:`\geq 4` and ``atdebug`` is not specified, the default *function* ``mdump`` is used. In some cases, extra arguments could be passed to the method. (default: ``fnil`` ). 
	 Example: ``atdebug = myaction`` .
	 
**info**
	 A *number* specifying the information level to control the verbosity of the output on the console. (default: ``nil``). 
	 Example: ``info = 2``.

**debug**
	 A *number* specifying the debug level to perform extra assertions and to control the verbosity of the output on the console. (default: ``nil``). 
	 Example: ``debug = 2``.

**usrdef** 
	 Any user defined data that will be attached to the tracked map flow, which is internally passed to the elements method :meth:`:survey` and to their underlying maps. (default: ``nil``). 
	 Example: ``usrdef = { myvar=somevalue }``.

**mflow** 
	 A *mflow* containing the current state of a ``survey`` command. If a map flow is provided, all attributes are discarded except ``nstep``, ``info`` and ``debug``, as the command was already set up upon its creation. (default: ``nil``). 
	 Example: ``mflow = mflow0``.


The ``survey`` command returns the following objects in this order:

**mtbl** 
	A *mtable* corresponding to the TFS table of the ``survey`` command.

**mflw** 
	A *mflow* corresponding to the map flow of the ``survey`` command.

**eidx**
	 An optional *number* corresponding to the last surveyed element index in the sequence when ``nstep`` was specified and stopped the command before the end of the ``range``.


Survey mtable
-------------
.. _sec.survey.mtable:

The ``survey`` command returns a *mtable* where the information described hereafter is the default list of fields written to the TFS files. [#f3]_ 



**name**
	 The name of the command that created the ``"survey"``.
**type**
	 The type of the ``"survey"``.
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
**range**
	 The value of the command attribute ``range``. [#f4]_ 
**__seq**
	 The *sequence* from the command attribute ``sequence``. [#f5]_



**name**
	 The name of the element.
**kind**
	 The kind of the element.
**s**
	 The :math:`s`-position at the end of the element slice.
**l**
	 The length from the start of the element to the end of the element slice.
**angle**
	 The angle from the start of the element to the end of the element slice.
**tilt**
	 The tilt of the element.
**x**
	 The global coordinate :math:`x` at the :math:`s`-position.
**y**
	 The global coordinate :math:`y` at the :math:`s`-position.
**z**
	 The global coordinate :math:`z` at the :math:`s`-position.
**theta**
	 The global angle :math:`\theta` at the :math:`s`-position.
**phi**
	 The global angle :math:`\phi` at the :math:`s`-position.
**psi**
	 The global angle :math:`\psi` at the :math:`s`-position.
**slc**
	 The slice number ranging from ``- 2`` to ``nslice``.
**turn**
	 The turn number.
**tdir**
	 The :math:`t`-direction of the tracking in the element.
**eidx**
	 The index of the element in the sequence.
**__map**
	 The orientation *matrix* at the :math:`s`-position. [#f5]_


Geometrical tracking
--------------------

:numref:`fig.survey.trkslc` presents the scheme of the geometrical tracking through an element sliced with ``nslice=3``. The actions ``atentry`` (index ``- 1``), ``atslice`` (indexes ``0..3``), and ``atexit`` (index ``- 2``) are reversed between the forward tracking (``dir=1`` with increasing :math:`s`-position) and the backward tracking (``dir=- 1`` with decreasing :math:`s`-position). By default, the action ``atsave`` is attached to the exit slice, and hence it is also reversed in the backward tracking.


.. _fig.survey.trkslc:
.. figure:: fig/dyna-trck-slice-crop.png

	Geometrical tracking with slices.

Slicing
"""""""



	#.	 A *number* of the form ``nslice=``:math:`N` that specifies the number of slices with indexes ``0..N``. This defines a uniform slicing with slice length :math:`l_{\text{slice}} = l_{\text{elem}}/N`.

	#.	 An *iterable* of the form ``nslice={lw_1,lw_2,..lw_N}`` with :math:`\sum_i lw_i=1` that specifies the fraction of length of each slice with indexes ``0..N`` where :math:`N`=\ ``#nslice``. This defines a non-uniform slicing with a slice length of :math:`l_i = lw_i\times l_{\text{elem}}`.

	#.	 A *callable* ``(elm, mflw, lw)`` returning one of the two previous forms of slicing. The arguments are in order, the current element, the tracked map flow, and the length weight of the step, which should allow to return a user-defined element-specific slicing. 


The surrounding ``P`` and ``P``\ :math:`^{-1}` maps represent the patches applied around the body of the element to change the frames, after the ``atentry`` and before the ``atexit`` actions:

	#.	 The misalignment of the element to move from the *global frame* to the *element frame* if the command attribute ``misalign`` is set to ``true``.

	#.	 The tilt of the element to move from the element frame to the *titled frame* if the element attribute ``tilt`` is non-zero. The ``atslice`` actions take place in this frame.



Sub-elements
""""""""""""

The ``survey`` command takes sub-elements into account, mainly for compatibility with the ``track`` command. In this case, the slicing specification is taken between sub-elements, e.g. 3 slices with 2 sub-elements gives a final count of 9 slices. It is possible to adjust the number of slices between sub-elements with the third form of slicing specifier, i.e. by using a callable where the length weight argument is between the current (or the end of the element) and the last sub-elements (or the start of the element).

Examples
--------



.. rubric:: Footnotes

.. [#f1] MAD-NG implements only two tracking codes denominated the *geometric* and *dynamic* tracking
.. [#f2] An orientation matrix can be obtained from the 3 angles with ``W=matrix(3):rotzxy(- phi,theta,psi)``
.. [#f3] The output of mtable in TFS files can be fully customized by the user.
.. [#f4] This field is not saved in the TFS table by default.
.. [#f5] Fields and columns starting with two underscores are protected data and never saved to TFS files.
