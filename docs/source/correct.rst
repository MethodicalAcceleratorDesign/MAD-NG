Correct
=======
.. _ch.cmd.correct:

The ``correct`` command (i.e. orbit correction) provides a simple interface to compute the orbit steering correction and setup the kickers of the sequences from the analysis of their ``track`` and ``twiss`` mtables.

.. code-block:: 
	:caption: Synopsis of the ``correct`` command with default setup.
	:name: fig-correct-synop

	mlst = correct { 
		sequence=nil,	-- sequence(s) (required) 
		range=nil,  	-- sequence(s) range(s) (or sequence.range) 
		title=nil,  	-- title of mtable (default seq.name) 
		model=nil,  	-- mtable(s) with twiss functions (required) 
		orbit=nil,  	-- mtable(s) with measured orbit(s), or use model 
		target=nil,  	-- mtable(s) with target orbit(s), or zero orbit 
		kind='ring',  	--\ 'line' or 'ring' 
		plane='xy',  	--\ 'x', 'y' or 'xy' 
		method='micado',--\ 'LSQ', 'SVD' or 'MICADO' 
		ncor=0,  	-- number of correctors to consider by method, 0=all 
		tol=1e-5,  	-- rms tolerance on the orbit 
		units=1,  	-- units in [m] of the orbit 
		corcnd=false,  	-- precond of correctors using 'svdcnd' or 'pcacnd' 
		corcut=0,  	-- value to theshold singular values in precond 
		cortol=0,  	-- value to theshold correctors in svdcnd 
		corset=true,  	-- update correctors correction strengths 
		monon=false,  	-- fraction (0<?<=1) of randomly available monitors 
		moncut=false,  	-- cut monitors above moncut sigmas 
		monerr=false,  	-- 1:use mrex and mrey alignment errors of monitors 
						-- 2:use msex and msey scaling errors of monitors 
		info=nil,  	-- information level (output on terminal) 
		debug=nil, 	-- debug information level (output on terminal) 
	}

.. _sec.correct.synop:

Command synopsis
----------------

The ``correct`` command format is summarized in :numref:`fig-correct-synop`, including the default setup of the attributes.
The ``correct`` command supports the following attributes:

.. _correct.attr:

	**sequence**
	 The *sequence* (or a list of *sequence*) to analyze. (no default, required). 
	 Example: ``sequence = lhcb1``.

	**range**
	 A *range* (or a list of *range*) specifying the span of the sequence to analyze. If no range is provided, the command looks for a range attached to the sequence, i.e. the attribute . (default: ``nil``). 
	 Example: ``range = "S.DS.L8.B1/E.DS.R8.B1"``.

	**title**
	 A *string* specifying the title of the *mtable*. If no title is provided, the command looks for the name of the sequence, i.e. the attribute ``seq.name``. (default: ``nil``). 
	 Example: ``title = "Correct orbit around IP5"``.

	**model** 
	 A *mtable* (or a list of *mtable*) providing ``twiss``-like information, e.g. elements, orbits and optical functions, of the corresponding sequences. (no default, required). 
	 Example: ``model = twmtbl``.

	**orbit**
	 A *mtable* (or a list of *mtable*) providing ``track``-like information, e.g. elements and measured orbits, of the corresponding sequences. If this attribute is ``nil``, the model orbit is used. (default: ``nil``). 
	 Example: ``orbit = tkmtbl``.

	**target** 
	 A *mtable* (or a list of *mtable*) providing ``track``-like information, e.g. elements and target orbits, of the corresponding sequences. If this attribute is ``nil``, the design orbit is used. (default: ``nil``). 
	 Example: ``target = tgmtbl``.

	**kind** 
	 A *string* specifying the kind of correction to apply among ``line`` or ``ring``. The kind ``line`` takes care of the causality between monitors, correctors and sequences directions, while the kind ``ring`` considers the system as periodic. (default: ). 
	 Example: ``kind = 'line'``.

	**plane**
	 A *string* specifying the plane to correct among ``x``, , ``y`` and ``xy``. (default: ``'xy'``). 
	 Example: ``plane = 'x'``.

	**method**
	 A *string* specifying the method to use for correcting the orbit among ``LSQ``, ``SVD`` or ``micado``. These methods correspond to the solver used from the :doc:`linear algebra <linalg>` module to find the orbit correction, namely ``solve``, ``ssolve`` or ``nsolve``. (default: ``'micado'``). 
	 Example: ``method = 'svd'``.

	**ncor**
	 A *number* specifying the number of correctors to consider with the method ``micado``, zero meaning all available correctors. (default: ``0``). 
	 Example: ``ncor = 4``.

	**tol** 
	 A *number* specifying the rms tolerance on the residuals for the orbit correction. (default: 1e-6). 
	 Example: ``tol = 1e- 6``.

	**unit**
	 A *number* specifying the unit of the ``orbit`` and ``target`` coordinates. (default: ``1`` [m]). 
	 Example: ``units = 1e- 3`` [m], i.e. [mm].

	**corcnd** 
	 A *log* or a *string* specifying the method to use among ``svdcnd`` and ``pcacnd`` from the :doc:`linear algebra <linalg>` module for the preconditioning of the system. A ``true`` value corresponds to . (default: ``false``). 
	 Example: ``corcnd = 'pcacnd'``.

	**corcut** 
	 A *number* specifying the thresholds for the singular values to pass to the ``svdcnd`` and ``pcacnd`` method for the preconditioning of the system. (default: ``0``). 
	 Example: ``cortol = 1e- 6``.

	**cortol**
	 A *number* specifying the thresholds for the correctors to pass to the ``svdcnd`` method for the preconditioning of the system. (default: ``0``). 
	 Example: ``cortol = 1e- 8``.

	**corset**
	 A *log* specifying to update the correctors strengths for the corrected orbit. (default: ``true``). 
	 Example: ``corset = false``.

	**monon**
	 A *number* specifying a fraction of available monitors selected from a uniform RNG. (default: ``false``). 
	 Example: ``monon = 0.8``, keep 80% of the monitors.

	**moncut**
	 A *number* specifying a threshold in number of sigma to cut monitor considered as outliers. (default: ``false``). 
	 Example: ``moncut = 2``, cut monitors above :math:`2\sigma`.

	**monerr**
	 A *number* in ``0..3`` specifying the type of monitor reading errors to consider: ``1`` use scaling errors ``msex`` and ``msey``, ``2`` use alignment errors ``mrex``, ``mrey`` and ``dpsi``, ``3`` use both. (default: ``false``). 
	 Example: ``monerr = 3``.

	**info**
	 A *number* specifying the information level to control the verbosity of the output on the console. (default: ``nil``). 
	 Example: ``info = 2``.

	**debug**
	 A *number*\ specifying the debug level to perform extra assertions and to control the verbosity of the output on the console. (default: ``nil``). 
	 Example: ``debug = 2``.


The ``correct`` command returns the following object:

	``mlst``
	 A *mtable* (or a list of *mtable*) corresponding to the TFS table of the ``correct`` command. It is a list when multiple sequences are corrected together.


Correct mtable
--------------
.. _sec.correct.mtable:

The ``correct`` command returns a *mtable* where the information described hereafter is the default list of fields written to the TFS files. [#f1]_ 



	**name**
	 The name of the command that created the ``"correct"``.
	**type**
	 The type of the ``"correct"``.
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
	**range**
	 The value of the command attribute ``range``. [#f2]_ 
	**__seq**
	 The *sequence* from the command attribute ``sequence``. [#f3]_ .. _ref.track.mtbl1}:



	**name**
	 The name of the element.
	**kind**
	 The kind of the element.
	**s**
	 The :math:`s`-position at the end of the element slice.
	**l**
	 The length from the start of the element to the end of the element slice.
	**x_old**
	 The local coordinate :math:`x` at the :math:`s`-position before correction.
	**y_old**
	 The local coordinate :math:`y` at the :math:`s`-position before correction.
	**x**
	 The predicted local coordinate :math:`x` at the :math:`s`-position after correction.
	**y**
	 The predicted local coordinate :math:`y` at the :math:`s`-position after correction.
	**rx**
	 The predicted local residual :math:`r_x` at the :math:`s`-position after correction.
	**ry**
	 The predicted local residual :math:`r_y` at the :math:`s`-position after correction.
	**hkick_old**
	 The local horizontal kick at the :math:`s`-position before correction.
	**vkick_old**
	 The local vertical kick at the :math:`s`-position before correction.
	**hkick**
	 The predicted local horizontal kick at the :math:`s`-position after correction.
	**vkick**
	 The predicted local vertical kick at the :math:`s`-position after correction.
	**shared**
	 A *log* indicating if the element is shared with another sequence.
	**eidx**
	 The index of the element in the sequence.

Note that ``correct`` does not take into account the particles and damaps ``id``s present in the (augmented) ``track`` *mtable*, hence the provided tables should contain single particle or damap information.

Examples
--------



.. rubric:: Footnotes

.. [#f1] The output of mtable in TFS files can be fully customized by the user.
.. [#f2] This field is not saved in the TFS table by default.
.. [#f3] Fields and columns starting with two underscores are protected data and never saved to TFS files.\label{ref:track:mtbl1
