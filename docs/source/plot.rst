Plot
====
.. _ch.cmd.plot:

The ``plot`` command provides a simple interface to the `Gnuplot <http://www.gnuplot.info>`_ application. The Gnuplot release 5.2 or higher must be installed and visible in the user  by MAD to be able to run this command.

Command synopsis
----------------
.. _sec.plot.synop:

.. code-block:: lua
	:name: fig-plot-synop
	:caption: Synopsis of the ``plot`` command with default setup.

	cmd = plot { 
		sid			= 1,	 -- stream id 1 <= n <= 25 (Gnuplot instances)
		data		= nil, 	 -- { x=tbl.x, y=vec } (precedence over table) 
		table		= nil, 	 -- mtable 
		tablerange	= nil, 	 -- mtable range (default table.range) 
		sequence	= nil, 	 -- seq | { seq1, seq2, ...,} | "keep" 
		range		= nil, 	 -- sequence range (default sequence.range) 
		name		= nil, 	 -- (default table.title) 
		date		= nil,   -- (default table.date) 
		time		= nil,   -- (default table.time) 
		output		= nil, 	 -- "filename" -> pdf | number -> wid 
		scrdump		= nil, 	 -- "filename" 
		survey-attributes,
		windows-attributes,
		layout-attributes,
		labels-attributes,
		axis-attributes,
		ranges-attributes,
		data-attributes,
		plots-attributes,
		custom-attributes,
		info		=nil,  	-- information level (output on terminal) 
		debug		=nil, 	-- debug information level (output on terminal) 
	}


The ``plot`` command format is summarized in :numref:`fig-plot-synop`, including the default setup of the attributes.
The ``plot`` command supports the following attributes:

.. _plot.attr:

	**info**
		A *number* specifying the information level to control the verbosity of the output on the console. (default: ``nil``). 
		Example: ``info = 2``.

	**debug**
		A *number* specifying the debug level to perform extra assertions and to control the verbosity of the output on the console. (default: ``nil``). 
		Example: ``debug = 2``.

The ``plot`` command returns itself.

.. rubric:: Footnotes

