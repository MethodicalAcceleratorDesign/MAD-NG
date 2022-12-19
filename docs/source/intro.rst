Introduction
============
.. _ch.gen.intro:

Presentation
------------

The Methodical Accelerator Design -- Next Generation application is an all-in-one standalone versatile tool for particle accelerator design, modeling, and optimization, and for beam dynamics and optics studies. Its general purpose scripting language is based on the simple yet powerful Lua programming language (with a few extensions) and embeds the state-of-art Just-In-Time compiler LuaJIT. Its physics is based on symplectic integration of differential maps made out of GTPSA (Generalized Truncated Power Series). The physics of the transport maps and the normal form analysis were both strongly inspired by the PTC/FPP library from E. Forest. MAD-NG development started in 2016 by the author as a side project of MAD-X, hence MAD-X users should quickly become familiar with its ecosystem, e.g. lattices definition.

MAD-NG is free open-source software, distributed under the GNU General Public License v3. The source code, units tests, integration tests, and examples are all available on its Github `repository <https://github.com/MethodicalAcceleratorDesign/MAD>`_, including the `documentation <https://github.com/MethodicalAcceleratorDesign/MADdocs>`_ and its LaTeX source. For convenience, the binaries and few examples are also made available from the `releases repository <http://cern.ch/mad/releases/madng/>`_ located on the AFS shared file system at CERN.

Installation
------------

Download the binary corresponding to your platform from the `releases repository`_ and install it in a local directory. Update (or check) that the  and install it in a local directory. Update (or check) that the ``PATH`` environment variable contains the path to your local directory or prefix  environment variable contains the path to your local directory or prefix ``mad`` with this path to run it. Rename the application from  with this path to run it. Rename the application from ``mad-arch-v.m.n`` to  to ``mad`` and make it executable with the command ' and make it executable with the command '``chmod u+x mad``' on Unix systems or add the  extension on Windows.

.. code-block:: console
	
	$ ./mad - h 
	usage: ./mad [options]... [script [args]...]. 
	Available options are: 
		- e chunk  	Execute string 'chunk'.
		- l name   	Require library 'name'.
		- b ...    	Save or list bytecode.
		- j cmd    	Perform JIT control command.
		- O[opt]   	Control JIT optimizations.
		- i        	Enter interactive mode after executing 'script'.
		- q        	Do not show version information.
		- M        	Do not load MAD environment.
		- Mt[=num] 	Set initial MAD trace level to 'num'.
		- MT[=num] 	Set initial MAD trace level to 'num' and location.
		- E        	Ignore environment variables.
		--        	 Stop handling options.
		-         	 Execute stdin and stop handling options.

Releases version
""""""""""""""""

MAD-NG releases are tagged on the Github repository and use mangled binary names on the releases repository, i.e. ``mad-arch-v.m.n`` where:

**arch**
	 is the platform architecture for binaries among ``linux``, ``macos`` and ``windows``.
**v**
	 is the **v**\ ersion number, ``0`` meaning beta-version under active development.
**m**
	 is the **m**\ ajor release number corresponding to features completeness.
**n**
	 is the mi\ **n**\ or release number corresponding to bug fixes.


Interactive Mode
----------------

To run MAD-NG in interactive mode, just typewrite its name on the Shell invite like any command-line tool. It is recommended to wrap MAD-NG with the `readline wrapper <http://github.com/hanslub42/rlwrap>`_  in interactive mode for easier use and commands history:

.. code-block:: console

  $ rlwrap ./mad
     ____  __   ______    ______     |   Methodical Accelerator Design
      /  \/  \   /  _  \   /  _  \   |   release: 0.9.0 (OSX 64)
     /  __   /  /  /_/ /  /  /_/ /   |   support: http://cern.ch/mad
    /__/  /_/  /__/ /_/  /_____ /    |   licence: GPL3 (C) CERN 2016+
                                     |   started: 2020-08-01 20:13:51


.. code-block::
	
	> print "hello world!"
	"hello world!"

Here the application is assumed to be installed in the current directory '``.``' and the character '``>``' is the prompt waiting for user input in interactive mode. If you write an incomplete statement, the interpreter waits for its completion by issuing a different prompt:

.. code-block::
	
	> print                -- 1st level prompt, incomplete statement
	>> "hello world!"      -- 2nd level prompt, complete the statement
	hello world!           -- execute

Typing the character '``=``' right after the 1st level prompt is equivalent to call the :literal:`print` function:

.. code-block::
	
	> = "hello world!"     -- 1st level prompt followed by =
	hello world!           -- execute print "hello world!"
	> = MAD.option.numfmt
	% -.10g


To quit the application typewrite ``Crtl+D`` to send ``EOF`` (end-of-file) on the input, [#f2]_ ``Crtl+\`` to send the ``SIGQUIT`` (quit) signal, or ``Crtl+C`` to send the stronger ``SIGINT`` (interrupt) signal. If the application is stalled or looping for ever, typewriting a single ``Crtl+\`` or ``Crtl+C`` twice will stop it:

.. code-block::
	
	> while true do end    -- loop forever, 1st Crtl+C doesn't stop it
	pending interruption (*in*) VM! (next will exit)         -- 2nd Crtl+C
	interrupted!           -- application stopped
	
	> while true do end    -- loop forever, a single Crtl+\ does stop it
	Quit: 3                -- Signal 3 caught, application stopped


In interactive mode, each line input is run in its own *chunk* [#f3]_, which also rules variables scopes. Hence ``local``, variables are not visible between chunks, i.e. input lines. The simple solutions are either to use global variables or to enclose local statements into the same chunk delimited by the ``do ... end`` keywords:

.. code-block::
	
	> local a = "hello"
	> print(a.." world!")
	  stdin:1: attempt to concatenate global 'a' (a nil value)
	  stack traceback:
	  stdin:1: (*in*) main chunk
	  [C]: at 0x01000325c0
	
	> do                   -- 1st level prompt, open the chunck
	>> local a = "hello"   -- 2nd level prompt, waiting for statement completion
	>> print(a.." world!") -- same chunk, local 'a' is visible
	>> end                 -- close and execute the chunk
	hello world!
	> print(a)             -- here 'a' is an unset global variable
	nil
	> a = "hello"          -- set global 'a'
	> print(a.." world!")  -- works but pollutes the global environment
	hello world!


Batch Mode
----------

To run MAD-NG in batch mode, just run it in the shell with files as arguments on the command line:

.. code-block:: console
	
	$ ./mad [mad options] myscript1.mad myscript2.mad ...


where the scripts contains programs written in the MAD-NG programming language (see :doc:`Scripting <script>`).

Online Help
-----------

MAD-NG is equipped with an online help system [#f4]_ useful in interactive mode to quickly search for information displayed in the ``man``-like Unix format :

.. code-block:: console


	> help()
    Related topics:
    MADX, aperture, beam, cmatrix, cofind, command, complex, constant, correct,
    ctpsa, cvector, dynmap, element, filesys, geomap, gfunc, gmath, gphys, gplot,
    gutil, hook, lfun, linspace, logrange, logspace, match, matrix, mflow,
    monomial, mtable, nlogrange, nrange, object, operator, plot, range, reflect,
    regex, sequence, strict, survey, symint, symintc, tostring, totable, tpsa,
    track, twiss, typeid, utest, utility, vector.

    > help "MADX"
    NAME
    MADX environment to emulate MAD-X workspace.

    SYNOPSIS
    local lhcb1 in MADX

    DESCRIPTION
    This module provide the function 'load' that read MADX sequence and optics
    files and load them in the MADX global variable. If it does not exist, it will
    create the global MADX variable as an object and load into it all elements,
    constants, and math functions compatible with MADX.

    RETURN VALUES
    The MADX global variable.

    EXAMPLES
    MADX:open()
    -- inline definition
    MADX:close()

    SEE ALSO
    element, object.


Complementary to the ``help`` function, the function ``show`` displays the type and value of variables, and if they have attributes, the list of their names in the lexicographic order:

.. code-block:: console
	
	> show "hello world!"
	:string: hello world!
	> show(MAD.option)
	:table: MAD.option
	colwidth           :number: 18
	hdrwidth           :number: 18
	intfmt             :string: % -10d
	madxenv            :boolean: false
	nocharge           :boolean: false
	numfmt             :string: % -.10g
	ptcmodel           :boolean: false
	strfmt             :string: % -25s


.. rubric:: Footnotes

.. [#f1] MAD-NG embeds the libraries  `FFTW <http://github.com/FFTW>`_ `NFFT <http://github.com/NFFT>`_ and `NLopt <http://github.com/stevengj/nlopt>`_ released under GNU (L)GPL too.
.. [#f5] MAD-NG has few thousands unit tests that do few millions checks, and it is constantly growing.
.. [#f2] Note that sending ``Crtl+D`` twice from MAD-NG invite will quit both MAD-NG and its parent Shell...
.. [#f3] A chunk is the unit of execution in Lua (see `Lua 5.2 <http://github.com/MethodicalAcceleratorDesign/MADdocs/blob/master/lua52-refman-madng.pdf>`_ §3.3.2).
.. [#f4] The online help is far incomplete and will be completed, updated and revised as the application evolves.
