MTables
=======
.. _ch.gen.mtbl:




The ``mtable`` object is the *root object* of the TFS tables that store information relative to tables.

The ``mtable`` module extends the :doc:`typeid <types>` module with the ``is_mtable`` function, which returns ``true`` if its argument is a ``mtable`` object, ``false`` otherwise.

Attributes
----------

The ``mtable`` object provides the following attributes:

**type**
	 A *string* specifying the type of the mtable (often) set to the name of the command that created it, like ``survey``, ``track`` or ``twiss``. (default: ``'user'``).

**title**
	 A *string* specifying the title of the mtable (often) set to the attribute ``title`` of the command that created it. (default: ``'no-title'``).

**origin**
	 A *string* specifying the origin of the mtable. (default: ``"MAD version os arch"``"}).

**date**
	 A *string* specifying the date of creation of the mtable. (default: ``day/month/year``"}).

**time**
	 A *string* specifying the time of creation of the mtable. (default: ``"hour:min:sec"``).

**refcol**
	 A *string* specifying the name of the reference column used to build the dictionary of the mtable, and to mangle values with counts. (default: ``nil``).

**header**
	 A *list* specifying the augmented attributes names (and their order) used by default for the header when writing the mtable to files. Augmented meaning that the *list* is concatenated to the *list* held by the parent mtable during initialization. 
	 (default: ``{'name', 'type', 'title', 'origin', 'date', 'time', 'refcol'}``).

**column**
	 A *list* specifying the augmented columns names (and their order) used by default for the columns when writing the mtable to files. Augmented meaning that the *list* is concatenated to the *list* held by the parent mtable during initialization. (default: ``nil``).

**novector**
	 A *logical* specifying to not convert (``novector == true``) columns containing only numbers to vectors during the insertion of the second row. The attribute can also be a *list* specifying the columns names to remove from the specialization. If the *list* is empty or ``novector ~= true``, all numeric columns will be converted to vectors, and support all methods and operations from the :doc:`linear algebra <linalg>` module. (default: ``nil``).

**owner**
	 A *logical* specifying if an *empty* mtable is a view with no data (``owner ~= true``), or a mtable holding data (``owner == true``). (default: ``nil``).

**reserve**
	 A *number* specifying an estimate of the maximum number of rows stored in the mtable. If the value is underestimated, the mtable will still expand on need. (default: ``8``).


**Warning**: the following private and read-only attributes are present in all mtables and should *never be used, set or changed*; breaking this rule would lead to an *undefined behavior*:

**__dat**
	 A *table* containing all the private data of mtables.

**__seq**
	 A *sequence* attached to the mtable by the ``survey`` and ``track`` commands and used by the methods receiving a *reference* to an element as argument. (default: ``nil``).

**__cycle**
	 A *reference* to the row registered with the ``:cycle`` method. (default: ``nil``).


Methods
-------

The ``mtable`` object provides the following methods:

**nrow**
	 A *method*	``()`` returning the *number* of rows in the mtable.

**ncol**
	 A *method*	``()`` returning the *number* of columns in the mtable.

**ngen**
	 A *method*	``()`` returning the *number* of columns generators in the mtable. The *number* of columns with data is given by ``:ncol() - :ngen()``.

**colname**
	 A *method*	``(idx)`` returning the *string* name of the ``idx``-th column in the mtable or ``nil``.

**colnames**
	 A *method*	``([lst])`` returning the *list* ``lst`` (default: ``{}``) filled with all the columns names of the mtable.

**index**
	 A *method*	``(idx)`` returning a positive index, or ``nil``. If ``idx`` is negative, it is reflected versus the size of the mtable, e.g. ``- 1`` becomes ``#self``, the index of the last row.

**name_of**
	 A *method*	``(idx, [ref])`` returning a *string* corresponding to the (mangled) *value* from the reference column of the row at the index ``idx``, or ``nil``. A row *value* appearing more than once in the reference column will be mangled with an absolute count, e.g. ``mq[3]``, or a relative count versus the reference row determined by ``:index_of(ref)``, e.g. ``mq{- 2}``.

**index_of**
	 A *method*	``(a, [ref], [dir])`` returning a *number* corresponding to the positive index of the row determined by the first argument or ``nil``. If ``a`` is a *number* (or a *string* representing a *number*), it is interpreted as the index of the row and returned as a second *number*. If ``a`` is a *string*, it is interpreted as the (mangled) *value* of the row in the reference column as returned by ``:name_of``. Finally, ``a`` can be a *reference* to an element to search for **if** the mtable has both, an attached sequence, and a column named ``'eidx'`` mapping the indexes of the elements to the attached sequence. [#f1]_ The argument ``ref`` (default: ``nil``) specifies the reference row determined by ``:index_of(ref)`` to use for relative indexes, for decoding mangled values with relative counts, or as the reference row to start searching from. The argument ``dir`` (default: ``1``) specifies the direction of the search with values ``1`` (forward), ``- 1`` (backward), or ``0`` (no direction), which correspond respectively to the rounding methods ``ceil``, ``floor`` and ``round`` from the lua math module.

**range_of**
	 A *method*	``([rng], [ref], [dir])`` returning three *number*\ s corresponding to the positive indexes *start* and *end* of the range and its direction *dir* (default: ``1``), or ``nil`` for an empty range. If ``rng`` is omitted, it returns ``1``, ``#self``, ``1``, or ``#self``, ``1``, ``- 1`` if ``dir`` is negative. If ``rng`` is a *number* or a *string* with no ``'/'`` separator, it is interpreted as *start* and *end*, both determined by ``:index_of``. If ``rng`` is a *string* containing the separator ``'/'``, it is split in two *string*\ s interpreted as *start* and *end*, both determined by ``:index_of``. If ``rng`` is a *list*, it will be interpreted as { *start*, *end*, ``[ref]``, ``[dir]`` }, both determined by ``:index_of``. The arguments ``ref`` and ``dir`` are forwarded to all invocations of ``:index_of`` with a higher precedence than ones in the *list* ``rng``, and a runtime error is raised if the method returns ``nil``, i.e. to disambiguate between a valid empty range and an invalid range.

**length_of**
	 A *method*	``([rng], [ntrn], [dir])`` returning a *number* specifying the length of the range optionally including ``ntrn`` extra turns (default: ``0``), and calculated from the indexes returned by ``:range_of([rng], nil, [dir])``.

**get**
	 A *method*	``(row, col, [cnt])`` returning the *value* stored in the mtable at the cell ``(row,col)``, or ``nil``. If ``row`` is a not a row index determined by ``:index(row)``, it is interpreted as a (mangled) *value* to search in the reference column, taking into account the count ``cnt`` (default: ``1``). If ``col`` is not a column index, it is interpreted as a column name.

**set**
	 A *method*	``(row, col, val, [cnt])`` returning the mtable itself after updating the cell ``(row,col)`` to the value ``val``, or raising an error if the cell does not exist. If ``row`` is a not a row index determined by ``:index(row)``, it is interpreted as a (mangled) *value* to search in the reference column, taking into account the count ``cnt`` (default: ``1``). If ``col`` is not a column index, it is interpreted as a column name.

**getcol**
	 A *method*	``(col)`` returning the column ``col``, or ``nil``. If ``col`` is not a column index, it is interpreted as a column name.

**setcol**
	 A *method*	``(col, val)`` returning the mtable itself after updating the column ``col`` with the values of ``val``, or raising an error if the column does not exist. If ``col`` is not a column index, it is interpreted as a column name. If the column is a generator, so must be ``val`` or an error will be raised. If the column is not a generator and ``val`` is a *callable* ``(ri)``, it will be invoked with the row index ``ri`` as its sole argument, using its returned value to update the column cell. Otherwise ``val`` must be an *iterable* or an error will be raised. If the column is already a specialized *vector*, the *iterable* must provide enough numbers to fill it entirely as ``nil`` is not a valid value.

**inscol**
	 A *method*	``([ref], col, val, [nvec])`` returning the mtable itself after inserting the column data ``val`` with the *string* name ``col`` at index ``ref`` (default: ``:ncol()+1``). If ``ref`` is not a column index, it is interpreted as a column name. If ``val`` is a *callable* ``(ri)``, it will be added as a column generator. Otherwise ``val`` must be an *iterable* or an error will be raised. The *iterable* will used to fill the new column that will be specialized to a *vector* if its first value is a *number* and ``nvec ~= true`` (default: ``nil``).

**addcol**
	 A *method*	``(col, val, [nvec])`` equivalent to ``:inscol(nil, col, val, [nvec])``.

**remcol**
	 A *method*	``(col)`` returning the mtable itself after removing the column ``col``, or raising an error if the column does not exist. If ``col`` is not a column index, it is interpreted as a column name.

**rencol**
	 A *method*	``(col, new)`` returning the mtable itself after renaming the column ``col`` to the *string* ``new``, or raising an error if the column does not exist. If ``col`` is not a column index, it is interpreted as a column name.

**getrow**
	 A *method*	``(row, [ref])`` returning the *mappable* (proxy) of the row determined by the method ``:index_of(row, [ref])``, or ``nil``.

**setrow**
	 A *method*	``(row, val, [ref])`` returning the mtable itself after updating the row at index determined by ``:index_of(row, [ref])`` using the values provided by the *mappable* ``val``, which can be a *list* iterated as pairs of (*index*, *value*) or a *set* iterated as pairs of (*key*, *value*) with *key* being the column names, or a combination of the two. An error is raised if the column does not exist.

**insrow**
	 A *method*	``(row, val, [ref])`` returning the mtable itself after inserting a new row at index determined by ``:index_of(row, [ref])`` and filled with the values provided by the *mappable* ``val``, which can be a *list* iterated as pairs of (*index*, *value*) or a *set* iterated as pairs of (*key*, *value*) with *key* being the column names or a combination of the two.

**addrow**
	 A *method*	``(val)`` equivalent to ``:insrow(#self+1, val)``.

**remrow**
	 A *method*	``(row, [ref])`` returning the mtable itself after removing the row determined by the method ``:index_of(row, [ref])``, or raising an error if the row does not exist.

**swprow**
	 A *method*	``(row1, row2, [ref1], [ref2])`` returning the mtable itself after swapping the content of the rows, both determined by the method ``:index_of(row, [ref])``, or raising an error if one of the row does not exist.

**clrrow**
	 A *method*	``(row, [ref])`` returning the mtable itself after clearing the row determined by the method ``:index_of(row, [ref])``, or raising an error if the row does not exist; where clearing the row means to set *vector* value to ``0`` and ``nil`` otherwise.

**clear**
	 A *method*	``()`` returning the mtable itself after clearing all the rows, i.e. ``#self == 0``, with an opportunity for new columns specialization.

**iter**
	 A *method*	``([rng], [ntrn], [dir])`` returning an iterator over the mtable rows. The optional range is determined by ``:range_of([rng], [dir])``, optionally including ``ntrn`` turns (default: ``0``). The optional direction ``dir`` specifies the forward ``1`` or the backward ``- 1`` direction of the iterator. If ``rng`` is not provided and the mtable is cycled, the *start* and *end* indexes are determined by ``:index_of(self.__cycle)``. When used with a generic :literal:`for` loop, the iterator returns at each rows the index and the row *mappable* (proxy).

**foreach**
	 A *method*	``(act, [rng], [sel], [not])`` returning the mtable itself after applying the action ``act`` on the selected rows. If ``act`` is a *set* representing the arguments in the packed form, the missing arguments will be extracted from the attributes ``action``, ``range``, ``select`` and ``default``. The action ``act`` must be a *callable* ``(row, idx)``  applied to a row passed as first argument and its index as second argument. The optional range is used to generate the loop iterator ``:iter([rng])``. The optional selector ``sel`` is a *callable* ``(row, idx)`` predicate selecting eligible rows for the action from the row itself passed as first argument and its index as second argument. The selector ``sel`` can be specified in other ways, see :doc:`row selections <numrange>` for details. The optional *logical* ``not`` (default: ``false``) indicates how to interpret default selection, as *all* or *none*, depending on the semantic of the action. [#f2]_ method needs remove all rows if no selector is provided.}

**select**
	 A *method*	``([rng], [sel], [not])`` returning the mtable itself after selecting the rows using ``:foreach(sel_act, [rng], [sel], [not])``. By default mtable have all their rows deselected, the selection being stored as *boolean* in the column at index ``0`` and named ``is_selected``.

**deselect**
	 A *method*	``([rng], [sel], [not])`` returning the mtable itself after deselecting the rows using ``:foreach(desel_act, [rng], [sel], [not])``. By default mtable have all their rows deselected, the selection being stored as *boolean* in the column at index ``0`` and named ``is_selected``.

**filter**
	 A *method*	``([rng], [sel], [not])`` returning a *list* containing the positive indexes of the rows determined by ``:foreach(filt_act, [rng], [sel], [not])``, and its size.

**insert**
	 A *method*	``(row, [rng], [sel])`` returning the mtable itself after inserting the rows in the *list* ``row`` at the indexes determined by ``:filter([rng], [sel], true)``. If the arguments are passed in the packed form, the extra attribute ``rows`` will be used as a replacement for the argument ``row``, and if the attribute ``where="after"`` is defined then the rows will be inserted after the selected indexes. The insertion scheme depends on the number :math:`R` of rows in the *list* ``row`` versus the number :math:`S` of rows selected by ``:filter``; :math:`1\times 1` (one row inserted at one index), :math:`R\times 1` (:math:`R` rows inserted at one index), :math:`1\times S` (one row inserted at :math:`S` indexes) and :math:`R\times S` (:math:`R` rows inserted at :math:`S` indexes). Hence, the insertion schemes insert respectively :math:`1`, :math:`R`, :math:`S`, and :math:`\min(R, S)` rows.

**remove**
	 A *method*	``([rng], [sel])`` returning the mtable itself after removing the rows determined by ``:filter([rng], [sel], true)``.

**sort**
	 A *method*	``(cmp, [rng], [sel])`` returning the mtable itself after sorting the rows at the indexes determined by ``:filter([rng], [sel], true)`` using the ordering *callable* ``cmp(row1, row2)``. The arguments ``row1`` and ``row2`` are *mappable* (proxies) referring to the current rows being compared and providing access to the columns values for the comparison. [#f3]_ The argument ``cmp`` can be specified in a compact ordering form as a *string* that will be converted to an ordering *callable* by the function ``str2cmp`` from the :doc:`utility <numrange>` module. For example, the *string* "-y,x" will be converted by the method to the following *lambda* :literal:`\r1,r2 -> r1.y > r2.y or r1.y == r2.y and r1.x < r2.x`, where ``y`` and ``x`` are the columns used to sort the mtable in descending (``-``) and ascending (``+``) order respectively. The compact ordering form is not limited in the number of columns and avoids making mistakes in the comparison logic when many columns are involved.

**cycle**
	 A *method*	``(a)`` returning the mtable itself after checking that ``a`` is a valid reference using ``:index_of(a)``, and storing it in the ``__cycle`` attribute, itself erased by the methods editing the mtable like ``:insert``, ``:remove`` or ``:sort``.

**copy**
	 A *method*	``([name], [owner])`` returning a new mtable from a copy of ``self``, with the optional ``name`` and the optional attribute ``owner`` set. If the mtable is a view, so will be the copy unless ``owner == true``.

**is_view**
	 A *method*	``()`` returning ``true`` if the mtable is a view over another mtable data, ``false`` otherwise.

**set_readonly**
	 Set the mtable as read-only, including the columns and the rows proxies.

**read**
	 A *method*	``([filname])`` returning a new instance of ``self`` filled with the data read from the file determined by ``openfile(filename, 'r', {'.tfs','.txt','.dat'})`` from the :doc:`utility <miscfuns>` module. This method can read columns containing the data types *nil*, *boolean*, *number*, *complex number*, (numerical) *range*, and (quoted) *string*. The header can also contain tables saved as *string* and decoded with *function* ``str2tbl`` from the :doc:`utility <miscfuns>` module.

**write**
	 A *method*	``([filname], [clst], [hlst], [rsel])`` returning the mtable itself after writing its content to the file determined by ``openfile(filename, 'w', {'.tfs', '.txt', '.dat'})`` from the :doc:`utility <miscfuns>` module. The columns to write and their order is determined by ``clst`` or ``self.column`` (default: ``nil`` :math:`\equiv` all columns). The attributes to write in the header and their order is determined by ``hlst`` or ``self.header``. The *logical* ``rsel`` indicates to save all rows or only rows selected by the ``:select`` method (``rsel == true``). This method can write columns containing the data types *nil*, *boolean*, *number*, *complex number*, (numerical) *range*, and (quoted) *string*. The header can also contain tables saved as *string* and encoded with *function* ``tbl2str`` from the :doc:`utility <miscfuns>` module.

**print**
	 A *method*	``([clst], [hlst], [rsel])`` equivalent to ``:write(nil, [clst], [hlst], [rsel])``.

**save_sel**
	 A *method*	``([sel])`` saving the rows selection to the optional *iterable* ``sel`` (default: ``{}``) and return it.

**restore_sel**
	 A *method*	``(sel)`` restoring the rows selection from the *iterable* ``sel``. The indexes of ``sel`` must match the indexes of the rows in the mtable.

**make_dict**
	 A *method*	``([col])`` returning the mtable itself after building the rows dictionnary from the values of the reference column determined by ``col`` (default: ``refcol``) for fast row access. If ``col`` is not a column index, it is interpreted as a column name except for the special name ``'none'`` that disables the rows dictionnary and reset ``refcol`` to ``nil``.

**check_mtbl**
	 A *method*	``()`` checking the integrity of the mtable and its dictionary (if any), for debugging purpose only.


Metamethods
-----------

The ``mtable`` object provides the following metamethods:


**__len**
	 A *metamethod*	``()`` called by the length operator ``#`` to return the number of rows in the mtable.

**__add**
	 A *metamethod*	``(val)`` called by the plus operator ``+`` returning the mtable itself after appending the row ``val`` at its end, similiar to the ``:addrow`` method.

**__index**
	 A *metamethod*	``(key)`` called by the indexing operator ``[key]`` to return the *value* of an attribute determined by *key*. The *key* is interpreted differently depending on its type with the following precedence:

		#. A *number* is interpreted as a row index and returns an *iterable* on the row (proxy) or ``nil``.
		#. Other *key* types are interpreted as *object* attributes subject to object model lookup.
		#. If the *value* associated with *key* is ``nil``, then *key* is interpreted as a column name and returns the column if it exists, otherwise...
		#. If *key* is not a column name, then *key* is interpreted as a value in the reference column and returns either an *iterable* on the  row (proxy) determined by this value or an *iterable* on the rows (proxies) holding this non-unique value. 
		#. Otherwise returns ``nil``.

**__newindex**
	 A *metamethod*	``(key, val)`` called by the assignment operator ``[key]=val`` to create new attributes for the pairs (*key*, *value*). If *key* is a *number* or a value specifying a row in the reference column or a *string* specifying a column name, the following error is raised:

.. code-block::
	
	"invalid mtable write access (use 'set' methods)"


**__init**
	 A *metamethod*	``()`` called by the constructor to build the mtable from the column names stored in its *list* part and some attributes, like ``owner``, ``reserve`` and ``novector``.

**__copy**
	 A *metamethod*	``()`` similar to the ``copy``.




**__mtbl**
	 A unique private *reference* that characterizes mtables.


MTables creation
----------------
.. _sec.tbl.create:


Any column name in the *list* that is enclosed by braces is designated as the refererence column for the dictionnary that provides fast row indexing, and the attribute ``refcol`` is set accordingly.

Some attributes are considered during the creation by the ``__init``, like ``owner``, ``reserve`` and ``novector``, and some others are initialized with defined values like ``type``, ``title``, ``origin``, ``date``, ``time``, and ``refcol``. The attributes ``header`` and ``column`` are concatenated with the the parent ones to build incrementing *list* of attributes names and columns names used by default when writing the mtable to files, and these lists are not provided as arguments.



.. code-block::
	
	local mtable in MAD
	local tbl = mtable 'mytable' {
	
	   {'name'}, 'x', 'y' } -- column 'name' is the refcol
	  + { 'p11', 1.1, 1.2 }
	  + { 'p12', 2.1, 2.2 }
	  + { 'p13', 2.1, 3.2 }
	  + { 'p11', 3.1, 4.2 }
	print(tbl.name, tbl.refcol, tbl:getcol'name')
	-- display: mytable  name   mtable reference column: 0x010154cd10

**Pitfall:** When a column is named ``'name'``, it must be explicitly accessed, e.g. with the ``:getcol`` method, as the indexing operator ``[]`` gives the precedence to object's attributes and methods. Hence, ``tbl.name`` returns the table name ``'mytable'``, not the column ``'name'``.

Rows selections
---------------
.. _sec.tbl.rowsel:

The row selection in mtable use predicates in combination with iterators. The mtable iterator manages the range of rows where to apply the selection, while the predicate says if a row in this range is illegible for the selection. In order to ease the use of methods based on the ``:foreach`` method, the selector predicate ``sel`` can be built from different types of information provided in a *set* with the following attributes:

**selected**
	 A *boolean* compared to the rows selection stored in column ``'is_selected'``.

**pattern**
	 A *string* interpreted as a pattern to match the *string* in the reference column, which must exist, using ``string.match`` from the standard library, see `Lua 5.2 <http://github.com/MethodicalAcceleratorDesign/MADdocs/blob/master/lua52-refman-madng.pdf>`_ §6.4 for details. If the reference column does not exist, it can be built using the method.

**list**
	 An *iterable* interpreted as a *list* used to build a *set* and select the rows by their name, i.e. the built predicate will use ``tbl[row.name]`` as a *logical*, meaning that column ``name`` must exists. An alternate column name can be provided through the key ``colname``, i.e. used as ``tbl[row[colname]]``. If the *iterable* is a single item, e.g. a *string*, it will be converted first to a *list*.

**table**
	 A *mappable* interpreted as a *set* used to select the rows by their name, i.e. the built predicate will use ``tbl[row.name]`` as a *logical*, meaning that column ``name`` must exists. If the *mappable* contains a *list* or is a single item, it will be converted first to a *list* and its *set* part will be discarded.

**kind**
	 An *iterable* interpreted as a *list* used to build a *set* and select the rows by their kind, i.e. the built predicate will use ``tbl[row.kind]`` as a *logical*, meaning that column ``kind`` must exists. If the *iterable* is a single item, e.g. a *string*, it will be converted first to a *list*. This case is equivalent to ``list`` with ``colname='kind'``.

**select**
	 A *callable* interpreted as the selector itself, which allows to build any kind of predicate or to complete the restrictions already built above.

All these attributes are used in the aforementioned order to incrementally build predicates that are combined with logical conjunctions, i.e. ``and``'ed, to give the final predicate used by the ``:foreach`` method. If only one of these attributes is needed, it is possible to pass it directly in ``sel``, not as an attribute in a *set*, and its type will be used to determine the kind of predicate to build. For example, ``tbl:foreach(act, "\POW MB")`` is equivalent to ``tbl:foreach{action=act, pattern="\POW MB"}``.

Indexes, names and counts
-------------------------

Indexing a mtable triggers a complex look up mechanism where the arguments will be interpreted in various ways as described in the metamethod ``__index``. A *number* will be interpreted as a relative row index in the list of rows, and a negative index will be considered as relative to the end of the mtable, i.e. ``- 1`` is the last row. Non-*number* will be interpreted first as an object key (can be anything), looking for mtable methods or attributes; then as a column name or as a row *value* in the reference column if nothing was found.

If a row exists but its *value* is not unique in the reference column, an *iterable* is returned. An *iterable* supports the length ``#`` operator to retrieve the number of rows with the same *value*, the indexing operator ``[]`` waiting for a count :math:`n` to retrieve the :math:`n`-th row from the start with that *value*, and the iterator ``ipairs`` to use with generic ``for`` loops.



**Note:** Compared to the sequence, the indexing operator ``[]`` and the method ``:index_of`` of the mtable always interprets a *number* as a (relative) row index. To find a row from a :math:`s`-position [m] in the mtable if the column exists, use the functions ``lsearch`` or ``bsearch`` (if they are monotonic) from the :doc:`utility <miscfuns>` module.



.. code-block::
	
	local mtable in MAD
	local tbl = mtable { {'name'}, 'x', 'y' } -- column 'name' is the refcol
	                   + { 'p11', 1.1, 1.2 }
	                   + { 'p12', 2.1, 2.2 }
	                   + { 'p13', 2.1, 3.2 }
	                   + { 'p11', 3.1, 4.2 }
	print(tbl[ 1].y) -- display: 1.2
	print(tbl[-1].y) -- display: 4.2
	
	print(#tbl.p11, tbl.p12.y, tbl.p11[2].y)            -- display: 2 2.2 4.2
	for _,r in ipairs(tbl.p11) do io.write(r.x," ") end -- display: 1.1 3.1
	for _,v in ipairs(tbl.p12) do io.write(v,  " ") end -- display: 'p12' 2.1 2.2
	
	-- print name of point with name p11 in absolute and relative to p13.
	print(tbl:name_of(4))       -- display: p11[2]  (2nd p11 from start)
	print(tbl:name_of(1, -2))   -- display: p11{-1} (1st p11 before p13)



Iterators and ranges
--------------------

Ranging a mtable triggers a complex look up mechanism where the arguments will be interpreted in various ways as described in the method ``:range_of``, itself based on the methods ``:index_of`` and ``:index``. The number of rows selected by a mtable range can be computed by the ``:length_of`` method, which accepts an extra *number* of turns to consider in the calculation.

The mtable iterators are created by the method ``:iter``, based on the method ``:range_of`` as mentioned in its description and includes an extra *number* of turns as for the method ``:length_of``, and a direction ``1`` (forward) or ``- 1`` (backward) for the iteration.

The method ``:foreach`` uses the iterator returned by ``:iter`` with a range as its sole argument to loop over the rows where to apply the predicate before executing the action. The methods ``:select``, ``:deselect``, ``:filter``, ``:insert``, and ``:remove`` are all based directly or indirectly on the ``:foreach`` method. Hence, to iterate backward over a mtable range, these methods have to use either its *list* form or a numerical range. For example the invocation :literal:`tbl:foreach(\r -> print(r.name), {-2, 2, nil, -1})` will iterate backward over the entire mtable excluding the first and last rows, equivalently to the invocation :literal:`tbl:foreach(\r -> print(r.name), -2..2..-1)`.

The following example shows how to access to the rows with the ``:foreach`` method:

.. code-block::
	
	local mtable in MAD
	local tbl = mtable { {'name'}, 'x', 'y' }
	                   + { 'p11', 1.1, 1.2 }
	                   + { 'p12', 2.1, 2.2 }
	                   + { 'p13', 2.1, 3.2 }
	                   + { 'p11', 3.1, 4.2 }
	
	local act = \r -> print(r.name, r.y)
	tbl:foreach(act, -2..2..-1)
	-- display:  p13   3.2
	!            p12   2.2
	tbl:foreach(act, "p11[1]/p11[2]")
	-- display:  p11   1.2
	!            p12   2.2
	!            p13   3.2
	!            p11   4.2
	tbl:foreach{action=act, range="p11[1]/p13"}
	-- display:  p11   1.2
	!            p12   2.2
	!            p13   3.2
	tbl:foreach{action=act, pattern="[^1](*\$*)"}
	-- display:  p12   2.2
	!            p13   3.2
	local act = \r -> print(r.name, r.y, r.is_selected)
	tbl:select{pattern="p.1"}:foreach{action=act, range="1/-1"}
	-- display:  p11   1.2   true
	!            p12   2.2   nil
	!            p13   3.2   nil
	!            p11   4.2   true


Examples
--------

Creating a MTable
"""""""""""""""""

The following example shows how the ``track`` command, i.e. ``self`` hereafter, creates its MTable:

.. code-block::
	
	local header = { -- extra attributes to save in track headers
	  'direction', 'observe', 'implicit', 'misalign', 'deltap', 'lost' }
	
	local function make_mtable (self, (*range*), nosave)
	  local title, dir, observe, implicit, misalign, deltap, savemap in self
	  local sequ, nrow = self.sequence, nosave and 0 or 16
	
	  return mtable(sequ.name, { -- keep column order!
	    (*type*)='track', title=title, header=header,
	    direction=dir, observe=observe, implicit=implicit, misalign=misalign,
	    deltap=deltap, lost=0, (*range*)=(*range*), reserve=nrow, __seq=sequ,
	    {'name'}, 'kind', 's', 'l', 'id', 'x', 'px', 'y', 'py', 't', 'pt',
	    'slc', 'turn', 'tdir', 'eidx', 'status', savemap and '__map' or nil })
	end


Extending a MTable
""""""""""""""""""

The following example shows how to extend the MTable created by a ``twiss`` command with the elements tilt, angle and integrated strengths from the attached sequence:

.. code-block::
	
	-- The prelude creating the sequence seq is omitted.
	local tws = twiss { sequence=seq, method=4, cofind=true }
	
	local is_integer in MAD.typeid
	tws:addcol('angle', \ri => -- add angle column
	      local idx = tws[ri].eidx
	      return is_integer(idx) and tws.__seq[idx].angle or 0 end)
	   :addcol('tilt', \ri => -- add tilt column
	      local idx = tws[ri].eidx
	      return is_integer(idx) and tws.__seq[idx].tilt or 0 end)
	
	for i=1,6 do -- add k(*\IT{i}*)l and k(*\IT{i}*)sl columns
	tws:addcol('k'..i-1..'l', \ri =>
	      local idx = tws[ri].eidx
	      if not is_integer(idx) then return 0 end -- implicit drift
	      local elm = tws.__seq[idx]
	      return (elm['k'..i-1] or 0)*elm.l + ((elm.knl or {})[i] or 0)
	    end)
	   :addcol('k'..i-1..'sl', \ri =>
	      local idx = tws[ri].eidx
	      if not is_integer(idx) then return 0 end -- implicit drift
	      local elm = tws.__seq[idx]
	      return (elm['k'..i-1..'s'] or 0)*elm.l + ((elm.ksl or {})[i] or 0)
	    end)
	end
	
	local cols = {'name', 'kind', 's', 'l', 'angle', 'tilt',
	    'x', 'px', 'y', 'py', 't', 'pt',
	    'beta11', 'beta22', 'alfa11', 'alfa22', 'mu1', 'mu2', 'dx', 'ddx',
	    'k1l', 'k2l', 'k3l', 'k4l', 'k1sl', 'k2sl', 'k3sl', 'k4sl'}
	
	tws:write("twiss", cols) -- write header and columns to file twiss.tfs

Hopefully, the :doc:`physics <gphys>` module provides the *function* ``melmcol(mtbl, cols)`` to achieve the same task easily:

.. code-block::
	
	-- The prelude creating the sequence seq is omitted.
	local tws = twiss { sequence=seq, method=4, cofind=true }
	
	-- Add element properties as columns
	local melmcol in MAD.gphys
	local melmcol(tws, {'angle', 'tilt', 'k1l' , 'k2l' , 'k3l' , 'k4l',
	                                     'k1sl', 'k2sl', 'k3sl', 'k4sl'})
	
	-- write TFS table
	tws:write("twiss", {
	    'name', 'kind', 's', 'l', 'angle', 'tilt',
	    'x', 'px', 'y', 'py', 't', 'pt',
	    'beta11', 'beta22', 'alfa11', 'alfa22', 'mu1', 'mu2', 'dx', 'ddx',
	    'k1l', 'k2l', 'k3l', 'k4l', 'k1sl', 'k2sl', 'k3sl', 'k4sl'})


.. rubric:: Footnotes

.. [#f1] These information are usually provided by the command creating the ``mtable``, like ``survey`` and ``track``.
.. [#f2] For example, the ``:remove`` method needs ``not=true`` to *not* remove all rows if no selector is provided.
.. [#f3] A *mappable* supports the length operator ``#``, the indexing operator ``[]``, and generic ``for`` loops with ``pairs``.
.. [#f4] An *iterable* supports the length operator ``#``, the indexing operator ``[]``, and generic ``for`` loops with ``ipairs``.
