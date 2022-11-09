.. index::
   Utility functions

***********************
Miscellaneous Functions
***********************

This chapter lists some useful functions from the module :mod:`MAD.utility` that are complementary to the standard library for manipulating files, strings, tables, and more.

Files Functions
===============

.. function:: openfile (filename_, mode_, extension_)

.. function:: filexists (filename)

.. function:: fileisnewer (filename1, filename2, timeattr_)

.. function:: filesplitname (filename)

.. object:: mockfile

Formating Functions
===================

.. function:: printf (str, ...)

.. function:: fprintf (file, str, ...)

.. function:: assertf (str, ...)

.. function:: errorf (str, ...)

Strings Functions
=================

.. function:: strinter (str, var, policy_)

.. function:: strtrim (str, ini_)

.. function:: strnum (str, ini_)

.. function:: strident (str, ini_)

.. function:: strquote (str, ini_)

.. function:: strbracket (str, ini_)

.. function:: strsplit (str, seps, ini_)

.. function:: strqsplit (str, seps, ini_)

.. function:: strqsplitall (str, seps, ini_, r_)

.. function:: is_identifier (str)

Tables Functions
================

.. function:: kpairs (tbl, n_)

.. function:: tblrep (val, n_, tbldst_)

.. function:: tblicpy (tblsrc, mtflag_, tbldst_)

.. function:: tblcpy (tblsrc, mtflag_, tbldst_)

.. function:: tbldeepcpy (tblsrc, mtflag_, xrefs_, tbldst_)

.. function:: tblcat (tblsrc1, tblsrc2, mtflag_, tbldst_)

.. function:: tblorder (tbl, key, n_)

Iterable Functions
==================

.. function:: rep (x, n_)

.. function:: clearidxs (a, i_, j_)

.. function:: setidxs (a, k_, i_, j_)

.. function:: bsearch (tbl, val, [cmp_,] low_, high_)

.. function:: lsearch (tbl, val, [cmp_,] low_, high_)

.. function:: monotonic (tbl, [strict_,] [cmp_,] low_, high_)

Mappable Functions
==================

.. function:: clearkeys (a, pred_)

.. function:: setkeys (a, k_, i_, j_)

.. function:: countkeys (a)

.. function:: keyscount (a, c_)

.. function:: val2keys (a)

Conversion Functions
====================

.. function:: log2num (log)

.. function:: num2log (num)

.. function:: num2str (num)

.. function:: int2str (int)

.. function:: str2str (str)

.. function:: str2cmp (str)

.. function:: tbl2str (tbl, sep_)

.. function:: str2tbl (str, match_, ini_)

.. function:: lst2tbl (lst, tbl_)

.. function:: tbl2lst (tbl, lst_)

Generic Functions
=================

.. function:: same (a, ...)

.. function:: copy (a, ...)

.. function:: tostring (a, ...)

.. function:: totable (a, ...)

.. function:: toboolean (a)

Special Functions
=================

.. function:: pause (msg_, val_)

.. function:: atexit (fun, uniq_)

.. function:: runonce (fun, ...)

.. function:: collectlocal (fun_, env_)


