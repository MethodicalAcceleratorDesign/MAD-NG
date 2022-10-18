.. index::
   Linear algebra
   Vector and matrix

**************
Linear Algebra
**************

This chapter describes the real :type:`matrix`, complex :type:`cmatrix` and integer :type:`imatrix` objects as supported by MAD-NG. The module for `Vector <https://en.wikipedia.org/wiki/Vector_(mathematics_and_physics)>`_ and `Matrix <https://en.wikipedia.org/wiki/Matrix_(mathematics)>`_ is not exposed, only the contructors are visible from the :mod:`MAD` environment and thus, matrices are handled directly by their methods or by the generic functions of the same name from the module :mod:`MAD.gmath`. The :type:`imatrix`, i.e. matrix of integers, are mainly used for indexing other types of matrix and therefore supports only a limited subset of the features. Column and row vectors are shortcuts for :math:`[n\times 1]` and :math:`[1\times n]` matrices respectively. Note that :type:`matrix`, :type:`cmatrix` and :type:`imatrix` are all defined as C structures containing their elements in `row-major order <https://en.wikipedia.org/wiki/Row-_and_column-major_order>`_ for direct compliance with the C API. 

Types promotion
===============

The matrix operations may involve other data types like real and complex numbers leading to many combinations of types. In order to simplify the descriptions, the generic names :var:`num`, :var:`cpx` and :var:`idx` (indexes) are used for real, complex and integer numbers respectively, :var:`vec`, :var:`cvec` and :var:`ivec` for real, complex and integer vectors respectively, and :var:`mat`, :var:`cmat` and :var:`imat` for real, complex and integer matrices respectively. For example, the sum of a complex number :var:`cpx` and a real matrix :var:`mat` gives a complex matrix :var:`cmat`. The case of :var:`idx` means that a :type:`number` will be interpreted as an index and automatically rounded if it does not hold an integer value. The following table summarizes all valid combinations of types for binary operations involving at least one matrix type:

=================  ==================  ===============
Left Operand Type  Right Operand Type  Result Type
=================  ==================  ===============
:type:`number`     :type:`imatrix`     :type:`imatrix`
:type:`imatrix`    :type:`number`      :type:`imatrix`  
:type:`imatrix`    :type:`imatrix`     :type:`imatrix`
                                       
:type:`number`     :type:`matrix`      :type:`matrix` 
:type:`matrix`     :type:`number`      :type:`matrix`  
:type:`matrix`     :type:`matrix`      :type:`matrix`  
                                       
:type:`number`     :type:`cmatrix`     :type:`cmatrix`
:type:`complex`    :type:`matrix`      :type:`cmatrix` 
:type:`complex`    :type:`cmatrix`     :type:`cmatrix`
:type:`matrix`     :type:`complex`     :type:`cmatrix`
:type:`matrix`     :type:`cmatrix`     :type:`cmatrix`
:type:`cmatrix`    :type:`number`      :type:`cmatrix`  
:type:`cmatrix`    :type:`complex`     :type:`cmatrix`
:type:`cmatrix`    :type:`matrix`      :type:`cmatrix`  
:type:`cmatrix`    :type:`cmatrix`     :type:`cmatrix`
=================  ==================  ===============

Constructors
============

The constructors for vectors and matrices are directly available from the :mod:`MAD` environment. Note that real, complex or integer matrix with zero size are not allowed, i.e. the smallest allowed matrix has sizes of :math:`[1\times 1]`.

.. function::  vector(nrow)
              cvector(nrow)
              ivector(nrow)

   Return a real, complex or integer column vector (i.e. a matrix of size :math:`[n_{\text{row}}\times 1]`) filled with zeros. If :var:`nrow` is a table, it is equivalent to :expr:`vector(#nrow):fill(nrow)`. 

.. function::  matrix(nrow, ncol_)
              cmatrix(nrow, ncol_)
              imatrix(nrow, ncol_)

   Return a real, complex or integer matrix of size :math:`[n_{\text{row}}\times n_{\text{col}}]` filled with zeros. If :var:`nrow` is a table, it is equivalent to :expr:`matrix(#nrow, #nrow[1] or 1):fill(nrow)`, and ignoring :var:`ncol`. Default: :expr:`ncol_ = rnow`. 

.. function:: linspace([start_,] stop, size_)

   Return a real or complex column vector of length :var:`size` filled with values equally spaced between :var:`start` and :var:`stop` on a linear scale. Note that numerical :type:`range` can generate the same *real* sequence of values in a more compact form. Default: :expr:`start_ = 0`, :expr:`size_ = 100`.

.. function:: logspace([start_,] stop, size_)

   Return a real or complex column vector of length :var:`size` filled with values equally spaced between :var:`start` and :var:`stop` on a logarithmic scale. Note that numerical :type:`logrange` can generate the same *real* sequence of values in a more compact form. Default: :expr:`start_ = 1`, :expr:`size_ = 100`.

Attributes
==========

.. constant:: mat.nrow

   The number of rows of the real, complex or integer matrix :var:`mat`.

.. constant:: mat.ncol

   The number of columns of the real, complex or integer matrix :var:`mat`.

Functions
=========

.. function:: is_vector (a)
              is_cvector (a)
              is_ivector (a)

   Return :const:`true` if :var:`a` is respectively a real, complex or integer matrix of size :math:`[n_{\text{row}}\times 1]` or :math:`[1\times n_{\text{row}}]`, :const:`false` otherwise. These functions are only available from the module :mod:`MAD.typeid`.

.. function:: isa_vector (a)

   Return :const:`true` if :var:`a` is a real or complex vector (i.e. is-a vector), :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

.. function:: isy_vector (a)

   Return :const:`true` if :var:`a` is a real, complex or integer vector (i.e. is-any vector), :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

.. function:: is_matrix (a)
              is_cmatrix (a)
              is_imatrix (a)

   Return :const:`true` if :var:`a` is respectively a real, complex or integer matrix, :const:`false` otherwise. These functions are only available from the module :mod:`MAD.typeid`.

.. function:: isa_matrix (a)

   Return :const:`true` if :var:`a` is a real or complex matrix (i.e. is-a matrix), :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

.. function:: isy_matrix (a)

   Return :const:`true` if :var:`a` is a real, complex or integer matrix (i.e. is-any matrix), :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

Special Constructors
--------------------

.. function:: mat:vec ()

   Return a vector of the same type as :var:`mat` filled with the values of the elements of the `vectorized <https://en.wikipedia.org/wiki/Vectorization_(mathematics)>`_ real, complex or integer matrix :var:`mat` equivalent to :expr:`mat:t():reshape(#mat,1)`.

.. function:: mat:vech ()

   Return a vector of the same type as :var:`mat` filled with the values of the elements of the `half vectorized <https://en.wikipedia.org/wiki/Vectorization_(mathematics)#Half-vectorization>`_ real, complex or integer *symmetric* matrix :var:`mat`. The symmetric property can be pre-checked bu the user with :func:`mat:is_symm()`.

.. function:: mat:diag (k_)

   Return a zero square matrix of the same type as :var:`mat` where :expr:`n = #mat + abs(k)`, and with its :math:`k`-th diagonal :expr:`-n <= k <= n` filled with the values of the elements of :var:`mat`, i.e. interpreting the matrix as the :math:`k`-th diagonal vector. Default: :expr:`k_ = 0`.

Sizes and Indexing
------------------

.. function:: mat:sizes ()

   Return the number of rows :var:`nrow` and columns :var:`ncol` of the real, complex or integer matrix :var:`mat`. Note that :expr:`#mat` returns the full size :expr:`nrow * ncol` of the matrix.

.. function:: mat:tsizes ()

   Return the number of columns :var:`ncol` and rows :var:`nrow` (i.e. transposed sizes) of the real, complex or integer matrix :var:`mat` equivalent to :expr:`swap(mat:sizes())`.

.. function:: mat:getij (ij_, ri_, rj_)

   Return two :type:`ivector` or :var:`ri` and :var:`rj` containing the indexes :expr:`(i,j)` extracted from the :type:`iterable` :var:`ij` for the real, complex or integer matrix :var:`mat`. If :var:`ij` is a number, the two returned items are also numbers. This method is the reverse method of :func:`mat:getidx()` to convert 1D indexes into 2D indexes for the given matrix sizes. Default: :expr:`ij_ = 1..#mat`.

.. function:: mat:getidx (ir_, jc_, rij_)

   Return an :type:`ivector` or :var:`rij` containing :expr:`#ir * #jc` vector indexes in row-major order given by the :type:`iterable` :var:`ir` and :var:`jc` of the real, complex or integer matrix :var:`mat`, followed by :var:`ir` and :var:`jc` potentially set from defaults. If both :var:`ir` and :var:`jc` are numbers, it returns a single number. This method is the reverse method of :func:`mat:getij()` to convert 2D indexes into 1D indexes for the given matrix sizes. Default: :expr:`ir_ = 1..nrow`, :expr:`jc_ = 1..ncol`.

.. function:: mat:getdidx (k_)

   Return an :type:`iterable` describing the indexes of the :math:`k`-th diagonal of the real, complex or integer matrix :var:`mat` where :expr:`-nrow <= k <= ncol`. This method is useful to build the 1D indexes of the :math:`k`-th diagonal for the given matrix sizes.

Getters and Setters
-------------------

.. function:: mat:get (i, j)

   Return the value of the element at the indexes :expr:`(i,j)` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= i <= nrow` and :expr:`1 <= j <= ncol`, :const:`nil` otherwise.

.. function:: mat:set (i, j, v)

   Assign the value :var:`v` to the element at the indexes :expr:`(i,j)` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= i <= nrow` and :expr:`1 <= j <= ncol` and return the matrix, otherwise raise an *"index out of bounds"* error.

.. function:: mat:geti (n)

   Return the value of the element at the vector index :var:`n` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= n <= #mat`, i.e. interpreting the matrix as a vector, :const:`nil` otherwise.

.. function:: mat:seti (n, v)

   Assign the value :var:`v` to the element at the vector index :var:`n` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= n <= #mat` and return the matrix, i.e. interpreting the matrix as a vector, otherwise raise an *"index out of bounds"* error.

.. function:: mat:getvec (ij, r_)

   Return a column vector or :var:`r` containing the values of the elements at the vector indexes given by the :type:`iterable` :var:`ij` of the real, complex or integer matrix :var:`mat`, i.e. interpreting the matrix as a vector.

.. function:: mat:setvec (ij, a, p_, s_)

   Return the real, complex or integer matrix :var:`mat` after filling the elements at the vector indexes given by the :type:`iterable` :var:`ij`, i.e. interpreting the matrix as a vector, with the values given by :var:`a` depending of its kind:

   - if :var:`a` is a :type:`scalar`, it is will be used repetitively.

   - if :var:`a` is an :type:`iterable` then the matrix will be filled with values from :var:`a[n]` for :expr:`1 <= n <= #a` and recycled repetitively if :expr:`#a < #ij`.

   - if :var:`a` is a :type:`callable`, then :var:`a` is considered as a *stateless iterator*, and the matrix will be filled with the values :var:`v` returned by iterating :expr:`s, v = a(p, s)`.

.. function:: mat:swpvec (ij, ij2)

   Return the real, complex or integer matrix :var:`mat` after swapping the values of the elements at the vector indexes given by the :type:`iterable` :var:`ij` and :var:`ij2`, i.e. interpreting the matrix as a vector.

.. function:: mat:remvec (ij)

   Return the real, complex or integer matrix :var:`mat` after removing the elements at the vector indexes given by the :type:`iterable` :var:`ij`, i.e. interpreting the matrix as a shrinking vector, and reshaped as a *column vector* of size :expr:`#mat - #ij`.

.. function:: mat:insvec (ij, a)

   Return the real, complex or integer matrix :var:`mat` after inserting the elements at the vector indexes given by the :type:`iterable` :var:`ij`, i.e. interpreting the matrix as a vector, with the values given by :var:`a` depending of its kind:
   
   - if :var:`a` is a :type:`scalar`, it is will be used repetitively.

   - if :var:`a` is an :type:`iterable` then the matrix will be filled with values from :var:`a[n]` for :expr:`1 <= n <= #a`.
   
   The elements after the inserted indexes are shifted toward the end of the matrix in row-major order and discarded if they go beyond the last index.

.. function:: mat:getsub (ir_, jc_, r_)

   Return a :math:`[` :expr:`#ir` :math:`\times` :expr:`#jc` :math:`]` matrix or :var:`r` containing the values of the elements at the indexes given by the :type:`iterable` :var:`ir` and :var:`jc` of the real, complex or integer matrix :var:`mat`. If :expr:`ir == nil`, :expr:`jc ~= nil` and :var:`r` is a 1D :type:`iterable`, then the latter is filled using column-major indexes. Default: as :func:`mat:getidx()`.

.. function:: mat:setsub (ir_, jc_, a, p_, s_)

   Return the real, complex or integer matrix :var:`mat` after filling the elements at the indexes given by the :type:`iterable` :var:`ir` and :var:`jc` with the values given by :var:`a` depending of its kind:

   - if :var:`a` is a :type:`scalar`, it is will be used repetitively.

   - if :var:`a` is an :type:`iterable` then the rows and columns will be filled with values from :var:`a[n]` for :expr:`1 <= n <= #a` and recycled repetitively if :expr:`#a < #ir * #ic`.

   - if :var:`a` is a :type:`callable`, then :var:`a` is considered as a *stateless iterator*, and the columns will be filled with the values :var:`v` returned by iterating :expr:`s, v = a(p, s)`.

   If :expr:`ir == nil`, :expr:`jc ~= nil` and :var:`a` is a 1D :type:`iterable`, then the latter is used to filled the matrix in the column-major order. Default: as :func:`mat:getidx()`.

.. function:: mat:swpsub (ir_, jc_, ir2_, jc2_)

   Return the real, complex or integer matrix :var:`mat` after swapping the elements at indexes given by the iterable :type:`iterable` :var:`ir` and :var:`jc` with the elements at indexes given by :type:`iterable` :var:`ir2` and :var:`jc2`. Default: as :func:`mat:getidx()`.

.. function:: mat:remsub (ir_, jc_)

   Return the real, complex or integer matrix :var:`mat` after removing the rows and columns at the indexes given by the :type:`iterable` :var:`ir` and :var:`jc` and reshaping the matrix accordingly. Default: as :func:`mat:getidx()`.
  
.. function:: mat:inssub (ir_, jc_, a)

   Return the real, complex or integer matrix :var:`mat` after inserting elements at the indexes :expr:`(i,j)` given by the :type:`iterable` :var:`ir` and :var:`jc` with the values given by :var:`a` depending of its kind:
   
   - if :var:`a` is a :type:`scalar`, it is will be used repetitively.

   - if :var:`a` is an :type:`iterable` then the matrix will be filled with values from :var:`a[n]` for :expr:`1 <= n <= #a`.
   
   The values after the inserted indexes are pushed toward the end of the matrix and discarded if they go beyond the last index. If :expr:`ir == nil`, :expr:`jc ~= nil` and :var:`a` is a 1D :type:`iterable`, then the latter is used to filled the matrix in the column-major order. Default: as :func:`mat:getidx()`.

.. function:: mat:getrow (ir, r_)

    Equivalent to :func:`mat:getsub()` with :expr:`jc = nil`.

.. function:: mat:setrow (ir, a, p_, s_)

   Equivalent to :func:`mat:setsub()` with :expr:`jc = nil`.

.. function:: mat:swprow (ir, ir2)

   Equivalent to :func:`mat:swpsub()` with :expr:`jc = nil` and :expr:`jc2 = nil`.

.. function:: mat:remrow (ir)

   Equivalent to :func:`mat:remsub()` with :expr:`jc = nil`.

.. function:: mat:insrow (ir, a)

   Equivalent to :func:`mat:inssub()` with :expr:`jc = nil`.

.. function:: mat:getcol (jc, r_)

   Equivalent to :func:`mat:getsub()` with :expr:`ir = nil`.

.. function:: mat:setcol (jc, a, p_, s_)

   Equivalent to :func:`mat:setsub()` with :expr:`ir = nil`.

.. function:: mat:swpcol (jc, jc2)

   Equivalent to :func:`mat:swpsub()` with :expr:`ir = nil` and :expr:`ir2 = nil`.

.. function:: mat:remcol (jc)

   Equivalent to :func:`mat:remsub()` with :expr:`ir = nil`.

.. function:: mat:inscol (jc, a)

   Equivalent to :func:`mat:inssub()` with :expr:`ir = nil`.

.. function:: mat:getdiag ([k_,] r_)

   Return a column vector of length :expr:`min(nrow, ncol)` or :var:`r` containing the values of the elements on the :math:`k`-th diagonal of the real, complex or integer matrix :var:`mat` using :expr:`mat:getvec()`. Default: as :func:`mat:getdidx()`.

.. function:: mat:setdiag (a, [k_,] p_, s_)

   Return the real, complex or integer matrix :var:`mat` after filling the elements on its :math:`k`-th diagonal with the values given by :var:`a` using :func:`mat:setvec()`. Default: as :func:`mat:getdidx()`.

Cloning and Reshaping
---------------------

.. function:: mat:same ([nr_, nc_,] v_)

   Return a matrix with elements of the type of :var:`v` and with :var:`nr` rows and :var:`nc` columns. Default: :expr:`v_ = mat[1]`, :expr:`nr_ = nrow`, :expr:`nc_ = ncol`.

.. function:: mat:copy (r_)

   Return a copy of the real, complex or integer matrix :var:`mat`.

.. function:: mat:reshape (nr_, nc_)

   Return the real, complex or integer matrix :var:`mat` reshaped to the new sizes :var:`nr` and :var:`nc` that must result into an equal or smaller number of elements, or it will raise an *invalid new sizes* error. Default: :expr:`nr_ = #mat`, :expr:`nc_ = 1`.

.. function:: mat:_reshapeto (nr, nc_)

   Same as :func:`mat:reshape()` except that :var:`nr` must be explicitly provided as this method allows for a larger size than :var:`mat` current size. A typical use of this method is to expand a vector after an explicit shrinkage, while keeping track of its original size, e.g. similar to :expr:`vector(100) :reshape(1):seti(1,1) :_reshapeto(2):seti(2,1)` that would raise an *"index out of bounds"* error without the call to :func:`_reshapeto()`. Default :expr:`nc_ = 1`.

   *WARNING: This method is unsafe and may crash MAD-NG with, e.g. a* `Segmentation fault <https://en.wikipedia.org/wiki/Segmentation_fault>`__ *, if wrongly used. It is the responsibility of the user to ensure that* :var:`mat` *was created with a size greater than or equal to the new size.* 

Matrix Properties
-----------------

.. function:: mat:is_const (tol_)

   Return true if all elements are equal within the tolerance :var:`tol`, false otherwise. Default: :expr:`tol_ = 0`. 

.. function:: mat:is_real (tol_)

   Return true if the imaginary part of all elements are equal to zero within the tolerance :var:`tol`, false otherwise. Default: :expr:`tol_ = 0`. 

.. function:: mat:is_imag (tol_)

   Return true if the real part of all elements are equal to zero within the tolerance :var:`tol`, false otherwise. Default: :expr:`tol_ = 0`. 

.. function:: mat:is_diag (tol_)

   Return true if all elements off the diagonal are equal to zero within the tolerance :var:`tol`, false otherwise. Default: :expr:`tol_ = 0`. 

.. function:: mat:is_symm ([tol_,] [sk_,] c_)

   Return true if :var:`mat` is a `symmetric matrix <https://en.wikipedia.org/wiki/Symplectic_matrix>`_, i.e. :math:`M = M^*` within the tolerance :var:`tol`, false otherwise. It checks for a `skew-symmetric matrix <https://en.wikipedia.org/wiki/Skew-symmetric_matrix>`_ if :expr:`sk == true`, and for a `Hermitian matrix <https://en.wikipedia.org/wiki/Hermitian_matrix>`_ if :expr:`c ~= false`, and a `skew-Hermitian matrix <https://en.wikipedia.org/wiki/Skew-Hermitian_matrix>`_ if both are combined. Default: :expr:`tol_ = 0`.

.. function:: mat:is_symp (tol_)

   Return true if :var:`mat` is a `symplectic matrix <https://en.wikipedia.org/wiki/Symplectic_matrix>`_, i.e. :math:`M^* S_{2n} M = S_{2n}`
   within the tolerance :var:`tol`, false otherwise. Default: :expr:`tol_ = eps`.

Filling and Moving
------------------

.. function:: mat:zeros ()

   Return the real, complex or integer matrix :var:`mat` filled with zeros.

.. function:: mat:ones (v_)

   Return the real, complex or integer matrix :var:`mat` filled with the value of :var:`v`. Default: :expr:`v_ = 1`.

.. function:: mat:eye (v_)

   Return the real, complex or integer matrix :var:`mat` filled with the value of :var:`v` on the diagonal and zeros elsewhere. The name of this method comes from the spelling of the `Identity matrix <https://en.wikipedia.org/wiki/Identity_matrix>`_ :math:`I`. Default: :expr:`v_ = 1`.

.. function:: mat:seq ([v_,] d_)

   Return the real, complex or integer matrix :var:`mat` filled with the indexes of the elements (i.e. starting at 1) and shifted by the value of :var:`v`. The matrix is filled in the column-major direction for :expr:`d == 'col'` and in the row-major direction otherwise. Default: :expr:`v_ = 0`, :expr:`d_ = 'row'`.

.. function:: mat:random (f_, ...)

   Return the real or complex matrix :var:`mat` filled with random values generated by :var:`f(...)`, and called twice for each element of :type:`cmatrix`. Default: :expr:`f_ = math.random`.

.. function:: mat:shuffle ()

   Return the real, complex or integer matrix :var:`mat` with its elements randomly swapped using the `Fisher–Yates or Knuth shuffle <https://en.wikipedia.org/wiki/Fisher–Yates_shuffle>`_ algorithm and :func:`math.random` as PRNG.

.. function:: mat:symp ()

   Return the real or complex matrix :var:`mat` filled with the block diagonal unitary `Symplectic matrix <https://en.wikipedia.org/wiki/Symplectic_matrix>`_ sometimes named :math:`J_{2n}` or :math:`S_{2n}`. The matrix :var:`mat` must be square with even number of rows and columns otherwise a *"2n square matrix expected"* error is raised.

.. function:: mat:circ (v)

   Return the real or complex matrix :var:`mat` filled as a `Circulant matrix <https://en.wikipedia.org/wiki/Circulant_matrix>`_ using the values from the :type:`iterable` :var:`v`, and rotating elements after each row or column depending on the shape of :var:`v`.

.. function:: mat:fill (a, p_, s_)

   Return the real, complex or integer matrix :var:`mat` filled with values provided by :var:`a` depending of its kind:

   - if :var:`a` is a :type:`scalar`, it is equivalent to :func:`mat:ones()`.

   - if :var:`a` is a :type:`callable`, then:

     - if :var:`p` and :var:`s` are provided, then :var:`a` is considered as a *stateless iterator*, and the matrix will be filled with the values :var:`v` returned by iterating :expr:`s, v = a(p, s)`.

     - otherwise :var:`a` is considered as a *generator*, and the matrix will be filled with values returned by calling :expr:`a(mat:get(i,j), i, j)`.

   - if :var:`a` is an :type:`iterable` then:
   
      - if :var:`a[1]` is also an :type:`iterable`, the matrix will be filled with the values from :var:`a[i][j]` for :expr:`1 <= i <= nrow` and :expr:`1 <= j <= ncol`, i.e. treated as a 2D container.

      - otherwise the matrix will be filled with values from :var:`a[n]` for :expr:`1 <= n <= #mat`, i.e. treated as a 1D container.

.. function:: mat:roll (nr_, nc_)

   Return the real, complex or integer matrix :var:`mat` after rolling its rows and columns by :var:`nr` and :var:`nc` respectively. Default: :expr:`nr_ = 0`, :expr:`nc_ = 0`.  

.. function:: mat:movev (i, j, k, r_)

   Return the real, complex or integer matrix :var:`r` after moving the elements in :expr:`mat[i..j]` to :expr:`r[k..k+j-i]` with :expr:`1 <= i <= j <= #mat` and :expr:`1 <= k <= k+j-i <= #r`. Default: :expr:`r_ = mat`.

.. function:: mat:shiftv (i, n_)

   Return the real, complex or integer matrix :var:`mat` after shifting the elements in :expr:`mat[i..]` to :expr:`mat[i+n..]` if :expr:`n > 0` and in the opposite direction if :expr:`n < 0`, i.e. it is equivalent to :expr:`mat:movev(i, #mat-n, i+n)` for :expr:`n > 0` and to :expr:`mat:movev(i-n, #mat+n, i)` for :expr:`n < 0`. Default: :expr:`n_ = 1`.

Mapping and Folding
-------------------

   This section lists the high-order functions `map <https://en.wikipedia.org/wiki/Map_(higher-order_function)>`_, `fold <https://en.wikipedia.org/wiki/Fold_(higher-order_function)>`_ and their variants useful in `functional programming <https://en.wikipedia.org/wiki/Functional_programming>`_ [#f1]_, followed by sections that list their direct application.

.. function:: mat:foreach ([ij_,] f, r_)

   Return the real, complex or integer matrix :var:`mat` after applying the :type:`callable` (or the operator string) :var:`f` to the elements at the indexes given by the :type:`iterable` :var:`ij` using :expr:`f(mat[n], n)`, i.e. interpreting the matrix as a vector. If :var:`r` is provided then it is filled with the values returned by :var:`f`. If :expr:`r = 'in'` then it is assigned :var:`mat`. Default: :expr:`ij_ = 1..#mat`.

.. function:: mat:filter ([ij_,] p, r_)

   Return a matrix or :var:`r` filled with the values of the elements of the real, complex or integer matrix :var:`mat` at the indexes given by the :type:`iterable` :var:`ij` if they are selected by the :type:`callable` `predicate <https://en.wikipedia.org/wiki/First-order_logic>`_ :var:`p` using :expr:`p(mat[n], n) == true`, i.e. interpreting the matrix as a vector. This method returns next to the matrix, a :type:`table` if :var:`r` is a table or a :type:`ivector` otherwise, containing the indexes of the selected elements. If :expr:`r = 'in'` then it is assigned :var:`mat`. Default: :expr:`ij_ = 1..#mat`.

.. function:: mat:filter_out ([ij_,] p, r_)

   Equivalent to :expr:`map:filter(ij_, compose(lnot,p), r_)`, where the functions :func:`compose()` and :func:`lnot()` are provided by the module :mod:`MAD.gfunc`.

.. function:: mat:map ([ij_,] f, r_)

   Return a matrix or :var:`r` filled with the values returned by the :type:`callable` (or the operator string) :var:`f` applied to the elements of the real, complex or integer matrix :var:`mat` at the indexes given by the :type:`iterable` :var:`ij` using :expr:`f(mat[n], n)`, i.e. interpreting the matrix as a vector. If :expr:`r = 'in'` or :expr:`r = nil` and expr:`ij ~= nil` then it is assigned :var:`mat`, i.e. map in place. If :expr:`r = nil`, the type of the returned matrix depends on the type of the first values returned by :expr:`f(mat[1], 1)`. Default: :expr:`ij_ = 1..#mat`.

.. function:: mat:map2 (y, [ij_,] f, r_)

   Equivalent to :func:`mat:map()` but with two arguments passed to :var:`f`, i.e. using :expr:`f(mat[n], y[n], n)`.

.. function:: mat:map3 (y, z, [ij_,] f, r_)

   Equivalent to :func:`mat:map()` but with three arguments passed to :var:`f`, i.e. using :expr:`f(mat[n], y[n], z[n], n)`.

.. function:: mat:foldl (f, [x0_,] [d_,] r_)

   Return a scalar, a vector or :var:`r` filled with the values returned by the :type:`callable` (or the operator string) :var:`f` applied iteratively to the elements of the real, complex or integer matrix :var:`mat` using the folding left (forward with increasing indexes) expression :expr:`v = f(v, mat[n])` in the direction depending on the :type:`string` :var:`d`:

   - If :expr:`d = 'vec'`, the folding left iteration runs on the entire matrix :var:`mat` interpreted as a vector and a scalar is returned.

   - If :expr:`d = 'diag'`, the folding left iteration runs on the diagonal of the matrix :var:`mat` and a scalar is returned.

   - If :expr:`d = 'row'`, the folding left iteration runs on the rows of the matrix :var:`mat` and a column vector is returned.

   - If :expr:`d = 'col'`, the folding left iteration runs on the columns of the matrix :var:`mat` and a row vector is returned.

   Default: :expr:`x0 = mat[1]` (or first row or column element), :expr:`d = 'vec'`.

.. function:: mat:foldr (f, [x0_,] [d_,] r_)

   Same as :func:`mat:foldl()` but the :type:`callable` (or the operator string) :var:`f` is applied iteratively using the folding right (backward with decreasing indexes) expression :expr:`v = f(mat[n], v)`. Default: :expr:`x0 = mat[#mat]` (or last row or column element), :expr:`d = 'vec'`.

.. function:: mat:scanl (f, [x0_,] [d_,] r_)

   Return a vector, a matrix or :var:`r` filled with the values returned by the :type:`callable` (or the operator string) :var:`f` applied iteratively to the elements of the real, complex or integer matrix :var:`mat` using the scanning left (forward with increasing indexes) expression :expr:`v = f(v, mat[n])` in the direction depending on the :type:`string` :var:`d`:

   - If :expr:`d = 'vec'`, the sanning left iteration runs on the entire matrix :var:`mat` interpreted as a vector and a vector is returned.

   - If :expr:`d = 'diag'`, the sanning left iteration runs on the diagonal of the matrix :var:`mat` and a vector is returned.

   - If :expr:`d = 'row'`, the sanning left iteration runs on the rows of the matrix :var:`mat` and a matrix is returned.

   - If :expr:`d = 'col'`, the sanning left iteration runs on the columns of the matrix :var:`mat` and a matrix is returned.

   Default: :expr:`x0 = mat[1]` (or first row or column element), :expr:`d = 'vec'`.

.. function:: mat:scanr (f, [x0_,] [d_,] r_)

   Same as :func:`mat:scanl()` but the :type:`callable` (or the operator string) :var:`f` is applied iteratively using the folding right (backward with decreasing indexes) expression :expr:`v = f(mat[n], v)`. Default: :expr:`x0 = mat[#mat]` (or last row or column element), :expr:`d = 'vec'`.

Mapping Real-like Methods
-------------------------

The following table lists the methods built from the application of :func:`mat:map()` and variants to the real-like functions from the module :mod:`MAD.gmath` for :type:`matrix` and :type:`cmatrix`. Only the mehods :func:`mat:abs()`, :func:`mat:sqr()` and :func:`mat:sign()` are available for :type:`imatrix`.

==========================  ===============================
Functions                   Equivalent Mapping
==========================  ===============================
:func:`mat:abs(r_)`         :expr:`mat:map(abs,r_)`
:func:`mat:acos(r_)`        :expr:`mat:map(acos,r_)`
:func:`mat:acosh(r_)`       :expr:`mat:map(acosh,r_)`
:func:`mat:acot(r_)`        :expr:`mat:map(acot,r_)`
:func:`mat:acoth(r_)`       :expr:`mat:map(acoth,r_)`
:func:`mat:asin(r_)`        :expr:`mat:map(asin,r_)`
:func:`mat:asinh(r_)`       :expr:`mat:map(asinh,r_)`
:func:`mat:asinc(r_)`       :expr:`mat:map(asinc,r_)`
:func:`mat:asinhc(r_)`      :expr:`mat:map(asinhc,r_)`
:func:`mat:atan(r_)`        :expr:`mat:map(atan,r_)`
:func:`mat:atan2(y,r_)`     :expr:`mat:map2(y,atan2,r_)`
:func:`mat:atanh(r_)`       :expr:`mat:map(atanh,r_)`
:func:`mat:ceil(r_)`        :expr:`mat:map(ceil,r_)`
:func:`mat:cos(r_)`         :expr:`mat:map(cos,r_)`
:func:`mat:cosh(r_)`        :expr:`mat:map(cosh,r_)`
:func:`mat:cot(r_)`         :expr:`mat:map(cot,r_)`
:func:`mat:coth(r_)`        :expr:`mat:map(coth,r_)`
:func:`mat:exp(r_)`         :expr:`mat:map(exp,r_)`
:func:`mat:floor(r_)`       :expr:`mat:map(floor,r_)`
:func:`mat:frac(r_)`        :expr:`mat:map(frac,r_)`
:func:`mat:hypot(y,r_)`     :expr:`mat:map2(y,hypot,r_)`
:func:`mat:hypot3(y,z,r_)`  :expr:`mat:map3(y,z,hypot3,r_)`
:func:`mat:log(r_)`         :expr:`mat:map(log,r_)`
:func:`mat:log10(r_)`       :expr:`mat:map(log10,r_)`
:func:`mat:round(r_)`       :expr:`mat:map(round,r_)`
:func:`mat:sign(r_)`        :expr:`mat:map(sign,r_)`
:func:`mat:sign1(r_)`       :expr:`mat:map(sign1,r_)`
:func:`mat:sin(r_)`         :expr:`mat:map(sin,r_)`
:func:`mat:sinc(r_)`        :expr:`mat:map(sinc,r_)`
:func:`mat:sinh(r_)`        :expr:`mat:map(sinh,r_)`
:func:`mat:sinhc(r_)`       :expr:`mat:map(sinhc,r_)`
:func:`mat:sqr(r_)`         :expr:`mat:map(sqr,r_)`
:func:`mat:sqrt(r_)`        :expr:`mat:map(sqrt,r_)`
:func:`mat:tan(r_)`         :expr:`mat:map(tan,r_)`
:func:`mat:tanh(r_)`        :expr:`mat:map(tanh,r_)`
:func:`mat:trunc(r_)`       :expr:`mat:map(trunc,r_)`
==========================  ===============================

Mapping Complex-like Methods
----------------------------

The following table lists the methods built from the application of :func:`mat:map()` to the the complex-like functions from the module :mod:`MAD.gmath` for :type:`matrix` and :type:`cmatrix`.

==========================  ===============================
Functions                   Equivalent Mapping
==========================  ===============================
:func:`mat:cabs(r_)`        :expr:`mat:map(cabs,r_)`
:func:`mat:carg(r_)`        :expr:`mat:map(carg,r_)`
:func:`mat:conj(r_)`        :expr:`mat:map(conj,r_)`
:func:`mat:cplx(im_,r_)`    :expr:`mat:map2(im_ or 0, cplx, r_)`
:func:`mat:fabs(r_)`        :expr:`mat:map(fabs,r_)`
:func:`mat:imag(r_)`        :expr:`mat:map(imag,r_)`
:func:`mat:polar(r_)`       :expr:`mat:map(polar,r_)`
:func:`mat:proj(r_)`        :expr:`mat:map(proj,r_)`
:func:`mat:real(r_)`        :expr:`mat:map(real,r_)`
:func:`mat:rect(r_)`        :expr:`mat:map(rect,r_)`
:func:`mat:reim(re_, im_)`  :expr:`re_ and mat:real(re_), im_ and mat:imag(im_)`
==========================  ===============================

Mapping Error-like Methods
--------------------------

The following table lists the methods built from the application of :func:`mat:map()` to the error-like functions from the module :mod:`MAD.gmath` for :type:`matrix` and :type:`cmatrix`.

==========================  ===============================
Functions                   Equivalent Mapping
==========================  ===============================
:func:`mat:erf(r_)`         :expr:`mat:map(erf,r_)`
:func:`mat:erfc(r_)`        :expr:`mat:map(erfc,r_)`
:func:`mat:erfcx(r_)`       :expr:`mat:map(erfcx,r_)`
:func:`mat:erfi(r_)`        :expr:`mat:map(erfi,r_)`
:func:`mat:wf(r_)`          :expr:`mat:map(wf,r_)`
==========================  ===============================

Mapping Vector-like Methods
---------------------------

The following table lists the methods built from the application of :func:`mat:map2()` to the vector-like functions from the module :mod:`MAD.gfunc` for :type:`matrix`, :type:`cmatrix`, and :type:`imatrix`.

==========================  ===============================
Functions                   Equivalent Mapping
==========================  ===============================
:func:`mat:emul(mat2,r_)`   :expr:`mat:map2(mat2,mul,r_)`
:func:`mat:ediv(mat2,r_)`   :expr:`mat:map2(mat2,div,r_)`
:func:`mat:emod(mat2,r_)`   :expr:`mat:map2(mat2,mod,r_)`
:func:`mat:epow(mat2,r_)`   :expr:`mat:map2(mat2,pow,r_)`
==========================  ===============================

Folding Methods
---------------

The following table lists the methods built from the application of :func:`mat:foldl()` to the functions from the module :mod:`MAD.gmath` for :type:`matrix`, :type:`cmatrix`, and :type:`imatrix`.

==========================  ===============================
Functions                   Equivalent Folding
==========================  ===============================
:func:`mat:all(p,d_,r_)`    :expr:`mat:foldl(all(p),false,d_,r_)`
:func:`mat:any(p,d_,r_)`    :expr:`mat:foldl(any(p),true,d_,r_)`
:func:`mat:min(d_,r_)`      :expr:`mat:foldl(min,nil,d_,r_)`
:func:`mat:max(d_,r_)`      :expr:`mat:foldl(max,nil,d_,r_)`
:func:`mat:sum(d_,r_)`      :expr:`mat:foldl(add,nil,d_,r_)`
:func:`mat:prod(d_,r_)`     :expr:`mat:foldl(mul,nil,d_,r_)`
:func:`mat:sumsqr(d_,r_)`   :expr:`mat:foldl(sumsqrl,0,d_,r_)`
:func:`mat:sumabs(d_,r_)`   :expr:`mat:foldl(sumabsl,0,d_,r_)`
:func:`mat:minabs(d_,r_)`   :expr:`mat:foldl(minabsl,inf,d_,r_)`
:func:`mat:maxabs(d_,r_)`   :expr:`mat:foldl(maxabsl,0,d_,r_)`
==========================  ===============================

Where :func:`any()` and :func:`all()` are functions that bind the predicate :var:`p` to the propagation of the logical AND and the logical OR respectively, that can be implemented like:

   - :expr:`all = \p -> \r,x -> lbool(land(r, p(x)))`
   - :expr:`any = \p -> \r,x -> lbool(lor (r, p(x)))`

Scanning Methods
----------------

The following table lists the methods built from the application of :func:`mat:scanl()` and :func:`mat:scanr()` to the functions from the module :mod:`MAD.gmath` for :type:`matrix` and :type:`cmatrix`.

=============================  ===============================
Functions                      Equivalent Scanning
=============================  ===============================
:func:`mat:accmin(d_,r_)`      :expr:`mat:scanl(min,nil,d_,r_)`
:func:`mat:accmax(d_,r_)`      :expr:`mat:scanl(max,nil,d_,r_)`
:func:`mat:accsum(d_,r_)`      :expr:`mat:scanl(add,nil,d_,r_)`
:func:`mat:accprod(d_,r_)`     :expr:`mat:scanl(mul,nil,d_,r_)`
:func:`mat:accsumsqr(d_,r_)`   :expr:`mat:scanl(sumsqrl,0,d_,r_)`
:func:`mat:accsumabs(d_,r_)`   :expr:`mat:scanl(sumabsl,0,d_,r_)`
:func:`mat:accminabs(d_,r_)`   :expr:`mat:scanl(minabsl,inf,d_,r_)`
:func:`mat:accmaxabs(d_,r_)`   :expr:`mat:scanl(maxabsl,0,d_,r_)`
:func:`mat:raccmin(d_,r_)`     :expr:`mat:scanr(min,nil,d_,r_)`
:func:`mat:raccmax(d_,r_)`     :expr:`mat:scanr(max,nil,d_,r_)`
:func:`mat:raccsum(d_,r_)`     :expr:`mat:scanr(add,nil,d_,r_)`
:func:`mat:raccprod(d_,r_)`    :expr:`mat:scanr(mul,nil,d_,r_)`
:func:`mat:raccsumsqr(d_,r_)`  :expr:`mat:scanr(sumsqrr,0,d_,r_)`
:func:`mat:raccsumabs(d_,r_)`  :expr:`mat:scanr(sumabsr,0,d_,r_)`
:func:`mat:raccminabs(d_,r_)`  :expr:`mat:scanr(minabsr,inf,d_,r_)`
:func:`mat:raccmaxabs(d_,r_)`  :expr:`mat:scanr(maxabsr,0,d_,r_)`
=============================  ===============================

The method :func:`mat:accumulate()` is also available as a more common name (i.e. an alias) for :func:`mat:accsum()`.

.. _matrix-functions:

Matrix Functions
----------------

The following table lists the methods built from the application of :func:`mat:mfun()` to the real-like functions from the module :mod:`MAD.gmath` for :type:`matrix` and :type:`cmatrix`.

==========================  ===============================
Functions                   Equivalent Matrix Function
==========================  ===============================
:func:`mat:macos()`         :expr:`mat:mfun(acos)`
:func:`mat:macosh()`        :expr:`mat:mfun(acosh)`
:func:`mat:macot()`         :expr:`mat:mfun(acot)`
:func:`mat:macoth()`        :expr:`mat:mfun(acoth)`
:func:`mat:masin()`         :expr:`mat:mfun(asin)`
:func:`mat:masinh()`        :expr:`mat:mfun(asinh)`
:func:`mat:masinc()`        :expr:`mat:mfun(asinc)`
:func:`mat:masinhc()`       :expr:`mat:mfun(asinhc)`
:func:`mat:matan()`         :expr:`mat:mfun(atan)`
:func:`mat:matanh()`        :expr:`mat:mfun(atanh)`
:func:`mat:mcos()`          :expr:`mat:mfun(cos)`
:func:`mat:mcosh()`         :expr:`mat:mfun(cosh)`
:func:`mat:mcot()`          :expr:`mat:mfun(cot)`
:func:`mat:mcoth()`         :expr:`mat:mfun(coth)`
:func:`mat:mexp()`          :expr:`mat:mfun(exp)`
:func:`mat:mlog()`          :expr:`mat:mfun(log)`
:func:`mat:mlog10()`        :expr:`mat:mfun(log10)`
:func:`mat:msin()`          :expr:`mat:mfun(sin)`
:func:`mat:msinc()`         :expr:`mat:mfun(sinc)`
:func:`mat:msinh()`         :expr:`mat:mfun(sinh)`
:func:`mat:msinhc()`        :expr:`mat:mfun(sinhc)`
:func:`mat:msqrt()`         :expr:`mat:mfun(sqrt)`
:func:`mat:mtan()`          :expr:`mat:mfun(tan)`
:func:`mat:mtanh()`         :expr:`mat:mfun(tanh)`
==========================  ===============================

Operator-like Methods
---------------------

.. function:: mat:unm (r_)

   Equivalent to :expr:`mat:map(unm,r_)`, where :func:`unm()` is from module :mod:`gmath`.

.. function:: mat:add (a, r_)

   Equivalent to :expr:`mat + a` with the possibility to place the result in :var:`r`.

.. function:: mat:sub (a, r_)

   Equivalent to :expr:`mat - a` with the possibility to place the result in :var:`r`.

.. function:: mat:mul (a, r_)

   Equivalent to :expr:`mat * a` with the possibility to place the result in :var:`r`.

.. function:: mat:div (a, r_, rcond_)

   Equivalent to :expr:`mat / a` with the possibility to place the result in :var:`r`, and to specify the conditional number :var:`rcond` used by the solver to determine the effective rank of non-square systems. Default: :expr:`rcond = eps`.

.. function:: mat:inv (r_, rcond_)

   Equivalent to :expr:`mat.div(1, mat, r_, rcond_)`. 

.. function:: mat:mod (a, r_)

   Equivalent to :expr:`mat % a` with the possibility to place the result in :var:`r`.

.. function:: mat:pow (n, r_)

   Equivalent to :expr:`mat ^ n` with the possibility to place the result in :var:`r`.

.. function:: mat:tmul (mat2, r_)

   Return a real or complex matrix or :var:`r` filled with the product of the transpose of :var:`mat` by :var:`mat2`, i.e. equivalent to :expr:`mat:t() * mat2`.
   
.. function:: mat:mult (mat2, r_)

   Return a real or complex matrix or :var:`r` filled with the product of :var:`mat` by the transpose of :var:`mat2`, i.e. equivalent to :expr:`mat * mat2:t()`.

.. function:: mat:eq (a, tol_)

   Return :const:`false` if :var:`a` is any matrix with incompatible sizes or if any element differ in a one-to-one comparison by more than :var:`tol`, :const:`true` otherwise. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix. Default: :expr:`tol_ = 0`.

.. function:: mat:concat (mat2, d_, r_)

   Return a real, complex or integer matrix resulting from concatenation of :var:`mat` and :var:`mat2` in the direction determined by :var:`d_`:
   
   - vector-oriented (appended) for :expr:`d = 'vec'`
   - row-oriented (horizontal) for :expr:`d = 'row'`
   - column-oriented (vectical) for :expr:`d = 'col'`
   
   The type of the returned matrix is given by the type promotion between :var:`mat` and the first element of :var:`mat2` except for :type:`imatrix`. Default: :var:`d_ = 'row'`.

Special Methods
---------------

.. function:: mat:conjugate (r_)
              mat:conj (r_)

   Equivalent to :expr:`mat:map(conj,r_)`.

.. function:: mat:transpose ([r_,] c_)
              mat:t ([r_,] c_)

   Return a real, complex or integer matrix or :var:`r` resulting from the conjugate transpose of the matrix :var:`mat` unless :expr:`c == false`. If :expr:`r = 'in'` then it is assigned :var:`mat`.

.. function:: mat:sympconj (r_)
              mat:bar (r_)

   Return a real or complex matrix or :var:`r` resulting from the symplectic conjugate of the matrix :var:`mat`, with :math:`\bar{M} = -S_{2n} M^* S_{2n}`, and :math:`M^{-1} = \bar{M}` if :math:`M` is symplectic. If :expr:`r = 'in'` then it is assigned :var:`mat`.

.. function:: mat:symperr (r_)

   Return the norm of the symplectic deviation matrix given by :math:`M^* S_{2n} M - S_{2n}` of the real or complex matrix :var:`mat`. If :var:`r` is provided, it is filled with the symplectic deviation matrix.

.. function:: mat:trace ()
              mat:tr()

   Return the `Trace <https://en.wikipedia.org/wiki/Trace_(linear_algebra)>`_ of the real or complex :var:`mat`.

.. function:: mat:inner (y)
              mat:dot (y)

   Return the `Inner Product <https://en.wikipedia.org/wiki/Dot_product>`_ of the two real or complex matrices :var:`mat` and :var:`y` with compatible sizes, i.e. return :math:`x^* . y` interpreting matrices as vectors. Note that multiple dot products, i.e. not interpreting matrices as vectors, can be achieved with :func:`mat:tmul()`.

.. function:: mat:outer (y, r_)

   Return the real or complex matrix resulting from the `Outer Product <https://en.wikipedia.org/wiki/Outer_product>`_ of the two real or complex matrices :var:`mat` and :var:`y`, i.e. return :math:`x . y^*` interpreting matrices as vectors.

.. function:: mat:cross (y, r_)

   Return the real or complex matrix resulting from the `Cross Product <https://en.wikipedia.org/wiki/Cross_product>`_ of the two real or complex matrices :var:`mat` and :var:`y` with compatible sizes, i.e. return :math:`x \times y` interpreting matrices as a list of :math:`[3 \times 1]` column vectors.

.. function:: mat:mixed (y, z, r_)

   Return the real or complex matrix resulting from the `Mixed Product <https://en.wikipedia.org/wiki/Triple_product>`_ of the three real or complex matrices :var:`mat`, :var:`y` and :var:`z` with compatible sizes, i.e. return :math:`x^* . (y \times z)` interpreting matrices as a list of :math:`[3 \times 1]` column vectors.

.. function:: mat:norm ()

   Return the `Frobenius norm <https://en.wikipedia.org/wiki/Matrix_norm#Frobenius_norm>`_ of the matrix :math:`\| M \|_2`. Other :math:`L_p` matrix norms and variants can be easily calculated using already provided methods, e.g. :math:`L_1` :expr:`= mat:sumabs'col':max()`, :math:`L_{\infty}` :expr:`= mat:sumabs'row':max()`, and :math:`L_2` :expr:`= mat:svd():max()`.

.. function:: mat:distance (y)

   Equivalent to :expr:`(mat - y):norm()`.

.. function:: mat:unit (r_)

   Equivalent to :expr:`mat:div(mat:norm(), r_)`.

.. function:: mat:center (d_, r_)

   Equivalent to :expr:`mat:sub(mat:mean(), r_)`. If :expr:`d = 'vec'`, :expr:`d = 'row'` or :expr:`d = 'col'` then centering will be vector-wise, row-wise or column-wise respectively. Default: :expr:`d_ = 'vec'`.

.. function:: mat:angle (y, n_)

   Return the angle between the two real or complex vectors :var:`mat` and :var:`y` using the method :func:`mat:inner()`. If :var:`n` is provided, the sign of :expr:`mat:mixed(y, n)` is used to define the angle in :math:`[-\pi,\pi]`, otherwise it is defined in :math:`[0,\pi]`.

.. function:: mat:minmax (abs_)

   Return the minimum and maximum values of the elements of the real, complex or integer matrix :var:`mat`. If :expr:`abs == true`, it returns the minimum and maximum absolute values of the elements. Default: :expr:`abs_ = false`.

.. function:: mat:iminmax (abs_)

   Return the two vector-like indexes of the minimum and maximum values of the elements of the real, complex or integer matrix :var:`mat`. If :expr:`abs == true`, it returns the indexes of the minimum and maximum absolute values of the elements. Default: :expr:`abs_ = false`.

.. function:: mat:mean ()

   Equivalent to :expr:`mat:sum()/#mat`, i.e. interpreting the matrix as a vector.

.. function:: mat:variance ()

   Equivalent to :expr:`(mat - mat:mean()):sumsqr()/(#mat-1)`, i.e. return the unbiased estimator of the variance with second order correction, interpreting the matrix as a vector.

.. function:: mat:ksum ()
              mat:kdot (y)

   Same as :func:`mat:sum()` and :func:`mat:dot()` respectively, except that they use the more accurate `Kahan Babushka Neumaier <https://en.wikipedia.org/wiki/Kahan_summation_algorithm>`_ algorithm for the summation, e.g. the sum of the elements of the vector :math:`[1,10^{100},1,-10^{100}]` should return :math:`0` with :func:`sum()` and the correct answer :math:`2` with :func:`ksum()`.

.. function:: mat:kadd (a, x)

   Return the real or complex matrix :var:`mat` filled with the linear combination of the compatible matrices stored in :var:`x` times the scalars stored in :var:`a`, i.e. :expr:`mat = a[1]*x[1] + a[2]*x[2] ...`

.. function:: mat:eval (x0)

   Return the evaluation of the real or complex matrix :var:`mat` at the value :var:`x0`, i.e. interpreting the matrix as a vector of polynomial coefficients of increasing orders in :var:`x` evaluated at :expr:`x = x0` using `Horner's method <https://en.wikipedia.org/wiki/Horner%27s_method>`_.

Solvers and Decompositions
--------------------------

   Except for :func:`nsolve()`, the solvers hereafter are wrappers around the library `Lapack <https://netlib.org/lapack/explore-html/index.html>`_ [#f2]_.

.. function:: mat:solve (b, rcond_)

   Return the real or complex :math:`[ n \times p ]` matrix :var:`x` as the minimum-norm solution of the linear least square problem :math:`\min \| A x - B \|` where :math:`A` is the real or complex :math:`[ m \times n ]` matrix :var:`mat` and :math:`B` is a :math:`[ m \times p ]` matrix of the same type as :math:`A`, using LU, QR or LQ factorisation depending on the shape of the system. The conditional number :var:`rcond` is used by the solver to determine the effective rank of non-square system. This method also returns the rank of the system. Default: :expr:`rcond_ = eps`.

.. function:: mat:ssolve (b, rcond_)

   Return the real or complex :math:`[ n \times p ]` matrix :var:`x` as the minimum-norm solution of the linear least square problem :math:`\min \| A x - B \|` where :math:`A` is the real or complex :math:`[ m \times n ]` matrix :var:`mat` and :math:`B` is a :math:`[ m \times p ]` matrix of the same type as :math:`A`, using SVD factorisation. The conditional number :var:`rcond` is used by the solver to determine the effective rank of the system. This method also returns the rank of the system followed by the real :math:`[ \min(m,n) \times 1 ]` vector of singluar values. Default: :expr:`rcond_ = eps`.

.. function:: mat:gsolve (b, c, d)

   Return the real or complex :math:`[ n \times 1 ]` vector :var:`x` as the minimum-norm solution of the linear least square problem :math:`\min \| A x - C \|` under the constraint :math:`B x = D` where :math:`A` is the real or complex :math:`[ m \times n ]` matrix :var:`mat`, :math:`B` is a :math:`[ p \times n ]` matrix, :math:`C` is a :math:`[ m \times 1 ]` vector and :math:`D` is a :math:`[ p \times 1 ]` vector, all of the same type as :math:`A`, using QR or LQ factorisation depending on the shape of the system.This method also returns the norm of the residues and the status :var:`info`.

.. function:: mat:gmsolve (b, d)

   Return the real or complex :math:`[ n \times 1 ]` vector :var:`x` and :math:`[ p \times 1 ]` matrix :var:`y` as the minimum-norm solution of the linear Gauss-Markov problem :math:`\min_x \| y \|` under the constraint :math:`A x + B y = D` where :math:`A` is the :math:`[ m \times n ]` real or complex matrix :var:`mat`, :math:`B` is a :math:`[ m \times p ]` matrix, and :math:`D` is a :math:`[ m \times 1 ]` vector, both of the same type as :math:`A`, using QR or LQ factorisation depending on the shape of the system. This method also returns the status :var:`info`.

.. function:: mat:nsolve (b, tol_, nc_)

   Return the real :math:`[ n \times 1 ]` vector :var:`x` (of correctors kicks) as the minimum-norm solution of the linear (best-kick) least square problem :math:`\min \| A x - B \|` where :math:`A` is the real :math:`[ m \times n ]` (response) matrix :var:`mat` and :math:`B` is a real :math:`[ m \times 1 ]` vector (of monitors readings), using the MICADO [#f3]_ algorithm based on the Householder-Golub method [MICADO]_. The argument :var:`tol` is a convergence threshold (on the residues) to stop the (orbit) correction if :math:`\| A x - B \| \leq m \times` :var:`tol`, and the argument :var:`nc` is the maximum number of correctors to use with :math:`0 < n_c \leq n`. This method also returns the updated number of correctors :math:`n_c` effectively used during the correction followed by the real :math:`[ m \times 1 ]` vector of residues. Default: :expr:`tol_ = eps`, :expr:`nc_ = ncol`, i.e. use all correctors.

.. function:: mat:pcacnd (ns_, rcond_)

   Return the integer column vector :var:`ic` containing the indexes of the columns to remove from the real or complex :math:`[ m \times n ]` matrix :var:`mat` using the Principal Component Analysis. The argument :var:`ns` is the maximum number of singular values to consider and :var:`rcond` is the conditionning number used to select the singular values versus the largest one, i.e. consider the :var:`ns` larger singular values :math:`\sigma_i` such that :math:`\sigma_i > \sigma_{\max}\times`:var:`rcond`. This method also returns the real :math:`[ \min(m,n) \times 1 ]` vector of singluar values. Default: :expr:`ns_ = min(nrow,ncol)`, :expr:`rcond_ = eps`.

.. function:: mat:svdcnd (ns_, rcond_, tol_)

   Return the integer column vector :var:`ic` containing the indexes of the columns to remove from the real or complex :math:`[ m \times n ]` matrix :var:`mat` based on the analysis of the right matrix :math:`V` from the SVD decomposition :math:`U S V`. The argument :var:`ns` is the maximum number of singular values to consider and :var:`rcond` is the conditionning number used to select the singular values versus the largest one, i.e. consider the :var:`ns` larger singular values :math:`\sigma_i` such that :math:`\sigma_i > \sigma_{\max}\times`:var:`rcond`. The argument :var:`tol` is a threshold similar to :var:`rcond` used to reject components in :math:`V` that have similar or opposite effect than components already encountered. This method also returns the real :math:`[ \min(m,n) \times 1 ]` vector of singluar values. Default: :expr:`ns_ = min(nrow,ncol)`, :expr:`rcond_ = eps`.

.. function:: mat:svd ()

   Return the real vector of the singular values and the two real or complex matrices resulting from the `SVD factorisation <https://en.wikipedia.org/wiki/Singular_value_decomposition>`_ of the real or complex matrix :var:`mat`, followed the status :var:`info`. The singular values are positive and sorted in decreasing order of values, i.e. largest first.

.. function:: mat:eigen (vr_, vl_)

   Return the complex vector filled with the eigenvalues followed by the by the status :var:`info` and the two optional real or complex matrices :var:`vr` and :var:`vl` containing the left and right eigenvectors resulting from the `Eigen Decomposition <https://en.wikipedia.org/wiki/Eigendecomposition_of_a_matrix>`_ of the real or complex square matrix :var:`mat`. The eigenvectors are normalized to have unit Euclidean norm and their largest component real, and satisfy :math:`A v_r = \lambda v_r` and :math:`v_l A = \lambda v_l`.

.. function:: mat:det ()

   Return the `Determinant <https://en.wikipedia.org/wiki/Determinant>`_ of the real or complex square matrix :var:`mat` using LU factorisation for better numerical stability, followed by the status :var:`info`.

.. function:: mat:mfun (fun)

   Return the real or complex matrix resulting from the matrix function :var:`fun` applyied to the real or complex matrix :var:`mat`. So far, :func:`mat:mfun()` uses the eigen decomposition of the matrix :var:`mat`, which must be `Diagonalizable <https://en.wikipedia.org/wiki/Diagonalizable_matrix>`_. See the section :ref:`matrix-functions` for the list of matrix functions already provided. Future versions of this method may be extended to use the more general Schur-Parlett algorithm [MATFUN]_, and other specialized versions for :func:`msqrt()`, :func:`mpow`, :func:`mexp`, and :func:`mlog` may be implemented too.

Fourier Transforms and Convolutions
-----------------------------------

The methods described is this section are based on the `FFTW <https://fftw.org>`_ and `NFFT <https://www-user.tu-chemnitz.de/~potts/nfft/>`_ libraries.

.. function:: mat:fft ([d_,] r_)

   Return the complex :math:`[n_r \times n_c]` vector, matrix or :var:`r` resulting from the 1D or 2D `Fourier Transform <https://en.wikipedia.org/wiki/Fourier_transform>`_ of the real or complex :math:`[n_r \times n_c]` vector or matrix :var:`mat` in the direction given by :var:`d`:

   - vector-oriented for :expr:`d = 'vec'`, return a 1D FFT of length :math:`n_r n_c`.
   - row-oriented for :expr:`d = 'row'`, return :math:`n_r` 1D FFTs of length :math:`n_c`.
   - column-oriented for :expr:`d = 'col'`, return :math:`n_c` 1D FFTs of length :math:`n_r`.
   - otherwise, return a 2D FFT of sizes :math:`[n_r \times n_c]`.

.. function:: mat:ifft ([d_,] r_)

   Return the complex :math:`[n_r \times n_c]` vector, matrix or :var:`r` resulting from the 1D or 2D inverse `Fourier Transform <https://en.wikipedia.org/wiki/Fourier_transform>`_ of the complex :math:`[n_r \times n_c]` vector or matrix :var:`mat`. See :func:`mat:fft()` for the direction :var:`d`.

.. function:: mat:rfft ([d_,] r_)

   Return the complex :math:`[n_r \times \lfloor n_c/2+1\rfloor]` vector, matrix or :var:`r` resulting from the 1D or 2D `Fourier Transform <https://en.wikipedia.org/wiki/Fourier_transform>`_ of the *real* :math:`[n_r \times n_c]` vector or matrix :var:`mat`. This method used an optimized version of the FFT for real data, which is about twice as fast and compact as the method :func:`mat:fft()`. See :func:`mat:fft()` for the direction :var:`d`.

.. function:: mat:irfft ([d_,] r)

   Return the *real* :math:`[n_r \times n_c]` vector, matrix or :var:`r` resulting from the 1D or 2D inverse `Fourier Transform <https://en.wikipedia.org/wiki/Fourier_transform>`_ of the complex :math:`[n_r \times \lfloor n_c/2+1\rfloor]` vector or matrix :var:`mat` as computed by the method :func:`mat:rfft()`. See :func:`mat:fft()` for the direction :var:`d`. Note that :var:`r` must be provided to specify the correct :math:`n_c` of the result.

.. function:: mat:nfft (p_, r_)

   Return the complex vector, matrix or :var:`r` resulting from the 1D or 2D *Nonequispaced* `Fourier Transform <https://en.wikipedia.org/wiki/Fourier_transform>`_ of the real or complex vector or matrix :var:`mat` respectively at :var:`p` time nodes. 

.. function:: mat:infft (p_, r_)

   Return the complex vector, matrix or :var:`r` resulting from the 1D or 2D *Nonequispaced* inverse `Fourier Transform <https://en.wikipedia.org/wiki/Fourier_transform>`_ of the real or complex vector or matrix :var:`mat` respectively at :var:`p` frequency nodes.

.. function:: mat:conv (y, [d_], r_)

   Return the real or complex vector, matrix or :var:`r` resulting from the 1D or 2D `Convolution <https://en.wikipedia.org/wiki/Convolution>`_ between the compatible real or complex vectors or matrices :var:`mat` and :var:`y` respectively. See :func:`mat:fft()` for the direction :var:`d`.

.. function:: mat:corr (y, [d_], r_)

   Return the real or complex vector, matrix or :var:`r` resulting from the 1D or 2D `Correlation <https://en.wikipedia.org/wiki/Cross-correlation>`_ between the compatible real or complex vectors or matrices :var:`mat` and :var:`y` respectively. See :func:`mat:fft()` for the direction :var:`d`.

.. function:: mat:covar (y, [d_,] r_)

   Return the real or complex vector, matrix or :var:`r` resulting from the 1D or 2D `Covariance <https://en.wikipedia.org/wiki/Covariance>`_ between the compatible real or complex vectors or matrices :var:`mat` and :var:`y` respectively. See :func:`mat:fft()` for the direction :var:`d`.

Rotations
---------

This section describe methods dealing with 2D and 3D rotations (see `Rotation Matrix <https://en.wikipedia.org/wiki/Rotation_matrix>`_) with angles in radians and trigonometric (counter-clockwise) direction for a right-handed frame, and where the following convention hold: :expr:`ax = -phi` is the *elevation* angle, :expr:`ay =  theta` is the *azimuthal* angle and :expr:`az =  psi` is the *roll/tilt* angle.

.. function:: mat:rot (a)

   Return the real :type:`matrix` :var:`mat` :math:`[2\times 2]` filled with a 2D rotation of angle :var:`a`.

.. function:: mat:rotx (a)
              mat:roty (a)
              mat:rotz (a)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation of angle :var:`a` around the x-axis, y-axis and z-axis respectively.

.. function:: mat:rotxy (ax, ay, inv_)
              mat:rotxz (ax, az, inv_)
              mat:rotyx (ay, ax, inv_)
              mat:rotyz (ay, az, inv_)
              mat:rotzx (az, ax, inv_)
              mat:rotzy (az, ay, inv_)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation of the first angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively *followed* by another 3D rotation of the second angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively of the frame rotated by the first rotation. If :var:`inv` is true, the returned matrix is the inverse rotation, i.e. the transposed matrix.

.. function:: mat:rotxyz (ax, ay, az, inv_)
              mat:rotxzy (ax, az, ay, inv_)
              mat:rotyxz (ay, ax, az, inv_)
              mat:rotyzx (ay, az, ax, inv_)
              mat:rotzxy (az, ax, ay, inv_)
              mat:rotzyx (az, ay, ax, inv_)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation of the first angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively *followed* by another 3D rotation of the second angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively of the frame rotated by the first rotation, and *followed* by a last 3D rotation of the third angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively of the frame already rotated by the two first rotations. If :var:`inv` is true, the returned matrix is the inverse rotations, i.e. the transposed matrix.

.. function:: mat:torotxyz (inv_)
              mat:torotxzy (inv_)
              mat:torotyxz (inv_)
              mat:torotyzx (inv_)
              mat:torotzxy (inv_)
              mat:torotzyx (inv_)

   Return three real :type:`number` representing the three angles :var:`ax`, :var:`ay` and :var:`az` (always in this order) of the 3D rotations stored in the real :type:`matrix` :var:`mat` :math:`[3\times 3]` by the methods with corresponding names. If :var:`inv` is true, the inverse rotations are returned, i.e. extracted from the transposed matrix.

.. function:: mat:rotv (v, av, inv_)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation of angle :var:`av` around the axis defined by the 3D vector-like :var:`v` (see `Axis-Angle representation <https://en.wikipedia.org/wiki/Axis–angle_representation>`_). If :var:`inv` is true, the returned matrix is the inverse rotation, i.e. the transposed matrix.

.. function:: mat:torotv (v_, inv_)

   Return a real :type:`number` representing the angle of the 3D rotation around the axis defined by a 3D vector as stored in the real :type:`matrix` :var:`mat` :math:`[3\times 3]` by the method :func:`mat:rotv()`. If the :type:`iterable`` :var:`v` is provided, it is filled with the components of the unit vector that defines the axis of the rotation.  If :var:`inv` is true, the inverse rotation is returned, i.e. extracted from the transposed matrix.

.. function:: mat:rotq (q, inv_)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation defined by the quaternion :var:`q` (see `Axis-Angle representation <https://en.wikipedia.org/wiki/Axis–angle_representation>`_). If :var:`inv` is true, the returned matrix is the inverse rotation, i.e. the transposed matrix.

.. function:: mat:torotq (q_, inv_)

   Return a quaternion representing the 3D rotation as stored in the real :type:`matrix` :var:`mat` :math:`[3\times 3]` by the method :func:`mat:rotq()`. If the :type:`iterable`` :var:`q` is provided, it is filled with the components of the quaternion otherwise the quaternion is returned in a :type:`list` of length 4.  If :var:`inv` is true, the inverse rotation is returned, i.e. extracted from the transposed matrix.

Conversions
-----------

.. function:: mat:tostring (sep_, lsep_)

   Return the string containing the real, complex or integer matrix converted to string. The argument :var:`sep` and :var:`lsep` are used as separator for columns and rows respectively. The elements values are formated using :func:`tostring()` that follows the :expr:`option.numfmt` string format for real numbers. Default: :expr:`sep = " "`, :expr:`lsep = "\n"`.

.. function:: mat:totable ([d_,] r_)

   Return the table or :var:`r` containing the real, complex or integer matrix converted to tables, i.e. one per row unless :var:`mat` is a vector or the direction :expr:`d = 'vec'`.  

Input and Output
----------------

.. function:: mat:write (filename_, name_, eps_, line_, nl_)

   Return the real, complex or integer matrix after writing it to the file :var:`filename` opened with :func:`MAD.utility.openfile()`. The content of the matrix :var:`mat` is preceded by a header containing enough information to read it back. If :var:`name` is provided, it is part of the header. If :expr:`line = 'line'`, the matrix is displayed on a single line with rows separated by a semicolumn, otherwise it is displayed on multiple lines separated by :var:`nl`. Elements with absolute value below :var:`eps` are displayed as zeros. The formats defined by :var:`MAD.option.numfmt` and :var:`MAD.option.intfmt` are used to format numbers of :type:`matrix`, :type:`cmatrix` and :type:`imatrix` respectively. Default: :expr:`filename_ = io.stdout`, :expr:`name_ = ''`, :expr:`eps_ = 0`, :expr:`line_ = nil`, :expr:`nl_ = '\\n'`.

.. function:: mat:print (name_, eps_, line_, nl_)

   Equivalent to :func:`mat:write(nil, name_, eps_, line_, nl_)`.

.. function:: mat:read (filename_)

   Return the real, complex or integer matrix read from the file :var:`filename` opened with :func:`MAD.utility.openfile()`. Note that the matrix :var:`mat` is only used to call the method :func:`:read()` and has no impact on the type and sizes of the returned matrix fully characterized by the content of the file. Default: :expr:`filename_ = io.stdin`.

Operators
=========

.. function:: #mat

   Return the size of the real, complex or integer matrix :var:`mat`, i.e. the number of elements interpreting the matrix as a vector.

.. function:: mat[n]

   Return the value of the element at index :var:`n` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= n <= #mat`, i.e. interpreting the matrix as a vector, :const:`nil` otherwise.

.. function:: mat[n] = v

   Assign the value :var:`v` to the element at index :var:`n` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= n <= #mat`, i.e. interpreting the matrix as a vector, otherwise raise an *"out of bounds"* error.

.. function:: -mat

   Return a real, complex or integer matrix resulting from the unary minus applied individually to all elements of the matrix :var:`mat`.

.. function:: num + mat
              mat + num
              mat + mat2

   Return a :type:`matrix` resulting from the sum of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: num + cmat
              cpx + mat
              cpx + cmat
              mat + cpx
              mat + cmat
              cmat + num
              cmat + cpx
              cmat + mat
              cmat + cmat2

   Return a :type:`cmatrix` resulting from the sum of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: idx + imat
              imat + idx
              imat + imat

   Return a :type:`imatrix` resulting from the sum of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: num - mat
              mat - num
              mat - mat2

   Return a :type:`matrix` resulting from the difference of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: num - cmat
              cpx - mat
              cpx - cmat
              mat - cpx
              mat - cmat
              cmat - num
              cmat - cpx
              cmat - mat
              cmat - cmat2

   Return a :type:`cmatrix` resulting from the difference of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: idx - imat
              imat - idx
              imat - imat

   Return a :type:`imatrix` resulting from the difference of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: num * mat
              mat * num
              mat * mat2

   Return a :type:`matrix` resulting from the product of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix. If the two operands are matrices, the mathematical `matrix multiplication <https://en.wikipedia.org/wiki/Matrix_multiplication>`_ is performed.

.. function:: num * cmat
              cpx * mat
              cpx * cmat
              mat * cpx
              mat * cmat
              cmat * num
              cmat * cpx
              cmat * mat
              cmat * cmat2

   Return a :type:`cmatrix` resulting from the product of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix. If the two operands are matrices, the mathematical `matrix multiplication <https://en.wikipedia.org/wiki/Matrix_multiplication>`_ is performed.

.. function:: idx * imat
              imat * idx

   Return a :type:`imatrix` resulting from the product of the left and right operands that must have compatible sizes. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: num / mat
              mat / num
              mat / mat2

   Return a :type:`matrix` resulting from the division of the left and right operands that must have compatible sizes. If the right operand is a scalar, the operator will be applied individually to all elements of the matrix. If the left operand is a scalar the operation :expr:`x/Y` is converted to :expr:`x (I/Y)` where :var:`I` is the identity matrix with compatible sizes. If the right operand is a matrix, the operation :expr:`X/Y` is performed using a system solver based on LU, QR or LQ factorisation depending on the shape of the system. 

.. function:: num / cmat
              cpx / mat
              cpx / cmat
              mat / cpx
              mat / cmat
              cmat / num
              cmat / cpx
              cmat / mat
              cmat / cmat2

   Return a :type:`cmatrix` resulting from the division of the left and right operands that must have compatible sizes. If the right operand is a scalar, the operator will be applied individually to all elements of the matrix. If the left operand is a scalar the operation :expr:`x/Y` is converted to :expr:`x (I/Y)` where :var:`I` is the identity matrix with compatible sizes. If the right operand is a matrix, the operation :expr:`X/Y` is performed using a system solver based on LU, QR or LQ factorisation depending on the shape of the system.

.. function:: imat / idx

   Return a :type:`imatrix` resulting from the division of the left and right operands, where the operator will be applied individually to all elements of the matrix.

.. function:: mat % num
              mat % mat

   Return a :type:`matrix` resulting from the modulo between the elements of the left and right operands that must have compatible sizes. If the right operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: cmat % num
              cmat % cpx
              cmat % mat
              cmat % cmat

   Return a :type:`cmatrix` resulting from the modulo between the elements of the left and right operands that must have compatible sizes. If the right operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: imat % idx
              imat % imat

   Return a :type:`imatrix` resulting from the modulo between the elements of the left and right operands that must have compatible sizes. If the right operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: mat ^ n
              cmat ^ n

   Return a :type:`matrix` or :type:`cmatrix` resulting from :var:`n` products of the square input matrix by itself. If :var:`n` is negative, the inverse of the matrix is used for the product.

.. function:: num == mat
              num == cmat
              num == imat
              cpx == mat
              cpx == cmat            
              mat == num
              mat == cpx
              mat == mat2
              mat == cmat
              cmat == num
              cmat == cpx
              cmat == mat
              cmat == cmat2
              imat == num
              imat == imat2

   Return :const:`false` if the left and right operands have incompatible sizes or if any element differ in a one-to-one comparison, :const:`true` otherwise. If one of the operand is a scalar, the operator will be applied individually to all elements of the matrix.

.. function:: mat .. mat2
              mat .. imat
              imat .. mat

   Return a :type:`matrix` resulting from the row-oriented (horizontal) concatenation of the left and right operands. If the first element of the right operand :var:`mat` (third case) is an integer, the resulting matrix will be a :type:`imatrix` instead.

.. function:: mat .. cmat
              imat .. cmat
              cmat .. mat
              cmat .. imat
              cmat .. cmat2

   Return a :type:`cmatrix` resulting from the row-oriented (horizontal) concatenation of the left and right operands.

.. function:: imat .. imat2

   Return a :type:`imatrix` resulting from the row-oriented (horizontal) concatenation of the left and right operands.

Iterators
=========

.. function:: ipairs(mat)
   :noindex:

   Return an :type:`ipairs` iterator suitable for generic :const:`for` loops. The generated values are those returned by :func:`mat[i]`. 

C API
=====

Real Vector
-----------

.. c:function:: void   mad_vec_zero   (                                           num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_seq    (                         num_t x        ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_fill   (                         num_t x        ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_shift  (       num_t x[],                                      ssz_t n, ssz_t d, int nshft)

.. c:function:: void   mad_vec_roll   (       num_t x[],                                      ssz_t n, ssz_t d, int nroll)

.. c:function:: void   mad_vec_copy   (const  num_t x[],                          num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_copyv  (const  num_t x[],                         cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_cvec   (const  num_t x[], const  num_t y[],       cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_minmax (const  num_t x[],       log_t abs       ,  idx_t r[2], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_center (const  num_t x[],                          num_t  r[], ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_abs    (const  num_t x[],                          num_t  r[], ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_eval   (const  num_t x[],        num_t x0,                     ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_sum    (const  num_t x[],                                      ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_ksum   (const  num_t x[],                                      ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_mean   (const  num_t x[],                                      ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_var    (const  num_t x[],                                      ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_norm   (const  num_t x[]                                     , ssz_t n, ssz_t d) 

.. c:function:: num_t  mad_vec_knorm  (const  num_t x[]                                     , ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_dist   (const  num_t x[], const  num_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_distv  (const  num_t x[], const cnum_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_dot    (const  num_t x[], const  num_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: num_t  mad_vec_kdot   (const  num_t x[], const  num_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: cnum_t mad_vec_dotv   (const  num_t x[], const cnum_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_dotv_r (const  num_t x[], const cnum_t y[]      , cnum_t *r  , ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_add    (const  num_t x[], const  num_t y[]      ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_addn   (const  num_t x[],        num_t y        ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_addc   (const  num_t x[],       cnum_t y        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_addc_r (const  num_t x[], num_t y_re, num_t y_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_kadd   (int k, const num_t a[], const num_t *x[],  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_sub    (const  num_t x[], const  num_t y[]      ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_subv   (const  num_t x[], const cnum_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_subn   (const  num_t y[],        num_t x        ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_subc   (const  num_t y[],       cnum_t x        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_subc_r (const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_mul    (const  num_t x[], const  num_t y[]      ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_muln   (const  num_t x[],        num_t y        ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_mulc   (const  num_t x[],       cnum_t y        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_mulc_r (const  num_t x[], num_t y_re, num_t y_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_div    (const  num_t x[], const  num_t y[]      ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_divv   (const  num_t x[], const cnum_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_divn   (const  num_t y[],        num_t x        ,  num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_divc   (const  num_t y[],       cnum_t x        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_divc_r (const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_vec_fft    (const  num_t x[],                         cnum_t  r[], ssz_t n)        

.. c:function:: void   mad_vec_rfft   (const  num_t x[],                         cnum_t  r[], ssz_t n)        

.. c:function:: void   mad_vec_nfft   (const  num_t x[], const num_t x_node[]  , cnum_t  r[], ssz_t n, ssz_t nr)

Complex Vector
--------------

.. c:function:: void   mad_cvec_zero  (                                          cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_seq   (                        cnum_t x        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_seq_r (                  num_t x_re, num_t x_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_fill  (                        cnum_t x        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_fill_r(                  num_t x_re, num_t x_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_shift (      cnum_t x[],                                      ssz_t n, ssz_t d, int nshft)

.. c:function:: void   mad_cvec_roll  (      cnum_t x[],                                      ssz_t n, ssz_t d, int nroll)

.. c:function:: void   mad_cvec_minmax(const cnum_t x[],                          idx_t r[2], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_center(const cnum_t x[],                         cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_copy  (const cnum_t x[],                         cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_vec   (const cnum_t x[],             num_t re[],  num_t ri[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_conj  (const cnum_t x[],                         cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: num_t  mad_cvec_abs   (const cnum_t x[],                          num_t  r[], ssz_t n, ssz_t d)

.. c:function:: cnum_t mad_cvec_eval  (const cnum_t x[],       cnum_t x0,                     ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_eval_r(const cnum_t x[],num_t x0_re,num_t x0_im, cnum_t *r  , ssz_t n, ssz_t d)

.. c:function:: cnum_t mad_cvec_sum   (const cnum_t x[],                                      ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_sum_r (const cnum_t x[],                         cnum_t *r  , ssz_t n, ssz_t d)

.. c:function:: cnum_t mad_cvec_mean  (const cnum_t x[],                                      ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_mean_r(const cnum_t x[],                         cnum_t *r  , ssz_t n, ssz_t d)

.. c:function:: cnum_t mad_cvec_var   (const cnum_t x[],                                      ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_var_r (const cnum_t x[],                         cnum_t *r  , ssz_t n, ssz_t d)

.. c:function:: num_t  mad_cvec_norm  (const cnum_t x[]                                     , ssz_t n, ssz_t d)

.. c:function:: num_t  mad_cvec_dist  (const cnum_t x[], const cnum_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: num_t  mad_cvec_distv (const cnum_t x[], const  num_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: cnum_t mad_cvec_dot   (const cnum_t x[], const cnum_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: cnum_t mad_cvec_dotv  (const cnum_t x[], const  num_t y[]                   , ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_dot_r (const cnum_t x[], const cnum_t y[]      , cnum_t *r  , ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_dotv_r(const cnum_t x[], const  num_t y[]      , cnum_t *r  , ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_add   (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_addv  (const cnum_t x[], const  num_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_addn  (const cnum_t x[],        num_t y        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_addc  (const cnum_t x[],       cnum_t y        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_addc_r(const cnum_t x[], num_t y_re, num_t y_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_kadd  (int k, const cnum_t a[],const cnum_t *x[],cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_sub   (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_subv  (const cnum_t x[], const  num_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_subn  (const cnum_t y[],        num_t x        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_subc  (const cnum_t y[],       cnum_t x        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_subc_r(const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_mul   (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_mulv  (const cnum_t x[], const  num_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_muln  (const cnum_t x[],        num_t y        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_mulc  (const cnum_t x[],       cnum_t y        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_mulc_r(const cnum_t x[], num_t y_re, num_t y_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_div   (const cnum_t x[], const cnum_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_divv  (const cnum_t x[], const  num_t y[]      , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_divn  (const cnum_t y[],        num_t x        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_divc  (const cnum_t y[],       cnum_t x        , cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_divc_r(const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_cvec_fft   (const cnum_t x[],                         cnum_t  r[], ssz_t n)        

.. c:function:: void   mad_cvec_nfft  (const cnum_t x[], const num_t x_node[]  , cnum_t  r[], ssz_t n, ssz_t nr)

.. c:function:: void   mad_cvec_ifft  (const cnum_t x[],                         cnum_t  r[], ssz_t n)          

.. c:function:: void   mad_cvec_irfft (const cnum_t x[],                          num_t  r[], ssz_t n)       

.. c:function:: void   mad_cvec_infft (const cnum_t x[], const num_t r_node[]  , cnum_t  r[], ssz_t n, ssz_t nx)

Integer Vector
--------------

.. c:function:: void   mad_ivec_zero  (                                           idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_seq   (                         idx_t x        ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_fill  (                         idx_t x        ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_shift (       idx_t x[],                                      ssz_t n, ssz_t d, int nshft)

.. c:function:: void   mad_ivec_roll  (       idx_t x[],                                      ssz_t n, ssz_t d, int nroll)

.. c:function:: void   mad_ivec_copy  (const  idx_t x[],                          idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_copyv (const  idx_t x[],                          num_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_minmax(const  idx_t x[],       log_t abs       ,  idx_t r[2], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_add   (const  idx_t x[], const  idx_t y[]      ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_addn  (const  idx_t x[],        idx_t y        ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_sub   (const  idx_t x[], const  idx_t y[]      ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_subn  (const  idx_t y[],        idx_t x        ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_mul   (const  idx_t x[], const  idx_t y[]      ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_muln  (const  idx_t x[],        idx_t y        ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_divn  (const  idx_t x[],        idx_t y        ,  idx_t  r[], ssz_t n, ssz_t d)

.. c:function:: void   mad_ivec_modn  (const  idx_t x[],        idx_t y        ,  idx_t  r[], ssz_t n, ssz_t d)

Real Matrix
-----------

.. c:function:: void   mad_mat_eye      (                         num_t x  ,        num_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_mat_seq      (                         num_t x  ,        num_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_mat_fill     (                         num_t x  ,        num_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_mat_roll     (       num_t x[],                                      ssz_t m, ssz_t n, int mroll, int nroll)

.. c:function:: void   mad_mat_copy     (const  num_t x[],                          num_t  r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)

.. c:function:: void   mad_mat_copym    (const  num_t x[],                         cnum_t  r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)

.. c:function:: void   mad_mat_trans    (const  num_t x[],                          num_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_mat_dot      (const  num_t x[], const  num_t y[],        num_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_mat_dotm     (const  num_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_mat_mul      (const  num_t x[], const  num_t y[],        num_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_mat_mulm     (const  num_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_mat_tmul     (const  num_t x[], const  num_t y[],        num_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_mat_tmulm    (const  num_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_mat_mult     (const  num_t x[], const  num_t y[],        num_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_mat_multm    (const  num_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: int    mad_mat_det      (const  num_t x[],                          num_t *r  ,          ssz_t n)                     

.. c:function:: int    mad_mat_invn     (const  num_t y[],        num_t x  ,        num_t  r[], ssz_t m, ssz_t n,          num_t rcond)

.. c:function:: int    mad_mat_invc     (const  num_t y[],       cnum_t x  ,       cnum_t  r[], ssz_t m, ssz_t n,          num_t rcond)

.. c:function:: int    mad_mat_invc_r   (const  num_t y[], num_t x_re, num_t x_im, cnum_t  r[], ssz_t m, ssz_t n,          num_t rcond)

.. c:function:: int    mad_mat_div      (const  num_t x[], const  num_t y[],        num_t  r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)

.. c:function:: int    mad_mat_divm     (const  num_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)

.. c:function:: int    mad_mat_solve    (const  num_t a[], const  num_t b[],        num_t  x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)

.. c:function:: int    mad_mat_nsolve   (const  num_t a[], const  num_t b[],        num_t  x[], ssz_t m, ssz_t n, ssz_t N, num_t rcond, num_t r_[])

.. c:function:: int    mad_mat_ssolve   (const  num_t a[], const  num_t b[],        num_t  x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond, num_t s_[])

.. c:function:: int    mad_mat_gsolve   (const  num_t a[], const  num_t b[], const  num_t  c[], const num_t d[], num_t  x[], ssz_t m, ssz_t n, ssz_t p, num_t *nrm_)

.. c:function:: int    mad_mat_gmsolve  (const  num_t a[], const  num_t b[], const  num_t  d[], num_t x[],        num_t  y[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: int    mad_mat_pcacnd   (const  num_t a[],        idx_t c[],                    ssz_t m, ssz_t n, ssz_t N, num_t cut, num_t s_[])

.. c:function:: int    mad_mat_svdcnd   (const  num_t a[],        idx_t c[],                    ssz_t m, ssz_t n, ssz_t N, num_t cut, num_t s_[], num_t tol)

.. c:function:: int    mad_mat_svd      (const  num_t x[], num_t u[], num_t s[],    num_t  v[], ssz_t m, ssz_t n)

.. c:function:: int    mad_mat_eigen    (const  num_t x[], cnum_t w[], num_t vl[],  num_t vr[],          ssz_t n)

.. c:function:: void   mad_mat_fft      (const  num_t x[],                         cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_mat_rfft     (const  num_t x[],                         cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_mat_nfft     (const  num_t x[], const num_t x_node[]  , cnum_t  r[], ssz_t m, ssz_t n, ssz_t nr)

.. c:function:: void   mad_mat_center   (const  num_t x[],                          num_t  r[], ssz_t m, ssz_t n, int d)

.. c:function:: void   mad_mat_sympconj (const  num_t x[],                          num_t  r[],          ssz_t n)

.. c:function:: num_t  mad_mat_symperr  (const  num_t x[],                          num_t  r[],          ssz_t n)

.. c:function:: num_t  mad_mat_vdot     (const  num_t x[], idx_t xs, const  num_t y[], idx_t ys,         ssz_t n)

Complex Matrix
--------------

.. c:function:: void   mad_cmat_eye     (                        cnum_t x  ,       cnum_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_cmat_eye_r   (                  num_t x_re, num_t x_im, cnum_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_cmat_seq     (                        cnum_t x  ,       cnum_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_cmat_seq_r   (                  num_t x_re, num_t x_im, cnum_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_cmat_fill    (                        cnum_t x  ,       cnum_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_cmat_fill_r  (                  num_t x_re, num_t x_im, cnum_t  r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_cmat_roll    (      cnum_t x[],                                      ssz_t m, ssz_t n, int mroll, int nroll)

.. c:function:: void   mad_cmat_copy    (const cnum_t x[],                         cnum_t  r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)

.. c:function:: void   mad_cmat_trans   (const cnum_t x[],                         cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_cmat_ctrans  (const cnum_t x[],                         cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_cmat_dot     (const cnum_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_cmat_dotm    (const cnum_t x[], const  num_t y[],       cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_cmat_mul     (const cnum_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_cmat_mulm    (const cnum_t x[], const  num_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_cmat_tmul    (const cnum_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_cmat_tmulm   (const cnum_t x[], const  num_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_cmat_mult    (const cnum_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: void   mad_cmat_multm   (const cnum_t x[], const  num_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: int    mad_cmat_det     (const cnum_t x[],                         cnum_t *r  ,          ssz_t n)

.. c:function:: int    mad_cmat_invn    (const cnum_t y[],        num_t x  ,       cnum_t  r[], ssz_t m, ssz_t n,          num_t rcond)

.. c:function:: int    mad_cmat_invc    (const cnum_t y[],       cnum_t x  ,       cnum_t  r[], ssz_t m, ssz_t n,          num_t rcond)

.. c:function:: int    mad_cmat_invc_r  (const cnum_t y[], num_t x_re, num_t x_im, cnum_t  r[], ssz_t m, ssz_t n,          num_t rcond)

.. c:function:: int    mad_cmat_div     (const cnum_t x[], const cnum_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)

.. c:function:: int    mad_cmat_divm    (const cnum_t x[], const  num_t y[],       cnum_t  r[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)

.. c:function:: int    mad_cmat_solve   (const cnum_t a[], const cnum_t b[],       cnum_t  x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond)

.. c:function:: int    mad_cmat_ssolve  (const cnum_t a[], const cnum_t b[],       cnum_t  x[], ssz_t m, ssz_t n, ssz_t p, num_t rcond, num_t s_[])

.. c:function:: int    mad_cmat_gsolve  (const cnum_t a[], const cnum_t b[], const cnum_t  c[], const cnum_t d[], cnum_t  x[], ssz_t m, ssz_t n, ssz_t p, num_t *nrm_)

.. c:function:: int    mad_cmat_gmsolve (const cnum_t a[], const cnum_t b[], const cnum_t  d[], cnum_t x[],       cnum_t  y[], ssz_t m, ssz_t n, ssz_t p)

.. c:function:: int    mad_cmat_pcacnd  (const cnum_t a[],        idx_t c[],                    ssz_t m, ssz_t n, ssz_t N, num_t cut, num_t s_[])

.. c:function:: int    mad_cmat_svd     (const cnum_t x[], cnum_t u[], num_t s[],  cnum_t  v[], ssz_t m, ssz_t n)

.. c:function:: int    mad_cmat_eigen   (const cnum_t x[], cnum_t w[], cnum_t vl[],cnum_t vr[],          ssz_t n)

.. c:function:: void   mad_cmat_fft     (const cnum_t x[],                         cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_cmat_nfft    (const cnum_t x[], const num_t x_node[]   ,cnum_t  r[], ssz_t m, ssz_t n, ssz_t nr)

.. c:function:: void   mad_cmat_ifft    (const cnum_t x[],                         cnum_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_cmat_irfft   (const cnum_t x[],                          num_t  r[], ssz_t m, ssz_t n)

.. c:function:: void   mad_cmat_infft   (const cnum_t x[], const num_t r_node[]   ,cnum_t  r[], ssz_t m, ssz_t n, ssz_t nx)

.. c:function:: void   mad_cmat_center  (const cnum_t x[],                         cnum_t  r[], ssz_t m, ssz_t n, int d)

.. c:function:: void   mad_cmat_sympconj(const cnum_t x[],                         cnum_t  r[],          ssz_t n)

.. c:function:: num_t  mad_cmat_symperr (const cnum_t x[],                         cnum_t  r[],          ssz_t n)

.. c:function:: cnum_t mad_cmat_vdot    (const cnum_t x[], idx_t xs, const cnum_t y[], idx_t ys,         ssz_t n)

.. c:function:: cnum_t mad_cmat_vdotm   (const cnum_t x[], idx_t xs, const  num_t y[], idx_t ys,         ssz_t n)

.. c:function:: void   mad_cmat_vdot_r  (const cnum_t x[], idx_t xs, const cnum_t y[], idx_t ys, cnum_t *r, ssz_t n)

.. c:function:: void   mad_cmat_vdotm_r (const cnum_t x[], idx_t xs, const  num_t y[], idx_t ys, cnum_t *r, ssz_t n)

Integer Matrix
--------------

.. c:function:: void   mad_imat_eye     (       idx_t x  ,                           idx_t r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_imat_seq     (       idx_t x  ,                           idx_t r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_imat_fill    (       idx_t x  ,                           idx_t r[], ssz_t m, ssz_t n,            ssz_t ldr)

.. c:function:: void   mad_imat_roll    (       idx_t x[],                                      ssz_t m, ssz_t n, int mroll, int nroll)

.. c:function:: void   mad_imat_copy    (const  idx_t x[],                           idx_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)

.. c:function:: void   mad_imat_copym   (const  idx_t x[],                           num_t r[], ssz_t m, ssz_t n, ssz_t ldx, ssz_t ldr)

.. c:function:: void   mad_imat_trans   (const  idx_t x[],                           idx_t r[], ssz_t m, ssz_t n)

Rotation Matrix
---------------

.. c:function:: void   mad_mat_rot      (      num_t x[2*2], num_t a)

.. c:function:: void   mad_mat_rotx     (      num_t x[3*3], num_t ax)

.. c:function:: void   mad_mat_roty     (      num_t x[3*3], num_t ay)

.. c:function:: void   mad_mat_rotz     (      num_t x[3*3], num_t az)

.. c:function:: void   mad_mat_rotxy    (      num_t x[3*3], num_t ax, num_t ay, log_t inv)

.. c:function:: void   mad_mat_rotxz    (      num_t x[3*3], num_t ax, num_t az, log_t inv)

.. c:function:: void   mad_mat_rotyz    (      num_t x[3*3], num_t ay, num_t az, log_t inv)

.. c:function:: void   mad_mat_rotxyz   (      num_t x[3*3], num_t ax, num_t ay, num_t az, log_t inv)

.. c:function:: void   mad_mat_rotxzy   (      num_t x[3*3], num_t ax, num_t ay, num_t az, log_t inv)

.. c:function:: void   mad_mat_rotyxz   (      num_t x[3*3], num_t ax, num_t ay, num_t az, log_t inv)

.. c:function:: void   mad_mat_torotxyz (const num_t x[3*3], num_t r[3]                  , log_t inv)

.. c:function:: void   mad_mat_torotxzy (const num_t x[3*3], num_t r[3]                  , log_t inv)

.. c:function:: void   mad_mat_torotyxz (const num_t x[3*3], num_t r[3]                  , log_t inv)

.. c:function:: void   mad_mat_rotv     (      num_t x[3*3], num_t v[3],         num_t av, log_t inv)

.. c:function:: num_t  mad_mat_torotv   (const num_t x[3*3], num_t v_[3]                 , log_t inv)

.. c:function:: void   mad_mat_rotq   (      num_t x[3*3], num_t q[4], log_t inv)

.. c:function:: void   mad_mat_torotq (const num_t x[3*3], num_t q[4], log_t inv)

Misalignments
-------------

.. c:function:: void   mad_mat_rtbar    (      num_t Rb[3*3],       num_t Tb[3], num_t el, num_t ang, num_t tlt,  const num_t R_[3*3], const num_t T [3])

Miscellaneous
-------------

.. c:function:: void   mad_fft_cleanup (void)

.. ------------------------------------------------------------

References
==========

.. [MICADO] B. Autin, and Y. Marti, *"Closed Orbit Correction of Alternating Gradient Machines using a Small Number of Magnets"*, CERN ISR-MA/73-17, Mar. 1973.

.. [MATFUN] N.J. Higham, and X. Liu, *"A Multiprecision Derivative-Free Schur–Parlett Algorithm for Computing Matrix Functions"*, SIAM J. Matrix Anal. Appl., Vol. 42, No. 3, pp. 1401-1422, 2021.

.. rubric:: Footnotes

.. [#f1] For *true* Functional Programming, see the module :mod:`MAD.lfun`, a binding of the `LuaFun <https://github.com/luafun/luafun>`_  library adapted to the ecosystem of MAD-NG.

.. [#f2] The solvers are based, among others, on the following Lapack drivers: 

   - :func:`dgesv()` and :func:`zgesv()` for LU factorization.
   - :func:`dgelsy()` and :func:`zgelsy()` for QR or LQ factorization.
   - :func:`dgelsd()` and :func:`zgelsd()` for SVD factorisation.
   - :func:`dgees()` and :func:`zgees()` for Schur factorisation.
   - :func:`dgglse()` and :func:`zgglse()` for equality-constrained linear Least Squares problems.
   - :func:`dggglm()` and :func:`zggglm()` for general Gauss-Markov linear model problems.

.. [#f3] MICADO stands for "Minimisation des CArrés des Distortions d'Orbite" in french. 
