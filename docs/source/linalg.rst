.. index::
   Linear algebra
   Vector and matrix

**************
Linear Algebra
**************

This chapter describes the real :type:`matrix`, complex :type:`cmatrix` and integer :type:`imatrix` objects as supported by MAD-NG. The module for `Matrix <https://en.wikipedia.org/wiki/Matrix_(mathematics)>`_ is not exposed, only the contructors are visible from the :mod:`MAD` environment and thus, matrices are handled directly by their methods or by the generic functions of the same name from the module :mod:`MAD.gmath`. The :type:`imatrix`, i.e. matrix of integers, are mainly used for indexing other types of matrix and therefore supports only a limited subset of the features. Column and row vectors are shortcuts for :math:`[n\times 1]` and :math:`[1\times n]` matrices respectively. Note that :type:`matrix`, :type:`cmatrix` and :type:`imatrix` are all defined as C structures for direct compliance with the C API. 

Types promotion
===============

The matrix operations may involve other data types like real and complex numbers leading to many combinations of types. In order to simplify the descriptions, the generic names :var:`num`, :var:`cpx` and :var:`idx` (indexes) are used for real, complex and integer-as-index numbers respectively, :var:`vec`, :var:`cvec` and :var:`ivec` for real, complex and integer vectors respectively, and :var:`mat`, :var:`cmat` and :var:`imat` for real, complex and integer matrices respectively. For example, the sum of a complex number :var:`cpx` and a real matrix :var:`mat` gives a complex matrix :var:`cmat`. The case of :var:`idx` means that a :type:`number` will be interpreted as an index and automatically rounded if it does not hold an integer value. The following table summarizes all valid combinations of types for binary operations involving at least one matrix type:

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

The constructors for vectors and matrices are directly available from the :mod:`MAD` environment. Note that real, complex or integer matrix with zero size are not allowed, i.e. the smallest allowed matrix has a size of :math:`[1\times 1]`.

.. function::  vector(nrow)
              cvector(nrow)
              ivector(nrow)

   Return a real, complex or integer column vector (i.e. a matrix of size :math:`[n_{\text{row}}\times 1]`) filled with zeros.

.. function::  matrix(nrow, ncol_)
              cmatrix(nrow, ncol_)
              imatrix(nrow, ncol_)

   Return a real, complex or integer matrix of size :math:`[n_{\text{row}}\times n_{\text{col}}]` filled with zeros. Default: :expr:`ncol_ = rnow`.

.. function:: linspace([start_,] stop, size_)

   Return a real or complex column vector of length :var:`size` filled with values equally spaced between :var:`start` and :var:`stop` on a linear scale. Default: :expr:`start_ = 0`, :expr:`size_ = 100`.

.. function:: logspace([start_,] stop, size_)

   Return a real or complex column vector of length :var:`size` filled with values equally spaced between :var:`start` and :var:`stop` on a logarithmic scale. Default: :expr:`start_ = 1`, :expr:`size_ = 100`.

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

   Return :const:`true` if :var:`a` is a real or complex vector, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

.. function:: is_matrix (a)
              is_cmatrix (a)
              is_imatrix (a)

   Return :const:`true` if :var:`a` is respectively a real, complex or integer matrix, :const:`false` otherwise. These functions are only available from the module :mod:`MAD.typeid`.

.. function:: isa_matrix (a)

   Return :const:`true` if :var:`a` is a real or complex matrix, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

Getters/Setters
---------------

.. function:: mat:sizes ()

   Return the number of rows and columns of the real, complex or integer matrix :var:`mat`.

.. function:: mat:tsizes ()

   Return the number of columns and rows (i.e. transposed) of the real, complex or integer matrix :var:`mat`.

.. function:: mat:get (i, j)

   Return the value of the element at indexes :var:`(i,j)` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= i <= nrow` and :expr:`1 <= j <= ncol`, :const:`nil` otherwise.

.. function:: mat:set (i, j, v)

   Assign the value :var:`v` to the element at indexes :var:`(i,j)` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= i <= nrow` and :expr:`1 <= j <= ncol` and return the matrix, otherwise raise an *"out of bounds"* error.

.. function:: mat:geti (n)

   Return the value of the element at index :var:`n` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= n <= #mat`, :const:`nil` otherwise.

.. function:: mat:seti (n, v)

   Assign the value :var:`v` to the element at index :var:`n` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= n <= #mat` and return the matrix, otherwise raise an *"out of bounds"* error.

Copy/Shape
----------

.. function:: mat:same ([nrow_, ncol_,] v_)

   Return a matrix with elements of the type of :var:`v` (ignored by :type:`imatrix`) and with :var:`nrow` rows and :var:`ncol` columns. Default: :expr:`v_ = mat[1]`, :expr:`nrow_ = mat.nrow`, :expr:`ncol_ = mat.ncol`.

.. function:: mat:copy (r_)

   Return a copy of the real, complex or integer matrix :var:`mat`.

.. function:: mat:reshape (nrow_, ncol_)

   Return the real, complex or integer matrix :var:`mat` reshaped to the new sizes :var:`nrow` and :var:`ncol` that must give to a smaller or an equal size, or raise an *invalid new sizes* error. Default: :expr:`nrow_ = #mat`, :expr:`ncol_ = 1`.

.. function:: mat:_reshapeto (nrow_, ncol_)

   Same as :func:`mat:reshape()` but allows for a new size larger than :var:`mat` current size.

   *WARNING: This method is unsafe and may crash MAD-NG, i.e. with a* `Segmentation fault <https://en.wikipedia.org/wiki/Segmentation_fault>`__ *, if wrongly used. It is the responsibility of the user to ensure that* :var:`mat` *contains enough allocated memory to fulfill the new sizes.* 

.. function:: vec:_appendto (v_)

   Return the real, complex or integer vector :var:`vec` with the value :var:`v` appended at its end and increments its number of rows or columns depending on the kind of vector. Default: :expr:`v_ = 0`.

   *WARNING: This method is unsafe and may crash MAD-NG, i.e. with a* `Segmentation fault <https://en.wikipedia.org/wiki/Segmentation_fault>`__ *, if wrongly used. It is the responsibility of the user to ensure that* :var:`vec` *is effectively a vector and contains enough allocated memory to append the value* :var:`v`. 

Filling/Moving
--------------

.. function:: mat:is_const (v_, tol_)

.. function:: mat:is_diag (v_, tol_)

.. function:: mat:zeros ()

.. function:: mat:ones (v_)

.. function:: mat:seq (v0_)

.. function:: mat:eye (v_)

.. function:: mat:random (f_, ...)

.. function:: mat:shuffle ()

.. function:: mat:symp ()

.. function:: mat:circ (v)

.. function:: mat:fill (a, p_, s_)

.. function:: mat:movev (i, j, yi, y_)

.. function:: mat:shiftv (ni, ns_)

.. function:: mat:roll (ns_, ms_)

Conversions
-----------

.. function:: mat:real (r_)

.. function:: mat:imag (r_)

.. function:: mat:rerim (re_, im_)

.. function:: mat:cplx (re_, im_, r_)

.. function:: mat:totable (v_, r_)

.. function:: mat:tostring (sep_, lsep_)

Operator-like Methods
---------------------

.. function:: mat:concat (mat2, dir_, r_)

   Return a real, complex or integer matrix resulting from concatenation of :var:`mat` and :var:`mat2` in the direction determined by :var:`dir_`:
   
   - row-oriented (horizontal) for :expr:`dir = 'row'`
   - column-oriented (vectical) for :expr:`dir = 'col'`
   - vector-oriented (appended) for :expr:`dir = 'vec'`
   
   The type of the returned matrix is given by the type promotion between :var:`mat` and the first element of :var:`mat2`. Default: :var:`dir_ = 'row'`.

Input/Output
------------

.. function:: mat:write (filename_, name_, eps_, line_, nl_)

   Return the real, complex or integer matrix after writing it to the file :var:`filename` opened with :func:`MAD.utility.openfile()`. The content of the matrix :var:`mat` is preceded by a header containing enough information to read it back. If :var:`name` is provided, it is part of the header. If :expr:`line = 'line'`, the matrix is displayed on a single line with rows separated by a semicolumn, otherwise it is displayed on multiple lines separated by :var:`nl`. Elements with absolute value below :var:`eps` are displayed as zeros. The formats defined by :var:`MAD.option.numfmt` and :var:`MAD.option.intfmt` are used to format numbers of :type:`matrix`, :type:`cmatrix` and :type:`imatrix` respectively. Default: :expr:`filename_ = io.stdout`, :expr:`name_ = ''`, :expr:`eps_ = 0`, :expr:`line_ = nil`, :expr:`nl_ = '\\n'`.

.. function:: mat:print (name_, eps_, line_, nl_)

   Equivalent to :func:`mat:write(nil, name_, eps_, line_, nl_)`.

.. function:: mat:read (filename_)

   Return the real, complex or integer matrix read from the file :var:`filename` opened with :func:`MAD.utility.openfile()`. Note that the matrix :var:`mat` is only used to call the method :func:`:read()` and has no impact on the type and sizes of the returned matrix fully characterized by the content of the file. Default: :expr:`filename_ = io.stdin`.

Rotations
---------

This section describe methods dealing with 2D and 3D rotations (see `Rotation Matrix <https://en.wikipedia.org/wiki/Rotation_matrix>`_) with angles in radians and trigonometric (counter-clockwise) direction for a right-handed frame, and where the following convention hold: :expr:`ax = -phi` is the *elevation* angle, :expr:`ay =  theta` is the *azimuthal* angle and :expr:`az =  psi` is the *roll/tilt* angle.

.. function:: mat:rot(a)

   Return the real :type:`matrix` :var:`mat` :math:`[2\times 2]` filled with a 2D rotation of angle :var:`a`.

.. function:: mat:rotx(a)
              mat:roty(a)
              mat:rotz(a)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation of angle :var:`a` around the x-axis, y-axis and z-axis respectively.

.. function:: mat:rotxy(ax, ay, inv_)
              mat:rotxz(ax, az, inv_)
              mat:rotyx(ay, ax, inv_)
              mat:rotyz(ay, az, inv_)
              mat:rotzx(az, ax, inv_)
              mat:rotzy(az, ay, inv_)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation of the first angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively *followed* by another 3D rotation of the second angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively of the frame rotated by the first rotation. If :var:`inv` is true, the returned matrix is the inverse rotation, i.e. the transposed matrix.

.. function:: mat:rotxyz(ax, ay, az, inv_)
              mat:rotxzy(ax, az, ay, inv_)
              mat:rotyxz(ay, ax, az, inv_)
              mat:rotyzx(ay, az, ax, inv_)
              mat:rotzxy(az, ax, ay, inv_)
              mat:rotzyx(az, ay, ax, inv_)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation of the first angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively *followed* by another 3D rotation of the second angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively of the frame rotated by the first rotation, and *followed* by a last 3D rotation of the third angle argument :var:`ax`, :var:`ay` or :var:`az` around the x-axis, y-axis or z-axis respectively of the frame already rotated by the two first rotations. If :var:`inv` is true, the returned matrix is the inverse rotations, i.e. the transposed matrix.

.. function:: mat:torotxyz(inv_)
              mat:torotxzy(inv_)
              mat:torotyxz(inv_)
              mat:torotyzx(inv_)
              mat:torotzxy(inv_)
              mat:torotzyx(inv_)

   Return three real :type:`number` representing the three angles :var:`ax`, :var:`ay` and :var:`az` (always in this order) of the 3D rotations stored in the real :type:`matrix` :var:`mat` :math:`[3\times 3]` by the methods with corresponding names. If :var:`inv` is true, the inverse rotations are returned, i.e. extracted from the transposed matrix.

.. function:: mat:rotv(v, av, inv_)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation of angle :var:`av` around the axis defined by the 3D vector-like :var:`v` (see `Axis-Angle representation <https://en.wikipedia.org/wiki/Axis–angle_representation>`_). If :var:`inv` is true, the returned matrix is the inverse rotation, i.e. the transposed matrix.

.. function:: mat:torotv(v_, inv_)

   Return a real :type:`number` representing the angle of the 3D rotation around the axis defined by a 3D vector as stored in the real :type:`matrix` :var:`mat` :math:`[3\times 3]` by the method :func:`mat:rotv()`. If the :type:`iterable`` :var:`v` is provided, it is filled with the components of the unit vector that defines the axis of the rotation.  If :var:`inv` is true, the inverse rotation is returned, i.e. extracted from the transposed matrix.

.. function:: mat:rotq(q, inv_)

   Return the real :type:`matrix` :var:`mat` :math:`[3\times 3]` filled with a 3D rotation defined by the quaternion :var:`q` (see `Axis-Angle representation <https://en.wikipedia.org/wiki/Axis–angle_representation>`_). If :var:`inv` is true, the returned matrix is the inverse rotation, i.e. the transposed matrix.

.. function:: mat:torotq(q_, inv_)

   Return a quaternion representing the 3D rotation as stored in the real :type:`matrix` :var:`mat` :math:`[3\times 3]` by the method :func:`mat:rotq()`. If the :type:`iterable`` :var:`q` is provided, it is filled with the components of the quaternion otherwise the quaternion is returned in a :type:`list` of length 4.  If :var:`inv` is true, the inverse rotation is returned, i.e. extracted from the transposed matrix.

Operators
=========

.. function:: #mat

   Return the size of the real, complex or integer matrix :var:`mat`, i.e. the number of elements interpreting the matrix as an array.

.. function:: mat[n]

   Return the value of the element at index :var:`n` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= n <= #mat`, i.e. interpreting the matrix as an array, :const:`nil` otherwise.

.. function:: mat[n] = v

   Assign the value :var:`v` to the element at index :var:`n` of the real, complex or integer matrix :var:`mat` for :expr:`1 <= n <= #mat`, i.e. interpreting the matrix as an array, and return the matrix, otherwise raise an *"out of bounds"* error.

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

