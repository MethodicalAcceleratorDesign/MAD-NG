.. index::
   Linear algebra
   Vector and matrix

**************
Linear Algebra
**************

This chapter describes the real :type:`matrix`, complex :type:`cmatrix` and integer :type:`imatrix` objects as supported by MAD-NG. The module for `Matrix <https://en.wikipedia.org/wiki/Matrix_(mathematics)>`_ is not exposed, only the contructors are visible from the :mod:`MAD` environment and thus, matrices are handled directly by their methods or by the generic functions of the same name from the module :mod:`MAD.gmath`. The :type:`imatrix`, i.e. matrix of integers, are mainly used for indexing other types of matrix with limited features. Column and row vectors are not separate types but only shortcuts for :math:`[n\times 1]` and :math:`[1\times n]` matrices. Note that :type:`matrix`, :type:`cmatrix` and :type:`imatrix` are all defined as C structures for direct compliance with the C API. 

Constructors
============

.. function::  vector(nrow)
              cvector(nrow)
              ivector(nrow)

   Return a real, complex or integer column vector (i.e. a matrix) of size :math:`[n_{\text{row}}\times 1]` filled with zeros.

.. function::  matrix(nrow, ncol_)
              cmatrix(nrow, ncol_)
              imatrix(nrow, ncol_)

   Return a real, complex or integer matrix of size :math:`[n_{\text{row}}\times n_{\text{col}}]` filled with zeros. Default: :expr:`ncol_ = rnow`.

.. function:: linspace(start_, stop, size_)

   Return a real or complex column vector of length :var:`size` filled with values equally spaced between :var:`start` and :var:`stop` on a linear scale. Default: :expr:`start_ = 0`, :expr:`size_ = 100`.

.. function:: logspace(start_, stop, size_)

   Return a real or complex column vector of length :var:`size` filled with values equally spaced between :var:`start` and :var:`stop` on a logarithmic scale. Default: :expr:`start_ = 1`, :expr:`size_ = 100`.

Functions
=========

.. function:: is_vector (a)
              is_cvector (a)
              is_ivector (a)

   Return :const:`true` if :var:`a` is respectively a real, complex or integer :type:`matrix` of size :math:`[n_{\text{row}}\times 1]` or :math:`[1\times n_{\text{row}}]`, :const:`false` otherwise. These functions are only available from the module :mod:`MAD.typeid`.

.. function:: isa_vector (a)

   Return :const:`true` if :var:`a` is a real or complex vector, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

.. function:: is_matrix (a)
              is_cmatrix (a)
              is_imatrix (a)

   Return :const:`true` if :var:`a` is respectively a real, complex or integer :type:`matrix`, :const:`false` otherwise. These functions are only available from the module :mod:`MAD.typeid`.

.. function:: isa_matrix (a)

   Return :const:`true` if :var:`a` is a real or complex matrix, :const:`false` otherwise. This function is only available from the module :mod:`MAD.typeid`.

Methods
=======

Operators
=========

In this section, :var:`num`, :var:`cpx`and :var:`idx` are generic names used for real, complex and integer numbers respectively, and :var:`mat`, :var:`cmat` and :var:`imat` are generic names used for real, complex and integer matrices respectively, unless otherwise stated.

.. function:: #mat

   Return the size of the real, complex or integer matrix :var:`mat` as computed by :func:`mat:size()`.

.. function:: mat[n]

   Return the value at index :var:`n` for :expr:`1 <= n <= #mat`, i.e. interpreting the real, complex or integer matrix :var:`mat` as an array, :const:`nil` otherwise.

.. function:: mat[n] = v

   Assign the value :var:`v` to index :var:`n` for :expr:`1 <= n <= #mat`, i.e. interpreting the real, complex or integer matrix :var:`mat` as an array, otherwise raise an *"out of bounds"* error.

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

