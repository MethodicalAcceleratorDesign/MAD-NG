local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  linalg -- linear algebra solvers

SYNOPSIS
  local linalg = require 'linalg'
  local X1 = linalg.solve_AX_B(A,B)
  local X2 = linalg.solve_XA_B(A,B:t())

DESCRIPTION
  The module linalg provides consistent front-end to lapack solvers:
  solve_AX_B, solve_XA_B

REMARK:
  By default, use_drivers is false (faster).

RETURN VALUES
  The module as table of solvers.

SEE ALSO
  math, gmath, complex, matrix, cmatrix, linalg
  http://www.netlib.org/lapack/lug (Lapack User's Guide)
]]

-- DEFS ------------------------------------------------------------------------

local ffi     = require 'ffi'
local lapack  = require 'lapack'
local matrix  = require 'matrix'
local cmatrix = require 'cmatrix'
local gmath   = require 'gmath'

-- FFI types

local int_ct  = ffi.typeof('int[1]')
local dbl_ct  = ffi.typeof('double[1]')
local cpx_ct  = ffi.typeof('complex[1]')

local iarr_ct = ffi.typeof('int[?]')
local darr_ct = ffi.typeof('double[?]')
local carr_ct = ffi.typeof('complex[?]')

-- local scalars
local n_      = int_ct()
local m_      = int_ct()
local p_      = int_ct()
local nrhs_   = int_ct()
local info_   = int_ct()
local rcond_  = dbl_ct()

local ismat, iscmat = gmath.is_matrix, gmath.is_cmatrix

-- solvers API

local function solve_nn_XA_B (A, B, err_)
  local nrhs, nr, nc = B:rows(), A:sizes() 
  assert(nr == nc and nr == B:cols(), "incompatible matrix sizes")

  local ipiv = iarr_ct(nr)

  n_[0], nrhs_[0] = nr, nrhs

  if ismat(A) and ismat(B) then
    lapack.dgesv_(n_, nrhs_, A.data, n_, ipiv, B.data, n_, info_)
  elseif iscmat(A) and iscmat(B) then
    lapack.zgesv_(n_, nrhs_, A.data, n_, ipiv, B.data, n_, info_)
  else
    error("invalid input arguments")
  end

  if err_ then
        if info_[0] < 0 then error("invalid input matrix")
    elseif info_[0] > 0 then error("singular matrix detected") end
  end

  return B, info_[0] -- C -> F90 + transpose
end

local function dsolve_nn_XA_B (A, B, err_)
  local nrhs, nr, nc = B:rows(), A:sizes()
  assert(nr == nc and nr == B:cols(), "incompatible matrix sizes")

  local trans = 'N'
  local ipiv  = iarr_ct(nr)

  local berr = darr_ct(nrhs)
  local C    = darr_ct(nr)
  local R    = darr_ct(nr)
  local ferr = matrix(nrhs)

  local AF, work, X, iwork, rwork
  if ismat(A) then
    iwork = iarr_ct(nr)
    AF    = darr_ct(nr*nr)
    work  = darr_ct(4*nr)
    X     = matrix(nrhs, nr)
  else
    rwork = darr_ct(2*nr)
    AF    = carr_ct(nr*nr)
    work  = carr_ct(2*nr)
    X     = cmatrix(nrhs, nr)
  end

  n_[0], nrhs_[0] = nr, nrhs

  if ismat(A) and ismat(B) then
    lapack.dgesvx_('E', 'N', n_, nrhs_, A.data, n_, AF, n_, ipiv, 'N', R, C,
                    B.data, n_, X.data, n_, rcond_, ferr.data, berr, work, iwork, info_)
  elseif iscmat(A) and iscmat(B) then
    lapack.zgesvx_('E', 'N', n_, nrhs_, A.data, n_, AF, n_, ipiv, 'N', R, C,
                    B.data, n_, X.data, n_, rcond_, ferr.data, berr, work, rwork, info_)
  else
    error("invalid input arguments")
  end

  if err_ then
        if info_[0] < 0 then error("invalid input matrix")
    elseif info_[0] > 0 then error("singular matrix detected") end
  end

  return X, info_[0], rcond_[0], ferr -- C -> F90 + transpose
end

-- Lua API

M.use_drivers = false -- default

function M.solve_AX_B (A, B, err_)
  -- XA = B -> X = B A^-1 = ((A^-1)^t B^t)^t
  assert(isamat(A) and isamat(B), "invalid argument to solver")

  if M.use_drivers ~= true then
    if A:rows() == A:cols() then
      local X, info = solve_nn_XA_B(A:trans(), B:trans(), err_)
      return X:trans(), info
    else
      error("NYI: solving non-square system")
    end
  else
    if A:rows() == A:cols() then
      local X, info, rcond, ferr = dsolve_nn_XA_B(A:trans(), B:trans(), err_)
      return X:trans(), info
    else
      error("NYI: solving non-square system")
    end    
  end

end

function M.solve_XA_B (A, B, err_)
  assert(isamat(A) and isamat(B), "invalid argument to solver")

  if M.use_drivers ~= true then
    if A:rows() == A:cols() then -- square system
      return solve_nn_XA_B(A:copy(), B:copy(), err_)
    else
      error("NYI: solving non-square system")
    end
  else
    if A:rows() == A:cols() then
      return dsolve_nn_XA_B(A:copy(), B:copy(), err_)
    else
      error("NYI: solving non-square system")
    end
  end
end

-- END -------------------------------------------------------------------------
return M
