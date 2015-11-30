local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  linalg -- linear algebra solvers

SYNOPSIS
  local linalg = require 'linalg'
  local X1 = linalg.solve_AX_B(A,B)
  local X2 = linalg.solve_XA_B(A,B:t())
  local X3 = linalg.inverse(A)

DESCRIPTION
  The module linalg provides consistent front-end to lapack solvers:
  solve_AX_B, solve_XA_B. The third argument X (dest) is optional.

REMARK:
  By default, use_expert_drivers is false (faster).
  By default, check_info is false (user's reponsibility).

RETURN VALUES
  The module as table of solvers.

SEE ALSO
  math, gmath, complex, matrix, cmatrix, linalg
  http://www.netlib.org/lapack/lug (Lapack User's Guide)
]]

-- DEFS ------------------------------------------------------------------------

local ffi    = require 'ffi'
local lapack = require 'lapack'
local linbas = require 'linbas'
local gmath  = require 'gmath'

-- local
local  matrix = linbas. matrix
local cmatrix = linbas.cmatrix

local isnum, iscpx, iscal, ismat, iscmat, isamat =
  gmath.is_number, gmath.is_complex, gmath.is_scalar,
  gmath.is_matrix, gmath.is_cmatrix, gmath.isa_matrix

local min, max = gmath.min, gmath.max

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
local nrhs_   = int_ct()
local info_   = int_ct()
local lwork_  = int_ct()
local rcond_  = dbl_ct()
local rsize_  = dbl_ct()
local csize_  = cpx_ct()

-- solvers API

-- Crosscheck with algos used by:
-- http://docs.scipy.org/doc/scipy/reference/linalg.html
-- http://github.com/numpy/numpy/blob/master/numpy/linalg/linalg.py
-- http://julia.readthedocs.org/en/latest/stdlib/linalg/

local function inverse_nn (A)
  local nr, nc = A:sizes()
  assert(nr == nc, "invalid matrix sizes")
  local ipiv = iarr_ct(nr)
  n_[0] = nr

  if ismat(A) then
    lapack.dgetrf_(n_, n_, A.data, n_, ipiv, info_);
    if info_[0] == 0 then
      lwork_[0], n_[0] = -1, nr -- retrieve optimal size
      lapack.dgetri_(n_, A.data, n_, ipiv, rsize_, lwork_, info_);
      lwork_[0] = rsize_[0]
      local rwork = darr_ct(lwork_[0])
      lapack.dgetri_(n_, A.data, n_, ipiv, rwork , lwork_, info_);
    end
  elseif iscmat(A) then
    lapack.zgetrf_(n_, n_, A.data, n_, ipiv, info_);
    if info_[0] == 0 then
      lwork_[0], n_[0] = -1, nr -- retrieve optimal size
      lapack.zgetri_(n_, A.data, n_, ipiv, csize_, lwork_, info_);
      lwork_[0] = csize_[0].re
      local cwork = carr_ct(lwork_[0])
      lapack.zgetri_(n_, A.data, n_, ipiv, cwork , lwork_, info_);
    end
  else
    error("invalid input arguments")
  end

  if M.check_info then
        if info_[0] < 0 then error("invalid input matrix")
    elseif info_[0] > 0 then error("singular matrix detected") end
  end

  return A, info_[0] -- C -> F90 + transpose
end

local function solve_mn_XA_B (A, B)
  local nrhs, nr, nc = B:rows(), A:sizes()
  assert(nr == B:cols(), "incompatible matrix sizes")

  n_[0], nrhs_[0], rcond_[0] = nr, nrhs, -1
  local S = darr_ct(min(nr,nc))

  if ismat(A) and ismat(B) then
    lwork_[0] = -1 -- retrieve optimal size
    lapack.dgelss_(m_, n_, nrhs_, A.data, n_, B.data, n_, S, rcond_, rank_, rsize_, lwork_, info_)
    lwork_[0] = rsize_[0]
    local rwork = darr_ct(lwork_[0])
    print('lwork_=', lwork_[0]) -- TODO!!!
    lapack.dgelss_(m_, n_, nrhs_, A.data, n_, B.data, n_, S, rcond_, rank_, rwork , lwork_, info_)
  elseif iscmat(A) and iscmat(B) then
    lwork_[0] = -1 -- retrieve optimal size
    lapack.zgelss_(m_, n_, nrhs_, A.data, n_, B.data, n_, S, rcond_, rank_, csize_, lwork_, rsize_, info_)
    lwork_[0] = csize_[0].re
    local cwork, rwork = carr_ct(lwork_[0]), darr_ct(5*min(nr,nc))
    lapack.zgelss_(m_, n_, nrhs_, A.data, n_, B.data, n_, S, rcond_, rank_, cwork , lwork_, rwork , info_)
  else
    error("invalid input arguments")
  end

  if M.check_info then
        if info_[0] < 0 then error("invalid input matrix")
    elseif info_[0] > 0 then error("singular matrix detected") end
  end

  return B, info_[0] -- C -> F90 + transpose
end

local function solve_nn_XA_B (A, B)
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

  if M.check_info then
        if info_[0] < 0 then error("invalid input matrix")
    elseif info_[0] > 0 then error("singular matrix detected") end
  end

  return B, info_[0] -- C -> F90 + transpose
end

local function dsolve_nn_XA_B (A, B)
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
    error("incompatible input arguments")
  end

  if M.check_info then
        if info_[0] < 0 then error("invalid input matrix")
    elseif info_[0] > 0 then error("singular matrix detected") end
  end

  return X, info_[0], rcond_[0], ferr -- C -> F90 + transpose
end

-- Lua API

M.check_info         = false -- default
M.use_expert_drivers = false -- default

function M.inverse (A)
  local nr, nc = A:sizes()

  if nr == nc then
    local X, info = inverse_nn(A:trans())
    return X:trans(), info
  else
    -- check if algo proposed in the links is better.
    -- http://math.stackexchange.com/questions/75789/what-is-step-by-step-logic-of-pinv-pseudoinverse
    -- http://go.helms-net.de/math/divers/ACheapimplementationforthePenrose.htm
    local n = max(nr,nc)
    local B = ismat(A) and matrix(n,n) or cmatrix(n,n)
    return solve_mn_XA_B(A:copy(),B:ones())
  end
end

function M.solve_AX_B (A, B)
  -- XA = B -> X = B A^-1 = ((A^-1)^t B^t)^t
  if iscal(B) then
    local X, info = M.inverse(A)
    return X:mul(B,X), info
  end

  if M.use_expert_drivers ~= true then
    if A:rows() == A:cols() then
      local X, info = solve_nn_XA_B(A:trans(), B:trans()) -- square system
      return X:trans(), info
    else
      local X, info = solve_mn_XA_B(A:trans(), B:trans()) -- non-square system
      return X:trans(), info
    end
  else
    if A:rows() == A:cols() then
      local X, info, rcond, ferr = dsolve_nn_XA_B(A:trans(), B:trans())
      return X:trans(), info, rcond, ferr
    else
      error("NYI: solving non-square system using expert driver")
    end    
  end
end

function M.solve_XA_B (A, B)
  -- case A [n x n] or [m x n] and m < n or m > n
  -- case B is_scalar or isa_matrix
  -- case X is provided or not
  if iscal(B) then
    local X, info = M.inverse(A)
    return X:mul(B,X), info
  end

  if M.use_expert_drivers ~= true then
    if A:rows() == A:cols() then 
      return solve_nn_XA_B(A:copy(), B:copy()) -- square system
    else
      return solve_mn_XA_B(A:copy(), B:copy()) -- non-square system
    end
  else
    if A:rows() == A:cols() then
      return dsolve_nn_XA_B(A:copy(), B:copy())
    else
      error("NYI: solving non-square system using expert driver")
    end
  end
end

-- END -------------------------------------------------------------------------
return M
