local M = { __author = 'ldeniau', __version = '2015.06', __help = {}, __test = {} }

-- MAD -------------------------------------------------------------------------

M.__help.self = [[
NAME
  matrix

SYNOPSIS
  local matrix = require 'matrix'

DESCRIPTION
  The module matrix implements the operators and math functions on matrices:
    (minus) -, +, -, *, /, %, ^, ==, #
    rows, cols, sizes,
    real, imag, conj, norm, angle,
    abs, arg, exp, log, pow, sqrt, proj,
    sin, cos, tan, sinh, cosh, tanh,
    asin, acos, atan, asinh, acosh, atanh,
    foldl, foldr, map, map2, maps, tostring.

RETURN VALUES
  The constructor of matrices

SEE ALSO
  math, complex, vector, generic
]]
 
-- DEFS ------------------------------------------------------------------------

local vector  = require 'vector'
local generic = require 'generic'
local tbl_new = require 'table.new'

local real, imag, conj, 
      abs, arg, exp, log, sqrt, proj,
      sin, cos, tan, sinh, cosh, tanh,
      asin, acos, atan, asinh, acosh, atanh,
      unm, add, sub, mul, div, mod, pow = 
      generic.real, generic.imag, generic.conj,
      generic.abs, generic.arg, generic.exp, generic.log, generic.sqrt, generic.proj,
      generic.sin, generic.cos, generic.tan, generic.sinh, generic.cosh, generic.tanh,
      generic.asin, generic.acos, generic.atan, generic.asinh, generic.acosh, generic.atanh,
      generic.unm, generic.add, generic.sub, generic.mul, generic.div, generic.mod, generic.pow

-- Lua API

local function ismat (x)
  return getmetatable(x) == M
end

local function matrix (rows, cols)
  if type(rows) == 'table' then
    local m = rows
    m[0] = #m
    for i=1,m[0] do
      m[i] = vector(m[i])
      assert(#m[i] == #m[1], "inhomegeneous sizes in table of tables")
    end
    return setmetatable(m, M)
  end

  if type(rows) == 'number' and rows > 0 and
     type(cols) == 'number' and cols > 0 then
    local m = tbl_new(rows,0)
    m[0] = rows
    for i=1,m[0] do m[i] = vector(cols) end
    return setmetatable(m, M)
  end

  error("invalid argument to matrix constructor, expecting (rows,cols) or table of tables")
end

function M.rows (x) return x[0] end
function M.cols (x) return #x[1] end
function M.sizes(x) return x:rows(), x:cols() end

function M.transpose (x)
  local rows, cols = x:sizes()
  local r = matrix(cols, rows)
  for i=1,rows do
    for j=1,cols do
      r[j][i] = x[i][j]
    end
  end
  return r
end

function M.norm (x)
  local r = vector(x:rows())
  for i=1,#r do r[i] = x[i]:norm() end
  return r
end

function M.foldl (x, r, f)
  for i=1,x:rows() do r[i] = x[i]:foldl(r[i], f) end
  return r
end

function M.foldr (x, r, f)
  for i=1,x:rows() do r[i] = x[i]:foldr(r[i], f) end
  return r
end

function M.map (x, f)
  local rows, cols = x:sizes()
  local r = matrix(rows, cols)
  for i=1,rows do
    for j=1,cols do
      r[i][j] = f(x[i][j])
    end
  end
  return r
end

function M.map2 (x, y, f)
  local rows, cols = x:sizes()

  if rows ~= y:rows() or cols ~= y:cols() then
    error("incompatible matrices sizes")
  end

  local r = matrix(rows, cols)
  for i=1,rows do
    for j=1,cols do
      r[i][j] = f(x[i][j], y[i][j])
    end
  end
  return r
end

function M.maps (x, y, f)
  if not ismat(x) then
    local rows, cols = y:sizes()
    local r = matrix(rows, cols)
    for i=1,rows do
      for j=1,cols do
        r[i][j] = f(x, y[i][j])
      end
    end
    return r
  end

  if not ismat(y) then
    local rows, cols = x:sizes()
    local r = matrix(rows, cols)
    for i=1,rows do
      for j=1,cols do
        r[i][j] = f(x[i][j], y)
      end
    end
    return r
  end

  return M.map2(x, y, f)
end

function M.real  (x)   return M.map (x,   real)  end
function M.imag  (x)   return M.map (x,   imag)  end
function M.conj  (x)   return M.map (x,   conj)  end
function M.abs   (x)   return M.map (x,   abs)   end
function M.arg   (x)   return M.map (x,   arg)   end
function M.exp   (x)   return M.map (x,   exp)   end
function M.log   (x)   return M.map (x,   log)   end
function M.sqrt  (x)   return M.map (x,   sqrt)  end
function M.proj  (x)   return M.map (x,   proj)  end
function M.sin   (x)   return M.map (x,   sin)   end
function M.cos   (x)   return M.map (x,   cos)   end
function M.tan   (x)   return M.map (x,   tan)   end
function M.sinh  (x)   return M.map (x,   sinh)  end
function M.cosh  (x)   return M.map (x,   cosh)  end
function M.tanh  (x)   return M.map (x,   tanh)  end
function M.asin  (x)   return M.map (x,   asin)  end
function M.acos  (x)   return M.map (x,   acos)  end
function M.atan  (x)   return M.map (x,   atan)  end
function M.asinh (x)   return M.map (x,   asinh) end
function M.acosh (x)   return M.map (x,   acosh) end
function M.atanh (x)   return M.map (x,   atanh) end
function M.__unm (x)   return M.map (x,   unm)   end
function M.__add (x,y) return M.maps(x,y, add)   end
function M.__sub (x,y) return M.maps(x,y, sub)   end
function M.__mul (x,y) return M.maps(x,y, mul)   end
function M.__div (x,y) return M.maps(x,y, div)   end
function M.__mod (x,y) return M.maps(x,y, mod)   end
function M.__pow (x,y) return M.maps(x,y, pow)   end

function M.__eq (x, y)
  local rows, cols = x:sizes()
  if rows ~= y:rows() or cols ~= y:cols() then return false end
  for i=1,rows do
    for j=1,cols do
      if x[i][j] ~= y[i][j] then return false end
    end
  end
  return true
end

function M.__tostring (x, sep, lsep)
  local rows, cols = x:sizes()
  local r = tbl_new(rows, 0)
  for i=1,rows do
    local c = tbl_new(cols, 0)
    for j=1,cols do
      c[j] = tostring(x[i][j])
    end
    r[i] = table.concat(c, sep or ' ')
  end
  return table.concat(r, lsep or '\n')
end

M.pow      = M.__pow
M.tostring = M.__tostring
M.__len    = function (x) return x:rows() * x:cols() end
M.__index  = M

-- END -------------------------------------------------------------------------
return matrix
