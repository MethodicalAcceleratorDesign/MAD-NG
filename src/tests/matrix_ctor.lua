local ffi = require 'ffi'

ffi.cdef[[
//typedef struct { int32_t nr, nc; double data[?]; } matrix_t;

void* malloc(size_t);
void  free  (void*);
void* mad_malloc(size_t);
void  mad_free  (void*);
]]
--[[
local matrix_ctor = ffi.typeof( 'matrix_t')

local function matrix(nr, nc)
  local nc  = nc or 1
  local len = nr*nc
  if len < 36 then
    return matrix_ctor(len, nr, nc)
  end

  local siz = ffi.sizeof('matrix_t', len)
  local ptr = ffi.gc(ffi.C.malloc(siz), ffi.C.free)
  local mat = ffi.cast('matrix_t*', ptr)
  ffi.fill(mat, siz)
  mat.nr, mat.nc = nr, nc
  return mat[0]
end

local mt = {
  __add = function (a,b) io.write('__add(a,b)\n'); return b ; end,
  __sub = function (a,b) io.write('__sub(a,b)\n'); return b ; end,
  __mul = function (a,b) io.write('__mul(a,b)\n'); return b ; end,
  __div = function (a,b) io.write('__div(a,b)\n'); return b ; end,
}

ffi.metatype( 'matrix_t', mt)

local a, b

-- use ffi ctor
a = matrix(5,5)
b = matrix(5,5)

print(a,b)
print(a+b, a-b, a*b, a/b)

-- use malloc/free
a = matrix(10,10)
b = matrix(10,10)

print(a,b)
print(a+b, a-b, a*b, a/b)
]]
local clib = ffi.C

for i=1,1e8 do
  clib.mad_free(clib.mad_malloc(10))
end