local gmath   = require 'gmath'
local matrix  = require 'matrix'
local cmatrix = require 'cmatrix'
local jit     = require('jit')
local jitv    = require('jit.v')
-- jit.opt.start("sizemcode=256","loopunroll=25")
jitv.start()

--

local tostring, min = gmath.tostring, gmath.min

gmath.format = "%18.15f"
gmath.format = "%+9.5f"

local x, y, r

-- 1D FFT/IFFT real -> complex -> real

x = matrix {{ 16,  2,  3, 13, 5, 11, 10,  8,  9,  7,  6, 12, 4, 14, 15,  1 }} :t()
io.write('x=\n', x:tostring(), '\n')

y = x:fft()
io.write('y=\n', y:tostring(), '\n')

r = y:ifft(matrix(x:sizes()))
io.write('r=\n', r:tostring(), '\n')

print("error1 =")
io.write((r-x):tostring(), '\n')

print("norm1 =")
print((r-x):norm())

x = x *1i
io.write('x=\n', x:tostring(), '\n')

y = x:fft()
io.write('y=\n', y:tostring(), '\n')

r = y:ifft(cmatrix(x:sizes()))
io.write('r=\n', r:tostring(), '\n')

print("error2 =")
io.write((r-x):tostring(), '\n')

print("norm2 =")
print((r-x):norm())

print("---------------------------")

x = matrix {
  {1, 2, 3}, 
  {3.5, 0.5, 8}, 
  {-1, 2, -3}, 
  {4, 9, 7}
}

io.write('x=\n', x:tostring(), '\n')

y = x:fft()
io.write('y=\n', y:tostring(), '\n')

r = y:ifft(matrix(x:sizes()))
io.write('r=\n', r:tostring(), '\n')

print("error3 =")
io.write((r-x):tostring(), '\n')

print("norm3 =")
print((r-x):norm())

x = x *1i
io.write('x=\n', x:tostring(), '\n')

y = x:fft()
io.write('y=\n', y:tostring(), '\n')

r = y:ifft(cmatrix(x:sizes()))
io.write('r=\n', r:tostring(), '\n')

print("error4 =")
io.write((r-x):tostring(), '\n')

print("norm4 =")
print((r-x):norm())
