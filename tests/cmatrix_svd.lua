local gmath   = require 'gmath'
local cmatrix = require 'cmatrix'
local jit     = require('jit')
local jitv    = require('jit.v')
-- jit.opt.start("sizemcode=256","loopunroll=25")
jitv.start()

--

local tostring, min = gmath.tostring, gmath.min

gmath.format = "%18.15f"
gmath.format = "%+9.5f"

local x, r, rs

-- SVD: square

x = cmatrix {
{ 16,  2,  3, 13 },
{  5, 11, 10,  8 },
{  9,  7,  6, 12 },
{  4, 14, 15,  1 }
} * 1i

io.write('x=\n', x:tostring(), '\n')

u, s, v = x:svd()
io.write('u=\n', u:tostring(), '\n')
io.write('s=\n', s:tostring(), '\n')
io.write('v=\n', v:tostring(), '\n')

n = min(x:sizes())
y = cmatrix(n,n):set_diag(s:t():totable()[1])
io.write('y=\n', y:tostring(), '\n')

r = x-u*y*v:t()
print("error1 =")
io.write(r:tostring(), '\n')

print("norm1 =")
print(r:norm())

-- SVD: non-square [4 x 2]

x = cmatrix {
  {1, 2}, 
  {3, 4}, 
  {5, 6}, 
  {7, 8}
} * 1i

io.write('x=\n', x:tostring(), '\n')

u, s, v = x:svd()
io.write('u=\n', u:tostring(), '\n')
io.write('s=\n', s:tostring(), '\n')
io.write('v=\n', v:tostring(), '\n')

y = cmatrix(x:sizes()):set_diag(s:t():totable()[1])
io.write('y=\n', y:tostring(), '\n')

r = x-u*y*v:t()
print("error2 =")
io.write(r:tostring(), '\n')

print("norm2 =")
print(r:norm())

-- SVD: non-square [3 x 4]

x = cmatrix {
  {1, 2, 3}, 
  {3.5, 0.5, 8}, 
  {-1, 2, -3}, 
  {4, 9, 7}
} : t() * 1i

io.write('x=\n', x:tostring(), '\n')

u, s, v = x:svd()
io.write('u=\n', u:tostring(), '\n')
io.write('s=\n', s:tostring(), '\n')
io.write('v=\n', v:tostring(), '\n')

y = cmatrix(x:sizes()):set_diag(s:t():totable()[1])
io.write('y=\n', y:tostring(), '\n')

r = x-u*y*v:t()
print("error3 =")
io.write(r:tostring(), '\n')

print("norm3 =")
print(r:norm())
