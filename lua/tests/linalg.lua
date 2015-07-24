-- time luajit -jv -Oloopunroll=50 -e "package.path = './?.lua;./lua/?.lua;./lib/?.lua;' .. package.path" lua/tests/matrix.lua

local ffi     = require 'ffi'
local gmath   = require 'gmath'
local complex = require 'complex'
local matrix  = require 'matrix'
local cmatrix = require 'cmatrix'
local linalg  = require 'linalg'

local A = matrix {
  {  6.80, -2.11,  5.66,  5.97,  8.23 },
  { -6.05, -3.30,  5.36, -4.44,  1.08 },
  { -0.45,  2.58, -2.70,  0.27,  9.04 },
  {  8.32,  2.71,  4.35, -7.17,  2.14 },
  { -9.67, -5.14, -7.26,  6.08, -6.87 } }

local B = matrix {
  {  4.02, -1.56,  9.81 },   
  {  6.19,  4.00, -4.09 },   
  { -8.22, -8.67, -4.57 },   
  { -7.57,  1.75, -8.61 },   
  { -3.03,  2.86,  8.99 } }
local Bt = B:trans()

gmath.format = "%- 6.2f"

io.write('A=\n', tostring(A), '\n')
io.write('B=\n', tostring(B), '\n')
io.write('Bt=\n', tostring(Bt), '\n')

gmath.format = "%- 10.5f"

-- solving AX=B and XA=B
local info, rcond, ferr

X1, info, rcond, err = linalg.solve_AX_B_eps(A,B )
print('info=', info, 'rcond=', rcond)
io.write('X1=\n', tostring(X1), '\n')
io.write('Err=\n', tostring(err), '\n')

for i=1,1e6 do
X2, info, rcond, err = linalg.solve_XA_B_eps(A,Bt)
end

print('info=', info, 'rcond=', rcond)
io.write('X2=\n', tostring(X2), '\n')
io.write('Err=\n', tostring(err), '\n')

-- end of solving AX=B and XA=B

--[=[ MAD
A=
 6.80  -2.11   5.66   5.97   8.23 
-6.05  -3.30   5.36  -4.44   1.08 
-0.45   2.58  -2.70   0.27   9.04 
 8.32   2.71   4.35  -7.17   2.14 
-9.67  -5.14  -7.26   6.08  -6.87 
B=
 4.02  -1.56   9.81 
 6.19   4.00  -4.09 
-8.22  -8.67  -4.57 
-7.57   1.75  -8.61 
-3.03   2.86   8.99 
X1= (solve AX=B)
-1.46576    0.26841    0.23346  
 2.71274   -1.04172   -0.62699  
 2.89937    0.22400    0.25025  
 1.85694   -0.36079    1.30327  
-0.94597   -0.57073   -0.27915  
X2= (solve XA=B:t())
-0.80071   -0.69524    0.59391    1.32173    0.56576  
-0.38962   -0.55443    0.84223   -0.10380    0.10571  
 0.95546    0.22066    1.90064    5.35766    4.04060  
]=]

--[=[ Matlab
A=[  6.80, -2.11,  5.66,  5.97,  8.23 ;
    -6.05, -3.30,  5.36, -4.44,  1.08 ;
    -0.45,  2.58, -2.70,  0.27,  9.04 ;
     8.32,  2.71,  4.35, -7.17,  2.14 ;
    -9.67, -5.14, -7.26,  6.08, -6.87 ]
B=[  4.02, -1.56,  9.81 ;   
     6.19,  4.00, -4.09 ;   
    -8.22, -8.67, -4.57 ;   
    -7.57,  1.75, -8.61 ;   
    -3.03,  2.86,  8.99 ]
A\B= (AX=B)
  -1.46576   0.26841   0.23346
   2.71274  -1.04172  -0.62699
   2.89937   0.22400   0.25025
   1.85694  -0.36079   1.30327
  -0.94597  -0.57073  -0.27915
B'/A= (XA=B')
  -0.80071  -0.69524   0.59391   1.32173   0.56576
  -0.38962  -0.55443   0.84223  -0.10380   0.10571
   0.95546   0.22066   1.90064   5.35766   4.04060
]=]
