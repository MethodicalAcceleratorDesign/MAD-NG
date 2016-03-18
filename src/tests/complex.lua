local lu  = require 'luaunit'
local C   = require 'complex'
local ref = require 'tests.complex_ref'

local isEqu    = lu.assertEquals
local isEquEps = lu.assertAlmostEquals

local pi  = math.pi
local nan = 0.0/0.0
local inf = 1.0/0.0

local eps = 1e-16

TestComplex = {}

function TestComplex:setUp()    end
function TestComplex:tearDown() end

function TestComplex:testExp()
  local exp = C.exp
  local z1  = ref.z1  
  local val, res

  for i=1,#z1 do 
    val = ref.expRef[i]
    res = exp( C(z1[i].re, z1[i].im) )

    isEquEps( res.re, val.re, eps )
    isEquEps( res.im, val.im, eps )
  end
end

--[[
function TestComplex:testSin()
  local sin = C.sin
  local exp = C.exp

  local z1     = ref.z1
  local sinRef = ref.sinRef
  local sinRes = {}
  local sinExp = {}

  for i = 1, #z1 do 
    sinRes[i] = sin( C( z1[i][1], z1[i][2]))

    isEquEps( sinRes[i].im, sinRef[i][2], eps )
    isEquEps( sinRes[i].re, sinRef[i][1], eps )
    
  --sin(C(x,y))=1/2*C(0,-1)*(exp(C(0,1)*C(x,y))-exp(C(0,1)*C(x,y)))
    sinExp[i] =1/2*C(0,1)*(exp(C(0,-1)*C( z1[i][1], z1[i][2]))-exp(C(0,1)*C( z1[i][1], z1[i][2])))

    isEquEps( sinRes[i].im, sinExp[i].im, eps )
    isEquEps( sinRes[i].re, sinExp[i].re, eps )
  end
end

function TestComplex:testCos()
  local cos = C.cos
  local exp = C.exp
  local z1  = ref.z1

  local cosRef = ref.cosRef
  local cosRes = {}

  --cos(C(x,y))=1/2*(exp(C(0,1)*C(x,y))+exp(C(0,-1)*C(x,y)))
  a=1/2*(exp(C(0,1)*C(1,1))+exp(C(0,-1)*C(1,1)))
  b=cos(C(1,1))

  isEquEps(a.re,b.re,eps)
  isEquEps(a.im,b.im,eps)
end

function TestComplex:testTan() -- to check sin/cos we need to check exp function
  local tan = C.tan
  local exp = C.exp

--tan(C(x,y))=(exp(C(0,1)*C(x,y))-1)/(C(0,1)*(exp(C(0,1)*C(x,y))+1))
  a = (exp(C(0,2)*C(1,1))-1)/(C(0,1)*(exp(C(0,2)*C(1,1))+1))
  b = tan(C(1,1))

  isEquEps(a.re,b.re,eps)
  isEquEps(a.im,b.im,eps)
  isEqu( tan( C(0,0))     ,  C(0,0) )
end

function TestComplex:testSinh() -- to check sin/cos we need to check exp function
  local sinh = C.sinh
  local exp = C.exp
  
  --sinh(C(x,y))=1/2*(exp(C(x,y))-exp((-1)*C(x,y)))
  a = 1/2*(exp(C(1,1))-exp((-1)*C(1,1)))
  b = sinh(C(1,1))
  
  isEquEps(a.re,b.re,eps)  
  isEquEps(a.im,b.im,eps)
  isEqu( sinh( C(0,0))     ,  C(0,0) )
end

function TestComplex:testCosh() -- to check sin/cos we need to check exp function
  local cosh = C.cosh
  local exp = C.exp

  --cosh(C(x,y))=1/2*(exp(C(x,y))+exp((-1)*C(x,y)))
  a = 1/2*(exp(C(1,1))+exp((-1)*C(1,1)))
  b = cosh(C(1,1))

  isEquEps(a.re,b.re,eps)    
  isEquEps(a.im,b.im,eps)
  isEqu( cosh( C(0,0))     ,  C(1,0) )
end

function TestComplex:testTanh() -- to check sin/cos we need to check exp function
  local tanh = C.tanh
  local exp = C.exp
  
  --tanh(C(x,y))=(1/2*(exp(C(x,y))-exp((-1)*C(x,y))))/(1/2*(exp(C(x,y))+exp((-1)*C(x,y))))
  a=(1/2*(exp(C(1,1))-exp((-1)*C(1,1))))/(1/2*(exp(C(1,1))+exp((-1)*C(1,1))))
  b = tanh(C(1,1))

  isEquEps(a.re,b.re,eps) 
  isEquEps(a.im,b.im,eps)
  isEqu( tanh( C(0,0))     ,  C(0,0) )
end
--]]

os.exit(lu.LuaUnit.run() )
