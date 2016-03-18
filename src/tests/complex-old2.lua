local lu = require 'luaunit'
local C  = require 'complex'

local pi  = math.pi
local nan = 0.0/0.0
local inf = 1.0/0.0
local eps = 1e-16

local isEqu    = lu.assertEquals
local isTrue   = lu.assertTrue
local isFalse  = lu.assertFalse
local isStr    = lu.assertIsString
local isEquEps = lu.assertAlmostEquals

TestComplex = {}

function TestComplex:setUp()
  -- nothing
end

function TestComplex:tearDown() -- moved close to setUp
  -- nothing
end

function TestComplex:testCtor()
  isEqu( C(1).re      ,    1 )
  isEqu( C(1).im      ,    0 )
  isEqu( C(0,1).re    ,    0 )
  isEqu( C(0,1).im    ,    1 )
  isEqu( C(-0.1,0).re , -0.1 )
  isEqu( C(-0.1,0).im ,    0 )
  isEqu( C(0,-0.1).re ,    0 )
  isEqu( C(0,-0.1).im , -0.1 )
end

function TestComplex:testReal() -- missing tests on pure real and/or imag
  local real = C.real
  isEqu( real(C(1,1))    ,    1 )
  isEqu( real(C(0,1))    ,    0 )
  isEqu( real(C(-0.1,1)) , -0.1 )
end

function TestComplex:testImag() -- missing tests on pure real and/or imag
  local imag = C.imag
  isEqu( imag(C(1,1))    ,    1 )
  isEqu( imag(C(1,0))    ,    0 )
  isEqu( imag(C(1,-0.1)) , -0.1 )
end

function TestComplex:testConj() -- missing tests on pure real and/or imag
  local conj = C.conj
  isEqu( conj(C(1,1))    , C(1,-1)    )
  isEqu( conj(C(0,1))    , C(0,-1)    )
  isEqu( conj(C(-0.1,1)) , C(-0.1,-1) )
end 

function TestComplex:testAbs() -- missing tests on pure real and/or imag
  local abs = C.abs
  isEqu( abs(C(0,-4)) , 4)
  isEqu( abs(C(0,0) ) , 0)
end

function TestComplex:testEq() -- missing tests on pure real and/or imag
  local equ = C.__eq -- suggest missing API entry, metamethod should not be used directly
  isTrue ( equ(C(1,1)   , C(1,1)   ) ) 
  isTrue ( equ(C(-0.1,1), C(-0.1,1)) )
  isFalse( equ(C(1,-1)  , C(2,1)   ) )  
end

function TestComplex:testUnm() -- missing tests on pure real and/or imag
  local unm = C.unm
  isEqu( unm(C(1,1)   ) , C(-1,-1)  ) 
  isEqu( unm(C(0,0)   ) , C(-0,-0)  )  
  isEqu( unm(C(-0.1,1)) , C(0.1,-1) ) 
end

function TestComplex:testAdd() -- missing tests on pure real and/or imag
  local add = C.add
  isEqu( add(C(1,1)   , C(1,1)) , C(2,2)   )
  isEqu( add(C(0,1)   , C(0,1)) , C(0,2)   )
  isEqu( add(C(1,0)   , C(1,0)) , C(2,0)   )
  isEqu( add(C(1,-1)  , C(2,1)) , C(3,0)   )
  isEqu( add(C(-0.1,1), C(1,1)) , C(0.9,2) )
end

function TestComplex:testSub() -- missing tests on pure real and/or imag
  local sub = C.sub
  isEqu( sub(C(1,1)   , C(1,1))   , C(0,0)    )
  isEqu( sub(C(0,1)   , C(0,1))   , C(0,0)    )
  isEqu( sub(C(1,0)   , C(1,0))   , C(0,0)    )
  isEqu( sub(C(1,-1)  , C(2,1))   , C(-1, -2) )
  isEqu( sub(C(-0.1,1), C(1,1))   , C(-1.1,0) )
end

function TestComplex:testMul() -- missing tests on pure real and/or imag
  local mul = C.mul
  isEqu( mul(C(1,1)   , C(1,1)) , C(0,2)      )
  isEqu( mul(C(1,0)   , C(1,0)) , C(1,0)      )
  isEqu( mul(C(1,1)   , C(2,0)) , C(2,2)      )
  isEqu( mul(C(0,1)   , C(0,1)) , C(-1,0)     )
  isEqu( mul(C(1,-1)  , C(2,1)) , C(3, -1)    )
  isEqu( mul(C(-0.1,1), C(1,1)) , C(-1.1,0.9) )
end

function TestComplex:testDiv() -- missing tests on pure real and/or imag + overflow rounding errors
  local div = C.div
  isEqu( div(C(0,1)   , C(0,1) )  , C(1,0)           )
  isEqu( div(C(6,3)   , C(7,-5))  , C(27/74, 51/74)  )
  isEqu( div(C(1,0)   , C(1,5) )  , C(1/26, (-5/26)) )
  isEqu( div(C(-0.1,1), C(1,1) )  , C(0.45,0.55)     )
  --isEqu(div(C(1,1)   , C(0,0))          , C(nan,nan)) 
end

function TestComplex:testTostring()
  local toStr = C.tostring
  isStr( toStr( C(1,1)    ) )
  isStr( toStr( C(-0.1,0) ) )
end
 
function TestComplex:testSqrt()
  local sqrt = C.sqrt
  isEqu( sqrt(C(0,0)), C(0,0) )
  isEqu( sqrt(C(0,2)), C(1,1) )
  isEqu( sqrt(C(1,0)), C(1,0) )

  local sqrtRe, sqrtIm = sqrt(C(2,2)).re, sqrt(C(2,2)).im
  isEquEps( sqrtRe, 1.55377397403003730, eps )   --only as an example
  isEquEps( sqrtIm, 0.64359425290558262, eps ) 
end

function TestComplex:testArg()   
  local arg = C.arg
  isEqu( arg(C(1,0) ) ,  0.0   * pi )
  isEqu( arg(C(-1,0)) ,  1.0   * pi )
  isEqu( arg(C(1,1) ) ,  1.0/4 * pi )
  isEqu( arg(C(0,1) ) ,  1.0/2 * pi )
  isEqu( arg(C(0,-1)) , -1.0/2 * pi ) 
end

function TestComplex:testLog() 
  local log = C.log
  isEqu( log(C(1,0) ) , C(0,0)           )
  isEqu( log(C(0,0) ) , C(-inf,0)        )
  isEqu( log(C(-1,0)) , C(0,pi)          )
  isEqu( log(C(0,1) ) , C(0, 1.0/2 * pi) )
  isEqu( log(C(0,-1)) , C(0,-1.0/2 * pi) )
end

os.exit( lu.LuaUnit.run() )