local lu  = require 'luaunit'
local C   = require 'complex'
local ref = require 'tests.complex_ref'

local isEqu    = lu.assertEquals
local isTrue   = lu.assertTrue
local isFalse  = lu.assertFalse
local isStr    = lu.assertIsString
local isEquEps = lu.assertAlmostEquals

local pi  = math.pi
local nan = 0.0/0.0
local inf = 1.0/0.0
local eps = 1e-16

TestComplex = {}

function TestComplex:setUp()    end
function TestComplex:tearDown() end

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

function TestComplex:testReal() 
  local real = C.real
  isEqu( real(C(1,1))    ,    1 )
  isEqu( real(C(0,1))    ,    0 )
  isEqu( real(C(1,0))    ,    1 )
  isEqu( real(C(-0.1,1)) , -0.1 )
end

function TestComplex:testImag()
  local imag = C.imag
  isEqu( imag(C(1,1))    ,    1 )
  isEqu( imag(C(1,0))    ,    0 )
  isEqu( imag(C(0,1))    ,    1 )
  isEqu( imag(C(1,-0.1)) , -0.1 )
end

function TestComplex:testConj() 
  local conj = C.conj
  isEqu( conj(C(1,1))    , C(1,-1)    )
  isEqu( conj(C(0,1))    , C(0,-1)    )
  isEqu( conj(C(1,0))    , C(1,-0)    )
  isEqu( conj(C(-0.1,1)) , C(-0.1,-1) )
end 

function TestComplex:testAbs() 
  local abs = C.abs
  isEqu( abs(C(0,-4)) , 4 )
  isEqu( abs(C(0,0) ) , 0 )
  isEqu( abs(C(1,0) ) , 1 )
  isEqu( abs(C(0,1) ) , 1 )
end

function TestComplex:testArg()   
  local arg = C.arg
  isEqu( arg(C(1,0) ) ,  0.0   * pi )
  isEqu( arg(C(-1,0)) ,  1.0   * pi )
  isEqu( arg(C(1,1) ) ,  1.0/4 * pi )
  isEqu( arg(C(0,1) ) ,  1.0/2 * pi )
  isEqu( arg(C(0,-1)) , -1.0/2 * pi ) 
end

function TestComplex:testEq() 
  local equ = C.__eq -- suggest missing API entry, metamethod should not be used directly
  isTrue ( equ(C(1,1)   , C(1,1)   ) ) 
  isTrue ( equ(C(0,1)   , C(0,1)   ) ) 
  isTrue ( equ(C(1,0)   , C(1,0)   ) ) 
  isTrue ( equ(C(-0.1,1), C(-0.1,1)) )
  isFalse( equ(C(1,-1)  , C(2,1)   ) )  
end

function TestComplex:testUnm() 
  local unm = C.unm
  isEqu( unm(C(1,1)   ) , C(-1,-1)  ) 
  isEqu( unm(C(0,1)   ) , C(-0,-1)  )
  isEqu( unm(C(1,0)   ) , C(-1,-0)  )
  isEqu( unm(C(0,0)   ) , C(-0,-0)  )  
  isEqu( unm(C(-0.1,1)) , C(0.1,-1) ) 
end

function TestComplex:testAdd()
  local add = C.add
  isEqu( add(C(1,1)   , C(1,1)) , C(2,2)   )
  isEqu( add(C(0,1)   , C(0,1)) , C(0,2)   )
  isEqu( add(C(1,0)   , C(1,0)) , C(2,0)   )
  isEqu( add(C(1,-1)  , C(2,1)) , C(3,0)   )
  isEqu( add(C(-0.1,1), C(1,1)) , C(0.9,2) )
end

function TestComplex:testSub() 
  local sub = C.sub
  isEqu( sub(C(1,1)   , C(1,1))   , C(0,0)    )
  isEqu( sub(C(0,1)   , C(0,1))   , C(0,0)    )
  isEqu( sub(C(1,0)   , C(1,0))   , C(0,0)    )
  isEqu( sub(C(1,-1)  , C(2,1))   , C(-1, -2) )
  isEqu( sub(C(-0.1,1), C(1,1))   , C(-1.1,0) )
end

function TestComplex:testMul() 
  local mul = C.mul
  isEqu( mul(C(1,1)   , C(1,1)) , C(0,2)      )
  isEqu( mul(C(1,0)   , C(1,0)) , C(1,0)      )
  isEqu( mul(C(1,1)   , C(2,0)) , C(2,2)      )
  isEqu( mul(C(0,1)   , C(0,1)) , C(-1,0)     )
  isEqu( mul(C(1,-1)  , C(2,1)) , C(3, -1)    )
  isEqu( mul(C(-0.1,1), C(1,1)) , C(-1.1,0.9) )
end

function TestComplex:testDiv() -- overflow rounding errors
  local div = C.div
  isEqu( div(C(0,1)   , C(0,1) )  , C(1,0)                  )
  isEqu( div(C(6,3)   , C(7,-5))  , C(27/74, 51/74)         )
  isEqu( div(C(1,0)   , C(1,5) )  , C(1/26, (-5/26))        )
  isEqu( div(C(-0.1,1), C(1,1) )  , C(0.45,0.55)            )
  isEqu( div(C(0,1)   , C(0,10))  , C(0.1,0)                )
  isEqu( div(C(0,10)  , C(0,3))   , C(3.3333333333333333,0) )
  --isEqu( C(1,1):div(C(0,0)):tostring()   , C(nan,nan):tostring() ) 
end

function TestComplex:testTostring()
  local toStr = C.tostring
  isStr( toStr( C(1,1)    ) )
  isStr( toStr( C(-0.1,0) ) )
end
 
function TestComplex:testSqrt()
  local sqrt = C.sqrt
  isEqu( sqrt(C(0,0)) , C(0,0) )
  isEqu( sqrt(C(0,2)) , C(1,1) )
  isEqu( sqrt(C(1,0)) , C(1,0) )
end 

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
  isEqu( exp(C(0,0)), C(1,0) )  
end

function TestComplex:testLog() 
  local log = C.log
  isEqu( log(C(1,0) ) , C(0,0)           )
  isEqu( log(C(0,0) ) , C(-inf,0)        )
  isEqu( log(C(-1,0)) , C(0,pi)          )
  isEqu( log(C(0,1) ) , C(0, 1.0/2 * pi) )
  isEqu( log(C(0,-1)) , C(0,-1.0/2 * pi) )
end

function TestComplex:testProj()
  local proj = C.proj
  isEqu( proj( C(1,2))     ,  C(1,2)     )
  isEqu( proj( C(inf,-1))  ,  C(inf,-1)  )
  isEqu( proj( C(0, -inf)) ,  C(0 ,-inf) )
end 
z
function TestComplex:testSin()
  local sin = C.sin
  local z1  = ref.z1
  local val, res

  for i = 1, #z1 do 
    val = ref.sinRef[i]
    res = sin( C(z1[i].re, z1[i].im) )

    isEquEps( res.re, val.re, eps )
    isEquEps( res.im, val.im, eps )
  end
  isEqu( sin(C(0,0)), C(0,0) )
end

function TestComplex:testCos()
  local cos = C.cos  
  local z1  = ref.z1
  local val, res

  for i = 1, #z1 do 
    val = ref.cosRef[i]
    res = cos( C(z1[i].re, z1[i].im) )

    isEquEps( res.re, val.re, eps )
    isEquEps( res.im, val.im, eps )
  end
  isEqu( cos(C(0,0)), C(1,0) )
end

function TestComplex:testTan()
  local tan = C.tan
  local z1  = ref.z1
  local val, res

  for i = 1, #z1 do 
    val = ref.tanRef[i]
    res = tan( C(z1[i].re, z1[i].im) )

    isEquEps( res.re, val.re, eps )
    isEquEps( res.im, val.im, eps )
  end
  isEqu( tan( C(0,0)) , C(0,0) )
end

function TestComplex:testSinh() 
  local sinh = C.sinh
  local z1  = ref.z1
  local val, res

  for i = 1, #z1 do 
    val = ref.sinhRef[i]
    res = sinh( C(z1[i].re, z1[i].im) )

    isEquEps( res.re, val.re, eps )
    isEquEps( res.im, val.im, eps )
  end
  isEqu( sinh( C(0,0)) , C(0,0) )
end

function TestComplex:testCosh() 
  local cosh = C.cosh
  local z1  = ref.z1
  local val, res

  for i = 1, #z1 do 
    val = ref.coshRef[i]
    res = cosh( C(z1[i].re, z1[i].im) )

    isEquEps( res.re, val.re, eps )
    isEquEps( res.im, val.im, eps )
  end
  isEqu( cosh( C(0,0)) , C(1,0) )
end

function TestComplex:testTanh()
  local tanh = C.tanh
  local z1  = ref.z1
  local val, res

  for i = 1, #z1 do 
    val = ref.tanhRef[i]
    res = tanh( C(z1[i].re, z1[i].im) )

    isEquEps( res.re, val.re, eps )
    isEquEps( res.im, val.im, eps )
  end
  isEqu( tanh( C(0,0)), C(0,0) )
end


function TestComplex:testAsin()  
  local asin = C.asin
  isEqu( asin( C(0,0) ) , C(0,0)         )
  isEqu( asin( C(1,0) ) , C(1/2 * pi,0)  )
  isEqu( asin( C(-1,0)) , C(-1/2 * pi,0) ) 
end

function TestComplex:testAcos() 
  local acos = C.acos
  isEqu( acos( C(1,0) ) , C(0,0)        )
  isEqu( acos( C(-1,0)) , C(pi,0)       )
  isEqu( acos( C(0,0) ) , C(1/2 * pi,0) )
end

function TestComplex:testAtan()  
  local atan = C.atan
  isEqu( atan( C(0,0)    ) , C(0,0)       )
  isEqu( atan( C(0,1)    ) , C(0,inf)     )
  isEqu( atan( C(0,-1)   ) , C(0,-inf)    )
  isEqu( atan( C(inf,0)  ) , C(1/2*pi,0)  )
  isEqu( atan( C(-inf,0) ) , C(-1/2*pi,0) )
end

function TestComplex:testAsinh()  
  local asinh = C.asinh
    isEqu( asinh( C(0,0)) , C(0,0) )

end

function TestComplex:testAtanh()  
  local atanh = C.atanh
    isEqu( atanh( C(0,0) ) , C(0,0) )
end

--[[
function TestComplex:testAcosh() 
  local acosh = C.acosh
end

function TestComplex:testSqrt()
  local sqrt = C.sqrt
  local z1  = ref.z1
  local val, res

  for i = 1, #z1 do 
    val = ref.sqrtRef[i]
    res = sqrt( C(z1[i].re, z1[i].im) )

    isEquEps( res.re, val.re, eps )
    isEquEps( res.im, val.im, eps )
  end

  isEqu( sqrt(C(0,0)), C(0,0) )
  isEqu( sqrt(C(0,2)), C(1,1) )
  isEqu( sqrt(C(1,0)), C(1,0) )
end

function TestComplex:testPow() end
function TestComplex:testMod() edn

--]]
os.exit(lu.LuaUnit.run() )
