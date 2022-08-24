from sympy import *#simplify, diff, symbols, Derivative, Symbol, ordered, log, sqrt, exp,  #pip install sympy
from generateTPSAFile import createFile, createStringOfValues, closeFile, multiVarDataGen, multiVarDataGenSlow
import time
if __name__ == '__main__':
    singleVarFunctions = {"inv": lambda x: 1/x, "invsqrt": lambda x: 1/sqrt(x), "sqrt": sqrt, "exp": exp, "log": log} 
    #This is a suggestion, add and remove as you like, however the keys must be a string and the string must be function in MAD

    ###########-------------------2 VARIABLE-------------------------######################
    multiVarFuncNames = ["axpbypc", "axypb"]
    functionEnding = "fun"
    filename = createFile(2, multiVarFuncNames, functionEnding)   
    x, y = symbols("x y")
    a = b = 1
    exprList = [
        a * x + b * y,
        a * x * y + b
        ]
    eps = [[1, 1], [6, 6], [1, 1], [26, 26], [1, 1]]
    trigVals = [Rational(1, 100), Rational(1, 10), Integer(100), Integer(10), pi*Rational(5,3)]
    multiVarDataGen(filename, singleVarFunctions, exprList,  eps, trigVals, 6, multiVarFuncNames, functionEnding) 

    #Trig functions
    trigFunctions = {"sin": sin, "cos": cos, "tan": tan}
    trigVals = [pi * Rational(1,3), pi * Rational(1,3), pi*Rational(2), pi, pi*Rational(5,3)]
    multiVarDataGen(filename, trigFunctions, exprList,  [[170, 170], [115, 115], [25, 25]], trigVals, 6, multiVarFuncNames, functionEnding) 
    closeFile(filename)

    # ###########-------------------3 VARIABLE-------------------------######################
    # #HINT: eps = [axy + bz + c; axyz + b; ax + by + cz]
    multiVarFuncNames = ["axypbzpc", "axyzpb", "axpbypcz"]
    filename = createFile(3, multiVarFuncNames, functionEnding)

    x, y, z = symbols("x y z")
    a = b = c = 1
    exprList = [
        a * x * y + b * z + c, 
        a * x * y * z + b * c,
        a * x + b* y + c * z, 
    ]
    eps = [[1, 4, 1], [1, 1, 1], [1, 2, 1], [18, 135, 9], [1, 8, 1]]
    values3VarPlus = [Rational(1, 10), Integer(6), pi*Rational(5,3)]
    multiVarDataGen(filename, singleVarFunctions, exprList,  eps, values3VarPlus, 5, multiVarFuncNames, functionEnding)

    #Trig functions
    trigVals = [pi * Rational(1,3), pi, pi*Rational(5,3)]
    multiVarDataGen(filename, trigFunctions, exprList,  [[100, 3700, 8], [160, 695, 8], [44, 345, 34]], trigVals, 6, multiVarFuncNames, functionEnding) 
    closeFile(filename)

    ###########-------------------4 VARIABLE-------------------------######################
    multiVarFuncNames = ["axypbvwpc", "avwxypb", "avpbwpcxpdy"]
    filename = createFile(4, multiVarFuncNames, functionEnding)

    v, w, x, y = symbols("v w x y")
    a = b = c = d = 1
    exprList = [
        a * v * w + b * x * y + c,
        a * v * w * x * y + b,
        a * v + b * w + c * x + d * y
    ]
    eps = [[1, 2, 4], [1, 4, 1], [1, 2, 1], [28, 64, 16], [1, 4, 1]]
    multiVarDataGen(filename, singleVarFunctions, exprList,  eps, values3VarPlus, 5, multiVarFuncNames, functionEnding)

    #Trig functions -Not done yet
    # multiVarDataGen(filename, trigFunctions, exprList,  [[1, 1, 1], [1, 1, 1], [1, 1, 1]], trigVals, 6, multiVarFuncNames, functionEnding) 
    closeFile(filename)

    ###########-------------------5 VARIABLE-------------------------######################
    multiVarFuncNames = ["avpbwpcxpdypez", "avwpbxypcz", "avwxyzpb"]
    filename = createFile(5, multiVarFuncNames, functionEnding)

    v, w, x, y, z = symbols("v w x y z")
    a = b = c = d = e = 1
    exprList = [
        a * v + b * w + c * x + d * y + e * z,
        a * v * w + b * x * y + c * z, 
        a * v * w * x * y * z + b
    ]
    eps = [[1, 2, 4], [1, 2, 10], [1, 1, 2], [22, 31, 132], [1, 1, 8]]
    multiVarDataGen(filename, singleVarFunctions, exprList,  eps, values3VarPlus, 5, multiVarFuncNames, functionEnding)

    #Trig functions -Not done yet
    # multiVarDataGen(filename, trigFunctions, exprList,  [[1, 1, 1], [1, 1, 1], [1, 1, 1]], trigVals, 6, multiVarFuncNames, functionEnding) 
    closeFile(filename)

    ###########-------------------6 VARIABLE-------------------------######################
    multiVarFuncNames = ["avpbwpcxpdypez", "avwpbxypcz", "avwxyzpb"]
    filename = createFile(6, multiVarFuncNames, functionEnding)

    t, v, w, x, y, z = symbols("t v w x y z")
    a = b = c = d = e = f = 1
    exprList = [
        a * t + b * v + c * w + d * x + e * y + f * z,
        a * t * v + b * w * x + c * y * z, 
        a * t * v * w * x * y * z + b
    ]
    tvwxyzValues = [Rational(1, 10), Integer(6)] #Three may take too long
    eps = [[1, 1, 2], [1, 1, 2], [1, 1, 3], [20, 41, 20], [1, 5, 4]]
    multiVarDataGen(filename, singleVarFunctions, exprList,  eps, tvwxyzValues, 5, multiVarFuncNames, functionEnding)

    #Trig functions -Not done yet
    # multiVarDataGen(filename, trigFunctions, exprList,  [[1, 1, 1], [1, 1, 1], [1, 1, 1]], trigVals, 6, multiVarFuncNames, functionEnding) 
    closeFile(filename)