from sympy import *#simplify, diff, symbols, Derivative, Symbol, ordered, log, sqrt, exp,  #pip install sympy
from generateTPSAFile import createFile, createStringOfValues, closeFile, multiVarDataGen, multiVarDataGenSlow
import time

singleVarFunctions = {"inv": lambda x: 1/x, "invsqrt": lambda x: 1/sqrt(x), "sqrt": sqrt, "exp": exp, "log": log} 
#This is a suggestion, add and remove as you like, however the keys must be a string and the string must be function in mad

###########-------------------2 VARIABLE-------------------------######################
# if __name__ == "__main__":
multiVarFuncNames = ["axpbypc", "axypb"]
functionEnding = "fun"
filename = createFile(7, multiVarFuncNames, functionEnding)   
x, y = symbols("x y")
a = b = 1
exprList = [
    a * x + b * y,
    a * x * y + b
    ]
eps = [1, 6, 1, 1, 26, 1]
xyValues = [Rational(1, 100), Rational(1, 10), Integer(100), Integer(10), pi*Rational(5,3)]
multiVarDataGen(filename, singleVarFunctions, exprList,  eps, xyValues, 6, multiVarFuncNames, functionEnding) 
trigFunctions = {"sin": sin, "cos": cos, "tan": tan}
xyValues = [pi * Rational(1,3), pi * Rational(1,3), pi*Rational(2), pi, pi*Rational(5,3)]
multiVarDataGen(filename, trigFunctions, exprList,  [170, 115, 25], xyValues, 6, multiVarFuncNames, functionEnding) 
closeFile(filename)

# ###########-------------------3 VARIABLE-------------------------######################
# #HINT: eps = [axy + bz + c; axyz + b; ax + by + cz]

# filename = createFile("3VarTPSA.dat")
# with open(r"data/" + filename, "a") as file:
#     file.write(f"M.axypbzpcfun={{}}\nM.axyzpbfun={{}}\nM.axpbypczfun={{}}\n\n")

# for name, func in singleVarFunctions.items():
#     eps = [1, 1, 1]
#     match name:
#         case "exp":
#             eps = [18, 135, 9]
#         case "inv":
#             eps[1] = 4
#         case "log":
#             eps[1] = 8
#         case "sqrt":
#             eps[1] = 2
            
#     xyzFileWrite(name, func, eps)
# xyzFileWrite("sin", sin, [100,3700, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# xyzFileWrite("cos", cos, [160, 695, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# xyzFileWrite("tan", tan, [44 , 345,34], [pi * Rational(1,3), pi, pi*Rational(5,3)])

# closeFile(filename)

###########-------------------4 VARIABLE-------------------------######################
# # #HINT: eps = []

# filename = createFile("4VarTPSA.dat")
# with open(r"data/" + filename, "a") as file:
#     file.write(f"M.axypbvwpcfun={{}}\nM.avwxypbfun={{}}\nM.avpbwpcxpdyfun={{}}\n\n")

# for name, func in singleVarFunctions.items():
#     eps = [1, 1, 1]
#     match name:
#         case "exp":
#             eps = [28, 64, 16]
#         case "inv":
#             eps = [1, 2, 4]
#         case "invsqrt":
#             eps[1] = 4
#         case "log":
#             eps[1] = 4
#         case "sqrt":
#             eps[1] = 2
            
#     vwxyFileWrite(name, func, eps)
# xyzFileWrite("sin", sin, [100,3700, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# xyzFileWrite("cos", cos, [160, 695, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# xyzFileWrite("tan", tan, [44 , 345,34], [pi * Rational(1,3), pi, pi*Rational(5,3)])

# closeFile(filename)

###########-------------------5 VARIABLE-------------------------######################
#HINT: eps = []

# filename = createFile("5VarTPSA.dat")
# with open(r"data/" + filename, "a") as file:
#     file.write(f"M.avpbwpcxpdypezfun={{}}\nM.avwpbxypczfun={{}}\nM.avwxyzpbfun={{}}\n\n")

# for name, func in singleVarFunctions.items():
    # eps = [1, 1, 1]
    # match name:
    #     case "exp":
    #         eps = [22, 31, 132]
    #     case "inv":
    #         eps = [1, 2, 4]
    #     case "log":
    #         eps[2] = 8
    #     case "sqrt":
    #         eps[2] = 2
    #     case "invsqrt":
    #         eps = [1, 2, 10]
    # # if name != "inv" or name != "invsqrt" :
    # vwxyzFileWrite(name, func, eps)
# xyzFileWrite("sin", sin, [100,3700, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# xyzFileWrite("cos", cos, [160, 695, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# xyzFileWrite("tan", tan, [44 , 345,34], [pi * Rational(1,3), pi, pi*Rational(5,3)])
# closeFile(filename)
# closeFile("5VarTPSA.dat")

###########-------------------6 VARIABLE-------------------------######################
# filename = createFile("6VarTPSA.dat")
# with open(r"data/" + filename, "a") as file:
#     file.write(f"M.atpbvpcwpdxpeypfzfun={{}}\nM.atvpbwxpcyzfun={{}}\nM.atvwxyzpbfun={{}}\n\n")  

# for name, func in singleVarFunctions.items():
#     eps = [1, 1, 1]
#     match name:
#         case "exp":
#             eps = [20, 41, 20]
#         case "inv":
#             eps[2] = 2
#         case "log":
#             eps = [1, 5, 4]
#         case "sqrt":
#             eps[2] = 3
#         case "invsqrt":
#             eps[2] = 2

#     tvwxyzFileWrite(name, func, eps)
# # xyzFileWrite("sin", sin, [100,3700, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# # xyzFileWrite("cos", cos, [160, 695, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# # xyzFileWrite("tan", tan, [44 , 345,34], [pi * Rational(1,3), pi, pi*Rational(5,3)])

# closeFile(filename)