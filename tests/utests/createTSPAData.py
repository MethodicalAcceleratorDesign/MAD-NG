from sympy import *#simplify, diff, symbols, Derivative, Symbol, ordered, log, sqrt, exp,  #pip install sympy
from generateTPSAFile import xyFileWrite, xyzFileWrite, vwxyFileWrite, vwxyzFileWrite, tvwxyzFileWrite

def createFile(filename: str = "nVarTPSA.dat"):
    with open(r"data/" + filename, "w") as file:
        file.write(f"local eps, nan, inf, pi\tin MAD.constant\nlocal M = {{}}\n")
    return filename

def closeFile(filename): #This will close any file
    with open(r"data/" + filename, "a") as file: #Could speed up process by passing file object?
        file.write("\nreturn M")

singleVarFunctions = {"inv": lambda x: 1/x, "invsqrt": lambda x: 1/sqrt(x), "sqrt": sqrt, "exp": exp, "log": log}

###########-------------------2 VARIABLE-------------------------######################
# filename = createFile("2VarTPSA.dat")
# with open(r"data/" + filename, "a") as file:
#     file.write(f"M.axpbypcfun={{}}\nM.axypbfun={{}}\n\n")
    
# for name, func in singleVarFunctions.items():
#     match name:
#         case "exp":
#             eps = 26
#         case "invsqrt":
#             eps = 6
#         case _:
#             eps = 1
#     xyFileWrite(name, func, eps)
# xyFileWrite("sin", sin, 170, [pi * Rational(1,2), pi * Rational(1,3), pi*Rational(2), pi, pi*Rational(5,3)])
# xyFileWrite("cos", cos, 115, [pi * Rational(1,2), pi * Rational(1,3), pi*Rational(2), pi, pi*Rational(5,3)])
# xyFileWrite("tan", tan, 25, [pi * Rational(1,3), pi * Rational(1,3), pi*Rational(2), pi, pi*Rational(5,3)])
# closeFile(filename)

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

filename = createFile("4VarTPSA.dat")
with open(r"data/" + filename, "a") as file:
    file.write(f"M.axypbvwpcfun={{}}\nM.avwxypbfun={{}}\nM.avpbwpcxpdyfun={{}}\n\n")

for name, func in singleVarFunctions.items():
    eps = [1, 1, 1]
    match name:
        case "exp":
            eps = [28, 64, 16]
        case "inv":
            eps = [1, 2, 4]
        case "invsqrt":
            eps[1] = 4
        case "log":
            eps[1] = 4
        case "sqrt":
            eps[1] = 2
            
    vwxyFileWrite(name, func, eps)
# xyzFileWrite("sin", sin, [100,3700, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# xyzFileWrite("cos", cos, [160, 695, 8], [pi * Rational(1,2), pi, pi*Rational(5,3)])
# xyzFileWrite("tan", tan, [44 , 345,34], [pi * Rational(1,3), pi, pi*Rational(5,3)])

closeFile(filename)

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