from taylorSeries import taylorSeries
from sympy import *
def createStringOfValues(singleVarFunc, multiVarFunc, varVals, order, eps):
    outputString = '{'
    _, values = taylorSeries(singleVarFunc(multiVarFunc), order, varVals) #For 2 thats up to 66 coefficients
    vars = list(ordered(multiVarFunc.atoms(Symbol)))
    if values["include"]:
        for i in range(len(vars)):
            outputString += f'{vars[i]}0 = {varVals[i]},'
        outputString += f"eps={eps}*eps,\n"
        
        valueList = list(values.values())
        for i in range(len(valueList) - 1):
            outputString += f'{valueList[i + 1]},'
            if (i+1) % 100 == 0:
                outputString += "\n"
        outputString += "\n},\n"
        return outputString
    else:
        return "\n"

def xyFileWrite(fName: str, f, eps, xyValues = [Rational(1, 100), Rational(1, 10), Integer(100), Integer(10), pi*Rational(5,3)]) -> None:
        axpbyStr = f"M.axpbypcfun.{fName} = {{\n"
        axypbStr = f"M.axypbfun.{fName} = {{\n"
        x, y = symbols("x y")
        a = b = 1
        axpby = a * x + b * y
        axypb = a * x * y + b
        for xVal in xyValues:
            for yVal in xyValues:
                axpbyStr += createStringOfValues(f, axpby, [xVal,yVal], 6, eps)
                # if not (fName == "exp" and (xVal == 100 or yVal == 100)): #To avoid e200+
                axypbStr += createStringOfValues(f, axypb, [xVal, yVal], 6, eps)

        axpbyStr += "}\n\n"
        axypbStr += "}\n\n"
        with open(r"data/2VarTPSA.dat", "a") as file:
            file.write(axpbyStr)
            file.write(axypbStr)

def xyzFileWrite(fName: str, f, eps, xyzValues = [Rational(1, 10), Integer(6), pi*Rational(5,3)]) -> None:
    stringList = [  
        f"M.axypbzpcfun.{fName} = {{\n", 
        f"M.axyzpbfun.{fName} = {{\n", 
        f"M.axpbypczfun.{fName} = {{\n"
    ]
    x, y, z = symbols("x y z")
    a = b = c = 1
    exprList = [
        a * x * y + b * z + c, 
        a * x * y * z + b * c,
        a * x + b* y + c * z, 
    ]
    for xVal in xyzValues:
        for yVal in xyzValues:
            for zVal in xyzValues:
                for i in range(len(stringList)):
                    stringList[i] += createStringOfValues(f, exprList[i], [xVal, yVal, zVal], 6, eps[i])

    with open(r"data/3VarTPSA.dat", "a") as file:
        for i in range(len(stringList)):
            stringList[i] += "}\n\n"
            file.write(stringList[i])

def vwxyFileWrite(fName: str, f, eps, vwxyValues = [Rational(1, 10), Integer(6), pi*Rational(5,3)]) -> None:
    stringList = [  
        f"M.axypbvwpcfun.{fName} = {{\n",
        f"M.avwxypbfun.{fName} = {{\n",
        f"M.avpbwpcxpdyfun.{fName} = {{\n"
    ]
    v, w, x, y = symbols("v w x y")
    a = b = c = d = 1
    exprList = [
        a * v * w + b * x * y + c,
        a * v * w * x * y + b,
        a * v + b * w + c * x + d * y
    ]
    for vVal in vwxyValues:
        for wVal in vwxyValues:
            for xVal in vwxyValues:
                for yVal in vwxyValues:
                    for i in range(len(stringList)):
                        stringList[i] += createStringOfValues(f, exprList[i], [vVal, wVal, xVal, yVal], 5, eps[i])
    
    with open(r"data/4VarTPSA.dat", "a") as file:
        for i in range(len(stringList)):
            stringList[i] += "}\n\n"
            file.write(stringList[i])

def vwxyzFileWrite(fName: str, f, eps, vwxyzValues = [Rational(1, 10), Integer(6), pi*Rational(5,3)]) -> None:
    stringList = [  
        f"M.avpbwpcxpdypezfun.{fName} = {{\n",
        f"M.avwpbxypczfun.{fName} = {{\n",
        f"M.avwxyzpbfun.{fName} = {{\n"
    ]
    v, w, x, y, z = symbols("v w x y z")
    a = b = c = d = e = 1
    exprList = [
        a * v + b * w + c * x + d * y + e * z,
        a * v * w + b * x * y + c * z, 
        a * v * w * x * y * z + b
    ]
    for vVal in vwxyzValues:
        for wVal in vwxyzValues:
            for xVal in vwxyzValues:
                for yVal in vwxyzValues:
                    for zVal in vwxyzValues:
                        for i in range(len(stringList)):
                            stringList[i] += createStringOfValues(f, exprList[i], [vVal, wVal, xVal, yVal, zVal], 5, eps[i])
    
    with open(r"data/5VarTPSA.dat", "a") as file:
        for i in range(len(stringList)):
            stringList[i] += "}\n\n"
            file.write(stringList[i])

def tvwxyzFileWrite(fName: str, func, eps, tvwxyzValues = [Rational(1, 10), Integer(6)]):#, pi*Rational(5,3)]) -> None:
    stringList = [  
        f"M.atpbvpcwpdxpeypfzfun.{fName} = {{\n",
        f"M.atvpbwxpcyzfun.{fName} = {{\n",
        f"M.atvwxyzpbfun.{fName} = {{\n"
    ]
    t, v, w, x, y, z = symbols("t v w x y z")
    a = b = c = d = e = f = 1
    exprList = [
        a * t + b * v + c * w + d * x + e * y + f * z,
        a * t * v + b * w * x + c * y * z, 
        a * t * v * w * x * y * z + b
    ]
    for tVal in tvwxyzValues:
        for vVal in tvwxyzValues:
            for wVal in tvwxyzValues:
                for xVal in tvwxyzValues:
                    for yVal in tvwxyzValues:
                        for zVal in tvwxyzValues:
                            for i in range(len(stringList)):
                                stringList[i] += createStringOfValues(func, exprList[i], [tVal, vVal, wVal, xVal, yVal, zVal], 5, eps[i])
    
    with open(r"data/6VarTPSA.dat", "a") as file:
        for i in range(len(stringList)):
            stringList[i] += "}\n\n"
            file.write(stringList[i])