from sympy import *
import math
def createExponents(nv, maxOrder):
    """nv: number of variables, maxOrder (mo): max order of differentiation"""
    exponentList = [[0] * nv]
    for order in range(1, maxOrder + 1):
        exponents = [0] * nv
        i = 0
        j = 1
        while exponents[-1] != order:
            exponents[0] += order - sum(exponents)
            exponentList.append(exponents.copy())
            if exponents[i] == order and i < nv - 1:
                exponents[i] = 0
                exponents[i + 1] += 1
                i += 1
                j = 1
            elif i < nv:
                if exponents[-1] == order - 1 and sum(exponents[:-2]) == 0:
                    exponents = [0] * nv
                    exponents[-1] = order
                elif sum(exponents[j:]) == order:
                    while exponents[j] == 0:
                        j += 1
                    exponents[j+1] += 1
                    exponents[j] = 0
                    j = 1
                else:
                    exponents[j-1] -= 1
                    exponents[j] += 1
        exponentList.append(exponents.copy())
    return exponentList

def taylorSeries(expr, maxOrder, varVals):
    """expr: expression you want the taylor series of;
maxorder: order of taylor series
args: list of variable values to expand around"""
    exprResult = {}
    evaluatedResult = {"include": True}
    vars = list(ordered(expr.atoms(Symbol)))
    nv = len(vars)
    exponentList = createExponents(nv, maxOrder)
    subsList = []
    for var in range(nv):
        subsList.append((vars[var], varVals[var]))
    for exponents in exponentList:
        newExpr = expr
        for var in range(len(exponents)):
            newExpr = diff(newExpr, vars[var], exponents[var])
        combinations = 1
        for exp in exponents:
            combinations *= factorial(exp)
        newExpr = newExpr / combinations
        evaledExpr = N(newExpr.subs(subsList), 20)
        if math.isinf(float(evaledExpr)):
            evaluatedResult["include"] = False
        exprResult[str(exponents)] = newExpr
        evaluatedResult[str(exponents)] = evaledExpr
    return exprResult, evaluatedResult

def createStringOfValues(singleVarFunc, multiVarFunc, varVals, order, eps):
    """"Calls taylor series and generates string for MAD interpretation"""
    outputString = '{'
    _, values = taylorSeries(singleVarFunc(multiVarFunc), order, varVals) #For 2 thats up to 66 coefficients
    vars = list(ordered(multiVarFunc.atoms(Symbol))) #Need to be ordered?
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