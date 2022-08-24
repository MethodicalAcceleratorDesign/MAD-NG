from taylorSeries import taylorSeries
from sympy import *
from itertools import permutations
from multiprocessing import Process, Manager


def createFile(numVar: int, multiVarFnams: list[str], end: str = "") -> str: #Typing is for docs and not enforced
    """numvar: number of variables in the TPSA data you are generating

multiVarFnams: The names of the multivariate functions, this will be the name of the function accessed from lua

end (optional): default (""), this is an addition to the end of the function names accessed from lua

returns filename of the file that has been created

Must be run from /MAD/tests/utests/data/tpsaTestResultGeneration folder
"""
    filename = str(numVar) + "VarTPSA.dat"
    with open(r"../" + filename, "w") as file:
        file.write(f"local eps, nan, inf, pi\tin MAD.constant\nlocal M = {{}}\n") #Do we need any other than eps, nan, inf and pi?
        for fnam in multiVarFnams:
            file.write(f"M.{fnam + end}={{}}\n")
        file.write("\n")
    return filename

def createStringOfValues(singleVarFunc, multiVarFunc, varVals, order, eps):
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

def closeFile(filename: str) -> None: #
    """Input the filename gained from creating the file and this will close off the file, allowing use with MAD
    
/MAD/tests/utests/data/tpsaTestResultGeneration"""
    with open(r"../" + filename, "a") as file: #Could speed up process by passing file object?
        file.write("\nreturn M")

def multiVarDataGenSlow(filename:str, singleVarFuncs: dict, multiVarExprs: list, epsList: list[int], expansionPoints: list, order: int, multiVarFnams: list[str], end: str = "") -> None:
    stringList = []
    keys = list(singleVarFuncs.keys())
    for multiVarFnam in multiVarFnams:
        for singleVarFnam in keys:
            stringList.append(f"M.{multiVarFnam + end}.{singleVarFnam} = {{\n")
    numVar = len(multiVarExprs[0].atoms(Symbol))
    #Generate permutations of values to loop through
    expansionPoints = [*expansionPoints] * numVar
    expansionPoints = list(set(permutations(expansionPoints, numVar)))

    for multiVarIndex in range(len(multiVarExprs)):
        for keyIndex in range(len(keys)): 
            for expansionPoint in expansionPoints:
                stringList[keyIndex + len(keys) * multiVarIndex] += createStringOfValues(singleVarFuncs[keys[keyIndex]], multiVarExprs[multiVarIndex], expansionPoint, order, epsList[keyIndex])
    
    with open(r"../" + filename, "a") as file:
        for i in range(len(stringList)):
            stringList[i] += "}\n\n"
            file.write(stringList[i])
def generateStringsPerMultiVarFunc(multiVarIndex, multiVarExprs, singleVarFuncs, keys, expansionPoints, order, epsList, returnList):
    stringList = [""] * len(keys)
    for keyIndex in range(len(keys)): 
        for expansionPoint in expansionPoints:
            stringList[keyIndex] += createStringOfValues(singleVarFuncs[keys[keyIndex]], multiVarExprs[multiVarIndex], expansionPoint, order, epsList[keyIndex][multiVarIndex])
    returnList[multiVarIndex] = stringList

def multiVarDataGen(filename:str, singleVarFuncs: dict, multiVarExprs: list, epsList: list[int], expansionPoints: list, order: int, multiVarFnams: list[str], end: str = "") -> None:
    stringList = []
    keys = list(singleVarFuncs.keys())
    for multiVarFnam in multiVarFnams:
        for singleVarFnam in keys:
            stringList.append(f"M.{multiVarFnam + end}.{singleVarFnam} = {{\n")

    #Generate permutations of values to loop through
    numVar = len(multiVarExprs[0].atoms(Symbol))
    expansionPoints = [*expansionPoints] * numVar
    expansionPoints = list(set(permutations(expansionPoints, numVar)))
    
    #Do each multivariate function as a different process
    with Manager() as manager:
        return_dict = manager.dict()
        jobs = []
        for i in range(len(multiVarExprs)):
            p = Process(target=generateStringsPerMultiVarFunc, args=(i, multiVarExprs, singleVarFuncs, keys, expansionPoints, order, epsList, return_dict))
            jobs.append(p)
            p.start()

        for proc in jobs: #Wait for processes to end
            proc.join() 

        for multiVarIndex in range(len(multiVarExprs)):
            for keyIndex in range(len(keys)): 
                    stringList[keyIndex + len(keys) * multiVarIndex] += return_dict.values()[multiVarIndex][keyIndex]
        for proc in jobs: #close the processes
            proc.close()
    
    with open(r"../" + filename, "a") as file:
        for i in range(len(stringList)):
            stringList[i] += "}\n\n"
            file.write(stringList[i])