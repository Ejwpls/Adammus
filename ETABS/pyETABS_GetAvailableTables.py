import pyETABS_Attach

# =================================================================

# =================================================================


# Name = "A"
NumberTables = 0
Tablekey = []
TableName = []
ImportType = []

#ret = pyETABS_Attach.SapModel.databasetables.GetAvailableTables(NumberTables,Tablekey, TableName,ImportType)
#print(ret)

# ret = pyETABS_Attach.SapModel.Results.Setup.SetCaseSelectedForOutput("AREA")

# ret = pyETABS_Attach.SapModel.Results.FrameForce(Name, eItemTypeElm, NumberResults, Obj, ObjSta, Elm, ElmSta, LoadCase,StepType,StepNum, P,V2,V3,T,M2,M3)

NumberSelecatedLoadCombinations = 0
LoadCombinationList = []

ret = pyETABS_Attach.SapModel.databasetables.GetLoadCombinationsSelectedForDisplay(NumberSelecatedLoadCombinations, LoadCombinationList)
print(ret)
