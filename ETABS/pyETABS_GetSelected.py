import pyETABS_Attach

# =================================================================
# Get selected object from ETABS model [cSelect.GetSelected Method]
# AS PER ETABS API DOCUMENTATION
# ObjType:
# 1. Point Object
# 2. Frame Object
# 3. Cable Object
# 4. Tendon Object
# 5. Area Object
# 6. Solid Object
# 7. Link Object
# =================================================================


# Name = "A"
# NumberItems = 0
ObjType = 2
ObjName = []

ret = pyETABS_Attach.SapModel.SelectObj.GetSelected(NumberItems, ObjType, ObjName)

# ret = pyETABS_Attach.SapModel.Results.Setup.SetCaseSelectedForOutput("AREA")

# ret = pyETABS_Attach.SapModel.Results.FrameForce(Name, eItemTypeElm, NumberResults, Obj, ObjSta, Elm, ElmSta, LoadCase,StepType,StepNum, P,V2,V3,T,M2,M3)


print(ret)

