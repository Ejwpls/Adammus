import pyETABS_Attach

def getfrGUID(label):
    # label = input("Column Label?")
    GUID = " "
    ret = pyETABS_Attach.SapModel.Frameobj.GetGUID(label, GUID)
    # print(ret[0])
    return ret[0]

def getfrForce(label):
    NumberResults = 0
    eItemTypeElm = 0

    Obj = []
    ObjSta = []
    Elm = []

    ElmSta = []
    LoadCase = []
    StepType = []
    StepNum = []

    P = []
    V2 = []
    V3 = []
    T = []
    M2 = []
    M3 = []



    ret = pyETABS_Attach.SapModel.Results.FrameForce(label, eItemTypeElm, NumberResults, Obj, ObjSta, Elm, ElmSta, LoadCase, StepType,StepNum, P,V2,V3,T,M2,M3)

    return ret

print((getfrForce("C40")))