import pyETABS_Attach
import pd as pandas


def getDrift():
    NumberResults = 0
    story = []
    StepType = []
    StepNum = []
    Direction = []
    Drift = []
    Label = []
    LoadCase = []
    X = []
    Y = []
    Z = []

    ret = pyETABS_Attach.SapModel.Results.Setup.SetCaseSelectedForOutput("SPEC XY u1")

    ret = pyETABS_Attach.SapModel.results.storydrift(NumberResults, story, LoadCase, StepType, StepNum, Direction,
                                                     Drift, Label, X, Y, Z)

    return ret

    print(getDrift())


